/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
/**
 * @file MLVariant.cpp
 * Implements class MLVariant.
 */

#include <algorithm>

#include "Lib/Backtrackable.hpp"
#include "Lib/BacktrackIterators.hpp"
#include "Lib/BinaryHeap.hpp"
#include "Lib/DArray.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Hash.hpp"
#include "Lib/Int.hpp"
#include "Lib/PairUtils.hpp"
#include "Lib/Stack.hpp"
#include "Lib/TriangularArray.hpp"

#include "Clause.hpp"
#include "Term.hpp"
#include "Matcher.hpp"
#include "MLVariant.hpp"

#if VDEBUG
#include <iostream>
using namespace std;
#endif

#define TRACE_LONG_MATCHING 0
#if TRACE_LONG_MATCHING
#include "Lib/Timer.hpp"
#endif

namespace Kernel
{

using namespace Lib;

typedef DHMap<unsigned,unsigned, IdentityHash, DefaultHash> UUMap;

namespace MLVariant_AUX
{

/**
 * Binder that stores bindings into a specified array. To be used
 * with MatchingUtils template methods.
 */
struct ArrayStoringBinder
{
  ArrayStoringBinder(unsigned* arr, UUMap& v2pos)
  : _arr(arr), _v2pos(v2pos) {}

  bool bind(unsigned var, TermList term)
  {
    ASS(term.isOrdinaryVar());
    _arr[_v2pos.get(var)]=term.var();
    return true;
  }

  void specVar(unsigned var, TermList term)
  { ASSERTION_VIOLATION; }
private:
  unsigned* _arr;
  UUMap& _v2pos;
};

bool createLiteralBindings(Literal* baseLit, LiteralList* alts, Clause* instCl,
    unsigned*& boundVarData, unsigned**& altBindingPtrs, unsigned*& altBindingData)
{
  static UUMap variablePositions;
  static BinaryHeap<unsigned,Int> varNums;
  variablePositions.reset();
  varNums.reset();

  VariableIterator bvit(baseLit);
  while(bvit.hasNext()) {
    unsigned var=bvit.next().var();
    varNums.insert(var);
  }

  unsigned nextPos=0;
  while(!varNums.isEmpty()) {
    unsigned var=varNums.pop();
    while(!varNums.isEmpty() && varNums.top()==var) {
      varNums.pop();
    }
    if(variablePositions.insert(var,nextPos)) {
      *(boundVarData++)=var;
      nextPos++;
    }
  }
  unsigned numVars=nextPos;

  LiteralList::Iterator ait(alts);
  while(ait.hasNext()) {
    //handle multiple matches in equality!
    Literal* alit=ait.next();
    if(alit->isEquality()) {
      //we must try both possibilities
      if(MatchingUtils::haveVariantArgs(baseLit,alit)) {
        ArrayStoringBinder binder(altBindingData, variablePositions);
        MatchingUtils::matchArgs(baseLit,alit,binder);
        *altBindingPtrs=altBindingData;
        altBindingPtrs++;
        altBindingData+=numVars;

        //add pointer to the literal at
        //the end of the binding sequence
        *(altBindingData++)=instCl->getLiteralPosition(alit);
      }
      if(MatchingUtils::haveReversedVariantArgs(baseLit, alit)) {
        ArrayStoringBinder binder(altBindingData, variablePositions);
        MatchingUtils::matchTerms(*baseLit->nthArgument(0),*alit->nthArgument(1),binder);
        MatchingUtils::matchTerms(*baseLit->nthArgument(1),*alit->nthArgument(0),binder);
        if(baseLit->isTwoVarEquality()){
          ASS(alit->isTwoVarEquality());
          MatchingUtils::matchTerms(baseLit->twoVarEqSort(),alit->twoVarEqSort(),binder);
        } //matchArgs automatically matches the sorts of literals if one is a twoVarEq literal
          //This is the reason for the difference between the two cases.
        *altBindingPtrs=altBindingData;
        altBindingPtrs++;
        altBindingData+=numVars;

        //add pointer to the literal at
        //the end of the binding sequence
        *(altBindingData++)=instCl->getLiteralPosition(alit);
      }

    } else {
      if(numVars) {
        ArrayStoringBinder binder(altBindingData, variablePositions);
        ALWAYS(MatchingUtils::matchArgs(baseLit,alit,binder));
      }

      *altBindingPtrs=altBindingData;
      altBindingPtrs++;
      altBindingData+=numVars;

      //add index of the literal in instance clause at
      //the end of the binding sequence
      *(altBindingData++) = (size_t)instCl->getLiteralPosition(alit);
    }
  }
  return true;
}

struct MatchingData {
  unsigned len;
  unsigned* varCnts;
  unsigned** boundVarNums;
  unsigned*** altBindings;
  TriangularArray<unsigned>* remaining;
  unsigned* nextAlts;

  Literal** bases;
  LiteralList** alts;
  Clause* instance;

  unsigned* boundVarNumStorage;
  unsigned** altBindingPtrStorage;
  unsigned* altBindingStorage;

  enum InitResult {
    OK,
    MUST_BACKTRACK,
    NO_ALTERNATIVE
  };

  unsigned getRemainingInCurrent(unsigned bi)
  {
    return remaining->get(bi,bi);
  }
  unsigned getAltRecordIndex(unsigned bi, unsigned alti)
  {
    return altBindings[bi][alti][varCnts[bi]];
  }

  /**
   * Return true if binding @b b1Index -th base literal that binds variables
   * to terms stored in @b i1Bindings is compatible to binding @b b2Index -th
   * base literal that binds variables to terms stored in
   * @b altBindings[b2Index][i2AltIndex] .
   */
  bool compatible(unsigned b1Index, unsigned* i1Bindings,
  	unsigned b2Index, unsigned i2AltIndex) const
  {
    if(varCnts[b1Index]==0 || varCnts[b2Index]==0) { return true; }

    //we'll create inverse binding to check that two variables
    //don't have the same binding
    static UUMap inverse;
    inverse.reset();

    //check that no variable is bound to two different ones
    unsigned* i2Bindings=altBindings[b2Index][i2AltIndex];
    unsigned* b1vn=boundVarNums[b1Index];
    unsigned* b2vn=boundVarNums[b2Index];
    unsigned* b1vnStop=boundVarNums[b1Index]+varCnts[b1Index];
    unsigned* b2vnStop=boundVarNums[b2Index]+varCnts[b2Index];
    while(true) {
      while(b1vn!=b1vnStop && *b1vn<*b2vn) {
        if(inverse.findOrInsert(*i1Bindings, *b1vn)!=*b1vn) { return false; }
        b1vn++; i1Bindings++;
      }
      if(b1vn==b1vnStop) { break; }
      while(b2vn!=b2vnStop && *b1vn>*b2vn) {
        if(inverse.findOrInsert(*i2Bindings, *b2vn)!=*b2vn) { return false; }
        b2vn++; i2Bindings++;
      }
      if(b2vn==b2vnStop) { break; }
      if(*b1vn==*b2vn) {
        if(*i1Bindings!=*i2Bindings) { 
          return false; 
        }
        if(inverse.findOrInsert(*i1Bindings, *b1vn)!=*b1vn) { return false; }
        b1vn++; i1Bindings++;
        b2vn++; i2Bindings++;
        if(b1vn==b1vnStop || b2vn==b2vnStop) { break; }
      }
    }
    while(b1vn!=b1vnStop) {
      if(inverse.findOrInsert(*i1Bindings, *b1vn)!=*b1vn) { return false; }
      b1vn++; i1Bindings++;
    }
    while(b2vn!=b2vnStop) {
      if(inverse.findOrInsert(*i2Bindings, *b2vn)!=*b2vn) { return false; }
      b2vn++; i2Bindings++;
    }
    return true;
  }

  bool bindAlt(unsigned bIndex, unsigned altIndex)
  {
    unsigned* curBindings=altBindings[bIndex][altIndex];
    for(unsigned i=bIndex+1; i<len; i++) {
      if(!isInitialized(i)) {
        break;
      }
      unsigned remAlts=remaining->get(i,bIndex);
      for(unsigned ai=0;ai<remAlts;ai++) {
        if(!compatible(bIndex,curBindings,i,ai)) {
          remAlts--;
          std::swap(altBindings[i][ai], altBindings[i][remAlts]);
          ai--;
        }
      }
      if(remAlts==0) {
        return false;
      }
      remaining->set(i,bIndex+1,remAlts);
    }
    return true;
  }

  bool isInitialized(unsigned bIndex) {
    return boundVarNums[bIndex];
  }

  InitResult ensureInit(unsigned bIndex)
  {
    if(!isInitialized(bIndex)) {
      boundVarNums[bIndex]=boundVarNumStorage;
      altBindings[bIndex]=altBindingPtrStorage;
      ALWAYS(createLiteralBindings(bases[bIndex], alts[bIndex], instance,
	  boundVarNumStorage, altBindingPtrStorage, altBindingStorage));
      varCnts[bIndex]=boundVarNumStorage-boundVarNums[bIndex];
      unsigned altCnt=altBindingPtrStorage-altBindings[bIndex];
      if(altCnt==0) {
        return NO_ALTERNATIVE;
      }
      remaining->set(bIndex, 0, altCnt);

      unsigned remAlts=0;
      for(unsigned pbi=0;pbi<bIndex;pbi++) { //pbi ~ previous base index
        unsigned* pbBindings=altBindings[pbi][nextAlts[pbi]-1];
        remAlts=remaining->get(bIndex, pbi);
        for(unsigned ai=0;ai<remAlts;ai++) {
          if(!compatible(pbi,pbBindings,bIndex,ai)) {
            remAlts--;
            std::swap(altBindings[bIndex][ai], altBindings[bIndex][remAlts]);
            ai--;
          }
        }
        remaining->set(bIndex,pbi+1,remAlts);
      }
      if(bIndex>0 && remAlts==0) {
        return MUST_BACKTRACK;
      }
    }
    return OK;
  }

};

MatchingData* getMatchingData(Literal* const * baseLits0, unsigned baseLen, Clause* instance, LiteralList** alts)
{
  static DArray<Literal*> baseLits(32);
  static DArray<LiteralList*> altsArr(32);
  baseLits.initFromArray(baseLen,baseLits0);
  altsArr.initFromArray(baseLen,alts);

  static DArray<unsigned> varCnts(32);
  static DArray<unsigned*> boundVarNums(32);
  static DArray<unsigned**> altPtrs(32);
  static TriangularArray<unsigned> remaining(32);
  static DArray<unsigned> nextAlts(32);

  static DArray<unsigned> boundVarNumData(64);
  //Referenced list of an alternative contains binding for each variable of
  //the base literal, and then a record identifying the alternative literal itself.
  static DArray<unsigned*> altBindingPtrs(128);
  static DArray<unsigned> altBindingsData(256);
  varCnts.ensure(baseLen);
  boundVarNums.init(baseLen,0);
  altPtrs.ensure(baseLen);
  remaining.setSide(baseLen);
  nextAlts.ensure(baseLen);

  unsigned zeroAlts=0;
  unsigned singleAlts=0;
  size_t baseLitVars=0;
  size_t altCnt=0;
  size_t altBindingsCnt=0;
  for(unsigned i=0;i<baseLen;i++) {
//    unsigned distVars=(*base)[i]->distinctVars();
    unsigned distVars=baseLits[i]->numVarOccs(); //an upper estimate is enough
    baseLitVars+=distVars;
    unsigned currAltCnt=0;
    LiteralList::Iterator ait(altsArr[i]);
    while(ait.hasNext()) {
      currAltCnt++;
      if(ait.next()->isEquality()) {
	currAltCnt++;
      }
    }
    altCnt+=currAltCnt+2;
    altBindingsCnt+=(distVars+1)*(currAltCnt+2);

    if(currAltCnt==0) {
      if(zeroAlts!=i) {
	std::swap(baseLits[i],baseLits[zeroAlts]);
	std::swap(altsArr[i],altsArr[zeroAlts]);
      }
      zeroAlts++;
      singleAlts++;
    } else if(currAltCnt==1) {
      if(singleAlts!=i) {
	std::swap(baseLits[i],baseLits[singleAlts]);
	std::swap(altsArr[i],altsArr[singleAlts]);
      }
      singleAlts++;
    }
  }
  boundVarNumData.ensure(baseLitVars);
  altBindingPtrs.ensure(altCnt);
  altBindingsData.ensure(altBindingsCnt);

  static MatchingData res;
  res.len=baseLen;
  res.varCnts=varCnts.array();
  res.boundVarNums=boundVarNums.array();
  res.altBindings=altPtrs.array();
  res.remaining=&remaining;
  res.nextAlts=nextAlts.array();

  res.bases=baseLits.array();
  res.alts=altsArr.array();
  res.instance=instance;

  res.boundVarNumStorage=boundVarNumData.array();
  res.altBindingPtrStorage=altBindingPtrs.array();
  res.altBindingStorage=altBindingsData.array();

  return &res;
}

} // MLVariant_AUX

using namespace MLVariant_AUX;

/**
 *
 */
bool MLVariant::isVariant(Literal* const * cl1Lits, Clause* cl2, LiteralList** alts, SimpleSubstitution *map)
{
  unsigned clen=cl2->length();

  MatchingData* md=getMatchingData(cl1Lits, clen, cl2, alts);
  if(!md) {
    return false;
  }
  static DArray<unsigned> matchRecord(32);
  unsigned matchRecordLen=clen;
  matchRecord.init(matchRecordLen,0xFFFFFFFF);


  unsigned matchedLen=md->len;

  md->nextAlts[0]=0;
  unsigned currBLit=0;

  while(true) {
    MatchingData::InitResult ires=md->ensureInit(currBLit);
    if(ires!=MatchingData::OK) {
      if(ires==MatchingData::MUST_BACKTRACK) {
        currBLit--;
        continue;
      } else {
        ASS_EQ(ires,MatchingData::NO_ALTERNATIVE);
        return false;
      }
    }

    unsigned maxAlt=md->getRemainingInCurrent(currBLit);
    while(md->nextAlts[currBLit]<maxAlt &&
	    ( matchRecord[md->getAltRecordIndex(currBLit, md->nextAlts[currBLit])]<currBLit ||
	    !md->bindAlt(currBLit,md->nextAlts[currBLit]) ) ) {
      md->nextAlts[currBLit]++;
    }
    if(md->nextAlts[currBLit] < maxAlt) {
      unsigned matchRecordIndex=md->getAltRecordIndex(currBLit, md->nextAlts[currBLit]);
      for(unsigned i=0;i<matchRecordLen;i++) {
        if(matchRecord[i]==currBLit) {
          matchRecord[i]=0xFFFFFFFF;
        }
      }
      if(matchRecord[matchRecordIndex]>currBLit) {
        matchRecord[matchRecordIndex]=currBLit;
      }
      md->nextAlts[currBLit]++;
      currBLit++;
      if(currBLit==matchedLen) { break; }
      md->nextAlts[currBLit]=0;
    } else {
      if(currBLit==0) { return false; }
      currBLit--;
    }
  }

  if(map)
    for(unsigned i = 0; i < cl2->length(); i++)
      ALWAYS(MatchingUtils::match(alts[i]->head(), cl1Lits[i], false, *map))

  return true;
}


bool MLVariant::isVariant(Literal* const * cl1Lits, Clause* cl2, bool complementary, SimpleSubstitution *map)
{
  bool fail=false;
  unsigned clen=cl2->length();
  static DArray<LiteralList*> alts(32);
  alts.init(clen, 0);

  for(unsigned i2=0;i2<clen;i2++) {
    bool cl2LitHasVariant=false;
    for(unsigned i1=0;i1<clen;i1++) {
      if(MatchingUtils::isVariant(cl1Lits[i1], (*cl2)[i2], complementary)) {
      cl2LitHasVariant=true;
      LiteralList::push((*cl2)[i2], alts[i1]);
      }
    }
    if(!cl2LitHasVariant) {
      fail=true;
      goto fin;
    }
  }
  for(unsigned i=0;i<clen;i++) {
    if(!alts[i]) {
      fail=true;
      goto fin;
    }
  }

  fail=!MLVariant::isVariant(cl1Lits,cl2,alts.array(),map);

fin:
  for(unsigned i=0;i<clen;i++) {
    LiteralList::destroy(alts[i]);
  }

  return !fail;
}

}
