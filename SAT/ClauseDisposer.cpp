/**
 * @file ClauseDisposer.cpp
 * Implements class ClauseDisposer.
 */

#include <algorithm>

#include "Lib/BinaryHeap.hpp"
#include "Lib/Int.hpp"

#include "TWLSolver.hpp"

#include "ClauseDisposer.hpp"

namespace SAT
{

using namespace std;

unsigned ClauseDisposer::varCnt() const
{
  CALL("ClauseDisposer::varCnt");

  return _solver._varCnt;
}

SATClauseStack& ClauseDisposer::getLearntStack()
{
  CALL("ClauseDisposer::getLearntStack");

  return _solver._learntClauses;
}

DArray<SATClauseStack>& ClauseDisposer::getWatchedStackArray()
{
  CALL("ClauseDisposer::getWatchedStackArray");

  return _solver._windex;
}

SATClause* ClauseDisposer::getAssignmentPremise(unsigned var)
{
  CALL("ClauseDisposer::getAssignmentPremise");

  return _solver._assignmentPremises[var];
}

/**
 * Clean the 'kept' flag of all learnt clauses that are not used
 * as a justification for some assignment.
 */
void ClauseDisposer::markAllRemovableUnkept()
{
  CALL("ClauseDisposer::markAllRemovableUnkept");

  SATClauseStack::Iterator lrnIt(getLearntStack());
  while(lrnIt.hasNext()) {
    SATClause* cl = lrnIt.next();
    cl->setKept(false);
  }

  unsigned vc = varCnt();
  for(unsigned i=0; i<vc; i++) {
    SATClause* cl = getAssignmentPremise(i);
    if(cl) {
      cl->setKept(true);
    }
  }
}

/**
 * We assume that all non-learnt clauses have their 'kept' flag
 * set to true.
 */
void ClauseDisposer::removeUnkept()
{
  CALL("ClauseDisposer::removeUnkept");

  unsigned watchCnt = varCnt()*2;
  DArray<SATClauseStack>& watches = getWatchedStackArray();

  for(unsigned i=0; i<watchCnt; i++) {
    SATClauseStack::Iterator wit(watches[i]);
    while(wit.hasNext()) {
      SATClause* cl = wit.next();
      if(!cl->kept()) {
	wit.del();
      }
    }
  }

  SATClauseStack::StableDelIterator lrnIt(getLearntStack());
  while(lrnIt.hasNext()) {
    SATClause* cl = lrnIt.next();
    if(!cl->kept()) {
      lrnIt.del();
      cl->destroy();
    }
  }
}

struct ClauseActivityComparator
{
  static Comparison compare(SATClause* c1, SATClause* c2)
  {
//    return Int::compare(c1->activity()*c2->length(), c2->activity()*c1->length());
    return Int::compare(c1->activity(), c2->activity());
  }
};

void ClauseDisposer::keepMostActive(size_t numberOfKept, ActivityType minActivity)
{
  CALL("ClauseDisposer::keepMostActive");

  static BinaryHeap<SATClause*, ClauseActivityComparator> mah; //most active heap
  mah.reset();

  SATClauseStack::Iterator lrnIt(getLearntStack());
  while(lrnIt.hasNext()) {
    SATClause* cl = lrnIt.next();
    if(cl->activity()<minActivity) {
      continue;
    }
    mah.insert(cl);
    if(mah.size()>numberOfKept) {
	mah.pop();
    }
  }
  while(!mah.isEmpty()) {
    mah.pop()->setKept(true);
  }
}

void ClauseDisposer::keepBinary()
{
  CALL("ClauseDisposer::keepBinary");

  size_t idx = 0;
  SATClauseStack::Iterator lrnIt(getLearntStack());
  while(lrnIt.hasNext()) {
    SATClause* cl = lrnIt.next();
    if(cl->size()<=2) {
      cl->setKept(true);
    }
    idx++;
  }
}

///////////////////////////
// DecayingClauseDisposer

/**
 * Decay the importance of former conflicts
 */
void DecayingClauseDisposer::onConflict()
{
  CALL("DecayingClauseDisposer::onConflict");

  _inc *= _decayFactor;
  if(_inc<1e30f) {
    return;
  }

  _inc /= 1e30f;
  SATClauseStack::StableDelIterator lrnIt(getLearntStack());
  while(lrnIt.hasNext()) {
    SATClause* cl = lrnIt.next();
    cl->activity() /= 1e30f;
  }
}


///////////////////////////
// MinisatClauseDisposer

void MinisatClauseDisposer::onRestart()
{
  CALL("MinisatClauseDisposer::onRestart");

  markAllRemovableUnkept();

  unsigned learntCnt = getLearntStack().size();

  keepMostActive(learntCnt/2, _inc/learntCnt);
  keepBinary();

  removeUnkept();
}


///////////////////////////
// GrowingClauseDisposer

void GrowingClauseDisposer::onNewInputClause(SATClause* cl)
{
  CALL("GrowingClauseDisposer::onNewInputClause");

  DecayingClauseDisposer::onNewInputClause(cl);
  if(cl->size()<3) {
    return;
  }
  _clauseCntAcc++;
  if(_clauseCntAcc>=4) {
    _survivorCnt += _clauseCntAcc/4;
    _clauseCntAcc += _clauseCntAcc%4;
  }
}

void GrowingClauseDisposer::onRestart()
{
  CALL("GrowingClauseDisposer::onRestart");

  markAllRemovableUnkept();

  keepMostActive(_survivorCnt, 0);
  keepBinary();

  removeUnkept();

  _phaseIdx++;
  if(_phaseIdx==_phaseLen) {
    _survivorCnt = _survivorCnt+max(_survivorCnt/20,static_cast<size_t>(1));
    _phaseIdx = 0;
    _phaseLen++;
  }
}


}
