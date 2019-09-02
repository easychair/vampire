
/*
 * File CombinatorDemodISE.cpp.
 *
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 *
 * In summary, you are allowed to use Vampire for non-commercial
 * purposes but not allowed to distribute, modify, copy, create derivatives,
 * or use in competitions. 
 * For other uses of Vampire please contact developers for a different
 * licence, which we will make an effort to provide. 
 */
/**
 * @file CombinatorDemodISE.cpp
 * Implements class CombinatorDemodISE.
 */

#include "Lib/Random.hpp"
#include "Lib/Environment.hpp"
#include "Lib/DArray.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/ApplicativeHelper.hpp"
#include "Kernel/TermIterators.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/SKIKBO.hpp"
#include "Kernel/SortHelper.hpp"
#include "Shell/Statistics.hpp"
#include "CombinatorDemodISE.hpp"

using namespace Lib;
using namespace Kernel;
using namespace Inferences;

typedef ApplicativeHelper AH; 

struct CombinatorDemodISE::HigherOrderTermInfo
{
public:
  HigherOrderTermInfo(TermList h, TermList hs, unsigned an){
    head = h;
    headSort = hs;
    argNum = an;
  }

  TermList head;
  TermList headSort;
  unsigned argNum;
};


Clause* CombinatorDemodISE::simplify(Clause* c)
{
  CALL("CombinatorDemodISE::simplify");

  Literal* newLit;
  LiteralStack litStack;
  bool modified = false;

  for(unsigned i = 0; i < c->length(); i++){
    Literal* lit = (*c)[i];
    ASS(lit->isEquality());
    TermList t0 = *lit->nthArgument(0);
    TermList t1 = *lit->nthArgument(1);
    
    TermList t0r = t0.isVar() ? t0 : reduce(t0);
    TermList t1r = t1.isVar() ? t1 : reduce(t1);      
    
    if((t0r != t0) || (t1r != t1)){
      modified = true;
      newLit = Literal::createEquality(lit->polarity(), TermList(t0r), TermList(t1r), SortHelper::getResultSort(t0.term()));
      litStack.push(newLit);
    } else {
      litStack.push(lit);
    }  
  }

  if(!modified){
    return c;
  }

  Inference* inference = new Inference1(Inference::COMBINATOR_DEMOD, c);
  Clause* newC = Clause::fromStack(litStack, c->inputType(), inference);
  //cout << "into CombinatorDemodISE " + c->toString() << endl;
  //cout << "out of CombinatorDemodISE " + newC->toString() << endl;
  if(!newC){ cout << "RETURNING NULL CLAUSE" << endl; }
  return newC;
}

TermList CombinatorDemodISE::reduce(TermList t)
{
  CALL("CombinatorDemodISE::reduce");
  
  ASS(!t.isVar());
    
  static Stack<Term*> terms(8);
  static Stack<HigherOrderTermInfo> infos(8);
  static Stack<bool> modified(8);
  static Stack<ApplicativeArgsIt*> argIts(8);
  static TermStack args;

  ASS(argIts.isEmpty());
  ASS(terms.isEmpty());
  ASS(infos.isEmpty());
  modified.reset();
  args.reset();

  headNormalForm(t);
  modified.push(false);
  argIts.push(new ApplicativeArgsIt(t, false));
  ApplicativeArgsIt* argsIt = argIts.top();
  infos.push(HigherOrderTermInfo(argsIt->head(), argsIt->headSort(), argsIt->argNum()));

  for (;;) {
    if (!argIts.top()->hasNext()) {
      delete argIts.pop();
      if (terms.isEmpty()) {
        //we're done, args stack contains modified arguments
        //of the literal.
        ASS(argIts.isEmpty());
        break;
      }
      Term* orig = terms.pop();
      HigherOrderTermInfo hoti=infos.pop();
      if (!modified.pop()) {
        args.truncate(args.length() - hoti.argNum);
        args.push(TermList(orig));
        continue;
      }
      //here we assume, that stack is an array with
      //second topmost element as &top()-1, third at
      //&top()-2, etc...
      TermList* argLst=&args.top() - (hoti.argNum - 1);
      args.truncate(args.length() - hoti.argNum);

      TermList trm = AH::createAppTerm(hoti.headSort, hoti.head, argLst, hoti.argNum);
      args.push(trm);
      modified.setTop(true);
      continue;
    }

    TermList tl= argIts.top()->next();
    bool reduced = headNormalForm(tl);
    if(reduced){
      modified.setTop(true);
    }
    if (tl.isVar()) {
      args.push(tl);
      continue;
    }
    ASS(tl.isTerm());
    Term* t=tl.term();
    terms.push(t);
    modified.push(false);
    argIts.push(new ApplicativeArgsIt(tl, false));
    argsIt = argIts.top();
    infos.push(HigherOrderTermInfo(argsIt->head(), argsIt->headSort(), argsIt->argNum()));
  }
  ASS(argIts.isEmpty());
  ASS(terms.isEmpty());
  ASS_EQ(modified.length(),1);
  ASS_EQ(infos.length(),1);
  HigherOrderTermInfo hoti=infos.pop();
  ASS_EQ(args.length(),hoti.argNum);

  if (!modified.pop()) {
    return t;
  }

  TermList* argLst=&args.top() - (hoti.argNum-1);
  ASS(!t.term()->isLiteral());
  return AH::createAppTerm(hoti.headSort, hoti.head, argLst, hoti.argNum);;
}

bool CombinatorDemodISE::headNormalForm(TermList& t)
{
  CALL("CombinatorDemodISE::headNormalForm");

  static TermStack args;
  TermList head;
  
  bool modified = false;
  
  for(;;){
    args.reset();
    AH::getHeadAndArgs(t, head, args);
    if(AH::isComb(head) && !AH::isUnderApplied(head, args.size())){
      modified = true;
      t = SKIKBO::reduce(args, head);
    } else {
      break; 
    }
  }
  return modified;
}



