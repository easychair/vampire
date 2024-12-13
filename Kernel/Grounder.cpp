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
 * @file Grounder.cpp
 * Implements class Grounder.
 */

#include "Kernel/Clause.hpp"
#include "Kernel/Renaming.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/TermIterators.hpp"

#include "SAT/SATLiteral.hpp"
#include "SAT/SATSolver.hpp"

#include "Grounder.hpp"

using namespace Kernel;

void GlobalSubsumptionGrounder::groundNonProp(Clause* cl, SATLiteralStack& acc)
{
  static DArray<Literal*> lits;

  unsigned clen = cl->length();

  lits.initFromArray(clen, *cl);
  Literal **normLits = lits.array();

  normalize(clen, normLits);

  for(unsigned i=0; i<clen; i++) {
    SATLiteral lit = groundNormalized(normLits[i]);
    acc.push(lit);
  }
}

SATLiteral GlobalSubsumptionGrounder::groundLiteral(Literal* lit)
{
  Literal* norm = lit;
  normalize(1, &norm);
  SATLiteral slit = groundNormalized(norm);
  return slit;
}

SATLiteral GlobalSubsumptionGrounder::groundNormalized(Literal* lit)
{
  bool isPos = lit->isPositive();
  Literal* posLit = Literal::positiveLiteral(lit);

  unsigned* pvar;
  if(_asgn.getValuePtr(posLit, pvar)) {
    *pvar = _satSolver.newVar();
  }
  return SATLiteral(*pvar, isPos);
}

struct GlobalSubsumptionGrounder::OrderNormalizingComparator
{
  Literal** _lits;
  OrderNormalizingComparator(Literal** lits) : _lits(lits) {}

  bool operator()(unsigned a, unsigned b) {
    Literal* la = _lits[a];
    Literal* lb = _lits[b];
    if(la==lb) { return false; }
    if(la->numVarOccs()!=lb->numVarOccs()) {
      //first, we want literals with less variables to appear in the
      //beginning as there is better chance to get some sharing across clauses
      return la->numVarOccs()<lb->numVarOccs();
    }
    if(la->weight()!=lb->weight()) {
      return la->weight()<lb->weight();
    }
    if(la->header()!=lb->header()) {
      return la->header()<lb->header();
    }
    //now get just some total deterministic order
    static DisagreementSetIterator dsit;
    dsit.reset(la, lb, false);
    ALWAYS(dsit.hasNext());
    auto da = dsit.next();
    if(da.first.isVar()!=da.second.isVar()) {
      return da.first.isVar();
    }
    if(da.first.isVar()) {
      ASS_NEQ(da.first.var(),da.second.var());
      return da.first.var()<da.second.var();
    }
    else {
      ASS_NEQ(da.first.term()->functor(),da.second.term()->functor());
      return da.first.term()->functor()<da.second.term()->functor();
    }
    return a<b; //if nothing else applies, we keep the original order
  }
};

void GlobalSubsumptionGrounder::normalize(unsigned cnt, Literal** lits)
{
  if(cnt==0) { return; }
  if(cnt==1) {
    lits[0] = Renaming::normalize(lits[0]);
  }

  static Stack<unsigned> litOrder;
  litOrder.reset();
  litOrder.loadFromIterator(getRangeIterator(0u,cnt));

  std::sort(litOrder.begin(), litOrder.end(), OrderNormalizingComparator(lits));

  static Renaming normalizer;
  normalizer.reset();

  for(unsigned i=0; i<cnt; i++) {
    unsigned li = litOrder[i];
    normalizer.normalizeVariables(lits[li]);
    lits[li] = normalizer.apply(lits[li]);
  }
}
