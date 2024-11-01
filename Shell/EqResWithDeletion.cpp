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
 * @file EqResWithDeletion.cpp
 * Implements class EqResWithDeletion.
 */

#include "Kernel/Clause.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/Unit.hpp"
#include "Kernel/FormulaVarIterator.hpp"

#include "EqResWithDeletion.hpp"

namespace Shell
{

using namespace Lib;
using namespace Kernel;

void EqResWithDeletion::apply(Problem& prb)
{
  if(apply(prb.units())) {
    prb.invalidateByRemoval();
  }
}

/**
 * Perform equality resolution with deletion and return
 * true iff some clause was modified.
 */
bool EqResWithDeletion::apply(UnitList*& units)
{
  bool modified = false;

  UnitList::DelIterator uit(units);
  while(uit.hasNext()) {
    Clause* cl=static_cast<Clause*>(uit.next());
    ASS(cl->isClause());
    Clause* cl2=apply(cl);
    if(cl!=cl2) {
      modified = true;
      uit.replace(cl2);
    }
  }
  return modified;
}

/**
 * @warning The application of this rule can currently be quadratic.
 *
 * The reason this is so is that "t1.containsSubterm(t0)" and "t0.containsSubterm(t1)" below
 * don't suffice when considering simultaneous substitution. E.g. X != f(Y) | Y = g(X) | rest ...
 */
Clause* EqResWithDeletion::apply(Clause* cl)
{
  //TODO: make the procedure linear time
start_applying:

  unsigned clen=cl->length();
  if (env.options->questionAnswering() == Options::QuestionAnsweringMode::SYNTHESIS) {
    _ansLit = cl->getAnswerLiteral();
  }

  _subst.reset();

  RStack<Literal*> resLits;

  bool foundResolvable=false;
  std::unordered_set<Literal *> resolved;
  for(unsigned i=0;i<clen;i++) {
    Literal* lit=(*cl)[i];
    if(!foundResolvable && scan(lit)) {
      foundResolvable=true;
      if(env.options->proofExtra() == Options::ProofExtra::FULL)
        resolved.insert(lit);
    } else {
      resLits->push(lit);
    }
  }
  if(!foundResolvable) {
    return cl;
  }

  for(unsigned i=0;i<resLits->size();i++) {
    (*resLits)[i] = SubstHelper::apply((*resLits)[i], *this);
  }

  cl = Clause::fromStack(*resLits,
      SimplifyingInference1(InferenceRule::EQUALITY_RESOLUTION_WITH_DELETION, cl));
  if(env.options->proofExtra() == Options::ProofExtra::FULL)
    env.proofExtra.insert(cl, new EqResWithDeletionExtra(std::move(resolved)));
  goto start_applying;
}

TermList EqResWithDeletion::apply(unsigned var)
{
  TermList res;
  if(_subst.find(var, res)) {
    return res;
  } else {
    return TermList(var, false);
  }
}

bool EqResWithDeletion::scan(Literal* lit)
{
  using Kernel::isFreeVariableOf;

  if(lit->isEquality() && lit->isNegative()) {
    TermList t0=*lit->nthArgument(0);
    TermList t1=*lit->nthArgument(1);
    if( t0.isVar() && !t1.containsSubterm(t0) && (!_ansLit || !t1.isTerm() || t1.term()->computableOrVar() || !isFreeVariableOf(_ansLit,t0.var()))) {
      if(_subst.insert(t0.var(), t1)) {
        return true;
      }
    }
    if( t1.isVar() && !t0.containsSubterm(t1) && (!_ansLit || !t0.isTerm() || t0.term()->computableOrVar() || !isFreeVariableOf(_ansLit,t1.var()))) {
      if(_subst.insert(t1.var(), t0)) {
        return true;
      }
    }
  }
  return false;
}

void EqResWithDeletionExtra::output(std::ostream &out) const {
  bool first = true;
  out << "resolved=[";
  for(Literal *l : resolved) {
    if(!first)
      out << ",";
    first = false;
    out << "(" << l->toString() << ")";
  }
  out << "]";

}

}
