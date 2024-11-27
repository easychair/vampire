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
 * @file SMTCheck.hpp
 * Routines for producing a ground SMT proof-check script
 */

#include <unordered_set>

#include "Inferences/BinaryResolution.hpp"
#include "SMTCheck.hpp"

#include "Kernel/RobSubstitution.hpp"
#include "Kernel/SortHelper.hpp"

using namespace Kernel;

// get N parents of a unit
// TODO merge with Dedukti version
template<unsigned N, typename T>
std::array<T *, N> getParents(T *unit) {
  std::array<T *, N> parents;
  UnitIterator it = unit->getParents();
  for(unsigned i = 0; i < N; i++) {
    ALWAYS(it.hasNext())
    parents[i] = static_cast<T *>(it.next());
  }
  ALWAYS(!it.hasNext())
  return parents;
}


static void outputName(std::ostream &out, const std::string &name) {
  out << '|' << name << '|';
}

static void outputVar(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, unsigned var) {
  if(conclSorts.findPtr(var))
    out << "v" << var;
  else
    out << "inhabit";
}

static void outputArgs(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, TermList *start) {
  ASS(start->isNonEmpty())

  Stack<TermList *> todo;
  TermList *current = start;
  while(true) {
    out << " ";
    if(current->isVar()) {
      outputVar(out, conclSorts, current->var());
      current = current->next();
    }
    else if(current->isTerm()) {
      Term *term = current->term();
      if(term->arity()) {
        out << "(";
        outputName(out, term->functionName());
        todo.push(current->next());
        current = term->args();
      }
      else {
        outputName(out, term->functionName());
        current = current->next();
      }
    }

    while(current->isEmpty()) {
      if(todo.isEmpty())
        return;

      current = todo.pop();
      out << ")";
    }
  }
}

static void outputTerm(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, TermList tl) {
  if(tl.isVar()) {
    outputVar(out, conclSorts, tl.var());
    return;
  }

  else if(tl.isTerm()) {
    Term *term = tl.term();
    if(term->arity())
      out << "(";
    outputName(out, term->functionName());
    if(term->arity())
      outputArgs(out, conclSorts, term->args());
    if(term->arity())
      out << ")";
  }
}

static void outputLiteral(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, Literal *literal) {
  if(!literal->polarity())
    out << "(not ";
  if(literal->arity())
    out << "(";
  outputName(out, literal->predicateName());
  if(literal->arity()) {
    outputArgs(out, conclSorts, literal->args());
    out << ")";
  }
  if(!literal->polarity())
    out << ")";
}

static void outputPremise(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, Clause *cl, const RobSubstitution &subst, int bank) {
  out << "(assert (or";
  for(Literal *l : cl->iterLits()) {
    out << " ";
    outputLiteral(out, conclSorts, subst.apply(l, bank));
  }
  out << "))\n";
}


static void outputConclusion(std::ostream &out, DHMap<unsigned, TermList> &conclSorts, Clause *cl) {
  for(Literal *l : cl->iterLits()) {
    out << "(assert ";
    outputLiteral(out, conclSorts, Literal::complementaryLiteral(l));
    out << ")\n";
  }
}

static DHMap<unsigned, TermList> outputSkolems(std::ostream &out, Unit *u) {
  DHMap<unsigned, TermList> sorts;
  SortHelper::collectVariableSorts(u, sorts);
  auto it = sorts.items();
  while(it.hasNext()) {
    auto [var, sort] = it.next();
    ASS_EQ(sort, AtomicSort::defaultSort())
    out << "(declare-const v" << var << " I)\n";
  }
  return sorts;
}


static void outputResolution(std::ostream &out, Clause *concl) {
  auto conclSorts = outputSkolems(out, concl);
  auto [left, right] = getParents<2>(concl);
  const auto &br = env.proofExtra.get<Inferences::BinaryResolutionExtra>(concl);

  // compute unifier for selected literals
  RobSubstitution subst;
  Literal *selectedLeft = br.selectedLiteral.selectedLiteral;
  Literal *selectedRight = br.otherLiteral;
  ASS_NEQ(selectedLeft->polarity(), selectedRight->polarity())
  ALWAYS(subst.unify(TermList(selectedLeft), 0, TermList(selectedRight), 1));

  for(unsigned i = 0; i < left->length(); i++)
    if((*left)[i] != selectedLeft)
      subst.apply((*left)[i], 0);
  for(unsigned i = 0; i < right->length(); i++)
    if((*right)[i] != selectedRight)
      subst.apply((*right)[i], 1);

  outputPremise(out, conclSorts, left->asClause(), subst, 0);
  outputPremise(out, conclSorts, right->asClause(), subst, 1);
  outputConclusion(out, conclSorts, concl->asClause());
}

namespace Shell {
namespace SMTCheck {

void outputTypeDecl(std::ostream &out, const std::string &name, OperatorType *type) {
  out << "(declare-fun ";
  outputName(out, name);

  TermList range = type->result();

  // we don't support polymorphism yet
  ASS_EQ(type->numTypeArguments(), 0)
  // we don't support many-sorted logic yet
  ASS(type->isAllDefault())
  // we don't support many-sorted logic yet
  ASS(range.isEmpty() || range == AtomicSort::defaultSort())

  out << " (";
  for(unsigned i = 0; i < type->arity(); i++)
    out << (i == 0 ? "" : " ") << "I";
  out << ") " << (range.isEmpty() ? "Bool" : "I") << ")\n";
}

void outputStep(std::ostream &out, Unit *u) {
  InferenceRule rule = u->inference().rule();
  if(rule == InferenceRule::INPUT || rule == InferenceRule::NEGATED_CONJECTURE)
    return;

  out << "\n; step " << u->number() << "\n";
  out << "(push)\n";
  switch(rule) {
  case InferenceRule::RESOLUTION:
    outputResolution(out, u->asClause());
    break;
  default:
    out << "(echo \"sorry: " << ruleName(rule) << "\")\n";
    out << "(assert false)\n";
  }
  out << "(set-info :status unsat)\n";
  out << "(check-sat)\n";
  out << "(pop)\n";
}

}
}
