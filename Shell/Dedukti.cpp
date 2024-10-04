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
 * @file Dedukti.cpp
 * Routines for Dedukti output
 */

#include "Dedukti.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/SortHelper.hpp"
#include "Parse/TPTP.hpp"
#include "Shell/UIHelper.hpp"

#include "Inferences/BinaryResolution.hpp"
using namespace Inferences;

const char *PRELUDE = R"((; Prop ;)
Prop : Type.
def Prf : (Prop -> Type).
true : Prop.
[] Prf true --> (r : Prop -> ((Prf r) -> (Prf r))).
false : Prop.
[] Prf false --> (r : Prop -> (Prf r)).
not : (Prop -> Prop).
[p] Prf (not p) --> ((Prf p) -> (r : Prop -> (Prf r))).
and : (Prop -> (Prop -> Prop)).
[p, q] Prf (and p q) --> (r : Prop -> (((Prf p) -> ((Prf q) -> (Prf r))) -> (Prf r))).
or : (Prop -> (Prop -> Prop)).
[p, q] Prf (or p q) --> (r : Prop -> (((Prf p) -> (Prf r)) -> (((Prf q) -> (Prf r)) -> (Prf r)))).
imp : (Prop -> (Prop -> Prop)).
[p, q] Prf (imp p q) --> ((Prf p) -> (Prf q)).

(; Set ;)
Set : Type.
injective El : (Set -> Type).
iota : Set.

(; Quant ;)
forall : (a : Set -> (((El a) -> Prop) -> Prop)).
[a, p] Prf (forall a p) --> (x : (El a) -> (Prf (p x))).
exists : (a : Set -> (((El a) -> Prop) -> Prop)).
[a, p] Prf (exists a p) --> (r : Prop -> ((x : (El a) -> ((Prf (p x)) -> (Prf r))) -> (Prf r))).

(; Classic ;)
def cPrf : (Prop -> Type) := (p : Prop => (Prf (not (not p)))).
def cand : (Prop -> (Prop -> Prop)) := (p : Prop => (q : Prop => (and (not (not p)) (not (not q))))).
def cor : (Prop -> (Prop -> Prop)) := (p : Prop => (q : Prop => (or (not (not p)) (not (not q))))).
def cimp : (Prop -> (Prop -> Prop)) := (p : Prop => (q : Prop => (imp (not (not p)) (not (not q))))).
def cforall : (a : Set -> (((El a) -> Prop) -> Prop)) := (a : Set => (p : ((El a) -> Prop) => (forall a (x : (El a) => (not (not (p x))))))).
def cexists : (a : Set -> (((El a) -> Prop) -> Prop)) := (a : Set => (p : ((El a) -> Prop) => (exists a (x : (El a) => (not (not (p x))))))).

(; Clauses ;)
def prop_clause : Type.
def ec : prop_clause.
def cons : (Prop -> (prop_clause -> prop_clause)).
def clause : Type.
def cl : (prop_clause -> clause).
def bind : (A : Set -> (((El A) -> clause) -> clause)).
def Prf_prop_clause : (prop_clause -> Type).

[] Prf_prop_clause ec --> (Prf false).
[p, c] Prf_prop_clause (cons p c) --> ((Prf p -> Prf false) -> (Prf_prop_clause c)).
def Prf_clause : (clause -> Type).
[c] Prf_clause (cl c) --> (Prf_prop_clause c).
[A, f] Prf_clause (bind A f) --> (x : (El A) -> (Prf_clause (f x))).

)";

using namespace Kernel;

namespace Shell {
namespace Dedukti {

void outputPrelude(std::ostream &out) {
  out << PRELUDE;
}

void outputTypeDecl(std::ostream &out, const char *name, OperatorType *type) {
  out << name << ": ";

  // we don't support polymorphism yet
  ASS_EQ(type->numTypeArguments(), 0)
  // we don't support many-sorted logic yet
  ASS(type->isAllDefault())

  for(unsigned i = 0; i < type->arity(); i++)
    out << "(El iota) -> ";

  TermList range = type->result();
  // we don't support many-sorted logic yet
  ASS(range.isEmpty() || range == AtomicSort::defaultSort())

  // predicate
  if(range.isEmpty())
    out << "Prop";
  // function
  else
    out << "(El iota)";

  out << "." << std::endl;
}

static void outputTermList(std::ostream &out, TermList *start) {
  ASS(start->isNonEmpty())

  Stack<TermList *> todo;
  TermList *current = start;
  while(true) {
    out << " ";
    if(current->isVar()) {
      out << "_" << current->var();
      current = current->next();
    }
    else if(current->isTerm()) {
      Term *term = current->term();
      if(term->arity()) {
        out << "(" << term->functionName();
        todo.push(current->next());
        current = term->args();
      }
      else {
        out << term->functionName();
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

void outputClause(std::ostream &out, Clause *clause) {

  if(clause->isEmpty()) {
    out << "Prf_clause (cl ec)";
    return;
  }

  out << "Prf_clause ";
  DHMap<unsigned, TermList> sorts;
  SortHelper::collectVariableSorts(clause, sorts);

  auto it = sorts.items();
  while(it.hasNext()) {
    std::pair<unsigned, TermList> next = it.next();
    unsigned var = next.first;
    TermList sort = next.second;
    // we don't support many-sorted logic yet
    ASS(sort == AtomicSort::defaultSort())
    out << "(bind iota (_" << var << " : (El iota) => ";
  }

  out << "(cl ";
  for(unsigned i = 0; i < clause->length(); i++) {
    Literal *literal = (*clause)[i];

    out << "(cons ";
    if(!literal->polarity())
      out << "(not ";
    if(literal->arity())
      out << "(";

    out << literal->predicateName();
    if(literal->arity())
      outputTermList(out, literal->args());

    if(literal->arity())
      out << ")";
    if(!literal->polarity())
      out << ")";
    out << " ";
  }
  out << "ec";

  for(unsigned i = 0; i < clause->length(); i++)
    out << ")";
  for(unsigned i = 0; i < sorts.size(); i++)
    out << "))";

  out << ")";
}

void outputAxiomName(std::ostream &out, Unit *axiom) {
  std::string recoveredName;
  out << "_axiom_";
  if(Parse::TPTP::findAxiomName(axiom, recoveredName))
    out << recoveredName;
  else
    out << axiom->number();
}

void outputAxiom(std::ostream &out, Unit *axiom) {
  // we don't support non-clause inputs yet
  ASS(axiom->isClause())
  // we don't support non-axiom inputs yet
  ASS(
    axiom->inputType() == UnitInputType::AXIOM ||
    axiom->inputType() == UnitInputType::ASSUMPTION ||
    axiom->inputType() == UnitInputType::NEGATED_CONJECTURE
  )

  outputAxiomName(out, axiom);
  out << ": ";
  Clause *clause = static_cast<Clause *>(axiom);
  outputClause(out, clause);
  out << "." << std::endl;
}

void outputDeductionPrefix(std::ostream &out, Unit *deduction) {
  // we don't support non-clause deductions yet
  ASS(deduction->isClause())

  out << "def _deduction_" << deduction->number() << ": ";
  Clause *clause = static_cast<Clause *>(deduction);
  outputClause(out, clause);
  out << " := ";
}

void outputSorry(std::ostream &out, Unit *admit) {
  // we don't support non-clause deductions yet
  ASS(admit->isClause())

  Clause *clause = static_cast<Clause *>(admit);
  out << "_sorry_" << admit->number() << ": ";

  UnitIterator parents = admit->getParents();
  while(parents.hasNext()) {
    Unit *parent = parents.next();
    ASS(parent->isClause());
    Clause *parent_clause = static_cast<Clause *>(parent);
    outputClause(out, parent_clause);
    out << " -> ";
  }
  outputClause(out, clause);
  out << "." << std::endl;

  outputDeductionPrefix(out, admit);
  out << "_sorry_" << admit->number();
  parents = admit->getParents();
  while(parents.hasNext())
    out << " _deduction_" << parents.next()->number();
  out << "." << std::endl;
}


void outputUnit(std::ostream &out, Unit *u) {
  ASS(u->isClause())
  InferenceRule rule = u->inference().rule();
  switch (rule)
  {
  case InferenceRule::INPUT:
  case InferenceRule::NEGATED_CONJECTURE:
    outputAxiom(out, u);
    outputDeductionPrefix(out, u);
    outputAxiomName(out, u);
    out << "." << std::endl;
    break;
  case InferenceRule::RESOLUTION: {
    const BinaryResolutionExtra &br = env.proofExtra.get<Inferences::BinaryResolutionExtra>(u);
  }
  default: outputSorry(out, u);
  }
}

void outputInput(std::ostream &out, Unit *input) {
  // TODO figure out what to do about literal selection
  outputAxiom(out, input);
  outputUnit(out, input);
  outputAxiomName(out, input);
  out << "." << std::endl;
}




// ProofPrinter::ProofPrinter(std::ostream &out, InferenceStore *store) : InferenceStore::ProofPrinter(out, store) {
//   outputPrelude(out);
//   UIHelper::outputSymbolDeclarations(out);
// }

// void ProofPrinter::printStep(Unit *unit) {
//   Inference &inference = unit->inference();
//   InferenceRule rule = inference.rule();

//   if(rule == InferenceRule::INPUT) {
//     outputInput(out, unit);
//     return;
//   }

//   UnitIterator parents = _is->getParents(unit);
//   switch(rule) {
//   // case InferenceRule::RESOLUTION:
//   // {
//   //   BinaryResolution *br = static_cast<BinaryResolution *>(datum);
//   //   Clause *left = static_cast<Clause *>(parents.next());
//   //   Clause *right = static_cast<Clause *>(parents.next());
//   //   out << "(; binary resolution ;)" << std::endl;
//   //   out << "(; " << left->toString() << " ;)" << std::endl;
//   //   out << "(; " << right->toString() << " ;)" << std::endl;
//   //   out << "(; ---------------------------- ;)" << std::endl;
//   //   out << "(; " << unit->toString() << " ;)" << std::endl;

//   //   TermList latom((*left)[br->leftIndex]);
//   //   TermList ratom(Literal::complementaryLiteral((*right)[br->rightIndex]));
//   //   RobSubstitution subst;
//   //   ALWAYS(subst.unify(latom, 0, ratom, 1));
//   //   out << "(; " << subst << " ;)" << std::endl;
//   //   outputOops(out, unit, _is);
//   //   break;
//   // }
//   default:
//     // should not register inferences for rules we don't handle
//     ASSERTION_VIOLATION;
//   }
// }

}
}
