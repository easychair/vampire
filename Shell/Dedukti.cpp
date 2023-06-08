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
#include "Kernel/SortHelper.hpp"

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

)";
const char *IOTA = "iota";
const char *ELEMENT = "El";
const char *PROPOSITION = "Prop";
const char *PROOF = "Prf";

using namespace Kernel;

namespace Shell {
namespace Dedukti {

void outputPrelude(std::ostream &out) {
  CALL("Dedukti::outputPrelude")
  out << PRELUDE;
}

void outputTypeDecl(std::ostream &out, const char *name, OperatorType *type) {
  CALL("Dedukti::outputTypeDecl")
  out << name << ": ";

  // we don't support polymorphism yet
  ASS_EQ(type->numTypeArguments(), 0)
  // we don't support many-sorted logic yet
  ASS(type->isAllDefault())

  for(unsigned i = 0; i < type->arity(); i++)
    out << "(" << ELEMENT << " " << IOTA << ")" << " -> ";

  TermList range = type->result();
  // we don't support many-sorted logic yet
  ASS(range.isEmpty() || range == AtomicSort::defaultSort())

  // predicate
  if(range.isEmpty())
    out << PROPOSITION;
  // function
  else
    out << "(" << ELEMENT << " " << IOTA << ")";

  out << "." << std::endl;
}

static void outputTermList(std::ostream &out, TermList *start) {
  CALL("Dedukti::outputTermList")
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

static void outputClause(std::ostream &out, Clause *clause) {
  CALL("Dedukti::outputClause")
  ASS(!clause->isEmpty())

  out << PROOF << " ";
  DHMap<unsigned, TermList> sorts;
  SortHelper::collectVariableSorts(clause, sorts);

  auto it = sorts.items();
  while(it.hasNext()) {
    std::pair<unsigned, TermList> next = it.next();
    unsigned var = next.first;
    TermList sort = next.second;
    // we don't support many-sorted logic yet
    ASS(sort == AtomicSort::defaultSort())
    out << "(forall iota (_" << var << " : (El iota) => ";
  }

  for(unsigned i = 0; i < clause->length(); i++) {
    bool last = i + 1 == clause->length();
    Literal *literal = (*clause)[i];

    if(!last)
      out << "(or ";
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
    if(!last)
      out << " ";
  }

  for(unsigned i = 1; i < clause->length(); i++)
    out << ")";
  for(unsigned i = 0; i < sorts.size(); i++)
    out << "))";
}

void outputAxiom(std::ostream &out, Unit *axiom) {
  CALL("Dedukti::outputAxiom")
  // we don't support non-clause inputs yet
  ASS(axiom->isClause())
  // we don't support non-axiom inputs yet
  ASS(
    axiom->inputType() == UnitInputType::AXIOM ||
    axiom->inputType() == UnitInputType::ASSUMPTION ||
    axiom->inputType() == UnitInputType::NEGATED_CONJECTURE
  )

  out << "_ax" << axiom->number() << ": ";

  Clause *clause = static_cast<Clause *>(axiom);
  outputClause(out, clause);
  out << "." << std::endl;

}

}
}