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

#include "Indexing/SubstitutionTree.hpp"
#include "Inferences/BinaryResolution.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/SortHelper.hpp"
#include "Parse/TPTP.hpp"
#include "Shell/UIHelper.hpp"


#include <set>

using namespace Inferences;

const char *PRELUDE = R"(
(; !Prop ;)
!Prop : Type.
def !Prf : (!Prop -> Type).
!true : !Prop.
[] !Prf !true --> (r : !Prop -> ((!Prf r) -> (!Prf r))).
!false : !Prop.
[] !Prf !false --> (r : !Prop -> (!Prf r)).
!not : (!Prop -> !Prop).
[p] !Prf (!not p) --> ((!Prf p) -> (r : !Prop -> (!Prf r))).
!and : (!Prop -> (!Prop -> !Prop)).
[p, q] !Prf (!and p q) --> (r : !Prop -> (((!Prf p) -> ((!Prf q) -> (!Prf r))) -> (!Prf r))).
!or : (!Prop -> (!Prop -> !Prop)).
[p, q] !Prf (!or p q) --> (r : !Prop -> (((!Prf p) -> (!Prf r)) -> (((!Prf q) -> (!Prf r)) -> (!Prf r)))).
!imp : (!Prop -> (!Prop -> !Prop)).
[p, q] !Prf (!imp p q) --> ((!Prf p) -> (!Prf q)).

(; !Set ;)
!Set : Type.
injective !El : (!Set -> Type).
!iota : !Set.
!inhabit : !El !iota.

(; Equality ;)
def !eq : (!El !iota) -> (!El !iota) -> Type.
[] !eq --> x : (!El !iota) => y : (!El !iota) => p : ((!El !iota) -> !Prop) -> (!Prf (p x) -> !Prf (p y)).
def !refl : x : (!El !iota) -> !eq x x.
[] !refl --> x : (!El !iota) => p : ((!El !iota) -> !Prop) => t : !Prf (p x) => t.

(; Quant ;)
!forall : (a : !Set -> (((!El a) -> !Prop) -> !Prop)).
[a, p] !Prf (!forall a p) --> (x : (!El a) -> (!Prf (p x))).
!exists : (a : !Set -> (((!El a) -> !Prop) -> !Prop)).
[a, p] !Prf (!exists a p) --> (r : !Prop -> ((x : (!El a) -> ((!Prf (p x)) -> (!Prf r))) -> (!Prf r))).

(; Classic ;)
def !cPrf : (!Prop -> Type) := (p : !Prop => (!Prf (!not (!not p)))).
def !cand : (!Prop -> (!Prop -> !Prop)) := (p : !Prop => (q : !Prop => (!and (!not (!not p)) (!not (!not q))))).
def !cor : (!Prop -> (!Prop -> !Prop)) := (p : !Prop => (q : !Prop => (!or (!not (!not p)) (!not (!not q))))).
def !cimp : (!Prop -> (!Prop -> !Prop)) := (p : !Prop => (q : !Prop => (!imp (!not (!not p)) (!not (!not q))))).
def !cforall : (a : !Set -> (((!El a) -> !Prop) -> !Prop)) := (a : !Set => (p : ((!El a) -> !Prop) => (!forall a (x : (!El a) => (!not (!not (p x))))))).
def !cexists : (a : !Set -> (((!El a) -> !Prop) -> !Prop)) := (a : !Set => (p : ((!El a) -> !Prop) => (!exists a (x : (!El a) => (!not (!not (p x))))))).

(; Clauses ;)
def !prop_clause : Type.
def !ec : !prop_clause.
def !cons : (!Prop -> (!prop_clause -> !prop_clause)).
def !clause : Type.
def !cl : (!prop_clause -> !clause).
def !bind : (A : !Set -> (((!El A) -> !clause) -> !clause)).
def !Prf_prop_clause : (!prop_clause -> Type).

[] !Prf_prop_clause !ec --> (!Prf !false).
[p, c] !Prf_prop_clause (!cons p c) --> ((!Prf p -> !Prf !false) -> (!Prf_prop_clause c)).
def !Prf_clause : (!clause -> Type).
[c] !Prf_clause (!cl c) --> (!Prf_prop_clause c).
[A, f] !Prf_clause (!bind A f) --> (x : (!El A) -> (!Prf_clause (f x))).

(;----------------------------------------------------------------------------;)

)";

using namespace Kernel;

static std::vector<Literal *> canonicalise(Clause *cl) {
  std::vector<Literal *> lits;
  for(unsigned i = 0; i < cl->length(); i++)
    lits.push_back((*cl)[i]);
  sort(lits.begin(), lits.end());
  return lits;
}

static std::set<unsigned> variables(Clause *cl) {
  std::set<unsigned> vars;
  auto it = cl->getVariableIterator();
  while(it.hasNext())
    vars.insert(it.next());
  return vars;
}

struct AlwaysCare {
  bool operator()(unsigned _) { return true; }
};

template<typename Care>
static void outputVar(std::ostream &out, unsigned var, Care care) {
  if(care(var))
    out << "!" << var;
  else
    out << "!inhabit";
}

template<typename Care>
static void outputArgs(std::ostream &out, TermList *start, Care care) {
  ASS(start->isNonEmpty())

  Stack<TermList *> todo;
  TermList *current = start;
  while(true) {
    out << " ";
    if(current->isVar()) {
      outputVar(out, current->var(), care);
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

template<typename Care>
static void outputTermList(std::ostream &out, TermList tl, Care care) {
  if(tl.isVar()) {
    outputVar(out, tl.var(), care);
    return;
  }

  else if(tl.isTerm()) {
    Term *term = tl.term();
    if(term->arity())
      out << "(";
    out << term->functionName();
    if(term->arity())
      outputArgs(out, term->args(), care);
    if(term->arity())
      out << ")";
  }
}

template<typename Care>
static void outputLiteral(std::ostream &out, Literal *literal, Care care) {
  if(!literal->polarity())
    out << "(!not ";
  if(literal->arity())
    out << "(";

  out << literal->predicateName();
  if(literal->arity())
    outputArgs(out, literal->args(), care);

  if(literal->arity())
    out << ")";
  if(!literal->polarity())
    out << ")";
}

static void outputClause(std::ostream &out, Clause *clause) {
  if(clause->isEmpty()) {
    out << "!Prf_clause (!cl !ec)";
    return;
  }

  out << "!Prf_clause ";
  std::set<unsigned> vars = variables(clause);
  for(unsigned var : vars)
    out << "(!bind !iota (!" << var << " : (!El !iota) => ";

  out << "(!cl ";
  auto canonical = canonicalise(clause);
  for(Literal *literal : canonical) {
    out << "(!cons ";
    outputLiteral(out, literal, AlwaysCare {});
    out << " ";
  }
  out << "!ec";

  for(unsigned i = 0; i < clause->length(); i++)
    out << ")";
  for(unsigned i = 0; i < vars.size(); i++)
    out << "))";

  out << ")";
}

static void outputAxiomName(std::ostream &out, Unit *axiom) {
  std::string recoveredName;
  out << "!axiom_";
  if(Parse::TPTP::findAxiomName(axiom, recoveredName))
    out << recoveredName;
  else
    out << axiom->number();
}

static void outputDeductionPrefix(std::ostream &out, Unit *deduction) {
  // we don't support non-clause deductions yet
  ASS(deduction->isClause())

  out << "def !deduction" << deduction->number() << ": ";
  Clause *clause = static_cast<Clause *>(deduction);
  outputClause(out, clause);
  out << " := ";
}

static void sorry(std::ostream &out, Unit *admit) {
  // we don't support non-clause deductions yet
  ASS(admit->isClause())

  Clause *clause = static_cast<Clause *>(admit);
  out << "!sorry" << admit->number() << ": ";

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
  out << "!sorry" << admit->number();
  parents = admit->getParents();
  while(parents.hasNext())
    out << " !deduction" << parents.next()->number();
  out << "." << std::endl;
}


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
    out << "(!El !iota) -> ";

  TermList range = type->result();
  // we don't support many-sorted logic yet
  ASS(range.isEmpty() || range == AtomicSort::defaultSort())

  // predicate
  if(range.isEmpty())
    out << "!Prop";
  // function
  else
    out << "(!El !iota)";

  out << "." << std::endl;
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

void outputDeduction(std::ostream &out, Unit *u) {
  ASS(u->isClause())
  Clause *derived = u->asClause();
  InferenceRule rule = u->inference().rule();
  UnitIterator parents= u->getParents();

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
    outputDeductionPrefix(out, u);
    ALWAYS(parents.hasNext())
    Clause *left = parents.next()->asClause();
    ALWAYS(parents.hasNext())
    Clause *right = parents.next()->asClause();
    const BinaryResolutionExtra &br = env.proofExtra.get<Inferences::BinaryResolutionExtra>(u);

    // compute unifier for selected literals
    RobSubstitution subst;
    Literal *leftSelected = br.selectedLiteral.selectedLiteral;
    TermList latom(leftSelected);
    Literal *rightSelected = br.otherLiteral;
    TermList ratom(rightSelected);
    ALWAYS(subst.unify(latom, 0, ratom, 1));

    // apply subst to all of the parent literals in the same order as BinaryResolution does it
    // this will, please God, ensure that the variables are mapped in the same way
    for(unsigned i = 0; i < left->length(); i++)
      if((*left)[i] != leftSelected)
        subst.apply((*left)[i], 0);
    for(unsigned i = 0; i < right->length(); i++)
      if((*right)[i] != rightSelected)
        subst.apply((*right)[i], 1);

    auto derivedLits = canonicalise(derived);
    auto derivedVars = variables(derived);
    auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

    // now also apply to the selected literals because we need it later
    Literal *leftSelectedSubst = subst.apply(leftSelected, 0);
    Literal *rightSelectedSubst = subst.apply(rightSelected, 1);

    auto leftVars = variables(left);
    auto rightVars = variables(right);
    // canonicalise parent literal order
    auto litsLeft = canonicalise(left);
    auto litsRight = canonicalise(right);

    // TODO there must be a nicer way of doing this
    unsigned leftSelectedIndex = 0, rightSelectedIndex = 0;
    while(litsLeft[leftSelectedIndex] != leftSelected) leftSelectedIndex++;
    while(litsRight[rightSelectedIndex] != rightSelected) rightSelectedIndex++;

    std::vector<Literal *> substLeft, substRight;
    for(Literal *l : litsLeft)
      substLeft.push_back(subst.apply(l, 0));
    for(Literal *r : litsRight)
      substRight.push_back(subst.apply(r, 1));

    for(unsigned v : derivedVars)
      out << " !" << v << " : (!El !iota) => ";
    for(unsigned i = 0; i < derivedLits.size(); i++) {
      Literal *l = derivedLits[i];
      out << "!" << l << " : (!Prf ";
      outputLiteral(out, l, care);
      out << " -> !Prf !false) => ";
    }
    out << "!deduction" << left->number();

    for(unsigned v : leftVars) {
      out << " ";
      outputTermList(out, subst.apply(TermList(v, false), 0), care);
    }
    for(unsigned i = 0; i < leftSelectedIndex; i++) {
      Literal *l = substLeft[i];
      out << " !" << l;
    }

    const char *tp = leftSelected->isPositive() ? "!tp" : "!tnp";
    const char *tnp = rightSelected->isPositive() ? "!tp" : "!tnp";
    out << " (" << tp << ": (!Prf ";
    outputLiteral(out, leftSelectedSubst, care);
    out << ") => " << "!deduction" << right->number();
    for(unsigned v : rightVars) {
      out << " ";
      outputTermList(out, subst.apply(TermList(v, false), 1), care);
    }
    for(unsigned i = 0; i < rightSelectedIndex; i++) {
      Literal *l = substRight[i];
      out << " !" << l;
    }
    out << " (" << tnp << ": !Prf ";
    outputLiteral(out, rightSelectedSubst, care);

    out << " => (!tnp !tp)";
    out << ")";
    for(unsigned i = rightSelectedIndex + 1; i < litsRight.size(); i++) {
      Literal *l = substRight[i];
      out << " !" << l;
    }
    out << ")";
    for(unsigned i = leftSelectedIndex + 1; i < litsLeft.size(); i++) {
      Literal *l = substLeft[i];
      out << " !" << l;
    }
    out << ".\n";
    break;
  }
  default: sorry(out, u);
  }
}
}
}
