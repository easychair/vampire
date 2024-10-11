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
#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/SortHelper.hpp"
#include "Parse/TPTP.hpp"
#include "Shell/UIHelper.hpp"

#include "Inferences/BinaryResolution.hpp"
#include "Inferences/Superposition.hpp"
#include "Inferences/EqualityResolution.hpp"

#include <array>
#include <set>

using namespace Inferences;

const char *PRELUDE = R"(
(; Prop ;)
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
inhabit : El iota.

(; Equality ;)
def eq : (El iota) -> (El iota) -> Prop.
[x, y] Prf (eq x y) --> p : ((El iota) -> Prop) -> (Prf (p x) -> Prf (p y)).
def refl : x : (El iota) -> Prf (eq x x).
[x] refl x --> p : ((El iota) -> Prop) => t : Prf (p x) => t.
def comm : x : (El iota) -> y : (El iota) -> Prf (eq x y) -> Prf (eq y x).
[x, y] comm x y --> e : (Prf (eq x y)) => p : ((El iota) -> Prop) => e (z : (El iota) => imp (p z) (p x)) (t : (Prf (p x)) => t).
def comml : x : (El iota) -> y : (El iota) -> (Prf (eq x y) -> Prf false) -> (Prf (eq y x) -> Prf false).
[x, y] comml x y --> l : (Prf (eq x y) -> Prf false) => e : Prf (eq y x) => l (comm y x e).

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
static void outputVar(std::ostream &out, unsigned var, bool special, Care care) {
  if(special)
    out << "z";
  else if(care(var))
    out << var;
  else
    out << "inhabit";
}

static void outputName(std::ostream &out, const std::string &name) {
  out << "{|" << name << "|}";
}

template<typename Care>
static void outputArgs(std::ostream &out, TermList *start, Care care) {
  ASS(start->isNonEmpty())

  Stack<TermList *> todo;
  TermList *current = start;
  while(true) {
    out << " ";
    if(current->isVar()) {
      outputVar(out, current->var(), current->isSpecialVar(), care);
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

template<typename Care>
static void outputTermList(std::ostream &out, TermList tl, Care care) {
  if(tl.isVar()) {
    outputVar(out, tl.var(), tl.isSpecialVar(), care);
    return;
  }

  else if(tl.isTerm()) {
    Term *term = tl.term();
    if(term->arity())
      out << "(";
    outputName(out, term->functionName());
    if(term->arity())
      outputArgs(out, term->args(), care);
    if(term->arity())
      out << ")";
  }
}

template<typename Care>
static void outputLiteral(std::ostream &out, Literal *literal, Care care) {
  if(!literal->polarity())
    out << "(not ";
  if(literal->arity())
    out << "(";

  if(literal->isEquality())
    out << "eq";
  else
    outputName(out, literal->predicateName());
  if(literal->arity())
    outputArgs(out, literal->args(), care);

  if(literal->arity())
    out << ")";
  if(!literal->polarity())
    out << ")";
}

template<typename Care>
static void outputSubstitutedLiteralPtr(std::ostream &out, RobSubstitution &subst, unsigned bank, Literal *literal, Care care) {
  if(!literal->isEquality()) {
    out << subst.apply(literal, bank);
    return;
  }

  TermList leftSubst = subst.apply(literal->termArg(0), bank);
  TermList rightSubst = subst.apply(literal->termArg(1), bank);
  Literal *after_subst = subst.apply(literal, bank);
  bool need_swap = after_subst->termArg(0) != leftSubst;
  if(need_swap) {
    out << "(comml ";
    outputTermList(out, rightSubst, care);
    out << " ";
    outputTermList(out, leftSubst, care);
    out << " ";
  }
  out << after_subst;
  if(need_swap)
    out << ")";
}

static void outputClause(std::ostream &out, Clause *clause) {
  if(clause->isEmpty()) {
    out << "Prf_clause (cl ec)";
    return;
  }

  out << "Prf_clause ";
  auto vars = variables(clause);
  for(unsigned var : vars)
    out << "(bind iota (" << var << " : El iota => ";

  out << "(cl ";
  auto canonical = canonicalise(clause);
  for(Literal *literal : canonical) {
    out << "(cons ";
    outputLiteral(out, literal, AlwaysCare {});
    out << " ";
  }
  out << "ec";

  for(unsigned i = 0; i < clause->length(); i++)
    out << ")";
  for(unsigned i = 0; i < vars.size(); i++)
    out << "))";

  out << ")";
}

static void outputAxiomName(std::ostream &out, Unit *axiom) {
  std::string recoveredName;
  out << "axiom_";
  if(Parse::TPTP::findAxiomName(axiom, recoveredName))
    out << recoveredName;
  else
    out << axiom->number();
}

static void outputDeductionPrefix(std::ostream &out, Unit *deduction) {
  // we don't support non-clause deductions yet
  ASS(deduction->isClause())

  out << "def deduction" << deduction->number() << ": ";
  Clause *clause = static_cast<Clause *>(deduction);
  outputClause(out, clause);
  out << " := ";
}

static void sorry(std::ostream &out, Unit *admit) {
  // we don't support non-clause deductions yet
  ASS(admit->isClause())

  Clause *clause = static_cast<Clause *>(admit);
  out << "sorry" << admit->number() << ": ";

  UnitIterator parents = admit->getParents();
  while(parents.hasNext()) {
    Unit *parent = parents.next();
    ASS(parent->isClause());
    Clause *parent_clause = static_cast<Clause *>(parent);
    outputClause(out, parent_clause);
    out << " -> ";
  }
  outputClause(out, clause);
  out << ".\n";

  outputDeductionPrefix(out, admit);
  out << "sorry" << admit->number();
  parents = admit->getParents();
  while(parents.hasNext())
    out << " deduction" << parents.next()->number();
}

// get N parents of a unit
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

static void outputResolution(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);

  auto [left, right] = getParents<2>(derived);
  const auto &br = env.proofExtra.get<Inferences::BinaryResolutionExtra>(derived);

  // compute unifier for selected literals
  RobSubstitution subst;
  Literal *selectedLeft = br.selectedLiteral.selectedLiteral;
  Literal *selectedRight = br.otherLiteral;
  ASS_NEQ(selectedLeft->polarity(), selectedRight->polarity())
  ALWAYS(subst.unify(TermList(selectedLeft), 0, TermList(selectedRight), 1));

  // apply subst to all of the parent literals in the same order as BinaryResolution does it
  // this will, I fervently hope, ensure that output variables are mapped in the same way
  for(unsigned i = 0; i < left->length(); i++)
    if((*left)[i] != selectedLeft)
      subst.apply((*left)[i], 0);
  for(unsigned i = 0; i < right->length(); i++)
    if((*right)[i] != selectedRight)
      subst.apply((*right)[i], 1);

  // now also apply subst to the selected literals because we need it later
  Literal *leftSelectedSubst = subst.apply(selectedLeft, 0);
  Literal *rightSelectedSubst = subst.apply(selectedRight, 1);

  // canonicalise order of literals in all clauses
  auto derivedLits = canonicalise(derived);
  auto litsLeft = canonicalise(left);
  auto litsRight = canonicalise(right);

  // variables in numerical order
  auto derivedVars = variables(derived);
  auto leftVars = variables(left);
  auto rightVars = variables(right);

  // for variables in the substitution that do not appear in the output
  // consider e.g. p(X) and ~p(Y): X -> Y, but output is $false and has no variables
  auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

  // bind variables present in the derived clause
  for(unsigned v : derivedVars)
    out << " " << v << " : El iota => ";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << "" << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) => ";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  out << "deduction" << left->number();

  for(unsigned v : leftVars) {
    out << " ";
    outputTermList(out, subst.apply(TermList(v, false), 0), care);
  }
  unsigned litLeft;
  for(litLeft = 0; litsLeft[litLeft] != selectedLeft; litLeft++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsLeft[litLeft], care);
  }

  const char *tp = "tp", *tnp = "tnp";
  if(selectedLeft->isNegative())
    std::swap(tp, tnp);

  out << " (" << tp << ": (Prf ";
  outputLiteral(out, leftSelectedSubst, care);
  out << ") => " << "deduction" << right->number();
  for(unsigned v : rightVars) {
    out << " ";
    outputTermList(out, subst.apply(TermList(v, false), 1), care);
  }
  unsigned litRight;
  for(litRight = 0; litsRight[litRight] != selectedRight; litRight++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 1, litsRight[litRight], care);
  }
  out << " (" << tnp << ": Prf ";
  outputLiteral(out, rightSelectedSubst, care);

  out << " => (tnp tp)";
  out << ")";
  for(litRight++; litRight < litsRight.size(); litRight++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 1, litsRight[litRight], care);
  }
  out << ")";
  for(litLeft++; litLeft < litsLeft.size(); litLeft++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsLeft[litLeft], care);
  }
}

static void outputSuperposition(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);

  auto [left, right] = getParents<2>(derived);
  const auto &sup = env.proofExtra.get<Inferences::SuperpositionExtra>(derived);

  // compute unifier for selected literals
  RobSubstitution subst;
  Literal *selectedLeft = sup.selected.selectedLiteral.selectedLiteral;
  Literal *selectedRight = sup.selected.otherLiteral;
  TermList from = sup.rewrite.lhs;
  TermList to = EqHelper::getOtherEqualitySide(selectedRight, from);
  bool fromisLHS = selectedRight->termArg(0) == from;
  TermList target = sup.rewrite.target;
  ASS(selectedRight->isEquality())
  ASS(selectedRight->isPositive())
  ASS(selectedRight->termArg(0) == from || selectedRight->termArg(1) == from)
  ALWAYS(subst.unify(target, 0, from, 1))

  // apply subst to all of the parent literals in the same order as Superposition does it
  // this will, I fervently hope, ensure that output variables are mapped in the same way
  for(unsigned i = 0; i < left->length(); i++)
    if((*left)[i] != selectedLeft)
      subst.apply((*left)[i], 0);
  for(unsigned i = 0; i < right->length(); i++)
    if((*right)[i] != selectedRight)
      subst.apply((*right)[i], 1);

  // now also apply subst to the selected literals because we need it later
  Literal *leftSelectedSubst = subst.apply(selectedLeft, 0);
  TermList fromSubst = subst.apply(from, 1);
  TermList toSubst = subst.apply(to, 1);
  Literal *leftSelectedSubstRewritten = EqHelper::replace(leftSelectedSubst, fromSubst, toSubst);

  // canonicalise order of literals in all clauses
  auto derivedLits = canonicalise(derived);
  auto litsLeft = canonicalise(left);
  auto litsRight = canonicalise(right);

  // variables in numerical order
  auto derivedVars = variables(derived);
  auto leftVars = variables(left);
  auto rightVars = variables(right);

  // for variables in the substitution that do not appear in the output
  // consider e.g. p(X) and ~p(Y): X -> Y, but output is $false and has no variables
  auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

  // bind variables present in the derived clause
  for(unsigned v : derivedVars)
    out << " " << v << " : El iota => ";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << "" << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) => ";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  out << "deduction" << left->number();

  for(unsigned v : leftVars) {
    out << " ";
    outputTermList(out, subst.apply(TermList(v, false), 0), care);
  }
  unsigned litLeft;
  for(litLeft = 0; litsLeft[litLeft] != selectedLeft; litLeft++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsLeft[litLeft], care);
  }

  // TODO q could also be one of the swappy things under the substitution
  // life is pain
  out << " (q : (Prf (";
  outputLiteral(out, leftSelectedSubst, care);
  out << ")) => deduction" << right->number();

  for(unsigned v : rightVars) {
    out << " ";
    outputTermList(out, subst.apply(TermList(v, false), 1), care);
  }
  unsigned litRight;
  for(litRight = 0; litsRight[litRight] != selectedRight; litRight++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 1, litsRight[litRight], care);
  }

  out << " (r : (Prf (eq ";
  outputTermList(out, subst.apply(selectedRight->termArg(0), 1), care);
  out << " ";
  outputTermList(out, subst.apply(selectedRight->termArg(1), 1), care);
  out << ")) => "  << leftSelectedSubstRewritten << " (";

  if(fromisLHS)
    out << "r";
  else {
    out << "(comm ";
    outputTermList(out, subst.apply(to, 1), care);
    out << " ";
    outputTermList(out, subst.apply(from, 1), care);
    out << " r)";
  }
  out << " (z : (El iota) => ";

  TermList z(0, true);
  Literal *replaced = EqHelper::replace(leftSelectedSubst, subst.apply(from, 1), z);
  outputLiteral(out, replaced, care);
  out << ") q))";

  for(litRight++; litRight < litsRight.size(); litRight++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 1, litsRight[litRight], care);
  }
  out << ")";
  for(litLeft++; litLeft < litsLeft.size(); litLeft++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsLeft[litLeft], care);
  }
}

static void outputEqualityResolution(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);

  auto [parent] = getParents<1>(derived);
  const auto &er = env.proofExtra.get<Inferences::EqualityResolutionExtra>(derived);

  // compute unifier for selected literal
  RobSubstitution subst;
  Literal *selected = er.selectedLiteral;
  ASS(selected->isEquality())
  ASS(selected->isNegative())
  TermList s = selected->termArg(0), t = selected->termArg(1);
  ALWAYS(subst.unify(s, 0, t, 0));

  // apply subst to all of the parent literals in the same order as EqualityResolution does it
  // this will, I fervently hope, ensure that output variables are mapped in the same way
  for(unsigned i = 0; i < parent->length(); i++)
    if((*parent)[i] != selected)
      subst.apply((*parent)[i], 0);

  // canonicalise order of literals in all clauses
  auto derivedLits = canonicalise(derived);
  auto litsParent = canonicalise(parent);

  // variables in numerical order
  auto derivedVars = variables(derived);
  auto parentVars = variables(parent);

  // for variables in the substitution that do not appear in the output
  // consider e.g. p(X) and ~p(Y): X -> Y, but output is $false and has no variables
  auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

  // bind variables present in the derived clause
  for(unsigned v : derivedVars)
    out << " " << v << " : El iota => ";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << "" << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) => ";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  out << "deduction" << parent->number();

  for(unsigned v : parentVars) {
    out << " ";
    outputTermList(out, subst.apply(TermList(v, false), 0), care);
  }
  unsigned lit;
  for(lit = 0; litsParent[lit] != selected; lit++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsParent[lit], care);
  }

  out << " (p : Prf (";
  outputLiteral(out, subst.apply(selected, 0), care);
  out << ") => p (refl ";
  outputTermList(out, subst.apply(s, 0), care);
  out << "))";

  for(lit++; lit < litsParent.size(); lit++) {
    out << " ";
    outputSubstitutedLiteralPtr(out, subst, 0, litsParent[lit], care);
  }
}


namespace Shell {
namespace Dedukti {

void outputPrelude(std::ostream &out) {
  out << PRELUDE;
}

void outputTypeDecl(std::ostream &out, const std::string &name, OperatorType *type) {
  outputName(out, name);
  out << ": ";

  // we don't support polymorphism yet
  ASS_EQ(type->numTypeArguments(), 0)
  // we don't support many-sorted logic yet
  ASS(type->isAllDefault())

  for(unsigned i = 0; i < type->arity(); i++)
    out << "El iota -> ";

  TermList range = type->result();
  // we don't support many-sorted logic yet
  ASS(range.isEmpty() || range == AtomicSort::defaultSort())

  // predicate
  if(range.isEmpty())
    out << "Prop";
  // function
  else
    out << "El iota";

  out << ".\n";
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
  out << ".\n";
}

void outputDeduction(std::ostream &out, Unit *u) {
  switch(u->inference().rule())
  {
  case InferenceRule::INPUT:
  case InferenceRule::NEGATED_CONJECTURE:
    outputAxiom(out, u);
    outputDeductionPrefix(out, u);
    outputAxiomName(out, u);
    break;
  case InferenceRule::RESOLUTION:
    outputResolution(out, u->asClause());
    break;
  case InferenceRule::EQUALITY_RESOLUTION:
    outputEqualityResolution(out, u->asClause());
    break;
  case InferenceRule::SUPERPOSITION:
    outputSuperposition(out, u->asClause());
    break;
  default:
    sorry(out, u);
  }
  out << ".\n";
}
}
}
