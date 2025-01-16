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

#include "Inferences/ProofExtra.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Kernel/Matcher.hpp"
#include "Kernel/MLVariant.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/TermIterators.hpp"
#include "Lib/SharedSet.hpp"
#include "Parse/TPTP.hpp"
#include "SATSubsumption/SATSubsumptionAndResolution.hpp"
#include "Saturation/SaturationAlgorithm.hpp"
#include "Saturation/Splitter.hpp"

#include "Inferences/BinaryResolution.hpp"
#include "Inferences/EqualityResolution.hpp"
#include "Inferences/Factoring.hpp"
#include "Inferences/Superposition.hpp"
#include "Shell/EqResWithDeletion.hpp"
#include "Shell/FunctionDefinition.hpp"

#include <array>
#include <unordered_map>
#include <set>
#include <sstream>

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
[p, q] Prf (or p q) --> (((Prf p) -> (Prf false)) -> (((Prf q) -> (Prf false)) -> (Prf false))).
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
def comml_not : x : (El iota) -> y : (El iota) -> (Prf (not (eq x y)) -> Prf false) -> (Prf (not (eq y x)) -> Prf false).
[x, y] comml_not x y --> l : ((Prf (eq x y) -> Prf false) -> Prf false) => ne : (Prf (eq y x) -> Prf false) => l (e : Prf (eq x y) => ne (comm x y e)).

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

def av_clause : Type.
def acl : clause -> av_clause.
def if : Prop -> av_clause -> av_clause.
def Prf_av_clause : av_clause -> Type.

[c] Prf_av_clause (acl c) --> Prf_clause c.
[sp, c] Prf_av_clause (if sp c) --> (Prf (not sp) -> Prf false) -> Prf_av_clause c.

(;----------------------------------------------------------------------------;)

)";

using namespace Kernel;

template<typename S, typename T>
struct Sequence;

template<typename Derived>
struct Transformation {
  template<typename Other>
  Sequence<Derived, Other> then(Other &&other) {
    return Sequence<Derived, Other> { std::move(*static_cast<Derived *>(this)), other };
  }
};

template<typename S, typename T>
struct Sequence : public Transformation<Sequence<S, T>> {
  S s;
  T t;
  Sequence(S s, T t) : s(s), t(t) {}
  Literal *operator()(Literal *literal) { return t(s(literal)); }
  TermList operator()(TermList term) { return t(s(term)); }
};

struct DoSubstitution : public Transformation<DoSubstitution> {
  SimpleSubstitution &substitution;
  DoSubstitution(SimpleSubstitution &substitution) : substitution(substitution) {}

  Literal *operator()(Literal *literal) { return SubstHelper::apply(literal, substitution); }
  TermList operator()(TermList term) { return SubstHelper::apply(term, substitution); }
};

template<unsigned bank>
struct DoRobSubstitution : public Transformation<DoRobSubstitution<bank>> {
  RobSubstitution &substitution;
  DoRobSubstitution(RobSubstitution &substitution) : substitution(substitution) {}

  Literal *operator()(Literal *literal) { return substitution.apply(literal, bank); }
  TermList operator()(TermList term) { return substitution.apply(term, bank); }
};
using DoRobSubstitution0 = DoRobSubstitution<0>;
using DoRobSubstitution1 = DoRobSubstitution<1>;

struct DoReplacement : public Transformation<DoReplacement> {
  TermList from, to;
  DoReplacement(TermList from, TermList to) : from(from), to(to) {}
  Literal *operator()(Literal *literal) { return EqHelper::replace(literal, from, to); }

  TermList operator()(TermList term) {
    if(term.isVar())
      return term;
    // EqHelper::replace only replaces _occurrences_
    if(term == from)
      return to;
    return TermList(EqHelper::replace(term.term(), from, to));
  }
};


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
static void outputTerm(std::ostream &out, TermList tl, Care care) {
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


template<typename Transform, typename Care>
static void outputLiteral(
  std::ostream &out,
  Literal *literal,
  Care care,
  Transform transform
) {
  if(!literal->isEquality())
    return outputLiteral(out, transform(literal), care);

  if(!literal->polarity())
    out << "(not ";
  out << "(eq ";
  outputTerm(out, transform(literal->termArg(0)), care);
  out << " ";
  outputTerm(out, transform(literal->termArg(1)), care);
  out << ")";
  if(!literal->polarity())
    out << ")";
}

template<typename Transform, typename Care>
static void outputLiteralName(
  std::ostream &out,
  Literal *literal,
  const char *name,
  const char *comm,
  Care care,
  Transform transform
) {
  if(!literal->isEquality()) {
    out << name;
    return;
  }

  Literal *after = transform(literal);
  TermList leftAfter = transform(literal->termArg(0));
  TermList rightAfter = transform(literal->termArg(1));
  bool need_commute = after->termArg(0) != leftAfter || after->termArg(1) != rightAfter;
  if(need_commute) {
    out << "(" << comm << " ";
    // TODO this is switched with respect to `outputLiteralPtr` - come up with some nice abstraction for this
    outputTerm(out, leftAfter, care);
    out << " ";
    outputTerm(out, rightAfter, care);
    out << " ";
  }
  out << name;
  if(need_commute)
    out << ")";
}

template<typename Care, typename Transform>
static void outputLiteralPtr(
  std::ostream &out,
  Literal *literal,
  Care care,
  Transform transform
) {
  if(!literal->isEquality()) {
    out << transform(literal);
    return;
  }

  Literal *after = transform(literal);
  TermList leftAfter = transform(literal->termArg(0));
  TermList rightAfter = transform(literal->termArg(1));
  bool need_commute = after->termArg(0) != leftAfter || after->termArg(1) != rightAfter;
  if(need_commute) {
    out << "(comml";
    if(literal->isNegative())
      out << "_not";
    out << " ";
    outputTerm(out, rightAfter, care);
    out << " ";
    outputTerm(out, leftAfter, care);
    out << " ";
  }
  out << after;
  if(need_commute)
    out << ")";
}

static void outputSplit(std::ostream &out, SATLiteral split, bool flip = false) {
  if(split.polarity() == flip)
    out << "(not ";
  out << "sp" << split.var();
  if(split.polarity() == flip)
    out << ")";
}

static void outputSplit(std::ostream &out, unsigned splitName, bool flip = false) {
  outputSplit(out, Splitter::getLiteralFromName(splitName), flip);
}

static void outputClause(std::ostream &out, Clause *clause) {
  unsigned close_parens = 0;
  if(clause->noSplits())
    out << "Prf_clause";
  else {
    out << "Prf_av_clause";
    SplitSet &splits = *clause->splits();
    for(unsigned i = 0; i < splits.size(); i++) {
      out << " (if";
      SATLiteral split = Splitter::getLiteralFromName(splits[i]);
      if(!split.polarity())
        out << " (not";
      out << " sp" << split.var();
      if(!split.polarity())
        out << ")";
    }
    out << " (acl";
    close_parens = splits.size() + 1;
  }
  auto vars = variables(clause);
  for(unsigned var : vars)
    out << " (bind iota (" << var << " : El iota =>";
  close_parens += 2 * vars.size();

  out << " (cl";
  close_parens++;
  auto canonical = canonicalise(clause);
  for(Literal *literal : canonical) {
    out << " (cons ";
    outputLiteral(out, literal, AlwaysCare {});
  }
  close_parens += canonical.size();
  out << " ec";

  for(unsigned i = 0; i < close_parens; i++)
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

static void outputSplitSet(std::ostream &out, SplitSet &splits) {
  for(unsigned i = 0; i < splits.size(); i++) {
    SATLiteral split = Splitter::getLiteralFromName(splits[i]);
    out << " nnsp" << split.var();
    out << " : (Prf";
    if(!split.polarity())
      out << " (not";
    out <<  " (not sp" << split.var();
    if(!split.polarity())
      out << ")";
    out << ") -> Prf false) =>";
  }
}

static void outputDeductionPrefix(std::ostream &out, Unit *deduction) {
  // we don't support non-clause deductions yet
  ASS(deduction->isClause())

  out << "def deduction" << deduction->number() << ": ";
  Clause *clause = static_cast<Clause *>(deduction);
  outputClause(out, clause);
  out << " := ";

  if(clause->noSplits())
    return;

  outputSplitSet(out, *clause->splits());
}

static void outputParentWithSplits(std::ostream &out, Clause *parent) {
  out << "deduction" << parent->number();
  if(parent->noSplits())
    return;
  SplitSet &splits = *parent->splits();
  for(unsigned i = 0; i < splits.size(); i++)
    out << " nnsp" << Splitter::getLiteralFromName(splits[i]).var();
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

  out << "#PRINT \"sorry: " << ruleName(admit->inference().rule()) << "\".\n";

  outputDeductionPrefix(out, admit);
  out << "sorry" << admit->number();
  parents = admit->getParents();
  while(parents.hasNext())
    out << " deduction" << parents.next()->number();
}

// get N parents of a unit
template<unsigned N>
std::array<Clause *, N> getParents(Unit *unit) {
  std::array<Clause *, N> parents;
  UnitIterator it = unit->getParents();
  for(unsigned i = 0; i < N; i++) {
    ALWAYS(it.hasNext())
    parents[i] = static_cast<Clause *>(it.next());
  }
  ALWAYS(!it.hasNext())
  return parents;
}

// get N parents of a unit
template<unsigned N>
std::pair<std::array<Clause *, N>, std::vector<Clause *>> getParentsVariadic(Unit *unit) {
  std::array<Clause *, N> parents;
  UnitIterator it = unit->getParents();
  for(unsigned i = 0; i < N; i++) {
    ALWAYS(it.hasNext())
    parents[i] = static_cast<Clause *>(it.next());
  }
  std::vector<Clause *> rest;
  while(it.hasNext())
    rest.push_back(static_cast<Clause *>(it.next()));
  return {parents, rest};
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
  outputParentWithSplits(out, left);

  for(unsigned v : leftVars) {
    out << " ";
    outputTerm(out, subst.apply(TermList(v, false), 0), care);
  }
  unsigned litLeft;
  for(litLeft = 0; litsLeft[litLeft] != selectedLeft; litLeft++) {
    out << " ";
    outputLiteralPtr(out, litsLeft[litLeft], care, DoRobSubstitution0(subst));
  }

  const char *tp = "tp", *tnp = "tnp";
  if(selectedLeft->isNegative())
    std::swap(tp, tnp);

  out << " (" << tp << ": (Prf ";
  outputLiteral(out, leftSelectedSubst, care);
  out << ") => ";
  outputParentWithSplits(out, right);
  for(unsigned v : rightVars) {
    out << " ";
    outputTerm(out, subst.apply(TermList(v, false), 1), care);
  }
  unsigned litRight;
  for(litRight = 0; litsRight[litRight] != selectedRight; litRight++) {
    out << " ";
    outputLiteralPtr(out, litsRight[litRight], care, DoRobSubstitution1(subst));
  }
  out << " (" << tnp << ": Prf ";
  outputLiteral(out, rightSelectedSubst, care);

  out << " => (tnp tp)";
  out << ")";
  for(litRight++; litRight < litsRight.size(); litRight++) {
    out << " ";
    outputLiteralPtr(out, litsRight[litRight], care, DoRobSubstitution1(subst));
  }
  out << ")";
  for(litLeft++; litLeft < litsLeft.size(); litLeft++) {
    out << " ";
    outputLiteralPtr(out, litsLeft[litLeft], care, DoRobSubstitution0(subst));
  }
}

static void outputSubsumptionResolution(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [left, right] = getParents<2>(derived);
  auto sr = env.proofExtra.get<Inferences::LiteralInferenceExtra>(derived);
  Literal *m = sr.selectedLiteral;

  // reconstruct match by calling into the SATSR code
  SATSubsumption::SATSubsumptionAndResolution satSR;
  ALWAYS(satSR.checkSubsumptionResolutionWithLiteral(right, left, left->getLiteralPosition(m)))
  auto subst = satSR.getBindingsForSubsumptionResolutionWithLiteral();

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
    out << " " << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) => ";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  //
  // TODO this should be more or less the 'resolution' term, but with multiple tnp/tp lambdas factored out
  // Anja is a smart cookie
  outputParentWithSplits(out, left);

  for(unsigned v : leftVars) {
    out << " ";
    outputVar(out, v, false, care);
  }

  // TODO can we make this case distinction less unpleasant?
  if(m->isPositive()) {
    for(Literal *l : litsLeft) {
      out << " ";
      if(l != m) {
        out << l;
        continue;
      }

      out << "(tp: (Prf ";
      outputLiteral(out, m, care);
      out << ") => ";
      outputParentWithSplits(out, right);
      for(unsigned v : rightVars) {
        out << " ";
        outputTerm(out, SubstHelper::apply(TermList(v, false), subst), care);
      }
      for(Literal *k : litsRight) {
        out << " ";
        Literal *ksubst = SubstHelper::apply(k, subst);
        if(Literal::complementaryLiteral(ksubst) != m) {
          outputLiteralPtr(out, k, care, DoSubstitution(subst));
          continue;
        }
        out << "(tnp: Prf ";
        outputLiteral(out, k, care, DoSubstitution(subst));
        out << " => (";
        outputLiteralName(out, k, "tnp", "comml", care, DoSubstitution(subst));
        out << " tp))";
      }
      out << ")";
    }
  }
  else {
    for(Literal *l : litsLeft) {
      out << " ";
      if(l != m) {
        out << l;
        continue;
      }

      out << "(tnp: (Prf ";
      outputLiteral(out, m, care);
      out << ") => ";
      outputParentWithSplits(out, right);
      for(unsigned v : rightVars) {
        out << " ";
        outputTerm(out, SubstHelper::apply(TermList(v, false), subst), care);
      }
      for(Literal *k : litsRight) {
        out << " ";
        Literal *ksubst = SubstHelper::apply(k, subst);
        if(Literal::complementaryLiteral(ksubst) != m) {
          outputLiteralPtr(out, k, care, DoSubstitution(subst));
          continue;
        }
        out << "(tp: Prf ";
        outputLiteral(out, k, care, DoSubstitution(subst));
        out << " => (tnp ";
        outputLiteralName(out, k, "tp", "comm", care, DoSubstitution(subst));
        out << "))";
      }
      out << ")";
    }
  }
}

static void outputSuperposition(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [left, right] = getParents<2>(derived);
  auto sup = env.proofExtra.get<Inferences::SuperpositionExtra>(derived);

  // compute unifier for selected literals
  RobSubstitution subst;
  Literal *leftSelected = sup.selected.selectedLiteral.selectedLiteral;
  Literal *rightSelected = sup.selected.otherLiteral;
  TermList from = sup.rewrite.lhs;
  TermList to = EqHelper::getOtherEqualitySide(rightSelected, from);
  bool fromisLHS = rightSelected->termArg(0) == from;
  TermList target = sup.rewrite.rewritten;
  ASS(rightSelected->isEquality())
  ASS(rightSelected->isPositive())
  ASS(rightSelected->termArg(0) == from || rightSelected->termArg(1) == from)
  ALWAYS(subst.unify(target, 0, from, 1))

  // apply subst to all of the parent literals in the same order as Superposition does it
  // this will, I fervently hope, ensure that output variables are mapped in the same way
  subst.apply(to, 1);
  subst.apply(leftSelected, 0);
  subst.apply(target, 0);
  subst.apply(from, 1);
  for(unsigned i = 0; i < left->length(); i++)
    if((*left)[i] != leftSelected)
      subst.apply((*left)[i], 0);
  for(unsigned i = 0; i < right->length(); i++)
    if((*right)[i] != rightSelected)
      subst.apply((*right)[i], 1);

  // now also apply subst to the selected literals because we need it later
  TermList fromSubst = subst.apply(from, 1);
  TermList toSubst = subst.apply(to, 1);

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
    out << " " << v << " : El iota =>";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << " " << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) =>";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  outputParentWithSplits(out, left);

  for(unsigned v : leftVars) {
    out << " ";
    outputTerm(out, subst.apply(TermList(v, false), 0), care);
  }
  for(Literal *litLeft : litsLeft) {
    out << " ";

    bool needs_rewrite = find(
      derivedLits.begin(),
      derivedLits.end(),
      subst.apply(litLeft, 0)
    ) == derivedLits.end();
    if(!needs_rewrite) {
      outputLiteralPtr(out, litLeft, care, DoRobSubstitution0(subst));
      continue;
    }

    out << " (q : (Prf (";
    outputLiteral(out, litLeft, care, DoRobSubstitution0(subst));
    out << ")) => ";
    outputParentWithSplits(out, right);

    for(unsigned v : rightVars) {
      out << " ";
      outputTerm(out, subst.apply(TermList(v, false), 1), care);
    }
    unsigned litRight;
    for(litRight = 0; litsRight[litRight] != rightSelected; litRight++) {
      out << " ";
      outputLiteralPtr(out, litsRight[litRight], care, DoRobSubstitution1(subst));
    }

    out << " (r : (Prf ";
    outputLiteral(out, rightSelected, care, DoRobSubstitution1(subst));
    out << ") => ";
    outputLiteralPtr(out, litLeft, care, DoRobSubstitution0(subst).then(DoReplacement(fromSubst, toSubst)));
    out << " (";

    if(!fromisLHS) {
      out << "comm ";
      outputTerm(out, toSubst, care);
      out << " ";
      outputTerm(out, fromSubst, care);
      out << " ";
    }
    out << "r (z : (El iota) => ";

    TermList z(0, true);
    outputLiteral(out, litLeft, care, DoRobSubstitution0(subst).then(DoReplacement(fromSubst, z)));
    out << ") q))";

    for(litRight++; litRight < litsRight.size(); litRight++) {
      out << " ";
      outputLiteralPtr(out, litsRight[litRight], care, DoRobSubstitution1(subst));
    }

    out << ")";
  }
}

static bool isDemodulatorFor(Literal *demodulator, TermList target) {
  ASS(demodulator->isEquality())
  ASS(demodulator->isPositive())

  SimpleSubstitution subst;
  if(!MatchingUtils::matchTerms((*demodulator)[0], target, subst))
    return false;

  Ordering &ordering = *Ordering::tryGetGlobalOrdering();
  return ordering.compare(
    SubstHelper::apply((*demodulator)[0], subst),
    SubstHelper::apply((*demodulator)[1], subst)
  ) == Ordering::GREATER;
}

static void outputDemodulation(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [left, right] = getParents<2>(derived);
  ASS_EQ(right->length(), 1)
  auto rw = env.proofExtra.get<Inferences::RewriteInferenceExtra>(derived);

  // compute unifier for selected literals
  SimpleSubstitution subst;
  Literal *rightLit = (*right)[0];
  TermList target = rw.rewritten;
  TermList from = isDemodulatorFor(rightLit, target)
    ? (*rightLit)[0]
    : (*rightLit)[1];
  TermList to = EqHelper::getOtherEqualitySide(rightLit, from);
  bool fromisLHS = rightLit->termArg(0) == from;
  ASS(rightLit->isEquality())
  ASS(rightLit->isPositive())
  ASS(rightLit->termArg(0) == from || rightLit->termArg(1) == from)
  ALWAYS(MatchingUtils::matchTerms(from, target, subst))

  // now also apply subst to the selected literals because we need it later
  TermList toSubst = SubstHelper::apply(to, subst);

  // canonicalise order of literals in all clauses
  auto derivedLits = canonicalise(derived);
  auto litsLeft = canonicalise(left);

  // variables in numerical order
  auto derivedVars = variables(derived);
  auto leftVars = variables(left);
  auto rightVars = variables(right);

  // for variables in the substitution that do not appear in the output
  // consider e.g. p(X) and ~p(Y): X -> Y, but output is $false and has no variables
  auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

  // bind variables present in the derived clause
  for(unsigned v : derivedVars)
    out << " " << v << " : El iota =>";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << " " << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) =>";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  outputParentWithSplits(out, left);

  for(unsigned v : leftVars) {
    out << " ";
    outputVar(out, v, false, care);
  }
  for(Literal *litLeft : litsLeft) {
    out << " ";

    bool needs_rewrite = find(
      derivedLits.begin(),
      derivedLits.end(),
      litLeft
    ) == derivedLits.end();
    if(!needs_rewrite) {
      out << litLeft;
      continue;
    }

    out << " (q : (Prf (";
    outputLiteral(out, litLeft, care);
    out << ")) => ";
    outputParentWithSplits(out, right);

    for(unsigned v : rightVars) {
      out << " ";
      outputTerm(out, SubstHelper::apply(TermList(v, false), subst), care);
    }
    out << " (r : (Prf ";
    outputLiteral(out, rightLit, care, DoSubstitution(subst));
    out << ") => ";
    outputLiteralPtr(out, litLeft, care, DoReplacement(target, toSubst));
    out << " (";

    if(!fromisLHS) {
      out << "comm ";
      outputTerm(out, toSubst, care);
      out << " ";
      outputTerm(out, target, care);
      out << " ";
    }
    out << "r (z : (El iota) => ";

    TermList z(0, true);
    outputLiteral(out, litLeft, care, DoReplacement(target, z));
    out << ") q)))";
  }
}

static void outputDefinitionUnfolding(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [parents, defs] = getParentsVariadic<1>(derived);
  auto [parent] = parents;
  const auto &rw = env.proofExtra.get<Inferences::FunctionDefinitionExtra>(derived);

  // canonicalise order of literals
  auto derivedLits = canonicalise(derived);
  auto parentLits = canonicalise(parent);
  // variables in numerical order
  auto derivedVars = variables(derived);
  auto parentVars = variables(parent);

  // for variables in the substitution that do not appear in the output
  // consider e.g. p(X) and ~p(Y): X -> Y, but output is $false and has no variables
  auto care = [&](unsigned var) -> bool { return derivedVars.count(var); };

  // bind variables present in the derived clause
  for(unsigned v : derivedVars)
    out << " " << v << " : El iota =>";
  // bind literals in the derived clause
  for(unsigned i = 0; i < derivedLits.size(); i++) {
    Literal *l = derivedLits[i];
    out << " " << l << " : (Prf ";
    outputLiteral(out, l, care);
    out << " -> Prf false) =>";
  }

  out << " deduction" << parent->number();
  for(unsigned v : derivedVars) {
    out << " ";
    outputVar(out, v, false, care);
  }

  std::unordered_map<unsigned, unsigned> rewrites;
  for(unsigned i = 0; i < defs.size(); i++)
    rewrites.insert({rw.lhs[i]->functor(), i});

  std::vector<std::string> after;
  for(Literal *l : parentLits) {
restart:
    out << " (t_" << l << " : Prf ";
    outputLiteral(out, l, care);
    out << " =>";
    NonVariableIterator subterms(l);
    while(subterms.hasNext()) {
      Term *candidate = subterms.next().term();
      if(auto found = rewrites.find(candidate->functor()); found != rewrites.end()) {
        auto [_, i] = *found;
        Clause *def = defs[i];
        Literal *eq = (*def)[0];
        TermList lhs(rw.lhs[i]);
        TermList rhs = EqHelper::getOtherEqualitySide(eq, lhs);
        SimpleSubstitution subst;
        ALWAYS(MatchingUtils::matchTerms(TermList(lhs), TermList(candidate), subst))
        TermList lhsSubst(SubstHelper::apply(lhs, subst));
        TermList rhsSubst = SubstHelper::apply(rhs, subst);
        out << " deduction" << def->number();

        auto defVars = variables(def);
        for(unsigned v : defVars) {
          out << " ";
          outputTerm(out, subst.apply(v), care);
        }
        out << " (e_" << eq << " : Prf ";
        outputLiteral(out, eq, care, DoSubstitution(subst));
        out << " => ";

        std::stringstream deferred;

        DoReplacement transform(lhsSubst, rhsSubst);
        Literal *lAfter = transform(l);

        // decide if needs commutation
        bool needs_comm = false;
        if(l->isEquality()) {
          TermList leftAfter = transform(l->termArg(0));
          TermList rightAfter = transform(l->termArg(1));
          needs_comm = lAfter->termArg(0) != leftAfter || lAfter->termArg(1) != rightAfter;

          if (needs_comm) {
            deferred << "(comm" << (lAfter->polarity() ? "" : "l") << " ";
            outputTerm(deferred, leftAfter, care);
            deferred << " ";
            outputTerm(deferred, rightAfter, care);
          }
        }

        deferred << " (";
        if(eq->termArg(0) == lhs)
          deferred << "e_" << eq;
        else {
          deferred << "(comm ";
          outputTerm(
            deferred,
            SubstHelper::apply(eq->termArg(0), subst),
            care
          );
          deferred << " ";
          outputTerm(
            deferred,
            SubstHelper::apply(eq->termArg(1), subst),
            care
          );
          deferred << " e_" << eq << ")";
        }
        deferred << " (z : El iota => ";
        outputLiteral(deferred, l, care, DoReplacement(lhsSubst, TermList(0, true)));
        deferred << ") t_" << l << ")))";
        if (needs_comm){
          deferred << ")";
        }
        l = lAfter;

        after.push_back(deferred.str());
        goto restart;
      }
    }
    out << " " << l << " t_" << l << ")";
    while(!after.empty()) {
      out << after.back();
      after.pop_back();
    }
  }
}

static void outputTrivialInequalityRemoval(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [parent] = getParents<1>(derived);

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
  outputParentWithSplits(out, parent);

  for(unsigned v : parentVars) {
    out << " ";
    outputVar(out, v, false, care);
  }
  for(Literal *l : litsParent) {
    out << " ";
    if(l->polarity() || !l->isEquality() || l->termArg(0) != l->termArg(1)) {
      out << l;
      continue;
    }
    out << "(p : Prf (";
    outputLiteral(out, l, care);
    out << ") => p (refl ";
    outputTerm(out, l->termArg(0), care);
    out << "))";
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
  outputParentWithSplits(out, parent);

  for(unsigned v : parentVars) {
    out << " ";
    outputTerm(out, subst.apply(TermList(v, false), 0), care);
  }
  unsigned lit;
  for(lit = 0; litsParent[lit] != selected; lit++) {
    out << " ";
    outputLiteralPtr(out, litsParent[lit], care, DoRobSubstitution0(subst));
  }

  out << " (p : Prf (";
  outputLiteral(out, subst.apply(selected, 0), care);
  out << ") => p (refl ";
  outputTerm(out, subst.apply(s, 0), care);
  out << "))";

  for(lit++; lit < litsParent.size(); lit++) {
    out << " ";
    outputLiteralPtr(out, litsParent[lit], care, DoRobSubstitution0(subst));
  }
}

static void outputEqualityResolutionWithDeletion(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);

  auto [parent] = getParents<1>(derived);
  const auto &er = env.proofExtra.get<Inferences::EqResWithDeletionExtra>(derived);

  // compute unifier for selected literal
  SimpleSubstitution subst;
  for(Literal *l : er.resolved) {
    ASS(l->isNegative())
    ASS(l->isEquality())

    TermList s = l->termArg(0), t = l->termArg(1);
    ASS(s.isVar() || t.isVar())
    bool lhs = s.isVar() && !t.containsSubterm(s);
    if(!lhs)
      std::swap(s, t);
    ASS(s.isVar() && !t.containsSubterm(s))
    ALWAYS(subst.bind(s.var(), t))
  }

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
    outputTerm(out, subst.apply(v), care);
  }

  for(Literal *l : litsParent) {
    out << " ";
    if(!er.resolved.count(l)) {
      outputLiteralPtr(out, l, care, DoSubstitution(subst));
      continue;
    }
    out << " (p : Prf (";
    outputLiteral(out, SubstHelper::apply(l, subst), care);
    out << ") => p (refl ";
    outputTerm(out, SubstHelper::apply(l->termArg(0), subst), care);
    out << "))";
  }
}


static void outputDuplicateLiteral(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);
  auto [parent] = getParents<1>(derived);

  // canonicalise order of literals in all clauses
  auto derivedLits = canonicalise(derived);
  auto litsParent = canonicalise(parent);

  // variables in numerical order
  auto vars = variables(derived);

  // bind variables present in the derived clause
  for(unsigned v : vars)
    out << " " << v << " : El iota => ";
  // bind literals in the derived clause
  for(Literal *l : derivedLits) {
    out << "" << l << " : (Prf ";
    outputLiteral(out, l, AlwaysCare {});
    out << " -> Prf false) => ";
  }

  // construct the proof term: refer to
  // "A Shallow Embedding of Resolution and Superposition Proofs into the λΠ-Calculus Modulo"
  // Guillaume Burel
  outputParentWithSplits(out, parent);

  for(unsigned v : vars) {
    out << " ";
    outputVar(out, v, false, AlwaysCare {});
  }
  for(Literal *l : litsParent)
    out << " " << l;
}

static void outputFactoring(std::ostream &out, Clause *derived) {
  outputDeductionPrefix(out, derived);

  auto [parent] = getParents<1>(derived);
  const auto &fact = env.proofExtra.get<Inferences::FactoringExtra>(derived);

  // compute unifier for selected literal
  RobSubstitution subst;
  Literal *selected = fact.selectedLiteral.selectedLiteral;
  Literal *other = fact.otherLiteral;
  ASS_EQ(selected->polarity(), other->polarity())
  ASS_EQ(selected->functor(), other->functor())
  ALWAYS(subst.unify(TermList(selected), 0, TermList(other), 0));

  // apply subst to all of the parent literals in the same order as Factoring does it
  // this will, I fervently hope, ensure that output variables are mapped in the same way
  for(Literal *l : parent->iterLits())
    if(l != other)
      subst.apply(l, 0);

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
  outputParentWithSplits(out, parent);

  for(unsigned v : parentVars) {
    out << " ";
    outputTerm(out, subst.apply(TermList(v, false), 0), care);
  }
  for(Literal *l : litsParent) {
    out << " ";
    outputLiteralPtr(out, l, care, DoRobSubstitution0(subst));
  }
}

static void outputAVATARDefinition(std::ostream &out, Unit *def) {
  Clause *component = env.proofExtra.get<SplitDefinitionExtra>(def).component;
  unsigned componentName = component->splits()->sval();
  SATLiteral split = Splitter::getLiteralFromName(componentName);
  out << "def sp" << split.var() << " : Prop :=";

  auto vars = variables(component);
  for(unsigned var : vars)
    out << " forall iota (" << var << " : El iota =>";
  auto literals = canonicalise(component);
  for(Literal *literal : literals) {
    out << " (imp (not ";
    outputLiteral(out, literal, AlwaysCare {});
    out << ")";
  }
  out << " false";
  for(unsigned i = 0; i < literals.size() + vars.size(); i++)
    out << ")";
}

static void outputAVATARComponent(std::ostream &out, Clause *component) {
  ASS(component->splits())
  ASS_EQ(component->splits()->size(), 1)
  outputDeductionPrefix(out, component);

  unsigned componentName = component->splits()->sval();
  SATLiteral split = Splitter::getLiteralFromName(componentName);
  auto literals = canonicalise(component);
  // variables in numerical order
  auto componentVars = variables(component);
  // bind variables present in the derived clause
  for(unsigned v : componentVars)
    out << " " << v << " : El iota => ";
  // bind literals in the derived clause
  for(Literal *l : literals) {
    out << "" << l << " : (Prf ";
    outputLiteral(out, l, AlwaysCare {});
    out << " -> Prf false) => ";
  }
  out << "nnsp" << split.var() << "(psp : Prf ";
  outputSplit(out, componentName);
  out << " => psp";

  for(unsigned v : componentVars)
    out << " " << v;
  for(Literal *l : literals)
    out << " " << l;
  out << ")";
}

static void outputAVATARSplitClause(std::ostream &out, Unit *derived) {
  UnitIterator parents = derived->getParents();
  ALWAYS(parents.hasNext())
  Clause *parent = parents.next()->asClause();
  SATClause *sat = env.proofExtra.get<SplitClauseExtra>(derived).clause;

  std::unordered_map<unsigned, Clause *> components;
  while(parents.hasNext()) {
    Unit *parent = parents.next();
    auto dex = env.proofExtra.get<SplitDefinitionExtra>(parent);
    ASS(dex.component->isComponent())
    unsigned component = dex.component->splits()->sval();
    components.insert({component, dex.component});
  }

  DHSet<SATLiteral> splitSet;
  if(!parent->noSplits()) {
    SplitSet &splits = *parent->splits();
    for(unsigned i = 0; i < splits.size(); i++)
      splitSet.insert(Splitter::getLiteralFromName(splits[i]));
  }
  std::vector<SATLiteral> existingSplits, newSplits;
  for(unsigned i = 0; i < sat->length(); i++) {
    if(splitSet.contains((*sat)[i]))
      existingSplits.push_back((*sat)[i]);
    else
      newSplits.push_back((*sat)[i]);
  }

  out << "def deduction" << derived->number() << " : Prf_av_clause";
  for(SATLiteral sl : existingSplits) {
    out << " (if ";
    outputSplit(out, sl);
  }
  out << " (acl (cl";
  for(SATLiteral sl : newSplits) {
    out << " (cons ";
    outputSplit(out, sl);
  }
  out << " ec";
  for(unsigned i = 0; i < sat->length(); i++)
    out << ")";
  out << ")) :=";

  Stack<LiteralStack> disjointLiterals;
  if(!Splitter::getComponents(parent, disjointLiterals)) {
    disjointLiterals.reset();
    LiteralStack component;
    for(Literal *l : parent->iterLits())
      component.push(l);
    disjointLiterals.push(std::move(component));
  }

  for(SATLiteral l : existingSplits) {
    out << " nnsp" << l.var() << " : ((Prf ";
    outputSplit(out, l);
    out << " -> Prf false) -> Prf false) =>";
  }
  for(SATLiteral l : newSplits) {
    out << " nnsp" << l.var() << " : (Prf ";
    outputSplit(out, l);
    out << " -> Prf false) =>";
  }

  unsigned closeParens = 0;
  decltype(disjointLiterals)::Iterator classes(disjointLiterals);
  // map variables in the parent to corresponding variables in the child splits
  SimpleSubstitution parent2SplitVars;
  while(classes.hasNext()) {
    LiteralStack klass = classes.next();
    Clause *found = nullptr;
    SimpleSubstitution subst;
    bool flip = false;
    for(auto [name, component] : components) {
      if(klass.size() == component->length()) {
        subst.reset();
        if(klass.size() == 1 && klass[0]->ground() && Literal::positiveLiteral(klass[0]) == Literal::positiveLiteral((*component)[0])) {
          found = component;
          flip = klass[0] != (*component)[0];
          break;
        }
        if(MLVariant::isVariant(klass.begin(), component, /*complementary=*/false, &subst)) {
          found = component;
          break;
        }
      }
    }
    ASS(found)

    SATLiteral split = Splitter::getLiteralFromName(found->splits()->sval());
    out << " nnsp" << split.var();
    for(unsigned foundVar : variables(found)) {
      out << " (" << subst.apply(foundVar).var() << " : El iota =>";
      ALWAYS(parent2SplitVars.bind(subst.apply(foundVar).var(), TermList(foundVar, false)))
      closeParens++;
    }
    for(Literal *l : canonicalise(found)) {
      if(flip)
        l = Literal::complementaryLiteral(l);
      out << " (" << SubstHelper::apply(l, subst) << " : (Prf ";
      outputLiteral(out, l, AlwaysCare {}, DoSubstitution(subst));
      out << " -> Prf false) =>";
      closeParens++;
    }
  }

  out << " deduction" << parent->number();
  // basically outputParentWithSplits() but NB comment below
  if(!parent->noSplits()) {
    SplitSet &splits = *parent->splits();
    for(unsigned i = 0; i < splits.size(); i++) {
      SATLiteral l = Splitter::getLiteralFromName(splits[i]);
      // if there is *already* a split in the parent (and it's negative),
      // need to double-negate it in order to match what the parent expects
      if(!l.polarity())
        out
          << "(nnnsp" << l.var() << " : (Prf (not sp" << l.var() << ") -> Prf false) =>"
          << "nnnsp" << l.var() << " nnsp" << l.var() << ")";
      else
        out << " nnsp" << l.var();
    }
  }

  for(unsigned parentVar : variables(parent))
    out << " " << parentVar; //parent2SplitVars.apply(parentVar).var();
  for(Literal *l : canonicalise(parent)) {
    out << ' ';
    Literal *lsubst = SubstHelper::apply(l, parent2SplitVars);
    bool needs_commute = l->isEquality() && (
         SubstHelper::apply((*l)[0], parent2SplitVars) != (*lsubst)[0]
      || SubstHelper::apply((*l)[1], parent2SplitVars) != (*lsubst)[1]);
    if(needs_commute) {
      out << "(comml ";
      outputTerm(out, (*l)[1], AlwaysCare {});
      out << ' ';
      outputTerm(out, (*l)[0], AlwaysCare {});
      out << ' ';
    }
    out << l;
    if(needs_commute)
      out << ")";
  }
  while(closeParens--)
    out << ')';
}

static void outputAVATARContradictionClause(std::ostream &out, Unit *contradiction) {
  auto [parent] = getParents<1>(contradiction);
  out << "def deduction" << contradiction->number() << " : ";
  outputClause(out, parent);
  out << " := deduction" << parent->number();
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
  case InferenceRule::DEFINITION_UNFOLDING:
    outputDefinitionUnfolding(out, u->asClause());
    break;
  case InferenceRule::RESOLUTION:
    outputResolution(out, u->asClause());
    break;
  case InferenceRule::REMOVE_DUPLICATE_LITERALS:
    outputDuplicateLiteral(out, u->asClause());
    break;
  case InferenceRule::FACTORING:
    outputFactoring(out, u->asClause());
    break;
  case InferenceRule::SUBSUMPTION_RESOLUTION:
    outputSubsumptionResolution(out, u->asClause());
    break;
  case InferenceRule::TRIVIAL_INEQUALITY_REMOVAL:
    outputTrivialInequalityRemoval(out, u->asClause());
    break;
  case InferenceRule::EQUALITY_RESOLUTION:
    outputEqualityResolution(out, u->asClause());
    break;
  case InferenceRule::EQUALITY_RESOLUTION_WITH_DELETION:
    outputEqualityResolutionWithDeletion(out, u->asClause());
    break;
  case InferenceRule::SUPERPOSITION:
    outputSuperposition(out, u->asClause());
    break;
  case InferenceRule::FORWARD_DEMODULATION:
  case InferenceRule::BACKWARD_DEMODULATION:
    outputDemodulation(out, u->asClause());
    break;
  case InferenceRule::AVATAR_DEFINITION:
    outputAVATARDefinition(out, u);
    break;
  case InferenceRule::AVATAR_COMPONENT:
    outputAVATARComponent(out, u->asClause());
    break;
  case InferenceRule::AVATAR_SPLIT_CLAUSE:
    outputAVATARSplitClause(out, u);
    break;
  case InferenceRule::AVATAR_CONTRADICTION_CLAUSE:
    outputAVATARContradictionClause(out, u);
    // TODO just the identity function
    break;
  case InferenceRule::AVATAR_REFUTATION:
    // TODO
    return;
  default:
    sorry(out, u);
  }
  out << ".\n";
}
}
}
