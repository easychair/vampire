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
 * @file TestUtils.cpp
 * Implements class TestUtils.
 */

#include <cstdarg>


#include "Lib/List.hpp"
#include "Kernel/TermIterators.hpp"

#include "Kernel/Formula.hpp"
#include "Kernel/Unit.hpp"

#include "SAT/SATClause.hpp"
#include "SAT/SATLiteral.hpp"

#include "TestUtils.hpp"

namespace Test {

using namespace Kernel;
using namespace Shell;

SAT::SATClause* TestUtils::buildSATClause(unsigned len,...)
{
  using namespace SAT;

  va_list litArgs;
  va_start(litArgs,len);

  static SATLiteralStack lits;
  lits.reset();
  for(unsigned i=0; i<len; ++i) {
    int litVal = va_arg(litArgs, int);
    bool pos = litVal>0;
    unsigned var = pos ? litVal : -litVal;
    lits.push(SATLiteral(var, pos));
  }

  va_end(litArgs);

  return SATClause::fromStack(lits);
}

std::ostream& printOp(std::ostream& out, const Term* t, const char* op, bool par = true) {
  auto l = *t->nthArgument(0);
  auto r = *t->nthArgument(1);
  if (par) out << "(";
  out << pretty(l) << " " << op << " " << pretty(r);
  if (par) out << ")";
  return out;
}

template<>
std::ostream& Pretty<Kernel::TermList>::prettyPrint(std::ostream& out) const
{ return out << _self; }

template<>
std::ostream& Pretty<Literal*>::prettyPrint(std::ostream& out) const
{ return out << pretty(*_self); }


template<>
std::ostream& Pretty<Clause>::prettyPrint(std::ostream& out) const
{ 
  if (_self.isEmpty()) {
    return out << "_|_";
  }
  auto iter = _self.iterLits();
    out << "{ ";
  if (iter.hasNext()) {
    out << pretty(*iter.next());
    while(iter.hasNext()) {
      out << " | " << pretty(*iter.next());
    }
  } else {
    // out << "[]";
  }
    out << " }";
  return out;
}

template<>
std::ostream& Pretty<Literal>::prettyPrint(std::ostream& out) const
{
  const Literal& lit = _self;
  auto print = [&]() -> std::ostream& {

    auto func = lit.functor();
    if(theory->isInterpretedPredicate(func)) {
      switch(theory->interpretPredicate(func)) {
#define NUM_CASE(oper) \
        case Kernel::Theory::INT_  ## oper: \
        case Kernel::Theory::REAL_ ## oper: \
        case Kernel::Theory::RAT_  ## oper

        NUM_CASE(GREATER):
          return printOp(out, &lit, ">", /* par = */ false);
        NUM_CASE(GREATER_EQUAL):
          return printOp(out, &lit, ">=", /* par = */ false);
        NUM_CASE(LESS_EQUAL):
          return printOp(out, &lit, "<=", /* par = */ false);
        case Kernel::Theory::EQUAL:
          return printOp(out, &lit, "=", /* par = */ false);
        default: 
        {
        }
#undef NUM_CASE
      }
    }
    Signature::Symbol* sym = env.signature->getPredicate(func);
    out << sym->name();
    if (sym->arity() > 0) {
      out << "(" << pretty(*lit.nthArgument(0));
      for (unsigned i = 1; i < sym->arity(); i++) {
        out << ", " << pretty(*lit.nthArgument(i));
      }
      out << ")";
    }
    return out;
  };


  if (!lit.polarity()) {
    out << "~(";
  }
  print();
  if (!lit.polarity()) {
    out << ")";
  }
  return out;
}


bool TestUtils::isAC(Term* t) 
{
  auto f = t->functor();
  if (t->isLiteral()) {
    return theory->isInterpretedPredicate(f) && isAC(theory->interpretPredicate(f));
  } else if (t->isSort()) {
    return false;
  } else {
    return theory->isInterpretedFunction(f) && isAC(theory->interpretFunction(f));
  }
}
bool TestUtils::isAC(Theory::Interpretation i) 
{
  switch (i) {
#define NUM_CASE(oper) \
    case Kernel::Theory::INT_  ## oper: \
    case Kernel::Theory::REAL_ ## oper: \
    case Kernel::Theory::RAT_  ## oper

    NUM_CASE(PLUS):
    NUM_CASE(MULTIPLY):
    case Kernel::Theory::EQUAL:
      return true;
    default: 
      return false;

#undef NUM_CASE
  }
}

bool TestUtils::eqModAC(const Kernel::Clause* lhs, const Kernel::Clause* rhs)
{ return permEq(*lhs, *rhs, [](Literal* l, Literal* r) -> bool { return TestUtils::eqModAC(l, r); }); }

bool TestUtils::eqModAC(Kernel::Literal* lhs, Kernel::Literal* rhs)
{ return lhs->isPositive() == rhs->isPositive() 
      && TestUtils::eqModAC(TermList(lhs), TermList(rhs)); }


void __collect(unsigned functor, Term* t, Stack<TermList>& out) {
  ASS_EQ(t->functor(), functor);
  for (unsigned i = 0; i < t->arity(); i++) {
    auto trm = t->nthArgument(i);
    if (trm->isVar()) {
      out.push(*trm);
    } else {
      ASS(trm->isTerm());
      if (trm->term()->functor() == functor) {
        __collect(functor, trm->term(), out);
      } else {
        out.push(*trm);
      }
    }
  }
}

Stack<TermList> collect(unsigned functor, Term* t) {
  Stack<TermList> out;
  __collect(functor, t, out);
  return out;
}


template<class Comparisons>
bool TestUtils::eqModAC_(TermList lhs, TermList rhs, Comparisons comp)
{
  if (lhs.isVar() && rhs.isVar()) {
    return comp.var(lhs.var(), rhs.var());
  } else if (lhs.isTerm() && rhs.isTerm()) {
    auto& l = *lhs.term();
    auto& r = *rhs.term();
    if ( l.functor() != r.functor() ) return false;
    auto fun = l.functor();
    if (isAC(&l)) {
      Stack<TermList> lstack = collect(fun, &l);
      Stack<TermList> rstack = collect(fun, &r);
      return permEq(lstack, rstack, [&](TermList l, TermList r) -> bool {
            return comp.subterm(l, r);
      });
    } else {
      for (unsigned i = 0; i < l.arity(); i++) {
        if (!comp.subterm(*l.nthArgument(i), *r.nthArgument(i))) {
          return false;
        }
      }
      return true;
    }
  } else {
    return false;
  }
}

bool TestUtils::eqModAC(TermList lhs, TermList rhs) 
{
  struct Comparisons {
    bool var(unsigned lhs, unsigned rhs) const 
    { return lhs == rhs; }

    bool subterm(TermList lhs, TermList rhs) const 
    { return eqModAC(lhs,rhs); }
  };
  Comparisons c {};

  return eqModAC_(lhs, rhs, c);
}

struct AcRectComp {
  Stack<TermList> const& vl;
  Stack<TermList> const& vr;
  DArray<unsigned> const& perm;

  bool var(unsigned lhs, unsigned rhs) const 
  { 
    auto l = iterTraits(getRangeIterator<unsigned>(0, vl.size()))
              .find([&](auto i) { return vl[i].var() == lhs; })
              .unwrap();
    auto r = iterTraits(getRangeIterator<unsigned>(0, vr.size()))
              .find([&](auto i) { return vr[perm[i]].var() == rhs; })
              .unwrap();
    return l == r;
  }

  bool subterm(TermList lhs, TermList rhs) const 
  { return TestUtils::eqModAC_(lhs, rhs, *this); }
};

bool TestUtils::eqModACRect(Kernel::TermList lhs, Kernel::TermList rhs)
{ 
  auto vl = iterTraits(vi(new VariableIterator(lhs))).collect<Stack>();
  vl.sort();
  vl.dedup();
  auto vr = iterTraits(vi(new VariableIterator(lhs))).collect<Stack>();
  vr.sort();
  vr.dedup();

  if (vl.size() != vr.size()) return false;

  return anyPerm(vl.size(), [&](auto& perm) {
    AcRectComp c {vl, vr, perm};
    return eqModAC_(lhs, rhs, c);
  });
}


template<class Lits>
bool TestUtils::_eqModACRect(Lits const& lhs, Lits const& rhs)
{ 
  auto vars = [](Lits const& lits) {
    auto vs = iterTraits(getRangeIterator(0u, (unsigned) lits.size()))
      .map([&](auto i) { return lits[i]; })
      .flatMap([](Literal* lit) { return vi(new VariableIterator(lit)); })
      .template collect<Stack>();
    vs.sort();
    vs.dedup();
    return vs;
  };
  auto vl = vars(lhs);
  auto vr = vars(rhs);

  if (vl.size() != vr.size()) return false;

  return anyPerm(vl.size(), [&](auto& perm) {
     AcRectComp c {vl, vr, perm};
      return permEq(lhs, rhs, [&](Literal* l, Literal* r) {
        return l->isPositive() == r->isPositive() 
            && eqModAC_(TermList(l), TermList(r), c); 
      });
  });
}


bool TestUtils::eqModACRect(const Clause* lhs, const Clause* rhs)
{ return _eqModACRect(*lhs, *rhs); }

bool TestUtils::eqModACRect(Stack<Literal*> const& lhs, Stack<Literal*> const& rhs)
{ return _eqModACRect(lhs, rhs); }


bool TestUtils::eqModACRect(Literal* lhs, Literal* rhs)
{ return lhs->isPositive() == rhs->isPositive() 
      && TestUtils::eqModACRect(TermList(lhs), TermList(rhs)); }


} // namespace Test
