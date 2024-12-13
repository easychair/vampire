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
 * @file UnificationWithAbstraction.cpp
 * Defines class AbstractionOracle.
 *
 */

#include <vector>

#include "Lib/Backtrackable.hpp"
#include "Lib/Coproduct.hpp"
#include "Shell/Options.hpp"
#include "Lib/Environment.hpp"


#include "Forwards.hpp"
#include "Signature.hpp"
#include "Term.hpp"
#include "RobSubstitution.hpp"

#include "UnificationWithAbstraction.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/TermIterators.hpp"
#include "NumTraits.hpp"
#include "Kernel/TermIterators.hpp"
#include "Lib/Output.hpp"
#include "Debug/Tracer.hpp"
#define DEBUG(...) // DBG(__VA_ARGS__)
#define DEBUG_FINALIZE(LVL, ...) if (LVL < 0) DBG(__VA_ARGS__)
#define DEBUG_UNIFY(LVL, ...) if (LVL < 0) DBG(__VA_ARGS__)

namespace Kernel
{

Shell::Options::UnificationWithAbstraction AbstractionOracle::create()
{
  if (env.options->unificationWithAbstraction()!=Options::UnificationWithAbstraction::OFF) {
    return env.options->unificationWithAbstraction();
  } else if (env.options->functionExtensionality() == Options::FunctionExtensionality::ABSTRACTION && env.getMainProblem()->getProperty()->higherOrder()) { 
    return Options::UnificationWithAbstraction::FUNC_EXT;
  } else {
    return Options::UnificationWithAbstraction::OFF;
  }
}

Shell::Options::UnificationWithAbstraction AbstractionOracle::createOnlyHigherOrder()
{
  if (env.options->functionExtensionality() == Options::FunctionExtensionality::ABSTRACTION && env.getMainProblem()->getProperty()->higherOrder()) { 
    return Options::UnificationWithAbstraction::FUNC_EXT;
  } else {
    return Options::UnificationWithAbstraction::OFF;
  }
}

bool AbstractionOracle::isInterpreted(unsigned functor) const 
{
  auto f = env.signature->getFunction(functor);
  return f->interpreted() || f->termAlgebraCons();
}


class AcIter {
  unsigned _function;
  Recycled<Stack<TermSpec>> _todo;
  RobSubstitution const* _subs;
public:

  AcIter(unsigned function, TermSpec t, RobSubstitution const* subs) : _function(function), _todo(), _subs(subs) 
  { _todo->push(std::move(t)); }

  DECL_ELEMENT_TYPE(TermSpec);

  bool hasNext() const { return !_todo->isEmpty(); }

  TermSpec next() {
    ASS(!_todo->isEmpty());
    auto t = _todo->pop();
    auto* dt = &_subs->derefBound(t);
    while (dt->isTerm() && dt->functor() == _function) {
      ASS_EQ(dt->nTermArgs(), 2);
      _todo->push(dt->termArg(1));
      t = dt->termArg(0);
      dt = &_subs->derefBound(t);
    }
    return *dt;
  }
};

bool AbstractionOracle::canAbstract(AbstractingUnifier* au, TermSpec const& t1, TermSpec const& t2) const 
{
  if(!(t1.isTerm() && t2.isTerm())) return false;
  if(t1.isSort() || t2.isSort()) return false;

  bool t1Interp = isInterpreted(t1.functor());
  bool t2Interp = isInterpreted(t2.functor());
  bool bothNumbers = t1.isNumeral() && t2.isNumeral();

  switch(_mode) {
    case Shell::Options::UnificationWithAbstraction::INTERP_ONLY:
      return (t1Interp && t2Interp && !bothNumbers);
    case Shell::Options::UnificationWithAbstraction::ONE_INTERP:
      return !bothNumbers && (t1Interp || t2Interp);
      break;
    case Shell::Options::UnificationWithAbstraction::CONSTANT:
      return !bothNumbers && (t2Interp || t2Interp)
            && (t1Interp || t1.nTermArgs())
            && (t2Interp || t2.nTermArgs());
    case Shell::Options::UnificationWithAbstraction::ALL:
    case Shell::Options::UnificationWithAbstraction::GROUND:
      return true;
    case Shell::Options::UnificationWithAbstraction::OFF:
      return false;
    case Shell::Options::UnificationWithAbstraction::AC1: 
    case Shell::Options::UnificationWithAbstraction::AC2: 
    case Shell::Options::UnificationWithAbstraction::FUNC_EXT: 
      ASSERTION_VIOLATION
  }
  ASSERTION_VIOLATION;
}

Option<AbstractionOracle::AbstractionResult> funcExt(
    AbstractingUnifier* au, 
    TermSpec const& t1, TermSpec const& t2)
{
  ASS(t1.isTerm() || t2.isTerm())

  auto isApp = [](auto& t) { return env.signature->isAppFun(t.functor()); };
  if ( (t1.isTerm() && t1.isSort()) 
    || (t2.isTerm() && t2.isSort()) ) return Option<AbstractionOracle::AbstractionResult>();
  if (t1.isTerm() && t2.isTerm()) {
    if (isApp(t1) && isApp(t2)) {
      auto argSort1 = au->subs().derefBound(t1.typeArg(0));
      auto argSort2 = au->subs().derefBound(t2.typeArg(0));
      if (t1.isVar() || t2.isVar()
       || argSort1.isVar() || argSort2.isVar()
       || env.signature->isArrowCon(argSort1.functor())
       || env.signature->isArrowCon(argSort2.functor())
       || env.signature->isBoolCon(argSort1.functor())
       || env.signature->isBoolCon(argSort2.functor())
       ) {
        auto& arg1 = au->subs().derefBound(t1.termArg(1));
        auto& arg2 = au->subs().derefBound(t2.termArg(1));
        auto out = AbstractionOracle::EqualIf()
              .unify (UnificationConstraint(t1.typeArg(0), t2.typeArg(0), TermSpec(AtomicSort::superSort(), 0)),
                      UnificationConstraint(t1.typeArg(1), t2.typeArg(1), TermSpec(AtomicSort::superSort(), 0)),
                      UnificationConstraint(t1.termArg(0), t2.termArg(0), t1.termArgSort(0)));

        auto argsEq = UnificationConstraint(arg1, arg2, t1.termArgSort(1));
        auto res = some(AbstractionOracle::AbstractionResult(
              // if both are variables we don't want to introduce a constraint
              arg1.isVar() && arg2.isVar()
                ? std::move(out).unify(std::move(argsEq))
                : std::move(out).constr(std::move(argsEq)))) ;
        return res;
      }
    }
  }
  return Option<AbstractionOracle::AbstractionResult>();
}


Option<AbstractionOracle::AbstractionResult> AbstractionOracle::tryAbstract(AbstractingUnifier* au, TermSpec const& t1, TermSpec const& t2) const
{
  using Uwa = Shell::Options::UnificationWithAbstraction;
  ASS(_mode != Uwa::OFF)

  auto intSort = []() { return TermSpec(IntTraits::sort(), 0); };

  if (_mode == Uwa::AC1 || _mode == Uwa::AC2) {
      if ( (t1.isTerm() && t1.isSort())
        || (t2.isTerm() && t2.isSort())
        || !(t1.isTerm() && theory->isInterpretedFunction(t1.functor(), IntTraits::addI))
        || !(t2.isTerm() && theory->isInterpretedFunction(t2.functor(), IntTraits::addI))
        ) {
        return Option<AbstractionResult>();
      }
      auto a1 = iterTraits(AcIter(IntTraits::addF(), t1, &au->subs())).template collect<Stack>();
      auto a2 = iterTraits(AcIter(IntTraits::addF(), t2, &au->subs())).template collect<Stack>();
      auto cmp = [&](TermSpec const& lhs, TermSpec const& rhs) { return TermSpec::compare(lhs, rhs, [&](auto& t) -> TermSpec const& { return au->subs().derefBound(t); }); };
      auto less = [&](TermSpec const& lhs, TermSpec const& rhs) { return cmp(lhs, rhs) < 0; };
      a1.sort(less);
      a2.sort(less);

      Recycled<Stack<TermSpec>> diff1_;
      Recycled<Stack<TermSpec>> diff2_;
      auto& diff1 = *diff1_;
      auto& diff2 = *diff2_;
      diff1.loadFromIterator(iterSortedDiff(arrayIter(a1), arrayIter(a2), cmp));
      diff2.loadFromIterator(iterSortedDiff(arrayIter(a2), arrayIter(a1), cmp));
      auto sum = [&](auto& diff) {
          return arrayIter(diff)
            .map([](TermSpec& t) -> TermSpec { return t; })
            .fold([&](TermSpec l, TermSpec r) -> TermSpec
              { return au->subs().createTerm(IntTraits::addF(), l, r); })
            .unwrap(); };
      auto diffConstr = [&]() 
      { return UnificationConstraint(sum(diff1), sum(diff2), intSort()); };

      auto functors = [](auto& diff) 
      { return arrayIter(diff).map([](auto& f) { return f.functor(); }); };

      if (diff1.size() == 0 && diff2.size() == 0) {
        return some(AbstractionResult(EqualIf()));

      } else if (( diff1.size() == 0 && diff2.size() != 0 )
              || ( diff2.size() == 0 && diff1.size() != 0 ) ) {
        return some(AbstractionResult(NeverEqual{}));

      } else if (_mode == Uwa::AC2 && diff1.size() == 1 && diff1[0].isVar()) {
        return some(AbstractionResult(EqualIf().unify(UnificationConstraint(std::move(diff1[0]), sum(diff2), intSort()))));

      } else if (_mode == Uwa::AC2 && diff2.size() == 1 && diff2[0].isVar()) {
        return some(AbstractionResult(EqualIf().unify(UnificationConstraint(std::move(diff2[0]), sum(diff1), intSort()))));

      } else if (concatIters(arrayIter(diff1), arrayIter(diff2)).any([](auto& x) { return x.isVar(); })) {
        return some(AbstractionResult(EqualIf().constr(diffConstr())));

      } else if (iterSortedDiff(functors(diff1), functors(diff2)).hasNext()
              || iterSortedDiff(functors(diff2), functors(diff1)).hasNext()) {
        return some(AbstractionResult(NeverEqual{}));

      } else {
        return some(AbstractionResult(EqualIf().constr(diffConstr())));
      }

  } else if (_mode == Shell::Options::UnificationWithAbstraction::FUNC_EXT) {
    return funcExt(au, t1, t2);

  } else {
    auto abs = canAbstract(au, t1, t2);
    DEBUG("canAbstract(", t1, ",", t2, ") = ", abs);
    return someIf(abs, [&](){
        return AbstractionResult(EqualIf().constr(UnificationConstraint(t1, t2, 
                t1.isTerm() ? TermSpec(SortHelper::getResultSort(t1.term.term()), t1.index)
                            : TermSpec(SortHelper::getResultSort(t2.term.term()), t2.index)
                )));
    });
  }
}

void UnificationConstraintStack::add(UnificationConstraint c, Option<BacktrackData&> bd)
{ 
  DEBUG("introduced constraint: ", c)
  _cont.push(c);
  if (bd) {
    bd->addClosure([this]() { _cont.pop(); });
  }
}


UnificationConstraint UnificationConstraintStack::pop(Option<BacktrackData&> bd)
{ 
  auto old = _cont.pop();
  if (bd) {
    bd->addClosure([this, old]() mutable { _cont.push(std::move(old)); });
  }
  return old;
}

Recycled<Stack<Literal*>> UnificationConstraintStack::literals(RobSubstitution& s)
{ 
  Recycled<Stack<Literal*>> out;
  out->reserve(_cont.size());
  out->loadFromIterator(literalIter(s));
  return out;
}


Option<Literal*> UnificationConstraint::toLiteral(RobSubstitution& s)
{ 
  auto t1 = _t1.toTerm(s);
  auto t2 = _t2.toTerm(s);
  auto sort = _sort.toTerm(s);
  return t1 == t2 
    ? Option<Literal*>()
    : Option<Literal*>(Literal::createEquality(false, t1, t2, sort));
}


}

bool AbstractingUnifier::fixedPointIteration()
{
  TIME_TRACE("uwa fixed point")
  Recycled<Stack<UnificationConstraint>> todo;
  while (!constr().isEmpty()) { 
    todo->push(constr().pop(bd()));
  }

  DEBUG_FINALIZE(1, "finalizing: ", todo)
  while (!todo->isEmpty()) {
    auto c = todo->pop();
    DEBUG_FINALIZE(2, "popped: ", c);
    bool progress;
    auto res = unify(c.lhs(), c.rhs(), progress);
    if (!res) {
      DEBUG_FINALIZE(1, "finalizing failed");
      return false;
    } else if (progress) {
      while (!constr().isEmpty()) { 
        todo->push(constr().pop(bd()));
      }
    } else {
      // no progress. we keep the constraints
    }
  }
  DEBUG_FINALIZE(1, "finalizing successful: ", *this);
  return true;
}

Option<Recycled<Stack<unsigned>>> AbstractingUnifier::unifiableSymbols(SymbolId fid)
{
  auto f = fid.functor;
  if (fid.kind == TermKind::SORT)
    // we don't perform UWA on sorts
    return some(recycledStack(f));

  ASS(fid.kind == TermKind::TERM) // not implemented for literals

  auto anything = []() -> Option<Recycled<Stack<unsigned>>> { return {}; };
  auto nothing  = []() -> Option<Recycled<Stack<unsigned>>> { return some(recycledStack<unsigned>()); };
  switch (_uwa._mode) {
    case Options::UnificationWithAbstraction::OFF: return some(recycledStack(f));
    case Options::UnificationWithAbstraction::INTERP_ONLY: return theory->isInterpretedFunction(f) ? anything() : some(recycledStack(f));
    case Options::UnificationWithAbstraction::ONE_INTERP: return anything();
    case Options::UnificationWithAbstraction::CONSTANT: return theory->isInterpretedConstant(f) ? anything() : nothing();
    case Options::UnificationWithAbstraction::ALL: return anything();
    case Options::UnificationWithAbstraction::GROUND: anything();
    case Options::UnificationWithAbstraction::FUNC_EXT: anything();
    case Options::UnificationWithAbstraction::AC1: return some(recycledStack(f));
    case Options::UnificationWithAbstraction::AC2: return some(recycledStack(f));
  }
  ASSERTION_VIOLATION
}

bool AbstractingUnifier::unify(TermList term1, unsigned bank1, TermList term2, unsigned bank2)
{
  if (_uwa._mode == Shell::Options::UnificationWithAbstraction::OFF) 
    return _subs->unify(term1, bank1, term2, bank2);

  bool progress;
  return unify(TermSpec(term1, bank1), TermSpec(term2, bank2), progress);
}

bool AbstractingUnifier::unify(TermSpec t1, TermSpec t2, bool& progress)
{
  TIME_TRACE("unification with abstraction")
  ASS_NEQ(_uwa._mode, Shell::Options::UnificationWithAbstraction::OFF) 
  DEBUG_UNIFY(1, *this, ".unify(", t1, ",", t2, ")")
  progress = false;

  if(t1 == t2) {
    progress = true;
    return true;
  }

  auto impl = [&]() -> bool {

    Recycled<Stack<std::pair<TermSpec, TermSpec>>> toDo;
    toDo->push(std::make_pair(t1, t2));
    
    // Save encountered unification pairs to avoid
    // recomputing their unification
    Recycled<DHSet<std::pair<TermSpec,TermSpec>>> encountered;

    Option<AbstractionOracle::AbstractionResult> absRes;
    auto doAbstract = [&](auto& l, auto& r) -> bool
    { 
      absRes = _uwa.tryAbstract(this, l, r);
      if (absRes) {
        DEBUG_UNIFY(2, "abstraction result: ", absRes)
      }
      return absRes.isSome();
    };

    auto pushTodo = [&](auto pair) {
        // we unify each subterm pair at most once, to avoid worst-case exponential runtimes
        // in order to safe memory we do ot do this for variables.
        // (Note by joe:  didn't make this decision, but just keeping the implemenntation 
        // working as before. i.e. as described in the paper "Comparing Unification 
        // Algorithms in First-Order Theorem Proving", by Krystof and Andrei)
        // TODO restore this branch?
        // if (pair.first.isVar() && isUnbound(pair.first.varSpec()) &&
        //     pair.second.isVar() && isUnbound(pair.second.varSpec())) {
        //   toDo.push(pair);
        // } else 
        if (!encountered->find(pair)) {
          encountered->insert(pair);
          toDo->push(std::move(pair));
        }
    };

    auto occurs = [this](auto& var, auto& term) {
      Recycled<Stack<TermSpec>> todo;
      todo->push(term);
      while (todo->isNonEmpty()) {
        auto t = todo->pop();
        auto& dt = subs().derefBound(t);
        if (dt.isVar()) {
          if (dt == var) {
            return true;
          }
        } else {
          todo->loadFromIterator(dt.allArgs());
        }
      }
      return false;
    };


    while (toDo->isNonEmpty()) {
      auto cur = toDo->pop();
      auto& dt1 = subs().derefBound(cur.first);
      auto& dt2 = subs().derefBound(cur.second);
      DEBUG_UNIFY(2, "popped: ", dt1, " = ", dt2)
      if (dt1.deepEqCheck(dt2)) {
        progress = true;

      } else if(dt1.isVar() && !occurs(dt1, dt2)) {
        progress = true;
        subs().bind(dt1.varSpec(), dt2);

      } else if(dt2.isVar() && !occurs(dt2, dt1)) {
        progress = true;
        subs().bind(dt2.varSpec(), dt1);

      } else if(doAbstract(dt1, dt2)) {

        ASS(absRes);
        if (absRes->is<AbstractionOracle::NeverEqual>()) {
          return false;

        } else {
          ASS(absRes->is<AbstractionOracle::EqualIf>())
          auto& conditions = absRes->unwrap<AbstractionOracle::EqualIf>();
          auto deepEqCheck = [](UnificationConstraint& c, TermSpec const& lhs, TermSpec const& rhs) 
          { 
            return (c.lhs().deepEqCheck(lhs) && c.rhs().deepEqCheck(rhs)) 
                || (c.lhs().deepEqCheck(rhs) && c.rhs().deepEqCheck(lhs)); };
          if (progress 
              || conditions.constr().size() != 1 
              || conditions.unify().size() != 0
              || !deepEqCheck(conditions.constr()[0], t1, t2)
              ) {
            progress = true;
          }
          for (auto& x : conditions.unify()) {
            auto pair = std::make_pair(x.lhs(), x.rhs());
            ASS_NEQ(pair, cur)
            pushTodo(pair);
          }
          for (auto& x: conditions.constr()) {
            _constr->add(std::move(x), bd());
          }
        }

      } else if(dt1.isTerm() && dt2.isTerm() && dt1.functor() == dt2.functor()) {
        
        for (auto p : dt1.allArgs().zip(dt2.allArgs())) {
          pushTodo(p);
        }

      } else {
        return false;
      }

    }
    return true;
  };

  BacktrackData localBD;
  bdRecord(localBD);
  bool success = impl();
  bdDone();

  if(!success) {
    localBD.backtrack();
  } else {
    if(subs().bdIsRecording()) {
      subs().bdCommit(localBD);
    }
    localBD.drop();
  }

  DEBUG_UNIFY(1, *this, " (", success ? "success" : "fail", ")")
  return success;
}
