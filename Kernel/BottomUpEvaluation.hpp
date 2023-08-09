/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __LIB__BOTTOM_UP_EVALUATION_HPP__
#define __LIB__BOTTOM_UP_EVALUATION_HPP__

#define DEBUG(...) // DBG(__VA_ARGS__)

/**
 * @file Kernel/BottomUpEvaluation.hpp
 *
 * This file contains mainly the function Lib::evaluateBottomUp, that can be used to evaluate an arbitrary
 * acyclic graph structure bottom up. It uses iteration instead of recursion and can be equipped with memoization structures
 * found in Lib::Memo.
 *
 * \see UnitTests/tBottomUpEvaluation.cpp
 * \see Lib::evaluateBottomUp
 */

#include "Lib/Stack.hpp"
#include "Lib/Recycled.hpp"
#include "Lib/Option.hpp"
#include "Lib/TypeList.hpp"
#include <utility>

namespace Lib {

namespace Memo {

  /** a mocked memoization that does not store any results */
  template<class Arg, class Result>
  struct None
  {
    Option<Result> get(Arg const&)
    { return Option<Result>(); }

    template<class Init> Result getOrInit(Arg const& orig, Init init)
    { return init(); }
  };

  /** a memoization realized as a hashmap */
  template<class Arg, class Result, class Hash = DefaultHash>
  class Hashed
  {
    Map<Arg, Result, Hash> _memo;

  public:
    Hashed() : _memo(decltype(_memo)()) {}

    template<class Init> Result getOrInit(Arg const& orig, Init init)
    { return _memo.getOrInit(Arg(orig), init); }

    Option<Result> get(const Arg& orig)
    { return _memo.tryGet(orig).toOwned(); }
  };

} // namespace Memo

/**
 * An iterator over the children of a node in a Directed Acyclic Graph (DAG).
 * The DAG is a structure to be evaluate bottom up using evaluateBottomUp.
 * A template-specialization of this struct with the methods like in this template
 * must be implemented for a class in order to bottom-up evaluatable.
 */
template<class A>
struct BottomUpChildIter
{
  /** constructs an iterator over the children of the current node */
  BottomUpChildIter(A a);

  /** returns the node this iterator was constructed with */
  A self();

  /** returns the next child of the node this this object was constructed with */
  A next();

  /** returns the next child of the current node in the structure to be traversed */
  bool hasNext();

  /** returns how many children this node has */
  unsigned nChildren();
};

template<class A> BottomUpChildIter<A> bottomUpChildIter(A a)
{ return BottomUpChildIter<A>(a); }

/**
 * Evaluates a term-like datastructure (i.e.: a Directed Acyclic Graph (DAG)), without using recursion.
 *
 * Optionly a memoization method (i.e. a class from Kernel::Memo) can be specified. The memo can be a static,
 * variable, in order to keep cached results for multiple runs of the funcion.
 *
 * The term-ish structure is evaluated according to the structure EvalFn. It is expected to have the following structure:
 * class EvalFn {
 *    using Arg    = ...; // <- the term-ish structure that will be evaluated.
 *                        //    A specialization template<> BottomUpChildIter<Arg> must exist
 *    using Result = ...; // <- the type the structure will be evaluated to
 *
 *    // The actual evaluation function. It will be called once for each node in the directed acyclic graph, together with
 *    // the already recursively evaluated children.
 *    Result operator()(Arg const& orig, Result* evaluatedChildren);
 * }
 *
 * The term to be evaluated will be traversed using a BottomUpChildIter<Arg>.
 */
template<class EvalFn, class Memo, class... Context>
typename EvalFn::Result evaluateBottomUpWithMemo(typename EvalFn::Arg const& term, EvalFn evaluateStep, Memo& memo, Context... cs)
{
  CALL("evaluateBottomUp(...)")
  using Result = typename EvalFn::Result;
  using Arg    = typename EvalFn::Arg;

  static_assert(std::is_same<ResultOf<EvalFn, Arg, Result*>, Result>::value, "evaluation function must have signature `Result eval(Arg term, Result* evaluatedArgs)`");


  /* recursion state. Contains a stack of items that are being recursed on. */
  Recycled<Stack<BottomUpChildIter<Arg>>> recState;
  Recycled<Stack<Result>> recResults;

  recState->push(BottomUpChildIter<Arg>(term, cs...));

  while (!recState->isEmpty()) {
    if (recState->top().hasNext(cs...)) {
      Arg t = recState->top().next(cs...);

      Option<Result> cached = memo.get(t);
      if (cached.isSome()) {
        recResults->push(std::move(cached).unwrap());
      } else {
        recState->push(BottomUpChildIter<Arg>(t, cs...));
      }

    } else {

      BottomUpChildIter<Arg> orig = recState->pop();
      Result eval = memo.getOrInit(orig.self(), [&](){
            CALL("evaluateBottomUp(..)::closure@1")
            Result* argLst = NULL;
            if (orig.nChildren(cs...) != 0) {
              ASS_GE(recResults->size(), orig.nChildren(cs...));
              argLst = static_cast<Result*>(&((*recResults)[recResults->size() - orig.nChildren(cs...)]));
            }
            return evaluateStep(orig.self(), argLst);
          });

      DEBUG("evaluated: ", orig.self(), " -> ", eval);
      recResults->pop(orig.nChildren(cs...));
      recResults->push(std::move(eval));
    }
  }
  ASS(recState->isEmpty())


  ASS(recResults->size() == 1);
  auto result = recResults->pop();
  DEBUG("eval result: ", term, " -> ", result);
  return result;
}

/** convenience wrapper for using evaluateBottomUp without a memo. */
template<class EvalFn, class... Context>
typename EvalFn::Result evaluateBottomUp(typename EvalFn::Arg const& term, EvalFn evaluateStep, Context... cs) 
{
  using namespace Memo;
  auto memo = None<typename EvalFn::Arg, typename EvalFn::Result>();
  return evaluateBottomUpWithMemo(term, evaluateStep, memo, std::move(cs)...);
}
}

#include "Kernel/Term.hpp"

namespace Lib {
// specialisation for TermList
// iterate up through TermLists, ignoring sort arguments
template<>
struct BottomUpChildIter<Kernel::TermList>
{
  Kernel::TermList _self;
  unsigned _idx;

  BottomUpChildIter(Kernel::TermList self) : _self(self), _idx(0)
  { }

  Kernel::TermList next()
  {
    ASS(hasNext());
    return _self.term()->termArg(_idx++);
  }

  bool hasNext() const
  { return _self.isTerm() && _idx < _self.term()->numTermArguments(); }

  unsigned nChildren() const
  { return _self.isVar() ? 0 : _self.term()->numTermArguments(); }

  Kernel::TermList self() const
  { return _self; }
};
}

#include "TypedTermList.hpp"

namespace Lib {
// specialisation for TypedTermList
template<>
struct BottomUpChildIter<Kernel::TypedTermList>
{
  Kernel::TypedTermList _self;
  unsigned      _idx;

  BottomUpChildIter(Kernel::TypedTermList self) : _self(self), _idx(0)
  {}

  Kernel::TypedTermList next(int);
  Kernel::TypedTermList next()
  {
    ASS(hasNext());
    auto cur = self().term();
    auto next = cur->termArg(_idx);
    auto sort = Kernel::SortHelper::getTermArgSort(cur, _idx);
    ASS_NEQ(sort, Kernel::AtomicSort::superSort())
    _idx++;
    return Kernel::TypedTermList(next, sort);
  }

  bool hasNext() const
  { return _self.isTerm() && _idx < _self.term()->numTermArguments(); }

  unsigned nChildren() const
  { return _self.isVar() ? 0 : _self.term()->numTermArguments(); }

  Kernel::TypedTermList self() const
  { return _self; }
};

template<class EvalFn, class Memo>
Kernel::Literal* evaluateLiteralBottomUp(Kernel::Literal* const& lit, EvalFn evaluateStep, Memo& memo)
{
  using namespace Kernel;
  Recycled<Stack<TermList>> args;
  for (unsigned i = 0; i < lit->arity(); i++) {
    args->push(evaluateBottomUpWithMemo(TypedTermList(*lit->nthArgument(i), SortHelper::getArgSort(lit, i)), evaluateStep, memo));
  }
  return Literal::create(lit, args->begin());
}


template<class EvalFn>
Kernel::Literal* evaluateLiteralBottomUp(Kernel::Literal* const& lit, EvalFn evaluateStep)
{
  using namespace Memo;
  auto memo = None<typename EvalFn::Arg, typename EvalFn::Result>();
  return evaluateLiteralBottomUp(lit, evaluateStep, memo);
}
}

#include "Polynomial.hpp"

namespace Lib {
// specialisation for PolyNf
template<>
struct BottomUpChildIter<Kernel::PolyNf>
{
  struct PolynomialBottomUpChildIter
  {
    Kernel::AnyPoly _self;
    unsigned _idx1;
    unsigned _idx2;
    unsigned _nChildren;

    PolynomialBottomUpChildIter(Kernel::AnyPoly self) : _self(self), _idx1(0), _idx2(0), _nChildren(0)
    {
      while (_idx1 < _self.nSummands() && _self.nFactors(_idx1) == 0) {
        _idx1++;
      }
      for (unsigned i = 0; i < _self.nSummands(); i++) {
        _nChildren += self.nFactors(i);
      }
    }

    bool hasNext() const
    { return _idx1 < _self.nSummands(); }

    Kernel::PolyNf next()
    {
      auto out = _self.termAt(_idx1, _idx2++);
      if (_idx2 >= _self.nFactors(_idx1)) {
        _idx1++;
        while (_idx1 < _self.nSummands() && _self.nFactors(_idx1) == 0) {
          _idx1++;
        }
        _idx2 = 0;
      }
      return out;
    }

    unsigned nChildren() const
    { return _nChildren; }

    friend ostream& operator<<(ostream& out, PolynomialBottomUpChildIter const& self)
    { return out << self._self << "@(" << self._idx1 << ", " << self._idx2 << ")"; }
  };

  struct FuncTermBottomUpChildIter
  {

    Perfect<Kernel::FuncTerm> _self;
    unsigned _idx;

    FuncTermBottomUpChildIter(Perfect<Kernel::FuncTerm> self) : _self(self), _idx(0) {}

    bool hasNext() const
    { return _idx < _self->numTermArguments(); }

    Kernel::PolyNf next()
    { return _self->arg(_idx++); }

    unsigned nChildren() const
    { return _self->numTermArguments(); }

    friend ostream& operator<<(ostream& out, FuncTermBottomUpChildIter const& self)
    { return out << self._self << "@" << self._idx; }
  };


  struct VariableBottomUpChildIter
  {
    Kernel::Variable _self;
    VariableBottomUpChildIter(Kernel::Variable self) : _self(self) {}

    bool hasNext() const
    { return false; }

    Kernel::PolyNf next()
    { ASSERTION_VIOLATION }

    unsigned nChildren() const
    { return 0; }

    friend ostream& operator<<(ostream& out, VariableBottomUpChildIter const& self)
    { return out << self._self; }
  };

  using Inner = Coproduct<FuncTermBottomUpChildIter, VariableBottomUpChildIter, PolynomialBottomUpChildIter>;
  Inner _self;

  BottomUpChildIter(Kernel::PolyNf self) : _self(self.match(
        [&](Perfect<Kernel::FuncTerm> self) { return Inner(FuncTermBottomUpChildIter( self ));            },
        [&](Kernel::Variable                  self) { return Inner(VariableBottomUpChildIter( self ));            },
        [&](Kernel::AnyPoly           self) { return Inner(PolynomialBottomUpChildIter(std::move(self))); }
      ))
  {}

  Kernel::PolyNf next()
  { ALWAYS(hasNext()); return _self.apply([](auto& x) -> Kernel::PolyNf { return x.next(); }); }

  bool hasNext() const
  { return _self.apply([](auto& x) { return x.hasNext(); }); }

  unsigned nChildren() const
  { return _self.apply([](auto& x) { return x.nChildren(); }); }

  Kernel::PolyNf self() const
  { return _self.apply([](auto& x) { return Kernel::PolyNf(x._self); }); }

  friend ostream& operator<<(ostream& out, BottomUpChildIter const& self)
  { return out << self._self; }
};

namespace TL = Lib::TypeList;

template<unsigned I, class Indexed>
struct MapTupleElem;

template<unsigned I, unsigned J, class A>
struct MapTupleElem<I, TL::Indexed<J, A>>
{
  template<class Tup, class F>
  inline static auto apply(Tup& bs, F& f) -> A
  { return std::get<J>(bs); }
};

template<unsigned I, class A>
struct MapTupleElem<I, TL::Indexed<I, A>>
{
  template<class Tup, class F>
  inline static auto apply(Tup& bs, F& f) -> decltype(auto) 
  { return std::move(f)(std::get<I>(bs)); }
};



template<unsigned N, class F, class Tup, class... Indexed>
auto __mapTupleElem(Tup tup, F f, TL::List<Indexed...>) -> decltype(auto) {
  return std::tuple<
    decltype(MapTupleElem<N, Indexed>::apply(tup, f))...
    >(MapTupleElem<N, Indexed>::apply(tup, f)...);
}

template<unsigned N, class F, class... As> 
auto mapTupleElem(std::tuple<As...> tup, F f) -> decltype(auto) 
{ return __mapTupleElem<N>(std::move(tup), f, TL::WithIndices<TL::List<As...>>{}); }

template<unsigned N, class B, class... As> 
auto replaceTupleElem(std::tuple<As...> tup, B b) -> decltype(auto) 
{ return mapTupleElem<N>(std::move(tup), [&](auto) -> B { return move_if_value<B>(b); }); }

template<class Type, Type value>
struct ReturnConst {
  template<class... As>
  constexpr Type operator()(As...) const { return value; }
};
template<class Type>
struct ReturnNone {
  template<class... As>
  constexpr Option<Type> operator()(As...) const { return Option<Type>(); }
};

// static constexpr auto returnNoneClosure = [](auto){
//   return None<();
// };

template<class Arg, class Result>
using NoMemo = Memo::None<Arg, Result>;

#define FOR_FIELD(MACRO)                                                                  \
  MACRO(0, Function , function      , (std::tuple<>())           )                        \
  MACRO(1, EvNonRec , evNonRec      , (ReturnNone<Result>{})     )                        \
  MACRO(2, Memo     , memo          , (NoMemo<Arg, Result>()))                            \
  MACRO(3, Context  , context       , (std::tuple<>())           )                        \
  MACRO(4, Int      , memoThreshold , (unsigned(0))              )                        \
  /*    ^  ^^^^^^^^^  ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^----> default value */
  /*    |      |         +--------------------------------------------> field name */
  /*    |      +------------------------------------------------------> type param name */
  /*    +-------------------------------------------------------------> index */

template< class Arg
        , class Result  
#                                       define foreach(idx, Type, name, defaultVal)       \
        , class Type = decltype(defaultVal) 
                                        FOR_FIELD(foreach)
#                                       undef foreach
  >
class BottomUpEvaluation {

  template<class A, class R
#                                       define foreach(idx, Type, name, defaultVal)       \
  , class Type ## _
                                        FOR_FIELD(foreach)
#                                       undef foreach
                                        > 
 friend class BottomUpEvaluation;

  std::tuple<> _dummy;
#                                       define foreach(idx, Type, name, defaultVal)       \
  Type _ ## name;
                                        FOR_FIELD(foreach)
#                                       undef foreach

  BottomUpEvaluation(
      std::tuple<
#                                       define foreach(idx, Type, name, defaultVal)       \
                Type,
                                        FOR_FIELD(foreach)
#                                       undef foreach
                std::tuple<>> elems)
    : _dummy()

#                                       define foreach(idx, Type, name, defaultVal)       \
    , _ ## name(std::get<idx>(elems))
                                        FOR_FIELD(foreach)
#                                       undef foreach
  { }

  template<
#                                       define foreach(idx, Type, name, defaultVal)       \
                class Type ## _,
                                        FOR_FIELD(foreach)
#                                       undef foreach
                class... Dummies
    >

  static auto fromTuple(
      std::tuple<
#                                       define foreach(idx, Type, name, defaultVal)       \
        Type ## _,
                                        FOR_FIELD(foreach)
#                                       undef foreach
        std::tuple<>
      > tup)
  { return BottomUpEvaluation< Arg
                             , Result
#                                       define foreach(idx, Type, name, defaultVal)       \
                             , Type ## _
                                        FOR_FIELD(foreach)
#                                       undef foreach
                             >(std::move(tup)); }

  auto intoTuple() &&
  { return std::tuple<
#                                       define foreach(idx, Type, name, defaultVal)       \
             Type,
                                        FOR_FIELD(foreach)
#                                       undef foreach
             std::tuple<>
    >(
#                                       define foreach(idx, Type, name, defaultVal)       \
             move_if_value<Type>(_ ## name),
                                        FOR_FIELD(foreach)
#                                       undef foreach
             make_tuple()); }

public:
  BottomUpEvaluation() 
    : BottomUpEvaluation(
        std::tuple<
#                                       define foreach(idx, Type, name, defaultVal)       \
             Type,
                                        FOR_FIELD(foreach)
#                                       undef foreach
             std::tuple<>
    >(
#                                       define foreach(idx, Type, name, defaultVal)       \
             defaultVal,
                                        FOR_FIELD(foreach)
#                                       undef foreach
             make_tuple()))
  {}


#                                       define foreach(idx, Type, name, defaultVal)       \
  template<class New>                                                                     \
  auto name(New val) &&                                                                   \
  { return fromTuple(replaceTupleElem<idx, New>(std::move(*this).intoTuple(), move_if_value<New>(val))); }\
                                                                                          \
  Type& name() { return _ ## name; }                                                      \


                                        FOR_FIELD(foreach)
#                                       undef foreach

  Result apply(Arg const& toEval) 
  {
    CALL("evaluateBottomUp(...)")

    // static_assert(std::is_same<ResultOf<Function, Arg, Result*>, Result>::value, "evaluation function must have signature `Result eval(Arg term, Result* evaluatedArgs)`");


    /* recursion state. Contains a stack of items that are being recursed on. */
    Recycled<Stack<BottomUpChildIter<Arg>>> recState;
    Recycled<Stack<Result>> recResults;

    recState->push(BottomUpChildIter<Arg>(toEval, _context));

    while (!recState->isEmpty()) {
      if (recState->top().hasNext(_context)) {
        Arg t = recState->top().next(_context);

        Option<Result> nonRec = _evNonRec(t);
        if (nonRec) {
          recResults->push(move_if_value<Result>(*nonRec));

        } else {
          Option<Result> cached = _memo.get(t);
          if (cached.isSome()) {
            recResults->push(std::move(cached).unwrap());
          } else {
            recState->push(BottomUpChildIter<Arg>(t, _context));
          }

        }

      } else {

        BottomUpChildIter<Arg> orig = recState->pop();

        Result* argLst = orig.nChildren(_context) == 0 
          ? nullptr 
          : static_cast<Result*>(&((*recResults)[recResults->size() - orig.nChildren(_context)]));

        Result eval = _memo.getOrInit(orig.self(), 
                        [&](){ return _function(orig.self(), argLst); });

        DEBUG("evaluated: ", orig.self(), " -> ", eval);
        recResults->pop(orig.nChildren(_context));
        recResults->push(std::move(eval));
      }
    }
    ASS(recState->isEmpty())


    ASS(recResults->size() == 1);
    auto result = recResults->pop();
    DEBUG("eval result: ", toEval, " -> ", result);
    return result;
  }
};

template<class R, class A, class F, class Memo, class... Context>
R evalBottomUpWithMemo(A const& term, F fun, Memo& memo, Context ... cs)
{
  struct Eval {
    using Result = R;
    using Arg = A;
    F& fun;
    Result operator()(Arg const& a, Result* rs)
    { return fun(a,rs); }
  };
  return evaluateBottomUpWithMemo(term, Eval{fun}, memo, std::move(cs)...);
}

template<class R, class A, class F, class... Context>
R evalBottomUp(A const& term, F fun, Context ... cs)
{
  Memo::None<A, R> memo;
  return evalBottomUpWithMemo<R>(term, std::move(fun), memo, std::move(cs)...);
}


} // namespace Lib

#undef DEBUG
#endif // __LIB__BOTTOM_UP_EVALUATION_HPP__

