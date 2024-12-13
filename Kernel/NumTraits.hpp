/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#ifndef __NUM_TRAITS_H__
#define __NUM_TRAITS_H__

#include "Forwards.hpp"
#include "Term.hpp"
#include "Theory.hpp"
#include "Signature.hpp"

namespace Kernel {

/** This struct provides a unified interface to the "number" theories. i.e. these 
 * are the theories of integers, rationals, and reals.
 *
 * For each ConstantType in {IntegerConstantType, RationalConstantType, RealConstantType},
 * there is a specialisation NumTraits<ConstantType> which provides functions for building 
 * terms and literals, accessing the functor/interpretation of some predicate or function. 
 *
 * There are associated constants for the related Sorts::DefaultSorts values, constexpr s 
 * to pattern match on interpretations and special constants (like zero and one), etc.
 *
 * =====
 *
 * To access the functor of some symbol "sym", there is a function 
 * static unsigned symF();
 *
 * e.g.: NumTraits<IntegerConstantType>::lessF()
 *       NumTraits<IntegerConstantType>::addF()
 *
 * =====
 *
 * To access the interpretation of some symbol "sym", there is a constant 
 * static Theory::Interpretation symI;
 *
 * e.g.: NumTraits<IntegerConstantType>::lessI // == Theory::Interpretation INT_LESS;
 *       NumTraits<IntegerConstantType>::addI  // == Theory::Interpretation INT_PLUS;
 *
 * =====
 *
 * To build a term from some function symbol "sym", there is a function
 * static TermList static sym(TermList...);
 *
 * e.g.: NumTraits<IntegerConstantType>::add(lhs, rhs) 
 *
 * =====
 *
 * To build a literal from some predicate symbol "sym", there is a function
 * static Literal* static sym(bool polarity, TermList...);
 *
 * e.g.: NumTraits<IntegerConstantType>::less(true, lhs, rhs) 
 *
 * =====
 *
 * For a complete picture build the doxygen documentaion.
 *
 */
template<class ConstantType>
struct NumTraits;

#define FOR_ARITY_RANGE_0(f)
#define FOR_ARITY_RANGE_1(f) f(0)
#define FOR_ARITY_RANGE_2(f) FOR_ARITY_RANGE_1(f), f(1)
#define FOR_ARITY_RANGE_3(f) FOR_ARITY_RANGE_2(f), f(2)
#define FOR_ARITY_RANGE_4(f) FOR_ARITY_RANGE_3(f), f(3)
#define FOR_ARITY_RANGE_5(f) FOR_ARITY_RANGE_4(f), f(4)
#define FOR_ARITY_RANGE_6(f) FOR_ARITY_RANGE_5(f), f(5)
#define FOR_ARITY_RANGE_7(f) FOR_ARITY_RANGE_6(f), f(6)
#define FOR_ARITY_RANGE_8(f) FOR_ARITY_RANGE_7(f), f(7)
#define FOR_ARITY_RANGE_9(f) FOR_ARITY_RANGE_8(f), f(8)
#define FOR_ARITY_RANGE_10(f) FOR_ARITY_RANGE_9(f), f(9)
#define FOR_ARITY_RANGE_11(f) FOR_ARITY_RANGE_10(f), f(10)

#define FOR_ARITY_RANGE(arity) FOR_ARITY_RANGE_ ## arity

#define IMPL_NUM_TRAITS__ARG_DECL_1(Type) Type a1
#define IMPL_NUM_TRAITS__ARG_DECL_2(Type) IMPL_NUM_TRAITS__ARG_DECL_1(Type), Type a2
#define IMPL_NUM_TRAITS__ARG_DECL_3(Type) IMPL_NUM_TRAITS__ARG_DECL_2(Type), Type a3
#define IMPL_NUM_TRAITS__ARG_DECL_4(Type) IMPL_NUM_TRAITS__ARG_DECL_3(Type), Type a4
#define IMPL_NUM_TRAITS__ARG_DECL_5(Type) IMPL_NUM_TRAITS__ARG_DECL_4(Type), Type a5
#define IMPL_NUM_TRAITS__ARG_DECL_6(Type) IMPL_NUM_TRAITS__ARG_DECL_5(Type), Type a6
#define IMPL_NUM_TRAITS__ARG_DECL_7(Type) IMPL_NUM_TRAITS__ARG_DECL_6(Type), Type a7
#define IMPL_NUM_TRAITS__ARG_DECL_8(Type) IMPL_NUM_TRAITS__ARG_DECL_7(Type), Type a8
#define IMPL_NUM_TRAITS__ARG_DECL_9(Type) IMPL_NUM_TRAITS__ARG_DECL_8(Type), Type a9
#define IMPL_NUM_TRAITS__ARG_DECL_10(Type) IMPL_NUM_TRAITS__ARG_DECL_9(Type), Type a10

#define IMPL_NUM_TRAITS__ARG_DECL(Type, arity) IMPL_NUM_TRAITS__ARG_DECL_ ## arity (Type)
#define IMPL_NUM_TRAITS__ARG_EXPR(arity) IMPL_NUM_TRAITS__ARG_DECL_ ## arity ()

#define IMPL_NUM_TRAITS__INTERPRETED_SYMBOL(name, Name, SORT_SHORT, _INTERPRETATION)      \
    static constexpr Theory::Interpretation name ## I = Theory::SORT_SHORT ## _INTERPRETATION;      \
                                                                                          \
    static unsigned name ## F() {                                                         \
      static const unsigned functor = env.signature->getInterpretingSymbol(name ## I);    \
      return functor;                                                                     \
    }                                                                                     \


#define IMPL_NUM_TRAITS__INTERPRETED_PRED(name, Name, SORT_SHORT, _INTERPRETATION, arity) \
    IMPL_NUM_TRAITS__INTERPRETED_SYMBOL(name, Name, SORT_SHORT, _INTERPRETATION)          \
                                                                                          \
    static bool is ## Name(unsigned f)                                                    \
    { return theory->isInterpretedPredicate(f, name ## I); }                              \
                                                                                          \
    template<class F>                                                                     \
    static auto if ## Name(Literal* t, F fun) {                                           \
      return someIf(is ## Name(t->functor()),                                             \
          [&]() { return fun(t->isPositive(), FOR_ARITY_RANGE(arity)(t->termArg)); });    \
    }                                                                                     \
                                                                                          \
    static Literal* name(bool polarity, IMPL_NUM_TRAITS__ARG_DECL(TermList, arity)) {     \
      return Literal::create(                                                             \
                  name##F(),                                                              \
                  polarity,                                                               \
                  { IMPL_NUM_TRAITS__ARG_EXPR( arity ) });                                \
    }                                                                                     \



#define IMPL_NUM_TRAITS__INTERPRETED_FUN(name, Name, SORT_SHORT, _INTERPRETATION, arity)  \
    IMPL_NUM_TRAITS__INTERPRETED_SYMBOL(name, Name, SORT_SHORT, _INTERPRETATION)          \
                                                                                          \
    static bool is ## Name(unsigned f)                                                    \
    { return theory->isInterpretedFunction(f, name ## I); }                               \
                                                                                          \
    static bool is ## Name(Term* t)                                                       \
    { return is ## Name(t->functor()); }                                                  \
                                                                                          \
    static bool is ## Name(TermList t)                                                    \
    { return t.isTerm() && is ## Name(t.term()); }                                        \
                                                                                          \
    template<class F>                                                                     \
    static auto if ## Name(Term* t, F fun) {                                              \
      return someIf(is ## Name(t->functor()),                                             \
          [&]() { return fun(FOR_ARITY_RANGE(arity)(t->termArg)); });                     \
    }                                                                                     \
    template<class F>                                                                     \
    static auto if ## Name(TermList t, F fun) {                                           \
      return someIf(t.isTerm(), [&]() { return if ## Name(t.term(), fun); }).flatten();   \
    }                                                                                     \
                                                                                          \
    static TermList name(IMPL_NUM_TRAITS__ARG_DECL(TermList, arity)) {                    \
      return TermList(                                                                    \
          Term::create(                                                                   \
            name##F(),                                                                    \
            { IMPL_NUM_TRAITS__ARG_EXPR(arity) }));                                       \
    }                                                                                     \

#define IMPL_NUM_TRAITS__SPECIAL_CONSTANT(name, Name, value)                              \
    static ConstantType name ## C() {                                                     \
      return ConstantType(value);                                                         \
    }                                                                                     \
    static Term* name ## T() {                                                            \
      static Term* trm = theory->representConstant(name ## C());                          \
      return trm;                                                                         \
    }                                                                                     \
    static unsigned name ## F() {                                                         \
      static unsigned f = name ## T()->functor();                                         \
      return f;                                                                           \
    }                                                                                     \
    static TermList name()                                                                \
    { return TermList(name ## T()); }                                                     \
                                                                                          \
    static bool is ## Name(const TermList& l)                                             \
    { return l == name(); }                                                               \

#define IMPL_NUM_TRAITS__QUOTIENT_REMAINDER(SHORT, X)                                     \
    IMPL_NUM_TRAITS__INTERPRETED_FUN( quotient ## X,  Quotient ## X, SHORT,  _QUOTIENT_ ## X, 2)    \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(remainder ## X, Remainder ## X, SHORT, _REMAINDER_ ## X, 2)    \
    

#define IMPL_NUM_TRAITS(CamelCase, lowerCase, LONG, SHORT)                                \
  template<> struct NumTraits<CamelCase ## ConstantType> {                                \
    /* dummy operator== to be able to compare Coproduct<IntTraits, RatTraits, ...> */     \
    friend bool operator==(NumTraits const& l, NumTraits const& r)                        \
    { return true; }                                                                      \
                                                                                          \
    friend bool operator!=(NumTraits const& l, NumTraits const& r)                        \
    { return !(l == r); }                                                                 \
                                                                                          \
    using ConstantType = CamelCase ## ConstantType;                                       \
    static TermList sort() { return AtomicSort::lowerCase ## Sort(); };                   \
                                                                                          \
    template<class I1, class I2, class... Is>                                             \
    static TermList sum(I1 i1, I2 i2, Is... is)                                           \
    { return sum(concatIters(i1, i2, is...)); };                                          \
                                                                                          \
    static TermList mulSimpl(ConstantType c, TermList t)                                  \
    { return c == ConstantType(1) ? t                                                     \
           : c == ConstantType(-1) ? (t == zero() ? t : minus(t))                         \
           : t == zero() ? zero()                                                         \
           : t == one() ? constantTl(c)                                                   \
           : NumTraits::mul(constantTl(c), t); }                                          \
                                                                                          \
    template<class Iter>                                                                  \
    static TermList sum(Iter iter) {                                                      \
      if (iter.hasNext()) {                                                               \
        auto out = iter.next();                                                           \
        while (iter.hasNext()) {                                                          \
          out = NumTraits::add(iter.next(), out);                                         \
        }                                                                                 \
        return out;                                                                       \
      } else {                                                                            \
        return NumTraits::zero();                                                         \
      }                                                                                   \
    };                                                                                    \
                                                                                          \
    template<class Iter>                                                                  \
    static TermList product(Iter iter) {                                                  \
      if (iter.hasNext()) {                                                               \
        auto out = iter.next();                                                           \
        while (iter.hasNext()) {                                                          \
          out = NumTraits::mul(iter.next(), out);                                         \
        }                                                                                 \
        return out;                                                                       \
      } else {                                                                            \
        return NumTraits::one();                                                         \
      }                                                                                   \
    };                                                                                    \
                                                                                          \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(less,    Less,    SHORT, _LESS,          2)         \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(leq,     Leq,     SHORT, _LESS_EQUAL,    2)         \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(greater, Greater, SHORT, _GREATER,       2)         \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(geq,     Geq,     SHORT, _GREATER_EQUAL, 2)         \
                                                                                          \
    static Literal* eq(bool polarity, TermList lhs, TermList rhs)                         \
    { return Literal::createEquality(polarity, lhs, rhs, sort()); }                       \
                                                                                          \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(toInt,   ToInt,   SHORT, _TO_INT       , 1)          \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(toRat,   ToRat,   SHORT, _TO_RAT       , 1)          \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(toReal,  ToReal,  SHORT, _TO_REAL      , 1)          \
                                                                                          \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(isInt,   IsInt,   SHORT, _IS_INT       , 1)         \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(isRat,   IsRat,   SHORT, _IS_RAT       , 1)         \
    IMPL_NUM_TRAITS__INTERPRETED_PRED(isReal,  IsReal,  SHORT, _IS_REAL      , 1)         \
                                                                                          \
    IMPL_NUM_TRAITS__QUOTIENT_REMAINDER(SHORT, E)                                         \
    IMPL_NUM_TRAITS__QUOTIENT_REMAINDER(SHORT, T)                                         \
    IMPL_NUM_TRAITS__QUOTIENT_REMAINDER(SHORT, F)                                         \
                                                                                          \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(minus, Minus, SHORT, _UNARY_MINUS, 1)                \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(binMinus, BinMinus, SHORT, _MINUS, 2)                \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(add  , Add  , SHORT, _PLUS       , 2)                \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(mul  , Mul  , SHORT, _MULTIPLY   , 2)                \
    IMPL_NUM_TRAITS__INTERPRETED_FUN(floor, Floor, SHORT, _FLOOR, 1)                      \
    __NUM_TRAITS_IF_FRAC(SHORT,                                                           \
        IMPL_NUM_TRAITS__INTERPRETED_FUN(div, Div, SHORT, _QUOTIENT, 2)                   \
        static ConstantType constant(int num, int den) { return ConstantType(num, den); } \
        static Term* constantT(int num, int den) { return theory->representConstant(constant(num, den)); }    \
        static TermList constantTl(int num, int den) { return TermList(constantT(num, den)); }      \
        static bool isFractional() { return true; }                                       \
    )                                                                                     \
                                                                                          \
    __NUM_TRAITS_IF_NOT_FRAC(SHORT,                                                       \
        static bool isFractional() { return false; }                                      \
        IMPL_NUM_TRAITS__INTERPRETED_PRED(divides, Divides, SHORT, _DIVIDES, 2)           \
    )                                                                                     \
                                                                                          \
    IMPL_NUM_TRAITS__SPECIAL_CONSTANT(one , One , 1)                                      \
    IMPL_NUM_TRAITS__SPECIAL_CONSTANT(zero, Zero, 0)                                      \
                                                                                          \
                                                                                          \
    static ConstantType constant(int i) { return ConstantType(i); }                       \
    static Term* constantT(int i) { return constantT(constant(i)); }                      \
    static Term* constantT(ConstantType i) { return theory->representConstant(i); }       \
    static TermList constantTl(int i) { return TermList(constantT(i)); }                  \
    static TermList constantTl(ConstantType i) { return TermList(constantT(i)); }         \
    template<class TermOrFunctor>                                                         \
    static Option<ConstantType> tryNumeral(TermOrFunctor t) {                             \
      ConstantType out;                                                                   \
      if (theory->tryInterpretConstant(t,out)) {                                          \
        return Option<ConstantType>(out);                                                 \
      } else {                                                                            \
        return Option<ConstantType>();                                                    \
      }                                                                                   \
    }                                                                                     \
    template<class TermOrFunctor>                                                         \
    static bool isNumeral(TermOrFunctor t) { return tryNumeral(t).isSome(); }             \
    template<class TermOrFunctor>                                                         \
    static bool isNumeral(TermOrFunctor t, ConstantType n) { return tryNumeral(t) == some(n); }     \
    template<class Term, class F>                                                         \
    static auto ifNumeral(Term t, F fun) -> Option<std::invoke_result_t<F, ConstantType&&>> \
    { return tryNumeral(t).map([&](ConstantType n) { return fun(std::move(n)); }); }      \
    static unsigned numeralF(ConstantType c) { return constantT(c)->functor(); }          \
                                                                                          \
    static const char* name() {return #CamelCase;}                                        \
  };                                                                                      \

#define __NUM_TRAITS_IF_FRAC(sort, ...) __NUM_TRAITS_IF_FRAC_ ## sort (__VA_ARGS__)
#define __NUM_TRAITS_IF_FRAC_INT(...) 
#define __NUM_TRAITS_IF_FRAC_REAL(...) __VA_ARGS__
#define __NUM_TRAITS_IF_FRAC_RAT(...) __VA_ARGS__

#define __NUM_TRAITS_IF_NOT_FRAC(sort, ...) __NUM_TRAITS_IF_NOT_FRAC_ ## sort (__VA_ARGS__)
#define __NUM_TRAITS_IF_NOT_FRAC_INT(...)  __VA_ARGS__
#define __NUM_TRAITS_IF_NOT_FRAC_REAL(...)
#define __NUM_TRAITS_IF_NOT_FRAC_RAT(...)

IMPL_NUM_TRAITS(Rational, rational, RATIONAL, RAT )
IMPL_NUM_TRAITS(Real    , real    , REAL    , REAL)
IMPL_NUM_TRAITS(Integer , int     , INTEGER , INT )

#define FOR_NUM_TRAITS(macro)                                                             \
  macro(Kernel::NumTraits<Kernel:: IntegerConstantType>)                                  \
  macro(Kernel::NumTraits<Kernel::    RealConstantType>)                                  \
  macro(Kernel::NumTraits<Kernel::RationalConstantType>)                                  \

using IntTraits  = NumTraits< IntegerConstantType>;
using RatTraits  = NumTraits<RationalConstantType>;
using RealTraits = NumTraits<    RealConstantType>;


template<class Clsr>
auto forAnyNumTraits(Clsr clsr) {
  return clsr( IntTraits{}) 
      || clsr( RatTraits{})
      || clsr(RealTraits{});
}

template<class Clsr>
auto forAllNumTraits(Clsr clsr) {
  return clsr( IntTraits{}) 
      && clsr( RatTraits{})
      && clsr(RealTraits{});
}

template<class Clsr>
auto numTraitsIter(Clsr clsr) {
  return getConcatenatedIterator( clsr( IntTraits{}),
         getConcatenatedIterator( clsr( RatTraits{}),
                                  clsr(RealTraits{})));
}

template<class Clsr>
auto forEachNumTraits(Clsr clsr) {
  clsr( IntTraits{});
  clsr( RatTraits{});
  clsr(RealTraits{});
}

template<class Clsr>
auto tryNumTraits(Clsr clsr) {
               return clsr( IntTraits{}) 
      || [&] { return clsr( RatTraits{}); }
      || [&] { return clsr(RealTraits{}); };
}


template<class Clsr>
auto tryFracNumTraits(Clsr clsr) {
                      clsr( RatTraits{})
      || [&] { return clsr(RealTraits{}); };
}

template<class IfInt, class Else>
auto ifIntTraits(IntTraits n, IfInt ifIntF, Else) 
{ return ifIntF(std::move(n)); }

template<class IfInt, class Else, class NumTraits>
auto ifIntTraits(NumTraits n, IfInt, Else elseF) 
{ return elseF(std::move(n)); }


// template<class T, class IfT, class IfNotT, class Val>
// auto ifOfType(Val val, IfT ifT, IfNotT ifNotT)
// { return ifNotT(std::move(val)); }
//
// template<class T, class IfT, class IfNotT>
// auto ifOfType<T, IfT, IfNotT, T>(T val, IfT ifT, IfNotT ifNotT)
// { return ifT(std::move(val)); }


template<class T, class Val, class IfT, class IfNotT>
struct IfOfType 
{
  auto operator()(Val val, IfT ifT, IfNotT ifNotT)
  { return ifNotT(std::move(val)); }
};

template<class T, class IfT, class IfNotT>
struct IfOfType<T, T, IfT, IfNotT>
{
  auto operator()(T val, IfT ifT, IfNotT ifNotT)
  { return ifT(std::move(val)); }
};

template<class T, class Val, class IfT, class IfNotT>
auto ifOfType(Val val, IfT ifT, IfNotT ifNotT)
{ return IfOfType<T, Val, IfT, IfNotT>{}(std::move(val), ifT, ifNotT); }


}


#endif // __NUM_TRAITS_H__
