/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#include "Inverters.hpp"

namespace Kernel {
namespace Rebalancing {
namespace Inverters {
#define DEBUG(...) //DBG(__VA_ARGS__)

template<class A> void __ignoreWarnUnusedLocalTypedefHack() {}

#define CASE_INVERT(sort, fun, expr)                                           \
  case NumTraits<sort>::fun##I: {                                              \
    using Number = NumTraits<sort>;                                            \
    __ignoreWarnUnusedLocalTypedefHack<Number>();                              \
    return expr;                                                               \
  }

#define CASE_INVERT_INT(fun, expr) CASE_INVERT(IntegerConstantType, fun, expr)

#define CASE_INVERT_FRAC(fun, expr)                                            \
  CASE_INVERT(RealConstantType, fun, expr) \
  CASE_INVERT(RationalConstantType, fun, expr)

template<class NumTraits> bool canInvertDiv(const InversionContext &ctxt);
template<class NumTraits> TermList doInvertDiv(const InversionContext &ctxt);
bool canInvertMulInt(const InversionContext &ctxt);
TermList doInvertMulInt(const InversionContext &ctxt);

template <class Number> bool nonZero(const TermList &t);

bool dtorIsPredicate(Signature::Symbol const& ctor, unsigned index) 
{ return ctor.fnType()->arg(index) == AtomicSort::boolSort(); }

bool NumberTheoryInverter::canInvertTop(const InversionContext &ctxt) {
  auto &t = ctxt.topTerm();
  auto fun = t.functor();

  DEBUG("canInvert ", ctxt.topTerm().toString(), "@", ctxt.topIdx())
  if (theory->isInterpretedFunction(fun)) {
    auto inter = theory->interpretFunction(fun);
    switch (inter) {
      CASE_INVERT_FRAC(add, true)
      CASE_INVERT_FRAC(minus, true)
      CASE_INVERT_FRAC(mul, nonZero<Number>(t[1 - ctxt.topIdx()]))
      CASE_INVERT_FRAC(div, canInvertDiv<Number>(ctxt))
      CASE_INVERT_INT(mul, canInvertMulInt(ctxt))
      CASE_INVERT_INT(add, true)
      CASE_INVERT_INT(minus, true)
    // case Theory::Interpretation::ARRAY_STORE:
    //   /* store(t, i, x) = s ==> x = select(s, i) */
    //   return ctxt.topIdx() == 2;
      default:;
    }
    // DBG("WARNING: unknown interpreted function: ", t.toString())
    return false;
  } else { 
    /* cannot invert uninterpreted functions */
    return false;
  }
}

#define CASE_DO_INVERT(sort, fun, expr)                                        \
  case NumTraits<sort>::fun##I: {                                              \
    using Number = NumTraits<sort>;                                            \
    __ignoreWarnUnusedLocalTypedefHack<Number>();                              \
    return expr;                                                               \
  }

#define CASE_DO_INVERT_FRAC(fun, expr)                                         \
  CASE_DO_INVERT(RealConstantType, fun, expr)                                  \
  CASE_DO_INVERT(RationalConstantType, fun, expr)

#define CASE_DO_INVERT_ALL(fun, expr)                                          \
  CASE_DO_INVERT_INT(fun, expr)                                                \
  CASE_DO_INVERT_FRAC(fun, expr)

#define CASE_DO_INVERT_INT(fun, expr)                                          \
  CASE_DO_INVERT(IntegerConstantType, fun, expr)

TermList NumberTheoryInverter::invertTop(const InversionContext &ctxt) {
  ASS(canInvertTop(ctxt))
  // DBG("inverting: ", ctxt.topTerm().toString())
  auto &t = ctxt.topTerm();
  auto index = ctxt.topIdx();
  auto toWrap = ctxt.toWrap();
  auto fun = t.functor();
  DEBUG("inverting ", ctxt.topTerm().toString())
  if(theory->isInterpretedFunction(fun)) {
    switch (theory->interpretFunction(fun)) {

      CASE_DO_INVERT_ALL(add, Number::add(toWrap, Number::minus(t[1 - index])))
      CASE_DO_INVERT_ALL(minus, Number::minus(toWrap))

      CASE_DO_INVERT_FRAC( mul, Number::mul(toWrap, Number::div(Number::one(), t[1 - index])))
      CASE_DO_INVERT_FRAC(div, doInvertDiv<Number>(ctxt))

      CASE_DO_INVERT_INT(mul, doInvertMulInt(ctxt))

    default:
      ASSERTION_VIOLATION;
    }
  } else {
    ASSERTION_VIOLATION
  }
};

template<class NumTraits>
bool tryInvertDiv(const InversionContext &ctxt, TermList &out) {
  /* auto s = ctxt.topTerm()[0]; */
  auto t = ctxt.topTerm()[1];
  auto u = ctxt.toWrap();

  if (ctxt.topIdx() == 0) {
    // s / t = u ===> s = u * t (if nonZero(t) )
    if (nonZero<NumTraits>(t)) {
      out = NumTraits::mul(u, t);
      return true;
    } else {
      return false;
    }
  } else {
    ASS_EQ(ctxt.topIdx(), 1)
    return false;
  }
}

template<class NumTraits>
TermList doInvertDiv(const InversionContext &ctxt) {
  TermList out;
  ALWAYS(tryInvertDiv<NumTraits>(ctxt, out)) 
  return out;
}

template<class NumTraits>
bool canInvertDiv(const InversionContext &ctxt) {
  TermList _inv;
  return tryInvertDiv<NumTraits>(ctxt, _inv);
}

bool tryInvertMulInt(const InversionContext &ctxt, TermList &out) {
  using Number = NumTraits<IntegerConstantType>;

  auto a_ = ctxt.topTerm()[1 - ctxt.topIdx()];
  IntegerConstantType a;
  if ( theory->tryInterpretConstant(a_, a)) {
    if (a == IntegerConstantType(1)) {
      out = ctxt.toWrap();
      return true;

    } else if (a == IntegerConstantType(-1)) {
      out = Number::mul(a_, ctxt.toWrap());
      return true;

    } else {
      return false;
    }
  } else {
    return false;
  }
}

TermList doInvertMulInt(const InversionContext &ctxt) {
  TermList out;
  ALWAYS(tryInvertMulInt(ctxt, out)) 
  return out;
}

bool canInvertMulInt(const InversionContext &ctxt) {
  TermList _inv;
  return tryInvertMulInt(ctxt, _inv);
}

template <class Number> bool nonZero(const TermList &t) {
  typename Number::ConstantType c;
  return theory->tryInterpretConstant(t, c) && Number::constant(0) != c;
}


} // namespace Inverters
} // namespace Rebalancing
} // namespace Kernel
