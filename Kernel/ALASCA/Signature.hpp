/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __ALASCA_Signature__
#define __ALASCA_Signature__

#include "Kernel/NumTraits.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/Signature.hpp"

namespace Kernel {

enum class AlascaPredicate {
  EQ,
  NEQ,
  GREATER,
  GREATER_EQ,
};

/** returns true iff the predicate is > or >= */
inline bool isInequality(AlascaPredicate const& self) 
{
  switch(self) {
    case AlascaPredicate::EQ: 
    case AlascaPredicate::NEQ: return false;
    case AlascaPredicate::GREATER: 
    case AlascaPredicate::GREATER_EQ: return true;
  }
  ASSERTION_VIOLATION
}

template<class NumTraits>
Literal* createLiteral(AlascaPredicate self, TermList t)
{
  auto zero = NumTraits::linMul(NumTraits::constant(0), NumTraits::one());
  switch(self) {
    case AlascaPredicate::EQ: return NumTraits::eq(true, t, zero);
    case AlascaPredicate::NEQ: return NumTraits::eq(false, t, zero);
    case AlascaPredicate::GREATER: return NumTraits::greater(true, t, zero);
    case AlascaPredicate::GREATER_EQ: return NumTraits::geq(true, t, zero);
  }
  ASSERTION_VIOLATION
}
bool isIsInt(AlascaPredicate const& self);

inline std::ostream& operator<<(std::ostream& out, AlascaPredicate const& self)
{ 
  switch(self) {
    case AlascaPredicate::EQ: return out << "==";
    case AlascaPredicate::NEQ: return out << "!=";
    case AlascaPredicate::GREATER: return out << ">";
    case AlascaPredicate::GREATER_EQ: return out << ">=";
  } 
  ASSERTION_VIOLATION
}


template<class NumTraits> 
struct AlascaSignature : public NumTraits {
  using Numeral = typename NumTraits::ConstantType;
  static Option<Numeral> oneN;
  static Option<TermList> oneT;
  static Option<TermList> sortT;
  static Option<TermList> zeroT;

  template<class T>
  static Option<Numeral const&> tryNumeral(T t) {
    if (t == AlascaSignature::one()) {
      return Option<Numeral const&>(oneN.unwrapOrInit([&]() { return NumTraits::constant(1); }));
    } else {
      return NumTraits::ifLinMul(t, [](auto& c, auto t) {
          return someIf(t == NumTraits::one(), [&]() -> auto& { return c; });
      }).flatten();
    }
  }

  template<class T> static void numeralF(T) = delete;
  template<class T> static void constantTl(T) = delete;

  static auto addF() { return NumTraits::addF(); }
  static auto minusF() { return NumTraits::linMulF(Numeral(-1)); }

  template<class T>
  static auto isUninterpreted(T t) 
  { return !NumTraits::isFloor(t) && !NumTraits::isAdd(t) && !NumTraits::isLinMul(t) && !AlascaSignature::isOne(t); }



  static auto isVar(TermList t) { return t.isVar(); }
  static auto isVar(Term* t) { return false; }

  // TODO check if any of the uses should include isOne as atomic
  template<class T>
  static auto isAtomic(T t) 
  { return !isVar(t) && AlascaSignature::isUninterpreted(t); }

  template<class T>
  static TermList linMul(Numeral const& c, T t) 
  { return c == 1 ? TermList(t) : NumTraits::linMul(c, t); }

  // TODO faster
  static TermList const& zero() { return zeroT.unwrapOrInit([&]() { return NumTraits::zero(); }); }
  static TermList const& one() { return oneT.unwrapOrInit([&]() { return NumTraits::one(); }); }
  static TermList const& sort() { return sortT.unwrapOrInit([&]() { return NumTraits::sort(); }); }

  static bool isOne(unsigned f) { return AlascaSignature::one().term()->functor() == f; }
  static bool isOne(TermList t) { return AlascaSignature::one() == t; }
  static bool isOne(Term* t) { return AlascaSignature::one() == TermList(t); }

  static Kernel::TermList numeralTl(int c) 
  { return numeralTl(NumTraits::constant(c)); }

  static Kernel::TermList numeralTl(typename NumTraits::ConstantType const& c) 
  { return TermList(NumTraits::linMul(c, NumTraits::one())); }

};

template<class NumTraits>
AlascaSignature<NumTraits> asig(NumTraits n) 
{ return AlascaSignature<NumTraits> {}; }

template<typename NumTraits> Option<typename AlascaSignature<NumTraits>::Numeral> AlascaSignature<NumTraits>::oneN;
template<typename NumTraits> Option<TermList> AlascaSignature<NumTraits>::oneT;
template<typename NumTraits> Option<TermList> AlascaSignature<NumTraits>::sortT;
template<typename NumTraits> Option<TermList> AlascaSignature<NumTraits>::zeroT;

// TODO rename
template<class NumTraits, class F>
Option<std::invoke_result_t<F, AlascaPredicate, TermList, unsigned>> ifAlascaLiteral(Literal* lit, F f) {
  // TODO assert normalized
  if (NumTraits::isGreater(lit->functor())) {
    ASS(lit->termArg(1) == NumTraits::constantTl(0))
    return some(f(AlascaPredicate::GREATER   , lit->termArg(0), 0));
  }
  if (NumTraits::isGeq(lit->functor())    ) {
    ASS(lit->termArg(1) == NumTraits::constantTl(0))
    return some(f(AlascaPredicate::GREATER_EQ, lit->termArg(0), 0));
  }
  if (lit->isEquality() && SortHelper::getEqualityArgumentSort(lit) == NumTraits::sort()) {
    auto i = 0;
    if (auto n = NumTraits::tryNumeral(lit->termArg(0))) {
      if (*n == 0) {
        i++;
      }
    }
    return some(f(lit->isPositive() ? AlascaPredicate::EQ : AlascaPredicate::NEQ, lit->termArg(i), i));
  }
  return {};
}

// TODO rename
template<class NumTraits, class T, class F>
[[deprecated("use AlascaSignature instead")]]
auto ifNumMul(T term, F f) {
  return NumTraits::ifMul(term, [&](auto l, auto r) {
      ASS(!NumTraits::isNumeral(r))
      return NumTraits::ifNumeral(l, [&](auto l) { return f(l, r, 1); });
  }).flatten()
  .orElse([&](){
      return NumTraits::ifMinus(term, [&](auto t) { return  f(NumTraits::constant(-1), t, 0); });
      });
}

// TODO rename
template<class NumTraits, class T>
auto isNumMul(T term) 
{ return ifNumMul<NumTraits>(term, [](auto...) { return 0; }).isSome(); }

// TODO rename
template<class NumTraits, class T>
auto isAlascaLiteral(T term) 
{ return ifAlascaLiteral<NumTraits>(term, [](auto...) { return 0; }).isSome(); }

// TODO rename
template<class NumTraits, class T>
auto isUninterpreted(T t) 
{ return !NumTraits::isFloor(t) && !NumTraits::isAdd(t) && !isNumMul<NumTraits>(t); }


} // namespace Kernel


#endif // __ALASCA_Signature__
