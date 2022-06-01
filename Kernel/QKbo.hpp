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
 * @file QKbo.hpp
 *
 * Defines a Knuth-Bendix Ordering for linear arithmetic. This ignores numeral multiplications,
 * (i.e. it considers terms t and nt where n is a numeral as equaivalent). Further it is not sensitive to
 * bracketing or permuting of multiplications and additions.
 */

#ifndef __QKBO__
#define __QKBO__

#include "Forwards.hpp"

#include "Lib/DArray.hpp"
#include "Kernel/IRC.hpp"

#include "Ordering.hpp"
#include "Lib/DArray.hpp"
#include "Kernel/LaLpo.hpp"
#include "Kernel/KBO.hpp"
#include "Kernel/OrderingUtils.hpp"
#include "Lib/Coproduct_.hpp"

namespace Kernel {

using namespace Lib;

// template<class... Is>
// class CoproductIter 
// {
//   Coproduct<Is...> _inner;
// public:
//   DECL_ELEMENT_TYPE(ELEMENT_TYPE(TypeList::Get<0, TypeList::List<Is...>>));
// 
//   template<class I>
//   CoproductIter(I i) : _inner(Coproduct<Is...>(std::move(i))) {}
// 
//   bool hasNext()
//   { 
//     using Co = Coproduct_<unsigned, unsigned>;
//     Co bla_(Co::variant<0>(0));
//     Co& bla = bla_;
//     unsigned lala_ = 0;
//     unsigned& lala = lala_;
//     auto clsr = [](auto & x) {
//       // static_assert(!std::is_const<std::remove_reference_t<decltype(x)>>::value, "fail");
//       static_assert(std::is_same<decltype(x), unsigned&>::value, "fail");
//     };
//     #define BLA(REF) \
//     static_assert(std::is_same< \
//       std::result_of_t<decltype(clsr)(unsigned REF)>, \
//       void \
//     >::value,"la" );
// 
//     BLA(&)
//     bla.apply(clsr);
//   ASSERTION_VIOLATION
//   }
//     // Coproduct<Is...>& inner = _inner;
//     // return inner.apply([](auto& x) { 
//     //   static_assert(!std::is_const<std::remove_reference_t<decltype(x)>>::value);
//     //   return x.hasNext();}); }
// 
//   OWN_ELEMENT_TYPE next()
//   { ASSERTION_VIOLATION }
//   // { return _inner.apply([](auto&& x) { return x.next();}); }
// 
//   bool knowsSize() const 
//   { return _inner.apply([](auto& x) { return x.knowsSize();}); }
// 
//   size_t size() const
//   { return _inner.apply([](auto& x) { return x.size();}); }
// };

template<class IfIter, class ElseIter>
static auto ifElseIter(bool cond, IfIter ifIter, ElseIter elseIter) 
{ return iterTraits(
         cond ? CoproductIter<Lib::ResultOf<IfIter>, Lib::ResultOf<ElseIter>>(ifIter())
              : CoproductIter<Lib::ResultOf<IfIter>, Lib::ResultOf<ElseIter>>(elseIter())); }


// TODO move to right place (IRC.hpp ?)
struct SignedTerm 
{
  Sign sign;
  TermList term;

  static SignedTerm pos(TermList t) 
  { return { .sign = Sign::Pos, .term = t, }; }

  static SignedTerm neg(TermList t) 
  { return { .sign = Sign::Neg, .term = t, }; }

  static SignedTerm zero(TermList t) 
  { return { .sign = Sign::Zero, .term = t, }; }

  friend std::ostream& operator<<(std::ostream& out, SignedTerm const& self)
  { return out << self.sign << self.term; }

  friend bool operator==(SignedTerm const& lhs, SignedTerm const& rhs)
  { return lhs.sign == rhs.sign && lhs.term == rhs.term; }

  friend bool operator!=(SignedTerm const& lhs, SignedTerm const& rhs)
  { return !(lhs == rhs); }

  friend bool operator<(SignedTerm const& l, SignedTerm const& r)
  { return std::tie(l.term, l.sign) < std::tie(r.term, r.sign); }

  friend bool operator> (SignedTerm const& l, SignedTerm const& r) { return r < l; }
  friend bool operator<=(SignedTerm const& l, SignedTerm const& r) { return l == r || l < r; }
  friend bool operator>=(SignedTerm const& l, SignedTerm const& r) { return l == r || l > r; }
};

// represents a tuple of a numeral (1/k), and a multiset { ti | i in I } of signed terms
// with the intended semantics that the term that has been normalized to this
// sigmaNf is equivalent to 1/k * ( t1 + t2 + ... + tn )
struct SigmaNf {
  unsigned k;
  MultiSet<SignedTerm> sum;
  SigmaNf(unsigned k, MultiSet<SignedTerm> sum) : k(k), sum(std::move(sum)) {}
  friend std::ostream& operator<<(std::ostream& out, SigmaNf const& self)
  { return out << "(1 / " << self.k << ") " << self.sum; }

  friend bool operator==(SigmaNf const& lhs, SigmaNf const& rhs)
  { return lhs.k == rhs.k && lhs.sum == rhs.sum; }

  friend bool operator!=(SigmaNf const& lhs, SigmaNf const& rhs)
  { return !(lhs == rhs); }
};


class QKbo 
  : public Ordering
{
public:
  CLASS_NAME(QKbo);
  USE_ALLOCATOR(QKbo);

  QKbo(QKbo&& kbo) = default;
  QKbo& operator=(QKbo&& kbo) = default;
  explicit QKbo(Precedence);

  QKbo(Problem& prb, const Options& opts) 
    : QKbo(Precedence(prb,opts)) {}

  virtual ~QKbo() {}

  Result compare(Literal* l1, Literal* l2) const override;
  Result compare(TermList, TermList) const final override;

  void show(ostream& out) const final override;

  Comparison compareFunctors(unsigned fun1, unsigned fun2) const override { ASSERTION_VIOLATION }

  void setState(std::shared_ptr<IrcState> s) { _shared = std::move(s); }

private:
  template<class NumTraits>
  SigmaNf rmNum(std::tuple<unsigned, Perfect<Polynom<NumTraits>>> t) const
  {
    auto counts =  std::get<1>(t)->iterSummands()
          .map([](auto s) {
            auto count  = Int::safeAbs(ifOfType<IntegerConstantType>(s.numeral, 
                           [&](IntegerConstantType num) { return num; },
                            /* decltype(num) in { RatTraits, RealTraits } */
                           [&](auto num) { 
                             ASS_EQ(num.denominator(), IntegerConstantType(1))
                             return num.numerator();
                           }).toInner());
            if (count == 0) {
              ASS(s.numeral.sign() == Sign::Zero)
              count = 1;
            }
            SignedTerm term = { 
              .sign = s.numeral.sign(),
              .term = s.factors->denormalize(),
            };
            return std::make_tuple(term, count);
          })
          .template collect<Stack>();
    std::sort(counts.begin(), counts.end(), [](auto& l, auto& r) { return std::get<0>(l) < std::get<0>(r); });
    return SigmaNf(std::get<0>(t), MultiSet<SignedTerm>::fromSortedStack(std::move(counts)));
  }

  std::tuple<unsigned, Perfect<Polynom<IntTraits>>> divNf(Perfect<Polynom<IntTraits>> t) const
  { return std::make_tuple(1, t); }

  template<class NumTraits>
  std::tuple<unsigned, Perfect<Polynom<NumTraits>>> divNf(Perfect<Polynom<NumTraits>> t) const
  {
    auto l = t->iterSummands()
      .map([](auto s) { return s.numeral.denominator(); })
      .fold(IntegerConstantType(1), [&](auto acc, auto next) 
               { return IntegerConstantType::lcm(acc, next); });
    return std::make_tuple(Int::safeAbs(l.toInner()), typename NumTraits::ConstantType(l, IntegerConstantType(1)) * t);
  }

  auto asClosure() const 
  { return [this](auto const& l, auto const& r) { return this->compare(l, r); }; }

#ifdef VDEBUG
public:
#endif
  template<class NumTraits>
  SigmaNf sigmaNf(TermList t) const
  {
    // TODO: without normalizing?
    auto norm = _shared->normalize(TypedTermList(t, NumTraits::sort())).template wrapPoly<NumTraits>();
    return rmNum(divNf(norm));
  }

  template<class NumTraits>
  MultiSet<TermList> absEq(Literal* l) const
  {
    ASS(l->isEquality())
    // TODO this code was never executed yet!!
    ASSERTION_VIOLATION
  }

 
  // template<class NumTraits> 
  // auto iterAtomsT(TermList t) const
  // {
  //   // TODO: without normalizing?
  //   auto nf = sigmaNf(t);
  //   return nf.sum.iter()
  //     .map([](auto ))
  //   // auto norm = _shared->normalize(TypedTermList(t, NumTraits::sort())) .template wrapPoly<NumTraits>();
  //   // return norm->iterSummands()
  //   //   .map([&](auto const& monom) -> TermList 
  //   //       { return monom.factors->denormalize(); });
  // }
 
  // template<class NumTraits> 
  // auto iterAtomsL(Literal* t) const
  // {
  //   return ifElseIter(t->isEquality(), 
  //       [&]() { return concatIters(iterAtomsT<NumTraits>(t->termArg(0)), iterAtomsT<NumTraits>(t->termArg(1))); }
  //     , [&]() { return iterAtomsT<NumTraits>(t->termArg(0)); }
  //   );
  // }


  static constexpr unsigned POS_EQ_LEVEL = 0;
  static constexpr unsigned NEG_EQ_LEVEL = 1;
  bool hasSubstitutionProperty(SigmaNf const& l) const;

  using AtomsStar = std::tuple<MultiSet<SignedTerm>, uint8_t>;

  template<class NumTraits>
  // the signs are not used here, but only there since we reuse code from sigmaNf and don't want to copy the datastructure.
  Option<AtomsStar> atomsStar(Literal* literal) const
  {
    auto mainIdx = literal->isEquality() && literal->termArg(0) == NumTraits::zero() ? 1 : 0;
    ASS_REP(NumTraits::zero() == literal->termArg(1 - mainIdx), literal);
    auto level = literal->isEquality() && literal->isPositive() ? POS_EQ_LEVEL : NEG_EQ_LEVEL;
    auto nf = sigmaNf<NumTraits>(literal->termArg(mainIdx));
    if (hasSubstitutionProperty(nf)) {
      return Option<AtomsStar>(std::make_tuple(std::move(nf.sum), level));
    } else {
      return Option<AtomsStar>();
    }
  }

 
  Option<AtomsStar> atomsStar(Literal* literal) const
  {
    using Out = Option<AtomsStar>;
    return tryNumTraits([&](auto numTraits) {
      using NumTraits = decltype(numTraits);
      if (NumTraits::sort() != SortHelper::getTermArgSort(literal, 0)) {
        return Option<Out>();
      } else {
        return Option<Out>(atomsStar<NumTraits>(literal));
      }
    }) || [&]() -> Out {
      ASS(literal->isEquality())
      auto level = literal->isPositive() ? POS_EQ_LEVEL : NEG_EQ_LEVEL;
      auto multiset = MultiSet<SignedTerm>({ 
        // the sign is only a dummy here to match the type of atomsStar<NumTraits>
          SignedTerm::zero(literal->termArg(0)), 
          SignedTerm::zero(literal->termArg(1)), 
        });
      return Option<AtomsStar>(std::make_tuple(std::move(multiset), level));
    };
  } 

private:
  Result compare(SigmaNf l, SigmaNf r) const;

  using FlatSum = Stack<std::tuple<Option<TermList>, RationalConstantType>>;
  Ordering::Result cmpSum(FlatSum const& l, FlatSum const& r) const;
  FlatSum flatWithCoeffs(TermList t) const;
  Result cmpNonAbstr(TermList, TermList) const;
  Option<TermList> abstr(TermList) const;

  Precedence _prec;
  std::shared_ptr<IrcState> _shared;
  KBO _kbo;
};

} // namespace Kernel

#endif // __QKBO__