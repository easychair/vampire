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
 * @file InequalityResolution.cpp
 * Implements class InequalityResolution.
 */

#include "InequalityResolution.hpp"
#include "Saturation/SaturationAlgorithm.hpp"
#include "Shell/Statistics.hpp"

#define DEBUG(...) // DBG(__VA_ARGS__)

namespace Inferences {
namespace IRC {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INDEXING STUFF
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void InequalityResolution::attach(SaturationAlgorithm* salg) 
{
  CALL("InequalityResolution::attach");
  GeneratingInferenceEngine::attach(salg);

  ASS(!_lhsIndex);
  ASS(!_rhsIndex);

  _lhsIndex = static_cast<decltype(_lhsIndex)>(_salg->getIndexManager()->request(IRC_INEQUALITY_RESOLUTION_LHS_SUBST_TREE));
  _rhsIndex = static_cast<decltype(_rhsIndex)>(_salg->getIndexManager()->request(IRC_INEQUALITY_RESOLUTION_RHS_SUBST_TREE));
  _rhsIndex->setShared(_shared);
  _lhsIndex->setShared(_shared);
}

void InequalityResolution::detach() 
{
  CALL("InequalityResolution::detach");
  ASS(_salg);
  GeneratingInferenceEngine::detach();

  // _index=0;
  // _salg->getIndexManager()->release(IRC_INEQUALITY_RESOLUTION_SUBST_TREE);
}

#if VDEBUG
void InequalityResolution::setTestIndices(Stack<Indexing::Index*> const& indices)
{
  _lhsIndex = (decltype(_lhsIndex)) indices[0]; 
  _rhsIndex = (decltype(_rhsIndex)) indices[1]; 
  _lhsIndex->setShared(_shared);
  _rhsIndex->setShared(_shared);
}
#endif

using Lhs = InequalityResolution::Lhs;
using Rhs = InequalityResolution::Rhs;

ClauseIterator InequalityResolution::generateClauses(Clause* premise) 
{
  CALL("InequalityResolution::generateClauses(Clause* premise)")
  ASS(_lhsIndex)
  ASS(_rhsIndex)
  ASS(_shared)
  Stack<Clause*> out;

  for (auto const& lhs : Lhs::iter(*_shared, premise)) {
    DEBUG("lhs: ", lhs)
    for (auto rhs_sigma : _rhsIndex->find(lhs.monom())) {
      auto& rhs   = std::get<0>(rhs_sigma);
      auto& sigma = std::get<1>(rhs_sigma);
      DEBUG("  rhs: ", rhs)
      auto res = applyRule(lhs, 0, rhs, 1, sigma);
      if (res.isSome()) {
        out.push(res.unwrap());
      }
    }
  }

  for (auto const& rhs : Rhs::iter(*_shared, premise)) {
    DEBUG("rhs: ", rhs)
    for (auto lhs_sigma : _lhsIndex->find(rhs.monom())) {
      auto& lhs   = std::get<0>(lhs_sigma);
      auto& sigma = std::get<1>(lhs_sigma);
      DEBUG("  lhs: ", lhs)
      auto res = applyRule(lhs, 1, rhs, 0, sigma);
      if (res.isSome()) {
        out.push(res.unwrap());
      }
    }
  }
  return pvi(ownedArrayishIterator(std::move(out)));
}

// Fourier Motzkin normal:
//
// C₁ \/ +j s₁ + t₁ >₁ 0         C₂ \/ -k s₂ + t₂ >₂ 0 
// --------------------------------------------------
//           (C₁ \/ C₂ \/ k t₁ + j t₂ > 0)σ \/ Cnst
//
// where 
// - (σ, Cnst) = uwa(s₁, s₂)
// - C₁σ /⪯ (+j s₁  + t₁ >₁ 0)σ
// - C₂σ /⪯ (-k s₂ + t₂ >₂ 0)σ
// - s₁, s₂ are not variables
// - {>} ⊆ {>₁,>₂} ⊆ {>,≥}
//
// TODO this condition is not in the theory yet:
// • ( -k s₂ + t₂ >₂ 0 )σ /⪯  ( +j s₁ + t₁ >₁ 0 )σ
//
// Fourier Motzkin tight:
//
// C₁ \/ +j s₁ + t₁ ≥ 0                 C₂ \/ -k s₂ + t₂ ≥ 0 
// --------------------------------------------------------
// (C₁ \/ C₂ \/ k t₁ + j t₂ > 0 \/ -k s₂ + t₂ ≈ 0)σ \/ Cnst
//
// where 
// • (σ, Cnst) = uwa(s₁, s₂)
// • (+j s₁ + t₁ >₁ 0)σ /⪯ C₁σ
// • (-k s₂ + t₂ >₂ 0)σ /⪯ C₂σ
// • s₁σ /⪯ t₁σ 
// • s₂σ /≺ t₂σ 
// • s₁, s₂ are not variables
//
// TODO this condition is not in the theory yet:
// • ( -k s₂ + t₂ >₂ 0 )σ /⪯  ( +j s₁ + t₁ >₁ 0 )σ
//
Option<Clause*> InequalityResolution::applyRule(
    Lhs const& lhs, unsigned lhsVarBank,
    Rhs const& rhs, unsigned rhsVarBank,
    UwaResult& uwa
    ) const 
{
  CALL("InequalityResolution::applyRule")

  return lhs.numTraits().apply([&](auto numTraits) {
    using NumTraits = decltype(numTraits);


    ASS(lhs.sign() == Sign::Pos)
    ASS(rhs.sign() == Sign::Neg)
    ASS(lhs.literal()->functor() == NumTraits::geqF()
     || lhs.literal()->functor() == NumTraits::greaterF())
    ASS(rhs.literal()->functor() == NumTraits::geqF()
     || rhs.literal()->functor() == NumTraits::greaterF())

    bool tight = lhs.literal()->functor() == NumTraits::geqF()
              && rhs.literal()->functor() == NumTraits::geqF();

    Stack<Literal*> out( lhs.clause()->size() - 1 // <- C1
                       + rhs.clause()->size() - 1 // <- C2
                       + 1                        // <- k t₁ + j t₂ > 0
                       + (tight ? 1 : 0)          // <- -k s₂ + t₂ ≈ 0
                       + uwa.cnst().size());      // Cnst


#define check_side_condition(cond, cond_code)                                                       \
    if (!(cond_code)) {                                                                             \
      DEBUG("side condition not fulfiled: " cond)                                                   \
      return Option<Clause*>();                                                                     \
    }                                                                                               \

    ASS(!lhs.monom().isVar() && !rhs.monom().isVar())

    // check_side_condition(
    //     "s₁, s₂ are not variables",
    //     !lhs.monom().isVar() && !rhs.monom().isVar())

    auto L1σ = uwa.sigma(lhs.literal(), rhsVarBank);
    check_side_condition( 
        "(+j s₁ + t₁ >₁ 0)σ /⪯ C₁σ",
        lhs.contextLiterals()
           .all([&](auto L) {
             auto Lσ = uwa.sigma(L, lhsVarBank);
             out.push(Lσ);
             return _shared->notLeq(L1σ, Lσ);
           }));

    auto L2σ = uwa.sigma(rhs.literal(), rhsVarBank);
    check_side_condition(
        "(-k s₂ + t₂ >₂ 0)σ /⪯ C₂σ",
        rhs.contextLiterals()
           .all([&](auto L) {
             auto Lσ = uwa.sigma(L, rhsVarBank);
             out.push(Lσ);
             return _shared->notLeq(L2σ, Lσ);
           }));


    auto s1σ = uwa.sigma(lhs.monom(), lhsVarBank);
    auto s2σ = uwa.sigma(rhs.monom(), rhsVarBank);
    // ASS_REP(_shared->equivalent(sσ.term(), s2σ().term()), make_pair(sσ, s2σ()))
    Stack<TermList> t1σ(rhs.nContextTerms());
    Stack<TermList> t2σ(lhs.nContextTerms());

    check_side_condition(
        "s₁σ /⪯ t₁σ",
        lhs.contextTerms<NumTraits>()
           .all([&](auto ti) {
             auto tiσ = uwa.sigma(ti.factors->denormalize(), lhsVarBank);
             t1σ.push(NumTraits::mulSimpl(ti.numeral, tiσ));
             return _shared->notLeq(s1σ, tiσ);
           }))

    check_side_condition(
        "s₂σ /≺ t₂σ ",
        rhs.contextTerms<NumTraits>()
           .all([&](auto ti) {
             auto tiσ = uwa.sigma(ti.factors->denormalize(), rhsVarBank);
             t2σ.push(NumTraits::mulSimpl(ti.numeral, tiσ));
             return _shared->notLess(s2σ, tiσ);
           }))

    check_side_condition(
        "( -k s₂ + t₂ >₂ 0 )σ /⪯  ( +j s₁ + t₁ >₁ 0 )σ",
        _shared->notLeq(L2σ, L1σ));

    auto j = lhs.numeral().unwrap<typename NumTraits::ConstantType>();
    auto k = rhs.numeral().unwrap<typename NumTraits::ConstantType>().abs();

    auto add = [](auto l, auto r) {
      return l == NumTraits::zero() ? r 
           : r == NumTraits::zero() ? l
           : NumTraits::add(l, r); };

    auto resolventTerm // -> (k t₁ + j t₂)σ
        = add( NumTraits::mulSimpl(k, NumTraits::sum(t1σ.iterFifo())),
               NumTraits::mulSimpl(j, NumTraits::sum(t2σ.iterFifo())));

    if (std::is_same<IntTraits, NumTraits>::value) {
      resolventTerm = add(resolventTerm, NumTraits::constantTl(-1));
    }

    out.push(NumTraits::greater(true, resolventTerm, NumTraits::zero()));

    if (tight) {
      auto rhsSum = // -> (-k s₂ + t₂)σ
        uwa.sigma(rhs.literal(), rhsVarBank)->termArg(0);
      out.push(NumTraits::eq(true, rhsSum, NumTraits::zero()));
    }

    out.loadFromIterator(uwa.cnstLiterals());

    Inference inf(GeneratingInference2(Kernel::InferenceRule::IRC_INEQUALITY_RESOLUTION, lhs.clause(), rhs.clause()));
    auto cl = Clause::fromStack(out, inf);
    DEBUG("out: ", *cl);
    return Option<Clause*>(cl);
  });
}

} // namespace IRC 
} // namespace Inferences 