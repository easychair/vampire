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
 * @file InequalityResolution.hpp
 * Defines class InequalityResolution
 *
 */

#ifndef __IRC_InequalityResolution__
#define __IRC_InequalityResolution__

#include "Forwards.hpp"

#include "Inferences/InferenceEngine.hpp"
#include "Kernel/Ordering.hpp"
#include "Shell/UnificationWithAbstractionConfig.hpp"
#include "Indexing/InequalityResolutionIndex.hpp"
#include "Shell/Options.hpp"

namespace Inferences {
namespace IRC {

using namespace Kernel;
using namespace Indexing;
using namespace Saturation;

class InequalityResolution
: public GeneratingInferenceEngine
{
public:
  CLASS_NAME(InequalityResolution);
  USE_ALLOCATOR(InequalityResolution);

  InequalityResolution(InequalityResolution&&) = default;
  InequalityResolution(shared_ptr<IrcState> shared) 
    : _shared(std::move(shared))
    , _index(nullptr)
  {  }

  class Lhs : public SelectedSummand { 
  public: 
    Lhs(SelectedSummand s) : SelectedSummand(std::move(s)) {} 
    Lhs(Lhs&&) = default;
    Lhs& operator=(Lhs&&) = default;
  };

  class Rhs : public SelectedSummand { 
  public: 
    Rhs(SelectedSummand s) : SelectedSummand(std::move(s)) {} 
    Rhs(Rhs&&) = default;
    Rhs& operator=(Rhs&&) = default;
  };

  auto iterLhs(Clause* cl) const 
  { return _shared->selectedSummands(cl, /* strictly max literal*/ true, /* strictly max term */ true)
            .filter([&](auto const& selected) { return selected.sign() == Sign::Pos; })
            .map([&]   (auto selected)        { return Lhs(std::move(selected)); }); }

  auto iterRhs(Clause* cl) const 
  { return _shared->selectedSummands(cl, /* strictly max literal*/ true, /* strictly max term */ false)
            .filter([&](auto const& selected) { return selected.sign() == Sign::Neg; })
            .map([&]   (auto selected)        { return Rhs(std::move(selected)); }); }

  void attach(SaturationAlgorithm* salg) final override;
  void detach() final override;

  ClauseIterator generateClauses(Clause* premise) final override;

#if VDEBUG
  virtual void setTestIndices(Stack<Indexing::Index*> const&) final override;
#endif

private:

  Option<Clause*> applyRule(
      Lhs const& lhs, unsigned lhsVarBank,
      Rhs const& rhs, unsigned rhsVarBank,
      UwaResult& uwa
      ) const;

  template<class NumTraits, class Subst, class CnstIter> Option<Clause*> applyRule(
      Clause* hyp1, Literal* lit1, IrcLiteral<NumTraits> l1, Monom<NumTraits> j_s1,
      Clause* hyp2, Literal* lit2, IrcLiteral<NumTraits> l2, Monom<NumTraits> k_s2,
      Subst sigma, CnstIter cnst, unsigned nCnst
      ) const;

  template<class NumTraits> ClauseIterator generateClauses(Clause* clause, Literal* lit, IrcLiteral<NumTraits> l1, Monom<NumTraits> j_s1) const;

  shared_ptr<IrcState> _shared;
  InequalityResolutionIndex* _index;
  InequalityResolutionIndex* _lhsIndex;
  InequalityResolutionIndex* _rhsIndex;
};

} // namespace IRC 
} // namespace Inferences 


#endif /*__IRC_InequalityResolution__*/