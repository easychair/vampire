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
 * @file ForwardDemodulation.cpp
 * Implements class ForwardDemodulation.
 */

#include "Debug/RuntimeStatistics.hpp"

#include "Lib/DHSet.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/Metaiterators.hpp"
#include "Debug/TimeProfiling.hpp"
#include "Lib/Timer.hpp"
#include "Lib/VirtualIterator.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Ordering.hpp"
#include "Kernel/Renaming.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/TermIterators.hpp"
#include "Kernel/ColorHelper.hpp"
#include "Kernel/RobSubstitution.hpp"

#include "Indexing/Index.hpp"
#include "Indexing/IndexManager.hpp"
#include "Indexing/TermIndex.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"

#include "DemodulationHelper.hpp"

#include "ForwardDemodulation.hpp"

namespace Inferences {

using namespace Lib;
using namespace Kernel;
using namespace Indexing;
using namespace Saturation;

namespace {

struct Applicator : SubstApplicator {
  Applicator(ResultSubstitution* subst) : subst(subst) {}
  TermList operator()(unsigned v) const override {
    return subst->applyToBoundResult(v);
  }
  ResultSubstitution* subst;
};

struct ApplicatorWithEqSort : SubstApplicator {
  ApplicatorWithEqSort(ResultSubstitution* subst, const RobSubstitution& vSubst) : subst(subst), vSubst(vSubst) {}
  TermList operator()(unsigned v) const override {
    return vSubst.apply(subst->applyToBoundResult(v), 0);
  }
  ResultSubstitution* subst;
  const RobSubstitution& vSubst;
};

} // end namespace

void ForwardDemodulation::attach(SaturationAlgorithm* salg)
{
  ForwardSimplificationEngine::attach(salg);
  _index=static_cast<DemodulationLHSIndex*>(
	  _salg->getIndexManager()->request(DEMODULATION_LHS_CODE_TREE) );

  auto opt = getOptions();
  _preorderedOnly = opt.forwardDemodulation()==Options::Demodulation::PREORDERED;
  _encompassing = opt.demodulationRedundancyCheck()==Options::DemodulationRedundancyCheck::ENCOMPASS;
  _precompiledComparison = opt.demodulationPrecompiledComparison();
  _skipNonequationalLiterals = opt.demodulationOnlyEquational();
  _helper = DemodulationHelper(opt, &_salg->getOrdering());
}

void ForwardDemodulation::detach()
{
  _index=0;
  _salg->getIndexManager()->release(DEMODULATION_LHS_CODE_TREE);
  ForwardSimplificationEngine::detach();
}

template <bool combinatorySupSupport>
bool ForwardDemodulationImpl<combinatorySupSupport>::perform(Clause* cl, Clause*& replacement, ClauseIterator& premises)
{
  TIME_TRACE("forward demodulation");

  Ordering& ordering = _salg->getOrdering();

  //Perhaps it might be a good idea to try to
  //replace subterms in some special order, like
  //the heaviest first...

  static DHSet<TermList> attempted;
  attempted.reset();

  unsigned cLen=cl->length();
  for(unsigned li=0;li<cLen;li++) {
    Literal* lit=(*cl)[li];
    if (lit->isAnswerLiteral()) {
      continue;
    }
    if (_skipNonequationalLiterals && !lit->isEquality()) {
      continue;
    }
    typename std::conditional<!combinatorySupSupport,
      NonVariableNonTypeIterator,
      FirstOrderSubtermIt>::type it(lit);
    while(it.hasNext()) {
      TypedTermList trm = it.next();
      if(!attempted.insert(trm)) {
        //We have already tried to demodulate the term @b trm and did not
        //succeed (otherwise we would have returned from the function).
        //If we have tried the term @b trm, we must have tried to
        //demodulate also its subterms, so we can skip them too.
        it.right();
        continue;
      }

      bool redundancyCheck = _helper.redundancyCheckNeededForPremise(cl, lit, trm);

      auto git = _index->getGeneralizations(trm.term(), /* retrieveSubstitutions */ true);
      while(git.hasNext()) {
        auto qr=git.next();
        ASS_EQ(qr.data->clause->length(),1);

        if(!ColorHelper::compatible(cl->color(), qr.data->clause->color())) {
          continue;
        }

        auto lhs = qr.data->term;

        // TODO:
        // to deal with polymorphic matching
        // Ideally, we would like to extend the substitution
        // returned by the index to carry out the sort match.
        // However, ForwardDemodulation uses a CodeTree as its
        // indexing mechanism, and it is not clear how to extend
        // the substitution returned by a code tree.
        static RobSubstitution eqSortSubs;
        if(lhs.isVar()){
          eqSortSubs.reset();
          TermList querySort = trm.sort();
          TermList eqSort = qr.data->term.sort();
          if(!eqSortSubs.match(eqSort, 0, querySort, 1)){
            continue;
          }
        }

        TermList rhs = qr.data->rhs;
        bool preordered = qr.data->preordered;

        auto subs = qr.unifier;
        ASS(subs->isIdentityOnQueryWhenResultBound());

        ApplicatorWithEqSort applWithEqSort(subs.ptr(), eqSortSubs);
        Applicator applWithoutEqSort(subs.ptr());
        auto appl = lhs.isVar() ? (SubstApplicator*)&applWithEqSort : (SubstApplicator*)&applWithoutEqSort;

        if (_precompiledComparison) {
          qr.data->comparator->reset();
          if (!preordered && (_preorderedOnly || !qr.data->comparator->next(appl))) {
#if DEBUG_ORDERING
            if (ordering.isGreaterOrEq(AppliedTerm(trm),AppliedTerm(rhs,appl,true))==Ordering::GREATER) {
              INVALID_OPERATION("greater");
            }
#endif
            continue;
          }
#if DEBUG_ORDERING
          if (ordering.isGreaterOrEq(AppliedTerm(trm),AppliedTerm(rhs,appl,true))!=Ordering::GREATER) {
            INVALID_OPERATION("not greater");
          }
#endif
        } else {
          if (!preordered && (_preorderedOnly || ordering.isGreaterOrEq(AppliedTerm(trm),AppliedTerm(rhs,appl,true))!=Ordering::GREATER)) {
            continue;
          }
        }

        // encompassing demodulation is fine when rewriting the smaller guy
        if (redundancyCheck && _encompassing) {
          // this will only run at most once;
          // could have been factored out of the getGeneralizations loop,
          // but then it would run exactly once there
          Ordering::Result litOrder = ordering.getEqualityArgumentOrder(lit);
          if ((trm==*lit->nthArgument(0) && litOrder == Ordering::LESS) ||
              (trm==*lit->nthArgument(1) && litOrder == Ordering::GREATER)) {
            redundancyCheck = false;
          }
        }

        TermList rhsS = subs->applyToBoundResult(rhs);
        if (lhs.isVar()) {
          rhsS = eqSortSubs.apply(rhsS, 0);
        }

        if (redundancyCheck && !_helper.isPremiseRedundant(cl, lit, trm, rhsS, lhs, appl)) {
          continue;
        }

        Literal* resLit = EqHelper::replace(lit,trm,rhsS);
        if(EqHelper::isEqTautology(resLit)) {
          env.statistics->forwardDemodulationsToEqTaut++;
          premises = pvi( getSingletonIterator(qr.data->clause));
          return true;
        }

        RStack<Literal*> resLits;
        resLits->push(resLit);

        for(unsigned i=0;i<cLen;i++) {
          Literal* curr=(*cl)[i];
          if(curr!=lit) {
            resLits->push(curr);
          }
        }

        env.statistics->forwardDemodulations++;

        premises = pvi( getSingletonIterator(qr.data->clause));
        replacement = Clause::fromStack(*resLits, SimplifyingInference2(InferenceRule::FORWARD_DEMODULATION, cl, qr.data->clause));
        if(env.options->proofExtra() == Options::ProofExtra::FULL)
          env.proofExtra.insert(replacement, new ForwardDemodulationExtra(lhs, trm));
        return true;
      }
    }
  }

  return false;
}

// This is necessary for templates defined in cpp files.
// We are happy to do it for ForwardDemodulationImpl, since it (at the moment) has only two specializations:
template class ForwardDemodulationImpl<false>;
template class ForwardDemodulationImpl<true>;

}
