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
 * @file BinaryResolution.cpp
 * Implements class BinaryResolution.
 */

#include "Debug/RuntimeStatistics.hpp"

#include "Indexing/ResultSubstitution.hpp"
#include "Kernel/UnificationWithAbstraction.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/PairUtils.hpp"
#include "Lib/VirtualIterator.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/ColorHelper.hpp"
#include "Kernel/Formula.hpp"
#include "Kernel/Unit.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/LiteralSelector.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/RobSubstitution.hpp"

#include "Indexing/Index.hpp"
#include "Indexing/LiteralIndex.hpp"
#include "Indexing/IndexManager.hpp"
#include "Indexing/SubstitutionTree.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

#include "Shell/AnswerLiteralManager.hpp"
#include "Shell/ConditionalRedundancyHandler.hpp"
#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"

#include "BinaryResolution.hpp"
#define DEBUG_RESOLUTION(lvl, ...) if (lvl < 0) { DBG("resolution: ", __VA_ARGS__) }

namespace Inferences
{

using namespace std;
using namespace Lib;
using namespace Kernel;
using namespace Indexing;
using namespace Saturation;

void BinaryResolution::attach(SaturationAlgorithm* salg)
{
  ASS(!_index);

  GeneratingInferenceEngine::attach(salg);
  _index=static_cast<BinaryResolutionIndex*> (
	  _salg->getIndexManager()->request(BINARY_RESOLUTION_SUBST_TREE) );
}

void BinaryResolution::detach()
{
  ASS(_salg);

  _index=0;
  _salg->getIndexManager()->release(BINARY_RESOLUTION_SUBST_TREE);
  GeneratingInferenceEngine::detach();
}

/**
 * Ordering aftercheck is performed iff ord is not 0,
 * in which case also ls is assumed to be not 0.
 */
template<class ComputeConstraints>
Clause* BinaryResolution::generateClause(Clause* queryCl, Literal* queryLit, Clause* resultCl, Literal* resultLit,
          ResultSubstitutionSP subs, ComputeConstraints computeConstraints, const Options& opts, bool afterCheck,
          PassiveClauseContainer* passiveClauseContainer, Ordering* ord, LiteralSelector* ls,
          ConditionalRedundancyHandler const* condRedHandler)
{
  DEBUG_RESOLUTION(0, "lhs: ", *queryLit, " (clause: ", queryCl->number(), ")")
  DEBUG_RESOLUTION(0, "rhs: ", *resultLit, " (clause: ", resultCl->number(), ")")
  DEBUG_RESOLUTION(0, "subs: ", *subs)
  ASS(resultCl->store()==Clause::ACTIVE);//Added to check that generation only uses active clauses

  if(!ColorHelper::compatible(queryCl->color(),resultCl->color()) ) {
    env.statistics->inferencesSkippedDueToColors++;
    if(opts.showBlocked()) {
      std::cout << "Blocked resolution of " << *queryCl << " and " << * resultCl << endl;
    }
    return 0;
  }

  unsigned clength = queryCl->length();
  unsigned dlength = resultCl->length();

  // LRS-specific optimization:
  // check whether we can conclude that the resulting clause will be discarded by LRS since it does not fulfil the age/weight limits (in which case we can discard the clause)
  // we already know the age here so we can immediately conclude whether the clause fulfils the age limit
  // since we have not built the clause yet we compute lower bounds on the weight of the clause after each step and recheck whether the weight-limit can still be fulfilled.
  unsigned wlb=0;//weight lower bound
  unsigned numPositiveLiteralsLowerBound = // lower bound on number of positive literals, don't know at this point whether duplicate positive literals will occur
      Int::max(queryLit->isPositive() ?  queryCl->numPositiveLiterals()-1 :  queryCl->numPositiveLiterals(),
              resultLit->isPositive() ? resultCl->numPositiveLiterals()-1 : resultCl->numPositiveLiterals());

  auto constraints = computeConstraints();
  auto nConstraints = constraints->size();
  Inference inf(GeneratingInference2(nConstraints == 0 ?  InferenceRule::RESOLUTION : InferenceRule::CONSTRAINED_RESOLUTION, queryCl, resultCl));
  Inference::Destroyer inf_destroyer(inf); // will call destroy on inf when coming out of scope unless disabled

  bool andThatsIt = false;
  bool hasAgeLimitStrike = passiveClauseContainer && passiveClauseContainer->mayBeAbleToDiscriminateClausesUnderConstructionOnLimits()
                                                  && passiveClauseContainer->exceedsAgeLimit(numPositiveLiteralsLowerBound, inf, andThatsIt);
  if (hasAgeLimitStrike && andThatsIt) { // we are dealing with purely age-limited container (no need for weight-related investigations)
    RSTAT_CTR_INC("binary resolutions skipped for (pure) age limit before building clause");
    env.statistics->discardedNonRedundantClauses++;
    return 0;
  }

  if(hasAgeLimitStrike) {
    for(unsigned i=0;i<clength;i++) {
      Literal* curr=(*queryCl)[i];
      if(curr!=queryLit) {
        wlb+=curr->weight();
      }
    }
    for(unsigned i=0;i<dlength;i++) {
      Literal* curr=(*resultCl)[i];
      if(curr!=resultLit) {
        wlb+=curr->weight();
      }
    }
    if(passiveClauseContainer->exceedsWeightLimit(wlb, numPositiveLiteralsLowerBound, inf)) {
      RSTAT_CTR_INC("binary resolutions skipped for weight limit before building clause");
      env.statistics->discardedNonRedundantClauses++;
      return 0;
    }
  }

  bool synthesis = (env.options->questionAnswering() == Options::QuestionAnsweringMode::SYNTHESIS);
  Literal* cAnsLit = synthesis ? queryCl->getAnswerLiteral() : nullptr;
  Literal* dAnsLit = synthesis ? resultCl->getAnswerLiteral() : nullptr;
  bool bothHaveAnsLit = (cAnsLit != nullptr) && (dAnsLit != nullptr);

  RStack<Literal*> resLits;

  Literal* queryLitAfter = 0;
  if (afterCheck && queryCl->numSelected() > 1) {
    TIME_TRACE(TimeTrace::LITERAL_ORDER_AFTERCHECK);
    queryLitAfter = subs->applyToQuery(queryLit);
  }

  resLits->loadFromIterator(constraints->iterFifo());
  for(unsigned i=0;i<clength;i++) {
    Literal* curr=(*queryCl)[i];
    if(curr!=queryLit && (!bothHaveAnsLit || curr!=cAnsLit)) {
      Literal* newLit = subs->applyToQuery(curr);
      if(hasAgeLimitStrike) {
        wlb+=newLit->weight() - curr->weight();
        if(passiveClauseContainer->exceedsWeightLimit(wlb, numPositiveLiteralsLowerBound, inf)) {
          RSTAT_CTR_INC("binary resolutions skipped for weight limit while building clause");
          env.statistics->discardedNonRedundantClauses++;
          return nullptr;
        }
      }
      if (queryLitAfter && i < queryCl->numSelected()) {
        TIME_TRACE(TimeTrace::LITERAL_ORDER_AFTERCHECK);

        Ordering::Result o = ord->compare(newLit,queryLitAfter);

        if (o == Ordering::GREATER ||
            (ls->isPositiveForSelection(newLit)    // strict maximimality for positive literals
                && o == Ordering::EQUAL)) {
          env.statistics->inferencesBlockedForOrderingAftercheck++;
          return nullptr;
        }
      }
      resLits->push(newLit);
    }
  }

  Literal* qrLitAfter = 0;
  if (afterCheck && resultCl->numSelected() > 1) {
    TIME_TRACE(TimeTrace::LITERAL_ORDER_AFTERCHECK);
    qrLitAfter = subs->applyToResult(resultLit);
  }

  for(unsigned i=0;i<dlength;i++) {
    Literal* curr=(*resultCl)[i];
    if(curr!=resultLit && (!bothHaveAnsLit || curr!=dAnsLit)) {
      Literal* newLit = subs->applyToResult(curr);
      if(hasAgeLimitStrike) {
        wlb+=newLit->weight() - curr->weight();
        if(passiveClauseContainer->exceedsWeightLimit(wlb, numPositiveLiteralsLowerBound, inf)) {
          RSTAT_CTR_INC("binary resolutions skipped for weight limit while building clause");
          env.statistics->discardedNonRedundantClauses++;
          return nullptr;
        }
      }
      if (qrLitAfter && i < resultCl->numSelected()) {
        TIME_TRACE(TimeTrace::LITERAL_ORDER_AFTERCHECK);

        Ordering::Result o = ord->compare(newLit,qrLitAfter);

        if (o == Ordering::GREATER ||
            (ls->isPositiveForSelection(newLit)   // strict maximimality for positive literals
                && o == Ordering::EQUAL)) {
          env.statistics->inferencesBlockedForOrderingAftercheck++;
          return nullptr;
        }
      }
      resLits->push(newLit);
    } 
  }

  if (nConstraints == 0 && condRedHandler) {
    if (!condRedHandler->handleResolution(queryCl, queryLit, resultCl, resultLit, subs.ptr())) {
      return 0;
    }
  }

   if (bothHaveAnsLit) {
     Literal* newLitC = subs->applyToQuery(cAnsLit);
     Literal* newLitD = subs->applyToResult(dAnsLit);
     bool cNeg = queryLit->isNegative();
     Literal* condLit = cNeg ? subs->applyToResult(resultLit) : subs->applyToQuery(queryLit);
     resLits->push(SynthesisALManager::getInstance()->makeITEAnswerLiteral(condLit, cNeg ? newLitC : newLitD, cNeg ? newLitD : newLitC));
   }

  if(nConstraints != 0){
    env.statistics->cResolution++;
  }
  else{
    env.statistics->resolution++;
  }

  inf_destroyer.disable(); // ownership passed to the the clause below
  Clause *cl = Clause::fromStack(*resLits, inf);
  if(env.options->proofExtra() == Options::ProofExtra::FULL)
    env.proofExtra.insert(cl, new BinaryResolutionExtra(queryLit, resultLit));

  return cl;

}


Clause* BinaryResolution::generateClause(Clause* queryCl, Literal* queryLit,
                                         Clause* resultCl, Literal* resultLit,
                                         AbstractingUnifier& uwa, const Options& opts, SaturationAlgorithm* salg) {
  // perform binary resolution on query results
  auto subs = ResultSubstitution::fromSubstitution(&uwa.subs(), RetrievalAlgorithms::DefaultVarBanks::query, RetrievalAlgorithms::DefaultVarBanks::internal);
  bool doAfterCheck = opts.literalMaximalityAftercheck() && salg->getLiteralSelector().isBGComplete();
  return BinaryResolution::generateClause(queryCl, queryLit, resultCl, resultLit, subs,
      [&](){ return uwa.computeConstraintLiterals(); },
      opts, doAfterCheck, salg->getPassiveClauseContainer(),
      &salg->getOrdering(), &salg->getLiteralSelector(), &salg->condRedHandler());
}

ClauseIterator BinaryResolution::generateClauses(Clause* premise)
{
  return pvi(TIME_TRACE_ITER("resolution",
      premise->getSelectedLiteralIterator()
      // TODO filter out >= in alasca
        .filter([](auto l) { return !l->isEquality(); })
        .flatMap([this,premise](auto lit) {
            // find query results for literal `lit`
            return iterTraits(_index->getUwa(lit, /* complementary */ true,
                                             env.options->unificationWithAbstraction(),
                                             env.options->unificationWithAbstractionFixedPointIteration()))
                     .map([this,lit,premise](auto qr) { return BinaryResolution::generateClause(premise, lit, qr.data->clause, qr.data->literal, *qr.unifier, this->getOptions(), _salg); });
        })
        .filter(NonzeroFn())
  ));
}

}
