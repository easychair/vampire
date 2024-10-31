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
 * @file CodeTreeForwardSubsumptionAndResolution.cpp
 * Implements class CodeTreeForwardSubsumptionAndResolution.
 */

#include "Saturation/SaturationAlgorithm.hpp"
#include "Shell/Statistics.hpp"

#include "ProofExtra.hpp"
#include "CodeTreeForwardSubsumptionAndResolution.hpp"

namespace Inferences {

void CodeTreeForwardSubsumptionAndResolution::attach(SaturationAlgorithm *salg)
{
  ForwardSimplificationEngine::attach(salg);
  auto index = static_cast<CodeTreeSubsumptionIndex*>(
    _salg->getIndexManager()->request(FW_SUBSUMPTION_CODE_TREE));
  _ct = index->getClauseCodeTree();
}

void CodeTreeForwardSubsumptionAndResolution::detach()
{
  _ct = nullptr;
  _salg->getIndexManager()->release(FW_SUBSUMPTION_CODE_TREE);
  ForwardSimplificationEngine::detach();
}

bool CodeTreeForwardSubsumptionAndResolution::perform(Clause *cl, Clause *&replacement, ClauseIterator &premises)
{
  if (_ct->isEmpty()) {
    return false;
  }

  static ClauseCodeTree::ClauseMatcher cm;

  cm.init(_ct, cl, _subsumptionResolution);

  Clause* premise;
  int resolvedQueryLit;

  while ((premise = cm.next(resolvedQueryLit))) {
    if (resolvedQueryLit == -1) {
      ASS(satSubs.checkSubsumption(premise, cl));
      premises = pvi(getSingletonIterator(premise));
      env.statistics->forwardSubsumed++;
      cm.reset();
      return true;
    }
    ASS(satSubs.checkSubsumptionResolutionWithLiteral(premise, cl, resolvedQueryLit));

    LiteralStack res;
    for (unsigned i = 0; i < cl->length(); i++) {
      if (i == (unsigned)resolvedQueryLit) {
        continue;
      }
      res.push((*cl)[i]);
    }
    replacement = Clause::fromStack(res, SimplifyingInference2(InferenceRule::SUBSUMPTION_RESOLUTION, cl, premise));
    if(env.options->proofExtra() == Options::ProofExtra::FULL)
      env.proofExtra.insert(replacement, new LiteralInferenceExtra((*cl)[resolvedQueryLit]));
    premises = pvi(getSingletonIterator(premise));
    env.statistics->forwardSubsumptionResolution++;
    cm.reset();
    return true;
  }

  cm.reset();
  return false;
}

} // namespace Inferences
