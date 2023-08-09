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
 * @file BinaryResolution.hpp
 * Defines class BinaryResolution
 *
 */

#ifndef __BinaryResolution__
#define __BinaryResolution__

#include "Forwards.hpp"

#include "InferenceEngine.hpp"
#include "Kernel/Ordering.hpp"
#include "Shell/UnificationWithAbstractionConfig.hpp"
#include "Indexing/TermIndex.hpp"

namespace Inferences
{

using namespace Kernel;
using namespace Indexing;
using namespace Saturation;

class BinaryResolution
: public GeneratingInferenceEngine
{
public:
  CLASS_NAME(BinaryResolution);
  USE_ALLOCATOR(BinaryResolution);

  BinaryResolution() 
    : _index(0),
    _unificationWithAbstraction(false)
  {  }

  void attach(SaturationAlgorithm* salg);
  void detach();

  static Clause* generateClause(Clause* queryCl, Literal* queryLit, SLQueryResult res, const Options& opts, const Ordering& ord, PassiveClauseContainer* passive=0, bool afterCheck=false, LiteralSelector* ls = 0, DemodulationLHSIndex* demLhsIndex = 0);
  ClauseIterator generateClauses(Clause* premise);

private:
  struct UnificationsFn;
  struct ResultFn;

  BinaryResolutionIndex* _index;
  bool _unificationWithAbstraction;
  DemodulationLHSIndex* _demLhsIndex;
};

};

#endif /*__BinaryResolution__*/
