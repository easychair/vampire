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
#include "ProofExtra.hpp"
#include "Kernel/Ordering.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Indexing/LiteralIndex.hpp"

namespace Inferences
{

using namespace Kernel;
using namespace Indexing;
using namespace Saturation;

class BinaryResolution
: public GeneratingInferenceEngine
{
public:
  BinaryResolution() 
    : _index(0)
  {  }

  void attach(SaturationAlgorithm* salg);
  void detach();

  ClauseIterator generateClauses(Clause* premise);

private:
  Clause* generateClause(
    Clause* queryCl, Literal* queryLit, Clause* resultCl, Literal* resultLit,
    ResultSubstitutionSP subs, AbstractingUnifier* absUnif);

  BinaryResolutionIndex* _index;
};

using BinaryResolutionExtra = TwoLiteralInferenceExtra;

};

#endif /*__BinaryResolution__*/
