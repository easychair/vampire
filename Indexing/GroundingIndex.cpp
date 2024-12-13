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
 * @file GroundingIndex.cpp
 * Implements class GroundingIndex.
 */

#include "GroundingIndex.hpp"

#include "Lib/SharedSet.hpp"

#include "Kernel/Grounder.hpp"
#include "Kernel/Inference.hpp"

#include "Shell/Options.hpp"

#include "SAT/MinisatInterfacing.hpp"
#include "SAT/BufferedSolver.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

namespace Indexing
{

using namespace std;

GroundingIndex::GroundingIndex(const Options& opt)
{
  _solver = new MinisatInterfacing(opt,true);
  _grounder = new Kernel::GlobalSubsumptionGrounder(*_solver);
}

void GroundingIndex::handleClause(Clause* c, bool adding)
{
  //We are adding clauses into the index when performing the subsumption check
}

}
