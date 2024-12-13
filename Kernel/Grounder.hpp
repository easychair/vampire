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
 * @file Grounder.hpp
 * Defines class Grounder.
 */

#ifndef __Grounder__
#define __Grounder__

#include "Forwards.hpp"

#include "Lib/DHMap.hpp"

namespace Kernel {

using namespace Lib;
using namespace SAT;

class GlobalSubsumptionGrounder {
  struct OrderNormalizingComparator;

public:
  GlobalSubsumptionGrounder(SATSolver &satSolver) : _satSolver(satSolver) {}

  /**
   * Return SATClause that is a result of grounding of the
   * non-propositional part of @c cl.
   *
   * The order of literals in @c cl is preserved.
   *
   * @param cl the clause
   * @param acc previously accumulated literals
   */
  void groundNonProp(Clause* cl, SATLiteralStack& acc);

private:
  /**
   * Normalize literals before grounding.
   *
   * The order of literals in @c lits must be preserved.
   */
  void normalize(unsigned cnt, Literal** lits);

  /**
   * Return SATLiteral corresponding to @c lit.
   */
  SATLiteral groundLiteral(Literal* lit);

  /**
   * Return SATLiteral corresponding to @c lit.
   */
  SATLiteral groundNormalized(Literal*);

  /** Map from positive literals to SAT variable numbers */
  DHMap<Literal*, unsigned> _asgn;

  /** Reference to a SATSolver instance for which the grounded clauses
   * are being prepared. Used to request new variables from the Solver.
   *
   * Also used to communicate source literals with IGGrounder. */
  SATSolver &_satSolver;

};

}

#endif // __Grounder__
