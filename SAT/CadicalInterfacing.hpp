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
 * @file CadicalInterfacing.hpp
 * Defines class CadicalInterfacing
 */
#ifndef __CadicalInterfacing__
#define __CadicalInterfacing__

#include "SATSolver.hpp"
#include "SATLiteral.hpp"
#include "SATClause.hpp"

#include <cadical.hpp>

namespace SAT{

class CadicalInterfacing : public PrimitiveProofRecordingSATSolver
{
public: 
	CadicalInterfacing(const Shell::Options& opts, bool generateProofs=false);

  /**
   * Can be called only when all assumptions are retracted
   *
   * A requirement is that in a clause, each variable occurs at most once.
   */
  virtual void addClause(SATClause* cl) override;
  
  /**
   * Opportunity to perform in-processing of the clause database.
   *
   * (Minisat deletes unconditionally satisfied clauses.)
   */
  virtual void simplify() override {
    _solver.simplify();
  }

  virtual Status solve(unsigned conflictCountLimit) override;
  
  /**
   * If status is @c SATISFIABLE, return assignment of variable @c var
   */
  virtual VarAssignment getAssignment(unsigned var) override;

  /**
   * If status is @c SATISFIABLE, return 0 if the assignment of @c var is
   * implied only by unit propagation (i.e. does not depend on any decisions)
   */
  virtual bool isZeroImplied(unsigned var) override;
  /**
   * Collect zero-implied literals.
   *
   * Can be used in SATISFIABLE and UNKNOWN state.
   *
   * @see isZeroImplied()
   */
  virtual void collectZeroImplied(SATLiteralStack& acc) override;
  /**
   * Return a valid clause that contains the zero-implied literal
   * and possibly the assumptions that implied it. Return 0 if @c var
   * was an assumption itself.
   * If called on a proof producing solver, the clause will have
   * a proper proof history.
   */
  virtual SATClause* getZeroImpliedCertificate(unsigned var) override;

  virtual void ensureVarCount(unsigned newVarCnt) override;

  virtual unsigned newVar() override;

  virtual void suggestPolarity(unsigned var, unsigned pol) override {
    _solver.phase(vampire2Cadical(pol, var));
  }

  /**
   * Add an assumption into the solver.
   */
  virtual void addAssumption(SATLiteral lit) override;

  virtual void retractAllAssumptions() override {
    _assumptions.clear();
    _status = Status::UNKNOWN;
  };
  
  virtual bool hasAssumptions() const override {
    return (_assumptions.size() > 0);
  };

  Status solveUnderAssumptions(const SATLiteralStack& assumps, unsigned conflictCountLimit) override;

  /**
   * Use minisat and solving under assumptions to minimize the given set of premises (= unsat core extraction).
   *
   * Assumes @b premises in conjunction with @b assumps unsat.
   * Returns a "small" subset of premises which is still unsat under assumps.
   */
  static SATClauseList* minimizePremiseList(SATClauseList* premises, SATLiteralStack& assumps);

  /**
   * Assuming that @b first together with @b second is inconsistent,
   * produce (in @b result) a set of clauses over the signature of @b first,
   * such that @b second |= @b result and
   * @b first together with @b result is also inconsistent.
   */
  static void interpolateViaAssumptions(unsigned maxVar, const SATClauseStack& first, const SATClauseStack& second, SATClauseStack& result);

protected:
  void solveModuloAssumptionsAndSetStatus(unsigned conflictCountLimit = UINT_MAX);

private:
  static int vampire2Cadical(bool polarity, unsigned atom) {
    ASS_NEQ(atom, 0)
    return polarity ? atom : -(int)(atom);
  }

  static int vampire2Cadical(SATLiteral vampire) {
    return vampire2Cadical(vampire.polarity(), vampire.var());
  }

  static SATLiteral cadical2Vampire(int cadical) {
    return SATLiteral(std::abs(cadical), cadical < 0);
  }

  Status _status;
  std::vector<int> _assumptions;
  CaDiCaL::Solver _solver;
};

}//end SAT namespace

 #endif /*CadicalInterfacing*/
