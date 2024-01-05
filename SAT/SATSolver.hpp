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
 * @file SATSolver.hpp
 * Defines class SATSolver.
 */

#ifndef __SATSolver__
#define __SATSolver__

#include "SATLiteral.hpp"
#include "SATInference.hpp"

#include <climits>

namespace SAT {

class SATSolver {
public:
  enum VarAssignment {
    TRUE,
    FALSE,
    DONT_CARE,  // to represent partial models
    NOT_KNOWN
  };
  enum Status {
    SATISFIABLE,
    UNSATISFIABLE,
    /**
     * Solving for just a bounded number of conflicts may return UNKNOWN.
     **/
    UNKNOWN
  };

  virtual ~SATSolver() {}
  
  /**
   * Add a clause to the solver.
   *
   * In the clause each variable occurs at most once. (No repeated literals, no tautologies.)
   *
   * Memory-wise, the clause is NOT owned by the solver. Yet it shouldn't be destroyed while the solver lives.
   * TODO: This is not ideal and should be addressed! (reference counting?)
   */
  virtual void addClause(SATClause* cl) = 0;

  void addClausesIter(SATClauseIterator cit) {
    while (cit.hasNext()) {
      addClause(cit.next());
    }
  }

  /**
   * When a solver supports partial models (via DONT_CARE values in the assignment),
   * a partial model P computed must satisfy all the clauses added to the solver
   * via addClause and there must be a full model extending P which also
   * satisfies clauses added via addClauseIgnoredInPartialModel.
   * 
   * This is a default implementation of addClauseIgnoredInPartialModel
   * for all the solvers which return total models
   * for which addClause and addClauseIgnoredInPartialModel
   * naturally coincide.
   */  
  virtual void addClauseIgnoredInPartialModel(SATClause* cl) { addClause(cl); }
  
  /**
   * Opportunity to perform in-processing of the clause database.
   */
  virtual void simplify() {}

  /**
   * Establish Status of the clause set inserted so far.
   * 
   * If conflictCountLimit==0,
   * do only unit propagation, if conflictCountLimit==UINT_MAX, do
   * full satisfiability check, and for values in between, restrict
   * the number of conflicts, and in case it is reached, stop with
   * solving and assign the status to UNKNOWN.
   */
  virtual Status solve(unsigned conflictCountLimit) = 0;
  
  Status solve(bool onlyPropagate=false) { return solve(onlyPropagate ? 0u : UINT_MAX); }
    
  /**
   * If status is @c SATISFIABLE, return assignment of variable @c var
   */
  virtual VarAssignment getAssignment(unsigned var) = 0;

  /**
   * If status is @c SATISFIABLE, return true if the assignment of @c var is
   * implied only by unit propagation (i.e. does not depend on any decisions)
   */
  virtual bool isZeroImplied(unsigned var) = 0;
  /**
   * Collect zero-implied literals.
   *
   * Can be used in SATISFIABLE and UNKNOWN state.
   *
   * @see isZeroImplied()
   */
  virtual void collectZeroImplied(SATLiteralStack& acc) = 0;
  /**
   * Return a valid clause that contains the zero-implied literal
   * and possibly the assumptions that implied it. Return 0 if @c var
   * was an assumption itself.
   * If called on a proof producing solver, the clause will have
   * a proper proof history.
   */
  virtual SATClause* getZeroImpliedCertificate(unsigned var) = 0;

  /**
   * Ensure that clauses mentioning variables 1..newVarCnt can be handled.
   * 
   * See also newVar for a different (and conceptually incompatible)
   * way for managing variables in the solver.
   */
  virtual void ensureVarCount(unsigned newVarCnt) {}
  
  /**
   * Allocate a slot for a new (previosly unused) variable in the solver
   * and return the variable. 
   * 
   * Variables start from 1 and keep increasing by 1.
   */
  virtual unsigned newVar() = 0;

  virtual void suggestPolarity(unsigned var, unsigned pol) = 0;
  virtual void setPolarity(unsigned var, unsigned pol) = 0;

  /**
   * Suggest random polarity for variables up to maxVar (inclusive),
   * so that the next call to solver will tend to produce
   * a different model (provided the status will be satisfiable).
   */
  virtual void randomizeForNextAssignment(unsigned maxVar) {
    for (unsigned var=1; var <= maxVar; var++) {
      suggestPolarity(var,Random::getBit());
    }
  }

  /**
   * Immediately after a call to solveXXX that returned UNSAT,
   * this method can be used to obtain the corresponding
   * empty SATClause as a root of a corresponding refutation tree.
   * 
   * (However, the empty clause may be invalidated later on.)
   */    
  virtual SATClause* getRefutation() = 0;

  /**
   * Under the same conditions as getRefutation
   * a solver may return a list of SAT clauses which
   * where shown unsatisfiable
   * (possibly under additional assumptions the caller keeps track of themselves).
   *
   * A solver may ignore to implement this function
   * and will return an empty list instead.
   */
  virtual SATClauseList* getRefutationPremiseList() { return SATClauseList::empty(); }

  /**
   * If status is @c SATISFIABLE, return assignment of variable @c var
   */
  bool trueInAssignment(SATLiteral lit)
  {
    VarAssignment asgn = getAssignment(lit.var());
    VarAssignment desired = lit.polarity() ? TRUE : FALSE;
    return asgn==desired;
  }

  /**
   * If status is @c SATISFIABLE, return assignment of variable @c var
   */
  bool falseInAssignment(SATLiteral lit)
  {
    VarAssignment asgn = getAssignment(lit.var());
    VarAssignment desired = lit.polarity() ? FALSE: TRUE;
    return asgn==desired;
  }  
};

inline std::ostream& operator<<(std::ostream& out, SATSolver::Status const& s)
{ 
  switch (s)  {
    case SATSolver::SATISFIABLE: return out << "SATISFIABLE";
    case SATSolver::UNSATISFIABLE: return out << "UNSATISFIABLE";
    case SATSolver::UNKNOWN: return out << "UNKNOWN";
    default: ASSERTION_VIOLATION; return  out << "INVALID STATUS";
  }
}

class SATSolverWithAssumptions: 
      public SATSolver {
public:

  // The first three functions represent the original TWLSolver-style of interface ...

  /**
   * Add assumption to the solver. Perhaps process partially,
   * but mainly ensure the assumption is considered during the next calls to solve()
   */
  virtual void addAssumption(SATLiteral lit) = 0;

  /**
   * Retract all the assumptions added so far.
   * The solver becomes assumption-free.
   *
   * Note: this may destroy the model computed during a previous call to solve
   * (as it currently does in TWL).
   */
  virtual void retractAllAssumptions() = 0;

  /**
   * Test whether any assumptions are currently registered by the solver.
   */
  virtual bool hasAssumptions() const = 0;

  // ... a better alternative could be the interface below.

  // It is currently implemented in terms of the one above
  // (but may be overridden in SATSolver implementations).
  // Note the the below interface allows the solver
  // to identify a subset of the given assumptions that were only needed for
  // previous contradiction to be derived.

  // Note also, however, that solveUnderAssumptions cannot guarantee
  // access to the model after it returns SATISFIABLE, as it uses retractAllAssumptions
  // to clean in the end.

  virtual Status solveUnderAssumptions(const SATLiteralStack& assumps, unsigned conflictCountLimit, bool onlyProperSubusets) {
    ASS(!hasAssumptions());
    _failedAssumptionBuffer.reset();

    Status res = solve(conflictCountLimit);
    if (res == UNSATISFIABLE) {
      return res;
    }

    SATLiteralStack::ConstIterator it(assumps);
    while (it.hasNext()) {
      SATLiteral lit = it.next();
      addAssumption(lit);
      _failedAssumptionBuffer.push(lit);

      res = solve(conflictCountLimit);
      if (res == UNSATISFIABLE) {
        break;
      }
    }

    retractAllAssumptions();
    return res;
  }

  /**
   * Solve under the given set of assumptions @b assumps.
   * If UNSATISFIABLE is returned, a subsequent call to
   * failedAssumptions() returns a subset of these
   * that are already sufficient for the unsatisfiability.
   *
   * @b onlyPropagate suggests that a limited (and potentially incomplete)
   * solving strategy should be employed which only performs unit propagation.
   *
   * If @b onlyProperSubusets, time can be saved by
   * skipping the case when all the given assumptions
   * would need to be considered to obtain unsatisfiability
   * and UNKOWN can be returned instead right away.
   */
  Status solveUnderAssumptions(const SATLiteralStack& assumps, bool onlyPropagate=false, bool onlyProperSubusets=false) {
    return solveUnderAssumptions(assumps,onlyPropagate ? 0u : UINT_MAX,onlyProperSubusets);
  }

  /**
   * When solveUnderAssumptions(assumps) returns UNSATISFIABLE,
   * failedAssumptions contain a subset of assumps that is sufficient
   * for this UNSATISFIABLE status.
   */
  virtual const SATLiteralStack& failedAssumptions() {
    return _failedAssumptionBuffer;
  }

  /**
   * Apply fixpoint minimization to already obtained failed assumption set
   * and return the result (as failedAssumptions).
   */
  const SATLiteralStack& explicitlyMinimizedFailedAssumptions(bool onlyPropagate=false, bool randomize = false) {
    return explicitlyMinimizedFailedAssumptions(onlyPropagate ? 0u : UINT_MAX, randomize);
  }

  virtual const SATLiteralStack& explicitlyMinimizedFailedAssumptions(unsigned conflictCountLimit, bool randomize) {
    // assumes solveUnderAssumptions(...,conflictCountLimit,...) just returned UNSAT and initialized _failedAssumptionBuffer

    ASS(!hasAssumptions());

    unsigned sz = _failedAssumptionBuffer.size();

    if (!sz) { // a special case. Because of the bloody unsigned (see below)!
      return _failedAssumptionBuffer;
    }

    if (randomize) {
      // randomly permute the content of _failedAssumptionBuffer
      // not to bias minimization from one side or another
      for(unsigned i=sz-1; i>0; i--) {
        unsigned tgtPos=Random::getInteger(i+1);
        std::swap(_failedAssumptionBuffer[i], _failedAssumptionBuffer[tgtPos]);
      }
    }

    unsigned i = 0;
    while (i < sz) {
      // load all but i-th
      for (unsigned j = 0; j < sz; j++) {
        if (j != i) {
          addAssumption(_failedAssumptionBuffer[j]);
        }
      }

      if (solve(conflictCountLimit) == UNSATISFIABLE) {
        // leave out forever by overwriting by the last one (buffer shrinks implicitly)
        _failedAssumptionBuffer[i] = _failedAssumptionBuffer[--sz];
      } else {
        // move on
        i++;
      }

      retractAllAssumptions();
    }

    _failedAssumptionBuffer.truncate(sz);
    return _failedAssumptionBuffer;
  }

protected:
  SATLiteralStack _failedAssumptionBuffer;
};

/**
 * A convenience class for solvers which do not track actual refutations
 * and so return the whole set of clauses added so far as refutations.
 * 
 * This need not necessarily inherit from SATSolverWithAssumptions,
 * but why bother with multiple inheritance if we know the only 
 * two descendants of this class will need it...
 */
class PrimitiveProofRecordingSATSolver : public SATSolverWithAssumptions {
public:
  PrimitiveProofRecordingSATSolver() :  
    _addedClauses(0), _refutation(new(0) SATClause(0)), _refutationInference(new PropInference(SATClauseList::empty()))
    {
      _refutation->setInference(_refutationInference);    
    }
  
  virtual ~PrimitiveProofRecordingSATSolver() {
    // cannot clear the list - some inferences may be keeping its suffices till proof printing phase ...
    // _addedClauses->destroy(); // we clear the list but not its content
  }

  virtual void addClause(SATClause* cl) override 
  {
    SATClauseList::push(cl,_addedClauses);
  }
  
  virtual SATClause* getRefutation() override
  {
    // connect the added clauses ... 
    SATClauseList* prems = _addedClauses;  

    // ... with the current assumptions

    // TODO: the assumption set will be empty after a call to solveUnderAssumptions()
    // This does not matter much since refutations are only ever passed to collectFOPremises
    // and there are no FO premises of assumption inferences

    // So the below is commented out to prevent useless leaking

    /*
    for (size_t i=0; i < _assumptions.size(); i++) {
      SATClause* unit = new(1) SATClause(1);
      (*unit)[0] = _assumptions[i];
      unit->setInference(new AssumptionInference());
      SATClauseList::push(unit,prems);
    }
    */

    _refutationInference->setPremises(prems);

    return _refutation; 
  }
  
  virtual SATClauseList* getRefutationPremiseList() override {
    return _addedClauses;
  }

private:
  // to be used for the premises of a refutation
  SATClauseList* _addedClauses;
  
  /**
   * Empty clause to be returned by the getRefutation call.
   * Recycled between consecutive getRefutation calls.
   */
  SATClause* _refutation;
  /**
   * The inference inside _refutation.
   */
  PropInference* _refutationInference;  
};


}

#endif // __SATSolver__
