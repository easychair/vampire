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
 * @file MinimizingSolver.hpp
 * 
 * Defines class MinimizingSolver which supports partial models.
 */

#ifndef __MinimizingSolver__
#define __MinimizingSolver__

#include "Forwards.hpp"

#include "Lib/DArray.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/DHSet.hpp"
#include "Lib/ScopedPtr.hpp"
#include "Lib/Stack.hpp"
#include "Lib/DynamicHeap.hpp"
#include "Lib/ArrayMap.hpp"

#include "SATSolver.hpp"



namespace SAT {

using namespace Lib;

class MinimizingSolver : public SATSolver {
public:
  MinimizingSolver(SATSolver* inner);

  virtual SATClause* getRefutation() override { return _inner->getRefutation(); }
  virtual SATClauseList* getRefutationPremiseList() override {
    return _inner->getRefutationPremiseList();
  }
  virtual void randomizeForNextAssignment(unsigned maxVar) override {
    _inner->randomizeForNextAssignment(maxVar); _assignmentValid = false;
  }

  virtual void addClause(SATClause* cl) override;
  virtual void addClauseIgnoredInPartialModel(SATClause* cl) override;
  virtual Status solve(unsigned conflictCountLimit) override;
  
  virtual VarAssignment getAssignment(unsigned var) override;
  virtual bool isZeroImplied(unsigned var) override;
  virtual void collectZeroImplied(SATLiteralStack& acc) override { _inner->collectZeroImplied(acc); }
  virtual SATClause* getZeroImpliedCertificate(unsigned var) override { return _inner->getZeroImpliedCertificate(var); }

  virtual void ensureVarCount(unsigned newVarCnt) override;

  virtual unsigned newVar() override {
    DEBUG_CODE(unsigned oldVC = _varCnt);
    ensureVarCount(_varCnt+1);
    ASS_EQ(_varCnt,oldVC+1);
    return _varCnt;
  }

  virtual void suggestPolarity(unsigned var, unsigned pol) override { _inner->suggestPolarity(var,pol); }
  virtual void setPolarity(unsigned var, unsigned pol) override {
    _inner->setPolarity(var, pol);
    _explicit[var] = true;
    _polarity[var] = pol;
  }

private:
  bool admitsDontcare(unsigned var) { 
    ASS_G(var,0); ASS_LE(var,_varCnt);

    return _watcher[var].isEmpty() && !_inner->isZeroImplied(var);
    
    /**
     * TODO: as an optimization, the _watcher stack for zero implied variables
     * could be reset. It will not be needed anymore.
     */
  }
  
  void selectVariable(unsigned var);

  bool tryPuttingToAnExistingWatch(SATClause* cl);
  void putIntoIndex(SATClause* cl);

  void processInnerAssignmentChanges();
  void processUnprocessedAndFillHeap();
  void updateAssignment();

  unsigned _varCnt;
  ScopedPtr<SATSolver> _inner;

  /**
   * If true, _asgn assignment corresponds to the assignment in
   * the inner solver
   */
  bool _assignmentValid;

  /**
   * Clauses of which we yet need to ensure they are satisfied
   *
   * Invariant: outside of this object when _assignmentValid, the stack is empty.
   */
  SATClauseStack _unprocessed;

  /**
   * A total extension of the current assignment. 
   * The current assignment will report don't-care for those variables
   * for which admitsDontcare is true.
   * 
   * We define a literal "corresponding to variable var in _asgn"
   * as SATLiteral(var,_asgn[var]). Used below.
   */
  DArray<bool> _asgn;

  // explicitly set by setPolarity()
  DArray<bool> _explicit;
  DArray<bool> _polarity;

  /**
   * Array of clauses made satisfied by giving up the don't-care value
   * for a particular variable and using the value dictated by _asgn instead.   
   *
   * The length of the array is _varCnt.
   */
  DArray<SATClauseStack> _watcher;

  typedef DArray<unsigned> CntArray;
  
  /**
   * Number of unsatisfied clauses for each literal
   * corresponding to a variable in _asgn.
   * 
   * Indexed by literal's var.
   *
   * The length of the array is _varCnt.
   */
  CntArray _unsClCnt;

  struct CntComparator
  {
    CntComparator(CntArray& ctr) : _ctr(ctr) {}

    Comparison compare(unsigned v1, unsigned v2)
    {
      // DynamicHeap is minimal and we want maximum, 
      // so we need to swap the arguments
      return Int::compare(_ctr[v2], _ctr[v1]);
    }
    CntArray& _ctr;
  };
  
  /**
   * Heap "on top of" _unsClCnt to facilitate fast "extract max" operation.
   * 
   * Heap is empty when _assignmentValid.
   */  
  DynamicHeap<unsigned, CntComparator, ArrayMap<size_t> > _heap;
  
  /**
   * Not yet satisfied clauses indexed by each variable
   * whose corresponding literal in _asgn can make the clause true.
   *
   * A clause appears once for every such literal.
   * All entries should be empty when _assignmentValid.
   *
   * The length of the array is _varCnt.
   */
  DArray<SATClauseStack> _clIdx;

  /**
   * A set of satisfied clauses. To correctly maintain
   * _unsClCnt, when there is more than one way to make clause
   * satisfied.
   */  
  DHSet<SATClause*> _satisfiedClauses;
  
};

}

#endif // __MinimizingSolver__
