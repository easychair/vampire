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
 * @file FallbackSolverWrapper.hpp
 * Defines class FallbackSolverWrapper.
 *
 * The idea is to run two solvers next to each other 'falling back' to one if the first isn't giving
 * us an answer. The intended setting is where the main solver is a SMT solver that can return Unknown
 * and the fallback solver is a SAT solver that cannot 
 *
 * @author Giles
 */

#ifndef __FallbackSolverWrapper__
#define __FallbackSolverWrapper__

#include "Forwards.hpp"

#include "Lib/DArray.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/DHSet.hpp"
#include "Lib/ScopedPtr.hpp"
#include "Lib/Stack.hpp"

#include "SATSolver.hpp"

namespace SAT {

using namespace Lib;

class FallbackSolverWrapper : public SATSolver {
public:
  FallbackSolverWrapper(SATSolver* inner,SATSolver* fallback);

  virtual SATClause* getRefutation() override { 
    if(_usingFallback){
      return _fallback->getRefutation();
    }
    return _inner->getRefutation(); 
  }
  virtual SATClauseList* getRefutationPremiseList() override {
    if(_usingFallback){
      return _fallback->getRefutationPremiseList();
    }
    return _inner->getRefutationPremiseList();
  }
  virtual void randomizeForNextAssignment(unsigned maxVar) override {
    _fallback->randomizeForNextAssignment(maxVar);
    _inner->randomizeForNextAssignment(maxVar);
  }

  virtual void addClause(SATClause* cl) override;
  virtual Status solve(unsigned conflictCountLimit) override;
  virtual VarAssignment getAssignment(unsigned var) override;

  virtual bool isZeroImplied(unsigned var) override {
    ASS_G(var,0); ASS_LE(var,_varCnt);

    if(_usingFallback){
      return _fallback->isZeroImplied(var);
    }

    // alternatively, we could directly refer to _inner, it must handle variables up to _varCnt as well
    return  _inner->isZeroImplied(var);
  }
  virtual void collectZeroImplied(SATLiteralStack& acc) override { 
    if(_usingFallback){
      _fallback->collectZeroImplied(acc);
      return;
    }
    _inner->collectZeroImplied(acc); 
  }

  virtual SATClause* getZeroImpliedCertificate(unsigned var) override { 
    if(_usingFallback){
      return _fallback->getZeroImpliedCertificate(var);
    }
    return _inner->getZeroImpliedCertificate(var); 
  }

  virtual void ensureVarCount(unsigned newVarCnt) override { 
    _inner->ensureVarCount(newVarCnt); 
    _fallback->ensureVarCount(newVarCnt); 
    _varCnt=std::max(_varCnt,newVarCnt); 
  }


  virtual unsigned newVar() override { 
    ALWAYS(_inner->newVar() == ++_varCnt);
    ALWAYS(_fallback->newVar() == _varCnt);
    return _varCnt;
  }
  
  virtual void suggestPolarity(unsigned var,unsigned pol) override { 
    _inner->suggestPolarity(var,pol); 
    _fallback->suggestPolarity(var,pol); 
  }

  virtual void setPolarity(unsigned var,unsigned pol) override {
    _inner->setPolarity(var,pol);
    _fallback->setPolarity(var,pol);
  }

private:

  ScopedPtr<SATSolver> _inner;
  ScopedPtr<SATSolver> _fallback;

  bool _usingFallback;

  unsigned _varCnt;

};

}

#endif // __FallbackSolverWrapper__
