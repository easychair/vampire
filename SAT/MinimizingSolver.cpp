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
 * @file MinimizingSolver.cpp
 * Implements class MinimizingSolver.
 */

#include "SAT/SATClause.hpp"

#include "MinimizingSolver.hpp"

namespace SAT
{

MinimizingSolver::MinimizingSolver(SATSolver* inner)
 : _varCnt(0), _inner(inner), _assignmentValid(false), _heap(CntComparator(_unsClCnt))
{
}

void MinimizingSolver::ensureVarCount(unsigned newVarCnt)
{
  if (newVarCnt<= _varCnt) {
    return;
  }
  _varCnt = newVarCnt;
  _inner->ensureVarCount(newVarCnt);
  _asgn.expand(newVarCnt+1);
  _explicit.expand(newVarCnt+1);
  _polarity.expand(newVarCnt+1);
  _watcher.expand(newVarCnt+1);
  _unsClCnt.expand(newVarCnt+1, 0);
  _heap.elMap().expand(newVarCnt+1);
  _clIdx.expand(newVarCnt+1);
  _assignmentValid = false;
}

void MinimizingSolver::addClause(SATClause* cl)
{
  // pass it to inner ...
  _inner->addClause(cl);
  _assignmentValid = false;
  
  // ... and also keep track for minimization
  if (cl->length()!=0) {
    //we need to filter out the empty clause -- it won't have any influence on our algorithm
    //(as it will make the problem unsat and we process only satisfiable assignment), but it
    //is a corner case that needs to be handled

    _unprocessed.push(cl);
  }
}

void MinimizingSolver::addClauseIgnoredInPartialModel(SATClause* cl)
{
  // just passing to _inner, but for minimization it will be ignored
  _inner->addClause(cl);
  _assignmentValid = false;
}

SATSolver::Status MinimizingSolver::solve(unsigned conflictCountLimit) 
{
  _assignmentValid = false;
  return _inner->solve(conflictCountLimit);
}

SATSolver::VarAssignment MinimizingSolver::getAssignment(unsigned var)
{
  ASS_G(var,0); ASS_LE(var,_varCnt);

  if(!_assignmentValid) {
    updateAssignment();
  }

  if(admitsDontcare(var)) {
    if(_explicit[var])
      return _polarity[var] ? SATSolver::TRUE : SATSolver::FALSE;
    return SATSolver::DONT_CARE;
  }
  return _asgn[var] ? SATSolver::TRUE : SATSolver::FALSE;
}

bool MinimizingSolver::isZeroImplied(unsigned var)
{
  ASS_G(var,0); ASS_LE(var,_varCnt);

  bool res = _inner->isZeroImplied(var);
  ASS(!res || getAssignment(var)!=DONT_CARE); //zero-implied variables will not become a don't care
  return res;
}

/**
 * Give a concrete value (as opposed to don't-care) to the given variable.
 */
void MinimizingSolver::selectVariable(unsigned var)
{
  ASS_G(var,0); ASS_LE(var,_varCnt);
  ASS_G(_unsClCnt[var],0);
  
  SATClauseStack& satisfied = _clIdx[var];
  SATClauseStack& watch = _watcher[var];
  while(satisfied.isNonEmpty()) {
    SATClause* cl = satisfied.pop();
    if(!_satisfiedClauses.insert(cl)) {
      continue;
    }
    watch.push(cl);
    auto cit = cl->iter();
    while(cit.hasNext()) {
      SATLiteral cl_lit = cit.next();
      unsigned cl_var = cl_lit.var(); 
      if (cl_lit.polarity() == _asgn[cl_var]) {
        ASS_G(_unsClCnt[cl_var], 0);
        _unsClCnt[cl_var]--;
        if (cl_var != var) { // var has been just popped
          _heap.notifyIncrease(cl_var); //It was an increase wrt max-heap
        }
      }
    }
  }  
}

void MinimizingSolver::putIntoIndex(SATClause* cl)
{
  auto cit = cl->iter();
  while(cit.hasNext()) {
    SATLiteral lit = cit.next();
    unsigned var = lit.var();

    if (lit.polarity() == _asgn[var]) {
      _clIdx[var].push(cl);
      _unsClCnt[var]++;
    }
  }
}

bool MinimizingSolver::tryPuttingToAnExistingWatch(SATClause* cl)
{
  auto cit = cl->iter();
  while(cit.hasNext()) {
    SATLiteral lit = cit.next();
    unsigned var = lit.var();
    if(_asgn[var]==lit.polarity() && !admitsDontcare(var)) {
      ALWAYS(_satisfiedClauses.insert(cl));
      _watcher[var].push(cl);
      return true;
    }
  }
  return false;
}

/**
 * Move satisfied unprocessed clauses into an appropriate watch, and
 * unsatisfied unprocessed clauses into _clIdx
 */
void MinimizingSolver::processUnprocessedAndFillHeap()
{
  while(_unprocessed.isNonEmpty()) {
    SATClause* cl = _unprocessed.pop();
    ASS_G(cl->length(),0)
    if(!tryPuttingToAnExistingWatch(cl)) {
      putIntoIndex(cl);
    }
  }
  
  for(unsigned var=1; var<=_varCnt; var++) {
    ASS(!_heap.contains(var));
    if(_unsClCnt[var]>0) {            
      _heap.addToEnd(var);      
    }
  }
  _heap.heapify();    
}

/**
 * Update the values in _asgn and move the clauses whose watch
 * became unsatisfied to _unprocessed.
 */
void MinimizingSolver::processInnerAssignmentChanges()
{
  for(unsigned v=1; v<=_varCnt; v++) {
    VarAssignment va = _inner->getAssignment(v);
    bool changed;
    switch(va) {
    case DONT_CARE:
      changed = false;
      break;
    case TRUE:
      changed = !_asgn[v];
      _asgn[v] = true;
      break;
    case FALSE:
      changed = _asgn[v];
      _asgn[v] = false;
      break;
    case NOT_KNOWN:
    default:
      ASSERTION_VIOLATION;
      break;
    }

    if(changed) {
      SATClauseStack& watch = _watcher[v];
      _unprocessed.loadFromIterator(SATClauseStack::Iterator(watch));
      _satisfiedClauses.removeIteratorElements(SATClauseStack::Iterator(watch));
      watch.reset();
    }
  }
}

void MinimizingSolver::updateAssignment()
{
  TIME_TRACE("minimizing solver time");
  
  processInnerAssignmentChanges();
  processUnprocessedAndFillHeap();

  while (!_heap.isEmpty()) {
    unsigned best_var = _heap.pop();
    if (_unsClCnt[best_var] > 0) {
      selectVariable(best_var);
      ASS_EQ(_unsClCnt[best_var],0);
      ASS(_clIdx[best_var].isEmpty());
    } else {
      _clIdx[best_var].reset();
    }
  }
  
  _assignmentValid = true;
}

}
