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
 * @file CadicalInterfacing.cpp
 * Implements class CadicalInterfacing
 */

#include "CadicalInterfacing.hpp"

namespace SAT
{

using namespace Shell;
using namespace Lib;

CadicalInterfacing::CadicalInterfacing(const Shell::Options& opts, bool generateProofs):
  _status(Status::SATISFIABLE)
{
  // TODO: consider tuning minisat's options to be set for _solver
  // (or even forwarding them to vampire's options)  
}
  
SATSolver::Status CadicalInterfacing::solveUnderAssumptions(const SATLiteralStack& assumps, unsigned conflictCountLimit)
{
  ASS(!hasAssumptions());

  // load assumptions:
  SATLiteralStack::ConstIterator it(assumps);
  while (it.hasNext()) {
    _assumptions.push_back(vampire2Cadical(it.next()));
  }

  solveModuloAssumptionsAndSetStatus(conflictCountLimit);

  if (_status == Status::UNSATISFIABLE) {
    // unload minisat's internal conflict clause to _failedAssumptionBuffer
    _failedAssumptionBuffer.reset();
    for(int i : _assumptions)
      if(_solver.failed(i))
        _failedAssumptionBuffer.push(cadical2Vampire(i).opposite());
  }

  _assumptions.clear();
  return _status;
}

/**
 * Solve modulo assumptions and set status.
 * @b conflictCountLimit as with addAssumption.
 */
void CadicalInterfacing::solveModuloAssumptionsAndSetStatus(unsigned conflictCountLimit)
{
  _solver.limit("conflicts", conflictCountLimit);
  for(int assumption : _assumptions)
    _solver.assume(assumption);
  int res = _solver.solve();

  if (res == CaDiCaL::SATISFIABLE) {
    _status = Status::SATISFIABLE;
  } else if (res == CaDiCaL::UNSATISFIABLE) {
    _status = Status::UNSATISFIABLE;
  } else {
    _status = Status::UNKNOWN;
  }
}

/**
 * Add clause into the solver.
 *
 */
void CadicalInterfacing::addClause(SATClause* cl)
{
  // store to later generate the refutation
  PrimitiveProofRecordingSATSolver::addClause(cl);
  // TODO: consider measuring time
  ASS_EQ(_assumptions.size(),0);

  unsigned clen=cl->length();
  for(unsigned i=0;i<clen;i++) {
    SATLiteral l = (*cl)[i];
    _solver.add(vampire2Cadical(l));
  }
  _solver.add(0);
}

/**
 * Perform solving and return status.
 */
SATSolver::Status CadicalInterfacing::solve(unsigned conflictCountLimit)
{
  solveModuloAssumptionsAndSetStatus(conflictCountLimit);
  return _status;
}

void CadicalInterfacing::addAssumption(SATLiteral lit)
{
  _assumptions.push_back(vampire2Cadical(lit));
}

SATSolver::VarAssignment CadicalInterfacing::getAssignment(unsigned var)
{
	ASS_EQ(_status, Status::SATISFIABLE);
	ASS_G(var,0); ASS_L(var,_next);

  if(var > _solver.vars())
    return VarAssignment::DONT_CARE;
  int phase = _solver.val(var);
  return phase > 0 ? VarAssignment::TRUE : VarAssignment::FALSE;
}

bool CadicalInterfacing::isZeroImplied(unsigned var)
{
  ASS_G(var,0); ASS_L(var, _next);
  return _solver.fixed(var);
}

void CadicalInterfacing::collectZeroImplied(SATLiteralStack& acc)
{
  for(int i = 1, val; i <= _solver.vars(); i++)
    if((val = _solver.fixed(i)))
      acc.push(cadical2Vampire(val > 0 ? i : -i));
}

SATClause* CadicalInterfacing::getZeroImpliedCertificate(unsigned)
{
  // Currently unused anyway.
  /* The whole SATSolver interface should be revised before
   implementing functions like this one properly */
  NOT_IMPLEMENTED;
}

SATClauseList* CadicalInterfacing::minimizePremiseList(SATClauseList* premises, SATLiteralStack& assumps)
{
  CaDiCaL::Solver solver;

  // TODO this should be a vector
  static DHMap<int,SATClause*> var2prem;
  var2prem.reset();

  static std::vector<int> ass; // assumptions for the final call
  ass.clear();

  int cl_no = 0;

  SATClauseList* it= premises;
  while(it) {
    // store the link for fast lookup
    var2prem.insert(cl_no,it->head());

    // corresponding assumption
    ass.push_back(cl_no + 1); // posive as the assumption

    // allocate the var for the clause
    solver.reserve(cl_no + 1);

    cl_no++;
    it=it->tail();
  }

  // from now on, offset will mark the translation of premises' original variables to the ones in solver here
  int offset = cl_no;

  // start counting from 0 and traversing from the beginning again
  cl_no = 0;
  it = premises;
  while(it) {
    SATClause* cl = it->head();

    // translate the clause to minisat's language (shift vars by offset)
    unsigned clen=cl->length();
    for(unsigned i=0;i<clen;i++) {
      SATLiteral l = (*cl)[i];
      int var = offset + l.var();

      // make sure vars are allocated
      solver.reserve(var);
      solver.add(vampire2Cadical(l.isPositive(), var));
    }

    solver.add(-(cl_no + 1));
    solver.add(0);

    cl_no++;
    it=it->tail();
  }

  // add assumptions from assumps
  SATLiteralStack::Iterator ait(assumps);
  while (ait.hasNext()) {
    SATLiteral l = ait.next();
    int var = offset + l.var();

    ass.push_back(vampire2Cadical(l.isPositive(), var));
  }

  // solve
  for(int assumption : ass)
    solver.assume(assumption);
  ALWAYS(solver.solve() == CaDiCaL::UNSATISFIABLE); // should be unsat

  SATClauseList* result = SATClauseList::empty();

  // extract the used ones
  for (int assumption : ass) {
    if(!solver.failed(assumption))
      continue;

    SATClause* cl;

    if (var2prem.find(assumption - 1,cl)) {
      SATClauseList::push(cl,result);
    } // it could also be one of the "assumps"
  }
  return result;
}

void CadicalInterfacing::interpolateViaAssumptions(unsigned maxVar, const SATClauseStack& first, const SATClauseStack& second, SATClauseStack& result)
{
  NOT_IMPLEMENTED;
  // TODO below could work, just don't know how to test it
  /*
  Minisat::Solver solver_first;
  Minisat::Solver solver_second;

  // TODO: this may be quite wasteful for at least two reasons:
  // 1) 1..maxVar may be a large superset of symb(first \cup second)
  // 2) symb(first \cup second) may be a large superset of symb(second)
  // However, making sure variables are not wasted would require
  // setting up various renamings to maintain correspondence
  // between literals in the solvers
  for(unsigned v = 0; v <= maxVar; v++) { // ... and variable 0 will not be used either
    solver_first.newVar();
    solver_second.newVar();
  }

  DArray<bool> varOfFirst;
  varOfFirst.expand(maxVar+1,false);

  vec<Lit> tmp;

  // load first into solver_first and mark in varOfFirst
  SATClauseStack::ConstIterator it1(first);
  while(it1.hasNext()) {
    SATClause* cl = it1.next();

    unsigned clen=cl->length();
    for(unsigned i=0;i<clen;i++) {
      SATLiteral l = (*cl)[i];
      varOfFirst[l.var()] = true;
      tmp.push(mkLit(l.var(),l.isNegative()));
    }

    solver_first.addClause(tmp);
    tmp.clear();
  }

  // load second into solver_second
  SATClauseStack::ConstIterator it2(second);
  while(it2.hasNext()) {
    SATClause* cl = it2.next();

    unsigned clen=cl->length();
    for(unsigned i=0;i<clen;i++) {
      SATLiteral l = (*cl)[i];
      tmp.push(mkLit(l.var(),l.isNegative()));
    }

    solver_second.addClause(tmp);
    tmp.clear();
  }

  SATLiteralStack vlits;

  // generate the interpolant clauses
  while (solver_first.solve()) {
    // turn model into assumptions for solver_second
    for (int i = 1; i <= (int)maxVar; i++) {
      if (varOfFirst[i]) {
        tmp.push(mkLit(i,solver_first.model[i]==l_False));
      }
    }

    NEVER(solver_second.solve(tmp));
    tmp.clear();

    // conflict is a new clause to put into result and solver_first to look for a different model
    LSet& conflict = solver_second.conflict;
    for (int i = 0; i < conflict.size(); i++) {
      Lit l = conflict[i];

      tmp.push(l);
      vlits.push(SATLiteral(var(l),sign(l) ? 0 : 1));
    }

    solver_first.addClause(tmp);
    tmp.clear();

    result.push(SATClause::fromStack(vlits));
    vlits.reset();
  }
*/
}


} // namespace SAT
