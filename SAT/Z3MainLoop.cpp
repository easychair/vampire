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
 * @file Z3MainLoop.cpp
 * Defines class Z3MainLoop.
 */

#if VZ3

#include "Forwards.hpp"

#include "Shell/Statistics.hpp"

#include "Z3MainLoop.hpp"

namespace SAT
{

using namespace Shell;
using namespace Lib;

Z3MainLoop::Z3MainLoop(Problem& prb, const Options& opt)
: MainLoop(prb,opt)
{
}

void Z3MainLoop::init()
{
}

MainLoopResult Z3MainLoop::runImpl()
{
  if(!_prb.getProperty()->allNonTheoryClausesGround()){
    return MainLoopResult(Statistics::INAPPROPRIATE);
  }

  SAT2FO s2f;
  Z3Interfacing solver(_opt,s2f, /* unsat core */ false, /* export smtlib problem */ "", Shell::Options::ProblemExportSyntax::SMTLIB);

 ClauseIterator cit(_prb.clauseIterator());
 while(cit.hasNext()){
   Clause* cl = cit.next();

   if(cl->varCnt() > 0){
     ASS(cl->isTheoryAxiom());
     continue;
   }

   unsigned len = cl->size();
   SATClause* sc = new(len) SATClause(len);
   unsigned i=0;
   for (auto l : cl->iterLits()) {
     SATLiteral sl = s2f.toSAT(l);
     (*sc)[i++] = sl;
   }
   solver.addClause(sc);
 }

 SATSolver::Status status = solver.solve(UINT_MAX);

 Statistics::TerminationReason reason = Statistics::UNKNOWN; 
                
 if(status == SATSolver::Status::UNSATISFIABLE){ reason = Statistics::REFUTATION; }
 if(status == SATSolver::Status::SATISFIABLE){ reason = Statistics::SATISFIABLE; }

 return MainLoopResult(reason);
}

}

#endif
