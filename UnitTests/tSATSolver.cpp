/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#include "Lib/List.hpp"
#include "Lib/Stack.hpp"
#include "Lib/Environment.hpp"

#include "SAT/SATClause.hpp"
#include "SAT/SATLiteral.hpp"
#include "SAT/SATInference.hpp"
#include "SAT/SATSolver.hpp"
#include "SAT/CadicalInterfacing.hpp"
#include "SAT/Z3Interfacing.hpp"

#include "Test/UnitTesting.hpp"

using namespace std;
using namespace Lib;
using namespace SAT;

/**
 * c is either lower or upper case letter, upper
 * case is positive literals, lower negative.
 *
 * SATSolver to use these literals should have varCnt at least 27.
 */
SATLiteral getLit(char c)
{
  bool pol;
  unsigned var;
  if(c>='a' && c<='z') {
    var = c-'a'+1;
    pol = false;
  }
  else {
    ASS_REP(c>='A' && c<='Z',c)
    var = c-'A'+1;
    pol = true;
  }
  return SATLiteral(var, pol);
}

/**
 * spec consists of lower and upper case letters, upper
 * case are positive literals, lower negative.
 *
 * SATSolver to use these clauses should have varCnt at least 27.
 */
SATClause* getClause(const char* spec)
{
  SATLiteralStack lits;
  while(*spec) {
    char c = *spec;
    lits.push(getLit(c));
    spec++;
  }
  return SATClause::fromStack(lits);
}
void ensurePrepared(SATSolver& s)
{
  s.ensureVarCount(27);
}

// TODO this test used to work with TWLSolver
// It doesn't work with Minisat but could
// see Minisat::getZeroImpliedCertificate
/*
void testZICert1(SATSolverWithAssumptions& s)
{
  ensurePrepared(s);
  s.addClause(getClause("ab"));
  s.addClause(getClause("c"));

  unsigned cVar = getLit('c').var();
  ASS(s.isZeroImplied(cVar));
  SATClause* ccert = s.getZeroImpliedCertificate(cVar);
  ASS_EQ(ccert->size(),1);
  ASS_EQ((*ccert)[0],getLit('c'));

  s.solve(); // just so that TWL propagates

  unsigned bVar = getLit('b').var();
  ASS(!s.isZeroImplied(bVar));
  s.addAssumption(getLit('A'));

  s.solve(); // just so that TWL propagates

  ASS(s.isZeroImplied(bVar));
  SATClause* bcert = s.getZeroImpliedCertificate(bVar);
  ASS_EQ(bcert->size(),2);
}

TEST_FUN(satSolverZeroImpliedCert)
{
  MinisatInterfacing s(*env.options,true);
  testZICert1(s);
}
*/

void testProofWithAssumptions(SATSolver& s)
{
  s.ensureVarCount(2);
  s.addClause(getClause("a"));
  s.addClause(getClause("A"));

  ASS_EQ(s.solve(),SATSolver::Status::UNSATISFIABLE);

  SATClause* refutation = s.getRefutation();
  PropInference* inf = static_cast<PropInference*>(refutation->inference());

  // cout << endl << "Refutation: " << refutation->toString() << endl;

  List<SATClause*>* prems = inf->getPremises();

  // cout << "Inference length: " << prems->length() << endl;

  while(prems){
    // cout << prems->head()->toString() << endl;
    prems = prems->tail();
  }

}

TEST_FUN(testProofWithAssums)
{
  CadicalInterfacing s(*env.options,true);
  testProofWithAssumptions(s);
}

void testInterface(SATSolverWithAssumptions &s) {
  ensurePrepared(s);

  ASS_EQ(s.solve(),SATSolver::Status::SATISFIABLE);

  unsigned a = getLit('a').var();
  unsigned b = getLit('b').var();
  unsigned c = getLit('c').var();
  unsigned d = getLit('d').var();

  /*
  s.suggestPolarity(a,0);
  ASS_EQ(s.solve(),SATSolver::SATISFIABLE);
  ASS(s.trueInAssignment(getLit('a')));
  s.suggestPolarity(a,1);
  ASS_EQ(s.solve(),SATSolver::SATISFIABLE);
  ASS(!s.trueInAssignment(getLit('a')));

  s.forcePolarity(a,0);
  ASS_EQ(s.solve(),SATSolver::SATISFIABLE);
  ASS(s.trueInAssignment(getLit('a')));
  s.forcePolarity(a,1);
  ASS_EQ(s.solve(),SATSolver::SATISFIABLE);
  ASS(!s.trueInAssignment(getLit('a')));
  */

  // no longer hard enough to get UNKNOWN
  s.addClause(getClause("ab"));
  ASS_EQ(s.solve(true),SATSolver::Status::SATISFIABLE);
  s.addClause(getClause("aB"));
  ASS_EQ(s.solve(true),SATSolver::Status::SATISFIABLE);
  s.addClause(getClause("Ab"));
  ASS_EQ(s.solve(true),SATSolver::Status::SATISFIABLE);
  s.addClause(getClause("C"));
  ASS_EQ(s.solve(),SATSolver::Status::SATISFIABLE);

  ASS(s.trueInAssignment(getLit('a')));
  ASS(s.trueInAssignment(getLit('b')));
  ASS(s.falseInAssignment(getLit('c')));

  // for a and b depends on learned clauses, which depend on decide polarity
  // but should be both at the same time, or none of the two
  ASS(s.isZeroImplied(a) == s.isZeroImplied(b));

  ASS(s.isZeroImplied(c));
  ASS(!s.isZeroImplied(d));

  cout << " Random: ";
  for (int i = 0; i < 10; i++) {
    s.randomizeForNextAssignment(27);
    s.solve();
    cout << s.trueInAssignment(getLit('d'));
  }
  cout << "  Fixed: ";
  for (int i = 0; i < 10; i++) {
    cout << s.trueInAssignment(getLit('d'));
  }
  cout << endl;

  s.addAssumption(getLit('d'));
  s.addAssumption(getLit('a'));
  ASS(s.hasAssumptions());
  ASS_EQ(s.solve(),SATSolver::Status::SATISFIABLE);
  s.retractAllAssumptions();
  ASS(!s.hasAssumptions());
  // ASS(s.getStatus() == SATSolver::UNKNOWN || s.getStatus() == SATSolver::SATISFIABLE);

  ASS(!s.hasAssumptions());
  s.addAssumption(getLit('A'));
  ASS(s.hasAssumptions());
  ASS_EQ(s.solve(),SATSolver::Status::UNSATISFIABLE);
  s.retractAllAssumptions();
  ASS(!s.hasAssumptions());
  // ASS(s.getStatus() == SATSolver::UNKNOWN || s.getStatus() == SATSolver::SATISFIABLE);

  s.addAssumption(getLit('a'));
  ASS(s.hasAssumptions());
  ASS_EQ(s.solve(),SATSolver::Status::SATISFIABLE);
  s.retractAllAssumptions();
  ASS(!s.hasAssumptions());
  // ASS(s.getStatus() == SATSolver::UNKNOWN || s.getStatus() == SATSolver::SATISFIABLE);

  ASS_EQ(s.solve(),SATSolver::Status::SATISFIABLE);
}

TEST_FUN(testSATSolverInterface)
{ 
  cout << endl << "CaDiCaL" << endl;
  CadicalInterfacing sMini(*env.options,true);
  testInterface(sMini);

  /* Not fully conforming - does not support zeroImplied and resource-limited solving
  cout << endl << "Z3" << endl;
  {
    SAT2FO sat2fo;
    Z3Interfacing sZ3(*env.options,sat2fo);
    testInterface(sZ3);
  }
  */
}

void testAssumptions(SATSolverWithAssumptions &s) {
  ensurePrepared(s);

  s.addClause(getClause("ab"));
  s.addClause(getClause("cde"));

  static SATLiteralStack assumps;
  assumps.reset();

  assumps.push(getLit('X'));
  assumps.push(getLit('A'));
  assumps.push(getLit('B'));
  assumps.push(getLit('C'));
  assumps.push(getLit('D'));
  assumps.push(getLit('E'));
  assumps.push(getLit('Y'));

  ASS_EQ(s.solveUnderAssumptions(assumps),SATSolver::Status::UNSATISFIABLE);

  const SATLiteralStack& failed = s.failedAssumptions();
  for (unsigned i = 0; i < failed.size(); i++) {
    SATLiteral lit = failed[i];
    if (lit.polarity()) {
      cout << (char)('A' + lit.var()-1);
    } else {
      cout << (char)('a' + lit.var()-1);
    }
  }
  cout << endl;

  const SATLiteralStack& minimized = s.explicitlyMinimizedFailedAssumptions();
  for (unsigned i = 0; i < minimized.size(); i++) {
    SATLiteral lit = minimized[i];
    if (lit.polarity()) {
      cout << (char)('A' + lit.var()-1);
    } else {
      cout << (char)('a' + lit.var()-1);
    }
  }
  cout << endl;
}

TEST_FUN(testSolvingUnderAssumptions)
{
  cout << endl << "CaDiCaL" << endl;
  CadicalInterfacing sMini(*env.options,true);
  testAssumptions(sMini);

  /*cout << endl << "Z3" << endl;
  {
    SAT2FO sat2fo;
    Z3Interfacing sZ3(*env.options,sat2fo,true);
    testAssumptions(sZ3);
  }*/
}
