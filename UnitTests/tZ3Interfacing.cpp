/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#include "Test/UnitTesting.hpp"
#include "Test/SyntaxSugar.hpp"
#include "Test/SyntaxSugar.hpp"
#include "SAT/Z3Interfacing.hpp"
#include <fstream>

#if VZ3

#define DBG_ON 1
// #if DBG_ON
// #define DEBUG(...)
// #else
// #define DEBUG(...) DBG(__VA_ARGS__)
// #endif

using namespace std;
using namespace Shell;
using namespace SAT;

using Sort = unsigned;
using FuncId = unsigned;

/////////////////////////////////////////////////////////////////////////////////////////
// 1) TEST SOLVING
/////////////////////////////////////////////////////////////////////////////////////////

/** runs z3 on a bunch of vampire literals as assumptions, and checks the status afterwards */
void checkStatus(SAT::Z3Interfacing& z3, SAT2FO& s2f, SATSolver::Status expected, Stack<Literal*> assumptions) 
{
  for (auto a : assumptions) {
    // Stack<SATLiteral> clause{s2f.toSAT(a)};
    // z3.addClause(SATClause::fromStack(clause));
    z3.addAssumption(s2f.toSAT(a));
  }

  auto status = z3.solve();
  if(status != expected) {
    cout << "[ input    ] " << endl;
    for (auto a : assumptions) {
      cout << "\t" << *a << endl;
    }
    cout << "[ expected ] " <<  expected << endl;
    cout << "[ is       ] " <<  status << endl;
    if (status == SATSolver::Status::SATISFIABLE) {
      cout << "[ model    ] " <<  z3.getModel() << endl;

    }
    exit(-1);
  }
  z3.retractAllAssumptions();
}

void checkStatus(SATSolver::Status expected, Stack<Literal*> assumptions)
{
  SAT2FO s2f;
  SAT::Z3Interfacing z3(s2f, /* show z3 */ DBG_ON == 1, /* unsat core */ false, /* export smtlib */ "", Shell::Options::ProblemExportSyntax::SMTLIB);
  checkStatus(z3, s2f, expected, assumptions);
}


////////////////////////////////////
// Real Tests
/////////////////

TEST_FUN(solve__real__simple_01) {
  NUMBER_SUGAR(Real)
  checkStatus(
      SATSolver::Status::UNSATISFIABLE,
      { num(3) == num(0) });
}

TEST_FUN(solve__real__simple_02) {
  NUMBER_SUGAR(Real)
  DECL_CONST(a, Real)
  checkStatus(
      SATSolver::Status::SATISFIABLE,
      { num(3) == a });
}


TEST_FUN(solve__rat__simple_03) {
  NUMBER_SUGAR(Real)
  checkStatus(
      SATSolver::Status::UNSATISFIABLE,
      { num(3) == num(3) + 2 * num(7) });
}

TEST_FUN(solve__rat__simple_04) {
  NUMBER_SUGAR(Real)
  checkStatus(
      SATSolver::Status::UNSATISFIABLE,
      { num(17) != num(3) + 2 * num(7) });
}

////////////////////////////////////
// FOOL Tests
/////////////////

TEST_FUN(solve__fool__simple_01) {
  DECL_VAR(x, 0);
  checkStatus(
      SATSolver::Status::SATISFIABLE,
      { fool(false) == x });
}

TEST_FUN(solve__fool__simple_02) {
  DECL_VAR(x, 0);
  checkStatus(
      SATSolver::Status::SATISFIABLE,
      { fool(true) == x });
}


TEST_FUN(solve__fool__simple_03) {
  checkStatus(
      SATSolver::Status::UNSATISFIABLE,
      { fool(true) == fool(false) });
}

////////////////////////////////////
// term algebra tests
/////////////////

#define DECL_LIST(sort)                                                                                       \
  DECL_SORT(list)                                                                                             \
                                                                                                              \
  DECL_CONST(nil, list)                                                                                       \
  DECL_FUNC(cons, { sort, list  }, list)                                                                      \
                                                                                                              \
  DECL_TERM_ALGEBRA(list, {nil, cons})                                                                        \
  __ALLOW_UNUSED(                                                                                             \
    auto head = cons.dtor(0);                                                                                 \
    auto tail = cons.dtor(1);                                                                                 \
  )                                                                                                           \



TEST_FUN(solve__dty__01) {
  DECL_SORT(alpha)
  DECL_LIST(alpha)

  DECL_CONST(a0, alpha)
  DECL_CONST(a1, alpha)

  checkStatus(SATSolver::Status::UNSATISFIABLE, { cons(a0, nil) == nil });
  checkStatus(SATSolver::Status::UNSATISFIABLE, { cons(a0, nil) == cons(a1, nil), a0 != a1 });
}


// data Even = Zero | SuccEven Odd
// data Odd = SuccOdd Even
#define DECL_EVEN_ODD                                                                                         \
  DECL_SORT(even)                                                                                             \
  DECL_SORT(odd)                                                                                              \
                                                                                                              \
  DECL_CONST(zero, even)                                                                                      \
  DECL_FUNC(succEven, { odd  }, even)                                                                         \
                                                                                                              \
  DECL_FUNC(succOdd, { even }, odd)                                                                           \
                                                                                                              \
  DECL_TERM_ALGEBRA(even, {zero, succEven})                                                                   \
  DECL_TERM_ALGEBRA(odd , {      succOdd })                                                                   \


TEST_FUN(solve__dty__02) {
  DECL_EVEN_ODD
  checkStatus(SATSolver::Status::UNSATISFIABLE, { succEven(succOdd(zero)) == zero });
}


TEST_FUN(solve__dty__03_01) {
  // we have a non-mutually recursive datatype that depends on a non-mutual but recursive datatype that 
  DECL_EVEN_ODD
  DECL_LIST(even)
  // request non-mutual first
  checkStatus(SATSolver::Status::UNSATISFIABLE, { cons(succEven(succOdd(zero)), nil) == cons(zero, nil) });
}

TEST_FUN(solve__dty__03_02) {
  // we have a non-mutually recursive datatype that depends on a non-mutual but recursive datatype that 
  DECL_EVEN_ODD
  DECL_LIST(even)
  // request mutual only
  checkStatus(SATSolver::Status::UNSATISFIABLE, { succEven(succOdd(zero)) == zero });
}

TEST_FUN(solve__dty__03_03) {
  // we have a non-mutually recursive datatype that depends on a non-mutual but recursive datatype that 
  DECL_EVEN_ODD
  DECL_LIST(even)
  // request mutual first
  checkStatus(SATSolver::Status::UNSATISFIABLE, { succEven(succOdd(zero)) == zero });
  checkStatus(SATSolver::Status::UNSATISFIABLE, { cons(succEven(succOdd(zero)), nil) == cons(zero, nil) });
}

/////////////////////////////////////////////////////////////////////////////////////////
// 2) TEST INSTANTIATION
/////////////////////////////////////////////////////////////////////////////////////////

void checkInstantiation(SAT::Z3Interfacing& z3, SAT2FO& s2f, Stack<Literal*> assumptions, TermList toInstantiate, TermList expected)
{
  for (auto a : assumptions) {
    z3.addAssumption(s2f.toSAT(a));
  }

  auto status = z3.solve();
  ASS_EQ(status, SATSolver::Status::SATISFIABLE);
  auto result = z3.evaluateInModel(toInstantiate.term());
  if (result != expected.term()) {
    cout << "[ input    ] " << endl;
    for (auto a : assumptions) {
      cout << "\t" << *a << endl;
    }
    cout << "[ toInstantiate ] " <<  toInstantiate << endl;
    cout << "[ expected      ] " <<  expected << endl;
    cout << "[ is            ] " <<  ( result == nullptr ? "null" : result->toString() ) << endl;
    cout << "[ model         ] " <<  z3.getModel() << endl;
    exit(-1);
  }
  z3.retractAllAssumptions();
}


/**
 * Runs z3 on a bunch of vampire literals as assumptions, that need to be satisfyable. 
 * Then  the term toInstantiate will be instantiated with the model. The instantiated 
 * term will be checked to be equal to the term expected.
 */
void checkInstantiation(Stack<Literal*> assumptions, TermList toInstantiate, TermList expected)
{
  SAT2FO s2f;
  SAT::Z3Interfacing z3(s2f, /* show z3 */ DBG_ON == 1, /* unsat core */ false, /* export smtlib */ "", Shell::Options::ProblemExportSyntax::SMTLIB);
  return checkInstantiation(z3, s2f, assumptions, toInstantiate, expected);
}


////////////////////////////////////
// Real Tests
/////////////////

TEST_FUN(instantiate__rat__simple_01) {
  NUMBER_SUGAR(Real)
  DECL_CONST(c, Real)
  checkInstantiation(
      { c * 3 == 9 },
      c, num(3)
      );
}

TEST_FUN(instantiate__rat__simple_02) {
  NUMBER_SUGAR(Real)
  DECL_CONST(c, Real)
  checkInstantiation(
      { c * c == 9, c < 0 },
      c, num(-3)
      );
}

////////////////////////////////////
// term algebra tests
/////////////////


TEST_FUN(instantiate__list_01) {
  NUMBER_SUGAR(Real)
  DECL_LIST(Real)

  DECL_CONST(c, Real)

  checkInstantiation(
      { cons(c, nil) == cons(num(3), nil) },
      c, num(3)
      );
}

TEST_FUN(instantiate__list_02) {
  NUMBER_SUGAR(Real)
  DECL_LIST(Real)

  DECL_CONST(l, list)
  DECL_CONST(h, Real)
  DECL_CONST(t, list)

  checkInstantiation(
      { tail(l) == cons(num(2), nil)
      , head(l) == 1
      , cons(h,t) == l  // <- necessary since selectors are partial
      },
      l, cons(1,cons(2,nil))
      );
}

TEST_FUN(segfault01) {
  Z3_config config = Z3_mk_config();
  Z3_context context = Z3_mk_context(config);

  Z3_symbol consName = Z3_mk_string_symbol(context, "e00");
  Z3_symbol sortNames = Z3_mk_string_symbol(context, "Enum");

  Z3_func_decl enumCtor;
  Z3_func_decl enumDiscr;
  Z3_sort sorts = Z3_mk_enumeration_sort(context,
      sortNames, 1,
      &consName, &enumCtor, &enumDiscr);

  Z3_symbol c1_sym = Z3_mk_string_symbol(context, "c1");
  Z3_symbol c2_sym = Z3_mk_string_symbol(context, "c1");

  Z3_ast c1 = Z3_mk_const(context, c1_sym, sorts);
  Z3_ast c2 = Z3_mk_const(context, c2_sym, sorts);

  Z3_solver solver = Z3_mk_solver(context);

  Z3_solver_push(context, solver);

    Z3_ast expr = Z3_mk_eq(context, c1, c2);
    Z3_solver_assert(context, solver, expr);
    Z3_solver_check(context, solver);
    Z3_solver_check(context, solver);

  Z3_solver_pop(context, solver, 1);

  // std::cout << "segfault occurs between here ..." << std::endl;
  Z3_solver_check(context, solver);
  // std::cout << "... and here" << std::endl;
}

TEST_FUN(segfault02) {
  DECL_SORT(Enum)

  DECL_CONST(e00, Enum)

  DECL_TERM_ALGEBRA(Enum, { e00 })

  DECL_CONST(inst159, Enum)
  DECL_CONST(inst160, Enum)


  SAT2FO s2f;
  SAT::Z3Interfacing z3(s2f, /* show z3 */ DBG_ON == 1, /* unsat core */ false, /* export smtlib */ "", Shell::Options::ProblemExportSyntax::SMTLIB);

  checkStatus(z3, s2f, SATSolver::Status::SATISFIABLE, { inst159 == inst160 });
  z3.solve();
}

#endif // VZ3
