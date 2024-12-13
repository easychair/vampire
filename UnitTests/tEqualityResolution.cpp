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
#include "Test/TestUtils.hpp"
#include "Test/GenerationTester.hpp"

#include "Inferences/EqualityResolution.hpp"

using namespace Test;

REGISTER_GEN_TESTER(Test::Generation::GenerationTester<Inferences::EqualityResolution>(EqualityResolution()))

/**
 * NECESSARY: We neet to tell the tester which syntax sugar to import for creating terms & clauses. 
 * See Test/SyntaxSugar.hpp for which kinds of syntax sugar are available
 */
#define MY_SYNTAX_SUGAR                                                                                       \
  DECL_DEFAULT_VARS                                                                                           \
  DECL_SORT(s)                                                                                                \
  DECL_FUNC(f, {s}, s)                                                                                        \
  DECL_FUNC(g, {s}, s)                                                                                        \
  DECL_CONST(a, s)                                                                                            \
  DECL_PRED (p, {s})                                                                                          \
  DECL_PRED (q, {s})                                                                                          \


/** Defines a test case. */
TEST_GENERATION(test_01,                                   // <- name
    Generation::AsymmetricTest()
      .input(     clause({  selected(x != f(a)), p(x)  })) // <- input clause
      .expected(exactly(                                   // <- a list of exactly which clauses are expected
            clause({  p(f(a))  })                          //    to be returned. Order matters!
      ))
      .premiseRedundant(false)                             // <- shall the premis be removed from the search 
                                                           //    space after the rule application ? 
                                                           //    (default value: false)
    )

TEST_GENERATION(test_02,
    Generation::AsymmetricTest()
      .input(     clause({  x != f(a), selected(p(x))  }))
      .expected( exactly())
    )

TEST_GENERATION(test_03,
    Generation::AsymmetricTest()
      .input(     clause({  selected(x != f(a)), selected(p(x))  }))
      .expected( exactly( clause({  p(f(a))                              })))
    )

TEST_GENERATION(test_04,
    Generation::AsymmetricTest()
      .input(     clause({  selected(g(x) != f(a)), p(x)  }))
      .expected( exactly())
    )

TEST_GENERATION(test_05,
    Generation::AsymmetricTest()
      .input(     clause({  selected(f(g(x)) != f(y))  }))
      .expected( exactly( clause({})))
    )

TEST_GENERATION(test_06,
    Generation::AsymmetricTest()
      .input(     clause({  selected(f(g(x)) != f(x))  }))
      .expected( exactly())
    )
