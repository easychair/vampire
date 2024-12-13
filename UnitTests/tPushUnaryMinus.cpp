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
#include "Indexing/TermSharing.hpp"
#include "Inferences/PushUnaryMinus.hpp"
#include "Inferences/InterpretedEvaluation.hpp"
#include "Kernel/Ordering.hpp"
#include "Inferences/PolynomialEvaluation.hpp"
#include "Inferences/Cancellation.hpp"

#include "Test/SyntaxSugar.hpp"
#include "Test/TestUtils.hpp"
#include "Lib/Coproduct.hpp"
#include "Test/SimplificationTester.hpp"
#include "Test/GenerationTester.hpp"
#include "Kernel/KBO.hpp"

using namespace std;
using namespace Kernel;
using namespace Inferences;
using namespace Test;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// TEST UNIT INITIALIZATION
/////////////////////////////////////

/** 
 * NECESSARY: We need a subclass of SimplificationTester
 */
class PumSimplTester : public Test::Simplification::SimplificationTester
{
public:

  virtual Kernel::Clause* simplify(Kernel::Clause* in) const override 
  {
    PushUnaryMinus pum;
    return pum.simplify(in);
  }

  virtual bool eq(Kernel::Clause* lhs, Kernel::Clause* rhs) const  override
  { return TestUtils::eqModAC(lhs, rhs); }
};

REGISTER_SIMPL_TESTER(PumSimplTester)

#define MY_SYNTAX_SUGAR                                                                                       \
  NUMBER_SUGAR(Real)                                                                                          \
  DECL_DEFAULT_VARS                                                                                           \
  DECL_CONST(a, Real)                                                                                         \
  DECL_CONST(b, Real)                                                                                         \
  DECL_FUNC(f, {Real}, Real)                                                                                  \
  DECL_PRED(p, {Real})                                                                                        \
  DECL_PRED(q, {Real})                                                                                        \

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// TEST CASES
/////////////////////////////////////

TEST_SIMPLIFY(test_1,
    Simplification::Success()
      .input(    clause({  p(-(a + b))  }))
      .expected( clause({  p(-a + -b)  }))
    )

TEST_SIMPLIFY(test_2,
    Simplification::NotApplicable()
      .input( clause({ p(-a + -b) }))
    )

TEST_SIMPLIFY(test_3,
    Simplification::NotApplicable()
      .input( clause({ p(a) }))
    )

TEST_SIMPLIFY(test_4,
    Simplification::Success()
      .input(    clause({  p(-(a + b)), p(a)  }))
      .expected( clause({  p(-a + -b), p(a)  }))
    )

TEST_SIMPLIFY(test_5,
    Simplification::Success()
      .input(    clause({  p(a), p(-(a + b))  }))
      .expected( clause({  p(a), p(-a + -b)  }))
    )

// TEST_SIMPLIFY(test_3,
//     Simplification::Success()
//       .input(    clause({  3 * x != 6, x < x  }))
//       .expected( clause({  /* 2 < 2 */  }))
//     )
//
//   // 2x + y = x + y ==> 0 = 2x + y - x - y ==> 0 = x
// TEST_SIMPLIFY(test_4,
//     Simplification::Success()
//       .input(    clause({  2 * x + y != x + y, p(x) }))
//       .expected( clause({  p(0)  }))
//     )
//
// TEST_SIMPLIFY(test_uninterpreted,
//     Simplification::Success()
//       .input(    clause({  3 * f(x) != y, x < y  }))
//       .expected( clause({  x < 3 * f(x)  }))
//     )
//
//   // x!=4 \/ x+y != 5 \/ C[x]
//   //         4+y != 5 \/ C[4]
//   //                     C[4]
// TEST_SIMPLIFY(test_multiplesteps_1,
//     Simplification::Success()
//       .input(    clause({  x != 4, x + y != 5, x < f(x)  }))
//       .expected( clause({  4 < f(4)  }))
//     )
//
//   // x!=4 \/ x+y != 5 \/ C[x,y]
//   //         4+y != 5 \/ C[4,y]
//   //                     C[4,1]
// TEST_SIMPLIFY(test_multiplesteps_2,
//     Simplification::Success()
//       .input(    clause({  x != 4, x + y !=  5, x < f(y)  }))
//       .expected( clause({  4 < f(1)  }))
//     )
//
// TEST_SIMPLIFY(test_div,
//     Simplification::Success()
//       .input(    clause({  x / 3 != 4, p(x)  }))
//       .expected( clause({  p(12)  }))
//     )
