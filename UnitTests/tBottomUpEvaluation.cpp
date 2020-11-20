#include "Test/UnitTesting.hpp"
#include "Test/SyntaxSugar.hpp"
#include "Kernel/BottomUpEvaluation.hpp"
#include "Kernel/BottomUpEvaluation/TermList.hpp"
#include "Kernel/Term.hpp"

using namespace Kernel;
using namespace Inferences;
using namespace Test;


#define UNIT_ID BottomUpEvaluation
UT_CREATE;

TEST_FUN(example_01__replace_all_vars_by_term) {
  /* syntax sugar imports */
  DECL_DEFAULT_VARS
  DECL_SORT(s)
  DECL_CONST(a, s)
  DECL_FUNC(f, {s}, s)
  DECL_FUNC(g, {s,s}, s)

  /* defines how to evaluate bottom up. 
   * all variables are being replaced by a constant term in this case */
  struct Eval {
    TermList replacement;

    using Arg    = TermList;
    using Result = TermList;

    TermList operator()(TermList toEval, TermList* evaluatedChildren) {
      if (toEval.isVar()) {
        return replacement;
      } else {
        return TermList(Term::create(toEval.term(), evaluatedChildren));
      }
    }
  };


  /* test specification */
  TermList input    = g(f(x), y);
  TermList expected = g(f(a), a);

  /* actual evaluation */
  TermList result =  evaluateBottomUp(input, Eval{a});

  ASS_EQ(result, expected)
}

TEST_FUN(example_02__compute_size) {
  /* syntax sugar imports */
  DECL_DEFAULT_VARS
  DECL_SORT(s)
  DECL_CONST(a, s)
  DECL_FUNC(f, {s}, s)
  DECL_FUNC(g, {s,s}, s)

  /* defines how to evaluate bottom up. 
   * computes the size of the term (number of function & variable symbols) */
  struct Eval {
    TermList replacement;

    using Arg    = TermList;
    using Result = unsigned;

    unsigned operator()(TermList toEval, unsigned* evaluatedChildren) {
      if (toEval.isVar()) {
        return 1;
      } else {
        unsigned out = 1;
        for (int i = 0; i < toEval.term()->arity(); i++) {
          out += evaluatedChildren[i];
        }
        return out;
      }
    }
  };


  /* test specification */
  TermList input    = g(f(x), f(f(x)));
  //                    ^^^^    ^^^^ size of this sub-term will only be evaluated once due to memo

  /* actual evaluation */
  Memo::Hashed<TermList, unsigned> memo{};
  auto size =  evaluateBottomUp(input, Eval{a}, memo);

  ASS_EQ(size, 6)
}
