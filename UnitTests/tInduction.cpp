/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */


#include <unordered_set>

#include "Test/UnitTesting.hpp"
#include "Test/SyntaxSugar.hpp"
#include "Test/TestUtils.hpp"
#include "Test/GenerationTester.hpp"

#include "Indexing/TermIndex.hpp"
#include "Indexing/TermSubstitutionTree.hpp"
#include "Indexing/LiteralIndex.hpp"
#include "Indexing/LiteralSubstitutionTree.hpp"
#include "Indexing/TermIndex.hpp"
#include "Indexing/TermSubstitutionTree.hpp"

#include "Kernel/FormulaUnit.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/FormulaVarIterator.hpp"

#include "Inferences/Induction.hpp"

using namespace std;
using namespace Test;
using namespace Test::Generation;

#define SKOLEM_VAR_MIN 100
#define DECL_SKOLEM_VAR(x, i) DECL_VAR(x, i+SKOLEM_VAR_MIN)

namespace InductionTest {

LiteralIndex<LiteralClause>* comparisonIndex() {
  return new UnitIntegerComparisonLiteralIndex(new LiteralSubstitutionTree<LiteralClause>());
}

TermIndex<TermLiteralClause>* intInductionIndex() {
  return new InductionTermIndex(new TermSubstitutionTree<TermLiteralClause>());
}

TermIndex<TermLiteralClause>* structInductionIndex() {
  return new StructInductionTermIndex(new TermSubstitutionTree<TermLiteralClause>());
}

auto getIndices() {
  return Stack<std::function<Index*()>> { 
    std::function<Index*()>([]() { return comparisonIndex(); }), 
    std::function<Index*()>([]() { return intInductionIndex(); }), 
    std::function<Index*()>([]() { return structInductionIndex(); }), 
  };
}

inline Clause* fromInduction(Clause* cl) {
  cl->inference().setInductionDepth(1);
  return cl;
}

InductionContext inductionContext(TermSugar t, std::initializer_list<Clause*> cls) {
  InductionContext res({ t.sugaredExpr().term() });
  for (const auto& cl : cls) {
    for (unsigned i = 0; i < cl->length(); i++) {
      res.insert(cl, (*cl)[i]);
    }
  }
  return res;
}

namespace Inferences {
std::ostream& operator<<(std::ostream& out, const InductionContext& context) {
  out << context.toString();
  return out;
}
}

void assertContextReplacement(ContextReplacement& cr, Stack<InductionContext> contexts) {
  Stack<InductionContext> res;
  while (cr.hasNext()) {
    res.push(cr.next());
  }
  ASS_EQ(res.size(), contexts.size());
  ASS(TestUtils::permEq(res, contexts, [](const InductionContext& lhs, const InductionContext& rhs) {
    return InductionFormulaIndex::represent(lhs) == InductionFormulaIndex::represent(rhs);
  }));
}

class GenerationTesterInduction
  : public GenerationTester<Induction>
{
public:
  GenerationTesterInduction()
    : GenerationTester<Induction>(Induction()), _subst()
  {}

  ~GenerationTesterInduction() {
    _btd.drop();
  }

  /**
   * Generated induction clauses are special in that they contain fresh
   * Skolem constants. In order to check these, we use variables instead
   * of the constants we cannot predefine, and require that these variables
   * are mapped bijectively to the new Skolem constants, hence this override.
   */
  bool eq(Kernel::Clause* lhs, Kernel::Clause* rhs) override
  {
    // there can be false positive matches which later (in a different literal
    // or clause) can turn out to be the wrong ones and we have to backtrack
    // TODO: are all of these backtracking calls necessary?
    _subst.bdRecord(_btd);
    if (!TestUtils::permEq(*lhs, *rhs, [this](Literal* l, Literal* r) -> bool {
      if (l->polarity() != r->polarity()) {
        _btd.backtrack();
        return false;
      }
      VList::Iterator vit(freeVariables(r));
      while (vit.hasNext()) {
        auto v = vit.next();
        if (!_varsMatched.count(v)) {
          _btd.addBacktrackObject(new MatchedVarBacktrackObject(_varsMatched, v));
          _varsMatched.insert(v);
        }
      }
      _subst.bdRecord(_btd);
      if (_subst.match(Kernel::TermList(r), 0, Kernel::TermList(l), 1)) {
        if (matchAftercheck()) {
          _subst.bdDone();
          return true;
        }
      }

      _subst.bdDone();
      _btd.backtrack();
      _subst.bdRecord(_btd);
      if (l->isEquality() && r->isEquality() &&
        _subst.match(*r->nthArgument(0), 0, *l->nthArgument(1), 1) &&
        _subst.match(*r->nthArgument(1), 0, *l->nthArgument(0), 1))
      {
        if (matchAftercheck()) {
          _subst.bdDone();
          return true;
        }
      }
      _subst.bdDone();
      _btd.backtrack();
      return false;
    })) {
      _subst.bdDone();
      _btd.backtrack();
      return false;
    }
    _subst.bdDone();
    return true;
  }

private:
  bool matchAftercheck() {
    DHMap<TermList, unsigned> inverse;
    for (const auto& i : _varsMatched) {
      auto t = _subst.apply(TermList(i,false),0);
      unsigned v;
      // we check that each variable encountered so far from
      // the expected set is bijectively mapped to something
      if (inverse.find(t,v)) {
        return false;
      } else {
        inverse.insert(t,i);
      }
      // "Skolem" variables should bind to Skolem constants
      if (i >= SKOLEM_VAR_MIN) {
        if (t.isVar() || !env.signature->getFunction(t.term()->functor())->skolem()) {
          return false;
        }
      // normal variables should bind to variables
      } else {
        if (t.isTerm()) {
          return false;
        }
      }
    }
    return true;
  }

  Kernel::RobSubstitution _subst;
  std::unordered_set<unsigned> _varsMatched;
  BacktrackData _btd;

  class MatchedVarBacktrackObject : public BacktrackObject {
  public:
    MatchedVarBacktrackObject(unordered_set<unsigned>& s, unsigned i) : _s(s), _i(i) {}
    void backtrack() override {
      _s.erase(_i);
    }
  private:
    unordered_set<unsigned>& _s;
    unsigned _i;
  };
};

#define TEST_GENERATION_INDUCTION(name, expr)                                                                 \
  TEST_FUN(name) {                                                                                            \
    GenerationTesterInduction tester;                                                                         \
    __ALLOW_UNUSED(MY_SYNTAX_SUGAR)                                                                           \
    auto test = expr;                                                                                         \
    test.run(tester);                                                                                         \
  }                                                                                                           \

/**
 * NECESSARY: We need to tell the tester which syntax sugar to import for creating terms & clauses.
 * See Test/SyntaxSugar.hpp for which kinds of syntax sugar are available
 */
#define MY_SYNTAX_SUGAR                                                                    \
  DECL_DEFAULT_VARS                                                                        \
  DECL_SKOLEM_VAR(skx0,0)                                                                  \
  DECL_SKOLEM_VAR(skx1,1)                                                                  \
  DECL_SKOLEM_VAR(skx2,2)                                                                  \
  DECL_SKOLEM_VAR(skx3,3)                                                                  \
  DECL_SKOLEM_VAR(skx4,4)                                                                  \
  DECL_SKOLEM_VAR(skx5,5)                                                                  \
  DECL_SKOLEM_VAR(skx6,6)                                                                  \
  DECL_SKOLEM_VAR(skx7,7)                                                                  \
  DECL_SKOLEM_VAR(skx8,8)                                                                  \
  DECL_SKOLEM_VAR(skx9,9)                                                                  \
  DECL_SKOLEM_VAR(skx10,10)                                                                \
  DECL_SKOLEM_VAR(skx11,11)                                                                \
  DECL_SKOLEM_VAR(skx12,12)                                                                \
  DECL_SKOLEM_VAR(skx13,13)                                                                \
  DECL_SKOLEM_VAR(skx14,14)                                                                \
  DECL_VAR(x3,3)                                                                           \
  DECL_VAR(x4,4)                                                                           \
  DECL_VAR(x5,5)                                                                           \
  DECL_SORT(s)                                                                             \
  DECL_SORT(u)                                                                             \
  DECL_SKOLEM_CONST(sK1, s)                                                                \
  DECL_FUN_DEF(def_s,sK1())                                                                \
  DECL_SKOLEM_CONST(sK2, s)                                                                \
  DECL_SKOLEM_CONST(sK3, s)                                                                \
  DECL_SKOLEM_CONST(sK4, s)                                                                \
  DECL_SKOLEM_CONST(sK5, u)                                                                \
  DECL_CONST(b, s)                                                                         \
  DECL_FUNC(r, {s}, s)                                                                     \
  DECL_TERM_ALGEBRA(s, {b, r})                                                             \
  __ALLOW_UNUSED(                                                                          \
    auto r0 = r.dtor(0);                                                                   \
    TermSugar ph_s(TermList(getPlaceholderForTerm({ sK1.sugaredExpr().term() }, 0)));      \
  )                                                                                        \
  DECL_CONST(b1, u)                                                                        \
  DECL_CONST(b2, u)                                                                        \
  DECL_FUNC(r1, {s, u, u}, u)                                                              \
  DECL_FUNC(r2, {u, s}, u)                                                                 \
  DECL_TERM_ALGEBRA(u, {b1, b2, r1, r2})                                                   \
  DECL_FUNC(f, {s, s}, s)                                                                  \
  DECL_FUNC(g, {s}, s)                                                                     \
  DECL_FUNC(h, {s, s}, s)                                                                  \
  DECL_PRED(p, {s})                                                                        \
  DECL_PRED_DEF(def_p, p(x))                                                               \
  DECL_PRED(p1, {s})                                                                       \
  DECL_PRED(q, {u})                                                                        \
  NUMBER_SUGAR(Int)                                                                        \
  DECL_PRED(pi, {Int})                                                                     \
  DECL_FUNC(fi, {Int, s}, Int)                                                             \
  DECL_FUNC(gi, {Int}, Int)                                                                \
  DECL_CONST(sK6, Int)                                                                     \
  DECL_CONST(sK7, Int)                                                                     \
  DECL_CONST(sK8, Int)                                                                     \
  DECL_CONST(bi, Int)

TEST_FUN(test_tester1) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  GenerationTesterInduction tester;
  // first literal is matched both ways but none of them works
  ASS(!tester.eq(
    clause({ r(sK1) == r(x), f(r(sK1),y) != z }),
    clause({ r(skx1) == r(x3), f(r(x3),x4) != x5 })));
}

TEST_FUN(test_tester2) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  GenerationTesterInduction tester;
  // second clause cannot be matched because of x4
  ASS(!tester.eq(
    clause({ r(sK1) == r(x), f(r(sK1),y) != z }),
    clause({ r(skx1) == r(x3), f(r(skx1),x4) != x4 })));
}

TEST_FUN(test_tester3) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  GenerationTesterInduction tester;
  // y is matched to both y4 and y5
  ASS(!tester.eq(
    clause({ r(sK1) == r(x), f(r(sK1),y) != y }),
    clause({ r(skx1) == r(x3), f(r(skx1),x4) != x5 })));
}

TEST_FUN(test_tester4) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  GenerationTesterInduction tester;
  // normal match
  ASS(tester.eq(
    clause({ r(sK1) == r(x), f(r(sK1),y) != z }),
    clause({ r(skx1) == r(x3), f(r(skx1),x4) != x5 })));
}

// positive literals are not considered 1
TEST_GENERATION_INDUCTION(test_01,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  p(f(sK1,sK2)) }))
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// positive literals are not considered 2
TEST_GENERATION_INDUCTION(test_02,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  f(sK1,sK2) == g(sK1) }))
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// non-ground literals are not considered
TEST_GENERATION_INDUCTION(test_03,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  f(sK1,skx0) != g(sK1) }))
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// normal case sik=one
TEST_GENERATION_INDUCTION(test_04,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  ~p(f(sK1,sK2)) }))
      .expected({
        clause({ ~p(f(b,sK2)), p(f(skx0,sK2)) }),
        clause({ ~p(f(b,sK2)), ~p(f(r(skx0),sK2)) }),
        clause({ ~p(f(sK1,b)), p(f(sK1,skx1)) }),
        clause({ ~p(f(sK1,b)), ~p(f(sK1,r(skx1))) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 2) })
    )

// normal case sik=two
TEST_GENERATION_INDUCTION(test_05,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" }, { "structural_induction_kind", "two" } })
      .indices(getIndices())
      .input( clause({  ~p(f(sK1,sK2)) }))
      .expected({
        clause({ skx0 != r(r0(skx0)), p(f(r0(skx0),sK2)) }),
        clause({ ~p(f(skx0,sK2)) }),
        clause({ skx1 != r(r0(skx1)), p(f(sK1,r0(skx1))) }),
        clause({ ~p(f(sK1,skx1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 2) })
    )

// TODO this case is a bit hard to test since new predicates are introduced,
// so we need to customize the test suite for this even more, checking certain
// properties of these new predicates and matching it with some literals.
// This induction mode is not that useful compared to other sik modes to make
// the effort worthwhile.
// // normal case sik=three
// TEST_GENERATION_INDUCTION(test_06,
//     Generation::AsymmetricTest()
//       .options({ { "induction", "struct" }, { "structural_induction_kind", "three" } })
//       .indices(getIndices())
//       .input( clause({  f(sK1,sK2) != g(sK1) }))
//       .expected({ })
//     )

// generalizations
TEST_GENERATION_INDUCTION(test_07,
    Generation::AsymmetricTest()
      .options({ { "induction_gen", "on" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,sK3))) }) )
      .expected({
        // sK1 100
        clause({ f(f(g(b),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(skx0),f(sK2,sK4)),sK1) == g(f(sK1,f(sK2,sK3))) }),
        clause({ f(f(g(b),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(r(skx0)),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,sK3))) }),

        // sK1 010
        clause({ f(f(g(sK1),f(sK2,sK4)),b) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),skx1) == g(f(sK1,f(sK2,sK3))) }),
        clause({ f(f(g(sK1),f(sK2,sK4)),b) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),r(skx1)) != g(f(sK1,f(sK2,sK3))) }),

        // sK1 001
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(b,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),sK1) == g(f(skx2,f(sK2,sK3))) }),
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(b,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(r(skx2),f(sK2,sK3))) }),

        // sK1 110
        clause({ f(f(g(b),f(sK2,sK4)),b) != g(f(sK1,f(sK2,sK3))), f(f(g(skx3),f(sK2,sK4)),skx3) == g(f(sK1,f(sK2,sK3))) }),
        clause({ f(f(g(b),f(sK2,sK4)),b) != g(f(sK1,f(sK2,sK3))), f(f(g(r(skx3)),f(sK2,sK4)),r(skx3)) != g(f(sK1,f(sK2,sK3))) }),

        // sK1 101
        clause({ f(f(g(b),f(sK2,sK4)),sK1) != g(f(b,f(sK2,sK3))), f(f(g(skx4),f(sK2,sK4)),sK1) == g(f(skx4,f(sK2,sK3))) }),
        clause({ f(f(g(b),f(sK2,sK4)),sK1) != g(f(b,f(sK2,sK3))), f(f(g(r(skx4)),f(sK2,sK4)),sK1) != g(f(r(skx4),f(sK2,sK3))) }),

        // sK1 011
        clause({ f(f(g(sK1),f(sK2,sK4)),b) != g(f(b,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),skx5) == g(f(skx5,f(sK2,sK3))) }),
        clause({ f(f(g(sK1),f(sK2,sK4)),b) != g(f(b,f(sK2,sK3))), f(f(g(sK1),f(sK2,sK4)),r(skx5)) != g(f(r(skx5),f(sK2,sK3))) }),

        // sK1 111
        clause({ f(f(g(b),f(sK2,sK4)),b) != g(f(b,f(sK2,sK3))), f(f(g(skx6),f(sK2,sK4)),skx6) == g(f(skx6,f(sK2,sK3))) }),
        clause({ f(f(g(b),f(sK2,sK4)),b) != g(f(b,f(sK2,sK3))), f(f(g(r(skx6)),f(sK2,sK4)),r(skx6)) != g(f(r(skx6),f(sK2,sK3))) }),

        // sK2 10
        clause({ f(f(g(sK1),f(b,sK4)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(skx7,sK4)),sK1) == g(f(sK1,f(sK2,sK3))) }),
        clause({ f(f(g(sK1),f(b,sK4)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(r(skx7),sK4)),sK1) != g(f(sK1,f(sK2,sK3))) }),

        // sK2 01
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(b,sK3))), f(f(g(sK1),f(sK2,sK4)),sK1) == g(f(sK1,f(skx8,sK3))) }),
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(b,sK3))), f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(r(skx8),sK3))) }),

        // sK2 11
        clause({ f(f(g(sK1),f(b,sK4)),sK1) != g(f(sK1,f(b,sK3))), f(f(g(sK1),f(skx9,sK4)),sK1) == g(f(sK1,f(skx9,sK3))) }),
        clause({ f(f(g(sK1),f(b,sK4)),sK1) != g(f(sK1,f(b,sK3))), f(f(g(sK1),f(r(skx9),sK4)),sK1) != g(f(sK1,f(r(skx9),sK3))) }),

        // sK3 1
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,b))), f(f(g(sK1),f(sK2,sK4)),sK1) == g(f(sK1,f(sK2,skx10))) }),
        clause({ f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,b))), f(f(g(sK1),f(sK2,sK4)),sK1) != g(f(sK1,f(sK2,r(skx10)))) }),

        // sK4 1
        clause({ f(f(g(sK1),f(sK2,b)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(sK2,skx11)),sK1) == g(f(sK1,f(sK2,sK3))) }),
        clause({ f(f(g(sK1),f(sK2,b)),sK1) != g(f(sK1,f(sK2,sK3))), f(f(g(sK1),f(sK2,r(skx11))),sK1) != g(f(sK1,f(sK2,sK3))) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->generalizedInductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 4),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 12),
                        TEST_FN_ASS_EQ(env.statistics->generalizedInductionApplication, 8) })
    )

// complex terms
TEST_GENERATION_INDUCTION(test_08,
    Generation::AsymmetricTest()
      .options({ { "induction_on_complex_terms", "on" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != g(f(sK1,f(sK2,g(sK1)))) }) )
      .expected({
        // sK1
        clause({ f(f(g(b),f(sK2,sK3)),b) != g(f(b,f(sK2,g(b)))), f(f(g(skx0),f(sK2,sK3)),skx0) == g(f(skx0,f(sK2,g(skx0)))) }),
        clause({ f(f(g(b),f(sK2,sK3)),b) != g(f(b,f(sK2,g(b)))), f(f(g(r(skx0)),f(sK2,sK3)),r(skx0)) != g(f(r(skx0),f(sK2,g(r(skx0))))) }),

        // sK2
        clause({ f(f(g(sK1),f(b,sK3)),sK1) != g(f(sK1,f(b,g(sK1)))), f(f(g(sK1),f(skx1,sK3)),sK1) == g(f(sK1,f(skx1,g(sK1)))) }),
        clause({ f(f(g(sK1),f(b,sK3)),sK1) != g(f(sK1,f(b,g(sK1)))), f(f(g(sK1),f(r(skx1),sK3)),sK1) != g(f(sK1,f(r(skx1),g(sK1)))) }),

        // sK3
        clause({ f(f(g(sK1),f(sK2,b)),sK1) != g(f(sK1,f(sK2,g(sK1)))), f(f(g(sK1),f(sK2,skx3)),sK1) == g(f(sK1,f(sK2,g(sK1)))) }),
        clause({ f(f(g(sK1),f(sK2,b)),sK1) != g(f(sK1,f(sK2,g(sK1)))), f(f(g(sK1),f(sK2,r(skx3))),sK1) != g(f(sK1,f(sK2,g(sK1)))) }),

        // g(sK1)
        clause({ f(f(b,f(sK2,sK3)),sK1) != g(f(sK1,f(sK2,b))), f(f(skx4,f(sK2,sK3)),sK1) == g(f(sK1,f(sK2,skx4))) }),
        clause({ f(f(b,f(sK2,sK3)),sK1) != g(f(sK1,f(sK2,b))), f(f(r(skx4),f(sK2,sK3)),sK1) != g(f(sK1,f(sK2,r(skx4)))) }),

        // f(sK2,sK3)
        clause({ f(f(g(sK1),b),sK1) != g(f(sK1,f(sK2,g(sK1)))), f(f(g(sK1),skx5),sK1) == g(f(sK1,f(sK2,g(sK1)))) }),
        clause({ f(f(g(sK1),b),sK1) != g(f(sK1,f(sK2,g(sK1)))), f(f(g(sK1),r(skx5)),sK1) != g(f(sK1,f(sK2,g(sK1)))) }),

        // f(g(sK1),f(sK2,sK3))
        clause({ f(b,sK1) != g(f(sK1,f(sK2,g(sK1)))), f(skx6,sK1) == g(f(sK1,f(sK2,g(sK1)))) }),
        clause({ f(b,sK1) != g(f(sK1,f(sK2,g(sK1)))), f(r(skx6),sK1) != g(f(sK1,f(sK2,g(sK1)))) }),

        // f(f(g(sK1),f(sK2,sK3)),sK1)
        clause({ b != g(f(sK1,f(sK2,g(sK1)))), skx7 == g(f(sK1,f(sK2,g(sK1)))) }),
        clause({ b != g(f(sK1,f(sK2,g(sK1)))), r(skx7) != g(f(sK1,f(sK2,g(sK1)))) }),

        // f(sK2,g(sK1))
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != g(f(sK1,b)), f(f(g(sK1),f(sK2,sK3)),sK1) == g(f(sK1,skx8)) }),
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != g(f(sK1,b)), f(f(g(sK1),f(sK2,sK3)),sK1) != g(f(sK1,r(skx8))) }),

        // f(sK1,f(sK2,g(sK1)))
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != g(b), f(f(g(sK1),f(sK2,sK3)),sK1) == g(skx9) }),
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != g(b), f(f(g(sK1),f(sK2,sK3)),sK1) != g(r(skx9)) }),

        // g(f(sK1,f(sK2,g(sK1))))
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != b, f(f(g(sK1),f(sK2,sK3)),sK1) == skx10 }),
        clause({ f(f(g(sK1),f(sK2,sK3)),sK1) != b, f(f(g(sK1),f(sK2,sK3)),sK1) != r(skx10) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 10),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 10) })
    )

// positive literals are considered 1
TEST_GENERATION_INDUCTION(test_09,
    Generation::AsymmetricTest()
      .options({ { "induction_neg_only", "off" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  p(sK1) }))
      .expected({
        clause({ p(b), ~p(skx0), }),
        clause({ p(b), p(r(skx0)), }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 1) })
    )

// positive literals are considered 2
TEST_GENERATION_INDUCTION(test_10,
    Generation::AsymmetricTest()
      .options({ { "induction_neg_only", "off" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  sK1 == g(sK1) }))
      .expected({
        clause({ b == g(b), skx0 != g(skx0), }),
        clause({ b == g(b), r(skx0) == g(r(skx0)), }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 1) })
    )

// non-unit clauses are considered
TEST_GENERATION_INDUCTION(test_11,
    Generation::AsymmetricTest()
      .options({ { "induction_unit_only", "off" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  sK1 != g(sK1), p(g(sK2)), ~p(f(sK3,sK4)) }))
      .expected({
        // 1. literal sK1
        clause({ b != g(b), skx0 == g(skx0), p(g(sK2)), ~p(f(sK3,sK4)) }),
        clause({ b != g(b), r(skx0) != g(r(skx0)), p(g(sK2)), ~p(f(sK3,sK4)) }),

        // 3. literal sK3
        clause({ ~p(f(b,sK4)), p(f(skx1,sK4)), p(g(sK2)), sK1 != g(sK1) }),
        clause({ ~p(f(b,sK4)), ~p(f(r(skx1),sK4)), p(g(sK2)), sK1 != g(sK1) }),

        // 3. literal sK4
        clause({ ~p(f(sK3,b)), p(f(sK3,skx2)), p(g(sK2)), sK1 != g(sK1) }),
        clause({ ~p(f(sK3,b)), ~p(f(sK3,r(skx2))), p(g(sK2)), sK1 != g(sK1) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 3),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 3) })
    )

// "same induction" (i.e. generalized literal is same) is not done twice
// but resulting clauses are resolved with both literals
//
// TODO: this should be done with two inputs rather than with a non-unit clause
TEST_GENERATION_INDUCTION(test_12,
    Generation::AsymmetricTest()
      .options({ { "induction_unit_only", "off" }, { "induction", "struct" } })
      .indices(getIndices())
      .input( clause({  sK1 != g(sK1), sK2 != g(sK2) }))
      .expected({
        clause({ b != g(b), skx0 == g(skx0), sK2 != g(sK2) }),
        clause({ b != g(b), r(skx0) != g(r(skx0)), sK2 != g(sK2) }),

        clause({ b != g(b), skx0 == g(skx0), sK1 != g(sK1) }),
        clause({ b != g(b), r(skx0) != g(r(skx0)), sK1 != g(sK1) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->structInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->structInduction, 1) })
    )

// upward infinite interval integer induction
TEST_GENERATION_INDUCTION(test_13,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~pi(sK6) }) )
      .expected({
        clause({ ~pi(1), ~(skx0 < num(1)) }),
        clause({ ~pi(1), pi(skx0) }),
        clause({ ~pi(1), ~pi(skx0+1) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1) })
    )

// use bounds for upward+downward infinite interval integer induction
TEST_GENERATION_INDUCTION(test_14,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" }, { "int_induction_interval", "infinite" } })
      .context({ clause({ ~(sK6 < num(1)) }), clause({ ~(bi < sK6) }) })
      .indices(getIndices())
      .input( clause({ ~pi(sK6) }) )
      .expected({
        // upward induction
        clause({ ~pi(1), ~(skx0 < num(1)) }),
        clause({ ~pi(1), pi(skx0) }),
        clause({ ~pi(1), ~pi(skx0+1) }),

        // downard induction
        clause({ ~pi(bi), ~(bi < skx1) }),
        clause({ ~pi(bi), pi(skx1) }),
        clause({ ~pi(bi), ~pi(skx1+num(-1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 1) })
    )

// use bounds for upward+downward finite interval integer induction
TEST_GENERATION_INDUCTION(test_15,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" }, { "int_induction_interval", "finite" } })
      .context({ clause({ ~(sK6 < num(1)) }), clause({ ~(bi < sK6) }) })
      .indices(getIndices())
      .input( clause({ ~pi(sK6) }) )
      .expected({
        // upward induction
        clause({ ~pi(1), ~(skx0 < num(1)) }),
        clause({ ~pi(1), skx0 < bi }),
        clause({ ~pi(1), pi(skx0) }),
        clause({ ~pi(1), ~pi(skx0+1) }),

        // downard induction
        clause({ ~pi(bi), num(1) < skx1 }),
        clause({ ~pi(bi), ~(bi < skx1) }),
        clause({ ~pi(bi), pi(skx1) }),
        clause({ ~pi(bi), ~pi(skx1+num(-1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intFinUpInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intFinDownInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->intFinUpInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intFinDownInduction, 1) })
    )

// use default bound for downward integer induction,
// but for upward use the bound from index
TEST_GENERATION_INDUCTION(test_16,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" },
                 { "int_induction_interval", "infinite" },
                 { "int_induction_default_bound", "on" } })
      .context({ clause({ ~(sK6 < num(0)) }) })
      .indices(getIndices())
      .input( clause({ ~pi(sK6) }) )
      .expected({
        // upward induction
        clause({ ~pi(0), ~(skx0 < num(0)) }),
        clause({ ~pi(0), pi(skx0) }),
        clause({ ~pi(0), ~pi(skx0+1) }),

        // upward induction with default bound
        clause({ ~pi(0), ~(skx1 < num(0)), sK6 < 0 }),
        clause({ ~pi(0), pi(skx1), sK6 < 0 }),
        clause({ ~pi(0), ~pi(skx1+1), sK6 < 0 }),

        // downward induction with default bound
        clause({ ~pi(0), ~(num(0) < skx2), 0 < sK6 }),
        clause({ ~pi(0), pi(skx2), 0 < sK6 }),
        clause({ ~pi(0), ~pi(skx2+num(-1)), 0 < sK6 }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intDBUpInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intDBDownInduction, 0), })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 3),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intDBUpInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intDBDownInduction, 1) })
    )

// upward infinite interval induction triggered by the comparison literal
TEST_GENERATION_INDUCTION(test_17,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~pi(sK6) }) })
      .indices(getIndices())
      .input( clause({ ~(sK6 < num(1)) }) )
      .expected({
        clause({ ~pi(1), ~(skx0 < num(1)) }),
        clause({ ~pi(1), pi(skx0) }),
        clause({ ~pi(1), ~pi(skx0+1) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1) })
    )

// infinite+finite downward interval induction triggered by the comparison literal
TEST_GENERATION_INDUCTION(test_18,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~pi(sK6) }), clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ sK6 < bi }) )
      .expected({
        // infinite induction
        clause({ ~pi(bi), ~(bi < skx0) }),
        clause({ ~pi(bi), pi(skx0) }),
        clause({ ~pi(bi), ~pi(skx0+num(-1)) }),

        // finite induction
        clause({ ~pi(bi), ~(bi < skx1) }),
        clause({ ~pi(bi), num(1) < skx1 }),
        clause({ ~pi(bi), pi(skx1) }),
        clause({ ~pi(bi), ~pi(skx1+num(-1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intFinDownInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intFinDownInduction, 1) })
    )

// given the default strictness, induction is not applied on an interpreted constant
// (any strictness with term strictness != none works the same)
TEST_GENERATION_INDUCTION(test_19,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~pi(1) }) )
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// given a suitable strictness, induction is applied on an interpreted constant
// (any strictness with term strictness = none works the same)
TEST_GENERATION_INDUCTION(test_20,
    Generation::AsymmetricTest()
      .options({
        { "induction", "int" },
        { "int_induction_strictness_eq",   "always" },
        { "int_induction_strictness_comp", "always" },
        { "int_induction_strictness_term", "none" }
      })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~pi(1) }) )
      .expected({
        clause({ ~pi(sK6), ~(sK6 < skx0) }),
        clause({ ~pi(sK6), pi(skx0) }),
        clause({ ~pi(sK6), ~pi(skx0+num(-1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 1) })
    )

// given a suitable strictness, induction is applied on a term occuring only
// as one of the top-level arguments of "<"
// (any strictness with comparison strictness = none, term strictness in {none, interpreted_constant} works the same)
TEST_GENERATION_INDUCTION(test_21,
    Generation::AsymmetricTest()
      .options({
        { "induction", "int" },
        { "int_induction_strictness_eq",   "always" },
        { "int_induction_strictness_comp", "none" },
      })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~(bi < sK6) }) )
      .expected({
        // input used as main literal
        clause({ ~(bi < num(1)), ~(skx0 < num(1)) }),
        clause({ ~(bi < num(1)), bi < skx0 }),
        clause({ ~(bi < num(1)), ~(bi < skx0+1) }),
        // context used as main literal
        clause({ ~(bi < num(1)), ~(bi < skx1) }),
        clause({ ~(bi < num(1)), skx1 < num(1) }),
        clause({ ~(bi < num(1)), ~(skx1+num(-1) < num(1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 2),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfDownInduction, 1) })
    )

// given the default strictness, induction is applied on a term occuring in only
// one of the arguments of "<", but not to a term occuring only as a top-level
// argument of "<" (the "sK6" in context)
// (any strictness with comparison strictness != none, term strictness in {none, interpreted_constant} works the same)
TEST_GENERATION_INDUCTION(test_22,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~(bi < gi(sK6)) }) )
      .expected({
        clause({ ~(bi < gi(1)), ~(skx0 < num(1)) }),
        clause({ ~(bi < gi(1)), bi < gi(skx0) }),
        clause({ ~(bi < gi(1)), ~(bi < gi(skx0+1)) }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1) })
    )

// given the default suitable strictness, no induction is applied on a term occuring only
// as one of the top-level arguments of "<"
// (any strictness with comparison strictness != none, term strictness in {none, interpreted_constant} works the same)
TEST_GENERATION_INDUCTION(test_23,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ ~(bi < sK6) }) )
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// given the default strictness, induction is applied on a term occuring only
// as one of the top-level arguments of "="
// (any strictness with equality strictness != none, term strictness in {none, interpreted_constant} works the same)
TEST_GENERATION_INDUCTION(test_24,
    Generation::AsymmetricTest()
      .options({ { "induction", "int" } })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ bi != sK6 }) )
      .expected({
        clause({ bi != num(1), ~(skx0 < num(1)) }),
        clause({ bi != num(1), bi == skx0 }),
        clause({ bi != num(1), bi != skx0+1 }),
      })
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0),
                       TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 1),
                        TEST_FN_ASS_EQ(env.statistics->intInfUpInduction, 1) })
    )

// given a suitable strictness, no induction is applied on a term occuring only
// as one of the top-level arguments of "="
// (any strictness with equality strictness != none works the same)
TEST_GENERATION_INDUCTION(test_25,
    Generation::AsymmetricTest()
      .options({
        { "induction", "int" },
        { "int_induction_strictness_eq",   "toplevel_not_in_other" },
        { "int_induction_strictness_comp", "none" },
        { "int_induction_strictness_term", "none" }
      })
      .context({ clause({ ~(sK6 < num(1)) }) })
      .indices(getIndices())
      .input( clause({ bi != sK6 }) )
      .expected(none())
      .preConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
      .postConditions({ TEST_FN_ASS_EQ(env.statistics->inductionApplication, 0) })
    )

// all skolems are replaced when the hypothesis strengthening options is on, sik=one
TEST_GENERATION_INDUCTION(test_26,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" },
                 { "induction_strengthen_hypothesis", "on" } })
      .indices(getIndices())
      .input( clause({ f(sK1,sK2) != g(sK3) }) )
      .expected({
        // sK1
        clause({ f(b,skx0) != g(skx1), f(skx2,x) == g(y) }),
        clause({ f(b,skx0) != g(skx1), f(r(skx2),skx3) != g(skx4) }),

        // sK2
        clause({ f(skx5,b) != g(skx6), f(z,skx7) == g(x3) }),
        clause({ f(skx5,b) != g(skx6), f(skx8,r(skx7)) != g(skx9) }),

        // sK3
        clause({ f(skx10,skx11) != g(b), f(x4,x5) == g(skx12) }),
        clause({ f(skx10,skx11) != g(b), f(skx13,skx14) != g(r(skx12)) }),
      })
    )

// all skolems are replaced when the hypothesis strengthening options is on, sik=two
TEST_GENERATION_INDUCTION(test_27,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" }, { "structural_induction_kind", "two" },
                 { "induction_strengthen_hypothesis", "on" } })
      .indices(getIndices())
      .input( clause({ f(sK1,sK2) != g(sK3) }) )
      .expected({
        // sK1
        clause({ skx0 != r(r0(skx0)), f(r0(skx0),x) == g(y) }),
        clause({ f(skx0,skx1) != g(skx2) }),

        // sK2
        clause({ skx3 != r(r0(skx3)), f(z,r0(skx3)) == g(x3) }),
        clause({ f(skx4,skx3) != g(skx5) }),

        // sK3
        clause({ skx6 != r(r0(skx6)), f(x4,x5) == g(r0(skx6)) }),
        clause({ f(skx7,skx8) != g(skx6) }),
      })
    )

// multi-clause use case 1 (induction depth 0)
TEST_GENERATION_INDUCTION(test_28,
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" }, { "non_unit_induction", "on" } })
      .context({ clause({ p(sK1) })})
      .indices(getIndices())
      .input( clause({ sK2 != g(f(sK1,sK1)) }))
      .expected({
        // sK2
        clause({ b != g(f(sK1,sK1)), skx0 == g(f(sK1,sK1)) }),
        clause({ b != g(f(sK1,sK1)), r(skx0) != g(f(sK1,sK1)) }),

        // sK1 multiple literals
        clause({ sK2 != g(f(b,b)), sK2 == g(f(skx1,skx1)), ~p(skx1) }),
        clause({ sK2 != g(f(b,b)), p(r(skx1)) }),
        clause({ sK2 != g(f(b,b)), sK2 != g(f(r(skx1),r(skx1))) }),
        clause({ p(b), sK2 == g(f(skx1,skx1)), ~p(skx1) }),
        clause({ p(b), p(r(skx1)) }),
        clause({ p(b), sK2 != g(f(r(skx1),r(skx1))) }),

        // sK1 single literal
        clause({ sK2 != g(f(b,b)), sK2 == g(f(skx2,skx2)) }),
        clause({ sK2 != g(f(b,b)), sK2 != g(f(r(skx2),r(skx2))) }),
      })
    )

// multi-clause use case 2 (induction depth non-0)
TEST_GENERATION_INDUCTION(test_29,
    Generation::AsymmetricTest()
      .options({
        { "induction_on_complex_terms", "on" },
        { "induction", "struct" },
        { "non_unit_induction", "on" }
      })
      .context({
        fromInduction(clause({ p(g(sK3)) })),
        fromInduction(clause({ ~p(f(sK4,sK3)) }))
      })
      .indices(getIndices())
      .input( fromInduction(clause({ ~p(f(g(sK3),f(sK4,sK3))) })) )
      .expected({
        // g(sK3) multiple literals
        clause({ ~p(f(b,f(sK4,sK3))), p(f(skx0,f(sK4,sK3))), ~p(skx0) }),
        clause({ ~p(f(b,f(sK4,sK3))), ~p(f(r(skx0),f(sK4,sK3))) }),
        clause({ ~p(f(b,f(sK4,sK3))), p(r(skx0)) }),
        clause({ p(b), p(f(skx0,f(sK4,sK3))), ~p(skx0) }),
        clause({ p(b), ~p(f(r(skx0),f(sK4,sK3))) }),
        clause({ p(b), p(r(skx0)) }),

        // g(sK3) single literal
        clause({ ~p(f(b,f(sK4,sK3))), p(f(skx1,f(sK4,sK3))) }),
        clause({ ~p(f(b,f(sK4,sK3))), ~p(f(r(skx1),f(sK4,sK3))) }),

        // sK3
        clause({ ~p(f(g(b),f(sK4,b))), p(f(g(skx2),f(sK4,skx2))) }),
        clause({ ~p(f(g(b),f(sK4,b))), ~p(f(g(r(skx2)),f(sK4,r(skx2)))) }),

        // sK4
        clause({ ~p(f(g(sK3),f(b,sK3))), p(f(g(sK3),f(skx3,sK3))) }),
        clause({ ~p(f(g(sK3),f(b,sK3))), ~p(f(g(sK3),f(r(skx3),sK3))) }),

        // f(sK4,sK3) multiple literals
        clause({ ~p(f(g(sK3),b)), p(f(g(sK3),skx4)), p(skx4) }),
        clause({ ~p(f(g(sK3),b)), ~p(f(g(sK3),r(skx4))) }),
        clause({ ~p(f(g(sK3),b)), ~p(r(skx4)) }),
        clause({ ~p(b), p(f(g(sK3),skx4)), p(skx4) }),
        clause({ ~p(b), ~p(f(g(sK3),r(skx4))) }),
        clause({ ~p(b), ~p(r(skx4)) }),

        // f(sK4,sK3) single literal
        clause({ ~p(f(g(sK3),b)), p(f(g(sK3),skx5)) }),
        clause({ ~p(f(g(sK3),b)), ~p(f(g(sK3),r(skx5))) }),

        // f(g(sK3),f(sK4,sK3))
        clause({ ~p(b), p(skx6) }),
        clause({ ~p(b), ~p(r(skx6)) }),
      })
    )

// multi-clause use case 2 (main literal is from index)
TEST_GENERATION_INDUCTION(test_30,
    Generation::AsymmetricTest()
      .options({
        { "induction_on_complex_terms", "on" },
        { "induction", "struct" },
        { "non_unit_induction", "on" },
      })
      .context({
        fromInduction(clause({ ~p(g(g(sK3))) }))
      })
      .indices(getIndices())
      .input( fromInduction(clause({ ~p(g(sK3)) })) )
      .expected({
        // g(sK3) given clause
        clause({ ~p(b), p(skx0) }),
        clause({ ~p(b), ~p(r(skx0)) }),

        // sK3 given clause
        clause({ ~p(g(b)), p(g(skx1)) }),
        clause({ ~p(g(b)), ~p(g(r(skx1))) }),

        // sK3 multi-clause
        clause({ ~p(b), ~p(g(r(skx2))) }),
        clause({ ~p(b), ~p(r(skx2)) }),
        clause({ ~p(g(b)), ~p(g(r(skx2))) }),
        clause({ ~p(g(b)), ~p(r(skx2)) }),
        clause({ ~p(b), p(skx2), p(g(skx2)) }),
        clause({ ~p(g(b)), p(skx2), p(g(skx2)) }),
      })
    )

// multi-clause use case 1 (induction depth 0), non-unit
TEST_GENERATION_INDUCTION(test_31,
    Generation::AsymmetricTest()
      .options({
        { "induction_unit_only", "off" },
        { "induction", "struct" },
        { "non_unit_induction", "on" }, })
      .indices(getIndices())
      .input( clause({ ~p(sK1), ~p1(sK1) }) )
      .expected({
        // sK1, first literal
        clause({ ~p(b), p(skx0), ~p1(sK1) }),
        clause({ ~p(b), ~p(r(skx0)), ~p1(sK1) }),

        // sK1, second literal
        clause({ ~p1(b), p1(skx1), ~p(sK1) }),
        clause({ ~p1(b), ~p1(r(skx1)), ~p(sK1) }),

        // sK1, both literals, triggered by ~p(sK1)
        clause({ ~p(b), ~p1(b), p(skx2) }),
        clause({ ~p(b), ~p1(b), p1(skx2) }),
        clause({ ~p(b), ~p1(b), ~p(r(skx2)), ~p1(r(skx2)) }),

        // sK1, both literals, triggered by ~p1(sK1) (same as above)
        clause({ ~p(b), ~p1(b), p(skx2) }),
        clause({ ~p(b), ~p1(b), p1(skx2) }),
        clause({ ~p(b), ~p1(b), ~p(r(skx2)), ~p1(r(skx2)) }),
      })
    )

// multi-clause does not add tautological clauses
TEST_GENERATION_INDUCTION(test_32,
    // TODO enable multi-clause option when there is one
    Generation::AsymmetricTest()
      .options({ { "induction", "struct" }, { "non_unit_induction", "on" } })
      .context({ clause({ p(sK1) }) })
      .indices(getIndices())
      .input( clause({ ~p(sK1) }) )
      .expected({
        // sK1, given clause
        clause({ ~p(b), p(skx0) }),
        clause({ ~p(b), ~p(r(skx0)) }),
      })
    )

// multi-clause generalized occurrences
TEST_GENERATION_INDUCTION(test_33,
    Generation::AsymmetricTest()
      .options({
        { "induction", "struct" },
        { "induction_gen", "on" },
        { "non_unit_induction", "on" },
      })
      .context({ clause({ sK1 == sK2 }) })
      .indices(getIndices())
      .input( clause({ ~p(f(sK1,sK1)) }) )
      .expected({
        // sK1, given clause 01
        clause({ ~p(f(sK1,b)), p(f(sK1,skx0)) }),
        clause({ ~p(f(sK1,b)), ~p(f(sK1,r(skx0))) }),

        // sK1, given clause 10
        clause({ ~p(f(b,sK1)), p(f(skx1,sK1)) }),
        clause({ ~p(f(b,sK1)), ~p(f(r(skx1),sK1)) }),

        // sK1, given clause 11
        clause({ ~p(f(b,b)), p(f(skx2,skx2)) }),
        clause({ ~p(f(b,b)), ~p(f(r(skx2),r(skx2))) }),

        // sK1, multi-clause 101
        clause({ b == sK2, p(f(sK1,skx3)), skx3 != sK2 }),
        clause({ b == sK2, ~p(f(sK1,r(skx3))) }),
        clause({ b == sK2, r(skx3) == sK2 }),
        clause({ ~p(f(sK1,b)), p(f(sK1,skx3)), skx3 != sK2 }),
        clause({ ~p(f(sK1,b)), ~p(f(sK1,r(skx3))) }),
        clause({ ~p(f(sK1,b)), r(skx3) == sK2 }),

        // sK1, multi-clause 110
        clause({ b == sK2, p(f(skx4,sK1)), skx4 != sK2 }),
        clause({ b == sK2, ~p(f(r(skx4),sK1)) }),
        clause({ b == sK2, r(skx4) == sK2 }),
        clause({ ~p(f(b,sK1)), p(f(skx4,sK1)), skx4 != sK2 }),
        clause({ ~p(f(b,sK1)), ~p(f(r(skx4),sK1)) }),
        clause({ ~p(f(b,sK1)), r(skx4) == sK2 }),

        // sK1, multi-clause 111
        clause({ b == sK2, p(f(skx5,skx5)), skx5 != sK2 }),
        clause({ b == sK2, ~p(f(r(skx5),r(skx5))) }),
        clause({ b == sK2, r(skx5) == sK2 }),
        clause({ ~p(f(b,b)), p(f(skx5,skx5)), skx5 != sK2 }),
        clause({ ~p(f(b,b)), ~p(f(r(skx5),r(skx5))) }),
        clause({ ~p(f(b,b)), r(skx5) == sK2 }),
      })
    )

//
// FUNCTION DEFINITION TESTS
//

auto setup = [](SaturationAlgorithm& salg) {
  salg.getFunctionDefinitionHandler().initAndPreprocessLate(salg.getProblem(),salg.getOptions());
};

ClauseStack fnDefContext() {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);

  return {
    clause({ def_s(f(b,y),    y) }),
    clause({ def_s(f(r(x),b), f(x,x)) }),
    clause({ def_s(f(r(x),r(y)), h(f(x,g(y)),f(r(x),y))) }),

    clause({ def_s(h(b,y), y) }),
    clause({ def_s(h(r(x),y), h(x,y)) }),

    clause({ def_p(b()) }),
    clause({ ~def_p(r(r(x))), p(x), x != b() }),
  };
}

TEST_GENERATION_INDUCTION(test_34,
    Generation::AsymmetricTest()
      .setup(setup)
      .context(fnDefContext())
      .options({
        { "induction", "struct" },
        { "structural_induction_kind", "recursion" },
      })
      .indices(getIndices())
      .input( clause({ ~p(sK1) }) )
      .expected({
        clause({ ~p(b), ~p(r(b)), p(skx0) }),
        clause({ ~p(b), ~p(r(b)), ~p(r(r(skx0))) }),
      })
    )

TEST_GENERATION_INDUCTION(test_35,
    Generation::AsymmetricTest()
      .setup(setup)
      .context(fnDefContext())
      .options({
        { "induction", "struct" },
        { "structural_induction_kind", "recursion" },
      })
      .indices(getIndices())
      .input( clause({ f(sK1,sK2) != g(sK2) }) )
      .expected({
        clause({ f(b,skx0) != g(skx0), f(skx1,skx1) == g(skx1), f(skx2,g(skx3)) == g(g(skx3)) }),
        clause({ f(b,skx0) != g(skx0), f(r(skx1),b) != g(b), f(skx2,g(skx3)) == g(g(skx3)) }),
        clause({ f(b,skx0) != g(skx0), f(skx1,skx1) == g(skx1), f(r(skx2),skx3) == g(skx3) }),
        clause({ f(b,skx0) != g(skx0), f(r(skx1),b) != g(b), f(r(skx2),skx3) == g(skx3) }),
        clause({ f(b,skx0) != g(skx0), f(skx1,skx1) == g(skx1), f(r(skx2),r(skx3)) != g(r(skx3)) }),
        clause({ f(b,skx0) != g(skx0), f(r(skx1),b) != g(b), f(r(skx2),r(skx3)) != g(r(skx3)) }),
      })
    )

// only schemes from active positions are picked up
// and there is one generalization only for active occurrences
TEST_GENERATION_INDUCTION(test_36,
    Generation::AsymmetricTest()
      .setup(setup)
      .context(fnDefContext())
      .options({
        { "induction", "struct" },
        { "induction_on_active_occurrences", "on" },
      })
      .indices(getIndices())
      .input( clause({ h(h(sK1,sK2),f(sK3,sK4)) != g(sK4) }) )
      .expected({
        clause({ h(h(b,sK2),f(sK3,sK4)) != g(sK4), h(h(skx0,sK2),f(sK3,sK4)) == g(sK4) }),
        clause({ h(h(b,sK2),f(sK3,sK4)) != g(sK4), h(h(r(skx0),sK2),f(sK3,sK4)) != g(sK4) }),

        clause({ h(h(sK1,sK2),f(sK3,b)) != g(b), h(h(sK1,sK2),f(sK3,skx1)) == g(skx1) }),
        clause({ h(h(sK1,sK2),f(sK3,b)) != g(b), h(h(sK1,sK2),f(sK3,r(skx1))) != g(r(skx1)) }),

        clause({ h(h(sK1,sK2),f(sK3,sK4)) != g(b), h(h(sK1,sK2),f(sK3,sK4)) == g(skx1) }),
        clause({ h(h(sK1,sK2),f(sK3,sK4)) != g(b), h(h(sK1,sK2),f(sK3,sK4)) != g(r(skx1)) }),
      })
    )

// active occurrences are always inducted on in generalizations
TEST_GENERATION_INDUCTION(test_37,
    Generation::AsymmetricTest()
      .setup(setup)
      .context(fnDefContext())
      .options({
        { "induction", "struct" },
        { "induction_on_active_occurrences", "on" },
        { "induction_gen", "on" },
      })
      .indices(getIndices())
      .input( clause({ h(h(sK1,sK1),sK1) != g(sK1) }) )
      .expected({
        clause({ h(h(b,sK1),sK1) != g(b), h(h(skx0,sK1),sK1) == g(skx0) }),
        clause({ h(h(b,sK1),sK1) != g(b), h(h(r(skx0),sK1),sK1) != g(r(skx0)) }),

        clause({ h(h(b,b),sK1) != g(b), h(h(skx1,skx1),sK1) == g(skx1) }),
        clause({ h(h(b,b),sK1) != g(b), h(h(r(skx1),r(skx1)),sK1) != g(r(skx1)) }),

        clause({ h(h(b,sK1),b) != g(b), h(h(skx2,sK1),skx2) == g(skx2) }),
        clause({ h(h(b,sK1),b) != g(b), h(h(r(skx2),sK1),r(skx2)) != g(r(skx2)) }),

        clause({ h(h(b,b),b) != g(b), h(h(skx3,skx3),skx3) == g(skx3) }),
        clause({ h(h(b,b),b) != g(b), h(h(r(skx3),r(skx3)),r(skx3)) != g(r(skx3)) }),
      })
    )

//
// GENERALIZATION TESTS
//

// no generalization
TEST_FUN(generalizations_01) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  ContextReplacement it(inductionContext(sK1, {
    clause({ p(sK1) }),
    clause({ sK1 == sK2 }),
  }));

  assertContextReplacement(it, {
    inductionContext(sK1, {
      clause({ p(ph_s) }),
      clause({ ph_s == sK2 }),
    }),
  });
}

// test maximum subset size and generalizations
TEST_FUN(generalizations_02) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  ContextSubsetReplacement it(inductionContext(sK1, {
    clause({ p(f(sK1, sK1)) }),
    clause({ sK1 == f(sK2,sK1) }),
  }), 2);

  assertContextReplacement(it, {
    // 1 occurrence
    inductionContext(sK1, {
      clause({ p(f(ph_s,sK1)) }),
    }),
    inductionContext(sK1, {
      clause({ p(f(sK1,ph_s)) }),
    }),
    inductionContext(sK1, {
      clause({ ph_s == f(sK2,sK1) }),
    }),
    inductionContext(sK1, {
      clause({ sK1 == f(sK2,ph_s) }),
    }),
    // 2 occurrences
    inductionContext(sK1, {
      clause({ p(f(ph_s,ph_s)) }),
    }),
    inductionContext(sK1, {
      clause({ p(f(ph_s, sK1)) }),
      clause({ ph_s == f(sK2,sK1) }),
    }),
    inductionContext(sK1, {
      clause({ p(f(ph_s, sK1)) }),
      clause({ sK1 == f(sK2,ph_s) }),
    }),
    inductionContext(sK1, {
      clause({ p(f(sK1, ph_s)) }),
      clause({ ph_s == f(sK2,sK1) }),
    }),
    inductionContext(sK1, {
      clause({ p(f(sK1, ph_s)) }),
      clause({ sK1 == f(sK2,ph_s) }),
    }),
    inductionContext(sK1, {
      clause({ ph_s == f(sK2,ph_s) }),
    }),
    // 3 occurrences are missing, since it's more than
    // the maximum subset size and not all occurrences

    // all occurrences
    inductionContext(sK1, {
      clause({ p(f(ph_s, ph_s)) }),
      clause({ ph_s == f(sK2,ph_s) }),
    }),
  });
}

// no generalization if there are too many occurrences (20 currently)
TEST_FUN(generalizations_03) {
  __ALLOW_UNUSED(MY_SYNTAX_SUGAR);
  ContextSubsetReplacement it(inductionContext(sK1, {
    clause({ p(f(f(sK1,sK2), f(sK1,f(f(sK1,sK1),g(sK2))))) }),
    clause({ p1(f(f(f(sK1,sK1),sK2), f(sK1,sK1))) }),
    clause({ sK2 == f(f(f(f(g(sK1),sK1),f(sK1,sK2)),f(f(f(sK1,sK1),sK2),f(sK1,sK1))),
                      f(f(f(sK1,sK2),f(sK1,sK1)),f(f(sK1,g(sK2)),f(f(sK1,sK2),sK1)))) }),
  }), 0);

  // structure is preserved and all sK1 occurrences are replaced
  assertContextReplacement(it, {
    inductionContext(sK1, {
      clause({ p(f(f(ph_s,sK2), f(ph_s,f(f(ph_s,ph_s),g(sK2))))) }),
      clause({ p1(f(f(f(ph_s,ph_s),sK2), f(ph_s,ph_s))) }),
      clause({ sK2 == f(f(f(f(g(ph_s),ph_s),f(ph_s,sK2)),f(f(f(ph_s,ph_s),sK2),f(ph_s,ph_s))),
                        f(f(f(ph_s,sK2),f(ph_s,ph_s)),f(f(ph_s,g(sK2)),f(f(ph_s,sK2),ph_s)))) }),
    }),
  });
}

}
