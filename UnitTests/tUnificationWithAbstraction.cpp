/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#include "Forwards.hpp"
#include "Indexing/SubstitutionTree.hpp"
#include "Lib/Environment.hpp"

#include "Shell/Options.hpp"
#include "Test/TestUtils.hpp"

#include "Kernel/Unit.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/OperatorType.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/UnificationWithAbstraction.hpp"

#include "Indexing/Index.hpp"
#include "Indexing/LiteralSubstitutionTree.hpp"
#include "Indexing/LiteralIndex.hpp"
#include "Indexing/TermSubstitutionTree.hpp"
#include "Indexing/TermIndex.hpp"

#include "Test/UnitTesting.hpp"
#include "Test/SyntaxSugar.hpp"
#if VZ3
#include "z3++.h"
#endif
#include <ios>

using namespace std;
using namespace Kernel;
using namespace Indexing;

Clause* unit(Literal* lit)
{ return clause({ lit }); }

static const auto tld = [](auto t) { return TermWithoutValue(TypedTermList(t)); };


unique_ptr<TermSubstitutionTree<TermWithoutValue>> getTermIndexHOL()
{ return std::make_unique<TermSubstitutionTree<TermWithoutValue>>(); }

unique_ptr<TermSubstitutionTree<TermWithoutValue>> getTermIndex()
{ return std::make_unique<TermSubstitutionTree<TermWithoutValue>>();
}

auto getLiteralIndex()
{ return std::make_unique<LiteralSubstitutionTree<LiteralClause>>(); }

template<class TermOrLit>
struct UnificationResultSpec {
  TermOrLit querySigma;
  TermOrLit resultSigma;
  Stack<Literal*> constraints;

  friend bool operator==(UnificationResultSpec const& l, UnificationResultSpec const& r)
  {
    return Test::TestUtils::eqModAC(l.querySigma, r.querySigma)
       &&  Test::TestUtils::eqModAC(l.resultSigma, r.resultSigma)
       &&  Test::TestUtils::permEq(l.constraints, r.constraints,
             [](auto& l, auto& r) { return Test::TestUtils::eqModAC(l,r); });
  }

  friend std::ostream& operator<<(std::ostream& out, UnificationResultSpec const& self)
  { 
    out << "{ querySigma = " << Test::pretty(self.querySigma) << ", resultSigma = " << Test::pretty(self.resultSigma) << ", cons = [ ";
    for (auto& c : self.constraints) {
      out << *c << ", ";
    }
    return out << "] }";
  }
};

using TermUnificationResultSpec    = UnificationResultSpec<TermList>;
using LiteralUnificationResultSpec = UnificationResultSpec<Literal*>;

void checkLiteralMatches(LiteralSubstitutionTree<LiteralClause>& index, Options::UnificationWithAbstraction uwa, bool fixedPointIteration, Literal* lit, Stack<LiteralUnificationResultSpec> expected)
{
  Stack<LiteralUnificationResultSpec> is;
  for (auto qr : iterTraits(index.getUwa(lit, /* complementary */ false, uwa, fixedPointIteration)) ) {

    is.push(LiteralUnificationResultSpec {
        .querySigma = qr.unifier->subs().apply(lit, /* result */ QUERY_BANK),
        .resultSigma = qr.unifier->subs().apply(qr.data->literal, /* result */ RESULT_BANK),
        .constraints = *qr.unifier->constr().literals(qr.unifier->subs()),
    });
  }
  if (Test::TestUtils::permEq(is, expected, [](auto& l, auto& r) { return l == r; })) {
    cout << "[  OK  ] " << *lit << endl;
  } else {
    cout << "[ FAIL ] " << *lit << endl;
    cout << "tree: " << Output::multiline(index, 1) << endl;
    cout << "query: " << *lit << endl;

    cout << "is:" << endl;
    for (auto& x : is)
      cout << "         " << x << endl;

    cout << "expected:" << endl;
    for (auto& x : expected)
      cout << "         " << x << endl;

    exit(-1);
  }
}

template<class F>
void checkTermMatchesWithUnifFun(TermSubstitutionTree<TermWithoutValue>& index, TypedTermList term, Stack<TermUnificationResultSpec> expected, F unifFun)
{
  Stack<TermUnificationResultSpec> is;
  for (auto qr : iterTraits(unifFun(index, term))) {
    is.push(TermUnificationResultSpec {
        .querySigma  = qr.unifier->subs().apply(term, /* result */ QUERY_BANK),
        .resultSigma = qr.unifier->subs().apply(qr.data->term, /* result */ RESULT_BANK),
        .constraints = *qr.unifier->constr().literals(qr.unifier->subs()),
    });
  }
  if (Test::TestUtils::permEq(is, expected, [](auto& l, auto& r) { return l == r; })) {
    cout << "[  OK  ] " << term << endl;
  } else {
    cout << "[ FAIL ] " << term << endl;
    cout << "tree: " << Output::multiline(index, 1) << endl;
    cout << "query: " << term << endl;

    cout << "is:" << endl;
    for (auto& x : is)
      cout << "         " << x << endl;

    cout << "expected:" << endl;
    for (auto& x : expected)
      cout << "         " << x << endl;

    exit(-1);
  }

}

void checkTermMatches(TermSubstitutionTree<TermWithoutValue>& index, Options::UnificationWithAbstraction uwa, bool fixedPointIteration, TypedTermList term, Stack<TermUnificationResultSpec> expected)
{
  return checkTermMatchesWithUnifFun(index, term, expected, 
      [&](auto& idx, auto t) { return idx.getUwa(term, uwa, fixedPointIteration); });
}


struct IndexTest {
  std::unique_ptr<TermSubstitutionTree<TermWithoutValue>> index;
  Options::UnificationWithAbstraction uwa;
  bool fixedPointIteration = false;
  Stack<TypedTermList> insert;
  TermSugar query;
  Stack<TermUnificationResultSpec> expected;

  void run() {
    for (auto x : this->insert) {
      index->insert(TermWithoutValue(x));
    }

    checkTermMatches(*this->index, uwa, fixedPointIteration, TypedTermList(this->query, this->query.sort()),this->expected);

  }
};


struct LiteralIndexTest {
  unique_ptr<LiteralSubstitutionTree<LiteralClause>> index;
  Options::UnificationWithAbstraction uwa;
  bool fixedPointIteration = false;
  Stack<Literal*> insert;
  Literal* query;
  Stack<LiteralUnificationResultSpec> expected;

  void run() {
    DECL_PRED(dummy, {})
    for (auto x : this->insert) {
      index->insert(LiteralClause{ x, unit(dummy()) });
    }

    checkLiteralMatches(*index, uwa, fixedPointIteration, query, expected);
  }
};


#define INT_SUGAR                                                                                   \
   __ALLOW_UNUSED(                                                                                  \
      DECL_DEFAULT_VARS                                                                             \
      DECL_VAR(x0, 0)                                                                               \
      DECL_VAR(x1, 1)                                                                               \
      DECL_VAR(x2, 2)                                                                               \
      DECL_VAR(x3, 3)                                                                               \
      NUMBER_SUGAR(Int)                                                                             \
      DECL_PRED(p, {Int})                                                                           \
      DECL_FUNC(f, {Int}, Int)                                                                      \
      DECL_FUNC(g, {Int}, Int)                                                                      \
      DECL_FUNC(f2, {Int, Int}, Int)                                                                \
      DECL_FUNC(g2, {Int, Int}, Int)                                                                \
      DECL_CONST(a, Int)                                                                            \
      DECL_CONST(b, Int)                                                                            \
      DECL_CONST(c, Int)                                                                            \
    )                                                                                               \
 

#define RUN_TEST(name, sugar, ...)                                                                  \
  TEST_FUN(name) {                                                                                  \
       __ALLOW_UNUSED(sugar)                                                                        \
       __VA_ARGS__.run();                                                                           \
  }                                                                                                 \

RUN_TEST(term_indexing_one_side_interp_01,
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        f(1 + num(1)),
        f(1 + a),
      },
      .query = f(x),
      .expected = { 

          TermUnificationResultSpec 
          { .querySigma  = f(1 + a),
            .resultSigma = f(1 + a),
            .constraints = Stack<Literal*>() }, 

          TermUnificationResultSpec 
          { .querySigma  = f(1 + num(1)),
            .resultSigma = f(1 + num(1)),
            .constraints = Stack<Literal*>() }, 

      },
    })


RUN_TEST(term_indexing_one_side_interp_02, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        f(1 + num(1)),
        f(1 + a),
      },
      .query = g(x),
      .expected = Stack<TermUnificationResultSpec>(),
    })
 
RUN_TEST(term_indexing_one_side_interp_03, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
      },
      .query = x.sort(Int),
      .expected = { 

        TermUnificationResultSpec 
        { .querySigma  = 1 + a,
          .resultSigma = 1 + a,
          .constraints = Stack<Literal*>() }, 

        TermUnificationResultSpec 
        { .querySigma  = 1 + num(1),
          .resultSigma = 1 + num(1),
          .constraints = Stack<Literal*>() }, 

      }
    })


RUN_TEST(term_indexing_one_side_interp_04, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
      },
      .query = b + 2,
      .expected = { 

        TermUnificationResultSpec 
        { .querySigma  = 2 + b,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != 2 + b, } },

        TermUnificationResultSpec 
        { .querySigma  = 2 + b,
          .resultSigma = 1 + num(1),
          .constraints = { 2 + b != 1 + num(1), } }, 

      }
    })



RUN_TEST(term_indexing_one_side_interp_04_b, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + a,
      },
      .query = 2 + a,
      .expected = { 

        TermUnificationResultSpec 
        { .querySigma  = 2 + a,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != 2 + a, } },


      }
    })


RUN_TEST(term_indexing_one_side_interp_04_c, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        f(1 + num(1)),
        f(1 + a),
      },
      .query = f( b + 2 ),
      .expected = { 

        TermUnificationResultSpec 
        { .querySigma  = f( 2 + b ),
          .resultSigma = f( 1 + a ),
          .constraints = { 1 + a != 2 + b, } },

        TermUnificationResultSpec 
        { .querySigma  = f( 2 + b ),
          .resultSigma = f( 1 + num(1) ),
          .constraints = { 2 + b != 1 + num(1), } }, 

      }
    })

RUN_TEST(term_indexing_one_side_interp_04_d, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        g(f(1 + num(1))),
        g(f(1 + a)),
      },
      .query = g(f( b + 2 )),
      .expected = { 

        TermUnificationResultSpec 
        { .querySigma  = g(f( 2 + b )),
          .resultSigma = g(f( 1 + a )),
          .constraints = { 1 + a != 2 + b, } },

        TermUnificationResultSpec 
        { .querySigma  = g(f( 2 + b )),
          .resultSigma = g(f( 1 + num(1) )),
          .constraints = { 2 + b != 1 + num(1), } }, 

      }
    })

RUN_TEST(term_indexing_one_side_interp_05, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
        a,
      },
      .query = b + 2, 
      .expected = {
        TermUnificationResultSpec 
        { .querySigma  = 2 + b,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != 2 + b, } },

        TermUnificationResultSpec 
        { .querySigma  = 2 + b,
          .resultSigma = 1 + num(1),
          .constraints = { 2 + b != 1 + num(1), } }, 

        TermUnificationResultSpec 
        { .querySigma  = 2 + b,
          .resultSigma = a,
          .constraints = { 2 + b != a, } }, 

      }
})


RUN_TEST(term_indexing_one_side_interp_06, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
        a,
      },
      .query = x.sort(Int),
      .expected = {
        TermUnificationResultSpec 
        { .querySigma  = 1 + a,
          .resultSigma = 1 + a,
          .constraints = Stack<Literal*>{} },

        TermUnificationResultSpec 
        { .querySigma  = 1 + num(1),
          .resultSigma = 1 + num(1),
          .constraints = Stack<Literal*>{} }, 

        TermUnificationResultSpec 
        { .querySigma  = a,
          .resultSigma = a,
          .constraints = Stack<Literal*>{} }, 

      }
})


RUN_TEST(term_indexing_one_side_interp_07, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
        a,
        f(x),
      },
      .query = f(a),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = f(a),
          .resultSigma = 1 + a,
          .constraints = { 1 + a != f(a), } },

        TermUnificationResultSpec 
        { .querySigma  = f(a),
          .resultSigma = 1 + num(1),
          .constraints = { f(a) != 1 + num(1), } }, 


        TermUnificationResultSpec 
        { .querySigma  = f(a),
          .resultSigma = f(a),
          .constraints = Stack<Literal*>{} }, 

      }
})

RUN_TEST(term_indexing_one_side_interp_08, 
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .fixedPointIteration = false,
      .insert = {
        1 + num(1),
        1 + a,
        a,
        f(x),
      },
      .query = 3 + a,
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = 3 + a,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != 3 + a, } },

        TermUnificationResultSpec 
        { .querySigma  = 3 + a,
          .resultSigma = 1 + num(1),
          .constraints = { 3 + a != 1 + num(1), } }, 

        TermUnificationResultSpec 
        { .querySigma  = 3 + a,
          .resultSigma = a,
          .constraints = { 3 + a != a, } }, 

        TermUnificationResultSpec 
        { .querySigma  = 3 + a,
          .resultSigma = f(x),
          .constraints = { 3 + a != f(x) } }, 

      }
})

TEST_FUN(term_indexing_poly_01)
{
  auto uwa = Options::UnificationWithAbstraction::ONE_INTERP;
  auto fixedPointIteration = false;
  auto index = getTermIndex();

  DECL_DEFAULT_VARS
  DECL_DEFAULT_SORT_VARS  
  NUMBER_SUGAR(Int)
  DECL_CONST(a, Int) 
  DECL_POLY_CONST(h, 1, alpha)
  DECL_SORT(A)

  index->insert(tld(1 + a ));
  index->insert(tld(h(Int)));

  checkTermMatches(*index, uwa, fixedPointIteration, h(alpha), Stack<TermUnificationResultSpec>{

        TermUnificationResultSpec 
        { .querySigma  = h(Int),
          .resultSigma = h(Int),
          .constraints = Stack<Literal*>{  } }, 

        TermUnificationResultSpec 
        { .querySigma  = h(Int),
          .resultSigma = 1 + a,
          .constraints = { 1 + a != h(Int), } }, 

      });

  checkTermMatches(*index, uwa, fixedPointIteration, h(A), Stack<TermUnificationResultSpec>{ });
}


#define POLY_INT_SUGAR                                                                              \
  DECL_DEFAULT_VARS                                                                                 \
  DECL_DEFAULT_SORT_VARS                                                                            \
  NUMBER_SUGAR(Int)                                                                                 \
  DECL_POLY_CONST(b, 1, alpha)                                                                      \
  DECL_POLY_CONST(a, 1, alpha)                                                                      \
  DECL_POLY_FUNC(f, 1, {alpha}, alpha)                                                              \
  DECL_SORT(A)                                                                                      \
  DECL_CONST(someA, A)                                                                              \


#define HOL_SUGAR(...)                                                                              \
  DECL_DEFAULT_VARS                                                                                 \
  DECL_DEFAULT_SORT_VARS                                                                            \
  NUMBER_SUGAR(Int)                                                                                 \
  DECL_SORT(srt)                                                                                    \
  __VA_ARGS__

 


RUN_TEST(hol_0101,
    HOL_SUGAR(
      DECL_FUNC(f3, {srt, srt, srt}, srt)
      DECL_CONST(f1, arrow(srt, srt))
      DECL_CONST(f2, arrow(srt, srt))
      DECL_CONST(h, arrow(arrow(srt, srt), srt))
    ),
    IndexTest {
      .index = getTermIndexHOL(),
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT,
      .insert = {
               f3(x          , x, ap(h, f1)),
      },
      .query = f3(ap(h, f2), y, y          ),
      .expected =  {

        TermUnificationResultSpec 
        { .querySigma  = f3(ap(h, f2), ap(h, f1), ap(h, f1)),
          .resultSigma = f3(ap(h, f1), ap(h, f1), ap(h, f1)),
          .constraints = { f1 != f2 } }, 

      }
    })


RUN_TEST(hol_0102,
    HOL_SUGAR(
      DECL_FUNC(f3, {srt, srt, srt}, srt)
      DECL_CONST(f1, arrow(srt, srt))
      DECL_CONST(f2, arrow(srt, srt))
      DECL_CONST(h, arrow(arrow(srt, srt), srt))
    ),
    IndexTest {
      .index = getTermIndexHOL(),
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT,
      .insert = {
               f3(ap(h, f2), y, y          ),
      },
      .query = f3(x          , x, ap(h, f1)),
      .expected =  {

        TermUnificationResultSpec 
        { .querySigma  = f3(ap(h, f1), ap(h, f1), ap(h, f1)),
          .resultSigma = f3(ap(h, f2), ap(h, f1), ap(h, f1)),
          .constraints = { f1 != f2 } }, 

      }
    })


RUN_TEST(hol_02,
    HOL_SUGAR(
      DECL_FUNC(f3, {srt, srt, srt}, srt)
      DECL_CONST(f1, arrow(srt, srt))
      DECL_CONST(f2, arrow(srt, srt))
      DECL_CONST(a, srt)
      DECL_CONST(h, arrow(arrow(srt, srt), srt))
      ),
    IndexTest {
      .index = getTermIndexHOL(),
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT,
      .insert = {
               f3(a          , x, ap(h, f1)),
               f3(x          , x, ap(h, f1)),
      },
      .query = f3(ap(h, f2), y, y          ),
      .expected =  {

        TermUnificationResultSpec 
        { .querySigma  = f3(ap(h, f2), ap(h, f1), ap(h, f1)),
          .resultSigma = f3(ap(h, f1), ap(h, f1), ap(h, f1)),
          .constraints = { f1 != f2 } }, 

      }
    })


RUN_TEST(hol_03,
    HOL_SUGAR(
      DECL_FUNC(f3, {srt, srt, srt}, srt)
      DECL_CONST(f1, arrow(srt, srt))
      DECL_CONST(f2, arrow(srt, srt))
      DECL_CONST(h1, arrow(arrow(srt, srt), srt))
      DECL_CONST(h2, arrow(arrow(srt, srt), srt))
    ),
    IndexTest {
      .index = getTermIndexHOL(),
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT,
      .insert = {
               ap(h1, f1),
               ap(h2, f1),
      },
      .query = ap(h1, f2),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = ap(h1, f2),
          .resultSigma = ap(h1, f1),
          .constraints = { f1 != f2 } }, 
      }
    })

#define RUN_TEST_hol_04(idx, ...)                                                                   \
  RUN_TEST(hol_04_ ## idx,                                                                          \
    HOL_SUGAR(                                                                                      \
      DECL_FUNC(f3, {srt, srt, srt}, srt)                                                           \
      DECL_POLY_CONST(c1, 1, alpha)                                                                 \
      DECL_POLY_CONST(c2, 1, alpha)                                                                 \
      DECL_POLY_CONST(h, 2, arrow(alpha, beta))                                                     \
    ),                                                                                              \
    IndexTest {                                                                                     \
      .index = getTermIndexHOL(),                                                \
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT, \
      .insert = {                                                                                   \
               ap(h(arrow(srt, srt), srt), c1(arrow(srt, srt))),                                    \
               ap(h(srt            , srt), c1(srt)),                                                \
      },                                                                                            \
      __VA_ARGS__                                                                                   \
    })


RUN_TEST_hol_04(01,
      .query = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .resultSigma = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .constraints = Stack<Literal*>{} }, 
      }
    )

RUN_TEST_hol_04(02,
      .query = ap(h(arrow(srt,srt), srt), c2(arrow(srt, srt))),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = ap(h(arrow(srt,srt), srt), c2(arrow(srt, srt))),
          .resultSigma = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .constraints = Stack<Literal*>{ c1(arrow(srt,srt)) != c2(arrow(srt,srt)) } }, 
      }
    )


#define RUN_TEST_hol_05(idx, ...)                                                                   \
  RUN_TEST(hol_05_ ## idx,                                                                          \
    HOL_SUGAR(                                                                                      \
      DECL_FUNC(f3, {srt, srt, srt}, srt)                                                           \
      DECL_POLY_CONST(c1, 1, alpha)                                                                 \
      DECL_POLY_CONST(c2, 1, alpha)                                                                 \
      DECL_POLY_CONST(h, 2, arrow(alpha, beta))                                                     \
    ),                                                                                              \
    IndexTest {                                                                                     \
      .index = getTermIndexHOL(),                                                \
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT, \
      .insert = {                                                                                   \
               ap(h(arrow(srt, srt), srt), c1(arrow(srt, srt))),                                    \
               ap(h(srt            , srt), c2(srt)),                                                \
      },                                                                                            \
      __VA_ARGS__                                                                                   \
    })


RUN_TEST_hol_05(01,
      .query = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .resultSigma = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .constraints = Stack<Literal*>{} }, 
      }
    )

RUN_TEST_hol_05(02,
      .query = ap(h(arrow(srt,srt), srt), c2(arrow(srt, srt))),
      .expected =  {
        TermUnificationResultSpec 
        { .querySigma  = ap(h(arrow(srt,srt), srt), c2(arrow(srt, srt))),
          .resultSigma = ap(h(arrow(srt,srt), srt), c1(arrow(srt, srt))),
          .constraints = Stack<Literal*>{ c1(arrow(srt,srt)) != c2(arrow(srt,srt)) } }, 
      }
    )

RUN_TEST(hol_06,
    HOL_SUGAR(
      DECL_SORT_BOOL;
      DECL_SORT(A)
      DECL_CONST(f, arrow(Bool, A))
      DECL_CONST(a, Bool)
      DECL_CONST(b, Bool)
    ),
    IndexTest {
      .index = getTermIndexHOL(),
      .uwa = Options::UnificationWithAbstraction::FUNC_EXT,
      .insert = {
               ap(f, a),
               ap(f, b),
               a,
               b
      },
      .query = ap(f,a),
      .expected =  {

        TermUnificationResultSpec 
        { .querySigma  = ap(f,a),
          .resultSigma = ap(f,a),
          .constraints = Stack<Literal*>{ } }, 

        TermUnificationResultSpec 
        { .querySigma  = ap(f,a),
          .resultSigma = ap(f,b),
          .constraints = { a != b } }, 

      }
    })


RUN_TEST(term_indexing_poly_uwa_01,
    POLY_INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::ONE_INTERP,
      .insert = {
        f(alpha, a(alpha)),
        f(alpha, b(alpha)),
        f(A, someA),
        f(A, a(A)),
      },
      .query = f(Int, a(Int) + x),
      .expected =  {

        TermUnificationResultSpec 
        { .querySigma  = f(Int, a(Int) + x),
          .resultSigma = f(Int, a(Int)),
          .constraints = { a(Int) != a(Int) + x } }, 

        TermUnificationResultSpec 
        { .querySigma  = f(Int, a(Int) + x),
          .resultSigma = f(Int, b(Int)),
          .constraints = { b(Int) != a(Int) + x } }, 

      }
    })


TEST_FUN(term_indexing_interp_only)
{
  auto uwa = Options::UnificationWithAbstraction::INTERP_ONLY;
  auto fixedPointIteration = false;
  auto index = getTermIndex();

  DECL_DEFAULT_VARS
  NUMBER_SUGAR(Int)

  DECL_CONST(a, Int) 
  DECL_CONST(b, Int) 

  index->insert(tld(num(1) + num(1)));
  index->insert(tld(1 + a          ));

  checkTermMatches(*index, uwa, fixedPointIteration, b + 2, {

        TermUnificationResultSpec 
        { .querySigma  = b + 2,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != b + 2, } }, 

        TermUnificationResultSpec 
        { .querySigma  = b + 2,
          .resultSigma = 1 + num(1),
          .constraints = { 1 + num(1) != b + 2, } }, 

      });

  index->insert(tld(a));

  checkTermMatches(*index, uwa, fixedPointIteration, b + 2 , {

        TermUnificationResultSpec 
        { .querySigma  = b + 2,
          .resultSigma = 1 + a,
          .constraints = { 1 + a != b + 2, } }, 

        TermUnificationResultSpec 
        { .querySigma  = b + 2,
          .resultSigma = 1 + num(1),
          .constraints = { 1 + num(1) != b + 2, } }, 

      });

}

TEST_FUN(literal_indexing)
{
  auto uwa = Options::UnificationWithAbstraction::ONE_INTERP;
  auto fixedPointIteration = false;
  auto index = getLiteralIndex();

  DECL_DEFAULT_VARS
  NUMBER_SUGAR(Int)
  DECL_PRED(p, {Int})

  DECL_CONST(a, Int) 
  DECL_CONST(b, Int) 

  Clause* dummy = unit(p(a));

  index->insert(LiteralClause{ p(num(1) + num(1)), dummy });
  index->insert(LiteralClause{ p(1 + a          ), dummy });

  checkLiteralMatches(*index, uwa, fixedPointIteration, p(b + 2), {

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(num(1) + 1),
        .constraints = { b + 2 != num(1) + 1 }, },

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(a + 1),
        .constraints = { b + 2 != a + 1 }, },

      });

  index->insert(LiteralClause{ p(b + 2),unit(p(b + 2)) });
  index->insert(LiteralClause{ p(2 + b),unit(p(2 + b)) });

  checkLiteralMatches(*index, uwa, fixedPointIteration, p(b + 2), {

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(num(1) + 1),
        .constraints = { b + 2 != num(1) + 1 }, },

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(a + 1),
        .constraints = { b + 2 != a + 1 }, },

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(b + 2),
        .constraints = Stack<Literal*>{  }, },

      LiteralUnificationResultSpec {
        .querySigma = p(b + 2),
        .resultSigma = p(b + 2),
        .constraints = Stack<Literal*>{ b + 2 != 2 + b }, },

      });


}

TEST_FUN(higher_order)
{

  DECL_DEFAULT_VARS
  DECL_DEFAULT_SORT_VARS  
  NUMBER_SUGAR(Int)
  DECL_SORT(srt) 
  DECL_CONST(a, arrow(srt, srt))
  DECL_CONST(b, arrow(srt, srt))
  DECL_CONST(c, srt)  
  DECL_CONST(f, arrow(arrow(srt, srt), srt))
  DECL_CONST(g, arrow(srt, arrow(srt, srt)))
  auto uwa = Options::UnificationWithAbstraction::FUNC_EXT;
  auto fixedPointIteration = false;
  auto index = getTermIndexHOL();

  index->insert(tld(ap(f,a)));

  checkTermMatches(*index, uwa, fixedPointIteration, ap(f,b), Stack<TermUnificationResultSpec>{

        TermUnificationResultSpec 
        { .querySigma  = ap(f,b),
          .resultSigma = ap(f,a),
          .constraints = { a != b, } }, 

      });

  index->insert(tld(ap(g,c)));
  index->insert(tld(g));

  checkTermMatches(*index, uwa, fixedPointIteration, TypedTermList(x, arrow(srt, srt)), Stack<TermUnificationResultSpec>{

        TermUnificationResultSpec 
        { .querySigma  = ap(g,c),
          .resultSigma = ap(g,c),
          .constraints = Stack<Literal*>{} }, 

      });

  checkTermMatches(*index, uwa, fixedPointIteration, ap(f,b), Stack<TermUnificationResultSpec>{

        TermUnificationResultSpec 
        { .querySigma  = ap(f,b),
          .resultSigma = ap(f,a),
          .constraints = { a != b, } }, 

      });



}

TEST_FUN(higher_order2)
{
  auto index = getTermIndexHOL();

  DECL_DEFAULT_VARS
  DECL_DEFAULT_SORT_VARS  
  NUMBER_SUGAR(Int)
  DECL_SORT(srt) 
  DECL_CONST(a, arrow(srt, srt))
  DECL_CONST(b, arrow(srt, srt))
  DECL_CONST(f, arrow({arrow(srt, srt), arrow(srt, srt)}, srt))

  index->insert(tld(ap(ap(f,a),b)));

}

Option<TermUnificationResultSpec> runRobUnify(TypedTermList a, TypedTermList b, Options::UnificationWithAbstraction opt, bool fixedPointIteration) {
  auto au = AbstractingUnifier::unify(a, 0, b, 0, AbstractionOracle(opt), fixedPointIteration);

  if (au) {
    return some(TermUnificationResultSpec { 
     .querySigma  = au->subs().apply(a, 0), 
     .resultSigma = au->subs().apply(b, 0), 
     .constraints = *au->computeConstraintLiterals(),
    });

  } else {
    return {};
  }


}

void checkRobUnify(TypedTermList a, TypedTermList b, Options::UnificationWithAbstraction opt, bool fixedPointIteration, TermUnificationResultSpec exp)
{
  auto is = runRobUnify(a,b,opt, fixedPointIteration);
  if (is.isSome() && is.unwrap() == exp) {
    cout << "[  OK  ] " << a << " unify " << b << endl;
  } else {
    cout << "[ FAIL ] " << a << " unify " << b << endl;
    cout << "is:       " << is << endl;
    cout << "expected: " << exp << endl;
    exit(-1);
  }
}


void checkRobUnifyFail(TypedTermList a, TypedTermList b, Options::UnificationWithAbstraction opt, bool fixedPointIteration)
{
  auto is = runRobUnify(a,b,opt, fixedPointIteration);
  if(is.isNone()) {
      cout << "[  OK  ] " << a << " unify " << b << endl;
  } else {
    cout << "[ FAIL ] " << a << " unify " << b << endl;
    cout << "is:       " << is << endl;
    cout << "expected: nothing" << endl;
    exit(-1);
  }
}

#define ROB_UNIFY_TEST(name, opt, fixedPointIteration, lhs, rhs, ...)                                          \
  TEST_FUN(name)                                                                                    \
  {                                                                                                 \
    INT_SUGAR                                                                                       \
    checkRobUnify(lhs, rhs, opt, fixedPointIteration, __VA_ARGS__ );                                           \
  }                                                                                                 \

#define ROB_UNIFY_TEST_FAIL(name, opt, fixedPointIteration, lhs, rhs)                                          \
  TEST_FUN(name)                                                                                    \
  {                                                                                                 \
    INT_SUGAR                                                                                       \
    checkRobUnifyFail(lhs, rhs, opt, fixedPointIteration);                                                     \
  }                                                                                                 \

ROB_UNIFY_TEST(rob_unif_test_01,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f(b + 2), 
    f(x + 2),
    TermUnificationResultSpec { 
      .querySigma = f(b + 2),
      .resultSigma = f(x + 2),
      .constraints = { x + 2 != b + 2 },
    })

ROB_UNIFY_TEST(rob_unif_test_02,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f(b + 2), 
    f(x + 2),
    TermUnificationResultSpec { 
      .querySigma = f(b + 2),
      .resultSigma = f(x + 2),
      .constraints = { x + 2 != b + 2 },
    })

ROB_UNIFY_TEST(rob_unif_test_03,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f(x + 2), 
    f(a),
    TermUnificationResultSpec { 
      .querySigma = f(x + 2),
      .resultSigma = f(a),
      .constraints = { x + 2 != a },
    })

ROB_UNIFY_TEST_FAIL(rob_unif_test_04,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f(a), g(1 + a))


ROB_UNIFY_TEST(rob_unif_test_05,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f(a + b), 
    f(x + y),
    TermUnificationResultSpec { 
      .querySigma = f(a + b),
      .resultSigma = f(x + y),
      .constraints = { x + y != a + b },
    })

ROB_UNIFY_TEST(rob_unif_test_06,
    Options::UnificationWithAbstraction::ONE_INTERP,
    /* withFixedPointIteration */ false,
    f2(x, x + 1), 
    f2(a, a),
    TermUnificationResultSpec { 
      .querySigma = f2(a, a + 1),
      .resultSigma = f2(a, a),
      .constraints = { a != a + 1 },
    })

ROB_UNIFY_TEST(over_approx_test_2_bad_AC1,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(x, a + x),
    f2(c, b + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(c, a + c),
      .resultSigma = f2(c, b + a),
      .constraints = { c != b },
    })

ROB_UNIFY_TEST_FAIL(over_approx_test_2_bad_AC1_fixedPointIteration,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ true,
    f2(x, a + x),
    f2(c, b + a)
    )

ROB_UNIFY_TEST_FAIL(over_approx_test_2_good_AC1,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(a + x, x),
    f2(b + a, c))

ROB_UNIFY_TEST(bottom_constraint_test_1_bad_AC1,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(f2(y, x), a + y + x),
    f2(f2(b, c), c + b + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(f2(b,c), a + b + c),
      .resultSigma = f2(f2(b,c), c + b + a),
      .constraints = Stack<Literal*>{ b + c != c + b },
    })

ROB_UNIFY_TEST(bottom_constraint_test_1_bad_AC1_fixedPointIteration,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ true,
    f2(f2(y, x), a + y + x),
    f2(f2(b, c), c + b + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(f2(b,c), a + b + c),
      .resultSigma = f2(f2(b,c), c + b + a),
      .constraints = Stack<Literal*>{  },
    })

ROB_UNIFY_TEST(bottom_constraint_test_1_good_AC1,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(a + x + y, f2(x, y)),
    f2(c + b + a, f2(b, c)),
    TermUnificationResultSpec { 
      .querySigma  = f2(a + b + c, f2(b,c)),
      .resultSigma = f2(c + b + a, f2(b,c)),
      .constraints = Stack<Literal*>{},
    })


ROB_UNIFY_TEST(ac_bug_01,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    a + b + c + a,
    a + b + x + y,
    TermUnificationResultSpec { 
      .querySigma  = a + b + c + a,
      .resultSigma = a + b + x + y,
      .constraints = { c + a != x + y },
    })

ROB_UNIFY_TEST(ac_test_01_AC1,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(b, a + b + c),
    f2(b, x + y + c),
    TermUnificationResultSpec { 
      .querySigma  = f2(b, a + b + c),
      .resultSigma = f2(b, x + y + c),
      .constraints = { a + b != x + y },
    })

ROB_UNIFY_TEST(ac_test_02_AC1_good,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(a + b + c, c),
    f2(x + y + z, z),
    TermUnificationResultSpec { 
      .querySigma  = f2(a + b + c, c),
      .resultSigma = f2(x + y + c, c),
      .constraints = { a + b != x + y },
    })

ROB_UNIFY_TEST(ac_test_02_AC1_bad,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ false,
    f2(c, a + b + c),
    f2(z, x + y + z),
    TermUnificationResultSpec { 
      .querySigma  = f2(c, a + b + c),
      .resultSigma = f2(c, x + y + c),
      .constraints = { a + b + c != x + y + c },
    })

ROB_UNIFY_TEST(ac_test_02_AC1_bad_fixedPointIteration,
    Options::UnificationWithAbstraction::AC1,
    /* withFixedPointIteration */ true,
    f2(c, a + b + c),
    f2(z, x + y + z),
    TermUnificationResultSpec { 
      .querySigma  = f2(c, a + b + c),
      .resultSigma = f2(c, x + y + c),
      .constraints = { a + b != x + y },
    })

ROB_UNIFY_TEST(ac2_test_01,
    Options::UnificationWithAbstraction::AC2,
    /* withFixedPointIteration */ false,
    f2(x, a + b + c),
    f2(x, x + b + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(c, a + b + c),
      .resultSigma = f2(c, c + b + a),
      .constraints = Stack<Literal*>{},
    })

ROB_UNIFY_TEST(ac2_test_02,
    Options::UnificationWithAbstraction::AC2,
    /* withFixedPointIteration */ false,
    f2(a + b + c, f2(x,b)),
    f2(x + y + a, f2(x,y)),
    TermUnificationResultSpec { 
      .querySigma  = f2(a + b + c, f2(c,b)),
      .resultSigma = f2(c + b + a, f2(c,b)),
      .constraints = Stack<Literal*>{},
    })

ROB_UNIFY_TEST(ac2_test_02_bad,
    Options::UnificationWithAbstraction::AC2,
    /* withFixedPointIteration */ false,
    f2(f2(x,b), a + b + c),
    f2(f2(x,y), x + y + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(f2(x,b), a + b + c),
      .resultSigma = f2(f2(x,b), x + b + a),
      .constraints = Stack<Literal*>{ b + c != x + b },
    })

ROB_UNIFY_TEST(ac2_test_02_bad_fixedPointIteration,
    Options::UnificationWithAbstraction::AC2,
    /* withFixedPointIteration */ true,
    f2(f2(x,b), a + b + c),
    f2(f2(x,y), x + y + a),
    TermUnificationResultSpec { 
      .querySigma  = f2(f2(c,b), a + b + c),
      .resultSigma = f2(f2(c,b), c + b + a),
      .constraints = Stack<Literal*>{ },
    })


ROB_UNIFY_TEST(top_level_constraints_1,
    Options::UnificationWithAbstraction::AC2,
    /* withFixedPointIteration */ false,
    a + y + x,
    a + b + c,
    TermUnificationResultSpec { 
      .querySigma  = a + y + x,
      .resultSigma = a + b + c,
      .constraints = Stack<Literal*>{ b + c != x + y },
    })

RUN_TEST(top_level_constraints_2_with_fixedPointIteration,
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::AC2,
      .fixedPointIteration = true,
      .insert = {
        a + b + c,
        b,
        a + b + f(a) + c,
        f(x),
        f(a),
      },
      .query = a + y + x,
      .expected = { 

          TermUnificationResultSpec 
          { .querySigma  = a + x0 + x1,
            .resultSigma = a + b + c,
            .constraints = Stack<Literal*>{ b + c != x1 + x0 } }, 

          TermUnificationResultSpec 
          { .querySigma  = a + x0 + x1,
            .resultSigma = a + b + f(a) + c,
            .constraints = Stack<Literal*>{ b + f(a) + c != x1 + x0 } }, 

      },
    })


RUN_TEST(top_level_constraints_2,
    INT_SUGAR,
    IndexTest {
      .index = getTermIndex(),
      .uwa = Options::UnificationWithAbstraction::AC2,
      .fixedPointIteration = false,
      .insert = {
        a + b + c,
        b,
        a + b + a + c,
        f(x),
        f(a),
      },
      .query = a + y + x,
      .expected = { 

          TermUnificationResultSpec 
          { .querySigma  = a + x0 + x1,
            .resultSigma = a + b + c,
            .constraints = Stack<Literal*>{ a + b + c != a + x1 + x0 } }, 

          TermUnificationResultSpec 
          { .querySigma  = a + x0 + x1,
            .resultSigma = a + b + a + c,
            .constraints = Stack<Literal*>{ a + b + a + c != a + x1 + x0 } }, 

      },
    })

RUN_TEST(literal_tree_test_01,
    INT_SUGAR,
    LiteralIndexTest {
      .index = getLiteralIndex(),
      .uwa = Options::UnificationWithAbstraction::AC2,
      .fixedPointIteration = true,
      .insert = {
        x < (a + (x + 1)),
        x < (x + a),
      },
      .query = x < 1 + (x + 1),
      .expected = Stack<LiteralUnificationResultSpec>{ 

      },
    })


RUN_TEST(literal_tree_test_02,
    INT_SUGAR,
    LiteralIndexTest {
      .index = getLiteralIndex(),
      .uwa = Options::UnificationWithAbstraction::AC1,
      .fixedPointIteration = true,
      .insert = {
        x < x,
        x < (y + 1),
        x < (-1 + (x + 1)),
        x < (x + -1),
        x < y,
        (x + 1) < x,
        (1  + (x + 1)) < x,
      },
      .query = x < 1 + (x + 1),
      .expected = Stack<LiteralUnificationResultSpec>{ 

          // x < y + 1
          LiteralUnificationResultSpec 
          { .querySigma  = x < 1 + (x + 1),
            .resultSigma = x < 1 + (x + 1),
            .constraints = Stack<Literal*>{} }, 

          // x < y
          LiteralUnificationResultSpec 
          { .querySigma  = x < 1 + (x + 1),
            .resultSigma = x < 1 + (x + 1),
            .constraints = Stack<Literal*>{} }, 

      },
    })
