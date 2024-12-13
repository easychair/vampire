/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __LIB_TYPELIST__H__
#define __LIB_TYPELIST__H__

#include <functional>

namespace Lib {

  template<class F, class... As> using ResultOf = typename std::invoke_result<F, As...>::type;

  template<unsigned v>
  struct Constant { static constexpr unsigned value = v; };

namespace TypeList {

  template<class A>
  struct Token { using inner = A; };

  template<class Token>
  using TokenType = typename Token::inner;

  template<class F>
  auto __forEach(F& f) { }

  template<class F, class A, class... As>
  auto __forEach(F& f) {
    f(Token<A>{});
    __forEach<F, As...>(f);
  }

  /* Meta level type constructor.
   *
   * List : [class] -> class
   */
  template<class... As> struct List
  {
    template<class F>
    static auto forEach(F f)
    { __forEach<F, As...>(f); }


    template<class F>
    static auto toTuple(F f)
    { return std::make_tuple(f(Token<As>{})...); }
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// FUNTIONS RETURNING TYPES
  ////////////////////////////////////////

  /* Meta level type function.
   * into : ([class] -> class) -> List [class] -> class
   *
   * Transforms a List<As...> into another type with the same variadic type arguments.
   * E.g. Into< Coproduct, List<As...>> ==>  Coproduct<As...>
   *      Into<std::tuple, List<As...>> ==> std::tuple<As...>
   */
  template<template<class...> class HKT, class A>
  struct IntoImpl;

  template<template<class...> class HKT, class... As>
  struct IntoImpl<HKT, List<As...> >
  { using type = HKT<As...>; };

  template<template<class...> class HKT, class A>
  using Into = typename IntoImpl<HKT, A>::type;

  /* Meta level type function.
   * concat : List [class] -> List [class] -> List [class]
   *
   * Concatenates two List<...>s.
   *
   * E.g. Concat<List<A, B>, List<A, B>> ==>  List<A,B,A,B>
   */
  template<class A, class B> struct ConcatImpl;

  template<class... As, class... Bs>
  struct ConcatImpl<List<As...>, List<Bs...>>
  { using type = List<As..., Bs...>; };

  template<class A, class B>
  using Concat = typename ConcatImpl<A,B>::type;

  /*
   * E.g. Repeat<4, A> ==>  List<A,A,A,A>
   */
  template<unsigned N, class A> struct RepeatImpl
  { using type = Concat<List<A>, typename RepeatImpl<N-1, A>::type>; };

  template<class A>
  struct RepeatImpl<0, A>
  { using type = List<>; };

  template<unsigned N, class A>
  using Repeat = typename RepeatImpl<N,A>::type;


  /* Meta level type function.
   * get : List [class] -> unsigned -> class
   *
   * returns the nth element of the type list
   *
   * E.g. Get<2, List<A,B,C,D>> ==> C
   */
  template<unsigned i, class A> struct GetImpl;
  template<unsigned i, class A> using Get = typename GetImpl<i,A>::type;

  template<class A, class... As> struct GetImpl<0, List<A, As...>>
  { using type = A; };

  template<unsigned i, class A, class... As> struct GetImpl<i, List<A, As...>>
  { using type = Get<i - 1, List<As...>>; };

  /* Meta level type function.
   * zip : List [class] -> List [class] ->  List [List [class]]
   *
   * Zipps two lists of types
   *
   * E.g. Zip<List<A, B, A>, List<B, A, B>> ==> List<List<A, B>, Indexed<B, A>, Indexed<A, B>>
   */
  template<class As, class Bs> struct ZipImpl;

  template<class As, class Bs> using Zip = typename ZipImpl<As, Bs>::type;

  template<> struct ZipImpl<List<>, List<>>
  {
    using type = List<>;
  };

  template<class A, class... As, class B, class... Bs> struct ZipImpl<List<A, As...>, List<B, Bs...>>
  {
    using type = Concat<List<List<A,B>>, Zip<List<As...>, List<Bs...>>>; //Concat<List<Indexed<acc, A>>, typename WithIndicesImpl<acc + 1, List<As...>>::type>;
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// FUNTIONS RETURNING VALUES
  ////////////////////////////////////////

  /* Meta level function.
   * size : List [class] -> unsigned
   *
   * returns the number of elements in the list
   *
   * E.g. Cnt<B, List<A,B,C,B>>::val ==> 4
   */
  template<class As> struct Size;

  template<> struct Size<List<>>
  { static unsigned constexpr val = 0; };

  template<class A, class... As> struct Size<List<A, As...>>
  { static unsigned constexpr val = 1 + Size<List<As...>>::val; };


  /* Meta level function.
   * size : List [class] -> bool
   *
   * returns the number of elements in the list
   *
   * E.g. All<std::is_trivially_copyable, List<int,std::vector<int>>>::val ==> false
   *      All<std::is_trivially_copyable, List<int,bool>            >::val ==> true
   */
  template<template<class> class Pred, class As> struct All;

  template<template<class> class Pred> struct All<Pred, List<>>
  { static bool constexpr val = true; };

  template<template<class> class Pred, class A, class... As>
  struct All<Pred, List<A, As...>>
  { static bool constexpr val = Pred<A>::value && All<Pred, List<As...>>::val; };


  /* Meta level function.
   * cnt : class -> List [class] -> unsigned
   *
   * returns the number of occurences of an element in a list
   *
   * E.g. Cnt<B, List<A,B,C,B>>::val ==> 2
   */
  template<class A, class As> struct Cnt;

  template<class A> struct Cnt<A, List<>>
  { static unsigned constexpr val = 0; };

  template<class A,          class... Bs> struct Cnt<A, List<A, Bs...>>
  { static unsigned constexpr val = 1 + Cnt<A, List<Bs...>>::val; };

  template<class A, class B, class... Bs> struct Cnt<A, List<B, Bs...>>
  { static unsigned constexpr val =     Cnt<A, List<Bs...>>::val; };


  /* Meta level function.
   * contains : class -> List [class] -> bool
   *
   * returns whether some type is in the list or not.
   *
   * E.g. Contains<B, List<A,B,C,D>>::val ==> true
   */
  template<class A, class As> struct Contains;

  template<class A, class... As> struct Contains<A, List<As...>>
  { static bool constexpr val = Cnt<A,List<As...>>::val > 0; };

  /* Meta level function.
   * fistIdxOf : class -> List [class] -> unsigned
   *
   * returns the index of the first occurrence of a type
   *
   * E.g. FirstIdxOf<B, List<A,B,C,B>>::val ==> 1
   */
  template<class A, class As> struct FirstIdxOf;

  template<class A, class... As> struct FirstIdxOf<A, List<A, As...>>
  { static unsigned constexpr val = 0; };

  template<class A, class B, class... Bs> struct FirstIdxOf<A, List<B, Bs...>>
  { static unsigned constexpr val = 1 + FirstIdxOf<A, List<Bs...>>::val; };


  /* Meta level function.
   * idxOf : class -> List [class] -> unsigned
   *
   * returns the index of the unique occurrence of a type in a list
   *
   * E.g. IdxOf<B, List<A,B,C,B>>::val ==> compile fail
   *      IdxOf<C, List<A,B,C,B>>::val ==> 2
   */
  template<class A, class As> struct IdxOf;

  template<class A, class... As> struct IdxOf<A, List<As...>>
  {
    static_assert(Cnt<A, List<As...>>::val == 1, "Type must have one unique occurence in the list.");
    static unsigned constexpr val = FirstIdxOf<A, List<As...>>::val;
  };


  template<unsigned idx, class A>
  struct Indexed
  { static unsigned constexpr index = idx; };

  /*
   * Zipps the list of types terms with indices.
   *
   * E.g. WithIndices<List<A, B, A>> ==> List<Indexed<0, A>, Indexed<1, B>, Indexed<2, A>>
   */
  template<unsigned acc, class As> struct WithIndicesImpl;

  template<class As> using WithIndices = typename WithIndicesImpl<0, As>::type;


  template<unsigned acc> struct WithIndicesImpl<acc, List<>>
  {
    using type = List<>;
  };

  template<unsigned acc, class A, class... As> struct WithIndicesImpl<acc, List<A, As...>>
  {
    using type = Concat<List<Indexed<acc, A>>, typename WithIndicesImpl<acc + 1, List<As...>>::type>;
  };



  /*
   * Returns the list of indices of a list
   *
   * E.g. Indices<List<A, B, A>> ==> List<Constant<0>, Constant<1>, Constant<2>>
   */
  template<unsigned acc, class A> struct IndicesImpl;

  template<unsigned acc> struct IndicesImpl<acc, List<>>
  {
    using type = List<>;
  };

  template<unsigned acc, class A, class... As> struct IndicesImpl<acc, List<A, As...>>
  {
    using type = Concat<List<Constant<acc>>, typename IndicesImpl<acc + 1, List<As...>>::type>;
  };

  template<class As> using Indices = typename IndicesImpl<0, As>::type;

  /*
   * E.g. Range<0, 3> ==>  List<Constant<0>, Constant<1>, Constant<2>>
   */
  template<unsigned From, unsigned ToExcl> struct RangeImpl
  { using type = Concat<List<Constant<From>>, typename RangeImpl<From + 1, ToExcl>::type>; };

  template<unsigned From> struct RangeImpl<From, From>
  { using type = List<>; };

  template<unsigned From, unsigned ToExcl>
  using Range = typename RangeImpl<From,ToExcl>::type;

   template<class F>
   struct Closure {
      template<class A> using apply = std::invoke_result_t<F, A>;
   };

  /*
   * Maps a list A using the function F
   * struct Function {
   *   template<class A>
   *   using apply = std::tuple<A, A>;
   * }
   * E.g. Map<F, List<A, B, A>> ==> List<std::tuple<A, A>, std::tuple<B, B>, std::tuple<A, A>>
   */
  template<class F, class A> struct MapImpl;

  template<class F> struct MapImpl<F, List<>>
  {
    using type = List<>;
  };

  template<class F, class A, class... As> struct MapImpl<F, List<A, As...>>
  {
    using type = Concat<List<typename F::template apply<A>>, typename MapImpl<F, List<As...>>::type>;
  };

  template<class F, class As> using Map = typename MapImpl<F, As>::type;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// COMPILE TIME TESTS
  ////////////////////////////////

  namespace StaticTest {
    struct A{};
    struct B{};
    struct C{};
    struct D{};

    template<class... As>
    struct Tuple {
      Tuple(As...);
    };

    struct MkTuple {
      template<class A, class B>
      Tuple<A, B> operator()(A,B) const;
    };


#   define STATIC_TEST_TYPE_EQ(...) \
       static_assert(std::is_same< __VA_ARGS__ >::value, "TEST FAIL: types don't match fail");

#   define STATIC_TEST_VAL_EQ(...) \
       static_assert(__VA_ARGS__, "TEST FAIL: values don't match fail");

    STATIC_TEST_TYPE_EQ(
        std::invoke_result<MkTuple, A, B>::type,
        Tuple<A, B>)

    STATIC_TEST_TYPE_EQ(
        Into<Tuple, List<A, B, C> >,
        Tuple<A, B, C>)

    STATIC_TEST_TYPE_EQ(
       Concat<List<A, B>, List<A, B>>,
       List<A,B,A,B>)

    STATIC_TEST_TYPE_EQ(
       Get<2, List<A,B,C,D>>,
       C)

    STATIC_TEST_VAL_EQ(
       Size<List<A, A, B, C, A>>::val == 5)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<>>::val == 0)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<B>>::val == 1)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<B, A>>::val == 1)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<A, B>>::val == 1)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<B, B>>::val == 2)

    STATIC_TEST_VAL_EQ(
       Cnt<B, List<A,B,C,B>>::val == 2)

    STATIC_TEST_VAL_EQ(
       Contains<B, List<A,B,C,D>>::val == true)

    STATIC_TEST_VAL_EQ(
       FirstIdxOf<B, List<A,B,C,B>>::val == 1)

    STATIC_TEST_VAL_EQ(
       IdxOf<C, List<A,B,C,B>>::val == 2)

    STATIC_TEST_TYPE_EQ(
       Indices<List<A, B, A>>,
       List<Constant<0>, Constant<1>, Constant<2>>)

    template<class A> class TestFunctor { };
    struct TestFunctorFunction {
      template<class A> using apply = TestFunctor<A>;
    };

    STATIC_TEST_TYPE_EQ(
       Map<TestFunctorFunction, List<A, B, A>>,
       List<TestFunctor<A>, TestFunctor<B>, TestFunctor<A>>)

    inline void test_fun() {
      auto tuple = std::make_tuple('a', 'b', 3);
      auto tupleGet = [&](auto N) { return std::get<N.value>(tuple); };
      STATIC_TEST_TYPE_EQ(
          Map<Closure<decltype(tupleGet)>, List<Constant<0>, Constant<1>, Constant<2> >> ,
          List<char, char, int> )
    }

    inline void test_fun2() {
      auto clsr = [&](auto N) { return Constant<N.value + 1>{}; };
      STATIC_TEST_TYPE_EQ(
          Map<Closure<decltype(clsr)>, List<Constant<0>, Constant<1>, Constant<2> >> ,
          List<Constant<1>, Constant<2>, Constant<3> > )
    }

    STATIC_TEST_TYPE_EQ(
       Map<TestFunctorFunction, List<A, B, A>>,
       List<TestFunctor<A>, TestFunctor<B>, TestFunctor<A>>)

    STATIC_TEST_TYPE_EQ(
       WithIndices<List<A, B, A>>,
       List<Indexed<0, A>, Indexed<1, B>, Indexed<2, A>>)

    STATIC_TEST_TYPE_EQ(
       Zip<List<A, B>, List<C, D>>  ,
       List<List<A, C>, List<B, D>> )


    STATIC_TEST_TYPE_EQ(
       Zip<List<A, B, A>, List<B, A, B>>,
       List<List<A, B>, List<B, A>, List<A, B>>)

    STATIC_TEST_TYPE_EQ(
       Repeat<5, A> ,
       List<A, A, A, A, A> )

    STATIC_TEST_TYPE_EQ(
       Range<5, 8> ,
       List<Constant<5>, Constant<6>, Constant<7>> )

    STATIC_TEST_TYPE_EQ(
       Repeat<1, A> ,
       List<A> )

    STATIC_TEST_TYPE_EQ(
       Repeat<0, A> ,
       List<> )

    struct NotTrivCopy {
      NotTrivCopy(NotTrivCopy const&) {}
    };

    STATIC_TEST_VAL_EQ(
       std::is_trivially_copyable<NotTrivCopy>::value == false)

    STATIC_TEST_VAL_EQ(
       std::is_trivially_copyable<int>::value == true)

    STATIC_TEST_VAL_EQ(
             All<std::is_trivially_copyable, List<int,NotTrivCopy>>::val == false)

    STATIC_TEST_VAL_EQ(
             All<std::is_trivially_copyable, List<int,bool       >>::val == true)

  }

} // namespace TypeList

} // namespace Lib


#endif // __LIB_TYPELIST__H__
