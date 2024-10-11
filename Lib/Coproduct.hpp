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
 * This file defines the class Coproduct, which is a generic tagged union.
 *
 * \see UnitTests/tCoproduct.cpp for a tutorial
 */
#ifndef __LIB_COPRODUCT__H__
#define __LIB_COPRODUCT__H__

#include <type_traits>
#define MACRO_EXPANSION true

#include "Debug/Assertion.hpp"
#include "Lib/Hash.hpp"
#include "Lib/Comparison.hpp"
#include "Lib/Sort.hpp"
#include "Lib/TypeList.hpp"
#include "Lib/Option.hpp"
#include <memory>
#include <functional>
#include <tuple>
#include <vector>

namespace Lib {

namespace TL = TypeList;

template<unsigned v>
struct Constant { static constexpr unsigned value = v; };

template <class... As>
class Coproduct;

#define USE_SWITCH 0

#if USE_SWITCH

template<unsigned maxExcl> struct SwitchImpl{};

#define SWITCH_CONT_0
#define SWITCH_CONT_1  SWITCH_CONT_0  case 0:  return f(Constant<0>{});
#define SWITCH_CONT_2  SWITCH_CONT_1  case 1:  return f(Constant<1>{});
#define SWITCH_CONT_3  SWITCH_CONT_2  case 2:  return f(Constant<2>{});
#define SWITCH_CONT_4  SWITCH_CONT_3  case 3:  return f(Constant<3>{});
#define SWITCH_CONT_5  SWITCH_CONT_4  case 4:  return f(Constant<4>{});
#define SWITCH_CONT_6  SWITCH_CONT_5  case 5:  return f(Constant<5>{});
#define SWITCH_CONT_7  SWITCH_CONT_6  case 6:  return f(Constant<6>{});
#define SWITCH_CONT_8  SWITCH_CONT_7  case 7:  return f(Constant<7>{});
#define SWITCH_CONT_9  SWITCH_CONT_8  case 8:  return f(Constant<8>{});
#define SWITCH_CONT_10 SWITCH_CONT_9  case 9:  return f(Constant<9>{});
#define SWITCH_CONT_11 SWITCH_CONT_10 case 10: return f(Constant<10>{});
#define SWITCH_CONT_12 SWITCH_CONT_11 case 11: return f(Constant<11>{});
#define SWITCH_CONT_13 SWITCH_CONT_12 case 12: return f(Constant<12>{});
#define SWITCH_CONT_14 SWITCH_CONT_13 case 13: return f(Constant<13>{});
#define SWITCH_CONT_15 SWITCH_CONT_14 case 14: return f(Constant<14>{});
#define SWITCH_CONT_16 SWITCH_CONT_15 case 15: return f(Constant<15>{});

#define DECL_SWITCH_STRUCT(N)                                                             \
  template<>                                                                              \
  struct SwitchImpl<N> {                                                                  \
    template<class F>                                                                     \
    static auto apply(unsigned tag, F f) -> decltype(auto) {                              \
      if (tag < N) {                                                                      \
        switch (tag) {                                                                    \
          SWITCH_CONT_ ## N                                                               \
        }                                                                                 \
      }                                                                                   \
      ASSERTION_VIOLATION                                                                 \
    }                                                                                     \
  };                                                                                      \

DECL_SWITCH_STRUCT(1)
DECL_SWITCH_STRUCT(2)
DECL_SWITCH_STRUCT(3)
DECL_SWITCH_STRUCT(4)
DECL_SWITCH_STRUCT(5)
DECL_SWITCH_STRUCT(6)
DECL_SWITCH_STRUCT(7)
DECL_SWITCH_STRUCT(8)
DECL_SWITCH_STRUCT(9)
DECL_SWITCH_STRUCT(10)
DECL_SWITCH_STRUCT(11)
DECL_SWITCH_STRUCT(12)
DECL_SWITCH_STRUCT(13)
DECL_SWITCH_STRUCT(14)
DECL_SWITCH_STRUCT(15)
DECL_SWITCH_STRUCT(16)


template<unsigned N, class F> 
ResultOf<F, Constant<0>> switchN(unsigned tag, F fun)
{ return SwitchImpl<N>::apply(tag, std::move(fun)); }

#else // !USE_SWITCH

template<unsigned I, unsigned N>
struct SwitchImpl
{
  static_assert(I < N, "out of bounds");
  template<class F>
  inline static ResultOf<F, Constant<0>> apply(unsigned tag, F f) {
    if (tag == I) {
      return f(Constant<I>{});
    }
    return SwitchImpl<I + 1, N>::apply(tag, std::move(f));
  }
};

template<unsigned N>
struct SwitchImpl<N, N> {
  template<class F>
  inline static ResultOf<F, Constant<0>> apply(unsigned tag, F f)
  {
    ASS_EQ(tag, N)
    return f(Constant<N>{});
  }
};


template<unsigned N, class F> 
inline ResultOf<F, Constant<0>> switchN(unsigned tag, F fun)
{ return SwitchImpl<0, N - 1>::apply(tag, std::move(fun)); }

#endif  // if(USE_SWITCH) else


constexpr unsigned neededBits(unsigned i)
{ return i <= 1       ? 0
       : (i & 1) == 1 ? neededBits(i + 1)
                      : 1 + neededBits(i >> 1); }

constexpr unsigned bitMask(unsigned i)
{ return ~(unsigned(-1) << neededBits(i)); }


template<class... As>
class Coproduct;

/** This namespace constains helper classes and functions to implement the coproduct */
namespace CoproductImpl {

  template<class... As>
  class RawCoproduct;

  namespace TrivialOperations {

    template<class Op, class T> using trivial     = typename Op::template trivial<T>;
    template<class Op, class T> using DefaultImpl = typename Op::template DefaultImpl<T>;

    template<class A> struct RawCoproductTypes;

    template<class... As>
    struct RawCoproductTypes<RawCoproduct<As...>>
    { using type = TL::List<As...>; };

    template<template<class> class W, class A>
    struct RawCoproductTypes<W<A>>
    { using type = typename RawCoproductTypes<A>::type; };

    template<class Union> using Ts       = typename RawCoproductTypes<Union>::type;

    template<class Op, class ToWrap>
    using DefaultImplIfNeeded =
      std::conditional_t<TL::All<Op::template trivial , Ts<ToWrap>>::val,                 ToWrap ,
      std::conditional_t<TL::All<Op::template possible, Ts<ToWrap>>::val, DefaultImpl<Op, ToWrap>,
                                                                             ToWrap >>;

    struct Nothing {};

    template<class Op, class Ts>
    using DisableIfNeeded =
      std::conditional_t<TL::All<Op::template trivial, Ts>::val, Nothing, typename Op::Disable>;

    struct Destr {

      template<class A> using possible = std::is_destructible<A>;
      template<class A> using trivial  = std::is_trivially_destructible<A>;

      struct Disable { Disable() {}; ~Disable() {} };

      template<class T>
      struct DefaultImpl : public T {
        DefaultImpl() : T() {}
        ~DefaultImpl()
        {
          this->switchN([&](auto N) {
              using A = TL::Get<N.value, typename T::Ts>;
              this->template cast<A>().~A();
          });
        }
      };
    };

#define MK_CONS(ConsClass, REF, MOVE, move_OR_copy, OTHER_REF)                            \
    struct ConsClass {                                                                    \
                                                                                          \
      template<class A> using possible                                                    \
        = std::is_          ## move_OR_copy ## _constructible<A>;                         \
      template<class A> using trivial                                                     \
       = std::is_trivially_ ## move_OR_copy ## _constructible<A>;                         \
                                                                                          \
      struct Disable { Disable() {}; Disable(Disable REF) = delete; };                    \
                                                                                          \
      template<class T>                                                                   \
      struct DefaultImpl : public T {                                                     \
        DefaultImpl() : T() {}                                                            \
        ~DefaultImpl() = default;                                                         \
        DefaultImpl(DefaultImpl OTHER_REF other) = default;                               \
        DefaultImpl(DefaultImpl       REF other)                                          \
          : T()                                                                           \
        {                                                                                 \
          this->assignTag(other.tag());                                                   \
          this->switchN([&](auto N) {                                                     \
              using A = TL::Get<N.value, typename T::Ts>;                                 \
              ::new(&this->template cast<A>())                                            \
                A(MOVE(other.template cast<A>()));                                        \
          });                                                                             \
        }                                                                                 \
                                                                                          \
        DefaultImpl& operator=(DefaultImpl OTHER_REF other) = default;                    \
        DefaultImpl& operator=(DefaultImpl       REF other)                               \
        {                                                                                 \
          if (this == &other) return *this;                                               \
          this->switchN([&](auto N) {                                                     \
              using A = TL::Get<N.value, typename T::Ts>;                                 \
              this->template cast<A>().~A();                                              \
          });                                                                             \
          ::new(this) DefaultImpl(MOVE(other));                                           \
          return *this;                                                                   \
        }                                                                                 \
      };                                                                                  \
    };


  MK_CONS(CopyCons, const&,          , copy,     &&)
  MK_CONS(MoveCons,     &&, std::move, move, const&)


  }


  template<class... Ts>
  struct MaxSize;

  template<>
  struct MaxSize<>
  { static constexpr unsigned value = 0; };

  template<class T, class... Ts>
  struct MaxSize<T, Ts...>
  { static constexpr unsigned value = std::max<unsigned>(sizeof(T), MaxSize<Ts...>::value); };


  template<class... As>
  class RawCoproduct {

    template<class> friend struct TrivialOperations::CopyCons::DefaultImpl;
    template<class> friend struct TrivialOperations::MoveCons::DefaultImpl;
    template<class> friend struct TrivialOperations::Destr::DefaultImpl;
    template<class... Bs> friend class Lib::Coproduct;

    /** a type-level list of all types of this Coproduct */
    using Ts = TL::List<As...>;

    /** the number of alternatives */
    static constexpr unsigned size = TL::Size<Ts>::val;

    static constexpr unsigned nTags =
#if VDEBUG
                                size + 1;
#else //!VDEBUG
                                size;
#endif // VDEBUG
    static constexpr unsigned bitMask = ::bitMask(nTags);

    static_assert(nTags == 0 || nTags - 1 == ((nTags - 1) & bitMask), "bug in function neededBits");

    using Bytes = char [MaxSize<As...>::value];
    unsigned _tag: neededBits(nTags);
    Bytes _content;


    TrivialOperations::DisableIfNeeded<TrivialOperations::CopyCons, Ts> _copyCons;
    TrivialOperations::DisableIfNeeded<TrivialOperations::MoveCons, Ts> _moveCons;
    TrivialOperations::DisableIfNeeded<TrivialOperations::Destr   , Ts> _destr;

#define __COPRODUCT_CONTENT_INIT 0

#if VDEBUG
    RawCoproduct()
      : _tag(size)
    {
#if __COPRODUCT_CONTENT_INIT 
      for (unsigned i = 0; i < sizeof(Bytes); i++) {
        _content[i] = 0xFF;
      }
#endif // __COPRODUCT_CONTENT_INIT
    }
#else // !VDEBUG
    RawCoproduct() = default;
#endif // VDEBUG

    template<class F>
    ResultOf<F, Constant<0>> switchN(F f) const
    { return Lib::switchN<size>(_tag, std::move(f)); }

#define CONST_POLYMORPIHIC(CONST)                                                         \
    template<class B>                                                                     \
    B CONST& cast() CONST                                                                 \
    {                                                                                     \
      static_assert(TL::Contains<B, TL::List<As...>>::val, "invalid cast");               \
      return *(B CONST*)_content;                                                         \
    }                                                                                     \
                                                                                          \

    CONST_POLYMORPIHIC(const)
    CONST_POLYMORPIHIC(     )
#undef CONST_POLYMORPIHIC

    unsigned tag() const
    {
      ASS_REP(_tag < size, "access to uninitialized Coproduct")
      return _tag;
    }

    template<unsigned tag>
    void assignTag()
    {
      static_assert(tag < size, "tag out of bounds");
      static_assert((tag & bitMask) == tag, "unexpected lib author error");
      _tag = tag;
    }

    void assignTag(unsigned tag)
    {
      ASS_REP(tag < size, "tag out of bounds");
      ASS_REP((tag & bitMask) == tag, "unexpected lib author error");
      _tag = tag;
    }
  };

  template<class A>
  class RawCoproduct<A> {

    template<class> friend struct TrivialOperations::CopyCons::DefaultImpl;
    template<class> friend struct TrivialOperations::MoveCons::DefaultImpl;
    template<class> friend struct TrivialOperations::Destr::DefaultImpl;

    template<class... Bs> friend class Lib::Coproduct;

    /** a type-level list of all types of this Coproduct */
    using Ts = TL::List<A>;

    /** the number of alternatives */
    static constexpr unsigned size = 1;
    A _content;

    template<unsigned tag>
    void assignTag()
    { static_assert(tag == 0, "tag out of bounds"); }

    void assignTag(unsigned tag)
    { ASS_REP(tag == 0, "tag out of bounds"); }


    template<class F>
    ResultOf<F, Constant<0>> switchN(F f) const
    { return f(Constant<0>{}); }

    constexpr unsigned tag() const { return 0; }
  };


  template<class... As>
  using RawWithDefaultImpls =
     TrivialOperations::DefaultImplIfNeeded<TrivialOperations::CopyCons,
     TrivialOperations::DefaultImplIfNeeded<TrivialOperations::MoveCons,
     TrivialOperations::DefaultImplIfNeeded<TrivialOperations::Destr,
       RawCoproduct<As...>
      >>>;

  static_assert( std::is_trivially_copyable<RawCoproduct<int, int>>::value, "test 01");
  static_assert(!std::is_trivially_copyable<std::vector<int>>::value, "test 02");
  static_assert(!TL::All<std::is_trivially_copyable, TL::List<std::vector<int>, int>>::val, "test 03");
  static_assert(!std::is_trivially_copyable<RawCoproduct<std::vector<int>, int>>::value, "test 04");

  static_assert( std::is_trivially_destructible<RawCoproduct<int, int>>::value, "test 01");
  static_assert(!std::is_trivially_destructible<std::vector<int>>::value, "test 02");
  static_assert(!TL::All<std::is_trivially_destructible, TL::List<std::vector<int>, int>>::val, "test 03");
  static_assert(!std::is_trivially_destructible<RawCoproduct<std::vector<int>, int>>::value, "test 04");


} // namespace CoproductImpl



template<unsigned i, class A>
class Variant {
  A _self;
  template<class...>
  friend class Coproduct;
public:
  Variant(A a) : _self(move_if_value<A>(a)) {}
};

template<unsigned i, class A>
Variant<i, A> variant(A a)
{ return Variant<i, A>(move_if_value<A>(a)); };




/**
 * The actual Coproduct class.
 * A coproduct, also called Sum type, is a union of types, tagged with indices. It can be constructed with
 * either of the type/index pairs, and in this implementation the index can be left away if all types in this
 * coproduct are distinct.
 *
 * It is implemented as a tagged union.
 *
 * \see UnitTests/tCoproduct.cpp for usage
 */
template <class... As>
class Coproduct
{
  CoproductImpl::RawWithDefaultImpls<As...> _inner;


  /** a type-level list of all types of this Coproduct */
  using Ts = TL::List<As...>;

  /** the number of alternatives */
  static constexpr unsigned size = TL::Size<Ts>::val;

  /** unsafe default constructor, content will be uninit */
  // TODO allow uninit constructor if all alternatives are uninit constructible
  Coproduct() {}

  inline unsigned tag() const { return _inner.tag(); }

public:

  /** Returns whether this coproduct is the variant idx */
  template<unsigned idx> bool is() const
  {
    static_assert(idx < size, "out of bounds");
    return tag() == idx;
  }

  /**
   * Returns whether this coproduct is the variant witht he given type.
   * \pre is exactly one occurence of the type B in this Coproduct's types (As...).
   */
  template <class B> bool is() const
  { return is<TL::IdxOf<B, Ts>::val>(); }
                                                                                          \
  /**
   * constructs a new Coproduct with the variant idx.
   * \pre B must occur exactly once in As...
   */
  template<class B, std::enable_if_t<TL::Contains<B, Ts>::val, int> = 0>
  explicit Coproduct(B b)
    : Coproduct(Variant<TL::IdxOf<B, Ts>::val, B>(move_if_value<B>(b)))
  { }

#define REF_POLYMORPIHIC(REF, MOVE)                                                       \
                                                                                          \
  /* Coproduct &operator=(Coproduct REF other) {                                          \
    this->~Coproduct();                                                                   \
    ::new(this) Coproduct(MOVE(other));                                                   \
    return *this;                                                                         \
  }  */                                                                                   \
                                                                                          \
   /**                                                                                    \
   * transforms all variants of this Coproduct to the same type and retuns the result     \
   *                                                                                      \
   * The arguments F... must all be function whichs argument type must match the type of  \
   * the corresponding * variant of this Coproduct. The output types of the functions must\
   * all be the same type, which will be the return type of this function.                \
   */                                                                                     \
  template <class... F>                                                                   \
  inline ResultOf<TL::Get<0, TL::List<F...>>, TL::Get<0, Ts> REF> match(F... fs) REF {    \
    auto fs_ = std::tie(fs...);                                                           \
    return _inner.switchN([&](auto N) -> decltype(auto) {                                 \
        auto& f = std::get<N.value>(fs_);                                                 \
        return f(unwrap<N.value>());                                                      \
    });                                                                                   \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * transforms all variants of this Coproduct to the same type and retuns the result     \
   *                                                                                      \
   * This function works basically in the same way as match, but takes one polymorphic function object that   \
   * can transform any variant instead of multiple functions per variant.                 \
   */                                                                                     \
  template <class F>                                                                      \
  inline auto apply(F f) REF -> decltype(auto) {                                          \
    return _inner.switchN([&](auto N) -> decltype(auto) {                                 \
        return f((TL::Get<N.value, Ts> REF)MOVE(unwrap<N.value>()));                      \
    });                                                                                   \
  }                                                                                       \
                                                                                          \
  template <class F>                                                                      \
  auto applyWithIdx(F f) REF -> decltype(auto) {                                          \
    return _inner.switchN([&](auto N) -> decltype(auto) {                                 \
        return f(MOVE(unwrap<N.value>()), N);                                             \
    });                                                                                   \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * returns the value of this Coproduct if its variant is of type B. If ifs variant is of another type       \
   * the result is undefined.                                                             \
   *                                                                                      \
   * \pre B must occur exactly once in As...                                              \
   */                                                                                     \
  template <class B> inline B REF unwrap() REF                                            \
  { return MOVE(unwrap<TL::IdxOf<B, Ts>::val>()); }                                       \
                                                                                          \
  /**                                                                                     \
   * returns the value of this Coproduct if its variant's index is idx. otherwise the result is undefined.    \
   *                                                                                      \
   * \pre idx must be less than the number of variants of this Coproduct                  \
   */                                                                                     \
  template <unsigned idx>                                                                 \
  inline TL::Get<idx, Ts> REF unwrap() REF {                                              \
    static_assert(idx < size, "out of bounds");                                           \
    ASS_EQ(idx, tag());                                                                   \
    return MOVE(_inner.template cast<TL::Get<idx, Ts>>());                                \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * returns the value of this Coproduct if its variant is of type B. If ifs variant is of another type       \
   * an empty Option is returned.                                                         \
   *                                                                                      \
   * \pre B must occur exactly once in As...                                              \
   */                                                                                     \
  template <class B> inline Option<B REF> as() REF                                        \
  { return as<TL::IdxOf<B, Ts>::val>(); }                                                 \
                                                                                          \
  /**                                                                                     \
   * returns the value of this Coproduct if its variant's index is idx. otherwise an empty Option is returned.\
   *                                                                                      \
   * \pre idx must be less than the number of variants of this Coproduct                  \
   */                                                                                     \
  template <unsigned idx>                                                                 \
  inline Option<TL::Get<idx, Ts> REF> as() REF                                            \
  {                                                                                       \
    using B = TL::Get<idx, Ts>;                                                           \
    return is<idx>() ? Option<B REF>(MOVE(unwrap<idx>()))                                 \
                     : Option<B REF>();                                                   \
  }                                                                                       \

  FOR_REF_QUALIFIER(REF_POLYMORPIHIC)
#undef REF_POLYMORPIHIC

  // TODO trivial one
  friend bool operator==(const Coproduct &lhs, const Coproduct &rhs)
  {
    return lhs.tag() == rhs.tag()
      && lhs.applyWithIdx([&](auto& lhs, auto N) -> bool {
          return lhs == rhs.template unwrap<N.value>();
      });
  }

  // TODO trivial one
  friend bool operator!=(const Coproduct &lhs, const Coproduct &rhs)
  { return !(lhs == rhs); }

  template <unsigned idx, class B>
  Coproduct(Variant<idx, B> value)
  {
    static_assert(TL::Contains<B, Ts>::val, "not a variant of Coproduct");
    static_assert(idx < size, "variant index out of bounds");
    static_assert(std::is_same<B, TL::Get<idx, Ts>>::value, "illegal index for variant");

    _inner.template assignTag<idx>();
    ::new(&_inner._content) B(move_if_value<B>(value._self));
  }

  /**
   * constructs a new Coproduct with the variant idx. The argument type must match
   * `idx`th type of this * corpoduct's variants types (As...).
   */
  template <unsigned idx>
  static Coproduct variant(TL::Get<idx, Ts> value)
  { return Coproduct(Variant<idx, TL::Get<idx, Ts>>(move_if_value<TL::Get<idx, Ts>>(value))); }

  friend std::ostream &operator<<(std::ostream &out, const Coproduct &self)
  { return self.apply([&](auto const& x)  -> std::ostream&
        { return out << "var" << self.tag() << "(" << x << ")"; }); }

  friend struct std::hash<Coproduct>;

  inline Lib::Comparison compare(Coproduct const& rhs) const
  {
    auto& lhs = *this;
    return  lexCompare(DefaultComparator::compare(lhs.tag(), rhs.tag()),
        [&](){
          return lhs._inner.switchN([&](auto N){
              return DefaultComparator::compare(
                  lhs.template unwrap<N.value>(),
                  rhs.template unwrap<N.value>());
           });
        });
  }


  IMPL_COMPARISONS_FROM_COMPARE(Coproduct)

  unsigned defaultHash() const
  { return Lib::HashUtils::combine( std::hash<unsigned>{}(tag()), apply([](auto const& x){ return x.defaultHash(); })); }

  unsigned defaultHash2() const
  { return Lib::HashUtils::combine( std::hash<unsigned>{}(tag()), apply([](auto const& x){ return x.defaultHash2(); })); }

  inline Coproduct clone() const { return apply([](auto& x){ return Coproduct(x.clone()); }); }
}; // class Coproduct<As...>



} // Lib

template<class... Ts> struct std::hash<Lib::Coproduct<Ts...>>
{
  size_t operator()(Lib::Coproduct<Ts...> const& self) const
  {
    return Lib::HashUtils::combine(
        std::hash<unsigned>{}(self.tag()),
        self.apply([](auto const& x){ return std::hash<std::remove_const_t<std::remove_reference_t<decltype(x)>>>{}(x); }));
  }
};


#endif // __LIB_COPRODUCT__H__

