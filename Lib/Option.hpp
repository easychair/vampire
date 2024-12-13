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
 * This file mainly defined the class Option, which can be thought of as a NULLable pointer, that is 
 * stack-allocated, with RAII semantics. 
 *
 * \see UnitTests/tOption.cpp for examples of the usage
 */
#ifndef __OPTIONAL_H__
#define __OPTIONAL_H__

#include <type_traits>
#include "Debug/Assertion.hpp"
#include <iostream>
#include "Lib/Reflection.hpp"


namespace Lib {

template<
  class T,
  typename std::enable_if<std::is_reference<T>::value, bool>::type = true
  > 
T move_if_value(std::remove_reference_t<T> && t) 
{ return t; }


template<
  class T,
  typename std::enable_if<std::is_reference<T>::value, bool>::type = true
  > 
T move_if_value(std::remove_reference_t<T>& t) 
{ return t; }


template<
  class T,
  typename std::enable_if< !std::is_reference<T>::value , bool>::type = true
  > 
T move_if_value(T& t) 
{ return std::move(t); }

template<
  class T,
  typename std::enable_if< !std::is_reference<T>::value , bool>::type = true
  > 
T move_if_value(T const& t) 
{ return std::move(t); }


template<
  class T,
  typename std::enable_if< !std::is_reference<T>::value , bool>::type = true
  > 
T move_if_value(T && t) 
{ return std::move(t); }

#define FOR_REF_QUALIFIER(macro)                                                          \
  macro(const &, ) macro(&, ) macro(&&, std::move)

template<class T>
struct MaybeUninit {
  union Value { 
    T init; int uninint[0]; 
    Value() : uninint{} {};
    ~Value() {};
  } _elem;

   MaybeUninit() : _elem() {}
  ~MaybeUninit() {}
// This macro will b expanded for (REF,MV) in { (`&&`, `std::move`), (`&`, ``), (`const &`, ``) }
#define methods(REF, MV)                                                                  \
  operator T REF() REF                                                                    \
  { return MV(_elem.init); }                                                              \
                                                                                          \
  void init(T REF content)                                                                \
  { ::new(&_elem)T(MV(content)); }                                                        \
                                                                                          \
  MaybeUninit& operator=(T REF content)                                                   \
  {                                                                                       \
    _elem.init = MV(content);                                                             \
    return *this;                                                                         \
  }                                                                                       \

  FOR_REF_QUALIFIER(methods)

#undef methods 
};

template<class A>
class OptionBase 
{

  bool _isSome;
  MaybeUninit<A> _elem;
public:

  constexpr OptionBase() : _isSome(false) {}

  ~OptionBase() 
  { 
    if (isSome()) { 
      unwrap().~A(); 
    }
  }

#define methods(REF, MV)                                                                  \
  constexpr explicit OptionBase(A REF content)                                            \
    : _isSome(true)                                                                       \
      , _elem()                                                                           \
  {                                                                                       \
    _elem.init(move_if_value<A>(content));                                                \
  }                                                                                       \
                                                                                          \
  A REF unwrap() REF                                                                      \
  {                                                                                       \
    ASS(_isSome);                                                                         \
    return MV(_elem);                                                                     \
  }                                                                                       \
                                                                                          \
  OptionBase(OptionBase REF a) : _isSome(a._isSome)                                       \
  {                                                                                       \
    if (isSome()) {                                                                       \
      _elem.init(MV(a).unwrap());                                                         \
    }                                                                                     \
  }                                                                                       \

  FOR_REF_QUALIFIER(methods)

#undef methods

  OptionBase& operator=(OptionBase&& other)
  {
    if (_isSome) {
      unwrap().~A();
    }
    if (other._isSome) {
      _elem.init(move_if_value<A>(other.unwrap()));
    }
    _isSome = other._isSome;
    return *this;
  }

  OptionBase& operator=(OptionBase const& other)
  {
    if (_isSome) {
      unwrap().~A();
    }
    if (other._isSome) {
      _elem.init(other.unwrap());
    }
    _isSome = other._isSome;
    return *this;
  }


  bool isSome() const { return _isSome;   }
  bool isNone() const { return !isSome(); }

  static OptionBase fromPtr(A* ptr) 
  { return ptr == nullptr ? OptionBase() : OptionBase(*ptr); }

  friend bool operator==(OptionBase const& lhs, OptionBase const& rhs) 
  { 
    if (lhs._isSome != rhs._isSome) return false;
    
    if (lhs._isSome) {
      return lhs.unwrap() == rhs.unwrap();
    } else {
      return true;
    }
  }
    
};



template<class A>
class OptionBaseRef
{

  A* _elem;
public:

  constexpr OptionBaseRef(                ) : _elem(nullptr) {  }
  constexpr OptionBaseRef(A const* content) : _elem(const_cast<A*>(content)) { }
  constexpr OptionBaseRef(A      * content) : _elem(content) { }

  bool isSome() const { return _elem != nullptr;   }

  A const& unwrap() const& { ASS(isSome()); return           *_elem ; }
  A     && unwrap()     && { ASS(isSome()); return std::move(*_elem); }
  A      & unwrap()      & { ASS(isSome()); return           *_elem ; }

  constexpr OptionBaseRef(OptionBaseRef      & a) = default;
  constexpr OptionBaseRef(OptionBaseRef     && a) = default;
  constexpr OptionBaseRef(OptionBaseRef const& a) = default;

  OptionBaseRef& operator=(OptionBaseRef      & a) = default;
  OptionBaseRef& operator=(OptionBaseRef     && a) = default;
  OptionBaseRef& operator=(OptionBaseRef const& a) = default;

  static OptionBaseRef fromPtr(A* ptr) 
  { return ptr == nullptr ? OptionBaseRef() : *ptr; }

  friend bool operator==(OptionBaseRef const& lhs, OptionBaseRef const& rhs) 
  { return (lhs._elem == nullptr && rhs._elem == nullptr)
        || (lhs._elem != nullptr && rhs._elem != nullptr && (lhs._elem == rhs._elem || *lhs._elem == *rhs._elem )); }

};

template<class A>
class OptionBase<A const&> : public OptionBaseRef<A>
{
public:
  constexpr OptionBase() : OptionBaseRef<A>() {}
  constexpr OptionBase(A const& item) : OptionBaseRef<A>(&item) {}
  constexpr OptionBase(OptionBase const& b) : OptionBaseRef<A>(b) {}
};

template<class A>
class OptionBase<A&> : public OptionBaseRef<A>
{
public:
  constexpr OptionBase() : OptionBaseRef<A>() {}
  constexpr OptionBase(A& item) : OptionBaseRef<A>(&item) {}
  constexpr OptionBase(OptionBase const& b) : OptionBaseRef<A>(b) {}
};


template<class A>
class OptionBase<A&&> : public OptionBaseRef<A>
{
public:
  constexpr OptionBase() : OptionBaseRef<A>() {}
  constexpr OptionBase(A&& item) : OptionBaseRef<A>(&item) {}
  constexpr OptionBase(OptionBase const& b) : OptionBaseRef<A>(b) {}
};

/** The actual Option class
 * An Option<A> is a class that holds either a value of type A, or is none/empty.
 * It can be thought of a nullable pointer, that has the advantage that does not need to be allocated 
 * in a separate structure, and does not expose any uninitialized memory to the user. Further it 
 * automatically calls the destructor when it goes out of scope.
 *
 * \see UnitTests/tOption.cpp for usage examples
 */
template<class A>
class Option : OptionBase<A> {

  explicit constexpr Option(OptionBase<A>&& base) : OptionBase<A>(std::move(base)) {  }
public:
  using Content = A;

  /** constructs an option from a value of type A&, A const&, or A&&. */
  using OptionBase<A>::OptionBase;

  /** checks whether the Option holds a value */
  using OptionBase<A>::isSome;

  /** returns the Options value if it holds one */
  using OptionBase<A>::unwrap;

  friend bool operator==(Option const& lhs, Option const& rhs) 
  { return static_cast<OptionBase<A>const&>(lhs) == static_cast<OptionBase<A>const&>(rhs); }

  friend bool operator!=(Option const& lhs, Option const& rhs) 
  { return !(lhs == rhs); }

  /** creates an Option<A&>, or Option<A const&> from a pointer A*. if the pointer is NULL the option will be empty */
  template<class C> static Option<A> fromPtr(C self) 
  { return Option(OptionBase<A>::fromPtr(self)); }

  /** checks whether the option is empty */
  bool isNone() const { return !this->isSome(); }

  operator bool() const { return isSome(); }
  template<class T, std::enable_if_t<!std::is_same_v<T, bool>, bool> = true>
  operator T() const = delete;

  A const& operator*() const { return unwrap(); }
  A      & operator*()       { return unwrap(); }

  std::remove_reference_t<A> const* operator->() const { return &unwrap(); }
  std::remove_reference_t<A>      * operator->()       { return &unwrap(); }

  std::remove_reference_t<A>      * asPtr()       { return isSome() ? &unwrap() : nullptr; }
  std::remove_reference_t<A> const* asPtr() const { return isSome() ? &unwrap() : nullptr; }


  Option take() 
  {
    Option out;
    std::swap(*this,out);
    return out;
  }
  /** 
   * returns the value held by this option if there is one, or calls the given function f without arguments, 
   * initializes the closuer with the returned value, and returns a reference to the value afterwards.
   */ 
  template<class Clsr>
  A& unwrapOrInit(Clsr f) { 
    if (isNone()) {
      ::new(this) Option(f());
    }
    return this->unwrap();
  }

#define ref_polymorphic(REF, MOVE)                                                        \
                                                                                          \
  /**                                                                                     \
   * applies the given function to the value of this option and returns an option of the return type. \
   * if the Option was None an empty option of the function's return type is returned.    \
   */                                                                                     \
  template<class Clsr>                                                                    \
  Option<typename std::invoke_result<Clsr, A REF>::type> map(Clsr clsr) REF {                 \
    using OptOut = Option<typename std::invoke_result<Clsr, A REF>::type>;                    \
    return this->isSome() ? OptOut(clsr(MOVE(unwrap())))                                  \
                          : OptOut();                                                     \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * if the Option holds a value the first function is applied to the value.              \
   * if the Option is none the second function is called without arguments and the result is returned.\
   * \pre both CaseSome and CaseNone must have the same return type                       \
   */                                                                                     \
  template<class CaseSome, class CaseNone>                                                \
  typename std::invoke_result<CaseSome, A REF>::type match(CaseSome present, CaseNone none) REF {      \
    if (this->isSome()) {                                                                 \
      return present(MOVE((*this)).unwrap());                                             \
    } else {                                                                              \
      return none();                                                                      \
    }                                                                                     \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * returns the value held by this option if there is one, or returns the value alt otherwise      \
   */                                                                                     \
  A REF unwrapOr(A REF alt) REF {                                                         \
    if (this->isSome()) {                                                                 \
      return MOVE(*this).unwrap();                                                        \
    } else {                                                                              \
      return MOVE(alt);                                                                   \
    }                                                                                     \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * returns the value held by this option if there is one, or calls the given function f without arguments   \
   * and returns the value otherwise.                                                     \
   */                                                                                     \
  template<class Clsr>                                                                    \
  A unwrapOrElse(Clsr f) REF {                                                            \
    if (this->isSome()) {                                                                 \
      return MOVE(*this).unwrap();                                                        \
    } else {                                                                              \
      return f();                                                                         \
    }                                                                                     \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * Returns this, if this is Some, or uses the closure to create an alternative option if this is None.      \
   */                                                                                     \
  template<class Clsr,                                                                    \
           typename std::enable_if<std::is_same< typename std::invoke_result<Clsr>::type    \
                                      , Option                                            \
                                      >::value                                            \
                         , bool                                                           \
                         >::type = true                                                   \
          >                                                                               \
  auto orElse(Clsr clsr) REF -> Option                                                    \
  { return this->isSome() ? MOVE(*this) : clsr(); }                                       \
                                                                                          \
  /** Returns the value of this, if this is Some, or uses the closure to create a value othewise. */\
  template<class Clsr,                                                                    \
           typename std::enable_if<std::is_same< typename std::invoke_result<Clsr>::type    \
                                      , A                                                 \
                                      >::value                                            \
                         , bool                                                           \
                         >::type = true                                                   \
          >                                                                               \
  auto orElse(Clsr clsr) REF -> A                                                         \
  { return this->isSome() ? MOVE(*this).unwrap() : clsr(); }                              \
                                                                                          \
   /**                                                                                    \
   * applies a function to the value of this closure if ther is one. the function is expected to return\
   * another option. the resulting Option<Option<Result>> will then be flattened to an Option<Result>.\
   *                                                                                      \
   * This function is the same as flatMap/andThen/(>>=)  in other programming languages with monads.\
   */                                                                                     \
  template<class Clsr>                                                                    \
  typename std::invoke_result<Clsr, A REF>::type andThen(Clsr clsr) REF {                     \
    using OptOut = typename std::invoke_result<Clsr, A REF>::type;                            \
    return this->isSome() ? clsr(MOVE(*this).unwrap())                                    \
                          : OptOut();                                                     \
  }                                                                                       \
                                                                                          \
  template<class Clsr> auto flatMap(Clsr clsr) REF { return andThen(clsr); }              \
                                                                                          \
  template<class Pred>                                                                    \
  Option filter(Pred p) REF {                                                             \
    return isSome() && p(unwrap())                                                        \
                ? MOVE(*this)                                                             \
                : Option();                                                               \
  }                                                                                       \
                                                                                          \
  /**                                                                                     \
   * turns an Option<A&>, Option<A const&>, or Option<A&&> into an Option<A> by calling the appropriate move  \
   * or copy constructor.                                                                 \
   */                                                                                     \
  Option<typename std::remove_const<typename std::remove_reference<A>::type>::type>  toOwned() REF  \
  {                                                                                       \
    using Out = typename std::remove_const<typename std::remove_reference<A>::type>::type;\
    return map([](A REF  elem) -> Out { return Out(MOVE(elem)); });                       \
  }                                                                                       \

  FOR_REF_QUALIFIER(ref_polymorphic)

#undef ref_polymorphic

  class OptionIter {
    Option _self;
    OptionIter(Option self) : _self(std::move(self)) {}
  public:
    friend class Option;
    DECL_ELEMENT_TYPE(A);
    inline bool hasNext() const { return _self.isSome(); }
    inline bool hasNext()       { return _self.isSome(); }
    inline A next() { return _self.take().unwrap(); }
  };
  
  Option<A const&> asRef() const { return someIf(isSome(), [&]() -> decltype(auto) { return unwrap(); }); }
  Option<A      &> asRef()       { return someIf(isSome(), [&]() -> decltype(auto) { return unwrap(); }); }

  OptionIter intoIter() &&
  { return OptionIter(std::move(*this)); }

  auto iter() const& { return asRef().intoIter(); }
  auto iter()      & { return asRef().intoIter(); }

  friend std::ostream& operator<<(std::ostream& out, Option const& self) 
  { return self.isSome() ?  out << self.unwrap() : out << "None"; }


  auto flatten()      & { return andThen([](auto x) { return  std::move(x); }); }
  auto flatten()     && { return std::move(*this).andThen([](auto x) { return  std::move(x); }); }
  auto flatten() const& { return andThen([](auto x) { return  std::move(x); }); }

  friend bool operator<(Option const& lhs, Option const& rhs) 
  { 
    if (lhs.isSome() < rhs.isSome()) {
      return true;
    } else if (lhs.isSome() > rhs.isSome()) {
      return false;
    } else if (lhs.isNone()) {
      ASS(rhs.isNone()) 
      return false;
    } else {
      ASS(rhs.isSome()) 
      ASS(lhs.isSome()) 
      return lhs.unwrap() < rhs.unwrap();
    }
  }
};


template<class F> 
auto someIf(bool condition, F create) -> Option<decltype(create())> 
{ return condition ? Option<decltype(create())>(create()) 
                   : Option<decltype(create())>(); }


template<class T> constexpr Option<T> some(T const& t) { return Option<T>(t);            }
template<class T> constexpr Option<T> some(T      & t) { return Option<T>(t);            }
template<class T> constexpr Option<T> some(T     && t) { return Option<T>(std::move(t)); }

template<class T> constexpr Option<T> none() { return Option<T>(); }

template<class T> constexpr Option<T> optionalFromPtr(T* t) { return Option<T>::fromPtr(t); }
template<class T> constexpr Option<T> optionFromPtr(T* t) { return Option<T>::fromPtr(t); }

template<class T>
T operator||(Option<T> t, T c)
{ return std::move(t).unwrapOr(std::move(c)); }

template<class T, class Clsr>
auto operator||(Option<T> t, Clsr f) -> decltype(f())
{ return std::move(t).orElse(f); }

template<class T>
Option<T> operator||(Option<T> t, Option<T> c)
{ return std::move(t).orElse([&](){ return std::move(c); }); }


template<class T, class Clsr>
Option<T> operator&&(Option<T> t, Clsr c)
{ return std::move(t).andThen(c); }

template<class T>
Option<T> operator&&(Option<T> t, Option<T> c)
{ return std::move(t).andThen([&](){ return std::move(c); }); }

} // namespace Lib

#endif // __OPTIONAL_H__
