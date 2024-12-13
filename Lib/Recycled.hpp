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
 * @file Recycled.hpp
 */

#ifndef __Recycled__
#define __Recycled__

#include "Forwards.hpp"

#include "Stack.hpp"
#include "DArray.hpp"
#include <initializer_list>
#include <memory>
#include "Lib/Reflection.hpp"
#include "Lib/Hash.hpp"

namespace Lib
{

struct DefaultReset
{ 
  template<class T> void operator()(T& t) { t.reset(); } 

  template<unsigned i, unsigned sz, class Tup> 
  struct __ResetTuple
  {
    static void apply(Tup& self)
    {
      std::get<i>(self).reset();
      __ResetTuple<i + 1, sz, Tup>::apply(self);
    }
  };

  template<unsigned i, class Tup> 
  struct __ResetTuple<i, i, Tup>  {
    static void apply(Tup& self)
    { std::get<i>(self).reset(); }
  };

  void operator()(std::tuple<>& t) {}

  template<class A, class... As>
  void operator()(std::tuple<A, As...>& t) 
  { __ResetTuple<0, std::tuple_size<std::tuple<A, As...>>::value - 1, std::tuple<A, As...>>::apply(t); }

};

struct DefaultKeepRecycled
{ 
  template<class T> bool operator()(T const& t) { return t.keepRecycled(); } 

  template<unsigned i, unsigned sz, class Tup> 
  struct __KeepRecycledTuple
  {
    static bool apply(Tup const& self)
    { return std::get<i>(self).keepRecycled() 
        || __KeepRecycledTuple<i + 1, sz, Tup>::apply(self); }
  };

  template<unsigned i, class Tup> 
  struct __KeepRecycledTuple<i, i, Tup>  {
    static bool apply(Tup const& self)
    { return std::get<i>(self).keepRecycled(); }
  };

  bool operator()(std::tuple<> const& t) { return false; }

  template<class A, class... As>
  bool operator()(std::tuple<A, As...> const& t) 
  { return __KeepRecycledTuple<0, std::tuple_size<std::tuple<A, As...>>::value - 1, std::tuple<A, As...>>::apply(t); }
};

struct NoReset
{ template<class T> void operator()(T& t) {  } };

template<class T>
struct MaybeAlive
{
  T _self;
  bool* _alive;
  MaybeAlive(T self, bool* alive)
    : _self(std::move(self)), _alive(alive)
  { *_alive = true; }
  ~MaybeAlive() { *_alive = false; }


  T const& operator* () const { return  _self; }
  T const* operator->() const { return &_self; }
  T& operator* () { return  _self; }
  T* operator->() { return &_self; }
};

#define USE_PTRS 1


#if USE_PTRS
#define IF_USE_PTRS(l, r) l
#else
#define IF_USE_PTRS(l, r) r
#endif

/** A smart that lets you keep memory allocated, and reuse it.
 * Constructs an object of type T on the heap. When the Recycled<T> is destroyed,
 * the object is not returned, but the object is reset using Reset::operator(),
 * and returned to a pool of pre-allocated objects. When the next Recycled<T> is
 * constructed an object from the pool is returned instead of allocating heap memory.
 *
 * This may be useful if profiling suggests it, but in general you should not need this.
 * If you think you do, there are the following tradeoffs:
 * - more construction overhead than both automatic and (slightly) static storage
 * - extra indirection: Recycled<T> internally holds a pointer to T
 * - shares allocations of type T between all instances of Recycled<T>
 * - safer than static storage due to the Reset mechanism
 * - Recycled<T> is reentrant, whereas static is not
 * - "leaks" memory, retaining the peak number of concurrently-live T objects
 * - implements move semantics where T may not
 * - "destruction" often very cheap because T is never destroyed
 *
 * To summarise:
 * 1. Use automatic storage unless you really have to do something else.
 * 2. If you must, use Recycled<T> by preference to static for safety reasons.
 * 3. If you really must and you are sure of correctness (reentrancy, reset logic),
 *    use static storage.
 */
template<class T, class Reset = DefaultReset, class Keep = DefaultKeepRecycled>
class Recycled
{
  using Self = IF_USE_PTRS(std::unique_ptr<T>,T);
  Self _self;
  Reset _reset;
  Keep _keep;

  static bool memAlive;
  static Stack<Self>& mem() {
    static MaybeAlive<Stack<Self>> mem(Stack<Self>(), &memAlive);
    return *mem;
  }
  Recycled(Recycled const&) = delete;
  
  T      & self()       { ASS(alive()) return IF_USE_PTRS(*, )_self; }
  T const& self() const { ASS(alive()) return IF_USE_PTRS(*, )_self; }
  struct EmptyConstructMarker {};

  // Recycled(Self self)  : _self(std::move(self)), _reset(), _keep() {}
public:

  Recycled()
    : _self(mem().isNonEmpty() ? mem().pop() : IF_USE_PTRS(std::make_unique<T>(), T()))
  { }

  // template<class A, class... As>
  // Recycled(A a, As... as) : Recycled()
  // { 
  //   self().~T();
  //   new(&self) T(a, as...); 
  // }


  template<class A, class... As>
  Recycled(A a, As... as) 
    : _self([&](){ 
        if (mem().isNonEmpty()) {
          auto elem = mem().pop();
          elem->init(a, as...);
          return elem;
        } else {
          return IF_USE_PTRS(std::make_unique<T>, T)(a, as...);
        }
    }())
  { }

  template<class Clone>
  Recycled clone(Clone cloneFn) const
  {
    Recycled c;
    T const& from = **this;
    T& to = *c;
    cloneFn(to, from);
    return c;
  }

  bool alive() const { return IF_USE_PTRS(bool(_self), true); }

  auto asTuple() const -> decltype(auto) { return std::make_tuple(someIf(alive(), [this]() -> decltype(auto) { return self(); })); }
  IMPL_COMPARISONS_FROM_TUPLE(Recycled);
  IMPL_HASH_FROM_TUPLE(Recycled);

  Recycled(Recycled&& other) = default;
  Recycled& operator=(Recycled&& other) = default;

  ~Recycled()
  {
    if (IF_USE_PTRS(_self, true) && _keep(self()) && memAlive) {
      _reset(self());
      mem().push(std::move(_self));
    }
  }

  T const& operator* () const { return  self(); }
  T const* operator->() const { return &self(); }
  T& operator* () { return  self(); }
  T* operator->() { return &self(); }

  auto size() const { return self().size(); }
  template<class Idx> decltype(auto) operator[](Idx idx) const { return self()[idx]; }
  template<class Idx> decltype(auto) operator[](Idx idx)       { return self()[idx]; }

  friend std::ostream& operator<<(std::ostream& out, Recycled const& self)
  { if (self.alive())return out << self.self(); else return out << "Recycled(NULL)"; }
};

template<class T, class Reset, class Keep>
bool Recycled<T, Reset, Keep>::memAlive = true;

};

template<class T>
using RStack = Recycled<Stack<T>>;

template<class T>
Recycled<Stack<T>> recycledStack() 
{ return Recycled<Stack<T>>(); }


template<class T, class... Ts>
Recycled<Stack<T>> recycledStack(T t, Ts... ts) 
{
  Recycled<Stack<T>> out;
  out->pushMany(std::move(t), std::move(ts)...);
  return out;
}

#endif /*__Recycled__*/
