/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __UNIQUE_SHARED_HPP__
#define __UNIQUE_SHARED_HPP__

#include <functional>
#include "Lib/Map.hpp"
#include "Lib/Coproduct.hpp"
#include "Lib/Reflection.hpp"
#include "Lib/Sort.hpp"

namespace Lib {

#define DEBUG(...) // DBG(__VA_ARGS__)
struct PerfectPtrComparison ;
struct PerfectIdComparison ;

/** 
 * Smart pointer for perfectly sharing objects.
 *
 * This means that all objects of type T that are structurally equal, will be represented by the same pointer Perfect<T>.
 * This makes equality comparisons, hashing, and copying constant time operations.
 *
 * The type parameter DfltComparison defines how two objects of this type should be compared in operator<, operator==, 
 * and std::hash. Available options are:
 * - PerfectIdComparison : deterministic, an id is addigned to each term
 * - PerfectPtrComparison: indeterministic, pointers are compared
 *
 * T is required to be comparable with `bool operator==(const T&, const T&)`, and hashable with `std::hash<T>`.
 */
template<class T, class DfltComparison = PerfectIdComparison>
class Perfect 
{
  using IdMap = Map<const T*, Perfect, DerefPtrHash<StlHash>>;

  unsigned _id;
  const T* _ptr;
  static IdMap _ids;

  Perfect(unsigned id, const T* ptr) : _id(id), _ptr(ptr) {}

public:
  /** 
   * If an equal object to elem exists, a pointer to that object is returned.
   * Otherwise elem is moved to the heap, and a pointer to that heap location is returned.
   */
  explicit Perfect(T elem) 
    : Perfect(_ids.tryGet(&elem).toOwned()
        .unwrapOrElse([&](){
            auto entry = Perfect(_ids.size(),  new T(std::move(elem)));
            _ids.insert(entry._ptr, entry);
            DEBUG(*elemPtr, " -> ", T::className(),"#",entry._id);
            return entry;
          })) 
    { }

  /** copy constructor. Constant time. */
  Perfect(Perfect const& t) : _id(t._id), _ptr(t._ptr) {  }

  /** default constructor. for this T must be default-constructible itself. */
  Perfect() : Perfect(T()) {}

  template<class U, class C> friend bool operator==(Perfect<U, C> const& l, Perfect<U, C> const& r);

  /** dereferencing the smart pointer */
  T const* operator->() const& { return _ptr; }
  T const& operator*() const& { return *_ptr; }

  friend std::ostream& operator<<(std::ostream& out, const Perfect& self) 
  { return out << *self; }

  friend struct std::hash<Perfect<T, DfltComparison>>;

  friend struct PerfectPtrComparison;
  friend struct PerfectIdComparison;

  Lib::Comparison compare(Perfect const& rhs) const 
  { return DfltComparison::compare(*this, rhs); }
  IMPL_COMPARISONS_FROM_COMPARE(Perfect);
  IMPL_EQ_FROM_COMPARE(Perfect);

}; // class Perfect


/** instantiating the cache */
template<class T, class Cmp> typename Perfect<T, Cmp>::IdMap Perfect<T, Cmp>::_ids;

struct PerfectPtrComparison 
{
  template<class T, class Cmp>
  static Lib::Comparison compare(const Perfect<T, Cmp>& lhs, const Perfect<T, Cmp>& rhs)
  { return DefaultComparator::compare((size_t)lhs._ptr, (size_t)rhs._ptr); }

  template<class T, class Cmp>
  static bool equals(const Perfect<T, Cmp>& lhs, const Perfect<T, Cmp>& rhs) 
  { return compare(lhs, rhs) == Comparison::EQUAL; }

  template<class T, class Cmp>
  static size_t hash(Lib::Perfect<T, Cmp> const& self) 
  { return std::hash<size_t>{}((size_t)self._ptr); }
};


struct PerfectIdComparison 
{

  template<class T, class Cmp>
  static Lib::Comparison compare(const Perfect<T, Cmp>& lhs, const Perfect<T, Cmp>& rhs)
  { return DefaultComparator::compare(lhs._id, rhs._id); }

  template<class T, class Cmp>
  static bool equals(const Perfect<T, Cmp>& lhs, const Perfect<T, Cmp>& rhs) 
  { return compare(lhs, rhs) == Comparison::EQUAL; }

  template<class T, class Cmp>
  static size_t hash(Lib::Perfect<T, Cmp> const& self) 
  { return std::hash<unsigned>{}(self._id); }
};


/** function to create a Perfect<T> ergonomically (with the help of type deduction) */
template<class T, class Cmp = PerfectIdComparison> 
Perfect<T, Cmp> perfect(T t) 
{ return Perfect<T, Cmp>(std::move(t)); } } // namespace Lib

template<class A, class B, class Cmp> 
auto operator*(A const& l, Perfect<B, Cmp> const& r) 
{ return perfect(l * (*r)); }

template<class A, class B> 
auto operator*(Perfect<A> const& l, B const& r) 
{ return perfect((*l) * r); }

template<class A, class B> 
auto operator*(Perfect<A> const& l, Perfect<B> const& r) 
{ return perfect((*l) * (*r)); }

template<class A> 
auto operator-(Perfect<A> const& x) 
{ return perfect(-(*x)); }


template<class T, class Cmp> struct std::hash<Lib::Perfect<T, Cmp>> 
{
  size_t operator()(Lib::Perfect<T, Cmp> const& self) const 
  { return Cmp::hash(self); }
};


template<class T, class Cmp> struct std::less<Lib::Perfect<T, Cmp>> 
{
  bool operator()(Lib::Perfect<T, Cmp> const& lhs, Lib::Perfect<T, Cmp> const& rhs) const 
  { return Cmp{}(lhs, rhs); }
};


#undef DEBUG
#endif // __UNIQUE_SHARED_HPP__
