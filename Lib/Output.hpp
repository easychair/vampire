/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#ifndef __Lib_Output_HPP__
#define __Lib_Output_HPP__

#include <iostream>
#include <utility>
#include <tuple>

namespace std {
template<class... As> 
std::ostream& operator<<(std::ostream& out, tuple<As...> const& self);
} // namespace std

namespace Lib {
namespace Output {

/** Newtype in order to nicely output a pointer.
 * Usage: `out << Output::ptr(ptr) << std::endl;` 
 */
template<class T>
struct Ptr { T self; };

template<class T, std::enable_if_t<!std::is_pointer_v<T>, bool> = true>
Ptr<T const&> ptr(T const& self) { return { .self = self, }; }

template<class T, std::enable_if_t<std::is_pointer_v<T>, bool> = true>
Ptr<T> ptr(T self) { return { .self = self, }; }

template<class T>
std::ostream& operator<<(std::ostream& out, Ptr<T*> const& self)
{ return self.self == nullptr ? out << "NULL" : out << *self.self; }

template<class T>
std::ostream& operator<<(std::ostream& out, Ptr<T> const& self)
{ return out << self.self; }

template<class T>
struct Repeat { T const& to_repeat; unsigned times; };

template<class T>
Repeat<T> repeat(T const& c, unsigned times)
{ return Repeat<T>{c, times}; }

template<class T>
std::ostream& operator<<(std::ostream& out, Repeat<T> const& self)
{ for (unsigned i = 0; i < self.times; i++) out << self.to_repeat; return out; };

/** Newtype for outputting a datatype that implements it in multiline format.
 * Usage: `out << multiline(substitutioTree) << std::endl;` 
 *
 * You can implement this for a `MyType` by implementing 
 * std::ostream& operator<<(std::ostream&, Output::Multiline<MyType>)
 */
template<class T>
struct Multiline { 
  T const& self; 
  unsigned indent; 

  static void outputIndent(std::ostream& out, unsigned indent);
};

template<class T>
Multiline<T> multiline(T const& self, unsigned indent = 0)
{ return { self, indent, }; }

template<class T>
void Multiline<T>::outputIndent(std::ostream& out, unsigned indent)
{ out << Output::repeat("    ", indent); };



template<class Sep, class Iter>
struct Interleaved { Sep const& sep; Iter iter; };

template<class Sep, class Iter>
struct Interleaved<Sep,Iter> interleaved(Sep const& s, Iter i)
{ return Interleaved<Sep, Iter>{s, std::move(i)}; }

template<class Sep, class Iter>
std::ostream& operator<<(std::ostream& out, Interleaved<Sep, Iter> self)
{
  if (self.iter.hasNext()) {
    out << self.iter.next();
    while (self.iter.hasNext()) {
      out << self.sep << self.iter.next();
    }
  }
  return out;
}


template<class Iter>
auto commaSep(Iter i) { return Output::interleaved(", ", std::move(i)); }

struct Nothing {
  friend std::ostream& operator<<(std::ostream& out, Nothing const& self)
  { return out << ""; }
};

struct Comma {
  friend std::ostream& operator<<(std::ostream& out, Comma const& self)
  { return out << ", "; }
};


template<unsigned i, unsigned sz, class Tup, class Sep> 
struct __TupleHelper
{
  static void apply(std::ostream& out, Tup const& self)
  {
   out << Sep{} << std::get<i>(self);
    __TupleHelper<i + 1, sz, Tup, Sep>::apply(out, self);
  }
};

template<class Sep> 
struct __TupleHelper<0, 0, std::tuple<>, Sep>  
{ static void apply(std::ostream& out, std::tuple<> const& self) { } };

template<unsigned sz, class Tup, class Sep> 
struct __TupleHelper<sz, sz, Tup, Sep>  
{ static void apply(std::ostream& out, Tup const& self) { } };


template<unsigned sz, class Tup, class Sep> 
struct __TupleHelper<0, sz, Tup, Sep>  
{
  static void apply(std::ostream& out, Tup const& self)
  {
    out << std::get<0>(self);
    __TupleHelper<1, sz, Tup, Sep>::apply(out, self);
  }
};


template<class... Elems>
struct OutputCat { std::tuple<Elems...> elems; };

template<class... Elems>
struct OutputCat<Elems...> catOwned(Elems&&... elems)
{ return OutputCat<Elems...>{std::tuple<Elems...>(std::move(elems)...)}; }

template<class... Elems>
struct OutputCat<Elems const&...> cat(Elems const&... elems)
{ return OutputCat<Elems const&...>{std::tie(elems...)}; }

template<class... As> 
std::ostream& operator<<(std::ostream& out, OutputCat<As...> const& self)
{ 
  __TupleHelper<0, std::tuple_size<std::tuple<As...>>::value, std::tuple<As...>, Nothing>::apply(out, self.elems);
  return out;
}

} // namespace Output
} // namespace Lib

namespace std {

template<class A, class B>
std::ostream& operator<<(std::ostream& out, pair<A,B> const& self)
{ return out << "(" << self.first << ", " << self.second << ")"; }


template<class... As> 
std::ostream& operator<<(std::ostream& out, tuple<As...> const& self)
{ 
  out << "(";
  Lib::Output::__TupleHelper<0, std::tuple_size<std::tuple<As...>>::value, std::tuple<As...>, Lib::Output::Comma>::apply(out, self);
  out << ")";
  return out;
}

} // namespace std

#endif // __Lib_Output_HPP__
