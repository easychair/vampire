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
 * @file DArray.hpp
 * Defines a class of self-deallocating arrays. They should be used instead
 * of Array when the size is known in advance
 *
 * @since 30/12/2007 Manchester
 */

#ifndef __DArray__
#define __DArray__

#include "Forwards.hpp"
#include "Debug/Assertion.hpp"

#include "Allocator.hpp"
#include "Comparison.hpp"
#include "Exception.hpp"
#include "Random.hpp"
#include "Reflection.hpp"
#include "VirtualIterator.hpp"

namespace Lib {

/**
 * Class of fixed-size self-deallocating generic arrays.
 * @param C the type of content
 * @since 30/12/2007 Manchester
 */
template<typename C>
class DArray
{
public:
  USE_ALLOCATOR(DArray<C>);

  class Iterator;

  DECL_ELEMENT_TYPE(C);
  DECL_ITERATOR_TYPE(Iterator);
  /**
   * Create an array having the given @b size
   * @since 30/12/2007 Manchester
   */
  inline
  DArray (size_t size=0)
    : _size(size), _capacity(size)
  {
    if(size>0) {
      void* mem = ALLOC_KNOWN(sizeof(C)*_capacity,"DArray<>");
      _array = array_new<C>(mem, _capacity);
    } else {
      _array=0;
    }
  }

  explicit DArray(const DArray& o)
    : _size(o.size()), _capacity(o.size())
  {
    if(_size==0) {
      _array=0;
      return;
    }
    void* mem = ALLOC_KNOWN(sizeof(C)*_capacity,"DArray<>");
    _array = static_cast<C*>(mem);
    for(size_t i=0; i<_size; i++) {
      ::new (&_array[i]) C(o[i]);
    }
  }

  DArray clone() const { return DArray(*this); }

  void swap(DArray& other) {
    std::swap(other._size, _size);
    std::swap(other._capacity, _capacity);
    std::swap(other._array, _array);
  }

  bool keepRecycled() const { return _capacity > 0; }

  DArray(DArray&& other) : DArray() { swap(other); }
  DArray& operator=(DArray&& other) { swap(other); return *this; }

  /** Delete array */
  inline ~DArray()
  {
    if(_array) {
      array_delete(_array, _capacity);
      DEALLOC_KNOWN(_array,sizeof(C)*_capacity,"DArray<>");
    }
  }

  /** Return a reference to the n-th element of the array */
  inline C& operator[] (size_t n)
  {
    ASS_L(n,_size);
    return _array[n];
  } // operator[]

  /** Return a reference to the n-th element of the array */
  inline const C& operator[](size_t n) const
  {
    ASS_L(n,_size);
    return _array[n];
  }
 /**
   * Set array's size to @b s. @b s must be smaller or equal to the
   * current array size.
   */
  inline void shrink(size_t s)
  {
    ASS_LE(s,_size);

    _size = s;
  }

  inline bool operator==(const DArray& o) const
  {
    if(size()!=o.size()) { return false; }
    size_t sz = size();
    for(size_t i=0; i<sz; i++) {
      if(!((*this)[i]==o[i])) { return false; }
    }
    return true;
  }

  inline bool operator!=(const DArray& o) const { return !(*this==o); }

  /** return the standard array represented by this DArray */
  inline C* array() { return _array; }
  inline const C* array() const { return _array; }

  inline C* begin() { return _array; }
  inline C* end() { return _array+_size; }

  /**
   * Set array's size to @b s and that its capacity is at least @b s.
   * Return true iff array was not extended.
   * If the capacity is smaller, the array will be extended.
   *
   * @warning upon extension no copying of elements will be done
   *   so the operation is safe only before processing the array.
   * @since 02/01/2008 Manchester
   */
  inline bool ensure(size_t s)
  {
    if (_capacity >= s) {
      _size = s;
      return true;
    }

    size_t newCapacity = std::max(s, _capacity*2);

    void* mem = ALLOC_KNOWN(sizeof(C)*newCapacity,"DArray<>");
    C* newArray=array_new<C>(mem, newCapacity);

    if(_array) {
      array_delete(_array, _capacity);
      DEALLOC_KNOWN(_array,sizeof(C)*_capacity,"DArray<>");
    }
    _size = s;
    _capacity = newCapacity;
    _array = newArray;
    return false;
  } // ensure


  void reset() { ensure(0); }

  /**
   * Set array's size to @b s and that its capacity is at least @b s.
   * If the capacity is smaller, the array will expand, and all old
   * elements will be copied to the new array.
   *
   */
  void expand(size_t s)
  {
    if (_capacity >= s) {
      _size = s;
      return;
    }

    size_t oldCapacity=_capacity;
    size_t newCapacity=std::max(s, oldCapacity*2);
    void* mem = ALLOC_KNOWN(sizeof(C)*newCapacity,"DArray<>");

    C* oldArr = _array;

    _capacity = newCapacity;
    _array = static_cast<C*>(mem);

    C* optr=oldArr;
    C* nptr=_array;
    C* firstEmpty=_array+_size;
    C* afterLast=_array+_capacity;

    while(nptr!=firstEmpty) {
      C *oldAddr = optr++, *newAddr = nptr++;
      ::new (newAddr) C(std::move(*oldAddr));
      oldAddr->~C();
    }
    while(nptr!=afterLast) {
      ::new (nptr++) C();
    }
    _size = s;

    if(oldArr) {
      DEALLOC_KNOWN(oldArr,sizeof(C)*oldCapacity,"DArray<>");
    }

  } // expand

  /**
   * Set array's size to @b s and that its capacity is at least @b s.
   * If the capacity is smaller, the array will expand, and all old
   * elements will be copied to the new array.
   *
   */
  void expand(size_t s, C defVal)
  {
    size_t oldSize = _size;
    expand(s);

    if(s<=oldSize) {
      return;
    }

    for(size_t i=oldSize; i<s; i++) {
      (*this)[i] = defVal;
    }
  } // expand

  /** Return the ensured size of the array */
  inline size_t size() const { return _size; }

  /** Creates a new array that is initialized with @b value on every position */
  static DArray initialized(size_t count, const C& value=C()) {
    DArray out(count);
    out.init(count, value);
    return out;
  }

  /** Ensure that the array's size is at least @b count and
   * initialize first @b count elements of the array to @b value. */
  void init(size_t count, const C& value=C()) {
    ensure(count);
    C* ptr=_array+count;
    while(ptr!=_array) {
      *(--ptr)=value;
    }
  }

  /**
   * Ensure that the array's size is at least @b count and
   * initialize first @b count elements with values from @b src.
   *
   * @b src has to support C operator[](size_t).
   */
  template<typename Arr>
  void initFromArray(size_t count, const Arr& src) {
    ensure(count);
    C* ptr=_array+count;
    while(count) {
      *(--ptr)=src[--count];
    }
  }
  void initFromArray(size_t count, const C* src) {
    ensure(count);
    C* ptr=_array+count;
    while(count) {
      *(--ptr)=src[--count];
    }
  }

  /**
   * Initialize array by elements of the iterator. Size of the array
   * will be equal to number of elements of the iterator.
   */
  template<class It>
  void initFromIterator(It it, size_t count=0) {
    if(count) {
      ensure(count);
      C* ptr=_array;
      while(it.hasNext()) {
	*(ptr++)=it.next();
      }
    } else {
      ensure(0);
      count=0;
      while(it.hasNext()) {
	expand(++count);
	(*this)[count-1]=it.next();
      }
    }
  }

  /**
   * Creates an array initialized with all the elements of the iterator it.
   */
  template<class It>
  static DArray fromIterator(It it, size_t count=0) {
    DArray out;
    if (count != 0) {
      out.initFromIterator(it, count);
    } else if (it.knowsSize()) {
      out.initFromIterator(it, it.size());
    } else {
      out.initFromIterator(it);
    }
    return out;
  }


  /**
   * Sort first @b count items using @b Comparator::compare
   * as comparator.
   */
  template<typename Comparator>
  inline void sort(Comparator comp)
  {
    sortGen<false>(comp);
  }

  template<typename Comparator>
  inline void sortInversed(Comparator comp)
  {
    sortGen<true>(comp);
  }

  template<bool Inversed, typename Comparator>
  void sortGen(Comparator comp)
  {
    if(_size <= 1) {
      return;
    }

    // array behaves as a stack of calls to quicksort
    static DArray<size_t> ft(32);

    size_t from = 0;
    size_t to=size()-1;
    ft.ensure(to);

    size_t p = 0; // pointer to the next element in ft
    for (;;) {
      ASS(from<size() && to<size()); //checking for underflows
      ASS(from<to);
      // invariant: from < to
      size_t m = from + Random::getInteger(to-from+1);
      C mid = (*this)[m];
      size_t l = from;
      size_t r = to;
      while (l < m) {
        switch ((Inversed?-1:1)*comp.compare((*this)[l],mid))
  	{
  	case EQUAL:
  	case LESS:
  	  l++;
  	  break;
  	case GREATER:
  	  if (m == r) {
  	    (*this)[m] = (*this)[l];
  	    (*this)[l] = (*this)[m-1];
  	    (*this)[m-1] = mid;
  	    m--;
  	    r--;
  	  }
  	  else {
  	    ASS(m < r);
  	    C aux = (*this)[l];
  	    (*this)[l] = (*this)[r];
  	    (*this)[r] = aux;
  	    r--;
  	  }
  	  break;
  	}
      }
      // l == m
      // now literals in lits[from ... m-1] are smaller than lits[m]
      // and literals in lits[r+1 ... to] are greater than lits[m]
      while (m < r) {
        switch ((Inversed?-1:1)*comp.compare(mid,(*this)[m+1]))
  	{
  	case LESS:
  	  {
  	    C aux = (*this)[r];
  	    (*this)[r] = (*this)[m+1];
  	    (*this)[m+1] = aux;
  	    r--;
  	  }
  	  break;
  	case EQUAL:
  	case GREATER:
  	  (*this)[m] = (*this)[m+1];
  	  (*this)[m+1] = mid;
  	  m++;
  	}
      }
      // now literals in lits[from ... m-1] are smaller than lits[m]
      // and all literals in lits[m+1 ... to] are greater than lits[m]
      if (m+1 < to) {
        ft[p++] = m+1;
        ft[p++] = to;
      }

      to = m-1;
      if (m!=0 && from < to) {
        continue;
      }
      if (p != 0) {
        p -= 2;
        ASS(p >= 0);
        from = ft[p];
        to = ft[p+1];
        continue;
      }
      return;
    }
  }

protected:
  /** current array's size */
  size_t _size;
  /** capacity of currently allocated piece of memory */
  size_t _capacity;
  /** array's content */
  C* _array;

public:
  class Iterator
  {
  public:
    DECL_ELEMENT_TYPE(C);
    inline Iterator() : _next(0), _afterLast(0) {}
    inline Iterator(DArray& arr) : _next(arr._array),
    _afterLast(arr._array+arr._size) {}
    inline bool hasNext() { return _next!=_afterLast; }
    inline C next() { return *(_next++); }
  private:
   C* _next;
   C* _afterLast;
  };
  class ConstIterator
  {
  public:
    DECL_ELEMENT_TYPE(C);
    inline ConstIterator() : _next(0), _afterLast(0) {}
    inline ConstIterator(const DArray& arr) : _next(arr._array),
    _afterLast(arr._array+arr._size) {}
    inline bool hasNext() { return _next!=_afterLast; }
    inline const C& next() { return *(_next++); }
  private:
   const C* _next;
   const C* _afterLast;
  };
  class ReversedIterator
  {
  public:
    DECL_ELEMENT_TYPE(C);
    inline ReversedIterator(DArray& arr) : _curr(arr._array+arr._size),
    _first(arr._array) {}
    inline bool hasNext() { return _curr!=_first; }
    inline C next() { return *(--_curr); }
  private:
    C* _curr;
    C* _first;
  };

  friend std::ostream& operator<<(std::ostream& out, DArray const& self) 
  {
    ConstIterator iter(self);
    out << "[ ";
    if (iter.hasNext()) {
      out << iter.next();
      while (iter.hasNext()) {
        out << ", ";
        out << iter.next();
      }
    }
    out << " ]";
    return out;
  }
}; // class DArray

template<typename T>
VirtualIterator<T> getContentIterator(DArray<T>& arr)
{
  return pvi( typename DArray<T>::Iterator(arr) );
}

} // namespace Lib

#endif
