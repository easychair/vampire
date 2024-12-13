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
 * @file BinaryHeap.hpp
 * Defines class BinaryHeap<T, Comparator> of binary heaps.
 */

#ifndef __BinaryHeap__
#define __BinaryHeap__

#include <algorithm>
#include <limits>

#include "Debug/Assertion.hpp"
#include "Allocator.hpp"
#include "Exception.hpp"
#include "Comparison.hpp"
#include "Backtrackable.hpp"
#include "Metaiterators.hpp"

namespace Lib {

/**
 * Class BinaryHeap implements a binary minimum heap using an array which expands,
 * when additional space is needed..
 *
 * @param T a type, that will be contained in the BinaryHeap object
 * @param Comparator class, that contains a static method, that can be
 * 	called as Comparator::compare(a, b), where a and b are of type T,
 * 	and returns Lib::Comparison enumeration member. Also, in order to
 * 	use backtrackable insert, Comparator should contain static method
 * 	max() that returns object T such that all other objects T are
 * 	smaller.
 */
template <typename T, class Comparator>
class BinaryHeap
{
public:
  /** Create a new BinaryHeap */
  BinaryHeap()
  : _size(0), _capacity(0), _data(0), _data1(0)
  {
  }

  /** Deallocate the BinaryHeap */
  ~BinaryHeap()
  {
    if(_data) {
      T* ep=_data+_size;
      while(ep!=_data1) {
	(--ep)->~T();
      }
      DEALLOC_KNOWN(_data,_capacity*sizeof(T),"BinaryHeap::T");
    }
  }

  /** Make this @b BinaryHeap empty */
  void reset()
  {
    T* ep=_data+_size;
    while(ep!=_data1) {
	(--ep)->~T();
    }
    _size=0;
  }


  /** Return number of items stored in this @b BinaryHeap */
  inline
  unsigned size() const
  {
    ASS(_size>=0);
    return _size;
  }

  /** Return true, iff there are no items in the heap */
  inline
  bool isEmpty() const
  {
    ASS(_size>=0);
    return _size==0;
  }

  /** Insert an item to the heap */
  inline
  void insert(T obj)
  {
    ensureAvaiablePosition();
    _size++;
    ::new (&_data1[_size]) T(obj);
    bubbleUp(_size);
  }

  /** Return a const reference to the smallest item in the heap */
  inline
  const T& top()
  {
    ASS(!isEmpty());
    return _data[0];
  }

  /** Remove the smallest item in the heap and return it */
  inline
  T pop()
  {
    ASS(!isEmpty());
    T res=_data[0];
    _size--;
    if(_size) {
      std::swap(_data[0],_data[_size]);
      bubbleDown(1);
    }
    _data[_size].~T();
    return res;
  }

  /** Remove the smallest item in the heap together with all
   * items equal to it, and return it */
  inline
  T popWithAllEqual()
  {
    T res=pop();
    while(!isEmpty() && Comparator::compare(res, top())==EQUAL) {
      pop();
    }
    return res;
  }

  T backtrackablePop(unsigned& lastBubbleIndex)
  {
    ASS(!isEmpty());
    T res=_data[0];
    _size--;
    if(_size) {
      std::swap(_data[0],_data[_size]);
      lastBubbleIndex=bubbleDown(1);
    } else {
      lastBubbleIndex=1;
    }
    _data[_size].~T();
    return res;
  }
  inline
  T backtrackablePop(BacktrackData& bd)
  {
    unsigned lastBubbleIndex;
    T res=backtrackablePop(lastBubbleIndex);
    bd.addBacktrackObject(
	new BHPopBacktrackObject(this, res, lastBubbleIndex));
    return res;
  }

  unsigned backtrackableInsert(T obj)
  {
    ensureAvaiablePosition();
    _size++;
    ::new (&_data1[_size]) T(obj);
    return bubbleUp(_size);
  }
  inline
  void backtrackableInsert(T obj, BacktrackData& bd)
  {
    unsigned lastBubbleIndex=backtrackableInsert(obj);
    bd.addBacktrackObject(
	    new BHInsertBacktrackObject(this, lastBubbleIndex));
  }

  void backtrackInsert(unsigned lastBubbleIndex)
  {
    //We replace the inserted element with maximal possible
    //element, so that we know for sure, that when we do
    //bubbleDown() on it, a maximal element will be at the
    //last position. Also from the way how bubbleDown works
    //we know, that the heap will be exactly the same as before
    //inserting.
    _data1[lastBubbleIndex]=Comparator::max();
    bubbleDown(lastBubbleIndex);
    ASS(_data1[_size]==Comparator::max());
    _data1[_size].~T();
    _size--;
  }

  void backtrackPop(T val, unsigned lastBubbleIndex)
  {
    //During insertion, the first item is swapped with the last,
    //removed from the end of the array, and then the item at
    //the first position bubbles down, until the heap condition
    //is fulfilled. Here we reverse the process provided that
    //_lastBubbleIndex is the current index of the formerly last
    //element.
    _size++;
    ::new (&_data1[_size]) T(val);
    std::swap(_data1[_size], _data1[lastBubbleIndex]);
    //Now at the position _lastBubbleIndex is the smallest element
    //of the heap, so we know that it will bubble up to the first
    //position[1]. (There's only one way to do that, so the heap will
    //be exactly the same as before the popping occured.)
    //
    //[1] or, to be precise, to such position, that all elements
    //above will be equal to it.
    bubbleUp(lastBubbleIndex);
  }

  /**
   * Iterator on elements in the heap. It yields elements
   * in no particular order.
   */
  auto iter() const 
  { return arrayIter(_data, _size); }

  friend std::ostream& operator<<(std::ostream& out, BinaryHeap const& self) 
  { 
    out << "[";
    auto iter = self.iter();
    if (iter.hasNext()) {
      out << iter.next();
      while (iter.hasNext()) {
        out << ", " << iter.next();
      }
    }
    return out << "]";
  }



private:
  class BHPopBacktrackObject
  : public BacktrackObject
  {
  public:
    BHPopBacktrackObject(BinaryHeap* bh, T v, unsigned lastBubbleIndex)
    :_bh(bh), _val(v), _lastBubbleIndex(lastBubbleIndex) {}
    void backtrack()
    {
      _bh->backtrackPop(_val,_lastBubbleIndex);
    }
    USE_ALLOCATOR(BHPopBacktrackObject);
  private:
    BinaryHeap* _bh;
    T _val;
    unsigned _lastBubbleIndex;
  };

  class BHInsertBacktrackObject
  : public BacktrackObject
  {
  public:
    BHInsertBacktrackObject(BinaryHeap* bh, unsigned lastBubbleIndex)
    :_bh(bh), _lastBubbleIndex(lastBubbleIndex) {}
    void backtrack()
    {
      _bh->backtrackInsert(_lastBubbleIndex);
    }
    USE_ALLOCATOR(BHInsertBacktrackObject);
  private:
    BinaryHeap* _bh;
    unsigned _lastBubbleIndex;
  };

  /** Copy constructor is private and without a body, because we don't want any. */
  BinaryHeap(const BinaryHeap& obj);
  /** operator= is private and without a body, because we don't want any. */
  BinaryHeap& operator=(const BinaryHeap& obj);

  /** Make sure the heap property is not violated by the element
   * at @b index wrt its ancestors, and return its new index. */
  unsigned bubbleUp(unsigned index)
  {
    ASS(index>0 && index<=_size);
    unsigned nextIndex=index>>1;
    while(nextIndex) {
      if(Comparator::compare(_data1[index], _data1[nextIndex])==LESS) {
	std::swap(_data1[index], _data1[nextIndex]);
      } else {
	return index;
      }
      index=nextIndex;
      nextIndex=index>>1;
    }
    return 1;
  }

  /** Make sure the heap property is not violated by the element
   * at @b index wrt its descendants, and return its new index. */
  unsigned bubbleDown(unsigned index)
  {
    ASS(index>0 && index<=_size);
    unsigned nextIndex=index<<1;
    while(nextIndex<=_size) {
      if(nextIndex!=_size && Comparator::compare(_data1[index], _data1[nextIndex|1])==GREATER) {
	if(Comparator::compare(_data1[nextIndex|1], _data1[nextIndex])==GREATER) {
	  std::swap(_data1[index], _data1[nextIndex]);
	} else {
	  std::swap(_data1[index], _data1[nextIndex|1]);
	  nextIndex|=1;
	}
      } else if(Comparator::compare(_data1[index], _data1[nextIndex])==GREATER) {
	std::swap(_data1[index], _data1[nextIndex]);
      } else {
	return index;
      }
      index=nextIndex;
      nextIndex=index<<1;
    }
    return index;
  }

  /** Ensure there is at least one unused position at the end of _data array */
  inline
  void ensureAvaiablePosition()
  {
    ASS(_capacity>=_size);
    if(_capacity==_size)
      expand();
  }

  /**
   * Expand BinaryHeap to double of its current size.
   *
   * Should be called only when _capacity==_size.
   */
  void expand()
  {
    ASS(_capacity==_size);

    unsigned oldCapacity=_capacity;
    T* oldData=_data;

    _capacity= _capacity ? _capacity*2 : 4;

    void* mem = ALLOC_KNOWN(_capacity*sizeof(T),"BinaryHeap::T");
    _data = static_cast<T*>(mem);
    _data1 = _data-1;

    if(_size) {
      T* otp = oldData+_size;
      T* ntp = _data+_size;
      do {
	::new (--ntp) T(*(--otp));
	//because oldCapacity==_size, we destroy all elements of oldData array here
	otp->~T();
      } while(ntp!=_data);
    }

    if(oldData) {
      DEALLOC_KNOWN(oldData,oldCapacity*sizeof(T),"BinaryHeap::T");
    }
  }

  /** Number of entries stored in this BinaryHeap */
  unsigned _size;
  /** Size of the _data array */
  unsigned _capacity;

  /** Array containing the heap tree */
  T* _data;

  /**
   * Pointer to the T before the start of the _data
   * (we can use it for one-based access to _data)
   */
  T* _data1;

}; // class BinaryHeap

};

#endif // __BinaryHeap__


