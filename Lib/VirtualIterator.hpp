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
 * @file VirtualIterator.hpp
 * Defines VirtualIterator which allows iterators over various
 * structures to be stored in the same object.
 */

#ifndef __VirtualIterator__
#define __VirtualIterator__

#include "Forwards.hpp"

#include "Debug/Assertion.hpp"

#include "Allocator.hpp"
#include "Exception.hpp"
#include "Reflection.hpp"

namespace Lib {

///@addtogroup Iterators
///@{

template<typename T>
  class VirtualIterator;

/**
 * Base class of objects that provide the "virtual" core of
 * @b VirtualIterator objects.
 *
 * @tparam T type returned by the iterator
 *
 * @b IteratorCore objects can be used as ordinary stack allocated
 * or static iterators as well, but in that case they must not be
 * passed to a @b VirtualIterator object as an inside.
 *
 * If used as an inside of a @b VirtualIterator object, updating
 * the reference counter @b _refCnt is done by the @b VirtualIterator
 * object, as well as calling the destructor when the counter reaches
 * zero.
 */
template<typename T>
class IteratorCore {
public:
  IteratorCore(const IteratorCore&) = delete;
  IteratorCore& operator=(const IteratorCore&) = delete;
  IteratorCore(IteratorCore&&) = default;
  IteratorCore& operator=(IteratorCore&&) = default;

  DECL_ELEMENT_TYPE(T);
  /** Create new IteratorCore object */
  IteratorCore() : _refCnt(0) {}
  /** Destroy IteratorCore object */
  virtual ~IteratorCore() { ASS(_refCnt==0); }
  /** Return true if there is a next element */
  virtual bool hasNext() = 0;
  /**
   * Return the next element
   *
   * @warning Before each call to this function, the function @b hasNext() MUST be
   * called and return true.
   */
  virtual T next() = 0;
  /** Return true if the function @b size() can be called */
  virtual bool knowsSize() const { return false; }
  /**
   * Return the total number of elements of this iterator
   *
   * The number of elements at the construction of the iterator object
   * is always returned (even when there are no more elements left).
   *
   * @warning This function can be called only if the function @b knowsSize()
   * returns true.
   */
  virtual size_t size() const { INVALID_OPERATION("This iterator cannot retrieve its size."); }

private:
  /**
   * Reference counter field used by the @b VirtualIterator object
   */
  mutable int _refCnt;

  friend class VirtualIterator<T>;
};

/**
 * Core object for @b VirtualIterator, that represents
 * an empty iterator.
 */
template<typename T>
class EmptyIterator
: public IteratorCore<T>
{
public:
  USE_ALLOCATOR(EmptyIterator);

  EmptyIterator() {}
  bool hasNext() { return false; };
  T next() { INVALID_OPERATION("next() called on EmptyIterator object"); };
  bool knowsSize() const { return true; }
  size_t size() const { return 0; }
};

/**
 * Template class of iterators that can encapsulate different implementations
 * of element retrieval through the polymorphism of @b IteratorCore objects
 *
 * @tparam T type returned by the iterator
 *
 * The @b VirtualIterator object performs reference counting on @b IteratorCore
 * objects and deletes them when the reference count drops to zero. The reference
 * count is kept in the @b IteratorCore::_refCnt field.
 *
 * @see IteratorCore
 */
template<typename T>
class VirtualIterator {
public:
  USE_ALLOCATOR(VirtualIterator);

  DECL_ELEMENT_TYPE(T);

  /** Return an empty iterator */
  static VirtualIterator getEmpty()
  {
    static VirtualIterator inst(new EmptyIterator<T>());
    return inst;
  }

  /** Return an invalid iterator */
  static VirtualIterator getInvalid()
  {
    return VirtualIterator();
  }

  /**
   * Create an uninitialized object
   *
   * When created with this constructor, the object must be assigned
   * an initialized VirtualIterator object through the @b operator=(),
   * before any of the @b hasNext(), @b next(), @b knowsSize() or @b size()
   * functions can be called.
   */
  inline
  VirtualIterator() : _core(0) {}

  /**
   * Create an object with @b core as its core.
   */
  inline
  explicit VirtualIterator(IteratorCore<T>* core) : _core(core) { _core->_refCnt++; }

  inline
  VirtualIterator(const VirtualIterator& obj) : _core(obj._core)
  {
    if(_core) {
      _core->_refCnt++;
    }
  }

  inline
  ~VirtualIterator()
  {
    if(_core) {
	_core->_refCnt--;
	if(!_core->_refCnt) {
	  delete _core;
	}
    }
  }
  VirtualIterator& operator=(const VirtualIterator& obj)
  {
    IteratorCore<T>* oldCore=_core;
    _core=obj._core;
    if(_core) {
      _core->_refCnt++;
    }
    if(oldCore) {
      oldCore->_refCnt--;
      if(!oldCore->_refCnt) {
	delete oldCore;
      }
    }
    return *this;
  }

  /**
   * Remove reference to the iterator core.
   * Return true iff the the iterator core does not exist
   * any more after return from this function.
   *
   * The returned value can be useful for asserting
   * that the iterator core (and all resources it
   * used) was indeed released.
   */
  inline
  bool drop()
  {
    if(_core) {
      _core->_refCnt--;
      if(_core->_refCnt) {
	_core=0;
	return false;
      }
      else {
	delete _core;
	_core=0;
      }
    }
    _core=0;
    return true;
  }

  /** Return true if there is a next element */
  inline
  bool hasNext()
  {
    ASS(_core);

    return _core->hasNext();
  }
  /**
   * Return the next element
   *
   * @warning Before each call to this function, the function @b hasNext() MUST be
   * called and return true.
   */
  inline
  T next()
  {
    ASS(_core);

    return _core->next();
  }

  /** Return true if the function @b size() can be called */
  bool knowsSize() const
  {
    ASS(_core);

    return _core->knowsSize();
  }

  /**
   * Return the total number of elements of this iterator
   *
   * The number of elements at the construction of the iterator object
   * is always returned (even when there are no more elements left).
   *
   * @warning This function can be called only if the function @b knowsSize()
   * returns true.
   */
  size_t size() const
  {
    ASS(_core);
    ASS(knowsSize());

    return _core->size();
  }

  /**
   * Return true if the object is invalid (i.e. uninitialized to any IteratorCore)
   */
  bool isInvalid() { return !_core; }
private:
  /** The polymorphous core of this @b VirtualIterator object */
  IteratorCore<T>* _core;
};

/**
 * Encapsulate pointer to an @b IteratorCore object @b core into a
 * @b VirtualIteratod object
 */
template<typename T>
inline
VirtualIterator<T> vi(IteratorCore<T>* core)
{
  return VirtualIterator<T>(core);
}

/**
 * Core object for virtual iterators, that can proxy any
 * object that supports has @b hasNext() and @b next() functions.
 */
template<typename T, class Inner>
class ProxyIterator
: public IteratorCore<T>
{
public:
  USE_ALLOCATOR(ProxyIterator);
  DEFAULT_CONSTRUCTORS(ProxyIterator)
  virtual ~ProxyIterator() override {}
  
  explicit ProxyIterator(Inner inn) : _inn(std::move(inn)) {}
  bool hasNext() override { return _inn.hasNext(); };
  T next() override { return _inn.next(); };
private:
  Inner _inn;
};

/**
 * Encapsulate an ordinary iterator object @b it into a VirtualIterator
 *
 * @tparam Inner type of the iterator to be encapsulated. It must be
 *   a class/struct that declares its element type through the
 *   DECL_ELEMENT_TYPE macro, and has functions @b hasNext() and
 *   @b next()
 *
 * @see DECL_ELEMENT_TYPE
 */
template<class Inner>
inline
VirtualIterator<ELEMENT_TYPE(Inner)> pvi(Inner it)
{
  return VirtualIterator<ELEMENT_TYPE(Inner)>(new ProxyIterator<ELEMENT_TYPE(Inner),Inner>(std::move(it)));
}




///@}�

}

#endif
