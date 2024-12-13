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
 * @file Metaiterators.hpp
 * Defines class Metaiterators.
 */


#ifndef __Metaiterators__
#define __Metaiterators__

#include <tuple>
#include <utility>
#include <functional>

#include "Forwards.hpp"

#include "Lib/Recycled.hpp"
#include "Lib/Reflection.hpp"
#include "List.hpp"
#include "DHSet.hpp"
#include "VirtualIterator.hpp"
#include "Lib/Option.hpp"
#include "Lib/Coproduct.hpp"
#include "Debug/TimeProfiling.hpp"

namespace Lib {

///@addtogroup Iterators
///@{


/**
 * Return iterator on the container C
 *
 * The getContentIterator method makes it possible to
 * iterate on arbitrary containers, provided the @b ITERATOR_TYPE
 * macro can obtain the iterator type of the container.
 *
 * For types where the use of the @b ITERATOR_TYPE macro is not
 * suitable, this method can be overloaded to obtain the
 * iterator by some other means.
 */
template<class C>
ITERATOR_TYPE(C) getContentIterator(C& c)
{
  return ITERATOR_TYPE(C)(c);
}

/**
 * Iterator returning a single element
 *
 * The single element is being passed to the constructor of the iterator.
 */
template<typename T>
class SingletonIterator
{
public:
  DECL_ELEMENT_TYPE(T);
  DEFAULT_CONSTRUCTORS(SingletonIterator)
  explicit SingletonIterator(T el) : _finished(false), _el(el) {}
  inline bool hasNext() { return !_finished; };
  inline T next() { ASS(!_finished); _finished=true; return _el; };
  inline bool knowsSize() const { return true; }
  inline size_t size() const { return 1; }
private:
  bool _finished;
  T _el;
};

/**
 * Return iterator returning @b el as a single element
 *
 * @see SingletonIterator
 */
template<typename T>
inline
SingletonIterator<T> getSingletonIterator(T el)
{
  return SingletonIterator<T>(el);
}

/**
 * A functor class that returns true if the argument is non-zero
 *
 * The nonzeroness is tested by @b x!=0 .
 */
struct NonzeroFn
{
  template<typename T>
  bool operator()(T obj)
  {
    return obj!=0;
  }
};

template<bool cond> struct __doDeletion
{ template<class Iter> void operator()(Iter& iter) { } };

template<> struct __doDeletion<true>
{ template<class Iter> void operator()(Iter& iter) { iter.del(); } };


/**
 * Iterator class that returns elements of the inner iterator
 * for which the functor returns true
 *
 * @tparam Inner type of the inner iterator
 * @tparam Functor type of the functor used for filtering the
 *   elements returned by the inner iterator
 */
template<class Inner, class Functor, bool deleteFilteredOut = false>
class FilteredIterator
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Inner));
  DEFAULT_CONSTRUCTORS(FilteredIterator)

  FilteredIterator(Inner inn, Functor func)
  : _func(std::move(func)), _inn(std::move(inn)), _next() {}


  bool hasNext()
  {
    if(_next.isSome()) {
      return true;
    }
    while(_inn.hasNext()) {
      OWN_ELEMENT_TYPE next = move_if_value<OWN_ELEMENT_TYPE>(_inn.next());
      if(_func(next)) {
        _next = Option<OWN_ELEMENT_TYPE>(move_if_value<OWN_ELEMENT_TYPE>(next));
        return true;
      } else {
        __doDeletion<deleteFilteredOut>{}(_inn);
      }
    }
    return false;
  }

  OWN_ELEMENT_TYPE next()
  {
    ALWAYS(hasNext());
    ASS(_next.isSome());
    OWN_ELEMENT_TYPE out = move_if_value<OWN_ELEMENT_TYPE>(_next.unwrap());
    _next = Option<OWN_ELEMENT_TYPE>();
    return out;
  }

  bool knowsSize() const { return false; }
  size_t size() const { ASSERTION_VIOLATION }

private:
  
  Functor _func;
  Inner _inn;
  Option<OWN_ELEMENT_TYPE> _next;
};

template<class Iter, class Pred>
class TakeWhileIter
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Iter));
  DEFAULT_CONSTRUCTORS(TakeWhileIter)

  TakeWhileIter(Iter iter, Pred pred)
  : _pred(std::move(pred)), _iter(std::move(iter)), _next(), _break(false) {}

  bool hasNext()
  {
    if (_break) return false;
    if(_next.isSome()) return true;
    
    while(_iter.hasNext()) {
      OWN_ELEMENT_TYPE next = move_if_value<OWN_ELEMENT_TYPE>(_iter.next());
      if(_pred(next)) {
        _next = Option<OWN_ELEMENT_TYPE>(move_if_value<OWN_ELEMENT_TYPE>(next));
        return true;
      } else {
        _break = true;
        return false;
      }
    }
    return false;
  }

  OWN_ELEMENT_TYPE next()
  {
    ASS(!_break)
    ALWAYS(hasNext());
    ASS(_next.isSome());
    OWN_ELEMENT_TYPE out = move_if_value<OWN_ELEMENT_TYPE>(_next.unwrap());
    _next = Option<OWN_ELEMENT_TYPE>();
    return out;
  }

private:
  
  Pred _pred;
  Iter _iter;
  Option<OWN_ELEMENT_TYPE> _next;
  bool _break;
};

template<class I1, class I2, class Cmp>
class SortedIterDiff {
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(I1));
  DEFAULT_CONSTRUCTORS(SortedIterDiff)

  SortedIterDiff(I1 i1, I2 i2, Cmp cmp)
    : _i1(std::move(i1))
    , _i2(std::move(i2))
    , _curr1()
    , _curr2(someIf(_i2.hasNext(), [&]() -> OWN_ELEMENT_TYPE { return move_if_value<OWN_ELEMENT_TYPE>(_i2.next()); })) 
    , _cmp(std::move(cmp))
    {}

  void moveToNext() 
  {
#   if VDEBUG
    Option<OWN_ELEMENT_TYPE> old1;
#   endif
    while (_curr1.isNone() && _i1.hasNext()) {
      _curr1 = Option<OWN_ELEMENT_TYPE>(move_if_value<OWN_ELEMENT_TYPE>(_i1.next()));
      ASS_REP(!_curr1.isSome() || !old1.isSome() || _cmp(*old1, *_curr1) <= 0, "iterator I1 must be sorted");
      while (_curr2.isSome() && _cmp(*_curr2, *_curr1) < 0) {
#       if VDEBUG
        Option<OWN_ELEMENT_TYPE> old2 = _curr2.take();
#       endif
        _curr2 = someIf(_i2.hasNext(), [&]() -> OWN_ELEMENT_TYPE 
            { return move_if_value<OWN_ELEMENT_TYPE>(_i2.next()); });
        ASS_REP(!_curr2.isSome() || !old2.isSome() || _cmp(*old2, *_curr2) <= 0, "iterator I2 must be sorted");
      }
      if (( _curr1.isSome() && _curr2.isSome() && _cmp(*_curr1, *_curr2) == 0 )
          || (_curr1.isNone() && _curr2.isNone())) {
#   if VDEBUG
    old1 = _curr1.take();
#   endif
        _curr1 = Option<OWN_ELEMENT_TYPE>();
        _curr2 = someIf(_i2.hasNext(), [&]() -> OWN_ELEMENT_TYPE 
            { return move_if_value<OWN_ELEMENT_TYPE>(_i2.next()); });
      }
    }
  }

  bool hasNext()
  {
    moveToNext();
    return _curr1.isSome();
  }

  OWN_ELEMENT_TYPE next()
  {
    moveToNext();
    return move_if_value<OWN_ELEMENT_TYPE>(_curr1.take().unwrap());
  }
private:
  I1 _i1;
  I2 _i2;
  Option<OWN_ELEMENT_TYPE> _curr1;
  Option<OWN_ELEMENT_TYPE> _curr2;
  Cmp _cmp;
};

/**
 * Iterator that maps the contents of another iterator by a function. Whenever the function retuns a non-empty Option
 * this iterator will return the corresponding value. 
 */
template<class Inner, class Functor>
class FilterMapIter
{
public:
  DECL_ELEMENT_TYPE(typename std::invoke_result<Functor, ELEMENT_TYPE(Inner)>::type::Content);
  DEFAULT_CONSTRUCTORS(FilterMapIter)

  FilterMapIter(Inner inn, Functor func)
  : _func(std::move(func)), _inn(std::move(inn)), _next() {}

  bool hasNext()
  {
    if(_next.isSome()) {
      return true;
    }
    while(_inn.hasNext()) {
      _next = _func(move_if_value<ELEMENT_TYPE(Inner)>(_inn.next()));
      if(_next.isSome()) {
        return true;
      }
    }
    return false;
  };

  OWN_ELEMENT_TYPE next()
  {
    ALWAYS(hasNext());
    ASS(_next.isSome());
    OWN_ELEMENT_TYPE out = move_if_value<OWN_ELEMENT_TYPE>(_next.unwrap());
    _next = Option<OWN_ELEMENT_TYPE>();
    return out;
  };

private:
  Functor _func;
  Inner _inn;
  Option<OWN_ELEMENT_TYPE> _next;
};

/**
 * Return an iterator object that returns elements of the @b inn iterator
 * for which the functor @b func returns true
 *
 * @see FilteredIterator
 */
template<class Inner, class Functor>
inline
FilteredIterator<Inner,Functor> getFilteredIterator(Inner inn, Functor func)
{
  return FilteredIterator<Inner,Functor>(std::move(inn), std::move(func));
}

template<class Inner, class Functor>
inline
auto getFilteredDelIterator(Inner inn, Functor func)
{ return FilteredIterator<Inner, Functor, /*deleteFilteredOut=*/true>(inn, func); }


/**
 * Iterator class that returns elements of an inner iterator
 * only until the specified functor returns false for some element
 * (this element is already not returned)
 */
template<class Inner, class Functor>
class WhileLimitedIterator
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Inner));
  DEFAULT_CONSTRUCTORS(WhileLimitedIterator)
  WhileLimitedIterator(Inner inn, Functor func)
  : _func(std::move(func)), _inn(std::move(inn)), _nextStored(false) {}
  bool hasNext()
  {
    if(!_nextStored) {
      if(!_inn.hasNext()) {
        return false;
      }
      _next=_inn.next();
      _nextStored=true;
    }
    return _func(_next);
  };
  OWN_ELEMENT_TYPE next()
  {
    if(!_nextStored) {
      ALWAYS(hasNext());
      ASS(_nextStored);
    }
    _nextStored=false;
    return _next;
  };
private:
  Functor _func;
  Inner _inn;
  OWN_ELEMENT_TYPE _next;
  bool _nextStored;
};

/**
 * Return iterator object that returns elements of an inner iterator
 * @b inn only until the functor @b func returns false for some element
 * (this element is already not returned)
 *
 * @see WhileLimitedIterator
 */
template<class Inner, class Functor>
inline
WhileLimitedIterator<Inner,Functor> getWhileLimitedIterator(Inner inn, Functor func)
{
  return WhileLimitedIterator<Inner,Functor>(inn, func);
}


/**
 * Iterator that concatenates two other iterators
 *
 * The @b knowsSize() and @b size() functions of this iterator can be
 * called only if both underlying iterators contain these functions.
 */
template<class It1,class It2>
class CatIterator
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(It1));
  DEFAULT_CONSTRUCTORS(CatIterator)

  CatIterator(It1 it1, It2 it2)
  	:_first(true), _it1(std::move(it1)), _it2(std::move(it2)) {}
  bool hasNext()
  {
    if(_first) {
      if(_it1.hasNext()) {
	return true;
      }
      _first=false;
    }
    return  _it2.hasNext();
  };
  /**
   * Return the next value
   * @warning hasNext() must have been called before
   */
  OWN_ELEMENT_TYPE next()
  {
    ALWAYS(hasNext())
    if(_first) {
      //_it1 contains the next value, as hasNext must have
      //been called before. (It would have updated the
      //_first value otherwise.)
      return _it1.next();
    }
    return  _it2.next();
  };

  /**
   * Return true the size of the iterator can be obtained
   *
   * This function can be called only if both underlying iterators contain
   * the @b knowsSize() function.
   */
  bool knowsSize() const { return _it1.knowsSize() && _it2.knowsSize(); }
  /**
   * Return the initial number of elements of this iterator
   *
   * This function can be called only if both underlying iterators contain
   * the @b size() function, and if the @b knowsSize() function returns true.
   */
  size_t size() const { return _it1.size()+_it2.size(); }
private:
  /** False if we have already iterated through the first iterator */
  bool _first;
  It1 _it1;
  It2 _it2;
};

/**
 * Iterator that transforms elements of its inner iterator by
 * a specified functor
 *
 * The @b knowsSize() and @b size() functions of this iterator can be
 * called only if the underlying iterator contains these functions.
 */
template<typename Inner, typename Functor, typename ResultType=std::invoke_result_t<Functor, ELEMENT_TYPE(Inner)>>
class MappingIterator
{
public:
  DECL_ELEMENT_TYPE(ResultType);
  DEFAULT_CONSTRUCTORS(MappingIterator)
  explicit MappingIterator(Inner inner, Functor func)
  : _func(std::move(func)), _inner(std::move(inner)) {}
  inline bool hasNext() { return _inner.hasNext(); };
  inline ResultType next() { return _func(move_if_value<ELEMENT_TYPE(Inner)>(_inner.next())); };

  /**
   * Return true the size of the iterator can be obtained
   *
   * This function can be called only if the underlying iterator contains
   * the @b knowsSize() function.
   */
  inline bool knowsSize() const { return _inner.knowsSize(); }
  /**
   * Return the initial number of elements of this iterator
   *
   * This function can be called only if the underlying iterator contains
   * the @b size() function, and if the @b knowsSize() function returns true.
   */
  inline size_t size() const { return _inner.size(); }

  auto reverse() && 
  { return MappingIterator<decltype(std::move(_inner).reverse()), Functor, ResultType>(std::move(_inner).reverse(), std::move(_func)); }
private:
  Functor _func;
  Inner _inner;
};

template<typename Iter>
class IterTraits;

template<typename Iter>
IterTraits<Iter> iterTraits(Iter);

template<typename Iter>
class IterFromTryNext
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Iter));
  DEFAULT_CONSTRUCTORS(IterFromTryNext)
  IterFromTryNext(Iter iter)
    : _iter(std::move(iter)) {}


  bool hasNext() {
    if (!_nextCached) { _nextCached = some(_iter.tryNext()); }
    return _nextCached->isSome();
  }

  OWN_ELEMENT_TYPE next() {
    if (!_nextCached) { _nextCached = some(_iter.tryNext()); }
    return _nextCached.take()->unwrap();
  }

private:
  Iter _iter;
  Option<Option<OWN_ELEMENT_TYPE>> _nextCached;
};

template<typename Iter>
IterFromTryNext<Iter> iterFromTryNext(Iter iter)
{ return IterFromTryNext<Iter>(std::move(iter)); }


/**
 * Return iterator that returns elements of @b it transformed by
 * the functor @b f
 *
 * @see MappingIterator
 */
template<typename Inner, typename Functor>
MappingIterator<Inner,Functor,ResultOf<Functor, ELEMENT_TYPE(Inner)>> getMappingIterator(Inner it, Functor f)
{
  return MappingIterator<Inner,Functor, ResultOf<Functor, ELEMENT_TYPE(Inner)>>(std::move(it), std::move(f));
}

/**
 * Return iterator that returns elements of @b it transformed by
 * the functor @b f
 *
 * @see MappingIterator
 */
template<typename ResultType, typename Inner, typename Functor>
MappingIterator<Inner,Functor,ResultType> getMappingIteratorKnownRes(Inner it, Functor f)
{
  return MappingIterator<Inner,Functor,ResultType>(it, f);
}


/**
 * Iterator that uses elements of its inner iterator as argments to
 * single-parameter constructor @b Constructor, and yields created
 * objects.
 */
template<typename Constructor, typename Inner>
class ConstructingIterator
{
public:
  DECL_ELEMENT_TYPE(Constructor*);
  DEFAULT_CONSTRUCTORS(ConstructingIterator)
  explicit ConstructingIterator(Inner inner)
  : _inner(inner) {}
  inline bool hasNext() { return _inner.hasNext(); };
  inline Constructor* next() { return new Constructor(_inner.next()); };
private:
  Inner _inner;
};

/**
 * Return iterator that uses elements of @b it as arguments to
 * the single-paramater constructor of type @b Constructor and
 * returns the created elements
 *
 * @see ConstructingIterator
 */
template<typename Constructor, typename Inner>
inline
ConstructingIterator<Constructor,Inner> getConstructingIterator(Inner it)
{
  return ConstructingIterator<Constructor,Inner>(it);
}


/**
 * Iterator that takes iterator over iterators as its argument and
 * flattens it, returning elements of the inner iterators.
 *
 * @tparam Master The outer iterator to be flattened. It must be
 *   an iterator over iterators, which also means that the macro
 *   @b ELEMENT_TYPE() must be applicable to the result of
 *   @b ELEMENT_TYPE(Master)
 */
template<typename Master>
class FlatteningIterator
{
public:
  using Inner = ELEMENT_TYPE(Master);
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Inner));
  DEFAULT_CONSTRUCTORS(FlatteningIterator)

  explicit FlatteningIterator(Master master)
  : _master(std::move(master))
  , _current(_master.hasNext() 
        ? Option<Inner>(move_if_value<Inner>(_master.next()))
        : Option<Inner>())
  { }

  bool hasNext()
  {
    while (_current.isSome()) {
      if (_current->hasNext()) {
        return true;
      } else {
        _current = _master.hasNext() 
          ? Option<Inner>(move_if_value<Inner>(_master.next())) 
          : Option<Inner>();
      }
    }
    return false;
  }

  inline
  ELEMENT_TYPE(FlatteningIterator) next()
  {
    ASS(_current.isSome());
    ASS(_current.unwrap().hasNext());
    return move_if_value<OWN_ELEMENT_TYPE>(_current.unwrap().next());
  }
private:
  Master _master;
  Option<Inner> _current;
};

/**
 * Return iterator that flattens the iterator over iterators @b it
 * into an iterator over elements the inner elements
 *
 * @b it must be an iterator over iterators
 *
 * @see FlatteningIterator, FlatteningIterator<VirtualIterator<VirtualIterator<T>>>
 */
template<typename T>
inline
FlatteningIterator<T> getFlattenedIterator(T it)
{ return FlatteningIterator<T>(std::move(it)); }

template<class Inner, class Functor>
using FlatMapIter = FlatteningIterator<MappingIterator<Inner,Functor>>;

/**
 * Return iterator that applies functor @b f to elements of the @b it
 * iterator, treats the result as iterators and flattens them
 *
 * This function is a combination of the @b getMappingIterator() and
 * @b getFlattenedIterator() functions.
 *
 * The functor @b f must return an iterator
 *
 * @see getMappingIterator(), getFlattenedIterator()
 */
template<typename Inner, typename Functor>
inline
FlatMapIter<Inner,Functor> getMapAndFlattenIterator(Inner it, Functor f)
{
  return FlatteningIterator<MappingIterator<Inner,Functor> >(
	  MappingIterator<Inner,Functor>(std::move(it), f) );
}


/**
 * Iterator that in its constructor stores elements of an inner iterator
 * and then returns these elements later in a deterministic order, skipping
 * the duplicate ones
 *
 * The iterator object does not contain the copy constructor or
 * the operator=. If this behavior is required, it should be created
 * on the heap and pointer to it put inside a VirtualIterator object.
 *
 * @see VirtualIterator
 */
template<class Inner>
class UniquePersistentIterator
: public IteratorCore<ELEMENT_TYPE(Inner)>
{
public:
  typedef ELEMENT_TYPE(Inner) T;
private:
  typedef List<T> ItemList;
public:

  explicit UniquePersistentIterator(Inner& inn)
  {
    _items=getUniqueItemList(inn, _size);
  }
  ~UniquePersistentIterator()
  {
    if(_items) {
      ItemList::destroy(_items);
    }
  }
  inline bool hasNext() { return _items; };
  inline T next()
  {
    return ItemList::pop(_items);
  };

  inline bool knowsSize() const { return true; }
  inline size_t size() const { return _size; }
private:
  typedef DHSet<T> ItemSet;

  static ItemList* getUniqueItemList(Inner& inn, size_t& sizeRef)
  {
    ItemList* res=0;
    Recycled<ItemSet> iset;

    sizeRef=0;
    while(inn.hasNext()) {
      T el=inn.next();
      if(iset->insert(el)) {
	ItemList::push(el, res);
	sizeRef++;
      }
    }

    return res;
  }

  ItemList* _items;
  size_t _size;
};

/**
 * Return iterator that stores values of @b it in its constructor,
 * and then yields them in a deterministic order, skipping duplicate values
 *
 * After the call to this function, the iterator @b it and any resources
 * it holds may be released, since the elements are stored independently
 * of it.
 *
 * @see UniquePersistentIterator
 */
template<class Inner>
inline
VirtualIterator<ELEMENT_TYPE(Inner)> getUniquePersistentIterator(Inner it)
{
  if(!it.hasNext()) {
    return VirtualIterator<ELEMENT_TYPE(Inner)>::getEmpty();
  }
  return vi( new UniquePersistentIterator<Inner>(it) );
}


/**
 * Return iterator that stores values of the iterator pointed to by @b it
 * in its constructor, and then yields them in a deterministic order,
 * skipping duplicate values
 *
 * After the call to this function, the iterator pointed to by @b it and
 * any resources it holds may be released, since the elements are stored
 * independently of it.
 *
 * @see UniquePersistentIterator
 */
template<class Inner>
inline
VirtualIterator<ELEMENT_TYPE(Inner)> getUniquePersistentIteratorFromPtr(Inner* it)
{
  if(!it->hasNext()) {
    return VirtualIterator<ELEMENT_TYPE(Inner)>::getEmpty();
  }
  return vi( new UniquePersistentIterator<Inner>(*it) );
}

/**
 * Remove duplicate elements from the container @c cont
 */
template<class Container>
void makeUnique(Container& cont)
{
  VirtualIterator<ELEMENT_TYPE(Container)> uniqueIt = pvi(
      getUniquePersistentIterator(ITERATOR_TYPE(Container)(cont)) );
}

/**
 * Return number of elements in iterator @c it
 */
template<class It>
size_t countIteratorElements(It it)
{
  size_t res = 0;
  while(it.hasNext()) {
    it.next();
    res++;
  }
  return res;
}


/**
 * Iterator that goes from object @b from to the object @b to using the
 * postfix @b operator++. (The objects are passed in the constructor.)
 * The object @b to is not returned.
 */
template<typename T>
class RangeIterator
{
public:
  DECL_ELEMENT_TYPE(T);
  DEFAULT_CONSTRUCTORS(RangeIterator)
  inline
  RangeIterator(T from, T to)
  : _next(from), _from(from), _to(to) {}
  inline bool hasNext() { return _next<_to; };
  inline T next() { return _next++; };
  inline bool knowsSize() const { return true; }
  inline size_t size() const { return (_to>_from) ? (_to-_from) : 0; }
  auto reverse() &&
  { 
    auto to = _to;
    return getMappingIterator(std::move(*this), [to](auto i) { return to - 1 - i; }); 
  }
private:
  T _next;
  T _from;
  T _to;
};

/**
 * Return iterator that goes from object @b from to the object @b to
 * using the postfix @b operator++; the object @b to is not returned
 *
 * @see RangeIterator
 */
template<typename T>
RangeIterator<T> getRangeIterator(T from, T to)
{
  return RangeIterator<T>(from, to);
}

template<typename T>
class CombinationIterator
{
public:
  DECL_ELEMENT_TYPE(std::pair<T,T>);
  DEFAULT_CONSTRUCTORS(CombinationIterator)
  CombinationIterator(T from, T to)
  : _first(from), _second(from), _afterLast(to)
  {
    ASS_LE(from,to);
    if(from!=to) {
      if(from+1==to) {
	_second=_afterLast;
      } else {
	moveToNext();
      }
    }
  }
  inline bool hasNext()
  { ASS_LE(_first,_afterLast); return _second!=_afterLast; }
  std::pair<T,T> next()
  {
    ASS(hasNext());
    std::pair<T,T> res=std::pair<T,T>(_first,_second);
    moveToNext();
    return res;
  }
private:
  void moveToNext()
  {
    _second++;
    ASS_LE(_second,_afterLast);
    if(_second==_afterLast) {
      _first++;
      _second=_first;
      _second++;
      //now, if _second==_afterLast, there's no combination left
    }
  }
  T _first;
  T _second;
  T _afterLast;
};

/**
 * Return iterator, that yields all unordered pairs from set {@b from,
 * (@b from)+1, (@b from)+2,..., (@b to)-1}. (The addition is performed
 * by the operator++.) For a singleton set, nothing is yielded.
 */
template<typename T>
inline
CombinationIterator<T> getCombinationIterator(T from, T to)
{
  return CombinationIterator<T>(from, to);
}

template<typename T>
class Combination2Iterator
{
public:
  DECL_ELEMENT_TYPE(std::pair<T,T>);
  DEFAULT_CONSTRUCTORS(Combination2Iterator)
  Combination2Iterator(T from, T to1, T to2)
  : _first(from), _second(from), _afterLast1(to1), _afterLast2(to2)
  {
    ASS_LE(from,to1);
    ASS_LE(to1,to2);
    if(from!=to1) {
      moveToNext();
    }
  }
  inline bool hasNext()
  { return _first!=_afterLast1 && _second!=_afterLast2; }
  std::pair<T,T> next()
  {
    ASS(hasNext());
    std::pair<T,T> res=std::pair<T,T>(_first,_second);
    ASS_LE(_first,_afterLast1);
    ASS_LE(_second,_afterLast2);
    moveToNext();
    return res;
  }
private:
  void moveToNext()
  {
    _second++;
    ASS_LE(_second,_afterLast2);
    if(_second==_afterLast2) {
      _first++;
      _second=_first;
      _second++;
      //now, if _second==_afterLast, there's no combination left
    }
  }
  T _first;
  T _second;
  T _afterLast1;
  T _afterLast2;
};

/**
 * Return iterator, that yields all unordered pairs from set {@b from,
 * (@b from)+1, (@b from)+2,..., (@b to2)-1} where one of the pair is
 * less than @b to1. (The addition is performed by the operator++.)
 * For a singleton set, nothing is yielded.
 *
 * Parameter @b to1 must be less than or equal to @b to2.
 */
template<typename T>
inline
Combination2Iterator<T> getCombinationIterator(T from, T to1, T to2)
{
  return Combination2Iterator<T>(from, to1, to2);
}


/**
 * Wraps a context around specified iterator.
 *
 * Context is an object of type @b Ctx with methods
 * bool enter(T)
 * void leave(T)
 * where @b T is the return type of inner iterator.
 * Method @b enter is called before an element of inner
 * iterator is yielded (with this element as a parameter).
 * If @b enter returns false, the element is skipped (@b leave is
 * not called for it).
 * @b leave is called when an element becomes no longer
 * needed (after the hasNext method is called next time,
 * or when the iterator is being destroyed).
 */
template<class Inner, class Ctx>
class ContextualIterator
{
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Inner));
  DEFAULT_CONSTRUCTORS(ContextualIterator)

  ContextualIterator(Inner iit, Ctx context)
  : _inContext(false), _used(true), _context(std::move(context)), _iit(std::move(iit)) {}

  ~ContextualIterator()
  {
    assureContextLeft();
  }
  bool hasNext()
  {
    if(!_used) {
      return true;
    }
    assureContextLeft();
    do {
      if(!_iit.hasNext()) {
	return false;
      }
      _current = Option<ELEMENT_TYPE(Inner)>(_iit.next());
    } while (!_context.enter(_current.unwrap()));
    _inContext=true;

    _used=false;
    return true;
  }
  inline
  ELEMENT_TYPE(Inner) next()
  {
    ASS(!_used);
    _used=true;
    return move_if_value<ELEMENT_TYPE(Inner)>(_current.unwrap());
  }
private:
  void assureContextLeft()
  {
    if(_inContext) {
      _context.leave(_current.unwrap());
      _inContext=false;
    }
  }

  bool _inContext;
  bool _used;
  Ctx _context;
  Inner _iit;
  Option<ELEMENT_TYPE(Inner)> _current;
};

template<class Inner, class Ctx>
inline
ContextualIterator<Inner,Ctx> getContextualIterator(Inner it, Ctx context)
{
  return ContextualIterator<Inner,Ctx>(it, context);
}

#if VTIME_PROFILING

template<class Iter>
class TimeTracedIter
{
  const char* _name;
  Iter _iter;
public:
  TimeTracedIter(const char* name, Iter iter) 
    : _name(name)
    , _iter(std::move(iter)) 
  {}

  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Iter));

  OWN_ELEMENT_TYPE next() { TIME_TRACE(_name); return _iter.next(); }
  bool hasNext() { TIME_TRACE(_name); return _iter.hasNext(); }

  bool knowsSize() const 
  { return _iter.knowsSize(); }

  size_t size() const
  { return _iter.size(); }

};

template<class Iter>
auto timeTraceIter(const char* name, Iter iter) 
{ return TimeTracedIter<Iter>(name, std::move(iter)); }

#define TIME_TRACE_ITER(name, iter) timeTraceIter(name, iter)

#else // !VTIME_PROFILING

#define TIME_TRACE_ITER(name, iter) iter

#endif // VTIME_PROFILING


/**
 * Return true iff @c it1 and it2 contain the same values in the same order
 */
template<class It1, class It2>
bool iteratorsEqual(It1 it1, It2 it2)
{
  while(it1.hasNext()) {
    if(!it2.hasNext()) {
      return false;
    }
    if(it1.next()!=it2.next()) {
      return false;
    }
  }
  return !it2.hasNext();
}

template<typename T>
static bool lessThan(T a, T b) { return a<b; }

/**
 * Return true iff @c it is sorted according to the default < ordering
 * (each element is less than its successive element).
 */
template<class It>
bool isSorted(It it)
{
  if(!it.hasNext()) { return true; }

  ELEMENT_TYPE(It) prev = it.next();

  while(it.hasNext()) {
    ELEMENT_TYPE(It) curr = it.next();
    if(!(prev<curr)) {
      return false;
    }
    prev = curr;
  }
  return true;
}

/**
 * Return true iff @c it is sorted according to ordering specified by lessThan
 * (each element is less than its successive element).
 *
 * Transitivity of @c lessThan is assumed.
 */
template<class It, typename Pred>
bool isSorted(It it, Pred lessThan)
{
  if(!it.hasNext()) { return true; }

  ELEMENT_TYPE(It) prev = it.next();

  while(it.hasNext()) {
    ELEMENT_TYPE(It) curr = it.next();
    if(!lessThan(prev, curr)) {
      return false;
    }
    prev = curr;
  }
  return true;
}

/**
 * Return true iff pred returns true for all elements of it
 */
template<class It, typename Pred>
bool forAll(It it, Pred pred)
{
  while(it.hasNext()) {
    if(!pred(it.next())) {
      return false;
    }
  }
  return true;
}

/**
 * Return first element for which pred returns true.
 * There must be an element for which true is returned in the iterator.
 */
template<class It, typename Pred>
ELEMENT_TYPE(It) getFirstTrue(It it, Pred pred)
{
  while(it.hasNext()) {
    ELEMENT_TYPE(It) el = it.next();
    if(pred(el)) {
      return el;
    }
  }
  ASSERTION_VIOLATION;
}

/**
 * Do folding on iterator it using fn which is a function taking
 * iterator element as first argument and the intermediate result
 * as the second.
 */
template<class It, typename Fun, typename Res>
Res fold(It it, Fun fn, Res init)
{
  Res res = init;
  while(it.hasNext()) {
    res = fn(it.next(), res);
  }
  return res;
}

/**
 * Do folding on iterator it using fn which is a function taking
 * iterator element as first argument and the intermediate result
 * as the second.
 * it must be a non-empty iterator.
 */
template<class It, typename Fun>
ELEMENT_TYPE(It) fold(It it, Fun fn)
{
  ALWAYS(it.hasNext());
  ELEMENT_TYPE(It) init = it.next();
  return fold(it,fn,init);
}

/** sum function, useful for fold */
template<typename T>
T sumFn(T a1, T a2) { return a1+a2; }

/** max function, useful for fold */
template<typename T>
T maxFn(T a1, T a2) { return std::max(a1,a2); }

/** min function, useful for fold */
template<typename T>
T minFn(T a1, T a2) { return std::min(a1,a2); }


template<class It>
struct StmJoinAuxStruct
{
  StmJoinAuxStruct(std::string glue, It it) : _glue(glue), _it(it) {}
  std::string _glue;
  It _it;
};

template<class It>
StmJoinAuxStruct<It> join(std::string glue, It it)
{
  return StmJoinAuxStruct<It>(glue, it);
}
template<typename It>
std::ostream& operator<< (std::ostream& out, const StmJoinAuxStruct<It>& info )
{
  It it = info._it;
  while(it.hasNext()) {
    out << it.next();
    if(it.hasNext()) {
      out << info._glue;
    }
  }
  return out;
}


/**
 * Split iterator @c it into two iterators in the element satisfying
 * the @c edge predicate. If there is no element satisfying @c edge,
 * return false, put the original iterator into res1 and make res2 empty.
 * The first element for which edge succeeded is not present in any
 * of the resulting iterators.
 */
template<class It, class Pred>
bool splitIterator(It it, Pred edge, VirtualIterator<ELEMENT_TYPE(It)>& res1, VirtualIterator<ELEMENT_TYPE(It)>& res2)
{
  typedef ELEMENT_TYPE(It) T;

  bool success = false;
  List<T>* firstPart = 0;
  while(it.hasNext()) {
    T itm = it.next();
    if(edge(itm)) {
      success = true;
      break;
    }
    List<T>::push(itm, firstPart);
  }
  firstPart = firstPart->reverse();
  res1 = fpvi(typename List<T>::DestructiveIterator(firstPart));
  res2 = pvi(it);
  return success;
}

template<typename Inner>
struct NegPred
{
  NegPred(const Inner& inner) : _inner(inner) {}
  template<typename Arg>
  bool operator()(Arg a) { return !_inner(a); }
private:
  Inner _inner;
};

template<typename Inner>
NegPred<Inner> negPred(Inner inner) {
  return NegPred<Inner>(inner);
}

template<typename T>
struct ConstEqPred
{
  ConstEqPred(const T& val) : _val(val) {}
  template<typename Arg>
  bool operator()(Arg a) { return a==_val; }
private:
  T _val;
};

template<typename T>
ConstEqPred<T> constEqPred(const T& val) {
  return ConstEqPred<T>(val);
}

template<typename OuterFn, typename InnerFn>
struct CompositionFn {
  CompositionFn(OuterFn outer, InnerFn inner)
   : _outer(outer), _inner(inner) { }

  template<typename Arg>
  std::invoke_result_t<OuterFn, std::invoke_result_t<InnerFn, Arg>> operator()(Arg a) {
    return _outer(_inner(a));
  }
private:
  OuterFn _outer;
  InnerFn _inner;
};

template<typename OuterFn, typename InnerFn>
CompositionFn<OuterFn,InnerFn> getCompositionFn(OuterFn outer, InnerFn inner)
{
  return CompositionFn<OuterFn,InnerFn>(outer,inner);
}

template<class P>
struct GetFirstOfPair {
  typename P::first_type operator()(P p) {
    return p.first;
  }
};

template<class P>
struct GetSecondOfPair {
  typename P::second_type operator()(P p) {
    return p.second;
  }
};


template<class... Is>
class CoproductIter 
{
  Coproduct<Is...> _inner;
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(TypeList::Get<0, TypeList::List<Is...>>));
  DEFAULT_CONSTRUCTORS(CoproductIter)

  CoproductIter(Coproduct<Is...> i) : _inner(Coproduct<Is...>(std::move(i))) {}

  bool hasNext() const
  { 
    Coproduct<Is...> const& inner = _inner;;
    return inner.apply([](auto& x) { return x.hasNext();}); }

  bool hasNext()
  { 
    Coproduct<Is...> & inner = _inner;;
    return inner.apply([](auto& x) { return x.hasNext();}); }

  OWN_ELEMENT_TYPE next()
  { 
    Coproduct<Is...> & inner = _inner;;
    return inner.apply([](auto&& x) { return x.next();}); }

  bool knowsSize() const 
  { 
    Coproduct<Is...> const& inner = _inner;;
    return inner.apply([](auto& x) { return x.knowsSize();}); }

  size_t size() const
  { 
    Coproduct<Is...> const& inner = _inner;;
    return inner.apply([](auto& x) { return x.size();}); }
};

template<class... Is>
auto coproductIter(Coproduct<Is...> is)
{ return iterTraits(CoproductIter<Is...>(std::move(is))); }

template<class Closure>
struct ApplyClosure {
   template<class A> 
   using apply = std::invoke_result_t<Closure, A>;
};

template<unsigned N>
struct Times {
  template<class C>
  using apply = Constant<C::value * N>;
};

template<unsigned N>
struct Plus {
  template<class C>
  using apply = Constant<C::value + N>;
};

#define LAZY(...)  [&]() { return __VA_ARGS__; }

template<class A
  , std::enable_if_t<!std::is_invocable_v<A>, bool> = true
  >
auto makeLazy(A a) 
{ return [a = std::move(a)]() mutable -> A { return std::move(a); }; }

template<class A
  , std::enable_if_t<std::is_invocable_v<A>, bool> = true
  >
auto makeLazy(A a) -> A
{ return a; }

template<class... Args>
static auto __ifElseIter(Args... args)
{ 
  auto tuple = std::tie(args...);
  constexpr unsigned total = TL::Size<TL::List<Args...>>::val;
  auto tupleGetApplied = [&](auto N) -> decltype(auto) { return std::get<N.value>(tuple)(); };
  static_assert(total % 2 == 1);

  using Out = TL::Into<Coproduct, 
        TL::Map<
          TL::Closure<decltype(tupleGetApplied)>,
          TL::Concat<
            TL::Map<Plus<1>, TL::Map<Times<2>, TL::Range<0, total / 2>>>,
            TL::List<Constant<total - 1>>
          >>>;

  Option<Out> out;
  TL::Range<0, total / 2>::forEach([&](auto token){
      constexpr unsigned i = TL::TokenType<decltype(token)>::value;
      if (out.isNone()) {
        if (tupleGetApplied(Constant<2*i>{})) {
          out = some(Out::template variant<i>(tupleGetApplied(Constant<2*i + 1>{})));
        }
      }
  });
  return coproductIter(out ? *out : Out::template variant<total/2>(tupleGetApplied(Constant<total - 1>{})));
}

template<class... Args>
static auto ifElseIter(Args... args) 
{ return __ifElseIter(makeLazy(std::move(args))...); }

#define ifElseIter2(i0, t0, e) ifElseIter(LAZY(i0), LAZY(t0), LAZY(e))
#define ifElseIter3(i0, t0, i1, t1, e) ifElseIter(LAZY(i0), LAZY(t0), LAZY(i1), LAZY(t1), LAZY(e))
#define ifElseIter4(i0, t0, i1, t1, i2, t2, e) ifElseIter(LAZY(i0), LAZY(t0), LAZY(i1), LAZY(t1), LAZY(i2), LAZY(t2), LAZY(e))

template<class IfIterCons>
static auto ifIter(bool cond, IfIterCons ifCons) 
{ return iterTraits(someIf(cond, std::move(ifCons)).intoIter()).flatten(); }

template<class T>
struct EmptyIter
{
  DECL_ELEMENT_TYPE(T);
  bool hasNext() { return false; }
  T next() { ASSERTION_VIOLATION }
  unsigned size() { return 0; }
  bool knowsSize() { return true; }
};

/** If Pointer references a type that can be used as an iterator, then this wrapper type makes Pointer an Iterator.
 * This is useful for example if you want to have IterTraits methods for some `shared_pointer<I> ptr` where `I` is an iterator.
 */
template<class Pointer>
class IterPointer
{
  Pointer _p;
public:
  IterPointer(Pointer p) : _p(std::move(p)) {}
  DEFAULT_CONSTRUCTORS(IterPointer);
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(std::remove_reference_t<decltype(*_p)>));
  bool hasNext()       { return (*_p).hasNext(); }
  bool hasNext() const { return (*_p).hasNext(); }
  OWN_ELEMENT_TYPE next() { return (*_p).next(); }
  unsigned size() { return (*_p).size(); }
  bool knowsSize() { return (*_p).knowsSize(); }
};

template<class Iter>
class IterTraits
{
  Iter _iter;
public:
  DECL_ELEMENT_TYPE(ELEMENT_TYPE(Iter));
  DEFAULT_CONSTRUCTORS(IterTraits)
  using Elem = ELEMENT_TYPE(Iter);

  explicit IterTraits(Iter iter) : _iter(std::move(iter)) {}

  Elem next() 
  { return move_if_value<Elem>(_iter.next()); }

  bool hasNext() 
  { return _iter.hasNext(); }

  bool knowsSize()       { return _iter.knowsSize(); }
  bool knowsSize() const { return _iter.knowsSize(); }

  auto size()       { return _iter.size(); }
  auto size() const { return _iter.size(); }

  Option<Elem> tryNext() 
  { 
    return _iter.hasNext() 
        ? Option<Elem>(move_if_value<Elem>(_iter.next()))
        : Option<Elem>();
  }


  template<class P>
  bool any(P f) 
  {
    while (hasNext()) {
      if (f(next())) return true;
    }
    return false;
  }

  template<class P>
  bool all(P f) 
  {
    while (hasNext()) {
      if (!f(next())) return false;
    }
    return true;
  }

  template<class F>
  void forEach(F f) 
  {
    while (hasNext()) {
      f(next());
    }
  }

  template<class P>
  Option<Elem> find(P p) 
  {
    while (hasNext()) {
      Elem x = next();
      if (p(x)) {
        return some<Elem>(x);
      }
    }
    return none<Elem>();
  }

  template<class P>
  Option<unsigned> findPosition(P p) 
  {
    unsigned i = 0;
    while (hasNext()) {
      Elem x = next();
      if (p(x)) {
        return some<unsigned>(i);
      }
      i++;
    }
    return none<unsigned>();
  }

  template<class... Is>
  auto concat(Is... is) 
  { return concatIters(std::move(_iter), std::move(is)...); }

  template<class F>
  IterTraits<MappingIterator<Iter, F>> map(F f)
  { return iterTraits(getMappingIterator<Iter, F>(std::move(_iter), std::move(f))); }

  auto eval()
  { return map([](auto f){ return f(); }); }

  template<class F>
  auto inspect(F f)
  { return map([f = std::move(f)](auto x) { f(x); return x; }); }

  template<class F>
  IterTraits<FilteredIterator<Iter, F>> filter(F f)
  { return iterTraits(getFilteredIterator<Iter, F>(std::move(_iter), std::move(f))); }

  template<class F>
  IterTraits<FilterMapIter<Iter, F>> filterMap(F f)
  { return iterTraits(FilterMapIter<Iter, F>(std::move(_iter), std::move(f))); }

  template<class F>
  IterTraits<FlatMapIter<Iter, F>> flatMap(F f)
  { return iterTraits(getFlattenedIterator(getMappingIterator(std::move(_iter), std::move(f)))); }

  auto flatten()
  { return iterTraits(getFlattenedIterator(std::move(_iter))); }

  template<class Pred>
  auto takeWhile(Pred p)
  { return iterTraits(TakeWhileIter<Iter, Pred>(std::move(_iter), std::move(p))); }

  auto unique()
  { 
    Map<OWN_ELEMENT_TYPE, std::tuple<>> found;
    return iterTraits(std::move(*this)
        .filterMap([found = std::move(found)](OWN_ELEMENT_TYPE next) mutable {
          if (found.tryGet(next).isSome()) {
            return Option<OWN_ELEMENT_TYPE>();
          } else {
            found.insert(next, std::make_tuple());
            return Option<OWN_ELEMENT_TYPE>(std::move(next));
          }
        })); 
  }

  auto persistent()
  { 
    auto stack = collect<Stack>();
    return iterTraits(arrayIter(std::move(stack)));
  }

  /** 
   * returns the first minimal element wrt the function `less` 
   * less takes two arguments of this iterators element type and 
   * returns the wheter the first is smaller than the second.
   * */
  template<class IsLess>
  Option<Elem> minBy(IsLess isLess)
  { 
    if (hasNext()) {
      Elem min = next();
      while (hasNext()) {
        Elem e = next();
        if (isLess(e, min)) {
          min = e;
        }
      }
      return some<Elem>(min);
    } else {
      return none<Elem>();
    }
  }

  unsigned count()
  { 
    unsigned i = 0;
    while (hasNext()) {
      i++;
      next();
    }
    return i;
  }

  Option<Elem> min()
  { return minBy(std::less<Elem>{}); }

  template<class IsLess>
  Option<Elem> maxBy(IsLess isLess)
  { return minBy([isLess = std::move(isLess)](Elem const& l, Elem const& r) { return isLess(r,l); }); }

  Option<Elem> max()
  { return maxBy(std::less<Elem>{}); }

  auto timeTraced(const char* name)
  { return iterTraits(TIME_TRACE_ITER(name, std::move(_iter))); }

  template<class Result, class F>
  auto fold(Result init, F f) -> Result
  { 
    Result accum = std::move(init);
    while (hasNext()) {
      accum = f(std::move(accum), next());
    }
    return accum;
  }

  template<class F> 
  auto fold(F fun) -> Option<Elem>
  { return someIf(hasNext(), [&]() { return fold(next(), std::move(fun)); }); }

  template<class OtherIter>
  auto zip(OtherIter other)
  { return map([other = std::move(other)](Elem x) mutable { return std::make_pair(std::move(x), other.next()); }); }

  auto zipWithIndex()
  { return map([idx = 0](Elem x) mutable { return std::make_pair(std::move(x), idx++); }); }

  template<class Val>
  auto store(Val v) 
  { return map([v = std::move(v)](Elem x) -> Elem { return x; }); }

  auto reverse() &&
  { return iterTraits(std::move(_iter).reverse()); }

  auto sum()
  { return fold([](Elem l, Elem r) { return l + r; }) || []() { return Elem(0); }; }

  template<unsigned n>
  auto collectTuple()
  { return TL::Range<0, n>::toTuple([&](auto i) { return next(); }); }
 
  template<class Container>
  auto collect()
  { return Container::fromIterator(*this); }
  
  template<template<class> class Container>
  Container<Elem> collect()
  { 
    return Container<Elem>::fromIterator(std::move(*this)); 
  }

  IterTraits clone() 
  { return *this; }
  
  /** This class is to be used in the context of a for (auto x : ...) loop only. */
  class StlIter 
  {
    Option<IterTraits&> _iter; // <- nothing here encodes that this == end()
    Option<Elem>  _cur;

  public:
    StlIter(IterTraits& iter)  : _iter(Option<IterTraits&>(iter)), _cur(iter.tryNext()) {}
    StlIter()  : _iter(), _cur() {}

    void operator++() 
    { _cur = _iter.unwrap().tryNext(); }

    Elem operator*() 
    { return move_if_value<Elem>(_cur.unwrap()); } 

    friend bool operator!=(StlIter const& lhs, StlIter const& rhs) 
    { return !(lhs == rhs); }

    friend bool operator==(StlIter const& lhs, StlIter const& rhs) 
    { 
      ASS(rhs._iter.isNone()); 
      ASS(lhs._iter.isSome()); 
      return lhs._cur.isNone(); 
    }
  };

public:
  StlIter begin() { return StlIter(*this); }
  StlIter end() { return StlIter(); }
};

template<class Iter> 
IterTraits<Iter> iterTraits(Iter i) { return IterTraits<Iter>(std::move(i)); }

static const auto range = [](auto from, auto to) 
  { return iterTraits(getRangeIterator<decltype(to)>(decltype(to)(from), to)); };

template<class I1>
static auto concatIters(I1 i1) 
{ return iterTraits(std::move(i1)); }

template<class I1, class I2, class... Is>
static auto concatIters(I1 i1, I2 i2, Is... is) 
{ return iterTraits(CatIterator<I1, decltype(concatIters(std::move(i2), std::move(is)...))>(std::move(i1), concatIters(std::move(i2), std::move(is)...))); }

template<class I1, class I2, class Cmp>
auto iterSortedDiff(I1 i1, I2 i2, Cmp cmp) 
{ return iterTraits(SortedIterDiff<I1,I2, Cmp>(std::move(i1), std::move(i2), std::move(cmp))); }

template<class I1, class I2>
auto iterSortedDiff(I1 i1, I2 i2) 
{ return iterSortedDiff(std::move(i1), std::move(i2), [&](auto& l, auto& r) { return l == r ? 0 : l < r ? -1 : 1;  }); }


template<class P>
auto iterPointer(P p) { return iterTraits(IterPointer<P>(std::move(p))); }

///@}

template<class Iterator>
auto dropElementType(Iterator iter) 
{ return iterTraits(std::move(iter)).map([](auto _) { return std::make_tuple(); }); }

template<class Array, class Size>
auto arrayIter(Array const& a, Size s) { return range(0, s).map([&](auto i) -> decltype(auto) { return a[i]; }); }

template<class Array, class Size>
auto arrayIter(Array      & a, Size s) { return range(0, s).map([&](auto i) -> decltype(auto) { return a[i]; }); }

template<class Array, class Size>
auto arrayIter(Array     && a, Size s) { return range(0, s).map([a = std::move(a)](auto i) { return std::move(a[i]); }); }

template<class Array> auto arrayIter(Array const& a) { return arrayIter(          a , a.size()); }
template<class Array> auto arrayIter(Array     && a) { return arrayIter(std::move(a), a.size()); }
template<class Array> auto arrayIter(Array      & a) { return arrayIter(          a , a.size()); }

template<class Item, class... Items>
auto iterItems(Item item, Items... items) 
{ return arrayIter(Stack<Item>{std::move(item), std::move(items)...}); }

template<class Item>
auto iterItems()
{ return arrayIter(Stack<Item>{}); }

} // namespace Lib

template <typename Iterator>
class STLIterator
{
  private:
    Iterator begin;
    Iterator end;

  public:
    using value_type = typename Iterator::value_type;
    DECL_ELEMENT_TYPE(value_type);

    STLIterator(Iterator begin, Iterator end)
      : begin(begin), end(end)
    { }

    bool hasNext() {
      return begin != end;
    }

    value_type next() {
      value_type x = *begin;
      ++begin;
      return x;
    }
};

template <typename Iterator>
STLIterator<Iterator> getSTLIterator(Iterator begin, Iterator end)
{
  return STLIterator<Iterator>(begin, end);
}

/**
 * Return iterator that stores values of @b it in its constructor,
 * and then yields them in the same order
 *
 * After the call to this function, the iterator @b it and any resources
 * it holds may be released, since the elements are stored independently
 * of it.
 *
 * @see PersistentIterator
 */
template<class Inner>
auto getPersistentIterator(Inner it)
{ return pvi(arrayIter(iterTraits(it).template collect<Stack>())); }

template<class Iter>
class IterContOps {
  Iter const _iter;

public:
  IterContOps(Iter iter) : _iter(std::move(iter)) {}

  auto defaultHash() const { return DefaultHash::hashIter(Iter(_iter).map([](ELEMENT_TYPE(Iter) x) -> unsigned { return DefaultHash::hash(x); })); }
  auto defaultHash2() const { return DefaultHash::hashIter(Iter(_iter).map([](ELEMENT_TYPE(Iter) x) -> unsigned { return DefaultHash2::hash(x); })); }

  static int cmp(IterContOps const& lhs, IterContOps const& rhs) {
    auto l = lhs._iter;
    auto r = rhs._iter;
    while (l.hasNext() && r.hasNext()) {
      auto ln = l.next();
      auto rn = r.next();
      if (ln < rn) {
        return -1;
      } else if (rn < ln) {
        return 1;
      }
    }
    return !l.hasNext() ? (r.hasNext() ? -1 : 0) : 1;
  }
  friend bool operator<(IterContOps const& lhs, IterContOps const& rhs) 
  { return cmp(lhs, rhs) < 0; }

  friend bool operator==(IterContOps const& lhs, IterContOps const& rhs)
  { return cmp(lhs, rhs) == 0; }

  friend bool operator!=(IterContOps const& lhs, IterContOps const& rhs) 
  { return !(lhs == rhs); }
};

template<class Iter>
auto iterContOps(Iter iter) { return IterContOps<Iter>(std::move(iter)); }

#endif /* __Metaiterators__ */
