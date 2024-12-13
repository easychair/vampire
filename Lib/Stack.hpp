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
 * @file Stack.hpp
 * Defines a class of flexible-sized stacks
 *
 * @since 04/06/2005 Manchester
 */

#ifndef __Stack__
#define __Stack__

#include <algorithm>
#include <cstdlib>
#include <ostream>

#include "Forwards.hpp"
#include "Lib/TypeList.hpp"

#include "Debug/Assertion.hpp"

#include "Allocator.hpp"
#include "Lib/Reflection.hpp"

namespace std
{
template<typename T>
void swap(Lib::Stack<T>& s1, Lib::Stack<T>& s2);
}

namespace Lib {

/**
 * Class of flexible-size generic stacks.
 * @since 11/03/2006 Bellevue
 * @since 14/03/2006 Bellevue reimplemented in a slightly more efficient way.
 * @since 19/05/2007 Manchester reimplemented back due to errors
 */
template<class C>
class Stack
{
public:
  template<typename U>
  friend void std::swap(Stack<U>&,Stack<U>&);

  class Iterator;
  class ConstIterator;
  class BottomFirstIterator;

  DECL_ELEMENT_TYPE(C);
  DECL_ITERATOR_TYPE(Iterator);

  USE_ALLOCATOR(Stack);


  /**
   * Create a stack having initialCapacity.
   */
  inline
  explicit Stack (size_t initialCapacity=0)
    : _capacity(initialCapacity)
  {
    if(_capacity) {
      void* mem = ALLOC_KNOWN(_capacity*sizeof(C),className());
      _stack = static_cast<C*>(mem);
    }
    else {
      _stack = nullptr;
    }
    _cursor = _stack;
    _end = _stack+_capacity;
  }

  inline
  void reserve(size_t capacity) 
  {
    if (_capacity >= capacity) {
      return;
    }
    C* mem = static_cast<C*>(ALLOC_KNOWN(capacity*sizeof(C),className()));
    if (_stack) {
      for (unsigned i = 0; i < size(); i++) {
        ::new(&mem[i]) C(std::move((*this)[i]));
      }
      DEALLOC_KNOWN(_stack,_capacity*sizeof(C),className());

      _cursor = mem + (_cursor - _stack);
      _capacity = capacity;
      _stack = mem;
      _end = _stack + _capacity;
    } else {
      _stack = mem;
      _cursor = mem;
      _capacity = capacity;
      _end = _stack + _capacity;
    }
  }

  bool keepRecycled() const { return _capacity > 0; }


  Stack(const Stack& s)
   : _capacity(s._capacity)
  {
    if(_capacity) {
      void* mem = ALLOC_KNOWN(_capacity*sizeof(C),className());
      _stack = static_cast<C*>(mem);
    }
    else {
      _stack = 0;
    }
    _cursor = _stack;
    _end = _stack+_capacity;

    loadFromIterator(BottomFirstIterator(const_cast<Stack&>(s)));
  }

  Stack(Stack&& s) noexcept
  {
    _capacity = 0;
    _stack = _cursor = _end = nullptr;

    std::swap(*this,s);
  }

  /** De-allocate the stack
   * @since 13/01/2008 Manchester
   */
  inline ~Stack()
  {
    //The while cycle is completely eliminated by compiler
    //in "-O6 -DVDEBUG=0" mode for types without destructor,
    //so this destructor is constant time.
    C* p=_cursor;
    while(p!=_stack) {
      (--p)->~C();
    }
    if(_stack) {
      DEALLOC_KNOWN(_stack,_capacity*sizeof(C),className());
    }
    else {
      ASS_EQ(_capacity,0);
    }
  }

  Stack& operator=(const Stack& s)
  {
    if(&s == this) {
      return *this;
    }
    reset();
    loadFromIterator(BottomFirstIterator(const_cast<Stack&>(s)));
    return *this;
  }

  Stack& operator=(Stack&& s) noexcept
  {
    std::swap(*this,s);
    return *this;
  }


  /**
   * Put all elements of an iterator onto the stack.
   */
  template<class It>
  void loadFromIterator(It it) {
    // TODO check iterator.size() or iterator.sizeHint()
    while(it.hasNext()) {
      push(it.next());
    }
  }

  /**
   * Create a new stack with the contents of the itererator.
   */
  template<class It>
  static Stack fromIterator(It it) {
    Stack out;
    out.loadFromIterator(std::move(it));
    return out;
  }

  /* a first-in-first-out iterator  */
  BottomFirstIterator iterFifo() const
  { return BottomFirstIterator(*this); }

  /**
   * Return a reference to the n-th element of the stack.
   */
  inline
  C& operator[](size_t n)
  {
    ASS(n >= 0);
    ASS(_stack+n < _cursor);

    return _stack[n];
  } // operator[]

  /** Return a const reference to the n-th element of the stack */
  inline
  const C& operator[](size_t n) const
  {
    ASS(n >= 0);
    ASS(_stack+n < _cursor);

    return _stack[n];
  }

  friend int cmp(const Stack& l, const Stack& r)
  {
    int sdiff = int(l.size()) - int(r.size());
    if(sdiff) return sdiff;
    
    auto i1 = arrayIter(l);
    auto i2 = arrayIter(r);
    while (i1.hasNext()) {
      auto& e1 = i1.next();
      auto& e2 = i2.next();
      if (e1 != e2) {
        if (e1 < e2) return -1;
        if (e1 > e2) return 1;
      }
    }
    ASS(!i2.hasNext())
    return 0;
  }

  friend bool operator< (const Stack& l, const Stack& r) { return cmp(l,r) <  0; }
  friend bool operator<=(const Stack& l, const Stack& r) { return cmp(l,r) <= 0; }
  friend bool operator> (const Stack& l, const Stack& r) { return cmp(l,r) >  0; }
  friend bool operator>=(const Stack& l, const Stack& r) { return cmp(l,r) >= 0; }


  bool operator==(const Stack& o) const
  {
    if(size()!=o.size()) {
      return false;
    }
    size_t sz = size();
    for(size_t i=0; i!=sz; ++i) {
      if((*this)[i]!=o[i]) {
	return false;
      }
    }
    return true;
  }

  bool operator!=(const Stack& o) const
  { return !((*this)==o); }

  /**
   * Return the top of the stack.
   * @since 11/03/2006 Bellevue
   */
  inline
  C& top() const
  {
    ASS(_cursor > _stack);
    ASS(_cursor <= _end);

    return _cursor[-1];
  } // Stack::top()

  /**
   * Return the top but one of the stack.
   */
  inline
  C& scnd() const
  {
    ASS(_cursor > _stack + 1);
    ASS(_cursor <= _end);

    return _cursor[-2];
  } // Stack::top()


  /**
   * Set top to a new value.
   * @since 14/03/2006 Bellevue
   */
  inline
  void setTop(C elem)
  {
    ASS(_cursor > _stack);
    ASS(_cursor <= _end);

    _cursor[-1] = elem;
  } // Stack::top()

  /**
   * True if the stack is empty.
   * @since 11/03/2006 Bellevue
   */
  inline
  bool isEmpty() const
  {
    return _cursor == _stack;
  } // Stack::isEmpty()

  /**
   * True if the stack is non-empty.
   * @since 11/03/2006 Bellevue
   */
  inline
  bool isNonEmpty() const
  {
    return _cursor != _stack;
  } // Stack::nonempty()

  /**
   * Push new element on the stack.
   * @since 11/03/2006 Bellevue
   */
  inline
  void push(C elem)
  {
    if (_cursor == _end) {
      expand();
    }
    ASS(_cursor < _end);
    ::new(_cursor) C(std::move(elem));
    _cursor++;
  } // Stack::push()

  /**
   * Pop the stack and return the popped element.
   * @since 11/03/2006 Bellevue
   */
  inline
  C pop()
  {
    ASS(_cursor > _stack);
    _cursor--;

    C res = std::move(*_cursor);
    _cursor->~C();

    return res;
  } // Stack::pop()


  /** removes consecutive duplicates. instead of the operator== the given predicate is used */
  template<class Equal = std::equal_to<C>>
  void dedup(Equal eq = std::equal_to<C>{})
  { 
    auto& self = *this;
    if (self.size() == 0) return;
    unsigned offs = 0;
    for (unsigned i = 1;  i < self.size(); i++) {
      if (eq(self[offs], self[i])) {
        /* skip */
      } else {
        self[offs++ + 1] = std::move(self[i]);
      }
    }
    self.pop(self.size() - (offs + 1));
  }

  /** like Stack::dedup but moves the content out of `this` and returns the resulting Stack instead of changing the contents of this  */
  template<class Equal = std::equal_to<C>>
  Stack deduped(Equal eq = std::equal_to<C>{})
  { dedup(); return std::move(*this); }

  template<class Less = std::less<C>>
  void sort(Less less = std::less<C>{})
  { std::sort(begin(), end(), std::move(less)); }

  /** like Stack::sort but moves the content out of `this` and returns the resulting Stack instead of changing the contents of this */
  template<class Less = std::less<C>>
  Stack sorted(Less less = std::less<C>{})
  { sort(std::move(less)); return std::move(*this); }

  inline
  void pop(unsigned cnt)
  {
    while (cnt-- != 0) 
      pop();
  } // Stack::pop(unsigned)

  /**
   * If the element @b el is present in the stack, remove it and return
   * true, otherwise return false.
   */
  bool remove(C el)
  {
    Iterator it(*this);
    while(it.hasNext()) {
      if(it.next()==el) {
        it.del();
        return true;
      }
    }
    return false;
  }

  /**
   * removes the element at the given index, replacing it by the last element in the stack and shrinking the stack.
   * constant time operation.
   * returns the removed element.
   */
  C swapRemove(unsigned idx)
  {
    ASS(idx < size())
    ASS(size() > 0)
    std::swap((*this)[idx], (*this)[size() - 1]);
    return pop();
  }

  /**
   * Return the element past the end of the stack, can be used together
   * with begin() for iterating over the elements of the stack.
   * @since 11/03/2006 Bellevue
   */
  inline
  C* end() const
  {
    return _cursor;
  }

  inline
  C* begin() const
  {
    return _stack;
  }

  /** Empties the stack. */
  inline
  void reset()
  {
    C* p=_cursor;
    while(p!=_stack) {
      (--p)->~C();
    }
    _cursor = _stack;
  }

  void init(std::initializer_list<C> elems)
  {
    reserve(elems.size());
    for (auto& x : elems) {
      push(std::move(x));
    }
  }
private:
  void __pushMany() {}

  template<class... As>
  void __pushMany(C item, As... rest) 
  { 
    push(std::move(item)); 
    __pushMany(std::move(rest)...);
  }


public:

  template<class... As>
  void pushMany(As... items) 
  { 
    reserve(size() + TypeList::Size<TypeList::List<As...>>::val);
    __pushMany(std::move(items)...);
  }


  /** Sets the length of the stack to @b len
   *  @since 27/12/2007 Manchester */
  inline
  void truncate(size_t len)
  {
    ASS_LE(len,length());
    C* p=_stack+len;
    while(p!=_cursor) {
      (p++)->~C();
    }
    _cursor = _stack+len;
  } // truncate

  /** Return the number of elements in the stack, same as size() */
  inline
  size_t length() const
  { return _cursor - _stack; }

  /** Return the number of elements in the stack, same as length() */
  inline
  size_t size() const
  { return _cursor - _stack; }

  bool find(const C& el) const
  {
    Iterator it(const_cast<Stack&>(*this));
    while(it.hasNext()) {
      if(it.next()==el) {
	return true;
      }
    }
    return false;
  }

  friend class RefIterator;

  /** Iterator iterates over the elements of a stack and can
   *  delete elements from the stack.
   *  @warning After deletion the order of elements in the stack
   *           may change
   *  @warning The contents of the stack should not be changed by
   *           other operations when a stack is traversed using an
   *           iterator
   *  @since 13/02/2008 Manchester
   */
  class RefIterator {
  public:
    DECL_ELEMENT_TYPE(C&);
    /** create an iterator for @b s */
    inline
    explicit RefIterator (Stack& s)
      : _pointer(s._cursor),
	_stack(s)
#if VDEBUG
      , _last(0)
#endif
    {
    }

    /** true if there exists the next element */
    inline
#if VDEBUG
    bool hasNext()
#else
    bool hasNext() const
#endif
    {
#if VDEBUG
      _last = 2;
#endif

      return _pointer != _stack._stack;
    }

    /** return the next element */
    inline
    C& next()
    {
      ASS(_pointer > _stack._stack);
      ASS(_last == 2);
#if VDEBUG
      _last = 1;
#endif

      _pointer--;
      return *_pointer;
    }

    /** Delete the last element returned by next() */
    inline
    void del()
    {
      ASS(_pointer < _stack._cursor);
      ASS(_pointer >= _stack._stack);
      ASS(_last == 1);
#if VDEBUG
      _last = 3;
#endif

      *_pointer = _stack.pop();
    }

    /** Replace the last element returned by next() */
    inline
    void replace(C val)
    {
      ASS(_pointer < _stack._cursor);
      ASS(_pointer >= _stack._stack);
      ASS(_last == 1);

      *_pointer = val;
    }
  private:
    /** pointer to the stack element returned by next() */
    C* _pointer;
    /** stack over which we iterate */
    Stack& _stack;
#if VDEBUG
    /** last operation: 0(none), 1(next), 2(hasNext), 3(del) */
    int _last;
#endif
  };

  class Iterator : public RefIterator {
  public:
    Iterator(Stack & s) : RefIterator(s) {}
    DECL_ELEMENT_TYPE(C);
    C next() { return RefIterator::next(); }
  };

  class ConstRefIterator {
  public:
    DECL_ELEMENT_TYPE(C const&);
    /** create an iterator for @b s */
    inline
    explicit ConstRefIterator (const Stack& s)
      : _pointer(s._cursor),
	_stack(s)
    {
    }

    /** true if there exists the next element */
    inline
    bool hasNext() const
    {
      return _pointer != _stack._stack;
    }

    /** return the next element */
    inline
    C const& next()
    {
      ASS(_pointer > _stack._stack);
      _pointer--;
      return *_pointer;
    }

  private:
    /** pointer to the stack element returned by next() */
    C* _pointer;
    /** stack over which we iterate */
    const Stack& _stack;
  };

  class ConstIterator : public ConstRefIterator {
  public:
    ConstIterator(Stack const& s) : ConstRefIterator(s) {}
    DECL_ELEMENT_TYPE(C);
    C next() { return ConstRefIterator::next(); }
  };

  ConstIterator iterCloned() const&
  { return ConstIterator(*this); }

  RefIterator iter() &
  { return RefIterator(*this); }

  ConstRefIterator iter() const&
  { return ConstRefIterator(*this); }



  typedef Iterator DelIterator;
  typedef ConstIterator TopFirstIterator;

  /**
   * An iterator object that for stack @b s first yields element s[0]
   * and the element s.top() is last.
   */
  class BottomFirstIterator {
  public:
    DECL_ELEMENT_TYPE(C);
    /** create an iterator for @b s */
    inline
    explicit BottomFirstIterator (const Stack& s)
      : _pointer(s._stack),
	_afterLast(s._cursor)
    {
    }

    /** true if there exists the next element */
    inline
    bool hasNext() const
    {
      ASS_LE(_pointer, _afterLast);
      return _pointer != _afterLast;
    }

    /** return the next element */
    inline
    const C& next()
    {
      ASS_L(_pointer, _afterLast);
      return *(_pointer++);
    }

  private:
    /** pointer to the stack element returned by next call to @b next() */
    C* _pointer;
    /** pointer to element after the last element on the stack */
    C*  _afterLast;
  };

  /**
   * Iterator iterates over the elements of a stack from s[0] to s.top()
   * and can delete elements from the stack without changing the order of
   * the remaining elements.
   *  @warning The contents of the stack should not be changed by
   *           other operations when a stack is traversed using an
   *           iterator
   */
  class StableDelIterator {
    StableDelIterator(const StableDelIterator&);
    StableDelIterator& operator=(const StableDelIterator&);
  public:
    DECL_ELEMENT_TYPE(C);
    /** create an iterator for @b s */
    inline
    explicit StableDelIterator (Stack& s)
      : _reader(s._stack),
        _writer(s._stack),
	_stack(s)
#if VDEBUG
      , _last(0)
#endif
    {
    }

    ~StableDelIterator() {
      if(_reader!=_writer) {
	//if we deleted something, we must go through the rest of the stack
	//to shift the remaining elements
	while(hasNext()) {
	  next();
	}
      }
    }

    /** true if there exists the next element */
    inline
    bool hasNext()
    {
#if VDEBUG
      _last = 2;
#endif

      if(_reader==_stack._cursor) {
	if(_reader!=_writer) {
	  _stack._cursor = _writer;
	  _reader = _writer; //this is to handle properly repeated calls to this function
	}
	return false;
      }
      ASS_L(_reader,_stack._cursor);
      return true;
    }

    /** return the next element */
    inline
    C next()
    {
      ASS(_reader < _stack._cursor);
      ASS(_last == 2);
#if VDEBUG
      _last = 1;
#endif
      if(_reader!=_writer) {
	ASS_L(_writer, _reader);
	*_writer = *_reader;
      }
      const C& res = *_reader;
      _reader++;
      _writer++;

      return res;
    }

    /** Delete the last element returned by next() */
    inline
    void del()
    {
      ASS(_writer <= _stack._cursor);
      ASS(_writer >= _stack._stack);
      ASS(_last == 1);
#if VDEBUG
      _last = 3;
#endif

      _writer--;
    }

    /** Replace the last element returned by next() */
    inline
    void replace(C val)
    {
      ASS(_writer < _stack._cursor);
      ASS(_writer >= _stack._stack);
      ASS(_last == 1);

      *_writer = val;
    }
  private:
    /** pointer to the stack element returned by next() */
    C* _reader;
    /** pointer to the stack element returned by next() */
    C* _writer;
    /** stack over which we iterate */
    Stack& _stack;
#if VDEBUG
    /** last operation: 0(none), 1(next), 2(hasNext), 3(del) */
    mutable int _last;
#endif
  };

protected:
  /** Capacity of the stack */
  size_t _capacity;
  /** the stack itself */
  C* _stack;
  /** the cursor, points at the element after the top of the stack */
  C* _cursor;
  /** points to after the last possible value for _cursor */
  C* _end;

  /**
   * Expand the stack. Note: the function heavily uses
   * the fact that the expansion happens exactly when _cursor=_end
   * @since 11/03/2006 Redmond
   */
  void expand ()
  {
    ASS(_cursor == _end);

    size_t newCapacity = _capacity ? (2 * _capacity) : 8;

    // allocate new stack and copy old stack's content to the new place
    void* mem = ALLOC_KNOWN(newCapacity*sizeof(C),className());

    C* newStack = static_cast<C*>(mem);
    if(_capacity) {
      for (size_t i = 0; i<_capacity; i++) {
        ::new(newStack+i) C(std::move(_stack[i]));
        _stack[i].~C();
      }
      // deallocate the old stack
      DEALLOC_KNOWN(_stack,_capacity*sizeof(C),className());
    }

    _stack = newStack;
    _cursor = _stack + _capacity;
    _end = _stack + newCapacity;
    _capacity = newCapacity;
  } // Stack::expand

public:


  friend std::ostream& operator<<(std::ostream& out, const Stack<C>& s) {
    out << "[";
    auto iter = s.begin();
    if(iter != s.end()) {
      out << " " << *iter++;
      while (iter != s.end()) {
        out << ", " << *iter++;
      }
    }
    out << " ]";
    return out;
  }

  Stack(std::initializer_list<C> cont)
   : Stack(cont.size())
  {
    for (auto const& x : cont) {
      push(x);
    }
  }

};

} // namespace Lib

namespace std
{

template<typename T>
void swap(Lib::Stack<T>& s1, Lib::Stack<T>& s2)
{
  using std::swap;//ADL
  swap(s1._capacity, s2._capacity);
  swap(s1._cursor, s2._cursor);
  swap(s1._end, s2._end);
  swap(s1._stack, s2._stack);
}

}

#endif
