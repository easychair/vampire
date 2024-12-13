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
 * @file List.hpp
 * Implements class List<C> for lists.
 * @since 05/05/1999 Manchester
 * @since 29/05/1999 Manchester: setHead, setTail and conc added
 * @since 06/08/1999 moved from KK to Sat
 * @since 25/09/1999 Uppsala: lots of modifications
 * @since 17/10/1999 Manchester: moved back to kk
 * @since 26/02/2000 Manchester, inserted into library
 * @since 14/05/2000 Dallas, function length added
 * @since 15/05/2000 Dallas, function member added
 * @since 29/05/2000 Manchester, isNonEmpty added
 * @since 01/06/2000 Manchester, second added
 * @since 04/06/2000 Manchester, remove added
 * @since 04/06/2000 Manchester, nth() added
 * @since 10/12/2000 Manchester, class Iterator added
 * @since 17/12/2000 Manchester, Iterator::reset added
 * @since 02/01/2001 Manchester, addLast added
 * @since 07/01/2001 Manchester, pop changed to return a value
 * @since 29/03/2001 Manchester, two small bugs fixed
 * @since 09/05/2001 Manchester, DelIterator introduced
 * @since 05/05/2001 Manchester, deleteNth ()
 * @since 30/05/2001 Manchester, split
 * @since 30/05/2001 Manchester, append
 * @since 04/06/2001 Manchester, DelIterator::restore
 * @since 02/12/2003, Manchester, allocation changed to use Allocator
 */

#ifndef __list__
#define __list__


#include "Forwards.hpp"

#include "Debug/Assertion.hpp"

#include "Allocator.hpp"
#include "VirtualIterator.hpp"

#if VDEBUG

#include <ostream>

#endif

namespace Lib {

/**
 * Class of generic of single-linked lists.
 */
template <class C>
class List
{
public:
  DECL_ELEMENT_TYPE(C);

  /** builds a single-element list */
  explicit List (C head)
    :
    _head(head),
    _tail(0)
  {
  } // List::List

  /** cons */
  List (C head, List* tail)
    :
    _head(head),
    _tail(tail)
  {
  } // List::List

  /** creates an empty list */
  static List* empty ()
  { return 0; }

  /** true if the list is empty */
  static bool isEmpty (const List* l)
  {
    return l == 0;
  } // List::isEmpty

  /** true if the list is non-empty */
  static bool isNonEmpty (const List* l)
  {
    return l != 0;
  } // List::isNonEmpty

  /** return the tail of the list */
  List* tail() const
  {
    return _tail;
  } // List::tail

  /** return a reference to the tail of the list */
  List*& tailReference ()
  {
    return _tail;
  } // List::tailReference

  /** return a reference to the tail of the list */
  List** tailPtr ()
  {
    return &_tail;
  } // List::tailReference

  /** return the head of the list */
  const C head () const
  {
    return _head;
  } // List::head

  C head ()
  {
    return _head;
  } // List::head

  C& headRef ()
  {
    return _head;
  } // List::head

  C* headPtr ()
  {
    return &(this->_head);
  }


  /** return the second element of the list */
  C second () const
  {
    return tail()->head();
  } // List::second

  /** Set the head of the list to hd */
  void setHead (C hd)
  {
    _head = hd;
  } // List::setHead


  /** Set the tail of the list to tl
   *  and return the original tail. */
  List* setTail(List* tl)
  {
    List* res = _tail;
    _tail = tl;
    return res;
  } // list::setTail

  /** Destroy the list */
  static void destroy (List *l)
  {
    if (isEmpty(l)) return;
    List* current = l;

    for (;;) {
      List* next = current->tail();
      delete current;
      if (isEmpty(next)) return;
      current = next;
    }
  } // List::destroy

  /**
   * Destroy the list and apply delete() to every element of the list.
   * Naturally, the list must be a list of pointers.
   * @since 01/04/2006 Bellevue
   */
  static void destroyWithDeletion(List *l)
  {
    if (isEmpty(l)) return;
    List* current = l;

    for (;;) {
      List* next = current->tail();
      delete current->head();
      delete current;
      if (isEmpty(next)) return;
      current = next;
    }
  } // List::destroyWithDeletion

  /** create a shallow copy of the list */
  static List* copy (const List* l)
  {
    if (isEmpty(l)) return empty();

    List* result = new List(l->head());
    List* previous = result;
    List* rest = l->tail();

    while (! isEmpty(rest)) {
      List* tmp = new List(rest->_head);
      previous->_tail = tmp;
      previous = tmp;
      rest = rest->_tail;
    }

    previous->_tail = empty();
    return result;
  }  // List::copy

  /** appends snd to a copy of fst */
  static List* append (List* fst, List* snd)
  {
    if (isEmpty(fst)) return snd;

    List* result = new List(fst->head());
    List* previous = result;
    List* rest = fst->tail();

    while(isNonEmpty(rest)) {
      List* tmp = new List(rest->head());
      previous->setTail(tmp);
      previous = tmp;
      rest = rest->tail();
    }

    previous->setTail(snd);
    return result;
  }  // List::append

  /** return the list obtained by adding elem as the first element
   *  of this list */
  static List* cons(C elem, List* l)
  {
    return new List(elem, l);
  } // List::cons

  /** return list with one element, the given elem */
  static List* singleton(C elem)
  {
    return new List(elem);
  }

  /** push elem to lst */
  static void push(C elem,List* &lst)
  {
    lst = cons(elem, lst);
  } // List::push

  /**
   * Push all elements from @b it iterator to the list.
   *
   * The first element in the iterator will end up as last in the list
   */
  template<class It>
  static void pushFromIterator(It it, List* &lst)
  {
    while(it.hasNext()) {
      push(it.next(), lst);
    }
  }

  /** pop the first element and return it */
  static C pop(List* &lst)
  {
    ASS_NEQ(lst,0);

    List* tail = lst->tail();
    C result = lst->head();
    delete lst;
    lst = tail;

    return result;
  } // List::pop

  /** destructive list concatenation
   *
   * The @b first list is traversed and the @b second is attached
   * to its end.*/
  static List* concat(List* first,List* second)
  {
    if (first == 0) return second;
    if (second == 0) return first;

    List* current = first;
    for (;;) {
      List* next = current->tail();
      if (! next) {
        current->setTail(second);
        return first;
      }
      current = next;
    }
  } // List::concat

  /** Destructive list reversal */
  static List* reverse(List* l)
  {
    if (isEmpty(l)) return empty();

    List* result = empty();

    while (isNonEmpty(l)) {
      List* tl = l->tail();
      l->setTail(result);
      result = l;
      l = tl;
    }

    return result;
  } // List::reverse

  /** return the length of the list */
  static unsigned length(const List *l)
  {
    unsigned len = 0;

    while (isNonEmpty(l)) {
      len ++;
      l = l->tail();
    }

    return len;
  } // List::length

  /** True if elem is a member of the list, the comparison is made using == */
  static bool member (C elem, const List* l)
  {
    while (isNonEmpty(l)) {
      if (l->head() == elem) return true;
      l = l->tail();
    }

    return false;
  } // List::member

  /** Destructively removes the first occurrence of elem
   * from the list and returns the resulting list. Does nothing
   * if elem is not a member of the list.
   */
  static List* remove (C elem, List* l)
  {
    if (isEmpty(l)) return empty();

    if (l->head() == elem) {
      List* result = l->tail();
      delete l;
      return result;
    }
    if (isEmpty(l->tail())) return l;

    List* current = l;
    List* next = l->tail();

    for (;;) {
      if (next->head() == elem) { // element found
        current->setTail(next->tail());
        delete next;
        return l;
      }
      current = next;
      next = next->tail();
      if (isEmpty(next)) return l;
    }
  } // List::remove

  /** Return the nth element, counting from 0 */
  static C nth(const List *l, unsigned n)
  {
    while (n != 0) {
      ASS_NEQ(l,0);

      l = l->tail();
      n--;
    }

    return l->head();
  } // nth

  /** delete the nth element and return it */
  static C deleteNth(List*& lst, int n)
  {
    // nth element, counting from 0
    ASS (n >= 0);

    C result;
    List* l = lst;
    ASS (isNonEmpty(lst));

    if (n == 0) {
      result = l->head();
      lst = l->tail();
      delete l;
      return result;
    }

    // n != 0
    List* next = l->tail();

    while (--n != 0) {
      l = next;
      next = next->tail();
      ASS (isNonEmpty(next));
    }
    //  now next must be deleted
    result = next->head();
    l->setTail(next->tail());
    delete next;

    return result;
  } // deleteNth

  /** Add  elem as the last element and return the resulting list */
  static List* addLast (List* l, C elem)
  {
    if (! l) return new List (elem);

    // nonempty, trying to find the end
    List* current;
    for (current = l; current->_tail; current = current->_tail) {
    }

    current->setTail(new List(elem));

    return l;
  } // List::addLast

  /** Split the list into two sublists, first of the length n. Return
   *  the first sublist and save the second sublist in the variable rest. */
  static List* split (List* l, int n, List*& rest)
  {
    if (! l) {
      ASS_EQ(n,0);
      rest = empty();
      return empty();
    }

    if (n == 0) {
      rest = empty();
      return l;
    }

    List* nth = l;
    while (--n > 0) {
      ASS_NEQ(nth,0);
      nth = nth->_tail;
    }

    ASS_NEQ(nth,0);
    rest = nth->_tail;
    nth->_tail = empty();
    return l;
  } // List::split

#if VDEBUG
// Only works if called on a List of elements with toString functions
  std::string toString(){
    std::string h = _head->toString();
    if(_tail){
      return h+","+_tail->toString();
    }
    else return h;
  }
#endif

  /** iterator over the list elements */
  class Iterator {
  public:
    USE_ALLOCATOR(List::Iterator);
    
    DECL_ELEMENT_TYPE(C);

    explicit Iterator(const List* l) : _lst (l) {}

    /** return the next element */
    C next()
    {
      ASS_NEQ(_lst,0);

      C result = _lst->head();
      _lst = _lst->tail();
      return result;
    }

    /**
     * Return the element that will be returned by next (therefore hasNexgt()
     * must have returned true), but do not advance to a further element.
     */
    C peekAtNext()
    {
      return _lst->head();
    }

    /** True if there is a next element. */
    bool hasNext() const
    {
      return isNonEmpty(_lst);
    }

    void reset(const List* l) { _lst = l; }

   private:
    /** the rest of the list */
    const List* _lst;
  };

  Iterator iter() const { return Iterator(this); }

  /** iterator over references to list elements */
  class RefIterator {
   public:
     USE_ALLOCATOR(List::RefIterator);
     
     DECL_ELEMENT_TYPE(C&);

    explicit RefIterator(List* l) : _lst (l) {}

    /** return the next element */
    C& next()
    {
      ASS_NEQ(_lst,0);

      C& result = _lst->_head;
      _lst = _lst->tail();
      return result;
    }

    /** True if there is a next element. */
    bool hasNext() const
    { return isNonEmpty(_lst); }

    void reset(List* l) { _lst = l; }

   private:
    /** the rest of the list */
    List* _lst;
  };

  RefIterator iter() { return RefIterator(this); }

  /** Iterator that allows one to delete the current element */
  class DelIterator {
   public:
     USE_ALLOCATOR(List::DelIterator);
     
    DECL_ELEMENT_TYPE(C);
    DelIterator (List*& l)
      :
      _lst(l),
      _prev(0),
      _cur(0)
    {}

    /** Reset the iterator to the beginning of the list */
    void reset()
    {
      _prev = 0;
      _cur = 0;
    } // reset

    /** return the next element */
    C next()
    {
      if (_cur) { // there was an element previously returned by next()
	_prev = _cur;
	_cur = _cur->tail();
	ASS_NEQ(_cur,0);
      }
      else { // result must be the first element of the list
	_cur = _lst;
      }
      return _cur->head();
    }

    /** True if there is a next element */
    bool hasNext()
    {
      if (_cur) { // there was an element previously returned by next()
	return _cur->tail() != 0;
      }
      return isNonEmpty(_lst);
    }

    /** Delete the current element */
    void del()
    {
      // we can only delete the element returned by next
      ASS_NEQ(_cur,0);
      // check that two delete requests in row did not occur
      ASS_NEQ(_cur,_prev);

      if (_cur == _lst) { // the first element must be deleted
	_lst = _lst->tail();
	delete _cur;
	_cur = 0;
	return;
      }

      ASS_NEQ(_prev,0);
      // not the first element
      _prev->setTail(_cur->tail());
      delete _cur;
      _cur = _prev;
    } // List::DelIterator::del

    /**
     * Replace the current element by elem.
     */
    void replace(C elem)
    {
      ASS_NEQ(_cur,0);
      _cur->setHead(elem);
    } // DelIterator::replace

    /**
     * Insert a list of elements before the cursor. If the iteration
     * continues, the new list will not be returned by the iterator.
     * @pre At least one element should have been returned
     *   by a previous call to next() so _cur != null
     * @pre The preceding operation must not have been del()
     * @since 27/12/2007 Manchester
     */
    void insert (List* lst)
    {
      ASS_NEQ(_cur, _prev);
      if (! lst) return;

      List* last = lst;
      // lst is non-empty, find it's last element
      while (last->tail()) {
	last = last->tail();
      }

      if (_prev)
	_prev->setTail(lst);
      else _lst = lst;

      last->setTail(_cur);
      _prev = last;
    } // List::DelIterator::insert

    /**
     * Insert an element before the cursor. If the iteration
     * continues, the new element will not be returned by the iterator.
     * @pre At least one element should have been returned
     *   by a previous call to next() so _cur != null
     * @pre The preceding operation must not have been del()
     * @since 27/12/2007 Manchester
     */
    void insert (C elem)
    {
      ASS_NEQ(_cur, _prev);
      List* lst = new List(elem,_cur);
      if (_prev)
	_prev->setTail(lst);
      else _lst = lst;

      lst->setTail(_cur);
      _prev = lst;
    } // List::DelIterator::insert

  private:
    /** The reference to the list over which the iteration is done */
    List*& _lst;
    /** The element previous to _cur. If _cur is the first element
     *  then _prev=null */
    List* _prev;
    /** _cur is the element returned by the last next,
     * if _cur=null then no next was called */
    List* _cur;
  };
  RefIterator delIter() { return DelIterator(this); }

  template<class Iter>
  static List* fromIterator(Iter iter) {
    List* result = List::empty();
    while (iter.hasNext()) {
      List::push(iter.next(), result);
    }
    return result;
  }

  /**
   * iterator over the list elements
   *
   * @warning All elements of this iterator must be yielded
   * or a memory leak will occur.
   */
  class DestructiveIterator {
  public:
    USE_ALLOCATOR(List::DestructiveIterator);
    
    DECL_ELEMENT_TYPE(C);

    explicit
    DestructiveIterator(List* l)
      : _lst (l)
    {}

    /** return the next element */
    C next()
    {
      ASS_NEQ(_lst,0);
      return List::pop(_lst);
    }

    /** True if there is a next element. */
    bool hasNext() const
    {
      return isNonEmpty(_lst);
    }
   private:
    /** the rest of the list */
    List* _lst;
  };

  // use allocator to (de)allocate objects of this class
  USE_ALLOCATOR(List);

  /**
   * Class that allows to create a list initially by pushing element both at the beginning (pushFront, the usual push)
   * and at the end of it (pushBack, the FIFO style).
   *
   * The interal list is not owned by the FIFO. In a typical use case, the list will be retrieved via list() call
   * and kept after FIFO passes out of scope.
   * @since 06/04/2006 Bellevue
   */
  class FIFO {
  public:
    explicit FIFO() : _first(0), _last(0) {}

    bool empty() const {
      return !_first;
    }

    /* If you only need pushFront, you probably don't need FIFO. */
    void pushFront(C elem)
    {
      _first = cons(elem, _first);
      if (!_last) {
        _last = _first;
      }
    }

    void pushBack(C elem)
    {
      List* newLast = new List(elem);
      if (_last) {
        _last->setTail(newLast);
      } else {
        _first = newLast;
      }
      _last = newLast;
    }

    List* list() const
    {
      return _first;
    }

    /**
     * @brief If last's tail is not null, set it to null.
     *
     * @return Whatever was last's tail pointing to.
     *
     * Note that this will only return non-null,
     * if the list inside was modified outside this objects interface.
     */
    List* clipAtLast() const
    {
      List* beyondLast = List::empty();
      if (_last) {
        std::swap(_last->_tail,beyondLast);
      }
      return beyondLast;
    }

  protected:
    List* _first;
    List* _last;
  }; // class List::FIFO

protected:  // structure

  /** head of the list */
  C _head;
  /** tail of the list */
  List* _tail;
};  // class List

///@addtogroup Iterators
///@{

template<typename T>
typename List<T>::Iterator getContentIterator(List<T>* lst)
{
  return typename List<T>::Iterator(lst);
}

///@}?


///@addtogroup Reflection
///@{

/** see Reflection.hpp */
template<typename T>
struct ElementTypeInfo<List<T>* >
{
  typedef T Type;
};

/** see Reflection.hpp */
template<typename T>
struct IteratorTypeInfo<List<T>* >
{
  typedef typename List<T>::Iterator Type;
};
template<typename T>
struct IteratorTypeInfo<const List<T>*>
{
  typedef typename List<T>::Iterator Type;
};

///@}?

#if VDEBUG

template<typename T>
std::ostream& operator<< (std::ostream& out, const List<T>& lstr )
{
  const List<T>* lst=&lstr;
  out<<'[';

  while(lst) {
    out<<lst->head();
    lst=lst->tail();
    if (lst) out << ", ";
  }

  return out<<']';
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const List<T*>& lstr )
{
  const List<T*>* lst=&lstr;
  out<<'[';

  while(lst) {
    out<<(*lst->head());
    lst=lst->tail();
    if (lst) out << ",\n";
  }

  return out<<']';
}





#endif

}

#endif


