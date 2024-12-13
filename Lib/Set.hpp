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
 * @file Set.hpp
 * Defines class Set<Type> of arbitrary sets, implemented in the same way
 * as Map.
 */

#ifndef __Set__
#define __Set__

#include <cstdlib>

#include "Forwards.hpp"

#include "Allocator.hpp"
#include "Hash.hpp"
#include "Reflection.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/Output.hpp"

namespace std {
template<typename T>
void swap(Lib::Set<T>& s1, Lib::Set<T>& s2);
}


namespace Lib {

/**
 * Defines class Set<Val> of arbitrary sets, implemented in the same way
 * as Map. Values are compared using Hash::equals.
 *
 * As defined in Forwards, Hash defaults to Lib::Hash
 * So, if you want to use default hash then either add it to Lib::Hash
 * or provide something in place of Hash
 */
template <typename Val,class Hash>
class Set
{
protected:
  class Cell
  {
  public:
    /** Create a new cell */
    inline Cell ()
      : code(0)
    {
    } // Set::Cell::Cell

    /** True if cell is empty */
    inline bool empty() const
    {
      return code == 0;
    } // Set::Cell::empty

    /** True if contains a deleted element */
    inline bool deleted() const
    {
      return code == 1;
    } // Set::Cell::deleted

    /** True if contains an element */
    inline bool occupied() const
    {
      return code > 1;
    } // Set::Cell::occupied

    /** declared but not defined, to prevent on-heap allocation */
    void* operator new (size_t);

    /** Hash code of the value, 0 and 1 are reserved for occupied and deleted elements */
    unsigned code;
    /** The value in this cell (if any) */
    Val value;
  }; // class Set::Cell

public:
  // use allocator to (de)allocate objects of this class
  USE_ALLOCATOR(Set);

  /** Create a new Set */
  Set()
    : _capacity(0),
      _nonemptyCells(0),
      _size(0),
      _entries(0)
  {
    expand();
  } // Set::Set

  template<typename U>
  friend void std::swap(Set<U>& lhs, Set<U>& rhs);
  

  Set(Set&& other) : Set()
  { std::swap(other, *this); }

  /** Deallocate the set */
  inline ~Set ()
  {
    if (_entries) {
      array_delete(_entries,_capacity);
      DEALLOC_KNOWN(_entries,_capacity*sizeof(Cell),"Set::Cell");
    }
  } // Set::~Set

  /**
   * If the set contains value equal to @b key, return true,
   * and assign the value to @b result
   *
   * Hash class has to contain methods
   * Hash::hash(Key)
   * Hash::equals(Val,Key)
   */
  template<typename Key>
  bool find(Key key, Val& result) const
  {
    unsigned code = Hash::hash(key);
    if (code < 2) {
      code = 2;
    }
    for (Cell* cell = firstCellForCode(code);
	 ! cell->empty();
	 cell = nextCell(cell)) {
      if (cell->deleted()) {
	continue;
      }
      if (cell->code == code &&
	  Hash::equals(cell->value,key)) {
	result=cell->value;
	return true;
      }
    }
    return false;
  } // Set::find

  /**
   * True if the set contains @b val.
   * @since 29/09/2002 Manchester
   */
  bool contains (Val val) const
  {
    unsigned code = Hash::hash(val);
    if (code < 2) {
      code = 2;
    }
    for (Cell* cell = firstCellForCode(code);
	 ! cell->empty();
	 cell = nextCell(cell)) {
      if (cell->deleted()) {
	continue;
      }
      if (cell->code == code &&
	  Hash::equals(cell->value,val)) {
	return true;
      }
    }
    return false;
  } // Set::contains

  /**
   * Checks whether a value with a given hashCode is in the map. 
   * If the value is not present, a new value will be inserted. The new value will be 
   * created using the closure `create`. If a new value has been inserted the bool `inserted` will be set to true, or to false otherwise.
   * When checking whether the is already in the map the closure `isCorrectVal` is used to 
   * compare the value in the map to the one to be inserted. This funciton can be used 
   * in order to avoid allocating a new value when it is already present in the map.
   * a pseudo-code use case for this:
   *
   * ```
   * Set<Stack<char>> set;
   * ...
   * set.rawFindOrInsert(
   *   [](){ return Stack<char>{'a','c'}; },  // <- allocates new stack
   *   computeHashCode({'a','c'}), 
   *   [](auto& stack other) { return other.size() == 2 && other[0] == 'a' && other[1] == 'c' });
   * ```
   */
  template<class Create, class IsCorrectVal>
  Val& rawFindOrInsert(Create create, unsigned hashCode, IsCorrectVal isCorrectVal, bool& inserted)
  {
    inserted = false;

    auto correctHash = [](unsigned hash) { return hash < 2 ? 2 : hash; };
    hashCode = correctHash(hashCode);
    if (_nonemptyCells >= _maxEntries) { // too many entries
      expand();
    }

    Cell* found = 0;
    Cell* cell = firstCellForCode(hashCode);
    while (! cell->empty()) {
      if (cell->deleted()) {
        if (! found) {
          found = cell;
        }
        cell = nextCell(cell);
        continue;
      }
      if (cell->code == hashCode && isCorrectVal(cell->value)) {
        return cell->value;
      }
      cell = nextCell(cell);
    }
    if (found) { // cell deleted
      cell = found;
    }
    else { // cell is empty
      _nonemptyCells++;
    }
    _size++;
    cell->value = create();
    cell->code = hashCode;
    inserted = true;
    ASS_REP(correctHash(Hash::hash(cell->value)) == hashCode, cell->value)
    return cell->value;
  } // Set::insert

  template<class Create, class IsCorrectVal>
  Val& rawFindOrInsert(Create create, unsigned hashCode, IsCorrectVal isCorrectVal)
  { bool b; return rawFindOrInsert(std::move(create), hashCode, std::move(isCorrectVal), b); }


  /**
   * If a value equal to @b val is not contained in the set, insert @b val
   * in the set.
   * Return the value equal to @b val from the set.
   * @since 29/09/2002 Manchester
   * @since 09/12/2006 Manchester, reimplemented
   */
  inline Val insert(const Val val)
  { return insert(val,Hash::hash(val)); } 

  /**
   * Insert a value with a given code in the set.
   * The set must have a sufficient capacity
   * @since 09/12/2006 Manchester, reimplemented
   */
  Val insert(Val val, unsigned code)
  { bool dummy; return rawFindOrInsert([&]() { return std::move(val); },code, [&](auto v) { return Hash::equals(v, val); }, dummy); } // Set::insert

  /** Insert all elements from @b it iterator in the set */
  template<class It>
  void insertFromIterator(It it)
  {
    while(it.hasNext()) {
      insert(it.next());
    }
  }
	
  /** Return the number of (non-deleted) elements */
  inline unsigned size() const
  {
    return _size;
  }

  /**
   * Remove a value from the set. Return true if the value is found
   * @since 23/08/2010 Torrevieja
   */
  bool remove(const Val val)
  {
    unsigned code = Hash::hash(val);
    if (code < 2) {
      code = 2;
    }

    Cell* cell = firstCellForCode(code);
    while (! cell->empty()) {
      if (cell->deleted()) {
	cell = nextCell(cell);
	continue;
      }
      if (cell->code == code &&
	  Hash::equals(cell->value,val)) {
	cell->code = 1; // deleted
	_size--;
	return true;
      }
      cell = nextCell(cell);
    }
    return false;
  } // Set::remove

  /**
   * Make the set empty
   *
   * Unlike reset function for some other containers, the cost
   * of this function is linear in the size of the set.
   */
  void reset()
  {
    Cell* ptr = _entries;
    while(ptr!=_afterLast) {
      ptr->code = 0;
      ptr++;
    }
    _size = 0;
    _nonemptyCells = 0;
  }

  /**
   * Delete all entries.
   * @warning Works only for sets where the value type is a pointer.
   */
  void deleteAll()
  {
    for (int i = _capacity-1;i >= 0;i--) {
      Cell& e = _entries[i];
      if (e.occupied()) {
        delete e.value;
      }
    }
  } // deleteAll

  void outputRaw(std::ostream& out) 
  {
    out << "[";
    for (auto i : range(0, _capacity)) {
      if (_entries[i].occupied()) {
        out << _entries[i].value;
      } else {
        out << "_";
      }
      out << " ";
    }
    out << "]";
  }
private:
  Set(const Set&); //private non-defined copy constructor to prevent copying

  /** the current capacity */
  int _capacity;
  /** the current number of cells */
  int _nonemptyCells;
  /** the current size */
  unsigned _size;
  /** the array of entries */
  Cell* _entries;
  /** the cell after the last one, required since the
   *  array of entries is logically a ring */
  Cell* _afterLast; // cell after the last one
  /** the maximal number of entries for this capacity */
  int _maxEntries;

  /**
   * Expand the set to the next available capacity
   * @throws Exception if map cannot be expanded (that is,
   *         maximal capacity exceeded)
   * @since 29/09/2002 Manchester
   */
  void expand()
  {
    size_t newCapacity = _capacity ? _capacity * 2 : 31;
    Cell* oldEntries = _entries;

    void* mem = ALLOC_KNOWN(newCapacity*sizeof(Cell),"Set::Cell");

    _entries = array_new<Cell>(mem, newCapacity);
    _afterLast = _entries + newCapacity;
    _maxEntries = (int)(newCapacity * 0.8);
    size_t oldCapacity = _capacity;
    _capacity = newCapacity;

    // experiments using (a) random numbers (b) consecutive numbers
    // and 30,000,000 allocations
    // 0.6 :  8.49 7.55
    // 0.7 :  9.19 7.71
    // 0.8 :  9.10 8.44
    // 0.9 :  9.34 7.36 (fewer allocations (21 vs. 22), fewer cache faults)
    // 0.95: 10.21 7.48
    // copy old entries
    Cell* current = oldEntries;
    int remaining = _size;
    _nonemptyCells = 0;
    _size = 0;
    while (remaining != 0) {
      // find first occupied cell
      while (! current->occupied()) {
	current++;
      }
      // now current is occupied
      insert(current->value,current->code);
      current ++;
      remaining --;
    }

    if (oldEntries) {
      array_delete(oldEntries,oldCapacity);
      DEALLOC_KNOWN(oldEntries,oldCapacity*sizeof(Cell),"Set::Cell");
    }
  } // Set::expand

  /**
   * Return the cell next to @b cell.
   * @since 09/12/2006 Manchester
   */
  inline Cell* nextCell(Cell* cell) const
  {
    cell ++;
    // check if the cell is a valid one
    return cell == _afterLast ? _entries : cell;
  } // nextCell

  /**
   * Return the first cell for @b code.
   * @since 09/12/2006 Manchester
   */
  inline Cell* firstCellForCode(unsigned code) const
  {
    return _entries + (code % _capacity);
  } // Set::firstCellForCode

public:
  /**
   * Class to allow iteration over values stored in the set.
   * @since 13/08/2005 Novotel, Moscow
   */
  class Iterator {
  public:
    DECL_ELEMENT_TYPE(Val);

    /** Create a new empty iterator */
    inline Iterator() : _next(0), _last(0) {}

    /** Create a new iterator */
    explicit inline Iterator(const Set& set)
      : _next(set._entries),
	_last(set._afterLast)
    {
    } // Set::Iterator

    /**
     * True if there exists next element
     * @since 13/08/2005 Novotel, Moscow
     */
    bool hasNext()
    {
      while (_next != _last) {
	if (_next->occupied()) {
	  return true;
	}
	_next++;
      }
      return false;
    } // Set::Iterator::hasNext

    /**
     * Return the next value
     * @since 13/08/2005 Novotel, Moscow
     * @warning hasNext() must have been called before
     */
    Val next()
    {
      ASS(_next != _last);
      ASS(_next->occupied());
      Val result = _next->value;
      _next++;
      return result;
    } // Set::Iterator::next

  private:
    /** iterator will look for the next occupied cell starting with this one */
    Cell* _next;
    /** iterator will stop looking for the next cell after reaching this one */
    Cell* _last;
  };
  DECL_ITERATOR_TYPE(Iterator);

  IterTraits<Iterator> iter() const
  { return iterTraits(Iterator(*this)); }

}; // class Set


template<class A, class B>
std::ostream& operator<<(std::ostream& out, Set<A, B> const& self)
{ 
  out << "{ ";
  auto iter = self.iter();
  if (iter.hasNext()) {
    out << iter.next();
    while (iter.hasNext()) 
      out << ", " << iter.next();
  }
  return out << " }";
}

} // namespace Lib

namespace std {

template<typename T>
void swap(Lib::Set<T>& lhs, Lib::Set<T>& rhs)
{
  std::swap(lhs._capacity, rhs._capacity);
  std::swap(lhs._nonemptyCells, rhs._nonemptyCells);
  std::swap(lhs._size, rhs._size);
  std::swap(lhs._entries, rhs._entries);
  std::swap(lhs._afterLast, rhs._afterLast);
  std::swap(lhs._maxEntries, rhs._maxEntries);
}

}

#endif // __Set__

