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
 * @file DHMap.hpp
 * Defines class DHMap<Key,Val,Hash1,Hash2> of maps, implemented as
 * double hashed hashtables.
 */

#ifndef __DHSet__
#define __DHSet__

#include "Forwards.hpp"

#include "DHMap.hpp"

namespace Lib {

/**
 * Class DHSet implements generic sets with values of a class Val.
 *
 * @param Val anything that can be hashed using Hash1 and Hash2,
 *        and compared using ==
 * @param Hash1 class containig the hash function for keys which
 *	  determines position of entry in hashtable when no collision
 *	  occurs.
 * @param Hash2 class containig the hash function for keys which
 *	  will be used when collision occurs. Otherwise it will not be
 *	  enumerated.
 */
template <typename Val, class Hash1, class Hash2>
class DHSet
{
public:
  USE_ALLOCATOR(DHSet);

  /** Empty the DHSet */
  void reset()
  {
    _map.reset();
  }

  DHSet() : _map() {}
  DHSet(DHSet&&) = default;
  DHSet& operator=(DHSet&&) = default;

  bool keepRecycled() const { return _map.keepRecycled(); }

  /**
   *  Return true iff @b val is in the set.
   */
  inline
  bool find(Val const& val) const
  {
    return _map.find(val);
  }

  /**
   *  Return true iff @b val is in the set.
   *
   *  (synomym for the @b find function)
   */
  inline
  bool contains(Val const& val) const
  {
    return find(val);
  }

  /**
   * If the @b val is not in the set, insert it and return true.
   * Otherwise, return false.
   */
  bool insert(Val val)
  {
    return _map.insert(std::move(val), EmptyStruct());
  }


  /**
   * If there is a value stored under the @b key, remove
   * it and return true. Otherwise, return false.
   */
  bool remove(Val const& val)
  {
    return _map.remove(val);
  }

  /** Return mumber of entries stored in this DHMap */
  inline
  unsigned size() const
  {
    return _map.size();
  }

  /** Return true iff there are any entries stored in this DHMap */
  inline
  bool isEmpty() const
  {
    return _map.isEmpty();
  }

  /** Return one arbitrary item that is present in the set */
  Val getOneKey()
  {
    return _map.getOneKey();
  }

  /**
   * Insert all elements of the iterator @b it
   */
  template<class It>
  void loadFromIterator(It it) {
    while(it.hasNext()) {
      insert(it.next());
    }
  }

  /**
   * Remove all elements of the iterator @b it
   *
   * The iterator elements must be present in the set
   */
  template<class It>
  void removeIteratorElements(It it) {
    while(it.hasNext()) {
      ALWAYS(remove(it.next()));
    }
  }

  VirtualIterator<Val> iterator() const
  {
    return _map.domain();
  }

  auto iter() const { return iterator(); }

  friend std::ostream& operator<<(std::ostream& out, DHSet const& self) 
  {
    auto iter = self._map.items();
    out << "{";
    if (iter.hasNext()) {
      out << iter.next().first;
      while (iter.hasNext()) {
        out << ", " << iter.next().first;
      }
    }
    return out << "}";
  }
private:
  /** operator= is private and without a body, because we don't want any. */
  DHSet& operator=(const DHSet& obj);

  typedef DHMap<Val,EmptyStruct,Hash1,Hash2> InnerMap;

  InnerMap _map;

public:
  class Iterator
  {
  public:
    Iterator(const DHSet& parent) : _mit(parent._map) {}

    bool hasNext() { return _mit.hasNext(); }
    Val next() { return _mit.nextKey(); }

  private:
    typename InnerMap::Iterator _mit;
  };
  class DelIterator
  {
  public:
    DelIterator(DHSet& parent) : _mit(parent._map) {}

    bool hasNext() { return _mit.hasNext(); }
    Val next() { return _mit.nextKey(); }
    void del() { _mit.del(); }

  private:
    typename InnerMap::DelIterator _mit;
  };
}; // class DHSet

}

#endif // __DHSet__

