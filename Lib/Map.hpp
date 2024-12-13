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
 * @file Map.hpp
 * Defines class Map<Key,Val> of arbitrary maps and its companion.
 */

#ifndef __Map__
#define __Map__

#include <cstdlib>

#include "Debug/Assertion.hpp"

#include "Allocator.hpp"
#include "Hash.hpp"
#include "Exception.hpp"
#include "Option.hpp"
#include "Metaiterators.hpp"

namespace Lib {

/**
 * Class Map implements generic maps with keys of a class Key
 * and values of a class Value If you implement a map with
 * a new class Key. Hash is the class containing a function
 * hash() mapping keys to unsigned integer values.
 *
 * @param Key a pointer or integral value (e.g., integer or long):
 *        anything that can be hashed to an unsigned integer
 *        and compared using ==
 * @param Val values, can be anything
 * @param Hash class containing the "hash" and "equals" functions for keys
 */
template <typename Key, typename Val,class Hash>
class Map
{
public:

  using HashFn = Hash;
  class Entry
  {
    friend class Map;
  public:
    /** Create a new entry */
    inline Entry ()
      : code(0)
    {
    } // Map::Entry::Entry

    inline ~Entry()
    {
      if (occupied()) {
        key().~Key();
        value().~Val();
      }
    }
    void reset() {
      if (occupied()) {
        key().~Key();
        value().~Val();
      }
      code = 0;
    }

  private:
    /** Create a new entry */
    explicit Entry(Entry const& other)
      : code(0)
    {
      if (other.occupied()) {
        init(Key(other.key()), Val(other.value()), other.code);
      }
    } // Map::Entry::Entry

  private:

    /** True if the cell is occupied */
    inline bool occupied () const
    { return code; } // Map::Entry::occupied

    /** declared but not defined, to prevent on-heap allocation */
    void* operator new (size_t);

    /** Hash code, 0 if not occupied */
    unsigned code;

    /** this wrapper is required in order to leave the storage realy unininitialized, which is 
     * 1) a performance boost, and
     * 2) required in order to make Map work with types that do not have a default-constructor 
     */
    template<class T> using MaybeUninit = typename std::aligned_storage<sizeof(T), alignof(T)>::type;

    /** The key of this cell (if any) */
    MaybeUninit<Key> _key;

    /** The value in this cell (if any) */
    MaybeUninit<Val> _value;

  public:
    /* unwrap the wrapper type */
    Val      & value()      & { ASS(code); return *reinterpret_cast<Val*>(&_value); }
    Val     && value()     && { ASS(code); return std::move(*reinterpret_cast<Val*>(&_value)); }
    Val const& value() const& { ASS(code); return *reinterpret_cast<Val const*>(&_value); }
    Key      &   key()      & { ASS(code); return *reinterpret_cast<Key*>(&_key);   }
    Key     &&   key()     && { ASS(code); return std::move(*reinterpret_cast<Key*>(&_key));   }
    Key const&   key() const& { ASS(code); return *reinterpret_cast<Key const*>(&_key);   }

    friend std::ostream& operator<<(std::ostream& out, Entry const& self) 
    { return self.occupied() ? out << self.key() << " -> " << self.value() : out << "<empty entry>";   } 

  private:

    /** initialize value underlying the wrapper type */
    void init(Key key, Val val, unsigned code)
    {
      ASS_REP(this->code == 0, this->code)
      ASS(code != 0)
      ::new(&_key  ) Key(std::move(key));
      ::new(&_value) Val(std::move(val));
      this->code = code;
    }
  }; // class Map::Entry

  /** Create a new map */
  Map ()
    : _capacity(0),
      _noOfEntries(0),
      _entries(0)
  {
    expand();
  } // Map::Map

  explicit Map (Map const& other)
    : _capacity(other._capacity),
      _noOfEntries(other._noOfEntries),
      _entries((Entry*)ALLOC_KNOWN(sizeof(Entry)*_capacity,"Map<>")),
      _afterLast  (_entries + (other._afterLast - other._entries)),
      _maxEntries (other._maxEntries)
  {

    for (int i = 0; i < _capacity; i++) {
      ::new(&_entries[i]) Entry(other._entries[i]);
    }
  }


  Map (Map && other)
    : _capacity   (other._capacity),
      _noOfEntries(other._noOfEntries),
      _entries    (other._entries),
      _afterLast  (other._afterLast),
      _maxEntries (other._maxEntries)
  {
    other._capacity    = 0;
    other._noOfEntries = 0;
    other._entries     = nullptr;
    other._afterLast   = nullptr;
    other._maxEntries  = 0;
  }

  Map& operator=(Map&& other) {
    _capacity    = other._capacity;
    _noOfEntries = other._noOfEntries;
    _entries     = other._entries;
    _afterLast = other._afterLast;
    _maxEntries = other._maxEntries;

    other._capacity    = 0;
    other._noOfEntries = 0;
    other._entries     = nullptr;
    other._afterLast   = nullptr;
    other._maxEntries  = 0;

    return *this;
  }

  /** Deallocate the map */
  inline ~Map ()
  {
    clear();
  } // Map::~Map

  /**
   * True if there is a value stored under this key.
   * @since 08/08/2008 Manchester
   */
  inline bool find(Key const& key) const
  { return tryGet(key).isSome(); }

  inline unsigned hashCode(Key const& key) const 
  { 
    auto code = Hash::hash(key);
    return code == 0 ? 1 : code;
  }

  /**
   * Find value by the key, and return it if it exists. Return an empty option otherwise
   * 
   * @since 25/08/2020 Manchester
   */
  Option<Val const&> tryGet(Key const& key) const
  {
    using Opt = Option<Val const&>;

    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(),key)) {
        return Opt(entry->value());
      }
    }

    return Opt();
  } // Map::find

  /**
   * Find value by the key. The result is true if a pair with this key
   * is in the map. If such a pair is found then its value is
   * returned in found.
   *
   * @since 29/09/2002 Manchester
   */
  bool find(Key key, Val& found) const
  {
    auto out = tryGet(key);
    if (out.isSome()) {
      found = out.unwrap();
      return true;
    } else {
      return false;
    }
  } // Map::find


  /**
   * Return the value associated with the key if it is present, or a nullptr otherwise.
   * (mutable version)
   *
   * @since 02/06/2020 Manchester
   * @author Johannes Schoisswohl
   */
  Val* getPtr(const Key& key) 
  {
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(),key)) {
        return &entry->value();
      }
    }
    return nullptr;
  }// Map::getPtr

  /**
   * Return the value associated with the key if it is present, or a nullptr otherwise.
   * (immutable version)
   *
   * @since 02/06/2020 Manchester
   * @author Johannes Schoisswohl
   */
  const Val* getPtr(const Key& key) const
  {
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(),key)) {
        return &entry->value();
      }
    }
    return nullptr;
  } // Map::getPtr

  /** returns the number of entries */
  int size() const 
  { return _noOfEntries; }

  /**
   * Return the value by the key. The value must be stored in the
   * map already.
   * @since 26/08/2010 Torrevieja
   */
  Val& get(Key key) const
  {
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); !Hash::equals(entry->key(),key); entry = nextEntry(entry)) {
      ASS(entry->occupied());
    }
    ASS(entry->occupied());
    return entry->value();
  } // Map::get

  /**
   * Return the first entry for @b code.
   * @since 09/12/2006 Manchester
   */
  inline Entry* firstEntryForCode(unsigned code) const
  {
    ASS(_entries)
    return _entries + (code % _capacity);
  } // Map::firstEntryForCode

  /**
   * Return the entry next to @b entry.
   * @since 09/12/2006 Manchester
   */
  inline Entry* nextEntry(Entry* entry) const
  {
    entry ++;
    // check if the entry is a valid one
    return entry == _afterLast ? _entries : entry;
  } // nextEntry

  /**
   * If no value is stored under key @b key, insert pair (key,value) in the map.
   * Return the value stored under @b key.
   * @since 29/09/2002 Manchester
   * @since 09/12/2006 Manchester, reimplemented
   * @since 23/12/2013 Manchester, documentation changed
   * @author Andrei Voronkov
   */
  inline Val& insert(Key key,Val val)
  {
    if (_noOfEntries >= _maxEntries) { // too many entries
      expand();
    }
    auto code = hashCode(key);
    return insert(std::move(key),std::move(val),code);
  } // Map::insert


private:
  /**
   * If no value is stored under key @b key, insert pair (key,value) in the map.
   * Return the value stored under @b key. @b code is the previously computed
   * hash code of the value.
   * The set must have a sufficient capacity
   * @since 09/12/2006 Manchester, reimplemented
   * @since 23/12/2013 Manchester, documentation changed
   * @since 02/06/2020 Manchester, refactor to work with non-copyable types (by Johannes Schoisswohl)
   * @author Andrei Voronkov
   */
  Val& insert(Key&& key, Val&& val,unsigned code)
  {
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(),key)) {
        return entry->value();
      }
    }
    // entry is not occupied
    _noOfEntries++;
    entry->init(std::move(key), std::move(val), code);
    return entry->value();
  } // Map::insert

public:

  /**
   * If no value under key @b key is not contained in the set, insert
   * pair (key,value) in the map. Otherwise replace the value stored
   * under @b key by @b val. Returns true iff the value is replaced
	 * and false if the value is new.
   * @since 29/09/2002 Manchester
   * @since 09/12/2006 Manchester, reimplemented
   */
  bool replaceOrInsert(Key key,Val val)
  {
    if (_noOfEntries >= _maxEntries) { // too many entries
      expand();
    }
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(), key)) {
        entry->value() = std::move(val);
        return true;
      }
    }
    // entry is not occupied
    _noOfEntries++;
    entry->init(std::move(key), std::move(val), code);
    return false;
  } // Map::replaceOrInsert


  /**
   * Replace the value stored under @b key by @b val.
   * (There must be a value stored under @b key).
   * @since 29/09/2002 Manchester
   * @since 09/12/2006 Manchester, reimplemented
   */
  void replace(const Key key,const Val val)
  {
    if (_noOfEntries >= _maxEntries) { // too many entries
      expand();
    }
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(), key)) {
        entry->value() = val;
        return;
      }
    }
#if VDEBUG
    ASSERTION_VIOLATION;
#endif
  } // Map::replace

 
  /**
   * Find the entry with key @b key, or initialize it with the function init otherwise.
   *
   * @b init must have the signature `Val init() {...}`
   */
  template<class InitFn>
  Val& getOrInit(Key key, InitFn init) 
  {
    return updateOrInit(std::move(key), [](Val v) { return std::move(v); }, init);
  } 
 
  /**
    * like `Val& getOrInit(Key key, InitFn init)`, but uses default constructor for initialization.
   */
  Val& getOrInit(Key key) 
  {
    return getOrInit(std::move(key), []() { return Val(); });
  } 


 
  /**
   * Find the entry with key @b key, and update it with UpdateFn, or initialize the value 
   * if it was not present before
   *
   * @b update must have the signature `Val init(Val) {...}`
   * @b init must have the signature `Val init() {...}`
   */
  template<class UpdateFn, class InitFn>
  Val& updateOrInit(Key key, UpdateFn update, InitFn init) 
  {
    if (_noOfEntries >= _maxEntries) { // too many entries
      expand();
    }
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(), key)) {
        ASS_NO_EXCEPT(
          entry->value() = update(std::move(entry->value()));
        )
        return entry->value();
      }
    }
    // entry is not occupied
    _noOfEntries++;
    ASS_NO_EXCEPT(
      entry->init(std::move(key), init(), code);
    )
    return entry->value();
  } 


  /**
   * Assign pointer to value stored under @b key into @b pval.
   * If nothing was previously stored under @b key, initialize
   * the value with @b initial, and return true. Otherwise,
   * return false.
   */
  bool getValuePtr(const Key& key, Val*& pval, const Val& initial)
  {
    if (_noOfEntries >= _maxEntries) { // too many entries
      expand();
    }
    auto code = hashCode(key);
    Entry* entry;
    for (entry = firstEntryForCode(code); entry->occupied(); entry = nextEntry(entry)) {
      if (entry->code == code && Hash::equals(entry->key(), key)) {
        pval = &entry->value();
        return false;
      }
    }
    // entry is not occupied
    _noOfEntries++;
    entry->init(Key(key), Val(initial), code);
    pval = &entry->value();
    return true;
  }
  

  void clear()
  {
    if (_entries) {
      array_delete(_entries, _capacity);
      DEALLOC_KNOWN(_entries,sizeof(Entry)*_capacity,"Map<>");
    }
    _capacity    = 0;
    _noOfEntries = 0;
    _entries     = nullptr;
    _afterLast   = nullptr;
    _maxEntries  = 0;
  }
  
  /**
   * resets every entry in the map keeping the memory of _entries allocated
   */
  void reset()
  {
    for (int i = _capacity-1;i >= 0;i--) {
      _entries[i].reset();
    }
    _noOfEntries = 0;
  } // reset

 
  /**
   * Delete all entries.
   * @since 07/08/2005 Redmond
   * @warning Works only for maps where the value type is a pointer.
   */
  void deleteAll()
  {
    for (int i = _capacity-1;i >= 0;i--) {
      Entry& e = _entries[i];
      if (e.occupied()) {
        delete e.value();
      }
    }
  } // deleteAll

  /** Return the number of elements */
  inline int numberOfElements()
  {
    return _noOfEntries;
  }

  /**
   * Destroy all entries by applying destroy() to them.
   * @since 03/12/2006 Haifa
   * @warning Works only for maps where the value type is a pointer
   *          and having method destroy()
   */
  void destroyAll()
  {
    for (int i = _capacity-1;i >= 0;i--) {
      Entry& e = _entries[i];
      if (e.occupied()) {
        e.value->destroy();
      }
    }
  } // destroyAll

  /** the current capacity */
  int _capacity;
  /** the current number of entries */
  int _noOfEntries;
  /** the array of entries */
  Entry* _entries;
  /** the entry after the last one, required since the
   *  array of entries is logically a ring */
  Entry* _afterLast; // entry after the last one
  /** the maximal number of entries for this capacity */
  int _maxEntries;

  void expand()
  {
    size_t oldCapacity = _capacity;
    _capacity = _capacity ? _capacity * 2 : 32;

    Entry* oldEntries = _entries;

    void* mem = ALLOC_KNOWN(sizeof(Entry)*_capacity,"Map<>");
    _entries = array_new<Entry>(mem, _capacity);

    _afterLast = _entries + _capacity;
    _maxEntries = (int)(_capacity * 0.8);
    // experiments using (a) random numbers (b) consecutive numbers
    // and (1) int->int 20M allocations (2) string->int 10M allocations
    // and 30,000,000 allocations
    // 0.6 : 6.80 4.87 20.8 14.9 32.6 14
    // 0.7 : 6.58 5.61 23.1 15.2 35.2 16.6
    // 0.8 : 6.36 5.77 24.0 15.4 36.0 17.4
    // 0.9 : 7.54 6.04 25.4 15.2 37.0 18.4
    // copy old entries
    Entry* current = oldEntries;
    int remaining = _noOfEntries;
    _noOfEntries = 0;
    while (remaining != 0) {
      // find first occupied entry
      while (! current->occupied()) {
        current ++;
      }
      // now current is occupied
      insert(std::move(current->key()),std::move(current->value()),current->code);
      current ++;
      remaining --;
    }
    if (oldEntries) {
      array_delete(oldEntries, oldCapacity);
      DEALLOC_KNOWN(oldEntries,sizeof(Entry)*oldCapacity,"Map<>");
    }
  } // Map::expand

public:

  /**
   * Class to allow iteration over keys and values stored in the map.
   * @since 13/08/2005 Novotel, Moscow
   */
  class Iterator {
  public:
    DECL_ELEMENT_TYPE(Entry&);

    /** Create a new iterator */
    inline Iterator(const Map& map)
      : _next(map._entries), _last(map._afterLast)
    { } // Map::Iterator

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
    }

    /**
     * Return the next value
     * @since 13/08/2005 Novotel, Moscow
     * @warning hasNext() must have been called before
     */
    Entry& next()
    {
      ASS(_next != _last);
      ASS(_next->occupied());
      Entry& result = *_next;
      _next++;
      return result;
    }

  private:
    /** iterator will look for the next occupied cell starting with this one */
    Entry* _next;
    /** iterator will stop looking for the next cell after reaching this one */
    Entry* _last;
  };


  /**
   * Class to allow iteration over keys and values stored in the map.
   * @since 28/08/2020 Johannes Schoisswohl, Manchester
   */
  class ConstIterator {
  public:
    /** Create a new iterator */
    inline ConstIterator(Map const& map)
      : _next(map._entries), _last(map._afterLast)
    { } // Map::ConstIterator

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
    }

    /**
     * Return the next entry
      * @since 28/08/2020 Johannes Schoisswohl, Manchester
     * @warning hasNext() must have been called before
     */
    Entry const& next()
    {
      ASS(_next != _last);
      ASS(_next->occupied());
      Entry& result = *_next;
      _next++;
      return result;
    }

  private:
    /** iterator will look for the next occupied cell starting with this one */
    Entry* _next;
    /** iterator will stop looking for the next cell after reaching this one */
    Entry* _last;
  };

  bool keepRecycled() const { return size() != 0; }


  Iterator iter() 
  { return Iterator(*this); }

  ConstIterator iter() const
  { return ConstIterator(*this); }

  friend std::ostream& operator<<(std::ostream& out, Map const& self) 
  { 
    out << "{";
    auto iter = self.iter();
    if (iter.hasNext()) {
      out << iter.next();
      while (iter.hasNext()) {
        out << ", " << iter.next();
      }
    }
    return out << "}";
  }
}; // class Map

} // namespace Lib

#endif // __Map__

