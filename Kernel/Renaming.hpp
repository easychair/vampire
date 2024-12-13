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
 * @file Renaming.hpp
 * Defines class Renaming
 *
 */

#ifndef __Renaming__
#define __Renaming__

#include "Lib/DHMap.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Lib/Metaiterators.hpp"
#include "Kernel/TypedTermList.hpp"

#include "Term.hpp"

namespace Kernel {

using namespace Lib;

class Renaming {
public:
  Renaming() :
    _nextVar(0), _identity(true) {
  }

  /*
   * Construct a renaming with names starting at firstVar
   */
  Renaming(unsigned firstVar) :
    _nextVar(firstVar), _identity(true) {
  }

  void reset()
  {
    _data.reset();
    _nextVar = 0;
    _identity = true;
  }
  bool keepRecycled() const { return _data.keepRecycled() > 0; }

  unsigned getOrBind(unsigned v)
  {
    unsigned res;
    if (_data.findOrInsert(v, res, _nextVar)) {
      _nextVar++;
      if(v!=res) {
	_identity = false;
      }
    }
    return res;
  }
  unsigned get(unsigned v) const
  { return _data.get(v); }
  bool contains(unsigned v)
  { return _data.find(v); }

  Literal* apply(Literal* l);
  Term* apply(Term* l);
  TermList apply(TermList l);
  bool identity() const;

  void normalizeVariables(const Literal* t);
  void normalizeVariables(const Term* t);
  void normalizeVariables(TermList t);
  template<class A, class B>
  void normalizeVariables(Coproduct<A, B> t) 
  { return t.apply([&](auto& t){ return normalizeVariables(t); }); }
  void normalizeVariables(TypedTermList t)
  { normalizeVariables(TermList(t)); normalizeVariables(t.sort()); }
  void makeInverse(const Renaming& orig);
  unsigned nextVar() const
  { return _nextVar; }


  template<class A, class B>
  static Coproduct<A,B> normalize(Coproduct<A, B> t)
  { return t.apply([&](auto& t){ return Coproduct<A,B>(normalize(t)); }); }

  static Literal* normalize(Literal* l);
  static TypedTermList normalize(TypedTermList l);
  static Term* normalize(Term* t);
  static TermList normalize(TermList t);
  friend std::ostream& operator<<(std::ostream& out, Renaming const& self)
  { return out << self._data; }

#if VDEBUG
  void assertValid() const;
  std::string toString() const;
#endif
private:
  class Applicator
  {
  public:
    Applicator(Renaming* parent) : _parent(parent) {}
    TermList apply(unsigned var)
    { return TermList(_parent->getOrBind(var), false); }
  private:
    Renaming* _parent;
  };

  typedef DHMap<unsigned, unsigned, IdentityHash, DefaultHash> VariableMap;
  VariableMap _data;
  unsigned _nextVar;
  bool _identity;
public:
  typedef VariableMap::Item Item;
  VirtualIterator<Item> items() const { return _data.items(); }

};

};

#endif /*__Renaming__*/

