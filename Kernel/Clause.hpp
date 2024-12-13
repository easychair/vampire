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
 * @file Clause.hpp
 * Defines class Clause for units consisting of clauses
 *
 * @since 09/05/2007 Manchester
 */

#ifndef __Clause__
#define __Clause__

#include <iosfwd>

#include "Forwards.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/Event.hpp"
#include "Lib/InverseLookup.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/Reflection.hpp"
#include "Lib/Stack.hpp"

#include "Unit.hpp"
#include "Kernel/Inference.hpp"

namespace Kernel {

using namespace Lib;

/**
 * Class to represent clauses.
 * @since 10/05/2007 Manchester
 *
 * When creating a clause object, several things usually need to be done
 * besides calling a constructor:
 * - Fill the Clause with Literals
 * - Increase a relevant counter in the env.statistics object
 */
class Clause
  : public Unit
{
private:
  /** Should never be used, declared just to get rid of compiler warning */
  ~Clause() { ASSERTION_VIOLATION; }
  /** Should never be used, just that compiler requires it */
  void operator delete(void* ptr) { ASSERTION_VIOLATION; }

  template<class VarIt>
  void collectVars2(DHSet<unsigned>& acc);
public:
  DECL_ELEMENT_TYPE(Literal*);

  /** Storage kind */
  enum Store {
    /** passive clause */
    PASSIVE = 0u,
    /** active clause */
    ACTIVE = 1u,
    /** queue of unprocessed clauses */
    UNPROCESSED = 2u,
    /** anything else */
    NONE = 3u,
    /** clause is selected from the passive container
     * and is not added to the active one yet */
    SELECTED = 4u
  };


private:
  Clause(Literal* const* lits, unsigned length, Inference inf);
  void* operator new(size_t,unsigned length);
public:
  void operator delete(void* ptr,unsigned length);

  static Clause* fromArray(Literal*const* lits, unsigned size, Inference inf)
  { return new(size) Clause(lits, size, std::move(inf)); }

  static Clause* fromLiterals(std::initializer_list<Literal*> lits, Inference inf)
  { return fromArray(std::data(lits), lits.size(), std::move(inf)); }

  static Clause* empty(Inference inf)
  { return fromLiterals({}, inf); }


  static Clause* fromStack(const Stack<Literal*>& lits, Inference inf)
  { return new(lits.size()) Clause(lits.begin(), lits.size(), std::move(inf)); }

  template<class Iter>
  static Clause* fromIterator(Iter litit, const Inference& inf)
  {
    static Stack<Literal*> st;
    st.reset();
    st.loadFromIterator(std::move(litit));
    return fromStack(st, inf);
  }

  static Clause* fromClause(Clause* c);

  /**
   * Return the (reference to) the nth literal
   *
   * Positions of literals in the clause are cached in the _literalPositions
   * object. In order to keep it in sync, content of the clause can be changed
   * only right after clause construction (before the first call to the
   * getLiteralPosition method), or during the literal selection (as the
   * _literalPositions object is updated in call to the setSelected method).
   */
  Literal*& operator[] (int n)
  { return _literals[n]; }
  /** Return the (reference to) the nth literal */
  Literal*const& operator[] (int n) const
  { return _literals[n]; }

  /** Return the length (number of literals) */
  unsigned length() const { return _length; }
  /** Alternative name for length to conform with other containers */
  unsigned size() const { return _length; }

  /** Return a pointer to the array of literals.
   * Caller should not manipulate literals, with the exception of
   * clause construction and literal selection. */
  Literal** literals() { return _literals; }

  /** True if the clause is empty */
  bool isEmpty() const { return _length == 0; }

  void destroy();
  void destroyExceptInferenceObject();
  std::string literalsOnlyToString() const;
  std::string toString() const;
  std::string toTPTPString() const;
  std::string toNiceString() const;

  friend std::ostream& operator<<(std::ostream& out, Clause const& self);

  /** Return the clause store */
  Store store() const { return _store; }
  void setStore(Store s);

  /** Return the age */
  unsigned age() const { return inference().age(); }
  /** Set the age to @b a */
  void setAge(unsigned a) { inference().setAge(a); }

  /** Return the number of selected literals */
  unsigned numSelected() const { return _numSelected; }
  /** Mark the first s literals as selected */
  void setSelected(unsigned s)
  {
    ASS(s >= 0);
    ASS(s <= _length);
    _numSelected = s;
    notifyLiteralReorder();
  }

  /** Return the weight = sum of literal weights (usually the number of symbols) */
  unsigned weight() const
  {
    if(!_weight) {
      _weight = computeWeight();
    }
    return _weight;
  }
  unsigned computeWeight() const;

  /**
   * weight used for clause selection
   */
  unsigned weightForClauseSelection(const Shell::Options& opt)
  {
    if(!_weightForClauseSelection) {
      _weightForClauseSelection = computeWeightForClauseSelection(opt);
    }
    return _weightForClauseSelection;
  }
  unsigned computeWeightForClauseSelection(const Shell::Options& opt) const;

  /*
   * single source of truth for computation of weightForClauseSelection
   */
  static unsigned computeWeightForClauseSelection(unsigned w, unsigned splitWeight, unsigned numeralWeight, bool derivedFromGoal, const Shell::Options& opt);

  /** Return the color of a clause */
  Color color() const
  {
    if(static_cast<Color>(_color)==COLOR_INVALID) {
      computeColor();
    }
    return static_cast<Color>(_color);
  }
  void computeColor() const;
  void updateColor(Color c) {
    _color = c;
  }

  bool isExtensionality() const { return _extensionality; }
  bool isTaggedExtensionality() const { return _extensionalityTag; }
  void setExtensionality(bool e) { _extensionality = e; }

  bool isComponent() const { return _component; }
  void setComponent(bool c) { _component = c; }

  bool skip() const;

  unsigned getLiteralPosition(Literal* lit);
  void notifyLiteralReorder();

  bool shouldBeDestroyed();
  void destroyIfUnnecessary();

  void incRefCnt() { _refCnt++; }
  void decRefCnt()
  {
    ASS_G(_refCnt,0);
    _refCnt--;
    destroyIfUnnecessary();
  }

  unsigned getReductionTimestamp() { return _reductionTimestamp; }
  void invalidateMyReductionRecords()
  {
    _reductionTimestamp++;
    if(_reductionTimestamp==0) {
      INVALID_OPERATION("Clause reduction timestamp overflow!");
    }
  }
  bool validReductionRecord(unsigned savedTimestamp) {
    return savedTimestamp == _reductionTimestamp;
  }

  auto getSelectedLiteralIterator() { return arrayIter(*this,numSelected()); }
  auto iterLits()                   { return arrayIter(*this,size()); }
  auto iterLits() const             { return arrayIter(*this,size()); }
  // TODO remove this
  auto getLiteralIterator()         { return arrayIter(*this,size()); }

  bool isGround();
  bool isPropositional();
  bool isHorn();

  VirtualIterator<unsigned> getVariableIterator() const;

  bool contains(Literal* lit);
#if VDEBUG
  void assertValid();
#endif

  SplitSet* splits() const { return _inference.splits(); }
  bool noSplits() const;

  /**
   * set splits
   * in order to keep all splitting-related functionality separate from Saturation,
   * the splits are not set during clause-construction but are added later by the Splitter-class.
   * we depend on the invariant that splits are set only once, and that splits are set before clause-weights are
   * computed and cached (which happens at the first call to weight())
   */
  void setSplits(SplitSet* splits) {
    ASS(_weight == 0);
    _inference.setSplits(splits);
  }

  int getNumActiveSplits() const { return _numActiveSplits; }
  void setNumActiveSplits(int newVal) { _numActiveSplits = newVal; }
  void incNumActiveSplits() { _numActiveSplits++; }
  void decNumActiveSplits() { _numActiveSplits--; }

  VirtualIterator<std::string> toSimpleClauseStrings();

  void setAux()
  {
    ASS(_auxInUse);
    _auxTimestamp=_auxCurrTimestamp;
  }

  /** Set auxiliary value of this clause. */
  void setAux(void* ptr)
  {
    ASS(_auxInUse);
    _auxTimestamp=_auxCurrTimestamp;
    _auxData=ptr;
  }
  /**
   * If there is an auxiliary value stored in this clause,
   * return true and assign it into @b ptr. Otherwise
   * return false.
   */
  template<typename T>
  bool tryGetAux(T*& ptr)
  {
    ASS(_auxInUse);
    if(_auxTimestamp==_auxCurrTimestamp) {
      ptr=static_cast<T*>(_auxData);
      return true;
    }
    return false;
  }
  /** Return auxiliary value stored in this clause. */
  template<typename T>
  T* getAux()
  {
    ASS(_auxInUse);
    ASS(_auxTimestamp==_auxCurrTimestamp);
    return static_cast<T*>(_auxData);
  }
  bool hasAux()
  {
    return _auxTimestamp==_auxCurrTimestamp;
  }

  /**
   * Request usage of the auxiliary value in clauses.
   * All aux. values stored in clauses before are guaranteed
   * to be discarded.
   */
  static void requestAux()
  {
#if VDEBUG
    ASS(!_auxInUse);
    _auxInUse=true;
#endif
    _auxCurrTimestamp++;
    if(_auxCurrTimestamp==0) {
      INVALID_OPERATION("Auxiliary clause value timestamp overflow!");
    }
  }
  /**
   * Announce that the auxiliary value in clauses is no longer
   * in use and can be used by someone else.
   */
  static void releaseAux()
  {
#if VDEBUG
    ASS(_auxInUse);
    _auxInUse=false;
#endif
  }

  unsigned splitWeight() const;
  unsigned getNumeralWeight() const;

  void collectVars(DHSet<unsigned>& acc);
  void collectUnstableVars(DHSet<unsigned>& acc);


  unsigned varCnt();
  unsigned maxVar(); // useful to create fresh variables w.r.t. the clause

  unsigned numPositiveLiterals(); // number of positive literals in the clause

  Literal* getAnswerLiteral();

  bool hasAnswerLiteral() {
    return getAnswerLiteral() != nullptr;
  }

  bool computable();

protected:
  /** number of literals */
  unsigned _length : 20;
  /** clause color, or COLOR_INVALID if not determined yet */
  mutable unsigned _color : 2;
  /** Clause was matched as extensionality and is tracked in the extensionality
    * clause container. The matching happens at activation. If the clause
    * becomes passive and is removed from the container, also this bit is unset.
    */
  unsigned _extensionality : 1;
  unsigned _extensionalityTag : 1;
  /** Clause is a splitting component. */
  unsigned _component : 1;

  /** storage class */
  Store _store : 3;
  /** number of selected literals */
  unsigned _numSelected : 20;

  /** weight */
  mutable unsigned _weight;
  /** weight for clause selection */
  unsigned _weightForClauseSelection;

  /** number of references to this clause */
  unsigned _refCnt;
  /** for splitting: timestamp marking when has the clause been reduced or restored by splitting */
  unsigned _reductionTimestamp;
  /** a map that translates Literal* to its index in the clause */
  InverseLookup<Literal>* _literalPositions;

  int _numActiveSplits;

  size_t _auxTimestamp;
  void* _auxData;

  static size_t _auxCurrTimestamp;
#if VDEBUG
  static bool _auxInUse;
#endif


  /** Array of literals of this unit */
  Literal* _literals[1];
}; // class Clause

std::ostream& operator<<(std::ostream& out, Clause::Store const& clause);
}

#endif
