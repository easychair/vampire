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
 * @file OperatorType.hpp
 * Defines class OperatorType.
 */

#ifndef __OperatorType__
#define __OperatorType__

#include "Forwards.hpp"

#include "Lib/DArray.hpp"
#include "Lib/Map.hpp"
#include "Lib/Set.hpp"
#include "Lib/DHSet.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/Stack.hpp"
#include "Lib/Vector.hpp"
#include "Lib/Allocator.hpp"

#include "Term.hpp"

namespace Kernel {

/**
 * The OperatorType class represents the predicate and function types
 * the only difference between the two is that a predicate type has return type
 * $o whilst a function type has any other return type.
 *
 * The class can be used to store polymorphic types which are of the form:
 * !>[alpha_0, ..., alpha_m](sort1 * ... * sortn) > return_sort where sorts can only contain variables
 * from  in {alpha_0,...,alpha_m}. View "A Polymorphic Vampire" for more details:
 *
 * https://link.springer.com/chapter/10.1007/978-3-030-51054-1_21
 *
 * The class stores data in a Vector<TermList>*, of length
 * num_of_arg_sorts + 1. The number of bound variables is stored in
 * a field _typeArgsArity. It is assumed that the bound variables range from
 * 0 to _typeArgsArity - 1. In order for this assumption to be valid all sorts must
 * normalised (view SortHelper::normaliseArgSorts) before being passed to a type forming
 * function.
 *
 * The objects of this class are perfectly shared (so that equal predicate / function types correspond to equal pointers)
 * and are obtained via static methods (to guarantee the sharing).
 */
class OperatorType
{
public:
  class TypeHash {
  public:
    static bool equals(OperatorType* t1, OperatorType* t2)
    { return (*t1) == (*t2); }

    static unsigned hash(OperatorType* ot)
    {
      OperatorKey& key = *ot->key();
      unsigned typeArgsArity = ot->numTypeArguments();
      return HashUtils::combine(
        DefaultHash::hash(key),
        DefaultHash::hash(typeArgsArity)
      );
    }
  };

private:
  typedef Vector<TermList> OperatorKey; // Vector of argument sorts together with "0" appended for predicates and resultSort appended for functions
  OperatorKey* _key;
  unsigned _typeArgsArity; /** number of quantified variables of this type */

  // constructors kept private
  OperatorType(OperatorKey* key, unsigned vLength) : _key(key), _typeArgsArity(vLength) {}

  /**
   * Convenience functions for creating a key
   */
  static OperatorKey* setupKey(unsigned arity, const TermList* sorts);
  static OperatorKey* setupKey(std::initializer_list<TermList> sorts);
  static OperatorKey* setupKeyUniformRange(unsigned arity, TermList argsSort);

  typedef Set<OperatorType*,TypeHash> OperatorTypes;
  static OperatorTypes& operatorTypes(); // just a wrapper around a static OperatorTypes object, to ensure a correct initialization order

  static OperatorType* getTypeFromKey(OperatorKey* key, unsigned taArity);

public:
  ~OperatorType() { _key->deallocate(); }

  inline bool operator==(const OperatorType& t) const
  { return  *_key==*t._key &&
             _typeArgsArity==t._typeArgsArity; }

  static OperatorType* getPredicateType(unsigned arity, const TermList* sorts=0, unsigned taArity = 0) {
    OperatorKey* key = setupKey(arity,sorts);
    (*key)[arity] = TermList::empty();
    return getTypeFromKey(key,taArity);
  }

  static OperatorType* getPredicateType(std::initializer_list<TermList> sorts, unsigned taArity = 0) {
    OperatorKey* key = setupKey(sorts);
    (*key)[sorts.size()] = TermList::empty();
    return getTypeFromKey(key,taArity);
  }

  static OperatorType* getPredicateTypeUniformRange(unsigned arity, TermList argsSort, unsigned taArity = 0) {
    OperatorKey* key = setupKeyUniformRange(arity,argsSort);
    (*key)[arity] = TermList::empty();
    return getTypeFromKey(key, taArity);
  }

  static OperatorType* getFunctionType(unsigned arity, const TermList* sorts, TermList resultSort, unsigned taArity = 0) {
    OperatorKey* key = setupKey(arity,sorts);
    (*key)[arity] = resultSort;
    return getTypeFromKey(key, taArity);
  }

  static OperatorType* getFunctionType(std::initializer_list<TermList> sorts, TermList resultSort, unsigned taArity = 0) {
    OperatorKey* key = setupKey(sorts);
    (*key)[sorts.size()] = resultSort;
    return getTypeFromKey(key,taArity);
  }

  static OperatorType* getFunctionTypeUniformRange(unsigned arity, TermList argsSort, TermList resultSort, unsigned taArity = 0) {
    OperatorKey* key = setupKeyUniformRange(arity,argsSort);
    (*key)[arity] = resultSort;
    return getTypeFromKey(key,taArity);
  }

  /**
   * Convenience function for creating OperatorType for constants (as symbols).
   * Constants are function symbols of 0 arity, so just provide the result sort.
   */
  static OperatorType* getConstantsType(TermList resultSort, unsigned taArity = 0) {
    return getFunctionType(0,nullptr,resultSort, taArity);
  }

  /**
   * Convenience function for creating OperatorType for type constructors.
   */
  static OperatorType* getTypeConType(unsigned arity) {
    return getFunctionTypeUniformRange(arity, AtomicSort::superSort(), AtomicSort::superSort());
  }

  OperatorKey* key() const { return _key; }
  unsigned numTypeArguments() const { return _typeArgsArity; }
  unsigned arity() const { return _typeArgsArity + _key->length()-1; }

  /**
   * These are the free variables of the sorts of the arguments (and the return type in case of functions) in the polymorphic case.
   *
   * For the non-polymorphic case, _typeArgsArity==0 and this function cannot be used meaninfully.
   */
  TermList quantifiedVar(unsigned idx) const
  {
    ASS(idx < _typeArgsArity);
    return TermList(idx, false);
  }

   /**
    * In the polymorhpic case, the first _typeArgsArity arguments of a predicate / function symbol are actually sorts, so their sort is "superSort"
    *
    * For idx >= _typeArgsArity (== 0 in the non-polymorphic case) this returns the actual sort of the idx-th argument
    */
  TermList arg(unsigned idx) const
  {
    if(idx < _typeArgsArity){
      return AtomicSort::superSort();
    }
    return (*_key)[idx - _typeArgsArity];
  }

  //TODO functions below do not hold for higher-order
  //In higher-order we have boolean functions

  bool isPredicateType() const { return (*_key)[arity() - numTypeArguments()].isEmpty(); };
  bool isFunctionType() const { return !(*_key)[arity() - numTypeArguments()].isEmpty(); };

  /**
   * The result sort of function types; or empty for predicates.
   */
  TermList result() const {
    return (*_key)[arity() - numTypeArguments()];
  }

  friend std::ostream& operator<<(std::ostream& out, OperatorType const& self)
  { return out << self.toString(); }

  std::string toString() const;

  bool isSingleSortType(TermList sort) const;
  bool isAllDefault() const { return isSingleSortType(AtomicSort::defaultSort()); }

private:
  std::string argsToString() const;
};

}

#endif // __OperatorType__
