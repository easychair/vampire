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
 * @file Indexing/Index.hpp
 * Defines abstract Index class and some other auxiliary classes.
 */

#ifndef __Indexing_Index__
#define __Indexing_Index__

#include "Forwards.hpp"
#include "Lib/Output.hpp"

#include "Lib/Event.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/Term.hpp"
#include "Lib/Exception.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Saturation/ClauseContainer.hpp"
#include "ResultSubstitution.hpp"
#include "Kernel/UnificationWithAbstraction.hpp"
#include "Lib/Allocator.hpp"
#include "Kernel/OrderingComparator.hpp"

/**
 * Indices are parametrized by a LeafData, i.e. the bit of data you want to store in the index.
 * Each leaf data must have a key to store at the leave. The Key can currently be either a Literal* or a TypedTermList.
 * A LeafData must have a function  `<Key> key() const;` returns the key, and must have comparison operators (<=,<,>=,>,!=,==) implemented.
 * See e.g. TermLiteralClause below for examples.
 */
namespace Indexing
{
using namespace Kernel;
using namespace Lib;
using namespace Saturation;

struct LiteralClause 
{
  Literal* const& key() const
  { return literal; }

private:
  std::tuple<unsigned,unsigned> asTuple() const
  { return std::make_tuple(clause->number(), literal->getId()); }
public:

  IMPL_COMPARISONS_FROM_TUPLE(LiteralClause)

  Literal* literal = nullptr;
  Clause* clause = nullptr;

  friend std::ostream& operator<<(std::ostream& out, LiteralClause const& self)
  { return out << "{ " << Output::ptr(self.clause) << ", " << Output::ptr(self.literal) << " }"; }
};

template<class Value>
struct TermWithValue {
  TypedTermList term;
  Value value;

  TermWithValue() {}

  TermWithValue(TypedTermList term, Value v)
    : term(term)
    , value(std::move(v))
  {}

  TypedTermList const& key() const { return term; }

  std::tuple<const TypedTermList &,const Value &> asTuple() const
  { return std::tie(term, value); }

  IMPL_COMPARISONS_FROM_TUPLE(TermWithValue)

  friend std::ostream& operator<<(std::ostream& out, TermWithValue const& self)
  { return out << self.asTuple(); }
};

class TermWithoutValue : public TermWithValue<std::tuple<>> 
{
public:
  TermWithoutValue(TypedTermList t) 
    : TermWithValue(t, std::make_tuple()) 
  { }
};

struct TermLiteralClause 
{
  TypedTermList term;
  Literal* literal = nullptr;
  Clause* clause = nullptr;

  TypedTermList const& key() const { return term; }

  auto  asTuple() const
  { return std::make_tuple(clause->number(), literal->getId(), term); }

  IMPL_COMPARISONS_FROM_TUPLE(TermLiteralClause)

  friend std::ostream& operator<<(std::ostream& out, TermLiteralClause const& self)
  { return out << "("
               << self.term << ", "
               << self.literal
               << Output::ptr(self.clause)
               << ")"; }
};

/** Custom leaf data for forward demodulation to store the demodulator
 * left- and right-hand side normalized and cache preorderedness. */
struct DemodulatorData
{
  DemodulatorData(TypedTermList term, TermList rhs, Clause* clause, bool preordered, const Ordering& ord)
    : term(term), rhs(rhs), clause(clause), preordered(preordered), comparator(ord.createComparator())
  {
    comparator->insert({ { term, rhs, Ordering::GREATER } }, (void*)0x1);
#if VDEBUG
    ASS(term.containsAllVariablesOf(rhs));
    ASS(!preordered || ord.compare(term,rhs)==Ordering::GREATER);
    Renaming r;
    r.normalizeVariables(term);
    ASS_EQ(term,r.apply(term));
    ASS_EQ(rhs,r.apply(rhs));
#endif
  }

  // lhs, the identifier is required to be `term` by CodeTree
  TypedTermList term;
  TermList rhs;
  Clause* clause;
  bool preordered; // whether term > rhs
  OrderingComparatorUP comparator; // owned comparator for whether term > rhs

  TypedTermList const& key() const { return term; }

  auto asTuple() const
  { return std::make_tuple(clause->number(), term, rhs); }

  IMPL_COMPARISONS_FROM_TUPLE(DemodulatorData)

  friend std::ostream& operator<<(std::ostream& out, DemodulatorData const& self)
  { return out << "(" << self.term << " = " << self.rhs << Output::ptr(self.clause) << ")"; }
};

template<class T>
struct is_indexed_data_normalized
{ static constexpr bool value = false; };

template<>
struct is_indexed_data_normalized<DemodulatorData>
{ static constexpr bool value = true; };

/**
 * Class of objects which contain results of term queries.
 */
template<class Unifier, class Data>
struct QueryRes
{
  Unifier unifier;
  Data const* data;

  QueryRes() {}
  QueryRes(Unifier unifier, Data const* data) 
    : unifier(std::move(unifier))
    , data(std::move(data)) {}

  friend std::ostream& operator<<(std::ostream& out, QueryRes const& self)
  { 
    return out 
      << "{ data: " << self.data()
      << ", unifier: " << self.unifier
      << "}";
  }
};

template<class Unifier, class Data>
QueryRes<Unifier, Data> queryRes(Unifier unifier, Data const* d) 
{ return QueryRes<Unifier, Data>(std::move(unifier), std::move(d)); }

class Index
{
public:
  virtual ~Index();

  void attachContainer(ClauseContainer* cc);
protected:
  Index() {}

  void onAddedToContainer(Clause* c)
  { handleClause(c, true); }
  void onRemovedFromContainer(Clause* c)
  { handleClause(c, false); }

  virtual void handleClause(Clause* c, bool adding) {}

  //TODO: postponing index modifications during iteration (methods isBeingIterated() etc...)

private:
  SubscriptionData _addedSD;
  SubscriptionData _removedSD;
};

};
#endif /*__Indexing_Index__*/
