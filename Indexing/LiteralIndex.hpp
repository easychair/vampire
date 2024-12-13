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
 * @file LiteralIndex.hpp
 * Defines class LiteralIndex.
 */


#ifndef __LiteralIndex__
#define __LiteralIndex__

#include "Lib/Output.hpp"
#include "Lib/DHMap.hpp"

#include "Index.hpp"
#include "LiteralIndexingStructure.hpp"


namespace Indexing {

template<class Data>
class LiteralIndex
: public Index
{
public:
  VirtualIterator<LiteralClause> getAll()
  { return _is->getAll(); }

  VirtualIterator<QueryRes<ResultSubstitutionSP, LiteralClause>> getUnifications(Literal* lit, bool complementary, bool retrieveSubstitutions = true)
  { return _is->getUnifications(lit, complementary, retrieveSubstitutions); }

  VirtualIterator<QueryRes<AbstractingUnifier*, Data>> getUwa(Literal* lit, bool complementary, Options::UnificationWithAbstraction uwa, bool fixedPointIteration)
  { return _is->getUwa(lit, complementary, uwa, fixedPointIteration); }

  VirtualIterator<QueryRes<ResultSubstitutionSP, LiteralClause>> getGeneralizations(Literal* lit, bool complementary, bool retrieveSubstitutions = true)
  { return _is->getGeneralizations(lit, complementary, retrieveSubstitutions); }

  VirtualIterator<QueryRes<ResultSubstitutionSP, LiteralClause>> getInstances(Literal* lit, bool complementary, bool retrieveSubstitutions = true)
  { return _is->getInstances(lit, complementary, retrieveSubstitutions); }

  size_t getUnificationCount(Literal* lit, bool complementary)
  { return _is->getUnificationCount(lit, complementary); }

  friend std::ostream& operator<<(std::ostream& out,                 LiteralIndex const& self) { return out << *self._is; }
  friend std::ostream& operator<<(std::ostream& out, Output::Multiline<LiteralIndex>const& self) { return out << Output::multiline(*self.self._is, self.indent); }

protected:
  LiteralIndex(LiteralIndexingStructure<Data>* is) : _is(is) {}

  void handle(Data data, bool add)
  { _is->handle(std::move(data), add); }

  std::unique_ptr<LiteralIndexingStructure<Data>> _is;
};

class BinaryResolutionIndex
: public LiteralIndex<LiteralClause>
{
public:
  BinaryResolutionIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex<LiteralClause>(is) {};
protected:
  void handleClause(Clause* c, bool adding);
};

class BackwardSubsumptionIndex
: public LiteralIndex<LiteralClause>
{
public:
  BackwardSubsumptionIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex<LiteralClause>(is) {};
protected:
  void handleClause(Clause* c, bool adding);
};

class FwSubsSimplifyingLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  FwSubsSimplifyingLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
    : LiteralIndex<LiteralClause>(is)
  { }

protected:
  void handleClause(Clause* c, bool adding) override;
};

class FSDLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  FSDLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
    : LiteralIndex<LiteralClause>(is)
  { }

protected:
  void handleClause(Clause* c, bool adding) override;
};

class UnitClauseLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  UnitClauseLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex<LiteralClause>(is) {};
protected:
  void handleClause(Clause* c, bool adding);
};

class UnitClauseWithALLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  UnitClauseWithALLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex(is) {};
protected:
  void handleClause(Clause* c, bool adding);
};

class NonUnitClauseLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  NonUnitClauseLiteralIndex(LiteralIndexingStructure<LiteralClause>* is, bool selectedOnly=false)
  : LiteralIndex<LiteralClause>(is), _selectedOnly(selectedOnly) {};
protected:
  void handleClause(Clause* c, bool adding);
private:
  bool _selectedOnly;
};

class NonUnitClauseWithALLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  NonUnitClauseWithALLiteralIndex(LiteralIndexingStructure<LiteralClause>* is, bool selectedOnly=false)
  : LiteralIndex(is), _selectedOnly(selectedOnly) {};
protected:
  void handleClause(Clause* c, bool adding);
private:
  bool _selectedOnly;
};

class RewriteRuleIndex
: public LiteralIndex<LiteralClause>
{
public:
  RewriteRuleIndex(LiteralIndexingStructure<LiteralClause>* is, Ordering& ordering);
  ~RewriteRuleIndex();

  Clause* getCounterpart(Clause* c) {
    return _counterparts.get(c);
  }
protected:
  void handleClause(Clause* c, bool adding);
  Literal* getGreater(Clause* c);

private:
  void handleEquivalence(Clause* c, Literal* cgr, Clause* d, Literal* dgr, bool adding);

  LiteralIndexingStructure<LiteralClause>* _partialIndex;
  DHMap<Clause*,Clause*> _counterparts;
  Ordering& _ordering;
};

class DismatchingLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  DismatchingLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex<LiteralClause>(is) {};
  void handleClause(Clause* c, bool adding);
  void addLiteral(Literal* c);
};

class UnitIntegerComparisonLiteralIndex
: public LiteralIndex<LiteralClause>
{
public:
  UnitIntegerComparisonLiteralIndex(LiteralIndexingStructure<LiteralClause>* is)
  : LiteralIndex<LiteralClause>(is) {}

protected:
  void handleClause(Clause* c, bool adding);
};

};

#endif /* __LiteralIndex__ */
