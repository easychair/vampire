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
 * @file LiteralIndexingStructure.hpp
 * Defines class LiteralIndexingStructure.
 */


#ifndef __LiteralIndexingStructure__
#define __LiteralIndexingStructure__

#include "Forwards.hpp"
#include "Index.hpp"
#include "Kernel/UnificationWithAbstraction.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Shell/Options.hpp"

namespace Indexing {

template<class LeafData_>
class LiteralIndexingStructure {
public:
  using LeafData = LeafData_;
  virtual ~LiteralIndexingStructure() {}

  virtual void handle(LeafData ld, bool insert) = 0;
  void insert(LeafData ld) { handle(std::move(ld), /* insert = */ true ); }
  void remove(LeafData ld) { handle(std::move(ld), /* insert = */ false); }

  virtual VirtualIterator<LeafData> getAll() { NOT_IMPLEMENTED; }
  virtual VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>> getUnifications(Literal* lit, bool complementary, bool retrieveSubstitutions = true) { NOT_IMPLEMENTED; }
  virtual VirtualIterator<QueryRes<AbstractingUnifier*, LeafData>> getUwa(Literal* lit, bool complementary, Options::UnificationWithAbstraction uwa, bool fixedPointIteration) = 0;
  virtual VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>> getGeneralizations(Literal* lit, bool complementary, bool retrieveSubstitutions = true) { NOT_IMPLEMENTED; }
  virtual VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>> getInstances(Literal* lit, bool complementary, bool retrieveSubstitutions = true) { NOT_IMPLEMENTED; }
  virtual VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>> getVariants(Literal* lit, bool complementary, bool retrieveSubstitutions = true) { NOT_IMPLEMENTED; }

  virtual size_t getUnificationCount(Literal* lit, bool complementary)
  {
    return countIteratorElements(getUnifications(lit, complementary, false));
  }

  virtual void output(std::ostream& out, Option<unsigned> multilineIndent) const = 0;

  friend std::ostream& operator<<(std::ostream& out,                 LiteralIndexingStructure const& self) {      self.output(out, {}               ); return out; }
  friend std::ostream& operator<<(std::ostream& out, Output::Multiline<LiteralIndexingStructure>const& self) { self.self.output(out, some(self.indent)); return out; }
};

};

#endif /* __LiteralIndexingStructure__ */
