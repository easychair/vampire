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
 * @file InferenceStore.hpp
 * Defines class InferenceStore.
 */


#ifndef __InferenceStore__
#define __InferenceStore__

#include <utility>
#include <ostream>

#include "Forwards.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/DHMultiset.hpp"
#include "Lib/Stack.hpp"

#include "Kernel/Signature.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/Inference.hpp"

namespace Kernel {

using namespace Lib;

class InferenceStore
{
public:
  static InferenceStore* instance();

  typedef List<int> IntList;

  struct FullInference
  {
    FullInference(unsigned premCnt) : csId(0), premCnt(premCnt) { }

    void* operator new(size_t,unsigned premCnt)
    {
      size_t size=sizeof(FullInference)+premCnt*sizeof(Unit*);
      size-=sizeof(Unit*);

      return ALLOC_KNOWN(size,"InferenceStore::FullInference");
    }

    size_t occupiedBytes()
    {
      size_t size=sizeof(FullInference)+premCnt*sizeof(Unit*);
      size-=sizeof(Unit*);
      return size;
    }

    void increasePremiseRefCounters();

    int csId;
    unsigned premCnt;
    InferenceRule rule;
    Unit* premises[1];
  };

  void recordSplittingNameLiteral(Unit* us, Literal* lit);
  void recordIntroducedSymbol(Unit* u, SymbolType st, unsigned number);
  void recordIntroducedSplitName(Unit* u, std::string name);

  void outputUnsatCore(std::ostream& out, Unit* refutation);
  void outputProof(std::ostream& out, Unit* refutation);
  void outputProof(std::ostream& out, UnitList* units);
  struct ProofPrinter;

private:
  struct TPTPProofPrinter;
  struct ProofCheckPrinter;
  struct ProofPropertyPrinter;
  struct DeduktiProofPrinter;

  ProofPrinter* createProofPrinter(std::ostream& out);

  DHMultiset<Clause*> _nextClIds;

  DHMap<Unit*, Literal*> _splittingNameLiterals;


  /** first records the type of the symbol (PRED,FUNC or TYPE_CON), second is symbol number */
  typedef std::pair<SymbolType,unsigned> SymbolId;
  typedef Stack<SymbolId> SymbolStack;
  DHMap<unsigned,SymbolStack> _introducedSymbols;
  DHMap<unsigned,std::string> _introducedSplitNames;

};

};

#endif /* __InferenceStore__ */
