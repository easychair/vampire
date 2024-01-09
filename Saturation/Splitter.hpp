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
 * @file Splitter.hpp
 * Defines class Splitter for SAT solver-driven splitting (= AVATAR).
 */

#ifndef __Splitter__
#define __Splitter__

#include "Forwards.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/ArrayMap.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/Hash.hpp"
#include "Lib/Stack.hpp"
#include "Lib/ScopedPtr.hpp"

#include "Shell/Options.hpp"

#include "Kernel/RCClauseStack.hpp"

#include "Indexing/ClauseVariantIndex.hpp"

#include "SAT/SAT2FO.hpp"
#include "SAT/SATLiteral.hpp"
#include "SAT/SATSolver.hpp"

#include "DP/DecisionProcedure.hpp"
#include "DP/SimpleCongruenceClosure.hpp"

#include "Lib/Allocator.hpp"

namespace Saturation {

using namespace Lib;
using namespace Kernel;
using namespace Shell;
using namespace SAT;
using namespace DP;
using namespace Indexing;

typedef Stack<SplitLevel> SplitLevelStack;

class Splitter;

/**
 * Object that decides which splitting components are to be selected
 */
class SplittingBranchSelector {
public:
  SplittingBranchSelector(Splitter& parent) : _ccModel(false), _parent(parent), _solverIsSMT(false)  {}
  ~SplittingBranchSelector(){
#if VZ3
_solver=0;
#endif
  }

  /** To be called from Splitter::init() */
  void init();

  void updateVarCnt();
  void considerPolarityAdvice(SATLiteral lit);
  void trySetTrue(SATLiteral lit);

  void addSatClauseToSolver(SATClause* cl, bool refutation);
  void recomputeModel(SplitLevelStack& addedComps, SplitLevelStack& removedComps, bool randomize = false);

  void flush(SplitLevelStack& addedComps, SplitLevelStack& removedComps);

private:
  friend class Splitter;

  SATSolver::Status processDPConflicts();
  SATSolver::VarAssignment getSolverAssimentConsideringCCModel(unsigned var);

  void handleSatRefutation();
  void updateSelection(unsigned satVar, SATSolver::VarAssignment asgn,
      SplitLevelStack& addedComps, SplitLevelStack& removedComps);

  int assertedGroundPositiveEqualityCompomentMaxAge();

  //options
  bool _eagerRemoval;
  Options::SplittingLiteralPolarityAdvice _literalPolarityAdvice;
  bool _ccMultipleCores;
  bool _minSCO; // minimize wrt splitting clauses only
  bool _ccModel;

  Splitter& _parent;

  bool _solverIsSMT;
  ScopedPtr<SATSolver> _solver;
  ScopedPtr<DecisionProcedure> _dp;
  // use a separate copy of the decision procedure for ccModel computations and fill it up only with equalities
  ScopedPtr<SimpleCongruenceClosure> _dpModel;
  
  /**
   * Contains selected component names (splitlevels)
   */
  ArraySet _selected;
  
  /**
   * Keeps track of positive ground equalities true in the last ccmodel.
   */
  ArraySet _trueInCCModel;

#if VDEBUG
  unsigned lastCheckedVar;
#endif
};


/**
 * SplitLevel meanings:
 * even -- positive ground literals and non-ground components
 * odd -- negative ground literals (and possibly skolemized non-ground components -- not implemented yet)
 */
class Splitter {
private:

/**
 * ReductionRecord - records information to do with a reduction
 *
 */

  struct ReductionRecord
  {
    ReductionRecord(Clause* clause) : clause(clause), 
        timestamp(clause->getReductionTimestamp()) {}
    Clause* clause;
    unsigned timestamp;
  };

/**
 *
 * SplitRecord - records the split information for the clause component
 *
 * Let's call the SplitLevel associated with a comp its "name"
 * A corresponding SplitRecord is added to _db[name] 
 *
 * children - Clauses that rely on name (of comp), should be thrown away "on backtracking"
 * reduced - The clauses that have been *conditionally* reduced by this clause (and are therefore frozen)
 * active - component currently true in the model
 *
 * Comment by Giles
 */
  struct SplitRecord
  {
    SplitRecord(Clause* comp)
     : component(comp), active(false)
    {
      component->incRefCnt();
    }

    ~SplitRecord();

    void addReduced(Clause* cl);

    Clause* component;
    RCClauseStack children;
    Stack<ReductionRecord> reduced;
    bool active;

    USE_ALLOCATOR(SplitRecord);
  };
  
public:
  Splitter();
  ~Splitter();

  const Options& getOptions() const;
  Ordering& getOrdering() const;
  
  void init(SaturationAlgorithm* sa);

  bool doSplitting(Clause* cl);

  void onClauseReduction(Clause* cl, ClauseIterator premises, Clause* replacement);
  void onNewClause(Clause* cl);
  void onAllProcessed();
  bool handleEmptyClause(Clause* cl);

  SplitLevel getNameFromLiteral(SATLiteral lit) const;
  Unit* getDefinitionFromName(SplitLevel compName) const;

  static vstring splitsToString(SplitSet* splits);
  static SATLiteral getLiteralFromName(SplitLevel compName);
  static vstring getFormulaStringFromName(SplitLevel compName, bool negated = false);

  bool isUsedName(SplitLevel name) const {
    ASS_L(name,_db.size());
    return (_db[name] != 0);
  }
  Clause* getComponentClause(SplitLevel name) const;

  SplitLevel splitLevelCnt() const { return _db.size(); }
  unsigned maxSatVar() const { return _sat2fo.maxSATVar(); }

  SAT2FO& satNaming() { return _sat2fo; }

  UnitList* preprendCurrentlyAssumedComponentClauses(UnitList* clauses);
  static bool getComponents(Clause* cl, Stack<LiteralStack>& acc, bool shuffle = false);

  /*
   * Clauses with answer literals cannot be split -- hence if we obtain a clause with
   * avatar assertions that has an answer literal, we have to un-split it.
   * This method replaces `C \/ ans(r) <- A1,...,An` with `C \/ ans(r) \/ ~A1 \/ ... \/ ~An`
   */
  Clause* reintroduceAvatarAssertions(Clause* cl);

private:
  friend class SplittingBranchSelector;
  
  SplitLevel getNameFromLiteralUnsafe(SATLiteral lit) const;

  bool shouldAddClauseForNonSplittable(Clause* cl, unsigned& compName, Clause*& compCl);
  bool handleNonSplittable(Clause* cl);
  bool tryGetExistingComponentName(unsigned size, Literal* const * lits, SplitLevel& comp, Clause*& compCl);

  void addComponents(const SplitLevelStack& toAdd);
  void removeComponents(const SplitLevelStack& toRemove);

  void collectDependenceLits(SplitSet* splits, SATLiteralStack& acc) const;

  SplitLevel addNonGroundComponent(unsigned size, Literal* const * lits, Clause* orig, Clause*& compCl);
  SplitLevel addGroundComponent(Literal* lit, Clause* orig, Clause*& compCl);

  Clause* buildAndInsertComponentClause(SplitLevel name, unsigned size, Literal* const * lits, Clause* orig=0);

  SplitLevel tryGetComponentNameOrAddNew(const LiteralStack& comp, Clause* orig, Clause*& compCl);
  SplitLevel tryGetComponentNameOrAddNew(unsigned size, Literal* const * lits, Clause* orig, Clause*& compCl);

  void addSatClauseToSolver(SATClause* cl, bool refutation);

  SplitSet* getNewClauseSplitSet(Clause* cl);
  void assignClauseSplitSet(Clause* cl, SplitSet* splits);

  bool allSplitLevelsActive(SplitSet* s);

  void conjecture(unsigned, Literal **);

  //settings
  bool _showSplitting;

  Options::SplittingAddComplementary _complBehavior;
  Options::SplittingNonsplittableComponents _nonsplComps;
  unsigned _flushPeriod;
  float _flushQuotient;
  Options::SplittingDeleteDeactivated _deleteDeactivated;
  Options::SplittingCongruenceClosure _congruenceClosure;
  bool _shuffleComponents;
#if VZ3
  bool hasSMTSolver;
#endif

  //utility objects
  SplittingBranchSelector _branchSelector;
  ScopedPtr<ClauseVariantIndex> _componentIdx;
  /**
   * Registers all the sat variables and keeps track
   * of associated ground literals for those variables
   * which have one.
   */  
  SAT2FO _sat2fo;  
  /**
   * Information about a split level. Can be null if a split level does
   * not contain any components (e.g. for negations of non-ground
   * components)
   *
   * Invariant: if there is a clause with a level in its splitting history,
   * the _db record of this level is non-null.
   */
  Stack<SplitRecord*> _db;

  /**
   * Definitions of ground components C and ~C are shared and placed at the slot of C.
   * (So the key here is never odd!)
   **/
  DHMap<SplitLevel,Unit*> _defs;
  
  //state variable used for flushing:  
  /** When this number of generated clauses is reached, it will cause flush */
  unsigned _flushThreshold;
  /** true if there was a clause added to the SAT solver since last call to onAllProcessed */
  bool _clausesAdded;
  /** true if there was a refutation added to the SAT solver */
  bool _haveBranchRefutation;

  /* as there can be both limits, it's hard to covert between them,
   * and we terminate at the earlier one, let's just keep checking both. */
  unsigned _stopSplittingAtTime; // time elapsed in milliseconds
#if VAMPIRE_PERF_EXISTS
  unsigned _stopSplittingAtInst; // mega-instructions elapsed
#endif

  bool _fastRestart; // option's value copy
  /**
   * We are postponing to consider these clauses for a split 
   * because a conflict clause has been derived
   * and will invariably change the SAT model.
   */
  RCClauseStack _fastClauses;

  // clauses we already added to the SAT solver
  // not just optimisation: also prevents the SAT solver oscillating between two models in some cases
  Set<SATClause *, DerefPtrHash<DefaultHash>> _already_added;

  SaturationAlgorithm* _sa;
public:
  static vstring splPrefix;

  // for observing the current model
  SplitLevel splitLevelBound() { return _db.size(); }
  bool splitLevelActive(SplitLevel lev) {
    ASS_REP(lev<_db.size(), lev);
    return (_db[lev]!=0 && _db[lev]->active);
  }
};

}

#endif // __Splitter__
