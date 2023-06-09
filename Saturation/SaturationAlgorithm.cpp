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
 * @file SaturationAlgorithm.cpp
 * Implementing SaturationAlgorithm class.
 */
#include "Debug/RuntimeStatistics.hpp"

#include "Lib/DHSet.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/SharedSet.hpp"
#include "Lib/Stack.hpp"
#include "Lib/Timer.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Lib/System.hpp"

#include "Indexing/LiteralIndexingStructure.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Kernel/FormulaUnit.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/InferenceStore.hpp"
#include "Kernel/KBO.hpp"
#include "Kernel/LiteralSelector.hpp"
#include "Kernel/MLVariant.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/SubformulaIterator.hpp"
#include "Kernel/Unit.hpp"

#include "Inferences/InterpretedEvaluation.hpp"
#include "Inferences/PolynomialEvaluation.hpp"
#include "Inferences/PushUnaryMinus.hpp"
#include "Inferences/Cancellation.hpp"
#include "Inferences/GaussianVariableElimination.hpp"
#include "Inferences/EquationalTautologyRemoval.hpp"
#include "Inferences/Condensation.hpp"
#include "Inferences/FastCondensation.hpp"
#include "Inferences/DistinctEqualitySimplifier.hpp"

#include "Inferences/InferenceEngine.hpp"
#include "Inferences/BackwardDemodulation.hpp"
#include "Inferences/BackwardSubsumptionAndResolution.hpp"
#include "Inferences/BackwardSubsumptionDemodulation.hpp"
#include "Inferences/BinaryResolution.hpp"
#include "Inferences/CodeTreeForwardSubsumptionAndResolution.hpp"
#include "Inferences/EqualityFactoring.hpp"
#include "Inferences/EqualityResolution.hpp"
#include "Inferences/BoolEqToDiseq.hpp"
#include "Inferences/ExtensionalityResolution.hpp"
#include "Inferences/FOOLParamodulation.hpp"
#include "Inferences/Injectivity.hpp"
#include "Inferences/Factoring.hpp"
#include "Inferences/FunctionDefinitionRewriting.hpp"
#include "Inferences/ForwardDemodulation.hpp"
#include "Inferences/ForwardLiteralRewriting.hpp"
#include "Inferences/ForwardSubsumptionAndResolution.hpp"
#include "Inferences/InvalidAnswerLiteralRemovals.hpp"
#include "Inferences/ForwardSubsumptionDemodulation.hpp"
#include "Inferences/GlobalSubsumption.hpp"
#include "Inferences/InnerRewriting.hpp"
#include "Inferences/TermAlgebraReasoning.hpp"
#include "Inferences/Superposition.hpp"
#include "Inferences/ArgCong.hpp"
#include "Inferences/NegativeExt.hpp"
#include "Inferences/Narrow.hpp"
#include "Inferences/PrimitiveInstantiation.hpp"
#include "Inferences/Choice.hpp"
#include "Inferences/ElimLeibniz.hpp"
#include "Inferences/SubVarSup.hpp"
#include "Inferences/CNFOnTheFly.hpp"
#include "Inferences/URResolution.hpp"
#include "Inferences/Instantiation.hpp"
#include "Inferences/TheoryInstAndSimp.hpp"
#include "Inferences/Induction.hpp"
#include "Inferences/ArithmeticSubtermGeneralization.hpp"
#include "Inferences/TautologyDeletionISE.hpp"
#include "Inferences/CombinatorDemodISE.hpp"
#include "Inferences/CombinatorNormalisationISE.hpp"
#include "Inferences/BoolSimp.hpp"
#include "Inferences/CasesSimp.hpp"
#include "Inferences/Cases.hpp"
#include "Inferences/DefinitionIntroduction.hpp"

#include "Saturation/ExtensionalityClauseContainer.hpp"

#include "Shell/AnswerLiteralManager.hpp"
#include "Shell/ConditionalRedundancyHandler.hpp"
#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/UIHelper.hpp"
#include "Shell/Shuffling.hpp"

#include "Splitter.hpp"

#include "ConsequenceFinder.hpp"
#include "LabelFinder.hpp"
#include "Splitter.hpp"
#include "SymElOutput.hpp"
#include "SaturationAlgorithm.hpp"
#include "ManCSPassiveClauseContainer.hpp"
#include "AWPassiveClauseContainer.hpp"
#include "PredicateSplitPassiveClauseContainer.hpp"
#include "Discount.hpp"
#include "LRS.hpp"
#include "Otter.hpp"

using namespace std;
using namespace Lib;
using namespace Kernel;
using namespace Shell;
using namespace Saturation;

/** Print information changes in clause containers */
#define REPORT_CONTAINERS 0
/** Print information about performed forward simplifications */
#define REPORT_FW_SIMPL 0
/** Print information about performed backward simplifications */
#define REPORT_BW_SIMPL 0

SaturationAlgorithm* SaturationAlgorithm::s_instance = 0;

std::unique_ptr<PassiveClauseContainer> makeLevel0(bool isOutermost, const Options& opt, std::string name)
{
  return std::make_unique<AWPassiveClauseContainer>(isOutermost, opt, name + "AWQ");
}

std::unique_ptr<PassiveClauseContainer> makeLevel1(bool isOutermost, const Options& opt, std::string name)
{
  if (opt.useTheorySplitQueues()) {
    std::vector<std::unique_ptr<PassiveClauseContainer>> queues;
    auto cutoffs = opt.theorySplitQueueCutoffs();
    for (unsigned i = 0; i < cutoffs.size(); i++) {
      auto queueName = name + "ThSQ" + Int::toString(cutoffs[i]) + ":";
      queues.push_back(makeLevel0(false, opt, queueName));
    }
    return std::make_unique<TheoryMultiSplitPassiveClauseContainer>(isOutermost, opt, name + "ThSQ", std::move(queues));
  }
  else {
    return makeLevel0(isOutermost, opt, name);
  }
}

std::unique_ptr<PassiveClauseContainer> makeLevel2(bool isOutermost, const Options& opt, std::string name)
{
  if (opt.useAvatarSplitQueues()) {
    std::vector<std::unique_ptr<PassiveClauseContainer>> queues;
    auto cutoffs = opt.avatarSplitQueueCutoffs();
    for (unsigned i = 0; i < cutoffs.size(); i++) {
      auto queueName = name + "AvSQ" + Int::toString(cutoffs[i]) + ":";
      queues.push_back(makeLevel1(false, opt, queueName));
    }
    return std::make_unique<AvatarMultiSplitPassiveClauseContainer>(isOutermost, opt, name + "AvSQ", std::move(queues));
  }
  else {
    return makeLevel1(isOutermost, opt, name);
  }
}

std::unique_ptr<PassiveClauseContainer> makeLevel3(bool isOutermost, const Options& opt, std::string name)
{
  if (opt.useSineLevelSplitQueues()) {
    std::vector<std::unique_ptr<PassiveClauseContainer>> queues;
    auto cutoffs = opt.sineLevelSplitQueueCutoffs();
    for (unsigned i = 0; i < cutoffs.size(); i++) {
      auto queueName = name + "SLSQ" + Int::toString(cutoffs[i]) + ":";
      queues.push_back(makeLevel2(false, opt, queueName));
    }
    return std::make_unique<SineLevelMultiSplitPassiveClauseContainer>(isOutermost, opt, name + "SLSQ", std::move(queues));
  }
  else {
    return makeLevel2(isOutermost, opt, name);
  }
}

std::unique_ptr<PassiveClauseContainer> makeLevel4(bool isOutermost, const Options& opt, std::string name)
{
  if (opt.usePositiveLiteralSplitQueues()) {
    std::vector<std::unique_ptr<PassiveClauseContainer>> queues;
    std::vector<float> cutoffs = opt.positiveLiteralSplitQueueCutoffs();
    for (unsigned i = 0; i < cutoffs.size(); i++) {
      auto queueName = name + "PLSQ" + Int::toString(cutoffs[i]) + ":";
      queues.push_back(makeLevel3(false, opt, queueName));
    }
    return std::make_unique<PositiveLiteralMultiSplitPassiveClauseContainer>(isOutermost, opt, name + "PLSQ", std::move(queues));
  }
  else {
    return makeLevel3(isOutermost, opt, name);
  }
}

/**
 * Create a SaturationAlgorithm object
 *
 * The @b passiveContainer object will be used as a passive clause container, and
 * @b selector object to select literals before clauses are activated.
 */
SaturationAlgorithm::SaturationAlgorithm(Problem& prb, const Options& opt)
  : MainLoop(prb, opt),
    _clauseActivationInProgress(false),
    _fwSimplifiers(0), _simplifiers(0), _bwSimplifiers(0), _splitter(0),
    _consFinder(0), _labelFinder(0), _symEl(0), _answerLiteralManager(0),
    _instantiation(0), _fnDefHandler(prb.getFunctionDefinitionHandler()),
    _generatedClauseCount(0),
    _activationLimit(0)
{
  ASS_EQ(s_instance, 0);  //there can be only one saturation algorithm at a time

  _activationLimit = opt.activationLimit();

  _ordering = OrderingSP(Ordering::create(prb, opt));
  if (!Ordering::trySetGlobalOrdering(_ordering)) {
    // this is not an error, it may just lead to lower performance (and most likely not significantly lower)
    cerr << "SaturationAlgorithm cannot set its ordering as global" << endl;
  }
  _selector = LiteralSelector::getSelector(*_ordering, opt, opt.selection());

  _completeOptionSettings = opt.complete(prb);

  _unprocessed = new UnprocessedClauseContainer();

  if (opt.useManualClauseSelection()) {
    _passive = std::make_unique<ManCSPassiveClauseContainer>(true, opt);
  }
  else {
    _passive = makeLevel4(true, opt, "");
  }
  _active = new ActiveClauseContainer(opt);

  _active->attach(this);
  _passive->attach(this);

  _active->addedEvent.subscribe(this, &SaturationAlgorithm::onActiveAdded);
  _active->removedEvent.subscribe(this, &SaturationAlgorithm::activeRemovedHandler);
  _passive->addedEvent.subscribe(this, &SaturationAlgorithm::onPassiveAdded);
  _passive->removedEvent.subscribe(this, &SaturationAlgorithm::passiveRemovedHandler);
  _passive->selectedEvent.subscribe(this, &SaturationAlgorithm::onPassiveSelected);
  _unprocessed->addedEvent.subscribe(this, &SaturationAlgorithm::onUnprocessedAdded);
  _unprocessed->removedEvent.subscribe(this, &SaturationAlgorithm::onUnprocessedRemoved);
  _unprocessed->selectedEvent.subscribe(this, &SaturationAlgorithm::onUnprocessedSelected);

  if (opt.extensionalityResolution() != Options::ExtensionalityResolution::OFF) {
    _extensionality = new ExtensionalityClauseContainer(opt);
    //_active->addedEvent.subscribe(_extensionality, &ExtensionalityClauseContainer::addIfExtensionality);
  }
  else {
    _extensionality = 0;
  }

  s_instance = this;
}

/**
 * Destroy the SaturationAlgorithm object
 */
SaturationAlgorithm::~SaturationAlgorithm()
{
  ASS_EQ(s_instance,this);

  s_instance = 0;

  if (_splitter) {
    delete _splitter;
  }
  if (_consFinder) {
    delete _consFinder;
  }
  if (_symEl) {
    delete _symEl;
  }

  _active->detach();
  _passive->detach();

  if (_generator) {
    _generator->detach();
  }
  if (_immediateSimplifier) {
    _immediateSimplifier->detach();
  }

  while (_fwSimplifiers) {
    ForwardSimplificationEngine* fse = FwSimplList::pop(_fwSimplifiers);
    fse->detach();
    delete fse;
  }
  while (_simplifiers) {
    SimplificationEngine* fse = SimplList::pop(_simplifiers);
    fse->detach();
    delete fse;
  }
  while (_bwSimplifiers) {
    BackwardSimplificationEngine* bse = BwSimplList::pop(_bwSimplifiers);
    bse->detach();
    delete bse;
  }

  delete _unprocessed;
  delete _active;
}

void SaturationAlgorithm::tryUpdateFinalClauseCount()
{
  SaturationAlgorithm* inst = tryGetInstance();
  if (!inst) {
    return;
  }
  env.statistics->finalActiveClauses = inst->_active->sizeEstimate();
  env.statistics->finalPassiveClauses = inst->_passive->sizeEstimate();
  if (inst->_extensionality != 0) {
    env.statistics->finalExtensionalityClauses = inst->_extensionality->size();
  }
}

/**
 * Return true if the run of the prover so far is complete
 */
bool SaturationAlgorithm::isComplete()
{
  return _completeOptionSettings && !env.statistics->inferencesSkippedDueToColors;
}

ClauseIterator SaturationAlgorithm::activeClauses()
{
  return _active->clauses();
}

/**
 * A function that is called when a clause is added to the active clause container.
 */
void SaturationAlgorithm::onActiveAdded(Clause* c)
{
  if (env.options->showActive()) {
    std::cout << "[SA] active: " << c->toString() << std::endl;
  }
}

/**
 * A function that is called when a clause is removed from the active clause container.
 */
void SaturationAlgorithm::onActiveRemoved(Clause* c)
{
  ASS(c->store()==Clause::ACTIVE);
  c->setStore(Clause::NONE);
  // at this point the c object may be deleted
}

void SaturationAlgorithm::onAllProcessed()
{
  ASS(clausesFlushed());

  if (_symEl) {
    _symEl->onAllProcessed();
  }

  if (_splitter) {
    _splitter->onAllProcessed();
  }

  if (_consFinder) {
    _consFinder->onAllProcessed();
  }
}

/**
 * A function that is called when a clause is added to the passive clause container.
 */
void SaturationAlgorithm::onPassiveAdded(Clause* c)
{
  if (env.options->showPassive()) {
    std::cout << "[SA] passive: " << c->toString() << std::endl;
  }

  //when a clause is added to the passive container,
  //we know it is not redundant
  onNonRedundantClause(c);
}

/**
 * A function that is called when a clause is removed from the active clause container.
 * The function is not called when a selected clause is removed from the passive container.
 * In this case the @b onPassiveSelected method is called.
 */
void SaturationAlgorithm::onPassiveRemoved(Clause* c)
{
  ASS(c->store()==Clause::PASSIVE);
  c->setStore(Clause::NONE);
  // at this point the c object can be deleted
}

/**
 * A function that is called when a clause is selected and removed from the passive
 * clause container to be activated.
 *
 * The clause @b c might not necessarily get to the activation, it can still be
 * removed by some simplification rule (in case of the Discount saturation algorithm).
 */
void SaturationAlgorithm::onPassiveSelected(Clause* c)
{
}

/**
 * A function that is called when a clause is added to the unprocessed clause container.
 */
void SaturationAlgorithm::onUnprocessedAdded(Clause* c)
{
}

/**
 * A function that is called when a clause is removed from the unprocessed clause container.
 */
void SaturationAlgorithm::onUnprocessedRemoved(Clause* c)
{
}

void SaturationAlgorithm::onUnprocessedSelected(Clause* c)
{
}

/**
 * A function that is called whenever a possibly new clause appears.
 */
void SaturationAlgorithm::onNewClause(Clause* cl)
{
  if (_splitter) {
    _splitter->onNewClause(cl);
  }

  if (env.options->showNew()) {
    std::cout << "[SA] new: " << cl->toString() << std::endl;
  }

  if (cl->isPropositional()) {
    onNewUsefulPropositionalClause(cl);
  }

  if (_answerLiteralManager) {
    _answerLiteralManager->onNewClause(cl);
  }
}

void SaturationAlgorithm::onNewUsefulPropositionalClause(Clause* c)
{
  ASS(c->isPropositional());

  if (env.options->showNewPropositional()) {
    std::cout << "[SA] new propositional: " << c->toString() << std::endl;
  }

  if (_consFinder) {
    _consFinder->onNewPropositionalClause(c);
  }
  if (_labelFinder) {
    _labelFinder->onNewPropositionalClause(c);
  }
}

/**
 * Called when a clause successfully passes the forward simplification
 */
void SaturationAlgorithm::onClauseRetained(Clause* cl)
{
  //cout << "[SA] retained " << cl->toString() << endl;

}

/**
 * Called whenever a clause is simplified or deleted at any point of the
 * saturation algorithm
 */
void SaturationAlgorithm::onClauseReduction(Clause* cl, Clause **replacements, unsigned numOfReplacements,
                                            Clause* premise, bool forward)
{
  ASS(cl);

  ClauseIterator premises;

  if (premise) {
    premises = pvi(getSingletonIterator(premise));
  }
  else {
    premises = ClauseIterator::getEmpty();
  }

  onClauseReduction(cl, replacements, numOfReplacements, premises, forward);
}

void SaturationAlgorithm::onClauseReduction(Clause* cl, Clause **replacements, unsigned numOfReplacements,
                                            ClauseIterator premises, bool forward)
{
  ASS(cl);

  static ClauseStack premStack;
  premStack.reset();
  premStack.loadFromIterator(premises);

  Clause *replacement = numOfReplacements ? *replacements : 0;

  if (env.options->showReductions()) {
    std::cout << "[SA] " << (forward ? "forward" : "backward") << " reduce: " << cl->toString() << endl;
    for(unsigned i = 0; i < numOfReplacements; i++){
      Clause* replacement = *replacements;
      if(replacement){ std::cout << "      replaced by " << replacement->toString() << endl; }
      replacements++;
    }
    ClauseStack::Iterator pit(premStack);
    while(pit.hasNext()){
      Clause* premise = pit.next();
      if(premise){ std::cout << "     using " << premise->toString() << endl; }
    }
  }

  if (_splitter) {
    _splitter->onClauseReduction(cl, pvi(ClauseStack::Iterator(premStack)), replacement);
  }

  if (replacement) {
    // Where an inference has multiple conclusions, onParenthood will only be run
    // for the final conclusion. This is unsafe when running with symbol elimination.
    // At the moment the only simplification rules that have multiple conclusions
    // are higher-order and it is assumed that we will not run higher-order along
    // with symbol elimination.
    // In the future if a first-order simplification rule is added with multiple
    // conclusions, this code should be updated.
    onParenthood(replacement, cl);
    while (premStack.isNonEmpty()) {
      onParenthood(replacement, premStack.pop());
    }
  }

  if(!replacement && env.options->proof() == Options::Proof::DEDUKTI)
    Dedukti::unregisterUnit(cl);
}

void SaturationAlgorithm::onNonRedundantClause(Clause* c)
{
  if (_symEl) {
    _symEl->onNonRedundantClause(c);
  }
}

/**
 * Called for clauses derived in the run of the saturation algorithm
 * for each pair clause-premise
 *
 * The propositional parts of clauses may not be set properly (the
 * clauses are always valid, however), also the function is not called
 * for clause merging (when the non-propositional parts would coincide).
 */
void SaturationAlgorithm::onParenthood(Clause* cl, Clause *parent)
{
  if (_symEl) {
    _symEl->onParenthood(cl, parent);
  }
}

/**
 * This function is subscribed to the remove event of the active container
 * instead of the @b onActiveRemoved function in the constructor, as the
 * @b onActiveRemoved function is virtual.
 */
void SaturationAlgorithm::activeRemovedHandler(Clause* cl)
{
  onActiveRemoved(cl);
}

/**
 * This function is subscribed to the remove event of the passive container
 * instead of the @b onPassiveRemoved function in the constructor, as the
 * @b onPassiveRemoved function is virtual.
 */
void SaturationAlgorithm::passiveRemovedHandler(Clause* cl)
{
  onPassiveRemoved(cl);
}

/**
 * Add input clause @b cl into the SaturationAlgorithm object
 *
 * The clause @b cl is added into the unprocessed container, unless the
 * set-of-support option is enabled and @b cl has input type equal to
 * @b Clause::AXIOM. In this case, @b cl is put into the active container.
 */
void SaturationAlgorithm::addInputClause(Clause* cl)
{
  ASS_LE(toNumber(cl->inputType()),toNumber(UnitInputType::CLAIM)); // larger input types should not appear in proof search

  if (_symEl) {
    _symEl->onInputClause(cl);
  }

  bool sosForAxioms = _opt.sos() == Options::Sos::ON || _opt.sos() == Options::Sos::ALL;
  sosForAxioms = sosForAxioms && cl->inputType() == UnitInputType::AXIOM;

  bool sosForTheory = _opt.sos() == Options::Sos::THEORY && _opt.sosTheoryLimit() == 0;

  if (_opt.sineToAge()) {
    unsigned level = cl->getSineLevel();
    // cout << "Adding " << cl->toString() << " level " << level;
    if (level == UINT_MAX) {
      level = env.maxSineLevel - 1; // as the next available (unused) value
      // cout << " -> " << level;
    }
    // cout << endl;
    cl->setAge(level);
  }

  if (sosForAxioms || (cl->isPureTheoryDescendant() && sosForTheory)) {
    addInputSOSClause(cl);
  }
  else {
    addNewClause(cl);
  }

  if (_instantiation) {
    _instantiation->registerClause(cl);
  }

  env.statistics->initialClauses++;
}

/**
 * Return literal selector that is to be used for set-of-support clauses
 */
LiteralSelector& SaturationAlgorithm::getSosLiteralSelector()
{
  if (_opt.sos() == Options::Sos::ALL || _opt.sos() == Options::Sos::THEORY) {
    if (!_sosLiteralSelector) {
      _sosLiteralSelector = new TotalLiteralSelector(getOrdering(), getOptions());
    }
    return *_sosLiteralSelector;
  }
  else {
    return *_selector;
  }
}

/**
 * Add an input set-of-support clause @b cl into the active container
 */
void SaturationAlgorithm::addInputSOSClause(Clause* cl)
{
  ASS_EQ(toNumber(cl->inputType()),toNumber(UnitInputType::AXIOM));

  // we add an extra reference until the clause is added to some container, so that
  // it won't get deleted during some code e.g. in the onNewClause handler
  cl->incRefCnt();

  onNewClause(cl);

simpl_start:

  Clause *simplCl = _immediateSimplifier->simplify(cl);
  if (simplCl != cl) {
    if (!simplCl) {
      onClauseReduction(cl, 0, 0, 0);
      goto fin;
    }

    simplCl->incRefCnt();
    cl->decRefCnt(); // now cl is referenced from simplCl, so after removing the extra reference, it won't be destroyed

    onNewClause(simplCl);
    onClauseReduction(cl, &simplCl, 1, 0);
    cl = simplCl;
    goto simpl_start;
  }

  if (cl->isEmpty()) {
    addNewClause(cl);
    goto fin;
  }

  ASS(!cl->numSelected());
  {
    LiteralSelector& sosSelector = getSosLiteralSelector();
    sosSelector.select(cl);
  }

  cl->setStore(Clause::ACTIVE);
  env.statistics->activeClauses++;
  _active->add(cl);

  onSOSClauseAdded(cl);

fin:
  cl->decRefCnt();
}

/**
 * Insert clauses of the problem into the SaturationAlgorithm object
 * and initialize some internal structures.
 */
void SaturationAlgorithm::init()
{
  ClauseIterator toAdd;

  if (env.options->randomTraversals()) {
    TIME_TRACE(TimeTrace::SHUFFLING);

    Stack<Clause *> aux;
    aux.loadFromIterator(_prb.clauseIterator());
    Shuffling::shuffleArray(aux,aux.size());
    toAdd = pvi(arrayIter(std::move(aux)));
  } else {
    toAdd = _prb.clauseIterator();
  }

  while (toAdd.hasNext()) {
    Clause* cl = toAdd.next();
    addInputClause(cl);
  }

  if (_splitter) {
    _splitter->init(this);
  }
  if (_consFinder) {
    _consFinder->init(this);
  }
  if (_symEl) {
    _symEl->init(this);
  }
}

Clause *SaturationAlgorithm::doImmediateSimplification(Clause* cl0)
{
  TIME_TRACE("immediate simplification");

  static bool sosTheoryLimit = _opt.sos() == Options::Sos::THEORY;
  static unsigned sosTheoryLimitAge = _opt.sosTheoryLimit();
  static ClauseStack repStack;
  repStack.reset();

  SplitSet *splitSet = 0;

  if (sosTheoryLimit && cl0->isPureTheoryDescendant() && cl0->age() > sosTheoryLimitAge) {
    return 0;
  }

  Clause* cl = cl0;

  Clause *simplCl = _immediateSimplifier->simplify(cl);
  if (simplCl != cl) {
    if (simplCl) {
      addNewClause(simplCl);
    }
    onClauseReduction(cl, &simplCl, 1, 0);
    return 0;
  }

  ClauseIterator cIt = _immediateSimplifier->simplifyMany(cl);
  if (cIt.hasNext()) {
    while (cIt.hasNext()) {
      Clause *simpedCl = cIt.next();
      if (!splitSet) {
        splitSet = simpedCl->splits();
      }
      else {
        ASS(splitSet->isSubsetOf(simpedCl->splits()));
        ASS(simpedCl->splits()->isSubsetOf(splitSet));
      }
      ASS(simpedCl != cl);
      repStack.push(simpedCl);
      addNewClause(simpedCl);
    }
    onClauseReduction(cl, repStack.begin(), repStack.size(), 0);
    return 0;
  }

  return cl;
}

/**
 * Add a new clause to the saturation algorithm run
 *
 * At some point of the algorithm loop the @b newClausesToUnprocessed
 * function is called and all new clauses are added to the
 * unprocessed container.
 */
void SaturationAlgorithm::addNewClause(Clause* cl)
{
  if (env.options->randomTraversals()) {
    TIME_TRACE(TimeTrace::SHUFFLING);

    Shuffling::shuffle(cl);
  }

  // we increase the reference counter here so that the clause wouldn't
  // get destroyed during handling in the onNewClause handler
  //(there the control flow goes out of the SaturationAlgorithm class,
  // so we'd better not assume on what's happening out there)
  cl->incRefCnt();
  onNewClause(cl);
  _newClauses.push(cl);
  // we can decrease the counter here -- it won't get deleted because
  // the _newClauses RC stack already took over the clause
  cl->decRefCnt();
}

void SaturationAlgorithm::newClausesToUnprocessed()
{
  if (env.options->randomTraversals()) {
    TIME_TRACE(TimeTrace::SHUFFLING);

    Shuffling::shuffleArray(_newClauses.naked().begin(), _newClauses.size());
  }

  while (_newClauses.isNonEmpty()) {
    Clause* cl = _newClauses.popWithoutDec();
    switch (cl->store()) {
      case Clause::UNPROCESSED:
        break;
      case Clause::PASSIVE:
        onNonRedundantClause(cl);
        break;
      case Clause::NONE:
        addUnprocessedClause(cl);
        break;
      case Clause::SELECTED:
      case Clause::ACTIVE:
#if VDEBUG
        cout << "FAIL: " << cl->toString() << endl;
        // such clauses should not appear as new ones
        cout << cl->toString() << endl;
#endif
        ASSERTION_VIOLATION_REP(cl->store());
    }
    cl->decRefCnt(); // belongs to _newClauses.popWithoutDec()
  }
}

/**
 * Return true iff there are no clauses left to be processed
 *
 * More precisely, true is returned iff the unprocessed clause
 * container and the new clause stack are empty.
 */
bool SaturationAlgorithm::clausesFlushed()
{
  return _unprocessed->isEmpty() && _newClauses.isEmpty();
}

/**
 * Perform immediate simplifications and splitting on clause @b cl and add it
 * to unprocessed.
 *
 * Forward demodulation is also being performed on @b cl.
 */
void SaturationAlgorithm::addUnprocessedClause(Clause* cl)
{
  _generatedClauseCount++;
  env.statistics->generatedClauses++;

  cl=doImmediateSimplification(cl);
  if (!cl) {
    return;
  }

  if (cl->isEmpty()) {
    handleEmptyClause(cl);
    return;
  }

  cl->setStore(Clause::UNPROCESSED);
  _unprocessed->add(cl);
}

/**
 * Deal with clause that has an empty non-propositional part.
 *
 * The function receives a clause @b cl that has empty non-propositional part,
 * and returns a contradiction (an empty clause with false propositional part)
 * if it can be derived from @b cl and previously derived empty clauses.
 * Otherwise it returns 0.
 */
void SaturationAlgorithm::handleEmptyClause(Clause* cl)
{
  ASS(cl->isEmpty());

  if (isRefutation(cl)) {
    onNonRedundantClause(cl);

    throw RefutationFoundException(cl);
  }
  // as Clauses no longer have prop parts the only reason for an empty
  // clause not being a refutation is if it has splits

  if (_splitter && _splitter->handleEmptyClause(cl)) {
    return;
  }

  // splitter should only return false if splits isEmpty, which it cannot be
  ASSERTION_VIOLATION;
  // removed some code that dealt with the case where a clause is empty
  // but as a non-empty bdd prop part
}

/**
 * Forward-simplify the clause @b cl, return true iff the clause
 * should be retained
 *
 * If a weight-limit is imposed on clauses, it is being checked
 * by this function as well.
 */
bool SaturationAlgorithm::forwardSimplify(Clause* cl)
{
  TIME_TRACE("forward simplification");

  if (!_passive->fulfilsAgeLimit(cl) && !_passive->fulfilsWeightLimit(cl)) {
    RSTAT_CTR_INC("clauses discarded by weight limit in forward simplification");
    env.statistics->discardedNonRedundantClauses++;
    return false;
  }

  FwSimplList::Iterator fsit(_fwSimplifiers);

  while (fsit.hasNext()) {
    ForwardSimplificationEngine *fse = fsit.next();

    {
      Clause *replacement = 0;
      ClauseIterator premises = ClauseIterator::getEmpty();

      if (fse->perform(cl, replacement, premises)) {
        if (replacement) {
          addNewClause(replacement);
        }
        onClauseReduction(cl, &replacement, 1, premises);

        return false;
      }
    }
  }

  static ClauseStack repStack;

  repStack.reset();
  SimplList::Iterator sit(_simplifiers);

  while (sit.hasNext()) {
    SimplificationEngine *se = sit.next();

    {
      ClauseIterator results = se->perform(cl);

      if (results.hasNext()) {
        while (results.hasNext()) {
          Clause *simpedCl = results.next();
          ASS(simpedCl != cl);
          repStack.push(simpedCl);
          addNewClause(simpedCl);
        }
        onClauseReduction(cl, repStack.begin(), repStack.size(), 0);
        return false;
      }
    }
  }

  bool synthesis = (env.options->questionAnswering() == Options::QuestionAnsweringMode::SYNTHESIS);

  if (synthesis) {
    ASS((_answerLiteralManager != nullptr));
    Clause* ansLitCl = cl;
    if (_splitter && cl->hasAnswerLiteral() && !cl->noSplits() && cl->computable()) {
      ansLitCl = _splitter->reintroduceAvatarAssertions(cl);
    }
    Clause* reduced = _answerLiteralManager->recordAnswerAndReduce(ansLitCl);
    if (reduced) {
      ansLitCl = reduced;
    }
    if (ansLitCl != cl) {
      addNewClause(ansLitCl);
      onClauseReduction(cl, &ansLitCl, 1, 0);
      return false;
    }
  }

  //TODO: hack that only clauses deleted by forward simplification can be destroyed (other destruction needs debugging)
  cl->incRefCnt();

  if (_splitter && !_opt.splitAtActivation()) {
    if (_splitter->doSplitting(cl)) {
      return false;
    }
  }

  return true;
}

/**
 * The the backward simplification with the clause @b cl.
 */
void SaturationAlgorithm::backwardSimplify(Clause* cl)
{
  TIME_TRACE("backward simplification");

  BwSimplList::Iterator bsit(_bwSimplifiers);
  while (bsit.hasNext()) {
    BackwardSimplificationEngine *bse = bsit.next();

    BwSimplificationRecordIterator simplifications;
    bse->perform(cl, simplifications);
    while (simplifications.hasNext()) {
      BwSimplificationRecord srec = simplifications.next();
      Clause *redundant = srec.toRemove;
      ASS_NEQ(redundant, cl);

      Clause *replacement = srec.replacement;

      if (replacement) {
        addNewClause(replacement);
      }
      onClauseReduction(redundant, &replacement, 1, cl, false);

      // we must remove the redundant clause before adding its replacement,
      // as otherwise the redundant one might demodulate the replacement into
      // a tautology

      redundant->incRefCnt(); // we don't want the clause deleted before we record the simplification

      removeActiveOrPassiveClause(redundant);

      redundant->decRefCnt();
    }
  }
}

/**
 * Remove either passive or active (or reactivated, which is both)
 * clause @b cl
 *
 * In case the removal is requested during clause activation, when some indexes
 * might be traversed (and so cannot be modified), the clause deletion is postponed
 * until the clause activation is over. This is done by pushing the clause on the
 * @b _postponedClauseRemovals stack, which is then checked at the end of the
 * @b activate function.
 */
void SaturationAlgorithm::removeActiveOrPassiveClause(Clause* cl)
{
  if (_clauseActivationInProgress) {
    // we cannot remove clause now, as there indexes might be traversed now,
    // and so we cannot modify them
    _postponedClauseRemovals.push(cl);
    return;
  }

  switch (cl->store()) {
    case Clause::PASSIVE: {
      TIME_TRACE(TimeTrace::PASSIVE_CONTAINER_MAINTENANCE);
      _passive->remove(cl);
      break;
    }
    case Clause::ACTIVE:
      _active->remove(cl);
      break;
    default:
      ASS_REP2(false, cl->store(), *cl);
  }
  // at this point the cl object can be already deleted
}

/**
 * Add clause @b c to the passive container
 */
void SaturationAlgorithm::addToPassive(Clause* cl)
{
  ASS_EQ(cl->store(), Clause::UNPROCESSED);

  cl->setStore(Clause::PASSIVE);
  env.statistics->passiveClauses++;

  {
    TIME_TRACE(TimeTrace::PASSIVE_CONTAINER_MAINTENANCE);
    _passive->add(cl);
  }
}

void SaturationAlgorithm::removeSelected(Clause* cl)
{
  ASS_EQ(cl->store(), Clause::SELECTED);
  beforeSelectedRemoved(cl);
  cl->setStore(Clause::NONE);
}

/**
 * Activate clause @b cl
 *
 * This means putting the clause into the active container, and
 * performing generating inferences with it (in this order, so that
 * inferences such as self-superposition can happen).
 *
 * During clause activation the @b _clauseActivationInProgress value
 * is set to @b true, and clause removals by the @b removeBackwardSimplifiedClause
 * function are postponed. During the clause activation, generalisation
 * indexes should not be modified.
 */
void SaturationAlgorithm::activate(Clause* cl)
{
      TIME_TRACE("activation")

  {
    TIME_TRACE("redundancy check")
    if (_consFinder && _consFinder->isRedundant(cl)) {
      return removeSelected(cl);
    }
  }

  {
    TIME_TRACE("splitting")
    if (_splitter && _opt.splitAtActivation()) {
      if (_splitter->doSplitting(cl)) {
        return removeSelected(cl);
      }
    }
  }

  _clauseActivationInProgress = true;

  if (!cl->numSelected()) {
    TIME_TRACE("clause selection")
    TIME_TRACE("literal selection");

    if (env.options->randomTraversals()) {
      TIME_TRACE(TimeTrace::SHUFFLING);

      Shuffling::shuffle(cl);
    }

    _selector->select(cl);
  }

  ASS_EQ(cl->store(), Clause::SELECTED);
  cl->setStore(Clause::ACTIVE);
  env.statistics->activeClauses++;
  _active->add(cl);

  _conditionalRedundancyHandler->checkEquations(cl);

  auto generated = TIME_TRACE_EXPR(TimeTrace::CLAUSE_GENERATION, _generator->generateSimplify(cl));
  auto toAdd = TIME_TRACE_ITER(TimeTrace::CLAUSE_GENERATION, generated.clauses);

  while (toAdd.hasNext()) {
    Clause *genCl = toAdd.next();
    addNewClause(genCl);

    Inference::Iterator iit = genCl->inference().iterator();
    while (genCl->inference().hasNext(iit)) {
      Unit *premUnit = genCl->inference().next(iit);
      // Now we can get generated clauses having parents that are not clauses
      // Indeed, from induction we can have generated clauses whose parents do
      // not include the activated clause
      if (premUnit->isClause()) {
        Clause *premCl = static_cast<Clause *>(premUnit);
        onParenthood(genCl, premCl);
      }
    }
  }

  _clauseActivationInProgress = false;

  // now we remove clauses that could not be removed during the clause activation process
  if (env.options->randomTraversals()) {
    TIME_TRACE(TimeTrace::SHUFFLING);

    Shuffling::shuffleArray(_postponedClauseRemovals.begin(), _postponedClauseRemovals.size());
  }
  while (_postponedClauseRemovals.isNonEmpty()) {
    Clause* cl = _postponedClauseRemovals.pop();
    if (cl->store() != Clause::ACTIVE && cl->store() != Clause::PASSIVE) {
      continue;
    }
    TIME_TRACE("clause removal")
    removeActiveOrPassiveClause(cl);
  }

  if (generated.premiseRedundant) {
    _active->remove(cl);
  }

  return;
}

/**
 * Perform the loop that puts clauses from the unprocessed to the passive container.
 */
void SaturationAlgorithm::doUnprocessedLoop()
{
start:

  newClausesToUnprocessed();

  while (!_unprocessed->isEmpty()) {
    Clause* c = _unprocessed->pop();
    ASS(!isRefutation(c));

    if (forwardSimplify(c)) {
      onClauseRetained(c);
      addToPassive(c);
      ASS_EQ(c->store(), Clause::PASSIVE);
    }
    else {
      ASS_EQ(c->store(), Clause::UNPROCESSED);
      c->setStore(Clause::NONE);
    }

    newClausesToUnprocessed();
  }

  ASS(clausesFlushed());
  onAllProcessed();
  if (!clausesFlushed()) {
    // there were some new clauses added, so let's process them
    goto start;
  }
}

/**
 * Return true if clause can be passed to activation
 *
 * If false is returned, disposing of the clause is responsibility of
 * this function.
 */
bool SaturationAlgorithm::handleClauseBeforeActivation(Clause* c)
{
  return true;
}

/**
 * This function should be called if (and only if) we will use
 * the @c doOneAlgorithmStep() function to run the saturation
 * algorithm, instead of the @c MailLoop::run() function.
 */
void SaturationAlgorithm::initAlgorithmRun()
{
  init();
}

UnitList *SaturationAlgorithm::collectSaturatedSet()
{
  UnitList* res = 0;
  ClauseIterator it = _active->clauses();
  while (it.hasNext()) {
    Clause* cl = it.next();
    cl->incRefCnt();
    UnitList::push(cl, res);
  }
  return res;
}

/**
 *
 * This function may throw RefutationFoundException and TimeLimitExceededException.
 */
void SaturationAlgorithm::doOneAlgorithmStep()
{
  doUnprocessedLoop();

  if (_passive->isEmpty()) {
    MainLoopResult::TerminationReason termReason =
        isComplete() ? Statistics::SATISFIABLE : Statistics::REFUTATION_NOT_FOUND;
    MainLoopResult res(termReason);

    // if (termReason == Statistics::REFUTATION_NOT_FOUND){
    //   Shell::UIHelper::outputSaturatedSet(cout, pvi(UnitList::Iterator(collectSaturatedSet())));
    // }

    if (termReason == Statistics::SATISFIABLE && getOptions().proof() != Options::Proof::OFF) {
      res.saturatedSet = collectSaturatedSet();

      if (_splitter) {
        res.saturatedSet = _splitter->preprendCurrentlyAssumedComponentClauses(res.saturatedSet);
      }
    }
    throw MainLoopFinishedException(res);
  }

  Clause* cl = nullptr;
  {
    TIME_TRACE(TimeTrace::PASSIVE_CONTAINER_MAINTENANCE);
    cl = _passive->popSelected();
  }
  ASS_EQ(cl->store(), Clause::PASSIVE);
  cl->setStore(Clause::SELECTED);

  if (!handleClauseBeforeActivation(cl)) {
    return;
  }

  activate(cl);
}

/**
 * Perform saturation on clauses that were added through
 * @b addInputClauses function
 */
MainLoopResult SaturationAlgorithm::runImpl()
{
  unsigned l = 0;

  // could be more precise, but we don't care too much
  unsigned startTime = Timer::elapsedDeciseconds();
  try {
    for (;; l++) {
      if (_activationLimit && l > _activationLimit) {
        throw ActivationLimitExceededException();
      }
      if(_softTimeLimit && Timer::elapsedDeciseconds() - startTime > _softTimeLimit)
        throw TimeLimitExceededException();

      doOneAlgorithmStep();
      env.statistics->activations = l;
    }
  }
  catch (ThrowableBase&) {
    tryUpdateFinalClauseCount();
    throw;
  }
}

/**
 * Assign an generating inference object @b generator to be used
 *
 * This object takes ownership of the @b generator object
 * and will be responsible for its deletion.
 *
 * To use multiple generating inferences, use the @b CompositeGIE
 * object.
 */
void SaturationAlgorithm::setGeneratingInferenceEngine(SimplifyingGeneratingInference *generator)
{
  ASS(!_generator);
  _generator = generator;
  _generator->attach(this);
}

/**
 * Assign an immediate simplifier object @b immediateSimplifier
 * to be used
 *
 * This object takes ownership of the @b immediateSimplifier object
 * and will be responsible for its deletion.
 *
 * For description of what an immediate simplifier is, see
 * @b ImmediateSimplificationEngine documentation.
 *
 * To use multiple immediate simplifiers, use the @b CompositeISE
 * object.
 */
void SaturationAlgorithm::setImmediateSimplificationEngine(ImmediateSimplificationEngine *immediateSimplifier)
{
  ASS(!_immediateSimplifier);
  _immediateSimplifier = immediateSimplifier;
  _immediateSimplifier->attach(this);
}

/**
 * Add a forward simplifier, so that it is applied before the
 * simplifiers that were added before it. The object takes ownership
 * of the forward simplifier and will take care of destroying it.
 *
 * Forward demodulation simplifier should be added by the
 * @b setFwDemodulator function, not by this one.
 */
void SaturationAlgorithm::addForwardSimplifierToFront(ForwardSimplificationEngine *fwSimplifier)
{
  FwSimplList::push(fwSimplifier, _fwSimplifiers);
  fwSimplifier->attach(this);
}

void SaturationAlgorithm::addSimplifierToFront(SimplificationEngine *simplifier)
{
  SimplList::push(simplifier, _simplifiers);
  simplifier->attach(this);
}

/**
 * Add a backward simplifier, so that it is applied before the
 * simplifiers that were added before it. The object takes ownership
 * of the backward simplifier and will take care of destroying it.
 */
void SaturationAlgorithm::addBackwardSimplifierToFront(BackwardSimplificationEngine *bwSimplifier)
{
  BwSimplList::push(bwSimplifier, _bwSimplifiers);
  bwSimplifier->attach(this);
}

/**
 * @since 05/05/2013 Manchester, splitting changed to new values
 * @author Andrei Voronkov
 */
SaturationAlgorithm *SaturationAlgorithm::createFromOptions(Problem& prb, const Options& opt, IndexManager *indexMgr)
{
  SaturationAlgorithm* res;
  switch(opt.saturationAlgorithm()) {
  case Shell::Options::SaturationAlgorithm::DISCOUNT:
    res=new Discount(prb, opt);
    break;
  case Shell::Options::SaturationAlgorithm::LRS:
    res=new LRS(prb, opt);
    break;
  case Shell::Options::SaturationAlgorithm::OTTER:
    res=new Otter(prb, opt);
    break;
  default:
    NOT_IMPLEMENTED;
  }
  if (indexMgr) {
    res->_imgr = SmartPtr<IndexManager>(indexMgr, true);
    indexMgr->setSaturationAlgorithm(res);
  }
  else {
    res->_imgr = SmartPtr<IndexManager>(new IndexManager(res));
  }

  if (opt.splitting()) {
    res->_splitter = new Splitter();
  }

  // create generating inference engine
  CompositeGIE *gie = new CompositeGIE();

  if(opt.functionDefinitionIntroduction()) {
    gie->addFront(new DefinitionIntroduction);
  }

  //TODO here induction is last, is that right?
  if(opt.induction()!=Options::Induction::NONE){
    gie->addFront(new Induction());
  }

  if (opt.instantiation() != Options::Instantiation::OFF) {
    res->_instantiation = new Instantiation();
    // res->_instantiation->init();
    gie->addFront(res->_instantiation);
  }

  if (prb.hasEquality()) {
    gie->addFront(new EqualityFactoring());
    gie->addFront(new EqualityResolution());
    if (env.options->superposition()) {
      gie->addFront(new Superposition());
    }
  }
  else if (opt.unificationWithAbstraction() != Options::UnificationWithAbstraction::OFF) {
    gie->addFront(new EqualityResolution());
  }

  if (opt.combinatorySup()) {
    gie->addFront(new ArgCong());
    gie->addFront(new NegativeExt()); // TODO add option
    if (opt.narrow() != Options::Narrow::OFF) {
      gie->addFront(new Narrow());
    }
    if (!opt.pragmatic()) {
      gie->addFront(new SubVarSup());
    }
  }

  if(prb.hasFOOL() &&
    prb.isHigherOrder() && env.options->booleanEqTrick()){
  //  gie->addFront(new ProxyElimination::NOTRemovalGIE());
    gie->addFront(new BoolEqToDiseq());
  }

  if(opt.complexBooleanReasoning() && prb.hasBoolVar() &&
     prb.isHigherOrder() && !opt.lambdaFreeHol()){
    gie->addFront(new PrimitiveInstantiation()); //TODO only add in some cases
    gie->addFront(new ElimLeibniz());
  }

  if (env.options->choiceReasoning()) {
    gie->addFront(new Choice());
  }

  gie->addFront(new Factoring());
  if (opt.binaryResolution()) {
    gie->addFront(new BinaryResolution());
  }
  if (opt.unitResultingResolution() != Options::URResolution::OFF) {
    gie->addFront(new URResolution(opt.unitResultingResolution() == Options::URResolution::FULL));
  }
  if (opt.extensionalityResolution() != Options::ExtensionalityResolution::OFF) {
    gie->addFront(new ExtensionalityResolution());
  }
  if (opt.FOOLParamodulation()) {
    gie->addFront(new FOOLParamodulation());
  }
  if (opt.cases() && prb.hasFOOL() && !opt.casesSimp()) {
    gie->addFront(new Cases());
  }

  if((prb.hasLogicalProxy() || prb.hasBoolVar() || prb.hasFOOL()) &&
      prb.isHigherOrder() && !prb.quantifiesOverPolymorphicVar()){
    if(env.options->cnfOnTheFly() != Options::CNFOnTheFly::EAGER &&
       env.options->cnfOnTheFly() != Options::CNFOnTheFly::OFF){
      gie->addFront(new LazyClausificationGIE());
    }
  }

  if (opt.injectivityReasoning()) {
    gie->addFront(new Injectivity());
  }
  if (prb.hasEquality() && env.signature->hasTermAlgebras()) {
    if (opt.termAlgebraCyclicityCheck() == Options::TACyclicityCheck::RULE) {
      gie->addFront(new AcyclicityGIE());
    }
    else if (opt.termAlgebraCyclicityCheck() == Options::TACyclicityCheck::RULELIGHT) {
      gie->addFront(new AcyclicityGIE1());
    }
    if (opt.termAlgebraInferences()) {
      gie->addFront(new InjectivityGIE());
    }
  }
  if (env.options->functionDefinitionRewriting()) {
    gie->addFront(new FunctionDefinitionRewriting());
    res->addForwardSimplifierToFront(new FunctionDefinitionDemodulation());
  }

  CompositeSGI *sgi = new CompositeSGI();
  sgi->push(gie);

  auto& ordering = res->getOrdering();

  if (opt.evaluationMode() == Options::EvaluationMode::POLYNOMIAL_CAUTIOUS) {
    sgi->push(new PolynomialEvaluation(ordering));
  }

  if (env.options->cancellation() == Options::ArithmeticSimplificationMode::CAUTIOUS) {
    sgi->push(new Cancellation(ordering));
  }

  if (env.options->gaussianVariableElimination() == Options::ArithmeticSimplificationMode::CAUTIOUS) {
    sgi->push(new LfpRule<GaussianVariableElimination>(GaussianVariableElimination()));
  }

  if (env.options->arithmeticSubtermGeneralizations() == Options::ArithmeticSimplificationMode::CAUTIOUS) {
    for (auto gen : allArithmeticSubtermGeneralizations()) {
      sgi->push(gen);
    }
  }

#if VZ3
  if (opt.theoryInstAndSimp() != Shell::Options::TheoryInstSimp::OFF) {
    sgi->push(new TheoryInstAndSimp());
  }
#endif

  res->setGeneratingInferenceEngine(sgi);

  res->setImmediateSimplificationEngine(createISE(prb, opt, res->getOrdering()));

  // create simplification engine

  if((prb.hasLogicalProxy() || prb.hasBoolVar() || prb.hasFOOL()) &&
      prb.isHigherOrder() && !prb.quantifiesOverPolymorphicVar()){
    if(env.options->cnfOnTheFly() != Options::CNFOnTheFly::EAGER &&
       env.options->cnfOnTheFly() != Options::CNFOnTheFly::OFF){
      res->addSimplifierToFront(new LazyClausification());
    }
  }

  // create forward simplification engine
  if (prb.hasEquality() && opt.innerRewriting()) {
    res->addForwardSimplifierToFront(new InnerRewriting());
  }
  if (opt.globalSubsumption()) {
    res->addForwardSimplifierToFront(new GlobalSubsumption(opt));
  }
  if (opt.forwardLiteralRewriting()) {
    res->addForwardSimplifierToFront(new ForwardLiteralRewriting());
  }
  if (prb.hasEquality()) {
    // NOTE:
    // fsd should be performed after forward subsumption,
    // because every successful forward subsumption will lead to a (useless) match in fsd.
    if (opt.forwardSubsumptionDemodulation()) {
      res->addForwardSimplifierToFront(new ForwardSubsumptionDemodulation(false));
    }
  }
  if (prb.hasEquality()) {
    switch (opt.forwardDemodulation()) {
      case Options::Demodulation::ALL:
      case Options::Demodulation::PREORDERED:
        if (opt.combinatorySup()) {
          res->addForwardSimplifierToFront(new ForwardDemodulationImpl<true>());
        }
        else {
          res->addForwardSimplifierToFront(new ForwardDemodulationImpl<false>());
        }
        break;
      case Options::Demodulation::OFF:
        break;
#if VDEBUG
      default:
        ASSERTION_VIOLATION;
#endif
    }
  }

  if (opt.forwardSubsumption()) {
    if (opt.codeTreeSubsumption()) {
      res->addForwardSimplifierToFront(new CodeTreeForwardSubsumptionAndResolution(opt.forwardSubsumptionResolution()));
    } else {
      res->addForwardSimplifierToFront(new ForwardSubsumptionAndResolution(opt.forwardSubsumptionResolution()));
    }
  }
  else if (opt.forwardSubsumptionResolution()) {
    USER_ERROR("Forward subsumption resolution requires forward subsumption to be enabled.");
  }

  // create backward simplification engine
  if (prb.hasEquality()) {
    switch (opt.backwardDemodulation()) {
      case Options::Demodulation::ALL:
      case Options::Demodulation::PREORDERED:
        res->addBackwardSimplifierToFront(new BackwardDemodulation());
        break;
      case Options::Demodulation::OFF:
        break;
#if VDEBUG
      default:
        ASSERTION_VIOLATION;
#endif
    }
  }
  if (prb.hasEquality() && opt.backwardSubsumptionDemodulation()) {
    res->addBackwardSimplifierToFront(new BackwardSubsumptionDemodulation());
  }

  bool backSubsumption = opt.backwardSubsumption() != Options::Subsumption::OFF;
  bool backSR = opt.backwardSubsumptionResolution() != Options::Subsumption::OFF;
  bool subsumptionUnitOnly = opt.backwardSubsumption() == Options::Subsumption::UNIT_ONLY;
  bool srUnitOnly = opt.backwardSubsumptionResolution() == Options::Subsumption::UNIT_ONLY;
  if (backSubsumption || backSR) {
    res->addBackwardSimplifierToFront(new BackwardSubsumptionAndResolution(backSubsumption, subsumptionUnitOnly, backSR, srUnitOnly));
  }

  if (opt.mode() == Options::Mode::CONSEQUENCE_ELIMINATION) {
    res->_consFinder = new ConsequenceFinder();
  }
  if (opt.showSymbolElimination()) {
    res->_symEl = new SymElOutput();
  }

  res->_conditionalRedundancyHandler.reset(ConditionalRedundancyHandler::create(opt, &ordering, res->_splitter));

  res->_answerLiteralManager = AnswerLiteralManager::getInstance(); // selects the right one, according to options!
  ASS(!res->_answerLiteralManager||opt.questionAnswering()!=Options::QuestionAnsweringMode::OFF);
  ASS( res->_answerLiteralManager||opt.questionAnswering()==Options::QuestionAnsweringMode::OFF);
  return res;
} // SaturationAlgorithm::createFromOptions

/**
 * Create local clause simplifier for problem @c prb according to options @c opt
 */
ImmediateSimplificationEngine *SaturationAlgorithm::createISE(Problem& prb, const Options& opt, Ordering& ordering)
{
  CompositeISE* res=new CompositeISE();

  if (prb.hasEquality() && opt.equationalTautologyRemoval()) {
    res->addFront(new EquationalTautologyRemoval());
  }

  switch (opt.condensation()) {
    case Options::Condensation::ON:
      res->addFront(new Condensation());
      break;
    case Options::Condensation::FAST:
      res->addFront(new FastCondensation());
      break;
    case Options::Condensation::OFF:
      break;
  }

  if (env.options->combinatorySup()) {
    res->addFront(new CombinatorDemodISE());
    res->addFront(new CombinatorNormalisationISE());
  }

  if (env.options->choiceReasoning()) {
    res->addFront(new ChoiceDefinitionISE());
  }

  if((prb.hasLogicalProxy() || prb.hasBoolVar() || prb.hasFOOL()) &&
      prb.isHigherOrder() && !env.options->addProxyAxioms()){
    if(env.options->cnfOnTheFly() == Options::CNFOnTheFly::EAGER){
      /*res->addFrontMany(new ProxyISE());
      res->addFront(new OrImpAndProxyISE());
      res->addFront(new NotProxyISE());
      res->addFront(new EqualsProxyISE());
      res->addFront(new PiSigmaProxyISE());*/
      res->addFrontMany(new EagerClausificationISE());
    }
    else {
      res->addFront(new IFFXORRewriterISE());
    }
    res->addFront(new BoolSimp());
  }

  if (prb.hasFOOL() && opt.casesSimp() && !opt.cases()) {
    res->addFrontMany(new CasesSimp());
  }

  // Only add if there are distinct groups
  if (prb.hasEquality() && env.signature->hasDistinctGroups()) {
    res->addFront(new DistinctEqualitySimplifier());
  }
  if (prb.hasEquality() && env.signature->hasTermAlgebras()) {
    if (opt.termAlgebraInferences()) {
      res->addFront(new DistinctnessISE());
      res->addFront(new InjectivityISE());
      res->addFront(new NegativeInjectivityISE());
    }
  }
  if (prb.hasInterpretedOperations() || prb.hasNumerals()) {
    if (env.options->arithmeticSubtermGeneralizations() == Options::ArithmeticSimplificationMode::FORCE) {
      for (auto gen : allArithmeticSubtermGeneralizations()) {
        res->addFront(&gen->asISE());
      }
    }

    if (env.options->gaussianVariableElimination() == Options::ArithmeticSimplificationMode::FORCE) {
      res->addFront(&(new GaussianVariableElimination())->asISE());
    }

    if (env.options->cancellation() == Options::ArithmeticSimplificationMode::FORCE) {
      res->addFront(&(new Cancellation(ordering))->asISE());
    }

    switch (env.options->evaluationMode()) {
      case Options::EvaluationMode::SIMPLE:
        res->addFront(new InterpretedEvaluation(env.options->inequalityNormalization(), ordering));
        break;
      case Options::EvaluationMode::POLYNOMIAL_FORCE:
        res->addFront(&(new PolynomialEvaluation(ordering))->asISE());
        break;
      case Options::EvaluationMode::POLYNOMIAL_CAUTIOUS:
        break;
    }

    if (env.options->pushUnaryMinus()) {
      res->addFront(new PushUnaryMinus());
    }
  }
  if (prb.hasEquality()) {
    res->addFront(new TrivialInequalitiesRemovalISE());
  }
  res->addFront(new TautologyDeletionISE());
  if (env.options->newTautologyDel()) {
    res->addFront(new TautologyDeletionISE2());
  }
  res->addFront(new DuplicateLiteralRemovalISE());

  if (env.options->questionAnswering() == Options::QuestionAnsweringMode::PLAIN) {
    res->addFront(new AnswerLiteralResolver());
    if (env.options->questionAnsweringAvoidThese() != "") {
      res->addFront(new UndesiredAnswerLiteralRemoval(env.options->questionAnsweringAvoidThese()));
    }
  } else if (env.options->questionAnswering() == Options::QuestionAnsweringMode::SYNTHESIS) {
    res->addFront(new UncomputableAnswerLiteralRemoval());
    res->addFront(new MultipleAnswerLiteralRemoval());
  }
  return res;
}
