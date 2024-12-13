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
 * @file Statistics.cpp
 * Implements proof-search statistics
 *
 * @since 02/01/2008 Manchester
 */

#include <iostream>

#include "Debug/RuntimeStatistics.hpp"

#include "Lib/Environment.hpp"
#include "Lib/Timer.hpp"
#include "SAT/Z3Interfacing.hpp"

#include "Shell/UIHelper.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

#include "Options.hpp"
#include "Statistics.hpp"


using namespace std;
using namespace Lib;
using namespace Saturation;
using namespace Shell;

/**
 * Initialise statistics.
 * @since 02/01/2008 Manchester
 */
Statistics::Statistics()
  : inputClauses(0),
    inputFormulas(0),
    formulaNames(0),
    initialClauses(0),
    splitInequalities(0),
    purePredicates(0),
    trivialPredicates(0),
    unusedPredicateDefinitions(0),
    functionDefinitions(0),
    selectedBySine(0),
    sineIterations(0),
    factoring(0),
    resolution(0),
    urResolution(0),
    cResolution(0),
    forwardSuperposition(0),
    backwardSuperposition(0),
    cSelfSuperposition(0),
    cForwardSuperposition(0),
    cBackwardSuperposition(0),
    selfSuperposition(0),
    equalityFactoring(0),
    equalityResolution(0),
    forwardExtensionalityResolution(0),
    backwardExtensionalityResolution(0),
    theoryInstSimp(0),
    theoryInstSimpCandidates(0),
    theoryInstSimpTautologies(0),
    theoryInstSimpLostSolution(0),
    theoryInstSimpEmptySubstitution(0),
    maxInductionDepth(0),
    structInduction(0),
    structInductionInProof(0),
    intInfInduction(0),
    intInfInductionInProof(0),
    intFinInduction(0),
    intFinInductionInProof(0),
    intDBInduction(0),
    intDBInductionInProof(0),
    intInfUpInduction(0),
    intInfUpInductionInProof(0),
    intFinUpInduction(0),
    intFinUpInductionInProof(0),
    intDBUpInduction(0),
    intDBUpInductionInProof(0),
    intInfDownInduction(0),
    intInfDownInductionInProof(0),
    intFinDownInduction(0),
    intFinDownInductionInProof(0),
    intDBDownInduction(0),
    intDBDownInductionInProof(0),
    inductionApplication(0),
    inductionApplicationInProof(0),
    generalizedInductionApplication(0),
    generalizedInductionApplicationInProof(0),
    argumentCongruence(0),
    narrow(0),
    forwardSubVarSup(0),
    backwardSubVarSup(0),
    selfSubVarSup(0),
    negativeExtensionality(0),
    primitiveInstantiations(0),
    choiceInstances(0),
    proxyEliminations(0),
    leibnizElims(0),
    booleanSimps(0),
    skippedSuperposition(0),
    skippedResolution(0),
    skippedEqualityResolution(0),
    skippedEqualityFactoring(0),
    skippedFactoring(0),
    skippedInferencesDueToOrderingConstraints(0),
    skippedInferencesDueToAvatarConstraints(0),
    skippedInferencesDueToLiteralConstraints(0),
    duplicateLiterals(0),
    trivialInequalities(0),
    forwardSubsumptionResolution(0),
    backwardSubsumptionResolution(0),
    forwardDemodulations(0),
    forwardDemodulationsToEqTaut(0),
    backwardDemodulations(0),
    backwardDemodulationsToEqTaut(0),
    forwardSubsumptionDemodulations(0),
    forwardSubsumptionDemodulationsToEqTaut(0),
    backwardSubsumptionDemodulations(0),
    backwardSubsumptionDemodulationsToEqTaut(0),
    forwardLiteralRewrites(0),
    condensations(0),
    globalSubsumption(0),
    interpretedSimplifications(0),

    asgViolations(0),
    asgCnt(0),

    gveViolations(0),
    gveCnt(0),

    evaluationIncomp(0),
    evaluationGreater(0),
    evaluationCnt(0),
    innerRewrites(0),
    innerRewritesToEqTaut(0),
    deepEquationalTautologies(0),
    simpleTautologies(0),
    equationalTautologies(0),
    forwardSubsumed(0),
    backwardSubsumed(0),
    taDistinctnessSimplifications(0),
    taDistinctnessTautologyDeletions(0),
    taInjectivitySimplifications(0),
    taNegativeInjectivitySimplifications(0),
    taAcyclicityGeneratedDisequalities(0),
    generatedClauses(0),
    passiveClauses(0),
    activeClauses(0),
    extensionalityClauses(0),
    discardedNonRedundantClauses(0),
    inferencesBlockedForOrderingAftercheck(0),
    smtReturnedUnknown(false),
    smtDidNotEvaluate(false),
    inferencesSkippedDueToColors(0),
    finalPassiveClauses(0),
    finalActiveClauses(0),
    finalExtensionalityClauses(0),
    splitClauses(0),
    splitComponents(0),
    uniqueComponents(0),
    satClauses(0),
    unitSatClauses(0),
    binarySatClauses(0),

    satSplits(0),
    satSplitRefutations(0),

    smtFallbacks(0),

    satPureVarsEliminated(0),
    terminationReason(UNKNOWN),
    refutation(0),
    saturatedSet(0),
    phase(INITIALIZATION)
{
} // Statistics::Statistics

void Statistics::explainRefutationNotFound(std::ostream& out)
{
  // should be a one-liner for each case!
  if (discardedNonRedundantClauses) {
    out << "Refutation not found, non-redundant clauses discarded";
  }
  else if (inferencesSkippedDueToColors) {
    out << "Refutation not found, inferences skipped due to colors\n";
  }
  else if(smtReturnedUnknown){
    out << "Refutation not found, SMT solver inside AVATAR returned Unknown";
  }
  else if (smtDidNotEvaluate) {
    out << "Refutation not found, SMT solver inside AVATAR failed to evaluate a literal\n";
  }
  else {
    out << "Refutation not found, incomplete strategy";
  }
}

void Statistics::print(std::ostream& out)
{
  if (env.options->statistics() != Options::Statistics::NONE) {

  SaturationAlgorithm::tryUpdateFinalClauseCount();

  bool separable=false;
#define HEADING(text,num) if (num) { addCommentSignForSZS(out); out << ">>> " << (text) << endl;}
#define COND_OUT(text, num) if (num) { addCommentSignForSZS(out); out << (text) << ": " << (num) << endl; separable = true; }
#define SEPARATOR if (separable) { addCommentSignForSZS(out); out << endl; separable = false; }

  addCommentSignForSZS(out);
  out << "------------------------------\n";
  addCommentSignForSZS(out);
  out << "Version: " << VERSION_STRING << endl;
#if VZ3
  addCommentSignForSZS(out);
  out << "Linked with Z3 " << Z3Interfacing::z3_full_version() << endl;
#endif

  addCommentSignForSZS(out);
  out << "Termination reason: ";
  switch(terminationReason) {
  case Statistics::REFUTATION:
    out << "Refutation";
    break;
  case Statistics::TIME_LIMIT:
    out << "Time limit";
    break;
  case Statistics::MEMORY_LIMIT:
    out << "Memory limit";
    break;
  case Statistics::ACTIVATION_LIMIT:
    out << "Activation limit";
    break;
  case Statistics::REFUTATION_NOT_FOUND:
    explainRefutationNotFound(out);
    break;
  case Statistics::SATISFIABLE:
    out << "Satisfiable";
    break;
  case Statistics::UNKNOWN:
    out << "Unknown";
    break;
  case Statistics::INAPPROPRIATE:
    out << "Inappropriate";
    break;
  default:
    ASSERTION_VIOLATION;
  }
  out << endl;
  if (phase!=FINALIZATION) {
    addCommentSignForSZS(out);
    out << "Termination phase: " << phaseToString(phase) << endl;
  }
  out << endl;

  if (env.options->statistics()==Options::Statistics::FULL) {

  HEADING("Input",inputClauses+inputFormulas);
  COND_OUT("Input clauses", inputClauses);
  COND_OUT("Input formulas", inputFormulas);

  HEADING("Preprocessing",formulaNames+purePredicates+trivialPredicates+
    unusedPredicateDefinitions+functionDefinitions+selectedBySine+
    sineIterations+splitInequalities);
  COND_OUT("Introduced names",formulaNames);
  COND_OUT("Introduced skolems",skolemFunctions);
  COND_OUT("Pure predicates", purePredicates);
  COND_OUT("Trivial predicates", trivialPredicates);
  COND_OUT("Unused predicate definitions", unusedPredicateDefinitions);
  COND_OUT("Function definitions", functionDefinitions);
  COND_OUT("Selected by SInE selection", selectedBySine);
  COND_OUT("SInE iterations", sineIterations);
  COND_OUT("Split inequalities", splitInequalities);
  SEPARATOR;

  HEADING("Saturation",activeClauses+passiveClauses+extensionalityClauses+
      generatedClauses+finalActiveClauses+finalPassiveClauses+finalExtensionalityClauses+
      discardedNonRedundantClauses+inferencesSkippedDueToColors+inferencesBlockedForOrderingAftercheck);
  COND_OUT("Initial clauses", initialClauses);
  COND_OUT("Generated clauses", generatedClauses);
  COND_OUT("Activations started", activations);
  COND_OUT("Active clauses", activeClauses);
  COND_OUT("Passive clauses", passiveClauses);
  COND_OUT("Extensionality clauses", extensionalityClauses);
  COND_OUT("Blocked clauses", blockedClauses);
  COND_OUT("Final active clauses", finalActiveClauses);
  COND_OUT("Final passive clauses", finalPassiveClauses);
  COND_OUT("Final extensionality clauses", finalExtensionalityClauses);
  COND_OUT("Discarded non-redundant clauses", discardedNonRedundantClauses);
  COND_OUT("Inferences skipped due to colors", inferencesSkippedDueToColors);
  COND_OUT("Inferences blocked due to ordering aftercheck", inferencesBlockedForOrderingAftercheck);
  SEPARATOR;


  HEADING("Simplifying Inferences",duplicateLiterals+trivialInequalities+
      forwardSubsumptionResolution+backwardSubsumptionResolution+proxyEliminations+
      forwardDemodulations+backwardDemodulations+forwardLiteralRewrites+
      forwardSubsumptionDemodulations+backwardSubsumptionDemodulations+
      condensations+globalSubsumption+evaluationCnt
      +( gveCnt - gveViolations)
      +( asgCnt - asgViolations)
      +( evaluationCnt - evaluationIncomp - evaluationGreater)
      +innerRewrites
      +booleanSimps
      );
  COND_OUT("Duplicate literals", duplicateLiterals);
  COND_OUT("Trivial inequalities", trivialInequalities);
  COND_OUT("Fw subsumption resolutions", forwardSubsumptionResolution);
  COND_OUT("Bw subsumption resolutions", backwardSubsumptionResolution);
  COND_OUT("Fw demodulations", forwardDemodulations);
  COND_OUT("Bw demodulations", backwardDemodulations);
  COND_OUT("Fw subsumption demodulations", forwardSubsumptionDemodulations);
  COND_OUT("Bw subsumption demodulations", backwardSubsumptionDemodulations);
  COND_OUT("Fw literal rewrites", forwardLiteralRewrites);
  COND_OUT("Inner rewrites", innerRewrites);
  COND_OUT("Condensations", condensations);
  COND_OUT("Global subsumptions", globalSubsumption);
  COND_OUT("Interpreted simplifications", interpretedSimplifications);

  COND_OUT("asg count", asgCnt);
  COND_OUT("asg results not smaller than the premis", asgViolations);

  COND_OUT("gve count", gveCnt);
  COND_OUT("gve results not smaller than the premis", gveViolations);

  COND_OUT("Evaluation count",         evaluationCnt);
  COND_OUT("Evaluation results greater than premise", evaluationGreater);
  COND_OUT("Evaluation results incomparable to premise", evaluationIncomp);
  COND_OUT("Logicial proxy rewrites", proxyEliminations);
  COND_OUT("Boolean simplifications", booleanSimps)
  //COND_OUT("Interpreted simplifications", interpretedSimplifications);
  SEPARATOR;

  HEADING("Deletion Inferences",simpleTautologies+equationalTautologies+
      forwardSubsumed+backwardSubsumed+forwardDemodulationsToEqTaut+
      forwardSubsumptionDemodulationsToEqTaut+backwardSubsumptionDemodulationsToEqTaut+
      backwardDemodulationsToEqTaut+innerRewritesToEqTaut);
  COND_OUT("Simple tautologies", simpleTautologies);
  COND_OUT("Equational tautologies", equationalTautologies);
  COND_OUT("Deep equational tautologies", deepEquationalTautologies);
  COND_OUT("Forward subsumptions", forwardSubsumed);
  COND_OUT("Backward subsumptions", backwardSubsumed);
  COND_OUT("Fw demodulations to eq. taut.", forwardDemodulationsToEqTaut);
  COND_OUT("Bw demodulations to eq. taut.", backwardDemodulationsToEqTaut);
  COND_OUT("Fw subsumption demodulations to eq. taut.", forwardSubsumptionDemodulationsToEqTaut);
  COND_OUT("Bw subsumption demodulations to eq. taut.", backwardSubsumptionDemodulationsToEqTaut);
  COND_OUT("Inner rewrites to eq. taut.", innerRewritesToEqTaut);
  SEPARATOR;

  HEADING("Generating Inferences",resolution+urResolution+cResolution+factoring+
      forwardSuperposition+backwardSuperposition+selfSuperposition+
      cForwardSuperposition+cBackwardSuperposition+cSelfSuperposition+leibnizElims+
      equalityFactoring+equalityResolution+forwardExtensionalityResolution+
      backwardExtensionalityResolution+argumentCongruence+negativeExtensionality+
      +primitiveInstantiations+choiceInstances+narrow+forwardSubVarSup+backwardSubVarSup+selfSubVarSup+
      theoryInstSimp+theoryInstSimpCandidates+theoryInstSimpTautologies+theoryInstSimpLostSolution+inductionApplication+generalizedInductionApplication);
  COND_OUT("Binary resolution", resolution);
  COND_OUT("Unit resulting resolution", urResolution);
  COND_OUT("Binary resolution with abstraction",cResolution);
  COND_OUT("Factoring", factoring);
  COND_OUT("Forward superposition", forwardSuperposition);
  COND_OUT("Backward superposition", backwardSuperposition);
  COND_OUT("Self superposition", selfSuperposition);
  COND_OUT("Forward superposition with abstraction", cForwardSuperposition);
  COND_OUT("Backward superposition with abstraction", cBackwardSuperposition);
  COND_OUT("Self superposition with abstraction", cSelfSuperposition);
  COND_OUT("Equality factoring", equalityFactoring);
  COND_OUT("Equality resolution", equalityResolution);
  COND_OUT("Fw extensionality resolution", forwardExtensionalityResolution);
  COND_OUT("Bw extensionality resolution", backwardExtensionalityResolution);
  COND_OUT("TheoryInstSimp",theoryInstSimp);
  COND_OUT("TheoryInstSimpCandidates",theoryInstSimpCandidates);
  COND_OUT("TheoryInstSimpTautologies",theoryInstSimpTautologies);
  COND_OUT("TheoryInstSimpLostSolution",theoryInstSimpLostSolution);
  COND_OUT("TheoryInstSimpEmptySubstitutions",theoryInstSimpEmptySubstitution);
  COND_OUT("MaxInductionDepth",maxInductionDepth);
  COND_OUT("StructuralInduction",structInduction);
  COND_OUT("StructuralInductionInProof",structInductionInProof);
  COND_OUT("IntegerInfiniteIntervalInduction",intInfInduction);
  COND_OUT("IntegerInfiniteIntervalInductionInProof",intInfInductionInProof);
  COND_OUT("IntegerFiniteIntervalInduction",intFinInduction);
  COND_OUT("IntegerFiniteIntervalInductionInProof",intFinInductionInProof);
  COND_OUT("IntegerDefaultBoundInduction",intDBInduction);
  COND_OUT("IntegerDefaultBoundInductionInProof",intDBInductionInProof);
  COND_OUT("IntegerInfiniteIntervalUpInduction",intInfUpInduction);
  COND_OUT("IntegerInfiniteIntervalUpInductionInProof",intInfUpInductionInProof);
  COND_OUT("IntegerFiniteIntervalUpInduction",intFinUpInduction);
  COND_OUT("IntegerFiniteIntervalUpInductionInProof",intFinUpInductionInProof);
  COND_OUT("IntegerDefaultBoundUpInduction",intDBUpInduction);
  COND_OUT("IntegerDefaultBoundUpInductionInProof",intDBUpInductionInProof);
  COND_OUT("IntegerInfiniteIntervalDownInduction",intInfDownInduction);
  COND_OUT("IntegerInfiniteIntervalDownInductionInProof",intInfDownInductionInProof);
  COND_OUT("IntegerFiniteIntervalDownInduction",intFinDownInduction);
  COND_OUT("IntegerFiniteIntervalDownInductionInProof",intFinDownInductionInProof);
  COND_OUT("IntegerDefaultBoundDownInduction",intDBDownInduction);
  COND_OUT("IntegerDefaultBoundDownInductionInProof",intDBDownInductionInProof);
  COND_OUT("InductionApplications",inductionApplication);
  COND_OUT("InductionApplicationsInProof",inductionApplicationInProof);
  COND_OUT("GeneralizedInductionApplications",generalizedInductionApplication);
  COND_OUT("GeneralizedInductionApplicationsInProof",generalizedInductionApplicationInProof);
  COND_OUT("Argument congruence", argumentCongruence);
  COND_OUT("Negative extensionality", negativeExtensionality);
  COND_OUT("Primitive substitutions", primitiveInstantiations);
  COND_OUT("Elimination of Leibniz equalities", leibnizElims);
  COND_OUT("Choice axiom instances creatded", choiceInstances);
  COND_OUT("Narrow", narrow);
  COND_OUT("Forward sub-variable superposition", forwardSubVarSup);
  COND_OUT("Backward sub-variable superposition", backwardSubVarSup);
  COND_OUT("Self sub-variable superposition", selfSubVarSup);
  SEPARATOR;

  HEADING("Redundant Inferences",
    skippedSuperposition+skippedResolution+skippedEqualityResolution+skippedEqualityFactoring+
    skippedFactoring+skippedInferencesDueToOrderingConstraints+
    skippedInferencesDueToAvatarConstraints+skippedInferencesDueToLiteralConstraints);
  COND_OUT("Skipped superposition", skippedSuperposition);
  COND_OUT("Skipped resolution", skippedResolution);
  COND_OUT("Skipped equality resolution", skippedEqualityResolution);
  COND_OUT("Skipped equality factoring", skippedEqualityFactoring);
  COND_OUT("Skipped factoring", skippedFactoring);
  COND_OUT("Skipped inferences due to ordering constraints", skippedInferencesDueToOrderingConstraints);
  COND_OUT("Skipped inferences due to AVATAR constraints", skippedInferencesDueToAvatarConstraints);
  COND_OUT("Skipped inferences due to literal constraints", skippedInferencesDueToLiteralConstraints);
  SEPARATOR;

  HEADING("Term algebra simplifications",taDistinctnessSimplifications+
      taDistinctnessTautologyDeletions+taInjectivitySimplifications+
      taAcyclicityGeneratedDisequalities+taNegativeInjectivitySimplifications);
  COND_OUT("Distinctness simplifications",taDistinctnessSimplifications);
  COND_OUT("Distinctness tautology deletions",taDistinctnessTautologyDeletions);
  COND_OUT("Injectivity simplifications",taInjectivitySimplifications);
  COND_OUT("Negative injectivity simplifications",taNegativeInjectivitySimplifications);
  COND_OUT("Disequalities generated from acyclicity",taAcyclicityGeneratedDisequalities);

  HEADING("AVATAR",splitClauses+splitComponents+uniqueComponents+satSplits+
        satSplitRefutations);
  COND_OUT("Split clauses", splitClauses);
  COND_OUT("Split components", splitComponents);
  COND_OUT("Unique components", uniqueComponents);
  //COND_OUT("Sat splits", satSplits); // same as split clauses
  COND_OUT("Sat splitting refutations", satSplitRefutations);
  COND_OUT("SMT fallbacks",smtFallbacks);
  SEPARATOR;

  //TODO record statistics for FMB

  //TODO record statistics for MiniSAT
  HEADING("SAT Solver Statistics",satClauses+unitSatClauses+binarySatClauses+satPureVarsEliminated);
  COND_OUT("SAT solver clauses", satClauses);
  COND_OUT("SAT solver unit clauses", unitSatClauses);
  COND_OUT("SAT solver binary clauses", binarySatClauses);
  COND_OUT("Pure propositional variables eliminated by SAT solver", satPureVarsEliminated);
  SEPARATOR;

  }

  addCommentSignForSZS(out);
  out << "Time elapsed: ";
  Timer::printMSString(out,Timer::elapsedMilliseconds());
  out << endl;

  Timer::updateInstructionCount();
  unsigned instr = Timer::elapsedMegaInstructions();
  if (instr) {
    addCommentSignForSZS(out);
    out << "Instructions burned: " << instr << " (million)";
    out << endl;
  }

  addCommentSignForSZS(out);
  out << "------------------------------\n";

  RSTAT_PRINT(out);
  addCommentSignForSZS(out);
  out << "------------------------------\n";

#undef SEPARATOR
#undef COND_OUT
  } // if (env.options->statistics()!=Options::Statistics::NONE)

#if VTIME_PROFILING
  if (env.options && env.options->timeStatistics()) {
    TimeTrace::instance().printPretty(out);
  }
#endif // VTIME_PROFILING
}

const char* Statistics::phaseToString(ExecutionPhase p)
{
  switch(p) {
  case INITIALIZATION:
    return "Initialization";
  case PARSING:
    return "Parsing";
  case PROPERTY_SCANNING:
    return "Property scanning";
  case NORMALIZATION:
    return "Normalization";
  case SHUFFLING:
    return "shuffling";
  case SINE_SELECTION:
    return "SInE selection";
  case INCLUDING_THEORY_AXIOMS:
    return "Including theory axioms";
  case PREPROCESS_1:
    return "Preprocessing 1";
  case PREDIACTE_DEFINITION_MERGING:
    return "Predicate definition merging";
  case PREDICATE_DEFINITION_INLINING:
    return "Predicate definition inlining";
  case UNUSED_PREDICATE_DEFINITION_REMOVAL:
    return "Unused predicate definition removal";
  case BLOCKED_CLAUSE_ELIMINATION:
    return "Blocked clause elimination";
  case TWEE:
    return "Twee Goal Transformation";
  case ANSWER_LITERAL: 
    return "Answer literal addition";
  case PREPROCESS_2:
    return "Preprocessing 2";
  case NEW_CNF:
    return "NewCNF";
  case NAMING:
    return "Naming";
  case PREPROCESS_3:
    return "Preprocessing 3";
  case CLAUSIFICATION:
    return "Clausification";
  case FUNCTION_DEFINITION_ELIMINATION:
    return "Function definition elimination";
  case INEQUALITY_SPLITTING:
    return "Inequality splitting";
  case EQUALITY_RESOLUTION_WITH_DELETION:
    return "Equality resolution with deletion";
  case EQUALITY_PROXY:
    return "Equality proxy";
  case GENERAL_SPLITTING:
    return "General splitting";
  case SATURATION:
    return "Saturation";
  case FINALIZATION:
    return "Finalization";
  case UNKNOWN_PHASE:
    return "Unknown";
  case PREPROCESSING: 
    return "BP Preprocessing ";
  case SOLVING:
    return "Solving with Conflict Resolution";
  case SAT_SOLVING:
	  return "SAT Solving";
  case FMB_PREPROCESSING:
          return "Finite model building preprocessing";
  case FMB_CONSTRAINT_GEN:
          return "Finite model building constraint generation";
  case FMB_SOLVING:
          return "Finite model building SAT solving";
  default:
    ASSERTION_VIOLATION;
    return "Invalid ExecutionPhase value";
  }
}


