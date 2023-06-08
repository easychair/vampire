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
 * @file UIHelper.cpp
 * Implements class UIHelper.
 */

#include <fstream>

#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <sstream>

#include "Forwards.hpp"

#include "Lib/Environment.hpp"
#include "Debug/TimeProfiling.hpp"
#include "Lib/Allocator.hpp"
#include "Lib/ScopedLet.hpp"
#include "Lib/Timer.hpp"

#include "Kernel/InferenceStore.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/FormulaUnit.hpp"

#include "Parse/SMTLIB2.hpp"
#include "Parse/TPTP.hpp"

#include "AnswerLiteralManager.hpp"
#include "Dedukti.hpp"
#include "InterpolantMinimizer.hpp"
#include "Interpolants.hpp"
#include "LaTeX.hpp"
#include "LispLexer.hpp"
#include "LispParser.hpp"
#include "Options.hpp"
#include "Statistics.hpp"
#include "TPTPPrinter.hpp"
#include "UIHelper.hpp"

#include "SAT/Z3Interfacing.hpp"

#include "Lib/List.hpp"
#include "Lib/ScopedPtr.hpp"

namespace Shell {

using namespace Lib;
using namespace Kernel;
using namespace Saturation;
using namespace std;

bool outputAllowed(bool debug)
{
#if VDEBUG
  if (debug) { return true; }
#endif

  // spider and smtcomp output modes are generally silent
  return !Lib::env.options ||
    (
     Lib::env.options->outputMode() != Shell::Options::Output::SPIDER
     && Lib::env.options->outputMode() != Shell::Options::Output::SMTCOMP
     && Lib::env.options->outputMode() != Shell::Options::Output::UCORE
     );
}

void reportSpiderFail()
{
  reportSpiderStatus('!');
}

void reportSpiderStatus(char status)
{
#if VZ3
  if (UIHelper::spiderOutputDone) {
    return;
  }
  if (Lib::env.options && Lib::env.options->outputMode() != Shell::Options::Output::SPIDER) {
    return;
  }

  UIHelper::spiderOutputDone = true;

  // compute Vampire Z3 version and commit
  std::string version = VERSION_STRING;
  size_t versionPosition = version.find("commit ") + strlen("commit ");
  size_t afterVersionPosition = version.find(" ",versionPosition + 1);
  std::string commitNumber = version.substr(versionPosition,afterVersionPosition - versionPosition);
  std::string z3Version = Z3Interfacing::z3_full_version();
  size_t spacePosition = z3Version.find(" ");
  if (spacePosition != string::npos) {
    z3Version = z3Version.substr(0,spacePosition);
  }

  std::string problemName = Lib::env.options->problemName();
  std::cout
    << status << " "
    << (problemName.length() == 0 ? "unknown" : problemName) << " "
    << Timer::elapsedDeciseconds() << " "
    << Timer::elapsedMegaInstructions() << " "
    << (Lib::env.options ? Lib::env.options->testId() : "unknown") << " "
    << commitNumber << ':' << z3Version << endl;
#endif
}

bool szsOutputMode() {
  return (Lib::env.options && Lib::env.options->outputMode() == Shell::Options::Output::SZS);
}

ostream& addCommentSignForSZS(ostream& out)
{
  if (szsOutputMode()) {
    out << "% ";
    if (Lib::env.options && Lib::env.options->multicore() != 1) {
      out << "(" << getpid() << ")";
    }
  }
  return out;
}

Stack<UIHelper::LoadedPiece> UIHelper::_loadedPieces = { UIHelper::LoadedPiece() }; // start initialized with a singleton

bool UIHelper::s_expecting_sat=false;
bool UIHelper::s_expecting_unsat=false;

bool UIHelper::portfolioParent=false;
bool UIHelper::satisfiableStatusWasAlreadyOutput=false;

bool UIHelper::spiderOutputDone = false;

void UIHelper::outputAllPremises(ostream& out, UnitList* units, std::string prefix)
{
#if 1
  InferenceStore::instance()->outputProof(cerr, units);
#else
  Stack<UnitSpec> prems;
  Stack<UnitSpec> toDo;
  DHSet<UnitSpec> seen;

  //get the units to start with
  UnitList::Iterator uit(units);
  while (uit.hasNext()) {
    Unit* u = uit.next();
    toDo.push(UnitSpec(u));
  }

  while (toDo.isNonEmpty()) {
    UnitSpec us = toDo.pop();
    UnitSpecIterator pars = InferenceStore::instance()->getParents(us);
    while (pars.hasNext()) {
      UnitSpec par = pars.next();
      if (seen.contains(par)) {
	continue;
      }
      prems.push(par);
      toDo.push(par);
      seen.insert(par);
    }
  }

  std::sort(prems.begin(), prems.end(), UIHelper::unitSpecNumberComparator);

  Stack<UnitSpec>::BottomFirstIterator premIt(prems);
  while (premIt.hasNext()) {
    UnitSpec prem = premIt.next();
    out << prefix << prem.toString() << endl;
  }
#endif
}

void UIHelper::outputSaturatedSet(ostream& out, UnitIterator uit)
{
  addCommentSignForSZS(out);
  out << "# SZS output start Saturation." << endl;

  while (uit.hasNext()) {
    Unit* cl = uit.next();
    out << TPTPPrinter::toString(cl) << endl;
  }

  addCommentSignForSZS(out);
  out << "# SZS output end Saturation." << endl;
} // outputSaturatedSet

// String utility function that probably belongs elsewhere
static bool hasEnding (std::string const &fullString, std::string const &ending) {
  if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
  } else {
      return false;
  }
}

void UIHelper::tryParseTPTP(istream& input)
{
  LoadedPiece& curPiece = _loadedPieces.top();
  Parse::TPTP parser(input,curPiece._units);
  try {
    parser.parse();
    curPiece._units = parser.unitBuffer();
    curPiece._hasConjecture |= parser.containsConjecture();
  } catch (ParsingRelatedException& exception) {
    UnitList::destroy(curPiece._units.clipAtLast()); // destroy units that perhaps got already parsed
    throw;
  }
}

void UIHelper::tryParseSMTLIB2(istream& input)
{
  LoadedPiece& curPiece = _loadedPieces.top();
  Parse::SMTLIB2 parser(curPiece._units);
  try {
    parser.parse(input);
    Unit::onParsingEnd(); // dubious in interactiveMetamode (influences SMT goal guessing and InferenceStore::ProofPropertyPrinter)
    curPiece._units = parser.formulaBuffer();
    curPiece._smtLibLogic = parser.getLogic();
  } catch (ParsingRelatedException& exception) {
    UnitList::destroy(curPiece._units.clipAtLast()); // destroy units that perhaps got already parsed
    throw;
  }

#if VDEBUG
  const std::string& expected_status = parser.getStatus();
  if (expected_status == "sat") {
    s_expecting_sat = true;
  } else if (expected_status == "unsat") {
    s_expecting_unsat = true;
  }
#endif
}

void UIHelper::parseSingleLine(const std::string& lineToParse, Options::InputSyntax inputSyntax)
{
  LoadedPiece newPiece = _loadedPieces.top();  // copy everything
  newPiece._id = lineToParse;
  _loadedPieces.push(std::move(newPiece));

  ScopedLet<Statistics::ExecutionPhase> localAssing(env.statistics->phase,Statistics::PARSING);

  std::istringstream stream(lineToParse);
  try {
    switch (inputSyntax) {
      case Options::InputSyntax::TPTP:
        tryParseTPTP(stream);
        break;
      case Options::InputSyntax::SMTLIB2:
        tryParseSMTLIB2(stream);
        break;
      case Options::InputSyntax::AUTO:
        ASSERTION_VIOLATION;
        break;
    }
  } catch (ParsingRelatedException& exception) {
    _loadedPieces.pop();
    throw;
  }
}

// Call this function to report a parsing attempt has failed and to reset the input
void resetParsing(ParsingRelatedException& exception, istream& input, std::string nowtry)
{
  if (env.options->mode()!=Options::Mode::SPIDER) {
    addCommentSignForSZS(std::cout);
    std::cout << "Failed with\n";
    addCommentSignForSZS(std::cout);
    exception.cry(std::cout);
    addCommentSignForSZS(std::cout);
    std::cout << "Trying " << nowtry  << endl;
  }

  input.clear();
  input.seekg(0);
}

void UIHelper::parseStream(std::istream& input, Options::InputSyntax inputSyntax, bool verbose, bool preferSMTonAuto)
{
  switch (inputSyntax) {
  case Options::InputSyntax::AUTO:
    if (preferSMTonAuto){
      if (verbose) {
        addCommentSignForSZS(std::cout);
        std::cout << "Running in auto input_syntax mode. Trying SMTLIB2\n";
      }
      try {
        tryParseSMTLIB2(input);
      } catch (ParsingRelatedException& exception) {
        resetParsing(exception,input,"TPTP");
        tryParseTPTP(input);
      }
    } else {
      if (verbose) {
        addCommentSignForSZS(std::cout);
        std::cout << "Running in auto input_syntax mode. Trying TPTP\n";
      }
      try {
        tryParseTPTP(input);
      } catch (ParsingRelatedException& exception) {
        resetParsing(exception,input,"SMTLIB2");
        tryParseSMTLIB2(input);
      }
    }
    break;
  case Options::InputSyntax::TPTP:
    tryParseTPTP(input);
    break;
  case Options::InputSyntax::SMTLIB2:
    tryParseSMTLIB2(input);
    break;
  }
}

void UIHelper::parseStandardInput(Options::InputSyntax inputSyntax)
{
  LoadedPiece newPiece = _loadedPieces.top();  // copy everything
  newPiece._id = "<cin>";
  _loadedPieces.push(std::move(newPiece));

  if (inputSyntax == Options::InputSyntax::AUTO) {
    addCommentSignForSZS(std::cout);
    std::cout << "input_syntax=auto not supported for standard input parsing, switching to tptp.\n";

    inputSyntax = Options::InputSyntax::TPTP;
  }
  try {
    parseStream(cin,inputSyntax,false,false);
  } catch (ParsingRelatedException& exception) {
    _loadedPieces.pop();
    throw;
  }
}

void UIHelper::parseFile(const std::string& inputFile, Options::InputSyntax inputSyntax, bool verbose)
{
  LoadedPiece newPiece = _loadedPieces.top();  // copy everything
  newPiece._id = inputFile;
  _loadedPieces.push(std::move(newPiece));

  TIME_TRACE(TimeTrace::PARSING);
  ScopedLet<Statistics::ExecutionPhase> localAssing(env.statistics->phase,Statistics::PARSING);

  ifstream input(inputFile.c_str());
  if (input.fail()) {
    USER_ERROR("Cannot open problem file: "+inputFile);
  }

  try {
    parseStream(input,inputSyntax,verbose,hasEnding(inputFile,"smt") || hasEnding(inputFile,"smt2"));
  } catch (ParsingRelatedException& exception) {
    _loadedPieces.pop();
    throw;
  }
}

/**
 * After a single call (or a series of calls) to parse* functions,
 * return a problem object with the obtained units.
 *
 * No preprocessing is performed on the units.
 *
 * The Options object should intentionally not be part of this game,
 * as any form of "conditional parsing" compromises the effective use of the correspoding conditioning options
 * as a part of strategy development and use in portfolios. In other words, if you need getInputProblem or the parse* functions
 * to depend on an option, think twice, and if really needed, make it an explicit argument of that function.
 */
Problem* UIHelper::getInputProblem()
{
  LoadedPiece& topPiece = _loadedPieces.top();
  Problem* res = new Problem(topPiece._units.list());

  // NB this must happen immediately, as the Property relies on it
  res->setSMTLIBLogic(topPiece._smtLibLogic);

  if(res->isHigherOrder())
    USER_ERROR(
      "This version of Vampire is not yet HOLy.\n\n"
      "Support for higher-order logic is currently on the ahmed-new-hol branch.\n"
      "HOL should be coming to mainline 'soon'."
    );

  env.setMainProblem(res);
  return res;
}

void UIHelper::listLoadedPieces(ostream& out)
{
  auto it = _loadedPieces.iterFifo();
  ALWAYS(it.next()._id.empty()); // skip the first, empty, entry
  while (it.hasNext()) {
    out << " " << it.next()._id << endl;
  }
}

void UIHelper::popLoadedPiece(int numPops)
{
  while (numPops-- > 0) {
    if (_loadedPieces.size() > 1) {
      _loadedPieces.pop();
      UnitList::destroy(_loadedPieces.top()._units.clipAtLast());
    }
  }
}

/**
 * Output result based on the content of
 * @b env.statistics->terminationReason
 *
 * If LaTeX output is enabled, it is output in this function.
 *
 * If interpolant output is enabled, it is output in this function.
 */
void UIHelper::outputResult(ostream& out)
{
  switch (env.statistics->terminationReason) {
  case Statistics::REFUTATION: {
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unsat" << endl;
      return;
    }
    if(env.options->outputMode() == Options::Output::UCORE){
      out << "unsat" << endl;
      InferenceStore::instance()->outputUnsatCore(out, env.statistics->refutation);
      return;
    }
    addCommentSignForSZS(out);
    out << "Refutation found. Thanks to " << env.options->thanks() << "!"
      << endl; // let's have this flushed, as we might now take a bit of time to SAT-minimize the proof

    Unit* refutation = env.statistics->refutation;

    /**
     * Making explicit what was previously only done during proof printing.
     *
     * However, we need to know the correct input type to decide whether to print "Theorem" vs "ContradictoryAxioms" below
     * (and this is not correctly set in inference by Global subsumtions for heuristic reasons).
     *
     * So it makes sense to call this even if minimizeSatProofs is false.
     *
     * Also induction statistics deserve to be correct even if we don't print a proof.
     */
    bool seenInputInference = refutation->minimizeAncestorsAndUpdateSelectedStats();
    // minimization might have cause inductionDepth to change (in fact, decrease)
    env.statistics->maxInductionDepth = refutation->inference().inductionDepth();

    if (szsOutputMode()) {
      out << "% SZS status " <<
        (UIHelper::haveConjecture() ? ( refutation->derivedFromGoal() ? "Theorem" : "ContradictoryAxioms" ) : "Unsatisfiable")
	      << " for " << env.options->problemName() << endl;
    }
    if (env.options->questionAnswering()!=Options::QuestionAnsweringMode::OFF) {
      ASS(refutation->isClause());
      AnswerLiteralManager::getInstance()->tryOutputAnswer(static_cast<Clause*>(env.statistics->refutation),std::cout);
    }
    if (env.options->proof() != Options::Proof::OFF) {
      if (szsOutputMode()) {
        out << "% SZS output start Proof for " << env.options->problemName() << endl;
      }
      InferenceStore::instance()->outputProof(out, refutation);
      if (szsOutputMode()) {
        out << "% SZS output end Proof for " << env.options->problemName() << endl << flush;
      }

    }
    if (env.options->showInterpolant()!=Options::InterpolantMode::OFF) {
      ASS(refutation->isClause());

      Interpolants::removeConjectureNodesFromRefutation(refutation);
      Unit* formulifiedRefutation = Interpolants::formulifyRefutation(refutation);

      Formula* interpolant = nullptr;

      switch(env.options->showInterpolant()) {
      // new interpolation methods described in master thesis of Bernhard Gleiss
      case Options::InterpolantMode::NEW_HEUR:
        Interpolants().removeTheoryInferences(formulifiedRefutation); // do this only once for each proof!
        interpolant = Interpolants().getInterpolant(formulifiedRefutation, Interpolants::UnitWeight::VAMPIRE);
        break;
#if VZ3
      case Options::InterpolantMode::NEW_OPT:

        Interpolants().removeTheoryInferences(formulifiedRefutation); // do this only once for each proof!
        interpolant = InterpolantMinimizer().getInterpolant(formulifiedRefutation, Interpolants::UnitWeight::VAMPIRE);
        break;
#endif
      case Options::InterpolantMode::OFF:
        ASSERTION_VIOLATION;
      }

      out << "Symbol-weight minimized interpolant: " << TPTPPrinter::toString(interpolant) << endl;
      out << "Actual weight: " << interpolant->weight() << endl;
      out<<endl;
    }

    if (env.options->latexOutput() != "off") {
      ofstream latexOut(env.options->latexOutput().c_str());

      LaTeX formatter;
      latexOut << formatter.refutationToString(refutation);
    }

    // the following two sanity checks are performed only after the proof printing, so we can also have a look at the proof, when we get a report back

    if(refutation->isPureTheoryDescendant()) {
      INVALID_OPERATION("A pure theory descendant is empty, which means theory axioms are inconsistent.");
      break;
    }

    if(!seenInputInference){
      INVALID_OPERATION("The proof does not contain any input formulas.");
      break;
    }

    ASS(!s_expecting_sat);
    break;
  }
  case Statistics::TIME_LIMIT:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unknown" << endl;
      return;
    }
    addCommentSignForSZS(out);
    out << "Time limit reached!\n";
    break;
  case Statistics::MEMORY_LIMIT:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unknown" << endl;
      return;
    }
    addCommentSignForSZS(out);
    out << "Memory limit exceeded!\n";
    break;
  case Statistics::ACTIVATION_LIMIT: {
    addCommentSignForSZS(out);
    out << "Activation limit reached!\n";

    // HERE ADD MORE

    break;
  }
  case Statistics::REFUTATION_NOT_FOUND:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unknown" << endl;
      return;
    }
    addCommentSignForSZS(out);
    env.statistics->explainRefutationNotFound(out);
    break;
  case Statistics::SATISFIABLE:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "sat" << endl;
      return;
    }
    outputSatisfiableResult(out);

    ASS(!s_expecting_unsat);

    break;
  case Statistics::INAPPROPRIATE:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unknown" << endl;
      return;
    }
    addCommentSignForSZS(out);
    out << "Terminated due to inappropriate strategy.\n";
    break;
  case Statistics::UNKNOWN:
    if(env.options->outputMode() == Options::Output::SMTCOMP){
      out << "unknown" << endl;
      return;
    }
    addCommentSignForSZS(out);
    out << "Unknown reason of termination!\n";
    break;
  default:
    ASSERTION_VIOLATION;
  }
  env.statistics->print(out);
}

void UIHelper::outputSatisfiableResult(ostream& out)
{
  //out << "Satisfiable!\n";
  if (szsOutputMode() && !satisfiableStatusWasAlreadyOutput) {
    out << "% SZS status " << ( UIHelper::haveConjecture() ? "CounterSatisfiable" : "Satisfiable" )
	  <<" for " << env.options->problemName() << endl;
  }
  if (!env.statistics->model.empty()) {
    if (szsOutputMode()) {
	out << "% SZS output start FiniteModel for " << env.options->problemName() << endl;
    }
    out << env.statistics->model;
    if (szsOutputMode()) {
	out << "% SZS output end FiniteModel for " << env.options->problemName() << endl;
    }
  }
  else //if (env.statistics->saturatedSet)
       /* -- MS: it's never incorrect to print the empty one, in fact this prevents us from losing
        * points at CASC when the input gets completely emptied, by e.g. preprocessing
        */
  {
    outputSaturatedSet(out, pvi(UnitList::Iterator(env.statistics->saturatedSet)));
  }
}

/**
 * Output to @b out all symbol declarations for the current signature.
 * Symbols having default types will not be output.
 * @author Andrei Voronkov
 * @since 03/07/2013 Manchester
 */
void UIHelper::outputSymbolDeclarations(ostream& out)
{
  Signature& sig = *env.signature;

  unsigned typeCons = sig.typeCons();
  for (unsigned i=0; i<typeCons; ++i) {
    outputSymbolTypeDeclarationIfNeeded(out, false, true, i);
  }
  unsigned funcs = sig.functions();
  for (unsigned i=0; i<funcs; ++i) {
    if (!env.options->showFOOL()) {
      if (env.signature->isFoolConstantSymbol(true,i) || env.signature->isFoolConstantSymbol(false,i)) {
        continue;
      }
    }
    outputSymbolTypeDeclarationIfNeeded(out, true, false, i);
  }
  unsigned preds = sig.predicates();
  for (unsigned i=0; i<preds; ++i) {
    outputSymbolTypeDeclarationIfNeeded(out, false, false, i);
  }
} // UIHelper::outputSymbolDeclarations

/**
 * Output to @b out a function or a predicate symbol declaration.
 * Symbols having default types will not be output.
 * @author Andrei Voronkov
 * @since 03/07/2013 Manchester
 */
void UIHelper::outputSymbolTypeDeclarationIfNeeded(ostream& out, bool function, bool typeCon, unsigned symNumber)
{
  Signature::Symbol* sym;

  if(function){
    sym = env.signature->getFunction(symNumber);
  } else if(typeCon){
    sym = env.signature->getTypeCon(symNumber);
  } else {
    sym = env.signature->getPredicate(symNumber);
  }

  if (typeCon && (env.signature->isArrayCon(symNumber) ||
                  env.signature->isTupleCon(symNumber))){
    return;
  }

  if(typeCon && env.signature->isDefaultSortCon(symNumber) &&
    (!env.signature->isBoolCon(symNumber) || !env.options->showFOOL())){
    return;
  }

  if (sym->interpreted()) {
    //there is no need to output type definitions for interpreted symbols
    return;
  }

  if (sym->overflownConstant()) {
    // don't output definitions of numbers; not even big ones
    return;
  }

  unsigned dummy;
  if (!typeCon && Theory::tuples()->findProjection(symNumber, !function, dummy)) {
    return;
  }

  if (function) {
    TermList sort = env.signature->getFunction(symNumber)->fnType()->result();
    if (sort.isTupleSort()) {
      return;
    }
  }

  OperatorType* type = function ? sym->fnType() :
               (typeCon ? sym->typeConType() : sym->predType());

  if (type->isAllDefault() && env.options->proof() != Options::Proof::DEDUKTI) {//TODO required
    return;
  }

  //out << "tff(" << (function ? "func" : "pred") << "_def_" << symNumber << ", type, "
  //    << sym->name() << ": ";

  std::string symName = sym->name();
  if(typeCon && env.signature->isBoolCon(symNumber)){
    ASS(env.options->showFOOL());
    symName = "$bool";
  }

  //don't output type of app. It is an internal Vampire thing
  if(!(function && env.signature->isAppFun(symNumber))){
    if(env.options->proof() == Options::Proof::DEDUKTI) {
      Dedukti::outputTypeDecl(env.out(), symName.c_str(), type);
      return;
    }

    out << (env.getMainProblem()->isHigherOrder() ? "thf(" : "tff(")
        << (function ? "func" : (typeCon ?  "type" : "pred"))
        << "_def_" << symNumber << ", type, "
        << symName << ": ";
    out << type->toString();
    out << ")." << endl;
  }
  //out << ")." << endl;
}

} // namespace Shell
