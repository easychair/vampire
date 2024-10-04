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
 * @file vampire.cpp. Implements the top-level procedures of Vampire.
 */

#include <iostream>
#include <ostream>
#include <fstream>

#if VZ3
#include "z3++.h"
#endif

#include "Debug/TimeProfiling.hpp"
#include "Lib/Exception.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Random.hpp"
#include "Lib/Timer.hpp"
#include "Lib/List.hpp"
#include "Lib/System.hpp"
#include "Lib/StringUtils.hpp"
#include "Lib/Sys/Multiprocessing.hpp"
#include "Lib/Int.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/Formula.hpp"
#include "Kernel/FormulaUnit.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/Term.hpp"

#include "Inferences/InferenceEngine.hpp"
#include "Inferences/TautologyDeletionISE.hpp"

#include "CASC/PortfolioMode.hpp"
#include "Shell/CommandLine.hpp"

#include "Shell/Dedukti.hpp"
#include "Shell/Normalisation.hpp"
#include "Shell/Options.hpp"
#include "Shell/Property.hpp"
#include "Saturation/ProvingHelper.hpp"
#include "Shell/Preprocess.hpp"
#include "Shell/TheoryFinder.hpp"
#include "Shell/TPTPPrinter.hpp"
#include "Parse/TPTP.hpp"
#include "Shell/FOOLElimination.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/UIHelper.hpp"
#include "Shell/LaTeX.hpp"
#include "Shell/SineUtils.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

#include "FMB/ModelCheck.hpp"

using namespace std;

/**
 * Return value is non-zero unless we were successful.
 *
 * Being successful for modes that involve proving means that we have
 * either found refutation or established satisfiability.
 *
 *
 * If Vampire was interrupted (e.g. SIGINT, SIGHUP), value VAMP_RESULT_STATUS_INTERRUPTED is returned,
 * and in case of other signal we return VAMP_RESULT_STATUS_OTHER_SIGNAL. For implementation
 * of these return values see Lib/System.hpp.
 *
 * In case of an unhandled exception or user error, we return value
 * VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION.
 *
 * In case Vampire was terminated by the timer, return value is
 * uncertain (but definitely not zero), probably it will be 134
 * (we terminate by a call to the @b abort() function in this case).
 */
int vampireReturnValue = VAMP_RESULT_STATUS_UNKNOWN;

/**
 * Preprocess the given problem (in dependence of env.options).
 *
 * The problem is modified destructively.
 */
[[nodiscard]]
Problem* preprocessProblem(Problem* prb)
{
  // Here officially starts preprocessing of vampireMode
  // and that's the moment we want to set the random seed (no randomness in parsing, for the peace of mind)
  // the main reason being that we want to stay in sync with what profolio mode will do
  // cf ProvingHelper::runVampire
  Lib::Random::setSeed(env.options->randomSeed());

  TIME_TRACE(TimeTrace::PREPROCESSING);

  // this will provide warning if options don't make sense for problem
  if (env.options->mode()!=Options::Mode::SPIDER) {
    env.options->checkProblemOptionConstraints(prb->getProperty(), /*before_preprocessing = */ true);
  }

  Shell::Preprocess prepro(*env.options);
  //phases for preprocessing are being set inside the preprocess method
  prepro.preprocess(*prb);

  return prb;
} // getPreprocessedProblem

void explainException(Exception& exception)
{
  exception.cry(std::cout);
} // explainException

[[nodiscard]]
Problem *doProving(Problem* problem)
{
  // a new strategy randomization mechanism
  if (!env.options->strategySamplerFilename().empty()) {
    env.options->sampleStrategy(env.options->strategySamplerFilename());
  }

  env.options->setForcedOptionValues();
  env.options->checkGlobalOptionConstraints();

  Problem *prb = preprocessProblem(problem);

  // this will provide warning if options don't make sense for problem
  if (env.options->mode()!=Options::Mode::SPIDER) {
    env.options->checkProblemOptionConstraints(prb->getProperty(), /*before_preprocessing = */ false);
  }

  ProvingHelper::runVampireSaturation(*prb, *env.options);
  return prb;
}

/**
 * Read a problem and output profiling information about it.
 * @since 03/08/2008 Torrevieja
 */
void profileMode(Problem* problem)
{
  ScopedPtr<Problem> prb(problem);

  /* CAREFUL: Make sure that the order
   * 1) getProperty, 2) normalise, 3) TheoryFinder::search
   * is the same as in PortfolioMode::searchForProof
   * also, cf. the beginning of Preprocessing::preprocess*/
  Property* property = prb->getProperty();
  Normalisation().normalise(*prb);
  TheoryFinder(prb->units(), property).search();

  std::cout << property->categoryString() << ' ' << property->props() << ' '
	  << property->atoms() << "\n";

  //we have succeeded with the profile mode, so we'll terminate with zero return value
  vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
} // profileMode

// prints Unit u at an index to latexOut using the LaTeX object
void outputUnitToLaTeX(LaTeX& latex, ofstream& latexOut, Unit* u,unsigned index)
{
    std::string stringform = latex.toString(u);
    latexOut << index++ << " & ";
    unsigned count = 0;
    for(const char* p = stringform.c_str();*p;p++){
      latexOut << *p;
      count++;
      if(count>80 && *p==' '){
        latexOut << "\\\\ \n & ~~~~~";
        count=0;
      }
    }
    latexOut << "\\\\" << endl;
}

// print the clauses of a problem to a LaTeX file
void outputClausesToLaTeX(Problem* prb)
{
  ASS(env.options->latexOutput()!="off");

  LaTeX latex;
  ofstream latexOut(env.options->latexOutput().c_str());
  latexOut << latex.header() << endl;
  latexOut << "\\section{Problem "<<env.options->problemName() << "}" << endl;
  //TODO output more header
  latexOut << "\\[\n\\begin{array}{ll}" << endl;

  CompositeISE simplifier;
  simplifier.addFront(new TrivialInequalitiesRemovalISE());
  simplifier.addFront(new TautologyDeletionISE());
  simplifier.addFront(new DuplicateLiteralRemovalISE());

  unsigned index=0;
  ClauseIterator cit = prb->clauseIterator();
  while (cit.hasNext()) {
    Clause* cl = cit.next();
    cl = simplifier.simplify(cl);
    if (!cl) {
      continue;
    }
    outputUnitToLaTeX(latex,latexOut,cl,index++);
  }
  latexOut  << "\\end{array}\n\\]" << latex.footer() << "\n";

  //close ofstream?
}

// print the formulas of a problem to a LaTeX file
void outputProblemToLaTeX(Problem* prb)
{
  ASS(env.options->latexOutput()!="off");

  LaTeX latex;
  ofstream latexOut(env.options->latexOutput().c_str());
  latexOut << latex.header() << endl;
  latexOut << "\\section{Problem "<<env.options->problemName() << "}" << endl;
  //TODO output more header
  latexOut << "\\[\n\\begin{array}{ll}" << endl;

  //TODO  get symbol and sort declarations into LaTeX

  UnitList::Iterator units(prb->units());

  unsigned index = 0;
  while (units.hasNext()) {
    Unit* u = units.next();
    outputUnitToLaTeX(latex,latexOut,u,index++);
  }
  latexOut  << "\\end{array}\n\\]" << latex.footer() << "\n";

  //close ofstream?
}

/**
 * This mode only preprocesses the input using the current preprocessing
 * options and outputs it to stdout. It is useful for either preprocessing
 * per se or also for converting one syntax to another. For the latter, the input
 * and the output syntaxes must be set to different values. Note that for
 * simply translating one syntax to another, output mode is the right one.
 *
 * @author Andrei Voronkov
 * @since 02/07/2013 Manchester
 */
void preprocessMode(Problem* problem, bool theory)
{
  ScopedPtr<Problem> prb(problem);

  TIME_TRACE(TimeTrace::PREPROCESSING);

  // preprocess without clausification
  Shell::Preprocess prepro(*env.options);
  prepro.turnClausifierOff();
  if(env.options->mode() == Options::Mode::PREPROCESS2){
    prepro.keepSimplifyStep();
  }
  prepro.preprocess(*prb);

  //outputSymbolDeclarations also deals with sorts for now
  //UIHelper::outputSortDeclarations(std::cout);
  UIHelper::outputSymbolDeclarations(std::cout);
  UnitList::Iterator units(prb->units());
  while (units.hasNext()) {
    Unit* u = units.next();
    if (!env.options->showFOOL()) {
      if (u->inference().rule() == InferenceRule::FOOL_AXIOM_TRUE_NEQ_FALSE || u->inference().rule() == InferenceRule::FOOL_AXIOM_ALL_IS_TRUE_OR_FALSE) {
        continue;
      }
    }

    if (theory) {
      Formula* f = u->getFormula();

      // CONJECTURE as inputType is evil, as it cannot occur multiple times
      if (u->inference().inputType() == UnitInputType::CONJECTURE) {
        u->inference().setInputType(UnitInputType::NEGATED_CONJECTURE);
      }

      FormulaUnit* fu = new FormulaUnit(f,u->inference()); // we are stealing u's inference which is not nice
      std::cout << TPTPPrinter::toString(fu) << "\n";
    } else {
      std::cout << TPTPPrinter::toString(u) << "\n";
    }
  }

  if(env.options->latexOutput()!="off"){ outputProblemToLaTeX(prb.ptr()); }

  //we have successfully output all clauses, so we'll terminate with zero return value
  vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
} // preprocessMode

/**
 *
 * @author Giles
 * @since 6/10/2015
 */
void modelCheckMode(Problem* problem)
{
  ScopedPtr<Problem> prb(problem);
  
  if(env.getMainProblem()->hasPolymorphicSym() || env.getMainProblem()->isHigherOrder()){
    USER_ERROR("Polymorphic Vampire is not yet compatible with theory reasoning");
  }

  FMB::ModelCheck::doCheck(prb->units());
} // modelCheckMode


/**
 * This mode only outputs the input problem. It is useful for converting
 * one syntax to another.
 * @author Laura Kovacs and Andrei Voronkov
 * @since 02/07/2013 Gothenburg and Manchester
 */
void outputMode(Problem* problem)
{
  ScopedPtr<Problem> prb(problem);

  bool dedukti = env.options->proof() == Options::Proof::DEDUKTI;
  if(dedukti)
    Dedukti::outputPrelude(std::cout);

  //outputSymbolDeclarations also deals with sorts for now
  //UIHelper::outputSortDeclarations(std::cout);
  UIHelper::outputSymbolDeclarations(std::cout);
  UnitList::Iterator units(prb->units());

  while (units.hasNext()) {
    Unit* u = units.next();
    if(dedukti)
      Dedukti::outputAxiom(std::cout, u);
    else
      std::cout << TPTPPrinter::toString(u) << "\n";
  }

  if(env.options->latexOutput()!="off"){ outputProblemToLaTeX(prb.ptr()); }

  //we have successfully output all clauses, so we'll terminate with zero return value
  vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
} // outputMode

void vampireMode(Problem* problem)
{
  if (env.options->mode() == Options::Mode::CONSEQUENCE_ELIMINATION) {
    env.options->setUnusedPredicateDefinitionRemoval(false);
  }

  ScopedPtr<Problem> prb(doProving(problem));

  UIHelper::outputResult(std::cout);

  if (env.statistics->terminationReason == Statistics::REFUTATION
      || env.statistics->terminationReason == Statistics::SATISFIABLE) {
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
  }
} // vampireMode

void spiderMode(Problem* problem)
{
  env.options->setBadOptionChoice(Options::BadOption::HARD);
  env.options->setOutputMode(Options::Output::SPIDER);
  env.options->setNormalize(true);

  Exception* exception = 0;
#if VZ3
  z3::exception* z3_exception = 0;
#endif

  bool exceptionRaised = false;
  ScopedPtr<Problem> prb;
  try {
    prb = doProving(problem);
  } catch (Exception& e) {
    exception = &e;
    exceptionRaised = true;
#if VZ3
  } catch(z3::exception& e){
    z3_exception = &e;
    exceptionRaised = true;
#endif
  } catch (...) {
    exceptionRaised = true;
  }

  if (!exceptionRaised) {
    switch (env.statistics->terminationReason) {
    case Statistics::REFUTATION:
      reportSpiderStatus('+');
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
      break;
    case Statistics::TIME_LIMIT:
      reportSpiderStatus('t');
      break;
    case Statistics::MEMORY_LIMIT:
      reportSpiderStatus('m');
      break;
    case Statistics::UNKNOWN:
    case Statistics::INAPPROPRIATE:
      reportSpiderStatus('u');
      break;
    case Statistics::REFUTATION_NOT_FOUND:
      if (env.statistics->discardedNonRedundantClauses > 0) {
        reportSpiderStatus('n');
      }
      else{
        reportSpiderStatus('i');
      }
      break;
    case Statistics::SATISFIABLE:
      reportSpiderStatus('-');
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
      break;
    default:
      ASSERTION_VIOLATION;
    }
    return;
  }

  // exception
  vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;

#if VZ3
  if (z3_exception) {
    if (strcmp(z3_exception->msg(),"out of memory\n")) {
      reportSpiderStatus('m');
    }
    else {
      reportSpiderFail();
    }

    return;
  }
#endif

  reportSpiderFail();

  ASS(exception);
  explainException(*exception);
} // spiderMode

void clausifyMode(Problem* problem, bool theory)
{
  CompositeISE simplifier;
  simplifier.addFront(new TrivialInequalitiesRemovalISE());
  simplifier.addFront(new TautologyDeletionISE());
  simplifier.addFront(new DuplicateLiteralRemovalISE());

  if (!env.options->strategySamplerFilename().empty()) {
    env.options->sampleStrategy(env.options->strategySamplerFilename());
  }

  ScopedPtr<Problem> prb(preprocessProblem(problem));

  //outputSymbolDeclarations deals with sorts as well for now
  //UIHelper::outputSortDeclarations(std::cout);
  UIHelper::outputSymbolDeclarations(std::cout);

  ClauseIterator cit = prb->clauseIterator();
  bool printed_conjecture = false;
  while (cit.hasNext()) {
    Clause* cl = cit.next();
    cl = simplifier.simplify(cl);
    if (!cl) {
      continue;
    }
    printed_conjecture |= cl->inputType() == UnitInputType::CONJECTURE || cl->inputType() == UnitInputType::NEGATED_CONJECTURE;
    if (theory) {
      Formula* f = Formula::fromClause(cl);

      // CONJECTURE as inputType is evil, as it cannot occur multiple times
      if (cl->inference().inputType() == UnitInputType::CONJECTURE) {
        cl->inference().setInputType(UnitInputType::NEGATED_CONJECTURE);
      }

      FormulaUnit* fu = new FormulaUnit(f,cl->inference()); // we are stealing cl's inference, which is not nice!
      fu->overwriteNumber(cl->number()); // we are also making sure it's number is the same as that of the original (for Kostya from Russia to CASC, with love, and back again)
      std::cout << TPTPPrinter::toString(fu) << "\n";
    } else {
      std::cout << TPTPPrinter::toString(cl) << "\n";
    }
  }
  if(!printed_conjecture && UIHelper::haveConjecture()){
    unsigned p = env.signature->addFreshPredicate(0,"p");
    auto c = Clause::fromLiterals({
        Literal::create(p, /* polarity */ true , {}),
        Literal::create(p, /* polarity */ false, {})
      }, 
      NonspecificInference0(UnitInputType::NEGATED_CONJECTURE,InferenceRule::INPUT));
    std::cout << TPTPPrinter::toString(c) << "\n";
  }

  if (env.options->latexOutput() != "off") { outputClausesToLaTeX(prb.ptr()); }

  //we have successfully output all clauses, so we'll terminate with zero return value
  vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
} // clausifyMode

void axiomSelectionMode(Problem* problem)
{
  ScopedPtr<Problem> prb(problem);

  env.options->setSineSelection(Options::SineSelection::AXIOMS);

  if (prb->hasFOOL()) {
    FOOLElimination().apply(*prb);
  }

  // reorder units
  if (env.options->normalize()) {
    env.statistics->phase = Statistics::NORMALIZATION;
    Normalisation norm;
    norm.normalise(*prb);
  }

  env.statistics->phase = Statistics::SINE_SELECTION;
  Shell::SineSelector(*env.options).perform(*prb);

  env.statistics->phase = Statistics::FINALIZATION;

  UnitList::Iterator uit(prb->units());
  while (uit.hasNext()) {
    Unit* u = uit.next();
    std::cout << TPTPPrinter::toString(u) << "\n";
  }

  //we have successfully output the selected units, so we'll terminate with zero return value
  vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
}

void dispatchByMode(Problem* problem)
{
  switch (env.options->mode())
  {
  case Options::Mode::AXIOM_SELECTION:
    axiomSelectionMode(problem);
    break;
  case Options::Mode::SPIDER:
    spiderMode(problem);
    break;
  case Options::Mode::CONSEQUENCE_ELIMINATION:
  case Options::Mode::VAMPIRE:
    vampireMode(problem);
    break;

  case Options::Mode::CASC:
    env.options->setIgnoreMissing(Options::IgnoreMissing::WARN);
    if (env.options->intent() == Options::Intent::UNSAT) {
      env.options->setSchedule(Options::Schedule::CASC);
    } else {
      env.options->setSchedule(Options::Schedule::CASC_SAT);
    }
    env.options->setInputSyntax(Options::InputSyntax::TPTP);
    env.options->setOutputMode(Options::Output::SZS);
    env.options->setProof(Options::Proof::TPTP);
    env.options->setOutputAxiomNames(true);
    env.options->setNormalize(true);
    env.options->setRandomizeSeedForPortfolioWorkers(false);

    if (CASC::PortfolioMode::perform(problem)) {
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
    }
    break;

  case Options::Mode::CASC_HOL: {
    env.options->setIgnoreMissing(Options::IgnoreMissing::WARN);
    env.options->setSchedule(Options::Schedule::CASC_HOL_2020);
    env.options->setInputSyntax(Options::InputSyntax::TPTP);
    env.options->setOutputMode(Options::Output::SZS);
    env.options->setProof(Options::Proof::TPTP);
    env.options->setOutputAxiomNames(true);

    if (CASC::PortfolioMode::perform(problem)) {
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
    }
    break;
  }

  case Options::Mode::SMTCOMP:
    env.options->setIgnoreMissing(Options::IgnoreMissing::OFF);
    env.options->setInputSyntax(Options::InputSyntax::SMTLIB2);
    if(env.options->outputMode() != Options::Output::UCORE){
      env.options->setOutputMode(Options::Output::SMTCOMP);
    }
    env.options->setSchedule(Options::Schedule::SMTCOMP);
    env.options->setProof(Options::Proof::OFF);
    env.options->setNormalize(true);
    env.options->setRandomizeSeedForPortfolioWorkers(false);

    env.options->setMulticore(0); // use all available cores
    env.options->setTimeLimitInSeconds(1800);
    env.options->setStatistics(Options::Statistics::NONE);

    //TODO needed?
    // to prevent from terminating by time limit
    env.options->setTimeLimitInSeconds(100000);

    if (CASC::PortfolioMode::perform(problem)){
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
    }
    else {
      cout << "unknown" << endl;
    }
    break;

  case Options::Mode::PORTFOLIO:
    env.options->setIgnoreMissing(Options::IgnoreMissing::WARN);

    if (CASC::PortfolioMode::perform(problem)) {
      vampireReturnValue = VAMP_RESULT_STATUS_SUCCESS;
    }
    break;
  case Options::Mode::MODEL_CHECK:
    modelCheckMode(problem);
    break;

  case Options::Mode::CLAUSIFY:
    clausifyMode(problem,false);
    break;

  case Options::Mode::TCLAUSIFY:
    clausifyMode(problem,true);
    break;

  case Options::Mode::OUTPUT:
    outputMode(problem);
    break;

  case Options::Mode::PROFILE:
    profileMode(problem);
    break;

  case Options::Mode::PREPROCESS:
  case Options::Mode::PREPROCESS2:
    preprocessMode(problem,false);
    break;

  case Options::Mode::TPREPROCESS:
    preprocessMode(problem,true);
    break;
  }
}

void interactiveMetamode()
{
  Options& opts = *env.options;
  opts.setInteractive(false); // so that we don't pass the interactivity on to the workers

  ScopedPtr<Problem> prb;
  if (!opts.inputFile().empty()) {
    UIHelper::parseFile(opts.inputFile(),opts.inputSyntax(),true);
    opts.resetInputFile();
  } // no parsing of the whole cin in interactiveMetamode
  prb = UIHelper::getInputProblem();

  while (true) {
    std::string line;
    if (!getline(cin, line) || line.rfind("exit",0) == 0) {
      cout << "Bye." << endl;
      break;
    } else if (line.rfind("run",0) == 0) {
      // the whole running happens in a child (don't modify our options, don't crash here when parsing option rubbish, etc.)
      pid_t process = Lib::Sys::Multiprocessing::instance()->fork();
      ASS_NEQ(process, -1);
      if(process == 0) {
        Timer::reinitialise(); // start our timer (in the child)
        UIHelper::unsetExpecting(); // probably garbage at this point

        Stack<std::string> pieces;
        StringUtils::splitStr(line.c_str(),' ',pieces);
        StringUtils::dropEmpty(pieces);
        Stack<const char*> argv(pieces.size());
        for(auto it = pieces.iterFifo(); it.hasNext();) {
          argv.push(it.next().c_str());
        }
        Shell::CommandLine cl(argv.size(), argv.begin());
        cl.interpret(opts);
        if (!opts.inputFile().empty()) {
          UIHelper::parseFile(opts.inputFile(),opts.inputSyntax(),true);
          prb = UIHelper::getInputProblem();
        }
        dispatchByMode(prb.ptr());
        exit(vampireReturnValue);
      }
    } else if (line.rfind("load",0) == 0) {
      Stack<std::string> pieces;
      StringUtils::splitStr(line.c_str(),' ',pieces);
      StringUtils::dropEmpty(pieces);
      auto it = pieces.iterFifo();
      ALWAYS(it.next() == "load");
      while (it.hasNext()) {
        UIHelper::parseFile(it.next(),opts.inputSyntax(),true);
      }
      prb = UIHelper::getInputProblem();
    } else if (line.rfind("tptp ",0) == 0) {
      try {
        UIHelper::parseSingleLine(line.substr(5),Options::InputSyntax::TPTP);
        prb = UIHelper::getInputProblem();
      } catch (ParsingRelatedException& exception) {
        explainException(exception);
      }
    } else if (line.rfind("smt2 ",0) == 0) {
      try {
        UIHelper::parseSingleLine(line.substr(5),Options::InputSyntax::SMTLIB2);
        prb = UIHelper::getInputProblem();
      } catch (ParsingRelatedException& exception) {
        explainException(exception);
      }
    } else if (line.rfind("list",0) == 0) {
      UIHelper::listLoadedPieces(cout);
    } else if (line.rfind("pop",0) == 0) {
      Stack<std::string> pieces;
      StringUtils::splitStr(line.c_str(),' ',pieces);
      StringUtils::dropEmpty(pieces);
      int numPops = 1;
      if (pieces.size() > 1) {
        Int::stringToInt(pieces[1],numPops);
      }
      UIHelper::popLoadedPiece(numPops);
      prb = UIHelper::getInputProblem();
    } else {
      cout << "Unreconginzed command! Try 'run [options] [filename_to_load]', 'load <filenames>', 'tptp <one_line_input_in_tptp>',\n"
              "'smt2 <one_line_input_in_smt2>' 'pop [how_many_levels] (one is default)', 'list', or 'exit'." << endl;
    }
  }
}

/**
 * The main function.
 * @since 03/12/2003 many changes related to logging
 *        and exception handling.
 * @since 10/09/2004, Manchester changed to use knowledge bases
 */
int main(int argc, char* argv[])
{
  System::setSignalHandlers();

  try {
    Options& opts = *env.options;

    // read the command line and interpret it
    Shell::CommandLine cl(argc, argv);
    cl.interpret(opts);

#if VDEBUG
    std::cerr << "% WARNING: debug build, do not use in anger\n";
#endif

    if(opts.encodeStrategy()){
      cout << opts.generateEncodedOptions() << "\n";
    }

#if VTIME_PROFILING
    TimeTrace::instance().setEnabled(opts.timeStatistics());
#endif

    // If any of these options are set then we just need to output and exit
    if (opts.showHelp() || opts.showOptions() || opts.showExperimentalOptions() ||
       !opts.explainOption().empty() || opts.printAllTheoryAxioms()) {
      opts.output(std::cout);
      exit(0);
    }

    Lib::setMemoryLimit(env.options->memoryLimit() * 1048576ul);

    if (opts.mode() == Options::Mode::MODEL_CHECK) {
      opts.setOutputAxiomNames(true);
    }

    if (opts.interactive()) {
      interactiveMetamode();
    } else {
      // can only happen after reading options as it relies on `env.options`
      Timer::reinitialise(); // start our timer, so that we also limit parsing

#if VAMPIRE_PERF_EXISTS
      unsigned saveInstrLimit = env.options->instructionLimit();
      if (env.options->parsingDoesNotCount()) {
        env.options->setInstructionLimit(0);
      }
#endif

      if (opts.inputFile().empty()) {
        UIHelper::parseStandardInput(opts.inputSyntax());
      } else {
        UIHelper::parseFile(opts.inputFile(),opts.inputSyntax(),
                            opts.mode() != Options::Mode::SPIDER && opts.mode() != Options::Mode::PROFILE);
      }

#if VAMPIRE_PERF_EXISTS
      if (env.options->parsingDoesNotCount()) {
        Timer::updateInstructionCount();
        unsigned burnedParsing = Timer::elapsedMegaInstructions();

        addCommentSignForSZS(std::cout);
        std::cout << "Instructions burned parsing: " << burnedParsing << " (million)" << endl;

        env.options->setInstructionLimit(saveInstrLimit+burnedParsing);
      }
#endif

      dispatchByMode(UIHelper::getInputProblem());
    }

#if CHECK_LEAKS
    delete env.signature;
    env.signature = 0;
#endif
  }
#if VZ3
  catch (z3::exception& exception) {
    vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;
    if (outputAllowed()) {
      cout << "Z3 exception:\n" << exception.msg() << endl;
    }
    reportSpiderFail();
  }
#endif
  catch (UserErrorException& exception) {
    vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;
    reportSpiderFail();
    explainException(exception);
  }
catch (Parse::TPTP::ParseErrorException& exception) {
    vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;
    reportSpiderFail();
    explainException(exception);
  }
  catch (Exception& exception) {
    vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;
    reportSpiderFail();
    explainException(exception);
  } catch (std::bad_alloc& _) {
    vampireReturnValue = VAMP_RESULT_STATUS_UNHANDLED_EXCEPTION;
    reportSpiderFail();
    std::cout << "Insufficient system memory" << '\n';
  }

  return vampireReturnValue;
} // main
