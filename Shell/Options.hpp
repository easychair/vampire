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
 * @file Options.hpp
 * Defines Vampire options.
 *
 * INSTRUCTIONS on Adding a new Option
 *
 * Firstly, the easiest thing to do is copy what's been done for an existing option
 *
 * In Options.hpp
 * - Add an OptionValue object (see NOTE on OptionValues below)
 * - Add enum for choices if ChoiceOptionValue
 * - Add getter for OptionValue
 * - Only if necessary (usually not), add setter for OptionValue
 *
 * In Options.cpp
 * - Initialise the OptionValue member, to do this you need to
 * -- Call the constructor with at least a long name, short name and default value
 * -- Provide a description
 * -- Insert the option into lookup (this is essential)
 * -- Tag the option, otherwise it will not appear nicely in showOptions
 * -- Add value constraints, they can be soft or hard (see NOTE on OptionValueConstraints below)
 * -- Add problem constraints (see NOTE on OptionProblemConstraints)
 *
 */

#ifndef __Options__
#define __Options__

#include <type_traits>
#include <cstring>
#include <memory>

#include "Forwards.hpp"

#include "Debug/Assertion.hpp"

#include "Lib/VirtualIterator.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/StringUtils.hpp"
#include "Lib/DArray.hpp"
#include "Lib/Stack.hpp"
#include "Lib/Int.hpp"
#include "Lib/Comparison.hpp"

#include "Property.hpp"

#define VAMPIRE_CLAUSE_TRACING VDEBUG

namespace Shell {

using namespace Lib;
using namespace Kernel;

class Property;

/**
 * Class that represents Vampire's options.
 * 11/11/2004 Shrigley Hall, completely reimplemented
 *
 * @since Sep 14 reimplemented by Giles
 */
class Options
{
public:

    Options ();
    // It is important that we can safely copy Options for use in CASC mode
    void init();
    void copyValuesFrom(const Options& that);
    Options(const Options& that);
    Options& operator=(const Options& that);

    // used to print help and options
    void output (std::ostream&) const;

    // Dealing with encoded options. Used by --decode option
    void readFromEncodedOptions (std::string testId);
    void readOptionsString (std::string testId,bool assign=true);
    std::string generateEncodedOptions() const;

    // deal with completeness
    bool complete(const Problem&) const;

    // deal with constraints
    void setForcedOptionValues(); // not currently used effectively
    bool checkGlobalOptionConstraints(bool fail_early=false);
    bool checkProblemOptionConstraints(Property*, bool before_preprocessing, bool fail_early=false);

    /**
     * Sample a random strategy from a distribution described by the given file.
     *
     * The format of the sampler file should be easy to understant (Look for examples samplerFOL.txt, samplerFNT.txt, samplerSMT.txt under vampire root).
     * The file describes a sequence of sampling rules (one on each line, barring empty lines and comment lines starting with a #),
     * which are executed in order, and each rule (provided its preconditions are satisfied) triggers
     * sampling of a value for a particular option from a specified distribution.
     * The most common sampler is for the categorical distribution (~cat), which is specified by a list of values with corresponding integer frequencies.
     * Other samplers include ratios, uniform floats and integers, and a (shifted) geometric distribution for potentially unbounded integers.
     * A notable feature is the ability to sample also fake (non-existent) options, recognized by a $-sign prefixed, whose value can later be reference in the conditions.
     *
     * Example:
     *
     * # naming
     * > $nm ~cat Z:1,NZ:5
     * $nm=Z > nm ~cat 0:1
     * $nm=NZ > nm ~sgd 0.07,2
     *
     * First samples a fake option $nm with either a value Z (1 out of 6) or Nz (5 out of 6).
     * Then, if $nm is set to Z samples the (actual) naming option with value 0 (nm ~cat 0:1 is simply an assignment to the effect of nm := 0),
     * and if $nm is set to NZ, samples from a shifted geometric distribution with p=0.07 and a shift=2. (So 2 gets selected with a probability p,
     * 3 with a probability p(1-p), ... and 2+i with a probability p(1-p)^i).
     */
    void sampleStrategy(const std::string& samplerFileName);

    /**
     * Return the problem name
     *
     * The problem name is computed from the input file name in
     * the @b setInputFile function. If the input file is not set,
     * the problem name is equal to "unknown". The problem name can
     * be set to a specific value using setProblemName().
     */
    const std::string& problemName () const { return _problemName.actualValue; }
    void setProblemName(std::string str) { _problemName.actualValue = str; }
    
    void setInputFile(const std::string& newVal){ _inputFile.set(newVal); }
    std::string includeFileName (const std::string& relativeName);

    // standard ways of creating options
    void set(const std::string& name, const std::string& value); // implicitly the long version used here
    void set(const char* name, const char* value, bool longOpt);

public:
  //==========================================================
  // The Enums for Option Values
  //==========================================================
  //
  // If you create a ChoiceOptionValue you will also need to create an enum


    /**
     * Possible tags to group options by
     * Update _tagNames at the end of Options constructor if you add a tag
     * @author Giles
     */
    enum class OptionTag: unsigned int {
        UNUSED,
        OTHER,
        DEVELOPMENT,
        OUTPUT,
        FMB,
        SAT,
        AVATAR,
        INFERENCES,
        INDUCTION,
        THEORIES,
        LRS,
        SATURATION,
        PREPROCESSING,
        INPUT,
        HELP,
        HIGHER_ORDER,
        LAST_TAG // Used for counting the number of tags
    };
    // update _tagNames at the end of Options constructor if you add a tag

  enum class TheoryInstSimp : unsigned int {
    OFF,
    ALL,    // select all interpreted
    STRONG, // select strong only
    NEG_EQ, // select only positive equalities
    OVERLAP,
    FULL,   // <-+- deprecated. only exists to not break portfolio modes. behaves exactly like `ALL` now
    NEW,    // <-+
  };
  enum class UnificationWithAbstraction : unsigned int {
    OFF,
    INTERP_ONLY,
    ONE_INTERP,
    CONSTANT,
    ALL,
    GROUND,
    FUNC_EXT,
    AC1,
    AC2,
  };

  enum class Induction : unsigned int {
    NONE,
    STRUCTURAL,
    INTEGER,
    BOTH
  };
  enum class StructuralInductionKind : unsigned int {
    ONE,
    TWO,
    THREE,
    RECURSION,
    ALL
  };
  enum class IntInductionKind : unsigned int {
    ONE,
    TWO
  };
  enum class IntegerInductionInterval : unsigned int {
    INFINITE,
    FINITE,
    BOTH
  };
  enum class IntegerInductionLiteralStrictness: unsigned int {
    NONE,
    TOPLEVEL_NOT_IN_OTHER,
    ONLY_ONE_OCCURRENCE,
    NOT_IN_BOTH,
    ALWAYS
  };
  enum class IntegerInductionTermStrictness: unsigned int {
    NONE,
    INTERPRETED_CONSTANT,
    NO_SKOLEMS
  };

  enum class PredicateSineLevels : unsigned int {
    NO,   // no means 1) the reverse of "on", 2) use with caution, it is predicted to be the worse value
    OFF,
    ON
  };


  enum class InductionChoice : unsigned int {
    ALL,
    GOAL,                     // only apply induction to goal constants
                              // a goal constant is one appearing in an explicit goal, or if gtg is used
                              // a constant that is used to lift a clause to a goal (uniqueness or Skolem)
    GOAL_PLUS,                // above plus skolem terms introduced in induction inferences
  };

  enum class DemodulationRedundancyCheck : unsigned int {
    OFF,       // no check
    ORDERING,  // solely ordering-based check
    ENCOMPASS, // (1) positive unit equations (PUE) are smaller than non-PUE clauses,
               // (2) more general PUE are smaller than less general PUE, and
               // (3) equally general PUEs are ordered based on term ordering.
  };

  enum class TheoryAxiomLevel : unsigned int {
    ON,  // all of them
    OFF, // none of them
    CHEAP
  };

  enum class ProofExtra : unsigned int {
    OFF,
    FREE,
    FULL
  };
  enum class FMBWidgetOrders : unsigned int {
    FUNCTION_FIRST, // f(1) f(2) f(3) ... g(1) g(2) ...
    ARGUMENT_FIRST, // f(1) g(1) h(1) ... f(2) g(2) ...
    DIAGONAL,       // f(1) g(2) h(3) f(2) g(3) h(1) f(3) g(1) h(2)
  };
  enum class FMBSymbolOrders : unsigned int {
    OCCURENCE,
    INPUT_USAGE,
    PREPROCESSED_USAGE
  };
  enum class FMBAdjustSorts : unsigned int {
    OFF,
    EXPAND,
    GROUP,
    PREDICATE,
    FUNCTION
  };
  enum class FMBEnumerationStrategy : unsigned int {
    SBMEAM,
#if VZ3
    SMT,
#endif
    CONTOUR
  };

  enum class BadOption : unsigned int {
    HARD,
    FORCED,
    OFF,
    SOFT
  };

  enum class IgnoreMissing : unsigned int {
    ON,
    OFF,
    WARN
  };

  /**
   * Possible values for function_definition_elimination.
   * @since 29/05/2004 Manchester
   */
  enum class FunctionDefinitionElimination : unsigned int {
    ALL = 0,
    NONE = 1,
    UNUSED = 2
  };

  /**
   *
   *
   */
  enum class Instantiation : unsigned int {
    OFF = 0,
    ON = 1
  };

  /**
   * Possible values for the input syntax
   * @since 26/08/2009 Redmond
   */
  enum class InputSyntax : unsigned int {
    SMTLIB2 = 0,
    /** syntax of the TPTP prover */
    TPTP = 1,
    AUTO = 2
    //HUMAN = 4,
    //MPS = 5,
    //NETLIB = 6
  };


  /**
   * Possible values for mode_name.
   * @since 06/05/2007 Manchester
   */
  enum class Mode : unsigned int {
    AXIOM_SELECTION,
    CASC,
    CLAUSIFY,
    CONSEQUENCE_ELIMINATION,
    MODEL_CHECK,
    /** this mode only outputs the input problem, without any preprocessing */
    OUTPUT,
    PORTFOLIO,
    PREPROCESS,
    PREPROCESS2,
    PROFILE,
    SMTCOMP,
    SPIDER,
    TCLAUSIFY,
    TPREPROCESS,
    VAMPIRE
  };

  enum class Intent : unsigned int {
    UNSAT, // preferentially look for refutations, proofs, arguments of unsatisfiability etc.
    SAT    // preferentially look for (finite) models, saturations, etc.
  };

  enum class Schedule : unsigned int {
    CASC,
    CASC_2024,
    CASC_SAT,
    CASC_SAT_2024,
    FILE,
    INDUCTION,
    INTEGER_INDUCTION,
    INTIND_OEIS,
    LTB_DEFAULT_2017,
    LTB_HH4_2017,
    LTB_HLL_2017,
    LTB_ISA_2017,
    LTB_MZR_2017,
    SMTCOMP,
    SMTCOMP_2018,
    SNAKE_TPTP_UNS,
    SNAKE_TPTP_SAT,
    STRUCT_INDUCTION,
    STRUCT_INDUCTION_TIP
  };

/* TODO: use an enum for Selection. The current issue is the way these values are manipulated as ints
 *
  enum class Selection : unsigned int {
    TOTAL,
    MAXIMAL,
    TWO,
    THREE,
    FOUR,
    TEN,
    LOOKAHEAD,
    BEST_TWO,
    BEST_THREE,
    BEST_FOUR,
    BEST_TEN,
    BEST_LOOKAHED
  }
*/

  /** Various options for the output of statistics in Vampire */
  enum class Statistics : unsigned int {
    /** changed by the option "--statistics brief" */
    BRIEF = 0,
    /** changed by the option "--statistics full */
    FULL = 1,
    /** changed by the option "--statistics off" */
    NONE = 2
  };

  /** how much we want vampire talking and in what language */
  enum class Output : unsigned int {
    SMTCOMP,
    SPIDER,
    SZS,
    VAMPIRE,
    UCORE
  };

  /** Possible values for sat_solver */
  enum class SatSolver : unsigned int {
     MINISAT = 0
#if VZ3
     ,Z3 = 1
#endif
  };

  /** Possible values for saturation_algorithm */
  enum class SaturationAlgorithm : unsigned int {
     DISCOUNT,
     FINITE_MODEL_BUILDING,
     LRS,
     OTTER,
     Z3
   };

  /** Possible values for activity of some inference rules */
  enum class RuleActivity : unsigned int {
    INPUT_ONLY = 0,
    OFF = 1,
    ON = 2
  };

  enum class QuestionAnsweringMode : unsigned int {
    PLAIN = 0,
    SYNTHESIS = 1,
    OFF = 2
  };

  enum class InterpolantMode : unsigned int {
    NEW_HEUR,
#if VZ3
    NEW_OPT,
#endif
    OFF,
  };

  enum class LiteralComparisonMode : unsigned int {
    PREDICATE = 0,
    REVERSE = 1,
    STANDARD = 2
  };

  enum class Condensation : unsigned int {
    FAST = 0,
    OFF = 1,
    ON = 2
  };

  enum class Demodulation : unsigned int {
    ALL = 0,
    OFF = 1,
    PREORDERED = 2
  };

  enum class Subsumption : unsigned int {
    OFF = 0,
    ON = 1,
    UNIT_ONLY = 2
  };

  enum class URResolution : unsigned int {
    EC_ONLY = 0,
    OFF = 1,
    ON = 2,
    FULL = 3
  };

  enum class TermOrdering : unsigned int {
    KBO = 0,
    LPO = 1
  };

  enum class SymbolPrecedence : unsigned int {
    ARITY = 0,
    OCCURRENCE = 1,
    REVERSE_ARITY = 2,
    UNARY_FIRST = 3,
    CONST_MAX = 4,
    CONST_MIN = 5,
    SCRAMBLE = 6,
    FREQUENCY = 7,
    UNARY_FREQ = 8,
    CONST_FREQ = 9,
    REVERSE_FREQUENCY = 10,
    WEIGHTED_FREQUENCY = 11,
    REVERSE_WEIGHTED_FREQUENCY = 12
  };
  enum class SymbolPrecedenceBoost : unsigned int {
    NONE = 0,
    GOAL = 1,
    UNIT = 2,
    GOAL_UNIT = 3,
    NON_INTRO = 4,
    INTRO = 5,
  };
  enum class IntroducedSymbolPrecedence : unsigned int {
    TOP = 0,
    BOTTOM = 1
  };

  enum class SineSelection : unsigned int {
    AXIOMS = 0,
    INCLUDED = 1,
    OFF = 2
  };

  enum class Proof : unsigned int {
    OFF = 0,
    ON = 1,
    PROOFCHECK = 2,
    TPTP = 3,
    PROPERTY = 4,
    DEDUKTI = 5
  };

  /** Values for --equality_proxy */
  enum class EqualityProxy : unsigned int {
    R = 0,
    RS = 1,
    RST = 2,
    RSTC = 3,
    OFF = 4,
  };

  /** Values for --extensionality_resolution */
  enum class ExtensionalityResolution : unsigned int {
    FILTER = 0,
    KNOWN = 1,
    TAGGED = 2,
    OFF = 3
  };

  enum class SplittingLiteralPolarityAdvice : unsigned int {
    FALSE,
    TRUE,
    NONE,
    RANDOM
  };

  enum class SplittingMinimizeModel : unsigned int {
    OFF = 0,
    SCO = 1,
    ALL = 2
  };

  enum class SplittingDeleteDeactivated : unsigned int {
    ON,
    LARGE_ONLY,
    OFF
  };

  enum class SplittingAddComplementary : unsigned int {
    GROUND = 0,
    NONE = 1
  };

  enum class SplittingCongruenceClosure : unsigned int {
    MODEL = 0,
    OFF = 1,
    ON = 2
  };

  enum class SplittingNonsplittableComponents : unsigned int {
    ALL = 0,
    ALL_DEPENDENT = 1,
    KNOWN = 2,
    NONE = 3
  };

  enum class TweeGoalTransformation : unsigned int {
    OFF = 0,
    GROUND = 1,
    FULL = 2
  };

  enum class CCUnsatCores : unsigned int {
    FIRST = 0,
    SMALL_ONES = 1,
    ALL = 2
  };

  enum class GlobalSubsumptionSatSolverPower : unsigned int {
    PROPAGATION_ONLY,
    FULL
  };

  enum class GlobalSubsumptionExplicitMinim : unsigned int {
    OFF,
    ON,
    RANDOMIZED
  };

  enum class GlobalSubsumptionAvatarAssumptions : unsigned int {
    OFF,
    FROM_CURRENT,
    FULL_MODEL
  };

  enum class Sos : unsigned int{
    ALL = 0,
    OFF = 1,
    ON = 2,
    THEORY = 3
  };

  enum class TARules : unsigned int {
    OFF = 0,
    INJECTGEN = 1,
    INJECTSIMPL = 2,
    INJECTOPT = 2,
    FULL = 3
  };

  enum class TACyclicityCheck : unsigned int {
    OFF = 0,
    AXIOM = 1,
    RULE = 2,
    RULELIGHT = 3
  };

  enum class GoalGuess : unsigned int {
    OFF = 0,
    ALL = 1,
    EXISTS_TOP = 2,
    EXISTS_ALL = 3,
    EXISTS_SYM = 4,
    POSITION = 5
  };

  enum class EvaluationMode : unsigned int {
    SIMPLE,
    POLYNOMIAL_FORCE,
    POLYNOMIAL_CAUTIOUS,
  };

  enum class ArithmeticSimplificationMode : unsigned int {
    FORCE,
    CAUTIOUS,
    OFF,
  };

  enum class KboWeightGenerationScheme : unsigned int {
    CONST = 0,
    RANDOM = 1,
    ARITY = 2,
    INV_ARITY = 3,
    ARITY_SQUARED = 4,
    INV_ARITY_SQUARED = 5,
    PRECEDENCE = 6,
    INV_PRECEDENCE = 7,
    FREQUENCY = 8,
    INV_FREQUENCY = 9,
  };

  enum class KboAdmissibilityCheck : unsigned int {
    ERROR = 0,
    WARNING = 1,
  };

  enum class FunctionExtensionality : unsigned int {
    OFF = 0,
    AXIOM = 1,
    ABSTRACTION = 2
  };

  enum class CNFOnTheFly : unsigned int {
    EAGER = 0,
    LAZY_GEN = 1,
    LAZY_SIMP = 2,
    LAZY_SIMP_NOT_GEN = 3,
    LAZY_SIMP_NOT_GEN_BOOL_EQ_OFF = 4,
    LAZY_SIMP_NOT_GEN_BOOL_EQ_GEN = 5,
    OFF = 6
  };

  enum class PISet : unsigned int {
    ALL = 0,
    ALL_EXCEPT_NOT_EQ = 1,
    FALSE_TRUE_NOT = 2,
    FALSE_TRUE_NOT_EQ_NOT_EQ = 3
  };

  enum class Narrow : unsigned int {
    ALL = 0,
    SK = 1,
    SKI = 2,
    OFF = 3
  };

  enum class ProblemExportSyntax : unsigned int {
    SMTLIB = 0,
    API_CALLS = 1,
  };

    //==========================================================
    // The Internals
    //==========================================================
    // Here I define the internal structures used to specify Options
    // Normally these are not modified, see below for getters and values
    //
    // The internals consist of
    // - OptionChoiceValues: to store the names of a option choice
    // - OptionValue: stores an options value and meta-data
    // - OptionValueConstraint: to give a constraint on an option
    // - OptionProblemConstraint: to give a constraint on an option wrt the problem
    //
    // The details are explained in comments below
private:
    // helper function of sampleStrategy
    void strategySamplingAssign(std::string optname, std::string value, DHMap<std::string,std::string>& fakes);
    std::string strategySamplingLookup(std::string optname, DHMap<std::string,std::string>& fakes);

    /**
     * These store the names of the choices for an option.
     * They can be declared using initializer lists i.e. {"on","off","half_on"}
     *
     * TODO: this uses a linear search, for alternative see NameArray
     *
     * @author Giles
     * @since 30/07/14
     */
    class OptionChoiceValues{
      void check_names_are_short() {
        for (auto x : _names) {
          ASS(x.size() < 70) // or else cannot be printed on a line
        }
      }
    public:
        OptionChoiceValues() : _names() { };
        OptionChoiceValues(Stack<std::string> names) : _names(std::move(names))  
        {
          check_names_are_short();
        }

        OptionChoiceValues(std::initializer_list<std::string> list) : _names(list)
        {
          check_names_are_short();
        }
        
        int find(std::string value) const {
            for(unsigned i=0;i<_names.length();i++){
                if(value.compare(_names[i])==0) return i;
            }
            return -1;
        }
        const int length() const { return _names.length(); }
        const std::string operator[](int i) const{ return _names[i];}

    private:
        Stack<std::string> _names;
    };

    // Declare constraints here so they can be referred to, but define them below
    template<typename T>
    struct OptionValueConstraint;
    template<typename T>
    using OptionValueConstraintUP = std::unique_ptr<OptionValueConstraint<T>>;
    struct AbstractWrappedConstraint;
    typedef std::unique_ptr<AbstractWrappedConstraint> AbstractWrappedConstraintUP;
    struct OptionProblemConstraint;
    typedef std::unique_ptr<OptionProblemConstraint> OptionProblemConstraintUP;

    /**
     * An AbstractOptionValue includes all the information and functionality that does not
     * depend on the type of the stored option. This is inherited by the templated OptionValue.
     *
     * The main purpose of the AbstractOptionValue is to allow us to have a collection of pointers
     * to OptionValue objects
     *
     * @author Giles
     */
    struct AbstractOptionValue {
        AbstractOptionValue(){}
        AbstractOptionValue(std::string l,std::string s) :
        longName(l), shortName(s), experimental(false), is_set(false),_should_copy(true), _tag(OptionTag::LAST_TAG), supress_problemconstraints(false) {}

        // Never copy an OptionValue... the Constraint system would break
        AbstractOptionValue(const AbstractOptionValue&) = delete;
        AbstractOptionValue& operator=(const AbstractOptionValue&) = delete;

        // however move-assigment is needed for all the assigns in Options::init()
        AbstractOptionValue(AbstractOptionValue&&) = default;
        AbstractOptionValue& operator= (AbstractOptionValue && ) = default;

        virtual ~AbstractOptionValue() = default;

        // This is the main method, it sets the value of the option using an input string
        // Returns false if we cannot set (will cause a UserError in Options::set)
        virtual bool setValue(const std::string& value) = 0;

        bool set(const std::string& value, bool dont_touch_if_defaulting = false) {
          bool okay = setValue(value);
          if (okay && (!dont_touch_if_defaulting || !isDefault())) {
            is_set=true;
          }
          return okay;
        }

        // Experimental options are not included in help
        void setExperimental(){experimental=true;}

        // Meta-data
        std::string longName;
        std::string shortName;
        std::string description;
        bool experimental;
        bool is_set;

        // Checking constraits
        virtual bool checkConstraints() = 0;
        virtual bool checkProblemConstraints(Property* prop) = 0;

        // Tagging: options can be filtered by mode and are organised by Tag in showOptions
        void tag(OptionTag tag){ ASS(_tag==OptionTag::LAST_TAG);_tag=tag; }
        void tag(Options::Mode mode){ _modes.push(mode); }

        OptionTag getTag(){ return _tag;}
        bool inMode(Options::Mode mode){
            if(_modes.isEmpty()) return true;
            else return _modes.find(mode);
        }

        // This allows us to get the actual value in string form
        virtual std::string getStringOfActual() const = 0;
        // Check if default value
        virtual bool isDefault() const = 0;

        // For use in showOptions and explainOption
        virtual void output(std::ostream& out,bool linewrap) const {
            out << "--" << longName;
            if(!shortName.empty()){ out << " (-"<<shortName<<")"; }
            out << std::endl;

            if (experimental) {
              out << "\t[experimental]" << std::endl;
            }


            if(!description.empty()){
                // Break a the description into lines where there have been at least 70 characters
                // on the line at the next space
                out << "\t";
                int count=0;
                for(const char* p = description.c_str();*p;p++){
                    out << *p;
                    count++;
                    if(linewrap && count>70 && *p==' '){
                        out << std::endl << '\t';
                        count=0;
                    }
                    if(*p=='\n'){ count=0; out << '\t'; }
                }
                out << std::endl;
            }
            else{ out << "\tno description provided!" << std::endl; }
        }

        // Used to determine wheter the value of an option should be copied when
        // the Options object is copied.
        bool _should_copy;
        bool shouldCopy() const { return _should_copy; }
       
        typedef std::unique_ptr<DArray<std::string>> stringDArrayUP;

        typedef std::pair<OptionProblemConstraintUP,stringDArrayUP> RandEntry;
 
    private:
        // Tag state
        OptionTag _tag;
        Lib::Stack<Options::Mode> _modes;

        stringDArrayUP toArray(std::initializer_list<std::string>& list){
          DArray<std::string>* array = new DArray<std::string>(list.size());
          unsigned index=0;
          for(typename std::initializer_list<std::string>::iterator it = list.begin();
           it!=list.end();++it){ (*array)[index++] =*it; }
          return stringDArrayUP(array);
        }
    protected:
        // Note has LIFO semantics so use BottomFirstIterator
        bool supress_problemconstraints;
    };

    struct AbstractOptionValueCompatator{
      Comparison compare(AbstractOptionValue* o1, AbstractOptionValue* o2)
      {
        int value = strcmp(o1->longName.c_str(),o2->longName.c_str());
        return value < 0 ? LESS : (value==0 ? EQUAL : GREATER);
      }
    };

    /**
     * The templated OptionValue is used to store default and actual values for options
     *
     * There are also type-related helper functions
     *
     * @author Giles
     */
    template<typename T>
    struct OptionValue : public AbstractOptionValue {
        // We need to include an empty constructor as all the OptionValue objects need to be initialized
        // with something when the Options object is created. They should then all be reconstructed
        // This is annoying but preferable to the alternative in my opinion
        OptionValue(){}
        OptionValue(std::string l, std::string s,T def) : AbstractOptionValue(l,s),
        defaultValue(def), actualValue(def){}

        // We store the defaultValue separately so that we can check if the actualValue is non-default
        T defaultValue;
        T actualValue;

        virtual bool isDefault() const { return defaultValue==actualValue;}

        // Getting the string versions of values, useful for output
        virtual std::string getStringOfValue(T value) const{ ASSERTION_VIOLATION;}
        virtual std::string getStringOfActual() const { return getStringOfValue(actualValue); }
        
        // Adding and checking constraints
        // By default constraints are soft and reaction to them is controlled by the bad_option option
        // But a constraint can be added as Hard, meaning that it always causes a UserError
        void addConstraint(OptionValueConstraintUP<T> c){ _constraints.push(std::move(c)); }
        void addHardConstraint(OptionValueConstraintUP<T> c){ c->setHard();addConstraint(std::move(c)); }

        // A onlyUsefulWith constraint gives a constraint that must be true if this option's value is set
        // For example, split_at_activation is only useful with splitting being on
        // These are defined for OptionValueConstraints and WrappedConstraints - see below for explanation
        void onlyUsefulWith(AbstractWrappedConstraintUP c){
            _constraints.push(If(hasBeenSet<T>()).then(unwrap<T>(c)));
        }
        void onlyUsefulWith(OptionValueConstraintUP<T> c){
            _constraints.push(If(hasBeenSet<T>()).then(std::move(c)));
        }

        // similar to onlyUsefulWith, except the trigger is a non-default value
        // (as opposed to the explicitly-set flag)
        // we use it for selection and awr which cannot be not set via the decode string
        void onlyUsefulWith2(AbstractWrappedConstraintUP c){
            _constraints.push(If(getNotDefault()).then(unwrap<T>(c)));
        }
        void onlyUsefulWith2(OptionValueConstraintUP<T> c){
            _constraints.push(If(getNotDefault()).then(std::move(c)));
        }

        virtual OptionValueConstraintUP<T> getNotDefault(){ return isNotDefault<T>(); }

        // similar to onlyUsefulWith2, except its a hard constraint,
        // so that the user is strongly aware of situations when changing the
        // respective option has no effect
        void reliesOn(AbstractWrappedConstraintUP c){
            OptionValueConstraintUP<T> tc = If(getNotDefault()).then(unwrap<T>(c));
            tc->setHard();
            _constraints.push(std::move(tc));
        }
        void reliesOn(OptionValueConstraintUP<T> c){
            OptionValueConstraintUP<T> tc = If(getNotDefault()).then(c);
            tc->setHard();
            _constraints.push(std::move(tc));
        }
        // This checks the constraints and may cause a UserError
        bool checkConstraints();

        // Produces a separate constraint object based on this option
        /// Useful for IfThen constraints and onlyUsefulWith i.e. _splitting.is(equal(true))
        AbstractWrappedConstraintUP is(OptionValueConstraintUP<T> c);

        // Problem constraints place a restriction on problem properties and option values
        void addProblemConstraint(OptionProblemConstraintUP c){ _prob_constraints.push(std::move(c)); }
        bool hasProblemConstraints(){
          return !supress_problemconstraints && !_prob_constraints.isEmpty();
        }
        virtual bool checkProblemConstraints(Property* prop);

        virtual void output(std::ostream& out, bool linewrap) const {
            AbstractOptionValue::output(out,linewrap);
            out << "\tdefault: " << getStringOfValue(defaultValue) << std::endl;
        }

    private:
        Lib::Stack<OptionValueConstraintUP<T>> _constraints;
        Lib::Stack<OptionProblemConstraintUP> _prob_constraints;
    };

    /**
     * We now define particular OptionValues, see NOTE on OptionValues for high level usage
     */

    /**
     * A ChoiceOptionValue is templated by an enum, which must be defined above
     *
     * It is then necessary to provide names for the enum values.
     * We do not check that those names have the same length as the enum but this is very important.
     * The names must also be in the same order!
     *
     * @author Giles
     */
    template<typename T >
    struct ChoiceOptionValue : public OptionValue<T> {
        ChoiceOptionValue(){}
        ChoiceOptionValue(std::string l, std::string s,T def,OptionChoiceValues c) :
        OptionValue<T>(l,s,def), choices(c) {}
        ChoiceOptionValue(std::string l, std::string s,T d) : ChoiceOptionValue(l,s,d, T::optionChoiceValues()) {}
        
        bool setValue(const std::string& value){
            // makes reasonable assumption about ordering of every enum
            int index = choices.find(value.c_str());
            if(index<0) return false;
            this->actualValue = static_cast<T>(index);
            return true;
        }

        virtual void output(std::ostream& out,bool linewrap) const {
            AbstractOptionValue::output(out,linewrap);
            out << "\tdefault: " << choices[static_cast<unsigned>(this->defaultValue)];
            out << std::endl;
            std::string values_header = "values: ";
            out << "\t" << values_header;
            // Again we restrict line length to 70 characters
            int count=0;
            for(int i=0;i<choices.length();i++){
                if(i==0){
                    out << choices[i];
                }
                else{
                    out << ",";
                    std::string next = choices[i];
                    if(linewrap && next.size()+count>60){ // next.size() will be <70, how big is a tab?
                        out << std::endl << "\t";
                        for(unsigned j=0;j<values_header.size();j++){out << " ";}
                        count = 0;
                    }
                    out << next;
                    count += next.size();
                }
            }
            out << std::endl;
        }
        
        std::string getStringOfValue(T value) const {
            unsigned i = static_cast<unsigned>(value);
            return choices[i];
        }

    private:
        OptionChoiceValues choices;
    };


    /**
     * For Booleans - we use on/off rather than true/false
     * @author Giles
     */
    struct BoolOptionValue : public OptionValue<bool> {
        BoolOptionValue(){}
        BoolOptionValue(std::string l,std::string s, bool d) : OptionValue(l,s,d){}
        bool setValue(const std::string& value){
            if (! value.compare("on") || ! value.compare("true")) {
                actualValue=true;

            }
            else if (! value.compare("off") || ! value.compare("false")) {
                actualValue=false;
            }
            else return false;

            return true;
        }
        
        std::string getStringOfValue(bool value) const { return (value ? "on" : "off"); }
    };

    struct IntOptionValue : public OptionValue<int> {
        IntOptionValue(){}
        IntOptionValue(std::string l,std::string s, int d) : OptionValue(l,s,d){}
        bool setValue(const std::string& value){
            return Int::stringToInt(value.c_str(),actualValue);
        }
        std::string getStringOfValue(int value) const{ return Lib::Int::toString(value); }
    };

    struct UnsignedOptionValue : public OptionValue<unsigned> {
        UnsignedOptionValue(){}
        UnsignedOptionValue(std::string l,std::string s, unsigned d) : OptionValue(l,s,d){}

        bool setValue(const std::string& value){
            return Int::stringToUnsignedInt(value.c_str(),actualValue);
        }
        std::string getStringOfValue(unsigned value) const{ return Lib::Int::toString(value); }
    };
    
    struct StringOptionValue : public OptionValue<std::string> {
        StringOptionValue(){}
        StringOptionValue(std::string l,std::string s, std::string d) : OptionValue(l,s,d){}
        bool setValue(const std::string& value){
            actualValue = (value=="<empty>") ? "" : value;
            return true;
        }
        std::string getStringOfValue(std::string value) const{
            if(value.empty()) return "<empty>";
            return value;
        }
    };

    struct LongOptionValue : public OptionValue<long> {
        LongOptionValue(){}
        LongOptionValue(std::string l,std::string s, long d) : OptionValue(l,s,d){}
        bool setValue(const std::string& value){
            return Int::stringToLong(value.c_str(),actualValue);
        }
        std::string getStringOfValue(long value) const{ return Lib::Int::toString(value); }
    };

struct FloatOptionValue : public OptionValue<float>{
FloatOptionValue(){}
FloatOptionValue(std::string l,std::string s, float d) : OptionValue(l,s,d){}
bool setValue(const std::string& value){
    return Int::stringToFloat(value.c_str(),actualValue);
}
std::string getStringOfValue(float value) const{ return Lib::Int::toString(value); }
};

/**
* Ratios have two actual values and two default values
* Therefore, we often need to tread them specially
* @author Giles
*/
struct RatioOptionValue : public OptionValue<int> {
RatioOptionValue(){}
RatioOptionValue(std::string l, std::string s, int def, int other, char sp=':') :
OptionValue(l,s,def), sep(sp), defaultOtherValue(other), otherValue(other) {};

virtual OptionValueConstraintUP<int> getNotDefault() override { return isNotDefaultRatio(); }

virtual bool isDefault() const override { return defaultValue * otherValue == actualValue * defaultOtherValue; }

void addConstraintIfNotDefault(AbstractWrappedConstraintUP c){
    addConstraint(If(isNotDefaultRatio()).then(unwrap<int>(c)));
}

bool readRatio(const char* val,char seperator);
bool setValue(const std::string& value) override {
    return readRatio(value.c_str(),sep);
}

char sep;
int defaultOtherValue;
int otherValue;

virtual void output(std::ostream& out,bool linewrap) const override {
    AbstractOptionValue::output(out,linewrap);
    out << "\tdefault left: " << defaultValue << std::endl;
    out << "\tdefault right: " << defaultOtherValue << std::endl;
}

virtual std::string getStringOfValue(int value) const override { ASSERTION_VIOLATION;}
virtual std::string getStringOfActual() const override {
    return Lib::Int::toString(actualValue)+sep+Lib::Int::toString(otherValue);
}

};

// We now have a number of option-specific values
// These are necessary when the option needs to be read in a special way

/**
* Oddly gets set with a float value and then creates a ratio of value*100/100
* @author Giles
*/
struct NonGoalWeightOptionValue : public OptionValue<float>{
NonGoalWeightOptionValue(){}
NonGoalWeightOptionValue(std::string l, std::string s, float def) :
OptionValue(l,s,def), numerator(1), denominator(1) {};

bool setValue(const std::string& value);

// output does not output numerator and denominator as they
// are produced from defaultValue
int numerator;
int denominator;

virtual std::string getStringOfValue(float value) const{ return Lib::Int::toString(value); }
};

/**
* Selection is defined by a set of integers (TODO: make enum)
* For now we need to check the integer is a valid one
* @author Giles
*/
struct SelectionOptionValue : public OptionValue<int>{
SelectionOptionValue(){}
SelectionOptionValue(std::string l,std::string s, int def):
OptionValue(l,s,def){};

bool setValue(const std::string& value);

virtual void output(std::ostream& out,bool linewrap) const {
    AbstractOptionValue::output(out,linewrap);
    out << "\tdefault: " << defaultValue << std::endl;;
}

virtual std::string getStringOfValue(int value) const{ return Lib::Int::toString(value); }

AbstractWrappedConstraintUP isLookAheadSelection(){
  return AbstractWrappedConstraintUP(new WrappedConstraint<int>(*this,OptionValueConstraintUP<int>(new isLookAheadSelectionConstraint())));
}
};

/**
* This also updates problemName
* @author Giles
*/
struct InputFileOptionValue : public OptionValue<std::string>{
InputFileOptionValue(){}
InputFileOptionValue(std::string l,std::string s, std::string def,Options* p):
OptionValue(l,s,def), parent(p){};

bool setValue(const std::string& value);

virtual void output(std::ostream& out,bool linewrap) const {
    AbstractOptionValue::output(out,linewrap);
    out << "\tdefault: " << defaultValue << std::endl;;
}
virtual std::string getStringOfValue(std::string value) const{ return value; }
private:
Options* parent;

};
/**
* We need to decode the encoded option string
* @author Giles
*/
struct DecodeOptionValue : public OptionValue<std::string>{
DecodeOptionValue(){ AbstractOptionValue::_should_copy=false;}
DecodeOptionValue(std::string l,std::string s,Options* p):
OptionValue(l,s,""), parent(p){ AbstractOptionValue::_should_copy=false;}

bool setValue(const std::string& value){
    parent->readFromEncodedOptions(value);
    return true;
}
virtual std::string getStringOfValue(std::string value) const{ return value; }

private:
Options* parent;

};
/**
* Need to read the time limit. By default it assumes seconds (and stores deciseconds) but you can give
* a multiplier i.e. d,s,m,h,D for deciseconds,seconds,minutes,hours,Days
* @author Giles
*/
struct TimeLimitOptionValue : public OptionValue<int>{
TimeLimitOptionValue(){}
TimeLimitOptionValue(std::string l, std::string s, float def) :
OptionValue(l,s,def) {};

bool setValue(const std::string& value);

virtual void output(std::ostream& out,bool linewrap) const {
    AbstractOptionValue::output(out,linewrap);
    out << "\tdefault: " << defaultValue << "d" << std::endl;
}
virtual std::string getStringOfValue(int value) const{ return Lib::Int::toString(value)+"d"; }
};

/**
* NOTE on OptionValueConstraints
*
* OptionValueConstraints are used to declare constraints on and between option values
* these are checked in checkGlobalOptionConstraints, which should be called after
* Options is updated
*
* As usual, see Options.cpp for examples.
*
* There are two kinds of ValueConstraints (see below for ProblemConstraints)
*
* - Unary constraints such as greaterThan, equals, ...
* - If-then constraints that capture dependencies
*
* In both cases an attempt has been made to make the declaration of constraints
* in Options.cpp as readable as possible. For example, an If-then constraint is
* written as follows
*
*  If(equals(0)).then(_otherOption.is(lessThan(5)))
*
* Note that the equals(0) will apply to the OptionValue that the constraint belongs to
*
* WrappedConstraints are produced by OptionValue.is and are used to provide constraints
* on other OptionValues, as seen in the example above. Most functions work with both
* OptionValueConstraint and WrappedConstraint but in some cases one of these options
* may need to be added. In this case see examples from AndWrapper below.
*
* MS: While OptionValueConstraints are expressions which wait for a concrete value to be evaluated against:
* as in λ value. expression(value),
* WrappedConstraints have already been "closed" by providing a concrete value:
* as in (λ value. expression(value))[concrete_value]
* Finally, we can at anytime "unwrap" a WrappedConstraint by providing a "fake" lambda again on top, to turn it into a OptionValueConstraints again:
* as in λ value. expression_ignoring_value
*
* The tricky part (C++-technology-wise) here is that unwrapping needs to get a type for the value
* and this type is independent form the expression_ignoring_value for obvious reasons.
* So various overloads of things are needed until we get to the point, where the type is known and can be supplied.
* (e.g. there needs to be a separate hierarchy of Wrapped expressions along the one for OptionValueConstraint ones).
*/

template<typename T>
struct OptionValueConstraint{
OptionValueConstraint() : _hard(false) {}

virtual ~OptionValueConstraint() {} // virtual methods present -> there should be virtual destructor

virtual bool check(const OptionValue<T>& value) = 0;
virtual std::string msg(const OptionValue<T>& value) = 0;

// By default cannot force constraint
virtual bool force(OptionValue<T>* value){ return false;}
// TODO - allow for hard constraints
bool isHard(){ return _hard; }
void setHard(){ _hard=true;}
bool _hard;
};

    // A Wrapped Constraint takes an OptionValue and a Constraint
    // It allows us to supply a constraint on another OptionValue in an If constraint for example
    struct AbstractWrappedConstraint {
      virtual bool check() = 0;
      virtual std::string msg() = 0;
      virtual ~AbstractWrappedConstraint() {};
    };

    template<typename T>
    struct WrappedConstraint : AbstractWrappedConstraint {
        WrappedConstraint(const OptionValue<T>& v, OptionValueConstraintUP<T> c) : value(v), con(std::move(c)) {}

        bool check() override {
            return con->check(value);
        }
        std::string msg() override {
            return con->msg(value);
        }

        const OptionValue<T>& value;
        OptionValueConstraintUP<T> con;
    };

    struct WrappedConstraintOrWrapper : public AbstractWrappedConstraint {
        WrappedConstraintOrWrapper(AbstractWrappedConstraintUP l, AbstractWrappedConstraintUP r) : left(std::move(l)),right(std::move(r)) {}
        bool check() override {
            return left->check() || right->check();
        }
        std::string msg() override { return left->msg() + " or " + right->msg(); }

        AbstractWrappedConstraintUP left;
        AbstractWrappedConstraintUP right;
    };

    struct WrappedConstraintAndWrapper : public AbstractWrappedConstraint {
        WrappedConstraintAndWrapper(AbstractWrappedConstraintUP l, AbstractWrappedConstraintUP r) : left(std::move(l)),right(std::move(r)) {}
        bool check() override {
            return left->check() && right->check();
        }
        std::string msg() override { return left->msg() + " and " + right->msg(); }

        AbstractWrappedConstraintUP left;
        AbstractWrappedConstraintUP right;
    };

    template<typename T>
    struct OptionValueConstraintOrWrapper : public OptionValueConstraint<T>{
        OptionValueConstraintOrWrapper(OptionValueConstraintUP<T> l, OptionValueConstraintUP<T> r) : left(std::move(l)),right(std::move(r)) {}
        bool check(const OptionValue<T>& value){
            return left->check(value) || right->check(value);
        }
        std::string msg(const OptionValue<T>& value){ return left->msg(value) + " or " + right->msg(value); }

        OptionValueConstraintUP<T> left;
        OptionValueConstraintUP<T> right;
    };

    template<typename T>
    struct OptionValueConstraintAndWrapper : public OptionValueConstraint<T>{
        OptionValueConstraintAndWrapper(OptionValueConstraintUP<T> l, OptionValueConstraintUP<T> r) : left(std::move(l)),right(std::move(r)) {}
        bool check(const OptionValue<T>& value){
            return left->check(value) && right->check(value);
        }
        std::string msg(const OptionValue<T>& value){ return left->msg(value) + " and " + right->msg(value); }

        OptionValueConstraintUP<T> left;
        OptionValueConstraintUP<T> right;
    };

    template<typename T>
    struct UnWrappedConstraint : public OptionValueConstraint<T>{
        UnWrappedConstraint(AbstractWrappedConstraintUP c) : con(std::move(c)) {}

        bool check(const OptionValue<T>&){ return con->check(); }
        std::string msg(const OptionValue<T>&){ return con->msg(); }
        
        AbstractWrappedConstraintUP con;
    };

    template <typename T>
    static OptionValueConstraintUP<T> maybe_unwrap(OptionValueConstraintUP<T> c) { return c; }

    template <typename T>
    static OptionValueConstraintUP<T> unwrap(AbstractWrappedConstraintUP& c) { return OptionValueConstraintUP<T>(new UnWrappedConstraint<T>(std::move(c))); }

    template <typename T>
    static OptionValueConstraintUP<T> maybe_unwrap(AbstractWrappedConstraintUP& c) { return unwrap<T>(c); }

    /*
     * To avoid too many cases a certain discipline is required from the user.
     * Namely, OptionValueConstraints need to precede WrappedConstraints in the arguments of Or and And
     **/

    // the base case (the unary Or)
    template <typename T>
    OptionValueConstraintUP<T> Or(OptionValueConstraintUP<T> a) { return a; }
    AbstractWrappedConstraintUP Or(AbstractWrappedConstraintUP a) { return a; }

    template<typename T, typename... Args>
    OptionValueConstraintUP<T> Or(OptionValueConstraintUP<T> a, Args... args)
    {
      OptionValueConstraintUP<T> r = maybe_unwrap<T>(Or(std::move(args)...));
      return OptionValueConstraintUP<T>(new OptionValueConstraintOrWrapper<T>(std::move(a),std::move(r)));
    }

    template<typename... Args>
    AbstractWrappedConstraintUP Or(AbstractWrappedConstraintUP a, Args... args)
    {
      AbstractWrappedConstraintUP r = Or(std::move(args)...);
      return AbstractWrappedConstraintUP(new WrappedConstraintOrWrapper(std::move(a),std::move(r)));
    }

    // the base case (the unary And)
    template <typename T>
    OptionValueConstraintUP<T> And(OptionValueConstraintUP<T> a) { return a; }
    AbstractWrappedConstraintUP And(AbstractWrappedConstraintUP a) { return a; }

    template<typename T, typename... Args>
    OptionValueConstraintUP<T> And(OptionValueConstraintUP<T> a, Args... args)
    {
      OptionValueConstraintUP<T> r = maybe_unwrap<T>(And(std::move(args)...));
      return OptionValueConstraintUP<T>(new OptionValueConstraintAndWrapper<T>(std::move(a),std::move(r)));
    }

    template<typename... Args>
    AbstractWrappedConstraintUP And(AbstractWrappedConstraintUP a, Args... args)
    {
      AbstractWrappedConstraintUP r = And(std::move(args)...);
      return AbstractWrappedConstraintUP(new WrappedConstraintAndWrapper(std::move(a),std::move(r)));
    }

    template<typename T>
    struct Equal : public OptionValueConstraint<T>{
        Equal(T gv) : _goodvalue(gv) {}
        bool check(const OptionValue<T>& value){
            return value.actualValue == _goodvalue;
        }
        std::string msg(const OptionValue<T>& value){
            return value.longName+"("+value.getStringOfActual()+") is equal to " + value.getStringOfValue(_goodvalue);
        }
        T _goodvalue;
    };
    template<typename T>
    static OptionValueConstraintUP<T> equal(T bv){
        return OptionValueConstraintUP<T>(new Equal<T>(bv));
    }

    template<typename T>
    struct NotEqual : public OptionValueConstraint<T>{
        NotEqual(T bv) : _badvalue(bv) {}
        bool check(const OptionValue<T>& value){
            return value.actualValue != _badvalue;
        }
        std::string msg(const OptionValue<T>& value){ return value.longName+"("+value.getStringOfActual()+") is not equal to " + value.getStringOfValue(_badvalue); }
        T _badvalue;
    };
    template<typename T>
    static OptionValueConstraintUP<T> notEqual(T bv){
        return OptionValueConstraintUP<T>(new NotEqual<T>(bv));
    }

    // Constraint that the value should be less than a given value
    // optionally we can allow it be equal to that value also
    template<typename T>
    struct LessThan : public OptionValueConstraint<T>{
        LessThan(T gv,bool eq=false) : _goodvalue(gv), _orequal(eq) {}
        bool check(const OptionValue<T>& value){
            return (value.actualValue < _goodvalue || (_orequal && value.actualValue==_goodvalue));
        }
        std::string msg(const OptionValue<T>& value){
            if(_orequal) return value.longName+"("+value.getStringOfActual()+") is less than or equal to " + value.getStringOfValue(_goodvalue);
            return value.longName+"("+value.getStringOfActual()+") is less than "+ value.getStringOfValue(_goodvalue);
        }

        T _goodvalue;
        bool _orequal;
    };
    template<typename T>
    static OptionValueConstraintUP<T> lessThan(T bv){
        return OptionValueConstraintUP<T>(new LessThan<T>(bv,false));
    }
    template<typename T>
    static OptionValueConstraintUP<T> lessThanEq(T bv){
        return OptionValueConstraintUP<T>(new LessThan<T>(bv,true));
    }

    // Constraint that the value should be greater than a given value
    // optionally we can allow it be equal to that value also
    template<typename T>
    struct GreaterThan : public OptionValueConstraint<T>{
        GreaterThan(T gv,bool eq=false) : _goodvalue(gv), _orequal(eq) {}
        bool check(const OptionValue<T>& value){
            return (value.actualValue > _goodvalue || (_orequal && value.actualValue==_goodvalue));
        }
        
        std::string msg(const OptionValue<T>& value){
            if(_orequal) return value.longName+"("+value.getStringOfActual()+") is greater than or equal to " + value.getStringOfValue(_goodvalue);
            return value.longName+"("+value.getStringOfActual()+") is greater than "+ value.getStringOfValue(_goodvalue);
        }

        T _goodvalue;
        bool _orequal;
    };
    template<typename T>
    static OptionValueConstraintUP<T> greaterThan(T bv){
        return OptionValueConstraintUP<T>(new GreaterThan<T>(bv,false));
    }
    template<typename T>
    static OptionValueConstraintUP<T> greaterThanEq(T bv){
        return OptionValueConstraintUP<T>(new GreaterThan<T>(bv,true));
    }

    // Constraint that the value should be smaller than a given value
    // optionally we can allow it be equal to that value also
    template<typename T>
    struct SmallerThan : public OptionValueConstraint<T>{
        SmallerThan(T gv,bool eq=false) : _goodvalue(gv), _orequal(eq) {}
        bool check(const OptionValue<T>& value){
            return (value.actualValue < _goodvalue || (_orequal && value.actualValue==_goodvalue));
        }

        std::string msg(const OptionValue<T>& value){
            if(_orequal) return value.longName+"("+value.getStringOfActual()+") is smaller than or equal to " + value.getStringOfValue(_goodvalue);
            return value.longName+"("+value.getStringOfActual()+") is smaller than "+ value.getStringOfValue(_goodvalue);
        }

        T _goodvalue;
        bool _orequal;
    };
    template<typename T>
    static OptionValueConstraintUP<T> smallerThan(T bv){
        return OptionValueConstraintUP<T>(new SmallerThan<T>(bv,false));
    }
    template<typename T>
    static OptionValueConstraintUP<T> smallerThanEq(T bv){
        return OptionValueConstraintUP<T>(new SmallerThan<T>(bv,true));
    }

    /**
     * If constraints
     */

    template<typename T>
    struct IfConstraint;

    template<typename T>
    struct IfThenConstraint : public OptionValueConstraint<T>{
        IfThenConstraint(OptionValueConstraintUP<T> ic, OptionValueConstraintUP<T> c) :
        if_con(std::move(ic)), then_con(std::move(c)) {}

        bool check(const OptionValue<T>& value){
            ASS(then_con);
            return !if_con->check(value) || then_con->check(value);
        }
        
        std::string msg(const OptionValue<T>& value){
            return "if "+if_con->msg(value)+" then "+ then_con->msg(value);
        }

        OptionValueConstraintUP<T> if_con;
        OptionValueConstraintUP<T> then_con;
    };

    template<typename T>
    struct IfConstraint {
        IfConstraint(OptionValueConstraintUP<T> c) :if_con(std::move(c)) {}

        OptionValueConstraintUP<T> then(OptionValueConstraintUP<T> c){
          return OptionValueConstraintUP<T>(new IfThenConstraint<T>(std::move(if_con),std::move(c)));
        }
        OptionValueConstraintUP<T> then(AbstractWrappedConstraintUP c){
          return OptionValueConstraintUP<T>(new IfThenConstraint<T>(std::move(if_con),unwrap<T>(c)));
        }

        OptionValueConstraintUP<T> if_con;
    };

    template<typename T>
    static IfConstraint<T> If(OptionValueConstraintUP<T> c){
        return IfConstraint<T>(std::move(c));
    }
    template<typename T>
    static IfConstraint<T> If(AbstractWrappedConstraintUP c){
        return IfConstraint<T>(unwrap<T>(c));
    }

    /**
     * Option-(explicitly)-set constraint
     */
    template<typename T>
    struct HasBeenSet : public OptionValueConstraint<T> {
      HasBeenSet() {}

        bool check(const OptionValue<T>& value) override {
            return value.is_set;
        }
        std::string msg(const OptionValue<T>& value) override { return value.longName+"("+value.getStringOfActual()+") has been set";}
    };

    template<typename T>
    static OptionValueConstraintUP<T> hasBeenSet(){
        return OptionValueConstraintUP<T>(new HasBeenSet<T>());
    }

    /**
     * Default Value constraints
     */
    template<typename T>
    struct NotDefaultConstraint : public OptionValueConstraint<T> {
        NotDefaultConstraint() {}

        bool check(const OptionValue<T>& value){
            return value.defaultValue != value.actualValue;
        }
        std::string msg(const OptionValue<T>& value) { return value.longName+"("+value.getStringOfActual()+") is not default("+value.getStringOfValue(value.defaultValue)+")";}
    };
    struct NotDefaultRatioConstraint : public OptionValueConstraint<int> {
        NotDefaultRatioConstraint() {}

        bool check(const OptionValue<int>& value){
            const RatioOptionValue& rvalue = static_cast<const RatioOptionValue&>(value);
            return (rvalue.defaultValue != rvalue.actualValue ||
                    rvalue.defaultOtherValue != rvalue.otherValue);
        }
        std::string msg(const OptionValue<int>& value) { return value.longName+"("+value.getStringOfActual()+") is not default";}
        
    };

    // You will need to provide the type, optionally use addConstraintIfNotDefault
    template<typename T>
    static OptionValueConstraintUP<T> isNotDefault(){
        return OptionValueConstraintUP<T>(new NotDefaultConstraint<T>());
    }
    // You will need to provide the type, optionally use addConstraintIfNotDefault
    static OptionValueConstraintUP<int> isNotDefaultRatio(){
        return OptionValueConstraintUP<int>(new NotDefaultRatioConstraint());
    }

    struct isLookAheadSelectionConstraint : public OptionValueConstraint<int>{
        isLookAheadSelectionConstraint() {}
        bool check(const OptionValue<int>& value){
            return value.actualValue == 11 || value.actualValue == 1011 || value.actualValue == -11 || value.actualValue == -1011;
        }
        std::string msg(const OptionValue<int>& value){
            return value.longName+"("+value.getStringOfActual()+") is not lookahead selection";
        }
    };


    /**
     * NOTE on OptionProblemConstraint
     *
     * OptionProblemConstraints are used to capture properties of a problem that
     * should be present when an option is used. The idea being that a warning will
     * be emitted if an option is used for an inappropriate problem.
     *
     * TODO - this element of Options is still under development
     */

    struct OptionProblemConstraint{
      virtual bool check(Property* p) = 0;
      virtual std::string msg() = 0;
      virtual ~OptionProblemConstraint() {};
    };

    struct CategoryCondition : OptionProblemConstraint{
      CategoryCondition(Property::Category c,bool h) : cat(c), has(h) {}
      bool check(Property*p){
          ASS(p);
          return has ? p->category()==cat : p->category()!=cat;
      }
      std::string msg(){
        std::string m =" not useful for property ";
        if(has) m+="not";
        return m+" in category "+Property::categoryToString(cat);
      }
      Property::Category cat;
      bool has;
    };

    struct UsesEquality : OptionProblemConstraint{
      bool check(Property*p){
        ASS(p)
        return (p->equalityAtoms() != 0) ||
          // theories may introduce equality at various places of the pipeline!
          HasTheories::actualCheck(p);
      }
      std::string msg(){ return " only useful with equality"; }
    };

    struct HasHigherOrder : OptionProblemConstraint{
      bool check(Property*p){
        ASS(p)
        return (p->higherOrder());
      }
      std::string msg(){ return " only useful with higher-order problems"; }
    };

    struct OnlyFirstOrder : OptionProblemConstraint{
      bool check(Property*p){
        ASS(p)
        return (!p->higherOrder());
      }
      std::string msg(){ return " not compatible with higher-order problems"; }
    };

    struct MayHaveNonUnits : OptionProblemConstraint{
      bool check(Property*p){
        return (p->formulas() > 0) // let's not try to guess what kind of clauses these will give rise to
          || (p->clauses() > p->unitClauses());
      }
      std::string msg(){ return " only useful with non-unit clauses"; }
    };

    struct NotJustEquality : OptionProblemConstraint{
      bool check(Property*p){
        return (p->category()!=Property::PEQ || p->category()!=Property::UEQ);
      }
      std::string msg(){ return " not useful with just equality"; }
    };

    struct AtomConstraint : OptionProblemConstraint{
      AtomConstraint(int a,bool g) : atoms(a),greater(g) {}
      int atoms;
      bool greater;
      bool check(Property*p){
        return greater ? p->atoms()>atoms : p->atoms()<atoms;
      }
          
      std::string msg(){ 
        std::string m = " not with ";
        if(greater){ m+="more";}else{m+="less";}
        return m+" than "+Lib::Int::toString(atoms)+" atoms";
      }
    };

    struct HasTheories : OptionProblemConstraint {
      static bool actualCheck(Property*p);

      bool check(Property*p);
      std::string msg(){ return " only useful with theories"; }
    };

    struct HasFormulas : OptionProblemConstraint {
      bool check(Property*p) {
        return p->hasFormulas();
      }
      std::string msg(){ return " only useful with (non-cnf) formulas"; }
    };

    struct HasGoal : OptionProblemConstraint {
      bool check(Property*p){
        return p->hasGoal();
      }
      std::string msg(){ return " only useful with a goal: (conjecture) formulas or (negated_conjecture) clauses"; }
    };

    // Factory methods
    static OptionProblemConstraintUP notWithCat(Property::Category c){
      return OptionProblemConstraintUP(new CategoryCondition(c,false));
    }
    static OptionProblemConstraintUP hasCat(Property::Category c){
      return OptionProblemConstraintUP(new CategoryCondition(c,true));
    }
    static OptionProblemConstraintUP hasEquality(){ return OptionProblemConstraintUP(new UsesEquality); }
    static OptionProblemConstraintUP hasHigherOrder(){ return OptionProblemConstraintUP(new HasHigherOrder); }
    static OptionProblemConstraintUP onlyFirstOrder(){ return OptionProblemConstraintUP(new OnlyFirstOrder); }
    static OptionProblemConstraintUP mayHaveNonUnits(){ return OptionProblemConstraintUP(new MayHaveNonUnits); }
    static OptionProblemConstraintUP notJustEquality(){ return OptionProblemConstraintUP(new NotJustEquality); }
    static OptionProblemConstraintUP atomsMoreThan(int a){
      return OptionProblemConstraintUP(new AtomConstraint(a,true));
    }
    static OptionProblemConstraintUP atomsLessThan(int a){
      return OptionProblemConstraintUP(new AtomConstraint(a,false));
    }
    static OptionProblemConstraintUP hasFormulas() { return OptionProblemConstraintUP(new HasFormulas); }
    static OptionProblemConstraintUP hasTheories() { return OptionProblemConstraintUP(new HasTheories); }
    static OptionProblemConstraintUP hasGoal() { return OptionProblemConstraintUP(new HasGoal); }

    //Cheating - we refer to env.options to ask about option values
    // There is an assumption that the option values used have been
    // set to their final values
    // These are used in randomisation where we guarantee a certain
    // set of options will not be randomized and some will be randomized first

    struct OptionHasValue : OptionProblemConstraint{
      OptionHasValue(std::string ov,std::string v) : option_value(ov),value(v) {}
      bool check(Property*p);
      std::string msg(){ return option_value+" has value "+value; } 
      std::string option_value;
      std::string value; 
    };

    struct ManyOptionProblemConstraints : OptionProblemConstraint {
      ManyOptionProblemConstraints(bool a) : is_and(a) {}

      bool check(Property*p){
        bool res = is_and;
        Stack<OptionProblemConstraintUP>::RefIterator it(cons);
        while(it.hasNext()){
          bool n=it.next()->check(p);res = is_and ? (res && n) : (res || n);}
        return res;
      }

      std::string msg(){
        std::string res="";
        Stack<OptionProblemConstraintUP>::RefIterator it(cons);
        if(it.hasNext()){ res=it.next()->msg();}
        while(it.hasNext()){ res+=",and\n"+it.next()->msg();}
        return res;
      }

      void add(OptionProblemConstraintUP& c){ cons.push(std::move(c));}
      Stack<OptionProblemConstraintUP> cons;
      bool is_and;
    };

    static OptionProblemConstraintUP And(OptionProblemConstraintUP left,
                                        OptionProblemConstraintUP right){
       ManyOptionProblemConstraints* c = new ManyOptionProblemConstraints(true);
       c->add(left);c->add(right);
       return OptionProblemConstraintUP(c);
    }
    static OptionProblemConstraintUP And(OptionProblemConstraintUP left,
                                        OptionProblemConstraintUP mid,
                                        OptionProblemConstraintUP right){
       ManyOptionProblemConstraints* c = new ManyOptionProblemConstraints(true);
       c->add(left);c->add(mid);c->add(right);
       return OptionProblemConstraintUP(c);
    }
    static OptionProblemConstraintUP Or(OptionProblemConstraintUP left,
                                        OptionProblemConstraintUP right){
       ManyOptionProblemConstraints* c = new ManyOptionProblemConstraints(false);
       c->add(left);c->add(right);
       return OptionProblemConstraintUP(c);
    }
    static OptionProblemConstraintUP Or(OptionProblemConstraintUP left,
                                        OptionProblemConstraintUP mid,
                                        OptionProblemConstraintUP right){
       ManyOptionProblemConstraints* c = new ManyOptionProblemConstraints(false);
       c->add(left);c->add(mid);c->add(right);
       return OptionProblemConstraintUP(c);
    }

  //==========================================================
  // Getter functions
  // -currently disabled all unnecessary setter functions
  //==========================================================
  //
  // This is how options are accessed so if you add a new option you should add a getter
public:
  bool encodeStrategy() const{ return _encode.actualValue;}
  BadOption getBadOptionChoice() const { return _badOption.actualValue; }
  void setBadOptionChoice(BadOption newVal) { _badOption.actualValue = newVal; }
  std::string forcedOptions() const { return _forcedOptions.actualValue; }
  std::string forbiddenOptions() const { return _forbiddenOptions.actualValue; }
  std::string testId() const { return _testId.actualValue; }
  std::string protectedPrefix() const { return _protectedPrefix.actualValue; }
  Statistics statistics() const { return _statistics.actualValue; }
  void setStatistics(Statistics newVal) { _statistics.actualValue=newVal; }
  Proof proof() const { return _proof.actualValue; }
  bool minimizeSatProofs() const { return _minimizeSatProofs.actualValue; }
  ProofExtra proofExtra() const { return _proofExtra.actualValue; }
  bool traceback() const { return _traceback.actualValue; }
  void setTraceback(bool traceback) { _traceback.actualValue = traceback; }
  std::string printProofToFile() const { return _printProofToFile.actualValue; }
  int naming() const { return _naming.actualValue; }

  bool fmbNonGroundDefs() const { return _fmbNonGroundDefs.actualValue; }
  unsigned fmbStartSize() const { return _fmbStartSize.actualValue;}
  float fmbSymmetryRatio() const { return _fmbSymmetryRatio.actualValue; }
  FMBWidgetOrders fmbSymmetryWidgetOrders() { return _fmbSymmetryWidgetOrders.actualValue;}
  FMBSymbolOrders fmbSymmetryOrderSymbols() const {return _fmbSymmetryOrderSymbols.actualValue; }
  FMBAdjustSorts fmbAdjustSorts() const {return _fmbAdjustSorts.actualValue; }
  bool fmbDetectSortBounds() const { return _fmbDetectSortBounds.actualValue; }
  unsigned fmbDetectSortBoundsTimeLimit() const { return _fmbDetectSortBoundsTimeLimit.actualValue; }
  unsigned fmbSizeWeightRatio() const { return _fmbSizeWeightRatio.actualValue; }
  FMBEnumerationStrategy fmbEnumerationStrategy() const { return _fmbEnumerationStrategy.actualValue; }
  bool keepSbeamGenerators() const { return _fmbKeepSbeamGenerators.actualValue; }

  bool flattenTopLevelConjunctions() const { return _flattenTopLevelConjunctions.actualValue; }
  Mode mode() const { return _mode.actualValue; }
  void setMode(Mode mode) { _mode.actualValue = mode; }
  Intent intent() const { return _intent.actualValue; }
  Schedule schedule() const { return _schedule.actualValue; }
  std::string scheduleName() const { return _schedule.getStringOfValue(_schedule.actualValue); }
  void setSchedule(Schedule newVal) {  _schedule.actualValue = newVal; }
  std::string scheduleFile() const { return _scheduleFile.actualValue; }
  unsigned multicore() const { return _multicore.actualValue; }
  void setMulticore(unsigned newVal) { _multicore.actualValue = newVal; }
  float slowness() const {return _slowness.actualValue; }
  InputSyntax inputSyntax() const { return _inputSyntax.actualValue; }
  void setInputSyntax(InputSyntax newVal) { _inputSyntax.actualValue = newVal; }
  bool normalize() const { return _normalize.actualValue; }
  void setNormalize(bool normalize) { _normalize.actualValue = normalize; }
  GoalGuess guessTheGoal() const { return _guessTheGoal.actualValue; }
  unsigned gtgLimit() const { return _guessTheGoalLimit.actualValue; }
  void setMaxXX(unsigned max) { _maximumXXNarrows.actualValue = max; }

  void setNaming(int n){ _naming.actualValue = n;} //TODO: ensure global constraints
  std::string include() const { return _include.actualValue; }
  void setInclude(std::string val) { _include.actualValue = val; }
  std::string inputFile() const { return _inputFile.actualValue; }
  void resetInputFile() { _inputFile.actualValue = ""; }
  int activationLimit() const { return _activationLimit.actualValue; }
  unsigned randomSeed() const { return _randomSeed.actualValue; }
  void setRandomSeed(unsigned seed) { _randomSeed.actualValue = seed; }
  const std::string& strategySamplerFilename() const { return _sampleStrategy.actualValue; }
  bool printClausifierPremises() const { return _printClausifierPremises.actualValue; }

  // IMPORTANT, if you add a showX command then include showAll
  bool showAll() const { return _showAll.actualValue; }
  bool showActive() const { return showAll() || _showActive.actualValue; }
  bool showBlocked() const { return showAll() || _showBlocked.actualValue; }
  bool showDefinitions() const { return showAll() || _showDefinitions.actualValue; }
  bool showNew() const { return showAll() || _showNew.actualValue; }
  bool sineToAge() const { return _sineToAge.actualValue; }
  PredicateSineLevels sineToPredLevels() const { return _sineToPredLevels.actualValue; }
  bool showSplitting() const { return showAll() || _showSplitting.actualValue; }
  bool showNewPropositional() const { return showAll() || _showNewPropositional.actualValue; }
  bool showPassive() const { return showAll() || _showPassive.actualValue; }
  bool showReductions() const { return showAll() || _showReductions.actualValue; }
  bool showPreprocessing() const { return showAll() || _showPreprocessing.actualValue; }
  bool showSkolemisations() const { return showAll() || _showSkolemisations.actualValue; }
  bool showSymbolElimination() const { return showAll() || _showSymbolElimination.actualValue; }
  bool showTheoryAxioms() const { return showAll() || _showTheoryAxioms.actualValue; }
  bool showFOOL() const { return showAll() || _showFOOL.actualValue; }
  bool showFMBsortInfo() const { return showAll() || _showFMBsortInfo.actualValue; }
  bool showInduction() const { return showAll() || _showInduction.actualValue; }
  bool showSimplOrdering() const { return showAll() || _showSimplOrdering.actualValue; }
#if VAMPIRE_CLAUSE_TRACING
  int traceBackward() { return _traceBackward.actualValue; }
  int traceForward() { return _traceForward.actualValue; }
#endif // VAMPIRE_CLAUSE_TRACING

#if VZ3
  bool showZ3() const { return showAll() || _showZ3.actualValue; }
  ProblemExportSyntax problemExportSyntax() const { return _problemExportSyntax.actualValue; }
  std::string const& exportAvatarProblem() const { return _exportAvatarProblem.actualValue; }
  std::string const& exportThiProblem() const { return _exportThiProblem.actualValue; }
#endif

  // end of show commands

  bool showNonconstantSkolemFunctionTrace() const { return _showNonconstantSkolemFunctionTrace.actualValue; }
  void setShowNonconstantSkolemFunctionTrace(bool newVal) { _showNonconstantSkolemFunctionTrace.actualValue = newVal; }
  InterpolantMode showInterpolant() const { return _showInterpolant.actualValue; }
  bool showOptions() const { return _showOptions.actualValue; }
  bool lineWrapInShowOptions() const { return _showOptionsLineWrap.actualValue; }
  bool showExperimentalOptions() const { return _showExperimentalOptions.actualValue; }
  bool showHelp() const { return _showHelp.actualValue; }
  std::string explainOption() const { return _explainOption.actualValue; }

  bool printAllTheoryAxioms() const { return _printAllTheoryAxioms.actualValue; }

#if VZ3
  bool satFallbackForSMT() const { return _satFallbackForSMT.actualValue; }
  bool smtForGround() const { return _smtForGround.actualValue; }
  TheoryInstSimp theoryInstAndSimp() const { return _theoryInstAndSimp.actualValue; }
  bool thiGeneralise() const { return _thiGeneralise.actualValue; }
  bool thiTautologyDeletion() const { return _thiTautologyDeletion.actualValue; }
#endif
  UnificationWithAbstraction unificationWithAbstraction() const { return _unificationWithAbstraction.actualValue; }
  bool unificationWithAbstractionFixedPointIteration() const { return _unificationWithAbstractionFixedPointIteration.actualValue; }
  void setUWA(UnificationWithAbstraction value){ _unificationWithAbstraction.actualValue = value; }
  bool fixUWA() const { return _fixUWA.actualValue; }
  bool useACeval() const { return _useACeval.actualValue;}

  bool unusedPredicateDefinitionRemoval() const { return _unusedPredicateDefinitionRemoval.actualValue; }
  bool blockedClauseElimination() const { return _blockedClauseElimination.actualValue; }
  unsigned distinctGroupExpansionLimit() const { return _distinctGroupExpansionLimit.actualValue; }
  void setUnusedPredicateDefinitionRemoval(bool newVal) { _unusedPredicateDefinitionRemoval.actualValue = newVal; }
  SatSolver satSolver() const { return _satSolver.actualValue; }
  //void setSatSolver(SatSolver newVal) { _satSolver = newVal; }
  SaturationAlgorithm saturationAlgorithm() const { return _saturationAlgorithm.actualValue; }
  void setSaturationAlgorithm(SaturationAlgorithm newVal) { _saturationAlgorithm.actualValue = newVal; }
  int selection() const { return _selection.actualValue; }
  void setSelection(int v) { _selection.actualValue=v;}
  std::string latexOutput() const { return _latexOutput.actualValue; }
  bool latexUseDefault() const { return _latexUseDefaultSymbols.actualValue; }
  LiteralComparisonMode literalComparisonMode() const { return _literalComparisonMode.actualValue; }
  bool forwardSubsumptionResolution() const { return _forwardSubsumptionResolution.actualValue; }
  //void setForwardSubsumptionResolution(bool newVal) { _forwardSubsumptionResolution = newVal; }
  bool forwardSubsumptionDemodulation() const { return _forwardSubsumptionDemodulation.actualValue; }
  unsigned forwardSubsumptionDemodulationMaxMatches() const { return _forwardSubsumptionDemodulationMaxMatches.actualValue; }
  Demodulation forwardDemodulation() const { return _forwardDemodulation.actualValue; }
  bool binaryResolution() const { return _binaryResolution.actualValue; }
  bool superposition() const {return _superposition.actualValue; }
  URResolution unitResultingResolution() const { return _unitResultingResolution.actualValue; }
  bool simulatenousSuperposition() const { return _simultaneousSuperposition.actualValue; }
  bool innerRewriting() const { return _innerRewriting.actualValue; }
  bool equationalTautologyRemoval() const { return _equationalTautologyRemoval.actualValue; }
  bool conditionalRedundancyCheck() const { return _conditionalRedundancyCheck.actualValue; }
  bool conditionalRedundancyOrderingConstraints() const { return _conditionalRedundancyOrderingConstraints.actualValue; }
  bool conditionalRedundancyAvatarConstraints() const { return _conditionalRedundancyAvatarConstraints.actualValue; }
  bool conditionalRedundancyLiteralConstraints() const { return _conditionalRedundancyLiteralConstraints.actualValue; }
  bool arityCheck() const { return _arityCheck.actualValue; }
  //void setArityCheck(bool newVal) { _arityCheck=newVal; }
  Demodulation backwardDemodulation() const { return _backwardDemodulation.actualValue; }
  DemodulationRedundancyCheck demodulationRedundancyCheck() const { return _demodulationRedundancyCheck.actualValue; }
  bool demodulationPrecompiledComparison() const { return _demodulationPrecompiledComparison.actualValue; }
  bool demodulationOnlyEquational() const { return _demodulationOnlyEquational.actualValue; }

  //void setBackwardDemodulation(Demodulation newVal) { _backwardDemodulation = newVal; }
  Subsumption backwardSubsumption() const { return _backwardSubsumption.actualValue; }
  //void setBackwardSubsumption(Subsumption newVal) { _backwardSubsumption = newVal; }
  Subsumption backwardSubsumptionResolution() const { return _backwardSubsumptionResolution.actualValue; }
  bool backwardSubsumptionDemodulation() const { return _backwardSubsumptionDemodulation.actualValue; }
  unsigned backwardSubsumptionDemodulationMaxMatches() const { return _backwardSubsumptionDemodulationMaxMatches.actualValue; }
  bool forwardSubsumption() const { return _forwardSubsumption.actualValue; }
  bool forwardLiteralRewriting() const { return _forwardLiteralRewriting.actualValue; }
  int lrsFirstTimeCheck() const { return _lrsFirstTimeCheck.actualValue; }
  int lrsWeightLimitOnly() const { return _lrsWeightLimitOnly.actualValue; }
  int lrsRetroactiveDeletes() const { return _lrsRetroactiveDeletes.actualValue; }
  int lookaheadDelay() const { return _lookaheadDelay.actualValue; }
  int simulatedTimeLimit() const { return _simulatedTimeLimit.actualValue; }
  void setSimulatedTimeLimit(int newVal) { _simulatedTimeLimit.actualValue = newVal; }
  float lrsEstimateCorrectionCoef() const { return _lrsEstimateCorrectionCoef.actualValue; }
  TermOrdering termOrdering() const { return _termOrdering.actualValue; }
  SymbolPrecedence symbolPrecedence() const { return _symbolPrecedence.actualValue; }
  SymbolPrecedenceBoost symbolPrecedenceBoost() const { return _symbolPrecedenceBoost.actualValue; }
  IntroducedSymbolPrecedence introducedSymbolPrecedence() const { return _introducedSymbolPrecedence.actualValue; }
  KboWeightGenerationScheme kboWeightGenerationScheme() const { return _kboWeightGenerationScheme.actualValue; }
  bool kboMaxZero() const { return _kboMaxZero.actualValue; }
  const KboAdmissibilityCheck kboAdmissabilityCheck() const { return _kboAdmissabilityCheck.actualValue; }
  const std::string& functionWeights() const { return _functionWeights.actualValue; }
  const std::string& predicateWeights() const { return _predicateWeights.actualValue; }
  const std::string& functionPrecedence() const { return _functionPrecedence.actualValue; }
  const std::string& typeConPrecedence() const { return _typeConPrecedence.actualValue; }
  const std::string& predicatePrecedence() const { return _predicatePrecedence.actualValue; }
  // Return time limit in deciseconds, or 0 if there is no time limit
  int timeLimitInDeciseconds() const { return _timeLimitInDeciseconds.actualValue; }
  size_t memoryLimit() const { return _memoryLimit.actualValue; }
  void setMemoryLimitOptionValue(size_t newVal) { _memoryLimit.actualValue = newVal; }
#if VAMPIRE_PERF_EXISTS
  unsigned instructionLimit() const { return _instructionLimit.actualValue; }
  void setInstructionLimit(unsigned newVal) { _instructionLimit.actualValue = newVal; }
  unsigned simulatedInstructionLimit() const { return _simulatedInstructionLimit.actualValue; }
  unsigned setSimulatedInstructionLimit() const { return _simulatedInstructionLimit.actualValue; }
  bool parsingDoesNotCount() const { return _parsingDoesNotCount.actualValue; }
#endif
  bool interactive() const { return _interactive.actualValue; }
  void setInteractive(bool v) { _interactive.actualValue = v; }
  int inequalitySplitting() const { return _inequalitySplitting.actualValue; }
  int ageRatio() const { return _ageWeightRatio.actualValue; }
  void setAgeRatio(int v){ _ageWeightRatio.actualValue = v; }
  int weightRatio() const { return _ageWeightRatio.otherValue; }
  bool useTheorySplitQueues() const { return _useTheorySplitQueues.actualValue; }
  std::vector<int> theorySplitQueueRatios() const;
  std::vector<float> theorySplitQueueCutoffs() const;
  int theorySplitQueueExpectedRatioDenom() const { return _theorySplitQueueExpectedRatioDenom.actualValue; }
  bool theorySplitQueueLayeredArrangement() const { return _theorySplitQueueLayeredArrangement.actualValue; }
  bool useAvatarSplitQueues() const { return _useAvatarSplitQueues.actualValue; }
  std::vector<int> avatarSplitQueueRatios() const;
  std::vector<float> avatarSplitQueueCutoffs() const;
  bool avatarSplitQueueLayeredArrangement() const { return _avatarSplitQueueLayeredArrangement.actualValue; }
  bool useSineLevelSplitQueues() const { return _useSineLevelSplitQueues.actualValue; }
  std::vector<int> sineLevelSplitQueueRatios() const;
  std::vector<float> sineLevelSplitQueueCutoffs() const;
  bool sineLevelSplitQueueLayeredArrangement() const { return _sineLevelSplitQueueLayeredArrangement.actualValue; }
  bool usePositiveLiteralSplitQueues() const { return _usePositiveLiteralSplitQueues.actualValue; }
  std::vector<int> positiveLiteralSplitQueueRatios() const;
  std::vector<float> positiveLiteralSplitQueueCutoffs() const;
  bool positiveLiteralSplitQueueLayeredArrangement() const { return _positiveLiteralSplitQueueLayeredArrangement.actualValue; }
  void setWeightRatio(int v){ _ageWeightRatio.otherValue = v; }
  bool literalMaximalityAftercheck() const { return _literalMaximalityAftercheck.actualValue; }
  bool superpositionFromVariables() const { return _superpositionFromVariables.actualValue; }
  EqualityProxy equalityProxy() const { return _equalityProxy.actualValue; }
  bool useMonoEqualityProxy() const { return _useMonoEqualityProxy.actualValue; }
  bool equalityResolutionWithDeletion() const { return _equalityResolutionWithDeletion.actualValue; }
  ExtensionalityResolution extensionalityResolution() const { return _extensionalityResolution.actualValue; }
  bool FOOLParamodulation() const { return _FOOLParamodulation.actualValue; }
  bool termAlgebraInferences() const { return _termAlgebraInferences.actualValue; }
  bool termAlgebraExhaustivenessAxiom() const { return _termAlgebraExhaustivenessAxiom.actualValue; }
  TACyclicityCheck termAlgebraCyclicityCheck() const { return _termAlgebraCyclicityCheck.actualValue; }
  unsigned extensionalityMaxLength() const { return _extensionalityMaxLength.actualValue; }
  bool extensionalityAllowPosEq() const { return _extensionalityAllowPosEq.actualValue; }
  unsigned nongoalWeightCoefficientNumerator() const { return _nonGoalWeightCoefficient.numerator; }
  unsigned nongoalWeightCoefficientDenominator() const { return _nonGoalWeightCoefficient.denominator; }
  bool restrictNWCtoGC() const { return _restrictNWCtoGC.actualValue; }
  Sos sos() const { return _sos.actualValue; }
  unsigned sosTheoryLimit() const { return _sosTheoryLimit.actualValue; }
  //void setSos(Sos newVal) { _sos = newVal; }

  bool shuffleInput() const { return _shuffleInput.actualValue; }
  bool randomPolarities() const { return _randomPolarities.actualValue; }
  bool randomAWR() const { return _randomAWR.actualValue; }
  bool randomTraversals() const { return _randomTraversals.actualValue; }
  bool randomizeSeedForPortfolioWorkers() const { return _randomizSeedForPortfolioWorkers.actualValue; }
  void setRandomizeSeedForPortfolioWorkers(bool val) { _randomizSeedForPortfolioWorkers.actualValue = val; }

  bool ignoreConjectureInPreprocessing() const {return _ignoreConjectureInPreprocessing.actualValue;}

  FunctionDefinitionElimination functionDefinitionElimination() const { return _functionDefinitionElimination.actualValue; }
  unsigned functionDefinitionIntroduction() const { return _functionDefinitionIntroduction.actualValue; }
  TweeGoalTransformation tweeGoalTransformation() const { return _tweeGoalTransformation.actualValue; }
  bool codeTreeSubsumption() const { return _codeTreeSubsumption.actualValue; }
  bool outputAxiomNames() const { return _outputAxiomNames.actualValue; }
  void setOutputAxiomNames(bool newVal) { _outputAxiomNames.actualValue = newVal; }
  QuestionAnsweringMode questionAnswering() const { return _questionAnswering.actualValue; }
  bool questionAnsweringGroundOnly() const { return _questionAnsweringGroundOnly.actualValue; }
  std::string questionAnsweringAvoidThese() const { return _questionAnsweringAvoidThese.actualValue; }
  Output outputMode() const { return _outputMode.actualValue; }
  void setOutputMode(Output newVal) { _outputMode.actualValue = newVal; }
  bool ignoreMissingInputsInUnsatCore() {  return _ignoreMissingInputsInUnsatCore.actualValue; }
  std::string thanks() const { return _thanks.actualValue; }
  void setQuestionAnswering(QuestionAnsweringMode newVal) { _questionAnswering.actualValue = newVal; }
  bool globalSubsumption() const { return _globalSubsumption.actualValue; }
  GlobalSubsumptionSatSolverPower globalSubsumptionSatSolverPower() const { return _globalSubsumptionSatSolverPower.actualValue; }
  GlobalSubsumptionExplicitMinim globalSubsumptionExplicitMinim() const { return _globalSubsumptionExplicitMinim.actualValue; }
  GlobalSubsumptionAvatarAssumptions globalSubsumptionAvatarAssumptions() const { return _globalSubsumptionAvatarAssumptions.actualValue; }

  /** true if calling set() on non-existing options does not result in a user error */
  IgnoreMissing ignoreMissing() const { return _ignoreMissing.actualValue; }
  void setIgnoreMissing(IgnoreMissing newVal) { _ignoreMissing.actualValue = newVal; }
  bool increasedNumeralWeight() const { return _increasedNumeralWeight.actualValue; }
  TheoryAxiomLevel theoryAxioms() const { return _theoryAxioms.actualValue; }
  //void setTheoryAxioms(bool newValue) { _theoryAxioms = newValue; }
  Condensation condensation() const { return _condensation.actualValue; }
  bool generalSplitting() const { return _generalSplitting.actualValue; }
#if VTIME_PROFILING
  bool timeStatistics() const { return _timeStatistics.actualValue; }
  std::string const& timeStatisticsFocus() const { return _timeStatisticsFocus.actualValue; }
#endif // VTIME_PROFILING
  bool splitting() const { return _splitting.actualValue; }
  void setSplitting(bool value){ _splitting.actualValue=value; }
  bool nonliteralsInClauseWeight() const { return _nonliteralsInClauseWeight.actualValue; }
  unsigned sineDepth() const { return _sineDepth.actualValue; }
  unsigned sineGeneralityThreshold() const { return _sineGeneralityThreshold.actualValue; }
  unsigned sineToAgeGeneralityThreshold() const { return _sineToAgeGeneralityThreshold.actualValue; }
  SineSelection sineSelection() const { return _sineSelection.actualValue; }
  void setSineSelection(SineSelection val) { _sineSelection.actualValue=val; }
  float sineTolerance() const { return _sineTolerance.actualValue; }
  float sineToAgeTolerance() const { return _sineToAgeTolerance.actualValue; }

  bool colorUnblocking() const { return _colorUnblocking.actualValue; }

  Instantiation instantiation() const { return _instantiation.actualValue; }
  bool theoryFlattening() const { return _theoryFlattening.actualValue; }
  bool ignoreUnrecognizedLogic() const { return _ignoreUnrecognizedLogic.actualValue; }

  Induction induction() const { return _induction.actualValue; }
  StructuralInductionKind structInduction() const { return _structInduction.actualValue; }
  IntInductionKind intInduction() const { return _intInduction.actualValue; }
  InductionChoice inductionChoice() const { return _inductionChoice.actualValue; }
  unsigned maxInductionDepth() const { return _maxInductionDepth.actualValue; }
  bool inductionNegOnly() const { return _inductionNegOnly.actualValue; }
  bool inductionUnitOnly() const { return _inductionUnitOnly.actualValue; }
  bool inductionGen() const { return _inductionGen.actualValue; }
  bool inductionGenHeur() const { return _inductionGenHeur.actualValue; }
  bool inductionStrengthenHypothesis() const { return _inductionStrengthenHypothesis.actualValue; }
  unsigned maxInductionGenSubsetSize() const { return _maxInductionGenSubsetSize.actualValue; }
  bool inductionOnComplexTerms() const {return _inductionOnComplexTerms.actualValue;}
  bool functionDefinitionRewriting() const { return _functionDefinitionRewriting.actualValue; }
  bool integerInductionDefaultBound() const { return _integerInductionDefaultBound.actualValue; }
  IntegerInductionInterval integerInductionInterval() const { return _integerInductionInterval.actualValue; }
  IntegerInductionLiteralStrictness integerInductionStrictnessEq() const {return _integerInductionStrictnessEq.actualValue; }
  IntegerInductionLiteralStrictness integerInductionStrictnessComp() const {return _integerInductionStrictnessComp.actualValue; }
  IntegerInductionTermStrictness integerInductionStrictnessTerm() const {return _integerInductionStrictnessTerm.actualValue; }
  bool nonUnitInduction() const { return _nonUnitInduction.actualValue; }
  bool inductionOnActiveOccurrences() const { return _inductionOnActiveOccurrences.actualValue; }

  bool useHashingVariantIndex() const { return _useHashingVariantIndex.actualValue; }

  void setTimeLimitInSeconds(int newVal) { _timeLimitInDeciseconds.actualValue = 10*newVal; }
  void setTimeLimitInDeciseconds(int newVal) { _timeLimitInDeciseconds.actualValue = newVal; }

  bool splitAtActivation() const{ return _splitAtActivation.actualValue; }
  SplittingNonsplittableComponents splittingNonsplittableComponents() const { return _splittingNonsplittableComponents.actualValue; }
  SplittingAddComplementary splittingAddComplementary() const { return _splittingAddComplementary.actualValue; }
  SplittingMinimizeModel splittingMinimizeModel() const { return _splittingMinimizeModel.actualValue; }
  SplittingLiteralPolarityAdvice splittingLiteralPolarityAdvice() const { return _splittingLiteralPolarityAdvice.actualValue; }
  SplittingDeleteDeactivated splittingDeleteDeactivated() const { return _splittingDeleteDeactivated.actualValue;}
  bool splittingFastRestart() const { return _splittingFastRestart.actualValue; }
  bool splittingBufferedSolver() const { return _splittingBufferedSolver.actualValue; }
  int splittingFlushPeriod() const { return _splittingFlushPeriod.actualValue; }
  float splittingFlushQuotient() const { return _splittingFlushQuotient.actualValue; }
  float splittingAvatimer() const { return _splittingAvatimer.actualValue; }
  bool splittingEagerRemoval() const { return _splittingEagerRemoval.actualValue; }
  SplittingCongruenceClosure splittingCongruenceClosure() const { return _splittingCongruenceClosure.actualValue; }
  CCUnsatCores ccUnsatCores() const { return _ccUnsatCores.actualValue; }

  void setProof(Proof p) { _proof.actualValue = p; }
  bool newCNF() const { return _newCNF.actualValue; }
  bool getIteInlineLet() const { return _inlineLet.actualValue; }

  bool useManualClauseSelection() const { return _manualClauseSelection.actualValue; }
  bool inequalityNormalization() const { return _inequalityNormalization.actualValue; }
  EvaluationMode evaluationMode() const { return _highSchool.actualValue ? EvaluationMode::POLYNOMIAL_CAUTIOUS : _evaluationMode.actualValue; }
  ArithmeticSimplificationMode gaussianVariableElimination() const { return _highSchool.actualValue ? ArithmeticSimplificationMode::CAUTIOUS : _gaussianVariableElimination.actualValue; }
  bool pushUnaryMinus() const { return _pushUnaryMinus.actualValue || _highSchool.actualValue; }
  ArithmeticSimplificationMode cancellation() const { return _highSchool.actualValue ? ArithmeticSimplificationMode::CAUTIOUS : _cancellation.actualValue; }
  ArithmeticSimplificationMode arithmeticSubtermGeneralizations() const { return  _highSchool.actualValue ? ArithmeticSimplificationMode::CAUTIOUS : _arithmeticSubtermGeneralizations.actualValue; }

  //Higher-order Options

  bool addCombAxioms() const { return _addCombAxioms.actualValue; }
  bool addProxyAxioms() const { return _addProxyAxioms.actualValue; }
  bool combinatorySup() const { return _combinatorySuperposition.actualValue; }
  bool choiceAxiom() const { return _choiceAxiom.actualValue; }
  bool injectivityReasoning() const { return _injectivity.actualValue; }
  bool pragmatic() const { return _pragmatic.actualValue; }
  bool choiceReasoning() const { return _choiceReasoning.actualValue; }
  bool prioritiseClausesProducedByLongReduction() const { return _priortyToLongReducts.actualValue; }
  int maxXXNarrows() const { return _maximumXXNarrows.actualValue; }
  FunctionExtensionality functionExtensionality() const { return _functionExtensionality.actualValue; }
  CNFOnTheFly cnfOnTheFly() const { return _clausificationOnTheFly.actualValue; }
  PISet piSet() const { return _piSet.actualValue; }
  Narrow narrow() const { return _narrow.actualValue; }
  bool equalityToEquivalence () const { return _equalityToEquivalence.actualValue; }
  bool complexBooleanReasoning () const { return _complexBooleanReasoning.actualValue; }
  bool booleanEqTrick() const { return _booleanEqTrick.actualValue; }
  bool casesSimp() const { return _casesSimp.actualValue; }
  bool cases() const { return _cases.actualValue; }
  bool newTautologyDel() const { return _newTautologyDel.actualValue; }
  bool lambdaFreeHol() const { return _lambdaFreeHol.actualValue; }
  bool complexVarCondition() const { return _complexVarCondition.actualValue; }
  // For unit testing
  void useCombSup() {
    _combinatorySuperposition.actualValue = true;
    _complexVarCondition.actualValue = true;
  }

private:

    /**
     * A LookupWrapper is used to wrap up two maps for long and short names and query them
     */
    struct LookupWrapper {

        LookupWrapper() {}

        private:
          LookupWrapper operator=(const LookupWrapper&){ NOT_IMPLEMENTED;}
        public:

        void insert(AbstractOptionValue* option_value){
            ASS(!option_value->longName.empty());
            bool new_long =  _longMap.insert(option_value->longName,option_value);
            bool new_short = true;
            if(!option_value->shortName.empty()){
                new_short = _shortMap.insert(option_value->shortName,option_value);
            }
            if(!new_long || !new_short){ std::cout << "Bad " << option_value->longName << std::endl; }
            ASS(new_long && new_short);
        }
        AbstractOptionValue* findLong(std::string longName) const{
            if(!_longMap.find(longName)){ throw ValueNotFoundException(); }
            return _longMap.get(longName);
        }
        AbstractOptionValue* findShort(std::string shortName) const{
            if(!_shortMap.find(shortName)){ throw ValueNotFoundException(); }
            return _shortMap.get(shortName);
        }

        VirtualIterator<AbstractOptionValue*> values() const {
            return _longMap.range();
        }

    private:
        DHMap<std::string,AbstractOptionValue*> _longMap;
        DHMap<std::string,AbstractOptionValue*> _shortMap;
    };

    LookupWrapper _lookup;

    // The const is a lie - we can alter the resulting OptionValue
    AbstractOptionValue* getOptionValueByName(std::string name) const{
        try{
          return _lookup.findLong(name);
        }
        catch(ValueNotFoundException&){
          try{
            return _lookup.findShort(name);
          }
          catch(ValueNotFoundException&){
            return 0;
          }
        }
    }
  
    Stack<std::string> getSimilarOptionNames(std::string name, bool is_short) const;

    //==========================================================
    // Variables holding option values
    //==========================================================

 /**
  * NOTE on OptionValues
  *
  * An OptionValue stores the value for an Option as well as all the meta-data
  * See the definitions of different OptionValue objects above for details
  * but the main OptionValuse are
  *  - BoolOptionValue
  *  - IntOptionValue, UnsignedOptionValue, FloatOptionValue, LongOptionValue
  *  - StringOptionValue
  *  - ChoiceOptionValue
  *  - RatioOptionValue
  *
  * ChoiceOptionValue requires you to define an enum for the choice values
  *
  * For examples of how the different OptionValues are used see Options.cpp
  *
  * If an OptionValue needs custom assignment you will need to create a custom
  *  OptionValue. See DecodeOptionValue and SelectionOptionValue for examples.
  *
  */

  DecodeOptionValue _decode;
  BoolOptionValue _encode;

  RatioOptionValue _ageWeightRatio;

  BoolOptionValue _useTheorySplitQueues;
  StringOptionValue _theorySplitQueueRatios;
  StringOptionValue _theorySplitQueueCutoffs;
  IntOptionValue _theorySplitQueueExpectedRatioDenom;
  BoolOptionValue _theorySplitQueueLayeredArrangement;
  BoolOptionValue _useAvatarSplitQueues;
  StringOptionValue _avatarSplitQueueRatios;
  StringOptionValue _avatarSplitQueueCutoffs;
  BoolOptionValue _avatarSplitQueueLayeredArrangement;
  BoolOptionValue _useSineLevelSplitQueues;
  StringOptionValue _sineLevelSplitQueueRatios;
  StringOptionValue _sineLevelSplitQueueCutoffs;
  BoolOptionValue _sineLevelSplitQueueLayeredArrangement;
  BoolOptionValue _usePositiveLiteralSplitQueues;
  StringOptionValue _positiveLiteralSplitQueueRatios;
  StringOptionValue _positiveLiteralSplitQueueCutoffs;
  BoolOptionValue _positiveLiteralSplitQueueLayeredArrangement;
	BoolOptionValue _randomAWR;
  BoolOptionValue _literalMaximalityAftercheck;
  BoolOptionValue _arityCheck;

  BoolOptionValue _randomTraversals;

  ChoiceOptionValue<BadOption> _badOption;
  ChoiceOptionValue<Demodulation> _backwardDemodulation;
  ChoiceOptionValue<Subsumption> _backwardSubsumption;
  ChoiceOptionValue<Subsumption> _backwardSubsumptionResolution;
  BoolOptionValue _backwardSubsumptionDemodulation;
  UnsignedOptionValue _backwardSubsumptionDemodulationMaxMatches;
  BoolOptionValue _binaryResolution;

  BoolOptionValue _colorUnblocking;
  ChoiceOptionValue<Condensation> _condensation;

  ChoiceOptionValue<DemodulationRedundancyCheck> _demodulationRedundancyCheck;
  BoolOptionValue _demodulationPrecompiledComparison;
  BoolOptionValue _demodulationOnlyEquational;

  ChoiceOptionValue<EqualityProxy> _equalityProxy;
  BoolOptionValue _useMonoEqualityProxy;
  BoolOptionValue _equalityResolutionWithDeletion;
  BoolOptionValue _equivalentVariableRemoval;
  ChoiceOptionValue<ExtensionalityResolution> _extensionalityResolution;
  UnsignedOptionValue _extensionalityMaxLength;
  BoolOptionValue _extensionalityAllowPosEq;

  BoolOptionValue _FOOLParamodulation;

  BoolOptionValue _termAlgebraInferences;
  ChoiceOptionValue<TACyclicityCheck> _termAlgebraCyclicityCheck;
  BoolOptionValue _termAlgebraExhaustivenessAxiom;

  BoolOptionValue _fmbNonGroundDefs;
  UnsignedOptionValue _fmbStartSize;
  FloatOptionValue _fmbSymmetryRatio;
  ChoiceOptionValue<FMBWidgetOrders> _fmbSymmetryWidgetOrders;
  ChoiceOptionValue<FMBSymbolOrders> _fmbSymmetryOrderSymbols;
  ChoiceOptionValue<FMBAdjustSorts> _fmbAdjustSorts;
  BoolOptionValue _fmbDetectSortBounds;
  TimeLimitOptionValue _fmbDetectSortBoundsTimeLimit;
  UnsignedOptionValue _fmbSizeWeightRatio;
  ChoiceOptionValue<FMBEnumerationStrategy> _fmbEnumerationStrategy;
  BoolOptionValue _fmbKeepSbeamGenerators;

  BoolOptionValue _flattenTopLevelConjunctions;
  StringOptionValue _forbiddenOptions;
  BoolOptionValue _forceIncompleteness;
  StringOptionValue _forcedOptions;
  ChoiceOptionValue<Demodulation> _forwardDemodulation;
  BoolOptionValue _forwardLiteralRewriting;
  BoolOptionValue _forwardSubsumption;
  BoolOptionValue _forwardSubsumptionResolution;
  BoolOptionValue _forwardSubsumptionDemodulation;
  UnsignedOptionValue _forwardSubsumptionDemodulationMaxMatches;
  ChoiceOptionValue<FunctionDefinitionElimination> _functionDefinitionElimination;
  UnsignedOptionValue _functionDefinitionIntroduction;
  ChoiceOptionValue<TweeGoalTransformation> _tweeGoalTransformation;
  BoolOptionValue _codeTreeSubsumption;

  BoolOptionValue _generalSplitting;
  BoolOptionValue _globalSubsumption;
  ChoiceOptionValue<GlobalSubsumptionSatSolverPower> _globalSubsumptionSatSolverPower;
  ChoiceOptionValue<GlobalSubsumptionExplicitMinim> _globalSubsumptionExplicitMinim;
  ChoiceOptionValue<GlobalSubsumptionAvatarAssumptions> _globalSubsumptionAvatarAssumptions;
  ChoiceOptionValue<GoalGuess> _guessTheGoal;
  UnsignedOptionValue _guessTheGoalLimit;

  BoolOptionValue _simultaneousSuperposition;
  BoolOptionValue _innerRewriting;
  BoolOptionValue _equationalTautologyRemoval;
  BoolOptionValue _conditionalRedundancyCheck;
  BoolOptionValue _conditionalRedundancyOrderingConstraints;
  BoolOptionValue _conditionalRedundancyAvatarConstraints;
  BoolOptionValue _conditionalRedundancyLiteralConstraints;

  /** if true, then calling set() on non-existing options will not result in a user error */
  ChoiceOptionValue<IgnoreMissing> _ignoreMissing;
  StringOptionValue _include;
  /** if this option is true, Vampire will add the numeral weight of a clause
   * to its weight. The weight is defined as the sum of binary sizes of all
   * integers occurring in this clause. This option has not been tested and
   * may be extensive, see Clause::getNumeralWeight()
   */
  BoolOptionValue _increasedNumeralWeight;

  BoolOptionValue _ignoreConjectureInPreprocessing;

  IntOptionValue _inequalitySplitting;
  ChoiceOptionValue<InputSyntax> _inputSyntax;
  ChoiceOptionValue<Instantiation> _instantiation;
  BoolOptionValue _useHashingVariantIndex;

  ChoiceOptionValue<Induction> _induction;
  ChoiceOptionValue<StructuralInductionKind> _structInduction;
  ChoiceOptionValue<IntInductionKind> _intInduction;
  ChoiceOptionValue<InductionChoice> _inductionChoice;
  UnsignedOptionValue _maxInductionDepth;
  BoolOptionValue _inductionNegOnly;
  BoolOptionValue _inductionUnitOnly;
  BoolOptionValue _inductionGen;
  BoolOptionValue _inductionGenHeur;
  BoolOptionValue _inductionStrengthenHypothesis;
  UnsignedOptionValue _maxInductionGenSubsetSize;
  BoolOptionValue _inductionOnComplexTerms;
  BoolOptionValue _functionDefinitionRewriting;
  BoolOptionValue _integerInductionDefaultBound;
  ChoiceOptionValue<IntegerInductionInterval> _integerInductionInterval;
  ChoiceOptionValue<IntegerInductionLiteralStrictness> _integerInductionStrictnessEq;
  ChoiceOptionValue<IntegerInductionLiteralStrictness> _integerInductionStrictnessComp;
  ChoiceOptionValue<IntegerInductionTermStrictness> _integerInductionStrictnessTerm;
  BoolOptionValue _nonUnitInduction;
  BoolOptionValue _inductionOnActiveOccurrences;

  StringOptionValue _latexOutput;
  BoolOptionValue _latexUseDefaultSymbols;

  ChoiceOptionValue<LiteralComparisonMode> _literalComparisonMode;
  IntOptionValue _lookaheadDelay;
  IntOptionValue _lrsFirstTimeCheck;
  BoolOptionValue _lrsWeightLimitOnly;
  BoolOptionValue _lrsRetroactiveDeletes;

#if VAMPIRE_PERF_EXISTS
  UnsignedOptionValue _instructionLimit;
  UnsignedOptionValue _simulatedInstructionLimit;
  BoolOptionValue _parsingDoesNotCount;
#endif

  UnsignedOptionValue _memoryLimit; // should be size_t, making an assumption

  BoolOptionValue _interactive;

  ChoiceOptionValue<Mode> _mode;
  ChoiceOptionValue<Intent> _intent;
  ChoiceOptionValue<Schedule> _schedule;
  StringOptionValue _scheduleFile;
  UnsignedOptionValue _multicore;
  FloatOptionValue _slowness;
  BoolOptionValue _randomizSeedForPortfolioWorkers;

  IntOptionValue _naming;
  BoolOptionValue _nonliteralsInClauseWeight;
  BoolOptionValue _normalize;
  BoolOptionValue _shuffleInput;
  BoolOptionValue _randomPolarities;

  BoolOptionValue _outputAxiomNames;

  StringOptionValue _printProofToFile;
  BoolOptionValue _printClausifierPremises;
  StringOptionValue _problemName;
  ChoiceOptionValue<Proof> _proof;
  BoolOptionValue _minimizeSatProofs;
  ChoiceOptionValue<ProofExtra> _proofExtra;
  BoolOptionValue _traceback;

  StringOptionValue _protectedPrefix;

  ChoiceOptionValue<QuestionAnsweringMode> _questionAnswering;
  BoolOptionValue _questionAnsweringGroundOnly;
  StringOptionValue _questionAnsweringAvoidThese;

  UnsignedOptionValue _randomSeed;

  StringOptionValue _sampleStrategy;

  IntOptionValue _activationLimit;

  ChoiceOptionValue<SatSolver> _satSolver;
  ChoiceOptionValue<SaturationAlgorithm> _saturationAlgorithm;
  BoolOptionValue _showAll;
  BoolOptionValue _showActive;
  BoolOptionValue _showBlocked;
  BoolOptionValue _showDefinitions;
  ChoiceOptionValue<InterpolantMode> _showInterpolant;
  BoolOptionValue _showNew;
  BoolOptionValue _sineToAge;
  ChoiceOptionValue<PredicateSineLevels> _sineToPredLevels;
  BoolOptionValue _showSplitting;
  BoolOptionValue _showNewPropositional;
  BoolOptionValue _showNonconstantSkolemFunctionTrace;
  BoolOptionValue _showOptions;
  BoolOptionValue _showOptionsLineWrap;
  BoolOptionValue _showExperimentalOptions;
  BoolOptionValue _showHelp;
  BoolOptionValue _printAllTheoryAxioms;
  StringOptionValue _explainOption;
  BoolOptionValue _showPassive;
  BoolOptionValue _showReductions;
  BoolOptionValue _showPreprocessing;
  BoolOptionValue _showSkolemisations;
  BoolOptionValue _showSymbolElimination;
  BoolOptionValue _showTheoryAxioms;
  BoolOptionValue _showFOOL;
  BoolOptionValue _showFMBsortInfo;
  BoolOptionValue _showInduction;
  BoolOptionValue _showSimplOrdering;
#if VAMPIRE_CLAUSE_TRACING
  // TODO make unsigned option value
  IntOptionValue _traceBackward;
  IntOptionValue _traceForward;
#endif // VAMPIRE_CLAUSE_TRACING
#if VZ3
  BoolOptionValue _showZ3;
  ChoiceOptionValue<ProblemExportSyntax> _problemExportSyntax;
  StringOptionValue _exportAvatarProblem;
  StringOptionValue _exportThiProblem;
  BoolOptionValue _satFallbackForSMT;
  BoolOptionValue _smtForGround;
  ChoiceOptionValue<TheoryInstSimp> _theoryInstAndSimp;
  BoolOptionValue _thiGeneralise;
  BoolOptionValue _thiTautologyDeletion;
#endif
  ChoiceOptionValue<UnificationWithAbstraction> _unificationWithAbstraction;
  BoolOptionValue _unificationWithAbstractionFixedPointIteration;
  BoolOptionValue _fixUWA;
  BoolOptionValue _useACeval;
  TimeLimitOptionValue _simulatedTimeLimit;
  FloatOptionValue _lrsEstimateCorrectionCoef;
  UnsignedOptionValue _sineDepth;
  UnsignedOptionValue _sineGeneralityThreshold;
  UnsignedOptionValue _sineToAgeGeneralityThreshold;
  ChoiceOptionValue<SineSelection> _sineSelection;
  FloatOptionValue _sineTolerance;
  FloatOptionValue _sineToAgeTolerance;
  ChoiceOptionValue<Sos> _sos;
  UnsignedOptionValue _sosTheoryLimit;
  BoolOptionValue _splitting;
  BoolOptionValue _splitAtActivation;
  ChoiceOptionValue<SplittingAddComplementary> _splittingAddComplementary;
  ChoiceOptionValue<SplittingCongruenceClosure> _splittingCongruenceClosure;
  ChoiceOptionValue<CCUnsatCores> _ccUnsatCores;
  BoolOptionValue _splittingEagerRemoval;
  UnsignedOptionValue _splittingFlushPeriod;
  FloatOptionValue _splittingFlushQuotient;
  FloatOptionValue _splittingAvatimer;
  ChoiceOptionValue<SplittingNonsplittableComponents> _splittingNonsplittableComponents;
  ChoiceOptionValue<SplittingMinimizeModel> _splittingMinimizeModel;
  ChoiceOptionValue<SplittingLiteralPolarityAdvice> _splittingLiteralPolarityAdvice;
  ChoiceOptionValue<SplittingDeleteDeactivated> _splittingDeleteDeactivated;
  BoolOptionValue _splittingFastRestart;
  BoolOptionValue _splittingBufferedSolver;

  ChoiceOptionValue<Statistics> _statistics;
  BoolOptionValue _superpositionFromVariables;
  ChoiceOptionValue<TermOrdering> _termOrdering;
  ChoiceOptionValue<SymbolPrecedence> _symbolPrecedence;
  ChoiceOptionValue<SymbolPrecedenceBoost> _symbolPrecedenceBoost;
  ChoiceOptionValue<IntroducedSymbolPrecedence> _introducedSymbolPrecedence;
  ChoiceOptionValue<EvaluationMode> _evaluationMode;
  ChoiceOptionValue<KboWeightGenerationScheme> _kboWeightGenerationScheme;
  BoolOptionValue _kboMaxZero;
  ChoiceOptionValue<KboAdmissibilityCheck> _kboAdmissabilityCheck;
  StringOptionValue _functionWeights;
  StringOptionValue _predicateWeights;
  StringOptionValue _typeConPrecedence;
  StringOptionValue _functionPrecedence;
  StringOptionValue _predicatePrecedence;

  StringOptionValue _testId;
  ChoiceOptionValue<Output> _outputMode;
  BoolOptionValue _ignoreMissingInputsInUnsatCore;
  StringOptionValue _thanks;
  ChoiceOptionValue<TheoryAxiomLevel> _theoryAxioms;
  BoolOptionValue _theoryFlattening;
  BoolOptionValue _ignoreUnrecognizedLogic;

  /** Time limit in deciseconds */
  TimeLimitOptionValue _timeLimitInDeciseconds;
#if VTIME_PROFILING
  BoolOptionValue _timeStatistics;
  StringOptionValue _timeStatisticsFocus;
#endif // VTIME_PROFILING

  ChoiceOptionValue<URResolution> _unitResultingResolution;
  BoolOptionValue _unusedPredicateDefinitionRemoval;
  BoolOptionValue _blockedClauseElimination;
  UnsignedOptionValue _distinctGroupExpansionLimit;

  OptionChoiceValues _tagNames;

  NonGoalWeightOptionValue _nonGoalWeightCoefficient;
  BoolOptionValue _restrictNWCtoGC;

  SelectionOptionValue _selection;

  InputFileOptionValue _inputFile;

  BoolOptionValue _newCNF;
  BoolOptionValue _inlineLet;

  BoolOptionValue _manualClauseSelection;
  // arithmeitc reasoning options
  BoolOptionValue _inequalityNormalization;
  BoolOptionValue _pushUnaryMinus;
  BoolOptionValue _highSchool;
  ChoiceOptionValue<ArithmeticSimplificationMode> _gaussianVariableElimination;
  ChoiceOptionValue<ArithmeticSimplificationMode> _cancellation;
  ChoiceOptionValue<ArithmeticSimplificationMode> _arithmeticSubtermGeneralizations;


  //Higher-order options
  BoolOptionValue _addCombAxioms;
  BoolOptionValue _addProxyAxioms;
  BoolOptionValue _combinatorySuperposition;
  BoolOptionValue _choiceAxiom;
  BoolOptionValue _injectivity;
  BoolOptionValue _pragmatic;
  BoolOptionValue _choiceReasoning;
  BoolOptionValue _priortyToLongReducts;
  IntOptionValue  _maximumXXNarrows;
  ChoiceOptionValue<FunctionExtensionality> _functionExtensionality;
  ChoiceOptionValue<CNFOnTheFly> _clausificationOnTheFly;
  ChoiceOptionValue<PISet> _piSet;
  ChoiceOptionValue<Narrow> _narrow;
  BoolOptionValue _equalityToEquivalence;
  BoolOptionValue _complexBooleanReasoning;
  BoolOptionValue _booleanEqTrick;
  BoolOptionValue _superposition;
  BoolOptionValue _casesSimp;
  BoolOptionValue _cases;
  BoolOptionValue _newTautologyDel;
  BoolOptionValue _lambdaFreeHol;
  BoolOptionValue _complexVarCondition;

}; // class Options

// Allow printing of enums
template<typename T,
         typename = typename std::enable_if<std::is_enum<T>::value>::type>
std::ostream& operator<< (std::ostream& str,const T& val)
{
  return str << static_cast<typename std::underlying_type<T>::type>(val);
}

}

#endif
