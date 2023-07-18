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
 * @file Parse/TPTP.cpp
 * Implements class TPTP for parsing TPTP files
 *
 * @since 08/04/2011 Manchester
 */

#include <fstream>

#include "Debug/Assertion.hpp"

#include "Lib/Int.hpp"
#include "Lib/Environment.hpp"

#include "Kernel/Signature.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/FormulaUnit.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/Theory.hpp"
#include "Kernel/RobSubstitution.hpp"
#include "Kernel/FormulaVarIterator.hpp"

#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/DistinctGroupExpansion.hpp"
#include "Shell/UIHelper.hpp"

#include "Indexing/TermSharing.hpp"

#include "Parse/TPTP.hpp"

using namespace std;
using namespace Lib;
using namespace Kernel;
using namespace Shell;
using namespace Parse;

#define DEBUG_SHOW_TOKENS 0
#define DEBUG_SHOW_UNITS 0
#define DEBUG_SOURCE 0
DHMap<unsigned, std::string> TPTP::_axiomNames;
DHMap<unsigned, Map<int,std::string>*> TPTP::_questionVariableNames;

//Numbers chosen to avoid clashing with connectives.
//Unlikely to ever have 100 connectives, so this should be ok.
const int TPTP::HOL_CONSTANTS_LOWER_BOUND = 99u;
/** operator lambda */
const int TPTP::LAMBDA = 100u;
/** application of any number of terms */
const int TPTP::APP = 101u;
/** Pi function for universal quantification */
const int TPTP::PI = 102u;
/** Sigma function for existential quantification */
const int TPTP::SIGMA = 103u;

/**
 * Create a parser, parse the input and return the parsed list of units.
 * @since 13/07/2011 Manchester
 */
UnitList* TPTP::parse(istream& input)
{
  Parse::TPTP parser(input);
  parser.parse();
  return parser.units();
}

Clause* TPTP::parseClauseFromString(const std::string& str)
{
  std::stringstream input(str+")."); // to fake endFOF, which creates the clause
  Parse::TPTP parser(input);
  parser._isFof = true;
  parser._lastInputType = UnitInputType::AXIOM;
  parser._bools.push(false);     // this is what cnf normally pushes (but we start "from the middle")
  parser._strings.push("dummy_name");
  parser._states.push(END_FOF);  // this is what does the clause building
  parser.parseImpl(FORMULA);
  return parser._units.list()->head()->asClause();
}

/**
 * Initialise a lexer.
 * @since 27/07/2004 Torrevieja
 */
TPTP::TPTP(istream& in, UnitList::FIFO unitBuffer)
  : _containsConjecture(false),
    _allowedNames(0),
    _in(&in),
    _includeDirectory(""),
    _units(unitBuffer),
    _isThf(false),
    _containsPolymorphism(false),
    _currentColor(COLOR_TRANSPARENT),
    _lastPushed(TM),
    _modelDefinition(false),
    _insideEqualityArgument(0),
    _unitSources(0),
    _filterReserved(false),
    _seenConjecture(false)
{
} // TPTP::TPTP

/**
 * The destructor, does nothing.
 * @since 09/07/2012 Manchester
 */
TPTP::~TPTP()
{
} // TPTP::~TPTP

void TPTP::parse()
{
  try {
    parseImpl();
  } catch (UserErrorException &e) {
    e.line = lineNumber();
    throw;
  }
}

/**
 * Read all tokens one by one 
 * @since 08/04/2011 Manchester
 */
void TPTP::parseImpl(State initialState)
{
  // bulding tokens one by one
  _gpos = 0;
  _cend = 0;
  _tend = 0;
  _lineNumber = 1;
  _states.push(initialState);
  while (!_states.isEmpty()) {
    State s = _states.pop();
#ifdef DEBUG_SHOW_STATE
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << toString(s) << endl;
#endif
    switch (s) {
    case UNIT_LIST:
      unitList();
      break;
    case FOF:
      _isFof = true;
      fof(true);
      break;
    case THF:
      _isThf = true;
    case TFF:
      _isFof = false;
      tff();
      break;
    case CNF:
      _isFof = true;
      fof(false);
      break;
    case FORMULA:
      formula();
      break;
    case FUN_APP:
      funApp();
      break;
    case ARGS:
      args();
      break;
    case TERM:
      term();
      break;
    case TERM_INFIX:
      termInfix();
      break;
    case END_TERM:
      endTerm();
      break;
    case END_ARGS:
      endArgs();
      break;
    case FORMULA_INFIX:
      formulaInfix();
      break;
    case END_EQ:
      endEquality();
      break;
    case MID_EQ:
      midEquality();
      break;
    case VAR_LIST:
      varList();
      break;
    case TAG:
      tag();
      break;
    case END_FOF:
      endFof();
      break;
    case SIMPLE_FORMULA:
      simpleFormula();
      break;
    case END_FORMULA:
      endFormula();
      break;
    case END_APP:
      endApp();
      break;
    case HOL_FORMULA:
      holFormula();
      break;
    case END_HOL_FORMULA:
      endHolFormula();
      break;
    case HOL_TERM:
      holTerm();
      break;
    case FORMULA_INSIDE_TERM:
      formulaInsideTerm();
      break;
    case END_FORMULA_INSIDE_TERM:
      endFormulaInsideTerm();
      break;
    case END_TERM_AS_FORMULA:
      endTermAsFormula();
      break;
    case INCLUDE:
      include();
      break;
    case TYPE:
      type();
      break;
    case SIMPLE_TYPE:
      simpleType();
      break;
    case END_TYPE:
      endType();
      break;
    case END_TFF:
      endTff();
      break;
    case UNBIND_VARIABLES:
      unbindVariables();
      break;
    case VAMPIRE:
      vampire();
      break;
    case END_ITE:
      endIte();
      break;
    case LET_TYPE:
      letType();
      break;
    case END_LET_TYPES:
      endLetTypes();
      break;
    case DEFINITION:
      definition();
      break;
    case MID_DEFINITION:
      midDefinition();
      break;
    case END_DEFINITION:
      endDefinition();
      break;
    case SYMBOL_DEFINITION:
      symbolDefinition();
      break;
    case TUPLE_DEFINITION:
      if(!env.options->newCNF()){ USER_ERROR("Set --newcnf on if using tuples"); }
      tupleDefinition();
      break;
    case END_LET:
      endLet();
      break; 
    case END_THEORY_FUNCTION:
      endTheoryFunction();
      break;
    case END_TUPLE:
      if(!env.options->newCNF()){ USER_ERROR("Set --newcnf on if using tuples"); }
      endTuple();
      break;
    default:
#if VDEBUG
      throw ParseErrorException(((std::string)"Don't know how to process state ")+toString(s),_lineNumber);
#else
      throw ParseErrorException("Don't know how to process state ",_lineNumber);
#endif
    }
#ifdef DEBUG_SHOW_STATE
    cout << "----------------------------------------" << endl;
    printStacks();
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
#endif
  }
} // TPTP::parse()

/**
 * Return either the content or the string for this token
 * @since 11/04/2011 Manchester
 */
std::string TPTP::Token::toString() const
{
  std::string str = TPTP::toString(tag);
  return str == "" ? content : str;
} // Token::toString

/**
 * Return the string representation of this tag or "" is the representation
 * is not fixed (e.g. for T_NAME)
 * @since 11/04/2011 Manchester
 */
std::string TPTP::toString(Tag tag)
{
  switch (tag) {
  case T_EOF:
    return "<eof>";
  case T_LPAR:
    return "(";
  case T_RPAR:
    return ")";
  case T_LBRA:
    return "[";
  case T_RBRA:
    return "]";
  case T_COMMA:
    return ",";
  case T_COLON:
    return ":";
  case T_NOT:
    return "~";
  case T_AND:
    return "&";
  case T_EQUAL:
    return "=";
  case T_NEQ:
    return "!=";
  case T_FORALL:
    return "!";
  case T_EXISTS:
    return "?";
  case T_PI:
    return "??";
  case T_SIGMA:
    return "!!";
  case T_IMPLY:
    return "=>";
  case T_XOR:
    return "<~>";
  case T_IFF:
    return "<=>";
  case T_REVERSE_IMP:
    return "<=";
  case T_DOT:
    return ".";
  case T_OR:
    return "|";
  case T_ASS:
    return ":=";
  case T_LAMBDA:
    return "^";
  case T_APP:
    return "@";
  case T_STAR:
    return "*";
  case T_UNION:
    return "+";
  case T_ARROW:
    return ">";
  case T_SUBTYPE:
    return "<<";
  case T_NOT_OR:
    return "~|";
  case T_NOT_AND:
    return "~&";
  case T_SEQUENT:
    return "-->";
  case T_TYPE_QUANT:
    return "!>";
  case T_CHOICE:
    return "@+";
  case T_POLY_CHOICE:
    return "@@+";
  case T_DEF_DESC:
    return "@-";
  case T_POLY_DEF_DESC:
    return "@@-";
  case T_THF_QUANT_SOME:
    return "?*";
  case T_TRUE:
    return "$true";
  case T_FALSE:
    return "$false";
  case T_TTYPE:
    return "$tType";
  case T_BOOL_TYPE:
    return "$o";
  case T_DEFAULT_TYPE:
    return "$i";
  case T_RATIONAL_TYPE:
    return "$rat";
  case T_REAL_TYPE:
    return "$real";
  case T_INTEGER_TYPE:
    return "$int";
  case T_TUPLE:
    return "$tuple";
  case T_THEORY_SORT:
    return "";
  case T_THEORY_FUNCTION:
    return "";
  case T_FOT:
    return "$fot";
  case T_FOF:
    return "$fof";
  case T_TFF:
    return "$tff";
  case T_THF:
    return "$thf";
  case T_ITE:
    return "$ite";
  case T_LET:
    return "$let";
  case T_NAME:
  case T_REAL:
  case T_RAT:
  case T_INT:
  case T_VAR:
  case T_DOLLARS:
  case T_STRING:
    return "";
  default:
    ASSERTION_VIOLATION
  }
} // toString(Tag)

/**
 * Read all tokens one by one
 * @since 08/04/2011 Manchester
 */
bool TPTP::readToken(Token& tok)
{
  skipWhiteSpacesAndComments();
  tok.start = _gpos;
  switch (getChar(0)) {
  case 0:
    tok.tag = T_EOF;
    return false;
  case 'a':
  case 'b':
  case 'c':
  case 'd':
  case 'e':
  case 'f':
  case 'g':
  case 'h':
  case 'i':
  case 'j':
  case 'k':
  case 'l':
  case 'm':
  case 'n':
  case 'o':
  case 'p':
  case 'q':
  case 'r':
  case 's':
  case 't':
  case 'u':
  case 'v':
  case 'w':
  case 'x':
  case 'y':
  case 'z':
    tok.tag = T_NAME;
    readName(tok);
    return true;
  case '$':
    readReserved(tok);
    return true;
  case 'A':
  case 'B':
  case 'C':
  case 'D':
  case 'E':
  case 'F':
  case 'G':
  case 'H':
  case 'I':
  case 'J':
  case 'K':
  case 'L':
  case 'M':
  case 'N':
  case 'O':
  case 'P':
  case 'Q':
  case 'R':
  case 'S':
  case 'T':
  case 'U':
  case 'V':
  case 'W':
  case 'X':
  case 'Y':
  case 'Z':
  case '_':
    tok.tag = T_VAR;
    readName(tok);
    return true;
  case '0':
  case '1':
  case '2':
  case '3':
  case '4':
  case '5':
  case '6':
  case '7':
  case '8':
  case '9':
    tok.tag = readNumber(tok);
    return true;
  case '"':
    tok.tag = T_STRING;
    readString(tok);
    return true;
  case '\'':
    tok.tag = T_NAME;
    readAtom(tok);
    return true;
  case '(':
    tok.tag = T_LPAR;
    resetChars();
    return true;
  case ')':
    tok.tag = T_RPAR;
    resetChars();
    return true;
  case '[':
    tok.tag = T_LBRA;
    resetChars();
    return true;
  case ']':
    tok.tag = T_RBRA;
    resetChars();
    return true;
  case ',':
    tok.tag = T_COMMA;
    resetChars();
    return true;
  case ':':
    if (getChar(1) == '=') {
      tok.tag = T_ASS;
      resetChars();
      return true;
    }
    tok.tag = T_COLON;
    shiftChars(1);
    return true;
  case '~':
    if (getChar(1) == '&') {
      tok.tag = T_NOT_AND;
      resetChars();
      return true;
    }
    if (getChar(1) == '|') {
      tok.tag = T_NOT_OR;
      resetChars();
      return true;
    }
    tok.tag = T_NOT;
    shiftChars(1);
    return true;
  case '=':
    if (getChar(1) == '>') {
      tok.tag = T_IMPLY;
      resetChars();
      return true;
    }
    tok.tag = T_EQUAL;
    shiftChars(1);
    return true;
  case '&':
    tok.tag = T_AND;
    resetChars();
    return true;
  case '^':
    tok.tag = T_LAMBDA;
    resetChars();
    return true;
  case '@':
    if (getChar(1) == '+') {
      tok.tag = T_CHOICE;
      resetChars();
      return true;
    }
    if (getChar(1) == '-') {
      tok.tag = T_DEF_DESC;
      resetChars();
      return true;
    }
    if (getChar(1) == '@' && getChar(2) == '+'){
      tok.tag = T_POLY_CHOICE;
      resetChars();
      return true;
    }
    if (getChar(1) == '@' && getChar(2) == '-'){
      tok.tag = T_POLY_DEF_DESC;
      resetChars();
      return true;
    }    
    tok.tag = T_APP;
    shiftChars(1);
    return true;
  case '*':
    tok.tag = T_STAR;
    resetChars();
    return true;
  case '>':
    tok.tag = T_ARROW;
    resetChars();
    return true;
  case '!':
    if (getChar(1) == '=') {
      tok.tag = T_NEQ;
      resetChars();
      return true;
    }
    if (getChar(1) == '>') {
      tok.tag = T_TYPE_QUANT;
      resetChars();
      return true;
    }
    if (getChar(1) == '!') {
      tok.tag = T_PI;
      resetChars();
      return true;
    }
    tok.tag = T_FORALL;
    shiftChars(1);
    return true;
  case '?':
    if (getChar(1) == '?') {
      tok.tag = T_SIGMA;
      resetChars();
      return true;
    }
    if (getChar(1) == '*') {
      tok.tag = T_THF_QUANT_SOME;
      resetChars();
      return true;
    }
    tok.tag = T_EXISTS;
    shiftChars(1);
    return true;
  case '<':
    if (getChar(1) == '<') {
      tok.tag = T_SUBTYPE;
      resetChars();
      return true;
    }
    if (getChar(1) == '~' && getChar(2) == '>') {
      tok.tag = T_XOR;
      resetChars();
      return true;
    }
    if (getChar(1) != '=') {
      PARSE_ERROR("unrecognized symbol",_gpos);
    }
    if (getChar(2) == '>') {
      tok.tag = T_IFF;
      resetChars();
      return true;
    }
    tok.tag = T_REVERSE_IMP;
    shiftChars(2);
    return true;
  case '.':
    tok.tag = T_DOT;
    resetChars();
    return true;
  case '|':
    tok.tag = T_OR;
    resetChars();
    return true;
  case '-':
    if (getChar(1) == '-' && getChar(2) == '>') {
      tok.tag = T_SEQUENT;
      resetChars();
      return true;
    }
    tok.tag = readNumber(tok);
    return true;
  case '+':
    if (getChar(1) < '0' || getChar(1) > '9') {
      tok.tag = T_UNION;
      shiftChars(1);
      return true;
    }
    tok.tag = readNumber(tok);
    return true;
  default:
    PARSE_ERROR("Bad character",_gpos);
  }
} // TPTP::readToken()

/**
 * Skip all white spaces and comments in the input file
 * @since 08/04/2011 Manchester
 */
void TPTP::skipWhiteSpacesAndComments()
{
  for (;;) {
    switch (getChar(0)) {
    case 0: // end-of-file
      return;

    case '\n':
    case '\r':
      _lineNumber++;
    case ' ':
    case '\t':
    case '\f':
      resetChars();
      break;

    case '%': // end-of-line comment
    resetChars();
    for (int n=0;;n++) {
      int c = getChar(n);
      if (c == 0) {
        resetChars();
        getChar(0);
	return;
      }
      if (c == '\n') {
        _lineNumber++;
#if VDEBUG
        // Only check for Status if in preamble before any units read (also only in the top level file, not in includes)
        if(_units.list() == 0 && _inputs.isEmpty()){
          _chars[n]='\0';
          std::string cline(_chars.content());
          if(cline.find("Status")!=std::string::npos){
             if(cline.find("Theorem")!=std::string::npos){ UIHelper::setExpectingUnsat(); }
             else if(cline.find("Unsatisfiable")!=std::string::npos){ UIHelper::setExpectingUnsat(); }
             else if(cline.find("ContradictoryAxioms")!=std::string::npos){ UIHelper::setExpectingUnsat(); }
             else if(cline.find("Satisfiable")!=std::string::npos){ UIHelper::setExpectingSat(); }
             else if(cline.find("CounterSatisfiable")!=std::string::npos){ UIHelper::setExpectingSat(); }
          }
        }
#endif
        resetChars();
	break;
      }
    }
    break;

    case '/': // potential comment
      if (getChar(1) != '*') {
	return;
      }
      resetChars();
      // search for the end of this comment
      for (;;) {
	int c = getChar(0);
        if( c == '\n' || c == '\r'){ _lineNumber++; }
	if (!c) {
	  return;
	}
	resetChars();
	if (c != '*') {
	  continue;
	}
	// c == '*'
	c = getChar(0);
	resetChars();
	if (c != '/') {
	  continue;
	}
	break;
      }
      break;

    // skip to the end of comment
    default:
      return;
    }
  }
} // TPTP::skipWhiteSpacesAndComments

/**
 * Read the name
 * @since 08/04/2011 Manchester
 */
void TPTP::readName(Token& tok)
{
  for (int n = 1;;n++) {
    switch (getChar(n)) {
    case 'A':
    case 'B':
    case 'C':
    case 'D':
    case 'E':
    case 'F':
    case 'G':
    case 'H':
    case 'I':
    case 'J':
    case 'K':
    case 'L':
    case 'M':
    case 'N':
    case 'O':
    case 'P':
    case 'Q':
    case 'R':
    case 'S':
    case 'T':
    case 'U':
    case 'V':
    case 'W':
    case 'X':
    case 'Y':
    case 'Z':
    case '_':
    case 'a':
    case 'b':
    case 'c':
    case 'd':
    case 'e':
    case 'f':
    case 'g':
    case 'h':
    case 'i':
    case 'j':
    case 'k':
    case 'l':
    case 'm':
    case 'n':
    case 'o':
    case 'p':
    case 'q':
    case 'r':
    case 's':
    case 't':
    case 'u':
    case 'v':
    case 'w':
    case 'x':
    case 'y':
    case 'z':
    case '$':
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      break;
    default:
      ASS(_chars.content()[0] != '$');
      tok.content.assign(_chars.content(),n);
      shiftChars(n);
      return;
    }
  }
} // readName

/**
 * Read a reserved name (starting with a $)
 * @since 10/07/2011 Manchester
 */
void TPTP::readReserved(Token& tok)
{
  int n = 1;
  for (;;n++) {
    switch (getChar(n)) {
    case 'A':
    case 'B':
    case 'C':
    case 'D':
    case 'E':
    case 'F':
    case 'G':
    case 'H':
    case 'I':
    case 'J':
    case 'K':
    case 'L':
    case 'M':
    case 'N':
    case 'O':
    case 'P':
    case 'Q':
    case 'R':
    case 'S':
    case 'T':
    case 'U':
    case 'V':
    case 'W':
    case 'X':
    case 'Y':
    case 'Z':
    case '_':
    case 'a':
    case 'b':
    case 'c':
    case 'd':
    case 'e':
    case 'f':
    case 'g':
    case 'h':
    case 'i':
    case 'j':
    case 'k':
    case 'l':
    case 'm':
    case 'n':
    case 'o':
    case 'p':
    case 'q':
    case 'r':
    case 's':
    case 't':
    case 'u':
    case 'v':
    case 'w':
    case 'x':
    case 'y':
    case 'z':
    case '$':
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      break;
    default:
      tok.content.assign(_chars.content(),n);
      //shiftChars(n);
      goto out;
    }
  }
 out:
  if (tok.content == "$true") {
    tok.tag = T_TRUE;
  }
  else if (tok.content == "$false") {
    tok.tag = T_FALSE;
  }
  else if (tok.content == "$ite_f" || tok.content == "$ite_t" || tok.content == "$ite") {
    tok.tag = T_ITE;
    // $ite_t and $ite_f are left for compatibility, $ite is a generalisation of them
    tok.content = "$ite";
  }
  else if (tok.content == "$let_tt" || tok.content == "$let_tf" || tok.content == "$let_ft" || tok.content == "$let_ff" || tok.content == "$let") {
    tok.tag = T_LET;
    // all tokens of the form $let_XY are left for compatibility, $let is a generalisation of them
    tok.content = "$let";
  }
  else if (tok.content == "$tType") {
    tok.tag = T_TTYPE;
  }
  else if (tok.content == "$o" || tok.content == "$oType") {
    tok.tag = T_BOOL_TYPE;
  }
  else if (tok.content == "$i" || tok.content == "$iType") {
    tok.tag = T_DEFAULT_TYPE;
  }
  else if (tok.content == "$int") {
    tok.tag = T_INTEGER_TYPE;
  }
  else if (tok.content == "$rat") {
    tok.tag = T_RATIONAL_TYPE;
  }
  else if (tok.content == "$real") {
    tok.tag = T_REAL_TYPE;
  }
  else if (tok.content == "$tuple") {
      tok.tag = T_TUPLE;
  }
  else if (isTheoryFunction(tok.content)) {
    tok.tag = T_THEORY_FUNCTION;
  }
  else if (isTheorySort(tok.content)) {
    tok.tag = T_THEORY_SORT;
  }
  else if (tok.content == "$fot") {
    tok.tag = T_FOT;
  }
  else if (tok.content == "$fof") {
    tok.tag = T_FOF;
  }
  else if (tok.content == "$tff") {
    tok.tag = T_TFF;
  }
  else if (tok.content == "$thf") {
    tok.tag = T_THF;
  }
  else if (tok.content.substr(0,2) == "$$" && !_filterReserved) {
      tok.tag = T_DOLLARS;
  }
  else {
      
      // If _filterReserved is on then filter "$" from content
      if(_filterReserved){
          unsigned c=0;
          for(;;c++){ if(getChar(c)!='$') break;}
          shiftChars(c);
          n=n-c;
          tok.content.assign(_chars.content(),n);
      }
      
      tok.tag = T_NAME;
  }
  // Moved from above so that _filterReserved works
  shiftChars(n);
} // readReserved

/**
 * Read a string
 * @since 08/04/2011 Manchester
 */
void TPTP::readString(Token& tok)
{
  for (int n = 1;;n++) {
    int c = getChar(n);
    if (!c) {
      PARSE_ERROR("non-terminated string",_gpos);
    }
    if (c == '\\') { // escape
      c = getChar(++n);
      if (!c) {
        PARSE_ERROR("non-terminated string",_gpos);
      }
      continue;
    }
    if (c == '"') {
      tok.content.assign(_chars.content()+1,n-1);
      resetChars();
      return;
    }
  }
} // readString

/**
 * Read a quoted atom
 * @since 08/04/2011 Manchester
 */
void TPTP::readAtom(Token& tok)
{
  for (int n = 1;;n++) {
    int c = getChar(n);
    if (!c) {
      PARSE_ERROR("non-terminated quoted atom",_gpos);
    }
    if (c == '\\') { // escape
      c = getChar(++n);
      if (!c) {
        PARSE_ERROR("non-terminated quoted atom",_gpos);
      }
      continue;
    }
    if (c == '\'') {
      tok.content.assign(_chars.content()+1,n-1);
      resetChars();
      return;
    }
  }
} // readAtom

TPTP::ParseErrorException::ParseErrorException(std::string message,int pos, unsigned ln) : _ln(ln)
{
  _message = message + " at position " + Int::toString(pos);
} // TPTP::ParseErrorException::ParseErrorException

TPTP::ParseErrorException::ParseErrorException(std::string message,Token& tok, unsigned ln) : _ln(ln)
{
  _message = message + " at position " + Int::toString(tok.start) + " (text: " + tok.toString() + ')'; 
} // TPTP::ParseErrorException::ParseErrorException

/**
 * Exception printing a message. Currently computing a position is simplified
 * @since 08/04/2011 Manchester
 */
void TPTP::ParseErrorException::cry(ostream& str) const
{
  str << "Parsing Error on line " << _ln << ": ";
  str << _message << "\n";
}

/**
 * Read a number
 * @since 08/04/2011 Manchester
 */
TPTP::Tag TPTP::readNumber(Token& tok)
{
  // skip the sign
  int c = getChar(0);
  ASS(c);
  int pos = decimal((c == '+' || c == '-') ? 1 : 0);
  switch (getChar(pos)) {
  case '/':
    pos = positiveDecimal(pos+1);
    tok.content.assign(_chars.content(),pos);
    shiftChars(pos);
    return T_RAT;
  case 'E':
  case 'e':
    {
      char c = getChar(pos+1);
      pos = decimal((c == '+' || c == '-') ? pos+2 : pos+1);
      tok.content.assign(_chars.content(),pos);
      shiftChars(pos);
    }
    return T_REAL;
  case '.':
    {
      int p = pos;
      do {
        c = getChar(++pos);
      }
      while (c >= '0' && c <= '9');
      if (pos == p+1) {
        // something like 12.
        PARSE_ERROR("wrong number format",_gpos);
      }
      c = getChar(pos);
      if (c == 'e' || c == 'E') {
        c = getChar(pos+1);
        pos = decimal((c == '+' || c == '-') ? pos+2 : pos+1);
      }
      tok.content.assign(_chars.content(),pos);
      shiftChars(pos);
    }
    return T_REAL;
  default:
    tok.content.assign(_chars.content(),pos);
    shiftChars(pos);
    return T_INT;
  }
} // readNumber

/**
 * Read a decimal starting at position pos (see the TPTP grammar),
 * return the position after the last character in the decimal.
 * @since 08/04/2011 Manchester
 */
int TPTP::decimal(int pos)
{
  switch (getChar(pos)) {
  case '0':
    return pos+1;
  case '1':
  case '2':
  case '3':
  case '4':
  case '5':
  case '6':
  case '7':
  case '8':
  case '9':
    break;
  default:
    PARSE_ERROR("wrong number format",_gpos);
  }

  int c;
  do {
    c = getChar(++pos);
  }
  while (c >= '0' && c <= '9');
  return pos;
} // decimal

/**
 * Read a positive decimal starting at position pos (see the TPTP grammar),
 * return the position after the last character in the decimal.
 * @since 08/04/2011 Manchester
 */
int TPTP::positiveDecimal(int pos)
{
  switch (getChar(pos)) {
  case '1':
  case '2':
  case '3':
  case '4':
  case '5':
  case '6':
  case '7':
  case '8':
  case '9':
    break;
  default:
    PARSE_ERROR("wrong number format",_gpos);
  }

  int c;
  do {
    c = getChar(++pos);
  }
  while (c >= '0' && c <= '9');
  return pos;
} // positiveDecimal

/**
 * Process unit list declaration. If end of file is reached, terminates. Otherwise,
 * pushes on the state state UNIT_LIST and one of CNF, FOF, VAMPIRE
 * @since 10/04/2011 Manchester
 */
void TPTP::unitList()
{
  Token& tok = getTok(0);
  if (tok.tag == T_EOF) {
    resetToks();
    if (_inputs.isEmpty()) {
      return;
    }
    resetChars();
    delete _in;
    _in = _inputs.pop();
    _includeDirectory = _includeDirectories.pop();
    delete _allowedNames;
    _allowedNames = _allowedNamesStack.pop();
    _states.push(UNIT_LIST);
    return;
  }
  if (tok.tag != T_NAME) {
    PARSE_ERROR("cnf(), fof(), vampire() or include() expected",tok);
  }
  std::string name(tok.content);
  _states.push(UNIT_LIST);
  if (name == "cnf") {
    _states.push(CNF);
    resetToks();
    return;
  }
  if (name == "fof") {
      _states.push(FOF);
    resetToks();
    return;
  }
  if (name == "tff") {
    _states.push(TFF);
    resetToks();
    return;
  }
  if (name == "thf") {
    _states.push(THF);
    resetToks();
    return;
  }
  if (name == "vampire") {
    _states.push(VAMPIRE);
    resetToks();
    return;
  }
  if (name == "include") {
    _states.push(INCLUDE);
    resetToks();
    return;
  }
  PARSE_ERROR("cnf(), fof(), vampire() or include() expected",tok);
}

/**
 * Process fof() or cnf() declaration. Does the following:
 * <ol>
 *  <li>add 0 to _formulas</li>
 *  <li>save the input type to _lastInputType</li>
 *  <li>add unit name to _strings</li>
 *  <li>add to _states END_FOF,FORMULA</li>
 *  <li>adds to _bools true, if fof and false, if cnf</li>
 * </ol>
 * @since 10/04/2011 Manchester
 */
void TPTP::fof(bool fo)
{
  _bools.push(fo);
  consumeToken(T_LPAR);
  // save the name of this unit
  Token& tok = getTok(0);
  switch(tok.tag) {
  case T_NAME:
    _strings.push(tok.content);
    resetToks();
    break;
  case T_INT:
    _strings.push(tok.content);
    resetToks();
    break;
  default:
    PARSE_ERROR("Unit name expected",tok);
  }

  consumeToken(T_COMMA);
  tok = getTok(0);
  int start = tok.start;
  std::string tp = name();

  _isQuestion = false;
  if(_modelDefinition){
    _lastInputType = UnitInputType::MODEL_DEFINITION;
  }
  else if (tp == "axiom" || tp == "plain") {
    _lastInputType = UnitInputType::AXIOM;
  }
  else if(tp == "extensionality"){
    // this will be transformed to just AXIOM after clausification
    _lastInputType = UnitInputType::EXTENSIONALITY_AXIOM;
  }
  else if (tp == "definition") {
    _lastInputType = UnitInputType::AXIOM;
  }
  else if (tp == "conjecture") {
    _containsConjecture = true;
    _lastInputType = UnitInputType::CONJECTURE;
  }
  else if (tp == "question") {
    _isQuestion = true;
    _containsConjecture = true;
    _lastInputType = UnitInputType::CONJECTURE;
  }
  else if (tp == "negated_conjecture") {
    _lastInputType = UnitInputType::NEGATED_CONJECTURE;
  }
  else if (tp == "hypothesis" || tp == "theorem" || tp == "lemma") {
    _lastInputType = UnitInputType::ASSUMPTION;
  }
  else if (tp == "claim") {
    _lastInputType = UnitInputType::CLAIM;
  }
  else if (tp == "assumption" || tp == "unknown") {
    // MS: we were silently dropping these until now. I wonder why...
    USER_ERROR("Unsupported unit type '", tp, "' found");
  }
  else {
    PARSE_ERROR((std::string)"unit type, such as axiom or definition expected but " + tp + " found",
        start);
  }
  consumeToken(T_COMMA);
  _states.push(END_FOF);
  _states.push(FORMULA);
} // fof()

/**
 * Process fof() or cnf() declaration. Does the following:
 * <ol>
 *  <li>add 0 to _formulas</li>
 *  <li>save the input type to _lastInputType</li>
 *  <li>add unit name to _strings</li>
 *  <li>add to _states END_FOF,FORMULA</li>
 *  <li>adds to _bools true, if fof and false, if cnf</li>
 * </ol>
 * @since 10/04/2011 Manchester
 * @author Andrei Voronkov
 */
void TPTP::tff()
{
  consumeToken(T_LPAR);
  // save the name of this unit
  Token& tok = getTok(0);
  switch(tok.tag) {
  case T_NAME:
  case T_INT:
    _strings.push(tok.content);
    resetToks();
    break;
  default:
    PARSE_ERROR("Unit name expected",tok);
  }

  consumeToken(T_COMMA);
  tok = getTok(0);
  int start = tok.start;
  std::string tp = name();
  if (tp == "type") {
    // Read a TPTP type declaration.
    consumeToken(T_COMMA);
    // TPTP syntax allows for an arbitrary number of parentheses around a type
    // declaration
    int lpars = 0;
    for (;;) {
      tok = getTok(0);
      if (tok.tag != T_LPAR) {
        break;
      }
      lpars++;
      resetToks();
    }
    std::string nm = name();
    consumeToken(T_COLON);
    if(_isThf){
      tok = getTok(0);
      if (tok.tag == T_TTYPE) {
        resetToks();
        unsigned arity = getConstructorArity();
        bool added = false;
        unsigned fun = arity == 0
            ? addUninterpretedConstant(nm, _overflow, added)
            : env.signature->addFunction(nm, arity, added);
        Signature::Symbol* symbol = env.signature->getFunction(fun);
        OperatorType* ot = OperatorType::getTypeConType(arity);
        if (!added) {
          if(symbol->fnType()!=ot){
            PARSE_ERROR("Type constructor declared with two different types",tok);
          }
        } else{
          symbol->setType(ot);  
          _typeConstructorArities.insert(nm, arity);
        }       
        //cout << "added type constuctor " + nm + " of type " + symbol->fnType()->toString() << endl;
        while (lpars--) {
          consumeToken(T_RPAR);
        }
        consumeToken(T_RPAR);
        consumeToken(T_DOT);
        return;
      }
    }
    // the matching number of rpars will be read
    _ints.push(lpars);
    // remember type name
    _strings.push(nm);
    _states.push(END_TFF);
    _states.push(TYPE);
    return;
  }

  _bools.push(true); // to denote that it is an FOF formula
  _isQuestion = false;
  if(_modelDefinition){
    _lastInputType = UnitInputType::MODEL_DEFINITION;
  }
  else if (tp == "axiom" || tp == "plain") {
    _lastInputType = UnitInputType::AXIOM;
  }
  else if (tp == "extensionality"){
    // this will be transformed to just AXIOM after clausification
    _lastInputType = UnitInputType::EXTENSIONALITY_AXIOM;
  }
  else if (tp == "definition") {
    _lastInputType = UnitInputType::AXIOM;
  }
  else if (tp == "conjecture") {
    _containsConjecture = true;
    _lastInputType = UnitInputType::CONJECTURE;
  }
  else if (tp == "question") {
    _isQuestion = true;
    _containsConjecture = true;
    _lastInputType = UnitInputType::CONJECTURE;
  }
  else if (tp == "negated_conjecture") {
    _lastInputType = UnitInputType::NEGATED_CONJECTURE;
  }
  else if (tp == "hypothesis" || tp == "theorem" || tp == "lemma") {
    _lastInputType = UnitInputType::ASSUMPTION;
  }
  else if (tp == "assumption" || tp == "unknown") {
    // MS: we were silently dropping these until now. I wonder why...
    USER_ERROR("Unsupported unit type '", tp);
  }
  else if (tp == "claim") {
    _lastInputType = UnitInputType::CLAIM;
  }
  else {
    PARSE_ERROR((std::string)"unit type, such as axiom or definition expected but " + tp + " found",
        start);
  }
  consumeToken(T_COMMA);
  _states.push(END_FOF);
  _states.push(FORMULA);
} // tff()

unsigned TPTP::getConstructorArity()
{
  unsigned arity = 0;
  Token tok = getTok(0);
  while(tok.tag == T_ARROW || tok.tag == T_TTYPE){
    arity += (tok.tag == T_TTYPE);
    resetToks();
    tok = getTok(0);
  }
  return arity;
}

/**
  * Reads a HOL term of type $o
  * @since 08/11/2017
  * @author Ahmed Bhayat
  */

void TPTP::holFormula()
{
  Token tok = getTok(0);
  
  switch (tok.tag) {
  case T_NOT:
    resetToks();
    _connectives.push(NOT);
    _states.push(HOL_FORMULA);
    return;

  case T_SIGMA:
    resetToks();
    readTypeArgs(1);
    _termLists.push(createFunctionApplication("vSIGMA", 1));    
    return;
  
  case T_PI:
    resetToks();
    readTypeArgs(1);
    _termLists.push(createFunctionApplication("vPI", 1));      
    return;

  case T_FORALL:
  case T_EXISTS:
   // _states.push(UNBIND_VARIABLES);
  case T_LAMBDA:
    resetToks();
    consumeToken(T_LBRA);
    _connectives.push(tok.tag == T_FORALL ? FORALL : (tok.tag == T_LAMBDA ? LAMBDA : EXISTS));
    _states.push(HOL_FORMULA);
    addTagState(T_COLON);
    addTagState(T_RBRA);
    _states.push(VAR_LIST);
    return;

  case T_LPAR:
    resetToks();
    addTagState(T_RPAR);
    _connectives.push(-1);
    _states.push(END_HOL_FORMULA);
    _states.push(HOL_FORMULA);
    return;
    
  //higher order syntax wierdly allows (~) @ (...)
  case T_RPAR: {
    ASS(_connectives.top() == NOT);
    _connectives.pop();
    _termLists.push(createFunctionApplication("vNOT", 0));
    return;
  }

  case T_CHOICE:
  case T_DEF_DESC:
  case T_POLY_CHOICE:
  case T_POLY_DEF_DESC:
  {
    USER_ERROR("At the moment Vampire HOL cannot parse definite and indefinite description operators");
  }

  case T_STRING:
  case T_INT:
  case T_RAT:
  case T_REAL://TODO update for HOL - AYB
    _states.push(END_EQ);
    _states.push(TERM);
    _states.push(MID_EQ);
    _states.push(TERM);
    return;
  case T_TRUE:
    resetToks();
    _formulas.push(new Formula(true));
    _lastPushed = FORM;
    return;
  case T_FALSE:
    resetToks();
    _formulas.push(new Formula(false));
    _lastPushed = FORM;
    return;
  case T_AND:
  case T_OR:
  case T_IMPLY:
  case T_IFF:
  case T_NAME:
  case T_VAR:
  case T_ITE:
  case T_THEORY_FUNCTION:
  case T_LET:
  case T_LBRA:
    //_states.push(END_HOL_TERM);
    _states.push(HOL_TERM);
    return;
  case T_APP:
    //higher-order syntax allows for ~ @ fomrula  
    if(_connectives.top() == NOT ||
       _connectives.top() == PI  ||
       _connectives.top() == SIGMA){
      resetToks();
      _states.push(HOL_FORMULA);
      return;
    }
  //AYB ADDED, TO BE MODIFIED
  default:
    PARSE_ERROR("formula or term expected",tok);
  }
}

/**
  * Reads a term of a higher order logic. Either a higher-order constant or a variable
  * @since 08/11/2017
  * @author Ahmed Bhayat
  */

void TPTP::holTerm()
{
  Token tok = getTok(0);
  resetToks();

  std::string name = tok.content;
  unsigned arity = _typeArities.find(name) ? _typeArities.get(name) : 0;

  switch (tok.tag) {

    case T_BOOL_TYPE:
    case T_DEFAULT_TYPE: {
      resetToks();
      switch (tok.tag) {
        case T_BOOL_TYPE:
          _termLists.push(AtomicSort::boolSort());
          break;
        case T_DEFAULT_TYPE:
          _termLists.push(AtomicSort::defaultSort());
          break;             
        default:
          ASSERTION_VIOLATION;
      }
      return;
    }

    case T_REAL_TYPE:
    case T_RATIONAL_TYPE: 
    case T_STRING:
    case T_INT:
    case T_REAL:
    case T_RAT: {
      USER_ERROR("Vampire higher-order is currently not compatible with theory reasoning");
    }  

    case T_AND:
    case T_OR:
    case T_IMPLY:  
    case T_IFF:
    case T_XOR:{
      ASS(arity == 0);
      name = convert(tok.tag);
      _termLists.push(createFunctionApplication(name, arity)); //TODO fix this
      break;
    }    
    case T_NAME:{
      //AYB must be a nicer way of dealing with this?
      if(name.at(0) == '$'){
        USER_ERROR("Vampire higher-order is currently not compatible with theory reasoning");    
      }
      readTypeArgs(arity);
      _termLists.push(createFunctionApplication(name, arity)); // arity
      break;
    }
    case T_VAR:{
      unsigned var = (unsigned)_vars.insert(name);
      _termLists.push(TermList(var, false)); // dummy arity to indicate a variable
      break;
    }
    default:
      PARSE_ERROR("unexpected token", tok);
  }

  _lastPushed = TM;

}
  
std::string TPTP::convert(Tag t)
{
  switch(t){
    case T_AND:
      return "vAND";
    case T_OR:
      return "vOR";
    case T_IMPLY:
      return "vIMP";
    case T_IFF:
      return "vIFF";
    case T_XOR:
      return "vXOR";
    case T_NOT:
      return "vNOT";
    case T_PI:
      return "vPI";
    case T_SIGMA:
      return "vSIGMA";
    default:
      ASSERTION_VIOLATION;
  }
}


/**
  * Process the end of a HOL function.
  * If the main connective is any that operates on formulas (!, ?, &, |, -2...) and
  * The top term on the termlist is not of type $o and error is raised. Otherwise, 
  * @since 05/11/2017 Manchester
  * @author Ahmed Bhayat
  */

void TPTP::endHolFormula()
{
  int con = _connectives.pop();

  if (con == -2){
    if(_termLists.size() == 1){
      endTermAsFormula();
    }
    return;
  }  
  
  if ((con < HOL_CONSTANTS_LOWER_BOUND) && (con != -1) && (_lastPushed == TM)){
    //At the moment, APP and LAMBDA are the only connectives that can take terms of type
    //Other than $o as arguments.
    endTermAsFormula();
  }


  Formula* f;
  TermList fun;
  bool conReverse = false;
  switch (con) {
  case IMP:
  case AND:
  case OR:
    conReverse = _bools.pop();
    break;
  case IFF:
  case XOR:
  case APP:
  case -2:
  case -1:
    break;
  case NOT:
    f = _formulas.pop();
    _formulas.push(new NegatedFormula(f));
    _lastPushed = FORM;
    _states.push(END_HOL_FORMULA);

    return;

  case FORALL:
  case EXISTS:
    f = _formulas.pop();
    _formulas.push(new QuantifiedFormula((Connective)con,_varLists.pop(),_sortLists.pop(),f));
    _lastPushed = FORM;
    _states.push(END_HOL_FORMULA);
    _states.push(UNBIND_VARIABLES);
    return;
  case LAMBDA:{
     if(_lastPushed == FORM){
       endFormulaInsideTerm();
     }
     fun = _termLists.pop();
     TermList ts(Term::createLambda(fun, _varLists.pop(), _sortLists.pop(), sortOf(fun)));
     _termLists.push(ts);
     _lastPushed = TM;
     _states.push(END_HOL_FORMULA);
     _states.push(UNBIND_VARIABLES);
     return; 
    }
  case LITERAL:
  default:
    throw ::Exception((std::string)"tell me how to handle connective " + Int::toString(con));
  }

  Token& tok = getTok(0);
  Tag tag = tok.tag;
  int c;
  bool cReverse = false;
switch (tag) {
  case T_AND:
    c = AND;
    break;
  case T_NOT_AND:
    cReverse = true;
    c = AND;
    break;
  case T_NOT_OR:
    cReverse = true;
    c = OR;
    break;
  case T_OR:
    c = OR;
    break;
  case T_XOR:
    c = XOR;
    break;
  case T_IFF:
    c = IFF;
    break;
  case T_IMPLY:
    c = IMP;
    break;
  case T_REVERSE_IMP:
    cReverse = true;
    c = IMP;
    break;
  case T_APP:
    c = APP;
    break;
  case T_EQUAL:
  case T_NEQ: {
    // not connectives, but we allow formulas to be arguments to = and !=
    _states.push(END_EQ);
    _connectives.push(-1);
    _states.push(END_HOL_FORMULA);
    _states.push(HOL_FORMULA);
    _states.push(MID_EQ);
    if(_lastPushed == FORM){
      endFormulaInsideTerm();
      //equality is evaluated between two terms
    }
    return;
  }
  default:
    // the formula does not end at a binary connective, build the formula and terminate
    switch (con) {
    case IMP:
      f = _formulas.pop();
      if (conReverse) {
        f = new BinaryFormula((Connective)con,f,_formulas.pop());
      }
      else {
        f = new BinaryFormula((Connective)con,_formulas.pop(),f);
      }
      _formulas.push(f);
        _lastPushed = FORM;
      _states.push(END_HOL_FORMULA);
      return;
    
    case APP:
      _states.push(END_HOL_FORMULA);
      _states.push(END_APP);
      return;
    case IFF:
    case XOR:
      f = _formulas.pop();
      f = new BinaryFormula((Connective)con,_formulas.pop(),f);
      _formulas.push(f);
        _lastPushed = FORM;
      _states.push(END_HOL_FORMULA);
      return;

    case AND:
    case OR:
      f = _formulas.pop();
      f = makeJunction((Connective)con,_formulas.pop(),f);
      if (conReverse) {
        f = new NegatedFormula(f);
      }
      _formulas.push(f);
      _lastPushed = FORM;
      _states.push(END_HOL_FORMULA);
      return;

    case -1:
      return;
    default:
      ASSERTION_VIOLATION;
    }
  }

  if ((c != APP) && (con == -1) && (_lastPushed == TM)){
    endTermAsFormula();
  }

  
  // con and c are binary connectives
  if (higherPrecedence(con,c)) {
    if (con == APP){
      _states.push(END_HOL_FORMULA);
      _states.push(END_APP);
      return;  
    }
    f = _formulas.pop(); 
    Formula* g = _formulas.pop();
    if (con == AND || con == OR) {
      f = makeJunction((Connective)con,g,f);
      if (conReverse) {
        f = new NegatedFormula(f);
      }
    }
    else if (con == IMP && conReverse) {
      f = new BinaryFormula((Connective)con,f,g); 
    }else {
      f = new BinaryFormula((Connective)con,g,f);
    }
    _formulas.push(f);
    _lastPushed = FORM;
    _states.push(END_HOL_FORMULA);
    return;
  }

  // c is a binary connective
  _connectives.push(con);
  if (con == IMP || con == AND || con == OR) {
    _bools.push(conReverse);
  }
  _connectives.push(c);
  if (c == IMP || c == AND || c == OR) {
    _bools.push(cReverse);
  }
  resetToks();
  _states.push(END_HOL_FORMULA);
  _states.push(HOL_FORMULA);
}


/**
  * Process the end of an @ term
  * @since 05/11/2017 Manchester
  * @author Ahmed Bhayat
  */
void TPTP::endApp()
{
  if(_lastPushed == FORM){
     endFormulaInsideTerm();     
  }
  TermStack args;
  TermList rhs = _termLists.pop();
  TermList lhs = _termLists.pop();
  TermList lhsSort = sortOf(lhs);
  ASS_REP2(lhsSort.isTerm() && lhsSort.term()->arity() == 2, lhs.toString(), lhsSort.toString());
  TermList s1 = *(lhsSort.term()->nthArgument(0));
  TermList s2 = *(lhsSort.term()->nthArgument(1));
  args.push(s1);
  args.push(s2);
  args.push(lhs);
  args.push(rhs);
  unsigned app = env.signature->getApp();
  
  _termLists.push(TermList(Term::create(app, 4, args.begin())));
  _lastPushed = TM;
}

/**
 * Process the end of the $ite expression
 * @since 27/07/2011 Manchester
 * @since 16/04/2015 Gothenburg, major changes to support FOOL
 */
void TPTP::endIte()
{
  TermList elseBranch = _termLists.pop();
  TermList thenBranch = _termLists.pop();
  Formula* condition = _formulas.pop();
  TermList thenSort = sortOf(thenBranch);
  TermList ts(Term::createITE(condition,thenBranch,elseBranch,thenSort));
  TermList elseSort = sortOf(elseBranch);
  if (thenSort != elseSort) {
    USER_ERROR("sort mismatch in the if-then-else expression: " +
               thenBranch.toString() + " has the sort " + thenSort.toString() + ", whereas " +
               elseBranch.toString() + " has the sort " + elseSort.toString());
  }
  _termLists.push(ts);
} // endIte

/**
 *
 */
void TPTP::endTheoryFunction() {
  /**
   * Things get a bit awkward with theories + FOOL, because theory function can
   * return $o in such case be a predicate symbol rather than a function symbol.
   * The current solution is the following -- we always treat application of
   * theory functions as a a term (a formula wrapped inside boolean term, if
   * needed). If later on we discover that we should've taken it as a formula,
   * we simply pull the formula out of the boolean term. This is done in
   * endTermAsFormula().
   */

  Theory::Interpretation itp;
  TermList args[3]; // all theory function use up to 3 arguments as for now
  TermList arraySort;

  TheoryFunction tf = _theoryFunctions.pop();
  switch (tf) {
    case TF_SELECT: {
      TermList index = _termLists.pop();
      TermList array = _termLists.pop();

      arraySort = sortOf(array);
      if (!arraySort.isArraySort()) {
        USER_ERROR("$select is being incorrectly used on a type of array " + arraySort.toString() + " that has not be defined");
      }

      TermList indexSort = SortHelper::getIndexSort(arraySort);
      if (sortOf(index) != indexSort) {
        USER_ERROR("sort of index is not the same as the index sort of the array");
      }

      args[0] = array;
      args[1] = index;

      if (SortHelper::getInnerSort(arraySort) == AtomicSort::boolSort()) {
        itp = Theory::Interpretation::ARRAY_BOOL_SELECT;
      } else {
        itp = Theory::Interpretation::ARRAY_SELECT;
      }
      break;
    }
    case TF_STORE: {
      TermList value = _termLists.pop();
      TermList index = _termLists.pop();
      TermList array = _termLists.pop();

      arraySort = sortOf(array);
      if (!arraySort.isArraySort()) {
        USER_ERROR("store is being incorrectly used on a type of array that has not be defined");
      }

      TermList indexSort = SortHelper::getIndexSort(arraySort);
      if (sortOf(index) != indexSort) {
        USER_ERROR("sort of index is not the same as the index sort of the array");
      }

      TermList innerSort = SortHelper::getInnerSort(arraySort);
      if (sortOf(value) != innerSort) {
        USER_ERROR("sort of value is not the same as the value sort of the array");
      }

      args[0] = array;
      args[1] = index;
      args[2] = value;

      itp = Theory::Interpretation::ARRAY_STORE;

      break;
    }
    default:
      ASSERTION_VIOLATION_REP(tf);
  }

  OperatorType* type = Theory::getArrayOperatorType(arraySort,itp);
  unsigned symbol = env.signature->getInterpretingSymbol(itp, type);
  unsigned arity = Theory::getArity(itp);

  if (Theory::isFunction(itp)) {
    Term* term = Term::create(symbol, arity, args);
    _termLists.push(TermList(term));
  } else {
    Literal* literal = Literal::create(symbol, arity, true, args);
    _formulas.push(new AtomicFormula(literal));
    _states.push(END_FORMULA_INSIDE_TERM);
  }
} // endTheoryFunction

/**
 * Process include() declaration
 * @since 07/07/2011 Manchester
 */
void TPTP::include()
{
  consumeToken(T_LPAR);
  Token& tok = getTok(0);
  if (tok.tag != T_NAME) {
    PARSE_ERROR((std::string)"file name expected",tok);
  }
  std::string relativeName=tok.content;
  resetToks();
  bool ignore = _forbiddenIncludes.contains(relativeName);
  if (!ignore) {
    _allowedNamesStack.push(_allowedNames);
    _allowedNames = 0;
    _inputs.push(_in);
    _includeDirectories.push(_includeDirectory);
  }

  tok = getTok(0);
  if (tok.tag == T_COMMA) {
    if (!ignore) {
      _allowedNames = new Set<std::string>;
    }
    resetToks();
    consumeToken(T_LBRA);
    for(;;) {
      tok = getTok(0);
      if (tok.tag != T_NAME) {
	PARSE_ERROR((std::string)"formula name expected",tok);
      }
      std::string axName=tok.content;
      resetToks();
      if (!ignore) {
	_allowedNames->insert(axName);
      }
      tok = getTok(0);
      if (tok.tag == T_RBRA) {
	resetToks();
	break;
      }
      consumeToken(T_COMMA);
    }
  }
  consumeToken(T_RPAR);
  consumeToken(T_DOT);

  if (ignore) {
    return;
  }
  // here should be a computation of the new include directory according to
  // the TPTP standard, so far we just set it to ""
  _includeDirectory = "";
  std::string fileName(env.options->includeFileName(relativeName));
  _in = new ifstream(fileName.c_str());
  if (!*_in) {
    USER_ERROR((std::string)"cannot open file " + fileName);
  }
} // include

/** add a file name to the list of forbidden includes */
void TPTP::addForbiddenInclude(std::string file)
{
  _forbiddenIncludes.insert(file);
}

/**
 * Read the next token that must be a name.
 * @since 10/04/2011 Manchester
 */
std::string TPTP::name()
{
  Token& tok = getTok(0);
  if (tok.tag != T_NAME) {
    PARSE_ERROR("name expected",tok);
  }
  std::string nm = tok.content;
  resetToks();
  return nm;
} // name

/**
 * Read the next token that must have a given name.
 * @since 10/04/2011 Manchester
 */
void TPTP::consumeToken(Tag t)
{
  Token& tok = getTok(0);
  if (tok.tag != t) {
    std::string expected = toString(t);
    PARSE_ERROR(expected + " expected",tok);
  }
  resetToks();
} // consumeToken

/**
 * Read a formula and save it on the stack of formulas.
 * Adds to _states END_SIMPLE_FORMULA,SIMPLE_FORMULA
 * @since 10/04/2011 Manchester
 */
void TPTP::formula()
{
  if(_isThf){
    _connectives.push(-2); //special connective for HOL funcs
    _connectives.push(-1);
    _states.push(END_HOL_FORMULA);
    _states.push(END_HOL_FORMULA);
    _states.push(HOL_FORMULA);
  }else{
    _connectives.push(-1);
    _states.push(END_FORMULA);
    _states.push(SIMPLE_FORMULA);
  }
} // formula

/**
 *
 * @since 26/03/2015 Manchester
 */
void TPTP::termInfix()
{
  Token tok = getTok(0);
  switch (tok.tag) {
    case T_EQUAL:
    case T_NEQ:
      _states.push(END_FORMULA_INSIDE_TERM);
      _states.push(FORMULA_INFIX);
      return;
    case T_COMMA:
    case T_RPAR:
    case T_RBRA:
    case T_ASS:
      _states.push(END_TERM);
      return;
    case T_AND:
    case T_NOT_AND:
    case T_NOT_OR:
    case T_OR:
    case T_XOR:
    case T_IFF:
    case T_IMPLY:
    case T_REVERSE_IMP:
      if (_insideEqualityArgument > 0) {
        _states.push(END_TERM);
        return;
      }
      _connectives.push(-1);
      _states.push(END_FORMULA_INSIDE_TERM);
      _states.push(END_FORMULA);
      _states.push(FORMULA_INFIX);
      return;
    default:
      PARSE_ERROR("term or formula expected", tok);
  }
} // termInfix

/**
 * Read a TPTP type expression
 * @since 10/04/2011 Manchester
 * @author Andrei Voronkov
 */
void TPTP::type()
{
  _typeTags.push(TT_ATOMIC);
  _states.push(END_TYPE);
  _states.push(SIMPLE_TYPE);
} // type

/**
 * Read a function application or a variable and save the resulting literal
 * @since 10/04/2011 Manchester
 */
void TPTP::funApp()
{
  Token tok = getTok(0);
  resetToks();

  if (tok.tag == T_LBRA) {
    _strings.push(toString(T_TUPLE));
  } else {
    _strings.push(tok.content);
  }

  switch (tok.tag) {
    case T_THEORY_FUNCTION:
      consumeToken(T_LPAR);
      addTagState(T_RPAR);
      switch (getTheoryFunction(tok)) {
        case TF_SELECT:
          _states.push(TERM);
          addTagState(T_COMMA);
          _states.push(TERM);
          break;
        case TF_STORE:
          _states.push(TERM);
          addTagState(T_COMMA);
          _states.push(TERM);
          addTagState(T_COMMA);
          _states.push(TERM);
          break;
        default:
          ASSERTION_VIOLATION_REP(tok.content);
      }
      return;

    case T_ITE:
      consumeToken(T_LPAR);
      addTagState(T_RPAR);
      _states.push(TERM);
      addTagState(T_COMMA);
      _states.push(TERM);
      addTagState(T_COMMA);
      _states.push(FORMULA);
      return;

    case T_LET: {
      consumeToken(T_LPAR);
      addTagState(T_RPAR);
      _states.push(TERM);
      addTagState(T_COMMA);
      _states.push(DEFINITION);
      _letDefinitions.push(LetDefinitions());
      addTagState(T_COMMA);

      bool multipleLetTypes = false;
      if (getTok(0).tag == T_LBRA) {
        resetToks();
        addTagState(T_RBRA);
        multipleLetTypes = true;
      }
      _bools.push(multipleLetTypes);

      _states.push(END_LET_TYPES);
      _states.push(LET_TYPE);
      _letTypedSymbols.push(LetSymbols());
      return;
    }

    case T_LBRA:
      _states.push(ARGS);
      _ints.push(1); // the arity of the function symbol is at least 1
      return;

    case T_VAR:
      _ints.push(-1); // dummy arity to indicate a variable
      return;

    case T_NAME:
      if (getTok(0).tag == T_LPAR) {
        resetToks();
        _states.push(ARGS);
        _ints.push(1); // the arity of the function symbol is at least 1
      } else {
        _ints.push(0); // arity
      }
      return;

    default:
      PARSE_ERROR("unexpected token", tok);
  }
} // TPTP::funApp

void TPTP::letType()
{
  // We cannot use this method in TPTP::tff() because type declarations in the
  // "type" role TFF units allow declarations of types ($tType), which are not
  // allowed inside $lets

  _states.push(TYPE);
  addTagState(T_COLON);
  _strings.push(name());
} // TPTP::letType

void TPTP::endLetTypes()
{
  std::string name = _strings.pop();
  Type* t = _types.pop();
  OperatorType* type = constructOperatorType(t);

  unsigned arity = type->arity();
  bool isPredicate = type->isPredicateType();

  unsigned functor = isPredicate
                  ? env.signature->addFreshPredicate(arity, name.c_str())
                  : env.signature->addFreshFunction(arity,  name.c_str());
  Signature::Symbol *symbol = isPredicate
    ? env.signature->getPredicate(functor)
    : env.signature->getFunction(functor);

  symbol->setType(type);
  symbol->markLetBound();

  LetSymbolName symbolName(name, arity);
  LetSymbolReference symbolReference(functor, isPredicate);

  LetSymbols scope = _letTypedSymbols.pop();

  if (findLetSymbol(symbolName, scope, symbolReference)) {
    USER_ERROR("The symbol " + name + " of arity " + Int::toString(arity) + " is defined twice in a $let-expression.");
  }

  scope.push(LetSymbol(symbolName, symbolReference));
  _letTypedSymbols.push(scope);

  bool multipleLetTypes = _bools.pop();
  if (multipleLetTypes && getTok(0).tag == T_COMMA) {
    resetToks();
    _bools.push(multipleLetTypes);
    _states.push(END_LET_TYPES);
    _states.push(LET_TYPE);
  } 
} // TPTP::endLetTypes

void TPTP::definition()
{
  // At this point we parse one or more simultaneous definitions.
  // Simultaneous definitions are of the form `[D1, ..., Dn]` and each
  // definition is either of a function/predicate symbol `f(X,Y,Z) := t`
  // or a tuple `[a, b, c] := t`.

  // If the next token is '[', then the definition could either be
  // a single tuple definition or a list of simultaneous definitions.
  // This ambiguity is resolved by the next two tokens:
  // if they are T_NAME and then T_COMMA, then it is a tuple definition,
  // otherwise it is a a simultaneous definition.

  // The challenge here is how to direct the parser while keeping it
  // looking only one token ahead, like in the rest of TPTP.cpp.
  // Essentially, the trick is to
  //   1) have a boolean flag in _bools that tells whether the current let
  //      definition is simultaneous or not;
  //   2) consume the T_NAME token of a symbol definition here and not in
  //      SYMBOL_DEFINITION;
  //   3) consume the sequence T_LBRA, T_NAME, T_COMMA of tokens here and
  //      not in TUPLE_DEFINITION.

  switch (getTok(0).tag) {
    case T_NAME:
      _bools.push(false); // not a simultaneous definition
      _strings.push(name());
      _states.push(SYMBOL_DEFINITION);
      return;

    case T_LBRA: {
      resetToks();
      switch (getTok(0).tag) {
        case T_NAME:
          _strings.push(name());
          switch (getTok(0).tag) {
            case T_ASS:
            case T_LPAR:
              _bools.push(true); // is a simultaneous definition
              addTagState(T_RBRA);
              _states.push(SYMBOL_DEFINITION);
              return;

            case T_COMMA:
              resetToks();
              _bools.push(false); // not a simultaneous definition
              _states.push(TUPLE_DEFINITION);
              return;

            default:
              PARSE_ERROR(toString(T_ASS) + " or " + toString(T_LPAR) + " or " + toString(T_COMMA) + " expected",
                          getTok(0));
          }
          return;

        case T_LBRA:
          resetToks();
          _bools.push(true); // is a simultaneous definition
          addTagState(T_RBRA);
          _states.push(TUPLE_DEFINITION);
          return;

        default:
          PARSE_ERROR("name or " + toString(T_LBRA) + " expected",getTok(0));
      }
    }

    default:
      PARSE_ERROR("name or " + toString(T_LBRA) + " expected",getTok(0));
  }
} // TPTP::definition

void TPTP::midDefinition()
{
  switch (getTok(0).tag) {
    case T_NAME:
      _strings.push(name());
      _states.push(SYMBOL_DEFINITION);
      break;

    case T_LBRA:
      resetToks();
      _states.push(TUPLE_DEFINITION);
      break;

    default:
      PARSE_ERROR("name or " + toString(T_LBRA) + " expected",getTok(0));
  }
} // TPTP::midDefinition

void TPTP::symbolDefinition()
{
  std::string nm = _strings.pop();
  unsigned arity = 0;
  VList* vs = VList::empty();

  Stack<unsigned> vars;
  if (getTok(0).tag == T_LPAR) {
    resetToks();
    for (;;) {
      if (getTok(0).tag == T_VAR) {
        int var = _vars.insert(getTok(0).content);
        vars.push(var);
        resetToks();
      } else {
        PARSE_ERROR("variable expected", getTok(0));
      }

      if (getTok(0).tag == T_COMMA) {
        resetToks();
        continue;
      }

      if (getTok(0).tag == T_RPAR) {
        resetToks();
        break;
      }

      PARSE_ERROR("comma or closing bracket expected", getTok(0));
    }
    arity = (unsigned)vars.size();
  }

  LetSymbolName name(nm, arity);
  LetSymbolReference ref;
  if (!findLetSymbol(name, _letTypedSymbols.top(), ref)) {
    USER_ERROR("Symbol " + nm + " with arity " + Int::toString(arity) + " is used in a let definition without a declared type");
  }

  unsigned symbol = SYMBOL(ref);
  bool isPredicate = IS_PREDICATE(ref);

  if (arity > 0) {
    OperatorType* type = isPredicate
                       ? env.signature->getPredicate(symbol)->predType()
                       : env.signature->getFunction(symbol)->fnType();

    unsigned index = 0;
    while (vars.isNonEmpty()) {
      unsigned var = vars.pop();
      TermList sort = type->arg(arity - 1 - index++);
      bindVariable(var, sort);
      VList::push(var, vs);
    }

    _bindLists.push(vs);
    _states.push(UNBIND_VARIABLES);
  }

  LetDefinitions definitions = _letDefinitions.pop();
  definitions.push(LetSymbolReference(symbol, isPredicate));
  _letDefinitions.push(definitions);

  _varLists.push(vs);

  _states.push(END_DEFINITION);
  consumeToken(T_ASS);
  _states.push(TERM);
} // TPTP::symbolDefinition 

/**
 * Read a non-empty sequence of constants and save the resulting
 * sequence of TermList and their number
 * @since 20/04/2016 Gothenburg
 */
void TPTP::tupleDefinition()
{
  Set<std::string> uniqueConstants;
  Stack<unsigned> symbols;
  TermStack sorts;

  std::string constant = _strings.pop();
  do {
    if (uniqueConstants.contains(constant)) {
      USER_ERROR("The symbol " + constant + " is defined twice in a tuple $let-expression.");
    } else {
      uniqueConstants.insert(constant);
    }

    LetSymbolName constantName(constant, 0);
    LetSymbolReference ref;
    if (!findLetSymbol(constantName, _letTypedSymbols.top(), ref)) {
      USER_ERROR("Constant " + constant + " is used in a tuple let definition without a declared sort");
    }

    unsigned symbol = SYMBOL(ref);
    bool isPredicate = IS_PREDICATE(ref);

    symbols.push(symbol);

    TermList sort = isPredicate
                  ? AtomicSort::boolSort()
                  : env.signature->getFunction(symbol)->fnType()->result();
    sorts.push(sort);

    if (getTok(0).tag == T_NAME) {
      constant = name();
      if (getTok(0).tag == T_COMMA) {
        resetToks();
      }
    } else {
      break;
    }
  } while (true);

  TermList tupleSort = AtomicSort::tupleSort(sorts.size(), sorts.begin());
  unsigned tupleFunctor = Theory::tuples()->getFunctor(tupleSort);

  LetDefinitions definitions = _letDefinitions.pop();
  definitions.push(LetSymbolReference(tupleFunctor, false));
  _letDefinitions.push(definitions);

  VList* constants = VList::empty();
  VList::pushFromIterator(Stack<unsigned>::Iterator(symbols), constants);
  _varLists.push(constants);

  _states.push(END_DEFINITION);
  _states.push(TERM);
  addTagState(T_ASS);
  addTagState(T_RBRA);
} // tupleDefinition

void TPTP::endDefinition() {
  LetSymbolReference ref = _letDefinitions.top().top();
  unsigned symbol = SYMBOL(ref);
  bool isPredicate = IS_PREDICATE(ref);

  TermList definition = _termLists.top();
  TermList definitionSort = sortOf(definition);

  TermList refSort = isPredicate
                     ? AtomicSort::boolSort()
                     : env.signature->getFunction(symbol)->fnType()->result();

  if (refSort != definitionSort) {
    std::string definitionSortName = definitionSort.toString();
    std::string refSymbolName = isPredicate
                            ? env.signature->predicateName(symbol)
                            : env.signature->functionName(symbol);
    OperatorType* type = isPredicate
                         ? env.signature->getPredicate(symbol)->predType()
                         : env.signature->getFunction(symbol)->fnType();
    USER_ERROR("The term " + definition.toString() + " of the sort " + definitionSortName +
               " is used as definition of the symbol " + refSymbolName +
               " of the type " + type->toString());
  }

  bool multipleDefinitions = _bools.pop();
  if (multipleDefinitions && getTok(0).tag == T_COMMA) {
    resetToks();
    _bools.push(multipleDefinitions);
    _states.push(MID_DEFINITION);
  } else {
    _letSymbols.push(_letTypedSymbols.pop());
  }
} // endDefinition

bool TPTP::findLetSymbol(LetSymbolName symbolName, LetSymbolReference& symbolReference) {
  Stack<LetSymbols>::TopFirstIterator scopes(_letSymbols);
  while (scopes.hasNext()) {
    LetSymbols scope = scopes.next();
    if (findLetSymbol(symbolName, scope, symbolReference)) {
      return true;
    }
  }
  return false;
} // findLetSymbol(LetSymbolName,LetSymbolReference)

bool TPTP::findLetSymbol(LetSymbolName symbolName, LetSymbols scope, LetSymbolReference& symbolReference) {
  LetSymbols::Iterator symbols(scope);
  while (symbols.hasNext()) {
    LetSymbol symbol = symbols.next();
    if (symbol.first == symbolName) {
      symbolReference = symbol.second;
      return true;
    }
  }
  return false;
} // findLetSymbol(LetSymbolName,LetSymbols,LetSymbolReference)


/**
 * Process the end of the $let expression
 * @since 27/07/2011 Manchester
 */
void TPTP::endLet()
{
  TermList let = _termLists.pop();
  TermList sort = sortOf(let);

  _letSymbols.pop();
  LetDefinitions scope = _letDefinitions.pop(); // TODO: inlining this crashes the program, WTF?
  LetDefinitions::TopFirstIterator definitions(scope);
  while (definitions.hasNext()) {
    LetSymbolReference ref = definitions.next();
    unsigned symbol = SYMBOL(ref);
    bool isPredicate = IS_PREDICATE(ref);

    VList* varList = _varLists.pop();
    TermList definition = _termLists.pop();

    bool isTuple = false;
    if (!isPredicate) {
      TermList resultSort = env.signature->getFunction(symbol)->fnType()->result();
      isTuple = resultSort.isTupleSort();
    }

    if (isTuple) {
      let = TermList(Term::createTupleLet(symbol, varList, definition, let, sort));
    } else {
      let = TermList(Term::createLet(symbol, varList, definition, let, sort));
    }
  }
  _termLists.push(let);
} // endLet

/**
 * Process the end of the tuple expression
 * @since 19/04/2016 Gothenburg
 */
void TPTP::endTuple()
{
  unsigned arity = (unsigned)_ints.pop();
  ASS_GE(_termLists.size(), arity);

  DArray<TermList> elements(arity);
  DArray<TermList> sorts(arity);

  for (int i = arity - 1; i >= 0; i--) {
    TermList ts = _termLists.pop();
    elements[i] = ts;
    sorts[i] = sortOf(ts);
  }

  Term* t = Term::createTuple(arity, sorts.begin(), elements.begin());
  _termLists.push(TermList(t));
} // endTuple

/**
 * Read a non-empty sequence of arguments, including the right parentheses
 * and save the resulting sequence of TermList and their number
 * @since 10/04/2011 Manchester
 */
void TPTP::args()
{
  _states.push(END_ARGS);
  _states.push(TERM);
} // args

/**
 * Read a list of arguments after a term
 * @since 27/07/2011 Manchester
 */
void TPTP::endArgs()
{
 // check if there is any other term in the argument list
  Token tok = getTok(0);
  switch (tok.tag) {
  case T_COMMA:
    resetToks();
    _ints.push(_ints.pop()+1);
    _states.push(END_ARGS);
    _states.push(TERM);
    return;
  case T_RPAR:
    resetToks();
    return;
  case T_RBRA:
    resetToks();
    return;
  default:
    PARSE_ERROR(", ) or ] expected after an end of a term",tok);
  }
} // endArgs

/**
 * Bind a variable to a sort
 * @since 22/04/2011 Manchester
 */
void TPTP::bindVariable(unsigned var,TermList sort)
{
  SList** definitions;
  // definitions will be a pointer to the list inside _variableSorts,
  // either the one that was there, or a freshly inserted empty one
  _variableSorts.getValuePtr(var,definitions,SList::empty());
  SList::push(sort,*definitions); // and this will modify that list
} // bindVariable

/**
 * Read a non-empty sequence of variable and save the resulting
 * sequence of TermList and their number
 * @since 07/07/2011 Manchester
 * @since 16/04/2015 Gothenburg, do not parse the closing ']'
 */
void TPTP::varList()
{
  Stack<int> vars;
  for (;;) {
    Token& tok = getTok(0);

    if (tok.tag != T_VAR) {
      PARSE_ERROR("variable expected",tok);
    }
    int var = _vars.insert(tok.content);
    if (_isQuestion) {
      _curQuestionVarNames.insert(var,tok.content);
    }
    vars.push(var);
    resetToks();
    bool sortDeclared = false;
  afterVar:
    tok = getTok(0);
    switch (tok.tag) {
    case T_COLON: // v: type
      if (sortDeclared) {
        PARSE_ERROR("two declarations of variable sort",tok);
      }
      resetToks();
      bindVariable(var,(_isThf ? readArrowSort() : readSort()));
      sortDeclared = true;
      goto afterVar;

    case T_COMMA:
      if (!sortDeclared) {
        bindVariable(var,AtomicSort::defaultSort());
      }
      resetToks();
      break;

    default:
      {
        if (!sortDeclared) {
          bindVariable(var,AtomicSort::defaultSort());
        }
        VList* vs = VList::empty();
        SList* ss = SList::empty();
        while (!vars.isEmpty()) {
          int v = vars.pop();
          VList::push(v,vs);
          SList::push(sortOf(TermList(v,false)),ss);
        }
        _varLists.push(vs);
        _sortLists.push(ss);
        _bindLists.push(vs);
        return;
      }
    }
  }
} // varList

/**
 * Read a term and save the resulting TermList
 * @since 10/04/2011 Manchester
 * @since 13/04/2015 Gothenburg, major changes to support FOOL
 */
void TPTP::term()
{
  Token tok = getTok(0);
  switch (tok.tag) {
    case T_NAME:
    case T_THEORY_FUNCTION:
    case T_VAR:
    case T_ITE:
    case T_LET:
    case T_LBRA:
      _states.push(TERM_INFIX);
      _states.push(FUN_APP);
      return;

    
    case T_INTEGER_TYPE:
    case T_REAL_TYPE:
    case T_RATIONAL_TYPE: 
    case T_BOOL_TYPE:
    case T_DEFAULT_TYPE: {
      resetToks();
      switch (tok.tag) {
        case T_INTEGER_TYPE:
          _termLists.push(AtomicSort::intSort());
          break;
        case T_REAL_TYPE:
          _termLists.push(AtomicSort::realSort());
          break;        
        case T_RATIONAL_TYPE:
          _termLists.push(AtomicSort::rationalSort());
          break;
        case T_BOOL_TYPE:
          _termLists.push(AtomicSort::boolSort());
          break;
        case T_DEFAULT_TYPE:
          _termLists.push(AtomicSort::defaultSort());
          break;             
        default:
          ASSERTION_VIOLATION;
      }
      return;
    }

    case T_STRING:
    case T_INT:
    case T_REAL:
    case T_RAT: {
      resetToks();
      unsigned number;
      switch (tok.tag) {
        case T_STRING:
          number = env.signature->addStringConstant(tok.content);
          // "distinct_object"s are _always_ of sort $i, even in typed contexts
          env.signature->getFunction(number)->setType(
            OperatorType::getConstantsType(AtomicSort::defaultSort())
          );
          break;
        case T_INT:
          number = addIntegerConstant(tok.content,_overflow,_isFof);
          break;
        case T_REAL:
          number = addRealConstant(tok.content,_overflow,_isFof);
          break;
        case T_RAT:
          number = addRationalConstant(tok.content,_overflow,_isFof);
          break;
        default:
          ASSERTION_VIOLATION;
      }
      _termLists.push(TermList(Term::createConstant(number)));
      return;
    }

    default:
      _states.push(FORMULA_INSIDE_TERM);
  }
} // term

/**
 * Build a term assembled by term()
 * @since 09/07/2011 Manchester
 * @since 14/04/2015 Gothenburg, major changes to support FOOL
 */
void TPTP::endTerm()
{
  std::string name = _strings.pop();

  if (name == toString(T_ITE)) {
    _states.push(END_ITE);
    return;
  }

  if (name == toString(T_LET)) {
    _states.push(END_LET);
    return;
  }

  if (name == toString(T_TUPLE)) {
    _states.push(END_TUPLE);
    return;
  }

  TheoryFunction tf;
  if (findTheoryFunction(name, tf)) {
    _theoryFunctions.push(tf);
    _states.push(END_THEORY_FUNCTION);
    return;
  }

  int arity = _ints.pop();

  if (arity == -1) {
    // it was a variable
    unsigned var = (unsigned)_vars.insert(name);
    _termLists.push(TermList(var, false));
    return;
  }

  LetSymbolReference ref;
  if (env.signature->predicateExists(name, arity) ||
      (findLetSymbol(LetSymbolName(name, arity), ref) && IS_PREDICATE(ref)) ||
      findInterpretedPredicate(name, arity)) {
    // if the function symbol is actually a predicate,
    // we need to construct a formula and wrap it inside a term
    _formulas.push(createPredicateApplication(name, arity));
    _states.push(END_FORMULA_INSIDE_TERM);
    return;
  }

  if(env.signature->typeConExists(name, arity)){
    _termLists.push(createTypeConApplication(name, arity));    
    return;
  }

  _termLists.push(createFunctionApplication(name, arity)); 
} // endTerm

/**
 * Read after an end of atom or after lhs of an equality or inequality
 * @since 10/04/2011 Manchester
 * @since 13/04/2015 Gothenburg, major changes to support FOOL
 */
void TPTP::formulaInfix()
{
  Token tok = getTok(0);

  if (tok.tag == T_EQUAL || tok.tag == T_NEQ) {
    _states.push(END_EQ);
    _states.push(TERM);
    _states.push(MID_EQ);
    _states.push(END_TERM);
    return;
  }

  std::string name = _strings.pop();

  if (name == toString(T_ITE)) {
    _states.push(END_TERM_AS_FORMULA);
    _states.push(END_ITE);
    return;
  }

  TheoryFunction tf;
  if (findTheoryFunction(name, tf)) {
    switch (tf) {
      case TF_STORE:
        USER_ERROR("$store expression cannot be used as formula");
        break;
      case TF_SELECT:
        _theoryFunctions.push(tf);
        _states.push(END_TERM_AS_FORMULA);
        _states.push(END_THEORY_FUNCTION);
        break;
      default:
        ASSERTION_VIOLATION_REP(name);
    }
    return;
  }

  if (name == toString(T_LET)) {
    _states.push(END_TERM_AS_FORMULA);
    _states.push(END_LET);
    return;
  }

  int arity = _ints.pop();

  if (arity == -1) {
    // that was a variable
    unsigned var = (unsigned)_vars.insert(name);
    _termLists.push(TermList(var, false));
    _states.push(END_TERM_AS_FORMULA);
    return;
  }

  if(env.signature->functionExists(name, arity)){
    //with polymorphism, we can have function symbols that are used as predicates
    _termLists.push(createFunctionApplication(name, arity));
    _states.push(END_TERM_AS_FORMULA);
    return;
  }

  _formulas.push(createPredicateApplication(name, arity));
} // formulaInfix

/**
 * Read after an end of equality or inequality and save the (in)equality formula.
 * @since 09/07/2011 Manchester
 */
void TPTP::endEquality()
{
  _insideEqualityArgument--;

  if((_isThf) && (_lastPushed == FORM)){
    endFormulaInsideTerm();
  }

  TermList rhs = _termLists.pop();
  TermList lhs = _termLists.pop();

  if (sortOf(rhs) != sortOf(lhs)) {
    TermList rsort = sortOf(rhs); 
    TermList lsort = sortOf(lhs);
    USER_ERROR("Cannot create equality between terms of different types.\n"+
      rhs.toString()+" is "+rsort.toString()+"\n"+
      lhs.toString()+" is "+lsort.toString()
    );
  }

  Literal* l = createEquality(_bools.pop(),lhs,rhs);
  _formulas.push(new AtomicFormula(l));
  _lastPushed = FORM;
} // endEquality

/**
 * Read
 * @since 09/07/2011 Manchester
 */
void TPTP::midEquality()
{
  _insideEqualityArgument++;

  Token tok = getTok(0);
  switch (tok.tag) {
  case T_EQUAL:
    _bools.push(true);
    break;
  case T_NEQ:
    _bools.push(false);
    break;
  default:
    PARSE_ERROR("either = or != expected",tok);
  }
  resetToks();
} // midEquality

/**
 * Creates an equality literal and takes care of its sort when it
 * is an equality between two variables.
 * @since 21/07/2011 Manchester
 * @since 03/05/2013 Train Manchester-London, bug fix
 */
Literal* TPTP::createEquality(bool polarity,TermList& lhs,TermList& rhs)
{
  TermList masterVar;
  TermList sort;
  if (!SortHelper::getResultSortOrMasterVariable(lhs, sort, masterVar)) {
    // Master variable is a variable whose sort determines the sort of a term.
    // If term is a variable, the master variable is the variable itself. The
    // trickier case is when we have an if-then-else expression with variable
    // arguments.
    SList* vs;
    if (_variableSorts.find(masterVar.var(),vs) && vs) {
      sort = vs->head();
    }
    else { // this may happen when free variables appear in the formula (or clause)
      sort = AtomicSort::defaultSort();
    }
  }
   
  return Literal::createEquality(polarity,lhs,rhs,sort);
} // TPTP::createEquality

/**
 * Creates a formula that is a predicate application literal from
 * provided predicate symbol name and arity. If arity is greater than zero,
 * the arguments are assumed to be on the _termLists stack.
 * @since 27/03/1015 Manchester
 */
Formula* TPTP::createPredicateApplication(std::string name, unsigned arity)
{
  ASS_GE(_termLists.size(), arity);

  int pred;
  LetSymbolReference ref;
  if (findLetSymbol(LetSymbolName(name, arity), ref) && IS_PREDICATE(ref)) {
    pred = (int)SYMBOL(ref);
  } else {
    if (arity > 0) {
      bool dummy;
      pred = addPredicate(name, arity, dummy, _termLists.top());
    } else {
      pred = env.signature->addPredicate(name, 0);
    }
  }
  if (pred == -1) { // equality
    TermList rhs = _termLists.pop();
    TermList lhs = _termLists.pop();
    return new AtomicFormula(createEquality(true,lhs,rhs));//TODO equality sort?
  }
  if (pred == -2){ // distinct
    // TODO check that we are top-level
    // If fewer than 5 things are distinct then we add the disequalities
    if(arity<5){
      static Stack<unsigned> distincts;
      distincts.reset();
      for(int i=arity-1;i >= 0; i--){
        TermList t = _termLists.pop();
        if(t.isVar() || t.term()->arity()!=0){
          USER_ERROR("$distinct can only be used with constants. Found "+t.toString());
        }
        distincts.push(t.term()->functor());
      }
      Formula* distinct_formula = DistinctGroupExpansion(0 /* zero means "always expand"*/).expand(distincts);
      return distinct_formula;
    }else{
      // Otherwise record them as being in a distinct group
      unsigned grpIdx = env.signature->createDistinctGroup(0);
      for(int i = arity-1;i >=0; i--){
        TermList ts = _termLists.pop();
        if(!ts.isTerm() || ts.term()->arity()!=0){
          USER_ERROR("$distinct can only be used with constants. Found "+ts.toString());
        }
        env.signature->addToDistinctGroup(ts.term()->functor(),grpIdx);
      }
      return new Formula(true); // we ignore it, it evaluates to true as we have recorded it elsewhere
    }
  }
  // not equality or distinct
  auto args = nLastTermLists(arity);
  OperatorType* type = env.signature->getPredicate(pred)->predType();
  for (auto i : range(0, arity)) {
    TermList sort = type->arg(i);
    TermList ts = args[i];
    TermList tsSort = sortOf(ts);
    if(i < type->numTypeArguments()){
      if(tsSort != AtomicSort::superSort()){
        USER_ERROR("The sort ", tsSort, " of type argument ", ts, " is not $ttype as mandated by TF1");
      }
    } else {
      _substScratchpad.reset();
      if(!_substScratchpad.match(sort, 0, tsSort, 1)) {
        USER_ERROR("Failed to create predicate application for ", name, " of type ", type->toString(), "\n",
                   "The sort ", tsSort, " of the intended term argument ", ts, " (at index ", i, ") "
                   "is not an instance of sort ", sort);
      }
    }
  }
  auto out = new AtomicFormula(Literal::create(pred, arity, /* polarity */ true, args));
  _termLists.pop(arity);
  return out;
} // createPredicateApplication


/**
 * Creates a term that is a function application from
 * provided function symbol name and arity. If arity is greater than zero,
 * the arguments are assumed to be on the _termLists stack.
 * @since 13/04/2015 Gothenburg, major changes to support FOOL
 */
TermList TPTP::createFunctionApplication(std::string name, unsigned arity)
{ //TODO update to deal with wierd /\ @ ... syntax
  ASS_GE(_termLists.size(), arity);

  unsigned fun;
  LetSymbolReference ref;
  if (findLetSymbol(LetSymbolName(name, arity), ref) && !IS_PREDICATE(ref)) {
    fun = SYMBOL(ref);
  } else {
    bool dummy;
    if (arity > 0) {
      fun = addFunction(name, arity, dummy, _termLists.top());
    } else {
      fun = addUninterpretedConstant(name, _overflow, dummy);
    }
  }

  OperatorType* type = env.signature->getFunction(fun)->fnType();
  auto args = nLastTermLists(arity);
  for (unsigned i : range(0, arity)) {
    TermList sort = type->arg(i);
    TermList ss = args[i];
    TermList ssSort = sortOf(ss);
    if(i < type->numTypeArguments()){
      if(ssSort != AtomicSort::superSort()){
        USER_ERROR("The sort " + ssSort.toString() + " of type argument " + ss.toString() + " "
                   "is not $tType as mandated by TF1");
      }
    } else {
      _substScratchpad.reset();
      if(!_substScratchpad.match(sort, 0, ssSort, 1)){
        USER_ERROR("Failed to create function application for " + name + " of type " + type->toString() + "\n" +
                   "The sort " + ssSort.toString() + " of the intended term argument " + ss.toString() + " (at index " + Int::toString(i) +") "
                   "is not an instance of sort " + sort.toString());
      }
    }
  }
  auto t = TermList(Term::create(fun, arity, args));
  _termLists.pop(arity);
  return t;
}

/**
 * Creates a term that is a function application from
 * provided function symbol name and arity. If arity is greater than zero,
 * the arguments are assumed to be on the _termLists stack.
 * @since 13/04/2015 Gothenburg, major changes to support FOOL
 */
TermList TPTP::createTypeConApplication(std::string name, unsigned arity)
{ 
  ASS_GE(_termLists.size(), arity);

  bool added = false;
  //TODO not checking for overflown constant. Is that OK?
  //seems to be done this way for predicates as well.
  unsigned typeCon = env.signature->addTypeCon(name,arity,added);
  if(added)
    USER_ERROR("Undeclared type constructor ", name, "/", arity);

  auto args = nLastTermLists(arity);
  for (auto i : range(0, arity)) {
    auto term = args[i];
    auto sort = sortOf(term);
    if (sort != AtomicSort::superSort()) 
        USER_ERROR("The sort ", sort, " of type argument ", term, " is not $tType as mandated by TF1");
  }
  auto s = TermList(AtomicSort::create(typeCon, arity, args));
  _termLists.pop(arity);
  return s;
}

/**
 * Build a formula from previousy built subformulas
 * @since 10/04/2011 Manchester
 */
void TPTP::endFormula()
{
  int con = _connectives.pop();
  Formula* f;
  bool conReverse = false;
  switch (con) {
  case IMP:
  case AND:
  case OR:
    conReverse = _bools.pop();
    break;
  case IFF:
  case XOR:
  case -1:
    break;
  case NOT:
    f = _formulas.pop();
    // This gets rid of the annoying step in proof output where ~(L) is flattend to (~L)
    if(f->connective()==LITERAL){
      Literal* oldLit = static_cast<AtomicFormula*>(f)->literal();
      Literal* newLit = Literal::create(oldLit,!oldLit->polarity());
      _formulas.push(new AtomicFormula(newLit));
    }
    else{
      _formulas.push(new NegatedFormula(f));
    }
    _states.push(END_FORMULA);
    return;
  case FORALL:
  case EXISTS:
    f = _formulas.pop();
    _formulas.push(new QuantifiedFormula((Connective)con,_varLists.pop(),_sortLists.pop(),f));
    _states.push(END_FORMULA);
    return;
  case LITERAL:
  default:
    throw ::Exception((std::string)"tell me how to handle connective " + Int::toString(con));
  }

  Token& tok = getTok(0);
  Tag tag = tok.tag;
  Connective c;
  bool cReverse = false;
  switch (tag) {
  case T_AND:
    c = AND;
    break;
  case T_NOT_AND:
    cReverse = true;
    c = AND;
    break;
  case T_NOT_OR:
    cReverse = true;
    c = OR;
    break;
  case T_OR:
    c = OR;
    break;
  case T_XOR:
    c = XOR;
    break;
  case T_IFF:
    c = IFF;
    break;
  case T_IMPLY:
    c = IMP;
    break;
  case T_REVERSE_IMP:
    cReverse = true;
    c = IMP;
    break;
  case T_EQUAL:
  case T_NEQ: {
    // not connectives, but we allow formulas to be arguments to = and !=
    _states.push(END_EQ);
    _states.push(TERM);
    _states.push(MID_EQ);
    _states.push(END_FORMULA_INSIDE_TERM);
    return;
  }
  default:
    // the formula does not end at a binary connective, build the formula and terminate
    switch (con) {
    case IMP:
      f = _formulas.pop();
      if (conReverse) {
	f = new BinaryFormula((Connective)con,f,_formulas.pop());
      }
      else {
	f = new BinaryFormula((Connective)con,_formulas.pop(),f);
      }
      _formulas.push(f);
      _states.push(END_FORMULA);
      return;

    case IFF:
    case XOR:
      f = _formulas.pop();
      f = new BinaryFormula((Connective)con,_formulas.pop(),f);
      _formulas.push(f);
      _states.push(END_FORMULA);
      return;

    case AND:
    case OR:
      f = _formulas.pop();
      f = makeJunction((Connective)con,_formulas.pop(),f);
      if (conReverse) {
	f = new NegatedFormula(f);
      }
      _formulas.push(f);
      _states.push(END_FORMULA);
      return;

    case -1:
      return;
    default:
      ASSERTION_VIOLATION;
    }
  }

  // con and c are binary connectives
  if (higherPrecedence(con,c)) {
    f = _formulas.pop();
    Formula* g = _formulas.pop();
    if (con == AND || con == OR) {
      f = makeJunction((Connective)con,g,f);
      if (conReverse) {
	f = new NegatedFormula(f);
      }
    }
    else if (con == IMP && conReverse) {
      f = new BinaryFormula((Connective)con,f,g);
    }
    else {
      f = new BinaryFormula((Connective)con,g,f);
    }
    _formulas.push(f);
    _states.push(END_FORMULA);
    return;
  }

  // c is a binary connective
  _connectives.push(con);
  if (con == IMP || con == AND || con == OR) {
    _bools.push(conReverse);
  }
  _connectives.push(c);
  if (c == IMP || c == AND || c == OR) {
    _bools.push(cReverse);
  }
  resetToks();
  _states.push(END_FORMULA);
  _states.push(SIMPLE_FORMULA);
} // endFormula

/**
 * Builds a term that really is a formula
 * @author Evgeny Kotelnikov
 * @since 27/03/2015 Manchester
 */
void TPTP::formulaInsideTerm()
{
  _states.push(END_FORMULA_INSIDE_TERM);
  _states.push(FORMULA);
} // formulaInsideTerm

/**
 * Wraps a formula inside a term
 * @author Evgeny Kotelnikov
 * @since 27/03/2015 Manchester
 */
void TPTP::endFormulaInsideTerm()
{
  Formula* f = _formulas.pop();
  TermList ts(Term::createFormula(f));
  _termLists.push(ts);
  _lastPushed = TM;
} // endFormulaInsideTerm

/**
 * Makes a boolean term a formula
 * @author Evgeny Kotelnikov
 * @since 27/03/2015 Manchester
 */
void TPTP::endTermAsFormula()
{
  TermList t = _termLists.pop();
  TermList tSort = sortOf(t);
  if (tSort != AtomicSort::boolSort()) {
    USER_ERROR("Non-boolean term " + t.toString() + " of sort " + tSort.toString() + " is used in a formula context");
  }
  if (t.isTerm() && t.term()->isFormula()) {
    _formulas.push(t.term()->getSpecialData()->getFormula());
    _lastPushed = FORM;
  } else {
    _formulas.push(new BoolTermFormula(t));
    _lastPushed = FORM;
  }
} // endTermAsFormula

/**
 * Build a type from previousy built types
 * @since 14/07/2011 Manchester
 */
void TPTP::endType()
{
  TypeTag tt = _typeTags.pop();
  Type* t = _types.pop();
  switch (tt) {
  case TT_ATOMIC:
    break;
  case TT_PRODUCT:
    t = new ProductType(_types.pop(),t);
    tt = _typeTags.pop();
    break;
  case TT_ARROW:
    t = new ArrowType(_types.pop(),t);
    tt = _typeTags.pop();
    break;
  case TT_QUANTIFIED:
    VList* vl = _varLists.pop();
    _sortLists.pop();
    t = new QuantifiedType(t, vl);
    tt = _typeTags.pop();
    break;    
  }
  ASS(tt == TT_ATOMIC);
  _types.push(t);

  Token tok = getTok(0);
  switch (tok.tag) {
  case T_STAR:
    _typeTags.push(tt);
    _typeTags.push(TT_PRODUCT);
    break;
  case T_ARROW:
    _typeTags.push(tt);
    _typeTags.push(TT_ARROW);
    break;
  default:
    return;
  }
  resetToks();
  _states.push(END_TYPE);
  _states.push(SIMPLE_TYPE);
} // endType

/**
 * Skip a tag.
 * @since 10/04/2011 Manchester
 */
void TPTP::tag()
{
  consumeToken(_tags.pop());
} // tag

/**
 * Process the end of the fof() definition and build the corresponding unit.
 * @since 10/04/2011 Manchester
 */
void TPTP::endFof()
{
  TPTP::SourceRecord* source = 0;

  // are we interested in collecting sources?
  if (_unitSources) {
    source = getSource();
  }
#if DEBUG_SOURCE
  else{
    // create fake map
    _unitSources = new DHMap<Unit*,SourceRecord*>();
    source = getSource();
  }
#endif

  skipToRPAR();
  consumeToken(T_DOT);

  bool isFof = _bools.pop();
  Formula* f = _formulas.pop();
  std::string nm = _strings.pop(); // unit name
  if (_allowedNames && !_allowedNames->contains(nm)) {
    return;
  }

  Unit* unit;
  if (isFof) { // fof() or tff()
    env.statistics->inputFormulas++;
    unit = new FormulaUnit(f,FromInput(_lastInputType));
    unit->setInheritedColor(_currentColor);
  }
  else { // cnf()
    env.statistics->inputClauses++;
    // convert the input formula f to a clause
    Stack<Formula*> forms;
    Stack<Literal*> lits;
    Formula* g = nullptr;
    forms.push(f);
    while (! forms.isEmpty()) {
      g = forms.pop();
      switch (g->connective()) {
      case OR:
	{
	  FormulaList::Iterator fs(static_cast<JunctionFormula*>(g)->getArgs());
	  while (fs.hasNext()) {
	    forms.push(fs.next());
	  }
	}
	break;

      case LITERAL:
      case NOT:
	{
	  bool positive = true;
	  while (g->connective() == NOT) {
	    g = static_cast<NegatedFormula*>(g)->subformula();
	    positive = !positive;
	  }
	  if (g->connective() != LITERAL) {
	    USER_ERROR((std::string)"input formula not in CNF: " + f->toString());
	  }
	  Literal* l = static_cast<AtomicFormula*>(g)->literal();
	  lits.push(positive ? l : Literal::complementaryLiteral(l));
	}
	break;

      case TRUE:
	return;
      case FALSE:
	break;
      default:
	USER_ERROR((std::string)"input formula not in CNF: " + f->toString());
      }
    }
    unit = Clause::fromStack(lits,FromInput(_lastInputType));
    unit->setInheritedColor(_currentColor);
  }

  if(source){ 
    ASS(_unitSources);
    _unitSources->insert(unit,source);
  }

  if (env.options->outputAxiomNames() || env.options->proof() == Options::Proof::DEDUKTI) {
    assignAxiomName(unit,nm);
  }
#if DEBUG_SHOW_UNITS
  cout << "Unit: " << unit->toString() << "\n";
#endif
  if (!_inputs.isEmpty()) {
    unit->inference().markIncluded();
  }

  switch (_lastInputType) {
  case UnitInputType::CONJECTURE:
    if(!isFof) USER_ERROR("conjecture is not allowed in cnf");
    if(_seenConjecture) USER_ERROR("Vampire only supports a single conjecture in a problem");
    _seenConjecture=true;
    if (_isQuestion && ((env.options->mode() == Options::Mode::CLAUSIFY) || (env.options->mode() == Options::Mode::TCLAUSIFY)) && f->connective() == EXISTS) {
      // create an answer predicate
      QuantifiedFormula* g = static_cast<QuantifiedFormula*>(f);
      unsigned arity = VList::length(g->vars());
      unsigned pred = env.signature->addPredicate("$$answer",arity);
      env.signature->getPredicate(pred)->markAnswerPredicate();
      Recycled<Stack<TermList>> args;
      VList::Iterator vs(g->vars());
      while (vs.hasNext()) {
        args->push(TermList::var(vs.next()));
      }
      Literal* a = Literal::create(pred, arity, /* polarity */ true, args->begin());
      f = new QuantifiedFormula(FORALL,
        g->vars(),
        g->sorts(),
        new BinaryFormula(IMP,g->subformula(),new AtomicFormula(a)));
        unit = new FormulaUnit(f,FormulaTransformation(InferenceRule::ANSWER_LITERAL_INJECTION,unit));
    }
    else {
      VList* vs = freeVariables(f);
      if (VList::isEmpty(vs)) {
        f = new NegatedFormula(f);
      }
      else {
        // TODO can we use sortOf to get the sorts of vs?
        f = new NegatedFormula(new QuantifiedFormula(FORALL,vs,0,f));
      }
      unit = new FormulaUnit(f,
			     FormulaTransformation(InferenceRule::NEGATED_CONJECTURE,unit));
      if (_isQuestion) {
        _questionVariableNames.insert(unit->number(),new Map<int,std::string>(std::move(_curQuestionVarNames)));
      }
    }
    break;

  case UnitInputType::CLAIM:
    unit = processClaimFormula(unit,f,nm);
    break;

  default:
    break;
  }
  _units.pushBack(unit);
} // tag

/*
* The given unit has already been parsed (and had the role CLAIM).
* It's actually a FormulaUnit with the formula f wrapped inside.
* nm is the name of the claim in TPTP, but can be any other string
* that will serve as the name of the predice we introduce:
*
* Now instead of returning it directly, we turn it into an equivalence
* with a fresh predicate symbol (of name nm) and return that one.
* The new symbo is marked not to be eliminated during preprocessing.
*/
Unit* TPTP::processClaimFormula(Unit* unit, Formula * f, const std::string& nm)
{
  bool added;
  unsigned pred = env.signature->addPredicate(nm,0,added);
  if (!added) {
    USER_ERROR("Names of claims must be unique: "+nm);
  }
  env.signature->getPredicate(pred)->markLabel();
  Formula* claim = new AtomicFormula(Literal::create(pred, /* polarity */ true, {}));
  VList* vs = freeVariables(f);
  if (VList::isNonEmpty(vs)) {
    //TODO can we use sortOf to get sorts of vs?
    f = new QuantifiedFormula(FORALL,vs,0,f);
  }
  f = new BinaryFormula(IFF,claim,f);
  return new FormulaUnit(f,
      FormulaTransformation(InferenceRule::CLAIM_DEFINITION,unit));
}

/**
 * Add a state just reading a tag and save the tag in _tags.
 * @since 28/07/2011 Manchester
 */
void TPTP::addTagState(Tag t)
{
  _states.push(TAG);
  _tags.push(t);
} // TPTP::addTagState

/**
 * Process the end of the tff() definition and build the corresponding unit.
 * @since 14/07/2011 Manchester
 */
void TPTP::endTff()
{
  int rpars= _ints.pop();
  while (rpars--) {
    consumeToken(T_RPAR);
  }
  skipToRPAR();
  consumeToken(T_DOT);

  // build a TPTP out of the parse type
  ASS(_typeTags.isEmpty());
  Type* t = _types.pop();
  ASS(_types.isEmpty());

  OperatorType* ot = constructOperatorType(t);
  std::string name = _strings.pop();

  unsigned arity = ot->arity();
  bool isPredicate = ot->isPredicateType() && !_isThf;
  bool isTypeCon = !isPredicate && (ot->result() == AtomicSort::superSort());

  bool added;
  Signature::Symbol* symbol;
  if (isPredicate) {
    unsigned pred = env.signature->addPredicate(name, arity, added);
    symbol = env.signature->getPredicate(pred);
    if (!added) {
      // GR: Multiple identical type declarations for a symbol are allowed
      if(symbol->predType() != ot){
        USER_ERROR("Predicate symbol type is declared after its use: " + name);
      }
    }
    else{
      if (arity != 0) {
        symbol->setType(ot);
      }
    }
  } else if (isTypeCon){
    unsigned typeCon = env.signature->addTypeCon(name, arity, added);
    symbol = env.signature->getTypeCon(typeCon);
    if (!added) {
      // GR: Multiple identical type declarations for a symbol are allowed
      if(symbol->typeConType() != ot){
        USER_ERROR("Type constructor type is declared after its use: " + name);
      }
    }
    else{
      symbol->setType(ot);
    }
  } else {
    unsigned fun = arity == 0
                   ? addUninterpretedConstant(name, _overflow, added)
                   : env.signature->addFunction(name, arity, added);
    symbol = env.signature->getFunction(fun);
    if (!added) {
      if(symbol->fnType() != ot){
        USER_ERROR("Function symbol type is declared after its use: " + name);
      }
    }
    else {   
      symbol->setType(ot);
      //TODO check whether the below is actually required or not.
      if(_isThf){
        if(!_typeArities.insert(name, ot->numTypeArguments())){
          USER_ERROR("Symbol " + name + " used with different type arities");
        }
      }
    }
    //cout << "added: " + symbol->name() + " of type " + ot->toString() + " and functor " << fun << endl;
  }
} // endTff


OperatorType* TPTP::constructOperatorType(Type* t, VList* vars)
{
  TermList resultSort;
  Stack<TermList> argumentSorts;

  switch (t->tag()) {
    case TT_PRODUCT:
      USER_ERROR("product types are not supported");

    case TT_ATOMIC: {
      // atomic types: 0-ary predicates (propositions) and constants (0-ary functions, eg. int constant, array1 constants)
      resultSort = static_cast<AtomicType*>(t)->sort();
      break;
    }

    case TT_ARROW: {
      // non-atomic types, i.e. with arrows
      ArrowType* at = static_cast<ArrowType*>(t);
      Type* rhs = at->returnType();
      if (rhs->tag() != TT_ATOMIC) {
        USER_ERROR("complex return types are not supported");
      }

      resultSort = static_cast<AtomicType*>(rhs)->sort();
      Stack<Type*> types;
      types.push(at->argumentType());
      while (!types.isEmpty()) {
        Type *tp = types.pop();
        switch (tp->tag()) {
          case TT_ARROW:
            USER_ERROR("higher-order types are not supported");

          case TT_ATOMIC: {
            TermList sort = static_cast<AtomicType*>(tp)->sort();
            argumentSorts.push(sort);
            break;
          }

          case TT_PRODUCT: {
            ProductType* pt = static_cast<ProductType*>(tp);
            types.push(pt->rhs());
            types.push(pt->lhs());
            break;
          }
          default:
            ASSERTION_VIOLATION;
        }
      }
      break;
    }

    case TT_QUANTIFIED: {
      QuantifiedType* qt = static_cast<QuantifiedType*>(t);
      OperatorType* ot = constructOperatorType(qt->qtype(), qt->vars());
      return ot;
      //TODO check that all free variables in ot are from quantifiedVars
    }

    default:
      ASSERTION_VIOLATION;
  }

  bool isPredicate = resultSort == AtomicSort::boolSort();
  unsigned arity = (unsigned)argumentSorts.size();

  if(_containsPolymorphism){
    SortHelper::normaliseArgSorts(vars, argumentSorts);
    SortHelper::normaliseSort(vars, resultSort);
  }

  if (isPredicate && !_isThf) { //in THF, we treat predicates and boolean terms the same
    return OperatorType::getPredicateType(arity, argumentSorts.begin(), VList::length(vars));
  } else {
    return OperatorType::getFunctionType(arity, argumentSorts.begin(), resultSort, VList::length(vars));
  }
} // constructOperatorType

/**
 *
 * @author Giles
 */
TPTP::SourceRecord* TPTP::getSource()
{
  if (getTok(0).tag != T_COMMA) { // if comma is not there, source was not provided
    return 0;
  }

  consumeToken(T_COMMA);

  //Either source is a file or an inference, otherwise we don't care about it!
  //  therefore failing will return 0
 
  Token& source_kind = getTok(0);
  if(source_kind.tag != T_NAME) return 0;

  resetToks();
  if (getTok(0).tag != T_LPAR) {
    return 0;
  } else {
    resetToks();
  }
  
  //file
  if(source_kind.content == "file"){
    std::string fileName = getTok(0).content;
    resetToks();
    consumeToken(T_COMMA);
    resetToks();
    std::string nameInFile = getTok(0).content;
    resetToks();

    // cout << "Creating file source record for " << fileName << " and " << nameInFile << endl;

    consumeToken(T_RPAR);
    return new FileSourceRecord(fileName,nameInFile);
  }
  // inference
  else if(source_kind.content == "inference" || source_kind.content == "introduced"){
    bool introduced = (source_kind.content == "introduced");
    std::string name = getTok(0).content;
    resetToks();

    // cout << "Creating inference source record for " << name <<  endl;

    InferenceSourceRecord* r = new InferenceSourceRecord(name);

    if(introduced){
      // then we don't expect names and we don't care about middle info 
      resetToks();
      skipToRPAR();
      return r;
    }

    // now skip this middle information that is between [ and ]
    consumeToken(T_COMMA);
    consumeToken(T_LBRA);
    skipToRBRA();
    consumeToken(T_COMMA);
    consumeToken(T_LBRA);

    // read comma separated list of names
    Token tok;
    while((tok=getTok(0)).tag != T_RBRA){
      resetToks();
      if(tok.tag == T_COMMA) continue;
   
      if (tok.tag != T_NAME && tok.tag != T_INT) {
        cout << "read token " << tok.tag << " with content " << tok.content << endl;

        // TODO: parse errors are nice, but maybe we just want to ignore any info which we cannot understand?

        PARSE_ERROR("Source unit name expected",tok);
      }

      std::string premise = tok.content;

      tok = getTok(0);
      if (tok.tag != T_COMMA && tok.tag != T_RBRA) {
        // if the next thing is neither comma not RBRA, it is an ugly info piece we want to skip
        resetToks();
        skipToRPAR();
      } else {
        r->premises.push(premise);
        // cout << "pushed premise " << premise << endl;
      }
    }
    resetToks();

    consumeToken(T_RPAR);
    return r;
  } else {
    
    skipToRPAR();
  }

  return 0;
}


/**
 * Skip any sequence tokens, including matching pairs of left parentheses,
 * until an unmatched right parenthesis is found. Consume this right parenthesis
 * and terminate.
 * @since 15/07/2011 Manchester
 */
void TPTP::skipToRPAR()
{
  int balance = 0;
  for (;;) {
    Token tok = getTok(0);
    switch (tok.tag) {
    case T_EOF:
      PARSE_ERROR(") not found",tok);
    case T_LPAR:
      resetToks();
      balance++;
      break;
    case T_RPAR:
      resetToks();
      balance--;
      if (balance == -1) {
	return;
      }
      break;
    default:
      resetToks();
      break;
    }
  }
} // skipToRPAR

/**
 * A copy of skipToRPAR but for BRA
 */
void TPTP::skipToRBRA()
{
  int balance = 0;
  for (;;) {
    Token tok = getTok(0);
    switch (tok.tag) {
    case T_EOF:
      PARSE_ERROR(") not found",tok);
    case T_LBRA:
      resetToks();
      balance++;
      break;
    case T_RBRA:
      resetToks();
      balance--;
      if (balance == -1) {
        return;
      }
      break;
    default:
      resetToks();
      break;
    }
  }
} // skipToRBRA

/**
 * Read a simple formula (quantified formula, negation,
 * formula in parentheses, true or false).
 * @since 10/04/2011 Manchester
 */
void TPTP::simpleFormula()
{
  Token tok = getTok(0);

  switch (tok.tag) {
  case T_NOT:
    resetToks();
    _connectives.push(NOT);
    _states.push(SIMPLE_FORMULA);
    return;

  case T_FORALL:
  case T_EXISTS:
    resetToks();
    consumeToken(T_LBRA);
    _connectives.push(tok.tag == T_FORALL ? FORALL : EXISTS);
    _states.push(UNBIND_VARIABLES);
    _states.push(SIMPLE_FORMULA);
    addTagState(T_COLON);
    addTagState(T_RBRA);
    _states.push(VAR_LIST);
    return;

  case T_LPAR:
    resetToks();
    addTagState(T_RPAR);
    _states.push(FORMULA);
    return;

  case T_STRING:
    _states.push(END_EQ);
    _states.push(TERM);
    _states.push(MID_EQ);
    _states.push(TERM);
    return;
  case T_INT:
  case T_RAT:
  case T_REAL:
    //PARSE_ERROR("Sorry, polymorphic Vampire does not yet support theories", tok);
    _states.push(END_EQ);
    _states.push(TERM);
    _states.push(MID_EQ);
    _states.push(TERM);
    return;
  case T_TRUE:
    resetToks();
    _formulas.push(new Formula(true));
    return;
  case T_FALSE:
    resetToks();
    _formulas.push(new Formula(false));
    return;
  case T_NAME:
  case T_VAR:
  case T_ITE:
  case T_THEORY_FUNCTION:
  case T_LET:
  case T_LBRA:
    _states.push(FORMULA_INFIX);
    _states.push(FUN_APP);
    return;
  default:
    PARSE_ERROR("formula or term expected",tok);
  }
} // simpleFormula

/**
 * Unbind variable sort definition.
 * @since 14/07/2011 Manchester
 */
void TPTP::unbindVariables()
{
  VList::Iterator vs(_bindLists.pop());
  while (vs.hasNext()) {
    unsigned var = vs.next();
    SList** sorts = _variableSorts.getPtr(var); // sorts is now a pointer to the list stored inside _variableSorts
    ALWAYS(sorts);
    SList::pop(*sorts); // so this will modify that stored list
  }
} // unbindVariables

/**
 * Read a simple type: name or type in parentheses
 * @since 14/07/2011 Manchester
 */
void TPTP::simpleType()
{
  Token& tok = getTok(0);

  if(tok.tag == T_TYPE_QUANT) {
    _containsPolymorphism = true;
    resetToks();
    _typeTags.push(TT_QUANTIFIED);
    consumeToken(T_LBRA);
    _states.push(UNBIND_VARIABLES);
    _states.push(TYPE);
    addTagState(T_COLON);
    addTagState(T_RBRA);
    _states.push(VAR_LIST);
    return;
  }

  if(_isThf){
    _types.push(new AtomicType(readArrowSort()));
    return;
  } 

  if (tok.tag == T_LPAR) {
    resetToks();
    addTagState(T_RPAR);
    _states.push(TYPE);
    return;
  }
  _types.push(new AtomicType(readSort()));
  
} // simpleType


/**
 * Read a HOL sort and return its number 
 * @since 10/11/2017 Leicester
 * @author Ahmed Bhayat
 */
 
TermList TPTP::readArrowSort()
{
  int inBrackets = 0;
  TermStack terms;
  Token tok = getTok(0);
  TermList sort;
  TermList dummy = TermList(0, true);
  while((tok.tag != T_COMMA) && (tok.tag != T_RBRA) && (tok.tag != T_APP)){
    switch(tok.tag){
      case T_LPAR: //This will need changing when we read tuple types - AYB
        terms.push(dummy);
        inBrackets += 1;
        break;
      case T_ARROW:
        break;
      case T_RPAR:
        inBrackets -= 1;
        if(inBrackets < 0){_gpos = 0; goto afterWhile;}
          foldl(&terms);
          break;
      default:{
        sort = readSort();
        terms.push(sort);
        if(!sort.isVar() && sort.term()->arity()){
          tok = getTok(0);
          continue; 
        }               
      }
    }
    resetToks();
    tok = getTok(0);
  }
afterWhile:
  if(terms.size() != 1){
    foldl(&terms);
  }
  ASS(terms.size() == 1);
  return terms.pop();   
}
 
void TPTP::foldl(TermStack* terms)
{
  TermList item1 = terms->pop();
  TermList item2 = terms->pop();
  while(!(terms->isEmpty()) && (!item2.isSpecialVar())){
    item1 = AtomicSort::arrowSort(item2, item1);
    item2 = terms->pop();
  }
  if (!item2.isSpecialVar()){
    item1 = AtomicSort::arrowSort(item2, item1);;
  }
  terms->push(item1);
}   

void TPTP::readTypeArgs(unsigned arity)
{
  for(unsigned i = 0; i < arity; i++){
    consumeToken(T_APP);
    Token tok = getTok(0);
    if(tok.tag == T_LPAR){
      resetToks();
      _termLists.push(readArrowSort());
      consumeToken(T_RPAR);
    } else {
      _termLists.push(readArrowSort());            
    }
  }
}

/**
 * Read a sort and return its number. If a sort is not built-in, then raise an
 * exception if it has been declared and newSortExpected, or it has not been
 * declared and newSortExpected is false.
 * @since 14/07/2011 Manchester
 */
TermList TPTP::readSort()
{
  Token tok = getTok(0); 
  resetToks();
  switch (tok.tag) {
  case T_NAME:
    {
      unsigned arity = 0;
      std::string fname = tok.content;
      if(_isThf){
        arity = _typeConstructorArities.find(fname) ? _typeConstructorArities.get(fname) : 0;
        readTypeArgs(arity);
      } else {
        int c = getChar(0);
        //Polymorphic sorts of are of the form 
        //type_con(sort_1, ..., sort_n)
        //the same as standard first-order terms.
        //Code below works, but does not fit the philosophy of
        //this parser. However, recursive calls to readSort are
        //used in for array sorts and tuple sorts, so polymorphism
        //isn't uniquely evil on this front!
        if(c == '('){
          consumeToken(T_LPAR);    
          for(;;){
            arity++;
            _termLists.push(readSort());
            tok = getTok(0);
            if(tok.tag == T_COMMA){
              consumeToken(T_COMMA);
            }else if(tok.tag == T_RPAR){
              consumeToken(T_RPAR);
              break;
            } else{
              ASSERTION_VIOLATION;
            }
          }
        }
      } 
      return createTypeConApplication(fname, arity);
    }
  case T_VAR:
    {
      std::string vname = tok.content;
      unsigned var = (unsigned)_vars.insert(vname);
      return  TermList(var, false);
    }

  case T_DEFAULT_TYPE:
    return AtomicSort::defaultSort();

  case T_BOOL_TYPE:
    return AtomicSort::boolSort();

  case T_INTEGER_TYPE:
    return AtomicSort::intSort();

  case T_RATIONAL_TYPE:
    return AtomicSort::rationalSort();

  case T_REAL_TYPE:
    return AtomicSort::realSort();

  case T_TTYPE:
    return AtomicSort::superSort();

  case T_LBRA:
  {
    TermStack sorts;
    for (;;) {
      TermList sort = readSort();
      sorts.push(sort);
      if (getTok(0).tag == T_COMMA) {
        resetToks();
      } else {
        consumeToken(T_RBRA);
        break;
      }
    }

    if (sorts.length() < 2) {
      USER_ERROR("Tuple sort with less than two arguments");
    }

    return AtomicSort::tupleSort((unsigned) sorts.length(), sorts.begin());
  }
  case T_THEORY_SORT: {
    TermList sort;
    consumeToken(T_LPAR);
    switch (getTheorySort(tok)) {
      case TS_ARRAY: {
        TermList indexSort = readSort();
        consumeToken(T_COMMA);
        TermList innerSort = readSort();
        sort = AtomicSort::arraySort(indexSort, innerSort);
        break;
      }
      default:
        ASSERTION_VIOLATION;
    }
    consumeToken(T_RPAR);
    return sort;
  }
  default:
    PARSE_ERROR("sort expected",tok);
  }
} // readSort

/**
 * True if c1 has a strictly higher priority than c2.
 * @since 07/07/2011 Manchester
 */
bool TPTP::higherPrecedence(int c1,int c2)
{
  if (c1 == APP) return true;
  if (c1 == c2) return false;
  if (c1 == -1) return false;
  if (c2 == IFF) return true;
  if (c1 == IFF) return false;
  if (c2 == XOR) return true;
  if (c1 == XOR) return false;
  if (c2 == IMP) return true;
  if (c1 == IMP) return false;
  if (c2 == OR) return true;
  if (c1 == OR) return false;
  ASSERTION_VIOLATION;
} // higherPriority

bool TPTP::findInterpretedPredicate(std::string name, unsigned arity) {
  if (name == "$evaleq" || name == "$equal" || name == "$distinct") {
    return true;
  }

  if (name == "$is_int" || name == "$is_rat") {
    return arity == 1;
  }

  if (name == "$less" || name == "$lesseq" || name == "$greater" || name == "$greatereq" || name == "$divides") {
    return arity == 2;
  }

  return false;
}

/**
 * Create an and- or or-formula flattening its lhs and rhs if necessary.
 * @since 07/07/2011 Manchester
 */
Formula* TPTP::makeJunction (Connective c,Formula* lhs,Formula* rhs)
{
  if (lhs->connective() == c) {
    FormulaList* largs = lhs->args();

    if (rhs->connective() == c) {
      FormulaList::concat(largs,rhs->args());
      delete static_cast<JunctionFormula*>(rhs);
      return lhs;
    }
    // only lhs has c as the main connective
    FormulaList::concat(largs,new FormulaList(rhs));
    return lhs;
  }
  // lhs' connective is not c
  if (rhs->connective() == c) {
    static_cast<JunctionFormula*>(rhs)->setArgs(new FormulaList(lhs,
								rhs->args()));
    return rhs;
  }
  // both connectives are not c
  return new JunctionFormula(c,
			     new FormulaList(lhs,
					     new FormulaList(rhs)));
} // makeJunction

/** Add a function to the signature
 * @param name the function name
 * @param arity the function arity
 * @param added if the function is new, will be assigned true, otherwise false
 * @param arg some argument of the function, require to resolve its type for overloaded
 *        built-in functions
 */
unsigned TPTP::addFunction(std::string name,int arity,bool& added,TermList& arg)
{
  if (name == "$sum") {
    return addOverloadedFunction(name,arity,2,added,arg,
				 Theory::INT_PLUS,
				 Theory::RAT_PLUS,
				 Theory::REAL_PLUS);
  }
  if (name == "$difference") {
    return addOverloadedFunction(name,arity,2,added,arg,
				 Theory::INT_MINUS,
				 Theory::RAT_MINUS,
				 Theory::REAL_MINUS);
  }
  if (name == "$product") {
    return addOverloadedFunction(name,arity,2,added,arg,
				 Theory::INT_MULTIPLY,
				 Theory::RAT_MULTIPLY,
				 Theory::REAL_MULTIPLY);
  }
  // An odd leftover, maps to the 'most natural' kind of division
  if (name == "$divide") {
    return addOverloadedFunction(name,arity,2,added,arg,
				 Theory::INT_QUOTIENT_E,
				 Theory::RAT_QUOTIENT,
				 Theory::REAL_QUOTIENT);
  }
  if (name == "$modulo"){
    if(sortOf(arg)!=AtomicSort::intSort()){
      USER_ERROR("$modulo can only be used with integer type");
    }
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_REMAINDER_E,  // $modulo is the always positive remainder, therefore INT_REMAINDER_E
                                 Theory::INT_REMAINDER_E,  // will not be used
                                 Theory::INT_REMAINDER_E); // will not be used
  }
  if (name == "$abs"){
    if(sortOf(arg)!=AtomicSort::intSort()){
      USER_ERROR("$abs can only be used with integer type");
    }
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_ABS,
                                 Theory::INT_ABS,  // will not be used
                                 Theory::INT_ABS); // will not be used
  }
  if (name == "$quotient") {
    if(sortOf(arg)==AtomicSort::intSort()){
      USER_ERROR("$quotient cannot be used with integer type");
    }
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_QUOTIENT_E,// this is a dummy
                                 Theory::RAT_QUOTIENT,
                                 Theory::REAL_QUOTIENT);
  }
  if (name == "$quotient_e") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_QUOTIENT_E,
                                 Theory::RAT_QUOTIENT_E,
                                 Theory::REAL_QUOTIENT_E);
  }
  if (name == "$quotient_t") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_QUOTIENT_T,
                                 Theory::RAT_QUOTIENT_T,
                                 Theory::REAL_QUOTIENT_T);
  }
  if (name == "$quotient_f") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_QUOTIENT_F,
                                 Theory::RAT_QUOTIENT_F,
                                 Theory::REAL_QUOTIENT_F);
  }
  if (name == "$remainder_e") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_REMAINDER_E,
                                 Theory::RAT_REMAINDER_E,
                                 Theory::REAL_REMAINDER_E);
  }
  if (name == "$remainder_t") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_REMAINDER_T,
                                 Theory::RAT_REMAINDER_T,
                                 Theory::REAL_REMAINDER_T);
  }
  if (name == "$remainder_f") {
    return addOverloadedFunction(name,arity,2,added,arg,
                                 Theory::INT_REMAINDER_F,
                                 Theory::RAT_REMAINDER_F,
                                 Theory::REAL_REMAINDER_F);
  }
  if (name == "$uminus") {
    return addOverloadedFunction(name,arity,1,added,arg,
				 Theory::INT_UNARY_MINUS,
				 Theory::RAT_UNARY_MINUS,
				 Theory::REAL_UNARY_MINUS);
  }
  if (name == "$successor"){
    if(sortOf(arg)!=AtomicSort::intSort()){
      USER_ERROR("$succ can only be used with integer type");
    }
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_SUCCESSOR,
                                 Theory::INT_SUCCESSOR,  // will not be used
                                 Theory::INT_SUCCESSOR); // will not be used
  }
  if (name == "$floor") {
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_FLOOR,
                                 Theory::RAT_FLOOR,
                                 Theory::REAL_FLOOR);
  }
  if (name == "$ceiling") {
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_CEILING,
                                 Theory::RAT_CEILING,
                                 Theory::REAL_CEILING);
  }
  if (name == "$truncate") {
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_TRUNCATE,
                                 Theory::RAT_TRUNCATE,
                                 Theory::REAL_TRUNCATE);
  }
  if (name == "$round") {
    return addOverloadedFunction(name,arity,1,added,arg,
                                 Theory::INT_ROUND,
                                 Theory::RAT_ROUND,
                                 Theory::REAL_ROUND);
  }
  if (name == "$to_int") {
    return addOverloadedFunction(name,arity,1,added,arg,
				 Theory::INT_TO_INT,
				 Theory::RAT_TO_INT,
				 Theory::REAL_TO_INT);
  }
  if (name == "$to_rat") {
    return addOverloadedFunction(name,arity,1,added,arg,
				 Theory::INT_TO_RAT,
				 Theory::RAT_TO_RAT,
				 Theory::REAL_TO_RAT);
  }
  if (name == "$to_real") {
    return addOverloadedFunction(name,arity,1,added,arg,
				 Theory::INT_TO_REAL,
				 Theory::RAT_TO_REAL,
				 Theory::REAL_TO_REAL);
  } 
  if (name == "vPI"  || name == "vSIGMA"){
    return env.signature->getPiSigmaProxy(name); 
  }
  if (arity > 0) {
    return env.signature->addFunction(name,arity,added);
  }
  return addUninterpretedConstant(name,_overflow,added);
} // addFunction

/** Add a predicate to the signature
 * @param name the predicate name
 * @param arity the predicate arity
 * @param added if the predicate is new, will be assigned true, otherwise false
 * @param arg some argument of the predicate, require to resolve its type for overloaded
 *        built-in predicates
 * @return the predicate number in the signature, or -1 if it is a different name for an equality
 *         predicate
 */
int TPTP::addPredicate(std::string name,int arity,bool& added,TermList& arg)
{
  if (name == "$evaleq" || name == "$equal") {
    return -1;
  }
  if (name == "$less") {
    return addOverloadedPredicate(name,arity,2,added,arg,
				  Theory::INT_LESS,
				  Theory::RAT_LESS,
				  Theory::REAL_LESS);
  }
  if (name == "$lesseq") {
    return addOverloadedPredicate(name,arity,2,added,arg,
				  Theory::INT_LESS_EQUAL,
				  Theory::RAT_LESS_EQUAL,
				  Theory::REAL_LESS_EQUAL);
  }
  if (name == "$greater") {
    return addOverloadedPredicate(name,arity,2,added,arg,
				  Theory::INT_GREATER,
				  Theory::RAT_GREATER,
				  Theory::REAL_GREATER);
  }
  if (name == "$greatereq") {
    return addOverloadedPredicate(name,arity,2,added,arg,
				  Theory::INT_GREATER_EQUAL,
				  Theory::RAT_GREATER_EQUAL,
				  Theory::REAL_GREATER_EQUAL);
  }
  if (name == "$is_int") {
    return addOverloadedPredicate(name,arity,1,added,arg,
				  Theory::INT_IS_INT,
				  Theory::RAT_IS_INT,
				  Theory::REAL_IS_INT);
  }
  if (name == "$divides"){
    if(sortOf(arg)!=AtomicSort::intSort()){
      USER_ERROR("$divides can only be used with integer type");
    }
    return addOverloadedPredicate(name,arity,2,added,arg,
                                  Theory::INT_DIVIDES,
                                  Theory::INT_DIVIDES,  // will not be used
                                  Theory::INT_DIVIDES); // will not be used
  }
  if (name == "$is_rat") {
    return addOverloadedPredicate(name,arity,1,added,arg,
				  Theory::INT_IS_RAT,
				  Theory::RAT_IS_RAT,
				  Theory::REAL_IS_RAT);
  } 
  if(name == "$distinct"){
    // special case for distinct, dealt with in formulaInfix
    return -2;
  }
  return env.signature->addPredicate(name,arity,added);
} // addPredicate


unsigned TPTP::addOverloadedFunction(std::string name,int arity,int symbolArity,bool& added,TermList& arg,
				     Theory::Interpretation integer,Theory::Interpretation rational,
				     Theory::Interpretation real)
{
  if (arity != symbolArity) {
    USER_ERROR(name + " is used with " + Int::toString(arity) + " argument(s) when there were "+Int::toString(symbolArity)+" expected");
  }
  TermList srt = sortOf(arg);
  TermList* n = arg.next();
  for(int i=1;i<arity;i++){
    if(sortOf(*n)!=srt){
      std::string msg = "The interpreted function symbol " + name + " is not used with a single sort.";
      msg += "\nArgument 0 is "+srt.toString()+" and argument "+Lib::Int::toString(i)+" is "+sortOf(*n).toString();
      if(_isFof){ msg += "\nCheck that you are using tff if you want numbers to be interpreted"; }
      USER_ERROR(msg);
    }
    n = n->next();
  }
  if (srt == AtomicSort::intSort()) {
    return env.signature->addInterpretedFunction(integer,name);
  }
  if (srt == AtomicSort::rationalSort()) {
    return env.signature->addInterpretedFunction(rational,name);
  }
  if (srt == AtomicSort::realSort()) {
    return env.signature->addInterpretedFunction(real,name);
  }
  USER_ERROR((std::string)"The symbol " + name + " is used with a non-numeric type");
} // addOverloadedFunction

unsigned TPTP::addOverloadedPredicate(std::string name,int arity,int symbolArity,bool& added,TermList& arg,
				     Theory::Interpretation integer,Theory::Interpretation rational,
				     Theory::Interpretation real)
{
  if (arity != symbolArity) {
    USER_ERROR(name + " is used with " + Int::toString(arity) + " argument(s) when there were "+Int::toString(symbolArity)+" expected");
  }
  TermList srt = sortOf(arg);
  TermList* n = arg.next();
  for(int i=1;i<arity;i++){
    if(sortOf(*n)!=srt){
      std::string msg = "The interpreted predicate symbol " + name + " is not used with a single sort.";
      msg += "\nArgument 0 is "+srt.toString()+" and argument "+Lib::Int::toString(i)+" is "+sortOf(*n).toString();
      if(_isFof){ msg += "Check that you are using tff if you want numbers to be interpreted"; }
      USER_ERROR(msg);
    }
    n = n->next(); 
  }
  
  if (srt == AtomicSort::intSort()) {
    return env.signature->addInterpretedPredicate(integer,name);
  }
  if (srt == AtomicSort::rationalSort()) {
    return env.signature->addInterpretedPredicate(rational,name);
  }
  if (srt == AtomicSort::realSort()) {
    return env.signature->addInterpretedPredicate(real,name);
  }
  USER_ERROR((std::string)"The symbol " + name + " is used with a non-numeric type");
} // addOverloadedPredicate

/**
 * Return the sort of the term.
 * @since 29/07/2011 Manchester
 * @since 03/05/2013 train Manchester-London bug fix
 * @author Andrei Voronkov
 */
TermList TPTP::sortOf(TermList t)
{
  for (;;) {
    if (t.isVar()) {
      SList* sorts;
      if (_variableSorts.find(t.var(),sorts) && SList::isNonEmpty(sorts)) {
        return sorts->head();
      }
      // there might be variables whose sort is undeclared,
      // in this case they have the default sort
      TermList def = AtomicSort::defaultSort();
      bindVariable(t.var(), def);
      return def;
    }
    TermList sort;
    TermList mvar;
    if (SortHelper::getResultSortOrMasterVariable(t.term(), sort, mvar)) { //TODO update
      return sort;
    } else {
      t = mvar;
    }
  }
} // sortOf

/**
 * Add an integer constant by reading it from the std::string name.
 * If it overflows, create an uninterpreted constant of the
 * integer type and the name 'name'. Check that the name of the constant
 * does not collide with user-introduced names of uninterpreted constants.
 * @since 22/07/2011 Manchester
 * @since 03/05/2013 train Manchester-London, bug fix: integers are treated
 *   as terms of the default sort when fof() or cnf() is used
 * @author Andrei Voronkov
 */
unsigned TPTP::addIntegerConstant(const std::string& name, Set<std::string>& overflow, bool defaultSort)
{
  try {
    return env.signature->addIntegerConstant(name,defaultSort);
  }
  catch (Kernel::ArithmeticException&) {
    bool added;
    unsigned fun = env.signature->addFunction(name,0,added,true /* overflown constant*/);
    if (added) {
      overflow.insert(name);
      Signature::Symbol* symbol = env.signature->getFunction(fun);
      symbol->setType(OperatorType::getConstantsType(defaultSort ? AtomicSort::defaultSort() : AtomicSort::intSort()));
    }
    else if (!overflow.contains(name)) {
      USER_ERROR((std::string)"Cannot use name '" + name + "' as an atom name since it collides with an integer number");
    }
    return fun;
  }
} // TPTP::addIntegerConstant

/**
 * Add an rational constant by reading it from the std::string name.
 * If it overflows, create an uninterpreted constant of the
 * rational type and the name 'name'. Check that the name of the constant
 * does not collide with user-introduced names of uninterpreted constants.
 * @since 22/07/2011 Manchester
 * @since 03/05/2013 train Manchester-London, fix to handle difference
 *    between treating rationals using fof() and tff()
 * @author Andrei Voronkov
 */
unsigned TPTP::addRationalConstant(const std::string& name, Set<std::string>& overflow, bool defaultSort)
{
  size_t i = name.find_first_of("/");
  ASS(i != std::string::npos);
  try {
    return env.signature->addRationalConstant(name.substr(0,i),
					      name.substr(i+1),
					      defaultSort);
  }
  catch(Kernel::ArithmeticException&) {
    bool added;
    unsigned fun = env.signature->addFunction(name,0,added,true /* overflown constant*/);
    if (added) {
      overflow.insert(name);
      Signature::Symbol* symbol = env.signature->getFunction(fun);
      symbol->setType(OperatorType::getConstantsType(defaultSort ? AtomicSort::defaultSort() : AtomicSort::rationalSort()));
    }
    else if (!overflow.contains(name)) {
      USER_ERROR((std::string)"Cannot use name '" + name + "' as an atom name since it collides with an rational number");
    }
    return fun;
  }
} // TPTP::addRationalConstant

/**
 * Add an real constant by reading it from the std::string name.
 * If it overflows, create an uninterpreted constant of the
 * real type and the name 'name'. Check that the name of the constant
 * does not collide with user-introduced names of uninterpreted constants.
 * @since 22/07/2011 Manchester
 * @since 03/05/2013 train Manchester-London, fix to handle difference
 *    between treating rationals using fof() and tff()
 * @author Andrei Voronkov
 */
unsigned TPTP::addRealConstant(const std::string& name, Set<std::string>& overflow, bool defaultSort)
{
  try {
    return env.signature->addRealConstant(name,defaultSort);
  }
  catch(Kernel::ArithmeticException&) {
    bool added;
    unsigned fun = env.signature->addFunction(name,0,added,true /* overflown constant*/);
    if (added) {
      overflow.insert(name);
      Signature::Symbol* symbol = env.signature->getFunction(fun);
      symbol->setType(OperatorType::getConstantsType(defaultSort ? AtomicSort::defaultSort() : AtomicSort::realSort()));
    }
    else if (!overflow.contains(name)) {
      USER_ERROR((std::string)"Cannot use name '" + name + "' as an atom name since it collides with an real number");
    }
    return fun;
  }
}  // TPTP::addRealConstant


/**
 * Add an uninterpreted constant by reading it from the std::string name.
 * Check that the name of the constant does not collide with uninterpreted constants
 * created by the parser from overflown input numbers.
 * @since 22/07/2011 Manchester
 */
unsigned TPTP::addUninterpretedConstant(const std::string& name, Set<std::string>& overflow, bool& added)
{
  if (overflow.contains(name)) {
    USER_ERROR((std::string)"Cannot use name '" + name + "' as an atom name since it collides with an integer number");
  }
  //TODO make sure Vampire internal names are unique to Vampire
  //and cannot occur in the input AYB
  if(name == "vAND" || name == "vOR" || name == "vIMP" ||
     name == "vIFF" || name == "vXOR"){
    return env.signature->getBinaryProxy(name);
  } else if (name == "vNOT"){
    return env.signature->getNotProxy();
  }
  return env.signature->addFunction(name,0,added);
} // TPTP::addUninterpretedConstant

/**
 * Associate name @b name with unit @b unit
 * Each formula can have its name assigned at most once
 */
void TPTP::assignAxiomName(const Unit* unit, std::string& name)
{
  ALWAYS(_axiomNames.insert(unit->number(), name));
} // TPTP::assignAxiomName

/**
 * If @b unit has a name associated, assign it into @b result,
 * and return true; otherwise return false
 */
bool TPTP::findAxiomName(const Unit* unit, std::string& result)
{
  return _axiomNames.find(unit->number(), result);
} // TPTP::findAxiomName

/**
 * Process vampire() declaration
 * @since 25/08/2009 Redmond
 */
void TPTP::vampire()
{
  consumeToken(T_LPAR);
  std::string nm = name();

  if (nm == "option") { // vampire(option,age_weight_ratio,3)
    consumeToken(T_COMMA);
    std::string opt = name();
    consumeToken(T_COMMA);
    Token tok = getTok(0);
    switch (tok.tag) {
    case T_INT:
    case T_REAL:
    case T_NAME:
      env.options->set(opt,tok.content);
      resetToks();
      break;
    default:
      PARSE_ERROR("either atom or number expected as a value of a Vampire option",tok);
    }
  }
  // Allows us to insert LaTeX templates for predicate and function symbols
  else if(nm == "latex"){
    consumeToken(T_COMMA);
    std::string kind = name();
    bool pred;
    if (kind == "predicate") {
      pred = true;
    }
    else if (kind == "function") {
      pred = false;
    }
    else {
      PARSE_ERROR("either 'predicate' or 'function' expected",getTok(0));
    }
    consumeToken(T_COMMA);
    std::string symb = name();
    consumeToken(T_COMMA);
    Token tok = getTok(0);
    if (tok.tag != T_INT) {
      PARSE_ERROR("a non-negative integer (denoting arity) expected",tok);
    }
    unsigned arity;
    if (!Int::stringToUnsignedInt(tok.content,arity)) {
      PARSE_ERROR("a number denoting arity expected",tok);
    }
    resetToks();
    consumeToken(T_COMMA);
    tok = getTok(0);
    if(tok.tag != T_STRING){
      PARSE_ERROR("a template string expected",tok);
    }
    std::string temp = tok.content;
    resetToks();
    if(pred){
      consumeToken(T_COMMA);
      std::string pol= name();
      bool polarity;
      if(pol=="true"){polarity=true;}else if(pol=="false"){polarity=false;}
      else{ PARSE_ERROR("polarity expected (true/false)",getTok(0)); }
      unsigned f = env.signature->addPredicate(symb,arity);
      theory->registerLaTeXPredName(f,polarity,temp);
    }
    else{
      unsigned f = env.signature->addFunction(symb,arity);
      theory->registerLaTeXFuncName(f,temp);
    }
  }
  else if (nm == "symbol") {
    consumeToken(T_COMMA);
    std::string kind = name();
    bool pred;
    if (kind == "predicate") {
      pred = true;
    }
    else if (kind == "function") {
      pred = false;
    }
    else {
      PARSE_ERROR("either 'predicate' or 'function' expected",getTok(0));
    }
    consumeToken(T_COMMA);
    std::string symb = name();
    consumeToken(T_COMMA);
    Token tok = getTok(0);
    if (tok.tag != T_INT) {
      PARSE_ERROR("a non-negative integer (denoting arity) expected",tok);
    }
    unsigned arity;
    if (!Int::stringToUnsignedInt(tok.content,arity)) {
      PARSE_ERROR("a number denoting arity expected",tok);
    }
    resetToks();
    consumeToken(T_COMMA);
    Color color = COLOR_INVALID;
    bool skip = false, uncomputable = false;
    std::string lr = name();
    if (lr == "left") {
      color=COLOR_LEFT;
    }
    else if (lr == "right") {
      color=COLOR_RIGHT;
    }
    else if (lr == "skip") {
      skip = true;
    }
    else if (lr == "uncomputable") {
      uncomputable = true;
    }
    else {
      PARSE_ERROR("'left', 'right', 'skip' or 'uncomputable' expected",getTok(0));
    }
    if (!uncomputable) {
      env.colorUsed = true;
    }
    Signature::Symbol* sym = pred
                             ? env.signature->getPredicate(env.signature->addPredicate(symb,arity))
                             : env.signature->getFunction(env.signature->addFunction(symb,arity));
    if (skip) {
      sym->markSkip();
    }
    else if (uncomputable) {
      sym->markUncomputable();
    }
    else {
      ASS_NEQ(color, COLOR_INVALID);
      sym->addColor(color);
    }
  }
  else if (nm == "left_formula") { // e.g. vampire(left_formula)
    _currentColor = COLOR_LEFT;
  }
  else if (nm == "right_formula") { // e.g. vampire(left_formula)
    _currentColor = COLOR_RIGHT;
  }
  else if (nm == "end_formula") { // e.g. vampire(left_formula)
    _currentColor = COLOR_TRANSPARENT;
  }
  else if (nm == "model_check"){
    consumeToken(T_COMMA);
    std::string command = name();
    if(command == "formulas_start"){
      _modelDefinition = false;
    }
    else if(command == "formulas_end"){
      // do nothing
    }
    else if(command == "model_start"){
      _modelDefinition = true;
    }
    else if(command == "model_end"){
      _modelDefinition = false;
    }
    else USER_ERROR("Unknown model_check command");
  }
  else {
    USER_ERROR((std::string)"Unknown vampire directive: "+nm);
  }
  consumeToken(T_RPAR);
  consumeToken(T_DOT);
} // vampire

#if VDEBUG
const char* TPTP::toString(State s)
{
  switch (s) {
  case UNIT_LIST:
    return "UNIT_LIST";
  case CNF:
    return "CNF";
  case FOF:
    return "FOF";
  case VAMPIRE:
    return "VAMPIRE";
  case FORMULA:
    return "FORMULA";
  case END_FOF:
    return "END_FOF";
  case SIMPLE_FORMULA:
    return "SIMPLE_FORMULA";
  case END_FORMULA:
    return "END_FORMULA";
  case FORMULA_INSIDE_TERM:
    return "FORMULA_INSIDE_TERM";
  case END_FORMULA_INSIDE_TERM:
    return "END_FORMULA_INSIDE_TERM";
  case END_TERM_AS_FORMULA:
    return "END_TERM_AS_FORMULA";
  case VAR_LIST:
    return "VAR_LIST";
  case FUN_APP:
    return "FUN_APP";
  case FORMULA_INFIX:
    return "FORMULA_INFIX";
  case ARGS:
    return "ARGS";
  case TERM:
    return "TERM";
  case TERM_INFIX:
    return "TERM_INFIX";
  case END_TERM:
    return "END_TERM";
  case TAG:
    return "TAG";
  case INCLUDE:
    return "INCLUDE";
  case END_EQ:
    return "END_EQ";
  case TFF:
    return "TFF";
  case THF:
    return "THF";
  case TYPE:
    return "TYPE";
  case END_TFF:
    return "END_TFF";
  case END_APP:
    return "END_APP";
  case HOL_FORMULA:
    return "HOL_FORMULA";
  case END_HOL_FORMULA:
    return "END_HOL_FORMULA";
  case HOL_TERM:
    return "HOL_TERM";
  case END_TYPE:
    return "END_TYPE";
  case SIMPLE_TYPE:
    return "SIMPLE_TYPE";
  case END_THEORY_FUNCTION:
    return "END_THEORY_FUNCTION";
  case END_ARGS:
    return "END_ARGS";
  case MID_EQ:
    return "MID_EQ";
  case LET_TYPE:
    return "LET_TYPE";
  case END_LET_TYPES:
    return "END_LET_TYPES";
  case DEFINITION:
    return "DEFINITION";
  case MID_DEFINITION:
    return "MID_DEFINITION";
  case END_DEFINITION:
    return "END_DEFINITION";
  case SYMBOL_DEFINITION:
    return "SYMBOL_DEFINITION";
  case TUPLE_DEFINITION:
    return "TUPLE_DEFINITION";
  case END_LET:
    return "END_LET";
  case UNBIND_VARIABLES:
    return "UNBIND_VARIABLES";
  case END_ITE:
    return "END_ITE";
  case END_TUPLE:
    return "END_TUPLE";
  default:
    cout << (int)s << "\n";
    ASS(false);
    break;
  }
}
#endif

#ifdef DEBUG_SHOW_STATE
void TPTP::printStacks() {

  Stack<State>::Iterator stit(_states);
  cout << "States:";
  if   (!stit.hasNext()) cout << " <empty>";
  while (stit.hasNext()) cout << " " << toString(stit.next());
  cout << endl;

  /*Stack<Type*>::Iterator tyit(_types);
  cout << "Types:";
  if   (!tyit.hasNext()) cout << " <empty>";
  while (tyit.hasNext()) cout << " " << tyit.next()->tag();
  cout << endl;

  Stack<int>::Iterator cit(_connectives);
  cout << "Connectives:";
  if   (!cit.hasNext()) cout << " <empty>";
  while (cit.hasNext()) cout << " " << cit.next();
  cout << endl; 

  Stack<std::string>::Iterator sit(_strings);
  cout << "Strings:";
  if   (!sit.hasNext()) cout << " <empty>";
  while (sit.hasNext()) cout << " " << sit.next();
  cout << endl;

  Stack<int>::Iterator iit(_ints);
  cout << "Ints:";
  if   (!iit.hasNext()) cout << " <empty>";
  while (iit.hasNext()) cout << " " << iit.next();
  cout << endl;

  Stack<bool>::Iterator bit(_bools);
  cout << "Bools:";
  if   (!bit.hasNext()) cout << " <empty>";
  while (bit.hasNext()) cout << " " << bit.next();
  cout << endl;
  */

  Stack<TermList>::Iterator tit(_termLists);
  cout << "Terms:";
  if   (!tit.hasNext()) cout << " <empty>";
  while (tit.hasNext()) cout << " " << tit.next().toString();
  cout << endl;

  Stack<Formula*>::Iterator fit(_formulas);
  cout << "Formulas:";
  if   (!fit.hasNext()) cout << " <empty>";
  while (fit.hasNext()) cout << " " << fit.next()->toString();
  cout << endl;

  /*
  Stack<Formula::VarList*>::Iterator vlit(_varLists);
  cout << "Var lists:";
  if   (!vlit.hasNext()) cout << " <empty>";
  while (vlit.hasNext()) {
    Formula::VarList::Iterator vit(vlit.next());
    if (!vit.hasNext()) {
      cout << " <empty>";
    } else {
      cout << " [";
      while (vit.hasNext()) {
        cout << vit.next();
        if (vit.hasNext()) cout << " ";
      };
      cout << "]";
    }
  }
  cout << endl;

  Map<int, SortList*>::Iterator vsit(_variableSorts);
  cout << "Variables sorts:";
  if   (!vsit.hasNext()) cout << "<empty>";
  int vsitKey;
  SortList* vsitVal;
  while (vsit.hasNext()) {
    vsit.next(vsitKey, vsitVal);
    cout << " {" << vsitKey << " ->";
    SortList::Iterator slit(vsitVal);
    if   (!slit.hasNext()) cout << " <empty>";
    while (slit.hasNext()) cout << " " << (slit.next()).toString();
    cout << "}";
  }
  cout << endl;

  Stack<SortList*>::Iterator slsit(_sortLists);
  cout << "Sort lists: ";
  if   (!slsit.hasNext()) cout << "<empty>";
  while (slsit.hasNext()) {
    SortList* sl = slsit.next();
    SortList::Iterator slit(sl);
    if   (!slit.hasNext()) cout << "<empty>";
    while (slit.hasNext()) cout << (slit.next()).toString() << " ";
    cout << ";";
  }
  cout << endl;

  Stack<TheoryFunction>::Iterator tfit(_theoryFunctions);
  cout << "Theory functions: ";
  if   (!tfit.hasNext()) cout << " <empty>";
  while (tfit.hasNext()) cout << " " << tfit.next();
  cout << endl;

  Stack<LetSymbols>::Iterator lfsit(_letSymbols);
  cout << "Let functions scopes: ";
  if (!lfsit.hasNext()) cout << "<empty>";
  while (lfsit.hasNext()) {
    LetSymbols lfs = lfsit.next();
    LetSymbols::Iterator sit(lfs);
    if (!sit.hasNext()) {
      cout << "<empty>";
    } else {
      unsigned i = lfs.length();
      while (sit.hasNext()) {
        LetSymbol f    = sit.next();
        std::string name     = f.first.first;
        unsigned arity   = f.first.second;
        unsigned symbol  = f.second.first;
        bool isPredicate = f.second.second;

        std::string symbolName = isPredicate ? env.signature->predicateName(symbol)
                                         : env.signature->functionName (symbol);

        cout << name << "/" << arity << " -> " << symbolName;
        if (--i > 0) {
          cout << ", ";
        }
      };
    }
  }
  cout << endl;

  Stack<LetSymbols>::Iterator clfsit(_letTypedSymbols);
  cout << "Current let functions scopes: ";
  if (!clfsit.hasNext()) cout << "<empty>";
  while (clfsit.hasNext()) {
    LetSymbols lfs = clfsit.next();
    LetSymbols::Iterator csit(lfs);
    if (!csit.hasNext()) {
      cout << "<empty>";
    } else {
      unsigned i = lfs.length();
      while (csit.hasNext()) {
        LetSymbol f    = csit.next();
        std::string name     = f.first.first;
        unsigned arity   = f.first.second;
        unsigned symbol  = f.second.first;
        bool isPredicate = f.second.second;

        std::string symbolName = isPredicate ? env.signature->predicateName(symbol)
                                         : env.signature->functionName (symbol);

        cout << name << "/" << arity << " -> " << symbolName;
        if (--i > 0) {
          cout << ", ";
        }
      };
    }
  }
  cout << endl;

  Stack<LetDefinitions>::Iterator lbsit(_letDefinitions);
  cout << "Let definitions: ";
  if (!lbsit.hasNext()) cout << "<empty>";
  while (lbsit.hasNext()) {
    LetDefinitions lbs = lbsit.next();
    LetDefinitions::Iterator lbit(lbs);
    if (!lbit.hasNext()) {
      cout << "<empty>";
    } else {
      cout << "[";
      unsigned i = (unsigned)lbs.length();
      while (lbit.hasNext()) {
        LetSymbolReference ref = lbit.next();
        unsigned symbol = SYMBOL(ref);
        bool isPredicate = IS_PREDICATE(ref);
        std::string symbolName = isPredicate ? env.signature->predicateName(symbol)
                                         : env.signature->functionName (symbol);
        cout << symbolName;
        if (--i > 0) {
          cout << ", ";
        }
      }
      cout << "]";
    }
  }*/
  cout << endl;
}
#endif
