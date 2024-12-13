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
 *  @file Formula.cpp
 *  Implements class Formula.
 */

#include "Formula.hpp"

#include "Clause.hpp"
#include "SubformulaIterator.hpp"
#include "FormulaVarIterator.hpp"

namespace Kernel {

using namespace std;

std::string Formula::DEFAULT_LABEL = "none";

/**
 * Destroy the content of the formula. The destruction depends on the type
 * of this formula.
 *
 * @since 11/12/2004 Manchester, true and false added
 * @since 02/06/2007 Manchester, rewritten for new data types
 */
void Formula::destroy ()
{
  switch ( connective() ) {
  case LITERAL:
    delete static_cast<AtomicFormula*>(this);
    return;

  case AND:
  case OR:
    delete static_cast<JunctionFormula*>(this);
    return;

  case IMP:
  case IFF:
  case XOR:
    delete static_cast<BinaryFormula*>(this);
    return;

  case NOT:
    delete static_cast<NegatedFormula*>(this);
    return;

  case FORALL:
  case EXISTS:
    delete static_cast<QuantifiedFormula*>(this);
    return;

  case BOOL_TERM:
    delete static_cast<BoolTermFormula*>(this);

  case TRUE:
  case FALSE:
    delete this;
    return;

  case NAME:
    delete static_cast<NamedFormula*>(this);
    return;

  case NOCONN:
    ASSERTION_VIOLATION;
  }
}

/**
 * Convert the connective to a std::string.
 * @since 02/01/2004 Manchester
 */
std::string Formula::toString (Connective c)
{
  static std::string names [] =
    { "", "&", "|", "=>", "<=>", "<~>", "~", "!", "?", "$var", "$false", "$true","",""};
  ASS_EQ(sizeof(names)/sizeof(std::string), NOCONN+1);

  return names[(int)c];
} // Formula::toString (Connective c)

/**
 * Convert the formula to a std::string
 *
 * @since 12/10/2002 Tbilisi, implemented as ostream output function
 * @since 09/12/2003 Manchester
 * @since 11/12/2004 Manchester, true and false added
 */
std::string Formula::toString () const
{
  std::string res;

  // render a connective if specified, and then a Formula (or ")" of formula is nullptr)
  typedef struct {
    bool wrapInParenthesis;
    Connective renderConnective; // NOCONN means ""
    const Formula* theFormula;   // nullptr means render ")" instead
  } Todo;

  Stack<Todo> stack;
  stack.push({false,NOCONN,this});

  while (stack.isNonEmpty()) {
    Todo todo = stack.pop();

    // in any case start by rendering the connective passed from "above"
    {
      std::string con = toString(todo.renderConnective);
      if (con != "") {
        res += " "+con+" ";
      }
    }

    const Formula* f = todo.theFormula;

    if (!f) {
      res += ")";
      continue;
    }

    if (todo.wrapInParenthesis) {
      res += "(";
      stack.push({false,NOCONN,nullptr}); // render the final closing bracket
    }

    Connective c = f->connective();
    switch (c) {
    case NAME:
      res += static_cast<const NamedFormula*>(f)->name();
      continue;
    case LITERAL:
      res += f->literal()->toString();
      continue;

    case AND:
    case OR:
      {
        // we will reverse the order
        // but that should not matter

        const FormulaList* fs = f->args();
        ASS (FormulaList::length(fs) >= 2);

        while (FormulaList::isNonEmpty(fs)) {
          const Formula* arg = fs->head();
          fs = fs->tail();
          // the last argument, which will be printed first, is the only one not preceded by a rendering of con
          stack.push({arg->parenthesesRequired(c),FormulaList::isNonEmpty(fs) ? c : NOCONN,arg});
        }

        continue;
      }

    case IMP:
    case IFF:
    case XOR:
      // here we can afford to keep the order right

      stack.push({f->right()->parenthesesRequired(c),c,f->right()});      // second argument with con
      stack.push({f->left()->parenthesesRequired(c),NOCONN,f->left()}); // first argument without con

      continue;

    case NOT:
      {
        res += toString(c);

        const Formula* arg = f->uarg();
        stack.push({arg->parenthesesRequired(c),NOCONN,arg});

        continue;
      }
    case FORALL:
    case EXISTS:
      {
        res += toString(c) + " [";
        VList::Iterator vs(f->vars());
        SList::Iterator ss(f->sorts());
        bool hasSorts = f->sorts();
        bool first=true;
        while (vs.hasNext()) {
          int var = vs.next();
          if (!first) {
            res += ",";
          }
          res += Term::variableToString(var);
          TermList t;
          if (hasSorts) {
            ASS(ss.hasNext());
            t = ss.next();
            if (t != AtomicSort::defaultSort()) {
              res += " : " + t.toString();
            }
          } else if (SortHelper::tryGetVariableSort(var, const_cast<Formula*>(f),t) && t != AtomicSort::defaultSort()) {
            res += " : " + t.toString();
          }
          first = false;
        }
        res += "] : ";

        const Formula* arg = f->qarg();
        stack.push({arg->parenthesesRequired(c),NOCONN,arg});

        continue;
      }

    case BOOL_TERM: {
      std::string term = f->getBooleanTerm().toString();
      res += env.options->showFOOL() ? "$formula{" + term + "}" : term;

      continue;
    }

    case TRUE:
    case FALSE:
      res += toString(c);
      continue;

    case NOCONN:
      ASSERTION_VIOLATION;
  }
  }

  return res;
} // Formula::toString

/**
 * True if the formula needs parentheses around itself
 * when in scope of outer.
 *
 * @since 21/09/2002 Manchester
 * @since 11/12/2004 Manchester, true and false added
 */
bool Formula::parenthesesRequired (Connective outer) const
{
  switch (connective())
    {
    case LITERAL:
    case NOT:
    case FORALL:
    case EXISTS:
    case BOOL_TERM:
    case TRUE:
    case FALSE:
    case NAME:
      return false;

    case OR:
    case AND:
    case IMP:
    case IFF:
    case XOR:
      return true;

    case NOCONN:
      ASSERTION_VIOLATION;
    }

  ASSERTION_VIOLATION;
} // Formula::parenthesesRequired

/**
 * Return the list of all bound variables of the formula
 *
 * If a variable is bound multiple times in the formula,
 * it appears in the list the same number of times as well.
 */
VList* Formula::boundVariables () const
{
  VList* res = VList::empty();
  SubformulaIterator sfit(const_cast<Formula*>(this));
  while(sfit.hasNext()) {
    Formula* sf = sfit.next();
    if(sf->connective() == FORALL || sf->connective() == EXISTS) {
      VList* qvars = sf->vars();
      VList* qvCopy = VList::copy(qvars);
      res = VList::concat(qvCopy, res);
    }
  }
  return res;
}

/**
 * Compute the weight of the formula: the number of connectives plus the
 * weight of all atoms.
 * @since 23/03/2008 Torrevieja
 */
unsigned Formula::weight() const
{
  unsigned result=0;

  SubformulaIterator fs(const_cast<Formula*>(this));
  while (fs.hasNext()) {
    const Formula* f = fs.next();
    switch (f->connective()) {
    case LITERAL:
      result += f->literal()->weight();
      break;
    default:
      result++;
      break;
    }
  }
  return result;
} // Formula::weight

Formula* JunctionFormula::generalJunction(Connective c, FormulaList* args)
{
  if(!args) {
    if(c==AND) {
      return new Formula(true);
    }
    else {
      ASS_EQ(c,OR);
      return new Formula(false);
    }
  }
  if(!args->tail()) {
    return FormulaList::pop(args);
  }
  return new JunctionFormula(c, args);
}

/**
 * Return color of the formula
 *
 * We do not store the color of the formula, so it gets
 * computed again each time the function is called.
 */
Color Formula::getColor()
{
  SubformulaIterator si(this);
  while(si.hasNext()) {
    Formula* f=si.next();
    if(f->connective()!=LITERAL) {
      continue;
    }

    if(f->literal()->color()!=COLOR_TRANSPARENT) {
      return f->literal()->color();
    }
  }
  return COLOR_TRANSPARENT;
}

/**
 * Return true iff the formula is Skip for the purpose of
 * symbol elimination
 */
bool Formula::getSkip()
{
  SubformulaIterator si(this);
  while(si.hasNext()) {
    Formula* f=si.next();
    if(f->connective()!=LITERAL) {
      continue;
    }

    if(!f->literal()->skip()) {
      return false;
    }
  }
  return true;
}

Formula* Formula::trueFormula()
{
  static Formula* res = new Formula(true);
  return res;
}

Formula* Formula::falseFormula()
{
  static Formula* res = new Formula(false);
  return res;
}

/**
 * Creates a formula of the form $ite(c, a, b), where a, b, c are formulas
 * @since 16/04/2015 Gothenburg
 */
Formula* Formula::createITE(Formula* condition, Formula* thenArg, Formula* elseArg)
{
  TermList thenTerm(Term::createFormula(thenArg));
  TermList elseTerm(Term::createFormula(elseArg));
  TermList iteTerm(Term::createITE(condition, thenTerm, elseTerm, AtomicSort::boolSort()));
  return new BoolTermFormula(iteTerm);
}

/**
 * Creates a formula of the form $let(lhs := rhs, body), where body is a formula
 * and lhs and rhs form a binding for a function
 * @since 16/04/2015 Gothenburg
 */
Formula* Formula::createLet(unsigned functor, VList* variables, TermList body, Formula* contents)
{
  TermList contentsTerm(Term::createFormula(contents));
  TermList letTerm(Term::createLet(functor, variables, body, contentsTerm, AtomicSort::boolSort()));
  return new BoolTermFormula(letTerm);
}

/**
 * Creates a formula of the form $let(lhs := rhs, body), where body is a formula
 * and lhs and rhs form a binding for a predicate
 * @since 16/04/2015 Gothenburg
 */
Formula* Formula::createLet(unsigned predicate, VList* variables, Formula* body, Formula* contents)
{
  TermList bodyTerm(Term::createFormula(body));
  TermList contentsTerm(Term::createFormula(contents));
  TermList letTerm(Term::createLet(predicate, variables, bodyTerm, contentsTerm, AtomicSort::boolSort()));
  return new BoolTermFormula(letTerm);
}

Formula* Formula::quantify(Formula* f)
{

  DHMap<unsigned,TermList> tMap;
  SortHelper::collectVariableSorts(f,tMap,/*ignoreBound=*/true);

  //we have to quantify the formula
  VList::FIFO quantifiedVars;
  SList::FIFO theirSorts;

  DHMap<unsigned,TermList>::Iterator tmit(tMap);
  while(tmit.hasNext()) {
    unsigned v;
    TermList s;
    tmit.next(v, s);
    if(s.isTerm() && s.term()->isSuper()){
      // type variable must appear at the start of the list
      quantifiedVars.pushFront(v);
      theirSorts.pushFront(s);
    } else {
      quantifiedVars.pushBack(v);
      theirSorts.pushBack(s);
    }
  }
  if(!quantifiedVars.empty()) {
    f = new QuantifiedFormula(FORALL, quantifiedVars.list(), theirSorts.list(), f);
  }
  return f;
}


/**
 * Return formula equal to @b cl
 * that has all variables quantified
 */
Formula* Formula::fromClause(Clause* cl)
{
  FormulaList* resLst=0;
  unsigned clen=cl->length();
  for(unsigned i=0;i<clen;i++) {
    Formula* lf=new AtomicFormula((*cl)[i]);
    FormulaList::push(lf, resLst);
  }

  Formula* res=JunctionFormula::generalJunction(OR, resLst);
  return Formula::quantify(res);
}

std::ostream& operator<< (std::ostream& out, const Formula& f)
{
  return out << f.toString();
}

std::ostream& operator<< (std::ostream& out, const Formula* f)
{
  return out << f->toString();
}

}
