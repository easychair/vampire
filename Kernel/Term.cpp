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
 * @file Term.cpp
 * Implements class Term.
 *
 * @since 18/04/2006 Bellevue
 * @since 06/05/2007 Manchester, changed into a single class instead of three
 */

#include "Lib/Output.hpp"
#include "Indexing/TermSharing.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/StringUtils.hpp"

#include "SubstHelper.hpp"
#include "TermIterators.hpp"
#include "RobSubstitution.hpp"
#include "Lib/Metaiterators.hpp"

#include "Term.hpp"

using namespace std;
using namespace Lib;
using namespace Kernel;

constexpr unsigned Term::SPECIAL_FUNCTOR_LOWER_BOUND;

void Term::setId(unsigned id)
{
  if (env.options->randomTraversals() &&
      Random::seed() != 1) { // not until a proper seed has been set (i.e. after parsing!)
      // (cf ProvingHelper::runVampire and getPreprocessedProblem in vampire.cpp)
    id += Random::getInteger(1 << 12) << 20; // the twelve most significant bits are randomized
  }
   _args[0]._setId(id);
}

/**
 * Allocate enough bytes to fit a term of a given arity.
 * @since 01/05/2006 Bellevue
 */
void* Term::operator new(size_t,unsigned arity, size_t preData)
{
  //preData must be a multiple of pointer size to maintain alignment
  ASS_EQ(preData%sizeof(size_t), 0);

  size_t sz = sizeof(Term)+arity*sizeof(TermList)+preData;
  void* mem = ALLOC_KNOWN(sz,"Term");
  mem = reinterpret_cast<void*>(reinterpret_cast<char*>(mem)+preData);
  return (Term*)mem;
} // Term::operator new


inline bool argSafeToShare(TermList t)
{ return t.isSafe() && !t.isSpecialVar(); }

/**
 * Destroy the term.
 * @since 01/05/2006 Bellevue
 * @since 07/06/2007 Manchester, changed to new data structures
 */
void Term::destroy ()
{
  ASS(CHECK_LEAKS || ! shared());

  size_t sz = sizeof(Term)+_arity*sizeof(TermList)+getPreDataSize();
  void* mem = this;
  mem = reinterpret_cast<void*>(reinterpret_cast<char*>(mem)-getPreDataSize());
  DEALLOC_KNOWN(mem,sz,"Term");
} // Term::destroy

/**
 * If the term is not shared, destroy it and all its nonshared subterms.
 */
void Term::destroyNonShared()
{
  if (shared()) {
    return;
  }
  TermList selfRef;
  selfRef.setTerm(this);
  TermList* ts=&selfRef;
  static Stack<TermList*> stack(4);
  static Stack<Term*> deletingStack(8);

  for(;;) {
    if (ts->tag()==REF && !ts->term()->shared()) {
      stack.push(ts->term()->args());
      deletingStack.push(ts->term());
    }
    if (stack.isEmpty()) {
      break;
    }
    ts=stack.pop();
    if (!ts->next()->isEmpty()) {
      stack.push(ts->next());
    }
  }
  while (!deletingStack.isEmpty()) {
    deletingStack.pop()->destroy();
  }
}

/**
 * Return true if the term does not contain any unshared proper term.
 *
 * Not containing an unshared term also means that there are no
 * if-then-else or let...in expressions.
 */
bool TermList::isSafe() const
{
  return isVar() || term()->shared();
}

bool TermList::ground() const 
{ return isTerm() && term()->ground(); }

/**
 * Return true if @b ss and @b tt have the same top symbols, that is,
 * either both are the same variable or both are complex terms with the
 * same function symbol.
 * @since 16/08/2008 flight Sydney-San Francisco
 */
bool TermList::sameTop(TermList ss,TermList tt)
{
  if (ss.isVar()) {
    return ss==tt;
  }
  if (tt.isVar()) {
    return false;
  }
  return ss.term()->functor() == tt.term()->functor();
}
void TermList::Top::output(std::ostream& out) const
{
  if (this->var()) {
    out << TermList::var(*this->var());
  } else {
    ASS(this->functor())
    auto f = *this->functor();
    switch (f.kind) {
      case TermKind::LITERAL: out << *env.signature->getPredicate(f.functor);
      case TermKind::TERM:    out << *env.signature->getFunction(f.functor);
      case TermKind::SORT:    out << *env.signature->getTypeCon(f.functor);
    }
  }
}

/**
 * Return true if @b ss and @b tt are both complex terms with the
 * same function symbol.
 */
bool TermList::sameTopFunctor(TermList ss, TermList tt)
{
  if (!ss.isTerm() || !tt.isTerm()) {
    return false;
  }
  return ss.term()->functor() == tt.term()->functor();
}

/**
 * Return true if @b ss and @b tt are both complex terms with the
 * same function symbol.
 */
bool TermList::equals(TermList t1, TermList t2)
{
  static Stack<TermList*> stack(8);
  ASS(stack.isEmpty());

  TermList* ss=&t1;
  TermList* tt=&t2;
  for(;;) {
    if (ss->isTerm() && tt->isTerm() && (!ss->term()->shared() || !tt->term()->shared())) {
      Term* s=ss->term();
      Term* t=tt->term();
      if (s->functor()!=t->functor()) {
        stack.reset();
        return false;
      }
      stack.push(s->args());
      stack.push(t->args());
    }
    else if (ss->content()!=tt->content()) {
      stack.reset();
      return false;
    }

    if (stack.isEmpty()) {
      break;
    }
    tt=stack.pop();
    ss=stack.pop();
    if (!tt->next()->isEmpty()) {
      stack.push(ss->next());
      stack.push(tt->next());
    }
  }
  return true;
}

/**
 * Return true if all proper terms in the @ args list are shared
 */
bool TermList::allShared(TermList* args)
{
  while (args->isNonEmpty()) {
    if (args->isTerm() && !args->term()->shared()) {
      return false;
    }
    args = args->next();
  }
  return true;
}

unsigned TermList::weight() const
{
  return isVar() ? 1 : term()->weight();
}

bool TermList::isArrowSort()
{
  return !isVar() && term()->isSort() &&
         static_cast<AtomicSort*>(term())->isArrowSort();
}

bool TermList::isBoolSort()
{
  return !isVar() && term()->isSort() &&
         static_cast<AtomicSort*>(term())->isBoolSort();
}

bool TermList::isArraySort()
{
  return !isVar() && term()->isSort() &&
         static_cast<AtomicSort*>(term())->isArraySort();
}

bool TermList::isTupleSort()
{
  return !isVar() && term()->isSort() &&
         static_cast<AtomicSort*>(term())->isTupleSort();
}

bool AtomicSort::isArrowSort() const {
  return env.signature->isArrowCon(_functor);
}

bool AtomicSort::isBoolSort() const {
  return env.signature->isBoolCon(_functor);
}

bool AtomicSort::isArraySort() const {
  return env.signature->isArrayCon(_functor);
}

bool AtomicSort::isTupleSort() const {
  return env.signature->isTupleCon(_functor);
}

bool TermList::isApplication() const {
  return !isVar() && term()->isApplication();
}

bool Term::isApplication() const {
  return !isSort() && !isLiteral() && env.signature->isAppFun(_functor);
}

unsigned Term::numTypeArguments() const {
  ASS(!isSort());

  return isSpecial()
    ? 0
    : isLiteral()
      ? env.signature->getPredicate(_functor)->numTypeArguments()
      : env.signature->getFunction(_functor)->numTypeArguments();
}

TermList* Term::termArgs()
{
  ASS(!isSort());

  return _args + (_arity - numTypeArguments());
}

const TermList* Term::typeArgs() const
{ return numTypeArguments() == 0 ? nullptr : args(); }

unsigned Term::numTermArguments() const
{
  if(isSuper() || isSort())
    return 0;

  ASS(_arity >= numTypeArguments())
  return _arity - numTypeArguments();
}

bool TermList::containsSubterm(TermList trm) const
{
  if (!isTerm()) {
    return trm==*this;
  }
  return term()->containsSubterm(trm);
}

bool Term::containsSubterm(TermList trm) const
{
  ASS(!trm.isTerm() || trm.term()->shared());
  ASS(shared());

  if (trm.isTerm() && trm.term()==this) {
    ASS(!isLiteral());
    return true;
  }
  if (arity()==0) {
    return false;
  }

  const TermList* ts=args();
  static Stack<const TermList*> stack(4);
  stack.reset();
  for(;;) {
    if (*ts==trm) {
      return true;
    }
    if (!ts->next()->isEmpty()) {
      stack.push(ts->next());
    }
    if (ts->isTerm()) {
      ASSERT_VALID(*ts->term());
      if (ts->term()->arity()) {
	stack.push(ts->term()->args());
      }
    }
    if (stack.isEmpty()) {
      return false;
    }
    ts=stack.pop();
  }
}

size_t Term::countSubtermOccurrences(TermList subterm) {
  size_t res = 0;

  unsigned stWeight = subterm.isTerm() ? subterm.term()->weight() : 1;
  SubtermIterator stit(this);
  while(stit.hasNext()) {
    TermList t = stit.next();
    if(t==subterm) {
      res++;
      stit.right();
    }
    else if(t.isTerm()) {
      if(t.term()->weight()<=stWeight) {
        stit.right();
      }
    }
  }
  return res;
}

bool TermList::containsAllVariablesOf(TermList t) const
{
  Set<TermList> vars;
  TermIterator oldVars=Term::getVariableIterator(*this);
  while (oldVars.hasNext()) {
    vars.insert(oldVars.next());
  }
  TermIterator newVars=Term::getVariableIterator(t);
  while (newVars.hasNext()) {
    if (!vars.contains(newVars.next())) {
      return false;
    }
  }
  return true;
}

bool Term::containsAllVariablesOf(Term* t)
{
  static DHSet<TermList> vars;
  vars.reset();

  static VariableIterator vit;

  //collect own vars
  vit.reset(this);
  while (vit.hasNext()) {
    vars.insert(vit.next());
  }

  //check t's vars are among collected
  vit.reset(t);
  while (vit.hasNext()) {
    if (!vars.contains(vit.next())) {
      return false;
    }
  }
  return true;
}

bool Term::isShallow() const
{
  const TermList* t = args();
  while (!t->isEmpty()) {
    if (t->isTerm() && t->term()->arity()>0) {
      return false;
    }
    t = t->next();
  }
  return true;
}

TermIterator Term::getVariableIterator(TermList tl)
{
  if (tl.isVar()) {
    return pvi( getSingletonIterator(tl) );
  }
  ASS(tl.isTerm());
  return vi( new VariableIterator(tl.term()) );
}


/**
 * Return the string representation of variable var.
 * @since 16/05/2007
 */
std::string Term::variableToString(unsigned var)
{
  return (std::string)"X" + Int::toString(var);
} // variableToString

/**
 * Return the string representation of variable term var.
 * @since 16/05/2007
 */
std::string Term::variableToString(TermList var)
{
  ASS(var.isVar());

  if (var.isOrdinaryVar()) {
    return (std::string)"X" + Int::toString(var.var());
  }
  else {
    return (std::string)"S" + Int::toString(var.var());
  }
} // variableToString

/**
 * Return the std::string representation of the terms "head"
 * i.e., the function / predicate symbol name or the special term head.
 * Special term prints also '(' and the following arguments which are not args() and a comma
 * Normal term prints "(" if there are any args to follow
 */
std::string Term::headToString() const
{
  if (isSpecial()) {
    const Term::SpecialTermData* sd = getSpecialData();

    switch(specialFunctor()) {
      case SpecialFunctor::FORMULA: {
        ASS_EQ(arity(), 0);
        std::string formula = sd->getFormula()->toString();
        return env.options->showFOOL() ? "$term{" + formula + "}" : formula;
      }
      case SpecialFunctor::LET: {
        ASS_EQ(arity(), 1);
        TermList binding = sd->getBinding();
        bool isPredicate = binding.isTerm() && binding.term()->isBoolean();
        std::string functor = isPredicate ? env.signature->predicateName(sd->getFunctor())
                                      : env.signature->functionName(sd->getFunctor());
        OperatorType* type = isPredicate ? env.signature->getPredicate(sd->getFunctor())->predType()
                                         : env.signature->getFunction(sd->getFunctor())->fnType();

        const VList* variables = sd->getVariables();
        std::string variablesList = "";
        for (unsigned i = 0; i < VList::length(variables); i++) {
          unsigned var = VList::nth(variables, i);
          variablesList += Term::variableToString(var);
          if (i < VList::length(variables) - 1) {
            variablesList += ", ";
          }
        }
        if (VList::length(variables)) {
          variablesList = "(" + variablesList + ")";
        }
        return "$let(" + functor + ": " + type->toString() + ", " + functor + variablesList + " := " + binding.toString() + ", ";
      }
      case SpecialFunctor::ITE: {
        ASS_EQ(arity(),2);
        return "$ite(" + sd->getCondition()->toString() + ", ";
      }
      case SpecialFunctor::TUPLE: {
        ASS_EQ(arity(), 0);
        Term* term = sd->getTupleTerm();
        std::string termList = "";
        Term::Iterator tit(term);
        unsigned i = term->arity();
        while (tit.hasNext()) {
          termList += tit.next().toString();
          if (--i > 0) {
            termList += ", ";
          }
        }
        return "[" + termList + "]";
      }
      case SpecialFunctor::LET_TUPLE: {
        ASS_EQ(arity(), 1);
        VList* symbols = sd->getTupleSymbols();
        unsigned tupleFunctor = sd->getFunctor();
        TermList binding = sd->getBinding();

        OperatorType* fnType = env.signature->getFunction(tupleFunctor)->fnType();

        std::string symbolsList = "";
        std::string typesList = "";
        for (unsigned i = 0; i < VList::length(symbols); i++) {
          Signature::Symbol* symbol = (fnType->arg(i) == AtomicSort::boolSort())
            ? env.signature->getPredicate(VList::nth(symbols, i))
            : env.signature->getFunction(VList::nth(symbols, i));
          symbolsList += symbol->name();
          typesList += symbol->name() + ": " + fnType->arg(i).toString();
          if (i != VList::length(symbols) - 1) {
            symbolsList += ", ";
            typesList += ", ";
          }
        }

        return "$let([" + typesList + "], [" + symbolsList + "] := " + binding.toString() + ", ";
      }
      case SpecialFunctor::LAMBDA: {
        VList* vars = sd->getLambdaVars();
        SList* sorts = sd->getLambdaVarSorts();
        TermList lambdaExp = sd->getLambdaExp();

        std::string varList = "[";

        VList::Iterator vs(vars);
        SList::Iterator ss(sorts);
        bool first = true;
        while(vs.hasNext()) {
          if (!first){
            varList += ", ";
          }else{ first = false; }
          varList += Term::variableToString(vs.next()) + " : ";
          varList += ss.next().toString();
        }
        varList += "]";
        return "(^" + varList + " : (" + lambdaExp.toString() + "))";
      }
      case SpecialFunctor::MATCH: {
        // we simply let the arguments be written out
        return "$match(";
      }
      default:
        ASSERTION_VIOLATION;
    }
  } else {
    unsigned proj;
    if (!isSort() && Theory::tuples()->findProjection(functor(), isLiteral(), proj)) {
      return "$proj(" + Int::toString(proj) + ", ";
    }
    bool print = (isLiteral() || isSort() ||
                 (env.signature->getFunction(_functor)->combinator() == Signature::NOT_COMB)) && arity();
    std::string name = "";
    if(isLiteral()) {
      name = static_cast<const Literal *>(this)->predicateName();
    } else if (isSort()) {
      const AtomicSort* asSort = static_cast<const AtomicSort *>(this);
      if(env.options->showFOOL() && asSort->isBoolSort()){
        name = "$bool";
      } else {
        name = asSort->typeConName();
      }
    } else {
      name = functionName();
    }
    return name + (print ? "(" : "");
  }
}

/**
 * In combination with Term::headToString prepares
 * std::string representation of a term.
 * (this) has to come from arguments of a term of non-zero arity,
 * possibly a special one.
 * Will close the term printed with ')'
 */
std::string TermList::asArgsToString() const
{
  std::string res;

  Stack<const TermList*> stack(64);

  stack.push(this);

  while (stack.isNonEmpty()) {
    const TermList* ts = stack.pop();
    if (! ts) { // comma
      res += ',';
      continue;
    }
    if (ts->isEmpty()) {
      res += ')';
      continue;
    }
    const TermList* tail = ts->next();
    stack.push(tail);
    if (! tail->isEmpty()) {
      stack.push(0);
    }
    if (ts->isVar()) {
      res += Term::variableToString(*ts);
      continue;
    }
    const Term* t = ts->term();

    if(t->isSort() && static_cast<AtomicSort*>(const_cast<Term*>(t))->isArrowSort()){
      res += t->toString();
      continue;
    }

    res += t->headToString();

    if (t->arity()) {
      stack.push(t->args());
    }
  }

  return res;
}

/**
 * Write as a std::string the head of the term list.
 * @since 27/02/2008 Manchester
 */
std::string TermList::toString(bool topLevel) const
{
  if (isEmpty()) {
    return "<empty TermList>";
  }
  if (isVar()) {
    return Term::variableToString(*this);
  }
  return term()->toString(topLevel);
} // TermList::toString


/**
 * Return the result of conversion of a term into a std::string.
 * @since 16/05/2007 Manchester
 */
std::string Term::toString(bool topLevel) const
{
  bool printArgs = true;

  if(isSuper()){
    return "$tType";
  }

  if(!isSpecial() && !isLiteral()){
    if(isSort() && static_cast<AtomicSort*>(const_cast<Term*>(this))->isArrowSort()){
      ASS(arity() == 2);
      std::string res;
      TermList arg1 = *(nthArgument(0));
      TermList arg2 = *(nthArgument(1));
      res += topLevel ? "" : "(";
      res += arg1.toString(false) + " > " + arg2.toString();
      res += topLevel ? "" : ")";
      return res;
    }

    printArgs = isSort() || env.signature->getFunction(_functor)->combinator() == Signature::NOT_COMB;
  }

#if NICE_THEORY_OUTPUT
  auto theoryTerm = Kernel::tryNumTraits([&](auto numTraits) {
    using NumTraits = decltype(numTraits);
    auto uminus = [&]()  {
    std::stringstream out;
      auto maybePar = [&](auto t) { 
        auto needsPar = t.isTerm() && NumTraits::isAdd(t.term()->functor());
        return t.toString(!needsPar);
      };
      out << "-" << maybePar(termArg(0));
      return Option<std::string>(out.str());
    };
    auto binary = [&](auto sym)  {
      auto needsPar = !topLevel;
      auto maybePar = [&](auto t) { 
        return t.toString(
            t.isTerm() && (
              t.term()->functor() == _functor 
              || (NumTraits::isAdd(_functor) && NumTraits::isMul(t.term()->functor()))
              )
            );
      };
      std::stringstream out;
      out << (needsPar ? "(" : "");
      out << maybePar(termArg(0)) << sym << maybePar(termArg(1));
      out << (needsPar ? ")" : "");
      return Option<std::string>(out.str());
    };
    if (isLiteral()) {
      if (NumTraits::isGreater(_functor)) {
        return binary(">");
      } else if (NumTraits::isGeq(_functor)) {
        return binary(">=");
      }
      /* nothing */
    } else if (isSort()) {
      /* nothing */
    } else {
      if (NumTraits::isAdd(_functor)) {
        return binary(" + ");
      } else if (NumTraits::isMul(_functor)) {
        return binary(" ");
      } else if (NumTraits::isFloor(_functor)) {
        return some(outputToString("⌊", termArg(0), "⌋"));
      } else if (NumTraits::isMinus(_functor)) {
        return uminus();
      }
    }
    /* this means we have an uninterpteted term which will be formatted as usual */
    return Option<std::string>();
  });

  if (theoryTerm.isSome()) {
    return theoryTerm.unwrap();
  }
#endif // NICE_THEORY_OUTPUT

  std::stringstream out;
  out << headToString();
  
  if (_arity && printArgs) {
    out << Output::interleaved(',', anyArgIter(this)) << ")";
  }
  return out.str();
} // Term::toString

TermList Literal::eqArgSort() const {
  ASS(isEquality())
  return SortHelper::getEqualityArgumentSort(this);
}

/**
 * Return the result of conversion of a literal into a std::string.
 * @since 16/05/2007 Manchester
 */
std::string Literal::toString() const
{
  if (isEquality()) {
    const TermList* lhs = args();
    std::string s = lhs->toString();
    if (isPositive()) {
      s += " = ";
    }
    else {
      s += " != ";
    }

    std::string res = s + lhs->next()->toString();
    if (env.getMainProblem() == nullptr || env.getMainProblem()->isHigherOrder() || 
       (SortHelper::getEqualityArgumentSort(this) == AtomicSort::boolSort())){
      res = "("+res+")";
    }
    /*if(isTwoVarEquality()){
      res += "___ sort: " + twoVarEqSort().toString();
    }*/

    return res;
  }

#if NICE_THEORY_OUTPUT
  auto theoryTerm = Kernel::tryNumTraits([&](auto numTraits) {
    auto binary = [&](auto sym)  {
    std::stringstream out;
      if (!polarity()) out << "~(" ;
      out << *nthArgument(0) << " " << sym << " " << *nthArgument(1) ;
      if (!polarity()) out << ")" ;
      return Option<std::string>(out.str());
    };
    using NumTraits = decltype(numTraits);
    if (_functor == NumTraits::greaterF()) {
      return binary(">");
    } else if (_functor == NumTraits::geqF()) {
      return binary(">=");
    }
    return Option<std::string>();
  });
  if (theoryTerm.isSome()) {
    return theoryTerm.unwrap();
  }
#endif // NICE_THEORY_OUTPUT


  Stack<const TermList*> stack(64);
  std::string s = polarity() ? "" : "~";
  unsigned proj;
  if (Theory::tuples()->findProjection(functor(), true, proj)) {
    return s + "$proj(" + Int::toString(proj) + ", " + args()->asArgsToString();
  }
  s += predicateName();

  //cerr << "predicate: "<< predicateName()<<endl;
  if (_arity) {
    s += '(' + args()->asArgsToString(); // will also print the ')'
  }
  return s;
} // Literal::toString


/**
 * Return the print name of the function symbol of this term.
 * @since 18/05/2007 Manchester
 */
const std::string& Term::functionName() const
{
#if VDEBUG
  static std::string nonexisting("<function does not exists>");
  if (_functor>=static_cast<unsigned>(env.signature->functions())) {
    return nonexisting;
  }
#endif

  return env.signature->functionName(_functor);
} // Term::functionName

/**
 * Return the print name of the type constructor symbol of this sort.
 */
const std::string& AtomicSort::typeConName() const
{
#if VDEBUG
  static std::string nonexisting("<type constructor does not exists>");
  if (_functor>=static_cast<unsigned>(env.signature->typeCons())) {
    return nonexisting;
  }
#endif

  return env.signature->typeConName(_functor);
} // Term::functionName

/**
 * Return the print name of the function symbol of this literal.
 * @since 18/05/2007 Manchester
 */
const std::string& Literal::predicateName() const
{
#if VDEBUG
  static std::string nonexisting("<predicate does not exists>");
  if (_functor>=static_cast<unsigned>(env.signature->predicates())) {
    return nonexisting;
  }
#endif

  return env.signature->predicateName(_functor);
} // Literal::predicateName


bool Literal::isAnswerLiteral() const {
  return isNegative() && env.signature->getPredicate(functor())->answerPredicate();
}


/**
 * Apply @b subst to the term and return the result.
 * @since 28/12/2007 Manchester
 */
Term* Term::apply(Substitution& subst)
{
  return SubstHelper::apply(this, subst);
} // Term::apply


/**
 * Apply @b subst to the literal and return the result.
 * @since 28/12/2007 Manchester
 */
Literal* Literal::apply(Substitution& subst)
{
  return SubstHelper::apply(this, subst);
} // Literal::apply

/**
 * Return literal opposite to @b l.
 */
Literal* Literal::complementaryLiteral(Literal* l)
{
  Literal* res=env.sharing->tryGetOpposite(l);
  if (!res) {
    res=create(l,!l->polarity());
  }
  return res;
}


/** Create a new complex term, copy from @b t its function symbol and
 *  from the array @b args its arguments. Insert it into the sharing
 *  structure if all arguments are shared.
 * @since 07/01/2008 Torrevieja
 */
Term* Term::create(Term* t,TermList* args)
{
  ASS_EQ(t->getPreDataSize(), 0);
  ASS(!t->isLiteral())
  ASS(!t->isSort())
  return Term::create(t->functor(), t->arity(), args);
}

/** Create a new complex term, and insert it into the sharing
 *  structure if all arguments are shared.
 */
Term* Term::create(unsigned function, unsigned arity, const TermList* args)
{
  ASS_EQ(env.signature->functionArity(function), arity);

  bool share = range(0, arity).all([&](auto i) { return argSafeToShare(args[i]); });

  auto allocTerm = [&]() {
    Term* s = new(arity) Term;
    s->makeSymbol(function,arity);
    for (auto i : range(0, arity)) {
      *s->nthArgument(i) = args[i];
    }
    return s;
  };

  if (share) {
    bool created = false;
    auto shared =
      env.sharing->_terms.rawFindOrInsert(allocTerm,
        Term::termHash(function, [&](auto i){ return args[i]; }, arity),
        [&](Term* t) { return t->functor() == function && range(0, arity).all([&](auto i) { return args[i] == *t->nthArgument(i); }); },
        created);
    if (created) {
      env.sharing->computeAndSetSharedTermData(shared);
    }
    return shared;
  } else {
    return allocTerm();
  }
}


/** Create a new constant and insert in into the sharing
 *  structure.
 */
Term* Term::createConstant(const std::string& name)
{
  unsigned symbolNumber = env.signature->addFunction(name,0);
  return createConstant(symbolNumber);
}

/** Create a new complex term, copy from @b t its function symbol and
 *  from the array @b args its arguments. Do not insert it into the sharing
 *  structure.
 * @since 07/01/2008 Torrevieja
 */
Term* Term::createNonShared(Term* t,TermList* args)
{
  int arity = t->arity();
  Term* s = new(arity) Term(*t);
  TermList* ss = s->args();
  for (int i = 0;i < arity;i++) {
    ASS(!args[i].isEmpty());
    *ss-- = args[i];
  }
  return s;
} // Term::createNonShared(const Term* t,Term* args)



/** Create a new complex term, and do not insert it into the sharing
 *  structure.
 */
Term* Term::createNonShared(unsigned function, unsigned arity, TermList* args)
{
  ASS_EQ(env.signature->functionArity(function), arity);

  Term* s = new(arity) Term;
  s->makeSymbol(function,arity);

  TermList* ss = s->args();

  TermList* curArg = args;
  TermList* argStopper = args+arity;
  while (curArg!=argStopper) {
    *ss = *curArg;
    --ss;
    ++curArg;
  }
  return s;
}

/**
 * Create a (condition ? thenBranch : elseBranch) expression
 * and return the resulting term
 */
Term* Term::createITE(Formula * condition, TermList thenBranch, TermList elseBranch, TermList branchSort)
{
  Term* s = new(2,sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::ITE), 2);
  TermList* ss = s->args();
  *ss = thenBranch;
  ss = ss->next();
  *ss = elseBranch;
  ASS(ss->next()->isEmpty());
  s->getSpecialData()->_iteData.condition = condition;
  s->getSpecialData()->_iteData.sort = branchSort;
  return s;
}

/**
 * Create (let lhs <- rhs in t) expression and return
 * the resulting term
 */
Term* Term::createLet(unsigned functor, VList* variables, TermList binding, TermList body, TermList bodySort)
{
#if VDEBUG
  Set<unsigned> distinctVars;
  VList::Iterator vit(variables);
  while (vit.hasNext()) {
    distinctVars.insert(vit.next());
  }
  ASS_EQ(distinctVars.size(), VList::length(variables));

  bool isPredicate = binding.isTerm() && binding.term()->isBoolean();
  const unsigned int arity = isPredicate ? env.signature->predicateArity(functor)
                                         : env.signature->functionArity(functor);
  ASS_EQ(arity, VList::length(variables));
#endif

  Term* s = new(1,sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::LET), 1);
  TermList* ss = s->args();
  *ss = body;
  ASS(ss->next()->isEmpty());
  s->getSpecialData()->_letData.functor = functor;
  s->getSpecialData()->_letData.variables = variables;
  s->getSpecialData()->_letData.sort = bodySort;
  s->getSpecialData()->_letData.binding = binding;
  return s;
}

/**
 * Create (let [a, b, c] <- rhs in t) expression and return
 * the resulting term
 */
Term* Term::createTupleLet(unsigned tupleFunctor, VList* symbols, TermList binding, TermList body, TermList bodySort)
{
#if VDEBUG
  Signature::Symbol* tupleSymbol = env.signature->getFunction(tupleFunctor);
  ASS_EQ(tupleSymbol->arity(), VList::length(symbols));
  ASS_REP(tupleSymbol->fnType()->result().isTupleSort(), tupleFunctor);

  Set<pair<int,bool> > distinctSymbols;
  VList::Iterator sit(symbols);
  unsigned arg = 0;
  while (sit.hasNext()) {
    unsigned symbol = sit.next();
    bool isPredicate = tupleSymbol->fnType()->arg(arg) == AtomicSort::boolSort();
    if (!distinctSymbols.contains(make_pair(symbol, isPredicate))) {
      distinctSymbols.insert(make_pair(symbol, isPredicate));
    } else {
      ASSERTION_VIOLATION_REP(symbol);
    }
    arg++;
  }
#endif

  Term* s = new(1,sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::LET_TUPLE), 1);
  TermList* ss = s->args();
  *ss = body;
  ASS(ss->next()->isEmpty());
  s->getSpecialData()->_letTupleData.functor = tupleFunctor;
  s->getSpecialData()->_letTupleData.symbols = symbols;
  s->getSpecialData()->_letTupleData.sort = bodySort;
  s->getSpecialData()->_letTupleData.binding = binding;
  return s;
}

/**
 * Create a formula expression and return
 * the resulting term
 */
Term* Term::createFormula(Formula* formula)
{
  Term* s = new(0,sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::FORMULA), 0);
  s->getSpecialData()->_formulaData.formula = formula;
  return s;
}


/**
 * Create a lambda term from a list of lambda vars and an
 * expression and returns the resulting term
 */
Term* Term::createLambda(TermList lambdaExp, VList* vars, SList* sorts, TermList expSort){
  Term* s = new(0, sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::LAMBDA), 0);
  //should store body of lambda in args
  s->getSpecialData()->_lambdaData.lambdaExp = lambdaExp;
  s->getSpecialData()->_lambdaData._vars = vars;
  s->getSpecialData()->_lambdaData._sorts = sorts;
  s->getSpecialData()->_lambdaData.expSort = expSort;
  SList::Iterator sit(sorts);
  Stack<TermList> revSorts;
  TermList lambdaTmSort = expSort;
  while(sit.hasNext()){
    revSorts.push(sit.next());
  }
  while(!revSorts.isEmpty()){
    TermList varSort = revSorts.pop();
    lambdaTmSort = AtomicSort::arrowSort(varSort, lambdaTmSort);
  }
  s->getSpecialData()->_lambdaData.sort = lambdaTmSort;
  return s;
}

Term* Term::createTuple(unsigned arity, TermList* sorts, TermList* elements) {
  unsigned tupleFunctor = Theory::tuples()->getFunctor(arity, sorts);
  Term* tupleTerm = Term::create(tupleFunctor, arity, elements);
  return createTuple(tupleTerm);
}

Term* Term::createTuple(Term* tupleTerm) {
  Term* s = new(0, sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::TUPLE), 0);
  s->getSpecialData()->_tupleData.term = tupleTerm;
  return s;
}

Term *Term::createMatch(TermList sort, TermList matchedSort, unsigned int arity, TermList *elements) {
  Term *s = new (arity, sizeof(SpecialTermData)) Term;
  s->makeSymbol(toNormalFunctor(SpecialFunctor::MATCH), arity);
  TermList *ss = s->args();
  s->getSpecialData()->_matchData.sort = sort;
  s->getSpecialData()->_matchData.matchedSort = matchedSort;

  for (unsigned i = 0; i < arity; i++) {
    ASS(!elements[i].isEmpty());
    *ss = elements[i];
    ss = ss->next();
  }
  ASS(ss->isEmpty());

  return s;
}

/** Create a new complex term, copy from @b t its function symbol and arity.
 *  Initialize its arguments by a dummy special variable.
 */
Term* Term::createNonShared(Term* t)
{
  int arity = t->arity();
  Term* s = new(arity) Term(*t);
  TermList* ss = s->args();
  for (int i = 0;i < arity;i++) {
    (*ss--).makeSpecialVar(0);
  }
  return s;
} // Term::createNonShared(const Term* t)

/** Create clone of complex term @b t. Do not insert it into the sharing
 *  structure.
 */
Term* Term::cloneNonShared(Term* t)
{
  int arity = t->arity();
  TermList* args = t->args();
  Term* s = new(arity) Term(*t);
  TermList* ss = s->args();
  for (int i = 0;i < arity;i++) {
    *ss-- = args[-i];
  }
  return s;
} // Term::cloneNonShared(const Term* t,Term* args)

Term* Term::create1(unsigned fn, TermList arg)
{ return Term::create(fn, { arg }); }

Term* Term::create2(unsigned fn, TermList arg1, TermList arg2)
{ return Term::create(fn, {arg1, arg2}); }


Term* Term::create(unsigned fn, std::initializer_list<TermList> args)
{ return Term::create(fn, args.size(), args.begin()); }

/**
 * Create singleton FOOL constants
 */
Term* Term::foolTrue(){
  static Term* _foolTrue = createConstant(env.signature->getFoolConstantSymbol(true));
  return _foolTrue;
}

Term* Term::foolFalse(){
  static Term* _foolFalse = createConstant(env.signature->getFoolConstantSymbol(false));
  return _foolFalse;
}

/*
 * NOTE: by design the term that represent $tType is not shared
 * and also is not linked to a symbol in the signature.
 */
TermList AtomicSort::superSort(){
  static AtomicSort* _super = createNonSharedConstant(0);
  return TermList(_super);
}

TermList AtomicSort::defaultSort(){
  static AtomicSort* _default = createConstant(env.signature->getDefaultSort());
  return TermList(_default);
}

TermList AtomicSort::boolSort(){
  static AtomicSort* _bool = createConstant(env.signature->getBoolSort());
  return TermList(_bool);
}

TermList AtomicSort::intSort(){
  static AtomicSort* _int = createConstant(env.signature->getIntSort());
  return TermList(_int);
}

TermList AtomicSort::realSort(){
  static AtomicSort* _real = createConstant(env.signature->getRealSort());
  return TermList(_real);
}

TermList AtomicSort::rationalSort(){
  static AtomicSort* _rat = createConstant(env.signature->getRatSort());
  return TermList(_rat);
}

TermList AtomicSort::arrowSort(TermList s1, TermList s2){
  unsigned arrow = env.signature->getArrowConstructor();
  return TermList(create2(arrow, s1, s2));
}

TermList AtomicSort::arrowSort(TermList s1, TermList s2, TermList s3){
  return arrowSort(s1, arrowSort(s2, s3));
}

TermList AtomicSort::arrowSort(TermStack& domSorts, TermList range)
{
  TermList res = range;

  for(unsigned i = 0; i < domSorts.size(); i++){
    res = arrowSort(domSorts[i], res);
  }
  return res;
}

AtomicSort* AtomicSort::createConstant(const std::string& name)
{
  bool added;
  unsigned newSort = env.signature->addTypeCon(name,0,added);
  if(added){
    OperatorType* ot = OperatorType::getConstantsType(superSort());
    env.signature->getTypeCon(newSort)->setType(ot);
  }
  return createConstant(newSort);
}

TermList AtomicSort::arraySort(TermList indexSort, TermList innerSort)
{
  unsigned array = env.signature->getArrayConstructor();
  TermList sort = TermList(create2(array, indexSort, innerSort));
  return sort;
}

TermList AtomicSort::tupleSort(unsigned arity, TermList* sorts)
{ return TermList(AtomicSort::create(env.signature->getTupleConstructor(arity), arity, sorts)); }

unsigned Term::computeDistinctVars() const
{
  Set<unsigned> vars;
  VariableIterator vit(this);
  while (vit.hasNext()) {
    vars.insert(vit.next().var());
  }
  return vars.size();
}

/**
 * True if each function and predicate symbols in this term or literal are
 * marked as skip for the purpose of symbol elimination.
 * @since 04/05/2013 Manchester, changed to use the new NonVariable Iterator
 * @author Andrei Voronkov
 */
bool Term::skip() const
{
  if (isLiteral()) {
    if (!env.signature->getPredicate(functor())->skip()) {
      return false;
    }
  }
  else {
    if (!env.signature->getFunction(functor())->skip()) {
      return false;
    }
  }
  NonVariableIterator nvi(const_cast<Term*>(this));
  while (nvi.hasNext()) {
    unsigned func=nvi.next().term()->functor();
    if (!env.signature->getFunction(func)->skip()) {
      return false;
    }
  }
  return true;
} // skip

bool Term::isBoolean() const {
  const Term* term = this;
  while (true) {
    if (env.signature->isFoolConstantSymbol(true, term->functor()) ||
        env.signature->isFoolConstantSymbol(false, term->functor())) return true;
    if (!term->isSpecial()){
      bool val = !term->isLiteral() &&
      env.signature->getFunction(term->functor())->fnType()->result() == AtomicSort::boolSort();
      return val;
    }
    switch (term->specialFunctor()) {
      case SpecialFunctor::FORMULA:
        return true;
      case SpecialFunctor::TUPLE:
      case SpecialFunctor::LAMBDA:
        return false;
      case SpecialFunctor::ITE:
      case SpecialFunctor::LET:
      case SpecialFunctor::LET_TUPLE: {
        const TermList *ts = term->nthArgument(0);
        if (!ts->isTerm()) {
          return false;
        } else {
          term = ts->term();
          break;
        }
      }
      case SpecialFunctor::MATCH: {
        const TermList *ts = term->nthArgument(2);
        if (!ts->isTerm()) {
          return false;
        } else {
          term = ts->term();
          break;
        }
      }
      default:
        ASSERTION_VIOLATION_REP(term->toString());
    }
  }
  return false;
} // isBoolean

bool Term::isSuper() const {
  return this == AtomicSort::superSort().term();
}

/** Create a new sort, and insert it into the sharing
 *  structure if all arguments are shared.
 */
AtomicSort* AtomicSort::create(unsigned typeCon, unsigned arity, const TermList* args)
{
  ASS_EQ(env.signature->typeConArity(typeCon), arity);

  bool share = range(0, arity).all([&](auto i) { return argSafeToShare(args[i]); });

  auto allocTerm = [&]() {
    AtomicSort* s = new(arity) AtomicSort(typeCon,arity);
    s->makeSymbol(typeCon,arity);
    for (auto i : range(0, arity)) {
      *s->nthArgument(i) = args[i];
    }
    return s;
  };

  if (share) {
    bool created = false;
    auto shared =
      env.sharing->_sorts.rawFindOrInsert(allocTerm,
        Term::termHash(typeCon, [&](auto i){ return args[i]; }, arity),
        [&](AtomicSort* t) { return t->functor() == typeCon && range(0, arity).all([&](auto i) { return args[i] == *t->nthArgument(i); }); },
        created
        );
    if (created) {
      env.sharing->computeAndSetSharedSortData(shared);
    }
    return shared;
  } else {
    return allocTerm();
  }
}

/** Create a new complex sort, copy from @b sort its function symbol and
 *  from the array @b args its arguments. Insert it into the sharing
 *  structure if all arguments are shared.
 * @since 07/01/2008 Torrevieja
 */
AtomicSort* AtomicSort::create(AtomicSort const* sort,TermList* args)
{
  return AtomicSort::create(sort->functor(), sort->arity(), args);
}


AtomicSort* AtomicSort::create2(unsigned tc, TermList arg1, TermList arg2)
{
  TermList args[] = {arg1, arg2};
  return AtomicSort::create(tc, 2, args);
}


/** Create a new complex sort, and do not insert it into the sharing
 *  structure.
 */
AtomicSort* AtomicSort::createNonShared(unsigned typeCon, unsigned arity, TermList* args)
{
  ASS_EQ(env.signature->typeConArity(typeCon), arity);

  AtomicSort* s = new(arity) AtomicSort(typeCon, arity);
  TermList* ss = s->args();

  TermList* curArg = args;
  TermList* argStopper = args+arity;
  while (curArg!=argStopper) {
    *ss = *curArg;
    --ss;
    ++curArg;
  }
  return s;
}

/**
 * Return true iff headers of literals match each other.
 */
bool Literal::headersMatch(Literal* l1, Literal* l2, bool complementary)
{
  if (l1->_functor!=l2->_functor || (complementary?1:0)!=(l1->polarity()!=l2->polarity())) {
    return false;
  }

  return true;
}

/** Create a new literal, and insert it into the sharing
 *  structure if all arguments are shared.
 */
template<class GetArg>
Literal* Literal::create(unsigned predicate, unsigned arity, bool polarity, GetArg getArg, Option<TermList> twoVarEqSort)
{
  ASS(!twoVarEqSort || (predicate == 0 && arity == 2 && getArg(0).isVar() && getArg(1).isVar()))
  ASS(predicate != 0 || arity == 2)
  ASS_EQ(env.signature->predicateArity(predicate), arity);

  bool share = range(0, arity).all([&](auto i) { return argSafeToShare(getArg(i)); });
  bool swapArgs = share && predicate == 0 && Indexing::TermSharing::argNormGt(getArg(0), getArg(1));
  auto normArg = [&](auto i) { return swapArgs ? getArg(1 - i) : getArg(i); };

  auto allocLiteral = [&]() {
    Literal* l = new(arity) Literal(predicate, arity, polarity);
    for (auto i : range(0, arity)) {
      *l->nthArgument(i) = normArg(i);
    }
    if (twoVarEqSort) {
      l->markTwoVarEquality();
      l->setTwoVarEqSort(*twoVarEqSort);
    }
    return l;
  };

  if (share) {
    bool created = false;
    auto shared =
      env.sharing->_literals.rawFindOrInsert(allocLiteral,
        Literal::literalHash(predicate, polarity, normArg, arity, twoVarEqSort),
        [&](Literal* t) { return Literal::literalEquals(t, predicate, polarity, normArg, arity, twoVarEqSort); },
        created);

    if (created) {
      if (twoVarEqSort)
        env.sharing->computeAndSetSharedVarEqData(shared, *twoVarEqSort);
      else
        env.sharing->computeAndSetSharedLiteralData(shared);
    }
    ASS(predicate != 0 || rightArgOrder(*shared->nthArgument(0), *shared->nthArgument(1)))
    return shared;
  } else {
    return allocLiteral();
  }
}

Literal* Literal::create(unsigned predicate, unsigned arity, bool polarity, TermList* args)
{ return create(predicate, arity, polarity, [&](auto i) { return args[i]; }); }

/** Create a new literal, copy from @b l its predicate symbol and
 *  its arguments, and set its polarity to @b polarity. Insert it
 *  into the sharing structure if all arguments are shared.
 * @since 07/01/2008 Torrevieja
 */
Literal* Literal::create(Literal* l,bool polarity)
{
  ASS_EQ(l->getPreDataSize(), 0);

  return l->isEquality()
    ? Literal::createEquality(polarity, *l->nthArgument(0), *l->nthArgument(1), SortHelper::getEqualityArgumentSort(l))
    : Literal::create(l->functor(), l->arity(), polarity, [&](auto i) { return *l->nthArgument(i); });
} // Literal::create

/** Create a new literal, copy from @b l its predicate symbol and
 *  from the array @b args its arguments. Insert it into the sharing
 *  structure if all arguments are shared.
 * @since 07/01/2008 Torrevieja
 */
Literal* Literal::create(Literal* l,TermList* args)
{
  return l->isEquality()
    ? Literal::createEquality(l->polarity(), args[0], args[1], SortHelper::getEqualityArgumentSort(l))
    : Literal::create(l->functor(), l->arity(), l->polarity(), [&](auto i) { return args[i]; });
} // Literal::create


/**
 * Return a new equality literal, with polarity @b polarity and
 * arguments @b arg1 and @b arg2. These arguments must be of sort @c sort
 * (or more specific, in the polymorphic case) unless they are variables.
 * Insert the new literal into the sharing structure if all arguments
 * are shared.
 *
 * The equality may be between two variables.
 */
Literal* Literal::createEquality (bool polarity, TermList arg1, TermList arg2, TermList sort)
{
#if VDEBUG
   TermList srt1, srt2;
   static RobSubstitution checkSortSubst;
   checkSortSubst.reset();

   if (!SortHelper::tryGetResultSort(arg1, srt1)) {
     if (!SortHelper::tryGetResultSort(arg2, srt2)) {
       ASS_REP(arg1.isVar(), arg1.toString());
       ASS_REP(arg2.isVar(), arg2.toString());
     } else{
       ASS(env.sharing->isWellSortednessCheckingDisabled() || checkSortSubst.match(sort, 0, srt2, 1));
     }
   }
   else {
    ASS_REP2(env.sharing->isWellSortednessCheckingDisabled() || checkSortSubst.match(sort, 0, srt1, 1), sort.toString(), srt1.toString());
     if (SortHelper::tryGetResultSort(arg2, srt2)) {
       checkSortSubst.reset();
       ASS_REP2(env.sharing->isWellSortednessCheckingDisabled() || checkSortSubst.match(sort, 0, srt2, 1), sort.toString(), arg2.toString() + " :  " + srt2.toString());
     }
   }
#endif // VDEBUG

   auto getArg = [&](auto i) { ASS_L(i, 2); return i == 0 ? arg1 : arg2; };
   return Literal::create(/* predicate */ 0, /* arity */ 2, polarity, getArg, someIf(arg1.isVar() && arg2.isVar(), [&](){ return sort; }));
}

Literal* Literal::create(unsigned predicate, bool polarity, std::initializer_list<TermList> args)
{ return Literal::create(predicate, args.size(), polarity, [&](auto i) { return args.begin()[i]; }); }

Literal* Literal::create1(unsigned predicate, bool polarity, TermList arg)
{ return Literal::create(predicate, polarity, { arg }); }

Literal* Literal::create2(unsigned predicate, bool polarity, TermList arg1, TermList arg2)
{ return Literal::create(predicate, polarity, { arg1, arg2 }); }

/** create a new term and copy from t the relevant part of t's content */
Term::Term(const Term& t) throw()
  : _functor(t._functor),
    _arity(t._arity),
    _color(COLOR_TRANSPARENT),
    _hasInterpretedConstants(0),
    _isTwoVarEquality(0),
    _weight(0),
    _kboWeight(-1),
#if VDEBUG
    _kboInstance(nullptr),
#endif
    _vars(0)
{
  ASS(!isSpecial()); //we do not copy special terms

  _args[0] = t._args[0];
  _args[0]._setShared(false);
  _args[0]._setOrder(AO_UNKNOWN);
  _args[0]._setDistinctVars(TERM_DIST_VAR_UNKNOWN);
} // Term::Term

/** create a new literal and copy from l its content */
Literal::Literal(const Literal& l) throw()
  : Term(l)
{
}

/** create a new AtomicSort and copy from l its content */
AtomicSort::AtomicSort(const AtomicSort& p) throw()
  : Term(p)
{
}

/** dummy term constructor */
Term::Term() throw()
  :_functor(0),
   _arity(0),
   _color(COLOR_TRANSPARENT),
   _hasInterpretedConstants(0),
   _isTwoVarEquality(0),
   _weight(0),
   _kboWeight(-1),
#if VDEBUG
   _kboInstance(nullptr),
#endif
   _maxRedLen(0),
   _vars(0)
{
  _args[0].setContent(0);
  _args[0]._setTag(FUN);
  _args[0]._setDistinctVars(TERM_DIST_VAR_UNKNOWN);
} // Term::Term

Literal::Literal()
{
}

bool Literal::computable() const {
  if (!env.signature->getPredicate(this->functor())->computable()) {
    return false;
  }
  for (unsigned i = 0; i < arity(); ++i) {
    const TermList* t = nthArgument(i);
    if (!t->isTerm() || !t->term()->computable()) {
      return false;
    }
  }
  return true;
}

bool Literal::computableOrVar() const {
  if (!env.signature->getPredicate(this->functor())->computable()) {
    return false;
  }
  for (unsigned i = 0; i < arity(); ++i) {
    const TermList* t = nthArgument(i);
    if (t->isTerm() && !t->term()->computableOrVar()) {
      return false;
    }
  }
  return true;
}

AtomicSort::AtomicSort()
{
}

#if VDEBUG
std::string Term::headerToString() const
{
  std::string s("functor: ");
  s += Int::toString(_functor) + ", arity: " + Int::toString(_arity)
    + ", weight: " + Int::toString(_weight)
    + ", vars: " + Int::toString(_vars)
    + ", polarity: " + Int::toString(_args[0]._polarity())
    + ", shared: " + Int::toString(_args[0]._shared())
    + ", literal: " + Int::toString(_args[0]._literal())
    + ", order: " + Int::toString(_args[0]._order())
    + ", tag: " + Int::toString(_args[0]._tag());
  return s;
}

void Term::assertValid() const
{
  ASS_ALLOC_TYPE(this, "Term");
  ASS_EQ(_args[0]._tag(), FUN);
}

void TermList::assertValid() const
{
  if (this->isTerm()) {
    ASS_ALLOC_TYPE(_term, "Term");
    ASS_EQ(_term()->_args[0]._tag(), FUN);
  }
}

#endif

std::ostream& Kernel::operator<<(std::ostream& out, TermList const& tl)
{
  if (tl.isEmpty()) {
    return out<<"<empty TermList>";
  }
  if (tl.isVar()) {
    return out<<Term::variableToString(tl);
  }
  return out << *tl.term();
}

std::ostream& Kernel::operator<<(std::ostream& out, const Term& t)
{
  return out<<t.toString();
}
std::ostream& Kernel::operator<<(std::ostream& out, const Literal& l)
{
  return out<<l.toString();
}

bool Kernel::operator<(const TermList& lhs, const TermList& rhs)
{
  auto cmp = lhs.isTerm() - rhs.isTerm();
  if (cmp != 0) return cmp < 0;
  if (lhs.isTerm()) {
    ASS(rhs.isTerm())
    return lhs.term()->getId() < rhs.term()->getId();
  } else if (lhs.isEmpty() || rhs.isEmpty()) {
    auto cmp = lhs.isEmpty() - rhs.isEmpty();
    if (cmp != 0) return cmp < 0;
    else return false;

  } else {
    ASS(lhs.isVar())
    ASS(rhs.isVar())
    return lhs.var() < rhs.var();
  }
}

bool Literal::rightArgOrder(TermList const& lhs, TermList const& rhs)
{ return !Indexing::TermSharing::argNormGt(lhs,rhs); }

bool Kernel::positionIn(TermList& subterm,TermList* term,std::string& position)
{
   //cout << "positionIn " << subterm.toString() << " in " << term->toString() << endl;

  if(!term->isTerm()){
    if(subterm.isTerm()) return false;
    if (term->var()==subterm.var()){
      position = "1";
      return true;
    }
    return false;
  }
  return positionIn(subterm,term->term(),position);
}

bool Kernel::positionIn(TermList& subterm,Term* term,std::string& position)
{
  //cout << "positionIn " << subterm.toString() << " in " << term->toString() << endl;

  if(subterm.isTerm() && subterm.term()==term){
    position = "1";
    return true;
  }
  if(term->arity()==0) return false;

  unsigned pos=1;
  TermList* ts = term->args();
  while(true){
    if(*ts==subterm){
      position=Lib::Int::toString(pos);
      return true;
    }
    if(positionIn(subterm,ts,position)){
      position = Lib::Int::toString(pos) + "." + position;
      return true;
    }
    pos++;
    ts = ts->next();
    if(ts->isEmpty()) break;
  }

  return false;
}

TermList Term::termArg(unsigned n) const
{
  ASS_LE(0, n)
  ASS_L(n, numTermArguments())
  return *nthArgument(n + numTypeArguments());
}

TermList Term::typeArg(unsigned n) const
{
  ASS_LE(0, n)
  ASS_L(n, numTypeArguments())
  return *nthArgument(n);
}

bool Term::computable() const {
  if (!env.signature->getFunction(this->functor())->computable()) {
    return false;
  }
  SubtermIterator sit(this);
  while (sit.hasNext()) {
    TermList t = sit.next();
    if (!t.isTerm() || !env.signature->getFunction(t.term()->functor())->computable()) {
      return false;
    }
  }
  return true;
}

bool Term::computableOrVar() const {
  if (!env.signature->getFunction(this->functor())->computable()) {
    return false;
  }
  SubtermIterator sit(this);
  while (sit.hasNext()) {
    TermList t = sit.next();
    if (t.isTerm() && !env.signature->getFunction(t.term()->functor())->computable()) {
      return false;
    }
  }
  return true;
}

std::ostream& Kernel::operator<<(std::ostream& out, SpecialFunctor const& self)
{
  switch (self) {
    case SpecialFunctor::ITE: return out << "ITE";
    case SpecialFunctor::LET: return out << "LET";
    case SpecialFunctor::FORMULA: return out << "FORMULA";
    case SpecialFunctor::TUPLE: return out << "TUPLE";
    case SpecialFunctor::LET_TUPLE: return out << "LET_TUPLE";
    case SpecialFunctor::LAMBDA: return out << "LAMBDA";
    case SpecialFunctor::MATCH: return out << "SPECIAL_FUNCTOR_LAST ";
  }
  ASSERTION_VIOLATION
}
