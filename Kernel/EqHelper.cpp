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
 * @file EqHelper.cpp
 * Implements class EqHelper.
 */

#include "Lib/Environment.hpp"

#include "Shell/Options.hpp"

#include "Ordering.hpp"
#include "SortHelper.hpp"
#include "TermIterators.hpp"
#include "ApplicativeHelper.hpp"
#include "Lib/Metaiterators.hpp"
#include "Matcher.hpp"

#include "EqHelper.hpp"

namespace Kernel {

using namespace Shell;

/*
 * Turns an iterator with TermList elemenst to an iterator of TypedTermList elements
 * with the sort of the provided equality literal.
 *
 * Computes the sort exactly once, so it's wastful to call if the the iterator is empty.
 * (Or if the sort is already known.)
 */
template<class TermListIter>
auto withEqualitySort(Literal* eq, TermListIter iter)
{
  ASS(iter.hasNext());
  auto sort = SortHelper::getEqualityArgumentSort(eq);
  return pvi(iterTraits(std::move(iter))
    .map([sort](TermList t) { return TypedTermList(t, sort); }));
}

/**
 * Return the other side of an equality @b eq than the @b lhs
 */
TermList EqHelper::getOtherEqualitySide(Literal* eq, TermList lhs)
{
  ASS(eq->isEquality());

  if (*eq->nthArgument(0) == lhs) {
    return *eq->nthArgument(1);
  }
  ASS(*eq->nthArgument(1) == lhs);
  return *eq->nthArgument(0);
} // getOtherEqualitySide

bool EqHelper::hasGreaterEqualitySide(Literal* eq, const Ordering& ord, TermList& lhs, TermList& rhs)
{
  ASS(eq->isEquality());

  switch(ord.getEqualityArgumentOrder(eq)) {
    case Ordering::INCOMPARABLE:
      return false;
    case Ordering::GREATER:
      lhs = *eq->nthArgument(0);
      rhs = *eq->nthArgument(1);
      return true;
    case Ordering::LESS:
      lhs = *eq->nthArgument(1);
      rhs = *eq->nthArgument(0);
      return true;
    //there should be no equality literals of equal terms
    case Ordering::EQUAL:
      ASSERTION_VIOLATION;
  }
  ASSERTION_VIOLATION;
}

Literal* EqHelper::replace(Literal* lit, TermList what, TermList by)
{
  return static_cast<Literal*>(replace(static_cast<Term*>(lit), what, by));
}

/**
 * Replace all occurences of the subterm  @b tSrc by @b tDest in the term/literal
 * @b lit, and return the result
 *
 * Cannot be used to replace a sort
 */
Term* EqHelper::replace(Term* trm0, TermList tSrc, TermList tDest)
{
  ASS(trm0->shared());
  ASS(!trm0->isSort());
  ASS(tSrc.isVar() || !tSrc.term()->isSort());
  ASS(tDest.isVar() || !tDest.term()->isSort());


  static Stack<TermList*> toDo(8);
  static Stack<Term*> terms(8);
  static Stack<bool> modified(8);
  static Stack<TermList> args(8);
  ASS(toDo.isEmpty());
  ASS(terms.isEmpty());
  modified.reset();
  args.reset();

  modified.push(false);
  toDo.push(trm0->args());

  for (;;) {
    TermList* tt=toDo.pop();
    if (tt->isEmpty()) {
      if (terms.isEmpty()) {
	//we're done, args stack contains modified arguments
	//of the literal.
	ASS(toDo.isEmpty());
	break;
      }
      Term* orig=terms.pop();
      if (!modified.pop()) {
	args.truncate(args.length() - orig->arity());
	args.push(TermList(orig));
	continue;
      }
      //here we assume, that stack is an array with
      //second topmost element as &top()-1, third at
      //&top()-2, etc...
      TermList* argLst=&args.top() - (orig->arity()-1);
      args.truncate(args.length() - orig->arity());

      args.push(TermList(Term::create(orig,argLst)));
      modified.setTop(true);
      continue;
    }
    toDo.push(tt->next());

    TermList tl=*tt;
    if (tl == tSrc) {
      args.push(tDest);
      modified.setTop(true);
      continue;
    }
    if (tl.isVar() || tl.term()->isSort()) {
      args.push(tl);
      continue;
    }
    ASS(tl.isTerm());
    Term* t=tl.term();
    terms.push(t);
    modified.push(false);
    toDo.push(t->args());
  }
  ASS(toDo.isEmpty());
  ASS(terms.isEmpty());
  ASS_EQ(modified.length(),1);
  ASS_EQ(args.length(),trm0->arity());

  if (!modified.pop()) {
    // we call replace in superposition only if we already know,
    // there is something to be replaced.
    // ASSERTION_VIOLATION; // MS: but there is now a new use in InnerRewriting which does not like this extra check
    return trm0;
  }

  // here we assume, that stack is an array with
  // second topmost element as &top()-1, third at
  // &top()-2, etc...
  TermList* argLst=&args.top() - (trm0->arity()-1);
  if (trm0->isLiteral()) {
    Literal* lit = static_cast<Literal*>(trm0);
    ASS_EQ(args.size(), lit->arity());
    return Literal::create(lit,argLst);
  }
  return Term::create(trm0,argLst);
}


VirtualIterator<Term*> EqHelper::getSubtermIterator(Literal* lit, const Ordering& ord)
{
  return getRewritableSubtermIterator<NonVariableNonTypeIterator>(lit, ord);
}

TermIterator EqHelper::getBooleanSubtermIterator(Literal* lit, const Ordering& ord)
{
  return getRewritableSubtermIterator<BooleanSubtermIt>(lit, ord);
}

VirtualIterator<Term*> EqHelper::getFoSubtermIterator(Literal* lit, const Ordering& ord)
{
  return getRewritableSubtermIterator<FirstOrderSubtermIt>(lit, ord);
}

TermIterator EqHelper::getNarrowableSubtermIterator(Literal* lit, const Ordering& ord)
{
  return getRewritableSubtermIterator<NarrowableSubtermIt>(lit, ord);
} 

/*
 * Function is used in the higher-order inference SubVarSup
 */
VirtualIterator<TypedTermList> EqHelper::getRewritableVarsIterator(DHSet<unsigned>* unstableVars, Literal* lit, const Ordering& ord)
{
  ASS(lit->isEquality());
    
  TermList sel;
  switch(ord.getEqualityArgumentOrder(lit)) {
  case Ordering::INCOMPARABLE: {
    return pvi(iterTraits(RewritableVarsIt(unstableVars, lit))
      .unique()
      .persistent());
  }
  case Ordering::EQUAL:
  case Ordering::GREATER:
    sel=*lit->nthArgument(0);
    break;
  case Ordering::LESS:
    sel=*lit->nthArgument(1);
    break;
#if VDEBUG
  default:
    ASSERTION_VIOLATION;
#endif
  }
  if (!sel.isTerm()) {
    return VirtualIterator<TypedTermList>::getEmpty();
  }
  return pvi(iterTraits(RewritableVarsIt(unstableVars, sel.term(), true)).unique().persistent());
} 


/**
 * Return iterator on subterms of a literal, that can be rewritten by
 * superposition.
 */
template<class SubtermIterator>
VirtualIterator<ELEMENT_TYPE(SubtermIterator)> EqHelper::getRewritableSubtermIterator(Literal* lit, const Ordering& ord)
{
  if (lit->isEquality()) {
    TermList sel;
    switch(ord.getEqualityArgumentOrder(lit)) {
    case Ordering::INCOMPARABLE: {
      SubtermIterator si(lit);
      return getUniquePersistentIteratorFromPtr(&si);
    }
    case Ordering::EQUAL:
    case Ordering::GREATER:
      sel=*lit->nthArgument(0);
      break;
    case Ordering::LESS:
      sel=*lit->nthArgument(1);
      break;
#if VDEBUG
    default:
      ASSERTION_VIOLATION;
#endif
    }
    if (!sel.isTerm()) {
      return VirtualIterator<ELEMENT_TYPE(SubtermIterator)>::getEmpty();
    }
    return getUniquePersistentIterator(vi(new SubtermIterator(sel.term(), true)));
  }

  SubtermIterator si(lit);
  return getUniquePersistentIteratorFromPtr(&si);

}


/**
 * Return iterator on sides of the equality @b lit that can be used as an LHS
 * for a rewriting inference (i.e. the other side of the equality is not greater)
 *
 * If the literal @b lit is not a positive equality, empty iterator is returned.
 */
VirtualIterator<TypedTermList> EqHelper::getLHSIterator(Literal* lit, const Ordering& ord)
{
  if (lit->isEquality()) {
    if (lit->isNegative()) {
      return VirtualIterator<TypedTermList>::getEmpty();
    }
    TermList t0=*lit->nthArgument(0);
    TermList t1=*lit->nthArgument(1);
    switch(ord.getEqualityArgumentOrder(lit))
    {
    case Ordering::INCOMPARABLE:
      return withEqualitySort(lit, iterItems(t0, t1) );
    case Ordering::GREATER:
      return withEqualitySort(lit, iterItems(t0) );
    case Ordering::LESS:
      return withEqualitySort(lit, getSingletonIterator(t1) );
    //there should be no equality literals of equal terms
    case Ordering::EQUAL:
      ASSERTION_VIOLATION;
    }
    return VirtualIterator<TypedTermList>::getEmpty();
  } else {
    return VirtualIterator<TypedTermList>::getEmpty();
  }
}

/**
 * A functor that returns true iff its argument is a non-variable term
 */
struct EqHelper::IsNonVariable
{
  bool operator()(TermList t)
  { return t.isTerm(); }
};

/**
 * Return iterator on sides of the equality @b lit that can be used as an LHS
 * for superposition
 *
 * If the literal @b lit is not a positive equality, empty iterator is returned.
 */
VirtualIterator<TypedTermList> EqHelper::getSuperpositionLHSIterator(Literal* lit, const Ordering& ord, const Options& opt)
{
  if (opt.superpositionFromVariables()) {
    return getLHSIterator(lit, ord);
  }
  else {
    return pvi( getFilteredIterator(getLHSIterator(lit, ord), IsNonVariable()) );
  }
}



VirtualIterator<TypedTermList> EqHelper::getSubVarSupLHSIterator(Literal* lit, const Ordering& ord)
{
  ASS(lit->isEquality());

  TermList eqSort = SortHelper::getEqualityArgumentSort(lit);

  if (eqSort.isVar() || eqSort.isArrowSort()) {
    if (lit->isNegative()) {
      return VirtualIterator<TypedTermList>::getEmpty();
    }

    TermList t0=*lit->nthArgument(0);
    TermList t1=*lit->nthArgument(1);
    TermList t0Head = ApplicativeHelper::getHead(t0);
    TermList t1Head = ApplicativeHelper::getHead(t1);
    bool t0hisVarOrComb = ApplicativeHelper::isComb(t0Head) || t0Head.isVar();
    bool t1hisVarOrComb = ApplicativeHelper::isComb(t1Head) || t1Head.isVar();

    switch(ord.getEqualityArgumentOrder(lit))
    {
    case Ordering::INCOMPARABLE:
      if(t0hisVarOrComb && t1hisVarOrComb){
        return pvi(iterItems(TypedTermList(t0,eqSort), TypedTermList(t1,eqSort)));
      } else if( t0hisVarOrComb ){
        return pvi(getSingletonIterator(TypedTermList(t1,eqSort)));
      } else if( t1hisVarOrComb ) {
        return pvi(getSingletonIterator(TypedTermList(t0,eqSort)));
      }
      break;
    case Ordering::GREATER:
      if(t1hisVarOrComb){
        return pvi(getSingletonIterator(TypedTermList(t0,eqSort)));
      }
      break;
    case Ordering::LESS:
      if(t0hisVarOrComb){
        return pvi(getSingletonIterator(TypedTermList(t1,eqSort)));
      }
      break;
    case Ordering::EQUAL:
      //there should be no equality literals of equal terms
    default:
      ASSERTION_VIOLATION;
    }
    return VirtualIterator<TypedTermList>::getEmpty();
  } else {
    return VirtualIterator<TypedTermList>::getEmpty();
  }
}

/**
 * Return pair of (i) iterator on sides of the equality @b lit that can be used
 * as an LHS for demodulation and (ii) a boolean representing whether the
 * demodulator(s) are preordered. Note that (ii) can only be true if (i)
 * contains only one element but it is not necessarily true in this case as
 * e.g. the set of variables on one side of an incomparable equation can
 * be a strict subset of the that of the other side.
 *
 * If the literal @b lit is not a positive equality, empty iterator and false are returned.
 */
std::pair<VirtualIterator<TypedTermList>,bool> EqHelper::getDemodulationLHSIterator(
  Literal* lit, bool onlyPreordered, const Ordering& ord)
{
  bool isPreordered = false;
  if (lit->isEquality()) {
    if (lit->isNegative()) {
      return { VirtualIterator<TypedTermList>::getEmpty(), isPreordered };
    }
    TermList t0=*lit->nthArgument(0);
    TermList t1=*lit->nthArgument(1);
    switch(ord.getEqualityArgumentOrder(lit))
    {
    case Ordering::INCOMPARABLE:
      if ( onlyPreordered ) {
        return { VirtualIterator<TypedTermList>::getEmpty(), isPreordered };
      }
      if (t0.containsAllVariablesOf(t1)) {
        if (t1.containsAllVariablesOf(t0)) {
          // If the equation is its own variant when oriented
          // reversed, there's no need to index both sides
          if (MatchingUtils::matchReversedArgs(lit, lit)) {
            return { withEqualitySort(lit, getSingletonIterator(t0) ), isPreordered };
          }
          return { withEqualitySort(lit, iterItems(t0, t1)), isPreordered };
        }
        return { withEqualitySort(lit, getSingletonIterator(t0) ), isPreordered };
      }
      if (t1.containsAllVariablesOf(t0)) {
        return { withEqualitySort(lit, getSingletonIterator(t1) ), isPreordered };
      }
      break;
    case Ordering::GREATER:
      ASS(t0.containsAllVariablesOf(t1));
      isPreordered = true;
      return { withEqualitySort(lit, getSingletonIterator(t0) ), isPreordered };
    case Ordering::LESS:
      ASS(t1.containsAllVariablesOf(t0));
      isPreordered = true;
      return { withEqualitySort(lit, getSingletonIterator(t1) ), isPreordered };
    //there should be no equality literals of equal terms
    case Ordering::EQUAL:
      ASSERTION_VIOLATION_REP(*lit);
    }
    return { VirtualIterator<TypedTermList>::getEmpty(), isPreordered };
  } else {
    return { VirtualIterator<TypedTermList>::getEmpty(), isPreordered };
  }
}

TermIterator EqHelper::getEqualityArgumentIterator(Literal* lit)
{
  ASS(lit->isEquality());

  return pvi( iterItems( *lit->nthArgument(0), *lit->nthArgument(1)) );
}


}
