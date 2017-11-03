/**
 * @file KBOForEPR.cpp
 * Implements class KBOForEPR for instances of the Knuth-Bendix ordering for EPR problems
 */

#include "Debug/Tracer.hpp"


#include "Lib/Environment.hpp"
#include "Lib/Comparison.hpp"
#include "Lib/Int.hpp"

#include "Indexing/TermSharing.hpp"

#include "Shell/Property.hpp"

#include "Problem.hpp"
#include "Term.hpp"
#include "Signature.hpp"

#include "KBOForEPR.hpp"

namespace Kernel {

using namespace Lib;


/**
 * Create a KBO object.
 *
 * Function and predicate preferences and predicate levels
 * must be initialized after calling this constructor and
 * before any comparisons using this object are being made.
 */
KBOForEPR::KBOForEPR(Problem& prb, const Options& opt)
 : PrecedenceOrdering(prb, opt)
{
  CALL("KBOForEPR::KBOForEPR");
  ASS_EQ(prb.getProperty()->maxFunArity(),0);
}

/**
 * Compare arguments of literals l1 and l2 and return the result
 * of the comparison.
 * @since 07/05/2008 flight Manchester-Brussels
 */
Ordering::Result KBOForEPR::compare(Literal* l1, Literal* l2) const
{
  CALL("KBOForEPR::compare(Literal*...)");
  ASS(l1->shared());
  ASS(l2->shared());

  if (l1 == l2) {
    return EQUAL;
  }

  unsigned p1 = l1->functor();
  unsigned p2 = l2->functor();

  if( (l1->isNegative() ^ l2->isNegative()) && (p1==p2) &&
	  l1->weight()==l2->weight() && l1->vars()==l2->vars() &&
	  l1==env.sharing->tryGetOpposite(l2)) {
    return l1->isNegative() ? LESS : GREATER;
  }

  Result res;
  if (p1 != p2) {
    Comparison levComp=Int::compare(predicateLevel(p1),predicateLevel(p2));
    if(levComp!=Lib::EQUAL) {
      res = fromComparison(levComp);
      goto fin;
    }
  }

  if(l1->isEquality()) {
    ASS(l2->isEquality());
    return compareEqualities(l1, l2);
  }
  ASS(!l1->isEquality());
  ASS(!l2->isEquality());

  if (p1 != p2) {
    Comparison arComp=Int::compare(l1->arity(),l2->arity());
    if(arComp!=Lib::EQUAL) {
      //since on the ground level each literal argument must be a constant,
      //and all symbols are of weight 1, the literal with higher arity is
      //heavier and therefore greater
      res = fromComparison(arComp);
      goto fin;
    }

    Comparison prComp=Int::compare(predicatePrecedence(p1),predicatePrecedence(p2));
    ASS_NEQ(prComp, Lib::EQUAL); //precedence ordering is total
    res = fromComparison(prComp);
    goto fin;
  }

  {
    TermList* t1=l1->args();
    TermList* t2=l2->args();

    ASS(!t1->isEmpty());
    while(*t1==*t2) {
      t1=t1->next();
      t2=t2->next();
      ASS(!t1->isEmpty()); //if we're at the end of the term, the literals would have been the same
    }
    res = compare(*t1, *t2);
    goto fin;
  }

fin:
  if(_reverseLCM && (l1->isNegative() || l2->isNegative()) ) {
    if(l1->isNegative() && l2->isNegative()) {
      res = reverse(res);
    }
    else {
      res = l1->isNegative() ? LESS : GREATER;
    }
  }
  return res;
}

Ordering::Result KBOForEPR::compare(TermList tl1, TermList tl2) const
{
  CALL("KBOForEPR::compare(TermList)");
  ASS(!tl1.isTerm() || tl1.term()->arity()==0)
  ASS(!tl2.isTerm() || tl2.term()->arity()==0)

  if(tl1==tl2) {
    return EQUAL;
  }

  if(tl1.isOrdinaryVar() || tl2.isOrdinaryVar()) {
    return INCOMPARABLE;
  }

  //We're dealing with constants here -- if the top symbols were the same,
  //the terms would be equal as well.
  ASS_NEQ(tl1.term()->functor(), tl2.term()->functor());

  return compareFunctionPrecedences(tl1.term()->functor(), tl2.term()->functor());
}

}
