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
 * @file RewritingByRule.cpp
 * Implements class RewritingByRule.
 */

#include "RewritingByRule.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Saturation/SaturationAlgorithm.hpp"
#include "Indexing/TermSubstitutionTree.hpp"

using namespace Kernel;
using namespace Inferences;

ClauseIterator SuperpositionByRule::generateClauses(Clause* c)
{
  CALL("SuperpositionByRule::generateClauses");

  auto cLen = c->length();
  auto& ord = _salg->getOrdering();
  const auto& rules = c->getRewriteRules();
  auto stats = env.statistics;

  if (_tis) {
    delete _tis;
  }
  _tis = new TermSubstitutionTree();

  iterTraits(rules.items())
    .forEach([this](pair<TermList,TermList> kv) {
      _tis->insert(kv.first, kv.second);
    });

  return pvi(iterTraits(c->getSelectedLiteralIterator())
    .flatMap([&ord](Literal* lit) {
      TermIterator it = env.options->combinatorySup() ? EqHelper::getFoSubtermIterator(lit, ord) :
                                                        EqHelper::getSubtermIterator(lit, ord);
      return pvi( pushPairIntoRightIterator(lit, it) );
    })
    .flatMap([this](pair<Literal*, TermList> kv) {
      // TODO I think this unification is giving too many inferences,
      // what we need would be taking both the query and the result with
      // variables from the same bank, since we are rewriting by rules
      // that contain variables from this clause. 
      return pvi(pushPairIntoRightIterator(kv, _tis->getUnifications(kv.second, true)));
    })
    .map([c,&rules,&ord,stats](pair<pair<Literal*, TermList>,TermQueryResult> arg) -> Clause* {
      auto eqLHS = arg.second.term;
      auto subst = arg.second.substitution;
      auto tgtTerm = rules.get(eqLHS);

      TermList eqLHSS = subst->apply(eqLHS, 1);
      TermList tgtTermS = subst->apply(tgtTerm, 1);
      Literal* rwLitS = subst->apply(arg.first.first, 0);

      //check that we're not rewriting smaller subterm with larger
      if (Ordering::isGorGEorE(ord.compare(tgtTermS,eqLHSS))) {
        return nullptr;
      }

      if (rwLitS->isEquality()) {
        // check that we're not rewriting only the smaller side of an equality
        TermList arg0 = rwLitS->termArg(0);
        TermList arg1 = rwLitS->termArg(1);

        if(!arg0.containsSubterm(eqLHSS)) {
          if(Ordering::isGorGEorE(ord.getEqualityArgumentOrder(rwLitS))) {
            return nullptr;
          }
        } else if(!arg1.containsSubterm(eqLHSS)) {
          if(Ordering::isGorGEorE(Ordering::reverse(ord.getEqualityArgumentOrder(rwLitS)))) {
            return nullptr;
          }
        }
      }

      Clause* res = new(c->length()) Clause(c->length(), GeneratingInference1(InferenceRule::SUPERPOSITION_BY_RULE, c));
      for (unsigned i = 0;i < c->length(); i++) {
        Literal* currAfter = EqHelper::replace(subst->apply((*c)[i], 0),eqLHSS,tgtTermS);
        if(EqHelper::isEqTautology(currAfter)) {
          res->destroy();
          return nullptr;
        }
        (*res)[i] = currAfter;
      }

      auto rwIt = c->getRewriteRules().items();
      while (rwIt.hasNext()) {
        auto kv = rwIt.next();
        res->addRewriteRule(
          subst->apply(kv.first, 0),
          subst->apply(kv.second, 0)
        );
      }
      auto rwBIt = c->getBlockedTerms();
      while (rwBIt.hasNext()) {
        res->addBlockedTerm(subst->apply(rwBIt.next(), 0));
      }
      // auto rwstit = RewriteableSubtermsFn(ordering)(rwLit);
      // while (rwstit.hasNext()) {
      //   auto st = subst->apply(rwstit.next().second, eqIsResult);
      //   if (ordering.compare(rwTermS, st)==Ordering::Result::GREATER) {
      //     // cout << "adding blocked " << st << " for " << rwTermS << endl;
      //     res->addBlockedTerm(st);
      //   }
      // }
      res->addRewriteRule(eqLHSS,tgtTermS);

      stats->superpositionByRule++;
      return res;
    })
    .filter(NonzeroFn()));
}

Clause* DemodulationByRule::simplify(Clause* c)
{
  CALL("DemodulationByRule::simplify");

  auto cLen = c->length();
  auto& ord = _salg->getOrdering();
  for (unsigned i = 0; i < cLen; i++) {
    auto lit = (*c)[i];
    auto it = c->getRewriteRules().items();
    while (it.hasNext()) {
      auto kv = it.next();
      if (lit->containsSubterm(kv.first)) {
        if (ord.compare(kv.first,kv.second)!=Ordering::GREATER) {
          continue;
        }
        if (lit->isEquality() && (lit->termArg(0) == kv.first || lit->termArg(1) == kv.first)) {
          // TODO: perform demodulation redundancy check
          auto other = lit->termArg(0) == kv.first ? lit->termArg(1) : lit->termArg(0);
          Ordering::Result tord=ord.compare(kv.second, other);
          if (tord !=Ordering::LESS && tord!=Ordering::LESS_EQ) {
            bool isMax = true;
            for (unsigned j = 0; j < cLen; j++) {
              if (lit == (*c)[j]) {
                continue;
              }
              if (ord.compare(lit, (*c)[j]) == Ordering::LESS) {
                isMax=false;
                break;
              }
            }
            if(isMax) {
              // cout << "rule " << kv.first << "->" << kv.second << ", cl: " << *c << endl;
              continue;
            }
          }
        }

        Clause* res = new(cLen) Clause(cLen,
          SimplifyingInference1(InferenceRule::DEMODULATION_BY_RULE, c));
        for(unsigned j=0;j<cLen;j++) {
          (*res)[j] = EqHelper::replace((*c)[j],kv.first,kv.second);
        }
        auto rwIt = c->getRewriteRules().items();
        while (rwIt.hasNext()) {
          auto kv = rwIt.next();
          res->addRewriteRule(kv.first,kv.second);
        }
        auto rwBIt = c->getBlockedTerms();
        while (rwBIt.hasNext()) {
          res->addBlockedTerm(rwBIt.next());
        }
        env.statistics->demodulationByRule++;
        return res;
      }
    }
  }
  return c;
}
