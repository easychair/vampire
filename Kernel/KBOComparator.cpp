/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#include "KBO.hpp"
#include "SubstHelper.hpp"
#include "Term.hpp"
#include "TermIterators.hpp"

#include "KBOComparator.hpp"

namespace Kernel {

using namespace std;
using namespace Lib;
using namespace Shell;

void KBOComparator::expandTermCase()
{
  const auto& kbo = static_cast<const KBO&>(_ord);

  // weight and variable balances first

  // we only care about the non-zero weights and counts
  bool varInbalance = false;
  auto state = kbo._state;
#if VDEBUG
  // we make sure kbo._state is not used while we're using it
  kbo._state = nullptr;
#endif
  auto w = state->_weightDiff;
  decltype(state->_varDiffs)::Iterator vit(state->_varDiffs);
  Stack<VarCoeffPair> nonzeros;
  while (vit.hasNext()) {
    unsigned v;
    int cnt;
    vit.next(v,cnt);
    if (cnt!=0) {
      nonzeros.push({ v, cnt });
      w-=cnt; // we have to remove the variable weights from w
    }
    if (cnt<0) {
      varInbalance = true;
    }
  }
#if VDEBUG
  kbo._state = state;
  state = nullptr;
#endif

  auto node = _curr->node();
  ASS(node->lhs.isTerm() && node->rhs.isTerm());
  auto lhs = node->lhs.term();
  auto rhs = node->rhs.term();

  auto eqBranch = node->eqBranch;
  auto gtBranch = node->gtBranch;
  auto ngeBranch = node->ngeBranch;

  auto curr = _curr;
  bool weightAdded = (w < 0 || varInbalance);
  if (weightAdded) {
    // we mutate the original node
    curr->node()->tag = T_POLY;
    curr->node()->poly = Polynomial::get(w, nonzeros);
    curr->node()->gtBranch = gtBranch;
    curr->node()->ngeBranch = ngeBranch;
    curr = &curr->node()->eqBranch;
  }

  Ordering::Result prec = lhs->isSort()
    ? kbo.compareTypeConPrecedences(lhs->functor(),rhs->functor())
    : kbo.compareFunctionPrecedences(lhs->functor(),rhs->functor());
  switch (prec)
  {
    case Ordering::LESS: {
      *curr = ngeBranch;
      break;
    }
    case Ordering::GREATER: {
      *curr = gtBranch;
      break;
    }
    case Ordering::EQUAL: {
      for (unsigned i = 0; i < lhs->arity(); i++) {
        auto lhsArg = *lhs->nthArgument(i);
        auto rhsArg = *rhs->nthArgument(i);
        // we mutate the original node in the first iteration
        if (!weightAdded && i==0) {
          curr->node()->lhs = lhsArg;
          curr->node()->rhs = rhsArg;
        } else {
          *curr = Branch(lhsArg,rhsArg);
          curr->node()->gtBranch = gtBranch;
          curr->node()->ngeBranch = ngeBranch;
        }
        curr = &curr->node()->eqBranch;
      }
      *curr = eqBranch;
      break;
    }
    default: {
      ASSERTION_VIOLATION;
    }
  }
}

}
