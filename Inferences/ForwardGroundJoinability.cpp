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
 * @file ForwardGroundJoinability.cpp
 * Implements class ForwardGroundJoinability.
 */

#include "Lib/DHSet.hpp"
#include "Lib/Environment.hpp"
#include "Lib/VirtualIterator.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/EqHelper.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Ordering.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/TermIterators.hpp"
#include "Kernel/ColorHelper.hpp"

#include "Indexing/IndexManager.hpp"

#include "Saturation/SaturationAlgorithm.hpp"

#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"

#include "DemodulationHelper.hpp"

#include "ForwardGroundJoinability.hpp"

namespace Inferences {

using namespace Lib;
using namespace Kernel;
using namespace Indexing;
using namespace Saturation;
using namespace std;

namespace {

struct Applicator : SubstApplicator {
  Applicator(ResultSubstitution* subst) : subst(subst) {}
  TermList operator()(unsigned v) const override {
    return subst->applyToBoundResult(v);
  }
  ResultSubstitution* subst;
};

} // end namespace

void ForwardGroundJoinability::attach(SaturationAlgorithm* salg)
{
  ForwardSimplificationEngine::attach(salg);
  _index=static_cast<DemodulationLHSIndex*>(
	  _salg->getIndexManager()->request(DEMODULATION_LHS_CODE_TREE) );
}

void ForwardGroundJoinability::detach()
{
  _index=0;
  _salg->getIndexManager()->release(DEMODULATION_LHS_CODE_TREE);
  ForwardSimplificationEngine::detach();
}

using Position = Stack<unsigned>;

bool toTheLeftStrict(const Position& p1, const Position& p2, bool& prefix)
{
  prefix = false;
  for (unsigned i = 0; i < p1.size(); i++) {
    if (p2.size() <= i) {
      return false;
    }
    if (p1[i] != p2[i]) {
      return p1[i] < p2[i];
    }
  }
  prefix = true;
  return false;
}

bool isUnderVariablePosition(const Position& p, TermList lhs)
{
  if (lhs.isVar()) {
    return true;
  }
  auto curr = lhs.term();
  for (unsigned i = 0; i < p.size(); i++) {
    ASS_L(p[i],curr->arity());
    auto next = *curr->nthArgument(i);
    if (next.isVar()) {
      return true;
    }
    curr = next.term();
  }
  return false;
}

string posToString(const Position& pos)
{
  string res;
  for (const auto& i : pos) {
    res += "." + Int::toString(i);
  }
  return res;
}

TermList replace(TermList t, TermList from, TermList to)
{
  if (t == from) {
    return to;
  }
  if (t.isVar()) {
    return t;
  }
  return TermList(EqHelper::replace(t.term(), from, to));
}

bool ForwardGroundJoinability::perform(Clause* cl, Clause*& replacement, ClauseIterator& premises)
{
  Ordering& ordering = _salg->getOrdering();

  static DHSet<TermList> attempted;

  if (cl->length()>1) {
    return false;
  }

  auto clit = (*cl)[0];
  if (!clit->isEquality() || clit->isNegative()) {
    return false;
  }
  DHSet<Clause*> premiseSet;
  struct State {
    TermList left;
    TermList right;
    bool L;
    bool R;
  };
  Stack<State*> states;
  auto curr = new State{
    clit->termArg(0), clit->termArg(1), true, true
  };
  states.push(curr);
  RedundancyCheck checker(ordering, curr);
  auto tpo = TermPartialOrdering::getEmpty(ordering);
  unsigned cnt = 0;

  auto cleanup = [&]() {
    while (states.isNonEmpty()) {
      delete states.pop();
    }
  };

  while (curr) {
    if (cnt++ > 500) {
      cleanup();
      return false;
    }
    attempted.reset();

    auto kv = checker.next({ { curr->left, curr->right, Ordering::EQUAL } }, nullptr);
    curr = static_cast<State*>(kv.first);
    tpo = kv.second;
    if (!curr) {
      break;
    }

    NonVariableNonTypeIterator it(curr->left, curr->right);
    while (it.hasNext()) {
      TypedTermList trm = it.next();
      if (!attempted.insert(trm)) {
        it.right();
        continue;
      }

      bool redundancyCheck =
        (curr->L && trm == curr->left) || 
        (curr->R && trm == curr->right);

      auto git = _index->getGeneralizations(trm.term(), /* retrieveSubstitutions */ true);
      while(git.hasNext()) {
        auto qr=git.next();
        ASS_EQ(qr.data->clause->length(),1);

        auto lhs = qr.data->term;
        if (lhs.isVar()) {
          // we are not interested in these for now
          continue;
        }

        TermList rhs = qr.data->rhs;

        auto subs = qr.unifier;
        ASS(subs->isIdentityOnQueryWhenResultBound());

        Applicator appl(subs.ptr());

        OrderingConstraints cons;

        TermList rhsS;

        Result comp = ordering.compare(AppliedTerm(lhs, &appl, true), AppliedTerm(rhs, &appl, true), tpo);
        if (comp == Ordering::LESS || comp == Ordering::EQUAL) {
          continue;
        } else if (comp == Ordering::INCOMPARABLE) {
          rhsS = subs->applyToBoundResult(rhs);
          cons.push({ trm, rhsS, Ordering::GREATER });
        }

        // encompassing demodulation is fine when rewriting the smaller guy
        if (redundancyCheck) {
          // this will only run at most once;
          // could have been factored out of the getGeneralizations loop,
          // but then it would run exactly once there
          Ordering::Result litOrder = ordering.compare(curr->left,curr->right);
          if ((trm==curr->left && litOrder == Ordering::LESS) ||
              (trm==curr->right && litOrder == Ordering::GREATER)) {
            redundancyCheck = false;
          }
        }

        if (redundancyCheck && DemodulationHelper::isRenamingOn(&appl,lhs)) {
          TermList other = trm == curr->left ? curr->right : curr->left;
          auto redComp = ordering.compare(AppliedTerm(rhs, &appl, true), AppliedTerm(other));
          if (redComp == Ordering::LESS || redComp == Ordering::EQUAL) {
            continue;
          } else if (redComp == Ordering::INCOMPARABLE) {
            if (rhsS.isEmpty()) {
              rhsS = subs->applyToBoundResult(rhs);
            }
            cons.push({ other, rhsS, Ordering::GREATER });
          }
        }

        if (rhsS.isEmpty()) {
          rhsS = subs->applyToBoundResult(rhs);
        }

        auto newLeft = replace(curr->left, trm, rhsS);
        auto newRight = replace(curr->right, trm, rhsS);

        std::pair<void*,const TermPartialOrdering*> next;
        if (newLeft == newRight) {
          next = checker.next(cons, nullptr);
        } else {
          bool L = curr->L && curr->left == newLeft;
          bool R = curr->R && curr->right == newRight;
          auto newState = new State{ newLeft, newRight, L, R };
          states.push(newState);
          next = checker.next(cons, newState);
        }
        // TODO can it be that we go to a different node in the tree and get the same result still?
        if (next.first != curr || next.second != tpo) {
          curr = static_cast<State*>(next.first);
          tpo = next.second;
          premiseSet.insert(qr.data->clause);
          goto LOOP_END;
        }
      }
    }
    cleanup();
    return false;
LOOP_END:
    continue;
  }
  premises = pvi(getPersistentIterator(premiseSet.iterator()));

  env.statistics->groundRedundantClauses++;
  cleanup();
  return true;
}

TermList ForwardGroundJoinability::normalize(TermList t, const TermPartialOrdering* tpo) const
{
  return t;
}

ForwardGroundJoinability::RedundancyCheck::RedundancyCheck(const Ordering& ord, void* data)
  : comp(ord.createComparator(/*onlyVars=*/true))
{
  comp->_source = Branch(data, comp->_sink);
  comp->_source.node()->ready = true;
  path.push(&comp->_source);
}

std::pair<void*,const TermPartialOrdering*> ForwardGroundJoinability::RedundancyCheck::next(OrderingConstraints ordCons, void* data)
{
  static Ordering::Result ordVals[] = { Ordering::EQUAL, Ordering::GREATER, Ordering::INCOMPARABLE };
  ASS(path.isNonEmpty());

  auto curr = path.top();
  ASS_EQ(curr->node()->tag, BranchTag::T_DATA);
  ASS(curr->node()->data);
  ASS(curr->node()->ready); 
  ASS_EQ(curr->node()->refcnt,1);

  // current node has to be processed again
  curr->node()->ready = false;

  // We replace (not modify) the current node
  // with a new subtree containing ordCons and data
  // and pointing to the original node otherwise.

  Branch origB(*curr);
  Branch newB = data ? Branch(data, comp->_sink) : comp->_sink;

  for (const auto& [lhs,rhs,rel] : ordCons) {
    *curr = Branch(lhs, rhs);
    for (unsigned i = 0; i < 3; i++) {
      if (ordVals[i] != rel) {
        curr->node()->getBranch(ordVals[i]) = origB;
      }
    }
    curr = &curr->node()->getBranch(rel);
  }
  *curr = newB;

  while (path.isNonEmpty()) {
    if (path.size()==1) {
      comp->_prev = nullptr;
    } else {
      comp->_prev = path[path.size()-2];
    }
    comp->_curr = path.top();
    comp->expand();

    auto node = comp->_curr->node();
    if (node->tag == BranchTag::T_DATA && !node->data) {
      pushNext();
      continue;
    }

    if (!node->trace || node->trace->hasIncomp()) {
      // invalid branch, continue
      pushNext();
      continue;
    }

    if (node->tag == BranchTag::T_DATA) {
      ASS(node->data);
      return { (Literal*)node->data, node->trace };
    }
    path.push(&node->gtBranch);
  }

  ASS(path.isEmpty());
  return { nullptr, nullptr };
}

void ForwardGroundJoinability::RedundancyCheck::pushNext()
{
  while (path.isNonEmpty()) {
    auto curr = path.pop();
    if (path.isEmpty()) {
      continue;
    }

    auto prev = path.top()->node();
    ASS(prev->tag == BranchTag::T_POLY || prev->tag == BranchTag::T_TERM);
    // if there is a previous node and we were either in the gt or eq
    // branches, just go to next branch in order, otherwise backtrack
    if (curr == &prev->gtBranch) {
      path.push(&prev->eqBranch);
      break;
    }
    if (curr == &prev->eqBranch) {
      path.push(&prev->ngeBranch);
      break;
    }
  }
}

}