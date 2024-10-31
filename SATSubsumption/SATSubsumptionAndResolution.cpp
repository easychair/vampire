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
  * @file SATSubsumptionAndResolution.cpp
  * Implements class SATSubsumptionAndResolution.
  */

  /**
   * THEORETICAL BACKGROUND
   *
   * The subsumption and subsumption resolution are described in the papers:
   * - 2022: "First-Order Subsumption via SAT Solving." by Jakob Rath, Armin Biere and Laura Kovács
   * - 2023: "SAT-Based Subsumption Resolution" by Robin Coutelier, Jakob Rath, Michael Rawson and
   *         Laura Kovács
   * - 2024: "SAT Solving for Variants of First-Order Subsumption" by Robin Coutelier, Jakob Rath,
   *         Michael Rawson, Armin Biere and Laura Kovács
   *
   * Note that in this implementation, we removed the direct encoding for subsumption resolution
   * because it does not bring a significant improvement.
   * The indirect encoding described in "SAT Solving for Variants of First-Order Subsumption" was
   * improved by ignoring the cⱼ variables for which there exists only b⁻ᵢⱼ variables. This reduces
   * the gap between the direct and indirect encoding on smaller instances, therefore reducing the
   * added value of the direct encoding.
   *
   *
   * TLDR on subsumption and subsumption resolution:
   * The SAT-based subsumption and subsumption resolution inference are defined as follows:
   *
   * ----- Subsumption: -----
   * Let L and M be two clauses considered as multisets. L subsumes M iff there exists a substitution σ
   * such that σ(L) ⊑ M, where ⊑ denotes the submultiset relation.
   * Subsumption can occurs iff the three following conditions are satisfied:
   * 1. Completeness : All literals of L have a substitution to M.
   *    There exists σ such that forall lᵢ in L, σ(lᵢ) in M
   *    ∀i ∃j. σ(lᵢ) = mⱼ
   * 2. Multiplicity conservation : For each literal lᵢ in L, there exists at most one
   *    literal mⱼ in M such that σ(lᵢ) = mⱼ.
   *    ∀iji'. (i ≠ i' ∧ σ(lᵢ) = mⱼ) ⇒ σ(lᵢ') ≠ mⱼ
   * 3. Substitution compatibility : The substitution σ is compatible with all the sub-substitutions.
   *
   * ----- Subsumption Resolution: -----
   * Let L and M be two clauses considered as sets. L and M are said to be the base and instance
   * of a subsumption resolution inference, respectively iif
   *    there exists a substitution σ,
   *                 a set of literal L' included in L
   *                 a literal m' in M
   *    such that σ(L') = {¬m'} and σ(L \ L') ⊆ M \ {¬m'}.
   * Subsumption resolution can occur iff the 5 following conditions are satisfied:
   * 1. Existence             : There exists a literal lᵢ in L and a literal mⱼ such that
   *    σ(lᵢ) = ¬mⱼ (m' exists).
   *    i.e. ∃ij. σ(lᵢ) = ¬mⱼ
   * 2. Uniqueness   : There is only one literal mⱼ such that there exists a literal lᵢ in L
   *    such that σ(lᵢ) = ¬mⱼ. (m' is unique)
   *    i.e. ∃j' ∀ij. (σ(lᵢ) = ¬mⱼ ⇒ j = j')
   * 3. Completeness          : All literals of L must either have a substitution to M - {m'} or {~m'}.
   *    Forall lᵢ in L, there exists mⱼ in M such that σ(lᵢ) = mⱼ or σ(lᵢ) = ¬mⱼ
   *    i.e. ∀i ∃j. (σ(lᵢ) = ¬mⱼ ∨ σ(lᵢ) = mⱼ)
   * 4. Coherence             : Literals in M cannot be mapped by both positive and negative substitutions.
   *    Forall mⱼ in M, forall lᵢ, lᵢ' != lᵢ in L, σ(lᵢ) = mⱼ ⇒ σ(lᵢ') != ¬mⱼ
   *    i.e. ∀j. (∃i σ(lᵢ) = mⱼ ⇒ ∀i σ(lᵢ) ≠ ¬mⱼ)
   * 5. Substitution compatibility : The substitution σ is compatible with all the sub-substitutions.
   */

#include "Kernel/Matcher.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Shell/Statistics.hpp"
#include "Debug/RuntimeStatistics.hpp"
#include <algorithm>
#include <iostream>
#include <type_traits>

#include "SATSubsumption/SATSubsumptionAndResolution.hpp"
#include "Shell/Statistics.hpp"
#include "SATSubsumptionAndResolution.hpp"

using namespace Indexing;
using namespace Kernel;
using namespace SATSubsumption;
using namespace std;

const unsigned INVALID = std::numeric_limits<unsigned>::max();

/// @brief If 1, prints the SAT clauses added to the solver on the standard output
#define PRINT_CLAUSES_SUBS 0
/// @brief If 1, prints some comments about the subsumption resolution process
#define PRINT_CLAUSE_COMMENTS_SUBS 0


/****************************************************************************/
/*               SATSubsumptionAndResolution::MatchSet                      */
/****************************************************************************/

void SATSubsumptionAndResolution::MatchSet::indexMatrix()
{
  if (_matchesByJ.size())
    return;

  ASS_EQ(_matchesByJ.size(), 0)
  ASS_EQ(_indexI.size(), 0)
  ASS_EQ(_indexJ.size(), 0)

  for (Match match : _matchesByI)
    _matchesByJ.push_back(match);

  std::sort(
    _matchesByJ.begin(),
    _matchesByJ.end(),
    [](Match left, Match right) { return left.j < right.j; });

  for (unsigned i = 0, idx = 0; i < _m; i++) {
    _indexI.push_back(idx);
    while (idx < _matchesByI.size() && _matchesByI[idx].i == i)
      idx++;
  }
  _indexI.push_back(_matchesByI.size());

  for (unsigned j = 0, idx = 0; j < _n; j++) {
    _indexJ.push_back(idx);
    while (idx < _matchesByJ.size() && _matchesByJ[idx].j == j)
      idx++;
  }
  _indexJ.push_back(_matchesByJ.size());
}

bool SATSubsumptionAndResolution::MatchSet::hasPositiveMatchJ(unsigned j)
{
  // the wizardry is explained in the header file
  ASS(j < _n)
  return (_jStates[j / 4] & (1 << (2 * (j % 4)))) != 0;
}

bool SATSubsumptionAndResolution::MatchSet::hasNegativeMatchJ(unsigned j)
{
  // the wizardry is explained in the header file
  ASS(j < _n)
  return (_jStates[j / 4] & (2 << (2 * (j % 4)))) != 0;
}

/****************************************************************************/
/*       SATSubsumptionAndResolution::SATSubsumptionAndResolution           */
/****************************************************************************/

void SATSubsumptionAndResolution::loadProblem(Clause* L,
                                              Clause* M)
{
  ASS(L)
  ASS(M)
#if VDEBUG
  // Check that two literals are not the same in L and M
  static DHSet<Literal*> lits;
  lits.reset();
  for (unsigned i = 0; i < L->length(); i++)
    if (!lits.insert((*L)[i]))
      ASS(false)
  lits.reset();
  for (unsigned i = 0; i < M->length(); i++)
    if (!lits.insert((*M)[i]))
      ASS(false)
#endif

#if PRINT_CLAUSES_SUBS
  cout << "----------------------------------------------" << endl;
  cout << "Setting up problem " << L->toString() << " " << M->toString() << endl;
#endif
  _L = L;
  _M = M;
  _m = L->length();
  _n = M->length();

  _matchSet.clear();
  _matchSet.resize(_m, _n);

  _subsumptionImpossible = false;
  _srImpossible = false;

  _solver.clear();
  _bindingsManager.clear();
} // SATSubsumptionAndResolution::loadProblem

/**
 * @brief Heuristically determine whether it is impossible to find a subsumption
 * between _L and _M.
 *
 * The idea is to check whether the multiset of predicates in _L is a subset of
 * the multiset of predicates in _M. If it is not, then it is impossible to find
 * a subsumption.
 *
 * @return true if subsumption is impossible, false if we don't know
 */
bool SATSubsumptionAndResolution::pruneSubsumption()
{
  ASS(_L)
  ASS(_M)

  if (_L->length() > _M->length()) {
    // if #(L) > #(M) than it is impossible that σ(L) is a sub-multiset of M
    _subsumptionImpossible = true;
    return true;
  }

  auto& headerMultiset = _pruneStorage;
  prune_t& timestamp = _pruneTimestamp;
  static_assert(
    std::is_same<std::remove_reference_t<decltype(timestamp)>,
    std::remove_reference_t<decltype(headerMultiset)>::value_type>::value,
    "timestamp and storage should be the same type");

  // multiset of signed predicates in M
  headerMultiset.resize(2 * env.signature->predicates(), 0);
  ASS(std::all_of(headerMultiset.begin(), headerMultiset.end(), [&](prune_t x) { return x <= timestamp; }))

  // Our relative zero for counting is the timestamp.
  // We need to reset the vector only if the counts could overflow.
  if (isAdditionOverflow<prune_t>(timestamp, _M->length())) {
    // when timestamp wraps around, we reset the vector to 0.
    std::fill(headerMultiset.begin(), headerMultiset.end(), 0);
    timestamp = 0;
  }

  prune_t const zero = timestamp;
  timestamp += _M->length();
  ASS(std::all_of(headerMultiset.begin(), headerMultiset.end(), [&](prune_t x) { return x <= zero; }))

  // fill in the multiset of functors in M
  for (unsigned i = 0; i < _M->length(); i++) {
    unsigned const hdr = (*_M)[i]->header();
    headerMultiset[hdr] = std::max(headerMultiset[hdr], zero) + 1;
  }

  // check if the multiset of functors in L is a subset of the multiset of functors in M
  for (unsigned j = 0; j < _L->length(); j++) {
    unsigned const hdr = (*_L)[j]->header();
    // we need to do the check before decrementing to avoid wraparound and keep the invariant valid
    if (headerMultiset[hdr] <= zero) {
      _subsumptionImpossible = true;
      return true;
    }
    headerMultiset[hdr]--;
  }

  // WARNING !!!
  // In the implementation, it is assumed that the pruning for subsumption resolution
  // is stronger than the pruning for subsumption. Therefore, if the pruning for subsumption
  // resolution occurs, then subsumption must be pruned as well.
  // The remark is left for whoever wants to change this property, the piece of code in
  // checkSubsumption() using the pruning methods must be changed as well.
  //    _srImpossible = pruneSubsumptionResolution();
  //    _subsumptionImpossible = _srImpossible || pruneSubsumption();
  ASS(!pruneSubsumptionResolution())
  return false;
} // SATSubsumptionAndResolution::pruneSubsumption

/**
 * @brief Heuristically determine whether it is impossible to apply subsumption
 * resolution between _L and _M.
 *
 * The idea is to check whether the set of predicates in _L is a subset of the
 * set of predicates in _M. If it is not, then it is impossible to find a
 * subsumption resolution.
 *
 * TODO: functorSet is initialized from _M which is the same for the whole forward loop.
 *       we could re-use it for multiple pruning checks instead of filling it everytime.
 *
 * @return true if subsumption resolution is impossible, false if we don't know
 */
bool SATSubsumptionAndResolution::pruneSubsumptionResolution()
{
  ASS(_L)
  ASS(_M)

  auto& functorSet = _pruneStorage;
  auto& timestamp = _pruneTimestamp;

  functorSet.resize(env.signature->predicates(), 0);
  ASS(std::all_of(functorSet.begin(), functorSet.end(), [&](prune_t x) { return x <= timestamp; }))

  timestamp++;
  if (timestamp == 0) {
    // when timestamp wraps around, we reset the vector to 0 and increment again.
    timestamp++;
    std::fill(functorSet.begin(), functorSet.end(), 0);
  }
  ASS(std::all_of(functorSet.begin(), functorSet.end(), [&](prune_t x) { return x < timestamp; }));

  for (unsigned i = 0; i < _M->length(); i++)
    functorSet[(*_M)[i]->functor()] = timestamp;

  for (unsigned j = 0; j < _L->length(); j++)
    if (functorSet[(*_L)[j]->functor()] != timestamp)
      return true;

  return false;
} // SATSubsumptionAndResolution::pruneSubsumptionResolution

void SATSubsumptionAndResolution::addBinding(BindingsManager::Binder* binder,
                                             unsigned i,
                                             unsigned j,
                                             bool polarity,
                                             bool isNullary)
{
  ASS(binder || isNullary)
  ASS(i < _m)
  ASS(j < _n)
  ASS_EQ((*_L)[i]->functor(), (*_M)[j]->functor())
  ASS_EQ((*_L)[i]->polarity() == (*_M)[j]->polarity(), polarity)
  subsat::Var satVar = _solver.new_variable();
#if PRINT_CLAUSES_SUBS
  cout << satVar << " -> (" << (*_L)[i]->toString() << " " << (*_M)[j]->toString() << " " << (polarity ? "+" : "-") << ")" << endl;
#endif
  _matchSet.addMatch(i, j, polarity, satVar);
  if (!isNullary)
    _bindingsManager.commit_bindings(*binder, satVar);
} // SATSubsumptionAndResolution::addBinding

bool SATSubsumptionAndResolution::checkAndAddMatch(Literal* l_i,
                                                   Literal* m_j,
                                                   unsigned i,
                                                   unsigned j,
                                                   bool polarity)
{
  ASS(l_i)
  ASS(m_j)
  ASS_EQ((*_L)[i], l_i)
  ASS_EQ((*_M)[j], m_j)
  ASS_EQ(l_i->functor(), m_j->functor())
  ASS_EQ(l_i->polarity() == m_j->polarity(), polarity)

  bool match = false;
  {
    auto binder = _bindingsManager.start_binder();
    if (MatchingUtils::matchArgs(l_i, m_j, binder)) {
      addBinding(&binder, i, j, polarity, false);
      match = true;
    }
  }
  if (l_i->isEquality()) {
    auto binder = _bindingsManager.start_binder();
    if (MatchingUtils::matchReversedArgs(l_i, m_j, binder)) {
      addBinding(&binder, i, j, polarity, false);
      match = true;
    }
  }
  return match;
} // SATSubsumptionAndResolution::checkAndAddMatch

bool SATSubsumptionAndResolution::fillMatchesS()
{
  ASS(_L)
  ASS(_M)
  ASS_G(_m, 0)
  ASS_G(_n, 0)
  ASS_EQ(_matchSet._m, _m)
  ASS_EQ(_matchSet._n, _n)

  Literal* l_i, * m_j;

  // number of matches found is equal to the number of variables in the SAT solver
  for (unsigned i = 0; i < _m; ++i) {
    l_i = _L->literals()[i];
    bool foundMatch = false;

    for (unsigned j = 0; j < _n; ++j) {
      m_j = _M->literals()[j];
      if (l_i->functor() != m_j->functor() || l_i->polarity() != m_j->polarity()) {
        continue;
      }
      if (l_i->arity() == 0) {
        ASS(m_j->arity() == 0)
        ASS(l_i->functor() == m_j->functor())
        addBinding(nullptr, i, j, true, true);
        foundMatch = true;
        continue;
      }
      // it is important that foundMatch is "or-ed" after calling the function. Otherwise the function might not be called.
      // foundMatch |= checkAndAddMatch(l_i, m_j, i, j, true); is NOT correct.
      foundMatch = checkAndAddMatch(l_i, m_j, i, j, true) || foundMatch;
    } // for (unsigned j = 0; j < _n; ++j)

    if (!foundMatch) {
      _subsumptionImpossible = true;
      return false;
    } // if (!foundPositiveMatch)
  }   // for (unsigned i = 0; i < _m; ++i)

  return true;
} // SATSubsumptionAndResolution::fillMatchesS()

void SATSubsumptionAndResolution::fillMatchesSR(unsigned litToRemove)
{
  ASS(_L)
  ASS(_M)
  ASS_G(_m, 0)
  ASS_G(_n, 0)
  ASS_EQ(_matchSet._m, _m)
  ASS_EQ(_matchSet._n, _n)

  // stores whether on all the literals in L there is a negative match in M
  bool clauseHasNegativeMatch = false;

  // the first literal in L that only has a negative match (no positive)
  Literal* firstOnlyNegativeMatch = nullptr;

  for (unsigned i = 0; i < _m; ++i) {
    Literal* l_i = _L->literals()[i];

    // does lᵢ have a positive match in M?
    bool literalHasPositiveMatch = false;
    // does lᵢ have a negative match in M?
    bool literalHasNegativeMatch = false;

    for (unsigned j = 0; j < _n; ++j) {
      Literal* m_j = _M->literals()[j];
      if (l_i->functor() != m_j->functor())
        continue;
      if (l_i->arity() == 0) {
        ASS(m_j->arity() == 0)
        ASS(l_i->functor() == m_j->functor())
        if (l_i->polarity() == m_j->polarity()) {
          addBinding(nullptr, i, j, true, true);
          literalHasPositiveMatch = true;
          continue;
        }
        if (litToRemove != 0xFFFFFFFF && j != litToRemove)
          continue;
        addBinding(nullptr, i, j, false, true);
        clauseHasNegativeMatch = true;
        literalHasNegativeMatch = true;
        continue;
      }

      if (l_i->polarity() == m_j->polarity()) {
        // it is important that foundPositiveMatch is "or-ed" after calling the function. Otherwise the function might not be called.
        // foundPositiveMatch |= checkAndAddMatch(l_i, m_j, i, j, true); is NOT correct.
        literalHasPositiveMatch = checkAndAddMatch(l_i, m_j, i, j, true) || literalHasPositiveMatch;
        continue;
      }
      // check negative polarity matches
      // same comment as above
      if (litToRemove != 0xFFFFFFFF && j != litToRemove)
          continue;
      literalHasNegativeMatch = checkAndAddMatch(l_i, m_j, i, j, false) || literalHasNegativeMatch;
      clauseHasNegativeMatch |= literalHasNegativeMatch;
    } // for (unsigned j = 0; j < _nInstanceLits; ++j)

    // Check whether subsumption and subsumption resolution are possible
    if (!literalHasPositiveMatch) {
      _subsumptionImpossible = true;
      if (!literalHasNegativeMatch) {
        // no positive nor negative matches found
        _srImpossible = true;
        return;
      }
      // TODO try matching?
      if (!firstOnlyNegativeMatch)
        firstOnlyNegativeMatch = l_i;
      else if (firstOnlyNegativeMatch->header() != l_i->header()) {
        // there are two literals in L with only negative polarity matches, but they have a different functor
        // therefore, subsumption resolution is impossible
        _srImpossible = true;
        return;
      }
    } // if (!foundPositiveMatch)
  }   // for (unsigned i = 0; i < _nBaseLits; ++i)

  // If there are no negative matches, then the SR is not possible
  if (!clauseHasNegativeMatch)
    _srImpossible = true;

} // SATSubsumptionAndResolution::fillMatchesSR()

bool SATSubsumptionAndResolution::cnfForSubsumption()
{
  ASS(_L)
  ASS(_M)
  ASS_GE(_matchSet.allMatches().size(), _L->length())
  ASS_G(_L->length(), 0)
  ASS(!_subsumptionImpossible)

  _matchSet.indexMatrix();

  Solver& solver = _solver;

  /**** Completeness ****/
  // Build clauses stating that lᵢ must be matched to at least one corresponding mⱼ.
  // ⋀ᵢ ⋁ⱼ b⁺ᵢⱼ
  for (unsigned i = 0; i < _m; ++i) {
    solver.constraint_start();
    for (Match match : _matchSet.getIMatches(i)) {
      if (match.polarity) {
        solver.constraint_push_literal(match.var);
      }
    }
    auto handle = solver.constraint_end();
    solver.add_clause_unsafe(handle);
  }

  /**** Multiplicity conservation ****/
  // Here we store the AtMostOne constraints saying that each instance literal may be matched at most once.
  // Each instance literal can be matched by at most 2 boolean vars per base literal (two orientations of equalities).
  // NOTE: instance constraints cannot be packed densely because we only know their shape at the end.
  // uint32_t const instance_constraint_maxsize = 2 * base_len;
  // ⋀ⱼ AMO({b⁺ᵢⱼ ∣ i ∈ {1,...,m}})
  for (unsigned j = 0; j < _n; ++j) {
    solver.constraint_start();
    for (Match match : _matchSet.getJMatches(j)) {
      if (match.polarity) {
        solver.constraint_push_literal(match.var);
      }
    }
    auto handle = solver.constraint_end();
    solver.add_atmostone_constraint_unsafe(handle);
  }

  return !solver.inconsistent();
} // SATSubsumptionAndResolution::cnfForSubsumption()


/// @brief a vector used to store the sat variables that are subjected to the at most one constraint (will hold the cⱼ).
/// The unsigned value is the index of the literal in the instance clause
static std::vector<pair<unsigned, subsat::Var>> atMostOneVars;
/**
 * Indirect encoding of the sat subsumption resolution
 *
 * Let the base clause be L and the instance clause be M.
 *  - We name the literals in L as L1, L2, ..., Ln
 *  - We name the literals in M as M1, M2, ..., Mk
 *
 * We introduce the variables
 *  - b⁺ᵢⱼ if lᵢ has a substitution σ such that σ(lᵢ) = mⱼ
 *    we will say that b⁺ᵢⱼ as a positive polarity
 *  - b⁻ᵢⱼ if lᵢ has a substitution σ such that σ(lᵢ) = ¬mⱼ
 *    we will say that b⁻ᵢⱼ as a negative polarity
 *  - cⱼ if mⱼ is matched by at least one negative polarity variable b⁻ᵢⱼ
 *    if only one b⁻ᵢⱼ exist for some j, then we do not introduce cⱼ,
 *    we use b⁻ᵢⱼ directly
 * b⁺ᵢⱼ and b⁻ᵢⱼ are mutually existentially exclusive, therefore, we only introduce one variable bᵢⱼ for the sat solver, and remember its polarity.
 * It may be that both b⁺ᵢⱼ and b⁻ᵢⱼ do not exist.
 *
 * In the case that there exists only one b⁻ᵢⱼ for some j, we do not introduce cⱼ. It would be useless to add the constraints cⱼ ⇔ b⁻ᵢⱼ.
 *
 * We introduce the following subsumption resolution clauses clauses:
 *  - Encoding of cⱼ
 *    cⱼ is true iff mⱼ is matched by at least one negative polarity variable b⁻ᵢⱼ
 *    ∀ j. cⱼ ⇔ (∃ i b⁻ᵢⱼ)
 *    ∀ j. (¬cⱼ ∨ b⁻₁ⱼ ∨ ... ∨ b⁻ₙⱼ)
 *       ∧ (cⱼ ∨ ¬b⁻₁ⱼ) ∧ ... ∧ (cⱼ ∨ ¬b⁻ₙⱼ)
 *  - Existence
 *      c₁ ∨ c₂ ∨ ... ∨ cₙ
 *  - Uniqueness
 *      AMO(c₁, c₂, ..., cₙ)
 *  - Completeness
 *      b₁₁ ∨ ... ∨ b₁ₙ
 *    ∧ b₂₁ ∨ ... ∨ b₂ₙ
 *    ∧ ...
 *    ∧ bₘ₁ ∨ ... ∨ bₘₙ
 *  - Coherence
 *    if cⱼ is true, then there is no positive polarity match b⁺ᵢⱼ
 *      ∀ ij. (¬cⱼ ∨ ¬b⁺ᵢⱼ)
 *  - Substitution validity
 *    bᵢⱼ implies a certain substitution is valid
 *      ∀ ij. bᵢⱼ ⇒ (σ(lᵢ) = mⱼ ∨ σ(lᵢ) = ¬mⱼ)
 *    (This rule is enforced by the match set in the SAT solver)
 */
bool SATSubsumptionAndResolution::cnfForSubsumptionResolution()
{
  ASS(_L)
  ASS(_M)
  // This should be pruned when filling the match set.
  ASS_GE(_matchSet.allMatches().size(), _L->length())

  atMostOneVars.clear();
  _matchSet.indexMatrix();
  Solver& solver = _solver;

  /**** Existence ****/
  // ⋁ⱼ cⱼ
#if PRINT_CLAUSE_COMMENTS_SUBS
  cout << "Existence" << endl;
#endif
#if PRINT_CLAUSES_SUBS
  string s = "";
#endif
  solver.constraint_start();
  for (unsigned j = 0; j < _n; ++j) {
    // Do not add useless variables in the solver
    // It would compromise correctness since a cⱼ without binding could be true and satisfy c₁ ∨ c₂ ∨ ... ∨ cₙ without making any other clause false.
    if (!_matchSet.hasNegativeMatchJ(j))
      continue;
    // Do not create a new SAT variable if there is only one negative match to mⱼ
    // match_found is true if we found a match, false if we found a different number than one
    // We also know that there is at least one match
    bool one_match_found = false;
    subsat::Var negative_match_var = subsat::Var(0xffffffff);
    for (Match match : _matchSet.getJMatches(j)) {
      if (!match.polarity) {
        one_match_found = !one_match_found;
        if (!one_match_found)
          break;
        negative_match_var = match.var;
      }
    }
    if (one_match_found) {
      ASS(!_matchSet.getMatchForVar(negative_match_var).polarity)
      // we have exactly one negative match. No need to create a new variable
      atMostOneVars.push_back(make_pair(j, negative_match_var));
      solver.constraint_push_literal(negative_match_var);
#if PRINT_CLAUSES_SUBS
      s += Int::toString(negative_match_var.index()) + " ∨ ";
#endif
      continue;
    }
    subsat::Var c_j = solver.new_variable();
    atMostOneVars.push_back(make_pair(j, c_j));
    solver.constraint_push_literal(c_j);
#if PRINT_CLAUSES_SUBS
    s += Int::toString(c_j.index()) + " ∨ ";
#endif
  } // for (unsigned j = 0; j < _n; ++j)

  ASS(!atMostOneVars.empty())
  auto build = solver.constraint_end();
  solver.add_clause_unsafe(build);
#if PRINT_CLAUSES_SUBS
  cout << s.substr(0, s.size() - 3) << endl;
#endif

  /**** Uniqueness ****/
  // AMO(cⱼ | j ∈ {1,...,n})
#if PRINT_CLAUSE_COMMENTS_SUBS
  cout << "Uniqueness" << endl;
#endif
#if PRINT_CLAUSES_SUBS
  cout << "AtMostOne(";
#endif
  solver.constraint_start();
  for (auto var : atMostOneVars) {
    solver.constraint_push_literal(var.second);
#if PRINT_CLAUSES_SUBS
    cout << var.second;
    if (var != atMostOneVars.back()) {
      cout << ", ";
    }
#endif
  } // for (auto var : atMostOneVars)
  build = solver.constraint_end();
  solver.add_atmostone_constraint(build);
#if PRINT_CLAUSES_SUBS
  cout << ")" << endl;
#endif

  /**** Completeness ****/
  // ⋀ᵢ ⋁ⱼ bᵢⱼ
#if PRINT_CLAUSE_COMMENTS_SUBS
  cout << "Completeness" << endl;
#endif
  for (unsigned i = 0; i < _m; ++i) {
    solver.constraint_start();
    Slice<Match> matches = _matchSet.getIMatches(i);
    for (Match match : matches) {
#if PRINT_CLAUSES_SUBS
      cout << match.var;
      if (match != matches.back()) {
        cout << " ∨ ";
      }
#endif
      solver.constraint_push_literal(match.var);
    }
    auto build = solver.constraint_end();
    solver.add_clause_unsafe(build);
#if PRINT_CLAUSES_SUBS
    cout << endl;
#endif
  } // for (unsigned i = 0; i < _m; ++i)

  /**** Coherence ****/
  // ⋀ⱼ ⋀ᵢ (¬cⱼ ∨ ¬b⁺ᵢⱼ)
#if PRINT_CLAUSE_COMMENTS_SUBS
  cout << "Coherence" << endl;
#endif
  for (auto var : atMostOneVars) {
    unsigned j = var.first;
    subsat::Var c_j = var.second;
    if (_matchSet.hasPositiveMatchJ(j)) {
      for (Match match : _matchSet.getJMatches(j)) {
        subsat::Var b_ij = match.var;
        if (match.polarity) {
          solver.constraint_start();
          solver.constraint_push_literal(~c_j);
          solver.constraint_push_literal(~b_ij);
          build = solver.constraint_end();
          solver.add_clause_unsafe(build);
#if PRINT_CLAUSES_SUBS
          cout << ~c_j << " ∨ " << ~b_ij << endl;
#endif
        } // if (match.polarity)
      } // for (Match match : _matchSet.getJMatches(j))
    } // if (_matchSet.hasPositiveMatchJ(j))
  } // for (auto var : atMostOneVars)

#if PRINT_CLAUSE_COMMENTS_SUBS
  cout << "Stucturality" << endl;
#endif
  /**** Encoding of cⱼ ****/
  // -> cⱼ is true iff mⱼ is matched by at least one negative polarity variable b⁻ᵢⱼ
  //    ⋀ⱼ [(¬cⱼ ∨ ⋁ᵢ b⁻ᵢⱼ)
  //        ∧ ⋀ᵢ (cⱼ ∨ ¬b⁻ᵢⱼ)]
  for (auto& pair : atMostOneVars) {
    unsigned j = pair.first;
    subsat::Var c_j = pair.second;
    if (_matchSet.isMatchVar(c_j))
      // This is not a cⱼ. We cannot encode it as a cⱼ
      continue;
    Slice<Match> matches = _matchSet.getJMatches(j);
    solver.constraint_start();
    //   ¬cⱼ ∨ b⁻₁ⱼ ∨ ... ∨ b⁻ₘⱼ
    solver.constraint_push_literal(~c_j);
#if PRINT_CLAUSES_SUBS
    cout << ~c_j;
#endif
    for (Match match : matches) {
      if (!match.polarity) {
        solver.constraint_push_literal(match.var);
#if PRINT_CLAUSES_SUBS
        cout << " ∨ " << match.var;
#endif
      }
    } // for (Match match : matches)
#if PRINT_CLAUSES_SUBS
    cout << endl;
#endif
    auto build = solver.constraint_end();
    solver.add_clause_unsafe(build);
    //   (cⱼ ∨ ¬b⁻ᵢⱼ) ∧ ... ∧ (cⱼ ∨ ¬b⁻ₘⱼ)
    for (Match match : matches) {
      if (!match.polarity) {
        solver.constraint_start();
        solver.constraint_push_literal(c_j);
        solver.constraint_push_literal(~match.var);
        auto build = solver.constraint_end();
        solver.add_clause_unsafe(build);
#if PRINT_CLAUSES_SUBS
        cout << c_j << " ∨ " << ~match.var << endl;
#endif
      } // if (!match.polarity)
    } // for (Match match : matches)
  } // for (auto &pair : atMostOneVars)
  return !solver.inconsistent();
} // cnfForSubsumptionResolution

Clause* SATSubsumptionAndResolution::getSubsumptionResolutionConclusion(Clause* M,
  Literal* m_j,
  Clause* L)
{
  RStack<Literal*> resLits;

  int mlen = M->length();
  for (int i = 0; i < mlen; i++) {
    Literal* curr = (*M)[i];
    if (curr == m_j)
      continue;
    resLits->push(curr);
  }

  return Clause::fromStack(*resLits,SimplifyingInference2(InferenceRule::SUBSUMPTION_RESOLUTION, M, L));
}

Clause* SATSubsumptionAndResolution::generateConclusion()
{
  ASS(_L)
  ASS(_M)
  ASS(_m > 0)
  ASS(_n > 0)
  ASS_EQ(_matchSet._m, _m)
  ASS_EQ(_matchSet._n, _n)
  ASS(_model.size() > 0)
  ASS_GE(_matchSet.allMatches().size(), _L->length())

  // Provided the solution of the sat solver, we can create the conclusion clause
#if VDEBUG
  unsigned j = INVALID;
  // Check that there is only one negative polarity match to j inside the model
  for (subsat::Lit lit : _model) {
    if (lit.is_positive()) {
      // matches can be null if the variable is a cⱼ
      if (_matchSet.isMatchVar(lit.var())) {
        Match match = _matchSet.getMatchForVar(lit.var());
        if (!match.polarity) {
          if (j == INVALID)
            j = match.j;
          else
            ASS(j == match.j)
        }
      }
    }
  }
#endif
  unsigned toRemove = INVALID;
  // find the negative polarity match to j inside the model
#if PRINT_CLAUSES_SUBS
  cout << "Model: ";
  for (subsat::Lit lit : _model) {
    cout << lit << " ";
  }
  cout << endl;
#endif
  for (subsat::Lit lit : _model) {
    if (lit.is_positive()) {
      // matches can be null if the variable is a cⱼ
      if (_matchSet.isMatchVar(lit.var())) {
        Match match = _matchSet.getMatchForVar(lit.var());
        if (!match.polarity) {
          toRemove = match.j;
          break;
        }
      }
    }
  }
  ASS_EQ(_n, _M->size())
  ASS(toRemove != INVALID)
  return SATSubsumptionAndResolution::getSubsumptionResolutionConclusion(_M, (*_M)[toRemove], _L);
} // SATSubsumptionResolution::generateConclusion

bool SATSubsumptionAndResolution::checkSubsumption(Clause* L,
                                                   Clause* M,
                                                   bool setSR)
{
  ASS(L)
  ASS(M)

  loadProblem(L, M);

  // Fill the matches
  if (setSR) {
    _srImpossible = pruneSubsumptionResolution();
    // WARNING!!! This assumes that the check for subsumption resolution is stronger than
    // the check for subsumption.
    _subsumptionImpossible = _srImpossible || pruneSubsumption();
    if (_srImpossible) {
      ASS(_subsumptionImpossible);
      return false;
    }
    else {
      ASS(!_srImpossible)
      fillMatchesSR();
    }
    if (_subsumptionImpossible)
      return false;
  } // if (setSR)
  else if (pruneSubsumption() || !fillMatchesS())
    return false;

  ASS_GE(_matchSet.allMatches().size(), _L->length())

  // Create the constraints for the sat solver
  if (!cnfForSubsumption())
    return false;

  // Solve the SAT problem
  _solver.theory().setBindings(&_bindingsManager);
  auto const result = _solver.solve();

  bool const subsumed = result == subsat::Result::Sat;

  return subsumed;
} // SATSubsumptionAndResolution::checkSubsumption

Clause* SATSubsumptionAndResolution::checkSubsumptionResolution(Clause* L,
                                                                Clause* M,
                                                                bool usePreviousSetUp)
{
  ASS(L)
  ASS(M)
  if (usePreviousSetUp) {
    ASS(_L == L)
    ASS(_M == M)
    if (_srImpossible) {
#if PRINT_CLAUSES_SUBS
      cout << "SR impossible" << endl;
#endif
      return nullptr;
    }
    ASS_GE(_matchSet.allMatches().size(), _L->length())
    // do not clear the variables and bindings
    _solver.clear_constraints();
  }
  else {
    loadProblem(L, M);
    if (pruneSubsumptionResolution()) {
#if PRINT_CLAUSES_SUBS
      cout << "SR pruned" << endl;
#endif
      return nullptr;
    }
    fillMatchesSR();

    if (_srImpossible) {
#if PRINT_CLAUSES_SUBS
      cout << "SR impossible" << endl;
#endif
      return nullptr;
    }
  }

  // set up the clauses
  if (!cnfForSubsumptionResolution()) {
#if PRINT_CLAUSES_SUBS
    cout << "CNF building failed" << endl;
#endif
    return nullptr;
  }

  // Solve the SAT problem
  if (_solver.theory().empty()) {
    // -> bᵢⱼ implies a certain substitution is valid
    //    for each i, j : bᵢⱼ ⇒ (σ(lᵢ) = mⱼ ∨ σ(lᵢ) = ¬mⱼ)
    // These constraints are created in the fillMatches() function by filling the _bindingsManager
    _solver.theory().setBindings(&_bindingsManager);
  }
  Clause* conclusion = nullptr;

  auto const result = _solver.solve();

  if (result == subsat::Result::Sat) {
#if PRINT_CLAUSES_SUBS
    cout << "SAT solver succeeded" << endl;
#endif
    _model.clear();
    _solver.get_model(_model);
    conclusion = generateConclusion();
  }
#if PRINT_CLAUSES_SUBS
  else
    cout << "SAT solver failed (" << result << ")" << endl;
#endif

  // If the problem is SAT, then generate the conclusion clause
  return conclusion;
} // SATSubsumptionAndResolution::checkSubsumptionResolution

bool SATSubsumption::SATSubsumptionAndResolution::checkSubsumptionResolutionWithLiteral(Kernel::Clause* L, Kernel::Clause* M, unsigned resolutionLiteral)
{
  loadProblem(L, M);
  if (pruneSubsumptionResolution()) {
    return false;
  }
  fillMatchesSR(resolutionLiteral);

  if (_srImpossible) {
    return false;
  }

  // set up the clauses
  if (!cnfForSubsumptionResolution()) {
    return false;
  }

  ASS(_solver.theory().empty())
  // since we cannot use the previous setup, the bindings should be cleared during loading of the problem
  _solver.theory().setBindings(&_bindingsManager);

  _model.clear();
  return (_solver.solve() == subsat::Result::Sat);
}

SimpleSubstitution SATSubsumption::SATSubsumptionAndResolution::getBindingsForSubsumptionResolutionWithLiteral()
{
  SimpleSubstitution subst;
  _solver.get_model(_model);
  for(auto lit : _model) {
    if(lit.is_negative())
      continue;
    // there could be non-match vars with the (only existing) indirect encoding
    if(!_matchSet.isMatchVar(lit.var()))
      continue;
    // extra guard for below
    if(lit.var().index() >= _bindingsManager.size())
      continue;
    auto b = _bindingsManager.get_bindings(lit.var());
    // TODO put a bit in the Match to indicate which way an equation is oriented
    for(auto binding_index = b.index; binding_index < b.end(); binding_index++) {
      auto [var, term] = _bindingsManager.get_single_binding(binding_index);
      ALWAYS(subst.bind(var, term))
    }
  }
  return subst;
}
