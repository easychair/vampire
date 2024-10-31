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
 * Defines class SATSubsumptionAndResolution.
 */
#ifndef SAT_SUBSUMPTION_RESOLUTION_HPP
#define SAT_SUBSUMPTION_RESOLUTION_HPP

#include "Kernel/Clause.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Lib/Slice.hpp"

#include "./subsat/subsat.hpp"

namespace SATSubsumption {

/**
 * Class implementing the simplifying rules of subsumption and subsumption resolution using a SAT solver.
 *
 * The subsumption and subsumption resolution are described in the papers:
 * - 2022: "First-Order Subsumption via SAT Solving." by Jakob Rath, Armin Biere and Laura Kovács
 * - 2023: "SAT-Based Subsumption Resolution" by Robin Coutelier, Jakob Rath, Michael Rawson and Laura Kovács
 * - 2024: "SAT Solving for Variants of First-Order Subsumption" by Robin Coutelier, Jakob Rath, Michael Rawson, Armin Biere and Laura Kovács
 */
class SATSubsumptionAndResolution {
#if VDEBUG
  // Make it public to allow unit testing
public:
#else
private:
#endif
  using Solver = subsat::Solver;
  using BindingsManager = typename Solver::BindingsManager;

  /**
   * A Match represents a binding between two literals l_i and m_j of different clauses L and M.
   * The binding can be either positive or negative
   */
  struct Match {
    /// @brief The index of the literal in L (base clause for subsumption resolution)
    unsigned i;
    // The index of the literal in M (instance clause for subsumption resolution)
    unsigned j;
    /// @brief The variable associated in the sat solver
    subsat::Var var;
    /// @brief The polarity of the match (true for positive, false for negative)
    bool polarity;

    Match() : i(0),
              j(0),
              var(subsat::Var(0)),
              polarity(true) {}

    Match(unsigned baseLitIndex,
          unsigned instanceLitIndex,
          bool isPositive,
          subsat::Var satVar) : i(baseLitIndex),
                                j(instanceLitIndex),
                                var(satVar),
                                polarity(isPositive) {}

    bool operator==(const Match &other) const
    {
      return i == other.i && j == other.j && polarity == other.polarity && var == other.var;
    }

    bool operator!=(const Match &other) const
    {
      return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& out, Match const& m)
    {
      return out << "Match(" << m.i << ", " << m.j << ", " << m.polarity << + ", " << m.var.index() << ")";
    }
  };

  /**
   * A Match set is a facility to store matches, and allows to enumerate them either according to the clause L or M indices.
   * The Match Set, when filled contains all the bindings possible between the clauses L and M
   * The match set can be abstracted as a sparse matrix holding the associated sat variable and polarity of the matches
   *
   * @example For example, here is a table of matches (not necessarily coherent with the subsumption resolution problem)
   *
   *      i=0  i=1  i=2  i=3  i=4
   * j=0 | 0+ |    |    | 3- |    |
   *     --------------------------
   * j=1 |    | 1+ |    | 2+ |    |
   *     --------------------------
   * j=2 | 4- |    |    |    |    |
   *     --------------------------
   * j=3 |    |    | 5+ |    | 6- |
   *     --------------------------
   * j=4 |    | 7+ |    | 8- |    |
   *
   * The Match set allows to get a line or column of the matrix
   */
  struct MatchSet {
    USE_ALLOCATOR(SATSubsumptionAndResolution::MatchSet);
    friend struct Match;

    /// @brief Metadata remembering whether some positive match or negative match was found for each literal in M
    /// @remark
    /// This information only needs 2 bits
    /// - 00 -> no match
    /// - 01 -> positive match
    /// - 10 -> negative match
    /// - 11 -> both positive and negative match
    /// Since we only need 2 bits, 4 values fit in one byte
    /// The interest value for i is therefore
    /// state[j / 4] & (0b11 << (2 * (j % 4))) >> (2 * (j % 4))
    /// @example For example, if we have 5 literals in M, the state will be
    /// 0b 00 10 01 11 | 0b 00
    /// and the interest value for j=2 will be state[0]
    /// 0b 00 10 01 11
    /// After passing the mask 0b 00 11 0 00 (= (0b11 << (2 * (2 % 4)))) )
    /// 0b 00 10 00 00
    /// After shifting 4 bits to the right (>> (2 * (j % 4)))
    /// 0b 00 00 00 10
    /// We know that m_2 does not have a positive match, but has a negative match
    std::vector<uint8_t> _jStates;

    /// @brief number of literals in L
    unsigned _m;
    /// @brief number of literals in M
    unsigned _n;

    /// @brief contains the list of all the the matches indexed in the match set
    /// Ordered by the index 'i', then by 'j'
    /// If the properties of AddMatch are respected, then _matchesByI[i] is the match with the sat variable i
    std::vector<Match> _matchesByI;

    /// @brief same as _matchesByI, but ordered by the index 'j', then by 'i'
    /// filled by fillMatchesByColumn based on _matchesByI
    std::vector<Match> _matchesByJ;

    /// @brief indices allowing access into _matchesByI
    /// if you wish to access row i, [_indexI[i], _indexJ[i + 1]) is the range required
    std::vector<unsigned> _indexI;

    /// @brief indices allowing access into _matchesByJ
    /// if you wish to access column j, [_indexJ[j], _indexJ[j + 1]) is the range required
    std::vector<unsigned> _indexJ;

    /// @brief all matches, ordered by row
    std::vector<Match> &allMatches()
    {
      return _matchesByI;
    }

    /**
     * @brief get the matches at row i
     * @param i the index of the row
     * @return a slice of the matches at row i
     */
    Slice<Match> getIMatches(unsigned i)
    {
      ASS_L(i, _indexI.size())
      return Slice<Match>(
          &_matchesByI[_indexI[i]],
          &_matchesByI[_indexI[i + 1]]);
    }

    /**
     * @brief get the matches at column j
     * @param j the index of the column
     * @return a slice of the matches at column j
     */
    Slice<Match> getJMatches(unsigned j)
    {
      ASS_L(j, _indexJ.size())
      return Slice<Match>(
          &_matchesByJ[_indexJ[j]],
          &_matchesByJ[_indexJ[j + 1]]);
    }

    /**
     * Resizes the match matrix to the given size and clears the matrix
     * Guarantees that the a match can be added at index ( @b m - 1, @b n - 1 )
     *
     * @param m the new length of L
     * @param n the new length of M
     *
     * @warning the allocated matches will remain accessible but will be no longer be reserved. It would therefore be preferable to not keep the matches after calling this function.
     */
    void resize(unsigned m, unsigned n)
    {
      ASS_G(m, 0)
      ASS_G(n, 0)
      _jStates.resize(n / 4 + 1, 0);
      _m = m;
      _n = n;
    }

    /**
     * Adds a new match to the set
     *
     * @pre the match must be valid (i.e. @b i < _m and @b j < _n)
     * @pre the sar variable @b var must follow the one that was added before, or be 0 if no match was added before
     *
     * @param i the index of the literal in L
     * @param j the index of the literal in M
     * @param polarity the polarity of the match
     * @param var the sat variable associated with the match
     * @return the newly added match
     */
    Match addMatch(unsigned i, unsigned j, bool polarity, subsat::Var var)
    {
      ASS(i < _m)
      ASS(j < _n)
      // Make sure that the variables are pushed in order.
      // Otherwise would break the property that _matches[i] is the match associated
      //    to the sat variable i
      ASS_EQ(var.index(), _matchesByI.size())

      Match match(i, j, polarity, var);
      _matchesByI.push_back(match);

      // update the match state
      // the wizardry is explained in the header file
      _jStates[j / 4] |= 1 << (2 * (j % 4) + !polarity);
      return match;
    }

    /**
     * Populate _matchesByJ and indices based on _matchesByI
     * Entries of _matchesByI are inserted in row-order by addMatch
     */
    void indexMatrix();

    /**
     * Returns true if the j-th literal in M has a positive match in the set
     *
     * @param j the index of the literal in M
     * @return whether m_j has a positive match in the set
     */
    bool hasPositiveMatchJ(unsigned j);

    /**
     * Returns true if the j-th literal in M has a positive match in the set
     *
     * @param j the index of the literal in M
     * @return whether m_j has a negative match in the set
     */
    bool hasNegativeMatchJ(unsigned j);

    /**
     * Checks whether the sat variable @b v is linked to a match
     *
     * @pre Assumes that the matches are linked to variables in an increasing order
     * (the first match being associated with variable 0, then 1, and so on)
     * No leaps are authorized
     * @param v the variable
     * @return true if the variable is a match variable
     */
    bool isMatchVar(subsat::Var v)
    {
      return v.index() < _matchesByI.size();
    }

    /**
     * Returns the match for a given sat variable
     *
     * @pre Assumes that the matches are linked to variables in an increasing order
     * (the first match being associated with variable 0, then 1, and so on)
     * No leaps are authorized
     * @param v the variable
     * @return the match for the variable, or nullptr if no match exists
     */
    Match getMatchForVar(subsat::Var v)
    {
      ASS(isMatchVar(v))
      return _matchesByI[v.index()];
    }

    /**
     * Clears the match set
     *
     * @warning the allocated matches will remain accessible but will be no longer be reserved. I would therefore be preferable to not keep the matches after calling this function.
     */
    void clear()
    {
      _matchesByI.clear();
      _matchesByJ.clear();
      _jStates.clear();
      _indexI.clear();
      _indexJ.clear();
    }
  };


  /* Variables */

  /// @brief the base clause (side premise)
  Kernel::Clause *_L;
  /// @brief the instance clause (main premise)
  Kernel::Clause *_M;
  /// @brief the number of literals in the base clause (side premise)
  unsigned _m;
  /// @brief the number of literals in the instance clause (main premise)
  unsigned _n;
  /// @brief the SAT solver
  Solver _solver;
  /// @brief the bindings manager
  BindingsManager _bindingsManager;
  /// @brief the match set used to store the matches between the base and instance clauses
  MatchSet _matchSet;
  /// @brief model of the SAT solver
  std::vector<subsat::Lit> _model;

  /// @brief remembers if the fillMatchesSR concluded that subsumption is impossible
  bool _subsumptionImpossible;
  /// remembers if the fillMatchesSR concluded that subsumption resolution is impossible
  bool _srImpossible;

  /// @brief temporary storage, used by pruneSubsumption and pruneSubsumptionResolution
  /// invariant: for all x in _tmpStorage, x <= _tmpTimestamp
  using prune_t = unsigned;
  std::vector<prune_t> _pruneStorage;
  prune_t _pruneTimestamp = 0;

  /* Methods */
  /**
   * Sets up the problem and cleans the match set and bindings
   *
   * @param L the base clause
   * @param M the instance clause
   */
  void loadProblem(Kernel::Clause *L,
                   Kernel::Clause *M);

  /**
   * Heuristically predicts whether subsumption or subsumption resolution will fail.
   * This method should be fast.
   *
   * @return true if subsumption is impossible
   */
  bool pruneSubsumption();

  /**
   * Heuristically predicts whether subsumption resolution will fail.
   * This method should be fast
   *
   * @return true if subsumption resolution is impossible
   */
  bool pruneSubsumptionResolution();

  /**
   * Adds one binding to the SAT solver and the match set
   *
   * @param binder the binder between the literals
   * @param varNumber the variable number of the binder
   * @param i the index of the literal in the base clause
   * @param j the index of the literal in the instance clause
   * @param polarity the polarity of the match
   * @param isNullary whether the literals l_i and m_j are nullary. If so, no binder is added to the SAT solver
   */
  void addBinding(BindingsManager::Binder *binder,
                  unsigned i,
                  unsigned j,
                  bool polarity,
                  bool isNullary);

  /**
   * Adds the clauses for the subsumption problem to the sat solver
   *
   * @pre _L and _M must be set in the checker
   * @pre the Match set is already filled
   * @return false if no solution is possible and true if there may exist a solution.
   */
  bool cnfForSubsumption();

  /**
   * Function type for encoding the subsumption resolution problem to the sat solver
   *
   * @param checker the SATSubsumptionAndResolution object
   * @return false if no solution is possible and true if there may exist a solution.
   *
   * @details The encoding must use the variables in the MatchSet. It may add new variables
   * to the SAT solver, but the variables in the MatchSet must be used to interpret the model
   * and build the conclusion.
   */

  /**
   * Adds the clauses for the subsumption resolution problem to the sat solver
   *
   * @remark The BindingsManager is not required to be set up in this method.
   * @pre _L and _M must be set in the checker
   * @pre the Match set is already filled
   * @return false if no solution is possible and true if there may exist a solution.
   */
  bool cnfForSubsumptionResolution();

  /**
   * Checks whether there exists a substitution from the literals @b l_i and @b m_j.
   * If there exists a substitution, it is added to the match set.
   *
   * @pre Assumes that the polarity and functors of the match was already checked.
   *
   * @param l_i one literal of the base clause
   * @param m_j one literal of the instance clause
   * @param i the index of the literal in the base clause
   * @param j the index of the literal in the instance clause
   * @param polarity the polarity of the match
   * @return true if the literals are matched
   */
  bool checkAndAddMatch(Kernel::Literal *l_i,
                        Kernel::Literal *m_j,
                        unsigned i,
                        unsigned j,
                        bool polarity);

  /**
   * Fills the match set and the bindings manager with all the possible positive bindings between the literals of L and M.
   *
   * @return false if no subsumption solution is possible, the number of sat variables allocated otherwise.
   */
  bool fillMatchesS();

  /**
   * Fills the match set and the bindings manager with all the possible positive and negative bindings between the literals of L and M.
   * If the argument litToRemove is set, it will not consider negative matches for the literals other than the one at index litToRemove.
   *
   * @return false if no subsumption resolution solution is possible, the number of sat variables allocated otherwise.
   */
  void fillMatchesSR(unsigned litToRemove = 0xFFFFFFFF);

  /**
   * Generates the conclusion clause based on the model provided by the sat solver
   *
   * @pre the match set must fill filled
   * @pre the model must be set by the sat solver
   * @return the conclusion clause
   */
  Kernel::Clause *generateConclusion();

public:
  using clock = std::chrono::steady_clock;

  SATSubsumptionAndResolution() : _L(nullptr), _M(nullptr),
                                 _m(0), _n(0)
  { }

  /**
   * Checks whether the instance clause is subsumed by the base clause
   *
   * @param L the base clause (side premise)
   * @param M the instance clause (main premise)
   * @param setSR if true, the Match set will be filled with negative matches as well. It will also check whether subsumption resolution is impossible.
   * @return true if M is subsumed by L
   *
   * @remark The @b setSR parameter is used to save time on the setup of subsumption resolution. If it is intended to check for subsumption resolution right after, it is better to set it to true and call checkSubsumptionResolution with the flag usePreviousMatchSet set to true.
   *
   * A clause L subsumes a clause M if there exists a substitution \f$\sigma\f$ such that \f$\sigma(L)\subseteq M\f$.
   * Where L and M are considered as multisets of literals.
   */
  bool checkSubsumption(Kernel::Clause *L,
                        Kernel::Clause *M,
                        bool setSR = false);

  /**
   * Checks whether a subsumption resolution can occur between the clauses @b L and @b M . If it is possible, returns the conclusion of the resolution, otherwise return NULL.
   *
   * @param L the base clause (side premise)
   * @param M the instance clause (main premise)
   * @param usePreviousMatchSet whether to use the previous match set or not. If false, the match set will be cleared and filled again.
   * @return the conclusion of the resolution, or NULL if no resolution is possible
   *
   * @warning if the @b usePreviousMatchSet flag is true, the last call must have been a call to checkSubsumption with the flag setSR set to true.
   *
   * L /\ M => L /\ M* where M* is the conclusion of the subsumption resolution
   *
   * Subsumption resolution is possible for two clauses (L* V L') and (m' V M*) if there exists a substitution "s" such that
   * s(L') = {~m'} /\ s(L*) \\subseteq M*
   * Where m' is a single literal. The conclusion of the subsumption resolution becomes M*.
   * @example
   * [p(x1, x2) V p(f(x2), x3)] /\ [~p2(f(y1), d) V p2(g(y1), c) V ~p2(f(c), e)]
   * ---------------------------------------------------------------------------
   *                  ~p2(f(y1),d) V p2(g(y1),c)
   * thank's to the substitution S = {x1 -> g(y1), x2 -> c, x3 -> e}
   * @remark Subsumption resolution may not be unique, this method has no guarantee on which solution will be generated.
   * @example
   * [p(x) V q(x) V r(x)] /\ [P(c) \/ Q(c) \/ ~R(c) \/ P(d) \/ ~Q(d) \/ R(d)]
   * may have several conclusion.
   */
  Kernel::Clause *checkSubsumptionResolution(Kernel::Clause *L,
                                             Kernel::Clause *M,
                                             bool usePreviousMatchSet = false);


  /**
   * @brief Checks whether a subsumption resolution can occur between the clauses @b L and @b M with the literal @b resolutionLiteral as the resolution literal.
   *
   * @param L the base clause (side premise)
   * @param M the instance clause (main premise)
   * @param resolutionLiteral the index of the resolution literal in the instance clause @b M
   *
   * @note The conclusion of subsumption resolution is not generated. The conclusion can be later generated using getSubsumptionResolutionConclusion( @b M, @b resolutionLiteral, @b L), provided that @b resolutionLiteral is a valid resolution literal.
   * @note If @b resolutionLiteral is 0xFFFFFFFF, the method will check if subsumption resolution is possible with any literal. But the conclusion will not be generated.
   */
  bool checkSubsumptionResolutionWithLiteral(Kernel::Clause *L,
                                             Kernel::Clause *M,
                                             unsigned resolutionLiteral);

  // return the substitution required for the last subsumption resolution
  // precondition: the last subsumption resolution succeeded
  SimpleSubstitution getBindingsForSubsumptionResolutionWithLiteral();

  /**
   * Creates a clause that is the subsumption resolution of @b M and @b L on @b m_j.
   * L V L' /\ M* V @b m_j => L V L' /\ M*
   *
   * @param M The clause to be subsumed after resolution.
   * @param m_j The literal on which the subsumption resolution is performed.
   * @param L The clause resolving with @b M.
   */
  static Kernel::Clause *getSubsumptionResolutionConclusion(Kernel::Clause *M,
                                                            Kernel::Literal *m_j,
                                                            Kernel::Clause *L);
  };

} // namespace SATSubsumption

#endif /* !SAT_SUBSUMPTION_RESOLUTION_HPP */
