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
 * @file Ordering.hpp
 * Defines (abstract) class Ordering for simplification orderings
 *
 * @since 30/04/2008 flight Brussels-Tel Aviv
 */

#ifndef __Ordering__
#define __Ordering__

#include "Forwards.hpp"

#include "Debug/Assertion.hpp"

#include "Lib/Comparison.hpp"
#include "Lib/SmartPtr.hpp"
#include "Lib/DArray.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/Portability.hpp"
#include "Kernel/SubstHelper.hpp"

#define DEBUG_ORDERING 0

namespace Kernel {

using namespace Shell;

/**
 * An abstract class for simplification orderings
 * @since 30/04/2008 flight Brussels-Tel Aviv
 */
class Ordering
{
public:
  /**
   * Represents the results of ordering comparisons
   *
   * Values of elements must be equal to values of corresponding elements
   * in the @c ArgumentOrderVals enum, so that one can convert between the
   * enums using static_cast.
   */
  enum [[nodiscard]] Result {
    GREATER=1,
    LESS=2,
    EQUAL=3,
    INCOMPARABLE=4
  };

  struct Constraint {
    TermList lhs;
    TermList rhs;
    Result rel;

    friend std::ostream& operator<<(std::ostream& out, const Constraint& self)
    { return out << self.lhs << " " << self.rhs << " " << resultToString(self.rel); }
  };

  virtual ~Ordering() = default;

  /** Return the result of comparing @b l1 and @b l2 */
  virtual Result compare(Literal* l1,Literal* l2) const = 0;

  /** Return the result of comparing terms (not term lists!)
   * @b t1 and @b t2 */
  virtual Result compare(TermList t1,TermList t2) const = 0;

  /** Same as @b compare, for applied (substituted) terms. */
  virtual Result compare(AppliedTerm lhs, AppliedTerm rhs) const
  { return compare(lhs.apply(), rhs.apply()); }

  /** Optimised function used for checking that @b t1 is greater than @b t2,
   * under some substitutions captured by @b AppliedTerm. */
  virtual Result isGreaterOrEq(AppliedTerm t1, AppliedTerm t2) const
  { return compare(t1, t2); }

  /** Creates optimised object for ordering checks. @see OrderingComparator. */
  virtual OrderingComparatorUP createComparator() const;

  virtual void show(std::ostream& out) const = 0;

  static bool isGreaterOrEqual(Result r) { return (r == GREATER || r == EQUAL); }

  void removeNonMaximal(LiteralList*& lits) const;

  static Result fromComparison(Comparison c);
  static Comparison intoComparison(Result c);

  static Result reverse(Result r)
  {
    switch(r) {
    case GREATER:
      return LESS;
    case LESS:
      return GREATER;
    case EQUAL:
    case INCOMPARABLE:
      return r;
    default:
      ASSERTION_VIOLATION;
    }
  }
  static const char* resultToString(Result r);

  static Ordering* create(Problem& prb, const Options& opt);

  static bool trySetGlobalOrdering(OrderingSP ordering);
  static Ordering* tryGetGlobalOrdering();

  Result getEqualityArgumentOrder(Literal* eq) const;
protected:

  Result compareEqualities(Literal* eq1, Literal* eq2) const;

private:
  /**
   * Helper methods for comparing literals s1=s2 and t1=t2.
   */
  Result compare_s1Gt1(TermList s1,TermList s2,TermList t1,TermList t2) const;
  Result compare_s1It1(TermList s1,TermList s2,TermList t1,TermList t2) const;
  Result compare_s1It1_s2It2(TermList s1,TermList s2,TermList t1,TermList t2) const;
  Result compare_s1Gt1_s2It2(TermList s1,TermList s2,TermList t1,TermList t2) const;
  Result compare_s1Gt1_s2Lt2(TermList s1,TermList s2,TermList t1,TermList t2) const;
  Result compare_s1Gt1_s1It2_s2It1(TermList s1,TermList s2,TermList t1,TermList t2) const;

  /**
   * We store orientation of equalities in this ordering inside
   * the term sharing structure. Setting an ordering to be global
   * does not change the behavior of Vampire, but may lead to
   * better performance, as the equality orientation will be cached
   * inside the sharing structure.
   */
  static OrderingSP s_globalOrdering;
}; // class Ordering

// orderings that rely on symbol precedence
class PrecedenceOrdering
: public Ordering
{
public:
  Result compare(Literal* l1, Literal* l2) const override;
  void show(std::ostream&) const override;
  virtual void showConcrete(std::ostream&) const = 0;

protected:
  // l1 and l2 are not equalities and have the same predicate
  virtual Result comparePredicates(Literal* l1,Literal* l2) const = 0;

  PrecedenceOrdering(const DArray<int>& funcPrec, const DArray<int>& typeConPrec,
                     const DArray<int>& predPrec, const DArray<int>& predLevels, bool reverseLCM);
  PrecedenceOrdering(Problem& prb, const Options& opt, const DArray<int>& predPrec);
  PrecedenceOrdering(Problem& prb, const Options& opt);


  static DArray<int> typeConPrecFromOpts(Problem& prb, const Options& opt);
  static DArray<int> funcPrecFromOpts(Problem& prb, const Options& opt);
  static DArray<int> predPrecFromOpts(Problem& prb, const Options& opt);
  static DArray<int> predLevelsFromOptsAndPrec(Problem& prb, const Options& opt, const DArray<int>& predicatePrecedences);

  Result compareFunctionPrecedences(unsigned fun1, unsigned fun2) const;
  Result compareTypeConPrecedences(unsigned tyc1, unsigned tyc2) const;

  int predicatePrecedence(unsigned pred) const;
  int predicateLevel(unsigned pred) const;

  /** number of predicates in the signature at the time the order was created */
  unsigned _predicates;
  /** number of functions in the signature at the time the order was created */
  unsigned _functions;
  /** Array of predicate levels */
  DArray<int> _predicateLevels;
  /** Array of predicate precedences */
  DArray<int> _predicatePrecedences;
  /** Array of function precedences */
  DArray<int> _functionPrecedences;
  /** Array of type con precedences */
  DArray<int> _typeConPrecedences;

  bool _reverseLCM;
};


inline std::ostream& operator<<(std::ostream& out, Ordering::Result const& r)
{
  switch (r) {
    case Ordering::Result::GREATER: return out << "GREATER";
    case Ordering::Result::LESS: return out << "LESS";
    case Ordering::Result::EQUAL: return out << "EQUAL";
    case Ordering::Result::INCOMPARABLE: return out << "INCOMPARABLE";
    default:
      return out << "UNKNOWN";
  }
  ASSERTION_VIOLATION
}

}

#endif
