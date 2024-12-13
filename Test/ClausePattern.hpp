/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __TEST__CLAUSE_PATTERN__
#define __TEST__CLAUSE_PATTERN__

#include "Kernel/Clause.hpp"
#include "Lib/Coproduct.hpp"
#include "Test/TestUtils.hpp"

namespace Test {

class ClausePattern;

/**
 * An alternative between two patterns. matches if either of the subpatterns matches.
 */
struct AnyOf 
{
  std::shared_ptr<ClausePattern> lhs;
  std::shared_ptr<ClausePattern> rhs;
};

/**
 * A pattern to match a clause against. 
 * Can be either a Clause, or a AnyOf which is a combination of two patterns.
 * A Clause matches a pattern Clause, if they are equal.
 * A Clause matches an AnyOf pattern if it matches both of the subpatterns.
 */
class ClausePattern : Coproduct<Kernel::Clause*, AnyOf>
{
  using Copro =  Coproduct<Kernel::Clause*, AnyOf>;
public:
  ClausePattern(Kernel::Clause* clause) 
    : Copro(clause) {}

  ClausePattern(ClausePattern l, ClausePattern r) : Copro(AnyOf {
        std::make_unique<ClausePattern>(std::move(l)),
        std::make_unique<ClausePattern>(std::move(r))
      }) {}

  template<class EqualityOperator>
  bool matches(EqualityOperator& equality, Kernel::Clause* result);
  friend std::ostream& operator<<(std::ostream& out, ClausePattern const& self);
};

inline std::ostream& operator<<(std::ostream& out, ClausePattern const& self) 
{
  return self.match(
      [&](Kernel::Clause* const& self) -> std::ostream&
      { return out << pretty(self); },

      [&](AnyOf const& self)  -> std::ostream&
      { return out << pretty(self.lhs) << " or " << pretty(self.rhs); });
}

template<class EqualityOperator>
bool ClausePattern::matches(EqualityOperator& equality, Kernel::Clause* result)
{
  return match(
      [&](Kernel::Clause* self) 
      { return equality.eq(result, self); },

      [&](AnyOf& self) 
      { return self.lhs->matches(equality, result) || self.rhs->matches(equality, result); });
}

inline ClausePattern anyOf(Kernel::Clause* lhs) 
{ return ClausePattern(lhs); }

template<class... As>
inline ClausePattern anyOf(Kernel::Clause* lhs, Kernel::Clause* rhs, As... rest) 
{ return ClausePattern(lhs, anyOf(rhs, rest...)); }


};

#endif // __TEST__CLAUSE_PATTERN__
