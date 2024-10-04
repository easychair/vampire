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
 * @file ProofExtra.hpp
 * Various objects that include "extra information" about an inference,
 * e.g. selected literals, unifier, etc.
 *
 * @since 09/09/2024 Oxford
 */

#ifndef __Inferences_ProofExtra__
#define __Inferences_ProofExtra__

#include "Kernel/Term.hpp"
#include "Lib/ProofExtra.hpp"

namespace Inferences {

// inferences that use one literal from their main premise
struct LiteralInferenceExtra : public InferenceExtra {
  LiteralInferenceExtra(Kernel::Literal *selected) : selectedLiteral(selected) {}

  virtual void output(std::ostream &out) const override;

  // the literal from the main premise
  Kernel::Literal *selectedLiteral;
};

// inferences that use one literal from their side premise
struct TwoLiteralInferenceExtra : public InferenceExtra {
  TwoLiteralInferenceExtra(Kernel::Literal *selected, Kernel::Literal *other)
    : selectedLiteral(selected), otherLiteral(other) {}

  virtual void output(std::ostream &out) const override;

  // selected literal
  LiteralInferenceExtra selectedLiteral;
  // the literal from the side premise
  Kernel::Literal *otherLiteral;
};

struct RewriteInferenceExtra : public InferenceExtra {
  RewriteInferenceExtra(Kernel::TermList lhs, Kernel::TermList target)
    : lhs(lhs), rewritten(target) {}

  virtual void output(std::ostream &out) const override;

  // the LHS used to rewrite with
  Kernel::TermList lhs;
  // the rewritten term
  Kernel::TermList rewritten;
};

struct TwoLiteralRewriteInferenceExtra : public InferenceExtra {
  TwoLiteralRewriteInferenceExtra(
    Kernel::Literal *selected,
    Kernel::Literal *other,
    Kernel::TermList lhs,
    Kernel::TermList rewritten)
    : selected(selected, other), rewrite(lhs, rewritten) {}
  virtual void output(std::ostream &out) const override;

  // selected literals
  TwoLiteralInferenceExtra selected;
  // rewrite information
  RewriteInferenceExtra rewrite;
};
}

#endif
