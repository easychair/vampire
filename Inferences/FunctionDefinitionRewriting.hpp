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
 * @file FunctionDefinitionRewriting.hpp
 * Defines class FunctionDefinitionRewriting.
 */

#ifndef __FunctionDefinitionRewriting__
#define __FunctionDefinitionRewriting__

#include "Forwards.hpp"

namespace Inferences {

using namespace Kernel;
using namespace Shell;

/**
 * Inference implementing function definition rewriting.
 * Function definitions are assumed to be available when saturation begins,
 * so there is only a forward version of the inference. Moreover, we use
 * a forward simplification variant to eagerly perform rewritings which
 * are also demodulations.
 */
class FunctionDefinitionRewriting
  : public GeneratingInferenceEngine
  , public ForwardSimplificationEngine
  {
public:
  USE_ALLOCATOR(FunctionDefinitionRewriting);

  ClauseIterator generateClauses(Clause *premise) override;
  bool perform(Clause* cl, Clause*& replacement, ClauseIterator& premises) override;

private:
  static Clause *perform(
      Clause *rwClause, Literal *rwLiteral, TermList rwTerm,
      Clause *eqClause, Literal *eqLiteral, TermList eqLHS,
      ResultSubstitutionSP subst, bool toplevelCheck,
      bool& isEqTautology, const Inference& inf, SaturationAlgorithm* salg = nullptr);
};

}; // namespace Inferences

#endif /* __FunctionDefinitionRewriting__ */