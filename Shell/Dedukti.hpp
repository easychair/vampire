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
 * @file Dedukti.hpp
 * Routines for Dedukti output
 */

#ifndef __Dedukti__
#define __Dedukti__

#include <ostream>

#include "Kernel/InferenceStore.hpp"
#include "Kernel/OperatorType.hpp"

namespace Shell {
namespace Dedukti {

void outputPrelude(std::ostream &out);
void outputTypeDecl(std::ostream &out, const char *name, Kernel::OperatorType *type);
void outputAxiom(std::ostream &out, Kernel::Unit *axiom);

struct Inference {};

struct BinaryResolutionInference: public Inference {
  CLASS_NAME(ResolutionProofInfo)
  USE_ALLOCATOR(BinaryResolutionInference)
  BinaryResolutionInference(unsigned left, unsigned right) : leftIndex(left), rightIndex(right) {}
  unsigned leftIndex, rightIndex;
};

void registerInference(Kernel::Unit *unit, Inference *info);

struct ProofPrinter: public InferenceStore::ProofPrinter {
  ProofPrinter(std::ostream &out, InferenceStore *store) : InferenceStore::ProofPrinter(out, store) {}
  void printStep(Unit* cs) override;
};

}
}

#endif
