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
void outputDeduction(std::ostream &out, Kernel::Unit *axiom);

}
}

#endif
