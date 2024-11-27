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
 * @file SMTCheck.hpp
 * Routines for producing a ground SMT proof-check script
 */

#include <ostream>

#include "Kernel/Unit.hpp"
#include "Shell/Dedukti.hpp"

namespace Shell {
namespace SMTCheck {

void outputTypeDecl(std::ostream &out, const std::string &name, Kernel::OperatorType *type);
void outputStep(std::ostream &out, Kernel::Unit *u);

}
}
