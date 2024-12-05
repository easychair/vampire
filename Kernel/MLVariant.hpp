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
 * @file MLVariant.hpp
 * Defines class MLVariant with methods
 * supporting multiliteral variant testing.
 */


#ifndef __MLVariant__
#define __MLVariant__

#include "Forwards.hpp"
#include "SubstHelper.hpp"

namespace Kernel {

using namespace Lib;

class MLVariant {
public:
  static bool isVariant(Literal* const * cl1Lits, Clause* cl2, LiteralList** alts, SimpleSubstitution *map = nullptr);
  static bool isVariant(Literal* const * cl1Lits, Clause* cl2, bool complementary=false, SimpleSubstitution *map = nullptr);
  static bool isVariant(Clause* cl1, Clause* cl2, bool complementary=false)
  {
    return cl1->length()==cl2->length() && isVariant(cl1->literals(), cl2, complementary);
  }

};


};

#endif /* __MLVariant__ */
