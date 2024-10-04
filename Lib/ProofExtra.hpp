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
 *
 * Include "extra information" about an inference, e.g. selected literals, unifier, etc.
 *
 * @since 09/09/2024 Oxford
 */

#ifndef __Lib_ProofExtra__
#define __Lib_ProofExtra__

#include "Lib/DHMap.hpp"

namespace Lib {

// base for extra inference information
struct InferenceExtra {
  virtual ~InferenceExtra() {};

  // for debug printing only, specific output formats should downcast and handle individually
  virtual void output(std::ostream &out) const = 0;

  std::string toString() const {
    std::stringstream ss;
    this->output(ss);
    return ss.str();
  }
};

inline std::ostream &operator<<(std::ostream &out, const InferenceExtra &extra) {
  extra.output(out);
  return out;
}

// container for extras
class ProofExtra {
  DHMap<Kernel::Unit *, std::unique_ptr<InferenceExtra>> extras;

public:
  // associate `extra` with `unit`, taking ownership of `extra`
  void insert(Kernel::Unit *unit, InferenceExtra *extra) {
    extras.insert(unit, std::unique_ptr<InferenceExtra>(extra));
  }

  // remove the extra information for this unit
  void remove(Kernel::Unit *unit) {
    extras.remove(unit);
  }

  // associated InferenceExtra if present, nullptr otherwise
  const InferenceExtra *find(const Kernel::Unit *unit) {
    auto *found = extras.findPtr(const_cast<Kernel::Unit *>(unit));
    if(!found)
      return nullptr;
    return found->get();
  }

  // associated InferenceExtra: must be present
  template<typename T>
  const T &get(const Kernel::Unit *unit) {
    auto &found = extras.get(const_cast<Kernel::Unit *>(unit));
    ASS(found)
    return static_cast<const T &>(*found);
  }
};
}

#endif
