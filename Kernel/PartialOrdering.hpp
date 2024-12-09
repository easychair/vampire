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
 * @file PartialOrdering.hpp
 * Defines class PartialOrdering.
 */

#ifndef __PartialOrdering__
#define __PartialOrdering__

#include <string>
#include <cstdint>

namespace Kernel {

enum class PoComp : uint8_t {
  UNKNOWN,
  GREATER,
  EQUAL,
  LESS,
  NGEQ,
  NLEQ,
  INCOMPARABLE,
};

bool checkCompatibility(PoComp old, PoComp curr, PoComp& res);
std::string po_to_infix(PoComp c);

class PartialOrdering
{
public:
  PoComp get(size_t x, size_t y) const;

  static const PartialOrdering* getEmpty();
  static const PartialOrdering* set(const PartialOrdering* po, size_t x, size_t y, PoComp v);
  static const PartialOrdering* extend(const PartialOrdering* po);

  // Returns if PO contains full incomparability yet.
  // Useful to discard branches when reasoning over ground terms.
  bool hasIncomp() const { return _hasIncomp; }

  std::string to_string() const;
  std::string all_to_string() const;

private:
  PartialOrdering();
  PartialOrdering(const PartialOrdering& other);
  ~PartialOrdering();
  PartialOrdering& operator=(const PartialOrdering&) = delete;

  void extend();

  PoComp get_unsafe(size_t x, size_t y) const;
  bool set_idx_of(size_t x, size_t y, PoComp v, bool& changed);
  bool set_idx_of_safe(size_t x, size_t y, PoComp v, bool& changed);

  bool set_inferred(size_t x, size_t y, PoComp result);
  bool set_inferred_loop(size_t x, size_t y, PoComp rel);
  bool set_inferred_loop_inc(size_t x, size_t y, PoComp wkn);
  bool set_inferred_loop_eq(size_t x, size_t y);

  size_t _size;
  PoComp* _array;
  bool _hasIncomp;
  const PartialOrdering* _prev;
};

};

#endif /* __PartialOrdering__ */