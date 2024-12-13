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
 * @file Unit.hpp
 * Defines class Unit for all kinds of proof units
 *
 * @since 08/05/2007 Manchester
 */

#ifndef __Unit__
#define __Unit__

#include "Forwards.hpp"

#include "Lib/List.hpp"
#include "Lib/Output.hpp"
#include "Kernel/Inference.hpp"

namespace Kernel {

using namespace Lib;

/**
 * Class to represent units of inference (such as clauses and formulas).
 * @since 08/05/2007 Manchester
 */
class Unit
{
protected:
  ~Unit() {}
public:
  /** Kind of unit. The integers should not be changed, they are used in
   *  Compare.
   */
  enum Kind {
    /** Clause unit */
    CLAUSE = 0,
    /** Formula unit */
    FORMULA = 1
  };

  void destroy();
  std::string toString() const;
  unsigned varCnt();

  std::string inferenceAsString() const;

  /** True if a clause unit */
  bool isClause() const
  { return _kind == CLAUSE; }

  Clause* asClause();

  /** Return the number of this unit */
  unsigned number() const { return _number; }
  /** Forcefully change the unit's number - use with care! - numbers should be unique across the whole board! */
  void overwriteNumber(unsigned newNumber) { _number = newNumber; }

  /** Return the inference of this unit */
  Inference& inference() { return _inference; }
  const Inference& inference() const { return _inference; }

  /** return the input type of the unit */
  UnitInputType inputType() const { return _inference.inputType(); }
  /** set the input type of the unit */
  void setInputType(UnitInputType it) { _inference.setInputType(it); }
  /** return true if inputType relates to a goal **/
  bool derivedFromGoal() const { return _inference.derivedFromGoal(); }
  /** see isPureTheoryDescendant in Inference.cpp */
  bool isPureTheoryDescendant() const { return _inference.isPureTheoryDescendant(); }
  /** see isCombAxiomsDescendant in Inference.cpp */
  bool isCombAxiomsDescendant() const { return _inference.isCombAxiomsDescendant(); }
  /** see isProxyAxiomsDescendant in Inference.cpp */
  bool isProxyAxiomsDescendant() const { return _inference.isProxyAxiomsDescendant(); }
  /** see isHolAxiomsDescendant in Inference.cpp */
  bool isHolAxiomsDescendant() const { return _inference.isHolAxiomsDescendant(); }
  /** see isTheoryAxiom in Inference.cpp */
  bool isTheoryAxiom() const { return _inference.isTheoryAxiom(); }

  unsigned char getSineLevel() const { return _inference.getSineLevel(); }
  /** true if the unit is read from a TPTP included file  */
  bool included() const { return _inference.included(); }

  /** This is just a more convenient (but also less efficient) way
   * for accessing the unit's parents than using its inferefence's iterator. */
  UnitIterator getParents() const;

  /**
   * Recursive minimization of ancestors for sat-based inferences (AVATAR_REFUTATION, GLOBAL_SUMSUPTION) triggered here.
   *
   * (Note that with minimizeSatProofs set to off, it's already too late to do this as we did not collect the necessary info,
   * however, it pays off to call this anyway, to correctly set any required stats in the inference objects.)
   *
   * Could be extended further if more stats provided by the inference tree should be found
   * to require some final touches before proof printing. Currently, we care about the InputType, the induction stats, and _isPureTheoryDescendant.
   *
   * returns true if the minimized tree traversal observed an INPUT inference (to facilitate a certain sanity check later).
  */
  bool minimizeAncestorsAndUpdateSelectedStats();

  /** Return the inherited color of the unit or COLOR_INVALID
   * if there isn't an inherited color.
   *
   * Inherited color is set by parser constructs left_formula,
   * right_formula and end_formula.
   */
  inline Color inheritedColor() const
  {
    return static_cast<Color>(_inheritedColor);
  }
  void setInheritedColor(Color color)
  {
    ASS_NEQ(color,COLOR_INVALID);
    _inheritedColor = color;
  } // setInheritedColor

  void invalidateInheritedColor() { _inheritedColor = COLOR_INVALID; }

  Color getColor();
  unsigned getWeight();

  Formula* getFormula();

  /**
   * Increase the number of references to the unit.
   *
   * Only clauses are reference-counted, for FormulaUnits nothing
   * is done.
   */
  void incRefCnt();
  /**
   * Decrease the number of references to the unit.
   *
   * Only clauses are reference-counted, for FormulaUnits nothing
   * is done.
   */
  void decRefCnt();

  /** Return true iff unit was created during preprocessing
   * (and not during the run of the saturation algorithm) */
  inline bool isFromPreprocessing()
  { return !_firstNonPreprocessingNumber || _number<_firstNonPreprocessingNumber; }

  void assertValid();

  static void onPreprocessingEnd();
  static void onParsingEnd(){ _lastParsingNumber = _lastNumber;}
  static unsigned getLastParsingNumber(){ return _lastParsingNumber;}

protected:
  /** Number of this unit, used for printing and statistics */
  unsigned _number;
  /** Kind  */
  unsigned _kind : 1;

  /** used in interpolation to denote parents of what color have been used */
  unsigned _inheritedColor : 2;

  /** inference used to obtain the unit */
  Inference _inference;

  Unit(Kind kind, Inference inf);

  /** Used to enumerate units */
  static unsigned _lastNumber;

  /** Used to determine which clauses come from preprocessing
   *
   * 0 means preprocessing is not over yet.*/
  static unsigned _firstNonPreprocessingNumber;

  static unsigned _lastParsingNumber;

private:
  template<class Pre, class Post>
  void _traverseParents(Unit& self, unsigned& depth, Pre& pre, Post& post) {
    pre(unsigned(depth), &self);
    depth++;
    auto infit = self.inference().iterator();
    while (self.inference().hasNext(infit)) {
      auto parent = self.inference().next(infit);
      _traverseParents(*parent, depth, pre, post);
    }
    depth--;
    post(unsigned(depth), &self);
  }


protected:
  template<class Pre, class Post>
  void traverseParents(Pre pre, Post post) {
    unsigned depth = 0;
    _traverseParents(*this, depth, pre, post);
  }

  template<class Pre>
  void traverseParentsPre (Pre pre ) { traverseParents(pre , [](auto...) {}); }

  template<class Post>
  void traverseParentsPost(Post post) { traverseParents([](auto...) {}, post); }

}; // class Unit

std::ostream& operator<< (std::ostream& out, const Unit& u );

}
#endif
