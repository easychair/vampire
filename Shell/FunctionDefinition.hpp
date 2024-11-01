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
 * @file FunctionDefinition.hpp
 * Defines class FunctionDefinition implementing work with definitions.
 *
 * @since 28/05/2004 Manchester
 * @note There is a lot of confusion in the definition handling
 *       mechanism of Vampire. There is much in common between
 *       this class and class Definition and they should be reworked.
 *       In addition, Preprocess eliminates unused definitions
 *       only from FO (not clausal) problems. Why? Options do
 *       not separate function and predicate definitions. This
 *       should be re-thought and re-implemented but it is too
 *       close to CASC now to make big changes.
 */

#ifndef __FunctionDefinition__
#define __FunctionDefinition__

#include "Forwards.hpp"

#include "Lib/DHMap.hpp"
#include "Lib/MultiCounter.hpp"
#include "Lib/ProofExtra.hpp"
#include "Lib/Stack.hpp"
#include "Kernel/Unit.hpp"

namespace Shell {

using namespace Lib;
using namespace Kernel;

/**
 * Class implementing procedures related to function
 * definitions. A function definition is of the form
 * f(x1,...,xn) = t, where
 * <ol>
 *  <li>all variables of t occur in {x1,...,xn}; </li>
 *  <li>f does not occur in t.</li>
 *  <li>if f is a constant, then t is a constant, too.
 * </ol>
 * @since 29/05/2004 Manchester
 */
class FunctionDefinition
{
public:
  struct Def;
  FunctionDefinition();
  ~FunctionDefinition();

  void removeAllDefinitions(Problem& prb, bool inHigherOrder);
  bool removeAllDefinitions(UnitList*& units, bool inHigherOrder);

  static void removeUnusedDefinitions(Problem& prb, bool inHigherOrder);
  static bool removeUnusedDefinitions(UnitList*& units, Problem* prb, bool inHigherOrder);


  static Def* isFunctionDefinition (Unit&, bool inHigherOrder);
  static void deleteDef(Def* def);

//   int removeAllDefinitions ();
//   static bool isFunctionDefinition (const Clause&, Term*& lhs, Term*& rhs);
//   static bool isFunctionDefinition (const Literal*, Term*& lhs, Term*& rhs);
//   static bool isFunctionDefinition (const Formula&, Term*& lhs, Term*& rhs);

private:
  static Def* isFunctionDefinition (Clause*, bool inHigherOrder);
  static Def* isFunctionDefinition (FormulaUnit&, bool inHigherOrder);
  static Def* isFunctionDefinition (Literal*, bool inHigherOrder);
  static Def* defines (Term* lhs, Term* rhs, bool inHigherOrder);
  static bool occurs (unsigned function, Term&);
  static void reverse(Def*);

  bool isDefined(Term* t);

  Term* applyDefinitions(Literal* t, Stack<Def*>* usedDefs);
  Clause* applyDefinitions(Clause* cl);

  void checkDefinitions(Def* t);
  void assignArgOccursData(Def* d);

//   void count (const Clause& c);
//   bool unfoldAllDefs ();
//   bool unfold (Def& def);
//   bool unfold (Term* t,UnitList* parents);
//   bool unfoldList (Term* ts,UnitList& parents);

//   void apply (LiteralList& ls,UnitList& parents);
//   void apply (Literal& l,UnitList& parents);
//   void apply (TermList& ls,UnitList& parents);
//   void apply (Term& l,UnitList& parents);

  typedef DHMap<int, Def*, IdentityHash, DefaultHash> Fn2DefMap;
  Fn2DefMap _defs;

  /** stack where definitions are put when they're marked as blocked */
  Stack<Def*> _blockedDefs;

  Stack<Def*> _safeDefs;
  /** Counters for occurrences of function symbols */
  MultiCounter _counter;
  /** The number of found definitions */
  int _found;
  /** The number of removed definitions */
  int _removed;

  Problem* _processedProblem;
}; // class FunctionDefinition

struct FunctionDefinitionExtra : public InferenceExtra {
  std::vector<Term *> lhs;
  FunctionDefinitionExtra(std::vector<Term *> lhs) : lhs(lhs) {}
  void output(std::ostream &out) const override;
};

}

#endif
