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
 * @file Unit.cpp
 * Defines class Unit for all kinds of proof units
 *
 * @since 09/05/2007 Manchester
 */

#include "Forwards.hpp"


#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/List.hpp"
#include "Lib/Set.hpp"
#include "Lib/DHMap.hpp"

#include "Shell/Statistics.hpp"

#include "Inference.hpp"
#include "InferenceStore.hpp"
#include "Clause.hpp"
#include "Formula.hpp"
#include "FormulaUnit.hpp"
#include "SubformulaIterator.hpp"

#include "Unit.hpp"

using namespace std;
using namespace Kernel;

unsigned Unit::_lastNumber = 0;
unsigned Unit::_firstNonPreprocessingNumber = 0;
unsigned Unit::_lastParsingNumber = 0;

/**
 * Should be called after the preprocessing and before the start
 * of the saturation algorithm.
 */
void Unit::onPreprocessingEnd()
{
  ASS(!_firstNonPreprocessingNumber);

  _firstNonPreprocessingNumber=_lastNumber+1;
}

/** New unit of a given kind */
Unit::Unit(Kind kind, Inference inf)
  : _number(++_lastNumber),
    _kind(kind),
    _inheritedColor(COLOR_INVALID),
    _inference(std::move(inf))
{
} // Unit::Unit

void Unit::incRefCnt()
{
  if(isClause()) {
    static_cast<Clause*>(this)->incRefCnt();
  }
}

void Unit::decRefCnt()
{
  if(isClause()) {
    static_cast<Clause*>(this)->decRefCnt();
  }
}

Clause* Unit::asClause() {
  ASS(isClause());
  return static_cast<Clause*>(this);
}


Color Unit::getColor()
{
  if(isClause()) {
    return static_cast<Clause*>(this)->color();
  }
  else {
    return static_cast<FormulaUnit*>(this)->getColor();
  }
}

unsigned Unit::getWeight()
{
  if(isClause()) {
    return static_cast<Clause*>(this)->weight();
  }
  else {
    return static_cast<FormulaUnit*>(this)->weight();
  }
}

void Unit::destroy()
{
  if(isClause()) {
    static_cast<Clause*>(this)->destroy();
  }
  else {
    static_cast<FormulaUnit*>(this)->destroy();
  }
}

std::string Unit::toString() const
{
  if(isClause()) {
    return static_cast<const Clause*>(this)->toString();
  }
  else {
    return static_cast<const FormulaUnit*>(this)->toString();
  }
}

unsigned Unit::varCnt()
{
  if(isClause()) {
    return static_cast<Clause*>(this)->varCnt();
  }
  else {
    return static_cast<FormulaUnit*>(this)->varCnt();
  }
}


/**
 * Return quantified formula equivalent to the unit.
 *
 * @since 16/01/14, removed BDDNode prop, Giles.
 */
Formula* Unit::getFormula()
{
  if(isClause()) {
    return Formula::fromClause(static_cast<Clause*>(this));//, prop);
  }
  else {
    return Formula::quantify(static_cast<FormulaUnit*>(this)->formula());
  }
}

/**
 * Print the inference as a std::string (used in printing units in
 * refutations).
 * @since 04/01/2008 Torrevieja
 */
std::string Unit::inferenceAsString() const
{
  const Inference& inf = inference();

  std::string result = (std::string)"[" + inf.name();
  bool first = true;

  auto it = inf.iterator();
  while (inf.hasNext(it)) {
    Unit* parent = inf.next(it);
    result += first ? ' ' : ',';
    first = false;
    result += Int::toString(parent->number());
  }

  // print extra if present
  if(env.options->proofExtra() == Options::ProofExtra::FULL) {
    auto *extra = env.proofExtra.find(this);
    if(extra) {
      if(!first)
        result += ',';
      result += extra->toString();
    }
  }

  return result + ']';
} // Unit::inferenceAsString()

void Unit::assertValid()
{
  if(isClause()) {
    ASS_ALLOC_TYPE(this,"Clause");
  }
  else {
    ASS_ALLOC_TYPE(this,"FormulaUnit");
  }
}

UnitIterator Unit::getParents() const
{
  // The unit itself stores the inference
  UnitList* res = 0;
  Inference::Iterator iit = _inference.iterator();
  while(_inference.hasNext(iit)) {
    Unit* premUnit = _inference.next(iit);
    UnitList::push(premUnit, res);
  }
  res = UnitList::reverse(res); // we want items in the same order
  return pvi(UnitList::DestructiveIterator(res));
}

bool Unit::minimizeAncestorsAndUpdateSelectedStats()
{
  Stack<std::pair<Unit*,bool>> todo;
  DHSet<Unit*> done;
  bool seenInputInference = false;

  todo.push(make_pair(this,false));
  while(!todo.isEmpty()) {
    Unit* current;
    bool collecting;
    std::tie(current,collecting) = todo.pop();
    if (collecting) {
      Inference& inf = current->inference();
      seenInputInference |= (inf.rule() == InferenceRule::INPUT);
      Inference::Iterator iit = inf.iterator();
      if (inf.hasNext(iit)) {
        UnitInputType uit = UnitInputType::AXIOM; // we compute a maximum, so start from the smallest element
        bool isPureTheoryDescendant = true; // see also Inference::initMany
        while(inf.hasNext(iit)) {
          Unit* premUnit = inf.next(iit);
          uit = getInputType(uit,premUnit->inputType());
          isPureTheoryDescendant &= premUnit->isPureTheoryDescendant();
        }
        current->setInputType(uit);
        inf.setPureTheoryDescendant(isPureTheoryDescendant);
      } else if (inf.rule() == InferenceRule::AVATAR_DEFINITION) {
        // don't touch _pureTheoryDescendant for AVATAR_DEFINITIONs (a split theory consequence is again a theory consequence)
        current->setInputType(UnitInputType::AXIOM); // AVATAR_DEFINITION might have inherited goaledness from its causal parent, which we want to reset
      } else {
        // no premises and not InferenceRule::AVATAR_DEFINITION
        // we leave input type unchanged and isPureTheoryDescendant is already correctly set to isTheoryAxiom
        ASS_EQ(inf.isPureTheoryDescendant(),inf.isTheoryAxiom());
      }
      inf.updateStatistics(); // in particular, update inductionDepth (which could have decreased, since we might have fewer parents after miniminization)

      switch (inf.rule()) {
        case InferenceRule::STRUCT_INDUCTION_AXIOM_ONE:
        case InferenceRule::STRUCT_INDUCTION_AXIOM_TWO:
        case InferenceRule::STRUCT_INDUCTION_AXIOM_THREE:
        case InferenceRule::STRUCT_INDUCTION_AXIOM_RECURSION:
          env.statistics->structInductionInProof++;
          break;
        case InferenceRule::INT_INF_UP_INDUCTION_AXIOM:
        case InferenceRule::INT_INF_DOWN_INDUCTION_AXIOM:
          env.statistics->intInfInductionInProof++;
          break;
        case InferenceRule::INT_FIN_UP_INDUCTION_AXIOM:
        case InferenceRule::INT_FIN_DOWN_INDUCTION_AXIOM:
          env.statistics->intFinInductionInProof++;
          break;
        case InferenceRule::INT_DB_UP_INDUCTION_AXIOM:
        case InferenceRule::INT_DB_DOWN_INDUCTION_AXIOM:
          env.statistics->intDBInductionInProof++;
          break;
        default:
          ;
      }
      switch (inf.rule()) {
        case InferenceRule::INT_INF_UP_INDUCTION_AXIOM:
          env.statistics->intInfUpInductionInProof++;
          break;
        case InferenceRule::INT_INF_DOWN_INDUCTION_AXIOM:
          env.statistics->intInfDownInductionInProof++;
          break;
        case InferenceRule::INT_FIN_UP_INDUCTION_AXIOM:
          env.statistics->intFinUpInductionInProof++;
          break;
        case InferenceRule::INT_FIN_DOWN_INDUCTION_AXIOM:
          env.statistics->intFinDownInductionInProof++;
          break;
        case InferenceRule::INT_DB_UP_INDUCTION_AXIOM:
          env.statistics->intDBUpInductionInProof++;
          break;
        case InferenceRule::INT_DB_DOWN_INDUCTION_AXIOM:
          env.statistics->intDBDownInductionInProof++;
          break;
        default:
          ;
      }

    } else {
      if (!done.insert(current)) {
        continue;
      }
      todo.push(make_pair(current,true)); // to collect stuff when children done
      Inference& inf = current->inference();
      inf.minimizePremises(); // here we do the minimization
      Inference::Iterator iit = inf.iterator();
      while(inf.hasNext(iit)) {
        Unit* premUnit = inf.next(iit);
        todo.push(make_pair(premUnit,false));
      }
    }
  }
  return seenInputInference;
}

std::ostream& Kernel::operator<<(std::ostream& out, const Unit& u)
{
  return out << u.toString();
}
