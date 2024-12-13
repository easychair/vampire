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
 * @file Clause.cpp
 * Implements class Clause for units consisting of clauses
 *
 * @since 18/05/2007 Manchester
 */

#include <ostream>

#include "Debug/RuntimeStatistics.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/DArray.hpp"
#include "Lib/Output.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/SharedSet.hpp"
#include "Lib/Stack.hpp"
#include "Lib/BitUtils.hpp"

#include "Saturation/ClauseContainer.hpp"
#include "Saturation/Splitter.hpp"

#include "SAT/SATClause.hpp"

#include "Shell/ConditionalRedundancyHandler.hpp"
#include "Shell/Options.hpp"

#include "Inference.hpp"
#include "Signature.hpp"
#include "Term.hpp"
#include "TermIterators.hpp"

#include <cmath>



#include "Clause.hpp"

#undef RSTAT_COLLECTION
#define RSTAT_COLLECTION 1

namespace Kernel
{

using namespace std;
using namespace Lib;
using namespace Saturation;
using namespace Shell;

size_t Clause::_auxCurrTimestamp = 0;
#if VDEBUG
bool Clause::_auxInUse = false;
#endif

/** New clause */
Clause::Clause(Literal* const* lits, unsigned length, Inference inf)
  : Unit(Unit::CLAUSE, std::move(inf)),
    _length(length),
    _color(COLOR_INVALID),
    _extensionality(false),
    _extensionalityTag(false),
    _component(false),
    _store(NONE),
    _numSelected(0),
    _weight(0),
    _weightForClauseSelection(0),
    _refCnt(0),
    _reductionTimestamp(0),
    _literalPositions(0),
    _numActiveSplits(0),
    _auxTimestamp(0)
{
  // MS: TODO: not sure if this belongs here and whether EXTENSIONALITY_AXIOM input types ever appear anywhere (as a vampire-extension TPTP formula role)
  if(inference().inputType() == UnitInputType::EXTENSIONALITY_AXIOM){
    //cout << "Setting extensionality" << endl;
    _extensionalityTag = true;
    inference().setInputType(UnitInputType::AXIOM);
  }

  for(unsigned i = 0; i < length; i++) {
    (*this)[i] = lits[i];
  }

#if VAMPIRE_CLAUSE_TRACING
  // TODO make unsigned
  if (env.options->traceBackward() && unsigned(env.options->traceBackward()) == number()) {
    traverseParentsPost(
        [&](unsigned depth, Unit* unit) {
          std::cout << "backward trace " <<  number() << ": " << Output::repeat("| ", depth) << unit->toString() << std::endl;
      });
  }

  // forward tracing
  // TODO make unsigned
  static int traceFwd = env.options->traceForward();
  if (traceFwd != -1) {

    bool doTrace = false;
    auto infit = inference().iterator();
    while (inference().hasNext(infit)) {
      if (inference().next(infit)->number() == unsigned(traceFwd)) {
        doTrace = true;
        break;
      }
    }
    if (doTrace) {
      std::cout << "forward trace " << traceFwd << ": " << toString() << std::endl;
    }
  }

#endif // VAMPIRE_CLAUSE_TRACING
}

/**
 * Allocate a clause having lits literals.
 * @since 18/05/2007 Manchester
 */
void* Clause::operator new(size_t sz, unsigned lits)
{
  ASS_EQ(sz,sizeof(Clause));

  RSTAT_CTR_INC("clauses created");

  //We have to get sizeof(Clause) + (_length-1)*sizeof(Literal*)
  //this way, because _length-1 wouldn't behave well for
  //_length==0 on x64 platform.
  size_t size = sizeof(Clause) + lits * sizeof(Literal*);
  size -= sizeof(Literal*);

  return ALLOC_KNOWN(size,"Clause");
}

void Clause::operator delete(void* ptr,unsigned length)
{
  RSTAT_CTR_INC("clauses deleted by delete operator");

  //We have to get sizeof(Clause) + (_length-1)*sizeof(Literal*)
  //this way, because _length-1 wouldn't behave well for
  //_length==0 on x64 platform.
  size_t size = sizeof(Clause) + length * sizeof(Literal*);
  size -= sizeof(Literal*);

  DEALLOC_KNOWN(ptr, size,"Clause");
}

void Clause::destroyExceptInferenceObject()
{
  if (_literalPositions) {
    delete _literalPositions;
  }

  ConditionalRedundancyHandler::destroyClauseData(this);

  RSTAT_CTR_INC("clauses deleted");

  //We have to get sizeof(Clause) + (_length-1)*sizeof(Literal*)
  //this way, because _length-1 wouldn't behave well for
  //_length==0 on x64 platform.
  size_t size = sizeof(Clause) + _length * sizeof(Literal*);
  size -= sizeof(Literal*);

  DEALLOC_KNOWN(this, size,"Clause");
}


/**
 * Create a clause with the same content as @c c. The inference of the
 * created clause refers to @c c using the REORDER_LITERALS inference.
 *
 * The age of @c c is used, however the selected literals are not kept.
 *
 * Splitting history from @c c is also copied into the new clause.
 */
Clause* Clause::fromClause(Clause* c)
{
  Clause* res = fromIterator(c->iterLits(), SimplifyingInference1(InferenceRule::REORDER_LITERALS, c));

  if (c->splits()) {
    res->setSplits(c->splits());
  }

  return res;
}

bool Clause::shouldBeDestroyed()
{
  return (_store == NONE) && _refCnt == 0 &&
    !isFromPreprocessing();
}

/**
 * If storage is set to NONE, there are no references to this clause,
 * an it is not an input clause, destroy it.
 */
void Clause::destroyIfUnnecessary()
{
  if (shouldBeDestroyed()) {
    destroy();
  }
}

/**
 * Destroy the clause by deleting the clause itself and all of its
 * literals.
 * @since 19/05/2007 Manchester
 */
void Clause::destroy()
{
  static Stack<Clause*> toDestroy(32);
  Clause* cl = this;
  for(;;) {
    if (env.options->proofExtra() == Options::ProofExtra::FULL) {
      env.proofExtra.remove(cl);
    }
    Inference::Iterator it = cl->_inference.iterator();
    while (cl->_inference.hasNext(it)) {
      Unit* refU = cl->_inference.next(it);
      if (!refU->isClause()) {
        continue;
      }
      Clause* refCl = static_cast<Clause*> (refU);
      refCl->_refCnt--;
      if (refCl->shouldBeDestroyed()) {
        toDestroy.push(refCl);
      }
    }
    cl->_inference.destroyDirectlyOwned();
    cl->destroyExceptInferenceObject();
    if (toDestroy.isEmpty()) {
      break;
    }
    cl = toDestroy.pop();
  }
} // Clause::destroy

/** Set the store to @b s
 *
 * Can lead to clause deletion if the store is @b NONE
 * and there Clause's reference counter is zero. */
void Clause::setStore(Store s)
{
#if VDEBUG
  //assure there is one selected clause
  static Clause* selected=0;
  if (_store==SELECTED) {
    ASS_EQ(selected, this);
    selected=0;
  }
  if (s==SELECTED) {
    ASS_EQ(selected, 0);
    selected=this;
  }
#endif
  _store = s;
  destroyIfUnnecessary();
}

/**
 * Return true iff clause contains no variables
 */
bool Clause::isGround()
{
  return iterLits().all([](auto l) { return l->ground(); });
}

/**
 * Return true iff clause contains no literals of non-zero arity
 */
bool Clause::isPropositional()
{
  return iterLits().all([](auto l) { return l->arity() == 0; });
}

/**
 * Return true iff clause is Horn
 */
bool Clause::isHorn()
{
  bool posFound=false;
  for (Literal* l : iterLits()) {
    if (l->isPositive()) {
      if (posFound) {
        return false;
      }
      else {
        posFound=true;
      }
    }
  }
  return true;
}

/**
 * Return iterator over clause variables
 */
VirtualIterator<unsigned> Clause::getVariableIterator() const
{
  return pvi( getUniquePersistentIterator(
      getMappingIterator(
	  getMapAndFlattenIterator(
	      iterLits(),
	      VariableIteratorFn()),
	  OrdVarNumberExtractorFn())));
}

/**
 * Return true if the clause does not depend on any splits
 * in the backtracking splitting.
 */
bool Clause::noSplits() const
{
  return !_inference.splits() || _inference.splits()->isEmpty();
}

/**
 * Convert non-propositional part of the clause to std::string.
 */
std::string Clause::literalsOnlyToString() const
{
  if (_length == 0) {
    return "$false";
  } else {
    std::string result;
    result += _literals[0]->toString();
    for(unsigned i = 1; i < _length; i++) {
      result += " | ";
      result += _literals[i]->toString();
    }
    return result;
  }
}

/**
 * Convert the clause to the TPTP-compatible std::string representation.
 *
 * The split history is omitted.
 */
std::string Clause::toTPTPString() const
{
  std::string result = literalsOnlyToString();

  return result;
}

/**
 * Convert the clause to easily readable std::string representation.
 */
std::string Clause::toNiceString() const
{
  std::string result = literalsOnlyToString();

  if (splits() && !splits()->isEmpty()) {
    result += std::string(" {") + splits()->toString() + "}";
  }

  return result;
}

std::ostream& operator<<(std::ostream& out, Clause const& self)
{
  if (self.size() == 0) {
    return out << "$false";
  } else {
    out << *self[0];
    for (unsigned i = 1; i < self.size(); i++){
      out << " | " << *self[i];
    }
    if (self.splits() && !self.splits()->isEmpty()) {
      out << "{" << *self.splits() << "}";
    }
  }
  return out;
}


/**
 * Convert the clause to the std::string representation
 * Includes splitting, age, weight, selected and inference
 */
std::string Clause::toString() const
{
  // print id and literals of clause
  std::string result = Int::toString(_number) + ". " + literalsOnlyToString();

  // print avatar components clause depends on
  if (splits() && !splits()->isEmpty()) {
    result += std::string(" <- (") + Splitter::splitsToString(splits()) + ")";
  }

  // print inference and ids of parent clauses
  result += " " + inferenceAsString();

  if(env.options->proofExtra() != Options::ProofExtra::OFF){
    // print statistics: each entry should have the form key:value
    result += std::string(" {");
      
    result += std::string("a:") + Int::toString(age());
    unsigned weight = (_weight ? _weight : computeWeight());
    result += std::string(",w:") + Int::toString(weight);
    
    unsigned weightForClauseSelection = (_weightForClauseSelection ? _weightForClauseSelection : computeWeightForClauseSelection(*env.options));
    if(weightForClauseSelection!=weight){
      result += std::string(",wCS:") + Int::toString(weightForClauseSelection);
    }

    if (numSelected()>0) {
      result += std::string(",nSel:") + Int::toString(numSelected());
    }

    if (env.colorUsed) {
      result += std::string(",col:") + Int::toString(color());
    }

    if(derivedFromGoal()){
      result += std::string(",goal:1");
    }
    if(env.maxSineLevel > 1) { // this is a cryptic way of saying "did we run Sine to compute sine levels?"
      result += std::string(",sine:")+Int::toString((unsigned)_inference.getSineLevel());
    }

    if(isPureTheoryDescendant()){
      result += std::string(",ptD:1");
    }

    if(env.options->induction() != Shell::Options::Induction::NONE){
      result += std::string(",inD:") + Int::toString(_inference.inductionDepth());
    }
    result += ",thAx:" + Int::toString((int)(_inference.th_ancestors));
    result += ",allAx:" + Int::toString((int)(_inference.all_ancestors));

    result += ",thDist:" + Int::toString( _inference.th_ancestors * env.options->theorySplitQueueExpectedRatioDenom() - _inference.all_ancestors);
    result += std::string("}");
  }

  return result;
}

/**
 * Convert the clause into sequence of strings, each containing
 * a proper clause
 */
VirtualIterator<std::string> Clause::toSimpleClauseStrings()
{
    return pvi(getSingletonIterator(literalsOnlyToString()));
}

/**
 * Return true iff the clause is skipped for the purpose
 * of symbol elimination reporting.
 */
bool Clause::skip() const
{
  unsigned clen = length();
  for(unsigned i = 0; i < clen; i++) {
    const Literal* lit = (*this)[i];
    if (!lit->skip()) {
      return false;
    }
  }
  return true;
}

/**
 * Compute the color of the clause and store it in @b _color
 * @pre All literals are shared, so their color is computed properly.
 */
void Clause::computeColor() const
{
  ASS_EQ(_color, COLOR_INVALID);

  Color color = COLOR_TRANSPARENT;

  if (env.colorUsed) {
    unsigned clen=length();
    for(unsigned i=0;i<clen;i++) {
      color = static_cast<Color>(color | (*this)[i]->color());
    }
    ASS_L(color, COLOR_INVALID);
  }

  _color=color;
}

/**
 * Compute the weight of the clause.
 * @pre All literals are shared, so their weight is computed properly.
 * @since 02/01/2008 Manchester.
 * @since 22/01/2015 include splitWeight in weight
 */
unsigned Clause::computeWeight() const
{
  unsigned result = 0;
  for (int i = _length-1; i >= 0; i--) {
    ASS(_literals[i]->shared());
    result += _literals[i]->weight();
  }

  return result;
} // Clause::computeWeight


/**
 * Return weight of the split part of the clause
 *
 * This weight is not included in the number returned by the
 * @b weight() function.
 */
unsigned Clause::splitWeight() const
{
  return splits() ? splits()->size() : 0;
}

/**
 * Returns the numeral weight of a clause. The weight is defined as the sum of
 * binary sizes of all integers occurring in this clause.
 * @warning Each call to this function recomputes the numeral weight, so the call may
 *          potentially result in traversing the whole clause
 * @since 04/05/2013 Manchester, updated to use new NonVariableIterator
 * @author Andrei Voronkov
 */

unsigned Clause::getNumeralWeight() const {
  unsigned res = 0;
  for (Literal* lit : iterLits()) {
    if (!lit->hasInterpretedConstants()) {
      continue;
    }
    NonVariableIterator nvi(lit);
    while (nvi.hasNext()) {
      const Term* t = nvi.next().term();
      if (t->arity() != 0 || t->isSort()) {
        continue;
      }
      IntegerConstantType intVal;

      if (theory->tryInterpretConstant(t, intVal)) {
        int w = BitUtils::log2(Int::safeAbs(intVal.toInner())) - 1;
        if (w > 0) {
          res += w;
        }
        continue;
      }
      RationalConstantType ratVal;
      RealConstantType realVal;
      bool haveRat = false;

      if (theory->tryInterpretConstant(t, ratVal)) {
        haveRat = true;
      } else if (theory->tryInterpretConstant(t, realVal)) {
        ratVal = RationalConstantType(realVal);
        haveRat = true;
      }
      if (!haveRat) {
        continue;
      }
      int wN = BitUtils::log2(Int::safeAbs(ratVal.numerator().toInner())) - 1;
      int wD = BitUtils::log2(Int::safeAbs(ratVal.denominator().toInner())) - 1;
      int v = wN + wD;
      if (v > 0) {
        res += v;
      }
    }
  }
  return res;
} // getNumeralWeight

/**
 * compute weight of the clause used by clause selection and cache it
 */
unsigned Clause::computeWeightForClauseSelection(const Options& opt) const
{
  unsigned w = 0;
  if (_weight) {
    w = _weight;
  } else {
    w = computeWeight();
  }

  unsigned splWeight = 0;
  if (opt.nonliteralsInClauseWeight()) {
    splWeight = splitWeight(); // no longer includes propWeight
  }

  // hack: computation of getNumeralWeight is potentially expensive, so we only compute it if
  // the option increasedNumeralWeight is set to true.
  unsigned numeralWeight = 0;
  if (opt.increasedNumeralWeight())
  {
    numeralWeight = getNumeralWeight();
  }

  bool derivedFromGoal = Unit::derivedFromGoal();
  if(derivedFromGoal && opt.restrictNWCtoGC()){
    bool found = false;
    for(unsigned i=0;i<_length;i++){
      NonVariableNonTypeIterator it(_literals[i]);
      while(it.hasNext()){
        found |= env.signature->getFunction(it.next()->functor())->inGoal();
      }
    }
    if(!found){ derivedFromGoal=false; }
  }

  return Clause::computeWeightForClauseSelection(w, splWeight, numeralWeight, derivedFromGoal, opt);
}

/*
 * note: we currently assume in Clause::computeWeightForClauseSelection(opt) that numeralWeight is only used here if
 * the option increasedNumeralWeight() is set to true.
 */
unsigned Clause::computeWeightForClauseSelection(unsigned w, unsigned splitWeight, unsigned numeralWeight, bool derivedFromGoal, const Shell::Options& opt)
{
  static unsigned nongoalWeightCoeffNum = opt.nongoalWeightCoefficientNumerator();
  static unsigned nongoalWeightCoefDenom = opt.nongoalWeightCoefficientDenominator();

  w += splitWeight;

  if (opt.increasedNumeralWeight()) {
    w = (2 * w + numeralWeight);
  }
  return w * ( !derivedFromGoal ? nongoalWeightCoeffNum : nongoalWeightCoefDenom);
}


void Clause::collectUnstableVars(DHSet<unsigned>& acc)
{
  collectVars2<UnstableVarIt>(acc);
}

void Clause::collectVars(DHSet<unsigned>& acc)
{
  collectVars2<VariableIterator>(acc);
}

template<class VarIt>
void Clause::collectVars2(DHSet<unsigned>& acc)
{
  for (Literal* lit : iterLits()) {
    VarIt vit(lit);
    while (vit.hasNext()) {
      TermList var = vit.next();
      ASS(var.isOrdinaryVar());
      acc.insert(var.var());
    }
  }
}

unsigned Clause::varCnt()
{
  static DHSet<unsigned> vars;
  vars.reset();
  collectVars(vars);
  return vars.size();
}

unsigned Clause::maxVar()
{
  unsigned max = 0;
  VirtualIterator<unsigned> it = getVariableIterator();

  while (it.hasNext()) {
    unsigned n = it.next();
    max = n > max ? n : max;
  }
  return max;
}

unsigned Clause::numPositiveLiterals()
{
  unsigned count = 0;
  for (int i = 0; i < _length; i++)
  {
    Literal *lit = (*this)[i];
    if (lit->isPositive())
    {
      count++;
    }
  }
  return count;
}

/**
 * Return index of @b lit in the clause
 *
 * @b lit has to be present in the clause
 */
unsigned Clause::getLiteralPosition(Literal* lit)
{
  switch(length()) {
  case 1:
    ASS_EQ(lit,(*this)[0]);
    return 0;
  case 2:
    if (lit==(*this)[0]) {
      return 0;
    } else {
      ASS_EQ(lit,(*this)[1]);
      return 1;
    }
  case 3:
    if (lit==(*this)[0]) {
      return 0;
    } else if (lit==(*this)[1]) {
      return 1;
    } else {
      ASS_EQ(lit,(*this)[2]);
      return 2;
    }
#if VDEBUG
  case 0:
    ASSERTION_VIOLATION;
#endif
  default:
    if (!_literalPositions) {
      _literalPositions=new InverseLookup<Literal>(_literals,length());
    }
    return static_cast<unsigned>(_literalPositions->get(lit));
  }
}

/**
 * This method should be called when literals of the clause are
 * reordered (e.g. after literal selection), so that the information
 * about literal positions can be updated.
 */
void Clause::notifyLiteralReorder()
{
  if (_literalPositions) {
    _literalPositions->update(_literals);
  }
}

#if VDEBUG

void Clause::assertValid()
{
  ASS_ALLOC_TYPE(this, "Clause");
  if (_literalPositions) {
    unsigned clen=length();
    for (unsigned i = 0; i<clen; i++) {
      ASS_EQ(getLiteralPosition((*this)[i]),i);
    }
  }
}

#endif

bool Clause::contains(Literal* lit)
{
  for (int i = _length-1; i >= 0; i--) {
    if (_literals[i]==lit) {
      return true;
    }
  }
  return false;
}

std::ostream& operator<<(std::ostream& out, Clause::Store const& store) 
{
  switch (store) {
    case Clause::PASSIVE:     return out << "PASSIVE";
    case Clause::ACTIVE:      return out << "ACTIVE";
    case Clause::UNPROCESSED: return out << "UNPROCESSED";
    case Clause::NONE:        return out << "NONE";
    case Clause::SELECTED:    return out << "SELECTED";
  }
  ASSERTION_VIOLATION;
}

Literal* Clause::getAnswerLiteral() {
  for (unsigned i = 0; i < _length; ++i) {
    if (_literals[i]->isAnswerLiteral()) {
      return _literals[i];
    }
  }
  return nullptr;
}

bool Clause::computable() {
  for (unsigned i = 0; i < length(); ++i) {
    if ((*this)[i]->isAnswerLiteral()) {
      continue;
    }
    if (!(*this)[i]->computable()) {
      return false;
    }
  }
  return true;
}

}
