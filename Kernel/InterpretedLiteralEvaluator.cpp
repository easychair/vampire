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
 * @file InterpretedLiteralEvaluator.cpp
 * Implements class InterpretedLiteralEvaluator.
 */

#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"

#include "Signature.hpp"
#include "SortHelper.hpp"
#include "OperatorType.hpp"
#include "TermIterators.hpp"
#include "Term.hpp"
#include "Theory.hpp"
#include "NumTraits.hpp"


#include "InterpretedLiteralEvaluator.hpp"

#if VDEBUG
#define _DEBUG(...) 
#define DEBUG(...) //DBG(__VA_ARGS__)
#else 
#define DEBUG(...)
#define _DEBUG(...)
#endif

namespace Kernel
{
using namespace std;
using namespace Lib;

struct PredEvalResult {
  enum status_t {
    Simplified,
    Trivial,
    Nop,
  } status;
  union {
    bool trivial_val;
    Literal* simplified_val;
  };
  static PredEvalResult nop() {return  PredEvalResult(Nop); }
  static PredEvalResult trivial(bool value) {return  PredEvalResult (value); }
  static PredEvalResult simplified(Literal* value) {return  PredEvalResult (value); }

private:
  explicit PredEvalResult(bool value) : status(Trivial), trivial_val(value) {}
  explicit PredEvalResult(Literal* value) : status(Simplified), simplified_val(value) {}
  explicit PredEvalResult(status_t stat) : status(stat) {}
};

/**
 * We use descendants of this class to evaluate various functions.
 *
 * One function can be evaluated only by one Evaluator object.
 */
class InterpretedLiteralEvaluator::Evaluator
{
public:
  virtual ~Evaluator() {}

  virtual bool canEvaluateFunc(unsigned func) { return false; }
  virtual bool canEvaluatePred(unsigned pred) { return false; }

  virtual bool tryEvaluateFunc(Term* trm, TermList& res) { return false; }
  virtual PredEvalResult tryEvaluatePred(Literal* trm)  { return PredEvalResult::nop(); }
};

  /* Traverses the given term. Only the subtree parts where pred is true are traversed.
   * For every traversed term where pred is *not* true, action is applied.
   *
   * The argument functions shall have the following signatures:
   * bool pred(const Term&)
   * void action(TermList)
   */
template<class Fn, class Predicate> 
void stackTraverseIf(TermList term, Predicate pred, Fn action) {

  Stack<TermList> todo;
  todo.push(term);

  while(!todo.isEmpty()){
    TermList t = todo.pop();
    if(t.isTerm()) {
      auto& trm = *t.term();
      if (pred(trm)) {
        for (unsigned i = 0; i < trm.numTermArguments(); i++) {
          todo.push(trm[i]);
        }
      } else {
        action(t);
      }
    } else {
      action(t);
    }
  }
}


/**
 * We want to simplify terms that are interpreted by abelian groups. e.g. (1+a)+1 -> 2 + a ...
 * the standard evaluation will not do this. 
 *
 * Additionally evaluator has a weekly normalizing behaviour. Namely it re-brackets sums, such that the lhs 
 * of the operation is always a non-operator term. Further all interpreted constants will be collapsed into 
 * the 'left-most' term. The left most term is omitted if it is the identity element.
 *
 * x + ( y + ( t + 4 ) + r ) + 5  ==> ( ( (9 + x) + y ) + t ) + r
 * x + ( y + 0 )                  ==> x + y
 *
 * (The name of this class comes from the Associative Commutative operation of the Group)
 *
 * @author Giles (refactorings by joe-hauns)
 * @since 06/12/18
 */
template<class AbelianGroup>
  class InterpretedLiteralEvaluator::ACFunEvaluator
   : public Evaluator
{
public:
  using ConstantType = typename AbelianGroup::ConstantType;

  ACFunEvaluator() : _fun(env.signature->getInterpretingSymbol(AbelianGroup::interpreation)) { }
  const unsigned _fun; 

  virtual bool canEvaluateFunc(unsigned func) { return func == _fun; }

  virtual bool tryEvaluateFunc(Term* trm, TermList& res) { 
    ASS_EQ(trm->functor(), _fun);
    ASS_EQ(trm->numTermArguments(),2);

    unsigned nums = 0;
    ConstantType acc = AbelianGroup::IDENTITY;
    Stack<TermList> keep;
    stackTraverseIf(TermList(trm), 
        /* we traverse only the parts with the same operation */
        [&](Term& t){ return t.functor() == _fun; },
        [&](TermList t) {
          ConstantType c;
          /* we eval constant parts */
          if (t.isTerm() && theory->tryInterpretConstant(t.term(), c)) {
            acc = AbelianGroup::groundEval(acc, c);
            nums++;
          } else {
            keep.push(t);
          }
        });
    if (nums == 0 || (nums == 1 && acc != AbelianGroup::IDENTITY)) {
      // only continue, if we updated acc more than once or at least saw a "_+0" somewhere
      return false;
    }

    if (acc != AbelianGroup::IDENTITY) {
      keep.push(TermList(theory->representConstant(acc)));
    }

    auto iter = Stack<TermList>::BottomFirstIterator(keep);
    if (!iter.hasNext()) {
      res = TermList(theory->representConstant(AbelianGroup::IDENTITY));
      return TermList(trm) != res;
    } else {
      TermList out = iter.next();
      while (iter.hasNext()) {
        auto t = iter.next();
        out = TermList(Term::create2(_fun, t, out));
      }
      res = out;
      return TermList(trm) != res;
    }
  }
};

template<class Inequality>
class InequalityNormalizer {

public:
  static Literal* normalize(Literal* in) {
    ASS(in->functor() == Inequality::functor());

    if (Inequality::isNormalized(in)) {
      /* nothing to do */
      return in;
    } else {
      return Inequality::normalizedLit(
          in->polarity(),
          in->termArg(0), 
          in->termArg(1));
    }
  }

};

template<class ConstantType>
class FracLess { 
  template<class Inequality> friend class InequalityNormalizer;

  using number = NumTraits<ConstantType>;
  inline static unsigned functor() { return number::lessF(); }

  static Literal* normalizedLit(bool polarity, TermList lhs, TermList rhs) {
    static auto zero = TermList(number::zeroT());
    return number::less(
              polarity,
              zero,
              number::add(rhs, number::minus(lhs)));
  }

  inline static bool isNormalized(Literal* in) {
    return number::isZero(in->termArg(0));
  }

}; 


class IntLess { 
  template<class Inequality> friend class InequalityNormalizer;

  using number = NumTraits<IntegerConstantType>;
  inline static unsigned functor() { return number::lessF(); }

  static Literal* normalizedLit(bool polarity, TermList lhs, TermList rhs) {
    static auto one  = TermList(number::oneT());
    static auto zero = TermList(number::zeroT());
    if (polarity) {
      return number::less(
              /* polarity */ true, 
              zero,
              number::add(rhs, number::minus(lhs)));
    } else {
      /* ~(l < r) ==> (r <= l) ==> r - 1 < l ==> 0 < l - r + 1 */
      return number::less(
              /* polarity */ true, 
              zero,
              number::add(number::add(lhs, one), number::minus(rhs)));
    }
  }

  inline static bool isNormalized(Literal* in) {
    return number::isZero(in->termArg(0)) && in->polarity();
  }

}; 


/**
 * Interpreted equality has to be treated specially. We do not have separate
 * predicate symbols for different kinds of equality so the sorts must be 
 * detected and the correct intepretation of constants carried out.
 *
 * Equality is only decided between constant terms.
 *
 * @author Giles
 * @since 10/11/14
 */
class InterpretedLiteralEvaluator::EqualityEvaluator
  : public Evaluator
{
  bool canEvaluatePred(unsigned pred) override {
    return Signature::isEqualityPredicate(pred);
  }

  template<typename T>
  bool checkEquality(Literal* lit, bool& res)
  {
    T arg1;
    if(!theory->tryInterpretConstant(lit->termArg(0).term(),arg1)){ 
      return false; 
    }
    T arg2;
    if(!theory->tryInterpretConstant(lit->termArg(1).term(),arg2)){ 
      return false;
    }

    res = (arg1 == arg2);

    return true;
  }

  /** is to be called when interpreted functions were already evaluated */
  PredEvalResult tryEvaluatePred(Literal* lit_) override
  {
    auto& lit = *lit_;
    ASS(lit.isEquality());
    TermList l = lit[0];
    TermList r = lit[1];
    if (l == r) {
      return PredEvalResult::trivial(lit.polarity());
    } else {
      // l != r
      // TODO lift to evaluate inductive datatypes as well
      if (theory->isInterpretedNumber(l) && theory->isInterpretedNumber(r)) {
        return PredEvalResult::trivial(!lit.polarity());
      } else {
        return PredEvalResult::nop();
      }
    }
  }
};

/**
 * An evaluator for dealing with conversions between sorts
 *
 */
class InterpretedLiteralEvaluator::ConversionEvaluator
  : public Evaluator
{
public:
  bool canEvaluateFunc(unsigned func) override
  {
    if (!theory->isInterpretedFunction(func)) {
      return false;
    }
    return theory->isConversionOperation(theory->interpretFunction(func));
  }

  virtual bool tryEvaluateFunc(Term* trm, TermList& res) override
  {
    ASS(theory->isInterpretedFunction(trm));

    try {
      Interpretation itp = theory->interpretFunction(trm);
      ASS(theory->isFunction(itp));
      ASS(theory->isConversionOperation(itp));
      ASS_EQ(theory->getArity(itp), 1);

      TermList argTrm = trm->termArg(0);
      switch(itp) {
      case Theory::INT_TO_RAT:
	{
	  IntegerConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) { 
	    return false;
	  }
	  RationalConstantType resNum(arg,1);
	  res = TermList(theory->representConstant(resNum));
	  return true;
	}
      case Theory::INT_TO_REAL:
	{
	  IntegerConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) {
	    return false;
	  }
	  RealConstantType resNum(RationalConstantType(arg,1));
	  res = TermList(theory->representConstant(resNum));
	  return true;
	}
      case Theory::RAT_TO_INT:
	{
	  RationalConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) {
	    return false;
	  }
	  IntegerConstantType resNum = IntegerConstantType::floor(arg);
	  res = TermList(theory->representConstant(resNum));
	  return true;
	}
      case Theory::RAT_TO_REAL:
	{
	  RationalConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) {
	    return false;
	  }
	  RealConstantType resNum(arg);
	  res = TermList(theory->representConstant(resNum));
	  return true;
	}
      case Theory::REAL_TO_INT:
	{
	  RealConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) {
	    return false;
	  }
	  IntegerConstantType resNum = IntegerConstantType::floor(RationalConstantType(arg));
	  res = TermList(theory->representConstant(resNum));
	  return true;
	}
      case Theory::REAL_TO_RAT:
	{
	  //this is correct only as long as we only represent rational real numbers
	  RealConstantType arg;
	  if (!theory->tryInterpretConstant(argTrm, arg)) {
	    return false;
	  }
	  RationalConstantType resNum(arg);
	  res = TermList(theory->representConstant(resNum));
	return true;
      }

      default:
	ASSERTION_VIOLATION;
      }
    }
    catch(ArithmeticException&)
    {
      return false;
    }
  }
};

/**
 * Evaluates constant theory expressions
 *
 * Evaluators for each sort implement tryEvaluate(Unary/Binary)(Func/Pred) 
 * 
 */
template<class T>
class InterpretedLiteralEvaluator::TypedEvaluator : public Evaluator
{
public:
  using Value = T;
  using number = NumTraits<Value>;

  TypedEvaluator() {}

  bool isZero(T arg) const { return number::constant(0) == arg; }
  TermList getZero() const {return number::zero(); }
  bool isOne(T arg) const { return number::constant(1) == arg; }
  bool isMinusOne(T arg) const { return typename number::ConstantType(-1) == arg; }
  TermList invert(TermList t) const { return number::minus(t); }
  bool isAddition(Interpretation interp) const { return interp == number::addI; }
  bool isProduct(Interpretation interp) const  { return interp == number::mulI; }
  virtual bool isDivision(Interpretation interp) const = 0;

  virtual bool canEvaluate(Interpretation interp)
  {
    // This is why we cannot evaluate Equality here... we cannot determine its sort
    if (!theory->hasSingleSort(interp)) { return false; } //To skip conversions and EQUAL

    if (theory->isPolymorphic(interp)) { return false; } // typed evaulator not for polymorphic stuff

    TermList opSort = theory->getOperationSort(interp);
    return opSort==T::getSort();
  }

  virtual bool tryEvaluateFunc(Term* trm, TermList& res) override
  {
    ASS(theory->isInterpretedFunction(trm));
    const auto num = NumTraits<Value>{};

    _DEBUG( "try evaluate ", trm->toString() );

    try {
      Interpretation itp = theory->interpretFunction(trm);
      ASS(theory->isFunction(itp));
      unsigned arity = theory->getArity(itp);

      if (arity!=1 && arity!=2) {
	INVALID_OPERATION("unsupported arity of interpreted operation: "+Int::toString(arity));
      }
      T resNum;
      TermList arg1Trm = trm->termArg(0);
      T arg1;
      if (arity==1) {
        if (theory->tryInterpretConstant(arg1Trm, arg1)){
          if (!tryEvaluateUnaryFunc(itp, arg1, resNum)) { return false;}
        } else if (itp == num.minusI){ 
          const unsigned umin = trm->functor();
          return trySimplifyUnaryMinus(umin, arg1Trm, res);
        } else{ 
          return false;
        }
      }
      else if(arity==2){


        // If one argument is not a constant and the other is zero, one or minus one then
        // we might have some special cases
        T arg2;
        TermList arg2Trm = trm->termArg(1);

        bool specialCase = true;
        T conArg;
        TermList nonConTerm;
        if (theory->tryInterpretConstant(arg1Trm, arg1) && (isZero(arg1) || isOne(arg1) || isMinusOne(arg1)) && 
            !theory->tryInterpretConstant(arg2Trm, arg2)) {
         conArg = arg1;
         nonConTerm = arg2Trm;
        }
        else if(theory->tryInterpretConstant(arg2Trm, arg2) && (isZero(arg2) || isOne(arg2) || isMinusOne(arg2)) && 
            !theory->tryInterpretConstant(arg1Trm, arg1)) {
         conArg = arg2;
         nonConTerm = arg1Trm;
        }
        else{
          specialCase = false;
        }
        if(specialCase){
	_DEBUG( "special case" );
 
          //Special case where itp is division and arg2 is '1'
          //   Important... this is a non-symmetric case!
          if(theory->tryInterpretConstant(arg2Trm, arg2) && isOne(arg2) && isDivision(itp)){
            res = arg1Trm;
            return true;
          }
          //Special case where itp is addition and conArg is '0'
          if(isZero(conArg) && isAddition(itp)){
            res = nonConTerm;
            return true;
          }
          //Special case where itp is multiplication and conArg  is '1'
          if(isOne(conArg) && isProduct(itp)){
            res = nonConTerm;
            return true;
          }
          //Special case where itp is multiplication and conArg  is '-1'
          if(isMinusOne(conArg) && isProduct(itp)){
            res = invert(nonConTerm); 
            return true;
          }
          //Special case where itp is multiplication and conArg is '0'
          if(isZero(conArg) && isProduct(itp)){
            res = getZero();
            return true;
          }
        }
        if(theory->tryInterpretConstant(arg1Trm, arg1) && theory->tryInterpretConstant(arg2Trm, arg2)){
	  if (!tryEvaluateBinaryFunc(itp, arg1, arg2, resNum)) { return false;}
        }
        else{ return false;}
      }
      // TODO FIXME: this point should not be reachable due to check for !(arity!=1 && arity!=2) before.
      res = TermList(theory->representConstant(resNum));
      return true;
    }
    catch(ArithmeticException&)
    {
       DEBUG( "ArithmeticException" );
      return false;
    }
  }

  virtual PredEvalResult tryEvaluatePred(Literal* lit) override
  {
    ASS(theory->isInterpretedPredicate(lit->functor()));
    bool res;

    try {
      Interpretation itp = theory->interpretPredicate(lit);
      ASS(!theory->isFunction(itp));
      unsigned arity = theory->getArity(itp);

      if (arity!=1 && arity!=2) {
	INVALID_OPERATION("unsupported arity of interpreted operation: "+Int::toString(arity));
      }
      TermList arg1Trm = lit->termArg(0);
      T arg1;
      if (!theory->tryInterpretConstant(arg1Trm, arg1)) { return PredEvalResult::nop(); }
      if (arity==1) {
	if (!tryEvaluateUnaryPred(itp, arg1, res)) { return PredEvalResult::nop();}
      }
      else {
	TermList arg2Trm = lit->termArg(1);
	T arg2;
	if (!theory->tryInterpretConstant(arg2Trm, arg2)) { return PredEvalResult::nop(); }
	if (!tryEvaluateBinaryPred(itp, arg1, arg2, res)) { return PredEvalResult::nop();}
      }
      if (lit->isNegative()) {
	res = !res;
      }
      return PredEvalResult::trivial(res);
    }
    catch(ArithmeticException&)
    {
      return PredEvalResult::nop();
    }

  }

  bool canEvaluateFunc(unsigned func) override
  {
    if (!theory->isInterpretedFunction(func)) {
      return false;
    }
    Interpretation interp = theory->interpretFunction(func);
    return canEvaluate(interp);
  }

  bool canEvaluatePred(unsigned pred) override
  {
    if (!theory->isInterpretedPredicate(pred)) {
      return false;
    }
    Interpretation interp = theory->interpretPredicate(pred);
    return canEvaluate(interp);
  }

protected:
  /** Tries to simplify a term involving unary minus. 
   * the term to be simplified is uminus(innner)
   * uminus functor argument is passed for the sake of perfomance
   */
  bool trySimplifyUnaryMinus(const unsigned& uminus_functor, const TermList& inner, TermList& result)
  { 
    DEBUG("trySimplifyUnaryMinus(uminus(", inner, "))")
    ASS_EQ(uminus_functor, env.signature->getInterpretingSymbol(number::minusI));
    if (inner.isTerm()) {
      /* complex term */
      auto& t = *inner.term();
      if (theory->isInterpretedFunction(t.functor())) {

        /* interpreted function */
        auto itp = theory->interpretFunction(t.functor());
        switch (itp) {
          case number::minusI:
            ASS_EQ(t.numTermArguments(), 1);
            result = t[0];
            return true;

          // case number::addI:
          //   ASS_EQ(t.numTermArguments(), 2);
          //   result = TermList(Term::create2(t.functor(), 
          //       simplifyUnaryMinus(uminus_functor, t[0]),
          //       simplifyUnaryMinus(uminus_functor, t[1])));
          //   return true; 

          default:
            /* interpreted function for which minus is not handled minus is not handled */
            return false;
        }

      } else {
        /* not an interpreted function */
        Value cons;
        if (theory->tryInterpretConstant(&t, cons)) {
          result = TermList(theory->representConstant(-cons));
          return true;
        } else {
          /* uninterpreted function */
          return false;
        }
      }
    } else {
      /* not a complex term */
      return false;
    }
  }

  TermList simplifyUnaryMinus(const unsigned& uminus_functor, const TermList& inner)
  {
    TermList out;
    if (trySimplifyUnaryMinus(uminus_functor, inner, out)) {
      return out;
    } else {
      return TermList(Term::create1(uminus_functor, inner));
    }
  }

  virtual bool tryEvaluateUnaryFunc(Interpretation op, const T& arg, T& res)
  { return false; }
  virtual bool tryEvaluateBinaryFunc(Interpretation op, const T& arg1, const T& arg2, T& res)
  { return false; }

  virtual bool tryEvaluateUnaryPred(Interpretation op, const T& arg1, bool& res)
  { return false; }
  virtual bool tryEvaluateBinaryPred(Interpretation op, const T& arg1, const T& arg2, bool& res)
  { return false; }
};

/**
 * Evaluates integer functions
 */
class InterpretedLiteralEvaluator::IntEvaluator : public TypedEvaluator<IntegerConstantType>
{
protected:

  virtual bool isDivision(Interpretation interp) const { 
    return interp==Theory::INT_QUOTIENT_E || interp==Theory::INT_QUOTIENT_T || 
           interp==Theory::INT_QUOTIENT_F; 
  }

  virtual bool tryEvaluateUnaryFunc(Interpretation op, const Value& arg, Value& res)
  {
    switch(op) {
    case Theory::INT_UNARY_MINUS:
      res = -arg;
      return true;
    case Theory::INT_ABS:
      if (arg < 0) {
        res = -arg;
      } else {
        res = arg;
      }
      return true;
    case Theory::INT_SUCCESSOR:
      res = arg+1;
      return true;
    case Theory::INT_FLOOR:
    case Theory::INT_CEILING:
    case Theory::INT_TRUNCATE:
    case Theory::INT_ROUND:
       // For integers these do nothing
      res = arg;
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryFunc(Interpretation op, const Value& arg1,
      const Value& arg2, Value& res)
  {
    switch(op) {
    case Theory::INT_PLUS:
      res = arg1+arg2;
      return true;
    case Theory::INT_MINUS:
      res = arg1-arg2;
      return true;
    case Theory::INT_MULTIPLY:
      res = arg1*arg2;
      return true;
    case Theory::INT_QUOTIENT_E:
      res = arg1.quotientE(arg2); // should be equivalent to arg1/arg2
      return true;
    case Theory::INT_QUOTIENT_T:
      res = arg1.quotientT(arg2);
      return true;
    case Theory::INT_QUOTIENT_F:
      res = arg1.quotientF(arg2);
      return true;
    // The remainder is left - (quotient * right)
    case Theory::INT_REMAINDER_E:
      res = arg1.remainderE(arg2);
      return true;
    case Theory::INT_REMAINDER_T:
      res = arg1 - (arg1.quotientT(arg2)*arg2);
      return true;
    case Theory::INT_REMAINDER_F:
      res = arg1 - (arg1.quotientF(arg2)*arg2);
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryPred(Interpretation op, const Value& arg1,
      const Value& arg2, bool& res)
  {
    switch(op) {
    case Theory::INT_GREATER:
      res = arg1>arg2;
      return true;
    case Theory::INT_GREATER_EQUAL:
      res = arg1>=arg2;
      return true;
    case Theory::INT_LESS:
      res = arg1<arg2;
      return true;
    case Theory::INT_LESS_EQUAL:
      res = arg1<=arg2;
      return true;
    case Theory::INT_DIVIDES:
      res = arg1.divides(arg2);
      return true;
    default:
      return false;
    }
  }
};

/**
 * Evaluations rational functions
 */
class InterpretedLiteralEvaluator::RatEvaluator : public TypedEvaluator<RationalConstantType>
{
protected:
  virtual bool isDivision(Interpretation interp) const { 
    return interp==Theory::RAT_QUOTIENT || interp==Theory::RAT_QUOTIENT_E || 
           interp==Theory::RAT_QUOTIENT_T || interp==Theory::RAT_QUOTIENT_F;
  }

  virtual bool tryEvaluateUnaryFunc(Interpretation op, const Value& arg, Value& res)
  {
    switch(op) {
    case Theory::RAT_UNARY_MINUS:
      res = -arg;
      return true;
    case Theory::RAT_FLOOR:
      res = arg.floor();
      return true;
    case Theory::RAT_CEILING:
      res = arg.ceiling();
      return true;
    case Theory::RAT_TRUNCATE:
      res = arg.truncate();
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryFunc(Interpretation op, const Value& arg1,
      const Value& arg2, Value& res)
  {
    switch(op) {
    case Theory::RAT_PLUS:
      res = arg1+arg2;
      return true;
    case Theory::RAT_MINUS:
      res = arg1-arg2;
      return true;
    case Theory::RAT_MULTIPLY:
      res = arg1*arg2;
      return true;
    case Theory::RAT_QUOTIENT:
      if (arg2 == RationalConstantType(0)) return false;
      res = arg1/arg2;
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryPred(Interpretation op, const Value& arg1,
      const Value& arg2, bool& res)
  {
    switch(op) {
    case Theory::RAT_GREATER:
      res = arg1>arg2;
      return true;
    case Theory::RAT_GREATER_EQUAL:
      res = arg1>=arg2;
      return true;
    case Theory::RAT_LESS:
      res = arg1<arg2;
      return true;
    case Theory::RAT_LESS_EQUAL:
      res = arg1<=arg2;
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateUnaryPred(Interpretation op, const Value& arg1,
      bool& res)
  {
    switch(op) {
    case Theory::RAT_IS_INT:
      res = arg1.isInt();
      return true;
    default:
      return false;
    }
  }
};

/**
 * Evaluates real functions. 
 * As reals are represented as rationals the operations are for reals.
 * See Kernel/Theory.hpp for how these operations are defined
 */
class InterpretedLiteralEvaluator::RealEvaluator : public TypedEvaluator<RealConstantType>
{
protected:
  virtual bool isDivision(Interpretation interp) const { 
    return interp==Theory::REAL_QUOTIENT || interp==Theory::REAL_QUOTIENT_E ||
           interp==Theory::REAL_QUOTIENT_T || interp==Theory::REAL_QUOTIENT_F;
  }

  virtual bool tryEvaluateUnaryFunc(Interpretation op, const Value& arg, Value& res)
  {
    switch(op) {
    case Theory::REAL_UNARY_MINUS:
      res = -arg;
      return true;
    case Theory::REAL_FLOOR:
      res = arg.floor();
      return true;
    case Theory::REAL_CEILING:
      res = arg.ceiling();
      return true;
    case Theory::REAL_TRUNCATE:
      res = arg.truncate();
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryFunc(Interpretation op, const Value& arg1,
      const Value& arg2, Value& res)
  {
    switch(op) {
    case Theory::REAL_PLUS:
      res = arg1+arg2;
      return true;
    case Theory::REAL_MINUS:
      res = arg1-arg2;
      return true;
    case Theory::REAL_MULTIPLY:
      res = arg1*arg2;
      return true;
    case Theory::REAL_QUOTIENT:
      if (arg2 == RealConstantType(0)) return false;
      res = arg1/arg2;
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateBinaryPred(Interpretation op, const Value& arg1,
      const Value& arg2, bool& res)
  {
    switch(op) {
    case Theory::REAL_GREATER:
      res = arg1>arg2;
      return true;
    case Theory::REAL_GREATER_EQUAL:
      res = arg1>=arg2;
      return true;
    case Theory::REAL_LESS:
      res = arg1<arg2;
      return true;
    case Theory::REAL_LESS_EQUAL:
      res = arg1<=arg2;
      return true;
    default:
      return false;
    }
  }

  virtual bool tryEvaluateUnaryPred(Interpretation op, const Value& arg1,
      bool& res)
  {
    switch(op) {
    case Theory::REAL_IS_INT:
      res = arg1.isInt();
      return true;
    case Theory::REAL_IS_RAT:
      //this is true as long as we can evaluate only rational reals.
      res = true;
      return true;
    default:
      return false;
    }
  }

};

template<Theory::Interpretation op>
struct AbelianGroup;

/** Creates an instance of struct AbelianGroup<oper>, for the use in ACFunEvaluator. */
#define IMPL_OPERATOR(oper, type, identity, eval) \
  template<> struct AbelianGroup<oper> { \
    const static Theory::Interpretation interpreation = oper; \
    using ConstantType = type; \
    const static type IDENTITY; \
    static type groundEval(type l, type r) { return eval; } \
    /*const static unsigned FUNCTOR;*/ \
  }; \
  const type     AbelianGroup<oper>::IDENTITY = identity; \

/* int opeators */
IMPL_OPERATOR(Theory::INT_MULTIPLY, IntegerConstantType, IntegerConstantType(1), l * r)
IMPL_OPERATOR(Theory::INT_PLUS, IntegerConstantType, IntegerConstantType(0), l + r)

/* rational opeators */
IMPL_OPERATOR(Theory::RAT_MULTIPLY, RationalConstantType, RationalConstantType(1), l * r)
IMPL_OPERATOR(Theory::RAT_PLUS, RationalConstantType, RationalConstantType(0), l + r)

/* real opeators */
IMPL_OPERATOR(Theory::REAL_MULTIPLY, RealConstantType, RealConstantType(RationalConstantType(1)), l * r)
IMPL_OPERATOR(Theory::REAL_PLUS, RealConstantType, RealConstantType(RationalConstantType(0)), l + r)

////////////////////////////////
// InterpretedLiteralEvaluator
//
// This is where the evaluators defined above are used.

InterpretedLiteralEvaluator::InterpretedLiteralEvaluator(bool doNormalize) : _normalize(doNormalize)
{
  // For an evaluator to be used it must be pushed onto _evals
  // We search this list, calling canEvaluate on each evaluator
  // An invariant we want to maintain is that for any literal only one
  //  Evaluator will return true for canEvaluate
  _evals.push(new IntEvaluator());
  _evals.push(new RatEvaluator());
  _evals.push(new RealEvaluator());
  _evals.push(new ConversionEvaluator());
  _evals.push(new EqualityEvaluator());

  if(env.options->useACeval()){

  // Special AC evaluators are added to be tried first for Plus and Multiply

  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::INT_PLUS>>()); 
  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::INT_MULTIPLY>>());
  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::RAT_PLUS>>());
  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::RAT_MULTIPLY>> ());
  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::REAL_PLUS>> ());
  _evals.push(new ACFunEvaluator<AbelianGroup<Theory::REAL_MULTIPLY>> ());

  }

  _funEvaluators.ensure(0);
  _predEvaluators.ensure(0);

}

InterpretedLiteralEvaluator::~InterpretedLiteralEvaluator()
{
  while (_evals.isNonEmpty()) {
    delete _evals.pop();
  }
}


/**
 * This checks if a literal is 'balancable' i.e. can be put into the form term=constant or term=var
 * 
 * This is still an experimental process and will be expanded/reworked later
 *
 * @author Giles
 * @since 11/11/14
 */
bool InterpretedLiteralEvaluator::balancable(Literal* lit)
{
  // Check that lit is compatible with this balancing operation
  // One thing that we cannot check, but assume is that it has already been simplified once
  // balance applies further checks

  // we can only rebalance for equality. 
  // (Inequalities would be possible as well but this would interfer with the normalization ( t < s ==> 0 < s - t ))
  if (!lit->isEquality()) return false;
  ASS_EQ(lit->numTermArguments(), 2)

  // one side must be a constant and the other interpretted
  // the other side can contain at most one variable or uninterpreted subterm 
  // but we do not check this second condition here, instead we detect it in balance
  TermList t1 = lit->termArg(0);
  TermList t2 = lit->termArg(1);

  bool t1Number = theory->isInterpretedNumber(t1);
  bool t2Number = theory->isInterpretedNumber(t2);

  if(!t1Number && !t2Number){ return false; } // cannot balance
  if(t1Number && t2Number){ return true; } // already balanced

  // so here exactly one of t1Number and t2Number is true

  if(t1Number){
    if(t2.isVar()){ return false;} // already balanced
    if(!theory->isInterpretedFunction(t2)){ return false;} // cannot balance
  }
  if(t2Number){
    if(t1.isVar()){ return false;} // already balanced
    if(!theory->isInterpretedFunction(t1)){ return false;} // cannot balance
  }

  return true;
}

/**
 * This attempts to 'balance' a literal i.e. put it into the form term=constant
 *
 * The current implementation is only applicable to a restricted set of cases.
 *
 * This is still an experimental process and will be expanded/reworked later
 *
 * @author Giles
 * @since 11/11/14
 */
bool InterpretedLiteralEvaluator::balance(Literal* lit,Literal*& resLit)
{
  ASS(balancable(lit));

  _DEBUG( "try balance ", lit->toString() );

  ASS(theory->isInterpretedPredicate(lit->functor()));

  // at the end this tells us if the predicate needs to be swapped, only applies if non-equality
  bool swap = false; 

  // only want lesseq and equality
  if(lit->numTermArguments()!=2) return false;

  TermList t1;
  TermList t2;
  // ensure that t1 is the constant
  if(theory->isInterpretedNumber(lit->termArg(0))){
    t1 = lit->termArg(0); t2 = lit->termArg(1);
  }else{
    t1 = lit->termArg(1); t2 = lit->termArg(0);
    swap=true;
  }
  // so we have t1 a constant and t2 something that has an interpreted function at the top

  Signature::Symbol* conSym = env.signature->getFunction(t1.term()->functor());
  TermList srt;
  if(conSym->integerConstant()) srt = AtomicSort::intSort();
  else if(conSym->rationalConstant()) srt = AtomicSort::rationalSort();
  else if(conSym->realConstant()) srt = AtomicSort::realSort();
  else{
     ASSERTION_VIOLATION_REP(t1);
    return false;// can't work out the sort, that's odd!
  }

  // unwrap t2, applying the wrappings to t1 until we get to something we can't unwrap
  // this updates t1 and t2 as we go

  // This relies on the fact that a simplified literal with a single non-constant
  // subterm will look like f(c,f(c,f(c,t)))=c
  // If we cannot invert an interpreted function f then we fail

  bool modified = false;

  while(theory->isInterpretedFunction(t2)){
    TermList* args = t2.term()->args();
    
    // find which arg of t2 is the non_constant bit, this is what we are unwrapping 
    TermList to_unwrap = TermList::empty();
    while(args->isNonEmpty()){
      if(!theory->isInterpretedNumber(*args)){
        if(!to_unwrap.isEmpty()){
          return false; // If there is more than one non-constant term this will not work
        }
        to_unwrap=*args;
      } 
      args= args->next();
    }
    //Should not happen if balancable passed and it was simplified
    if(to_unwrap.isEmpty()){ return false;} 
    
    // Now we do a case on the functor of t2
    Term* t2term = t2.term();
    Interpretation t2interp = theory->interpretFunction(t2term->functor());
    TermList result;
    bool okay=true;
    switch(t2interp){
      case Theory::INT_PLUS:
        okay=balancePlus(Theory::INT_PLUS,Theory::INT_UNARY_MINUS,t2term,to_unwrap,t1,result);
        break;
      case Theory::RAT_PLUS:
        okay=balancePlus(Theory::RAT_PLUS,Theory::RAT_UNARY_MINUS,t2term,to_unwrap,t1,result);
        break;
      case Theory::REAL_PLUS:
        okay=balancePlus(Theory::REAL_PLUS,Theory::REAL_UNARY_MINUS,t2term,to_unwrap,t1,result);
        break;

      case Theory::INT_MULTIPLY: 
      {
        okay=balanceIntegerMultiply(t2term,to_unwrap,t1,result,swap);
        break;
      }
      case Theory::RAT_MULTIPLY:
      {
        RationalConstantType zero(0,1);
        okay=balanceMultiply(Theory::RAT_QUOTIENT,zero,t2term,to_unwrap,t1,result,swap);
        break;
      }
      case Theory::REAL_MULTIPLY:
      {
        RealConstantType zero(RationalConstantType(0, 1));
        okay=balanceMultiply(Theory::REAL_QUOTIENT,zero,t2term,to_unwrap,t1,result,swap);
        break;
       }

      case Theory::RAT_QUOTIENT:
        okay=balanceDivide(Theory::RAT_MULTIPLY,t2term,to_unwrap,t1,result,swap);
        break;
      case Theory::REAL_QUOTIENT:
        okay=balanceDivide(Theory::REAL_MULTIPLY,t2term,to_unwrap,t1,result,swap);
        break;

      default:
        okay=false;
        break;
    }
    if(!okay){
        // cannot invert this one, finish here, either by giving up or going to the end
        if(!modified) return false;
        goto endOfUnwrapping; 
    }

    // update t1
    t1=result;
    // set t2 to the non-constant argument
    t2 = to_unwrap;
    modified = true;
  }
endOfUnwrapping:

  //Evaluate t1
  // We have rearranged things so that t2 is a non-constant term and t1 is a number
  // of interprted functions applied to some constants. By evaluating t1 we will
  //  get a constant (unless evaluation is not possible)

  // don't swap equality
  if(lit->functor()==0){
   resLit = BottomUpTermTransformer::transformLiteral(Literal::createEquality(lit->polarity(),t2,t1,srt));
  }
  else{
    // important, need to preserve the ordering of t1 and t2 in the original!
    if(swap){
      resLit = BottomUpTermTransformer::transformLiteral(Literal::create2(lit->functor(),lit->polarity(),t2,t1));
    }else{
      resLit = BottomUpTermTransformer::transformLiteral(Literal::create2(lit->functor(),lit->polarity(),t1,t2));
    }
  }
  return true;
}


bool InterpretedLiteralEvaluator::balancePlus(Interpretation plus, Interpretation unaryMinus, 
                                              Term* AplusB, TermList A, TermList C, TermList& result)
{
    unsigned um = env.signature->getInterpretingSymbol(unaryMinus);
    unsigned ip = env.signature->getInterpretingSymbol(plus);
    TermList B;
    if(AplusB->termArg(0)==A){
      B = AplusB->termArg(1);
    }
    else{
      ASS(AplusB->termArg(1)==A);
      B = AplusB->termArg(0);
    }

    TermList mB(Term::create1(um,B));
    result = TermList(Term::create2(ip,C,mB));
    return true;
}

template<typename ConstantType>
bool InterpretedLiteralEvaluator::balanceMultiply(Interpretation divide,ConstantType zero, 
                                                  Term* AmultiplyB, TermList A, TermList C, TermList& result,
                                                  bool& swap)
{
#if VDEBUG
    TermList srt = theory->getOperationSort(divide); 
    ASS(srt == AtomicSort::realSort() || srt == AtomicSort::rationalSort()); 
#endif

    unsigned div = env.signature->getInterpretingSymbol(divide);
    TermList B;
    if(AmultiplyB->termArg(0)==A){
      B = AmultiplyB->termArg(1);
    }
    else{
      ASS(AmultiplyB->termArg(1)==A);
      B = AmultiplyB->termArg(0);
    }
    result = TermList(Term::create2(div,C,B));

    ConstantType bcon;
    if(theory->tryInterpretConstant(B,bcon)){
      if(bcon.isZero()) return false;
      if(bcon.isNegative()){ swap=!swap; } // switch the polarity of an inequality if we're under one
      return true;
    }
    // Unsure exactly what the best thing to do here, so for now give up
    // This means we only balance when we have a constant on the variable side
    return false;
}

bool InterpretedLiteralEvaluator::balanceIntegerMultiply(
                                                  Term* AmultiplyB, TermList A, TermList C, TermList& result,
                                                  bool& swap)
{
    // only works if we in the end divid a number by a number
    IntegerConstantType ccon;
    if(!theory->tryInterpretConstant(C,ccon)){ return false; }

    // we are going to use rounding division but ensure that it is non-rounding
    unsigned div = env.signature->getInterpretingSymbol(Theory::INT_QUOTIENT_E);
    TermList B;
    if(AmultiplyB->termArg(0)==A){
      B = AmultiplyB->termArg(1);
    }
    else{
      ASS(AmultiplyB->termArg(1)==A);
      B = AmultiplyB->termArg(0);
    }
    result = TermList(Term::create2(div,C,B));

    IntegerConstantType bcon;
    if(theory->tryInterpretConstant(B,bcon)){
      if(bcon.isZero()){ return false; }
      if(!bcon.divides(ccon)){ return false;}
      if(bcon.isNegative()){ swap=!swap; } // switch the polarity of an inequality if we're under one
      return true;
    }
    return false;
}

bool InterpretedLiteralEvaluator::balanceDivide(Interpretation multiply, 
                       Term* AoverB, TermList A, TermList C, TermList& result, bool& swap)
{
#if VDEBUG
    TermList srt = theory->getOperationSort(multiply); 
    ASS(srt == AtomicSort::realSort() || srt == AtomicSort::rationalSort());
#endif

    unsigned mul = env.signature->getInterpretingSymbol(multiply);
    if(AoverB->termArg(0)!=A)return false;

    auto B = AoverB->termArg(1);

    result = TermList(Term::create2(mul,C,B));

    RationalConstantType bcon;
    if(theory->tryInterpretConstant(B,bcon)){
      ASS(!bcon.isZero());
      if(bcon.isNegative()){ swap=!swap; } // switch the polarity of an inequality if we're under one
      return true;
    }
    // Unsure exactly what the best thing to do here, so for now give up
    // This means we only balance when we have a constant on the variable side
    return false;    
}

/*
 // TODO document me
 Literal& balance(Literal& in)
 {
   // we only rebalance equalities
   if (!in.isEquality()) {
     return in;

   } else {
     ASS(in.numTermArguments() == 2);
     unsigned sort;
     if (!SortHelper::tryGetResultSort(in[0], sort) &&
         !SortHelper::tryGetResultSort(in[1], sort)) {
       return in;
     } else {

       Literal* out;
       switch (sort) {
 #define _CASE(SRT, ConstantType) \
         case Sorts::SRT: \
           if (!balance<ConstantType>(in, out)){ \
             return in; \
           }
         _CASE(SRT_REAL    ,    RealConstantType)
         _CASE(SRT_INTEGER , IntegerConstantType)
         _CASE(SRT_RATIONAL,RationalConstantType)
 #undef _CASE
         default: return in;
       }
       return *out;
     }
   }
 }
*/

class LiteralNormalizer 
{
  // static bool canNormalizePred(Theory::Interpretation pred);
  // static bool canNormalize(const Literal& l)
  // {
  //   return canNormalize(env.signature->getPredicate(l.getPredicate))
  // }

public:
  static Literal* normalize(Literal* in) {
    auto functor = in->functor();

    if (theory->isInterpretedPredicate(functor)) {
      auto i = theory->interpretPredicate(functor);

      switch (i) {
        case Interpretation::INT_LESS: 
          return InequalityNormalizer<IntLess>::normalize(in);

        case Interpretation::RAT_LESS: 
          return InequalityNormalizer<FracLess<RationalConstantType>>::normalize(in);

        case Interpretation::REAL_LESS: 
          return InequalityNormalizer<FracLess<RealConstantType>>::normalize(in);

        default: 
          return in;
      }
    } else {
      /* uninterpreted predicate */
      return in;
    }
  }
};
TermList InterpretedLiteralEvaluator::evaluate(TermList t) {
  if (t.isTerm())
    t = TermList(BottomUpTermTransformer::transform(t.term()));
  return InterpretedLiteralEvaluator::transformSubterm(t);
}

/**
 * Used to evaluate a literal, setting isConstant, resLit and resConst in the process
 *
 * Returns true if it has been evaluated, in which case resLit is set 
 * isConstant is true if the literal predicate evaluates to a constant value
 * resConst is set iff isConstant and gives the constant value (true/false) of resLit 
 */
bool InterpretedLiteralEvaluator::evaluate(Literal* lit, bool& isConstant, Literal*& resLit, bool& resConst)
{
  DEBUG( "evaluate ", lit->toString() );

  // This tries to transform each subterm using tryEvaluateFunc (see transform Subterm below)

  resLit = _normalize ? LiteralNormalizer::normalize(lit)
                      : lit;
  DEBUG( "\t0 ==> ", resLit->toString() );

  resLit = BottomUpTermTransformer::transformLiteral( resLit);
  DEBUG( "\t1 ==> ", resLit->toString() );

  unsigned pred = resLit->functor();

  // Now we try and evaluate the predicate
  Evaluator* predEv = getPredEvaluator(pred);
  if (predEv) {
    auto r = predEv->tryEvaluatePred(resLit);

    switch (r.status) {
      case PredEvalResult::Nop: 
        break;

      case PredEvalResult::Simplified: 
        // resLit = r.simplified_val;
        resLit = BottomUpTermTransformer::transformLiteral(r.simplified_val);
        break;

      case PredEvalResult::Trivial: 
        isConstant = true;
        resConst = r.trivial_val;
        DEBUG( "\t3 ==> ", resConst );
        return true;

    }
  }
  isConstant = false;
  auto out = resLit != lit;
  DEBUG( "\t3 ==> ", resLit->toString(), "(did evaluate: ", out, ")" );
  return out;
}




/**
 * This attempts to evaluate each subterm.
 * See Kernel/BottomUpTermTransformer for how it is used.
 * Terms are evaluated bottom-up
 */
TermList InterpretedLiteralEvaluator::transformSubterm(TermList trm)
{
  // Debug::Tracer::printStack(cout);

  // DEBUG( "transformSubterm for ", trm.toString() );

  //Nothing to evaluate in a sort
  if (!trm.isTerm() || trm.term()->isSort()) { return trm; }
  Term* t = trm.term();
  unsigned func = t->functor();

  Evaluator* funcEv = getFuncEvaluator(func);
  if (funcEv) {
    TermList res;
    if (funcEv->tryEvaluateFunc(t, res)) {
	return res;
    }
  } else {
    DEBUG("no transformer")
  }
  return trm;
}

/**
 * This searches for an Evaluator for a function
 */
template<class Fn>
InterpretedLiteralEvaluator::Evaluator* InterpretedLiteralEvaluator::getEvaluator(unsigned functor, DArray<Evaluator*>& evaluators, Fn canEval) 
{
  if (functor >= evaluators.size()) {
    unsigned oldSz = evaluators.size();
    unsigned newSz = functor + 1 ;
    evaluators.expand(newSz);
    for (unsigned i=oldSz; i<newSz; i++) {
        EvalStack::Iterator evit(_evals);
	while (evit.hasNext()) {
          Evaluator* ev = evit.next();
          // we only set the evaluator if it has not yet been set
	  if (canEval(ev, i)) {
	    evaluators[i] = ev;
            goto break_inner;
	  }
	}
break_inner:;
    }
  }
  return evaluators[functor];
}

/**
 * This searches for an Evaluator for a function
 */
InterpretedLiteralEvaluator::Evaluator* InterpretedLiteralEvaluator::getFuncEvaluator(unsigned func)
{
   return getEvaluator(func, 
       this->_funEvaluators, 
       [] (Evaluator* ev, unsigned i) {return ev->canEvaluateFunc(i); });
}

/**
 * This searches for an Evaluator for a predicate
 */
InterpretedLiteralEvaluator::Evaluator* InterpretedLiteralEvaluator::getPredEvaluator(unsigned pred)
{
   return getEvaluator(pred, this->_predEvaluators,
       [] (Evaluator* ev, unsigned i) {return ev->canEvaluatePred(i); });
}

}
