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
 * @file Theory.hpp
 * Defines class Theory.
 */

#ifndef __Theory__
#define __Theory__

#include <cmath>

#include "Forwards.hpp"

#include "Lib/DHMap.hpp"
#include "Lib/Exception.hpp"

#include "Shell/TermAlgebra.hpp"

#include "OperatorType.hpp"
#include "Term.hpp"

namespace Kernel {

class IntegerConstantType;
struct RationalConstantType;
class RealConstantType;

/**
 * Exception to be thrown when the requested operation cannot be performed,
 * e.g. because of overflow of a native type.
 */
class ArithmeticException : public Exception {
protected:
  ArithmeticException(const char* msg) : Exception(msg) {}
};

class MachineArithmeticException : public ArithmeticException 
{ 
public:
  MachineArithmeticException() : ArithmeticException("machine arithmetic exception"){} 
};

class DivByZeroException         : public ArithmeticException 
{ 
public:
  DivByZeroException() : ArithmeticException("divided by zero"){} 
};

enum class Sign : uint8_t {
  Zero = 0,
  Pos = 1,
  Neg = 2,
};

inline std::ostream& operator<<(std::ostream& out, Sign const& self)
{
  switch(self) {
    case Sign::Zero: return out << "0";
    case Sign::Pos: return out << "+";
    case Sign::Neg: return out << "-";
  }
  ASSERTION_VIOLATION
}

class IntegerConstantType
{
public:
  static TermList getSort() { return AtomicSort::intSort(); }

  typedef int InnerType;

  IntegerConstantType() {}
  IntegerConstantType(IntegerConstantType&&) = default;
  IntegerConstantType(const IntegerConstantType&) = default;
  IntegerConstantType& operator=(const IntegerConstantType&) = default;
  constexpr IntegerConstantType(InnerType v) : _val(v) {}
  explicit IntegerConstantType(const std::string& str);

  IntegerConstantType operator+(const IntegerConstantType& num) const;
  IntegerConstantType operator-(const IntegerConstantType& num) const;
  IntegerConstantType operator-() const;
  IntegerConstantType operator*(const IntegerConstantType& num) const;

  // true if this divides num
  bool divides(const IntegerConstantType& num) const ;
  float realDivide(const IntegerConstantType& num) const { 
    if(num._val==0) throw DivByZeroException();
    return ((float)_val)/num._val; 
  }
  int intDivide(const IntegerConstantType& num) const ;  
  IntegerConstantType remainderE(const IntegerConstantType& num) const; 
  IntegerConstantType quotientE(const IntegerConstantType& num) const; 
  IntegerConstantType quotientT(const IntegerConstantType& num) const;
  IntegerConstantType quotientF(const IntegerConstantType& num) const; 

  IntegerConstantType remainderT(const IntegerConstantType& num) const
  { return (*this) - num * quotientT(num); }
  IntegerConstantType remainderF(const IntegerConstantType& num) const
  { return (*this) - num * quotientF(num); } 

  bool operator==(const IntegerConstantType& num) const;
  bool operator>(const IntegerConstantType& num) const;

  bool operator!=(const IntegerConstantType& num) const { return !((*this)==num); }
  bool operator<(const IntegerConstantType& o) const { return o>(*this); }
  bool operator>=(const IntegerConstantType& o) const { return !(o>(*this)); }
  bool operator<=(const IntegerConstantType& o) const { return !((*this)>o); }

  InnerType toInner() const { return _val; }

  Sign sign() const 
  { return _val > 0 ? Sign::Pos
     : _val < 0 ? Sign::Neg
     :            Sign::Zero; }

  bool isZero(){ return _val==0; }
  bool isNegative(){ return _val<0; }

  static IntegerConstantType floor(RationalConstantType rat);
  static IntegerConstantType floor(IntegerConstantType rat);

  static IntegerConstantType ceiling(RationalConstantType rat);
  static IntegerConstantType ceiling(IntegerConstantType rat);
  IntegerConstantType abs() const;

  static Comparison comparePrecedence(IntegerConstantType n1, IntegerConstantType n2);
  size_t hash() const;

  std::string toString() const;
private:
  InnerType _val;
  IntegerConstantType operator/(const IntegerConstantType& num) const;
  IntegerConstantType operator%(const IntegerConstantType& num) const;
};

inline
std::ostream& operator<< (std::ostream& out, const IntegerConstantType& val) {
  return out << val.toInner();
}

/**
 * A class for representing rational numbers
 *
 * The class uses IntegerConstantType to store the numerator and denominator.
 * This ensures that if there is an overflow in the operations, an exception
 * will be raised by the IntegerConstantType methods.
 */
struct RationalConstantType {
  typedef IntegerConstantType InnerType;

  static TermList getSort() { return AtomicSort::rationalSort(); }

  RationalConstantType() {}
  RationalConstantType(RationalConstantType&&) = default;
  RationalConstantType(const RationalConstantType&) = default;
  RationalConstantType& operator=(const RationalConstantType&) = default;

  RationalConstantType(InnerType num, InnerType den);
  RationalConstantType(const std::string& num, const std::string& den);
  constexpr RationalConstantType(InnerType num) : _num(num), _den(1) {} //assuming den=1

  RationalConstantType operator+(const RationalConstantType& num) const;
  RationalConstantType operator-(const RationalConstantType& num) const;
  RationalConstantType operator-() const;
  RationalConstantType operator*(const RationalConstantType& num) const;
  RationalConstantType operator/(const RationalConstantType& num) const;

  RationalConstantType floor() const { 
    return RationalConstantType(IntegerConstantType::floor(*this));
  }
  RationalConstantType ceiling() const { 
    return RationalConstantType(IntegerConstantType::ceiling(*this));
  }
  RationalConstantType truncate() const { 
    return RationalConstantType(_num.quotientT(_den));
  }

  bool isInt() const;

  bool operator==(const RationalConstantType& num) const;
  bool operator>(const RationalConstantType& num) const;

  bool operator!=(const RationalConstantType& num) const { return !((*this)==num); }
  bool operator<(const RationalConstantType& o) const { return o>(*this); }
  bool operator>=(const RationalConstantType& o) const { return !(o>(*this)); }
  bool operator<=(const RationalConstantType& o) const { return !((*this)>o); }

  bool isZero(){ return _num.toInner()==0; } 
  // relies on the fact that cannonize ensures that _den>=0
  bool isNegative() const { ASS(_den>=0); return _num.toInner() < 0; }
  bool isPositive() const { ASS(_den>=0); return _num.toInner() > 0; }

  RationalConstantType abs() const;

  std::string toString() const;

  const InnerType& numerator() const { return _num; }
  const InnerType& denominator() const { return _den; }
  size_t hash() const;

  static Comparison comparePrecedence(RationalConstantType n1, RationalConstantType n2);

protected:
  void init(InnerType num, InnerType den);

private:
  void cannonize();

  InnerType _num;
  InnerType _den;
};

inline
std::ostream& operator<< (std::ostream& out, const RationalConstantType& val) {
  return out << val.toString();
}


class RealConstantType : public RationalConstantType
{
public:
  static TermList getSort() { return AtomicSort::realSort(); }

  RealConstantType() {}
  RealConstantType(RealConstantType&&) = default;
  RealConstantType(const RealConstantType&) = default;
  RealConstantType& operator=(const RealConstantType&) = default;

  explicit RealConstantType(const std::string& number);
  explicit constexpr RealConstantType(const RationalConstantType& rat) : RationalConstantType(rat) {}
  RealConstantType(int num, int den) : RationalConstantType(num, den) {}
  explicit constexpr RealConstantType(typename IntegerConstantType::InnerType number) : RealConstantType(RationalConstantType(number)) {}

  RealConstantType operator+(const RealConstantType& num) const
  { return RealConstantType(RationalConstantType::operator+(num)); }
  RealConstantType operator-(const RealConstantType& num) const
  { return RealConstantType(RationalConstantType::operator-(num)); }
  RealConstantType operator-() const
  { return RealConstantType(RationalConstantType::operator-()); }
  RealConstantType operator*(const RealConstantType& num) const
  { return RealConstantType(RationalConstantType::operator*(num)); }
  RealConstantType operator/(const RealConstantType& num) const
  { return RealConstantType(RationalConstantType::operator/(num)); }

  RealConstantType floor() const { return RealConstantType(RationalConstantType::floor()); }
  RealConstantType truncate() const { return RealConstantType(RationalConstantType::truncate()); }
  RealConstantType ceiling() const { return RealConstantType(RationalConstantType::ceiling()); }

  RealConstantType abs() const;

  std::string toNiceString() const;

  size_t hash() const;
  static Comparison comparePrecedence(RealConstantType n1, RealConstantType n2);

  /** 
   * returns the internal represenation of this RealConstantType. 
   * 
   * Currently we represent Reals as Rationals. We might
   * change this representation in the future in order to represent numerals other algebraic numbers (e.g.  sqrt(2)). 
   * In order to make this future proof this function is called in places where we rely on the representation of reals,
   * so we get a compiler error if we change the underlying datatype.
   */
  RationalConstantType representation() const;
private:
  static bool parseDouble(const std::string& num, RationalConstantType& res);

};

inline bool operator<(const RealConstantType& lhs ,const RealConstantType& rhs) { 
  return static_cast<const RationalConstantType&>(lhs) < static_cast<const RationalConstantType&>(rhs);
}
inline bool operator>(const RealConstantType& lhs, const RealConstantType& rhs) {
  return static_cast<const RationalConstantType&>(lhs) > static_cast<const RationalConstantType&>(rhs);
}
inline bool operator<=(const RealConstantType& lhs, const RealConstantType& rhs) {
  return static_cast<const RationalConstantType&>(lhs) <= static_cast<const RationalConstantType&>(rhs);
}
inline bool operator>=(const RealConstantType& lhs, const RealConstantType& rhs) {
  return static_cast<const RationalConstantType&>(lhs) >= static_cast<const RationalConstantType&>(rhs);
}

inline
std::ostream& operator<< (std::ostream& out, const RealConstantType& val) {
  return out << val.toString();
}


/**
 * A singleton class handling tasks related to theory symbols in Vampire
 */
class Theory
{
public:
  /**
   * Interpreted symbols and predicates
   *
   * If interpreted_evaluation is enabled, predicates GREATER_EQUAL,
   * LESS and LESS_EQUAL should not appear in the run of the
   * SaturationAlgorithm (they'll be immediately simplified by the
   * InterpretedEvaluation simplification).
   */
  enum Interpretation
  {
    //predicates
    EQUAL,

    INT_IS_INT,
    INT_IS_RAT,
    INT_IS_REAL,
    INT_GREATER,
    INT_GREATER_EQUAL,
    INT_LESS,
    INT_LESS_EQUAL,
    INT_DIVIDES,

    RAT_IS_INT,
    RAT_IS_RAT,
    RAT_IS_REAL,
    RAT_GREATER,
    RAT_GREATER_EQUAL,
    RAT_LESS,
    RAT_LESS_EQUAL,

    REAL_IS_INT,
    REAL_IS_RAT,
    REAL_IS_REAL,
    REAL_GREATER,
    REAL_GREATER_EQUAL,
    REAL_LESS,
    REAL_LESS_EQUAL,

    //numeric functions

    INT_SUCCESSOR,
    INT_UNARY_MINUS,
    INT_PLUS,  // sum in TPTP
    INT_MINUS, // difference in TPTP
    INT_MULTIPLY,
    INT_QUOTIENT_E,
    INT_QUOTIENT_T,
    INT_QUOTIENT_F,
    INT_REMAINDER_E,
    INT_REMAINDER_T,
    INT_REMAINDER_F,
    INT_FLOOR,
    INT_CEILING,
    INT_TRUNCATE,
    INT_ROUND,
    INT_ABS,

    RAT_UNARY_MINUS,
    RAT_PLUS, // sum in TPTP
    RAT_MINUS,// difference in TPTP
    RAT_MULTIPLY,
    RAT_QUOTIENT,
    RAT_QUOTIENT_E,
    RAT_QUOTIENT_T,
    RAT_QUOTIENT_F,
    RAT_REMAINDER_E,
    RAT_REMAINDER_T,
    RAT_REMAINDER_F,
    RAT_FLOOR,
    RAT_CEILING,
    RAT_TRUNCATE,
    RAT_ROUND,

    REAL_UNARY_MINUS,
    REAL_PLUS,  // plus in TPTP
    REAL_MINUS, // difference in TPTP
    REAL_MULTIPLY,
    REAL_QUOTIENT,
    REAL_QUOTIENT_E,
    REAL_QUOTIENT_T,
    REAL_QUOTIENT_F,
    REAL_REMAINDER_E,
    REAL_REMAINDER_T,
    REAL_REMAINDER_F,
    REAL_FLOOR,
    REAL_CEILING,
    REAL_TRUNCATE,
    REAL_ROUND,

    //conversion functions
    INT_TO_INT,
    INT_TO_RAT,
    INT_TO_REAL,
    RAT_TO_INT,
    RAT_TO_RAT,
    RAT_TO_REAL,
    REAL_TO_INT,
    REAL_TO_RAT,
    REAL_TO_REAL,

    // array functions
    ARRAY_SELECT,
    ARRAY_BOOL_SELECT,
    ARRAY_STORE,

    INVALID_INTERPRETATION // LEAVE THIS AS THE LAST ELEMENT OF THE ENUM
  };

  enum IndexedInterpretation {
    FOR_NOW_EMPTY
  };

  typedef std::pair<IndexedInterpretation, unsigned> ConcreteIndexedInterpretation;

  /**
   * Interpretations represent the abstract concept of an interpreted operation vampire knows about.
   *
   * Some of them are polymorphic, such the ones for ARRAYs, and only become a concrete operation when supplied with OperatorType*.
   * To identify these, MonomorphisedInterpretation (see below) can be used. However, notice that the appropriate Symbol always carries
   * an Interpretation (if interpreted) and an OperatorType*.
   *
   * Other operations might be, in fact, indexed families of operations, and need an additional index (unsigned) to be specified.
   * To keep the Symbol structure from growing for their sake, these IndexedInterpretations are instantiated to a concrete member of a family on demand
   * and we keep track of this instantiation in _indexedInterpretations (see below). Members of _indexedInterpretations
   * implicitly "inhabit" a range of values in Interpretation after INVALID_INTERPRETATION, so that again an
   * Interpretation (this time >= INVALID_INTERPRETATION) and an OperatorType* are enough to identify a member of an indexed family.
   */

  typedef std::pair<Interpretation, OperatorType*> MonomorphisedInterpretation;

private:
  DHMap<ConcreteIndexedInterpretation,Interpretation> _indexedInterpretations;

public:

  static unsigned numberOfFixedInterpretations() {
    return INVALID_INTERPRETATION;
  }

  Interpretation interpretationFromIndexedInterpretation(IndexedInterpretation ii, unsigned index)
  {
    ConcreteIndexedInterpretation cii = std::make_pair(ii,index);

    Interpretation res;
    if (!_indexedInterpretations.find(cii, res)) {
      res = static_cast<Interpretation>(numberOfFixedInterpretations() + _indexedInterpretations.size());
      _indexedInterpretations.insert(cii, res);
    }
    return res;
  }

  static bool isPlus(Interpretation i){
    return i == INT_PLUS || i == RAT_PLUS || i == REAL_PLUS;
  }

  static std::string getInterpretationName(Interpretation i);
  static unsigned getArity(Interpretation i);
  static bool isFunction(Interpretation i);
  static bool isInequality(Interpretation i);
  static OperatorType* getNonpolymorphicOperatorType(Interpretation i);

  static OperatorType* getArrayOperatorType(TermList arraySort, Interpretation i);

  static bool hasSingleSort(Interpretation i);
  static TermList getOperationSort(Interpretation i);
  static bool isConversionOperation(Interpretation i);
  static bool isLinearOperation(Interpretation i);
  static bool isNonLinearOperation(Interpretation i);
  bool isPartiallyInterpretedFunction(Term* t);
  bool partiallyDefinedFunctionUndefinedForArgs(Term* t);
  // static bool isPartialFunction(Interpretation i);

  static bool isPolymorphic(Interpretation i);

  unsigned getArrayExtSkolemFunction(TermList sort);

  static Theory theory_obj;
  static Theory* instance();

  void defineTupleTermAlgebra(unsigned arity, TermList* sorts);

  /** Returns true if the argument is an interpreted constant
   */
  bool isInterpretedConstant(unsigned func);
  bool isInterpretedConstant(Term* t);
  bool isInterpretedConstant(TermList t);
  /** Returns true if the argument is an interpreted number
   */
  bool isInterpretedNumber(Term* t);
  bool isInterpretedNumber(TermList t);

  /** Returns false if pred is equality.  Returns true if the argument is any
      other interpreted predicate.
   */
  bool isInterpretedPredicate(unsigned pred);

  bool isInterpretedEquality(Literal* lit);
  bool isInterpretedPredicate(Literal* lit, Interpretation itp);
  bool isInterpretedPredicate(unsigned f, Interpretation itp);
  bool isInterpretedPredicate(Literal* lit);

  bool isInterpretedFunction(unsigned func);
  bool isInterpretedFunction(Term* t);
  bool isInterpretedFunction(TermList t);
  bool isInterpretedFunction(unsigned func, Interpretation itp);
  bool isInterpretedFunction(Term* t, Interpretation itp);
  bool isInterpretedFunction(TermList t, Interpretation itp);

  // bool isInterpretedPartialFunction(unsigned func);
  bool isZero(TermList t);

  Interpretation interpretFunction(unsigned func);
  Interpretation interpretFunction(Term* t);
  Interpretation interpretFunction(TermList t);
  Interpretation interpretPredicate(unsigned pred);
  Interpretation interpretPredicate(Literal* t);

  void registerLaTeXPredName(unsigned func, bool polarity, std::string temp);
  void registerLaTeXFuncName(unsigned func, std::string temp);
  std::string tryGetInterpretedLaTeXName(unsigned func, bool pred,bool polarity=true);

private:
  // For recording the templates for predicate and function symbols
  DHMap<unsigned,std::string> _predLaTeXnamesPos;
  DHMap<unsigned,std::string> _predLaTeXnamesNeg;
  DHMap<unsigned,std::string> _funcLaTeXnames;

public:

  /**
   * Try to interpret the term list as an integer constant. If it is an
   * integer constant, return true and save the constant in @c res, otherwise
   * return false.
   */
  bool tryInterpretConstant(TermList trm, IntegerConstantType& res)
  {
    if (!trm.isTerm()) {
      return false;
    }
    return tryInterpretConstant(trm.term(),res);
  }
  bool tryInterpretConstant(const Term* t, IntegerConstantType& res);
  bool tryInterpretConstant(unsigned functor, IntegerConstantType& res);
  /**
   * Try to interpret the term list as an rational constant. If it is an
   * rational constant, return true and save the constant in @c res, otherwise
   * return false.
   */
  bool tryInterpretConstant(TermList trm, RationalConstantType& res)
  {
    if (!trm.isTerm()) {
      return false;
    }
    return tryInterpretConstant(trm.term(),res);
  }
  bool tryInterpretConstant(const Term* t, RationalConstantType& res);
  bool tryInterpretConstant(unsigned functor, RationalConstantType& res);
  /**
   * Try to interpret the term list as an real constant. If it is an
   * real constant, return true and save the constant in @c res, otherwise
   * return false.
   */
  bool tryInterpretConstant(TermList trm, RealConstantType& res)
  {
    if (!trm.isTerm()) {
      return false;
    }
    return tryInterpretConstant(trm.term(),res);
  }
  bool tryInterpretConstant(const Term* t, RealConstantType& res);
  bool tryInterpretConstant(unsigned functor, RealConstantType& res);

  Term* representConstant(const IntegerConstantType& num);
  Term* representConstant(const RationalConstantType& num);
  Term* representConstant(const RealConstantType& num);

  Term* representIntegerConstant(std::string str);
  Term* representRealConstant(std::string str);
private:
  Theory();
  static OperatorType* getConversionOperationType(Interpretation i);

  DHMap<TermList,unsigned> _arraySkolemFunctions;

public:
  class Tuples {
  public:
    bool isFunctor(unsigned functor);
    unsigned getFunctor(unsigned arity, TermList sorts[]);
    unsigned getFunctor(TermList tupleSort);
    unsigned getProjectionFunctor(unsigned proj, TermList tupleSort);
    bool findProjection(unsigned projFunctor, bool isPredicate, unsigned &proj);
  };

  static Theory::Tuples tuples_obj;
  static Theory::Tuples* tuples();
};

#define ANY_INTERPRETED_PREDICATE                                                                             \
         Kernel::Theory::EQUAL:                                                                               \
    case Kernel::Theory::INT_IS_INT:                                                                          \
    case Kernel::Theory::INT_IS_RAT:                                                                          \
    case Kernel::Theory::INT_IS_REAL:                                                                         \
    case Kernel::Theory::INT_GREATER:                                                                         \
    case Kernel::Theory::INT_GREATER_EQUAL:                                                                   \
    case Kernel::Theory::INT_LESS:                                                                            \
    case Kernel::Theory::INT_LESS_EQUAL:                                                                      \
    case Kernel::Theory::INT_DIVIDES:                                                                         \
    case Kernel::Theory::RAT_IS_INT:                                                                          \
    case Kernel::Theory::RAT_IS_RAT:                                                                          \
    case Kernel::Theory::RAT_IS_REAL:                                                                         \
    case Kernel::Theory::RAT_GREATER:                                                                         \
    case Kernel::Theory::RAT_GREATER_EQUAL:                                                                   \
    case Kernel::Theory::RAT_LESS:                                                                            \
    case Kernel::Theory::RAT_LESS_EQUAL:                                                                      \
    case Kernel::Theory::REAL_IS_INT:                                                                         \
    case Kernel::Theory::REAL_IS_RAT:                                                                         \
    case Kernel::Theory::REAL_IS_REAL:                                                                        \
    case Kernel::Theory::REAL_GREATER:                                                                        \
    case Kernel::Theory::REAL_GREATER_EQUAL:                                                                  \
    case Kernel::Theory::REAL_LESS:                                                                           \
    case Kernel::Theory::ARRAY_BOOL_SELECT:                                                                   \
    case Kernel::Theory::REAL_LESS_EQUAL

#define ANY_INTERPRETED_FUNCTION                                                                              \
         Kernel::Theory::INT_SUCCESSOR:                                                                       \
    case Kernel::Theory::INT_UNARY_MINUS:                                                                     \
    case Kernel::Theory::INT_PLUS:                                                                            \
    case Kernel::Theory::INT_MINUS:                                                                           \
    case Kernel::Theory::INT_MULTIPLY:                                                                        \
    case Kernel::Theory::INT_QUOTIENT_E:                                                                      \
    case Kernel::Theory::INT_QUOTIENT_T:                                                                      \
    case Kernel::Theory::INT_QUOTIENT_F:                                                                      \
    case Kernel::Theory::INT_REMAINDER_E:                                                                     \
    case Kernel::Theory::INT_REMAINDER_T:                                                                     \
    case Kernel::Theory::INT_REMAINDER_F:                                                                     \
    case Kernel::Theory::INT_FLOOR:                                                                           \
    case Kernel::Theory::INT_CEILING:                                                                         \
    case Kernel::Theory::INT_TRUNCATE:                                                                        \
    case Kernel::Theory::INT_ROUND:                                                                           \
    case Kernel::Theory::INT_ABS:                                                                             \
    case Kernel::Theory::RAT_UNARY_MINUS:                                                                     \
    case Kernel::Theory::RAT_PLUS:                                                                            \
    case Kernel::Theory::RAT_MINUS:                                                                           \
    case Kernel::Theory::RAT_MULTIPLY:                                                                        \
    case Kernel::Theory::RAT_QUOTIENT:                                                                        \
    case Kernel::Theory::RAT_QUOTIENT_E:                                                                      \
    case Kernel::Theory::RAT_QUOTIENT_T:                                                                      \
    case Kernel::Theory::RAT_QUOTIENT_F:                                                                      \
    case Kernel::Theory::RAT_REMAINDER_E:                                                                     \
    case Kernel::Theory::RAT_REMAINDER_T:                                                                     \
    case Kernel::Theory::RAT_REMAINDER_F:                                                                     \
    case Kernel::Theory::RAT_FLOOR:                                                                           \
    case Kernel::Theory::RAT_CEILING:                                                                         \
    case Kernel::Theory::RAT_TRUNCATE:                                                                        \
    case Kernel::Theory::RAT_ROUND:                                                                           \
    case Kernel::Theory::REAL_UNARY_MINUS:                                                                    \
    case Kernel::Theory::REAL_PLUS:                                                                           \
    case Kernel::Theory::REAL_MINUS:                                                                          \
    case Kernel::Theory::REAL_MULTIPLY:                                                                       \
    case Kernel::Theory::REAL_QUOTIENT:                                                                       \
    case Kernel::Theory::REAL_QUOTIENT_E:                                                                     \
    case Kernel::Theory::REAL_QUOTIENT_T:                                                                     \
    case Kernel::Theory::REAL_QUOTIENT_F:                                                                     \
    case Kernel::Theory::REAL_REMAINDER_E:                                                                    \
    case Kernel::Theory::REAL_REMAINDER_T:                                                                    \
    case Kernel::Theory::REAL_REMAINDER_F:                                                                    \
    case Kernel::Theory::REAL_FLOOR:                                                                          \
    case Kernel::Theory::REAL_CEILING:                                                                        \
    case Kernel::Theory::REAL_TRUNCATE:                                                                       \
    case Kernel::Theory::REAL_ROUND:                                                                          \
    case Kernel::Theory::INT_TO_INT:                                                                          \
    case Kernel::Theory::INT_TO_RAT:                                                                          \
    case Kernel::Theory::INT_TO_REAL:                                                                         \
    case Kernel::Theory::RAT_TO_INT:                                                                          \
    case Kernel::Theory::RAT_TO_RAT:                                                                          \
    case Kernel::Theory::RAT_TO_REAL:                                                                         \
    case Kernel::Theory::REAL_TO_INT:                                                                         \
    case Kernel::Theory::REAL_TO_RAT:                                                                         \
    case Kernel::Theory::REAL_TO_REAL:                                                                        \
    case Kernel::Theory::ARRAY_SELECT:                                                                        \
    case Kernel::Theory::ARRAY_STORE

typedef Theory::Interpretation Interpretation;

/**
 * Pointer to the singleton Theory instance
 */
extern Theory* theory;

std::ostream& operator<<(std::ostream& out, Kernel::Theory::Interpretation const& self);

}

template<>
struct std::hash<Kernel::IntegerConstantType>
{
  size_t operator()(Kernel::IntegerConstantType const& self) const noexcept 
  { return self.hash(); }
};

template<>
struct std::hash<Kernel::RationalConstantType>
{
  size_t operator()(Kernel::RationalConstantType const& self) const noexcept 
  { return self.hash(); }
};


template<>
struct std::hash<Kernel::RealConstantType>
{
  size_t operator()(Kernel::RealConstantType const& self) const noexcept 
  { return self.hash(); }
};

#endif // __Theory__
