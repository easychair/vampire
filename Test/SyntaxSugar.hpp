/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
/**! This file contains macros to provide syntax sugar for building formulas,
 * terms, etc. for test cases.
 *
 * Macros that are not meant to be used from outside of this file are prefixed 
 * with two underscores.
 *
 * @autor Johannes Schoisswohl
 * @date 2020-04-29
 */

#ifndef __TEST__SYNTAX_SUGAR__H__
#define __TEST__SYNTAX_SUGAR__H__

#include<functional>

#include "Forwards.hpp"
#include "Lib/Environment.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/NumTraits.hpp"
#include "Kernel/ApplicativeHelper.hpp"

#include "Indexing/TermSharing.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/TermIterators.hpp"
#include "Kernel/OperatorType.hpp"
#include "Shell/TermAlgebra.hpp"
#include "Shell/FunctionDefinitionHandler.hpp"

#define __TO_SORT_RAT RationalConstantType::getSort()
#define __TO_SORT_INT IntegerConstantType::getSort()
#define __TO_SORT_REAL RealConstantType::getSort()

#define __CONSTANT_TYPE_INT  IntegerConstantType
#define __CONSTANT_TYPE_REAL RealConstantType
#define __CONSTANT_TYPE_RAT  RationalConstantType


#define __ARGS_DECL(Type, arity) __ARGS_DECL_ ## arity(Type)
#define __ARGS_DECL_1(Type) Type arg0_ 
#define __ARGS_DECL_2(Type) Type arg0_ , Type arg1_ 

#define __ARGS_EXPR(Type, arity) __ARGS_EXPR_ ## arity(Type)
#define __ARGS_EXPR_1(Type) arg0_
#define __ARGS_EXPR_2(Type) arg0_, arg1_

#define __REPEAT_1(sort) sort
#define __REPEAT_2(sort) sort, __REPEAT_1(sort)
#define __REPEAT_3(sort) sort, __REPEAT_2(sort)
#define __REPEAT_4(sort) sort, __REPEAT_3(sort)
#define __REPEAT_5(sort) sort, __REPEAT_4(sort)
#define __REPEAT_6(sort) sort, __REPEAT_5(sort)
#define __REPEAT_7(sort) sort, __REPEAT_6(sort)
#define __REPEAT_8(sort) sort, __REPEAT_7(sort)
#define __REPEAT_9(sort) sort, __REPEAT_8(sort)
#define __REPEAT_10(sort) sort, __REPEAT_9(sort)
#define __REPEAT(arity, sort) __REPEAT_ ## arity(sort)

#define DECL_CONST(f, sort) auto f = ConstSugar(#f, sort);
#define DECL_SKOLEM_CONST(f, sort) auto f = ConstSugar(#f, sort, true);
#define DECL_FUNC(f, ...)   auto f = FuncSugar(#f, __VA_ARGS__);
#define DECL_POLY_FUNC(f, i, ...)   auto f = FuncSugar(#f, __VA_ARGS__, i); 
#define DECL_POLY_CONST(f, i, sort)   auto f = FuncSugar(#f, {}, sort, i);    
#define DECL_PRED(f, ...)   auto f = PredSugar(#f, __VA_ARGS__);
#define DECL_TYPE_CON(f, arity) auto f = TypeConSugar(#f, arity);    
#define DECL_SORT(s)        auto s = TypeConstSugar(#s);
#define DECL_SORT_BOOL      auto Bool = SortSugar(AtomicSort::boolSort());
#define DECL_VAR(x, i) auto x = TermSugar(TermList::var(i));
#define DECL_SORT_VAR(x, i) auto x = SortSugar(TermList::var(i));    
#define DECL_VAR_SORTED(x, i, s) auto x = TermSugar(TermList::var(i), s);
#define DECL_I_COMB(i) auto i = FuncSugar(env.signature->getCombinator(Signature::I_COMB));
#define DECL_K_COMB(k) auto k = FuncSugar(env.signature->getCombinator(Signature::K_COMB));
#define DECL_B_COMB(b) auto b = FuncSugar(env.signature->getCombinator(Signature::B_COMB));
#define DECL_C_COMB(c) auto c = FuncSugar(env.signature->getCombinator(Signature::C_COMB));
#define DECL_S_COMB(s) auto s = FuncSugar(env.signature->getCombinator(Signature::S_COMB));
#define DECL_FUN_DEF(d, t)  auto d = PredSugar(env.signature->getFnDef(t.sugaredExpr().term()->functor()));
#define DECL_PRED_DEF(d, t) auto d = PredSugar(env.signature->getBoolDef(((Literal*)t)->functor()));

#define DECL_DEFAULT_VARS                                                                 \
  __ALLOW_UNUSED(                                                                         \
    DECL_VAR(x, 0)                                                                        \
    DECL_VAR(y, 1)                                                                        \
    DECL_VAR(z, 2)                                                                        \
  )                                                                                       \

// Sort variables are just variables 
// numbers chosen to make a clash with term variables unlikely
// Only these vars should be used when creating polymorphic sorts
// if more are required this macro should be modified. 
#define DECL_DEFAULT_SORT_VARS                                                            \
  __ALLOW_UNUSED(                                                                         \
    DECL_SORT_VAR(alpha, 101)                                                             \
    DECL_SORT_VAR(beta, 102)                                                              \
    DECL_SORT_VAR(gamma, 103)                                                             \
  )

#define DECL_COMBINATORS                                                                  \
  __ALLOW_UNUSED(                                                                         \
    DECL_I_COMB(I)                                                                        \
    DECL_K_COMB(K)                                                                        \
    DECL_B_COMB(B)                                                                        \
    DECL_C_COMB(C)                                                                        \
    DECL_S_COMB(S)                                                                        \
  )


/** tldr: For examples on usage see UnitTesting/tSyntaxSugar.cpp
 *
 * This macro provides syntax sugar for building terms and clauses in a one-sorted theory of numbers. 
 *
 * i.e. the theories supported are:
 * NUMBER_SUGAR(Int) ==> integer arithmetic 
 * NUMBER_SUGAR(Rat) ==> rational arithmetic 
 * NUMBER_SUGAR(Real) ==> rational arithmetic 
 *
 * The macro is meant to be called within a TEST_FUN(...){...} block to import the syntax sugar for 
 * the corresponding test case. It provides a class TermSugar that implicitly converts number literals 
 * to terms in the corresponding sort, and operators to build complex terms and literals. Further it sets 
 * some (C++) variables to variable terms and some to constant terms.
 *
 * Additionally the macros `DECL_FUNC`, `DECL_CONST`, and `DECL_PRED` can be 
 * used to declare uninterpreted functions/predicates/constants over the sort.
 *
 * These are the automatically defined operators and variables
 *
 * local variables:
 * x ... variable X0
 * y ... variable X1
 * z ... variable X2
 *
 * Operators for creating complex terms:
 * operator* ... interpreted multiplication
 * operator+ ... interpreted addition
 * operator- ... interpreted unary minus
 *
 * Operators for creating literals:
 * operator>   ... interpreted greater relation
 * operator>=  ... interpreted greater or equal relation
 * operator<   ... interpreted less relation
 * operator<=  ... interpreted less or equal relation
 * operator==  ... interpreted equality relation
 * operator!=  ... interpreted equality relation
 * operator~   ... flipping a literal's polarity
 *
 * Other closures:
 * frac(int,int) ... creates a fractional interpreted constant (only for REAL, and RAT)
 * num(int)      ... explicity converts a number to an interpreted constant 
 *                   this can be needed in order to prevent the compiler from pre-evaluating integer expressions.
 *                   e.g. {
 *                      Literal* l1 = (a == (3 * 2));
 *                                          ^^^^^^^ this will be evaluated to 6 before transforming it into a term
 *                                                  in order to prevent this we can write:
 *                      Literal* l2 = (a == (num(3) * 2));
 *                   }
 *
 * For examples see UnitTesting/tSyntaxSugar.cpp.
 */
#define NUMBER_SUGAR(Sort)                                                                \
  __ALLOW_UNUSED(                                                                         \
    using NumTraits = Sort##Traits;                                                       \
    syntaxSugarGlobals().setNumTraits(NumTraits{});                                       \
    auto add = FuncSugar(NumTraits::addF());                                              \
    auto mul = FuncSugar(NumTraits::mulF());                                              \
    auto minus = FuncSugar(NumTraits::minusF());                                          \
    auto floor = FuncSugar(NumTraits::floorF());                                          \
    auto ceil = [&](auto t) { return minus(floor(minus(t))); };                           \
    auto Sort = SortSugar(NumTraits::sort());                                             \
  )                                                                                       \

#define DECL_TERM_ALGEBRA(...) createTermAlgebra(__VA_ARGS__);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SyntaxSugar {
public:
  static void reset() {
    env.signature = new Signature();
  }
};

class TermSugar;

class SyntaxSugarGlobals 
{
  static SyntaxSugarGlobals _instance;

  template<class NumTraits>
  void setAllNumTraits() 
  {
    createNumeral = [](int i) {return NumTraits::constantTl(i);};

    add = NumTraits::add;
    mul = NumTraits::mul;

    minus = NumTraits::minus;

    less = NumTraits::less;
    leq  = NumTraits::leq;
    greater = NumTraits::greater;
    geq  = NumTraits::geq;

    isInt = NumTraits::isInt;
    isRat = NumTraits::isRat;
    isReal = NumTraits::isReal;

  }


  template<class NumTraits>
  void setFracTraits() 
  { 
    setAllNumTraits<NumTraits>();
    div = NumTraits::div; 
    createFraction = [](int a, int b) {return NumTraits::constantTl(a,b);};
  }

public:
  static SyntaxSugarGlobals& instance() 
  { return _instance; }


  void setNumTraits(IntTraits)
  {
    setAllNumTraits<IntTraits>();

    quotientT = IntTraits::quotientT;
    remainderT = IntTraits::remainderT;
    quotientF = IntTraits::quotientF;
    remainderF = IntTraits::remainderF;
    quotientE = IntTraits::quotientE;
    remainderE = IntTraits::remainderE;
  }

  void setNumTraits(RatTraits)
  { setFracTraits<RatTraits>(); }


  void setNumTraits(RealTraits)
  { setFracTraits<RealTraits>(); }

  std::function<TermList(int, int)> createFraction;
  std::function<TermList(int)> createNumeral;

  std::function<TermList(TermList, TermList)> add;
  std::function<TermList(TermList, TermList)> mul;
  std::function<TermList(TermList, TermList)> div;

  std::function<TermList(TermList, TermList)> quotientT;
  std::function<TermList(TermList, TermList)> remainderT;
  std::function<TermList(TermList, TermList)> quotientF;
  std::function<TermList(TermList, TermList)> remainderF;
  std::function<TermList(TermList, TermList)> quotientE;
  std::function<TermList(TermList, TermList)> remainderE;

  std::function<TermList(TermList)> minus;


  std::function<Literal*(bool, TermList, TermList)> less;
  std::function<Literal*(bool, TermList, TermList)> leq ;
  std::function<Literal*(bool, TermList, TermList)> greater;
  std::function<Literal*(bool, TermList, TermList)> geq ;

  std::function<Literal*(bool, TermList)> isInt;
  std::function<Literal*(bool, TermList)> isRat;
  std::function<Literal*(bool, TermList)> isReal;

};

inline SyntaxSugarGlobals& syntaxSugarGlobals() 
{ return SyntaxSugarGlobals::instance(); }

class ExpressionSugar
{
public: 
  // TODO get rid of this default constructor. we never want to create uninitialized instances
  ExpressionSugar() 
   : _sugaredExpr(TermList::empty())
  { }

  ExpressionSugar(TermList sugaredExpr) : 
    _sugaredExpr(sugaredExpr){}

  /** explicit conversion */ 
  // TODO rename to what it actually is (an explicit concersion) toTermList
  TermList sugaredExpr() const { return _sugaredExpr;} 

  /** implicit conversion */ 
  operator TermList() const {return _sugaredExpr;}

protected:
  TermList _sugaredExpr;
};


struct SortSugar : public ExpressionSugar
{
public:

  SortSugar(SortId srt) : ExpressionSugar(srt) {
    ASS(srt.isVar() || srt.term()->isSort());    
  }

  SortSugar(const char* name, Stack<SortSugar> as_) 
  {
    if(as_.isEmpty()){
      _sugaredExpr = TermList(AtomicSort::createConstant(name));
    } else {
      Stack<SortId> as;     
      for (auto a : as_){ as.push(a.sugaredExpr()); }
      _sugaredExpr = AtomicSort::arrowSort(as, as.pop());      
    }
  }

  SortSugar(const char* name) 
    : SortSugar(TermList(AtomicSort::createConstant(name))) 
  {  }
};

class TermSugar : public ExpressionSugar
{
  SortId _srt;

public:
  TermSugar(bool foolConst) 
    : TermSugar(TermList(foolConst ? Term::foolTrue() : Term::foolFalse()))
  {}

  TermSugar(int trm) 
    : TermSugar(TermList(syntaxSugarGlobals().createNumeral(trm)))
  {}

  TermSugar(TermList trm) 
    : ExpressionSugar(trm)
  { 
    ASS_REP(!_sugaredExpr.isEmpty(), _sugaredExpr);
    if (_sugaredExpr.isVar()) {
      _srt = TermList::empty();
    } else {
      if (_sugaredExpr.term()->isLiteral()) {
        _srt = AtomicSort::boolSort();
      } else {
        _srt = SortHelper::getResultSort(_sugaredExpr.term());
      }
    }
  }

  TermSugar(TermList trm, SortSugar sort)
    : TermSugar(trm)
  {
    ASS(_sugaredExpr.isVar());
    _srt = sort.sugaredExpr();
  }

  SortId sort() const { return _srt; }

  TermSugar sort(SortId s) { _srt = s; return *this; }

  static TermSugar createConstant(const char* name, SortSugar s, bool skolem) {
    unsigned f = env.signature->addFunction(name,0);                                                                

    env.signature->getFunction(f)->setType(OperatorType::getFunctionType({}, s.sugaredExpr()));
    if (skolem) {
      env.signature->getFunction(f)->markSkolem();
    }
    return TermSugar(TermList(Term::createConstant(f)));                                                          
  }                                                                                                                 

  operator TypedTermList() const { return TypedTermList(TermList(*this), sort()); }
};

class SortedTermSugar : public TermSugar
{
  SortSugar _sort;
public:
  SortedTermSugar(TermList term, SortSugar sort) : TermSugar(term), _sort(sort) 
  {}

  SortSugar sort() const { return _sort; }
};


class Lit
{
  Literal* _lit;
  bool _selected;
public:
  Lit(Literal* lit) : _lit(lit), _selected(false) {}

  operator Literal*() const 
  { return _lit; }

  bool selected() const 
  { return _selected; }

  friend Lit selected(Lit l) 
  {
    l._selected = true;
    return l;
  }

  TermSugar wrapInTerm()
  {
    return TermSugar(TermList(_lit));
  }
};

inline SortedTermSugar sorted(TermList var, SortSugar sort) 
{ return SortedTermSugar(var, sort); }

inline TermSugar frac(int a, int b) 
{ return syntaxSugarGlobals().createFraction(a,b); }

inline TermSugar num(int a)
{ return syntaxSugarGlobals().createNumeral(a); }

inline TermSugar fool(bool b)
{ return TermSugar(b); }

////////////////////////// operators to create terms ////////////////////////// 

inline TermSugar operator-(TermSugar x) { return syntaxSugarGlobals().minus(x); }

inline TermSugar ap(SortSugar sort, TermSugar lhs, TermSugar rhs) 
{ return ApplicativeHelper::createAppTerm(sort, lhs, rhs); }  

inline TermSugar ap(TermSugar lhs, TermSugar rhs) 
{ return ap(lhs.sort(), lhs, rhs); }  

inline TermSugar operator+(TermSugar lhs, TermSugar rhs)  { return syntaxSugarGlobals().add(lhs, rhs); }  
inline TermSugar operator-(TermSugar lhs, TermSugar rhs)  { return lhs + -rhs; }  
inline TermSugar operator*(TermSugar lhs, TermSugar rhs)  { return syntaxSugarGlobals().mul(lhs, rhs); }  
inline TermSugar operator/(TermSugar lhs, TermSugar rhs)  { return syntaxSugarGlobals().div(lhs, rhs); }  

#define __IMPL_NUMBER_BIN_FUN(op, result_t)                                               \
  inline result_t op(int lhs, TermSugar rhs) { return op(TermSugar(lhs), rhs); }          \
  inline result_t op(TermSugar lhs, int rhs) { return op(lhs, TermSugar(rhs)); }          \

__IMPL_NUMBER_BIN_FUN(operator+, TermSugar)
__IMPL_NUMBER_BIN_FUN(operator*, TermSugar)
__IMPL_NUMBER_BIN_FUN(operator/, TermSugar)

#define __BIN_FUNC_QUOTIENT_REMAINDER(X)                                                  \
  inline TermSugar  quotient##X(TermSugar lhs, TermSugar rhs){ return syntaxSugarGlobals(). quotient##X(lhs, rhs); }   \
  inline TermSugar remainder##X(TermSugar lhs, TermSugar rhs){ return syntaxSugarGlobals().remainder##X(lhs, rhs); }   \
                                                                                          \
  __IMPL_NUMBER_BIN_FUN( quotient##X, TermSugar)                                          \
  __IMPL_NUMBER_BIN_FUN(remainder##X, TermSugar)                                          \

__BIN_FUNC_QUOTIENT_REMAINDER(E)
__BIN_FUNC_QUOTIENT_REMAINDER(T)
__BIN_FUNC_QUOTIENT_REMAINDER(F)

#undef __BIN_FUNC_QUOTIENT_REMAINDER



////////////////////////// operators to create literals ////////////////////////// 


inline Lit operator==(SortedTermSugar lhs, TermSugar rhs) 
{ return Literal::createEquality(true, lhs, rhs, lhs.sort().sugaredExpr()); }

inline Lit operator==(TermSugar lhs, SortedTermSugar rhs) 
{ return Literal::createEquality(true, lhs, rhs, rhs.sort().sugaredExpr()); }

inline Lit operator==(TermSugar lhs, TermSugar rhs) 
{
  SortId sort = lhs.sort().isNonEmpty() ? lhs.sort() : rhs.sort();
  return Literal::createEquality(true, lhs, rhs, sort);
}


inline Lit operator< (TermSugar lhs, TermSugar rhs) { return syntaxSugarGlobals().less   (true, lhs, rhs); }
inline Lit operator<=(TermSugar lhs, TermSugar rhs) { return syntaxSugarGlobals().leq    (true, lhs, rhs); }
inline Lit operator> (TermSugar lhs, TermSugar rhs) { return syntaxSugarGlobals().greater(true, lhs, rhs); }
inline Lit operator>=(TermSugar lhs, TermSugar rhs) { return syntaxSugarGlobals().geq    (true, lhs, rhs); }

inline Lit isInt (TermSugar trm) { return syntaxSugarGlobals().isInt (true, trm); }
inline Lit isRat (TermSugar trm) { return syntaxSugarGlobals().isRat (true, trm); }
inline Lit isReal(TermSugar trm) { return syntaxSugarGlobals().isReal(true, trm); }

inline Lit operator~(Lit lit) 
{
  Literal* l = lit;
  return Literal::create(l, !l->polarity());
}

inline Lit operator!=(SortedTermSugar lhs, TermSugar rhs)  { return ~(lhs == rhs); }
inline Lit operator!=(TermSugar lhs, SortedTermSugar rhs)  { return ~(lhs == rhs); }
inline Lit operator!=(TermSugar lhs, TermSugar rhs) { return ~(lhs == rhs); }

__IMPL_NUMBER_BIN_FUN(operator==, Lit)
__IMPL_NUMBER_BIN_FUN(operator!=, Lit)
__IMPL_NUMBER_BIN_FUN(operator< , Lit)
__IMPL_NUMBER_BIN_FUN(operator<=, Lit)
__IMPL_NUMBER_BIN_FUN(operator> , Lit)
__IMPL_NUMBER_BIN_FUN(operator>=, Lit)


inline SortSugar arrow(TermList args, TermList res) 
{ return AtomicSort::arrowSort({ args }, res);      }

inline SortSugar arrow(Stack<TermList> args, TermList res) 
{ return AtomicSort::arrowSort(args, res);      }
  
class FuncSugar {
  unsigned _functor;
  unsigned _arity;

public:
  explicit FuncSugar(unsigned functor) 
    : _functor(functor)
    , _arity(env.signature->getFunction(functor)->arity()) {}

  FuncSugar(std::string const& name, std::initializer_list<SortSugar> as_, 
    ExpressionSugar result, unsigned taArity = 0) 
  {
    Stack<SortId> as;
    for (auto a : as_) 
      as.push(a.sugaredExpr());

    bool added = false;
    _functor = env.signature->addFunction(name, as.size() + taArity, added);
    _arity = as.size();
    if (added){
      TermList res = result.sugaredExpr();

      if(taArity){
        TermStack vars = {TermList(101, false), TermList(102, false), TermList(103, false)};
        SortHelper::normaliseArgSorts(vars, as);
        SortHelper::normaliseSort(vars, res);
      }

      env.signature
        ->getFunction(_functor)
        ->setType(OperatorType::getFunctionType(as.size(), as.begin(), res, taArity));    
    }
  }

  FuncSugar dtor(unsigned i) const {
    ASS_L(i, arity())
    ASS (symbol()->termAlgebraCons()) 
    return FuncSugar(
        env.signature->getTermAlgebraConstructor(functor())
          ->destructorFunctor(i));
  }

  auto result()        const { return symbol()->fnType()->result(); }
  auto arg(unsigned i) const { return symbol()->fnType()->arg(i); }

  template<class... As>
  TermSugar operator()(As... args) const {
    Stack<TermList> as { TermSugar(args).sugaredExpr()... };
    return TermList(Term::create(_functor, 
        as.size(), 
        as.begin()));
  }
  unsigned functor() const { return _functor; }
  unsigned arity() const { return _arity; }
  Signature::Symbol* symbol() const { return env.signature->getFunction(functor()); }

  friend std::ostream& operator<<(std::ostream& out, FuncSugar const& self) 
  { return out << self.symbol()->name(); }
};

class ConstSugar : public TermSugar, public FuncSugar
{
public:

  ConstSugar(const char* name, SortSugar s, bool skolem = false)
    : TermSugar(TermSugar::createConstant(name, s, skolem).sugaredExpr())
    , FuncSugar(functor())
  { }
  unsigned functor() const { return this->sugaredExpr().term()->functor(); }
};

class TypeConSugar {
  unsigned _functor;

public:
  TypeConSugar(const char* name, unsigned arity)
  {
    bool added = false;
    _functor = env.signature->addTypeCon(name, arity, added);
    if (added)
      env.signature
        ->getTypeCon(_functor)
        ->setType(OperatorType::getTypeConType(arity));   
  }

  template<class... As>
  SortSugar operator()(As... args) const {
    Stack<TermList> as { SortSugar(args).sugaredExpr()... };
    return TermList(AtomicSort::create(_functor, 
        as.size(), 
        as.begin() ));
  }
  unsigned functor() const { return _functor; }
};

class TypeConstSugar : public SortSugar, public TypeConSugar
{
public:

  TypeConstSugar(const char* name) 
    : SortSugar(name) 
    , TypeConSugar(name, 0)
  { }
  unsigned functor() const { return this->sugaredExpr().term()->functor(); }
};

class PredSugar {
  unsigned _functor;

public:
  PredSugar(unsigned functor)  : _functor(functor) {}

  PredSugar(const char* name, std::initializer_list<SortSugar> args, unsigned taArity = 0)
  {
    Stack<SortId> as;
    for (auto a : args) {
      as.push(a.sugaredExpr());
    }

    if(taArity){
      // TODO don't haredcode these varible numbers?!
      TermStack vars = {TermList(101, false), TermList(102, false), TermList(103, false)};      
      SortHelper::normaliseArgSorts(vars, as);
    }

    _functor = env.signature->addPredicate(name, as.size() + taArity);
    env.signature
      ->getPredicate(_functor)
      ->setType(OperatorType::getPredicateType(as.size(), as.begin(), taArity));    
  }

  template<class... As>
  Lit operator()(As... args) const {
    Stack<TermList> as { TermSugar(args).sugaredExpr()... };
    return Literal::create(_functor,
        as.size(),
        /* polarity */ true,
        as.begin() );
  }
  unsigned functor() const { return _functor; }
};

inline Clause* clause(Stack<Lit> ls) { 
  static Inference testInf = Kernel::NonspecificInference0(UnitInputType::ASSUMPTION, InferenceRule::INPUT); 

  std::stable_sort(ls.begin(), ls.end(), [](Lit const& l1, Lit const& l2){ return l1.selected() > l2.selected(); });
  auto nSelected = iterTraits(ls.iterFifo())
    .findPosition([](Lit const& l) 
        { return !l.selected(); })
    .unwrapOrElse( [&]() {return ls.size(); });

  Clause& out = *Clause::fromIterator(arrayIter(ls)
      .map([](Lit l) -> Literal* { return l; }), testInf);

  out.setSelected(nSelected);
  return &out; 
}

inline Clause* clause(std::initializer_list<Lit> ls) 
{ return clause(Stack<Lit>(ls)); }

inline Stack<Clause*> clauses(std::initializer_list<std::initializer_list<Lit>> cls) { 
  auto out = Stack<Clause*>();
  for (auto cl : cls) {
    out.push(clause(cl));
  }
  return out;
}

inline void createTermAlgebra(SortSugar sort, std::initializer_list<FuncSugar> fs) {
  // avoid redeclaration
  if (env.signature->isTermAlgebraSort(sort.sugaredExpr())) {
    return;
  }

  using namespace Shell;

  Stack<FuncSugar> funcs = fs;
  Stack<TermAlgebraConstructor*> cons;

  for (auto f : funcs) {
    env.signature->getFunction(f.functor())
      ->markTermAlgebraCons();

    auto dtor = [&](unsigned i) {
      std::stringstream name;
      name << f << "@" << i;
      auto d = FuncSugar(name.str(), { f.result() }, f.arg(i));
      env.signature->getFunction(d.functor())
        ->markTermAlgebraDest();
      return d;
    };

    Array<unsigned> dtors(f.arity()); 
    for (unsigned i = 0; i < f.arity(); i++) {
      dtors[i] = dtor(i).functor();
    }

    cons.push(new TermAlgebraConstructor(f.functor(), dtors));
  }
  env.signature->addTermAlgebra(new TermAlgebra(sort.sugaredExpr(), cons.size(), cons.begin()));
}

#endif // __TEST__SYNTAX_SUGAR__H__
