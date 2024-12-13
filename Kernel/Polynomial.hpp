/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */

#ifndef __POLYNOMIAL__H__
#define __POLYNOMIAL__H__

/**
 * @file Kernel/Polynomial.hpp
 * This file mainly provides POLYnomial Normal Forms for terms PolyNf. In this normal form a term is defined as follows:
 *
 * PolyNf ::= FuncId [PolyNf]  // a normal function term
 *          | Variable         // a normal variable
 *          | AnyPoly          // a dynamically typed polynom
 *
 * AnyPoly ::= Polynom<IntTraits > 
 *           | Polynom<RatTraits >
 *           | Polynom<RealTraits>
 *
 * Polynom<Number> ::= [Monom<Number>]              // a sorted list of monoms
 * Monom<Number>   ::= Number [MonomFactor<Number>] // a numeral as factor, and a sorted list of other factors
 * MonomFactor     ::= PolyNf int                   // the term of the factor, and its power
 */

#include "Kernel/NumTraits.hpp"
#include <cassert>
#include "Lib/Coproduct.hpp"
#include "Lib/Option.hpp"
#include "Lib/Map.hpp"
#include "Kernel/Theory.hpp"
#include "Lib/Perfect.hpp"
#include "Kernel/NumTraits.hpp"
#include "Kernel/Ordering.hpp"
#include "Kernel/TypedTermList.hpp"
#include "Lib/Reflection.hpp"
#include <type_traits>

#define DEBUG(...) // DBG(__VA_ARGS__)

namespace Kernel {

// TODO use this newtype in Term.hpp
/** newtype for wrapping varible ids */
class Variable 
{
  unsigned _num;
public: 
  Variable();
  explicit Variable(unsigned num);
  unsigned id() const;

  friend struct std::hash<Variable>;
  friend bool operator==(Variable lhs, Variable rhs);
  friend bool operator!=(Variable lhs, Variable rhs);
  friend bool operator<(Variable const& lhs, Variable const& rhs);
  friend std::ostream& operator<<(std::ostream& out, const Variable& self);
};

} // namespace Kernel

template<> struct std::hash<Kernel::Variable>
{
  size_t operator()(Kernel::Variable const& self)
  { return std::hash<unsigned>{}(self._num); }
};

namespace Kernel {

// TODO use this newtype in Term.hpp
/** newtype for wrapping function ids (also often called functors in vampire) */
class FuncId 
{
  unsigned _num;
  const TermList* _typeArgs;
  
public: 
  explicit FuncId(unsigned num, const TermList* typeArgs);
  static FuncId symbolOf(Term* term);
  unsigned numTermArguments();
  TermList typeArg(unsigned i) const { return *(_typeArgs - i); }
  unsigned numTypeArgs() const { return env.signature->getFunction(_num)->numTypeArguments(); }

  friend std::ostream& operator<<(std::ostream& out, const FuncId& self);
  auto iterTypeArgs() const 
  { return range(0, numTypeArgs()).map([&](auto i) { return typeArg(i); }); }

  Signature::Symbol* symbol() const;

  unsigned id() const;
  Theory::Interpretation interpretation() const;
  bool isInterpreted() const;
  Option<Theory::Interpretation> tryInterpret() const;

  template<class Number>
  Option<typename Number::ConstantType> tryNumeral() const;
  
  auto asTuple() const { return std::tuple(_num, iterContOps(iterTypeArgs())); }
  IMPL_COMPARISONS_FROM_TUPLE(FuncId)
  IMPL_HASH_FROM_TUPLE(FuncId)
};

} // namespace Kernel


/////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations, needed to define PolyNf structure
/////////////////////////////////////////////////////////

namespace Kernel {

class PolyNf;
template<class Number> class Polynom;
template<class Number> class MonomFactors;
class AnyPoly;

////////////////////////////////////////////////////////////////////////////////////////////////////
// class declartions for PolyNf
/////////////////////////////////////////////////////////

/** 
 * Represents a summand in a polynom of the number type Number 
 * e.g. a term like 3 * (a*x) would be represented as 
 * Monom { 3, MonomFactors { a, x }}
 */
template<class Number> 
struct Monom 
{
  USE_ALLOCATOR(Monom)

  using Numeral = typename Number::ConstantType;

  Numeral numeral;
  Perfect<MonomFactors<Number>> factors;

  Monom(Numeral numeral, Perfect<MonomFactors<Number>> factors);

  static Monom zero();

  Option<Variable> tryVar() const;

  /** performs an integrity check on the datastructure, only has an effect in debug mode */
  void integrity() const;
};


/**
 * Represents an ordenary complex term, in the PolyNF term tree.
 */
class FuncTerm 
{
  FuncId _fun;
  Stack<PolyNf> _args;
public:
  USE_ALLOCATOR(FuncTerm)

  FuncTerm(FuncId f, Stack<PolyNf>&& args);
  FuncTerm(FuncId f, PolyNf* args);

  unsigned numTermArguments() const;
  FuncId function() const;
  PolyNf const& arg(unsigned i) const;

  template<class Number> 
  Option<typename Number::ConstantType> tryNumeral() const;

  friend std::ostream& operator<<(std::ostream& out, const FuncTerm& self);
  friend bool operator==(FuncTerm const& lhs, FuncTerm const& rhs);
  friend bool operator!=(FuncTerm const& lhs, FuncTerm const& rhs);
  friend struct std::hash<FuncTerm>;
};

using AnyPolySuper = Coproduct< 
    Perfect<Polynom< IntTraits>>
  , Perfect<Polynom< RatTraits>>
  , Perfect<Polynom<RealTraits>>
  >;

class AnyPoly : public AnyPolySuper
{
public:
  /** creates a new dynamically typed polynom from a statically typed one */
  template<class NumTraits> AnyPoly(Perfect<Polynom<NumTraits>> x);

  /** tries to turn this polynom into a polynom of the given NumTraits. */
  template<class NumTraits> Option<Perfect<Polynom<NumTraits>> const&> downcast() const&;

  /** returns wether this is a Polynom of the given NumTraits. */
  template<class NumTraits> bool isType() const;

  /** if this polynom has the right sort, and consist of a single summand that is a numeral, then this numeral
   * is returned. otherwise an empty Option is returned. */
  template<class NumTraits> Option<typename NumTraits::ConstantType> tryNumeral() const&;

  /** \see template<class N> Polynom<N>::replaceTerms */
  AnyPoly replaceTerms(PolyNf* newTs) const;

  /** \see template<class N> Polynom<N>::nSummands */
  unsigned nSummands() const;

  /** \see template<class N> Polynom<N>::nFactors */
  unsigned nFactors(unsigned i) const;

  /** \see template<class N> Polynom<N>::termAt */
  PolyNf const& termAt(unsigned summand, unsigned factor) const;

  /** \see template<class N> Polynom<N>::denormalize */
  TermList denormalize(TermList* results) const;

  friend std::ostream& operator<<(std::ostream& out, const AnyPoly& self);
  friend struct std::hash<AnyPoly>;
};

using PolyNfSuper = Lib::Coproduct<Perfect<FuncTerm>, Variable, AnyPoly>;

/**
 * Represents the polynomial normal form of a term, that is used for performing several simplifications and evaluations.
 *
 * See the file-level documentation for how this datatype is composed.
 */
class PolyNf : public PolyNfSuper
{
public:

  PolyNf(Perfect<FuncTerm> t);
  PolyNf(Variable               t);
  PolyNf(AnyPoly                t);

  static PolyNf normalize(TypedTermList t);

  /** 
   * If this term is a polynomial of sort NumTraits, it is downcasted to that sort,
   * otherwise an empty Option is returned
   */
  template<class NumTraits>
  Option<Perfect<Polynom<NumTraits>>> downcast() const;

  /** turns the normal form term PolyNf into an ordenary vampire term TermList. */
  TermList denormalize() const;

  /** 
   * Turns this PolyNf term into a typed polynom of sort Number.
   * this must have the same sort as Number. 
   * If this is already a polynom it will just be downcasted, 
   * otherwise (when it is a Variable, or a FuncTerm) it will be 
   * wrapped in a polynom.
   */
  template<class Number> 
  Perfect<Polynom<Number>> wrapPoly() const;
  

  /** if this PolyNf is a numeral, the numeral is returned */
  template<class Number>
  Option<typename Number::ConstantType> tryNumeral() const;

  /** if this PolyNf is a Variable, the variable is returned */
  Option<Variable> tryVar() const;

  /** an iterator over all PolyNf s that are subterms of this one */
  class SubtermIter;

  /** returns an iterator over all PolyNf s that are subterms of this one */
  IterTraits<SubtermIter> iterSubterms() const;

  friend struct std::hash<PolyNf>;
  friend bool operator==(PolyNf const& lhs, PolyNf const& rhs);
  friend bool operator!=(PolyNf const& lhs, PolyNf const& rhs);
  friend bool operator<(const PolyNf& lhs, const PolyNf& rhs);
  friend bool operator<=(const PolyNf& lhs, const PolyNf& rhs);
  friend std::ostream& operator<<(std::ostream& out, const PolyNf& self);
};


/** 
 * Represents a factor in a monom. Each unique term contained in the monom is stored 
 * together with the number of occurences of the term within that monom.
 *
 * e.g. a term like (x * x * a * x) is represented as 
 * MonomFactors { MonomFactor(x, 3), MonomFactor(a, 1) }
 */
template<class Number> 
struct MonomFactor 
{
  PolyNf term;
  int power;

  MonomFactor(PolyNf term, int power);
   
  /** if this monomfactor is a Variable and has power one it is turned into a variable */
  Option<Variable> tryVar() const;
};



/** 
 * Represents the non-numeral part of a monom. this means it is a sorted list of MonomFactor s.
 */
template<class Number>
class MonomFactors 
{
  using MonomFactor = Kernel::MonomFactor<Number>;
  using Monom       = Kernel::Monom<Number>;
  using Polynom     = Kernel::Polynom<Number>;
  using Numeral = typename Number::ConstantType;
  Stack<MonomFactor> _factors;
  friend struct std::hash<MonomFactors>;

public:
  USE_ALLOCATOR(MonomFactors)

  /** 
   * constructs a new MonomFactors. 
   * \pre factors must be sorted
   */
  MonomFactors(Stack<MonomFactor>&& factors);

  /** constructs an empty product, which corresponds to the numeral one.  */
  MonomFactors();

  /** constructs a singleton product */
  MonomFactors(PolyNf t);

  /** constructs the product of t1 and t2. There is no precondigion on the ordering of t1 and t2. */
  MonomFactors(PolyNf t1, PolyNf t2);

  /** returns the number of factors */
  unsigned nFactors() const;

  /** returns the nth factor */
  MonomFactor      & factorAt(unsigned i);

  /** returns the nth factor */
  MonomFactor const& factorAt(unsigned i) const;

  /** returns the number of factors */
  PolyNf const& termAt(unsigned i) const;

  /** returns whether this monom is a polynom, i.e. if its only factor is a polynom */
  bool isPolynom() const;

  /** turns this monom into a polynom. 
   * \pre isPolynom  must be true
   */
  Perfect<Polynom> asPolynom() const;

  /** returns the (empty) product one */
  static MonomFactors one();

  /** returns whether this is an empty product */
  bool isOne() const;


  /** if this MonomFactors consist of a single variable if will be returnd  */
  Option<Variable> tryVar() const;

  /** performs an integrity check on the datastructure, only has an effect in debug mode */
  void integrity() const;

  /** replaces all the factors, by new ones, keeping the power of each term the same  */
  MonomFactors replaceTerms(PolyNf* simplifiedTerms) const;


  /** returns an iterator over all factors */
  auto iter() const&
  { return arrayIter(_factors); }

  explicit MonomFactors(const MonomFactors&) = default;
  explicit MonomFactors(MonomFactors&) = default;

  MonomFactors& operator=(MonomFactors&&) = default;
  MonomFactors(MonomFactors&&) = default;

  template<class N> friend std::ostream& operator<<(std::ostream& out, const MonomFactors<N>& self);
  template<class N> friend bool operator==(const MonomFactors<N>& l, const MonomFactors<N>& r);

  /** helper function for PolyNf::denormalize() */
  TermList denormalize(TermList* results)  const;

  Stack<MonomFactor>& raw();
};

template<class Number>
class Polynom 
{
  friend struct std::hash<Polynom>;

  using Numeral      = typename Number::ConstantType;
  using MonomFactors = Kernel::MonomFactors<Number>;
  using Monom        = Kernel::Monom<Number>;

  Stack<Monom> _summands;

public:
  USE_ALLOCATOR(Polynom)

  /** 
   * constructs a new Polynom with a list of summands 
   * \pre summands must be sorted
   */
  explicit Polynom(Stack<Monom>&& summands);

  /** creates a Polynom that consists of only one summand */
  explicit Polynom(Monom m);

  /** creates a Polynom that consists of only one summand, that is a product of numeral, and term */
  explicit Polynom(Numeral numeral, PolyNf term);

  /** creates a Polynom that consists of only one summand */
  explicit Polynom(PolyNf t);

  /** creates a Polynom that consists of only one summand */
  explicit Polynom(Numeral constant);

  /** creates the empty Polynom that is interpreted as zero */
  static Polynom zero();

  /** if this Polynom consists only of one summand that is a numeral the numeral is returned */
  Option<Numeral> toNumber() const&;

  /** turns this Polynom into a numeral if it consists only of one summand that is a numeral, or throws an assertion violation 
   * \pre isNumeral is true*
   */
  Numeral unwrapNumber() const&;

  /** returns whether this Monom consists of only one summand that is a numeral */
  bool isNumber() const&;

  /** turns this polynom into a term */
  TermList denormalize() const;

  /** helper function for denormalize() */
  TermList denormalize(TermList* results) const;

  /**
   * replaces all subterms of this poly, by given array of subterms. the resulting polynom will be sorted correctly. 
   * here a monom does not count as a subterm, but all the subterms of the monom themselves do:
   *      Polynom(Monom(3, { f(x), y }), Monom(1, { z }))
   *         .replaceTerms({a,b,c})
   * ===> Polynom(Monom(3, { a   , b }), Monom(1, { c }))
   */
  Polynom replaceTerms(PolyNf* simplifiedTerms) const;

  /** returns the number of summands of this polynom */
  unsigned nSummands() const;

  /** returns the number of factors of the summand with the given index */
  unsigned nFactors(unsigned summand) const;

  /** returns the summand with the given index */
  Monom const& summandAt(unsigned summand) const;

  /** returns the summand with the given index */
  Monom      & summandAt(unsigned summand);

  /** integrity check of the data structure. does noly have an effect in debug mode */
  void integrity() const;

  /** returns iterator over all summands of this Polyom */
  auto iterSummands() const&
  { return arrayIter(_summands); }

  Stack<Monom>& raw();

  template<class N> friend bool operator==(const Polynom<N>& lhs, const Polynom<N>& rhs);
  template<class N> friend std::ostream& operator<<(std::ostream& out, const Polynom<N>& self);
};  


/** an iterator over a literal's arguments. The arguments are mapped to their corresponding PolNf s */
class IterArgsPnf
{
  Literal* _lit;
  unsigned _idx;
public:
  DECL_ELEMENT_TYPE(PolyNf);

  IterArgsPnf(Literal* lit);

  bool hasNext() const;
  PolyNf next();
};

/** convienent constructor for IterArgsPnf */
IterTraits<IterArgsPnf> iterArgsPnf(Literal* lit);

} // namespace Kernel

// include needs to go here, since we need the specialization BottomUpChildIter<PolyNf> to declare Iter
#include "Kernel/BottomUpEvaluation.hpp"

namespace Kernel {

class PolyNf::SubtermIter {
  Stack<BottomUpChildIter<PolyNf>> _stack;
public:
  DECL_ELEMENT_TYPE(PolyNf);

  SubtermIter(SubtermIter&&) = default;
  SubtermIter& operator=(SubtermIter&&) = default;

  SubtermIter(PolyNf p);

  PolyNf next();

  bool hasNext() const;
};

} // namespace Kernel

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// TEMPLATE STUFF IMPLEMENTATIONS
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
// impl Monom template stuff
////////////////////////////
namespace Kernel {

template<class Number>
Monom<Number>::Monom(Monom<Number>::Numeral numeral, Perfect<MonomFactors<Number>> factors) 
  : numeral(numeral), factors(factors)
{}

template<class Number>
Monom<Number> Monom<Number>::zero() 
{ 
  static Monom p = Monom(Numeral(0), perfect(MonomFactors<Number>()));
  return p; 
}

template<class Number>
Option<Variable> Monom<Number>::tryVar() const 
{
  using Opt = Option<Variable>;
  if (numeral == Numeral(1)) {
    return  factors->tryVar();
  } else {
    return  Opt();
  }
}

template<class Number>
void Monom<Number>::integrity() const 
{
#if VDEBUG
  this->factors->integrity();
#endif // VDEBUG
}


template<class Number>
bool operator<(Monom<Number> const& l, Monom<Number> const& r)
{ return std::tie(l.factors, l.numeral) < std::tie(r.factors, r.numeral); }

template<class Number>
bool operator==(Monom<Number> const& l, Monom<Number> const& r)
{ return std::tie(l.factors, l.numeral) == std::tie(r.factors, r.numeral); }

template<class Number>
bool operator!=(Monom<Number> const& l, Monom<Number> const& r)
{ return !(l == r); }

template<class Number>
std::ostream& operator<<(std::ostream& out, const Monom<Number>& self)
{ 
  if (self.numeral != typename Number::ConstantType(1)) {
    out << self.numeral;
  }
  return out << self.factors; 
}




} // namespace Kernel


/////////////////////////////////////////////////////////
// impl FuncId template stuff
////////////////////////////

namespace Kernel {

template<class Number>
Option<typename Number::ConstantType> FuncId::tryNumeral() const
{ 
  using Numeral = typename Number::ConstantType;
  Numeral out;
  if (theory->tryInterpretConstant(_num, out)) {
    return Option<Numeral>(out);
  } else {
    return Option<Numeral>();
  }
}

} // namespace Kernel


/////////////////////////////////////////////////////////
// impl FuncTerm template stuff
////////////////////////////


namespace Kernel {

template<class Number>
Option<typename Number::ConstantType> FuncTerm::tryNumeral() const
{ return _fun.template tryNumeral<Number>(); }

} // namespace Kernel


template<> struct std::hash<Kernel::FuncTerm> 
{
  size_t operator()(Kernel::FuncTerm const& f) const 
  { return Lib::HashUtils::combine(f._fun.defaultHash(), std::hash<Stack<Kernel::PolyNf>>{}(f._args));  }
};

/////////////////////////////////////////////////////////
// impl AnyPoly  template stuff
////////////////////////////

namespace Kernel {

template<class NumTraits>
AnyPoly::AnyPoly(Perfect<Polynom<NumTraits>> x) : Coproduct(std::move(x)) {  }

template<class NumTraits> 
Option<Perfect<Polynom<NumTraits>> const&>  AnyPoly::downcast() const& 
{ return as<Perfect<Polynom<NumTraits>>>(); }

template<class NumTraits> 
bool AnyPoly::isType() const 
{ return is<Perfect<Polynom<NumTraits>>>(); }

/* helper function for AnyPoly::tryNumeral */
template<class NumIn, class NumOut>
struct __ToNumeralHelper 
{
  Option<typename NumOut::ConstantType> operator()(Perfect<Polynom<NumIn>>) const
  { return Option<typename NumOut::ConstantType>(); }
};

/* helper function for AnyPoly::tryNumeral */
template<class Num>
struct __ToNumeralHelper<Num,Num>
{
  Option<typename Num::ConstantType> operator()(Perfect<Polynom<Num>> p) const
  { return p->toNumber(); }
};

template<class NumOut>  
struct PolymorphicToNumeral 
{
  template<class NumIn>
  Option<typename NumOut::ConstantType> operator()(Perfect<Polynom<NumIn>> const& p) const
  { return __ToNumeralHelper<NumIn, NumOut>{}(p); }
};


template<class NumTraits>
Option<typename NumTraits::ConstantType> AnyPoly::tryNumeral() const&
{ return apply(PolymorphicToNumeral<NumTraits>{}); }

} // namespace Kernel

template<> struct std::hash<Kernel::AnyPoly> 
{
  size_t operator()(Kernel::AnyPoly const& self) const 
  { return std::hash<Kernel::AnyPolySuper>{}(self); }
};



/////////////////////////////////////////////////////////
// impl PolyNf  template stuff
////////////////////////////

namespace Kernel {


template<class NumTraits>
Option<Perfect<Polynom<NumTraits>>> PolyNf::downcast()  const
{
  using Result = Perfect<Polynom<NumTraits>>;
  return as<AnyPoly>()
    .andThen([](AnyPoly const& p) { return p.as<Result>(); })
    .map([](Result const& p) -> Result { return p; });
}


template<class Number> 
Perfect<Polynom<Number>> PolyNf::wrapPoly() const
{
  if (this->is<AnyPoly>()) {
    return this->unwrap<AnyPoly>()
            .unwrap<Perfect<Polynom<Number>>>();
  } else {
    return perfect(Polynom<Number>(*this));
  }
}

template<class Number>
Option<typename Number::ConstantType> PolyNf::tryNumeral() const
{ 
  using Numeral = typename Number::ConstantType;
  return match(
      [](Perfect<FuncTerm> t) { return (*t).tryNumeral<Number>(); },
      [](Variable               t) { return Option<Numeral>();              },
      [](AnyPoly                t) { return t.template tryNumeral<Number>(); }
    );
}

} // namespace Kernel

template<> struct std::hash<Kernel::PolyNf> 
{
  size_t operator()(Kernel::PolyNf const& f) const 
  { return std::hash<Kernel::PolyNfSuper>{}(f); }
};

/////////////////////////////////////////////////////////
// impl MonomFactor  tempalte stuff
////////////////////////////

namespace Kernel {

template<class Number>
MonomFactor<Number>::MonomFactor(PolyNf term, int power) 
  : term(term)
  , power(power)
{}

template<class Number>
bool operator<(MonomFactor<Number> const& l, MonomFactor<Number> const& r)
{ return std::tie(l.term, l.power) < std::tie(r.term, r.power); }

template<class Number>
bool operator==(MonomFactor<Number> const& l, MonomFactor<Number> const& r)
{ return std::tie(l.term, l.power) == std::tie(r.term, r.power); }

template<class Number>
bool operator!=(MonomFactor<Number> const& l, MonomFactor<Number> const& r)
{ return !(l == r); }

template<class Number>
std::ostream& operator<<(std::ostream& out, const MonomFactor<Number>& self) {
  out << self.term; 
  if (self.power != 1) 
    out << "^" << self.power;
  return out;
}

template<class Number>
Option<Variable> MonomFactor<Number>::tryVar() const 
{ return power == 1 ? term.tryVar() : none<Variable>(); }

} // namespace Kernel

template<class NumTraits>
struct std::hash<Kernel::MonomFactor<NumTraits>> 
{
  size_t operator()(Kernel::MonomFactor<NumTraits> const& x) const noexcept 
  {
    return HashUtils::combine(
      StlHash::hash(x.term),
      StlHash::hash(x.power)
    );
  }
};


/////////////////////////////////////////////////////////
// impl MonomFactors  tempalte stuff
////////////////////////////

namespace Kernel {

template<class Number>
MonomFactors<Number>::MonomFactors(Stack<MonomFactor>&& factors) 
  : _factors(std::move(factors)) { }

template<class Number>
MonomFactors<Number>::MonomFactors() 
  : MonomFactors(decltype(_factors)()) { }

template<class Number>
MonomFactors<Number>::MonomFactors(PolyNf t) 
  : _factors { MonomFactor ( t, 1 ) }  { }

template<class Number>
MonomFactors<Number>::MonomFactors(PolyNf t1, PolyNf t2) 
  : _factors(t1 == t2 ? decltype(_factors) ({MonomFactor ( t1, 2 )  }) : 
             t1 <  t2 ? decltype(_factors) ({ MonomFactor ( t1, 1 ), MonomFactor ( t2, 1 ) }) 
                      : decltype(_factors) ({ MonomFactor ( t2, 1 ), MonomFactor ( t1, 1 ) }) 
                        )  { }

template<class Number>
unsigned MonomFactors<Number>::nFactors() const 
{ return _factors.size(); }

template<class Number>
MonomFactor<Number> & MonomFactors<Number>::factorAt(unsigned i) 
{ return _factors[i]; }

template<class Number>
MonomFactor<Number> const& MonomFactors<Number>::factorAt(unsigned i) const
{ return _factors[i]; }

template<class Number>
PolyNf const& MonomFactors<Number>::termAt(unsigned i) const
{ return _factors[i].term; }

template<class Number>
Stack<MonomFactor<Number>> & MonomFactors<Number>::raw()
{ return _factors; }

template<class Number>
bool MonomFactors<Number>::isPolynom() const
{ return nFactors() == 1 
    && _factors[0].power == 1 
    && _factors[0].term.template is<AnyPoly>(); }

template<class Number>
Perfect<Polynom<Number>> MonomFactors<Number>::asPolynom() const
{ 
  ASS(isPolynom());
  return _factors[0].term
    .template unwrap<AnyPoly>()
    .template unwrap<Perfect<Polynom>>(); 
}

template<class Number>
MonomFactors<Number> MonomFactors<Number>::one()
{ return MonomFactors(); }

template<class Number>
bool MonomFactors<Number>::isOne() const 
{ return _factors.begin() == _factors.end(); }

template<class Number>
std::ostream& operator<<(std::ostream& out, const MonomFactors<Number>& self) 
{
  out << "(";
  if (self._factors.size() == 0) {
    out << "MonomFactors()";
  } else {
    auto iter  = self._factors.begin();
    out << *iter;
    iter++;
    for (; iter != self._factors.end(); iter++) {
      out << " " << *iter;
    }
  }
  out << ")";
  return out;
}


template<class Number>
bool operator==(const MonomFactors<Number>& l, const MonomFactors<Number>& r) {
  return l._factors == r._factors;
}

template<class Number>
TermList MonomFactors<Number>::denormalize(TermList* results)  const
{
  if (_factors.size() == 0) {
    return Number::one();
  } else {

    auto powerTerm = [](TermList t, int pow) -> TermList {
      ASS(pow > 0)
      TermList out = t;
      for (int i = 1; i < pow; i++) {
        out = Number::mul(t,out);
      }
      return out;
    };

    TermList out = powerTerm(results[0], _factors[0].power);

    for (unsigned i = 1; i < nFactors(); i++) {
      out = Number::mul(out, powerTerm(results[i], _factors[i].power));
    }

    return out;
  }
}

template<class Number>
Option<Variable> MonomFactors<Number>::tryVar() const 
{
  using Opt = Option<Variable>;
  if (nFactors() == 1 ) {

    return _factors[0].tryVar();
  } else {
    return  Opt();
  }
}

template<class Number>
void MonomFactors<Number>::integrity() const 
{
#if VDEBUG
  if (_factors.size() > 0) {
    auto iter = this->_factors.begin();
    auto last = iter++;
    while (iter != _factors.end()) {
      ASS_REP(last->term <= iter->term, *this);
      last = iter++;
    }
  }
#endif
}

template<class Number>
MonomFactors<Number> MonomFactors<Number>::replaceTerms(PolyNf* simplifiedTerms) const 
{
  int offs = 0;
  MonomFactors out;
  out._factors.reserve(nFactors());

  for (auto& fac : _factors) {
    out._factors.push(MonomFactor(simplifiedTerms[offs++], fac.power));
  }

  return out;
}

} // namespace Kernel

template<class NumTraits>
struct std::hash<Kernel::MonomFactors<NumTraits>> 
{
  size_t operator()(Kernel::MonomFactors<NumTraits> const& x) const noexcept 
  {
    return StackHash<StlHash>::hash(x._factors);
  }
};


/////////////////////////////////////////////////////////
// impl Polynom tempalte stuff
////////////////////////////

namespace Kernel {


template<class Number>
Polynom<Number>::Polynom(Stack<Monom>&& summands) 
  : _summands(
      summands.isEmpty() 
        ? Stack<Monom>{Monom::zero()} 
        : std::move(summands)) 
{ }

template<class Number>
Polynom<Number>::Polynom(Monom m) 
  : Polynom(
      Stack<Monom>{m})
{  }

template<class Number>
Polynom<Number>::Polynom(Numeral numeral, PolyNf term) 
  : Polynom(Monom(numeral, perfect(MonomFactors(term))))
{  }

template<class Number>
Polynom<Number>::Polynom(PolyNf t) 
  : Polynom(Numeral(1), t) 
{  }

template<class Number>
Polynom<Number>::Polynom(Numeral constant) 
  : Polynom(Monom(constant, perfect(MonomFactors::one()))) 
{  }


template<class Number>
bool operator==(const Polynom<Number>& lhs, const Polynom<Number>& rhs)
{ return std::tie(lhs._summands) == std::tie(rhs._summands); }

template<class Number>
std::ostream& operator<<(std::ostream& out, const Polynom<Number>& self) {
  auto iter = self._summands.begin();
  out << "Poly(";
  if ( iter == self._summands.end() ) {
    out << "<empty>";
  } else {
    out << *iter;
    iter++;
    for (; iter != self._summands.end(); iter++) {
      out << " + " << *iter;
    }
  }
  out << ")";
  return out;
}



template<class Number>
Polynom<Number> Polynom<Number>::zero() 
{ 
  auto out = Polynom(Stack<Monom>{Monom::zero()}); 
  out.integrity();
  return std::move(out);
}

template<class Number>
Option<typename Number::ConstantType> Polynom<Number>::toNumber() const& 
{ 
  if ( _summands.size() == 0) {
    return Option<Numeral>(Numeral(0));

  } else if (_summands.size() == 1 && _summands[0].factors->nFactors() == 0) {
    return Option<Numeral>(_summands[0].numeral);

  } else {
    return Option<Numeral>();
  }
}

template<class Number>
bool Polynom<Number>::isNumber() const& 
{ 
  return toNumber().isSome();
}

template<class Number>
typename Number::ConstantType Polynom<Number>::unwrapNumber() const& 
{ return toNumber().unwrap(); }

template<class Number>
TermList Polynom<Number>::denormalize(TermList* results) const
{
  auto monomToTerm = [](Monom const& monom, TermList* t) -> TermList {
    auto c = TermList(theory->representConstant(monom.numeral));
    if (monom.factors->isOne()) {
      return c;
    } else {
      auto mon = monom.factors->denormalize(t);
      if (monom.numeral == Number::constant(1)) {
        return mon;
      } else if (monom.numeral == Number::constant(-1)) {
        return Number::minus(mon);
      } else {
        return Number::mul(c, mon);
      }
    }
  };

  if (_summands.size() == 0) {
    return Number::zero();
  } else {

    TermList out = monomToTerm(_summands[0], results);
    auto flatIdx = _summands[0].factors->nFactors();

    for (unsigned i = 1; i < nSummands(); i++) {
      auto& monom = _summands[i];
      out = Number::add(monomToTerm(monom, &results[flatIdx]), out);
      flatIdx += monom.factors->nFactors();
    }
    return out;
  }
}

template<class Number>
Stack<Monom<Number>>& Polynom<Number>::raw()
{ return _summands; }

template<class Number>
Polynom<Number> Polynom<Number>::replaceTerms(PolyNf* simplifiedTerms) const 
{
  int offs = 0;
  Stack<Monom> out;
  out.reserve(nSummands());

  for (auto& monom : _summands) {
    out.push(Monom(
          monom.numeral, 
          perfect(monom.factors->replaceTerms(&simplifiedTerms[offs]))));
    offs += monom.factors->nFactors();
  }
  return Polynom(std::move(out));
}


template<class Number>
unsigned Polynom<Number>::nSummands() const
{ return _summands.size(); }

template<class Number>
unsigned Polynom<Number>::nFactors(unsigned summand) const
{ return _summands[summand].factors->nFactors(); }

template<class Number>
Monom<Number> const& Polynom<Number>::summandAt(unsigned summand) const
{ return _summands[summand]; }

template<class Number>
Monom<Number>      & Polynom<Number>::summandAt(unsigned summand)
{ return _summands[summand]; }

template<class Number>
void Polynom<Number>::integrity() const {
#if VDEBUG
  ASS(_summands.size() > 0)
  if (_summands.size() > 0) {
    auto iter = this->_summands.begin();
    auto last = iter++;
    while (iter != _summands.end()) {
      // ASS_REP(std::less<Perfect<MonomFactors>>{}(last->factors, iter->factors), *this);
      ASS_REP(last->factors <= iter->factors, *this);
      iter->integrity();
      last = iter++;
    }
  } 
#endif
}


template<class Number> 
TermList Polynom<Number>::denormalize()  const
{ return PolyNf(AnyPoly(perfect(Polynom(*this)))).denormalize(); }

} // namespace Kernel

template<class NumTraits>
struct std::hash<Kernel::Polynom<NumTraits>> 
{
  size_t operator()(Kernel::Polynom<NumTraits> const& x) const noexcept 
  {
    using namespace Lib;
    using namespace Kernel;

    unsigned out = HashUtils::combine(0,0);
    for (auto c : x._summands) {
      out = HashUtils::combine(
        StlHash::hash(c.factors),
        StlHash::hash(c.numeral),
        out
      );
    }
    return out;
  }
};


#undef DEBUG
#endif // __POLYNOMIAL__H__
