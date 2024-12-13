/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */


#include "Kernel/Polynomial.hpp"
#include "Kernel/PolynomialNormalizer.hpp"
#include "Lib/Output.hpp"

#define DEBUG(...) // DBG(__VA_ARGS__)

namespace Kernel {

PolyNf PolyNf::normalize(TypedTermList t)  
{ return normalizeTerm(t); }

////////////////////////////////////////////////////////////////////////////////////////////////////
// impl Variable
/////////////////////////////////////////////////////////

Variable::Variable() : _num() {}

Variable::Variable(unsigned num) : _num(num) {}

unsigned Variable::id() const 
{ return _num; }


bool operator==(Variable lhs, Variable rhs) 
{ return lhs._num == rhs._num; }

bool operator!=(Variable lhs, Variable rhs) 
{ return !(lhs == rhs); }

bool operator<(Variable const& lhs, Variable const& rhs)
{ return lhs._num < rhs._num; }

std::ostream& operator<<(std::ostream& out, const Variable& self) 
{ return out << "X" << self._num; }

////////////////////////////////////////////////////////////////////////////////////////////////////
// impl FuncId
/////////////////////////////////////////////////////////

FuncId::FuncId(unsigned num, const TermList* typeArgs) : _num(num), _typeArgs(typeArgs) {}

FuncId FuncId::symbolOf(Term* term) 
{ return FuncId(term->functor(), term->typeArgs()); }

unsigned FuncId::numTermArguments() 
{ return symbol()->numTermArguments(); }

std::ostream& operator<<(std::ostream& out, const FuncId& self) 
{ 
  if (self.numTypeArgs() == 0) {
    return out << self.symbol()->name(); 
  } else {
    return out << self.symbol()->name() << "<" << Output::interleaved(", ", self.iterTypeArgs())  << ">"; 
  }
}

Signature::Symbol* FuncId::symbol() const 
{ return env.signature->getFunction(_num); }

unsigned FuncId::id() const 
{ return _num; }

Theory::Interpretation FuncId::interpretation() const 
{ return theory->interpretFunction(_num); }

bool FuncId::isInterpreted() const
{ return theory->isInterpretedFunction(_num); }

Option<Theory::Interpretation> FuncId::tryInterpret() const
{ 
  return isInterpreted() ? some<Theory::Interpretation>(interpretation())
                         : none<Theory::Interpretation>();
}
} // namespace Kernel

////////////////////////////////////////////////////////////////////////////////////////////////////
// impl FuncTerm
/////////////////////////////////////////////////////////

namespace Kernel {

FuncTerm::FuncTerm(FuncId f, Stack<PolyNf>&& args) 
  : _fun(f)
  , _args(std::move(args)) 
{ }

FuncTerm::FuncTerm(FuncId f, PolyNf* args) 
  : _fun(f)
  , _args(Stack<PolyNf>::fromIterator(arrayIter(args, f.numTermArguments()))) 
{ }

bool operator==(FuncTerm const& lhs, FuncTerm const& rhs) 
{ return lhs._fun == rhs._fun && lhs._args == rhs._args; }

bool operator!=(FuncTerm const& lhs, FuncTerm const& rhs) 
{ return !(lhs == rhs); }

unsigned FuncTerm::numTermArguments() const 
{ return _args.size(); }

FuncId FuncTerm::function() const 
{ return _fun; }

PolyNf const& FuncTerm::arg(unsigned i) const 
{ return _args[i]; }

std::ostream& operator<<(std::ostream& out, const FuncTerm& self) 
{ 
  out << self._fun;
  auto& stack = self._args;
  auto iter = stack.iterFifo();

  if (iter.hasNext()) {
    out << "(" << iter.next();
    while (iter.hasNext()) {
      out << ", " << iter.next();
    }
    out << ")";
  }

  return out;
}


/////////////////////////////////////////////////////////
// impl AnyPoly 
////////////////////////////


AnyPoly AnyPoly::replaceTerms(PolyNf* newTs) const 
{ return apply([&](auto& t) -> decltype(auto) { return AnyPoly(perfect(t->replaceTerms(newTs))); }); }

TermList AnyPoly::denormalize(TermList* results) const
{ return apply([&](auto& t) -> decltype(auto) { return t->denormalize(results); }); }

PolyNf const& AnyPoly::termAt(unsigned summand, unsigned factor)  const
{ return apply([&](auto& t) -> decltype(auto) { return t->summandAt(summand).factors->termAt(factor); }); }

unsigned AnyPoly::nSummands() const 
{ return apply([&](auto& t) -> decltype(auto) { return t->nSummands(); }); }

unsigned AnyPoly::nFactors(unsigned i) const 
{ return apply([&](auto& t) -> decltype(auto) { return t->nFactors(i); }); }

std::ostream& operator<<(std::ostream& out, const AnyPoly& self) 
{ return self.apply([&](auto& t) -> decltype(auto) { return out << *t; }); }


/////////////////////////////////////////////////////////
// impl PolyNf 
////////////////////////////

PolyNf::PolyNf(Perfect<FuncTerm> t) : Coproduct(t) {}
PolyNf::PolyNf(Variable          t) : Coproduct(t) {}
PolyNf::PolyNf(AnyPoly           t) : Coproduct(t) {}


bool operator==(PolyNf const& lhs, PolyNf const& rhs) 
{ return static_cast<PolyNfSuper const&>(lhs) == static_cast<PolyNfSuper const&>(rhs); }

bool operator!=(PolyNf const& lhs, PolyNf const& rhs) 
{ return !(lhs == rhs); }

std::ostream& operator<<(std::ostream& out, const PolyNf& self)
{ return self.apply([&](auto& t) -> decltype(auto) { return out << t; }); }

Option<Variable> PolyNf::tryVar() const 
{ return as<Variable>().toOwned(); }

IterTraits<PolyNf::SubtermIter> PolyNf::iterSubterms() const 
{ return iterTraits(SubtermIter(*this)); }

bool operator<(const PolyNf& lhs, const PolyNf& rhs) 
{ return std::less<PolyNfSuper>{}(lhs,rhs); }

bool operator<=(const PolyNf& lhs, const PolyNf& rhs) 
{ return lhs < rhs || lhs == rhs; }

PolyNf::SubtermIter::SubtermIter(PolyNf p) 
  : _stack(decltype(_stack){ BottomUpChildIter<PolyNf>(p) }) 
{  }

PolyNf PolyNf::SubtermIter::next() 
{
  ASS(_stack.size() != 0)
  while(_stack.top().hasNext()) {
    ASS(_stack.size() != 0)
    _stack.push(BottomUpChildIter<PolyNf>(_stack.top().next()));
  }
  ASS(_stack.size() != 0)
  return _stack.pop().self();
}

bool PolyNf::SubtermIter::hasNext() const 
{ 
  return !_stack.isEmpty(); 
}


/////////////////////////////////////////////////////////
// impl IterArgsPnf
////////////////////////////


IterArgsPnf::IterArgsPnf(Literal* lit) 
  : _lit(lit)
  , _idx(0) 
{ }

bool IterArgsPnf::hasNext() const  
{ return _idx < _lit->numTermArguments();  }

PolyNf IterArgsPnf::next()
{ 
  auto out = PolyNf::normalize(TypedTermList(_lit->termArg(_idx), SortHelper::getTermArgSort(_lit, _idx)));
  _idx++;
  return out;
}

IterTraits<IterArgsPnf> iterArgsPnf(Literal* lit) 
{ return iterTraits(IterArgsPnf(lit)); }

} // namespace Kernel
