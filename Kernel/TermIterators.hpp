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
 * @file TermIterators.hpp
 * Defines several iteratorn over terms.
 */

#ifndef __TermIterators__
#define __TermIterators__

#include "Forwards.hpp"

#include "Lib/Recycled.hpp"
#include "Lib/Stack.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/DHMultiset.hpp"

#include "Term.hpp"
#include "SortHelper.hpp"
#include "ApplicativeHelper.hpp"

namespace Kernel {

/**
 * Iterator that yields variables of specified
 * @b term in DFS left to right order.
 *
 * - This iterator returns sort variables
 * - If the sort of the returned variables is required, please
 *   use VariableIterator2 below, having read its documentation.
 */
class VariableIterator
: public IteratorCore<TermList>
{
public:
  DECL_ELEMENT_TYPE(TermList);
  VariableIterator() : _stack(8), _used(false) {}

  VariableIterator(const Term* term) : _stack(8), _used(false)
  {
    if(term->isLiteral() && static_cast<const Literal*>(term)->isTwoVarEquality()){
      _aux[0] = TermList::empty();
      _aux[1]=static_cast<const Literal*>(term)->twoVarEqSort();
      _stack.push(&_aux[1]);      
    }
    if(!term->shared() || !term->ground()) {
      _stack.push(term->args());
    }
  }

  VariableIterator(TermList t) : _stack(8), _used(false)
  {
    if(t.isVar()) {
      _aux[0] = TermList::empty();
      _aux[1]=t;
      _stack.push(&_aux[1]);
    }
    else {
      Term* term=t.term();
      if(!term->shared() || !term->ground()) {
	_stack.push(term->args());
      }
    }
  }

  void reset(const Term* term)
  {
    _stack.reset();
    _used = false;
    if(term->isLiteral() && static_cast<const Literal*>(term)->isTwoVarEquality()){
      _aux[0] = TermList::empty();
      _aux[1]=static_cast<const Literal*>(term)->twoVarEqSort();
      _stack.push(&_aux[1]);      
    }
    if(!term->shared() || !term->ground()) {
      _stack.push(term->args());
    }
  }

  void reset(TermList t)
  {
    _stack.reset();
    _used = false;
    if(t.isVar()) {
      _aux[0] = TermList::empty();
      _aux[1]=t;
      _stack.push(&_aux[1]);
    }
    else {
      Term* term=t.term();
      if(!term->shared() || !term->ground()) {
	_stack.push(term->args());
      }
    }
  }


  bool hasNext();
  /** Return the next variable
   * @warning hasNext() must have been called before */
  TermList next()
  {
    ASS(!_used);
    ASS(_stack.top()->isVar());
    _used=true;
    return *_stack.top();
  }
private:
  Stack<const TermList*> _stack;
  bool _used;
  TermList _aux[2];
};

struct VariableIteratorFn
{
  VirtualIterator<TermList> operator()(Term* t)
  {
    return vi( new VariableIterator(t) );
  }
  VirtualIterator<TermList> operator()(TermList t)
  {
    if(t.isVar()) {
      return pvi( getSingletonIterator(t) );
    }
    else {
      return (*this)(t.term());
    }
  }
};

struct OrdVarNumberExtractorFn
{
  unsigned operator()(TermList t)
  {
    ASS(t.isOrdinaryVar());

    return t.var();
  }
};


/**
 * Iterator that yields variables of @b term along with the 
 * types of the variables. Notes on the use off this iterator:
 *
 * - The iterator is NOT compatible with special terms
 * - The itertaor is NOT compatible with literals
 * - For situations where it can be guaranteed that a term is not special
 *   this iterator should be prefered to VariableIterator used in conjunction
 *   with SortHelper::collectVariableSorts as it is more efficient.
 */
class VariableWithSortIterator
: public IteratorCore<std::pair<TermList,TermList>>
{
public:

  VariableWithSortIterator(const Term* term) : _stack(8), _terms(8), _used(false)
  {
    ASS(!term->isLiteral());
    if(!term->shared() || !term->ground()) {
      _terms.push(term);
      _argNums.push(0);
      _stack.push(term->args());
    }
  }

  bool hasNext();
  /** Return the next variable
   * @warning hasNext() must have been called before */
  std::pair<TermList, TermList> next()
  {
    ASS(!_used);
    _used=true;
    return std::make_pair(*_stack.top(),  SortHelper::getArgSort(const_cast<Term*>(_terms.top()), _argNums.top()));
  }
private:
  Stack<const TermList*> _stack;
  Stack<const Term*> _terms;
  Stack<unsigned> _argNums;
  bool _used;
};

/**
 * Iterator that yields proper subterms
 * of @b term in DFS left to right order.
 */
class SubtermIterator
  : public IteratorCore<TermList>
{
public:
  SubtermIterator(const Term* term) : _used(false)
  {
    pushNext(term->args());
  }

  bool hasNext();
  /** Return next subterm
   * @warning hasNext() must have been called before */
  TermList next()
  {
    ASS(!_used && !_stack->isEmpty());
    _used=true;
    return *_stack->top();
  }

  /**
   * Do not iterate subterms of the recently returned term (i.e. go right
   * rather than down in the term).
   */
  void right();
protected:
  SubtermIterator() : _used(false)
  { }

  inline
  void pushNext(const TermList* t)
  {
    if(!t->isEmpty()) {
      _stack->push(t);
    }
  }

  Recycled<Stack<const TermList*>> _stack;
  bool _used;
};

//////////////////////////////////////////////////////////////////////////
///                                                                    ///
///            ITERATORS REQUIRED FOR HIGHER-ORDER REASONING           ///
///                                                                    ///
//////////////////////////////////////////////////////////////////////////

/**
 * Iterator that yields the arguments of @b
 * applicative term in left-right order. Rather non-
 * standardly, the iterator also has a function head() 
 * that can be called at any time to return the head 
 * of @b applicative term
 */
class ApplicativeArgsIt
  : public IteratorCore<TermList>
{
public:
  ApplicativeArgsIt(const TermList term, bool returnTypeArgs = true)
  {
    if(returnTypeArgs){
      ApplicativeHelper::getHeadAndAllArgs(term, _head, _stack);
    } else {
      ApplicativeHelper::getHeadSortAndArgs(term, _head, _headSort, _stack);      
    }
    _argNum = _stack.size();
  }

  bool hasNext(){
    return !_stack.isEmpty();
  }
  /** Return next arg of _head
   * @warning hasNext() must have been called before */
  TermList next()
  {
    ASS(!_stack.isEmpty());
    return _stack.pop();
  }

  TermList head(){
    return _head ;
  }

  TermList headSort(){
    return _headSort;
  }

  unsigned argNum(){
    return _argNum;
  }

  bool isVar(){
    return _head.isVar() && !_argNum;
  }

  bool varHead(){
    return _head.isVar();
  }

  unsigned headNum(){
    if(_head.isVar()){
      return _head.var();
    }
    return _head.term()->functor();
  }

protected:
  TermStack _stack;
  TermList _head;
  TermList _headSort;
  unsigned _argNum;
};

class TopLevelVarLikeTermIterator
  : public IteratorCore<Term*>
{
public:
  TopLevelVarLikeTermIterator(Term* term)
  {
    _next = 0;
    if(term->isApplication() && !term->ground()){
      _stack.push(term);
    }
  }

  bool hasNext();

  Term* next()
  {
    ASS(_next);
    Term* res = _next;
    _next = 0;
    return res;
  }

private:
  Stack<Term*> _stack;
  Term* _next;
};

class TopLevelVarIterator
  : public IteratorCore<TermList>
{
public:
  TopLevelVarIterator(TermList t);
  
  bool hasNext();

  TermList next()
  {
    ASS(!_next.isEmpty());
    TermList res = _next;
    _next = TermList::empty();
    return res;
  }

private:
  TermStack _stack;
  TermList _next;
};


class RewritableVarsIt
{
public: //includeSelf for compatibility
  DECL_ELEMENT_TYPE(TypedTermList);
  RewritableVarsIt(DHSet<unsigned>* unstableVars, Term* t, bool includeSelf = false) :  _next(), _stack(8)
  {
    _unstableVars = unstableVars;
    if(t->isLiteral()){
      TermList t0 = *t->nthArgument(0);
      TermList t1 = *t->nthArgument(1);
      if(!t0.isVar()){ 
        _stack.push(t0);
        _sorts.push(SortHelper::getResultSort(t0.term()));
      }
      if(!t1.isVar()){ 
        _stack.push(t1); 
        _sorts.push(SortHelper::getResultSort(t1.term()));
      }      
      return;      
    }     
    _stack.push(TermList(t));
    _sorts.push(SortHelper::getResultSort(t));
  }

  bool hasNext();
  TypedTermList next(){
    ASS(_next.isSome());
    ASS(_next->isVar());
    return *_next.take();
  }
private:
  Option<TypedTermList> _next;
  Stack<TermList> _stack;
  Stack<TermList> _sorts;
  DHSet<unsigned>* _unstableVars;
};

class UnstableVarIt
  : public IteratorCore<TermList>
{
public: 
  UnstableVarIt(Term* t) : _stable(8), _stack(8)
  {
    _next = TermList::empty();
    if(t->isLiteral()){
      _stack.push(*t->nthArgument(0));
      _stack.push(*t->nthArgument(1));
      _stable.push(true);
      _stable.push(true);   
      return;      
    }
    _stable.push(true);
    _stack.push(TermList(t)); 
  }
  
  bool hasNext();
  TermList next()
  {
    ASS(!_next.isEmpty());
    TermList res = _next;
    _next = TermList::empty();
    return res;
  }

private:
  TermList _next;
  Stack<bool> _stable;
  Stack<TermList> _stack;
};

class FirstOrderSubtermIt
: public IteratorCore<Term*>
{
public:
  FirstOrderSubtermIt(Term* term, bool includeSelf=false) 
  : _stack(8), _added(0)
  {
    if(term->isLiteral()){
      TermList t0 = *term->nthArgument(0);
      TermList t1 = *term->nthArgument(1);
      if(!t0.isVar()){ _stack.push(t0.term()); }
      if(!t1.isVar()){ _stack.push(t1.term()); }      
      return;      
    } 
    _stack.push(term);
    if (!includeSelf) {
      FirstOrderSubtermIt::next();
    }
  }

  bool hasNext(){ return !_stack.isEmpty(); }
  Term* next();
  void right();

private:
  Stack<Term*> _stack;
  int _added;
};


class NarrowableSubtermIt
: public IteratorCore<TermList>
{
public:
  NarrowableSubtermIt(Term* term, bool includeSelf=false) 
  : _used(true), _stack(8)
  {
    if(term->isLiteral()){
      TermList t0 = *term->nthArgument(0);
      TermList t1 = *term->nthArgument(1);
      if(!t0.isVar()){ _stack.push(t0.term()); }
      if(!t1.isVar()){ _stack.push(t1.term()); }      
      return;      
    } 
    _stack.push(term);
    //TODO
  }

  bool hasNext();
  TermList next(){
    ASS(!_used);
    _used = true;
    return _next;
  }

private:
  bool _used;
  TermList _next;
  Stack<Term*> _stack;
};

/*
 *  Returns Boolean subterms of a term.
 */
class BooleanSubtermIt
: public IteratorCore<TermList>
{
public:
  BooleanSubtermIt(Term* term, bool includeSelf=false) 
  : _used(true), _stack(8)
  {
    if(term->isLiteral()){
      TermList t0 = *term->nthArgument(0);
      TermList t1 = *term->nthArgument(1);
      if(!t0.isVar()){ _stack.push(t0.term()); }
      if(!t1.isVar()){ _stack.push(t1.term()); }      
      return;      
    } 
    _stack.push(term);
  }

  bool hasNext();
  TermList next(){
    ASS(!_used);
    _used = true;
    return _next;
  }

private:
  bool _used;
  TermList _next;
  Stack<Term*> _stack;
};

//////////////////////////////////////////////////////////////////////////
///                                                                    ///
///                    END OF HIGHER-ORDER ITERATORS                   ///
///                                                                    ///
//////////////////////////////////////////////////////////////////////////

/**
 * Iterator that yields proper subterms of commutative
 * literal @b lit in DFS left to right order with the
 * arguments of the literal reversed.
 */
class ReversedCommutativeSubtermIterator
: public SubtermIterator
{
public:
  ReversedCommutativeSubtermIterator(const Term* trm)
  {
    ASS(trm->commutative());
    ASS_EQ(trm->arity(),2);

    aux[0]=*trm->nthArgument(1);
    aux[1] = TermList::empty();
    aux[2]=*trm->nthArgument(0);
    aux[3] = TermList::empty();

    _stack->push(&aux[0]);
    _stack->push(&aux[2]);
  }
private:
  TermList aux[4];
};

/**
 * Iterator that yields proper subterms
 * of specified @b term, so that for each function it first yields
 * its arguments left to right, and then the function itself.
 */
class PolishSubtermIterator
: public IteratorCore<TermList>
{
public:
  PolishSubtermIterator(const Term* term) : _stack(8), _used(false)
  {
    pushNext(term->args());
  }

  bool hasNext();
  /** Return next subterm
   * @warning hasNext() must have been called before */
  TermList next()
  {
    ASS(!_used && !_stack.isEmpty());
    _used=true;
    return *_stack.top();
  }
private:
  inline
  void pushNext(const TermList* t)
  {
    while(!t->isEmpty()) {
      _stack.push(t);
      if(!t->isTerm()) {
	return;
      }
      t=t->term()->args();
    }
  }
  Stack<const TermList*> _stack;
  bool _used;
};

/**
 * Iterator that yields non-variable subterms
 * of specified @b term in DFS left to right order.
 *
 * - For polymorphic terms, this iterator returns both type and term subterms
 * - For monomorphic terms, this iterator and the one below behave identically
 * - This iterator should be used in circumstances where all non-variable subterms
 *   are required. 
 */
class NonVariableIterator
  : public IteratorCore<TermList>
{
public:
  NonVariableIterator(const NonVariableIterator&);
  /**
   * Create an iterator. If @c includeSelf is false, then only proper subterms
   * of @c term will be included.
   * @since 04/05/2013 Manchester, argument includeSelf added
   * @author Andrei Voronkov
   */
  NonVariableIterator(Term* term,bool includeSelf=false)
  : _stack(8),
    _added(0)
  {
    _stack.push(term);
    if (!includeSelf) {
      NonVariableIterator::next();
    }
  }
  // NonVariableIterator(TermList ts);

  /** true if there exists at least one subterm */
  bool hasNext() { return !_stack.isEmpty(); }
  TermList next();
  void right();
private:
  /** available non-variable subterms */
  Stack<Term*> _stack;
  /** the number of non-variable subterms added at the last iteration, used by right() */
  int _added;
}; // NonVariableIterator


/**
 * Iterator that yields non-variable subterms that are not type arguments
 * or subterms of type arguments of @b term in DFS left to right order.
 *
 * - This iterator should be used when only subterms of TERM arguments are required.
 *   Typical usecases are:
 *   - Subterms to be rewritten by an inference. See for example superposition,
 *     backward and forward demodulation.
 *   - Iterating through subterms to see whether a term occurs as a subterm of 
 *     another
 */
class NonVariableNonTypeIterator
  : public IteratorCore<Term*>
{
public:
  NonVariableNonTypeIterator(const NonVariableNonTypeIterator&);
  /**
   * If @c includeSelf is false, then only proper subterms of @c term will be included.
   */
  NonVariableNonTypeIterator(Term* term, bool includeSelf=false)
  : _stack(8),
    _added(0)
  {
    _stack.push(term);
    if (!includeSelf) {
      NonVariableNonTypeIterator::next();
    }
  }
  // NonVariableIterator(TermList ts);

  /** true if there exists at least one subterm */
  bool hasNext() { return !_stack.isEmpty(); }
  Term* next();
  void right();
private:
  /** available non-variable subterms */
  Stack<Term*> _stack;
  /** the number of non-variable subterms added at the last iteration, used by right() */
  int _added;
}; // NonVariableIterator


class TracedNonVariableNonTypeIterator
  : public IteratorCore<Term*>
{
public:
  TracedNonVariableNonTypeIterator(const TracedNonVariableNonTypeIterator&);
  /**
   * If @c includeSelf is false, then only proper subterms of @c term will be included.
   */
  TracedNonVariableNonTypeIterator(Literal* lit)
  : _stack(8),
    _ready(false)
  {
    for (unsigned i = 0; i <= 1; i++) {
      auto t = lit->termArg(i);
      if (t.isTerm()) {
        _stack.push(std::make_pair(t.term(),0));
      }
    }
    _ready = true;
  }

  /** true if there exists at least one subterm */
  bool hasNext();
  Term* next();
  void right();
  const Stack<Term*>& getTrace() const { return _trace; }
private:
  /** available non-variable subterms */
  Stack<std::pair<Term*,unsigned>> _stack;
  Stack<Term*> _trace;
  /** the number of non-variable subterms added at the last iteration, used by right() */
  bool _ready;
}; // TracedNonVariableNonTypeIterator

/**
 * Iterator that iterator over disagreement set of two terms
 * or literals in DFS left to right order.
 */
class DisagreementSetIterator
: public IteratorCore<std::pair<TermList, TermList> >
{
public:
  /**
   * Create an empty disagreement iterator
   *
   * In order to be used, it must be reset by the @b reset function
   */
  DisagreementSetIterator()
  {
    _arg1 = TermList::empty();
  }

  /**
   * Create an iterator over the disagreement set of two terms
   */
  DisagreementSetIterator(TermList t1, TermList t2, bool disjunctVariables=true)
  : _stack(8)
  {
    reset(t1, t2, disjunctVariables);
  }
  /**
   * Create an iterator over the disagreement set of two terms/literals
   * with the same top functor
   */
  DisagreementSetIterator(Term* t1, Term* t2, bool disjunctVariables=true)
  : _stack(8), _disjunctVariables(disjunctVariables)
  {
    reset(t1,t2,disjunctVariables);
  }

  void reset(TermList t1, TermList t2, bool disjunctVariables=true)
  {
    ASS(!t1.isEmpty());
    ASS(!t2.isEmpty());

    _stack.reset();
    _disjunctVariables=disjunctVariables;
    if(!TermList::sameTop(t1,t2) || (t1.isVar() && disjunctVariables)) {
      _arg1=t1;
      _arg2=t2;
      return;
    }
    _arg1 = TermList::empty();
    if(t1.isTerm() && t1.term()->arity()>0) {
      _stack.push(t1.term()->args());
      _stack.push(t2.term()->args());
    }
  }

  void reset(Term* t1, Term* t2, bool disjunctVariables=true)
  {
    ASS_EQ(t1->functor(), t2->functor());

    _stack.reset();
    _disjunctVariables=disjunctVariables;

    _arg1 = TermList::empty();
    if((t1->isLiteral() && static_cast<Literal*>(t1)->isTwoVarEquality()) ||
       (t2->isLiteral() && static_cast<Literal*>(t2)->isTwoVarEquality())){
      TermList s1 = SortHelper::getEqualityArgumentSort(static_cast<Literal*>(t1));
      TermList s2 = SortHelper::getEqualityArgumentSort(static_cast<Literal*>(t2));
      if(!TermList::sameTop(s1,s2) || (s1.isVar() && disjunctVariables)) {
        _arg1=s1;
        _arg2=s2;
      } else if(s1.isTerm() && s1.term()->arity()>0) {
        _stack.push(s1.term()->args());
        _stack.push(s2.term()->args());
      }
    }
    if(t1->arity()>0) {
      _stack.push(t1->args());
      _stack.push(t2->args());
    }
  }

  bool hasNext();

  /** Return next subterm
   * @warning hasNext() must have been called before */
  std::pair<TermList, TermList> next()
  {
    std::pair<TermList, TermList> res(_arg1,_arg2);
    _arg1 = TermList::empty();
    return res;
  }
private:
  Stack<TermList*> _stack;
  bool _disjunctVariables;
  TermList _arg1;
  TermList _arg2;
};



/**
 * Implements an iterator over function symbols of a term.
 *
 * Functions are yielded before their subterms.
 * @since 26/05/2007 Manchester, made from class TermVarIterator
 */
class TermFunIterator
: public IteratorCore<unsigned>
{
public:
  TermFunIterator (const Term*);

  bool hasNext();
  unsigned next();
private:
  /** True if the next symbol is found */
  bool _hasNext;
  /** next symbol, previously found */
  unsigned _next;
  /** Stack of term lists (not terms!) */
  Stack<const TermList*> _stack;
}; // class TermFunIterator


/**
 * Implements an iterator over variables of a term term list, or atom.
 * @since 06/01/2004 Manchester
 * @since 26/05/2007 Manchester, reimplemented for different data structures
 */
class TermVarIterator
: public IteratorCore<unsigned>
{
public:
  TermVarIterator (const Term*);
  TermVarIterator (const TermList*);

  bool hasNext ();
  unsigned next ();
private:
  /** True if the next variable is found */
  // bool _hasNext; // MS: unused
  /** next variable, previously found */
  unsigned _next;
  /** Stack of term lists (not terms!) */
  Stack<const TermList*> _stack;
}; // class TermVarIterator


class LiteralArgIterator 
{
  Literal* _lit;
  unsigned _idx;
public:
  DECL_ELEMENT_TYPE(TypedTermList);

  LiteralArgIterator(Literal* lit) : _lit(lit), _idx(0) {}

  inline bool hasNext() const { return _idx < _lit->arity(); }
  inline TermList next() { return TypedTermList(*_lit->nthArgument(_idx), SortHelper::getArgSort(_lit, _idx)); _idx++; }
  unsigned size() const { return _lit->arity(); }
};


/** iterator over all term arguments of @code term */
static const auto termArgIter = [](Term* term) 
  { return iterTraits(getRangeIterator<unsigned>(0, term->numTermArguments()))
      .map([=](auto i)
           { return term->termArg(i); }); };

/** iterator over all type arguments of @code term */
static const auto typeArgIter = [](Term* term) 
  { return iterTraits(getRangeIterator<unsigned>(0, term->numTypeArguments()))
      .map([=](auto i)
           { return term->typeArg(i); }); };

/** iterator over all type and term arguments of @code term */
static const auto anyArgIter = [](Term* term) 
  { return iterTraits(getRangeIterator<unsigned>(0, term->arity()))
      .map([=](auto i)
           { return *term->nthArgument(i); }); };


} // namespace Kernel

#endif // __TermIterators__
