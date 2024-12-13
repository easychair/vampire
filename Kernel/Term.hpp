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
 * @file Term.hpp
 * Defines class Term (also serving as term arguments)
 *
 * The way terms are laid out in memory is partly historical and certainly non-trivial.
 * Here are a few salient points to help you navigate:
 * - a "Term" represents a function (unsigned, see Kernel::Signature) applied to some number of arguments
 * - usually Terms are "perfectly shared" (see Indexing::TermSharing)
 * - the arguments are "TermList"s, i.e. a variable or a Term*
 * - TermList is a tagged union that relies on Term* being aligned (!) to achieve pointer tagging
 * - TermTag::REF == 0 because this does not change the value of an aligned Term*
 * - the arguments of a Term are laid out in reverse order (!)
 * - the last argument (i.e. the closest to the Term) is a sentinel with TermTag::FUN
 *   and some metadata about the enclosing Term
 *
 * @since 18/04/2006 Bellevue
 * @since 06/05/2007 Manchester, changed into a single class instead of three
 */

#ifndef __Term__
#define __Term__

#include "Lib/Output.hpp"
#include "Forwards.hpp"
#include "Debug/Assertion.hpp"

#include "Lib/Allocator.hpp"
#include "Lib/BitUtils.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/Portability.hpp"
#include "Lib/Comparison.hpp"
#include "Lib/Reflection.hpp"
#include "Lib/Sort.hpp"
#include "Lib/Stack.hpp"
#include "Lib/Hash.hpp"
#include "Lib/Coproduct.hpp"
#include "Lib/Recycled.hpp"

// the number of bits used for "TermList::_info::distinctVars"
#define TERM_DIST_VAR_BITS 22
// maximum number that fits in a TERM_DIST_VAR_BITS-bit unsigned integer
#define TERM_DIST_VAR_UNKNOWN ((1 << TERM_DIST_VAR_BITS)-1)

namespace Kernel {
  std::ostream& operator<<(std::ostream& out, Term const& self);
  std::ostream& operator<<(std::ostream& out, TermList const& self);
  std::ostream& operator<<(std::ostream& out, Literal const& self);
  bool operator<(TermList const&,TermList const&);

using namespace Lib;

/** Tag denoting the kind of this term
 * @since 19/02/2008 Manchester, moved outside of the Term class
 */
enum TermTag {
  /** reference to another term */
  REF = 0u,
  /** ordinary variable */
  ORD_VAR = 1u,
  /** (function) symbol */
  FUN = 2u,
  /** special variable */
  SPEC_VAR = 3u,
};

enum ArgumentOrderVals {
  /**
   * Values representing order of arguments in equality,
   * to be stores in the term sharing structure.
   *
   * The important thing is that the UNKNOWN value is
   * equal to 0, as this will be the default value inside
   * the term objects
   *
   * Values of elements must be equal to values of corresponding elements
   * in the @c Result enum, so that one can convert between the
   * enums using static_cast.
   */
  AO_UNKNOWN=0,
  AO_GREATER=1,
  AO_LESS=2,
  AO_EQUAL=3,
  AO_INCOMPARABLE=4,
};

inline std::ostream& operator<<(std::ostream& out, ArgumentOrderVals const& self)
{ 
  switch(self) {
    case AO_UNKNOWN: return out << "UNKNOWN";
    case AO_GREATER: return out << "GREATER";
    case AO_LESS: return out << "LESS";
    case AO_EQUAL: return out << "EQUAL";
    case AO_INCOMPARABLE: return out << "INCOMPARABLE";
  }
  ASSERTION_VIOLATION
}

enum class TermKind : unsigned {
  LITERAL,
  TERM,
  SORT,
};

/* a function symbol of a composite term. in addition to the function symbol id (the functor in vampire terminology) we store what kind of term (i.e. term, literal or sort) it is. */
struct SymbolId {
  unsigned functor;
  TermKind kind;
  friend bool operator==(SymbolId const& l, SymbolId const& r) 
  { return std::tie(l.functor, l.kind) == std::tie(r.functor, r.kind); }
  friend bool operator<(SymbolId const& l, SymbolId const& r) 
  { return std::tie(l.functor, l.kind) < std::tie(r.functor, r.kind); }
};

/**
 * Class containing either a pointer to a compound term or
 * a variable number or a functor.
 */
class TermList {
public:
  /* default constructor, satisfying isEmpty() */
  TermList() : _content(FUN) {}
  /** creates a term list containing a pointer to a term */
  explicit TermList(Term* t) : _content(0) {
    // NB we also zero-initialise _content so that the spare bits are zero on 32-bit platforms
    // dead-store eliminated on 64-bit
    _setTerm(t);
    ASS_EQ(tag(), REF);
  }
  /** creates a term list containing a variable. If @b special is true, then the variable
   * will be "special". Special variables are also variables, with a difference that a
   * special variables and ordinary variables have an empty intersection */
  TermList(unsigned var, bool special)
  {
    if (special) {
      makeSpecialVar(var);
    }
    else {
      makeVar(var);
    }
  }

  /** the tag */
  inline TermTag tag() const { return static_cast<TermTag>(_tag()); }
  /** the term list is empty */
  inline bool isEmpty() const
  { return tag() == FUN; }
  /** the term list is non-empty */
  inline bool isNonEmpty() const
  { return tag() != FUN; }
  /** next term in this list */
  inline TermList* next()
  { return this-1; }
  /** next term in this list */
  inline const TermList* next() const
  { return this-1; }
  /** the term contains a variable as its head */
  inline bool isVar() const { return tag() == ORD_VAR || tag() == SPEC_VAR; }
  /** the term contains an ordinary variable as its head */
  inline bool isOrdinaryVar() const { return tag() == ORD_VAR; }
  /** the term contains a special variable as its head */
  inline bool isSpecialVar() const { return tag() == SPEC_VAR; }

  /** return the variable number */
  inline unsigned var() const
  { ASS(isVar()); return _content / 4; }
  /** the term contains reference to Term class */
  inline bool isTerm() const
  { return tag() == REF; }
  inline const Term* term() const
  { ASS(isTerm()); return _term(); }
  inline Term* term()
  { ASS(isTerm()); return _term(); }
  /** True of the terms have the same content. Useful for comparing
   * arguments of shared terms. */
  inline bool sameContent(const TermList* t) const
  { return _content == t->_content ; }
  inline bool sameContent(const TermList& t) const
  { return sameContent(&t); }
  /** return the content, useful for e.g., term argument comparison */
  inline uint64_t content() const { return _content; }
  /** set the content manually - hazardous, such terms should then only be used as integers */
  void setContent(uint64_t content) { _content = content; }
  /** default hash is to hash the content */
  unsigned defaultHash() const { return DefaultHash::hash(content()); }
  unsigned defaultHash2() const { return content(); }

  std::string toString(bool needsPar = false) const;

  friend std::ostream& operator<<(std::ostream& out, Kernel::TermList const& tl);
  /** make the term into an ordinary variable with a given number */
  inline void makeVar(unsigned vnumber)
  { _content = vnumber * 4 + ORD_VAR; }
  /** make the term into a special variable with a given number */
  inline void makeSpecialVar(unsigned vnumber)
  { _content = vnumber * 4 + SPEC_VAR; }
  /** create an term empty (so that isEmpty() returns true)
   *  (can just be the default constructor now)
   */
  inline static TermList empty()
  { return TermList(); }

  /** the top of a term is either a function symbol or a variable id. this class is model this */
  class Top {
    using Inner = Coproduct<unsigned, SymbolId>;
    static constexpr unsigned VAR = 0;
    static constexpr unsigned FUN = 1;
    Inner _inner;
    
    Top(Inner self) : _inner(self) {}
  public:
    static Top var    (unsigned v) { return Top(Inner::variant<VAR>(v)); }
    // static Top functor(unsigned f) { return Top(Inner::variant<FUN>(f)); }
    template<class T>
    static Top functor(T const* t) { return Top(Inner::variant<FUN>({ t->functor(), t->kind(), })); }
    Option<unsigned> var()     const { return _inner.as<VAR>().toOwned(); }
    Option<SymbolId> functor() const { return _inner.as<FUN>().toOwned(); }
    Lib::Comparison compare(Top const& other) const 
    { return _inner.compare(other._inner); }
    IMPL_COMPARISONS_FROM_COMPARE(Top);
    friend bool operator==(Top const& l, Top const& r) { return l._inner == r._inner; }
    friend bool operator!=(Top const& l, Top const& r) { return      !(l == r);       }
    void output(std::ostream& out) const;

    friend std::ostream& operator<<(std::ostream& out, Kernel::TermList::Top const& self)
    { self.output(out); return out; }
  };

  /* returns the Top of a function (a variable id, or a function symbol depending on whether the term is a variable or a complex term) */
  Top top() const
  { return isTerm() ? TermList::Top::functor(term()) 
                    : TermList::Top::var(var());            }

  /** make the term into a reference */
  inline void setTerm(Term* t) {
    // NB we also zero-initialise _content so that the spare bits are zero on 32-bit platforms
    // dead-store eliminated on 64-bit
    _content = 0;
    _setTerm(t);
    ASS_EQ(tag(), REF);
  }
  static bool sameTop(TermList ss, TermList tt);
  static bool sameTopFunctor(TermList ss, TermList tt);
  static bool equals(TermList t1, TermList t2);
  static bool allShared(TermList* args);
  static TermList var(unsigned var, bool special = false) { return TermList(var, special); }
  /** if not var, the inner term must be shared */
  unsigned weight() const;
  /** returns true if this termList is wrapping a higher-order "arrow" sort */
  bool isArrowSort();
  bool isBoolSort();
  bool isArraySort();
  bool isTupleSort();
  bool isApplication() const;
  bool containsSubterm(TermList v) const;
  bool containsAllVariablesOf(TermList t) const;
  bool ground() const;
  bool isSafe() const;

#if VDEBUG
  void assertValid() const;
#endif

  inline bool operator==(const TermList& t) const
  { return _content==t._content; }
  inline bool operator!=(const TermList& t) const
  { return _content!=t._content; }

  friend bool operator<(const TermList& lhs, const TermList& rhs);

private:
  std::string asArgsToString() const;

  // the actual content of a TermList
  // this packs several things in:
  // 1. a Term *
  // 2. metadata (see below) such that _tag() is the lowest two bits of (1)
  // 3. "other", rarely used and handled specially
  uint64_t _content;

  // metadata used to be defined as this bitfield:
#if 0
  struct {
    /** a TermTag indicating what is stored here */
    unsigned tag : 2;
    /** polarity, used only for literals */
    unsigned polarity : 1;
    /** true if commutative/symmetric */
    unsigned commutative : 1;
    /** true if shared */
    unsigned shared : 1;
    /** true if literal */
    unsigned literal : 1;
    /** true if atomic sort */
    unsigned sort : 1;
    /** true if term contains at least one term var */
    unsigned hasTermVar : 1;
    /** Ordering comparison result for commutative term arguments, one of
     * 0 (unknown) 1 (less), 2 (equal), 3 (greater), 4 (incomparable)
     * @see Term::ArgumentOrder */
    unsigned order : 3;
    /** Number of distinct variables in the term, equal
     * to TERM_DIST_VAR_UNKNOWN if the number has not been
     * computed yet. */

    mutable unsigned distinctVars : TERM_DIST_VAR_BITS;
    /** term id hiding in this _info */
    // this should not be removed without care,
    // otherwise the bitfield layout might shift, resulting in broken pointer tagging
    unsigned id : 32;
    } _info;
#endif
  // but it was *not* portable because the layout of the bitfield is not guaranteed
  //
  // now we use a manual bitfield, as follows
  static constexpr unsigned
    TAG_BITS_START = 0,
    TAG_BITS_END = TAG_BITS_START + 2,
    POLARITY_BITS_START = TAG_BITS_END,
    POLARITY_BITS_END = POLARITY_BITS_START + 1,
    SHARED_BITS_START = POLARITY_BITS_END,
    SHARED_BITS_END = SHARED_BITS_START + 1,
    LITERAL_BITS_START = SHARED_BITS_END,
    LITERAL_BITS_END = LITERAL_BITS_START + 1,
    SORT_BITS_START = LITERAL_BITS_END,
    SORT_BITS_END = SORT_BITS_START + 1,
    HAS_TERM_VAR_BITS_START = SORT_BITS_END,
    HAS_TERM_VAR_BITS_END = HAS_TERM_VAR_BITS_START + 1,
    ORDER_BITS_START = HAS_TERM_VAR_BITS_END,
    ORDER_BITS_END = ORDER_BITS_START + 3,
    DISTINCT_VAR_BITS_START = ORDER_BITS_END,
    DISTINCT_VAR_BITS_END = DISTINCT_VAR_BITS_START + TERM_DIST_VAR_BITS,
    ID_BITS_START = DISTINCT_VAR_BITS_END,
    ID_BITS_END = ID_BITS_START + 32,
    TERM_BITS_START = 0,
    TERM_BITS_END = CHAR_BIT * sizeof(Term *);

  // various properties we want to check
  static_assert(TAG_BITS_START == 0, "tag must be the least significant bits");
  static_assert(TERM_BITS_START == 0, "term must be the least significant bits");
  static_assert(ID_BITS_END == 64, "whole thing must fit 64 bits exactly");
  static_assert(sizeof(void *) <= sizeof(uint64_t), "must be able to fit a pointer into a 64-bit integer");
  static_assert(AO_INCOMPARABLE < 8, "must be able to squash orderings into 3 bits");

  // getters and setters
  BITFIELD64_GET_AND_SET(unsigned, tag, Tag, TAG)
  BITFIELD64_GET_AND_SET(bool, polarity, Polarity, POLARITY)
  BITFIELD64_GET_AND_SET(bool, shared, Shared, SHARED)
  BITFIELD64_GET_AND_SET(bool, literal, Literal, LITERAL)
  BITFIELD64_GET_AND_SET(bool, sort, Sort, SORT)
  BITFIELD64_GET_AND_SET(bool, hasTermVar, HasTermVar, HAS_TERM_VAR)
  BITFIELD64_GET_AND_SET(unsigned, order, Order, ORDER)
  BITFIELD64_GET_AND_SET(uint32_t, distinctVars, DistinctVars, DISTINCT_VAR)
  BITFIELD64_GET_AND_SET(uint32_t, id, Id, ID)
  BITFIELD64_GET_AND_SET_PTR(Term*, term, Term, TERM)
  // end bitfield

  friend class Indexing::TermSharing;
  friend class Term;
  friend class Literal;
  friend class AtomicSort;
}; // class TermList
static_assert(sizeof(TermList) == 8, "size of TermList must be exactly 64 bits");

//special functor values
enum class SpecialFunctor {
  ITE,
  LET,
  FORMULA,
  TUPLE,
  LET_TUPLE,
  LAMBDA,
  MATCH, // <- keep this one the last, or modify SPECIAL_FUNCTOR_LAST accordingly
};
static constexpr SpecialFunctor SPECIAL_FUNCTOR_LAST = SpecialFunctor::MATCH;
std::ostream& operator<<(std::ostream& out, SpecialFunctor const& self);

/**
 * Class to represent terms and lists of terms.
 * @since 19/02/2008 Manchester, changed to use class TermList
 */
class alignas(8) Term
{
public:

  static constexpr unsigned SPECIAL_FUNCTOR_LOWER_BOUND  =  std::numeric_limits<unsigned>::max() - unsigned(SPECIAL_FUNCTOR_LAST);
  static SpecialFunctor toSpecialFunctor(unsigned f) {
    ASS_GE(f, SPECIAL_FUNCTOR_LOWER_BOUND);
    unsigned result = std::numeric_limits<unsigned>::max() - unsigned(f);
    ASS_LE(result, unsigned(SPECIAL_FUNCTOR_LAST))
    return SpecialFunctor(result);
  }
  static unsigned toNormalFunctor(SpecialFunctor f) 
  { return std::numeric_limits<unsigned>::max() - unsigned(f); }

  class SpecialTermData
  {
    friend class Term;
  private:
    union {
      struct {
        Formula * condition;
        TermList sort;
      } _iteData;
      struct {
        unsigned functor;
        VList* variables;
        TermList binding;
        TermList sort;
      } _letData;
      struct {
        Formula * formula;
      } _formulaData;
      struct {
        Term* term;
      } _tupleData;
      struct {
        unsigned functor;
        VList* symbols;
        TermList binding;
        TermList sort;
      } _letTupleData;
      struct {
        TermList lambdaExp;
        VList* _vars;
        SList* _sorts;  
        TermList sort; 
        TermList expSort;//TODO is this needed?
      } _lambdaData;
      struct {
        TermList sort;
        TermList matchedSort;
      } _matchData;
    };
    /** Return pointer to the term to which this object is attached */
    const Term* getTerm() const { return reinterpret_cast<const Term*>(this+1); }
  public:
    SpecialFunctor specialFunctor() const
    { return getTerm()->specialFunctor(); }

    Formula* getCondition() const { ASS_EQ(specialFunctor(), SpecialFunctor::ITE); return _iteData.condition; }
    unsigned getFunctor() const {
      ASS_REP(specialFunctor() == SpecialFunctor::LET || specialFunctor() == SpecialFunctor::LET_TUPLE, specialFunctor());
      return specialFunctor() == SpecialFunctor::LET ? _letData.functor : _letTupleData.functor;
    }
    VList* getLambdaVars() const { ASS_EQ(specialFunctor(), SpecialFunctor::LAMBDA); return _lambdaData._vars; }
    SList* getLambdaVarSorts() const { ASS_EQ(specialFunctor(), SpecialFunctor::LAMBDA); return _lambdaData._sorts; }
    TermList getLambdaExp() const { ASS_EQ(specialFunctor(), SpecialFunctor::LAMBDA); return _lambdaData.lambdaExp; }
    VList* getVariables() const { ASS_EQ(specialFunctor(), SpecialFunctor::LET); return _letData.variables; }
    VList* getTupleSymbols() const { return _letTupleData.symbols; }
    TermList getBinding() const {
      ASS_REP(specialFunctor() == SpecialFunctor::LET || specialFunctor() == SpecialFunctor::LET_TUPLE, specialFunctor());
      return TermList(specialFunctor() == SpecialFunctor::LET ? _letData.binding : _letTupleData.binding);
    }
    TermList getLambdaExpSort() const { ASS_EQ(specialFunctor(), SpecialFunctor::LAMBDA); return _lambdaData.expSort; }
    TermList getSort() const {
      switch (specialFunctor()) {
        case SpecialFunctor::ITE:
          return _iteData.sort;
        case SpecialFunctor::LET:
          return _letData.sort;
        case SpecialFunctor::LET_TUPLE:
          return _letTupleData.sort;
        case SpecialFunctor::LAMBDA:
          return _lambdaData.sort;
        case SpecialFunctor::MATCH:
          return _matchData.sort;
        default:
          ASSERTION_VIOLATION_REP(specialFunctor());
      }
    }
    Formula* getFormula() const { ASS_EQ(specialFunctor(), SpecialFunctor::FORMULA); return _formulaData.formula; }
    Term* getTupleTerm() const { return _tupleData.term; }
    TermList getMatchedSort() const { return _matchData.matchedSort; }
  };


  Term() throw();
  explicit Term(const Term& t) throw();
  static Term* create(unsigned function, unsigned arity, const TermList* args);
  static Term* create(unsigned fn, std::initializer_list<TermList> args);
  static Term* create(unsigned fn, Stack<TermList> const& args) { return Term::create(fn, args.length(), args.begin()); }
  template<class Iter>
  static Term* createFromIter(unsigned fn, Iter args) 
  { 
    Recycled<Stack<TermList>> stack;
    stack->loadFromIterator(args);
    return Term::create(fn, *stack); 
  }
  static Term* create(Term* t,TermList* args);
  static Term* createNonShared(unsigned function, unsigned arity, TermList* arg);
  static Term* createNonShared(Term* t,TermList* args);
  static Term* createNonShared(Term* t);
  static Term* cloneNonShared(Term* t);

  static Term* createConstant(const std::string& name);
  /** Create a new constant and insert in into the sharing structure */
  static Term* createConstant(unsigned symbolNumber) { return create(symbolNumber,0,0); }
  static Term* createITE(Formula * condition, TermList thenBranch, TermList elseBranch, TermList branchSort);
  static Term* createLet(unsigned functor, VList* variables, TermList binding, TermList body, TermList bodySort);
  static Term* createLambda(TermList lambdaExp, VList* vars, SList* sorts, TermList expSort);
  static Term* createTupleLet(unsigned functor, VList* symbols, TermList binding, TermList body, TermList bodySort);
  static Term* createFormula(Formula* formula);
  static Term* createTuple(unsigned arity, TermList* sorts, TermList* elements);
  static Term* createTuple(Term* tupleTerm);
  static Term* createMatch(TermList sort, TermList matchedSort, unsigned int arity, TermList* elements);
  static Term* create1(unsigned fn, TermList arg);
  static Term* create2(unsigned fn, TermList arg1, TermList arg2);

  //** fool constants
  static Term* foolTrue(); 
  static Term* foolFalse(); 

  VList* freeVariables() const;

  /** Return number of bytes before the start of the term that belong to it */
  size_t getPreDataSize() { return isSpecial() ? sizeof(SpecialTermData) : 0; }

  /** Function or predicate symbol of a term */
  const unsigned functor() const { return _functor; }


  SpecialFunctor specialFunctor() const 
  { return toSpecialFunctor(functor()); }
  std::string toString(bool topLevel = true) const;
  friend std::ostream& operator<<(std::ostream& out, Kernel::Term const& tl);
  static std::string variableToString(unsigned var);
  static std::string variableToString(TermList var);

  /** return the arguments 
   *
   *  WARNING: this function returns a pointer to the first argument
   *  which could be a sort when dealing with a polymorphic problem!
   * 
   *  Use with care! Consider whether the termArgs() function may be more
   *  suited to your needs before using this.
   */
  const TermList* args() const
  { return _args + _arity; }
  /** @see nthArguement(int) */ 
  const TermList* nthArgument(int n) const
  {
    ASS(n >= 0);
    ASS((unsigned)n < _arity);

    return _args + (_arity - n);
  }
  /** return the nth argument (counting from 0) 
   *
   *  Note that the arguments may be sort arguments as well as term arguments.
   *  i.e. nthArgument(n) will return 
   *    - a sort, for 0 <= n < numTypeArguemnts()
   *    - a term, for numTypeArguments() <= n < arity()
   *
   *  If you want to access a specific term or type argument use typeArg(int) or termArg(int) instead.
   */ 
  TermList* nthArgument(int n)
  {
    ASS(n >= 0);
    ASS((unsigned)n < _arity);

    return _args + (_arity - n);
  }

  /** returns the nth term argument. for 0 <= n <= numTermArguments  */
  TermList termArg(unsigned n) const;

  /** returns the nth type argument. for 0 <= n <= numTypeArguments  */
  TermList typeArg(unsigned n) const;

  /**
   * Return the number of type arguments for a polymorphic term (or 0 if monomorphic).
   */
  unsigned numTypeArguments() const;

  /**
   * Return the number of term arguments for a term (equal to _arity if monomorphic).
   */  
  unsigned numTermArguments() const;

  /** Return the 1st term argument for a polymorphic term.
    * Call hasTermArgs before calling this or test the result for
    * non-emptiness
    * In the monomorphic case, the same as args()
    */
  TermList* termArgs();

  /** Return the 1st type argument for a polymorphic term.
    * returns a nullpointer if the term not polymorphic
    * This is technically almost the same thing as calling args(), 
    * but can be used to increase readability of code.
    */
  const TermList* typeArgs() const;

  /** Indexing operator for accessing arguments */
  const TermList operator[](int i) const {
    return *nthArgument(i);
  }
  TermList operator[](int i) {
    return *nthArgument(i);
  }

  /** return the arguments 
   *
   *  WARNING: this function returns a pointer to the first argument
   *  which could be a sort when dealing with a polymorphic problem!
   * 
   *  Use with care! Consider whether the termArgs() function may be more
   *  suited to your needs before using this.
   */  
  TermList* args()
  { return _args + _arity; }


  template<class GetArg>
  static unsigned termHash(unsigned functor, GetArg getArg, unsigned arity) {
    return DefaultHash::hashIter(
        range(0, arity).map([&](auto i) {
          TermList t = getArg(i);
          return DefaultHash::hashBytes(
              reinterpret_cast<const unsigned char*>(&t),
              sizeof(TermList)
              );
          }),
        DefaultHash::hash(functor));
  }

  /**
   * Return the hash function of the top-level of a complex term.
   * @pre The term must be non-variable
   * @since 28/12/2007 Manchester
   */
  unsigned hash() const 
  { return termHash(_functor, [&](auto i) { return *nthArgument(i); }, _arity); }

  /** return the arity */
  unsigned arity() const
  { return _arity; }
  static void* operator new(size_t,unsigned arity,size_t preData=0);
  /** make the term into a symbol term */
  void makeSymbol(unsigned number,unsigned arity)
  {
    _functor = number;
    _arity = arity;
  }
  void destroy();
  void destroyNonShared();
  Term* apply(Substitution& subst);

  /** True iff all immediate arguments are variables */
  bool allArgumentsAreVariables() const
  {
    for(unsigned i = 0; i < arity(); i++)
      if(!nthArgument(i)->isVar())
        return false;

    return true;
  }

  /** True if the term is ground. Only applicable to shared terms */
  bool ground() const
  {
    ASS(_args[0]._shared());
    return numVarOccs() == 0;
  } // ground

  /** True if the term contains a term variable (type variables don't count)
   *  Only applicable to shared terms */
  bool hasTermVar() const
  {
    ASS(shared());
    return _args[0]._hasTermVar();
  } // ground

  /** True if the term is shared */
  bool shared() const
  { return _args[0]._shared(); } // shared

  /** Return the weight. Applicable only to shared terms */
  unsigned weight() const
  {
    ASS(shared());
    return _weight;
  }

  int maxRedLength() const
  {
    ASS(shared());
    return _maxRedLen;    
  }

  int kboWeight(const void* kboInstance) const
  {
    ASS(_kboInstance || _kboWeight == -1);
    ASS(!_kboInstance || _kboInstance == kboInstance);
    return _kboWeight;
  }

  void setKboWeight(int w, const void* kboInstance)
  {
#if VDEBUG
    ASS(!_kboInstance);
    _kboInstance = kboInstance;
#endif
    _kboWeight = w;
  }

  /** Mark term as shared */
  void markShared()
  {
    ASS(! shared());
    _args[0]._setShared(true);
  } // markShared

  /** Set term weight */
  void setWeight(unsigned w)
  {
    _weight = w;
  } // setWeight

  /** Set term id */
  void setId(unsigned id);

  /** Set (shared) term's id */
  unsigned getId() const
  {
    ASS(shared());
    return _args[0]._id();
  }
  
  void setMaxRedLen(int rl)
  {
    _maxRedLen = rl;
  } // setWeight

  /** Set the number of variable _occurrences_ */
  void setNumVarOccs(unsigned v)
  {
    if(_isTwoVarEquality) {
      ASS_EQ(v,2);
      return;
    }
    _vars = v;
  } // setVars

  void setHasTermVar(bool b)
  {
    ASS(shared() && !isSort());
    _args[0]._setHasTermVar(b);
  }

  /** Return the number of variable _occurrences_ */
  unsigned numVarOccs() const
  {
    ASS(shared());
    if(_isTwoVarEquality) {
      return _sort.isVar() ? 3 : 2 + _sort.term()->numVarOccs();
    }
    return _vars;
  } // vars()

  /**
   * Return true iff the object is an equality between two variables.
   *
   * This value is set during insertion into the term sharing structure or
   * for terms with special subterms during construction.
   * (I.e. can be used basically anywhere).
   */
  bool isTwoVarEquality() const
  {
    return _isTwoVarEquality;
  }

  const std::string& functionName() const;

  /** True if the term is, in fact, a literal */
  bool isLiteral() const { return _args[0]._literal(); }
  /** True if the term is, in fact, a sort */
  bool isSort() const { return _args[0]._sort(); }
  TermKind kind() const { return isSort() ? TermKind::SORT 
                               : isLiteral() ? TermKind::LITERAL
                               : TermKind::TERM; }
  /** true if the term is an application */
  bool isApplication() const;

  /** Return an index of the argument to which @b arg points */
  unsigned getArgumentIndex(TermList* arg)
  {
    unsigned res=arity()-(arg-_args);
    ASS_L(res,arity());
    return res;
  }

#if VDEBUG
  std::string headerToString() const;
  void assertValid() const;
#endif


  static TermIterator getVariableIterator(TermList tl);

  // the number of _distinct_ variables within the term
  unsigned getDistinctVars()
  {
    if(_args[0]._distinctVars()==TERM_DIST_VAR_UNKNOWN) {
      unsigned res=computeDistinctVars();
      if(res<TERM_DIST_VAR_UNKNOWN) {
        _args[0]._setDistinctVars(res);
      }
      return res;
    } else {
      ASS_L(_args[0]._distinctVars(),0x100000);
      return _args[0]._distinctVars();
    }
  }

  bool couldBeInstanceOf(Term* t)
  {
    ASS(shared());
    ASS(t->shared());
    if(t->functor()!=functor()) {
      return false;
    }
    return true;
  }

  bool containsSubterm(TermList v) const;
  bool containsAllVariablesOf(Term* t);
  size_t countSubtermOccurrences(TermList subterm);

  /** Return true if term has no non-constant functions as subterms */
  bool isShallow() const;

  /** set the colour of the term */
  void setColor(Color color)
  {
    ASS(_color == static_cast<unsigned>(COLOR_TRANSPARENT) || _color == static_cast<unsigned>(color));
    _color = color;
  } // setColor
  /** return the colour of the term */
  Color color() const { return static_cast<Color>(_color); }

  bool skip() const;

  /** Return true if term is an interpreted constant or contains one as its subterm */
  bool hasInterpretedConstants() const { return _hasInterpretedConstants; }
  /** Assign value that will be returned by the hasInterpretedConstants() function */
  void setInterpretedConstantsPresence(bool value) { _hasInterpretedConstants=value; }

  /** Return true if term is either an if-then-else or a let...in expression */
  bool isSpecial() const { return functor() >= SPECIAL_FUNCTOR_LOWER_BOUND; }

  bool isITE()      const { return functor() == toNormalFunctor(SpecialFunctor::ITE); }
  bool isLet()      const { return functor() == toNormalFunctor(SpecialFunctor::LET); }
  bool isTupleLet() const { return functor() == toNormalFunctor(SpecialFunctor::LET_TUPLE); }
  bool isTuple()    const { return functor() == toNormalFunctor(SpecialFunctor::TUPLE); }
  bool isFormula()  const { return functor() == toNormalFunctor(SpecialFunctor::FORMULA); }
  bool isLambda()   const { return functor() == toNormalFunctor(SpecialFunctor::LAMBDA); }
  bool isMatch()    const { return functor() == toNormalFunctor(SpecialFunctor::MATCH); }
  bool isBoolean() const;
  bool isSuper() const; 
  
  /** Return pointer to structure containing extra data for special terms such as
   * if-then-else or let...in */
  const SpecialTermData* getSpecialData() const { return const_cast<Term*>(this)->getSpecialData(); }
  /** Return pointer to structure containing extra data for special terms such as
   * if-then-else or let...in */
  SpecialTermData* getSpecialData() {
    ASS(isSpecial());
    return reinterpret_cast<SpecialTermData*>(this)-1;
  }

  virtual bool computable() const;
  virtual bool computableOrVar() const;

protected:
  std::string headToString() const;

  unsigned computeDistinctVars() const;

  /**
   * Return argument order value stored in term.
   *
   * The default value (which is returned if no value was set using the
   * @c setArgumentOrder() function) is zero.
   *
   * Currently, this function is used only by @c Ordering::getEqualityArgumentOrder().
   */
  ArgumentOrderVals getArgumentOrderValue() const
  {
    return static_cast<ArgumentOrderVals>(_args[0]._order());
  }

  /**
   * Store argument order value in term.
   *
   * The value must be non-negative and less than 8.
   *
   * Currently, this function is used only by @c Ordering::getEqualityArgumentOrder().
   */
  void setArgumentOrderValue(ArgumentOrderVals val)
  {
    ASS_GE(val,AO_UNKNOWN);
    ASS_LE(val,AO_INCOMPARABLE);

    _args[0]._setOrder(val);
  }

  /** The number of this symbol in a signature */
  unsigned _functor;
  /** Arity of the symbol */
  unsigned _arity : 28;
  /** colour, used in interpolation and symbol elimination */
  unsigned _color : 2;
  /** Equal to 1 if the term/literal contains any interpreted constants */
  unsigned _hasInterpretedConstants : 1;
  /** If true, the object is an equality literal between two variables */
  unsigned _isTwoVarEquality : 1;
  /** Weight of the symbol, i.e. sum of symbol and variable occurrences. */
  unsigned _weight;
  /** Cached weight of the term for KBO, otherwise -1 and invalid. Note that
   * KBO symbol weights are not necessarily 1, so this can differ from @b _weight. */
  int _kboWeight;
#if VDEBUG
  /** KBO instance that uses the cached value @b _kboWeight. */
  const void* _kboInstance;
#endif
  /** length of maximum reduction length */
  int _maxRedLen;
  union {
    /** If _isTwoVarEquality is false, this value is valid and contains
     * number of occurrences of variables */
    unsigned _vars;
    /** If _isTwoVarEquality is true, this value is valid and contains
     * the sort of the top-level variables */
    TermList _sort;
  };

  /** The list of arguments of size type arity + term arity + 1. The first
   *  argument stores the term weight and the mask (the last two bits are 0).
   */
  TermList _args[1];

  friend class TermList;
  friend class Indexing::TermSharing;
  friend class Ordering;

public:
  /**
   * Iterator returning arguments of a term left-to-right.
   */
  class Iterator
  {
  public:
    DECL_ELEMENT_TYPE(TermList);
    Iterator(Term* t) : _next(t->args()) {}
    bool hasNext() const { return _next->isNonEmpty(); }
    TermList next()
    {
      ASS(hasNext());
      TermList res = *_next;
      _next = _next->next();
      return res;
    }
  private:
    TermList* _next;
  }; // Term::Iterator
}; // class Term


/**
 * Class of AtomicSort.
 */
class AtomicSort
  : public Term
{
public:
  AtomicSort();
  explicit AtomicSort(const AtomicSort& t) throw();

  AtomicSort(unsigned functor,unsigned arity) throw()
  {
    _functor = functor;
    _arity = arity;
    _args[0]._setLiteral(false);
    _args[0]._setSort(true);
  }

  static AtomicSort* create(unsigned typeCon, unsigned arity, const TermList* args);
  static AtomicSort* create2(unsigned tc, TermList arg1, TermList arg2);
  static AtomicSort* create(AtomicSort const* t,TermList* args);
  static AtomicSort* createConstant(unsigned typeCon) { return create(typeCon,0,0); }
  static AtomicSort* createConstant(const std::string& name); 

  /** True if the sort is a higher-order arrow sort */
  bool isArrowSort() const;
  /** True if the sort $o */
  bool isBoolSort() const;
  /** true if sort is the sort of an array */
  bool isArraySort() const;
  /** true if sort is the sort of an tuple */
  bool isTupleSort() const;

  const std::string& typeConName() const;  
  
  static TermList arrowSort(TermStack& domSorts, TermList range);
  static TermList arrowSort(TermList s1, TermList s2);
  static TermList arrowSort(TermList s1, TermList s2, TermList s3);
  static TermList arraySort(TermList indexSort, TermList innerSort);
  static TermList tupleSort(unsigned arity, TermList* sorts);
  static TermList defaultSort();
  static TermList superSort();
  static TermList boolSort();
  static TermList intSort();
  static TermList realSort();
  static TermList rationalSort();

private:

  static AtomicSort* createNonShared(unsigned typeCon, unsigned arity, TermList* arg);
  static AtomicSort* createNonSharedConstant(unsigned typeCon) { return createNonShared(typeCon,0,0); }
};

/**
 * Class of literals.
 * @since 06/05/2007 Manchester
 */
class Literal
  : public Term
{
public:
  /** True if equality literal */
  bool isEquality() const
  { return functor() == 0; }

  Literal();
  explicit Literal(const Literal& l) throw();

  /**
   * Create a literal.
   * @since 16/05/2007 Manchester
   */
  Literal(unsigned functor,unsigned arity,bool polarity) throw()
  {
    _functor = functor;
    _arity = arity;
    _args[0]._setPolarity(polarity);
    _args[0]._setSort(false);
    _args[0]._setLiteral(true);
  }

  /**
   * A unique header, 2*p is negative and 2*p+1 if positive where p is
   * the number of the predicate symbol.
   */
  unsigned header() const
  { return 2*_functor + polarity(); }
  /**
   * Header of the complementary literal, 2*p+1 is negative and 2*p
   * if positive where p is the number of the predicate symbol.
   */
  unsigned complementaryHeader() const
  { return 2*_functor + 1 - polarity(); }

  static bool headersMatch(Literal* l1, Literal* l2, bool complementary);
  /** set polarity to true or false */
  void setPolarity(bool positive)
  { _args[0]._setPolarity(positive); }

  TermList eqArgSort() const;
  
  // prevent bugs through implicit bool <-> unsigned conversions
  template<class Iter> static Literal* createFromIter(unsigned predicate, unsigned polarity, Iter iter) = delete;
  template<class Iter> static Literal* createFromIter(    bool predicate, unsigned polarity, Iter iter) = delete;
  template<class Iter> static Literal* createFromIter(    bool predicate,     bool polarity, Iter iter) = delete;

  template<class Iter>
  static Literal* createFromIter(unsigned predicate, bool polarity, Iter iter) {
    RStack<TermList> args;
    while (iter.hasNext()) {
      args->push(iter.next());
    }
    return Literal::create(predicate, args->size(), polarity, args->begin());
  }

  template<class Iter>
  static Literal* createFromIter(Literal* lit, Iter iter) {
    if (lit->isEquality()) {
      return  Literal::createEquality(lit->polarity(), iter.tryNext().unwrap(), iter.tryNext().unwrap(), lit->eqArgSort());
    } else {
      return Literal::createFromIter(lit->functor(), bool(lit->polarity()), std::move(iter));
    }
  }

  static Literal* create(unsigned predicate, unsigned arity, bool polarity, TermList* args);
  static Literal* create(unsigned predicate, bool polarity, std::initializer_list<TermList>);
  static Literal* create(Literal* l,bool polarity);
  static Literal* create(Literal* l,TermList* args);
  static Literal* createEquality(bool polarity, TermList arg1, TermList arg2, TermList sort);
  static Literal* create1(unsigned predicate, bool polarity, TermList arg);
  static Literal* create2(unsigned predicate, bool polarity, TermList arg1, TermList arg2);

  /**
   * Return the hash function of the top-level of a literal.
   * @since 30/03/2008 Flight Murcia-Manchester
   */
  template<bool flip = false>
  unsigned hash() const
  {
    return Literal::literalHash(functor(), polarity() ^ flip,
        [&](auto i) -> TermList const& { return *nthArgument(i); }, arity(),
        someIf(isTwoVarEquality(), [&](){ return twoVarEqSort(); }));
  }

  template<class GetArg>
  static unsigned literalEquals(const Literal* lit, unsigned functor, bool polarity, GetArg getArg, unsigned arity, Option<TermList> twoVarEqSort) {
    if (functor != lit->functor() || polarity != lit->polarity()) return false;

    if (functor == 0) { // i.e., isEquality
      ASS_EQ(arity, 2)
      ASS(rightArgOrder(getArg(0), getArg(1)))
      ASS(rightArgOrder(*lit->nthArgument(0), *lit->nthArgument(1)))

      if (someIf(lit->isTwoVarEquality(), [&](){ return lit->twoVarEqSort(); }) != twoVarEqSort) {
        return false;
      }
      return std::make_tuple(*lit->nthArgument(0), *lit->nthArgument(1)) == std::make_tuple(getArg(0), getArg(1));

    } else {
      ASS(twoVarEqSort.isNone())
      return range(0, arity).all([&](auto i) { return *lit->nthArgument(i) == getArg(i); });
    }
  }

  static bool rightArgOrder(TermList const& lhs, TermList const& rhs);

  template<class GetArg>
  static unsigned literalHash(unsigned functor, bool polarity, GetArg getArg, unsigned arity, Option<TermList> twoVarEqSort) {
    if (functor == 0) { // i.e., isEquality
      ASS_EQ(arity, 2)
      ASS(rightArgOrder(getArg(0), getArg(1)))
      return HashUtils::combine(
          DefaultHash::hash(polarity),
          DefaultHash::hash(functor),
          DefaultHash::hash(twoVarEqSort),
          getArg(0).defaultHash(),
          getArg(1).defaultHash());
    } else {
      ASS(twoVarEqSort.isNone())
      return HashUtils::combine(
          DefaultHash::hash(polarity),
          Term::termHash(functor, getArg, arity));
    }
  }



  static Literal* complementaryLiteral(Literal* l);
  /** If l is positive, return l; otherwise return its complementary literal. */
  static Literal* positiveLiteral(Literal* l) {
    return l->isPositive() ? l : complementaryLiteral(l);
  }

  // destructively swap arguments of an equation
  // the term is assumed to be non-shared
  void argSwap() {
    ASS(isEquality() && !shared());
    ASS(arity() == 2);

    TermList* ts1 = args();
    TermList* ts2 = ts1->next();
    using std::swap;//ADL
    swap(ts1->_content, ts2->_content);
  }

  /** true if positive */
  bool isPositive() const
  {
    return polarity();
  } // isPositive

  /** true if negative */
  bool isNegative() const
  {
    return !polarity();
  } // isNegative

  /** return polarity, 1 if positive and 0 if negative */
  int polarity() const
  {
    return _args[0]._polarity();
  } // polarity

  /**
   * Mark this object as an equality between two variables.
   */
  void markTwoVarEquality()
  {
    ASS(!shared());
    ASS(isEquality());
    ASS(nthArgument(0)->isVar() || !nthArgument(0)->term()->shared());
    ASS(nthArgument(1)->isVar() || !nthArgument(1)->term()->shared());

    _isTwoVarEquality = true;
  }


  /** Return sort of the variables in an equality between two variables.
   * This value is set during insertion into the term sharing structure
   */
  TermList twoVarEqSort() const
  {
    ASS(isTwoVarEquality());

    return _sort;
  }

  /** Assign sort of the variables in an equality between two variables. */
  void setTwoVarEqSort(TermList sort)
  {
    ASS(isTwoVarEquality());

    _sort = sort;
  }

//   /** Applied @b subst to the literal and return the result */
  Literal* apply(Substitution& subst);


  inline bool couldBeInstanceOf(Literal* lit, bool complementary)
  {
    ASS(shared());
    ASS(lit->shared());
    return headersMatch(this, lit, complementary);
  }

  bool isAnswerLiteral() const;

  friend std::ostream& operator<<(std::ostream& out, Kernel::Literal const& tl);
  std::string toString() const;

  const std::string& predicateName() const;

  virtual bool computable() const;
  virtual bool computableOrVar() const;

private:
  template<class GetArg>
  static Literal* create(unsigned predicate, unsigned arity, bool polarity, GetArg args, Option<TermList> twoVarEqSort = Option<TermList>());
}; // class Literal

// TODO used in some proofExtra output
//      find a better place for this?
bool positionIn(TermList& subterm,TermList* term, std::string& position);
bool positionIn(TermList& subterm,Term* term, std::string& position);

} // namespace Kernel

template<>
struct std::hash<Kernel::TermList> {
  size_t operator()(Kernel::TermList const& t) const 
  { return t.defaultHash(); }
};

#endif
