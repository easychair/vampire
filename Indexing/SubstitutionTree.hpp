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
 * @file SubstitutionTree.hpp
 * Defines class SubstitutionTree.
 *
 * @since 16/08/2008 flight Sydney-San Francisco
 */

#ifndef __SubstitutionTree__
#define __SubstitutionTree__

#include <utility>

#include "Forwards.hpp"

#include "Kernel/UnificationWithAbstraction.hpp"
#include "Lib/Exception.hpp"
#include "Lib/Reflection.hpp"
#include "Lib/VirtualIterator.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/Comparison.hpp"
#include "Lib/Int.hpp"
#include "Lib/Stack.hpp"
#include "Lib/List.hpp"
#include "Lib/SkipList.hpp"
#include "Lib/BinaryHeap.hpp"
#include "Lib/Backtrackable.hpp"
#include "Lib/ArrayMap.hpp"
#include "Lib/Array.hpp"
#include "Lib/BiMap.hpp"
#include "Kernel/ApplicativeHelper.hpp"
#include "Lib/Recycled.hpp"

#include "Kernel/RobSubstitution.hpp"
#include "Kernel/Renaming.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/OperatorType.hpp"
#include "Debug/Tracer.hpp"
#include "Lib/Option.hpp"
#include "Kernel/Signature.hpp"

#include "Lib/Allocator.hpp"

#include "Index.hpp"

#if VDEBUG
#include <iostream>
#endif


static constexpr int QUERY_BANK=0;
static constexpr int RESULT_BANK=1;
static constexpr int NORM_RESULT_BANK=3;
using namespace Lib;
using namespace Kernel;

#define UARR_INTERMEDIATE_NODE_MAX_SIZE 4

#define REORDERING 1

namespace Indexing {


  /** a counter that is compiled away in release mode */
  struct Cntr {
#if VDEBUG
    Cntr() : self(0) {}
    int self;
    operator int() const { return self; }
#endif 
  };

template<class LeafData_>
class SubstitutionTree;
    
template<class LD> std::ostream& operator<<(std::ostream& out, SubstitutionTree<LD> const& self);
template<class LD> std::ostream& operator<<(std::ostream& out, Output::Multiline<SubstitutionTree<LD>> const& self);


  /** a reference to a Cntr that increments the counter when it is created and decrements it when it goes out of scope
   * This can be used to count the number of instances when an object of this type is added as a member field to the class 
   * that should be counted */
  class InstanceCntr {
  public:
#if VDEBUG
    Cntr* _cntr = nullptr;

    InstanceCntr(InstanceCntr&& other)
      : InstanceCntr()
    { std::swap(other._cntr, _cntr); }

    InstanceCntr(InstanceCntr const& other)
      : InstanceCntr()
    { 
      _cntr = other._cntr;
      if (_cntr) _cntr->self++; 
    }

    InstanceCntr& operator=(InstanceCntr&& other) 
    { std::swap(other._cntr, _cntr); return *this; }

    InstanceCntr& operator=(InstanceCntr const& other)
    { 
      if (_cntr) _cntr->self--; 
      _cntr = other._cntr; 
      if (_cntr) _cntr->self++;
      return *this; 
    }

    InstanceCntr() : _cntr(nullptr) 
    { }

    InstanceCntr(Cntr& cntr) : _cntr(&cntr) 
    { _cntr->self++; }

    ~InstanceCntr() 
    { if (_cntr) _cntr->self--; }

    void reset() { 
      if (_cntr) _cntr->self--; 
      _cntr = nullptr; 
    }
#else // VDEBUG
   InstanceCntr(){ }
   InstanceCntr(Cntr& parent) {}
   void reset() { }
#endif 
};

/**
 * Class of substitution trees. 
 *
 * We can either store typed terms, or literals in a subtitution tree.
 * Classically we'd think of inserting/removing only one term t into a substitution tree. 
 * This can be understood as inserting the substitution { S0 -> t } into the tree.
 *
 * In general we can insertt a substitution with more than just this one binding. 
 * This is what we do in order to store the sort of variables, and in order to insert all the arguments of a literal:
 * - For a term t of sort s we insert { S0 -> t; S1 -> s }
 * - For literals (~)P(t0..tn) we insert { S0 -> t0 .. Sn -> tn }.
 * (Note that we do not check the predicate or the polarity of literals here. This happens in LiteralSubstitutionTree)
 */
template<class LeafData_>
class SubstitutionTree final
{

public:
  using LeafData = LeafData_;

  static constexpr int QRS_QUERY_BANK = 0;
  static constexpr int QRS_RESULT_BANK = 1;

  SubstitutionTree(SubstitutionTree const&) = delete;
  SubstitutionTree& operator=(SubstitutionTree const& other) = delete;
  static void swap(SubstitutionTree& self, SubstitutionTree& other) {
    std::swap(self._nextVar, other._nextVar);
    std::swap(self._root,    other._root);
  }
  SubstitutionTree& operator=(SubstitutionTree && other) { swap(*this,other); return *this; }
  SubstitutionTree(SubstitutionTree&& other) : SubstitutionTree() { swap(*this, other); }

  SubstitutionTree() : _nextVar(0), _root(nullptr) {}

  ~SubstitutionTree()
  {
    ASS_EQ(_iterCnt,0);
    delete _root;
  }

#define VERBOSE_OUTPUT_OPERATORS 0
  friend std::ostream& operator<<(std::ostream& out, SubstitutionTree const& self)
  {
#if VERBOSE_OUTPUT_OPERATORS
    out << "{ nextVar: S" << self._nextVar << ", root: (";
#endif // VERBOSE_OUTPUT_OPERATORS
    if (self._root) {
      out << *self._root;
    } else {
      out << "<empty tree>";
    }
#if VERBOSE_OUTPUT_OPERATORS
    out << ") }";
#endif // VERBOSE_OUTPUT_OPERATORS
    return out;
  }
#undef VERBOSE_OUTPUT_OPERATORS
  friend std::ostream& operator<<(std::ostream& out, Output::Multiline<SubstitutionTree> const& self)
  {
    if (self.self._root) {
      self.self._root->output(out, true, /* indent */ 0);
    } else {
      out << "<empty tree>";
    }
    return out;
  }


  template<class T>
  friend std::ostream& operator<<(std::ostream& out, SubstitutionTree<T> const& self);
  template<class T>
  friend std::ostream& operator<<(std::ostream& out, Output::Multiline<SubstitutionTree<T>> const& self);
  typedef VirtualIterator<LeafData*> LDIterator;

  template<class I> using QueryResultIter = VirtualIterator<QueryRes<LeafData, typename I::Unifier>>;
  template<class I, class TermOrLit, class... Args>
  auto iterator(TermOrLit query, bool retrieveSubstitutions, bool reversed, Args... args)
  { return isEmpty() ? VirtualIterator<ELEMENT_TYPE(I)>::getEmpty()
                     : pvi(iterPointer(Recycled<I>(this, _root, query, retrieveSubstitutions, reversed, std::move(args)...)));
  }

  class LDComparator
  {
  public:
    template<class LD>
    static Comparison compare(const LD& ld1, const LD& ld2)
    {
      return ld1 < ld2 ? Comparison::LESS 
           : ld1 > ld2 ? Comparison::GREATER
           : Comparison::EQUAL;
    }
  };

  static bool isGround(Literal* literal) { return literal->ground(); }
  static bool isGround(TermList term) { return term.ground(); }

  enum NodeAlgorithm
  {
    UNSORTED_LIST=1,
    SKIP_LIST=2,
    SET=3
  };

  class Node {

    /** term at this node */
    TermList _term;
#define CACHE_FUNCTOR 0
#if CACHE_FUNCTOR
    /** if _term is a Term* we cache the functor so we are faster than calling _term->term(); */
    unsigned _functor;
#endif // CACHE_FUNCTOR

  public:
    friend std::ostream& operator<<(std::ostream& out, Output::Multiline<Node> const& self) 
    { self.self.output(out, /* multiline = */ true, self.indent); return out; }
    friend std::ostream& operator<<(std::ostream& out, Node const& self) 
    { self.output(out, /* multiline = */ false, /* indent */ 0); return out; }
    inline Node() : _term(TermList::empty()) {}
    inline Node(TermList ts) : Node() { setTerm(ts); }
    virtual ~Node();
    /** True if a leaf node */
    virtual bool isLeaf() const = 0;
    virtual bool isEmpty() const = 0;
    /**
     * Return number of elements held in the node.
     *
     * Descendant classes should override this method.
     */
    virtual int size() const { NOT_IMPLEMENTED; }
    virtual NodeAlgorithm algorithm() const = 0;

    /**
     * Remove all referenced structures without destroying them.
     *
     * This is used when the implementation of a node is being changed.
     * The current node will be deleted, but we don't want to destroy
     * structures, that are taken over by the new node implementation.
     */
    virtual void makeEmpty() { _term = TermList::empty(); }
    static void split(Node** pnode, TermList* where, int var);

    void setTerm(TermList t) { 
      _term = t; 
#if CACHE_FUNCTOR
      if (t.isTerm()) { _functor =  t.term()->functor(); } 
#endif
    }
    inline TermList::Top top() { 
#if CACHE_FUNCTOR
      auto out = _term.isTerm() ? TermList::Top::functor(_functor) : _term.top(); 
      ASS_EQ(_term.top(), out); 
      return out; 
#else // !CACHE_FUNCTOR
      return _term.top(); 
#endif
    }
    // if there is a write to this reference, a term with the same top() as the current one has to be written. otherwise use setTerm instead
    TermList& term() { return _term; }
    TermList const& term() const { return _term; }
    virtual void output(std::ostream& out, bool multiline, int indent) const = 0;
  };


    typedef VirtualIterator<Node**> NodeIterator;

    class IntermediateNode
          : public Node
    {
    public:
      /** Build a new intermediate node which will serve as the root*/
      inline
      IntermediateNode(unsigned childVar) : childVar(childVar) {}

      /** Build a new intermediate node */
      inline
      IntermediateNode(TermList ts, unsigned childVar) : Node(ts), childVar(childVar) {}

      inline
      bool isLeaf() const final override { return false; };

      virtual NodeIterator allChildren() = 0;
      virtual NodeIterator variableChildren() = 0;
      /**
       * Return pointer to pointer to child node with top symbol
       * of @b t. This pointer to node can be changed.
       *
       * If canCreate is true and such child node does
       * not exist, pointer to null pointer is returned, and it's
       * assumed, that pointer to newly created node with given
       * top symbol will be put there.
       *
       * If canCreate is false, null pointer is returned in case
       * suitable child does not exist.
       */
      virtual Node** childByTop(TermList::Top t, bool canCreate) = 0;


      /**
       * Remove child which points to node with top symbol of @b t.
       * This node has to still exist in time of the call to remove method.
       */
      virtual void remove(TermList::Top t) = 0;
      /**
       * Remove all children of the node without destroying them.
       */
      virtual void removeAllChildren() = 0;

      void destroyChildren();

      void makeEmpty() final override
      {
        Node::makeEmpty();
        removeAllChildren();
      }

      void loadChildren(NodeIterator children);

      const unsigned childVar;

      virtual void output(std::ostream& out, bool multiline, int indent) const override;
    }; // class SubstitutionTree::IntermediateNode

    class Leaf
    : public Node
    {
    public:
      /** Build a new leaf which will serve as the root */
      inline
      Leaf()
      {}
      /** Build a new leaf */
      inline
      Leaf(TermList ts) : Node(ts) {}

      inline
      bool isLeaf() const final override { return true; };
      virtual LDIterator allChildren() = 0;
      virtual void insert(LeafData ld) = 0;
      virtual void remove(LeafData ld) = 0;
      void loadChildren(LDIterator children);
      virtual void output(std::ostream& out, bool multiline, int indent) const override;
    };

    //These classes and methods are defined in SubstitutionTree_Nodes.cpp
    class UListLeaf;
    class SListIntermediateNode;
    class SListLeaf;
    class SetLeaf;
    static Leaf* createLeaf();
    static Leaf* createLeaf(TermList ts);
    static void ensureLeafEfficiency(Leaf** l);
    static IntermediateNode* createIntermediateNode(unsigned childVar);
    static IntermediateNode* createIntermediateNode(TermList ts, unsigned childVar);
    static void ensureIntermediateNodeEfficiency(IntermediateNode** inode);

    class UArrIntermediateNode
    : public IntermediateNode
    {
    public:
      inline
      UArrIntermediateNode(unsigned childVar) : IntermediateNode(childVar), _size(0)
      {
        _nodes[0]=0;
      }
      inline
      UArrIntermediateNode(TermList ts, unsigned childVar) : IntermediateNode(ts, childVar), _size(0)
      {
        _nodes[0]=0;
      }

      ~UArrIntermediateNode()
      {
        if(!isEmpty()) {
          IntermediateNode::destroyChildren();
        }
      }

      void removeAllChildren()
      {
        _size=0;
        _nodes[0]=0;
      }

      NodeAlgorithm algorithm() const { return UNSORTED_LIST; }
      bool isEmpty() const { return !_size; }
      int size() const { return _size; }
      NodeIterator allChildren()
      { return pvi( arrayIter(_nodes,_size).map([](Node *& n) { return &n; }) ); }

      NodeIterator variableChildren()
      {
        return pvi( arrayIter(_nodes, _size)
              .filter([](auto& n) { return n->term().isVar(); })
              .map([](Node *& n) { return &n; }));
      }
      virtual Node** childByTop(TermList::Top t, bool canCreate);
      void remove(TermList::Top t);

      USE_ALLOCATOR(UArrIntermediateNode);

      int _size;
      Node* _nodes[UARR_INTERMEDIATE_NODE_MAX_SIZE+1];
    };

    class SListIntermediateNode
    : public IntermediateNode
    {
    public:
      SListIntermediateNode(unsigned childVar) : IntermediateNode(childVar) {}
      SListIntermediateNode(TermList ts, unsigned childVar) : IntermediateNode(ts, childVar) {}

      ~SListIntermediateNode()
      {
        if(!isEmpty()) {
          IntermediateNode::destroyChildren();
        }
      }

      void removeAllChildren()
      {
        while(!_nodes.isEmpty()) {
          _nodes.pop();
        }
      }

      static IntermediateNode* assimilate(IntermediateNode* orig);

      inline
      NodeAlgorithm algorithm() const { return SKIP_LIST; }
      inline
      bool isEmpty() const { return _nodes.isEmpty(); }
      int size() const { return _nodes.size(); }
      inline
      NodeIterator allChildren()
      {
        return pvi(_nodes.ptrIter());
      }
      inline
      NodeIterator variableChildren()
      {
        return pvi( getWhileLimitedIterator(
                      _nodes.ptrIter(),
                      [](auto* n) {return (*n)->term().isVar(); }));
      }
      virtual Node** childByTop(TermList::Top t, bool canCreate)
      {
        Node** res;
        bool found=_nodes.getPosition(t,res,canCreate);
        if(!found) {
          if(canCreate) {
            *res=0;
          } else {
            res=0;
          }
        }
        return res;
      }

      inline void remove(TermList::Top t)
      { _nodes.remove(t); }

      USE_ALLOCATOR(SListIntermediateNode);

      class NodePtrComparator
      {
      public:
        inline static Comparison compare(TermList::Top& l, Node* r)
        { 
          if(l.var()) {
            return r->term().isVar() ? Int::compare(*l.var(), r->term().var())
                                     : LESS;
          } else {
            return r->term().isVar() ? GREATER
                                     : Int::compare(l.functor()->functor, r->term().term()->functor());
          }
        }
      };
      typedef SkipList<Node*,NodePtrComparator> NodeSkipList;
      NodeSkipList _nodes;
    };


    class Binding {
    public:
      /** Number of the variable at this node */
      unsigned var;
      /** term at this node */
      TermList term;
      /** Create new binding */
      Binding(int v,TermList t) : var(v), term(t) {}
    }; // class SubstitutionTree::Binding

    typedef DHMap<unsigned,TermList,IdentityHash,DefaultHash> BindingMap;
    typedef Stack<unsigned> VarStack;

    Leaf* findLeaf(BindingMap& svBindings)
    { return _root ? findLeaf(_root, svBindings) : nullptr; }

    Leaf* findLeaf(Node* root, BindingMap& svBindings);

    void setSort(TypedTermList const& term, LeafData& ld)
    {
      ASS_EQ(ld.term, term)
      ld.sort = term.sort();
    }

    void setSort(TermList const& term, LeafData& ld)
    {
      ASS_EQ(ld.term, term)
      if (term.isTerm()) {
        ld.sort = SortHelper::getResultSort(term.term());
      }
    }


    void setSort(Literal* literal, LeafData &ld)
    { 
      ASS_EQ(ld.literal, literal); 
      if (literal->isEquality()) {
        ld.sort = SortHelper::getEqualityArgumentSort(literal);
      }
    }


    void handle(LeafData ld, bool doInsert)
    {
      auto norm = Renaming::normalize(ld.key());
      Recycled<BindingMap> bindings;
      createBindings(norm, /* reversed */ false,
          [&](int var, auto term) { 
            _nextVar = std::max(_nextVar, var + 1);
            bindings->insert(var, term);
          });
      if (doInsert) insert(*bindings, ld);
      else          remove(*bindings, ld);
    }

  private:
    void insert(BindingMap& binding,LeafData ld);
    void remove(BindingMap& binding,LeafData ld);

    /** Number of the next variable */
    int _nextVar = 0;
    Node* _root = nullptr;
    Cntr _iterCnt;

  public:

    class RenamingSubstitution 
    : public ResultSubstitution 
    {
    public:
      Recycled<Renaming> _query;
      Recycled<Renaming> _result;
      RenamingSubstitution(): _query(), _result() {}
      virtual ~RenamingSubstitution() override {}
      virtual TermList applyToQuery(TermList t) final override { return _query->apply(t); }
      virtual Literal* applyToQuery(Literal* l) final override { return _query->apply(l); }
      virtual TermList applyToResult(TermList t) final override { return _result->apply(t); }
      virtual Literal* applyToResult(Literal* l) final override { return _result->apply(l); }

      virtual TermList applyTo(TermList t, unsigned index) final override { ASSERTION_VIOLATION; }
      virtual Literal* applyTo(Literal* l, unsigned index) final override { NOT_IMPLEMENTED; }

      virtual size_t getQueryApplicationWeight(TermList t) final override { return t.weight(); }
      virtual size_t getQueryApplicationWeight(Literal* l) final override  { return l->weight(); }
      virtual size_t getResultApplicationWeight(TermList t) final override { return t.weight(); }
      virtual size_t getResultApplicationWeight(Literal* l) final override { return l->weight(); }

      void output(std::ostream& out) const final override
      { out << "{ _query: " << _query << ", _result: " << _result << " }"; }
    };

    template<class Query>
    bool generalizationExists(Query query)
    {
      return _root == nullptr 
        ? false
        : FastGeneralizationsIterator(this, _root, query, /* retrieveSubstitutions */ false, /* reversed */ false).hasNext();
    }

    template<class Query>
    VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>> getVariants(Query query, bool retrieveSubstitutions)
    {
      auto renaming = retrieveSubstitutions ? std::make_unique<RenamingSubstitution>() : std::unique_ptr<RenamingSubstitution>(nullptr);
      ResultSubstitutionSP resultSubst = retrieveSubstitutions ? ResultSubstitutionSP(&*renaming) : ResultSubstitutionSP();

      Query normQuery;
      if (retrieveSubstitutions) {
        renaming->_query->normalizeVariables(query);
        normQuery = renaming->_query->apply(query);
      } else {
        normQuery = Renaming::normalize(query);
      }

      Recycled<BindingMap> svBindings;
      createBindings(normQuery, /* reversed */ false,
          [&](auto v, auto t) { {
            svBindings->insert(v, t);
          } });
      Leaf* leaf = findLeaf(*svBindings);
      if(leaf==0) {
        return VirtualIterator<QueryRes<ResultSubstitutionSP, LeafData>>::getEmpty();
      } else {
        return pvi(iterTraits(leaf->allChildren())
          .map([retrieveSubstitutions, renaming = std::move(renaming), resultSubst](LeafData* ld) 
            {
              ResultSubstitutionSP subs;
              if (retrieveSubstitutions) {
                renaming->_result->reset();
                renaming->_result->normalizeVariables(ld->key());
                subs = resultSubst;
              }
              return QueryRes(subs, ld);
            }));
      }
    }

    class LeafIterator
    {
    public:
      LeafIterator(LeafIterator&&) = default;
      LeafIterator& operator=(LeafIterator&&) = default;
      DECL_ELEMENT_TYPE(Leaf*);
      LeafIterator(SubstitutionTree* st);
      bool hasNext()
      { return _curr != nullptr; }
      Leaf* next();
    private:
      void skipToNextLeaf();
      Node* _curr;
      Stack<NodeIterator> _nodeIterators;
    };



     /**
     * Class that supports matching operations required by
     * retrieval of generalizations in substitution trees.
     */
    class GenMatcher
    {
      static unsigned weight(Literal* l) { return l->weight(); }
      static unsigned weight(TermList t) { return  t.weight(); }
    public:
      GenMatcher(GenMatcher&&) = default;
      GenMatcher& operator=(GenMatcher&&) = default;



      template<class TermOrLit>
      void init(TermOrLit query, unsigned nextSpecVar)
      {
        _boundVars.reset();
        // _specVars->reset(); <- does not need to be reset as we only write before we read 
        _bindings.reset();

        _maxVar = weight(query) - 1;
        if(_specVars.size()<nextSpecVar) {
          //_specVars can get really big, but it was introduced instead of hash table
          //during optimizations, as it raised performance by abour 5%.
          _specVars.ensure(std::max(static_cast<unsigned>(_specVars.size()*2), nextSpecVar));
        }
        _bindings.ensure(weight(query));
      }

      /**
       * @b nextSpecVar Number higher than any special variable present in the tree.
       * 	It's used to determine size of the array that stores bindings of
       * 	special variables.
       */
      template<class TermOrLit>
      GenMatcher(TermOrLit query, unsigned nextSpecVar)
      { init(query,nextSpecVar); }

      USE_ALLOCATOR(GenMatcher);

      /**
       * Bind special variable @b var to @b term. This method
       * should be called only before any calls to @b matchNext()
       * and @b backtrack().
       */
      void bindSpecialVar(unsigned var, TermList term)
      {
        (_specVars)[var]=term;
      }
      /**
       * Return term bound to special variable @b specVar
       */
      TermList getSpecVarBinding(unsigned specVar)
      { return (_specVars)[specVar]; }

      bool matchNext(unsigned specVar, TermList nodeTerm, bool separate=true);
      bool matchNextAux(TermList queryTerm, TermList nodeTerm, bool separate=true);
      void backtrack();

      ResultSubstitutionSP getSubstitution(Renaming* resultNormalizer);

      int getBSCnt()
      {
        int res=0;
        VarStack::Iterator vsit(_boundVars);
        while(vsit.hasNext()) {
      if(vsit.next()==BACKTRACK_SEPARATOR) {
        res++;
      }
        }
        return res;
      }

    protected:
      static const unsigned BACKTRACK_SEPARATOR=0xFFFFFFFF;

      struct Binder;
      struct Applicator;
      class Substitution;

      VarStack _boundVars;
      DArray<TermList> _specVars;
      /**
       * Inheritors must assign the maximal possible number of an ordinary
       * variable that can be bound during the retrievall process.
       */
      unsigned _maxVar;

      /**
       * Inheritors must ensure that the size of this map will
       * be at least @b _maxVar+1
       */
      ArrayMap<TermList> _bindings;
    };

    /**
     * creates the bindings that need to be set for querying the given `term` from a substitution tree.
     * This means the root speical variables of the substitution tree need to be set to the right values.
     * In the case of a term this means { S0 -> term, S0 -> sortOfTerm  }
     * In the case of a literal this means { S0 -> arg0, ..., SN -> argN  }
     */
    template<class BindingFunction>
    void createBindings(TypedTermList term, bool reversed, BindingFunction bindSpecialVar)
    {
      bindSpecialVar(0, term);
      bindSpecialVar(1, term.sort());
    }

    /** see createBindings(TypedTermList,...) */
    template<class BindingFunction>
    void createBindings(Literal* lit, bool reversed, BindingFunction bindSpecialVar)
    {
      if (lit->isEquality()) {

        if (reversed) {
          bindSpecialVar(1,*lit->nthArgument(0));
          bindSpecialVar(0,*lit->nthArgument(1));
        } else {
          bindSpecialVar(0,*lit->nthArgument(0));
          bindSpecialVar(1,*lit->nthArgument(1));
        }

        bindSpecialVar(2, SortHelper::getEqualityArgumentSort(lit));
      } else {

        TermList* args=lit->args();
        int nextVar = 0;
        while (! args->isEmpty()) {
          unsigned var = nextVar++;
          bindSpecialVar(var,*args);
          args = args->next();
        }
      }
    }

    /**
     * Iterator, that yields generalizations of given term/literal.
     */
    class FastGeneralizationsIterator
    {
    public:

      FastGeneralizationsIterator(FastGeneralizationsIterator&&) = default;
      FastGeneralizationsIterator& operator=(FastGeneralizationsIterator&&) = default;
      DECL_ELEMENT_TYPE(QueryRes<ResultSubstitutionSP, LeafData>);
      using Unifier = ResultSubstitutionSP;

      void reset() {
        _iterCntr.reset();
        _resultNormalizer.reset();
        _alternatives.reset();
        _specVarNumbers.reset();
        _nodeTypes.reset();
      }

      template<class TermOrLit>
      void init(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed) {
        _retrieveSubstitution = retrieveSubstitution;
        _inLeaf = root->isLeaf();
        _subst.init(query, parent->_nextVar);
        _ldIterator = _inLeaf ? static_cast<Leaf*>(root)->allChildren() : LDIterator::getEmpty();
        _root = root;

        _iterCntr = parent->_iterCnt;
        ASS(root);

        parent->createBindings(query, reversed,
            [&](unsigned var, TermList t) { _subst.bindSpecialVar(var, t); });
      }

      /**
       * If @b reversed If true, parameters of supplied binary literal are
       * 	reversed. (useful for retrieval commutative terms)
       */
      template<class TermOrLit>
      FastGeneralizationsIterator(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed)
      : _subst(query, parent->_nextVar)
      { init(parent, root, query, retrieveSubstitution, reversed); }

      QueryRes<ResultSubstitutionSP, LeafData> next();
      bool hasNext();
    protected:

      bool findNextLeaf();
      bool enterNode(Node*& node);

      /** We should include substitutions in the results */
      bool _retrieveSubstitution;
      /** The iterator is currently in a leaf
       *
       * This is false in the beginning when it is in the root */
      bool _inLeaf;

      GenMatcher _subst;

      LDIterator _ldIterator;

      Renaming _resultNormalizer;

      Node* _root;

      Stack<void*> _alternatives;
      Stack<unsigned> _specVarNumbers;
      Stack<NodeAlgorithm> _nodeTypes;
      InstanceCntr _iterCntr;
    public:
      bool keepRecycled() const 
      { return _resultNormalizer.keepRecycled() || _alternatives.keepRecycled() || _specVarNumbers.keepRecycled() || _nodeTypes.keepRecycled(); }
    };


    /**
     * Class that supports matching operations required by
     * retrieval of generalizations in substitution trees.
     */
    class InstMatcher 
    {
    public:

      USE_ALLOCATOR(InstMatcher);

      void reset() {
        _boundVars.reset();
        _bindings.reset();
        _derefBindings.reset();
      }

      struct TermSpec
      {
        TermSpec() : q(false) {
        #if VDEBUG
          t = TermList::empty();
        #endif
        }
        TermSpec(bool q, TermList t)
        : q(q), t(t)
        {
          //query does not contain special vars
          ASS(!q || !t.isTerm() || t.term()->shared());
          ASS(!q || !t.isSpecialVar());
        }

        std::string toString()
        { return (q ? "q|" : "n|")+t.toString(); }

        /**
         * Return true if the @b t field can be use as a binding for a query
         * term variable in the retrieved substitution
         */
        bool isFinal()
        {
          //the fact that a term is shared means it does not contain any special variables
          return q
        ? (t.isTerm() && t.term()->ground())
        : (t.isOrdinaryVar() || (t.isTerm() && t.term()->shared()) );
        }

        bool q;
        TermList t;
      };

      /**
       * Bind special variable @b var to @b term
       *
       * This method should be called only before any calls to @b matchNext()
       * and @b backtrack().
       */
      void bindSpecialVar(unsigned var, TermList term)
      {
        ASS_EQ(getBSCnt(), 0);

        ALWAYS(_bindings.insert(TermList(var,true),TermSpec(true,term)));
      }

      bool isSpecVarBound(unsigned specVar)
      {
        return _bindings.find(TermList(specVar,true));
      }

      /** Return term bound to special variable @b specVar */
      TermSpec getSpecVarBinding(unsigned specVar)
      { return _bindings.get(TermList(specVar,true)); }

      bool findSpecVarBinding(unsigned specVar, TermSpec& res)
      {
        return _bindings.find(TermList(specVar,true), res);
      }

      bool matchNext(unsigned specVar, TermList nodeTerm, bool separate=true);
      bool matchNextAux(TermList queryTerm, TermList nodeTerm, bool separate=true);

      void backtrack();
      ResultSubstitutionSP getSubstitution(Renaming* resultDenormalizer);

      int getBSCnt()
      {
        int res=0;
        TermStack::Iterator vsit(_boundVars);
        while(vsit.hasNext()) {
          if(vsit.next().isEmpty()) {
      res++;
          }
        }
        return res;
      }

      void onLeafEntered()
      {
        _derefBindings.reset();
      }

    private:

      class Substitution;

      TermList derefQueryBinding(unsigned var);

      bool isBound(TermList var)
      {
        ASS(var.isVar());

        return _bindings.find(var);
      }
      void bind(TermList var, TermSpec trm)
      {
        ASS(!var.isOrdinaryVar() || !trm.q); //we do not bind ordinary vars to query terms

        ALWAYS(_bindings.insert(var, trm));
        _boundVars.push(var);
      }

      TermSpec deref(TermList var);

      typedef DHMap<TermList, TermSpec> BindingMap;
      typedef Stack<TermList> TermStack;

      /** Stacks of bindings made on each backtrack level. Backtrack
       * levels are separated by empty terms. */
      TermStack _boundVars;

      BindingMap _bindings;

      /**
       * A cache for bindings of variables to result terms
       *
       * The map is reset whenever we enter a new leaf
       */
      DHMap<TermList,TermList> _derefBindings;

      struct DerefTask
      {
        DerefTask(TermList var) : var(var) { trm.t = TermList::empty(); }
        DerefTask(TermList var, TermSpec trm) : var(var), trm(trm) {}
        TermList var;
        TermSpec trm;
        bool buildDerefTerm() { return trm.t.isNonEmpty(); };
      };

      struct DerefApplicator
      {
        DerefApplicator(InstMatcher* im, bool query) : query(query), im(im) {}
        TermList apply(unsigned var)
        {
          if(query) {
            return im->_derefBindings.get(TermList(var, false));
          }
          else {
      return TermList(var, false);
          }
        }
        TermList applyToSpecVar(unsigned specVar)
        {
          ASS(!query);

          return im->_derefBindings.get(TermList(specVar, true));
        }
      private:
        bool query;
        InstMatcher* im;
      };
    };

    /**
     * Iterator, that yields generalizations of given term/literal.
     */
    class FastInstancesIterator
    {
    public:
      FastInstancesIterator(FastInstancesIterator&&) = default;
      FastInstancesIterator& operator=(FastInstancesIterator&&) = default;
      DECL_ELEMENT_TYPE(QueryRes<ResultSubstitutionSP, LeafData>);
      using Unifier = ResultSubstitutionSP;

      void reset() {
        _iterCntr.reset();
        _alternatives.reset();
        _specVarNumbers.reset();
        _nodeTypes.reset();
      }

      template<class TermOrLit>
      void init(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed) {
        _retrieveSubstitution = retrieveSubstitution;
        _inLeaf = root->isLeaf();
        _ldIterator = _inLeaf ? static_cast<Leaf*>(root)->allChildren() : LDIterator::getEmpty();
        _root = root;
        _iterCntr = parent->_iterCnt;
        _subst.reset();

        ASS(root);

        parent->createBindings(query, reversed,
            [&](unsigned var, TermList t) { _subst.bindSpecialVar(var, t); });

        if (_inLeaf) {
          _subst.onLeafEntered(); //we reset the bindings cache
        }
      }

      /**
       * If @b reversed If true, parameters of supplied binary literal are
       * 	reversed. (useful for retrieval commutative terms)
       */
      template<class TermOrLit>
      FastInstancesIterator(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed)
      { init(parent,root,query,retrieveSubstitution,reversed); }


      bool hasNext();
      QueryRes<ResultSubstitutionSP, LeafData> next();
    protected:
      bool findNextLeaf();

      bool enterNode(Node*& node);

    private:

      bool _retrieveSubstitution;
      bool _inLeaf;
      LDIterator _ldIterator;

      InstMatcher _subst;

      Renaming _resultDenormalizer;
      Node* _root;

      Stack<void*> _alternatives;
      Stack<unsigned> _specVarNumbers;
      Stack<NodeAlgorithm> _nodeTypes;
      InstanceCntr _iterCntr;

    public:
      bool keepRecycled() const 
      { return _alternatives.keepRecycled() || _specVarNumbers.keepRecycled() || _nodeTypes.keepRecycled(); }
    };

    /** * A generic iterator over a substitution tree that can be used to retrieve elements stored in the tree that match a certain retrieval condition R.
     * The most simple retrieval condition would be `R(x) <=> x unifies with query`.
     * Similarly we could have retrieval conditions that, find instances, generalizations, unification with 
     * abstraction etc. 
     * 
     * The retrieval condition is computed by objects of the type RetrievalAlgorithm. All commonly used ones can be 
     * found in Indexing::RetrievalAlgorithms, which also documents the interface of these objects. 
     *
     * Notes:
     * - currently instantiation and generalization don't use this generic approach, but the optimized iterators
     *   `Fast*Iterator`, which we hopefully can refactor away in the future without any loss in performance.
     * - We do not use subtyping but parametric polymorphism for them, as subtyping polymorphsim would require us to 
     *   have the same element type for all of them, which is not what we want.
     */
    template<class RetrievalAlgorithm>
    class Iterator final
    {
    public:
      Iterator(Iterator&&) = default;
      Iterator& operator=(Iterator&&) = default;
      using Unifier = typename RetrievalAlgorithm::Unifier;
      DECL_ELEMENT_TYPE(QueryRes<Unifier, LeafData>);

      void reset() {
        _iterCntr.reset();
        _svStack.reset();
        _nodeIterators.reset();
        _bdStack.reset();
        if(_normalizationRecording) {
          _algo.bdDone();
          _normalizationRecording=false;
          _normalizationBacktrackData.backtrack();
        }
        while(_bdStack.isNonEmpty()) {
          _bdStack.pop().backtrack();
        }
      }


      template<class TermOrLit, class...AlgoArgs>
      void init(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed, AlgoArgs... args) {
        _algo.init(args...);
        _retrieveSubstitution = retrieveSubstitution;
        _leafData = {};
        _normalizationRecording = false;
        _iterCntr = InstanceCntr(parent->_iterCnt);

#define DEBUG_QUERY(lvl, ...) if (lvl < 0) DBG(__VA_ARGS__)
        if(!root) {
          return;
        }

        parent->createBindings(query, reversed, 
            [&](unsigned var, TermList t) { _algo.bindQuerySpecialVar(var, t); });
        DEBUG_QUERY(1, "query: ", _algo)


        prepareChildren(root, /* backtrackable */ false);
      }

      template<class TermOrLit, class...AlgoArgs>
      Iterator(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed, AlgoArgs... args)
       : _algo(args...)
      { init(parent, root, query, retrieveSubstitution, reversed, args...); }


      ~Iterator()
      { reset(); }

      bool hasLeafData() { return _leafData.isSome() && _leafData->hasNext(); };

      bool hasNext()
      {
        if(_normalizationRecording) {
          _algo.bdDone();
          _normalizationRecording=false;
          _normalizationBacktrackData.backtrack();
        }

        while(!hasLeafData() && findNextLeaf()) {}
        return hasLeafData();
      }

      QueryRes<Unifier, LeafData> next()
      {
        while(!hasLeafData() && findNextLeaf()) {}
        ASS(hasLeafData());

        ASS(!_normalizationRecording);

        auto ld = _leafData->next();
        if (_retrieveSubstitution) {
            Renaming normalizer;
            normalizer.normalizeVariables(ld->key());

            ASS(_normalizationBacktrackData.isEmpty());
            _algo.bdRecord(_normalizationBacktrackData);
            _normalizationRecording=true;

            _algo.denormalize(normalizer);
        }

        DEBUG_QUERY(1, "leaf data: ", *ld)
        return QueryRes(_algo.unifier(), ld);
      }

    private:

      template<class F>
      bool runRecording(F f) 
      {
        _algo.bdRecord(_bdStack.top());
        bool success = f();
        _algo.bdDone();
        return success;
      }

      bool inLeaf() const { return _leafData.isSome(); }

      bool findNextLeaf()
      {
        if(_nodeIterators.isEmpty()) {
          //There are no node iterators in the stack, so there's nowhere
          //to look for the next leaf.
          //This shouldn't hapen during the regular retrieval process, but it
          //can happen when there are no literals inserted for a predicate,
          //or when predicates with zero arity are encountered.
          ASS(_bdStack.isEmpty());
          return false;
        }

        auto leaveLeaf = [&]() {
            ASS(!_normalizationRecording);
            _bdStack.pop().backtrack();
            _leafData = {};
        };

        if(_leafData.isSome()) {
          leaveLeaf();
        }

        ASS(!_normalizationRecording);
        ASS(_bdStack.length()+1==_nodeIterators.length());

        do {
          while (!_nodeIterators.top().hasNext() && !_bdStack.isEmpty()) {
            _bdStack.pop().backtrack();
          }
          if(!_nodeIterators.top().hasNext()) {
            return false;
          }
          Node* n=*_nodeIterators.top().next();
          DEBUG_QUERY(1, "trying S", _svStack.top(), " -> ", n->term())

          _bdStack.push(BacktrackData());

          if (runRecording([&]() { return _algo.associate(_svStack.top(), n->term());})) {
            prepareChildren(n, /* backtrackable */ true);
            if (_leafData.isSome() && !runRecording([&](){ return _algo.doFinalLeafCheck(); })) {
              leaveLeaf();
              continue;
            }
          } else {
            _bdStack.pop().backtrack();
            continue;
          }
        } while(_leafData.isNone());
        ASS(_leafData.isSome())
        ASS(_bdStack.size() != 0)
        return true;
      }

      /** if `n` is a leaf _ldIterator is prepared 
       * if `n` is internal, the next special variable is put on svStack and the children it should be unified with are being put on _nodeIterators
       */
      void prepareChildren(Node* n, bool backtrackable) {
        if(n->isLeaf()) {
          _leafData = some(static_cast<Leaf*>(n)->allChildren());
        } else {
          IntermediateNode* inode=static_cast<IntermediateNode*>(n);
          _svStack.push(inode->childVar);
          _leafData = {};
          DEBUG_QUERY(1, "entering node: S", _svStack.top())
          
          _nodeIterators.push(_algo.template selectPotentiallyUnifiableChildren<LeafData>(inode));
          if (backtrackable) {
            _bdStack.top().addClosure([&]() { 
                DEBUG_CODE(auto var = )_svStack.pop();
                DEBUG_QUERY(1, "backtracking node: S", var)
                _nodeIterators.pop(); 
            });
          }
        }
      }

      RetrievalAlgorithm _algo;
      VarStack _svStack;
      bool _retrieveSubstitution;
      Option<LDIterator> _leafData;
      Stack<NodeIterator> _nodeIterators;
      Stack<BacktrackData> _bdStack;
      bool _normalizationRecording;
      BacktrackData _normalizationBacktrackData;
      InstanceCntr _iterCntr;

    public:
      bool keepRecycled() const 
      { return _svStack.keepRecycled() || _nodeIterators.keepRecycled() || _bdStack.keepRecycled(); }
    };


  public:
    bool maybeEmpty() const { return _root == nullptr; }
    bool isEmpty() const { return _root == nullptr || _root->isEmpty(); }
  }; // class SubstiutionTree

  /* This namespace defines classes to be used as type parameter for SubstitutionTree::Iterator. 
   * 
   * They implement different retrieval operations, which determine which elements from the tree shall 
   * be returned.  The most classic retrieval operations are retrieving all key terms that unify with 
   * a query term, or that a generalization/instance of a query term.
   *
   * The only non-classic one implemented currently is UnificaitonWithAbstraction, but one could also think 
   * of implementing finding variable variants, equal terms (modulo some theory), etc using the same interface.
   *
   * The interface all RetrievalAlgorithms must conform to is documented at the example of class RobUnification.
   * Retrieval from a SubstitutionTree is preformed incrementally. We start first by inserting some query terms. 
   * This is done by Algorithm::bindQuerySpecialVar.
   *
   */ 
  namespace RetrievalAlgorithms {

      class RobUnification { 
        RobSubstitution _subs;
      public:
        RobUnification() : _subs() {}

        void init() { _subs.reset(); }

        /** a witness that the returned term matches the retrieval condition 
         * (could be a substitution, variable renaming, etc.) 
         */
        using Unifier = ResultSubstitutionSP;  
        
        /** before starting to retrieve terms from the tree we insert some query terms, which we are going to 
         *  match the terms in the tree with. This is done using this function. */
        void bindQuerySpecialVar(unsigned var, TermList term)
        { _subs.bindSpecialVar(var, term, QUERY_BANK); }

        /** we intrementally traverse the tree, and at every code we call this retrieval algorithm to check 
         * whether it is okay to bind a new special variable to some term in the tree.
         * This function returns true if the retrieval condition (like in this case unifyability) can still 
         * be achieved, or false if not. Depending on that the iterator will backtrack or continue to traverse 
         * deeper into the tree.
         *
         * On insert into a substitution tree the inserted terms are first normalized (the names of the variables)
         * Therefore the namespace of the variables passed here is different from the ones of the actually inserted 
         * terms. 
         * Matching them up again is done by the function denormalize.
         */
        bool associate(unsigned specialVar, TermList node)
        { return _subs.unify(TermList(specialVar, /* special */ true), QUERY_BANK, node, NORM_RESULT_BANK); }


        /** @see associate */
        void denormalize(Renaming& norm)
        { _subs.denormalize(norm, NORM_RESULT_BANK,RESULT_BANK); }

        /** whenever we arrive at a leave we return the currrent witness for the current leave term to unify
         * with the query term. The unifier is queried using this function.  */
        Unifier unifier() { return ResultSubstitution::fromSubstitution(&_subs, QUERY_BANK, RESULT_BANK); }

        /** same as in @Backtrackable */
        void bdRecord(BacktrackData& bd) { _subs.bdRecord(bd); }

        /** same as in @Backtrackable */
        void bdDone() { _subs.bdDone(); }

        /** This function is called once when the iterator arrives at a leaf. 
         * The function can do a final check whether the current state of the retrieved witness (e.g. substitution) 
         * is really unifying or not. 
         * If it returns true the leaf is returned, if it returns false the leaf is filtered out.
         * This is useful in the case of unificaiton with abstraction, where we overapproximate the 
         * set of potential unifiers. With this function we can filter out unnecessary unifiers that would be 
         * sound but are not needed. For examples and a bit more of an explanation have a look at the paper
         * Refining Unification with Abstraction from LPAR2023
         */
        bool doFinalLeafCheck() { return true; }

        /** 
         * Returns an iterator over all child nodes of n that should be attempted for unification.
         * This is only an optimization. One could allways return n->allChildren(), but in the case 
         * of syntactic unificaiton we can use SkipList (i.e. n->childByTop(...)) in order to skip 
         * unnecessary unification attempts.
         */
        template<class LD>
        typename SubstitutionTree<LD>::NodeIterator selectPotentiallyUnifiableChildren(typename SubstitutionTree<LD>::IntermediateNode* n)
        { return _selectPotentiallyUnifiableChildren<LD>(n, _subs); }

        template<class LD>
        static typename SubstitutionTree<LD>::NodeIterator _selectPotentiallyUnifiableChildren(typename SubstitutionTree<LD>::IntermediateNode* n, RobSubstitution& subs)
        {
          unsigned specVar=n->childVar;
          auto top = subs.getSpecialVarTop(specVar);
          if(top.var()) {
            return n->allChildren();
          } else {
            auto** match = n->childByTop(top, /* canCreate */ false);
            if(match) {
              return pvi(concatIters(
                           getSingletonIterator(match),
                           n->variableChildren()));
            } else {
              return n->variableChildren();
            }
          }
        }
        friend std::ostream& operator<<(std::ostream& out, RobUnification const& self)
        { return out << self._subs; }

      };

      class UnificationWithAbstraction { 
        AbstractingUnifier _unif;
        bool _fixedPointIteration;
      public:
        UnificationWithAbstraction(AbstractionOracle ao, bool fixedPointIteration) 
          : _unif(AbstractingUnifier::empty(ao)) 
          , _fixedPointIteration(fixedPointIteration) 
        {}

        void init(AbstractionOracle ao, bool fixedPointIteration) { 
          _unif.init(ao);
          _fixedPointIteration = fixedPointIteration;
        }

        using Unifier = AbstractingUnifier*;

        bool associate(unsigned specialVar, TermList node)
        { return _unif.unify(TermList(specialVar, /* special */ true), QUERY_BANK, node, NORM_RESULT_BANK); }

        Unifier unifier()
        { return &_unif; }

        void bindQuerySpecialVar(unsigned var, TermList term)
        { _unif.subs().bindSpecialVar(var, term, QUERY_BANK); }

        void bdRecord(BacktrackData& bd)
        { _unif.subs().bdRecord(bd); }

        void bdDone()
        { _unif.subs().bdDone(); }

        void denormalize(Renaming& norm)
        { _unif.subs().denormalize(norm, NORM_RESULT_BANK,RESULT_BANK); }

        bool doFinalLeafCheck()
        { return !_fixedPointIteration || _unif.fixedPointIteration(); }

        template<class LD>
        static typename SubstitutionTree<LD>::NodeIterator _selectPotentiallyUnifiableChildren(typename SubstitutionTree<LD>::IntermediateNode* n, AbstractingUnifier& unif)
        {
          if (unif.usesUwa()) {
            unsigned specVar = n->childVar;
            auto top = unif.subs().getSpecialVarTop(specVar);

            if(top.var()) {
              return n->allChildren();
            } else {
              auto syms = unif.unifiableSymbols(*top.functor());
              if (syms) {
                return pvi(concatIters(
                      arrayIter(std::move(*syms))
                        .map   ([=](auto f) { return n->childByTop(top, /* canCreate */ false); })
                        .filter([ ](auto n) { return n != nullptr; }),
                      n->variableChildren()
                      ));
              } else {
                return n->allChildren();
              }
            }
          } else {
            return RobUnification::template _selectPotentiallyUnifiableChildren<LD>(n, unif.subs());
          }
        }

        template<class LD>
        typename SubstitutionTree<LD>::NodeIterator selectPotentiallyUnifiableChildren(typename SubstitutionTree<LD>::IntermediateNode* n)
        { return _selectPotentiallyUnifiableChildren<LD>(n, _unif); }
        friend std::ostream& operator<<(std::ostream& out, UnificationWithAbstraction const& self)
        { return out << self._unif; }
      };
    };

  } // namespace Indexing


#include "Indexing/SubstitutionTree.cpp"
#include "Indexing/SubstitutionTree_Nodes.cpp"
#include "Indexing/SubstitutionTree_FastGen.cpp"
#include "Indexing/SubstitutionTree_FastInst.cpp"

#endif // __SubstitutionTree__
