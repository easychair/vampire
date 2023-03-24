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
#include "Lib/Recycled.hpp"
#include "Kernel/BottomUpEvaluation/TypedTermList.hpp"

#include "Kernel/RobSubstitution.hpp"
#include "Kernel/Renaming.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/OperatorType.hpp"
#include "Lib/Option.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/ApplicativeHelper.hpp"

#include "Lib/Allocator.hpp"

#include "Index.hpp"

#if VDEBUG
#include <iostream>
#endif

using namespace std;
using namespace Lib;
using namespace Kernel;

#define UARR_INTERMEDIATE_NODE_MAX_SIZE 4

#define REORDERING 1

namespace Indexing {

  namespace UnificationAlgoritms {
    struct Unification { };
    struct Generalization { };
    struct Instantiation { };
  }

class SubstitutionTree;
std::ostream& operator<<(std::ostream& out, SubstitutionTree const& self);
std::ostream& operator<<(std::ostream& out, OutputMultiline<SubstitutionTree> const& self);

template<class Key> struct SubtitutionTreeConfig;

/** a counter that is compiled away in release mode */
struct Cntr {
#if VDEBUG
  Cntr() : self(0) {}
  int self;
  operator int() const { return self; }
#endif 
};

/** a reference to a Cntr that increments the counter when it is created and decrements it when it goes out of scope
 * This can be used to count the number of instances when an object of this type is added as a member field to the class 
 * that should be counted */
class InstanceCntr {
public:
#if VDEBUG
  Cntr& _cntr;

  InstanceCntr& operator=(InstanceCntr&& other) 
  { swap(other._cntr, _cntr); return *this; }

  InstanceCntr(InstanceCntr&& other) 
    : _cntr(other._cntr)
  { other._cntr.self++; }

  InstanceCntr(Cntr& cntr) : _cntr(cntr) 
  { _cntr.self++; }
  ~InstanceCntr() 
  { _cntr.self--; }
#else // VDEBUG
  InstanceCntr(Cntr& parent) {}
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
class SubstitutionTree
{

public:
  static constexpr int QRS_QUERY_BANK = 0;
  static constexpr int QRS_RESULT_BANK = 1;
  CLASS_NAME(SubstitutionTree);
  USE_ALLOCATOR(SubstitutionTree);

  SubstitutionTree(bool useC, bool rfSubs);

  SubstitutionTree& operator=(SubstitutionTree const& other) = delete;
  SubstitutionTree(SubstitutionTree&& other)
  : SubstitutionTree(other._useC, /* rfSubst */ false)
  {
    std::swap(_nextVar, other._nextVar);
    std::swap(_useC, other._useC);
    std::swap(_functionalSubtermMap, other._functionalSubtermMap);
    std::swap(_root, other._root);
#if VDEBUG
    std::swap(_tag, other._tag);
#endif
  }


  virtual ~SubstitutionTree();

  friend std::ostream& operator<<(std::ostream& out, SubstitutionTree const& self);
  friend std::ostream& operator<<(std::ostream& out, OutputMultiline<SubstitutionTree> const& self);

  struct LeafData {
    LeafData() {}

    LeafData(Clause* cls, Literal* literal, TypedTermList term, TermList extraTerm)
    : clause(cls), literal(literal), term(term), sort(term.sort()), extraTerm(extraTerm) {}
    LeafData(Clause* cls, Literal* literal, TypedTermList term)
    : clause(cls), literal(literal), term(term), sort(term.sort()) { extraTerm.makeEmpty();}

    LeafData(Clause* cls, Literal* literal, TermList term, TermList extraTerm)
    : clause(cls), literal(literal), term(term), extraTerm(extraTerm) { sort.makeEmpty(); }
    LeafData(Clause* cls, Literal* literal, TermList term)
    : clause(cls), literal(literal), term(term) {  sort.makeEmpty(); extraTerm.makeEmpty();}

    LeafData(Clause* cls, Literal* literal)
    : clause(cls), literal(literal) { term.makeEmpty(); sort.makeEmpty(), extraTerm.makeEmpty(); }
    inline
    bool operator==(const LeafData& o)
    { return clause==o.clause && literal==o.literal && term==o.term; }

    Clause* clause;
    Literal* literal;
    TermList term;
    TermList sort;
    // In some higher-order use cases, we want to store a different term 
    // in the leaf to the indexed term. extraTerm is used for this purpose.
    // In all other situations it is empty
    TermList extraTerm;

  };
  typedef VirtualIterator<LeafData&> LDIterator;

  struct QueryResult {
    LeafData const* data; 
    ResultSubstitutionSP subst;
    UnificationConstraintStackSP constr;

    QueryResult(LeafData const& ld) : data(&ld), subst(), constr() {};
    QueryResult(LeafData const& ld, ResultSubstitutionSP subst, UnificationConstraintStackSP constr) : data(&ld), subst(subst), constr(constr) {}
  };

  using QueryResultIterator = VirtualIterator<QueryResult>;
  // TODO make const function
  template<class Iterator, class TermOrLit> 
  QueryResultIterator iterator(TermOrLit query, bool retrieveSubstitutions, bool withConstraints, bool reversed = false)
  {
    CALL("TermSubstitutionTree::iterator");
    return _root == nullptr 
      ? QueryResultIterator::getEmpty()
      : pvi(iterTraits(Iterator(this, _root, query, retrieveSubstitutions, reversed, withConstraints, _functionalSubtermMap.asPtr() )));
  }

  class LDComparator
  {
  public:
    inline
    static Comparison compare(const LeafData& ld1, const LeafData& ld2)
    {
      CALL("SubstitutionTree::LDComparator::compare");

      if(ld1.clause && ld2.clause && ld1.clause!=ld2.clause) {
        ASS_NEQ(ld1.clause->number(), ld2.clause->number());
        return (ld1.clause->number()<ld2.clause->number()) ? LESS : GREATER;
      }
      Comparison res;
      if(ld1.literal && ld2.literal && ld1.literal!=ld2.literal) {
        res = (ld1.literal->getId()<ld2.literal->getId())? LESS : GREATER;
      } else {
        ASS_EQ(ld1.clause,ld2.clause);
        ASS_EQ(ld1.literal,ld2.literal);

        if (ld1.term.isEmpty()) {
          ASS(ld2.term.isEmpty());
          res = EQUAL;
        } else {
          if (ld1.term.isVar()) {
            if (ld2.term.isVar()) {
              unsigned var1 = ld1.term.var();
              unsigned var2 = ld2.term.var();
              res=(var1<var2)? LESS : (var1>var2)? GREATER : EQUAL;
            }
            else{
              res = LESS;
            }
          } else {
            if (ld2.term.isVar()) {
              res = GREATER;
            } else {
              unsigned id1 = ld1.term.term()->getId();
              unsigned id2 = ld2.term.term()->getId();
              res=(id1<id2)? LESS : (id1>id2)? GREATER : EQUAL;
            }
          }
        }
      }
      return res;
    }
  };

  enum NodeAlgorithm
  {
    UNSORTED_LIST=1,
    SKIP_LIST=2,
    SET=3
  };

  class Node {
  public:
    friend std::ostream& operator<<(ostream& out, OutputMultiline<Node> const& self) 
    { self.self.output(out, /* multiline = */ true, /* indent */ 0); return out; }
    friend std::ostream& operator<<(ostream& out, Node const& self) 
    { self.output(out, /* multiline = */ false, /* indent */ 0); return out; }
    inline
    Node() { term.makeEmpty(); }
    inline
    Node(TermList ts) : term(ts) { }
    virtual ~Node();
    /** True if a leaf node */
    virtual bool isLeaf() const = 0;
    virtual bool isEmpty() const = 0;
    virtual bool withSorts(){ return false; }
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
    virtual void makeEmpty() { term.makeEmpty(); }
    static void split(Node** pnode, TermList* where, int var);

#if VDEBUG
    virtual void assertValid() const {};
#endif

    /** term at this node */
    TermList term;

    virtual void output(std::ostream& out, bool multiline, int indent) const = 0;
  };


  typedef VirtualIterator<Node**> NodeIterator;
  class IntermediateNode;
    
    //We can remove this class once we deal with UWA uniformly for
    //for theories and HOL AYB
    class ChildBySortHelper
    {
    public:
        
        CLASS_NAME(SubstitutionTree::ChildBySortHelper);
        USE_ALLOCATOR(ChildBySortHelper);
        
        ChildBySortHelper(IntermediateNode* p):  _parent(p)
        {
            bySort.ensure(Signature::FIRST_USER_CON);
            bySortTerms.ensure(Signature::FIRST_USER_CON);
        }

        void loadFrom(ChildBySortHelper* other){
          ASS(other->bySort.size() == other->bySortTerms.size());
          for(unsigned i=0;i<other->bySort.size();i++){
            DHSet<unsigned>::Iterator it1(other->bySort[i]);
            bySort[i].loadFromIterator(it1);
            Stack<TermList>::Iterator it2(other->bySortTerms[i]);
            bySortTerms[i].loadFromIterator(it2);
          }
        } 
        
        /**
         * Return an iterator of child nodes whose top term has the same sort
         * as Termlist t. Only consider interpreted sorts.
         *
         */
        NodeIterator childBySort(TermList t)
        {
          CALL("SubstitutionTree::ChildBySortHelper::childBySort");
          TermList srt;
          // only consider interpreted sorts
          if(SortHelper::tryGetResultSort(t,srt) && !srt.isVar()){
            ASS(srt.isTerm());
            unsigned con = srt.term()->functor(); 
            if(!srt.term()->isSuper() && env.signature->isInterpretedNonDefault(con)){
              unsigned top = t.term()->functor();
              TermStack::Iterator fit(bySortTerms[con]);
              auto withoutThisTop = getFilteredIterator(fit,NotTop(top));
              auto nodes = getMappingIterator(withoutThisTop,ByTopFn(this));
              return pvi(getFilteredIterator(nodes,NonzeroFn()));
             }
          }
          return NodeIterator::getEmpty();
        } 
        
        DArray<DHSet<unsigned>> bySort;
        DArray<Stack<TermList>> bySortTerms;
        
        IntermediateNode* _parent; 
        /*
         * This is used for recording terms that might
         */
        void mightExistAsTop(TermList t)
        {
          CALL("SubstitutionTree::ChildBySortHelper::mightExistAsTop");
          if(!t.isTerm()){ return; }
          TermList srt;
          if(SortHelper::tryGetResultSort(t,srt) &&  !srt.isVar() && 
             !srt.term()->isSuper()){
            unsigned con = srt.term()->functor();
            if(env.signature->isInterpretedNonDefault(con)){
              unsigned f = t.term()->functor();
              if(bySort[con].insert(f)){
                bySortTerms[con].push(t);
              }
            }
          }
        }
        void remove(TermList t)
        {
          CALL("SubstitutionTree::ChildBySortHelper::remove");
          if(!t.isTerm()){ return;}
          TermList srt;
          if(SortHelper::tryGetResultSort(t,srt) && !srt.isVar() &&  
             !srt.term()->isSuper()){
            unsigned con = srt.term()->functor();
            if(env.signature->isInterpretedNonDefault(con)){
              unsigned f = t.term()->functor();
              if(bySort[con].remove(f)){
                bySortTerms[con].remove(t);
              }
            }
          }
        }
        
    };// class SubstitutionTree::ChildBySortHelper
    
    

  class IntermediateNode
    	: public Node
  {
  public:
    /** Build a new intermediate node which will serve as the root*/
    inline
    IntermediateNode(unsigned childVar) : childVar(childVar),_childBySortHelper(0) {}

    /** Build a new intermediate node */
    inline
    IntermediateNode(TermList ts, unsigned childVar) : Node(ts), childVar(childVar),_childBySortHelper(0) {}

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
    virtual Node** childByTop(TermList t, bool canCreate) = 0;


    /**
     * Remove child which points to node with top symbol of @b t.
     * This node has to still exist in time of the call to remove method.
     */
    virtual void remove(TermList t) = 0;

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

    
    virtual NodeIterator childBySort(TermList t)
    {
        if(!_childBySortHelper) return NodeIterator::getEmpty();
        return _childBySortHelper->childBySort(t);
    }
    virtual void mightExistAsTop(TermList t) {
        if(_childBySortHelper){
          _childBySortHelper->mightExistAsTop(t);
        }
    }

    void loadChildren(NodeIterator children);

    const unsigned childVar;
    ChildBySortHelper* _childBySortHelper;

    virtual void output(std::ostream& out, bool multiline, int indent) const override;
  }; // class SubstitutionTree::IntermediateNode

    struct ByTopFn
    {
        ByTopFn(ChildBySortHelper* n) : node(n) {};
        Node** operator()(TermList t){
            return node->_parent->childByTop(t,false);
        }
    private:
        ChildBySortHelper* node;
    }; 
    struct NotTop
    {
        NotTop(unsigned t) : top(t) {};
        bool operator()(TermList t){
            return t.term()->functor()!=top;
        }
    private:
        unsigned top;
    };
    

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
  static IntermediateNode* createIntermediateNode(unsigned childVar,bool constraints);
  static IntermediateNode* createIntermediateNode(TermList ts, unsigned childVar,bool constraints);
  static void ensureIntermediateNodeEfficiency(IntermediateNode** inode);

  struct IsPtrToVarNodeFn
  {
    bool operator()(Node** n)
    {
      return (*n)->term.isVar();
    }
  };

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
	destroyChildren();
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
    { return pvi( PointerPtrIterator<Node*>(&_nodes[0],&_nodes[_size]) ); }

    NodeIterator variableChildren()
    {
      return pvi( getFilteredIterator(PointerPtrIterator<Node*>(&_nodes[0],&_nodes[_size]),
  	    IsPtrToVarNodeFn()) );
    }
    virtual Node** childByTop(TermList t, bool canCreate);
    void remove(TermList t);

#if VDEBUG
    virtual void assertValid() const
    {
      ASS_ALLOC_TYPE(this,"SubstitutionTree::UArrIntermediateNode");
    }
#endif

    CLASS_NAME(SubstitutionTree::UArrIntermediateNode);
    USE_ALLOCATOR(UArrIntermediateNode);

    int _size;
    Node* _nodes[UARR_INTERMEDIATE_NODE_MAX_SIZE+1];
  };

  class UArrIntermediateNodeWithSorts
  : public UArrIntermediateNode
  {
  public:
   UArrIntermediateNodeWithSorts(unsigned childVar) : UArrIntermediateNode(childVar) {
     _childBySortHelper = new ChildBySortHelper(this);
   }
   UArrIntermediateNodeWithSorts(TermList ts, unsigned childVar) : UArrIntermediateNode(ts, childVar) {
     _childBySortHelper = new ChildBySortHelper(this);
   }
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
	destroyChildren();
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
#if VDEBUG
    virtual void assertValid() const
    {
      ASS_ALLOC_TYPE(this,"SubstitutionTree::SListIntermediateNode");
    }
#endif
    inline
    NodeIterator allChildren()
    {
      return pvi( NodeSkipList::PtrIterator(_nodes) );
    }
    inline
    NodeIterator variableChildren()
    {
      return pvi( getWhileLimitedIterator(
  		    NodeSkipList::PtrIterator(_nodes),
  		    IsPtrToVarNodeFn()) );
    }
    virtual Node** childByTop(TermList t, bool canCreate)
    {
      CALL("SubstitutionTree::SListIntermediateNode::childByTop");

      Node** res;
      bool found=_nodes.getPosition(t,res,canCreate);
      if(!found) {
        if(canCreate) {
          mightExistAsTop(t);
          *res=0;
        } else {
          res=0;
        }
      }
      return res;
    }
    inline
    void remove(TermList t)
    {
      _nodes.remove(t);
      if(_childBySortHelper){
        _childBySortHelper->remove(t);
      }
    }

    CLASS_NAME(SubstitutionTree::SListIntermediateNode);
    USE_ALLOCATOR(SListIntermediateNode);

    class NodePtrComparator
    {
    public:
      static Comparison compare(TermList t1,TermList t2)
      {
	CALL("SubstitutionTree::SListIntermediateNode::NodePtrComparator::compare");

	if(t1.isVar()) {
	  if(t2.isVar()) {
	    return Int::compare(t1.var(), t2.var());
	  }
	  return LESS;
	}
	if(t2.isVar()) {
	  return GREATER;
	}
	return Int::compare(t1.term()->functor(), t2.term()->functor());
      }

      static Comparison compare(Node* n1, Node* n2)
      { return compare(n1->term, n2->term); }
      static Comparison compare(TermList t1, Node* n2)
      { return compare(t1, n2->term); }
    };
    typedef SkipList<Node*,NodePtrComparator> NodeSkipList;
    NodeSkipList _nodes;
  };


  
  class SListIntermediateNodeWithSorts
  : public SListIntermediateNode
  {
   public:
   SListIntermediateNodeWithSorts(unsigned childVar) : SListIntermediateNode(childVar) {
       _childBySortHelper = new ChildBySortHelper(this);
   }
   SListIntermediateNodeWithSorts(TermList ts, unsigned childVar) : SListIntermediateNode(ts, childVar) {
       _childBySortHelper = new ChildBySortHelper(this);
   }
  };

  class Binding {
  public:
    /** Number of the variable at this node */
    unsigned var;
    /** term at this node */
    TermList term;
    /** Create new binding */
    Binding(int v,TermList t) : var(v), term(t) {}

    struct Comparator
    {
      inline
      static Comparison compare(Binding& b1, Binding& b2)
      {
    	return Int::compare(b2.var, b1.var);
      }
    };
  }; // class SubstitutionTree::Binding

  struct SpecVarComparator
  {
    inline
    static Comparison compare(unsigned v1, unsigned v2)
    { return Int::compare(v2, v1); }
    inline
    static unsigned max()
    { return 0u; }
  };

  typedef DHMap<unsigned,TermList,IdentityHash,DefaultHash> BindingMap;
  //Using BinaryHeap as a BindingQueue leads to about 30% faster insertion,
  //that when SkipList is used.
  typedef BinaryHeap<Binding,Binding::Comparator> BindingQueue;
  typedef BinaryHeap<unsigned,SpecVarComparator> SpecVarQueue;
  typedef Stack<unsigned> VarStack;

  void getBindingsArgBindings(Term* t, BindingMap& binding);

  Leaf* findLeaf(BindingMap& svBindings)
  { ASS(!_root || !_root->isLeaf() )
    return _root ? findLeaf(_root, svBindings) : nullptr; }

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


  template<class Key>
  void handle(Key const& key, LeafData ld, bool doInsert)
  {
    auto norm = Renaming::normalize(key);
    if (_functionalSubtermMap.isSome()) {
      norm = ApplicativeHelper::replaceFunctionalAndBooleanSubterms(norm, &_functionalSubtermMap.unwrap());
    }

    Recycled<BindingMap> bindings;
    setSort(key, ld);
    createBindings(norm, /* reversed */ false,
        [&](auto var, auto term) { 
          bindings->insert(var, term);
          _nextVar = max(_nextVar, (int)var + 1);
        });
    if (doInsert) insert(*bindings, ld);
    else          remove(*bindings, ld);
  }

private:
  void insert(BindingMap& binding,LeafData ld);
  void remove(BindingMap& binding,LeafData ld);

  /** Number of the next variable */
  int _nextVar;
  /** Array of nodes */
  /** enable searching with constraints for this tree */
  bool _useC;
  /** functional subterms of a term are replaced by extra sepcial
      variables before being inserted into the tree */
  Option<FuncSubtermMap> _functionalSubtermMap;
  Node* _root;
#if VDEBUG
  bool _tag;
#endif
public:
#if VDEBUG
  // Tags are used as a debug tool to turn debugging on for a particular instance
  virtual void markTagged(){ _tag=true;}
#endif

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
  QueryResultIterator getVariants(Query query, bool retrieveSubstitutions)
  {
    CALL("LiteralSubstitutionTree::getVariants");


    RenamingSubstitution* renaming = retrieveSubstitutions ? new RenamingSubstitution() : nullptr;
    ResultSubstitutionSP resultSubst = retrieveSubstitutions ? ResultSubstitutionSP(renaming) : ResultSubstitutionSP();

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
          _nextVar = max<int>(_nextVar, v + 1); // TODO do we need this line?
          svBindings->insert(v, t);
        } });
    Leaf* leaf = findLeaf(*svBindings);
    if(leaf==0) {
      return QueryResultIterator::getEmpty();
    } else {
      return pvi(iterTraits(leaf->allChildren())
        .map([retrieveSubstitutions, renaming, resultSubst](LeafData& ld) 
          {
            ResultSubstitutionSP subs;
            if (retrieveSubstitutions) {
              renaming->_result->reset();
              renaming->_result->normalizeVariables(SubtitutionTreeConfig<Query>::getKey(ld));
              subs = resultSubst;
            }
            return QueryResult(ld, subs, UnificationConstraintStackSP());
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
    bool hasNext();
    Leaf* next();
  private:
    void skipToNextLeaf();
    Node* _curr;
    Stack<NodeIterator> _nodeIterators;
  };



  // TODO document
  template<class BindingFunction>
  void createBindings(TypedTermList term, bool reversed, BindingFunction bindSpecialVar)
  {
    bindSpecialVar(0, term);
    bindSpecialVar(1, term.sort());
  }

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

    } else if(reversed) {
      ASS(lit->commutative());
      ASS_EQ(lit->arity(),2);

      bindSpecialVar(1,*lit->nthArgument(0));
      bindSpecialVar(0,*lit->nthArgument(1));

    } else if (lit->arity() == 0) {
      // insert a dummy term
      bindSpecialVar(0, TermList::var(0));

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

  class SubstitutionTreeMismatchHandler : public UWAMismatchHandler 
  {
  public:
    SubstitutionTreeMismatchHandler(Stack<UnificationConstraint>& c, BacktrackData& bd) : 
      UWAMismatchHandler(c), _constraints(c), _bd(bd) {}
    //virtual bool handle(RobSubstitution* subst, TermList query, unsigned index1, TermList node, unsigned index2);
  private:
    virtual bool introduceConstraint(TermList t1,unsigned index1, TermList t2,unsigned index2);
    Stack<UnificationConstraint>& _constraints;
    BacktrackData& _bd;
  };

  class STHOMismatchHandler : public HOMismatchHandler 
  {
  public:
    STHOMismatchHandler(Stack<UnificationConstraint>& c, BacktrackData& bd) : 
      HOMismatchHandler(c), _constraints(c), _bd(bd) {}
    virtual bool handle(RobSubstitution* subst, TermList query, unsigned index1, TermList node, unsigned index2);
  private:
    Stack<UnificationConstraint>& _constraints;
    BacktrackData& _bd;
  };  

  template<class UnificationAlgorithm>
  class UnificationsIterator final
  {
  public:
    UnificationsIterator(UnificationsIterator&&) = default;
    UnificationsIterator& operator=(UnificationsIterator&&) = default;
    DECL_ELEMENT_TYPE(QueryResult);

    template<class TermOrLit>
    UnificationsIterator(SubstitutionTree* parent, Node* root, TermOrLit query, bool retrieveSubstitution, bool reversed, bool useC, FuncSubtermMap* funcSubtermMap)
      : _subst()
      , _svStack()
      , _literalRetrieval(std::is_same<TermOrLit, Literal*>::value)
      , _retrieveSubstitution(retrieveSubstitution)
      , _inLeaf(false)
      , _ldIterator(LDIterator::getEmpty())
      , _nodeIterators()
      , _bdStack()
      , _clientBDRecording(false)
      , _useUWAConstraints(useC)
      , _useHOConstraints(funcSubtermMap)
      , _constraints()
      , _iterCntr(parent->_iterCnt)
#if VDEBUG
      , _tag(parent->_tag)
#endif
    {
#define DEBUG_QUERY(...) // DBG(__VA_ARGS__)
      CALL("SubstitutionTree::UnificationsIterator::UnificationsIterator");

      ASS(!_useUWAConstraints || retrieveSubstitution);
      ASS(!_useUWAConstraints || parent->_useC);

      if(!root) {
        return;
      }

      if(funcSubtermMap){
        _subst->setMap(funcSubtermMap);
        query = ApplicativeHelper::replaceFunctionalAndBooleanSubterms(query, funcSubtermMap);
      }

      parent->createBindings(query, reversed, 
          [&](unsigned var, TermList t) { _subst->bindSpecialVar(var, t, QUERY_BANK); });
      DEBUG_QUERY("query: ", subst)


      BacktrackData bd;
      enter(root, bd);
      bd.drop();
    }


    ~UnificationsIterator()
    {
      if(_clientBDRecording) {
        _subst->bdDone();
        _clientBDRecording=false;
        _clientBacktrackData.backtrack();
      }
      if (_bdStack) 
        while(_bdStack->isNonEmpty()) {
          _bdStack->pop().backtrack();
        }
    }

    bool hasNext()
    {
      CALL("SubstitutionTree::UnificationsIterator::hasNext");

      if(_clientBDRecording) {
        _subst->bdDone();
        _clientBDRecording=false;
        _clientBacktrackData.backtrack();
      }

      while(!_ldIterator.hasNext() && findNextLeaf()) {}
      return _ldIterator.hasNext();
    }


    QueryResult next()
    {
      CALL("SubstitutionTree::UnificationsIterator::next");

      while(!_ldIterator.hasNext() && findNextLeaf()) {}
      ASS(_ldIterator.hasNext());

      ASS(!_clientBDRecording);

      LeafData& ld=_ldIterator.next();

      if(_retrieveSubstitution) {
        Renaming normalizer;
        if(_literalRetrieval) {
          normalizer.normalizeVariables(ld.literal);
        } else {
          normalizer.normalizeVariables(ld.term);
          if (ld.sort.isNonEmpty()) {
            normalizer.normalizeVariables(ld.sort);
          }
        }

        ASS(_clientBacktrackData.isEmpty());
        _subst->bdRecord(_clientBacktrackData);
        _clientBDRecording=true;

        _subst->denormalize(normalizer,NORM_RESULT_BANK,RESULT_BANK);

        return QueryResult(ld, ResultSubstitution::fromSubstitution( &*_subst, QUERY_BANK, RESULT_BANK),
            // TODO do we really wanna copy the whole constraints stack here?
                UnificationConstraintStackSP(new UnificationConstraintStack(*_constraints))); 
      } else {
        return QueryResult(ld);
      }
    }

  protected:
    bool associate(TermList query, TermList node, BacktrackData& bd)
    {
      CALL("SubstitutionTree::UnificationsIterator::associate");

      //The ordering of the if statements is important here. Higher-order problems
      //should never require theory resoning (at the moment, theories cannot be parsed in HOL)
      //However, a user can still set UWA option on. We don't wan't that to result in 
      //the wrong handler being used.
      if(_useHOConstraints){
        STHOMismatchHandler hndlr(*_constraints,bd);
        return _subst->unify(query,QUERY_BANK,node,NORM_RESULT_BANK,&hndlr);    
      }
      if(_useUWAConstraints){ 
        SubstitutionTreeMismatchHandler hndlr(*_constraints,bd);
        return _subst->unify(query,QUERY_BANK,node,NORM_RESULT_BANK,&hndlr);
      } 
      return _subst->unify(query,QUERY_BANK,node,NORM_RESULT_BANK);
    }

    NodeIterator getNodeIterator(IntermediateNode* n)
    {
      CALL("SubstitutionTree::UnificationsIterator::getNodeIterator");

      unsigned specVar=n->childVar;
      TermList qt=_subst->getSpecialVarTop(specVar);
      if(qt.isVar()) {
        return n->allChildren();
      } else {
        Node** match=n->childByTop(qt, false);
        if(match) {
          return pvi( 
            getConcatenatedIterator(
         getSingletonIterator(match),
         n->variableChildren() 
           ));
        } else {
          return n->variableChildren();
        }
      }
    }

    bool findNextLeaf()
    {
      CALL("SubstitutionTree::UnificationsIterator::findNextLeaf");

      if(_nodeIterators->isEmpty()) {
        //There are no node iterators in the stack, so there's nowhere
        //to look for the next leaf.
        //This shouldn't hapen during the regular retrieval process, but it
        //can happen when there are no literals inserted for a predicate,
        //or when predicates with zero arity are encountered.
        ASS(_bdStack->isEmpty());
        return false;
      }

      if(_inLeaf) {
        ASS(!_clientBDRecording);
        //Leave the current leaf
        _bdStack->pop().backtrack();
        _inLeaf=false;
      }

      ASS(!_clientBDRecording);
      ASS(_bdStack->length()+1==_nodeIterators->length());

      do {
        while(!_nodeIterators->top().hasNext() && !_bdStack->isEmpty()) {
          //backtrack undos everything that enter(...) method has done,
          //so it also pops one item out of the nodeIterators stack
          _bdStack->pop().backtrack();
          _svStack->pop();
        }
        if(!_nodeIterators->top().hasNext()) {
          return false;
        }
        Node* n=*_nodeIterators->top().next();

        BacktrackData bd;
        bool success=enter(n,bd);
        if(!success) {
          bd.backtrack();
          continue;
        } else {
          _bdStack->push(bd);
        }
      } while(!_inLeaf);
      return true;
    }
    bool enter(Node* n, BacktrackData& bd)
    {
      CALL("SubstitutionTree::UnificationsIterator::enter");

#if VDEBUG
      if(_tag){
        cout << "=========================================" << endl;
        cout << "entering..." << *n << endl;
        cout << "subst is " << _subst << endl;
        cout << "svstack is " << _svStack << endl;
        cout << "=========================================" << endl;
      } 
#endif

      bool success=true;
      bool recording=false;
      if(!n->term.isEmpty()) {
        //n is proper node, not a root

        TermList qt(_svStack->top(), true);

        recording=true;
        _subst->bdRecord(bd);
        success=associate(qt,n->term,bd);
      }
      if(success) {
        if(n->isLeaf()) {
          _ldIterator=static_cast<Leaf*>(n)->allChildren();
          _inLeaf=true;
        } else {
          IntermediateNode* inode=static_cast<IntermediateNode*>(n);
          _svStack->push(inode->childVar);
          NodeIterator nit=getNodeIterator(inode);
          if(_useUWAConstraints){
            TermList qt = _subst->getSpecialVarTop(inode->childVar);
            NodeIterator enit = pvi(getConcatenatedIterator(inode->childBySort(qt),nit));
            _nodeIterators->backtrackablePush(enit,bd);
          }
          else{
            _nodeIterators->backtrackablePush(nit, bd);
          }
        }
      }
      if(recording) {
        _subst->bdDone();
      }
      return success;
    }


    static const int QUERY_BANK=0;
    static const int RESULT_BANK=1;
    static const int NORM_RESULT_BANK=3;

    Recycled<RobSubstitution> _subst;
    Recycled<VarStack> _svStack;

  private:
    bool _literalRetrieval;
    bool _retrieveSubstitution;
    bool _inLeaf;
    LDIterator _ldIterator;
    Recycled<Stack<NodeIterator>> _nodeIterators;
    Recycled<Stack<BacktrackData>> _bdStack;
    bool _clientBDRecording;
    BacktrackData _clientBacktrackData;
    bool _useUWAConstraints;
    bool _useHOConstraints;
    Recycled<UnificationConstraintStack> _constraints;
    InstanceCntr _iterCntr;
#if VDEBUG
    bool _tag;
#endif
  };
  using FastInstancesIterator       = UnificationsIterator<UnificationAlgoritms::Instantiation>;
  using FastGeneralizationsIterator = UnificationsIterator<UnificationAlgoritms::Generalization>;


  template<class Query>
  bool generalizationExists(Query query)
  {
    return _root == nullptr 
      ? false
      : FastGeneralizationsIterator(this, _root, query, /* retrieveSubstitutions */ false, /* reversed */ false, /* useC */ false, /* funcSubtermMap */ nullptr).hasNext();
  }

#if VDEBUG
public:
  bool isEmpty() const;
#endif

  Cntr _iterCnt;
}; // class SubstiutionTree

template<> 
struct SubtitutionTreeConfig<Literal*> 
{
  static Literal* const& getKey(SubstitutionTree::LeafData const& ld)
  { return ld.literal;  }
};


template<> 
struct SubtitutionTreeConfig<TermList> 
{
  static TermList const& getKey(SubstitutionTree::LeafData const& ld)
  { return ld.term;  }
};





} // namespace Indexing

#endif
