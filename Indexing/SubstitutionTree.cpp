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
 * @file SubstitutionTree.cpp
 * Implements class SubstitutionTree.
 *
 * @since 16/08/2008 flight Sydney-San Francisco
 */

#define DEBUG_INSERT(lvl, ...) if (lvl < 0) DBG(__VA_ARGS__)
#define DEBUG_REMOVE(lvl, ...) if (lvl < 0) DBG(__VA_ARGS__)
#include <utility>

#include "Shell/Options.hpp"

#include "Kernel/Matcher.hpp"
#include "Kernel/Renaming.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/ApplicativeHelper.hpp"

#include "Lib/BinaryHeap.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Recycled.hpp"
#include "Lib/DHMultiset.hpp"

#include "TermSharing.hpp"

#include <iostream>
#include "Debug/Tracer.hpp"
#if VDEBUG
#include "Kernel/Signature.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#endif

#include "SubstitutionTree.hpp"

namespace Indexing {
struct UnresolvedSplitRecord
{
  UnresolvedSplitRecord() {}
  UnresolvedSplitRecord(unsigned var, TermList original)
  : var(var), original(original) {}

  unsigned var;
  TermList original;
};

template<class LeafData_>
struct BindingComparator
{
  static Comparison compare(const UnresolvedSplitRecord& r1, const UnresolvedSplitRecord& r2)
  {
    bool r1HasSpecVars=r1.original.isTerm() && !r1.original.term()->shared();
    bool r2HasSpecVars=r2.original.isTerm() && !r2.original.term()->shared();
    if( r1HasSpecVars && !r2HasSpecVars ) {
      return GREATER;
    }
    if( r2HasSpecVars && !r1HasSpecVars ) {
      return LESS;
    }
    return Int::compare(r2.var,r1.var);
  }
  static Comparison compare(const typename SubstitutionTree<LeafData_>::Binding& b1, const typename SubstitutionTree<LeafData_>::Binding& b2)
  {
#if REORDERING
    return Int::compare(b2.var,b1.var);
#else
    return Int::compare(b1.var,b2.var);
#endif
  }
};


/**
 * Insert an entry to the substitution tree.
 *
 * @b pnode is pointer to root of tree corresponding to
 * top symbol of the term/literal being inserted, and
 * @b bh contains its arguments.
 */
template<class LeafData_>
void SubstitutionTree<LeafData_>::insert(BindingMap& svBindings, LeafData ld)
{
  ASS_EQ(_iterCnt,0);
  auto pnode = &_root;
  DEBUG_INSERT(0, "insert: ", svBindings, " into ", *this)

  if(*pnode == 0) {
    if (svBindings.isEmpty()) {
      auto leaf = createLeaf();
      leaf->insert(std::move(ld));
      *pnode = leaf;
      DEBUG_INSERT(0, "out: ", *this);
      return;
    } else {
      *pnode=createIntermediateNode(svBindings.getOneKey());
    }
  }
  if(svBindings.isEmpty()) {
    ASS((*pnode)->isLeaf());
    ensureLeafEfficiency(reinterpret_cast<Leaf**>(pnode));
    static_cast<Leaf*>(*pnode)->insert(ld);
    DEBUG_INSERT(0, "out: ", *this);
    return;
  }

  typedef BinaryHeap<UnresolvedSplitRecord, BindingComparator<LeafData_>> SplitRecordHeap;
  static SplitRecordHeap unresolvedSplits;
  unresolvedSplits.reset();

  ASS((*pnode));
  ASS(!(*pnode)->isLeaf());

start:

#if REORDERING
  ASS(!(*pnode)->isLeaf() || !unresolvedSplits.isEmpty());
  bool canPostponeSplits=false;
  if((*pnode)->isLeaf() || (*pnode)->algorithm()!=UNSORTED_LIST) {
    canPostponeSplits=false;
  } else {
    UArrIntermediateNode* inode = static_cast<UArrIntermediateNode*>(*pnode);
    canPostponeSplits = inode->size()==1;
    if(canPostponeSplits) {
      unsigned boundVar=inode->childVar;
      Node* child=inode->_nodes[0];
      bool removeProblematicNode=false;
      if(svBindings.find(boundVar)) {
	TermList term=svBindings.get(boundVar);
	bool wouldDescendIntoChild = inode->childByTop(term.top(),false)!=0;
	ASS_EQ(wouldDescendIntoChild, TermList::sameTop(term, child->term()));
	if(!wouldDescendIntoChild) {
	  //if we'd have to perform all postponed splitting due to
	  //node with a single child, we rather remove that node
	  //from the tree and deal with the binding, it represented,
	  //later.
	  removeProblematicNode=true;
	}
      } else if(!child->term().isTerm() || child->term().term()->shared()) {
	//We can remove nodes binding to special variables undefined in our branch
	//of the tree, as long as we're sure, that during split resolving we put these
	//binding nodes below nodes that define spec. variables they bind.
	removeProblematicNode=true;
      } else {
	canPostponeSplits = false;
      }
      if(removeProblematicNode) {
	unresolvedSplits.insert(UnresolvedSplitRecord(inode->childVar, child->term()));
	child->setTerm(inode->term());
	*pnode=child;
	inode->makeEmpty();
	delete inode;
	goto start;
      }
    }
  }
  canPostponeSplits|=unresolvedSplits.isEmpty();
  if(!canPostponeSplits) {

    while(!unresolvedSplits.isEmpty()) {
      UnresolvedSplitRecord urr=unresolvedSplits.pop();

      Node* node=*pnode;
      IntermediateNode* newNode = createIntermediateNode(node->term(), urr.var);
      node->setTerm(urr.original);

      *pnode=newNode;

      Node** nodePosition=newNode->childByTop(node->top(), true);
      ASS(!*nodePosition);
      *nodePosition=node;
    }
  }
#endif
  ASS(!(*pnode)->isLeaf());

  IntermediateNode* inode = static_cast<IntermediateNode*>(*pnode);
  ASS(inode);

  unsigned boundVar=inode->childVar;
  TermList term=svBindings.get(boundVar);
  svBindings.remove(boundVar);

  //Into pparent we store the node, we might be inserting into.
  //So in the case we do insert, we might check whether this node
  //needs expansion.
  Node** pparent=pnode;
  pnode=inode->childByTop(term.top(),true);

  if (*pnode == 0) {
    BindingMap::Iterator svit(svBindings);
    BinaryHeap<Binding, BindingComparator<LeafData_>> remainingBindings;
    while (svit.hasNext()) {
      unsigned var;
      TermList term;
      svit.next(var, term);
      remainingBindings.insert(Binding(var, term));
    }
    while (!remainingBindings.isEmpty()) {
      Binding b=remainingBindings.pop();
      IntermediateNode* inode = createIntermediateNode(term, b.var);
      term=b.term;

      *pnode = inode;
      pnode = inode->childByTop(term.top(),true);
    }
    Leaf* lnode=createLeaf(term);
    *pnode=lnode;
    lnode->insert(ld);

    ensureIntermediateNodeEfficiency(reinterpret_cast<IntermediateNode**>(pparent));
    DEBUG_INSERT(0, "out: ", *this);
    return;
  }


  TermList* tt = &term;
  TermList* ss = &(*pnode)->term();

  ASS(TermList::sameTop(*ss, *tt));


  // ss is the term in node, tt is the term to be inserted
  // ss and tt have the same top symbols but are not equal
  // create the common subterm of ss,tt and an alternative node
  Stack<TermList*> subterms(64);
  for (;;) {
    if (*tt!=*ss && TermList::sameTop(*ss,*tt)) {
      // ss and tt have the same tops and are different, so must be non-variables
      ASS(! ss->isVar());
      ASS(! tt->isVar());

      Term* s = ss->term();
      Term* t = tt->term();

      ASS(s->arity() > 0);
      ASS(s->functor() == t->functor());

      if (s->shared()) {
        // create a shallow copy of s
        s = Term::cloneNonShared(s);
        ss->setTerm(s);
      }

      ss = s->args();
      tt = t->args();
      if (ss->next()->isEmpty()) {
        continue;
      }
      subterms.push(ss->next());
      subterms.push(tt->next());
    } else {
      if (! TermList::sameTop(*ss,*tt)) {
        unsigned x;
        if(!ss->isSpecialVar()) {
          x = _nextVar++;
        #if REORDERING
          unresolvedSplits.insert(UnresolvedSplitRecord(x,*ss));
          ss->makeSpecialVar(x);
        #else
          Node::split(pnode,ss,x);
        #endif
        } else {
          x=ss->var();
        }
        svBindings.set(x,*tt);
      }

      if (subterms.isEmpty()) {
        break;
      }
      tt = subterms.pop();
      ss = subterms.pop();
      if (! ss->next()->isEmpty()) {
        subterms.push(ss->next());
        subterms.push(tt->next());
      }
    }
  }

  if (svBindings.isEmpty()) {
    ASS((*pnode)->isLeaf());
    ensureLeafEfficiency(reinterpret_cast<Leaf**>(pnode));
    Leaf* leaf = static_cast<Leaf*>(*pnode);
    leaf->insert(ld);
    DEBUG_INSERT(0, "out: ", *this);
    return;
  }

  goto start;
} // // SubstitutionTree<LeafData_>::insert

/*
 * Remove an entry from the substitution tree.
 *
 * @b pnode is pointer to root of tree corresponding to
 * top symbol of the term/literal being removed, and
 * @b bh contains its arguments.
 *
 * If the removal results in a chain of nodes containing
 * no terms/literals, all those nodes are removed as well.
 */
template<class LeafData_>
void SubstitutionTree<LeafData_>::remove(BindingMap& svBindings, LeafData ld)
{
  ASS_EQ(_iterCnt,0);
  auto pnode = &_root;
  DEBUG_REMOVE(0, "remove: ", svBindings, " from ", *this)

  ASS(*pnode);

  Recycled<Stack<Node**>> history;

  while (! (*pnode)->isLeaf()) {
    history->push(pnode);

    IntermediateNode* inode=static_cast<IntermediateNode*>(*pnode);

    unsigned boundVar=inode->childVar;
    TermList t = svBindings.get(boundVar);

    pnode=inode->childByTop(t.top(),false);
    ASS(pnode);


    TermList* s = &(*pnode)->term();
    ASS(TermList::sameTop(*s,t));

    if(*s==t) {
      continue;
    }

    ASS(! s->isVar());
    TermList* ss = s->term()->args();
    ASS(!ss->isEmpty());

    // computing the disagreement set of the two terms
    Stack<TermList*> subterms(120);

    subterms.push(ss);
    subterms.push(t.term()->args());
    while (! subterms.isEmpty()) {
      TermList* tt = subterms.pop();
      ss = subterms.pop();
      if (tt->next()->isEmpty()) {
	ASS(ss->next()->isEmpty());
      }
      else {
	subterms.push(ss->next());
	subterms.push(tt->next());
      }
      if (*ss==*tt) {
	continue;
      }
      if (ss->isVar()) {
	ASS(ss->isSpecialVar());
	svBindings.set(ss->var(),*tt);
	continue;
      }
      ASS(! tt->isVar());
      ASS(ss->term()->functor() == tt->term()->functor());
      ss = ss->term()->args();
      if (! ss->isEmpty()) {
	ASS(! tt->term()->args()->isEmpty());
	subterms.push(ss);
	subterms.push(tt->term()->args());
      }
    }
  }

  ASS ((*pnode)->isLeaf());


  Leaf* lnode = static_cast<Leaf*>(*pnode);
  lnode->remove(ld);
  ensureLeafEfficiency(reinterpret_cast<Leaf**>(pnode));

  while( (*pnode)->isEmpty() ) {
    TermList term=(*pnode)->term();
    if(history->isEmpty()) {
      delete *pnode;
      *pnode=0;
      DEBUG_REMOVE(0, "out: ", *this);
      return;
    } else {
      Node* node=*pnode;
      IntermediateNode* parent=static_cast<IntermediateNode*>(*history->top());
      parent->remove(term.top());
      delete node;
      pnode = history->pop();
      ensureIntermediateNodeEfficiency(reinterpret_cast<IntermediateNode**>(pnode));
    }
  }
  DEBUG_REMOVE(0, "out: ", *this);
} // SubstitutionTree<LeafData_>::remove

/**
 * Return a pointer to the leaf that contains term specified by @b svBindings.
 * If no such leaf exists, return 0.
 */
template<class LeafData_>
typename SubstitutionTree<LeafData_>::Leaf* SubstitutionTree<LeafData_>::findLeaf(Node* root, BindingMap& svBindings)
{
  ASS(root);

  Node* node=root;

  while (! node->isLeaf()) {
    IntermediateNode* inode=static_cast<IntermediateNode*>(node);

    unsigned boundVar=inode->childVar;
    TermList t = svBindings.get(boundVar);

    Node** child=inode->childByTop(t.top(),false);
    if(!child) {
      return 0;
    }
    node=*child;


    TermList s = node->term();
    ASS(TermList::sameTop(s,t));

    if(s==t) {
      continue;
    }

    ASS(! s.isVar());
    TermList* ss = s.term()->args();
    ASS(!ss->isEmpty());

    // computing the disagreement set of the two terms
    Stack<TermList*> subterms(120);

    subterms.push(ss);
    subterms.push(t.term()->args());
    while (! subterms.isEmpty()) {
      TermList* tt = subterms.pop();
      ss = subterms.pop();
      if (tt->next()->isEmpty()) {
	ASS(ss->next()->isEmpty());
      }
      else {
	subterms.push(ss->next());
	subterms.push(tt->next());
      }
      if (*ss==*tt) {
	continue;
      }
      if (ss->isSpecialVar()) {
	svBindings.set(ss->var(),*tt);
	continue;
      }
      if(ss->isVar() || tt->isVar() || ss->term()->functor()!=tt->term()->functor()) {
	return 0;
      }
      ss = ss->term()->args();
      if (! ss->isEmpty()) {
	ASS(! tt->term()->args()->isEmpty());
	subterms.push(ss);
	subterms.push(tt->term()->args());
      }
    }
  }
  ASS(node->isLeaf());
  return static_cast<Leaf*>(node);
}


template<class LeafData_>
SubstitutionTree<LeafData_>::Node::~Node()
{
  if(term().isTerm()) {
    term().term()->destroyNonShared();
  }
}


template<class LeafData_>
void SubstitutionTree<LeafData_>::Node::split(Node** pnode, TermList* where, int var)
{
  Node* node=*pnode;

  IntermediateNode* newNode = createIntermediateNode(node->term(), var);
  node->setTerm(*where);
  *pnode=newNode;

  where->makeSpecialVar(var);

  Node** nodePosition=newNode->childByTop(node->top(), true);
  ASS(!*nodePosition);
  *nodePosition=node;
}

template<class LeafData_>
void SubstitutionTree<LeafData_>::IntermediateNode::loadChildren(NodeIterator children)
{
  while(children.hasNext()) {
    Node* ext=*children.next();
    Node** own=childByTop(ext->top(), true);
    ASS(! *own);
    *own=ext;
  }
}

template<class LeafData_>
void SubstitutionTree<LeafData_>::Leaf::loadChildren(LDIterator children)
{
  while(children.hasNext()) {
    // TODO move instead of copying here
    insert(*children.next());
  }
}
template<class LeafData_>
SubstitutionTree<LeafData_>::LeafIterator::LeafIterator(SubstitutionTree* st)
  : _curr()
  , _nodeIterators()
{
  if (st->_root->isLeaf()) {
    _curr = st->_root;
  } else {
    _curr = nullptr;
    _nodeIterators.push(static_cast<IntermediateNode*>(st->_root)->allChildren());
  }
}

template<class LeafData_>
typename SubstitutionTree<LeafData_>::Leaf* SubstitutionTree<LeafData_>::LeafIterator::next()
{
  ASS(_curr->isLeaf());
  auto out = _curr;
  skipToNextLeaf();
  return static_cast<Leaf*>(out);
}


template<class LeafData_>
void SubstitutionTree<LeafData_>::LeafIterator::skipToNextLeaf()
{
  for (;;) {
    while(!_nodeIterators.isEmpty() && !_nodeIterators.top().hasNext()) {
      _nodeIterators.pop();
      _curr = nullptr;
    }
    if (_nodeIterators.isEmpty()) {
      ASS_EQ(_curr,0)
      return;
    } else {
      _curr = *_nodeIterators.top().next();
      if (_curr->isLeaf()) {
        return;
      } else {
        _nodeIterators.push(static_cast<IntermediateNode*>(_curr)->allChildren());
        _curr = nullptr;
      }
    }
  }
}

template<class LeafData_>
void SubstitutionTree<LeafData_>::Leaf::output(std::ostream& out, bool multiline, int indent) const 
{ out << this->term(); }

template<class LeafData_>
void SubstitutionTree<LeafData_>::IntermediateNode::output(std::ostream& out, bool multiline, int indent) const 
{
  // TODO const version of allChildren
  auto childIter = iterTraits(((IntermediateNode*)this)->allChildren());
  if (!this->term().isEmpty()) {
    out << this->term()
        << " ; ";
  }
  out << "S" << this->childVar << " -> ";

  auto first = childIter.next();
  auto brackets = childIter.hasNext();
  if (brackets) {


    if (multiline) {
      auto outp = [&](Node** x) { 
        out << std::endl; 
        Output::Multiline<int>::outputIndent(out, indent + 1);
        out << "| ";
        (*x)->output(out, multiline, indent + 1);
      };
      out << "[";
      outp(first);
      while (childIter.hasNext()) {
        outp(childIter.next());
      }
      out << std::endl; 
      Output::Multiline<int>::outputIndent(out, indent + 1);
      out << "]";

    } else {
      out << "[ ";
      out << **first;
      while (childIter.hasNext()) {
        out <<  " | " << **childIter.next();
      } 
      out << " ]";
    }


  } else {
    (*first)->output(out, multiline, indent);
  }

}

} // namespace Indexing
