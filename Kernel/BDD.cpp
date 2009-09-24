/**
 * @file Clause.cpp
 * Implements class BDD for binary decision diagrams
 */

#include <utility>

#include "../Lib/Environment.hpp"
#include "../Lib/Exception.hpp"
#include "../Lib/DHMap.hpp"
#include "../Lib/Int.hpp"
#include "../Lib/List.hpp"
#include "../Lib/Stack.hpp"
#include "../Lib/Timer.hpp"
#include "../Lib/TimeCounter.hpp"

#include "../Kernel/Signature.hpp"

#include "../SAT/ClauseSharing.hpp"
#include "../SAT/Preprocess.hpp"
#include "../SAT/SATClause.hpp"
#include "../SAT/SATLiteral.hpp"
#include "../SAT/SingleWatchSAT.hpp"

#include "../Shell/Options.hpp"

#include "BDD.hpp"


namespace Kernel {

using namespace Lib;
using namespace SAT;

/**
 * Return the singleton instance of the BDD class
 */
BDD* BDD::instance()
{
  CALL("BDD::instance");

  static BDD* inst=0;
  if(!inst) {
    inst=new BDD();
  }
  return inst;
}

/**
 * Create a new BDD object
 */
BDD::BDD()
: _trueNode(-1,0,0), _falseNode(-1,0,0), _newVar(1)
{
}

/**
 * Destroy a BDD object
 */
BDD::~BDD()
{
  CALL("BDD::~BDD");

  NodeSet::Iterator nit(_nodes);
  while(nit.hasNext()) {
    BDDNode* node=nit.next();
    delete node;
  }
}

/**
 * Return a new BDD variable that will represent propositional
 * predicate symbol @b pred
 */
int BDD::getNewVar(unsigned pred)
{
  CALL("BDD::getNewVar(unsigned)");
  ASS_EQ(env.signature->predicateArity(pred), 0);

  int res=getNewVar();
  _predicateSymbols.insert(res, pred);
  return res;
}

/**
 * Return a propositional predicate name that can be used to represent
 * BDD variable @b var
 */
string BDD::getPropositionalPredicateName(int var)
{
  CALL("BDD::getPropositionalPredicateName");

  string prefix(BDD_PREDICATE_PREFIX);
  prefix+=env.options->namePrefix();
  string name = prefix + Int::toString(var);

  //We do not want the predicate to be already present!
  //(But we also don't want to insert it into signature,
  //as it would grow too much.)
  ASS(!env.signature->isPredicateName(name, 0));

  return name;
}

/**
 * If BDD variable has a corresponding propositional predicate symbol,
 * return true, and assign its name to @b res. Otherwise return false.
 */
bool BDD::getNiceName(int var, string& res)
{
  CALL("BDD::getNiceName");

  unsigned pred;
  bool found=_predicateSymbols.find(var, pred);
  if(found) {
    res=env.signature->predicateName(pred);
  }
  return found;
}

/**
 * If BDD variable has a corresponding propositional predicate symbol,
 * return pointer to its Signature::Symbol object. Otherwise return 0.
 */
Signature::Symbol* BDD::getSymbol(int var)
{
  CALL("BDD::getSymbol");

  unsigned pred;
  if(_predicateSymbols.find(var, pred)) {
    return env.signature->getPredicate(pred);
  }
  return 0;
}

/**
 * Return BDD representing an atomic
 */
BDDNode* BDD::getAtomic(int varNum, bool positive)
{
  CALL("BDD::getAtomic");
  ASS_G(varNum,0);

  if(varNum>=_newVar) {
    _newVar=varNum+1;
  }

  if(positive) {
    return getNode(varNum, getTrue(), getFalse());
  } else {
    return getNode(varNum, getFalse(), getTrue());
  }
}

/**
 * Return conjunction of @b n1 and @b n2
 */
BDDNode* BDD::conjunction(BDDNode* n1, BDDNode* n2)
{
  CALL("BDD::conjunction");
  return getBinaryFnResult(n1,n2, ConjunctionFn(this));
}

/**
 * Return disjunction of @b n1 and @b n2
 */
BDDNode* BDD::disjunction(BDDNode* n1, BDDNode* n2)
{
  CALL("BDD::disjunction");
  return getBinaryFnResult(n1,n2, DisjunctionFn(this));
}

/**
 * Return result of @b x | ~ @b y
 */
BDDNode* BDD::xOrNonY(BDDNode* x, BDDNode* y)
{
  CALL("BDD::xOrNonY");
  return getBinaryFnResult(x,y, XOrNonYFn(this));
}


/**
 * Return true iff @b x | ~ @b y is a constant formula with truth value
 * equal to @b resValue.
 */
bool BDD::isXOrNonYConstant(BDDNode* x, BDDNode* y, bool resValue)
{
  CALL("BDD::isXOrNonYConstant");
  return hasConstantResult(x,y, resValue, XOrNonYFn(this));
}

/**
 * Return result of applying the binary function specified by the
 * BinBoolFn functor to @b n1 and @b n2. The binary functor BinBoolFn
 * must allow to be called as
 *
 * @code
 * BDDNode* BinBoolFn(BDDNode* m1, BDDNode* m2)
 * @endcode
 *
 * and return either the result of applying the represented binary
 * function to @b m1 and @b m2, or 0 in case the result cannot
 * be determined locally by examining just the BDDNode objects pointed
 * by @b m1 and @b m2. It must not return 0 if both arguments are
 * either true or false formulas.
 */
template<class BinBoolFn>
BDDNode* BDD::getBinaryFnResult(BDDNode* n1, BDDNode* n2, BinBoolFn fn)
{
  CALL("BDD::getBinaryFnResult");
  ASS(n1);
  ASS(n2);

  TimeCounter tc(TC_BDD);

  int counter=0;

  static Stack<BDDNode*> toDo(8);
  //Results stack contains zeroes and proper pointers standing for
  //intermediate results.
  //It can be viewed as a prefix of an expression in prefix notation
  //with 0 being a binary function and non-zeroes constants.
  //The expression is being simplified every time a well formed
  //subexpression (i.e. zero followed by two non-zeroes) appears.
  //
  //For the purpose of caching, each 0 is preceeded by two pointers that
  //will be used as a key in the cache to store the result when it is
  //assembled.
  static Stack<BDDNode*> results(8);
  //Variable numbers to be used for building intermediate results in
  //the results stack.
  static Stack<int> vars(8);

  static DHMap<pair<BDDNode*,BDDNode*>, BDDNode*, PtrPairSimpleHash > cache;
  //if the cache was not reset, too much memory would be consumed
  cache.reset();

  for(;;) {
    counter++;
    if(counter==50000) {
      counter=0;
      //time to check if we aren't over the time limit
      if(env.timeLimitReached()) {
	throw TimeLimitExceededException();
      }
    }
    BDDNode* res=fn(n1,n2);
    if(res || cache.find(make_pair(n1, n2), res)) {
      while(results.isNonEmpty() && results.top()!=0) {
	BDDNode* pos=results.pop();
	BDDNode* neg=res;
	unsigned var=vars.pop();
	if(pos==neg) {
	  res=pos;
	} else {
	  res=getNode(var, pos, neg);
	}
	ASS_EQ(results.top(),0);
	results.pop();
	BDDNode* arg1=results.pop();
	BDDNode* arg2=results.pop();
	if( !(counter%4) ) {
	  cache.insert(make_pair(arg1, arg2), res);
	}
      }
      results.push(res);
    } else {
      //we split at variables with higher numbers first
      int splitVar=max(n1->_var, n2->_var);
      ASS_GE(splitVar,0);
      toDo.push((n2->_var==splitVar) ? n2->_neg : n2);
      toDo.push((n1->_var==splitVar) ? n1->_neg : n1);
      toDo.push((n2->_var==splitVar) ? n2->_pos : n2);
      toDo.push((n1->_var==splitVar) ? n1->_pos : n1);
      results.push(n2);
      results.push(n1);
      results.push(0);
      //now push arguments at the stack, so that we know store the
      //answer into the cache
      vars.push(splitVar);
    }

    if(toDo.isEmpty()) {
      break;
    }
    n1=toDo.pop();
    n2=toDo.pop();
  }

  ASS(toDo.isEmpty());
  ASS_EQ(results.length(),1);
  return results.pop();
}

/**
 * Return true iff the result of applying the binary function specified
 * by the BinBoolFn functor to @b n1 and @b n2 is a constant formula with truth
 * value equal to @b resValue.
 *
 * The binary functor BinBoolFn must allow to be called as
 *
 * @code
 * BDDNode* BinBoolFn(BDDNode* m1, BDDNode* m2)
 * @endcode
 *
 * and return either the result of applying the represented binary
 * function to @b m1 and @b m2, or 0 in case the result cannot
 * be determined locally by examining just the BDDNode objects pointed
 * by @b m1 and @b m2. It must not return 0 if both arguments are
 * either true or false formulas.
 */
template<class BinBoolFn>
bool BDD::hasConstantResult(BDDNode* n1, BDDNode* n2, bool resValue, BinBoolFn fn)
{
  CALL("BDD::hasConstantResult");
  ASS(n1);
  ASS(n2);

  TimeCounter tc(TC_BDD);

  int counter=0;

  static Stack<BDDNode*> toDo(8);
  toDo.reset();

  static DHMap<pair<BDDNode*,BDDNode*>, bool, PtrPairSimpleHash > cache;
  //if the cache was not reset, too much memory would be consumed
  cache.reset();

  for(;;) {
    counter++;
    if(counter==50000) {
      counter=0;
      //time to check if we aren't over the time limit
      if(env.timeLimitReached()) {
	throw TimeLimitExceededException();
      }
    }
    BDDNode* res=fn(n1,n2);
    if(res) {
      if( (resValue && !isTrue(res)) ||
	      (!resValue && !isFalse(res))) {
	return false;
      }
    } else {
      if(!cache.find(make_pair(n1, n2)))
      {
	//we split at variables with higher numbers first
	int splitVar=max(n1->_var, n2->_var);
	ASS_GE(splitVar,0);
	toDo.push((n2->_var==splitVar) ? n2->_neg : n2);
	toDo.push((n1->_var==splitVar) ? n1->_neg : n1);
	toDo.push((n2->_var==splitVar) ? n2->_pos : n2);
	toDo.push((n1->_var==splitVar) ? n1->_pos : n1);

	if( !(counter%4) ) {
	  cache.insert(make_pair(n1, n2), false);
	}
      }
    }

    if(toDo.isEmpty()) {
      break;
    }
    n1=toDo.pop();
    n2=toDo.pop();
  }

  return true;
}

/**
 * If it is possible to locally determine the result of the conjunction
 * of @b n1 and @b n2, return the result, otherwise return 0.
 */
BDDNode* BDD::ConjunctionFn::operator()(BDDNode* n1, BDDNode* n2)
{
  if(_parent->isFalse(n1) || _parent->isFalse(n2)) {
    return _parent->getFalse();
  }
  if(_parent->isTrue(n1)) {
    return n2;
  }
  if(_parent->isTrue(n2)) {
    return n1;
  }
  return 0;
}

/**
 * If it is possible to locally determine the result of the disjunction
 * of @b n1 and @b n2, return the result, otherwise return 0.
 */
BDDNode* BDD::DisjunctionFn::operator()(BDDNode* n1, BDDNode* n2)
{
  if(_parent->isTrue(n1) || _parent->isTrue(n2)) {
    return _parent->getTrue();
  }
  if(_parent->isFalse(n1)) {
    return n2;
  }
  if(_parent->isFalse(n2)) {
    return n1;
  }
  return 0;
}

/**
 * If it is possible to locally determine the result of the operation
 * @b n1 | ~ @b n2, return the result, otherwise return 0.
 */
BDDNode* BDD::XOrNonYFn::operator()(BDDNode* n1, BDDNode* n2)
{
  if(_parent->isTrue(n1) || _parent->isFalse(n2)) {
    return _parent->getTrue();
  }
  if(_parent->isTrue(n2)) {
    return n1;
  }
  return 0;
}


/**
 * Return a BDD node containing variable @b varNum that points
 * positively to @b pos and negatively to @b neg BDD node.
 */
BDDNode* BDD::getNode(int varNum, BDDNode* pos, BDDNode* neg)
{
  CALL("BDD::getNode");
  ASS_GE(varNum,0);
  ASS_L(varNum,_newVar);
  ASS_NEQ(pos,neg);

  //newNode contains either 0 or pointer to a BDDNode that
  //hasn't been used yet.
  static BDDNode* newNode=0;

  if(newNode==0) {
    newNode=new BDDNode();
  }

  newNode->_var=varNum;
  newNode->_pos=pos;
  newNode->_neg=neg;

  BDDNode* res=_nodes.insert(newNode);
  if(res==newNode) {
    newNode=0;
  }
  return res;
}


/**
 * Return a string representation of the formula represented by @b node.
 */
string BDD::toString(BDDNode* node)
{
  string res="";
  static Stack<BDDNode*> nodes(8);
  nodes.push(node);
  while(nodes.isNonEmpty()) {
    BDDNode* n=nodes.pop();
    bool canPrintSeparator=true;
    if(n==0) {
      res+=") ";
    } else if(isTrue(n)) {
      res+="$true ";
    } else if(isFalse(n)) {
      res+="$false ";
    } else {
      res+=string("( ")+Int::toString(n->_var)+" ? ";
      nodes.push(0);
      nodes.push(n->_neg);
      nodes.push(n->_pos);
      canPrintSeparator=false;
    }
    if(canPrintSeparator && nodes.isNonEmpty() && nodes.top()) {
      res+=": ";
    }
  }
  return res;
}

/**
 * Return the formula represented by @b node in a TPTP compatible format.
 * The @b bddPrefix string will be added as a prefix to the each BDD
 * variable number to form a predicate symbol name.
 *
 * @warning A recursion is used in this methos, which can lead to
 *   problems with very large BDDs.
 */
string BDD::toTPTPString(BDDNode* node, string bddPrefix)
{
  if(isTrue(node)) {
    return "$true";
  } else if(isFalse(node)) {
    return "$false";
  } else {
    return string("( ( ")+bddPrefix+Int::toString(node->_var)+" => "+toTPTPString(node->_pos, bddPrefix)+
      ") & ( ~"+bddPrefix+Int::toString(node->_var)+" => "+toTPTPString(node->_neg, bddPrefix)+" ) )";
  }
}

/**
 * Return the formula represented by @b node in a TPTP compatible format.
 *
 * @warning A recursion is used in this methos, which can lead to
 *   problems with very large BDDs.
 */
string BDD::toTPTPString(BDDNode* node)
{
  if(isTrue(node)) {
    return "$true";
  } else if(isFalse(node)) {
    return "$false";
  } else {
    return string("( ( ")+getPropositionalPredicateName(node->_var)+" => "+toTPTPString(node->_pos)+
      ") & ( ~"+getPropositionalPredicateName(node->_var)+" => "+toTPTPString(node->_neg)+" ) )";
  }
}


/**
 * Check whether two BDDNode objects are equal
 */
bool BDD::equals(const BDDNode* n1,const BDDNode* n2)
{
  return n1->_var==n2->_var && n1->_pos==n2->_pos && n1->_neg==n2->_neg;
}
/**
 * Return hash value of a BDDNode object
 */
unsigned BDD::hash(const BDDNode* n)
{
  CALL("BDD::hash");

  unsigned res=Hash::hash(n->_var);
  res=Hash::hash(n->_pos, res);
  res=Hash::hash(n->_neg, res);
  return res;
}

/**
 * Convert a BDDNode into a list of propositional clauses.
 */
SATClauseList* BDD::toCNF(BDDNode* node)
{
  CALL("BDD::toCNF");

  SATClauseList* res=0;
  static Stack<pair<BDDNode*,bool> > stack;
  stack.reset();

  for(;;) {
    while(!isConstant(node)) {
      stack.push(make_pair(node, false));
      node=node->_neg;
    }
    if(isFalse(node)) {
      unsigned clen=stack.size();
      SATClause* cl=new(clen) SATClause(clen, true);

      for(unsigned i=0;i<clen;i++) {
	(*cl)[i]=SATLiteral(stack[i].first->_var, !stack[i].second);
      }

      SATClauseList::push(cl, res);
    }
    while(stack.isNonEmpty() && stack.top().second==true) {
      stack.pop();
    }
    if(stack.isEmpty()) {
      return res;
    }
    ASS(!stack.top().second);
    stack.top().second=true;
    node=stack.top().first->_pos;
  }
}

/**
 * Add the formula represented by @b n to the conjunction represented
 * by this object
 */
void BDDConjunction::addNode(BDDNode* n)
{
  CALL("BDDConjunction::addNode");

  if(_isFalse) {
    return;
  }

  BDD* bdd=BDD::instance();

  if(bdd->isConstant(n)) {
    if(bdd->isFalse(n)) {
      _isFalse=true;
    } else {
      ASS(bdd->isTrue(n));
    }
    return;
  }

  TimeCounter tc(TC_SAT_SOLVER);

  if(n->_var > _maxVar) {
    _maxVar=n->_var;
  }


  unsigned varCnt=_maxVar+1;

  SATClauseList* newClLst=bdd->toCNF(n);

  _solver.ensureVarCnt(varCnt);
  _solver.addClauses(pvi( SATClauseList::DestructiveIterator(newClLst) ));

  if(_solver.getStatus()==TWLSolver::UNSATISFIABLE) {
    _isFalse=true;
  }

  return;
}


}
