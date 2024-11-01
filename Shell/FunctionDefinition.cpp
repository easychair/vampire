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
 * @file FunctionDefinition.cpp
 * Implements class FunctionDefinition implementing work with definitions.
 *
 * @since 28/05/2004 Manchester
 */


#include "Lib/Allocator.hpp"
#include "Lib/DHMultiset.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/ScopedLet.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/Formula.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/FormulaUnit.hpp"
#include "Kernel/Problem.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/TermIterators.hpp"

#include "Shell/Options.hpp"

#include "Statistics.hpp"

#include "FunctionDefinition.hpp"

#include <algorithm>
#if VDEBUG
#include <iostream>
#endif

namespace Shell {

using namespace std;
using namespace Lib;
using namespace Kernel;

/**
 * Auxiliary structure to represent information about function definitions.
 * @since 29/05/2004 Manchester
 */
struct FunctionDefinition::Def
{
  enum Mark {
    UNTOUCHED,
    SAFE,
    LOOP,
    BLOCKED,
    UNFOLDED,
    REMOVED
  };
  /** Unit containing a definition */
  Clause* defCl;
  /** The defined function symbol number */

  int fun;
  /** The lhs of the definition */
  Term* lhs;
  /** The rhs of the definition */
  Term* rhs;

  Mark mark;
  /** is the definition linear? */
  bool linear;
  /** strict means that all lhs variables occur in rhs */
  bool strict;

  bool twoConstDef;

  /** first defined function that is used in @b rhs, or -1 if there isn't such */
  int containedFn;

  int examinedArg;

  IntList* dependentFns;

  bool lhsIsBool(){
    return (env.signature->isFoolConstantSymbol(true , fun) ||
            env.signature->isFoolConstantSymbol(false, fun));
  }

  /**
   * If @b mark==SAFE or @b mark==UNFOLDED, contains @b bool array such that
   * @b argOccurs[i] is true iff i-th argument of @b lhs occurs in @b rhs after
   * all definition unfolding is performed.
   */
  bool* argOccurs;

  Def(Term* l, Term* r, bool lin, bool str)
    : fun(l->functor()),
      lhs(l),
      rhs(r),
      mark(UNTOUCHED),
      linear(lin),
      strict(str),
      twoConstDef(0),
      containedFn(-1),
      dependentFns(0),
      argOccurs(0)
  {
  }

  ~Def()
  {
    if(argOccurs) {
      DEALLOC_KNOWN(argOccurs, lhs->arity()*sizeof(bool), "FunctionDefinition::Def::argOccurs");
    }
  }
}; // class FunctionDefintion::Def


/**
 * Initialise a function definition.
 * Move all definitions from the problem to _definitions and _map.
 * @since 29/05/2004 Manchester
 */
FunctionDefinition::FunctionDefinition ()
  :
//    _defStore(32),
    _found(0),
    _removed(0),
    _processedProblem(0)
{
} // FunctionDefinition::FunctionDefinition

void FunctionDefinition::removeUnusedDefinitions(Problem& prb, bool inHigherOrder)
{
  if(removeUnusedDefinitions(prb.units(), &prb, inHigherOrder)) {
    prb.invalidateByRemoval();
  }
}

/**
 * From @b units remove clauses that contain function definitions
 * which are never used (i.e. the function occurs only in it's definition)
 * Return true iff set of clauses in @c units is modified. If the argument @c prb is
 * non-zero, information about removed function definitions is stored in it.
 *
 * The removal is performed iteratively, so that if one function is used
 * only in definition of an other one which is removed, the first definition
 * is removed as well.
 */
bool FunctionDefinition::removeUnusedDefinitions(UnitList*& units, Problem* prb, bool inHigherOrder)
{
  unsigned funs=env.signature->functions();

  Stack<Def*> defStack;
  DArray<Def*> def;
  DArray<unsigned> occCounter;
  def.init(funs, 0);
  occCounter.init(funs, 0);

  UnitList::DelIterator scanIterator(units);
  while(scanIterator.hasNext()) {
    Clause* cl=static_cast<Clause*>(scanIterator.next());
    unsigned clen=cl->length();
    ASS(cl->isClause());
    Def* d=isFunctionDefinition(cl,inHigherOrder);
    if(d) {
      d->defCl=cl;
      if(!def[d->fun]) {
        defStack.push(d);
        def[d->fun]=d;
        scanIterator.del();
      } else {
        delete d;
      }
    }
    for(unsigned i=0;i<clen;i++) {
      NonVariableNonTypeIterator nvit((*cl)[i]);
      while(nvit.hasNext()) {
        unsigned fn=nvit.next()->functor();
        occCounter[fn]++;
      }
    }
  }

  Stack<Def*> toDo;
  Stack<Def*>::Iterator dit(defStack);
  while(dit.hasNext()) {
    Def* d=dit.next();
    unsigned fn=d->fun;
    ASS_GE(occCounter[fn],1);
    if(occCounter[fn]==1) {
      toDo.push(d);
    }
  }

  while(toDo.isNonEmpty()) {
    Def* d=toDo.pop();
    d->mark=Def::REMOVED;
    ASS_EQ(d->defCl->length(), 1);
    ASS_EQ(occCounter[d->fun], 1);
    NonVariableNonTypeIterator nvit((*d->defCl)[0]);
    while(nvit.hasNext()) {
      unsigned fn=nvit.next()->functor();
      occCounter[fn]--;
      if(occCounter[fn]==1 && def[fn]) {
	toDo.push(def[fn]);
      }
    }
    ASS_EQ(occCounter[d->fun], 0);
  }

  bool modified = false;
  while(defStack.isNonEmpty()) {
    Def* d=defStack.pop();
    if(d->mark==Def::REMOVED) {
      modified = true;
      if(prb) {
	ASS_EQ(d->defCl->length(), 1);
	prb->addEliminatedFunction(d->fun, (*d->defCl)[0]);
      }
    }
    else {
      ASS_EQ(d->mark, Def::UNTOUCHED);
      UnitList::push(d->defCl, units);
    }
    delete d;
  }
  return modified;
}

void FunctionDefinition::removeAllDefinitions(Problem& prb, bool inHigherOrder)
{
  ScopedLet<Problem*> prbLet(_processedProblem, &prb);
  if(removeAllDefinitions(prb.units(),inHigherOrder)) {
    prb.invalidateByRemoval();
  }
}

void FunctionDefinition::reverse(Def* def){
  ASS(def->twoConstDef);
  Term* temp = def->lhs;
  def->lhs = def->rhs;
  def->rhs = temp;
  def->fun = def->lhs->functor();
}

/**
 * When possible, unfold function definitions in @b units and remove them
 * Return true iff the list of units was modified.
 */
bool FunctionDefinition::removeAllDefinitions(UnitList*& units, bool inHigherOrder)
{
  UnitList::DelIterator scanIterator(units);
  while(scanIterator.hasNext()) {
    Clause* cl=static_cast<Clause*>(scanIterator.next());
    ASS(cl->isClause());
    Def* d=isFunctionDefinition(cl,inHigherOrder);
    if(d) {
      d->defCl=cl;
      bool inserted = false;
      if(_defs.insert(d->fun, d)) {
        //cout<<"Found: "<<(*(*d->defCl)[0])<<endl;
        inserted = true;
        scanIterator.del();
      } else if(_defs.get(d->fun)->twoConstDef){
        Def* d2;
        ALWAYS(_defs.pop(d->fun, d2));
        reverse(d2);
        if(!d2->lhsIsBool() && _defs.insert(d2->fun, d2)){
          _defs.insert(d->fun, d);
          inserted = true;
          scanIterator.del();
        } else {
          reverse(d2); //back to original orientation
          ALWAYS(_defs.insert(d2->fun, d2));
        }
      } else if(d->twoConstDef){
         reverse(d);
         if(!d->lhsIsBool() && _defs.insert(d->fun, d)) {
           inserted = true;
           scanIterator.del();
         }
      }
      if(!inserted){
        delete d;
      }
    }
  }

  if(!_defs.size()) {
    return false;
  }

  Fn2DefMap::Iterator dit(_defs);
  while(dit.hasNext()) {
    Def* d=dit.next();

    if(d->mark==Def::SAFE || d->mark==Def::BLOCKED) {
      continue;
    }
    ASS(d->mark==Def::UNTOUCHED);
    checkDefinitions(d);
  }

  while(_blockedDefs.isNonEmpty()) {
    Def* d=_blockedDefs.pop();
    ASS_EQ(d->mark, Def::BLOCKED);
//    cout<<"Blocked: "<<(*(*d->defCl)[0])<<endl;

    UnitList::push(d->defCl, units);
    _defs.remove(d->fun);
    delete d;
  }

  if(_defs.isEmpty()) {
    return false;
  }

  ASS_EQ(_defs.size(), _safeDefs.size());
  //_safeDefs contains definitions in topologically ordered,
  //so that _safeDefs[i] uses only definitions up to
  //_safeDefs[i-1].
  for(unsigned i=0;i<_safeDefs.size(); i++) {
    Def* d=_safeDefs[i];
    ASS_EQ(d->mark, Def::SAFE);
//    cout<<"Safe: "<<(*(*d->defCl)[0]);

    //we temporarily block the definition, so that we can rewrite
    //the definition clause without rewriting the lhs
    d->mark=Def::BLOCKED;
//    Clause* oldCl=d->defCl;
    d->defCl=applyDefinitions(d->defCl);
//    cout<<" unfolded into "<<(*(*d->defCl)[0])<<endl;

    //update d->rhs with the right hand side of the equality
    Literal* defEq=(*d->defCl)[0];
    if( defEq->nthArgument(0)->term()==d->lhs ) {
      d->rhs=defEq->nthArgument(1)->term();
    } else {
      ASS_EQ(defEq->nthArgument(1)->term(),d->lhs);
      d->rhs=defEq->nthArgument(0)->term();
    }

    d->mark=Def::UNFOLDED;
    if(_processedProblem) {
      ASS_EQ(d->defCl->length(), 1);
      _processedProblem->addEliminatedFunction(d->fun, (*d->defCl)[0]);
    }

    if (env.options->showPreprocessing()) {
      std::cout << "[PP] fn def discovered: "<<(*d->defCl)<<"\n  unfolded: "<<(*d->rhs) << std::endl;
    }
    env.statistics->functionDefinitions++;
  }

  UnitList::DelIterator unfoldIterator(units);
  while(unfoldIterator.hasNext()) {
    Clause* cl=static_cast<Clause*>(unfoldIterator.next());
    ASS(cl->isClause());
    if(cl->isProxyAxiomsDescendant()){ continue; }
    Clause* newCl=applyDefinitions(cl);
    if(cl!=newCl) {
//      cout<<"D- "<<(*cl)<<endl;
//      cout<<"D+ "<<(*newCl)<<endl;
      unfoldIterator.replace(newCl);
    }
  }

  _safeDefs.reset();
  return true;
}

void FunctionDefinition::checkDefinitions(Def* def0)
{
  TermList t=TermList(def0->lhs);

  //Next argument of the current-level term to be processed.
  //An empty term means we've processed all arguments of the term.
  static Stack<TermList*> stack(4);
  //Definition whose rhs is the current-level term (or zero if none).
  static Stack<Def*> defCheckingStack(4);
  //Definition whose lhs is the current-level term (or zero if none).
  static Stack<Def*> defArgStack(4);
  //Current-level term.
  static Stack<Term*> termArgStack(4);

  //loop invariant: the four above stacks contain the same number of elements
  for(;;) {
    Def* d;
    if(t.isEmpty()) {
      d=defCheckingStack.pop();
      defArgStack.pop();
      termArgStack.pop();
      ASS(!d || d->mark==Def::LOOP);
      if(d && d->mark==Def::LOOP) {
        //the definition is safe (i.e. doesn't contain cycle of non-blocked definitions)
        assignArgOccursData(d);
        _safeDefs.push(d);
        d->mark=Def::SAFE;
      }
    } else if(t.isTerm() && !t.term()->isSort()) {
      //the sort check above is defensive programming
      Term* trm=t.term();
      Def* checkedDef=0;
    toplevel_def:
      if(!_defs.find(trm->functor(), d) || d->mark==Def::BLOCKED) {
        d=0;
      }
      if(trm->numTermArguments() > 0 || checkedDef) {
        stack.push(trm->termArgs());
        defCheckingStack.push(checkedDef);
        defArgStack.push(d);
        termArgStack.push(trm);
      }
      if(d) {
        if(d->mark==Def::UNTOUCHED) {
          //enter the definition
          d->mark=Def::LOOP;
          trm=d->rhs;
          checkedDef=d;
          goto toplevel_def;
        } else if(d->mark==Def::LOOP) {
          //unroll stacks until the point when the current
          //definition was entered
          do{
            stack.pop();

            defArgStack.pop();
            termArgStack.pop();
            d=defCheckingStack.pop();
          } while(!d);
          ASS_EQ(d->mark, Def::LOOP);
          d->mark=Def::BLOCKED;
          defArgStack.setTop(0);
          _blockedDefs.push(d);
        } else {
          ASS_EQ(d->mark, Def::SAFE);
        }
      }      
    }
    if(stack.isEmpty()) {
      break;
    }
    TermList* ts=stack.pop();
    if(ts->isNonEmpty()) {
      Def* argDef=defArgStack.top();
      if(argDef) {
        ASS_EQ(argDef->mark,Def::SAFE);
        Term* parentTerm=termArgStack.top();
        while(ts->isNonEmpty() && !argDef->argOccurs[parentTerm->getArgumentIndex(ts)]) {
          ts=ts->next();
        }
        if(ts->isNonEmpty()) {
          stack.push(ts->next());
        }
      } else {
        stack.push(ts->next());
      }
    }
    t=*ts;
  }
  ASS(defCheckingStack.isEmpty());
  ASS(defArgStack.isEmpty());
}

/**
 * Assign a @b bool array into @b updDef->argOccurs, so that
 * @b updDef->argOccurs[i] is true iff i-th argument of the
 * defined function appears also on the rhs of the definition.
 */
void FunctionDefinition::assignArgOccursData(Def* updDef)
{
  ASS(!updDef->argOccurs);

  if(!updDef->lhs->arity()) {
    return;
  }

  updDef->argOccurs=reinterpret_cast<bool*>(ALLOC_KNOWN(updDef->lhs->arity()*sizeof(bool),
	    "FunctionDefinition::Def::argOccurs"));
  std::memset(updDef->argOccurs, 0, updDef->lhs->arity() * sizeof(bool));

  static DHMap<unsigned, unsigned, IdentityHash, DefaultHash> var2argIndex;
  var2argIndex.reset();
  int argIndex=0;
  for (TermList* ts = updDef->lhs->args(); ts->isNonEmpty(); ts=ts->next()) {
    int w = ts->var();
    var2argIndex.insert(w, argIndex);
    argIndex++;
  }

  TermList t=TermList(updDef->rhs);
  static Stack<TermList*> stack(4);
  static Stack<Def*> defArgStack(4);
  static Stack<Term*> termArgStack(4);
  for(;;) {
    Def* d;
    if(t.isEmpty()) {
      defArgStack.pop();
      termArgStack.pop();
    } else if(t.isTerm()) {
      Term* trm=t.term();
      if(trm->arity()) {
        if(trm->isSort() || !_defs.find(trm->functor(), d) || d->mark==Def::BLOCKED) {
          d=0;
        }
        ASS(!d || d->mark==Def::SAFE);
        stack.push(trm->args());
        defArgStack.push(d);
        termArgStack.push(trm);
      }
    } else {
      ASS(t.isOrdinaryVar());
      updDef->argOccurs[var2argIndex.get(t.var())]=true;
    }
    if(stack.isEmpty()) {
      break;
    }
    TermList* ts=stack.pop();
    if(!ts->isEmpty()) {
      Def* argDef=defArgStack.top();
      if(argDef) {
        Term* parentTerm=termArgStack.top();
        while(ts->isNonEmpty() && !argDef->argOccurs[parentTerm->getArgumentIndex(ts)]) {
          ts=ts->next();
        }
        if(ts->isNonEmpty()) {
          stack.push(ts->next());
        }
      } else {
        stack.push(ts->next());
      }
    }
    t=*ts;
  }
}



typedef pair<unsigned,unsigned> BindingSpec;
typedef DHMap<BindingSpec, TermList> BindingMap;
typedef DHMap<BindingSpec, bool> UnfoldedSet;

Term* FunctionDefinition::applyDefinitions(Literal* lit, Stack<Def*>* usedDefs)
{
  //cout << "applying definitions to " + lit->toString() << endl;

  if (env.options->showPreprocessing()) {
    std::cout << "[PP] applying function definitions to literal "<<(*lit) << std::endl;
  }
  BindingMap bindings;
  UnfoldedSet unfolded;
  unsigned nextDefIndex=1;
  Stack<BindingSpec> bindingsBeingUnfolded;
  Stack<unsigned> defIndexes;
  //Terms that to be unfolded. When zero element appears on the
  //stack, a task from @b tasks stack should be finished.
  Stack<TermList*> toDo;
  Stack<Term*> terms;
  Stack<bool> modified;
  Stack<TermList> args;

  bindings.reset();
  unfolded.reset();

  defIndexes.reset();
  toDo.reset();
  terms.reset();
  modified.reset();
  args.reset();

  defIndexes.push(0);
  modified.push(false);
  toDo.push(lit->args());

  for(;;) {
    TermList* tt=toDo.pop();
    if(!tt) {
      BindingSpec spec=bindingsBeingUnfolded.pop();
      bindings.set(spec, args.top());
      unfolded.insert(spec, true);
      continue;
    }
    if(tt->isEmpty()) {
      if(terms.isEmpty()) {
        //we're done, args stack contains modified arguments
        //of the argument term.
        ASS(toDo.isEmpty());
        break;
      }
      defIndexes.pop();
      Term* orig=terms.pop();
      if(!modified.pop()) {
        args.truncate(args.length() - orig->arity());
        args.push(TermList(orig));
        continue;
      }
      //here we assume, that stack is an array with
      //second topmost element at &top()-1, third at
      //&top()-2, etc...
      TermList* argLst=&args.top() - (orig->arity()-1);

      Term* newTrm;
      if(orig->isSort()){
        newTrm=AtomicSort::create(static_cast<AtomicSort*>(orig), argLst);
      } else {
        newTrm=Term::create(orig,argLst);
      }
      args.truncate(args.length() - orig->arity());
      args.push(TermList(newTrm));

      modified.setTop(true);
      continue;
    }
    toDo.push(tt->next());

    TermList tl=*tt;
    unsigned defIndex=defIndexes.top();
    Term* t;

    if(tl.isVar()) {
      ASS(tl.isOrdinaryVar());

      if(defIndexes.top()) {
        modified.setTop(true);
        BindingSpec spec=make_pair(defIndexes.top(), tl.var());
        TermList bound=bindings.get(spec);
        if(bound.isVar() || unfolded.find(spec)) {
          args.push(bound);
          continue;
        } else {
          bindingsBeingUnfolded.push(spec);
          toDo.push(0);

          defIndex=0;
          t=bound.term();
        }
      } else {
        args.push(tl);
        continue;
      }

    } else {
      ASS(tl.isTerm());
      t=tl.term();
    }

    Def* d;
    //sorts can never contain definitions
    if(!t->isSort() && !defIndex && _defs.find(t->functor(), d) && d->mark!=Def::BLOCKED) {
      ASS_EQ(d->mark, Def::UNFOLDED);
      usedDefs->push(d);
      if (env.options->showPreprocessing()) {
        std::cout << "[PP] definition of "<<(*t)<<"\n  expanded to "<<(*d->rhs) << std::endl;
      }

      defIndex=nextDefIndex++;

      //bind arguments of definition lhs
      TermList* dargs=d->lhs->args();
      TermList* targs=t->args();
      while(dargs->isNonEmpty()) {
        ASS(targs->isNonEmpty());
        ALWAYS(bindings.insert(make_pair(defIndex, dargs->var()), *targs));
        dargs=dargs->next();
        targs=targs->next();
      }

      //use rhs of the definition instead of the original subterm
      t=d->rhs;
      modified.setTop(true);
    }
    defIndexes.push(defIndex);
    terms.push(t);
    modified.push(false);
    toDo.push(t->args());
  }
  ASS(toDo.isEmpty());
  ASS(terms.isEmpty());
  ASS_EQ(modified.length(),1);
  ASS_EQ(args.length(),lit->arity());

  if(!modified.pop()) {
    return lit;
  }

  //here we assume, that stack is an array with
  //second topmost element as &top()-1, third at
  //&top()-2, etc...
  TermList* argLst=&args.top() - (lit->arity()-1);
  return Literal::create(static_cast<Literal*>(lit),argLst);
}


Clause* FunctionDefinition::applyDefinitions(Clause* cl)
{
  unsigned clen=cl->length();

  static Stack<Def*> usedDefs(8);
  RStack<Literal*> resLits;
  ASS(usedDefs.isEmpty());

  bool modified=false;
  for(unsigned i=0;i<clen;i++) {
    Literal* lit=(*cl)[i];
    Literal* rlit=static_cast<Literal*>(applyDefinitions(lit, &usedDefs));
    resLits->push(rlit);
    modified|= rlit!=lit;
  }
  if(!modified) {
    ASS(usedDefs.isEmpty());
    return cl;
  }

  UnitList* premises=0;
  std::vector<Term *> extra;
  while(usedDefs.isNonEmpty()) {
    Def *def = usedDefs.pop();
    Clause* defCl=def->defCl;
    UnitList::push(defCl, premises);
    if(env.options->proofExtra() == Options::ProofExtra::FULL)
      extra.push_back(def->lhs);
  }
  std::reverse(extra.begin(), extra.end());
  UnitList::push(cl, premises);
  auto res = Clause::fromStack(*resLits, NonspecificInferenceMany(InferenceRule::DEFINITION_UNFOLDING, premises));
  if(env.options->proofExtra() == Options::ProofExtra::FULL)
    env.proofExtra.insert(res, new FunctionDefinitionExtra(std::move(extra)));
  res->setAge(cl->age()); // TODO isn't this dones automatically?
  return res;
}


/**
 * Destroy by deallocating all Defs.
 * @since 26/05/2007 Manchester
 */
FunctionDefinition::~FunctionDefinition ()
{
  Fn2DefMap::Iterator dit(_defs);
  while(dit.hasNext()) {
    delete dit.next();
  }
} // FunctionDefinition::~FunctionDefinition

/**
 * If the the unit is a function definition f(x1,...,xn) = t,
 * return the Def structure representing information about the definition.
 * @since 26/05/2007 Manchester
 */
FunctionDefinition::Def*
FunctionDefinition::isFunctionDefinition (Unit& unit, bool inHigherOrder)
{
  if(unit.derivedFromGoal() && env.options->ignoreConjectureInPreprocessing()){
    return 0;
  }

  if (unit.isClause()) {
    return isFunctionDefinition(static_cast<Clause*>(&unit), inHigherOrder);
  }
  return isFunctionDefinition(static_cast<FormulaUnit&>(unit), inHigherOrder);
} // Definition::isFunctionDefinition (const Clause& c)


/**
 * If the the clause if a function definition f(x1,...,xn) = t,
 * return the Def structure representing information about the definition.
 * @since 05/06/2004 Manchester
 * @since 26/05/2007 Manchester modified for new data structures
 */
FunctionDefinition::Def*
FunctionDefinition::isFunctionDefinition (Clause* clause, bool inHigherOrder)
{
  if (clause->length() != 1) {
    return 0;
  }
  return isFunctionDefinition((*clause)[0],inHigherOrder);
} // Definition::isFunctionDefinition (Clause* c)

/**
 * If the literal if a function definition f(x1,...,xn) = t.
 * return the Def structure representing information about the definition.
 */
FunctionDefinition::Def*
FunctionDefinition::isFunctionDefinition (Literal* lit, bool inHigherOrder)
{
  if (! lit->isPositive() ||
      ! lit->isEquality() ||
      ! lit->shared()) {
    return 0;
  }

  // the atom is an equality
  TermList* args = lit->args();
  if (args->isVar()) {
    return 0;
  }
  Term* l = args->term();
  args = args->next();
  if (args->isVar()) {
    return 0;
  }
  Term* r = args->term();
  Def* def = defines(l,r,inHigherOrder);
  if (def) {
    return def;
  }
  def = defines(r,l,inHigherOrder);
  if (def) {
    return def;
  }

  return 0;
} // FunctionDefinition::isFunctionDefinition (const Atom&...)


/**
 * If lhs=rhs if a function definition f(x1,...,xn) = t.
 * return the Def structure representing information about the definition.
  *
  * @since 09/05/2002, Manchester
  * @since 24/10/2002 Manchester, bug fixed
  * @since 05/01/2004, Manchester, completely reimplemented
  * @since 29/05/2004, Manchester, changed to check that the
  *   variables of t is a subset of {x1,...,xn} (before it was
  *   a check for equality.
  * @since 26/05/2007, Manchester, reimplemented for new datastructures
  *   and check for equality added
  */
FunctionDefinition::Def*
FunctionDefinition::defines (Term* lhs, Term* rhs, bool inHigherOrder)
{
  if(!lhs->shared() || !rhs->shared()) {
    return 0;
  }
  unsigned f = lhs->functor();
  if(env.signature->getFunction(f)->protectedSymbol()) {
    return 0;
  }
  if(env.signature->getFunction(f)->distinctGroups()!=0) {
    return 0;
  }
  if(lhs->color()==COLOR_TRANSPARENT && rhs->color()!=COLOR_TRANSPARENT) {
    return 0;
  }

  if (occurs(f,*rhs)) {
    return 0;
  }
  if (!lhs->arity()) {
    if(env.signature->isFoolConstantSymbol(true , f) ||
       env.signature->isFoolConstantSymbol(false, f)){
      return 0;
    }
    //Higher-order often contains definitions of the form f = ^x^y...
    if (rhs->arity() && !inHigherOrder) { // c = f(...)
      return 0;
    }
    if (rhs->functor() == f) {
      return 0;
    }
    if(!inHigherOrder){
      return new Def(lhs,rhs,true,true);
    }
  }

  int vars = 0; // counter of variables occurring in the lhs

  // First, iterate subterms in lhs and check that all of them are variables
  // and each of them occurs exactly once. counter will contain variables
  // occurring in lhs
  MultiCounter counter;
  for (const TermList* ts = lhs->args(); ts->isNonEmpty(); ts=ts->next()) {
    if (! ts->isVar()) {
      return 0;
    }
    int w = ts->var();
    if (counter.get(w) != 0) { // more than one occurrence
      return 0;
    }
    counter.inc(w);
    vars++;
  }

  bool linear = true;
  // now check that rhs contains only variables in the counter
  // Iterate over variables in rhs and check that they
  // are counted as 1 or 2.
  // found will be increased by the number of different variables
  // marked 1.
  TermVarIterator vs(rhs->args());
  while (vs.hasNext()) {
    int v = vs.next();
    switch (counter.get(v)) {
    case 0: // v does not occur in lhs
      return 0;

    case 1: // v occurs in lhs, it is first occurrence in rhs
      counter.inc(v);
      vars--;
      break;

    default: // v occurs in lhs, it is not its first occurrence in rhs
      linear = false;
      break;
    }
  }

  Def* res=new Def(lhs,rhs,linear,!vars);

  if(!lhs->arity() && !rhs->arity()){
    res->twoConstDef = true;
  }
  
  return res;
} // FunctionDefinition::defines


/**
 * True if function symbol f occurs in the term.
 * @since 28/09/2002 Manchester, changed to use numeric terms
 * @since 06/01/2004 Manchester, changed to use iterator
 */
bool FunctionDefinition::occurs (unsigned f, Term& t)
{
  TermFunIterator funs(&t);
  while (funs.hasNext()) {
    if (f == funs.next()) {
      return true;
    }
  }
  return false;
} // FunctionDefinition::occurs (const FunctionSymbol& f, const Term& t)


/**
 * If lhs=rhs if a function definition f(x1,...,xn) = t.
 * return the Def structure representing information about the definition.
 *
 * @since 26/07/2002 Bolzano, changed
 * @since 26/05/2007 Manchester, reimplemented using new datastructures
 */
FunctionDefinition::Def*
FunctionDefinition::isFunctionDefinition (FormulaUnit& unit, bool inHigherOrder)
{
  Formula* f = unit.formula();
  // skip all universal quantifiers in front of the formula
  while (f->connective() == FORALL) {
    f = f->qarg();
  }

  if (f->connective() != LITERAL) {
    return 0;
  }
  return isFunctionDefinition(f->literal(), inHigherOrder);
} // FunctionDefinition::isFunctionDefinition

/**
 * Delete the definition def.
 * @since 27/05/2007 flight Frankfurt-Lisbon
 */
void FunctionDefinition::deleteDef (Def* def)
{
  delete def;
} // FunctionDefinition::deleteDef

void FunctionDefinitionExtra::output(std::ostream &out) const {
  bool first = true;
  out << "inlined=[";
  for(Term *t : lhs) {
    if(!first)
      out << ",";
    first = false;
    out << t->toString();
  }
  out << "]";
}

}
