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
 * @file ModelCheck.hpp
 * Defines checking of finite models
 *
 * @since 6/10/2015 Manchester
 * @author Giles
 */

#ifndef __ModelCheck__
#define __ModelCheck__

#include "Lib/DHMap.hpp"
#include "Lib/Set.hpp"
#include "Lib/Environment.hpp"
#include "Lib/StringUtils.hpp"

#include "Kernel/Problem.hpp"
#include "Kernel/Unit.hpp"
#include "Kernel/Clause.hpp"
#include "Kernel/Term.hpp"

#include "FiniteModelMultiSorted.hpp"

namespace FMB{

using namespace Lib;
using namespace Kernel;

using std::cout;
using std::endl;

class ModelCheck{

public:
static void doCheck(UnitList* units)
{
  // find model size
  // looking for a domain axiom with a name starting 'finite_domain' (TODO search for something of the right shape)

  DHMap<unsigned,unsigned> sortSizes;
  DHMap<unsigned,std::unique_ptr<Set<Term*>>> domainConstantsPerSort;

  // first just search for finite_domain axiom (TODO: do this for every sort!)
  {
    UnitList::Iterator uit(units);
    while(uit.hasNext()){
      Unit* u = uit.next();
      if(u->inputType()!= UnitInputType::MODEL_DEFINITION) continue;
      std::string name;
      ALWAYS(Parse::TPTP::findAxiomName(u,name));
      if(StringUtils::starts_with(name,"finite_domain")){
        std::cout << "Finite domain axiom found:" << std::endl << u->toString() << std::endl;
        // Set model size and domainConstants
        // And check it is a finite domain axiom
        if(u->isClause()) USER_ERROR("finite_domain in cnf (use fof/tff instead)");

        Formula* formula = u->getFormula();
        if(formula->connective()!=Connective::FORALL) USER_ERROR("finite_domain is not a domain axiom");
        Formula* subformula = formula->qarg();

        if (VList::length(formula->vars()) != 1) USER_ERROR("finite_domain must start with forall over a single variable");
        unsigned curSort = (formula->sorts()) ? formula->sorts()->head().term()->functor() : env.signature->getDefaultSort();

        unsigned curModelSize = 0;
        auto curDomainConstants = std::make_unique<Set<Term*>>();

        int single_var = -1;
        if (subformula->connective()==Connective::OR) {
          FormulaList* args = subformula->args();
          FormulaList::Iterator fit(args);
          while(fit.hasNext()){
            curModelSize++;
            Formula* arg = fit.next();
            if(arg->connective()!=Connective::LITERAL) USER_ERROR("finite_domain is not a domain axiom");
            Literal* l = arg->literal();
            checkIsDomainLiteral(l,single_var,*curDomainConstants);
          }
        } else if (subformula->connective()==Connective::LITERAL) {
          curModelSize = 1;
          checkIsDomainLiteral(subformula->literal(),single_var,*curDomainConstants);
        } else USER_ERROR("finite_domain is not a domain axiom");

        sortSizes.set(curSort,curModelSize);
        domainConstantsPerSort.set(curSort,std::move(curDomainConstants));
      }
    }
  }

  if (sortSizes.size() == 0) {
    USER_ERROR("No domain definition found!\n Use vampire(model_check,model_start). MODEL_DEFINITION_FORMULAS. vampire(model_check,model_end).");
  }

  std::cout << "Loading model..." << std::endl;
  DArray<unsigned> sortSizesArray(env.signature->typeCons());
  {
    auto it = sortSizes.items();
    while (it.hasNext()) {
      auto [sort,curModelSize] = it.next();
      sortSizesArray[sort] = curModelSize;
    }
  }
  // TODO can we pass a reference here instead of clone()ing?
  FiniteModelMultiSorted model(sortSizesArray.clone());

  Set<Term*> domainConstants; // union of all the perSort ones
  DHMap<Term*,unsigned> domainConstantNumber;

  std::cout << "Detected model with " << sortSizes.size() << " sorts." << std::endl;
  auto it = sortSizes.items();
  while (it.hasNext()) {
    auto [sort,curModelSize] = it.next();
    auto& curDomainConstants = *domainConstantsPerSort.get(sort);
    ASS_EQ(curModelSize,(unsigned)curDomainConstants.size());

    std::cout << "  domain for " + env.signature->getTypeCon(sort)->name() + " of size " << curModelSize;
    std::cout << ", domain elements are:" << std::endl;

    // number the domain constants
    Set<Term*>::Iterator dit(curDomainConstants);
    unsigned count=1;
    while(dit.hasNext()){
      Term* con = dit.next();
      std::cout << "    " << con->toString() << std::endl;
      domainConstants.insert(con);
      domainConstantNumber.insert(con,count);
      model.addFunctionDefinition(con->functor(),DArray<unsigned>(0),count);
      count++;
    }
  }

  {
    UnitList::Iterator uit(units);
    while(uit.hasNext()){
      Unit* u = uit.next();
      if(u->inputType()!= UnitInputType::MODEL_DEFINITION) continue;
      std::string name;
      ALWAYS(Parse::TPTP::findAxiomName(u,name));
      if(StringUtils::starts_with(name,"finite_domain") || StringUtils::starts_with(name,"distinct_domain")) continue;

      // All model formulas should be conjunctions of definitions
      // TODO allow for unit clause definitions in the future
      if(u->isClause()) USER_ERROR("Expecting model to use formulas i.e. fof/tff");
      Formula* formula = u->getFormula();

      if(formula->connective()==Connective::NOT){
        Formula* inner = formula->uarg();
        if(inner->connective()==Connective::LITERAL){
          addDefinition(model,inner->literal(),true,domainConstants,domainConstantNumber);
        }
        else USER_ERROR("Unexpected negation in "+formula->toString());
      }
      else if(formula->connective()==Connective::LITERAL){
        addDefinition(model,formula->literal(),false,domainConstants,domainConstantNumber);
      }
      else{
        if(formula->connective()!=Connective::AND)
          USER_ERROR("Expecting conjunction of definitions in model:\n"+formula->toString());
        FormulaList* defs = formula->args();
        FormulaList::Iterator fit(defs);
        while(fit.hasNext()){
          unsigned negated = false;
          Formula* def = fit.next();
          if(def->connective()==Connective::NOT){
            def = def->uarg();
            negated=true;
          }
          if(def->connective()!=Connective::LITERAL) USER_ERROR("Badly formed definition");
          Literal* lit = def->literal();
          addDefinition(model,lit,negated,domainConstants,domainConstantNumber);
        }
      }
    }
    std::cout << "Model loaded" << std::endl;
    std::cout << "Checking formulas..." << std::endl;
    {
      bool allTrue = true;
      UnitList::Iterator uit(units);
      while(uit.hasNext()){
        Unit* u = uit.next();
        if(u->inputType()== UnitInputType::MODEL_DEFINITION) continue;

        std::cout << "Checking " << u->toString() << "..." << std::endl;
        bool res = model.evaluate(u);
        allTrue &= res;
        std::cout << "Evaluates to " << (res ? "True" : "False") << std::endl;
      }
      if (allTrue) {
        std::cout << "All formulas evaluated to True!" << std::endl;
      } else {
        std::cout << "There was a false formula!" << std::endl;
      }
    }
  }
}

private:

static void checkIsDomainLiteral(Literal* l, int& single_var, Set<Term*>& domainConstants)
{
  if(!l->isEquality()) USER_ERROR("finite_domain is not a domain axiom");

  // put var in left and constant in right
  TermList* left = l->nthArgument(0);
  TermList* right = l->nthArgument(1);
  if(right->isVar()){
    TermList* temp = left;left=right;right=temp;
  }
  if(right->isVar()) USER_ERROR("finite_domain is not a domain axiom");

  // store and check the single variable used
  if(single_var<0) single_var=left->var();
  if(left->var()!=(unsigned)single_var) USER_ERROR("finite_domain is not a domain axiom");

  // store and check the ground constant used
  Term* constant = right->term();
  unsigned f = constant->functor();
  if(env.signature->functionArity(f)!=0) USER_ERROR("finite_domain is not a domain axiom");
  if(domainConstants.contains(constant)) USER_ERROR("finite_domain is not a domain axiom");

  domainConstants.insert(constant);
}

static void addDefinition(FiniteModelMultiSorted& model,Literal* lit,bool negated,
                          Set<Term*>& domainConstants,
                          DHMap<Term*,unsigned>& domainConstantNumber)
{
  if(lit->isEquality()){
    if(!lit->polarity() || negated) USER_ERROR("Cannot have negated function definition");
    // Defining a function or constant
    TermList* left = lit->nthArgument(0);
    TermList* right = lit->nthArgument(1);
    if(domainConstants.contains(left->term())){
      TermList* temp = left; left=right;right=temp;
    }

    if(domainConstants.contains(left->term()))
      USER_ERROR("Cannot have equality between domain elements:\n"+lit->toString());
    unsigned res = domainConstantNumber.get(right->term());
    if(left->isVar()) USER_ERROR("Expect term on left of definition");
    Term* fun = left->term();
    unsigned f = fun->functor();
    unsigned arity = env.signature->functionArity(f);
    DArray<unsigned> args(arity);
    for(unsigned i=0;i<arity;i++){
      TermList* arg = fun->nthArgument(i);
      if(arg->isVar() || !domainConstants.contains(arg->term()))
        USER_ERROR("Expect term on left of definition to be grounded with domain constants");
      args[i] = domainConstantNumber.get(arg->term());
    }
    model.addFunctionDefinition(f,args,res);
  }else{
    // not sure this makes sense but...
    if(!lit->polarity()) negated=!negated;
    // Defining a predicate or proposition
    unsigned p = lit->functor();
    unsigned arity = env.signature->predicateArity(p);

    DArray<unsigned> args(arity);
    for(unsigned i=0;i<arity;i++){
      TermList* arg = lit->nthArgument(i);
      if(arg->isVar() || !domainConstants.contains(arg->term()))
        USER_ERROR("Expect term on left of definition to be grounded with domain constants");
      args[i] = domainConstantNumber.get(arg->term());
    }
    model.addPredicateDefinition(p,args,!negated);
  }
}

};

} // namespace FMB
#endif
