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
 * @file TheoryInstAndSimp.cpp
 * Implements class TheoryInstAndSimp.
 */

#if VZ3

#define DPRINT 0

#include "Debug/RuntimeStatistics.hpp"

#include "Lib/Environment.hpp"
#include "Lib/Metaiterators.hpp"
#include "Lib/VirtualIterator.hpp"

#include "Kernel/Clause.hpp"
#include "Kernel/Unit.hpp"
#include "Kernel/Inference.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/Term.hpp"
#include "Kernel/Substitution.hpp"
#include "Kernel/TermIterators.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/OperatorType.hpp"
#include "Kernel/Theory.hpp"
#include "Kernel/FormulaVarIterator.hpp"

#include "Saturation/SaturationAlgorithm.hpp"
#include "Saturation/Splitter.hpp"

#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/TheoryFlattening.hpp"
#include "Shell/UIHelper.hpp"

#include "SAT/SATLiteral.hpp"
#include "SAT/SAT2FO.hpp"
#include "SAT/Z3Interfacing.hpp"

#include "TheoryInstAndSimp.hpp"
#include "Kernel/NumTraits.hpp"
#include "Kernel/TermIterators.hpp"

#define DEBUG(...)  //DBG(__VA_ARGS__)

namespace Inferences
{

using namespace std;
using namespace Lib;
using namespace Kernel;
using namespace Shell;
using namespace Saturation;
using namespace SAT;

using SortId = SAT::Z3Interfacing::SortId;

TheoryInstAndSimp::TheoryInstAndSimp(Options& opts) : TheoryInstAndSimp(
    opts.theoryInstAndSimp(), 
    opts.thiTautologyDeletion(), 
    opts.showZ3(),  
    opts.thiGeneralise(),
    opts.exportThiProblem(),
    opts.problemExportSyntax()
    ) {}


Options::TheoryInstSimp manageDeprecations(Options::TheoryInstSimp mode) 
{
  switch (mode) {
    case Options::TheoryInstSimp::FULL:
    case Options::TheoryInstSimp::NEW:
      if(outputAllowed()) {
        std::cout << "WARNING: the modes full & new are deprecated for theory instantiation. using all instead." << std::endl;
      }
      return Options::TheoryInstSimp::ALL;
    default:
      return mode;
  }
}

TheoryInstAndSimp::TheoryInstAndSimp(Options::TheoryInstSimp mode, bool thiTautologyDeletion, bool showZ3, bool generalisation, std::string const& exportSmtlib, Shell::Options::ProblemExportSyntax problemExportSyntax) 
  : _splitter(0)
  , _mode(manageDeprecations(mode))
  , _thiTautologyDeletion(thiTautologyDeletion)
  , _naming()
  , _solver(new Z3Interfacing(_naming, showZ3, /* unsatCoresForAssumptions = */ generalisation, exportSmtlib, problemExportSyntax))
  , _generalisation(generalisation)
  , _instantiationConstants ("$inst")
  , _generalizationConstants("$inst$gen")
{ }


void TheoryInstAndSimp::attach(SaturationAlgorithm* salg)
{
  SimplifyingGeneratingInference::attach(salg);
  _splitter = salg->getSplitter();
}



bool TheoryInstAndSimp::calcIsSupportedSort(SortId sort)
{
  //TODO: extend for more sorts (arrays, datatypes)
  if(   sort == IntTraits::sort() 
     || sort == RatTraits::sort() 
     || sort == RealTraits::sort() ){
    return true;
  } else if (env.signature->isTermAlgebraSort(sort)) {
    return TermAlgebra::subSorts(sort).iter()
      .all([&](auto s){ return env.signature->isTermAlgebraSort(s) || calcIsSupportedSort(s); });
  } else {
    return false;
  }
}

bool TheoryInstAndSimp::isSupportedSort(SortId sort) 
{ return _supportedSorts.getOrInit(sort, [&](){ return calcIsSupportedSort(sort); }); }

/**
  Wraps around interpretePredicate to support interpreted equality
 */
bool TheoryInstAndSimp::isSupportedLiteral(Literal* lit) {
  //check equality spearately (X=Y needs special handling)
  if (lit->isEquality()) {
    return isSupportedSort(SortHelper::getEqualityArgumentSort(lit));
  }

  //check if predicate is interpreted
  if (! theory->isInterpretedPredicate(lit->functor())){
    return false;
  }

  //check if arguments of predicate are supported
  for (unsigned i=0; i<lit->numTermArguments(); i++) {
    TermList sort = SortHelper::getTermArgSort(lit,i);
    if (! isSupportedSort(sort))
      return false;
  }

  //check if this is an interpreted predicate that is not supported by Z3Interfacing
  switch(theory->interpretPredicate(lit->functor())){
    case Theory::INT_IS_RAT:
    case Theory::INT_IS_REAL:
    case Theory::RAT_IS_RAT:
    case Theory::RAT_IS_REAL:
    case Theory::REAL_IS_RAT:
    case Theory::REAL_IS_REAL:
      return false;
    default:;
  }

  return true;
}

bool TheoryInstAndSimp::isSupportedFunction(Term* trm) {
  auto sym = env.signature->getFunction(trm->functor());
  return !(theory->isInterpretedConstant(trm) 
      || (theory->isInterpretedFunction(trm) && isSupportedFunction(theory->interpretFunction(trm)) )
      || (sym->termAlgebraCons() && isSupportedSort(SortHelper::getResultSort(trm)))
      || (sym->termAlgebraDest() && isSupportedSort(SortHelper::getTermArgSort(trm, 0)))
      );
}


bool TheoryInstAndSimp::isSupportedFunction(Theory::Interpretation itp) {
  switch (itp) {
    case Theory::ARRAY_BOOL_SELECT:
    case Theory::ARRAY_SELECT:
    case Theory::ARRAY_STORE:
    case Theory::INT_SUCCESSOR:
    case Theory::INT_TO_INT:
    case Theory::RAT_TO_RAT:
    case Theory::REAL_TO_RAT:
    case Theory::REAL_TO_REAL:
      return false;
    default: return true;
  }
}



bool TheoryInstAndSimp::isPure(Literal* lit) {
  if (lit->isSpecial()) /* TODO: extend for let .. in / if then else */ {
#if DPRINT
    cout << "special lit " << lit -> toString() << endl;
#endif
    return false;
  }

  //check if the predicate is a theory predicate
  if (! isSupportedLiteral(lit) ) {
    return false;
  }
  //check all (proper) subterms
  NonVariableNonTypeIterator nvi(lit);
  while( nvi.hasNext() ) {
    Term* term = nvi.next();

    //we can stop if we found an uninterpreted function / constant
    if (isSupportedFunction(term)){
      return false;
    }
    //check if return value of term is supported
    if (! isSupportedSort(SortHelper::getResultSort(term))){
      return false;
    }
    //check if arguments of term are supported. covers e.g. f(X) = 0 where
    // f could map uninterpreted sorts to integer. when iterating over X
    // itself, its sort cannot be checked.
    for (unsigned i=0; i<term->numTermArguments(); i++) {
      TermList sort = SortHelper::getTermArgSort(term,i);
      if (! isSupportedSort(sort))
        return false;
    }
  }

#if DPRINT
  cout << "found pure literal: " << lit->toString() << endl;
#endif
  return true;
}

bool TheoryInstAndSimp::isXeqTerm(TermList left, TermList right) {
  bool r = left.isVar() &&
    right.isTerm() &&
    !VList::member(left.var(), freeVariables(right.term()));
  return r;
}

unsigned TheoryInstAndSimp::varOfXeqTerm(const Literal* lit,bool flip) {
  ASS(lit->isEquality());
  ASS(! lit->isPositive());
  //add assertion
  if (lit->isEquality()) {
    TermList left = lit->termArg(0);
    TermList right = lit->termArg(1);
    if (isXeqTerm(left,right)){ return left.var();}
    if (isXeqTerm(right,left)){ return right.var();}
    ASS(lit->isTwoVarEquality());
    if(flip){
      return left.var();
    }else{
      return right.var();
    }
  }
  ASSERTION_VIOLATION ;
  return -1; //TODO: do something proper to prevent compilation warnings
}

/** checks if variable v is contained in literal lit */
bool TheoryInstAndSimp::literalContainsVar(const Literal* lit, unsigned v) {
  SubtermIterator it(lit);
  while (it.hasNext()) {
    const TermList t = it.next();
    if ((t.isVar()) && (t.var() == v)){
      return true;
    }
  }
  return false;
}


/**
 * Scans through a clause C and selects the largest set T s.t. all literals in
 * T are trivial. A literal L is trivial in C if:
 *   1. L is of the form X != s where X does not occur in s
 *   2. L is pure
 *   3. for all literals L' in C that X (different from L) either
 *      + L' is not pure
 *      + L' is trivial in C
 *
 * some observations:
 *   - consider X != Y + 1 | Y != X - 1 | p(X,Y): then {}, as well as {X != Y+1, Y != X-1} are sets of trivial literals
 *   - we can partition the clause into pure and impure literals
 *   - trivial literals are always a subset of the pure literals
 *   - a literal that violates condition is pure and not trivial
 * the algorithm is as follows:
 *   - find the set of trivial candidates TC that fulfill conditions 1 and 2
 *   - define the set of certainly non-trivial pure literals NT as
 *     { X in C | X is pure, X not in TC}
 *   - move all X from TC to NT that do not fulfill criterion 3
 *     (by checking against all elements of NT)
 *   - repeat this step until no element was removed or TC is empty
 * the algorithm can be optimized by only checking the freshly removed elements
 **/
Stack<Literal*> TheoryInstAndSimp::selectTrivialLiterals(Clause* cl)
{
#if DPRINT
  cout << "selecting trivial literals in " << cl->toString() << endl ;
#endif
  /* find trivial candidates of the form x != t (x not occurring in t) */
  auto it = cl->iterLits();
  /* invariants:
       triv_candidates \cup nontriv_pure \cup impure = cl
       triv_candidates \cap nontriv_pure = 0
       triv_candidates \cap impure = 0
       nontriv_pure \cap impure = 0 */
  Stack<Literal*> triv_candidates;
  Stack<Literal*> nontriv_pure;
  Stack<Literal*> impure;
  while( it.hasNext() ) {
    Literal* c = it.next();
    if (isPure(c)) {
      //a liteal X != s is possibly trivial
      if (c->isNegative()
          && c->isEquality()) {
#if DPRINT
        cout << "checking " << c->toString() << endl;
#endif
        TermList left = c->termArg(0);
        TermList right = c->termArg(1);
        /* distinguish between X = s where s not a variable, X = Y and X = X */
        if (TheoryInstAndSimp::isXeqTerm(left, right) ||
            TheoryInstAndSimp::isXeqTerm(right, left) ) {
          triv_candidates.push(c);
        } else {
          // X=Y case
          if( left.isVar()
              && right.isVar()) {
            if (left.var() != right.var()) {
              triv_candidates.push(c);
            } else {
              //this is required by the definition, but making X=X trivial would
              //make more sense
              nontriv_pure.push(c);
            }
          } else { //term = term case
            nontriv_pure.push(c);
          }
        }
      } else {
        //mark as nontrivial pure
#if DPRINT
        cout << "non trivial pure found " << c->toString() << endl;
#endif
        nontriv_pure.push(c);
      }
    } else { // !isPure(c)
      impure.push(c);
    }
  }

#if DPRINT
  cout << "Found " << triv_candidates.length() << " trivial candidates." << endl;
  cout << "Found " << nontriv_pure.length() << " nontrivial pure literals." << endl;
  cout << "Found " << impure.length() << " impure literals." << endl;
#endif
  /* remove all candidates where the variable occurs in other pure
     non-trivial lits  */
  Stack<Literal*> nt_pure_tocheck(nontriv_pure);
  Stack<Literal*> nt_new;

  while( ! (nt_pure_tocheck.isEmpty() || triv_candidates.isEmpty()) ) {
    //for each candidate X=s, check if any literal in nt_pure_tocheck contains X
    //if yes, put it onto the removal list

    Stack<Literal*>::Iterator cand_it(triv_candidates);
    while(cand_it.hasNext() ) {
      Literal* cand = cand_it.next();
      Stack<Literal*>::Iterator tocheck_it(nt_pure_tocheck);
      while (tocheck_it.hasNext()) {
        Literal* checklit = tocheck_it.next();
        if (literalContainsVar(checklit, varOfXeqTerm(cand))) {
          nt_new.push(cand);
        }
        if(cand->isTwoVarEquality() && literalContainsVar(checklit,varOfXeqTerm(cand,true))){
          nt_new.push(cand);
        }
      } // ! nt_pure_tocheck.hasNext()
    }   // ! cand_it.hasNext()
    //remove nt_new from candidates, replace tocheck by nt_new
    Stack<Literal*>::Iterator nt_new_it(nt_new);
    while(nt_new_it.hasNext()) {
      triv_candidates.remove(nt_new_it.next());
    }
    nt_pure_tocheck = nt_new;
  }

#if DPRINT
  cout << "Found " << triv_candidates.length() << " trivial literals." << endl;
#endif
  
  return triv_candidates;
}


/** 
 * Selects the theory literals to be used for instantiation. These are all non-trivial pure theory literals.
 */
Stack<Literal*> TheoryInstAndSimp::selectTheoryLiterals(Clause* cl) {
#if DPRINT
  cout << "selectTheoryLiterals in " << cl->toString() << endl;
#endif

  ASS_NEQ(_mode, Shell::Options::TheoryInstSimp::OFF);

  Stack<Literal*> trivial_lits = selectTrivialLiterals(cl);
  Stack<Literal*> out;

  for (auto lit : cl->iterLits()) {
    // TODO this is O(n^2) runtime
    if (isPure(lit) && !trivial_lits.find(lit))
      out.push(lit);
  }
  return out;
}


void TheoryInstAndSimp::filterUninterpretedPartialFunctionDeep(Stack<Literal*>& theoryLits, Stack<Literal*>& filteredLits) {
#if DPRINT
  cout << "div zero filtering checking!" << endl;
#endif
  Stack<Literal*>::BottomFirstIterator it(theoryLits);
  while(it.hasNext()) {
    Literal* lit = it.next();
#if DPRINT
    cout << "div zero filtering checking: " << lit->toString() << endl;
#endif
    bool keep_lit = true;
    SubtermIterator sit(lit);
    while(sit.hasNext() && keep_lit){
      auto ts = sit.next();
      if(ts.isTerm()){
        Term* t = ts.term();
        if( theory->isPartiallyInterpretedFunction(t)
         && theory->partiallyDefinedFunctionUndefinedForArgs(t)
            ){
          // treat this literal as uninterpreted
          keep_lit = false;
#if DPRINT
          cout << "division by zero removed: " << lit->toString() << endl;
#endif
        }
      }
    }

    if (keep_lit) {
      filteredLits.push(lit);
    }
  }
}

void TheoryInstAndSimp::ConstantCache::reset()
{ for (auto& x : iterTraits(_inner.iter())) x.value().reset(); }

Term* TheoryInstAndSimp::ConstantCache::freshConstant(SortId sort) 
{ 
  auto& cache = _inner.getOrInit(sort, [](){ 
      DEBUG("new constant cache for sort ", _inner.size());
      return SortedConstantCache(); 
    });
  return cache.freshConstant(_prefix, sort);
}

void TheoryInstAndSimp::ConstantCache::SortedConstantCache::reset() 
{ _used = 0; }

Term* TheoryInstAndSimp::ConstantCache::SortedConstantCache::freshConstant(const char* prefix, SortId sort) 
{ 
  if (_constants.size() == _used)  {
    unsigned sym = env.signature->addFreshFunction(0, prefix);
    env.signature->getFunction(sym)
                 ->setType(OperatorType::getConstantsType(sort));
    DEBUG("new constant for sort ", sort, ": ", *env.signature->getFunction(sym));
    _constants.push(Term::createConstant(sym));
  }
  return _constants[_used++];
}

class TheoryInstAndSimp::GeneralisationTree {
  TermList _introduced;
  unsigned _functor;
  Stack<GeneralisationTree> _args;
public:
  GeneralisationTree(Term* name, TermList toAbstract, ConstantCache& cache) 
    : _introduced(TermList(name))
    , _functor(toAbstract.term()->functor())
    , _args(toAbstract.term()->numTermArguments())
  {
    for (unsigned i = 0; i < toAbstract.term()->numTermArguments(); i++) {
      auto arg  = toAbstract.term()->termArg(i);
      auto sort = SortHelper::getResultSort(arg.term());
      _args.push(GeneralisationTree(cache.freshConstant(sort), arg, cache));
    }
  }

  template<class F>
  void foreachDef(F f) 
  {
    Stack<TermList> args(_args.size());
    for (auto& a : _args) {
      args.push(a._introduced);
      a.foreachDef(f);
    }
    auto definition = TermList(Term::create(_functor, args.size(), args.begin()));
    f(*this, Literal::createEquality(true, _introduced, definition, SortHelper::getResultSort(_introduced.term())));
  }

  TermList buildGeneralTerm(Set<TermList> const& usedDefs, unsigned& freshVar)
  {
    if (usedDefs.contains(_introduced)) {
      Stack<TermList> args(_args.size());
      for (auto& a : _args) {
        args.push(a.buildGeneralTerm(usedDefs, freshVar));
      }
      return TermList(Term::create(_functor, args.size(), args.begin()));

    } else {
      return TermList::var(freshVar++);
    }
  }

  TermList constant() { return _introduced; }
};
 
Option<Substitution> TheoryInstAndSimp::instantiateGeneralised(
    SkolemizedLiterals skolem, unsigned freshVar)
{
  auto negatedClause = [](Stack<SATLiteral> lits) -> SATClause*
  { 
    for (auto& lit : lits) {
      lit = lit.opposite();
    }
    return SATClause::fromStack(lits);
  };

  return _solver->scoped([&]() {
    _solver->addClause(negatedClause(skolem.lits));

    Stack<SATLiteral> theoryLits;

    _generalizationConstants.reset();
    Map<SATLiteral, TermList> definitionLiterals;
    Stack<GeneralisationTree> gens;
    // unsigned freshVar = 0;
    for (auto v : skolem.vars) {
      ASS(v < freshVar);
      // freshVar = v >= freshVar ? v + 1 : freshVar;
      auto sk = skolem.subst.apply(v);
      auto val = _solver->evaluateInModel(sk.term());
      if (!val) {
        // Failed to obtain a value; could be an algebraic number or some other currently unhandled beast...
        env.statistics->theoryInstSimpLostSolution++;
        return Option<Substitution>();
      }

      auto gen = GeneralisationTree(sk.term(), TermList(val), _generalizationConstants);
      gens.push(gen);
      gen.foreachDef([&](GeneralisationTree gen, Literal* l){
          auto named = _naming.toSAT(l);
          theoryLits.push(named);
          definitionLiterals.insert(named, gen.constant());
      });
    }

    DEBUG_CODE(auto res =) _solver->solveUnderAssumptions(theoryLits, 0);
    ASS_EQ(res, SATSolver::Status::UNSATISFIABLE)

    Set<TermList> usedDefs;
    for (auto& x : _solver->failedAssumptions()) {
      definitionLiterals
        .tryGet(x)
        .andThen([&](TermList t)
            { usedDefs.insert(t); });
    }

    for (unsigned i = 0; i < skolem.vars.size(); i++) {
      skolem.subst.rebind(skolem.vars[i], gens[i].buildGeneralTerm(usedDefs, freshVar));
    }
    return Option<Substitution>(std::move(skolem.subst));
  });
};


Option<Substitution> TheoryInstAndSimp::instantiateWithModel(SkolemizedLiterals skolem)
{
  for (auto var : skolem.vars) {
    auto ev = _solver->evaluateInModel(skolem.subst.apply(var).term());
    if (ev) {
      skolem.subst.rebind(var, ev);
    } else {
      // Failed to obtain a value; could be an algebraic number or some other currently unhandled beast...
      env.statistics->theoryInstSimpLostSolution++;
      return Option<Substitution>();
    }
  }
  return Option<Substitution>(std::move(skolem.subst));
};

template<class IterLits> TheoryInstAndSimp::SkolemizedLiterals TheoryInstAndSimp::skolemize(IterLits lits) 
{
  // Currently we just get the single solution from Z3


  // Firstly, we need to consistently replace variables by constants (i.e. Skolemize)
  // Secondly, we take the complement of each literal and consider the conjunction
  // This subst is for the consistent replacement
  Substitution subst;

  Stack<SATLiteral> skolemized;
  Stack<unsigned> vars;
  _instantiationConstants.reset();
  for (auto lit : lits) {
    // replace variables consistently by fresh constants
    DHMap<unsigned, SortId> srtMap;
    SortHelper::collectVariableSorts(lit,srtMap);
    TermVarIterator vit(lit);
    while(vit.hasNext()){
      unsigned var = vit.next();
      auto sort = srtMap.get(var);
      TermList fc;
      if(!subst.findBinding(var,fc)){
        Term* fc = _instantiationConstants.freshConstant(sort);
        ASS_EQ(SortHelper::getResultSort(fc), sort);
        subst.bind(var,fc);
        vars.push(var);
      }
    }

    lit = SubstHelper::apply(lit,subst);

    skolemized.push(_naming.toSAT(lit));
  }

  return SkolemizedLiterals {
      .lits = std::move(skolemized),
      .vars = std::move(vars),
      .subst = std::move(subst),
  };
}



VirtualIterator<Solution> TheoryInstAndSimp::getSolutions(Stack<Literal*> const& theoryLiterals, Stack<Literal*> const& guards, unsigned freshVar) {
  auto skolemized = skolemize(concatIters(
        theoryLiterals.iterFifo(),
        guards.iterFifo()
        ));
  DEBUG("skolemized: ", iterTraits(skolemized.lits.iterFifo()).map([&](SATLiteral l){ return _naming.toFO(l)->toString(); }).collect<Stack>())

  // now we can call the solver
  SATSolver::Status status = _solver->solveUnderAssumptions(skolemized.lits, 0);

  if(status == SATSolver::Status::UNSATISFIABLE) {
    DEBUG("unsat")
    return pvi(getSingletonIterator(Solution::unsat()));

  } else if(status == SATSolver::Status::SATISFIABLE) {
    DEBUG("found model: ", _solver->getModel())
    auto subst = _generalisation ? instantiateGeneralised(skolemized, freshVar)
                                 : instantiateWithModel(skolemized);
    if (subst.isSome()) {
      return pvi(getSingletonIterator(Solution(std::move(subst).unwrap())));
    } else {
      DEBUG("could not build substituion from model.")
    }
  } else {
    // SMT solving was incomplete
    DEBUG("no solution.")
  }
  return VirtualIterator<Solution>::getEmpty();
}

Clause* instantiate(Clause* original, Substitution& subst, Stack<Literal*> const& theoryLits, Splitter* splitter)
{

  RStack<Literal*> instLits;
  RStack<Literal*> resLits;

  for (Literal* lit : original->iterLits()) {
    ASS_REP(SortHelper::areSortsValid(lit), *lit);
    Literal* lit_inst = SubstHelper::apply(lit,subst);
    SubtermIterator iter(lit_inst);
    while (iter.hasNext()) {
      DEBUG_CODE(auto t =) iter.next();
      ASS_REP(t.isVar() || SortHelper::areSortsValid(t.term()), t);
    }
    ASS_REP(SortHelper::areSortsValid(lit_inst), *lit_inst);
    // ASS()
    instLits->push(lit_inst);
    // we implicitly remove all theoryLits as the solution makes their combination false
    if(!theoryLits.find(lit)){
      resLits->push(lit_inst);
    }
  }
  Clause* inst = Clause::fromStack(*instLits, GeneratingInference1(InferenceRule::INSTANTIATION,original));
  if(splitter){
    splitter->onNewClause(inst);
  }
  return Clause::fromStack(*resLits, SimplifyingInference1(InferenceRule::INTERPRETED_SIMPLIFICATION,inst));
}


struct InstanceFn
{
  Clause* operator()(Solution sol, Clause* original, 
      Stack<Literal*> const& theoryLits, 
      Stack<Literal*> const& invertedLits,
      Stack<Literal*> const& guards, 
      Splitter* splitter,
      TheoryInstAndSimp* parent, 
      bool& redundant
    )
  {
    // We delete cl as it's a theory-tautology
    if(!sol.sat) {
      // now we run SMT solver again without guarding
      if(guards.isEmpty()){
        redundant = true;
      } else {
        auto skolem = parent->skolemize(iterTraits(invertedLits.iterFifo() /* without guards !! */));
        auto status = parent->_solver->solveUnderAssumptions(skolem.lits, 0);
        // we have an unsat solution without guards
        redundant = status == SATSolver::Status::UNSATISFIABLE;
      }

      if (redundant) {
        env.statistics->theoryInstSimpTautologies++;
      }

      DEBUG("tautology")
      return nullptr;
    }

    // If the solution is empty (for any reason) there is no point performing instantiation
    if(sol.subst.isEmpty()){
      env.statistics->theoryInstSimpEmptySubstitution++;
    }
    auto res = instantiate(original, sol.subst, theoryLits, splitter);
    env.statistics->theoryInstSimp++;
    return res;
  }
};

Stack<Literal*> computeGuards(Stack<Literal*> const& lits) 
{
  /* finds the constructor for a given distructor */
  auto findConstructor = [](TermAlgebra* ta, unsigned destructor, bool predicate) -> TermAlgebraConstructor* 
  {
    // TODO get rid of this wasteful search for the right constructor, and use some sort of hashing instead
    for (auto ctor : ta->iterCons()) {
      for (unsigned i = 0; i < ctor->arity(); i++) {
        auto p = ctor->argSort(i) == AtomicSort::boolSort();
        auto d = ctor->destructorFunctor(i);
        if(destructor == d && predicate == p) 
          return ctor;
      }
    }
    ASSERTION_VIOLATION
  };

  auto destructorGuard = [&findConstructor](Term* destr, SortId sort, bool predicate) -> Literal*
  {
      auto ctor = findConstructor(env.signature->getTermAlgebraOfSort(sort), destr->functor(), predicate);
      auto discr = ctor->discriminator();
      TermStack args;
      for (unsigned i = 0; i < destr->numTypeArguments(); i++) {
        args.push(destr->typeArg(i));
      }
      args.push(destr->termArg(0));
      // asserts e.g. isCons(l) for a term that contains the subterm head(l) for lists
      return Literal::create(discr, args.size(), /* polarity */ true, args.begin());
  };


  Stack<Literal*> out;
  for (auto lit : lits) {

    /* guards for predicates */
    auto predSym = env.signature->getPredicate(lit->functor());
    if (predSym->termAlgebraDest()) {
      out.push(destructorGuard(lit, predSym->predType()->arg(0), /* predicate */ true));
    }

    /* guards for subterms */
    SubtermIterator it(lit);
    for (auto t = it.next(); it.hasNext(); t = it.next()) {
      ASS_REP(t.isVar() || !t.term()->isLiteral(), t);
      if (t.isTerm()) {
        auto term = t.term();
        auto sym = env.signature->getFunction(t.term()->functor());
        auto fun = term->functor();
        if (theory->isInterpretedNumber(term)) {
          /* no guard */
        } else if (theory->isInterpretedFunction(fun) || theory->isInterpretedConstant(fun)) {

          switch (theory->interpretFunction(fun)) {
            case Theory::REAL_QUOTIENT:
            case Theory::REAL_REMAINDER_E:
            case Theory::REAL_QUOTIENT_F:
            case Theory::REAL_QUOTIENT_T:
            case Theory::REAL_REMAINDER_T:
            case Theory::REAL_REMAINDER_F:
              out.push(Literal::createEquality(false, RealTraits::zero(), term->termArg(1), RealTraits::sort()));
              break;

            case Theory::RAT_QUOTIENT:
            case Theory::RAT_QUOTIENT_T:
            case Theory::RAT_REMAINDER_T:
            case Theory::RAT_QUOTIENT_F:
            case Theory::RAT_REMAINDER_F:
            case Theory::RAT_REMAINDER_E:
              out.push(Literal::createEquality(false, RatTraits::zero(), term->termArg(1), RatTraits::sort()));
              break;

            case Theory::INT_QUOTIENT_F:
            case Theory::INT_REMAINDER_F:
            case Theory::INT_QUOTIENT_E:
            case Theory::INT_QUOTIENT_T:
            case Theory::INT_REMAINDER_T:
            case Theory::INT_REMAINDER_E:
              out.push(Literal::createEquality(false, IntTraits::zero(), term->termArg(1), IntTraits::sort()));
              break;

            default:; /* no guard */
          }
        } else if (sym->termAlgebraDest()) {
          out.push(destructorGuard(term, sym->fnType()->arg(0), /* predicate */ false));
        }
      }
    }
  }
  return out;
}

Stack<Literal*> filterLiterals(Stack<Literal*> lits, Options::TheoryInstSimp mode)
{
  auto isStrong = [](Literal* lit)  -> bool
  { return ( lit->isEquality() && lit->isNegative())
        || (!lit->isEquality() && theory->isInterpretedPredicate(lit->functor())); };

  auto freeVars = [](Literal* lit)
  { return iterTraits(VList::Iterator(freeVariables(lit))); };

  switch(mode) {
    case Options::TheoryInstSimp::ALL:
      return lits;

    case Options::TheoryInstSimp::STRONG:
      return iterTraits(lits.iterFifo())
                            .filter(isStrong)
                            .collect<Stack>();

    case Options::TheoryInstSimp::NEG_EQ:
      return iterTraits(lits.iterFifo())
                            .filter([](Literal* lit)
                               { return lit->isEquality() && lit->isNegative(); } )
                            .collect<Stack>();

    case Options::TheoryInstSimp::OVERLAP:
      {
        Set<unsigned> strongVars;

        for (auto l : lits) {
          if (isStrong(l)) {
            for (auto v : freeVars(l)) {
              strongVars.insert(v);
            }
          }
        }

        return iterTraits(lits.iterFifo())
                            .filter([&](Literal* lit) 
                                { return freeVars(lit)
                                          .any([&](unsigned var) 
                                              { return strongVars.contains(var); }); })
                            .collect<Stack>();
      }

    case Options::TheoryInstSimp::FULL:
    case Options::TheoryInstSimp::NEW:
    case Options::TheoryInstSimp::OFF:
      ASSERTION_VIOLATION
  }
  ASSERTION_VIOLATION
}

unsigned getFreshVar(Clause& clause) 
{
  unsigned freshVar = 0;
  for (unsigned i = 0; i < clause.size(); i++) {
    VariableIterator iter(clause[i]);
    while(iter.hasNext()) {
      auto var = iter.next();
      if (freshVar <= var.var()) {
        freshVar = var.var() + 1;
      }
    }
  }
  return freshVar;
}

SimplifyingGeneratingInference::ClauseGenerationResult TheoryInstAndSimp::generateSimplify(Clause* premise)
{
  auto empty = ClauseGenerationResult {
    .clauses          = ClauseIterator::getEmpty(),
    .premiseRedundant = false,
  };

  if(premise->isPureTheoryDescendant()){ 
    return empty;
  }


  Stack<Literal*> selectedLiterals = selectTheoryLiterals(premise);
  selectedLiterals = filterLiterals(std::move(selectedLiterals), _mode);

  // if there are no eligable theory literals selected then there is nothing to do
  if(selectedLiterals.isEmpty()){
    return empty;
  }

  // we have an eligable candidate
  env.statistics->theoryInstSimpCandidates++;
#if VTIME_PROFILING
  static const char* THEORY_INST_SIMP = "theory instantiation";
#endif // VTIME_PROFILING

  auto guards = computeGuards(selectedLiterals);

  DEBUG("input:             ", *premise);
  DEBUG("selected literals: ", iterTraits(selectedLiterals.iterFifo())
                                 .map([](Literal* l) -> std::string { return l->toString(); })
                                 .collect<Stack>())
  DEBUG("guards:            ", iterTraits(guards.iterFifo())
                                 .map([](Literal* l) -> std::string { return l->toString(); })
                                 .collect<Stack>())
  TIME_TRACE(THEORY_INST_SIMP);

  auto invertedLiterals = iterTraits(selectedLiterals.iterFifo())
    .map(Literal::complementaryLiteral)
    .collect<Stack>();

  bool premiseRedundant = false;

  auto it1 = iterTraits(getSolutions(invertedLiterals, guards, getFreshVar(*premise)))
    .map([&](Solution s)  { 
        DEBUG("found solution: ", s); 
        return InstanceFn{}(s, premise, selectedLiterals, invertedLiterals, guards, _splitter, this, premiseRedundant);
    })
    .filter([](Clause* cl) { return cl != nullptr; });

  auto it2 = TIME_TRACE_ITER(THEORY_INST_SIMP, it1);

  // we need to strictily evaluate the iterator to 
  auto clauses =  getPersistentIterator(it2);

  if (premiseRedundant && _thiTautologyDeletion) {
    return ClauseGenerationResult {
      .clauses          = ClauseIterator::getEmpty(),
      .premiseRedundant = true,
    };
  } else {
    return ClauseGenerationResult {
      .clauses          = clauses,
      .premiseRedundant = false,
    };
  }
}

std::ostream& operator<<(std::ostream& out, Solution const& self) 
{ return out << "Solution(" << (self.sat ? "sat" : "unsat") << ", " << self.subst << ")"; }

}

#endif
