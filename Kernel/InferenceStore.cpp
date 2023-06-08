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
 * @file InferenceStore.cpp
 * Implements class InferenceStore.
 */

#include "Lib/Allocator.hpp"
#include "Lib/DHSet.hpp"
#include "Lib/Environment.hpp"
#include "Lib/Int.hpp"
#include "Lib/ScopedPtr.hpp"
#include "Lib/SharedSet.hpp"
#include "Lib/Stack.hpp"
#include "Lib/StringUtils.hpp"
#include "Lib/ScopedPtr.hpp"

#include "Shell/Dedukti.hpp"
#include "Shell/LaTeX.hpp"
#include "Shell/Options.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/UIHelper.hpp"

#include "Parse/TPTP.hpp"

#include "Saturation/Splitter.hpp"

#include "Signature.hpp"
#include "Clause.hpp"
#include "Formula.hpp"
#include "FormulaUnit.hpp"
#include "FormulaVarIterator.hpp"
#include "Inference.hpp"
#include "Term.hpp"
#include "TermIterators.hpp"
#include "SortHelper.hpp"

#include "InferenceStore.hpp"

//TODO: when we delete clause, we should also delete all its records from the inference store

namespace Kernel
{

using namespace std;
using namespace Lib;
using namespace Shell;

void InferenceStore::FullInference::increasePremiseRefCounters()
{
  for(unsigned i=0;i<premCnt;i++) {
    if (premises[i]->isClause()) {
      premises[i]->incRefCnt();
    }
  }
}

/**
 * Records informations needed for outputting proofs of general splitting
 */
void InferenceStore::recordSplittingNameLiteral(Unit* us, Literal* lit)
{
  //each clause is result of a splitting only once
  ALWAYS(_splittingNameLiterals.insert(us, lit));
}


/**
 * Record the introduction of a new symbol
 */
void InferenceStore::recordIntroducedSymbol(Unit* u, SymbolType st, unsigned number)
{
  SymbolStack* pStack;
  _introducedSymbols.getValuePtr(u->number(),pStack);
  pStack->push(SymbolId(st,number));
}

/**
 * Record the introduction of a split name
 */
void InferenceStore::recordIntroducedSplitName(Unit* u, std::string name)
{
  ALWAYS(_introducedSplitNames.insert(u->number(),name));
}

/**
 * Return @b inner quantified over variables in @b vars
 *
 * It is caller's responsibility to ensure that variables in @b vars are unique.
 */
template<typename VarContainer>
std::string getQuantifiedStr(const VarContainer& vars, std::string inner, DHMap<unsigned,TermList>& t_map, bool innerParentheses=true){
  VirtualIterator<unsigned> vit=pvi( getContentIterator(vars) );
  std::string varStr;
  bool first=true;
  while(vit.hasNext()) {
    unsigned var =vit.next();
    std::string ty="";
    TermList t;

    if(t_map.find(var,t) && env.getMainProblem()->hasNonDefaultSorts()){
      //hasNonDefaultSorts is true if the problem contains a sort
      //that is not $i and not a variable
      ty=" : " + t.toString();
    }
    if(ty == " : $tType"){
      if (!first) { varStr = "," + varStr; }
      varStr=std::string("X")+Int::toString(var)+ty + varStr;
    } else {
      if (!first) { varStr+=","; }
      varStr+=std::string("X")+Int::toString(var)+ty;
    }
    first=false;
  }

  if (first) {
    //we didn't quantify any variable
    return inner;
  }

  if (innerParentheses) {
    return "( ! ["+varStr+"] : ("+inner+") )";
  }
  else {
    return "( ! ["+varStr+"] : "+inner+" )";
  }
}

/**
 * Return @b inner quentified over variables in @b vars
 *
 * It is caller's responsibility to ensure that variables in @b vars are unique.
 */
template<typename VarContainer>
std::string getQuantifiedStr(const VarContainer& vars, std::string inner, bool innerParentheses=true)
{
  static DHMap<unsigned,TermList> d;
  return getQuantifiedStr(vars,inner,d,innerParentheses);
}

/**
 * Return std::string containing quantified unit @b u.
 */
std::string getQuantifiedStr(Unit* u, List<unsigned>* nonQuantified=0)
{
  Set<unsigned> vars;
  std::string res;
  DHMap<unsigned,TermList> t_map;
  SortHelper::collectVariableSorts(u,t_map);
  if (u->isClause()) {
    Clause* cl=static_cast<Clause*>(u);
    unsigned clen=cl->length();
    for(unsigned i=0;i<clen;i++) {
      TermVarIterator vit( (*cl)[i] ); //TODO update iterator for two var lits?
      while(vit.hasNext()) {
        unsigned var=vit.next();
        if (List<unsigned>::member(var, nonQuantified)) {
          continue;
        }
        vars.insert(var);
      }
    }
    res=cl->literalsOnlyToString();
  } else {
    Formula* formula=static_cast<FormulaUnit*>(u)->formula();
    FormulaVarIterator fvit( formula );
    while(fvit.hasNext()) {
      unsigned var=fvit.next();
      if (List<unsigned>::member(var, nonQuantified)) {
        continue;
      }
      vars.insert(var);
    }
    res=formula->toString();
  }

  return getQuantifiedStr(vars, res, t_map);
}

struct InferenceStore::ProofPrinter
{
  ProofPrinter(ostream& out, InferenceStore* is)
  : _is(is), out(out)
  {
    outputAxiomNames=env.options->outputAxiomNames();
    delayPrinting=true;
  }

  void scheduleForPrinting(Unit* us)
  {
    outKernel.push(us);
    handledKernel.insert(us);
  }

  virtual ~ProofPrinter() {}

  virtual void print()
  {
    while(outKernel.isNonEmpty()) {
      Unit* cs=outKernel.pop();
      handleStep(cs);
    }
    out << fu->formula()->toString() << ' ';

    out <<"["<<Kernel::ruleName(rule);

  virtual bool hideProofStep(InferenceRule rule)
  {
    return false;
  }

  void requestProofStep(Unit* prem)
  {
    if (!handledKernel.contains(prem)) {
      handledKernel.insert(prem);
      outKernel.push(prem);
    }
  }

  virtual void printStep(Unit* cs)
  {
    if (cs->isClause()) {
      Clause* cl=cs->asClause();
      out << cl->toString() << "\n";
    } else {
      InferenceRule rule = cs->inference().rule();
      UnitIterator parents= cs->getParents();

      out << cs->number() << ". ";
      FormulaUnit* fu=static_cast<FormulaUnit*>(cs);
      if (env.colorUsed && fu->inheritedColor() != COLOR_INVALID) {
        out << " IC" << fu->inheritedColor() << " ";
      }
      out << fu->formula()->toString() << ' ';

      out <<"["<<cs->inference().name();

      if (outputAxiomNames && rule==InferenceRule::INPUT) {
        ASS(!parents.hasNext()); //input clauses don't have parents
        std::string name;
        if (Parse::TPTP::findAxiomName(cs, name)) {
          out << " " << name;
        }
      }

      bool first=true;
      while(parents.hasNext()) {
        Unit* prem=parents.next();
        out << (first ? ' ' : ',');
        out << prem->number();
        first=false;
      }
      out << "]" << endl;
    }
  }

  void handleStep(Unit* cs)
  {
    InferenceRule rule = cs->inference().rule();
    UnitIterator parents= cs->getParents();

    bool first=true;
    while(parents.hasNext()) {
      Unit* prem=parents.next();
      out << (first ? ' ' : ',');
      out << _is->getUnitIdStr(prem);
      first=false;
    }
    out << "]" << endl;
  }
}

void InferenceStore::ProofPrinter::handleStep(Unit* cs)
{
  CALL("InferenceStore::ProofPrinter::handleStep");
  InferenceRule rule;
  UnitIterator parents=_is->getParents(cs, rule);

  while(parents.hasNext()) {
    Unit* prem=parents.next();
    ASS(prem!=cs);
    requestProofStep(prem);
  }

  void printDelayed()
  {
    // Sort
    sort(
      delayed.begin(), delayed.end(),
      [](Unit *u1, Unit *u2) -> bool { return u1->number() < u2->number(); }
    );

  // Sort
  sort<UnitNumberComparator>(delayed.begin(),delayed.end());

  // Print
  for(unsigned i=0;i<delayed.size();i++){
    printStep(delayed[i]);
  }

}

struct InferenceStore::ProofPropertyPrinter
: public InferenceStore::ProofPrinter
{
  ProofPropertyPrinter(ostream& out, InferenceStore* is) : ProofPrinter(out,is)
  {
    max_theory_clause_depth = 0;
    for(unsigned i=0;i<11;i++){ buckets.push(0); }
    last_one = false;
  }

  void print()
  {
    ProofPrinter::print();
    for(unsigned i=0;i<11;i++){ out << buckets[i] << " ";}
    out << endl;
    if(last_one){ out << "yes" << endl; }
    else{ out << "no" << endl; }
  }

protected:

  void printStep(Unit* us)
  {
    static unsigned lastP = Unit::getLastParsingNumber();
    static float chunk = lastP / 10.0;
    if(us->number() <= lastP){
      if(us->number() == lastP){ 
        last_one = true;
      }
      unsigned bucket = (unsigned)(us->number() / chunk);
      buckets[bucket]++;
    }

    // TODO we could make clauses track this information, but I am not sure that that's worth it
    if(us->isClause() && us->isPureTheoryDescendant()){
      //cout << "HERE with " << us->toString() << endl;
      Inference* inf = &us->inference();
      while(inf->rule() == InferenceRule::EVALUATION){
        Inference::Iterator piit = inf->iterator();
        inf = &inf->next(piit)->inference();
      }
      Stack<Inference*> current;
      current.push(inf);
      unsigned level = 0;
      while(!current.isEmpty()){
        //cout << current.size() << endl;
        Stack<Inference*> next;
        Stack<Inference*>::Iterator it(current);
        while(it.hasNext()){
          Inference* inf = it.next();
          Inference::Iterator iit=inf->iterator();
          while(inf->hasNext(iit)) {
            Unit* premUnit=inf->next(iit);
            Inference* premInf = &premUnit->inference();
            while(premInf->rule() == InferenceRule::EVALUATION){
              Inference::Iterator piit = premInf->iterator();
              premUnit = premInf->next(piit);
              premInf = &premUnit->inference();
            }

//for(unsigned i=0;i<level;i++){ cout << ">";}; cout << premUnit->toString() << endl;
            next.push(premInf);
          }
        }
        level++;
        current = next;
      }
      level--;
      //cout << "level is " << level << endl;
      
      if(level > max_theory_clause_depth){
        max_theory_clause_depth=level;
      }
    }
  }

  unsigned max_theory_clause_depth;
  bool last_one;
  Stack<unsigned> buckets;

};


struct InferenceStore::TPTPProofPrinter
: public InferenceStore::ProofPrinter
{
  TPTPProofPrinter(ostream& out, InferenceStore* is)
  : ProofPrinter(out, is) {
    splitPrefix = Saturation::Splitter::splPrefix;
    // Don't delay printing in TPTP proof mode
    delayPrinting = false;
  }

  void print()
  {
    //outputSymbolDeclarations also deals with sorts for now
    //UIHelper::outputSortDeclarations(out);
    UIHelper::outputSymbolDeclarations(out);
    ProofPrinter::print();
  }

protected:
  std::string splitPrefix;

  std::string getRole(InferenceRule rule, UnitInputType origin)
  {
    switch(rule) {
    case InferenceRule::INPUT:
      if (origin==UnitInputType::CONJECTURE) {
        return "conjecture";
      } else if (origin==UnitInputType::NEGATED_CONJECTURE) {
        return "negated_conjecture";
      } else {
        return "axiom";
      }
    case InferenceRule::NEGATED_CONJECTURE:
      return "negated_conjecture";
    default:
      return "plain";
    }
  }

  std::string tptpRuleName(InferenceRule rule)
  {
    return StringUtils::replaceChar(ruleName(rule), ' ', '_');
  }

  std::string unitIdToTptp(std::string unitId)
  {
    return "f"+unitId;
  }

  std::string tptpUnitId(Unit* us)
  {
    return unitIdToTptp(Int::toString(us->number()));
  }

  std::string tptpDefId(Unit* us)
  {
    return unitIdToTptp(Int::toString(us->number())+"_D");
  }

  std::string splitsToString(SplitSet* splits)
  {
    ASS_G(splits->size(),0);

    if (splits->size()==1) {
      return Saturation::Splitter::getFormulaStringFromName(splits->sval(),true /*negated*/);
    }
    auto sit = splits->iter();
    std::string res("(");
    while(sit.hasNext()) {
      res+= Saturation::Splitter::getFormulaStringFromName(sit.next(),true /*negated*/);
      if (sit.hasNext()) {
	res+=" | ";
      }
    }
    res+=")";
    return res;
  }

  std::string quoteAxiomName(std::string n)
  {
    static std::string allowedFirst("0123456789abcdefghijklmnopqrstuvwxyz");
    const char* allowed="_ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz";

    if (n.size()==0 || allowedFirst.find(n[0])==std::string::npos ||
	n.find_first_not_of(allowed)!=std::string::npos) {
      n='\''+n+'\'';
    }
    return n;
  }

  std::string getFofString(std::string id, std::string formula, std::string inference, InferenceRule rule, UnitInputType origin=UnitInputType::AXIOM)
  {
    std::string kind = "fof";
    if(env.getMainProblem()->hasNonDefaultSorts()){ kind="tff"; }
    if(env.getMainProblem()->isHigherOrder()){ kind="thf"; }

    return kind+"("+id+","+getRole(rule,origin)+",("+"\n"
	+"  "+formula+"),\n"
	+"  "+inference+").";
  }

  std::string getFormulaString(Unit* us)
  {
    std::string formulaStr;
    if (us->isClause()) {
      Clause* cl=us->asClause();
      formulaStr=getQuantifiedStr(cl);
      if (cl->splits() && !cl->splits()->isEmpty()) {
	formulaStr+=" | "+splitsToString(cl->splits());
      }
    }
    else {
      FormulaUnit* fu=static_cast<FormulaUnit*>(us);
      formulaStr=getQuantifiedStr(fu);
    }
    return formulaStr;
  }

  bool hasNewSymbols(Unit* u) {
    bool res = _is->_introducedSymbols.find(u->number());
    ASS(!res || _is->_introducedSymbols.get(u->number()).isNonEmpty());
    if(!res){
      res = _is->_introducedSplitNames.find(u->number());
    }
    return res;
  }
  std::string getNewSymbols(std::string origin, std::string symStr) {
    return "new_symbols(" + origin + ",[" +symStr + "])";
  }
  /** It is an iterator over SymbolId */
  template<class It>
  std::string getNewSymbols(std::string origin, It symIt) {
    std::ostringstream symsStr;
    while(symIt.hasNext()) {
      SymbolId sym = symIt.next();
      if (sym.first == SymbolType::FUNC ) {
        symsStr << env.signature->functionName(sym.second);
      } else if (sym.first == SymbolType::PRED){
        symsStr << env.signature->predicateName(sym.second);
      } else {
        symsStr << env.signature->typeConName(sym.second);
      }
      if (symIt.hasNext()) {
        symsStr << ',';
      }
    }
    return getNewSymbols(origin, symsStr.str());
  }
  std::string getNewSymbols(std::string origin, Unit* u) {
    ASS(hasNewSymbols(u));

    if(_is->_introducedSplitNames.find(u->number())){
      return getNewSymbols(origin,_is->_introducedSplitNames.get(u->number()));
    }

    SymbolStack& syms = _is->_introducedSymbols.get(u->number());
    return getNewSymbols(origin, SymbolStack::ConstIterator(syms));
  }

  void printStep(Unit* us)
  {
    InferenceRule rule = us->inference().rule();
    UnitIterator parents= us->getParents();

    switch(rule) {
    //case Inference::AVATAR_COMPONENT:
    //  printSplittingComponentIntroduction(us);
    //  return;
    case InferenceRule::GENERAL_SPLITTING_COMPONENT:
      printGeneralSplittingComponent(us);
      return;
    case InferenceRule::GENERAL_SPLITTING:
      printSplitting(us);
      return;
    default: ;
    }

    //get std::string representing the formula

    std::string formulaStr=getFormulaString(us);

    //get inference std::string

    std::string inferenceStr;
    if (rule==InferenceRule::INPUT) {
      std::string fileName;
      if (env.options->inputFile()=="") {
	      fileName="unknown";
      }
      else {
	      fileName="'"+env.options->inputFile()+"'";
      }
      std::string axiomName;
      if (!outputAxiomNames || !Parse::TPTP::findAxiomName(us, axiomName)) {
	      axiomName="unknown";
      }
      inferenceStr="file("+fileName+","+quoteAxiomName(axiomName)+")";
    }
    else if (!parents.hasNext()) {
      std::string newSymbolInfo;
      if (hasNewSymbols(us)) {
        std::string newSymbOrigin;
        if (rule == InferenceRule::FUNCTION_DEFINITION ||
          rule == InferenceRule::FOOL_ITE_DEFINITION || rule == InferenceRule::FOOL_LET_DEFINITION ||
          rule == InferenceRule::FOOL_FORMULA_DEFINITION || rule == InferenceRule::FOOL_MATCH_DEFINITION) {
          newSymbOrigin = "definition";
        } else {
          newSymbOrigin = "naming";
        }
	      newSymbolInfo = getNewSymbols(newSymbOrigin,us);
      }
      inferenceStr="introduced("+tptpRuleName(rule)+",["+newSymbolInfo+"])";
    }
    else {
      ASS(parents.hasNext());
      std::string statusStr;
      if (rule==InferenceRule::SKOLEMIZE) {
	      statusStr="status(esa),"+getNewSymbols("skolem",us);
      }

      inferenceStr="inference("+tptpRuleName(rule);

      inferenceStr+=",["+statusStr+"],[";
      bool first=true;
      while(parents.hasNext()) {
        Unit* prem=parents.next();
        if (!first) {
          inferenceStr+=',';
        }
        inferenceStr+=tptpUnitId(prem);
        first=false;
      }
      inferenceStr+="])";
    }

    out<<getFofString(tptpUnitId(us), formulaStr, inferenceStr, rule, us->inputType())<<endl;
  }

  void printSplitting(Unit* us)
  {
    ASS(us->isClause());

    InferenceRule rule = us->inference().rule();
    UnitIterator parents= us->getParents();
    ASS(rule==InferenceRule::GENERAL_SPLITTING);

    std::string inferenceStr="inference("+tptpRuleName(rule)+",[],[";

    //here we rely on the fact that the base premise is always put as the first premise in
    //GeneralSplitting::apply

    ALWAYS(parents.hasNext());
    Unit* base=parents.next();
    inferenceStr+=tptpUnitId(base);

    ASS(parents.hasNext()); //we always split off at least one component
    while(parents.hasNext()) {
      Unit* comp=parents.next();
      ASS(_is->_splittingNameLiterals.find(comp));
      inferenceStr+=","+tptpDefId(comp);
    }
    inferenceStr+="])";

    out<<getFofString(tptpUnitId(us), getFormulaString(us), inferenceStr, rule)<<endl;
  }

  void printGeneralSplittingComponent(Unit* us)
  {
    ASS(us->isClause());

    InferenceRule rule = us->inference().rule();
    UnitIterator parents= us->getParents();
    ASS(!parents.hasNext());

    Literal* nameLit=_is->_splittingNameLiterals.get(us); //the name literal must always be stored

    std::string defId=tptpDefId(us);

    out<<getFofString(tptpUnitId(us), getFormulaString(us),
	    "inference("+tptpRuleName(InferenceRule::CLAUSIFY)+",[],["+defId+"])", InferenceRule::CLAUSIFY)<<endl;


    List<unsigned>* nameVars=0;
    VariableIterator vit(nameLit);
    while(vit.hasNext()) {
      unsigned var=vit.next().var();
      ASS(!List<unsigned>::member(var, nameVars)); //each variable appears only once in the naming literal
      List<unsigned>::push(var,nameVars);
    }

    std::string compStr;
    List<unsigned>* compOnlyVars=0;
    bool first=true;
    bool multiple=false;
    for (Literal* lit : us->asClause()->iterLits()) {
      if (lit==nameLit) {
	      continue;
      }
      if (first) {
	      first=false;
      }
      else {
	      multiple=true;
	      compStr+=" | ";
      }
      compStr+=lit->toString();

      VariableIterator lvit(lit);
      while(lvit.hasNext()) {
        unsigned var=lvit.next().var();
        if (!List<unsigned>::member(var, nameVars) && !List<unsigned>::member(var, compOnlyVars)) {
          List<unsigned>::push(var,compOnlyVars);
        }
      }
    }
    ASS(!first);

    compStr=getQuantifiedStr(compOnlyVars, compStr, multiple);
    List<unsigned>::destroy(compOnlyVars);

    std::string defStr=compStr+" <=> "+Literal::complementaryLiteral(nameLit)->toString();
    defStr=getQuantifiedStr(nameVars, defStr);
    List<unsigned>::destroy(nameVars);

    SymbolId nameSymbol = SymbolId(SymbolType::PRED,nameLit->functor());
    std::ostringstream originStm;
    originStm << "introduced(" << tptpRuleName(rule)
	      << ",[" << getNewSymbols("naming",getSingletonIterator(nameSymbol))
	      << "])";

    out<<getFofString(defId, defStr, originStm.str(), rule)<<endl;
  }

  void printSplittingComponentIntroduction(Unit* us)
  {
    ASS(us->isClause());

    Clause* cl=us->asClause();
    ASS(cl->splits());
    ASS_EQ(cl->splits()->size(),1);

    InferenceRule rule=InferenceRule::AVATAR_COMPONENT;

    std::string defId=tptpDefId(us);
    std::string splitPred = splitsToString(cl->splits());
    std::string defStr=getQuantifiedStr(cl)+" <=> ~"+splitPred;

    out<<getFofString(tptpUnitId(us), getFormulaString(us),
      "inference("+tptpRuleName(InferenceRule::CLAUSIFY)+",[],["+defId+"])", InferenceRule::CLAUSIFY)<<endl;

    std::stringstream originStm;
    originStm << "introduced(" << tptpRuleName(rule)
        << ",[" << getNewSymbols("naming",splitPred)
        << "])";

    out<<getFofString(defId, defStr, originStm.str(), rule)<<endl;
  }

};

struct InferenceStore::ProofCheckPrinter
: public InferenceStore::ProofPrinter
{
  ProofCheckPrinter(ostream& out, InferenceStore* is)
  : ProofPrinter(out, is) {}

protected:
  void printStep(Unit* cs)
  {
    InferenceRule rule = cs->inference().rule();
    UnitIterator parents= cs->getParents();

    //outputSymbolDeclarations also deals with sorts for now
    //UIHelper::outputSortDeclarations(out);
    UIHelper::outputSymbolDeclarations(out);

    std::string kind = "fof";
    if(env.getMainProblem()->hasNonDefaultSorts()){ kind="tff"; }
    if(env.getMainProblem()->isHigherOrder()){ kind="thf"; }

    out << kind
        << "(r"<< cs->number()
    	<< ",conjecture, "
    	<< getQuantifiedStr(cs)
    	<< " ). %"<<ruleName(rule)<<"\n";

    while(parents.hasNext()) {
      Unit* prem=parents.next();
      out << kind
        << "(pr"<<prem->number()
  	<< ",axiom, "
  	<< getQuantifiedStr(prem);
      out << " ).\n";
    }
    out << "%#\n";
  }


  bool hideProofStep(InferenceRule rule)
  {
    switch(rule) {
    case InferenceRule::INPUT:
    case InferenceRule::INEQUALITY_SPLITTING_NAME_INTRODUCTION:
    case InferenceRule::INEQUALITY_SPLITTING:
    case InferenceRule::SKOLEMIZE:
    case InferenceRule::EQUALITY_PROXY_REPLACEMENT:
    case InferenceRule::EQUALITY_PROXY_AXIOM1:
    case InferenceRule::EQUALITY_PROXY_AXIOM2:
    case InferenceRule::NEGATED_CONJECTURE:
    case InferenceRule::RECTIFY:
    case InferenceRule::FLATTEN:
    case InferenceRule::ENNF:
    case InferenceRule::NNF:
    case InferenceRule::CLAUSIFY:
    case InferenceRule::AVATAR_DEFINITION:
    case InferenceRule::AVATAR_COMPONENT:
    case InferenceRule::AVATAR_REFUTATION:
    case InferenceRule::AVATAR_SPLIT_CLAUSE:
    case InferenceRule::AVATAR_CONTRADICTION_CLAUSE:
    case InferenceRule::FOOL_LET_DEFINITION:
    case InferenceRule::FOOL_ITE_DEFINITION:
    case InferenceRule::FOOL_ELIMINATION:
    case InferenceRule::BOOLEAN_TERM_ENCODING:
    case InferenceRule::CHOICE_AXIOM:
    case InferenceRule::PREDICATE_DEFINITION:
      return true;
    default:
      return false;
    }
  }

  void print()
  {
    ProofPrinter::print();
    out << "%#\n";
  }
};

InferenceStore::ProofPrinter* InferenceStore::createProofPrinter(ostream& out)
{
  switch(env.options->proof()) {
  case Options::Proof::ON:
    return new ProofPrinter(out, this);
  case Options::Proof::PROOFCHECK:
    return new ProofCheckPrinter(out, this);
  case Options::Proof::TPTP:
    return new TPTPProofPrinter(out, this);
  case Options::Proof::PROPERTY:
    return new ProofPropertyPrinter(out,this);
  case Options::Proof::DEDUKTI:
    return new Dedukti::ProofPrinter(out, this);
  case Options::Proof::OFF:
    return nullptr;
  }
  ASSERTION_VIOLATION;
}

/**
 * Output a proof of refutation to out
 *
 *
 */
void InferenceStore::outputUnsatCore(ostream& out, Unit* refutation)
{
  out << "(" << endl;

  Stack<Unit*> todo;
  todo.push(refutation);
  Set<Unit*> visited;
  while(!todo.isEmpty()){

    Unit* u = todo.pop();
    visited.insert(u);

    if(u->inference().rule() ==  InferenceRule::INPUT){
      if(!u->isClause()){
        if(u->getFormula()->hasLabel()){
          std::string label =  u->getFormula()->getLabel();
          out << label << endl;
        }
        else{
          ASS(env.options->ignoreMissingInputsInUnsatCore() || u->getFormula()->hasLabel());
          if(!(env.options->ignoreMissingInputsInUnsatCore() || u->getFormula()->hasLabel())){
            cout << "ERROR: There is a problem with the unsat core. There is an input formula in the proof" <<  endl;
            cout << "that does not have a label. We expect all  input formulas to have labels as this  is what" << endl;
            cout << "smtcomp does. If you don't want this then use the ignore_missing_inputs_in_unsat_core option" << endl;
            cout << "The unlabelled  input formula is " << endl;
            cout << u->toString() << endl;
          }
        }
      }
      else{
        //Currently ignore clauses as they cannot come from SMT-LIB as input formulas
      }
    }
    else{
      UnitIterator parents = u->getParents();
      while(parents.hasNext()){
        Unit* parent = parents.next();
        if(!visited.contains(parent)){
          todo.push(parent);
        }
      }
    }
  }

  out << ")" << endl;
}



/**
 * Output a proof of refutation to out
 *
 *
 */
void InferenceStore::outputProof(ostream& out, Unit* refutation)
{
  ProofPrinter* p = createProofPrinter(out);
  if (!p) {
    return;
  }
  ScopedPtr<ProofPrinter> pp(p);
  pp->scheduleForPrinting(refutation);
  pp->print();
}

/**
 * Output a proof of units to out
 *
 */
void InferenceStore::outputProof(ostream& out, UnitList* units)
{
  ProofPrinter* p = createProofPrinter(out);
  if (!p) {
    return;
  }
  ScopedPtr<ProofPrinter> pp(p);
  UnitList::Iterator uit(units);
  while(uit.hasNext()) {
    Unit* u = uit.next();
    pp->scheduleForPrinting(u);
  }
  pp->print();
}

InferenceStore* InferenceStore::instance()
{
  static ScopedPtr<InferenceStore> inst(new InferenceStore());

  return inst.ptr();
}

}
