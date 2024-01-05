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
 * @file Z3Interfacing.hpp
 * Defines class Z3Interfacing
 */

#ifndef __Z3Interfacing__
#define __Z3Interfacing__

#if VZ3

/* an (imperfect and under development) version of tracing the Z3 interface
 *  so that vampire can be "factored-out" of runs which cause particular Z3
 *  behaviour. Should be useful for producing MWEs for the Z3 people.
 */
#define PRINT_CPP(X) // std::cout << X << std::endl;

#include <fstream>

#include "Lib/DHMap.hpp"
#include "Lib/Option.hpp"
#include "Lib/BiMap.hpp"
#include "Lib/Set.hpp"

#include "SATSolver.hpp"
#include "SATLiteral.hpp"
#include "SATClause.hpp"
#include "SATInference.hpp"
#include "SAT2FO.hpp"
#include "Lib/Option.hpp"
#include "Lib/Coproduct.hpp"

#include "Kernel/Signature.hpp"

#define __EXCEPTIONS 1
#include "z3++.h"
#include "z3_api.h"

namespace SAT{

  struct UninterpretedForZ3Exception : public ThrowableBase
  {
    UninterpretedForZ3Exception()
    {
    }
  };

class Z3Interfacing : public PrimitiveProofRecordingSATSolver
{
public:
  Z3Interfacing(const Shell::Options& opts, SAT2FO& s2f, bool unsatCoresForAssumptions, vstring const& exportSmtlib);
  Z3Interfacing(SAT2FO& s2f, bool showZ3, bool unsatCoresForAssumptions, vstring const& exportSmtlib);
  ~Z3Interfacing();

  static char const* z3_full_version();

  void addClause(SATClause* cl) override;

  Status solve();
  virtual Status solve(unsigned conflictCountLimit) override { return solve(); };
  /**
   * If status is @c SATISFIABLE, return assignment of variable @c var
   */
  virtual VarAssignment getAssignment(unsigned var) override;

  /**
   * If status is @c SATISFIABLE, return 0 if the assignment of @c var is
   * implied only by unit propagation (i.e. does not depend on any decisions)
   */
  virtual bool isZeroImplied(unsigned var) override;
  /**
   * Collect zero-implied literals.
   *
   * Can be used in SATISFIABLE and UNKNOWN state.
   *
   * @see isZeroImplied()
   */
  virtual void collectZeroImplied(SATLiteralStack& acc) override;
  /**
   * Return a valid clause that contains the zero-implied literal
   * and possibly the assumptions that implied it. Return 0 if @c var
   * was an assumption itself.
   * If called on a proof producing solver, the clause will have
   * a proper proof history.
   */
  virtual SATClause* getZeroImpliedCertificate(unsigned var) override;

  void ensureVarCount(unsigned newVarCnt) override {
    while (_varCnt < newVarCnt) {
      newVar();
    }
  }


  unsigned newVar() override;

  // Currently not implemented for Z3
  virtual void suggestPolarity(unsigned var, unsigned pol) override {}
  virtual void setPolarity(unsigned var, unsigned pol) override {}

  virtual void addAssumption(SATLiteral lit) override;
  virtual void retractAllAssumptions() override;
  virtual bool hasAssumptions() const override { return !_assumptions.isEmpty(); }

  virtual Status solveUnderAssumptions(const SATLiteralStack& assumps, unsigned conflictCountLimit, bool onlyProperSubusets) override;

  /**
   * The set of inserted clauses may not be propositionally UNSAT
   * due to theory reasoning inside Z3.
   * We cannot later minimize this set with minisat.
   *
   * TODO: think of extracting true refutation from Z3 instead.
   */
  SATClauseList* getRefutationPremiseList() override{ return 0; }

  SATClause* getRefutation() override;

  template<class F>
  auto scoped(F f)  -> decltype(f())
  {
    _solver.push();
    auto result = f();
    _solver.pop();
    return result;
  }

  using FuncId = unsigned;
  using PredId = unsigned;
  using SortId = TermList;

  struct FuncOrPredId
  {
    explicit FuncOrPredId(unsigned id, bool isPredicate, Term *forSorts = nullptr) : id(id), isPredicate(isPredicate), forSorts(forSorts) {}
    explicit FuncOrPredId(Term* term) :
      FuncOrPredId(
        term->functor(),
        term->isLiteral(),
        term->numTypeArguments() == 0 ? nullptr : term
      )
    {}
    static FuncOrPredId monomorphicFunction(FuncId id) { return FuncOrPredId (id, false); }
    static FuncOrPredId monomorphicPredicate(PredId id) { return FuncOrPredId (id, true); }
    unsigned id;
    bool isPredicate;
    /**
     * polymorphic symbol application: treat e.g. f(<sorts>, ...) as f<sorts>(...) for Z3
     * in the monomorphic case, nullptr
     * in the polymorphic case, some term of the form f(<sorts>, ...) which we use only for sort information
     */
    Term *forSorts;

    friend struct std::hash<FuncOrPredId> ;
    friend bool operator==(FuncOrPredId const& l, FuncOrPredId const& r)
    {
      if(l.id != r.id || l.isPredicate != r.isPredicate)
        return false;
      if(!l.forSorts)
        return true;
      ASS(r.forSorts != nullptr);

      // compare sort arguments
      for(unsigned i = 0; i < l.forSorts->numTypeArguments(); i++)
        // sorts are perfectly shared
        if(!l.forSorts->typeArg(i).sameContent(r.forSorts->typeArg(i)))
          return false;

      return true;
    }
    friend std::ostream& operator<<(std::ostream& out, FuncOrPredId const& self)
    {
      out << (self.isPredicate ? "pred " : "func ");
      out << (
        self.isPredicate
          ? env.signature->getPredicate(self.id)->name()
          : env.signature->getFunction(self.id)->name()
      );
      if(self.forSorts)
        for(unsigned i = 0; i < self.forSorts->numTypeArguments(); i++)
          out << " " << self.forSorts->typeArg(i).toString();
      return out;
    }
  };

private:

  Map<SortId, z3::sort> _sorts;
  struct Z3Hash {
    static unsigned hash(z3::func_decl const& c) { return c.hash(); }
    static unsigned hash(z3::expr const& c) { return c.hash(); }
    static bool equals(z3::func_decl const& l, z3::func_decl const& r) { return z3::eq(l,r); }
    static bool equals(z3::expr const& l, z3::expr const& r) { return z3::eq(l,r); }
  };
  Map<z3::func_decl, FuncOrPredId , Z3Hash > _fromZ3;
  Map<FuncOrPredId,  z3::func_decl, StlHash> _toZ3;
  Set<SortId> _createdTermAlgebras;

  z3::func_decl const& findConstructor(Term* t);
  void createTermAlgebra(TermList sort);

  z3::sort getz3sort(SortId s);

  z3::func_decl z3Function(FuncOrPredId function);

  friend struct ToZ3Expr;
  friend struct EvaluateInModel;
public:
  Term* evaluateInModel(Term* trm);
#if VDEBUG
  z3::model& getModel() { return _model; }
#endif

private:

  struct Representation
  {
    Representation(z3::expr expr, Stack<z3::expr> defs) : expr(expr), defs(defs) {}
    Representation(Representation&&) = default;
    z3::expr expr;
    Stack<z3::expr> defs;
  };

  Representation getRepresentation(Term* trm);
  Representation getRepresentation(SATLiteral lit);
  Representation getRepresentation(SATClause* cl);

  // arrays are a bit fragile in Z3, so we need to do things differently for them
  bool _hasSeenArrays;

  unsigned _varCnt; // just to conform to the interface
  SAT2FO& _sat2fo; // Memory belongs to Splitter

  Status _status;
  z3::config _config;
  z3::context _context;
  z3::solver _solver;
  z3::model _model;
  Stack<z3::expr> _assumptions;
  BiMap<SATLiteral, z3::expr, DefaultHash, Z3Hash> _assumptionLookup;
  const bool _showZ3;
  const bool _unsatCore;
  Option<std::ofstream> _out;
  Map<unsigned, z3::expr> _varNames;
  Map<TermList, z3::expr> _termIndexedConstants;
  Map<Signature::Symbol*, z3::expr> _constantNames;

  bool     isNamedExpr(unsigned var) const;
  z3::expr getNameExpr(unsigned var);

  z3::expr getNamingConstantFor(TermList name, z3::sort sort);
  z3::expr getConst(Signature::Symbol* symb, z3::sort srt);

                                 void __output(std::ostream& out               ) {                                 }
  template<class A, class... As> void __output(std::ostream& out, A a, As... as) { out << a; __output(out, as...); }

  template<class... As>
  void _output(bool endl, As... as)
  {
    if (_out.isSome()) {
      __output(_out.unwrap(), as...);
      if (endl)
        _out.unwrap() << std::endl;
    }
  }

  template<class... As> void outputln(As... as) { _output(true , as...); }
  template<class... As> void output  (As... as) { _output(false, as...); }
};

}//end SAT namespace
namespace std {
    template<>
    struct hash<SAT::Z3Interfacing::FuncOrPredId> {
      size_t operator()(SAT::Z3Interfacing::FuncOrPredId const& self) {
        unsigned hash = Lib::HashUtils::combine(self.id, self.isPredicate);
        if(self.forSorts)
          for(unsigned i = 0; i < self.forSorts->numTypeArguments(); i++)
            hash = Lib::HashUtils::combine(hash, self.forSorts->typeArg(i).content());
        return hash;
      }
    };
}

#endif /* if VZ3 */
#endif /*Z3Interfacing*/
