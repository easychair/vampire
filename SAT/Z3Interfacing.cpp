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
 * @file Z3Interfacing.cpp
 * Implements class Z3Interfacing
 */

#include "Lib/Allocator.hpp"
#if VZ3
#define UNIMPLEMENTED ASSERTION_VIOLATION
#define MODEL_COMPLETION true

#include "Forwards.hpp"
#include "Kernel/Term.hpp"
#include "Lib/StringUtils.hpp"
#include "z3.h"

#include "SATSolver.hpp"
#include "SATLiteral.hpp"
#include "SATClause.hpp"
#include "SATInference.hpp"

#include "Lib/Environment.hpp"
#include "Lib/System.hpp"

#include "Kernel/NumTraits.hpp"
#include "Kernel/Signature.hpp"
#include "Kernel/OperatorType.hpp"
#include "Kernel/SortHelper.hpp"
#include "Kernel/SubstHelper.hpp"
#include "Kernel/BottomUpEvaluation.hpp"
#include "Lib/Coproduct.hpp"

#include "Shell/UIHelper.hpp"
#include "Indexing/TermSharing.hpp"
#include "Z3Interfacing.hpp"

#define DEBUG(...) // DBG(__VA_ARGS__)

#define TRACE_Z3 0
#define INSTANTIATE_EXPRESSIONS 0
#define ENABLE_Z3_PROOF_GENERATION 0

namespace Lib {
using SortId = TermList;

template<>
struct BottomUpChildIter<z3::expr>
{
  unsigned _idx;
  z3::expr _self;

  /** constructs an iterator over the children of the current node */
  BottomUpChildIter(z3::expr a, EmptyContext = EmptyContext()) : _idx(0), _self(a) {}

  /** returns the node this iterator was constructed with */
  z3::expr self(EmptyContext = EmptyContext()) { return _self; }

  /** returns the next child of the node this this object was constructed with */
  z3::expr next(EmptyContext = EmptyContext()) { return _self.arg(_idx++); }

  /** returns the next child of the current node in the structure to be traversed */
  bool hasNext(EmptyContext = EmptyContext()) { return _self.is_app() && _idx < _self.num_args(); }

  /** returns how many children this node has */
  unsigned nChildren(EmptyContext = EmptyContext()) { return _self.is_app() ? _self.num_args() : 0; }
};

} // namespace Lib


auto quotient0_name(char c, z3::sort s) 
{ return Output::toString("quot-0-", c, "-", s); }

auto remainder0_name(char c, z3::sort s) 
{ return Output::toString("rem-0-", c, "-", s); }

template<class UInt64ToExpr>
z3::expr int_to_z3_expr(IntegerConstantType const& val, UInt64ToExpr toExpr) {
    auto sign = val.sign();
    auto abs = val.abs();

    RStack<uint64_t> digits;
    /* this is how we translate arbitrary big IntegerConstantType numbers to  z3::expr numbers */
    z3::expr base = // <- == 2^64
      toExpr(std::numeric_limits<uint64_t>::max()) + toExpr(1);
    while(!abs.fits<unsigned long>()) {
      auto ui = abs.truncate<unsigned long>();
      using ui_t = decltype(ui);
      static_assert(sizeof(ui_t) == sizeof(uint64_t), "unexpected number size");
      static_assert(sizeof(ui_t) == 64 / 8, "unexpected number size");
      static_assert(std::numeric_limits<ui_t>::max() == std::numeric_limits<uint64_t>::max(), "unexpected number size");
      digits->push(uint64_t(ui));
      abs.rshiftBits(64);
    }
    z3::expr res = toExpr(uint64_t(abs.template truncate<unsigned int>()));
    while(digits->isNonEmpty()) {
      res = toExpr(digits->pop()) + (res * base);
    }

    return sign == Sign::Neg ? -res : res;
};

namespace SAT
{

using namespace std;
using namespace Shell;
using namespace Lib;
using ProblemExportSyntax = Shell::Options::ProblemExportSyntax;

//using namespace z3;

Z3Interfacing::Z3Interfacing(const Shell::Options& opts, SAT2FO& s2f, bool unsatCore, std::string const& exportFile,Shell::Options::ProblemExportSyntax exportSyntax):
  Z3Interfacing(s2f, opts.showZ3(), /* unsatCore */ unsatCore, exportFile, exportSyntax)
{ }

const char* errToString(Z3_error_code code)
{
  switch (code) {
    case Z3_OK: return "Z3_OK";
    case Z3_SORT_ERROR: return "Z3_SORT_ERROR";
    case Z3_IOB: return "Z3_IOB";
    case Z3_INVALID_ARG: return "Z3_INVALID_ARG";
    case Z3_PARSER_ERROR: return "Z3_PARSER_ERROR";
    case Z3_NO_PARSER: return "Z3_NO_PARSER";
    case Z3_INVALID_PATTERN: return "Z3_INVALID_PATTERN";
    case Z3_MEMOUT_FAIL: return "Z3_MEMOUT_FAIL";
    case Z3_FILE_ACCESS_ERROR: return "Z3_FILE_ACCESS_ERROR";
    case Z3_INTERNAL_FATAL: return "Z3_INTERNAL_FATAL";
    case Z3_INVALID_USAGE: return "Z3_INVALID_USAGE";
    case Z3_DEC_REF_ERROR: return "Z3_DEC_REF_ERROR";
    case Z3_EXCEPTION: return "Z3_EXCEPTION";
  }
  ASSERTION_VIOLATION; return "UNKNOWN ERROR";
}

struct Z3MkConstructorCall {
  Z3_context c;
  Z3_symbol name;
  Z3_symbol tester;
  Stack<Z3_symbol> field_names;
  Stack<Z3_sort> sorts;
  Stack<unsigned> sort_refs;

  unsigned arity() { return field_names.size(); }
  Z3_constructor operator()() {
    return Z3_mk_constructor(
        c,
        name,
        tester,
        arity(),
        field_names.begin(),
        sorts.begin(),
        sort_refs.begin()
    );
  }
};

struct Z3Constructor
{
  z3::func_decl func;
  z3::func_decl tester;
  Stack<z3::func_decl> args;
};

struct Z3Datatype
{
  z3::sort sort;
  Stack<Z3Constructor> ctors;
};

struct Z3MkDatatypesCall
{
  std::unique_ptr<z3::context>& _context;
  Stack<Z3_symbol> sortNames;               // <- needed for Z3_mk_datatypes(...)
  Stack<Stack<Z3MkConstructorCall>> mkConstrs;


  Z3MkDatatypesCall(decltype(_context) context, Stack<TermAlgebra*> const& tas)
    : _context(context)
    , sortNames(tas.size())
    , mkConstrs(tas.size())
  { }

  unsigned nDtys() { return sortNames.size(); }

  Stack<Z3Datatype> operator()(){
    Array<Z3_sort> sorts;                     // <- needed for Z3_mk_datatypes(...)

    /* re-arranging datat for Z3_mk_datatypes call */
    Stack<Z3_constructor_list> z3_ctor_lists(nDtys()); // <- needed for Z3_mk_datatypes(...)
    Stack<Stack<Z3_constructor>> ctorss;      // <- needed for Z3_query_constructor(..)
    for (auto& mks : mkConstrs) {
      Stack<Z3_constructor> ctors(mks.size());
      for (auto& mkConstr : mks) {
        ctors.push(mkConstr());
      }
      z3_ctor_lists.push(Z3_mk_constructor_list(*_context, ctors.size(),  ctors.begin()));
      ctorss.push(std::move(ctors));
    }

    Z3_mk_datatypes(*_context, nDtys(), sortNames.begin(), sorts.begin(), z3_ctor_lists.begin());

    /* querying result of Z3_mk_datatypes call */
    Stack<Z3Datatype> out(nDtys());

    for (unsigned i = 0; i < ctorss.size(); i++) {
      auto sort = z3::sort(*_context, sorts[i]);
      Stack<Z3Constructor> ctors_res(ctorss[i].size());
      for (unsigned j = 0; j < ctorss[i].size(); j++) {
        Z3_func_decl func;
        Z3_func_decl tester;
        DArray<Z3_func_decl> args(mkConstrs[i][j].arity());

        Z3_query_constructor(*_context, ctorss[i][j], mkConstrs[i][j].arity(), &func, &tester, args.begin());

        ctors_res.push(Z3Constructor {
            .func   = z3::func_decl(*_context, func),
            .tester = z3::func_decl(*_context, tester),
            .args   = arrayIter(args)
                         .map([&](auto arg) { return z3::func_decl(*_context, arg); })
                         .template collect<Stack>(),
          });
      }
      out.push(Z3Datatype { .sort = sort, .ctors = std::move(ctors_res), });
    }

    /* clean up */

    for (auto& lst : z3_ctor_lists) {
      Z3_del_constructor_list(*_context, lst);
    }

    for (auto& ctors : ctorss) {
      for (auto& ctor : ctors) {
        Z3_del_constructor(*_context, ctor);
      }
    }

    return out;
  }
};


void handleZ3Error(Z3_context ctxt, Z3_error_code code)
{
  DEBUG(errToString(code))
  throw z3::exception(errToString(code));
}

Z3Interfacing::Z3Interfacing(SAT2FO& s2f, bool showZ3, bool unsatCore, std::string const& exportSmtlib, Options::ProblemExportSyntax exportSyntax):
  _hasSeenArrays(false),
  _varCnt(0),
  _sat2fo(s2f),
  _outSyntax(exportSyntax),
  _status(Status::SATISFIABLE),
  _context(new z3::context()),
  _solver(*_context),
  _model((_solver.check(), _solver.get_model())),
  _assumptions([]() {
      if (ENABLE_Z3_PROOF_GENERATION) {
        // needs to be called before _context is initialized, therefore we call it in
        // a closure that must be evaluated before _context is initialized
        z3::set_param("proof", true);
      }
      return decltype(_assumptions)();
      }()),
  _showZ3(showZ3),
  _unsatCore(unsatCore),
  _exporter([&](){
      using namespace ProblemExport;
      if (exportSmtlib == "") {
        return decltype(_exporter)(NoExport{});
      } else {
        std::ofstream file(exportSmtlib.c_str());
        if (file.fail())
          throw UserErrorException("Failed to open file: ", exportSmtlib);
        switch (exportSyntax) {
        case Shell::Options::ProblemExportSyntax::SMTLIB:    return decltype(_exporter)(Smtlib  (std::move(file), *_context));
        case Shell::Options::ProblemExportSyntax::API_CALLS: return decltype(_exporter)(ApiCalls(std::move(file), *_context));
        }
        ASSERTION_VIOLATION
      }
    }())
{
  _out = exportSmtlib == "" ? Option<std::ofstream>()
                            : Option<std::ofstream>(std::ofstream(exportSmtlib.c_str())) ;
  if (_out.isSome() && _out.unwrap().fail()) {
    throw UserErrorException("Failed to open file: ", exportSmtlib);
  }

  _exporter.apply([&](auto& e) { e.initialize(); });

  z3_set_param("rewriter.expand_store_eq", true);
  z3_set_param("model.completion", MODEL_COMPLETION);
  z3_set_param("model.compact", true); // keeps z3 from compressing its model. ~50% of the runtime of get_model is spent doing that otherwise

  if (_unsatCore) {
    z3_set_param(":unsat-core", true);
  }
  // Z3_set_error_handler(*_context, handleZ3Error); // MS: a handled error only reveals Z3_error_code, a propragated z3::exception is typically more informative

#if TRACE_Z3
  z3_enable_trace("arith");
#endif // TRACE_Z3

  for (auto c : { 'f', 't' }) {
    for (auto s : { _context->real_sort(), _context->int_sort() }) {
      // we need these auxilary variables to make $quotient_t and friends
      // uninterpreted functions for a zero divisor. i.e. we need to make
      // sure that they are completely freely interpreted, and that there
      // is for example no relationship between $quotient_t(2, 0),
      // $remainder_t(2, 0). in previous definitions they functionally
      // deptendet on the result of 2/0, which is not sound.
      // We make sure that they are freely interpreted by introducing
      // an uninterpreted function $quotient_t0, and defining
      // $quotient_t(x, y) = if(y == 0) $quotient_t0(y)
      //                     else <actual definition >
      //
      // This is done later in ToZ3Expr
      z3::sort_vector dom(*_context);
      dom.push_back(s);
      z3_declare_fun(quotient0_name(c, s), dom, s);
      z3_declare_fun(remainder0_name(c,s), dom, s);
    }
  }
}

void ProblemExport::Smtlib::initialize()                               {                                                        }
void ProblemExport::Smtlib::terminate()                                { out <<                                      std::endl; }
void ProblemExport::Smtlib::declareSort(z3::sort sort)                 { out << "(declare-sort " << sort << " 0)" << std::endl; }
void ProblemExport::Smtlib::eval(z3::expr const& x)                    { out << "(get-value (" << x << "))"       << std::endl; }
void ProblemExport::Smtlib::unsatCore()                                { out << "(get-unsat-core)"                << std::endl; }
void ProblemExport::Smtlib::addAssert(z3::expr const& x)               { out << "(assert " << x << ")"            << std::endl; }
void ProblemExport::Smtlib::get_model()                                { out << "(get-model)"                     << std::endl; }
void ProblemExport::Smtlib::reset()                                    { out << "(reset)"                         << std::endl; }

void ProblemExport::Smtlib::declare_const(std::string const& name, z3::sort codomain)
{ return declare_fun(name, z3::sort_vector(codomain.ctx()), codomain); }

void ProblemExport::Smtlib::declare_fun(std::string const& name, z3::sort_vector domain, z3::sort codomain) {
  out << "(declare-fun " << name << " (";
  for (auto s : domain)
    out << " "  << s;
  out << " ) " << codomain << ")" << std::endl;
}

void ProblemExport::Smtlib::check(Stack<z3::expr> const& assumptions)  {
  out << "(check-sat-assuming (";
  for (auto const& a : assumptions)
    out << " " << a;
  out << " ))" << std::endl;
}

void ProblemExport::Smtlib::instantiate_expression(z3::expr const&) { }

void ProblemExport::Smtlib::declare_array_sort(z3::sort array, z3::sort index, z3::sort result) { }

template<class Value>
void ProblemExport::Smtlib::set_param(const char* k, Value const& v)
{ out << ";- setting z3 parameter: " << k << "=" << v << std::endl; }

void ProblemExport::Smtlib::Z3_mk_datatypes(Z3MkDatatypesCall const& call) {
  auto quote = [&](auto x){
    std::stringstream s;
    s << x;
    auto str = s.str();
    if (str[0] == '\'') {
      return "|" + str + "|";
    } else {
      return str;
    }
  };

  out << "(declare-datatypes (" << std::endl;
  for (auto& s : call.sortNames) {
    out << " (" << quote(z3::symbol(_ctxt, s)) << " 0)";
  }
  out << " ) (" << std::endl;

  for (unsigned i = 0; i < call.sortNames.size(); i++) {
    out << "    ( ;-- datatype " << z3::symbol(_ctxt, call.sortNames[i]) << std::endl;
    for (auto& ctor : call.mkConstrs[i]) {
      out << "        ( " << quote(z3::symbol(_ctxt, ctor.name));
      for (unsigned j = 0; j < ctor.field_names.size(); j++) {
        out << " ( " << quote(z3::symbol(_ctxt, ctor.field_names[j])) << " ";
        if (ctor.sorts[j] == nullptr) out << quote(z3::symbol(_ctxt, call.sortNames[ctor.sort_refs[j]]));
        else                          out << z3::sort(_ctxt, ctor.sorts[j]);
        out << " )";
      }
      out << " )" << std::endl;
    }
    out << "    )";
  }
  out << "))" << std::endl;
}

std::string ProblemExport::ApiCalls::escapeVarName(z3::sort const& sym)
{
  if (sym.is_array()) {
    // Array sorts have argments. Hence we need to escape the arguments as well, not only the sort name
    return _escapeVarName(sym);
  } else {
    return Z3_ast(sym) == nullptr 
      ? std::string("nullptr")
      : _escapeVarName(sym.name());
  }
}

std::string ProblemExport::ApiCalls::escapeVarName(z3::symbol const& sym)
{ return _escapeVarName(sym); }


template<class Outputable>
std::string ProblemExport::ApiCalls::_escapeVarName(Outputable const& sym) {
  std::stringstream cvar;
  auto generatePrefix = [&](std::string const& toEscape) -> std::string {
    unsigned iter = 0;
    while (iter < toEscape.length()) {
      if (std::isalnum(toEscape[iter]) || toEscape[iter] == '_') break;
      else iter++;
    }
    if (iter == toEscape.length()) {
      /* none alphanumeric or _ */
      cvar << "_";
    } else {
      if ('0' <= toEscape[iter] && toEscape[iter] <= '9')
        cvar << "_";

      while (iter != toEscape.length()) {
        // we replace every letter that is not alphanumeric by '_'
        if (std::isalnum(toEscape[iter]) || toEscape[iter] == '_') {
          cvar << toEscape[iter];
        } else {
          cvar << '_';
        }
        iter++;
      }
    }
    return std::string(cvar.str());
  };


  auto origName = Output::toString(sym);
  return _escapedNames.getOrInit(origName, [&](){
    auto& ids = _escapePrefixes.getOrInit(generatePrefix(origName));
    auto nextId = ids.size();
    auto id = ids.getOrInit(origName, [&](){ return nextId; });
    if (id != 0)
      cvar << "_" << id;

    // DBG(sym, " -> ", cvar, " -> ", cvar.str())
    return cvar.str();
  });
}

void ProblemExport::ApiCalls::enableTrace(const char* name)
{
  out << "Z3_enable_trace(\"" << name << "\");" << std::endl;
}

void ProblemExport::ApiCalls::instantiate_expression(z3::expr const& expr)
{
#if INSTANTIATE_EXPRESSIONS
  out << "  (void) " << serialize(expr) << ";" << std::endl;
#endif
}

void ProblemExport::ApiCalls::initialize() {
  out << R"(
#include <z3++.h>
#include <z3_api.h>
#include <iostream>
#include <vector>

int main() {)" << std::endl;
#if ENABLE_Z3_PROOF_GENERATION
    out << "  z3::set_param(\"proof\", true);" << std::endl;
#endif
  out << R"(
  z3::context ctx;
  z3::solver solver(ctx);
  z3::model  model(ctx);
  auto sort_vec = [&](std::initializer_list<z3::sort> xs) {
    z3::sort_vector vec(ctx);
    for (auto s : xs) vec.push_back(s);
    return vec;
  };
  auto expr_vec = [&](std::initializer_list<z3::expr> xs) {
    z3::expr_vector vec(ctx);
    for (auto s : xs) vec.push_back(std::move(s));
    return vec;
  };

  auto query_constructor = [&](Z3_constructor& ctor,
                               z3::func_decl* name,
                               z3::func_decl* tester,
                               std::vector<z3::func_decl*> accessors) {

    Z3_func_decl _name;
    Z3_func_decl _tester;
    std::vector<Z3_func_decl> _accessors;
    for (auto a : accessors)
      _accessors.push_back(Z3_func_decl(a));
    Z3_query_constructor(ctx, ctor, accessors.size(), &_name, &_tester, accessors.size() == 0 ? nullptr : &_accessors[0]);
    *name   = z3::func_decl(ctx, _name);
    *tester = z3::func_decl(ctx, _tester);
    for (auto i = 0; i < accessors.size(); i++) {
      *accessors[i] = z3::func_decl(ctx, _accessors[i]);
    }
  };

  auto mk_constructor = [&](Z3_symbol name,
                            Z3_symbol tester,
                            std::vector<Z3_symbol> argNames,
                            std::vector<Z3_sort>   sorts,
                            std::vector<unsigned>  sortRefs) {
    return Z3_mk_constructor(ctx, name, tester, argNames.size(),
                             argNames.size() == 0 ? nullptr : &argNames[0],
                                sorts.size() == 0 ? nullptr :    &sorts[0],
                             sortRefs.size() == 0 ? nullptr : &sortRefs[0]);
  };


  auto mk_constructor_list = [&](std::vector<Z3_constructor> ctors) {
    return Z3_mk_constructor_list(ctx, ctors.size(), ctors.size() == 0 ? nullptr : &ctors[0]);
  };

)" << std::endl;
  #define ADD_BUILTIN_SORT(name, Name) \
    out << "  z3::sort " << escapeVarName(_ctxt.str_symbol(Name))  \
        << " = ctx." << name << "_sort();" << std::endl;
    ADD_BUILTIN_SORT("bool", "Bool")
    ADD_BUILTIN_SORT("int", "Int")
    ADD_BUILTIN_SORT("real", "Real")
  #undef ADD_BUILTIN_SORT
  out << endl;
}

void ProblemExport::ApiCalls::declare_array_sort(z3::sort array, z3::sort index, z3::sort result)
{
  out << "  z3::sort " << escapeVarName(array)
      << " = ctx.array_sort("
      << escapeVarName(index) << ", "
      << escapeVarName(result) << ");" << std::endl;
}

void ProblemExport::ApiCalls::terminate()
{
  out << "} // int main();" << std::endl;
}

struct ProblemExport::ApiCalls::EscapeString {
  std::string s;
  EscapeString(std::string s) : s(s) {}
  EscapeString(z3::expr const& x) : EscapeString(Output::toString(x)) {}
  friend std::ostream& operator<<(std::ostream& out, EscapeString const& self)
  { return out << "R\"(" << self.s << ")\""; }// TODO mask occurences of )"
};

std::ostream& ProblemExport::operator<<(std::ostream& out, ProblemExport::ApiCalls::Serialize<std::string> const& self)
{ return out << ProblemExport::ApiCalls::EscapeString{self.inner}; }

std::ostream& ProblemExport::operator<<(std::ostream& out, ProblemExport::ApiCalls::Serialize<bool> const& self)
{ return out << ( self.inner ? "true" : "false" ); }

std::ostream& ProblemExport::operator<<(std::ostream& out, ProblemExport::ApiCalls::Serialize<z3::expr> const& self)
{
  auto& x = self.inner;
  auto& state = self.state;
  #define ARG(idx) state.serialize(x.arg(idx))
  auto vec_func = [&](auto f) -> std::ostream& {
    out << f << "(expr_vec({";
    if (x.num_args() > 0) {
      out << ARG(0);
      for (unsigned i = 1; i < x.num_args(); i++)
        out << ", " << ARG(i);
    }
    return out << "}))";
  };
  auto func = [&](auto f) -> std::ostream& {
    if (x.num_args() > 4)
      return vec_func(f);
    else {
      if (x.num_args() == 0 && state._predeclaredConstants.contains(f)) {
        return out << f;
      } else {
        out << f << "(";
        if (x.num_args() > 0) {
          out << ARG(0);
          for (unsigned i = 1; i < x.num_args(); i++)
            out << ", " << ARG(i);
        }
        return out << ")";
      }
    }
  };

  auto bin = [&](auto op) -> std::ostream& {
    ASS_EQ(x.num_args(), 2)
    return out << "(" << ARG(0) << " " << op << " " << ARG(1) << ")";
  };

       if (x.is_eq())       return bin("==");
  else if (x.is_and())      return bin("&&");
  else if (x.is_or())       return bin("||");
  else if (x.is_not())      return func("!");
  else if (x.is_ite())      return func("z3::ite");
  else if (x.is_distinct()) return vec_func("z3::distinct");
  else if (x.is_implies())  return func("z3::implies");
  else if (x.is_true())     return out << "ctx.bool_val(true)";
  else if (x.is_false())    return out << "ctx.bool_val(false)";
  else if (x.is_numeral())  return out << "ctx.int_val(\""  << x << "\")";
  else if (x.is_app()) {
    auto f = x.decl();
    if (f.name().kind() == Z3_STRING_SYMBOL) {
      if (f.name().str() == "/" ) return bin("/");
      if (f.name().str() == "*" ) return bin("*");
      if (f.name().str() == "+" ) return bin("+");
      if (f.name().str() == "-" ) return x.num_args() == 1 ? func("-") : bin("-");
      if (f.name().str() == "<=") return bin("<=");
      if (f.name().str() == "<" ) return bin("<" );
      if (f.name().str() == ">=") return bin(">=");
      if (f.name().str() == ">" ) return bin(">" );
    }
    return func(self.state.escapeVarName(f.name()));
  } else  {
    ASSERTION_VIOLATION
  }
  #undef ARG
  // return out << "ctx.parse_string(" << EscapeString(self.inner) << ")[0]";
}

template<class A>
std::ostream& ProblemExport::operator<<(std::ostream& out, ProblemExport::ApiCalls::Serialize<A> const& self)
{ return out << self.inner; }

std::ostream& ProblemExport::operator<<(std::ostream& out, ProblemExport::ApiCalls::Serialize<z3::symbol> const& self)
{
  if (self.inner.kind() == Z3_INT_SYMBOL) {
    return out << "ctx.int_symbol(" << self.inner.to_int() << ")";
  } else  {
    auto str = Output::toString(self.inner);
    return out << "ctx.str_symbol(" << ProblemExport::ApiCalls::EscapeString(str) << ")";
  }
}

void ProblemExport::ApiCalls::declareSort(z3::sort sort) {
  out << "  z3::sort " << escapeVarName(sort)
      << " = ctx.uninterpreted_sort(" <<  serialize(sort.name()) << ");" << std::endl;
}

void ProblemExport::ApiCalls::eval(z3::expr const& x)
{ out << "  std::cout << \"model.eval(" << serialize(x) << ") = \" << model.eval(" << serialize(x) << " , " << MODEL_COMPLETION << ") << std::endl;" << std::endl; }

void ProblemExport::ApiCalls::unsatCore()
{
  out << "  std::cout << \"===== start solver.unsat_core() ====\" << std::endl;" << std::endl
      << "  std::cout << solver.unsat_core()                      << std::endl;" << std::endl
      << "  std::cout << \"=====   end solver.unsat_core() ====\" << std::endl;" << std::endl;
}


void ProblemExport::ApiCalls::addAssert(z3::expr const& x)
{
  // out << "  /* " << x <<  " */" << std::endl;
  out << "  solver.add(" << serialize(x) << ");" << std::endl;
}

void ProblemExport::ApiCalls::check(Stack<z3::expr> const& xs)
{
  out << std::endl;
  out << std::endl << "  std::cout << \"solver.check(..) = \" << solver.check(expr_vec({";
  for (auto& x : xs) {
    out << serialize(x) << ", ";
  }
  out << "})) << std::endl;" << std::endl;
  out << std::endl;
}

void ProblemExport::ApiCalls::get_model() {
  out << std::endl;
  out << "  model = solver.get_model();" << std::endl;
  out << "  std::cout                               << std::endl;" << std::endl;
  out << "  std::cout << \"===== start model ====\" << std::endl;" << std::endl;
  out << "  std::cout << model                      << std::endl;" << std::endl;
  out << "  std::cout << \"=====   end model ====\" << std::endl;" << std::endl;
  out << "  std::cout                               << std::endl;" << std::endl;
  out << std::endl;
}
void ProblemExport::ApiCalls::reset()     { out << "  std::cout << solver.reset() << std::endl;"     << std::endl; }

template<class Value>
void ProblemExport::ApiCalls::set_param(const char* k, Value const& v)
{ out << "  solver.set(" << EscapeString{k} << "," << serialize(v) << ");" << std::endl; }

template<class A, class F> struct InitList { A const& inner; F output; };
template<class A, class F> InitList<A, F> initList(A const& a, F f) { return InitList<A,F> { a, std::move(f), }; }
template<class A, class F> std::ostream& operator<<(std::ostream& out, InitList<A, F> const& self)
{
  out << "{ ";
  for (auto& x : self.inner) {
    self.output(out, x);
    out << ", ";
  }
  return out << "}";
}

void ProblemExport::ApiCalls::Z3_mk_datatypes(Z3MkDatatypesCall const& call) {

  out << std::endl << "  // datatypes:";
  for (auto s : call.sortNames) {
    out << " " << z3::symbol(_ctxt, s);
  }
  out << std::endl;

  for (auto& cs : call.mkConstrs) {
    for (auto& c : cs)  {
      out << "  z3::func_decl " << escapeVarName(z3::symbol(_ctxt, c.name))   << "(ctx);" << std::endl;
      out << "  z3::func_decl " << escapeVarName(z3::symbol(_ctxt, c.tester)) << "(ctx);" << std::endl;
      for (auto f : c.field_names)
        out << "  z3::func_decl " << escapeVarName(z3::symbol(_ctxt, f)) << "(ctx);" << std::endl;
    }
  }

  for (auto s : call.sortNames)
    out << "  z3::sort " << escapeVarName(z3::symbol(_ctxt,s)) << "(ctx);" << std::endl;
  // for (auto s : call.sortNames)
  //   out << "  z3::sort " << escapeVarName(z3::symbol(_ctxt,s)) << "(ctx);" << std::endl;

  out << "  {" << std::endl
      << "    Z3_symbol sort_names[] = " << initList(call.sortNames, [&](auto& out, auto& s) { z3::symbol sym(_ctxt, s); out << serialize(sym); }) << ";" << std::endl
      << "    Z3_sort sorts[] = "        << initList(call.sortNames, [&](auto& out, auto&  ) { out << "nullptr"; }) << ";" << std::endl;

  auto ctor_name = [&](auto& c) { return "ctor_" + escapeVarName(z3::symbol(_ctxt, c.name)); };
  for (auto& cs : call.mkConstrs) {
    for (auto& c : cs) {
      z3::symbol name(_ctxt, c.name);
      z3::symbol tester(_ctxt, c.tester);
      out << "    auto " << ctor_name(c) << " = mk_constructor("
          << serialize(name) << ", "
          << serialize(tester) << ", "
          << initList(c.field_names, [&](auto& out, auto& f){ out << serialize(z3::symbol(_ctxt,f)); }) << ", "
          << initList(c.sorts,       [&](auto& out, auto& s){ out << escapeVarName(z3::sort(_ctxt,s)); }) << ", "
          << initList(c.sort_refs,   [&](auto& out, auto& s){ out << serialize(s); }) << ");" << std::endl;
    }
  }

  out << "    Z3_constructor_list constructor_lists[] = {" << std::endl;
  for (auto& cs : call.mkConstrs) {
    out << "      mk_constructor_list({";
    for (auto c : cs)
      out << ctor_name(c) << ", ";
    out << "      })," << std::endl;
  }
  out << "    };" << std::endl;

  out << "    Z3_mk_datatypes(ctx, " << call.sortNames.size() << ", sort_names, sorts, constructor_lists);" << std::endl;
  int i = 0;
  for (auto s : call.sortNames)
    out << "    " << escapeVarName(z3::symbol(_ctxt,s)) << " = z3::sort(ctx, sorts[" << i++ <<"]);" << std::endl;


  for (auto& cs : call.mkConstrs) {
    for (auto c : cs){
      out << "      query_constructor("
          << ctor_name(c) << ", "
          << "&" << escapeVarName(z3::symbol(_ctxt, c.name)) << ", "
          << "&" << escapeVarName(z3::symbol(_ctxt, c.tester)) << ", "
          << "{";
      for (auto f : c.field_names) {
        out << "&" << escapeVarName(z3::symbol(_ctxt, f)) << ", ";
      }
      out << "});" << std::endl;
    }
  }

  // TODO z3::func_decl for ctors
  out << "  }" << std::endl << std::endl;
}

void ProblemExport::ApiCalls::declare_fun(std::string const& name, z3::sort_vector domain, z3::sort codomain) {
  out << "  z3::func_decl " << escapeVarName(_ctxt.str_symbol(name.c_str())) << " = ctx.function(" << EscapeString{name} << ", sort_vec({";
  for (auto s : domain)
    out << escapeVarName(s) << ", ";
  out << "}), " << escapeVarName(codomain) << " );" << std::endl;
}

void ProblemExport::ApiCalls::declare_const(std::string const& name, z3::sort codomain) {
  auto varName = escapeVarName(_ctxt.str_symbol(name.c_str()));
  out << "  z3::expr " << varName
      << " = ctx.constant(" << EscapeString{name} << ", " << escapeVarName(codomain) << " );" << std::endl;
  _predeclaredConstants.insert(std::move(varName));
}

z3::sort Z3Interfacing::z3_array_sort(z3::sort const& index_sort, z3::sort const& value_sort)
{
  auto z3_sort = _context->array_sort(index_sort,value_sort);
  _exporter.apply([&](auto& e) { e.declare_array_sort(z3_sort, index_sort, value_sort); });
  return z3_sort;
}

void Z3Interfacing::z3_enable_trace(const char* name) {
  Z3_enable_trace(name);
  _exporter.apply([&](auto& e) { e.enableTrace(name); });
}


z3::sort Z3Interfacing::z3_declare_sort(std::string const& name) {
  auto sort = _context->uninterpreted_sort(_context->str_symbol(name.c_str()));
  _exporter.apply([&](auto& e) { e.declareSort(sort); });
  return sort;
}

z3::expr Z3Interfacing::z3_eval(z3::expr const& x) {
  _exporter.apply([&](auto& e) { e.eval(x); });
  return _model.eval(x, MODEL_COMPLETION);
}

z3::expr_vector Z3Interfacing::z3_unsat_core() {
  _exporter.apply([&](auto& e) { e.unsatCore(); });
  return _solver.unsat_core();
}

void Z3Interfacing::z3_add(z3::expr const& x) {
  _exporter.apply([&](auto& e) { e.addAssert(x); });
  _solver.add(x);
}

z3::check_result Z3Interfacing::z3_check() {
  _exporter.apply([&](auto& e) { e.check(_assumptions); });
  return _solver.check(_assumptions.size(), _assumptions.begin());
}

z3::model Z3Interfacing::z3_get_model() {
  _exporter.apply([&](auto& e) { e.get_model(); });
  return _solver.get_model();
}

// void Z3Interfacing::z3_reset() {
//   _exporter.apply([&](auto& e) { e.reset(); });
//   _solver.reset();
// }

z3::expr Z3Interfacing::z3_declare_const(std::string const& name, z3::sort sort) {
  _exporter.apply([&](auto& e) { e.declare_const(name, sort); });
  return _context->function(name.c_str(), z3::sort_vector(*_context), sort)();
}


z3::func_decl Z3Interfacing::z3_declare_fun(std::string const& name, z3::sort_vector domain, z3::sort codomain) {
  _exporter.apply([&](auto& e) { e.declare_fun(name, domain, codomain); });
  return _context->function(name.c_str(), domain, codomain);
}

template<class Value>
void Z3Interfacing::z3_set_param(const char* k, Value const& v)
{
  _exporter.apply([&](auto& e) { e.set_param(k, v); });
  _solver.set(k, v);
}

char const* Z3Interfacing::z3_full_version()
{
  return Z3_get_full_version();
}

unsigned Z3Interfacing::newVar()
{
  ++_varCnt;

  // to make sure all the literals we will ask about later have allocated counterparts internally
  auto rep = getRepresentation(SATLiteral(_varCnt,1));
  _exporter.apply([&](auto& exp){ exp.instantiate_expression(rep.expr); });

  return _varCnt;
}

void Z3Interfacing::addClause(SATClause* cl)
{
  ASS(cl);

  // store to later generate the refutation
  PrimitiveProofRecordingSATSolver::addClause(cl);

  auto z3clause = getRepresentation(cl);

  if(_showZ3){
    std::cout << "[Z3] add (clause): " << z3clause.expr << std::endl;
  }

  for (auto def : z3clause.defs)  {
    DEBUG("adding def: ", def)
    z3_add(def);
  }

  z3_add(z3clause.expr);
  DEBUG("adding expr: ", z3clause.expr)
}

void Z3Interfacing::retractAllAssumptions()
{
  _assumptionLookup.clear();
  _assumptions.truncate(0);
}

void Z3Interfacing::addAssumption(SATLiteral lit)
{
  auto pushAssumption = [&](SATLiteral lit) -> z3::expr
  {
    auto repr = getRepresentation(lit);
    for (auto& def : repr.defs)
      _assumptions.push(def);

    _assumptions.push(repr.expr);
    return repr.expr;
  };

  if (_unsatCore) {
    _assumptionLookup.getOrInit(lit, [&]() { return pushAssumption(lit); });
  } else {
    pushAssumption(lit);
  }
}

Z3Interfacing::Representation Z3Interfacing::getRepresentation(SATClause* cl)
{

  z3::expr z3clause = _context->bool_val(false);

  Stack<z3::expr> defs;

  unsigned clen=cl->length();
  for(unsigned i=0;i<clen;i++){
    SATLiteral l = (*cl)[i];
    auto repr = getRepresentation(l);
    _exporter.apply([&](auto& exp){ exp.instantiate_expression(repr.expr); });

    defs.loadFromIterator(repr.defs.iterFifo());

    z3clause = z3clause || repr.expr;
  }

  return Representation(std::move(z3clause), std::move(defs));
}

SATSolver::Status Z3Interfacing::solve()
{
  DEBUG("assumptions: ", _assumptions);

  /* The purpose of this class is to conditionally disable variable elimination inside Z3's _solver.check,
   * which results in some literals not being evaluated to either true and false, that we need for AVATAR.
   * Why a class? To be able to rely on RAII for the call to pop() (via the destructor) and thus not forget about it.
   * Why conditional? Because push/pop slightly decreases z3's performance and so we want to do it only in
   * the cases where the problem has been observed - namely, when arrays are involved.
  */
  class ScopedPushAndPop {
    z3::solver& _s;
    bool _dpp;
  public:
    ScopedPushAndPop(z3::solver& s, bool doPushPop) : _s(s), _dpp(doPushPop) { if (_dpp) {_s.push();} }
    ~ScopedPushAndPop() { if (_dpp) {_s.pop();} }
  } _maybePushAndPop(_solver,_hasSeenArrays);


  auto result = z3_check();

  if(_showZ3){
    std::cout << "[Z3] solve result: " << result << std::endl;
  }

  if (_unsatCore) {
    auto core = z3_unsat_core();
    for (auto phi : core) {
      _assumptionLookup
             .tryGet(phi)
             .andThen([this](SATLiteral l)
                 { _failedAssumptionBuffer.push(l); });
    }
  }

  switch (result) {
    case z3::check_result::unsat:
      _status = Status::UNSATISFIABLE;
      break;
    case z3::check_result::sat:
      _status = Status::SATISFIABLE;
      _model = z3_get_model();
      break;
    case z3::check_result::unknown:
      _status = Status::UNKNOWN;
      break;
    default: ASSERTION_VIOLATION;
  }

  return _status;
}

SATSolver::Status Z3Interfacing::solveUnderAssumptions(const SATLiteralStack& assumps, unsigned conflictCountLimit)
{
  if (!_unsatCore) {
    return SATSolverWithAssumptions::solveUnderAssumptions(assumps,conflictCountLimit);
  }

  ASS(!hasAssumptions());

  for (auto a: assumps) {
    addAssumption(a);
  }
  auto result = solve();
  retractAllAssumptions();
  return result;
}

SATSolver::VarAssignment Z3Interfacing::getAssignment(unsigned var)
{
  ASS_EQ(_status,Status::SATISFIABLE);
  z3::expr rep = isNamedExpr(var) ? getNameExpr(var) : getRepresentation(SATLiteral(var,1)).expr;
  _exporter.apply([&](auto& exp){ exp.instantiate_expression(rep); });
  ASS(isNamedExpr(var) || getRepresentation(SATLiteral(var,1)).defs.isEmpty())
  auto assignment = z3_eval(rep);

  if(assignment.bool_value()==Z3_L_TRUE){
    return VarAssignment::TRUE;
  } else if(assignment.bool_value()==Z3_L_FALSE){
    return VarAssignment::FALSE;
  } else {
#if VDEBUG
    std::cout << std::endl;
    std::cout << "===== start _model ====" << std::endl;
    std::cout << _model << std::endl;
    std::cout << "=====   end _model ====" << std::endl;
    std::cout << std::endl;
    std::cout << rep << std::endl;
    ASSERTION_VIOLATION_REP(assignment);
#endif
    return VarAssignment::NOT_KNOWN;
  }
}

OperatorType* operatorType(Z3Interfacing::FuncOrPredId f)
{
  return f.isPredicate
    ? env.signature->getPredicate(f.id)->predType()
    : env.signature->getFunction (f.id)->fnType();
}


// TODO does this correctly work with polymorphism?
Term* createTermOrPred(Z3Interfacing::FuncOrPredId f, unsigned arity, TermList* ts)
{
  return f.isPredicate
    ? Literal::create(f.id, arity, true, ts)
    : Term::create(f.id, arity, ts);
}

Term* Z3Interfacing::evaluateInModel(Term* trm)
{
  DEBUG("in: ", *trm)
  DEBUG("model: \n", _model)
  ASS(!trm->isLiteral());

  auto ev = z3_eval(getRepresentation(trm).expr);
  ASS(getRepresentation(trm).defs.isEmpty())
  SortId sort = SortHelper::getResultSort(trm);

  DEBUG("z3 expr: ", ev)
  using Copro = Coproduct<Term*, RationalConstantType, IntegerConstantType>;
  using Result = Option<Copro>;
  auto result = BottomUpEvaluation<z3::expr, Option<Copro>>()
    .function(
      [&](z3::expr const& expr, Result* evaluatedArgs) -> Result
      {
        auto toTerm = [&](Copro const& co, SortId sort) -> Term* {
          return co.match(
                  [&](Term* t)
                  { return t; },

                  [&](RationalConstantType c)
                  {
                    return sort == RealTraits::sort()
                      ? theory->representConstant(RealConstantType(c))
                      : theory->representConstant(c);
                  },

                  [&](IntegerConstantType c)
                  { return theory->representConstant(c); }
                  );
        };


        DEBUG("in: ", expr)
        auto intVal = [](z3::expr e) -> Option<int> {
          int val;
          return e.is_numeral_i(val)
            ? Option<int>(val)
            : Option<int>();
        };

        if (expr.is_int()) {
          return intVal(expr)
            .map([](int i) { return Copro(IntTraits::constantT(i)); });

        } else if(expr.is_real()) {
          if (!expr.is_numeral()) {
            // non-numeral reals are, e.g., the algebraic numbers such as (root-obj (+ (^ x 2) (- 2)) 2),
            // which we currently cannot handle
            return Result();
          }      

          auto toFrac = [&](int l, int r)  { return Copro(RatTraits::constant(l,r)); };

          auto nonFractional = intVal(expr).map([&](int i) { return toFrac(i,1); });
          if (nonFractional.isSome()) {
            return nonFractional;
          } else {
            auto num = intVal(expr.numerator());
            auto den = intVal(expr.denominator());
            if (num.isSome() && den.isSome()) {
              return Result(Copro(toFrac(num.unwrap(), den.unwrap())));
            } else {
              return Result();
            }
          }

        } else if (expr.is_app()) {
          auto f = expr.decl();
          auto vfunc = _fromZ3.get(f);
          unsigned arity = f.arity();
          ASS(arity == 0 || evaluatedArgs != nullptr)
          Stack<TermList> args(arity);
          for (unsigned i = 0; i < arity; i++) {
            if (evaluatedArgs[i].isNone()) {
              // evaluation failed somewhere in a recursive call
              return Result();
            } else {
              auto argSort = operatorType(vfunc)->arg(i);
              auto t = TermList(toTerm(evaluatedArgs[i].unwrap(), argSort));
              args.push(t);
            }
          }
          return Result(Copro(createTermOrPred(vfunc, args.size(), args.begin())));
        } else {
          return Result();
        }
      })
    .apply(ev)
    .map([&](Copro co) {
        return co.match(
            [&](Term* t)
            { return t; },

            [&](RationalConstantType c)
            {
              return sort == RealTraits::sort()
                ? theory->representConstant(RealConstantType(c))
                : theory->representConstant(c);
            },

            [&](IntegerConstantType c)
            { return theory->representConstant(c); }
            );
      })
    .unwrapOrElse([](){ return nullptr; });
  DEBUG("vampire expr: ", ev)
  return result;

}

bool Z3Interfacing::isZeroImplied(unsigned var)
{
  // Safe. TODO consider getting zero-implied
  return false;
}

void Z3Interfacing::collectZeroImplied(SATLiteralStack& acc)
{
  NOT_IMPLEMENTED;
}

SATClause* Z3Interfacing::getZeroImpliedCertificate(unsigned)
{
  NOT_IMPLEMENTED;

  return 0;
}

z3::sort Z3Interfacing::getz3sort(SortId s)
{
  auto srt = _sorts.tryGet(s);
  if (srt.isSome()) {
    return srt.unwrap();
  } else {
    auto insert = [&](z3::sort x) { _sorts.insert(s, x); };
    // TODO what about built-in tuples?

    // Deal with known sorts differently
         if(s == AtomicSort::boolSort()) insert(_context->bool_sort());
    else if(s ==  IntTraits::sort()) insert( _context->int_sort());
    else if(s == RealTraits::sort()) insert(_context->real_sort());
    else if(s ==  RatTraits::sort()) insert(_context->real_sort()); // Drops notion of rationality
    // TODO: are we really allowed to do this ???                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    else if(s.isArraySort()) {
      _hasSeenArrays = true;
      insert(z3_array_sort(
            getz3sort(SortHelper::getIndexSort(s)),
            getz3sort(SortHelper::getInnerSort(s))
            ));

    } else if (env.signature->isTermAlgebraSort(s)) {
      createTermAlgebra(s);

    } else {
      insert(z3_declare_sort(s.toString()));
    }
  }
  return _sorts.get(s);
}

template<class A>
std::string to_string(A const& a)
{
  std::stringstream out;
  out << a;
  return out.str();
}

void Z3Interfacing::createTermAlgebra(TermList sort)
{
  ASS(sort.isTerm() && sort.term()->isSort());
  if (_createdTermAlgebras.contains(sort)) return;

  Stack<TermList> taSorts;        // <- stack of term algebra sorts
  Map<SortId, unsigned> recSorts; // <- mapping term algeba -> index

  auto subsorts = TermAlgebra::subSorts(sort);
  for (auto s : subsorts.iter()) {
    if (env.signature->isTermAlgebraSort(s)
        && !_createdTermAlgebras.contains(s)) {
      auto idx = taSorts.size();
      taSorts.push(s);
      recSorts.insert(s, idx);
    }
  }

  auto new_string_symbol = [this](unsigned f, const std::string& typePostfix, bool pred = false)
  { return Z3_mk_string_symbol(*_context,
    ((pred ? env.signature->getPredicate(f) : env.signature->getFunction(f))->name()+'_'+typePostfix).c_str()); };

  // create the data needed for Z3_mk_datatypes(...)
  Stack<Stack<Z3_constructor>> ctorss(taSorts.size());
  Stack<Z3_constructor_list> ctorss_z3(taSorts.size());
  Stack<Z3_symbol> sortNames(taSorts.size());
  // create the data needed to serialize the declare-datatypes to smtlib
  struct SerDtor { z3::symbol name; SortId sort; };
  struct SerCtor { z3::symbol name; Stack<SerDtor> dtors; };
  struct SerDtype { SortId name; Stack<SerCtor> ctors; };
  Stack<SerDtype> toSerialize;

  DEBUG("creating constructors: ");
  for (auto taSort : taSorts) {
    _createdTermAlgebras.insert(taSort);
    ASS(taSort.isTerm());
    auto taSortT = taSort.term();
    auto ta = env.signature->getTermAlgebraOfSort(taSort);
    Substitution typeSubst;
    std::string typePostfix = "$";
    ta->getTypeSub(taSortT, typeSubst);
    for(unsigned i = 0; i < taSortT->arity(); i++) {
      // MS: Is it OK not to use any separator here?
      typePostfix += taSortT->nthArgument(i)->toString();
    }
    Stack<Z3_constructor> ctors(ta->nConstructors());
    Stack<SerCtor> serCtors;

    for (auto cons : ta->iterCons()) {

      // data needed for the  Z3_mk_constructor call
      Stack<SerDtor> serDtors;
      auto consTermArity = cons->arity()-cons->numTypeArguments();
      Stack<Z3_sort> argSorts(consTermArity);
      Stack<unsigned> argSortRefs(consTermArity);
      Stack<Z3_symbol> argNames(consTermArity);

      for (unsigned i = cons->numTypeArguments(); i < cons->arity(); i++) {
        auto argSort = SubstHelper::apply(cons->argSort(i),typeSubst);
        auto dtorName = new_string_symbol(cons->destructorFunctor(i-cons->numTypeArguments()),typePostfix);
        if (_out.isSome())
          serDtors.push(SerDtor {
              .name = z3::symbol(*_context, dtorName),
              .sort = argSort,
          });
        argNames.push(dtorName);
        recSorts.tryGet(argSort)
          .match([&](unsigned idx) {
                // for sorts that are to be generated with the call of Z3_mk_datatypes we need to push their index, and a nullptr
                argSortRefs.push(idx);
                argSorts.push(nullptr);
              },
              [&]() {
                // for other sorts, we need to push the sort, and an arbitrary index
                argSortRefs.push(0);  // <- 0 will never be read
                argSorts.push(getz3sort(argSort));
              });
      }

      DEBUG("\t", taSort.toString(), "::", env.signature->getFunction(cons->functor())->name(), ": ");//, SubstHelper::apply(cons->sort(), typeSubst));

      Z3_symbol ctorName = Z3_mk_string_symbol(*_context, env.signature->getFunction(cons->functor())->name().c_str());

      ASS_EQ(argSortRefs.size(), consTermArity);
      ASS_EQ(   argSorts.size(), consTermArity);
      ASS_EQ(   argNames.size(), consTermArity);
      if (_out.isSome())
        serCtors.push(SerCtor{
            .name = z3::symbol(*_context, ctorName),
            .dtors = std::move(serDtors),
        });
      ctors.push(Z3_mk_constructor(*_context,
          ctorName,
          new_string_symbol(cons->discriminator(), typePostfix, true/*predicate*/),
          consTermArity,
          consTermArity == 0 ? nullptr : argNames.begin(),
          consTermArity == 0 ? nullptr : argSorts.begin(),
          consTermArity == 0 ? nullptr : argSortRefs.begin()
      ));

    }
    ASS_EQ(ctors.size(), ta->nConstructors());

    ctorss.push(std::move(ctors));
    ASS_EQ(ctorss.top().size(), ta->nConstructors());
    ctorss_z3.push(Z3_mk_constructor_list(*_context, ctorss.top().size(),  ctorss.top().begin()));
    sortNames.push(Z3_mk_string_symbol(*_context, taSort.toString().c_str()));
    if (_out.isSome())
      toSerialize.push(SerDtype {
        .name = taSort,
        .ctors = std::move(serCtors),
      });
    }
  //   mkDatatypes.sortNames.push(new_string_symbol(ta->sort().toString()));

  ASS_EQ(sortNames.size(), taSorts.size())
  ASS_EQ(ctorss.size()   , taSorts.size())
  ASS_EQ(ctorss_z3.size(), taSorts.size())

  Array<Z3_sort> sorts(taSorts.size());

  // actually created the datatypes
  Z3_mk_datatypes(*_context, taSorts.size(), sortNames.begin(), sorts.begin(), ctorss_z3.begin());

  // register the `z3::func_decl`s created by `Z3_mk_datatypes` in indices to be queried when needed
  for (unsigned iSort = 0; iSort < sorts.size(); iSort++) {
    _sorts.insert(taSorts[iSort], z3::sort(*_context, sorts[iSort]));
    auto ta = env.signature->getTermAlgebraOfSort(taSorts[iSort]);
    auto& ctors = ctorss[iSort];
    for (unsigned iCons = 0; iCons < ta->nConstructors(); iCons++) {
      auto ctor_v = ta->constructor(iCons);
      auto ctorTermArity = ctor_v->arity() - ctor_v->numTypeArguments();

      Z3_func_decl ctor_z3;
      Z3_func_decl discr_z3;
      Array<Z3_func_decl> destr(ctorTermArity);

      Z3_query_constructor(*_context,
                           ctors[iCons],
                           ctorTermArity,
                           &ctor_z3,
                           &discr_z3,
                           destr.begin());

      auto registerDecl = [&](auto vampireId, auto z3Id) {
        _toZ3.insert(vampireId, z3Id);
        _fromZ3.insert(z3Id, vampireId);
      };
      registerDecl(FuncOrPredId::monomorphicFunction(ctor_v->functor()), z3::func_decl(*_context, ctor_z3));
      if (ctor_v->hasDiscriminator()) {
        registerDecl(FuncOrPredId::monomorphicPredicate(ctor_v->discriminator()), z3::func_decl(*_context, discr_z3));
      }
      for (unsigned iDestr = 0; iDestr < ctorTermArity; iDestr++)  {
        auto dtor = z3::func_decl(*_context, destr[iDestr]);

        // careful: datatypes can have boolean fields!
        // NB: polymorphic destructors instantiated with booleans are still functions (?)
        bool is_predicate = ctor_v->argSort(iDestr).isBoolSort();
        auto id = FuncOrPredId(ctor_v->destructorFunctor(iDestr), is_predicate);

        registerDecl(id, dtor);
      }
    }

  }
}

z3::func_decl const& Z3Interfacing::findConstructor(Term* t)
{
  auto id = FuncOrPredId::monomorphicFunction(t->functor());
  auto f = _toZ3.tryGet(id);
  if (f.isSome()) {
    return f.unwrap();
  } else {
    createTermAlgebra(SortHelper::getResultSort(t));
    return _toZ3.get(id);
  }
}


z3::expr to_int(z3::expr e)
{ return z3::expr(e.ctx(), Z3_mk_real2int(e.ctx(), e)); }

namespace tptp {

  z3::expr floor(z3::expr e)
  { return to_real(to_int(e)); }

  z3::expr ceiling(z3::expr x)
  { return -tptp::floor(-x); }

  z3::expr truncate(z3::expr x)
  { return ite(x >= 0, tptp::floor(x), tptp::ceiling(x)); }

  z3::expr quotient0(char kind, z3::expr x)
  {
    std::string fname = quotient0_name(kind, x.get_sort());
      // uninterpreted remainder for zero division
      auto quotient0 = x.ctx().function(fname.c_str(), x.get_sort(), x.get_sort());
      return quotient0(x);
  }

  z3::expr remainder0(char kind, z3::expr x)
  {
    std::string fname = remainder0_name(kind, x.get_sort());
      // uninterpreted remainder for zero division
      auto remainder0 = x.ctx().function(fname.c_str(), x.get_sort(), x.get_sort());
      return remainder0(x);
  }

  z3::expr quotient_e(z3::expr l, z3::expr r)
  { return l / r; }

  z3::expr remainder_e(z3::expr l, z3::expr r)
  { return z3::mod(l, r); }

  z3::expr quotient_t(z3::expr l, z3::expr r)
  { return ite(r == 0, tptp::quotient0('t', r)
                     , tptp::truncate(l / r)); }

  z3::expr quotient_f(z3::expr l, z3::expr r)
  { return ite(r == 0, tptp::quotient0('f', l / r)
                     , tptp::floor(l / r)); }

  template<class F>
  struct LiftInt
  {
    F bin_real_func;

    z3::expr operator()(z3::expr l, z3::expr r)
    { return to_int(bin_real_func(to_real(l), to_real(r))); }
  };
  template<class F> LiftInt<F> liftInt(F f) { return LiftInt<F>{ f }; }

  template<class F>
  struct RemainderOp
  {
    char kind;
    F quotient;

    z3::expr operator()(z3::expr l, z3::expr r)
    { return ite(r == 0, remainder0(kind, l)
                       , l - r * quotient(l,r)); }
  };
  template<class F> RemainderOp<F> remainder(char kind, F f) { return RemainderOp<F>{ kind, f }; }
}


z3::func_decl Z3Interfacing::z3Function(FuncOrPredId functor)
{
  auto& self = *this;

  auto found = self._toZ3.tryGet(functor);
  if (found.isSome()) {
    return found.unwrap();
  } else {
    // function does not yet exist, create it
    auto symb = functor.isPredicate ? env.signature->getPredicate(functor.id)
                                    : env.signature->getFunction(functor.id);
    auto type = functor.isPredicate ? symb->predType() : symb->fnType();

    // polymorphic symbol application: treat f(<sorts>, ...) as f<sorts>(...) for Z3
    std::string namebuf = symb->name();
    Substitution typeSubst;
    if(functor.forSorts) {
      SortHelper::getTypeSub(functor.forSorts, typeSubst);
      namebuf += '$';
      for(unsigned i = 0; i < functor.forSorts->numTypeArguments(); i++)
        namebuf += functor.forSorts->typeArg(i).toString();
    }

    z3::sort_vector domain_sorts = z3::sort_vector(*self._context);
    for (unsigned i=type->numTypeArguments(); i<type->arity(); i++) {
      TermList arg = SubstHelper::apply(type->arg(i), typeSubst);
      domain_sorts.push_back(self.getz3sort(arg));
    }
    auto codomain = functor.isPredicate ? self._context->bool_sort() : self.getz3sort(type->result());
    auto decl = self.z3_declare_fun(namebuf, domain_sorts, codomain);
    self._toZ3.insert(functor, decl);
    return decl;
  }
}

/**
 * Translate a Vampire term into a Z3 term
 * - Assumes term is ground
 * - Translates the ground structure
 * - Some interpreted functions/predicates are handled
 */
Z3Interfacing::Representation Z3Interfacing::getRepresentation(Term* trm)
{
  Stack<z3::expr> defs;
  auto expr = BottomUpEvaluation<TermList, z3::expr>()
    .function(
      [&](TermList toEval, z3::expr* args) -> z3::expr
      {
        ASS(toEval.isTerm())
        auto trm = toEval.term();
        bool isLit = trm->isLiteral();

        Signature::Symbol* symb;
        SortId range_sort;
        if (isLit) {
          symb = env.signature->getPredicate(trm->functor());
          range_sort = AtomicSort::boolSort();
          // check for equality
          if( trm->functor()==0 || symb->equalityProxy()){
            ASS(trm->numTermArguments()==2);
            // both equality and equality proxy translated as z3 equality
            if (symb->wasFlipped()) {
              // equality proxy could have been flipped (by random_polarities)
              return args[0] != args[1];
            } else {
              return args[0] == args[1];
            }
          }
        } else {
          auto actualSort = SortHelper::getResultSort(trm);
          symb = env.signature->getFunction(trm->functor());
          OperatorType* ftype = symb->fnType();
          range_sort = ftype->result();
          if (env.signature->isTermAlgebraSort(actualSort) &&  !_createdTermAlgebras.contains(actualSort) ) {
            createTermAlgebra(actualSort);
          }
        }

        //if constant treat specially
        if(trm->arity() == 0) {
          if(symb->integerConstant()){
            return int_to_z3_expr(symb->integerValue(), [&](uint64_t i) { return _context->int_val(i); });
          }
          if(symb->realConstant()) {
            RealConstantType value = symb->realValue();
            auto num = int_to_z3_expr(value.numerator()  , [&](uint64_t i) { return _context->real_val(i); });
            auto den = int_to_z3_expr(value.denominator(), [&](uint64_t i) { return _context->real_val(i); });
            return num / den;
          }
          if(symb->rationalConstant()) {
            RationalConstantType value = symb->rationalValue();
            auto num = int_to_z3_expr(value.numerator()  , [&](uint64_t i) { return _context->real_val(i); });
            auto den = int_to_z3_expr(value.denominator(), [&](uint64_t i) { return _context->real_val(i); });
            return num / den;
          }
          if(!isLit && env.signature->isFoolConstantSymbol(true,trm->functor())) {
            return _context->bool_val(true);
          }
          if(!isLit && env.signature->isFoolConstantSymbol(false,trm->functor())) {
            return _context->bool_val(false);
          }
          if(symb->termAlgebraCons()) {
            auto ctor = findConstructor(trm);
            return ctor();
          }
          // TODO do we really have overflownConstants ?? not in evaluation(s) at least
          if (symb->overflownConstant()) {
            // too large for native representation, but z3 should cope
            auto s = symb->fnType()->result();
            if (s == IntTraits::sort()) {
              return _context->int_val(symb->name().c_str());
            } else if (s == RatTraits::sort()) {
              return _context->real_val(symb->name().c_str());
            } else if (s == RealTraits::sort()) {
              return _context->real_val(symb->name().c_str());
            } else {
              ; // intentional fallthrough; the input is fof (and not tff), so let's just treat this as a constant
            }
          }

          // If not value then create constant symbol
          return getConst(symb, getz3sort(range_sort));
        }
        ASS_G(trm->arity(), 0);

        // Currently do not deal with all intepreted operations, should extend
        // - constants dealt with above
        // - unary funs/preds like is_rat interpretation unclear
        if(symb->interpreted()){
          Interpretation interp = static_cast<Signature::InterpretedSymbol*>(symb)->getInterpretation();

          if (Theory::isPolymorphic(interp)) {
            switch(interp){
              case Theory::ARRAY_SELECT:
              case Theory::ARRAY_BOOL_SELECT:
                // select(array,index)
                return select(args[0],args[1]);

              case Theory::ARRAY_STORE:
                // store(array,index,value)
                return store(args[0],args[1],args[2]);

              default:
                {}//skip it and treat the function as uninterpretted
            }

          } else {
            auto int_zero = _context->int_val(0);
            auto real_zero = _context->real_val(0);

            switch(interp){
            // Numerical operations
            case Theory::INT_DIVIDES:
              {
              auto k = getNamingConstantFor(toEval, _context->int_sort());
              // a divides b <-> k * a ==  b
              return k * args[0] == args[1];
              }

            case Theory::INT_UNARY_MINUS:
            case Theory::RAT_UNARY_MINUS:
            case Theory::REAL_UNARY_MINUS:
              return -args[0];

            case Theory::INT_PLUS:
            case Theory::RAT_PLUS:
            case Theory::REAL_PLUS:
              return args[0] + args[1];

            // Don't really need as it's preprocessed away
            case Theory::INT_MINUS:
            case Theory::RAT_MINUS:
            case Theory::REAL_MINUS:
              return args[0] - args[1];

            case Theory::INT_MULTIPLY:
            case Theory::RAT_MULTIPLY:
            case Theory::REAL_MULTIPLY:
              return args[0] * args[1];

            case Theory::RAT_QUOTIENT:
            case Theory::REAL_QUOTIENT:
              return args[0] / args[1];

            /** TPTP's ${quotient,remainder}_e */
            case Theory::INT_QUOTIENT_E:  return args[0] / args[1];          /* <--- same semantics of tptp and smtlib2 for int */
            case Theory::INT_REMAINDER_E: return z3::mod(args[0], args[1]);  /* <---                                            */
            case Theory::RAT_QUOTIENT_E:
            case Theory::REAL_QUOTIENT_E:  return                 tptp::quotient_e (args[0], args[1]);
            case Theory::RAT_REMAINDER_E:
            case Theory::REAL_REMAINDER_E: return tptp::remainder('e', tptp::quotient_e)(args[0], args[1]);

             /** {quotient,remainder}_t */
            case Theory::INT_QUOTIENT_T:  return tptp::liftInt(                tptp::quotient_t )(args[0],args[1]);
            case Theory::INT_REMAINDER_T: return tptp::liftInt(tptp::remainder('t', tptp::quotient_t))(args[0],args[1]);
            case Theory::RAT_QUOTIENT_T:
            case Theory::REAL_QUOTIENT_T: return                 tptp::quotient_t (args[0], args[1]);
            case Theory::REAL_REMAINDER_T:
            case Theory::RAT_REMAINDER_T: return tptp::remainder('t', tptp::quotient_t)(args[0], args[1]);

            /** {quotient,remainder}_f */
            case Theory::INT_QUOTIENT_F:  return tptp::liftInt(                tptp::quotient_f )(args[0], args[1]);
            case Theory::INT_REMAINDER_F: return tptp::liftInt(tptp::remainder('f', tptp::quotient_f))(args[0],args[1]);
            case Theory::RAT_QUOTIENT_F:
            case Theory::REAL_QUOTIENT_F: return                 tptp::quotient_f (args[0], args[1]);
            case Theory::REAL_REMAINDER_F:
            case Theory::RAT_REMAINDER_F: return tptp::remainder('f', tptp::quotient_f)(args[0], args[1]);


            case Theory::RAT_TO_INT:
            case Theory::REAL_TO_INT:
            case Theory::INT_FLOOR:
            case Theory::RAT_FLOOR:
            case Theory::REAL_FLOOR:
              return to_real(to_int(args[0]));

            case Theory::RAT_TO_REAL:
              return args[0];

            case Theory::INT_TO_REAL:
            case Theory::INT_TO_RAT:
              return to_real(args[0]);

            case Theory::INT_CEILING:
            case Theory::RAT_CEILING:
            case Theory::REAL_CEILING:
              return tptp::ceiling(args[0]);

            case Theory::INT_TRUNCATE:
            case Theory::RAT_TRUNCATE:
            case Theory::REAL_TRUNCATE:
              return tptp::truncate(args[0]);

            case Theory::INT_ROUND:
            case Theory::RAT_ROUND:
            case Theory::REAL_ROUND: {
                z3::expr t = args[0];
                z3::expr i = to_int(t);
                z3::expr i2 = i + _context->real_val(1,2);
                return ite(t > i2, i+1, ite(t==i2, ite(z3::mod(i, 2),i ,i+1 ),i));
              }

            case Theory::INT_ABS: {
                z3::expr t = args[0];
                return ite(t > 0, t, -t);
              }

            case Theory::INT_IS_INT:
            case Theory::RAT_IS_INT:
            case Theory::REAL_IS_INT:
              return z3::is_int(args[0]);

            case Theory::INT_LESS:
            case Theory::RAT_LESS:
            case Theory::REAL_LESS:
              return args[0] < args[1];

            case Theory::INT_GREATER:
            case Theory::RAT_GREATER:
            case Theory::REAL_GREATER:
              return args[0] > args[1];

            case Theory::INT_LESS_EQUAL:
            case Theory::RAT_LESS_EQUAL:
            case Theory::REAL_LESS_EQUAL:
              return args[0] <= args[1];

            case Theory::INT_GREATER_EQUAL:
            case Theory::RAT_GREATER_EQUAL:
            case Theory::REAL_GREATER_EQUAL:
              return args[0] >= args[1];

            default:
              {}//skip it and treat the function as uninterpretted
            }
          }
        }

        // uninterpretd function
        auto f = z3Function(Z3Interfacing::FuncOrPredId(trm));
        return f(f.arity(), args);
      })
      .apply(TermList(trm));
  return Representation(expr, std::move(defs));
}

Z3Interfacing::Representation Z3Interfacing::getRepresentation(SATLiteral slit)
{
  //First, does this represent a ground literal
  Literal* lit = _sat2fo.toFO(slit);
  if(lit && lit->ground()){
    //cout << "getRepresentation of " << lit->toString() << endl;
    // Now translate it into an SMT object
    try{
      auto repr = getRepresentation(lit);
      _exporter.apply([&expr = repr.expr](auto& exp){ exp.instantiate_expression(expr); });

      /* we name all literals in order to make z3 cache their truth values.
       * this gives a massive performance boost in many cases.              */

      z3::expr bname = getNameExpr(slit.var());
      _exporter.apply([&](auto& exp){ exp.instantiate_expression(bname); });
      z3::expr naming = (bname == repr.expr);
      _exporter.apply([&](auto& exp){ exp.instantiate_expression(naming); });
      repr.defs.push(naming);
      repr.expr = bname;

      if(_showZ3){
        std::cout << "[Z3] add (naming): " << naming << std::endl;
      }

      if(slit.isNegative()) {
        repr.expr = !repr.expr;
        _exporter.apply([&expr = repr.expr](auto& exp){ exp.instantiate_expression(expr); });
      }


      return repr;
    }catch(z3::exception& exception){
     reportSpiderFail();
     cout << "Z3 exception:\n" << exception.msg() << endl;
     ASSERTION_VIOLATION_REP("Failed to create Z3 rep for " + lit->toString());
    }
  } else {
    //if non ground then just create a propositional variable
    z3::expr e = getNameExpr(slit.var());
    e = slit.isPositive() ? e : !e;
    _exporter.apply([&](auto& exp){ exp.instantiate_expression(e); });
    return Representation(e, Stack<z3::expr>());
  }
}

SATClause* Z3Interfacing::getRefutation()
{
  return PrimitiveProofRecordingSATSolver::getRefutation();

  // TODO: optionally, we could try getting an unsat core from Z3 (could be smaller than all the added clauses so far)
  // NOTE: this will not (necessarily) be the same option as _unsatCore, which takes care of minimization of added assumptions
  // also ':core.minimize' might need to be set to get some effect
}

Z3Interfacing::~Z3Interfacing()
{
  _sorts.clear();
  _toZ3.clear();
  _fromZ3.clear();
  _exporter.apply([&](auto& e) { e.terminate(); });
}



bool Z3Interfacing::isNamedExpr(unsigned var) const
{ return _varNames.find(var); }

z3::expr Z3Interfacing::getNameExpr(unsigned var)
{
      // this method is called very often in runs with a lot of avatar reasoning. Cache the constants to avoid that z3 has to search for the string name in its function index
  return _varNames.getOrInit(var, [&]()
      { return z3_declare_const("v"+Lib::Int::toString(var), _context->bool_sort()); });
}


z3::expr Z3Interfacing::getNamingConstantFor(TermList toName, z3::sort sort)
{
  return _termIndexedConstants.getOrInit(toName, [&]()
    { return z3_declare_const("n" + toName.toString(), sort); });
}

z3::expr Z3Interfacing::getConst(Signature::Symbol* symb, z3::sort sort)
{
  return _constantNames.getOrInit(symb, [&]()
    // careful: keep native constants' names distinct from the above ones (hence the "c"-prefix below)
    { return z3_declare_const("c" + symb->name(), sort); });
}

} // namespace SAT

#endif /** if VZ3 **/
