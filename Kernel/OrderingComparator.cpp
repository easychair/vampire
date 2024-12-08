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
 * @file OrderingComparator.cpp
 * Implements class OrderingComparator.
 */

#include "Lib/Stack.hpp"
#include "Lib/Environment.hpp"

#include "KBO.hpp"
#include "SubstHelper.hpp"

#include <deque>

#include "OrderingComparator.hpp"

using namespace std;

namespace Kernel {

std::ostream& operator<<(std::ostream& out, const OrderingComparator::BranchTag& t)
{
  switch (t) {
    case OrderingComparator::BranchTag::T_RESULT:
      out << "res";
      break;
    case OrderingComparator::BranchTag::T_TERM:
      out << "term";
      break;
    case OrderingComparator::BranchTag::T_POLY:
      out << "poly";
      break;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const OrderingComparator::Node& node)
{
  out << (OrderingComparator::BranchTag)node.tag << (node.ready?" ":"? ");
  switch (node.tag) {
    case OrderingComparator::BranchTag::T_RESULT: {
      out << node.result;
      break;
    }
    case OrderingComparator::BranchTag::T_POLY: {
      out << node.w;
      for (unsigned i = 0; i < node.varCoeffPairs->size(); i++) {
        const auto& [var, coeff] = (*node.varCoeffPairs)[i];
        ASS_NEQ(coeff,0);
        // output sign
        out << (coeff<0 ? " - " : " + ");
        // output coefficient
        if (abs(coeff) != 1) {
          out << abs(coeff) << " * ";
        }
        // output variable
        out << "X" << var;
      }
      break;
    }
    case OrderingComparator::BranchTag::T_TERM: {
      out << node.lhs << " " << node.rhs;
      break;
    }
  }
  return out;
}

std::ostream& operator<<(std::ostream& str, const OrderingComparator& comp)
{
  Stack<std::pair<const OrderingComparator::Branch*, unsigned>> todo;
  todo.push(std::make_pair(&comp._source,0));
  // Note: using this set we get a more compact representation
  DHSet<OrderingComparator::Node*> seen;

  while (todo.isNonEmpty()) {
    auto kv = todo.pop();
    for (unsigned i = 0; i < kv.second; i++) {
      if (i+1 == kv.second) {
        str << "  |--";
      } else {
        str << "  |  ";
      }
    }
    str << *kv.first->node() << std::endl;
    if (seen.insert(kv.first->node())) {
      if (kv.first->node()->tag!=OrderingComparator::BranchTag::T_RESULT) {
        todo.push(std::make_pair(&kv.first->node()->ngeBranch,kv.second+1));
        todo.push(std::make_pair(&kv.first->node()->eqBranch,kv.second+1));
        todo.push(std::make_pair(&kv.first->node()->gtBranch,kv.second+1));
      }
    }
  }
  return str;
}

std::string to_tikz_term(TermList t)
{
  if (t.isVar()) {
    switch (t.var()) {
      case 0:
        return "x";
      case 1:
        return "y";
      case 2:
        return "z";
      default:
        return "x_"+Int::toString(t.var());
    }
  }
  auto tt = t.term();
  auto res = tt->functionName() + "(";
  for (unsigned i = 0; i < tt->arity(); i++) {
    res += to_tikz_term(*tt->nthArgument(i)) + ",";
  }
  res += ")";
  return res;
}

std::string OrderingComparator::to_dot() const
{
  std::deque<const OrderingComparator::Node*> todo;
  todo.push_back(_source.node());
  // Note: using this set we get a more compact representation
  DHMap<const OrderingComparator::Node*,unsigned> seen;

  seen.insert(_source.node(), 0);
  unsigned cnt = 1;

  auto getId = [&todo,&seen,&cnt](auto n) {
    unsigned* ptr;
    if (seen.getValuePtr(n,ptr)) {
      todo.push_back(n);
      *ptr = cnt++;
    }
    return *ptr;
  };

  enum class EdgeTag {
    GT,
    EQ,
    INC,
    ALT,
  };
  auto getBranch = [&getId](unsigned from, auto to_node, EdgeTag tag) -> std::string
  {
    auto res = "n" + Int::toString(from) + " -> n" + Int::toString(getId(to_node));
    switch (tag) {
      case EdgeTag::GT:
        res += " [style = gtedge]";
        break;
      case EdgeTag::EQ:
        res += " [style = eqedge]";
        break;
      case EdgeTag::INC:
        res += " [style = ngeedge]";
        break;
      case EdgeTag::ALT:
        res += " [style = trueedge]";
        break;
    }
    return res + ";\n";
  };

  std::string nodes = "source [style = invisible, label = \"\"];\n";
  std::string edges = "source -> n0 [style = trueedge];\n";

  while (!todo.empty()) {
    auto n = todo.front();
    todo.pop_front();
    auto id = seen.get(n);

    std::string style = "";
    std::string label = "";
    if (n->ready) {
      style += "processed,";
    }
    switch (n->tag) {
      case BranchTag::T_RESULT: {
        if (!n->result) {
          // do not output anything for the fail node
          style += "sinknode,";
          label += "";
          break;
        }
        style += "datanode,";
        label += Int::toString((unsigned long)n->result);
        // edges += getBranch(id, alt, EdgeTag::ALT);
        break;
      }
      case BranchTag::T_TERM: {
        // nodes += "termnode] (n" + Int::toString(id) + ") at (" + Int::toString(x) + "," + Int::toString(y) + ") {$"
        //   + to_tikz_term(n->lhs) + " \\comp " + to_tikz_term(n->rhs) + "$};\n";
        style += "termnode,";
        label += "$" + to_tikz_term(n->lhs) + " \\comp " + to_tikz_term(n->rhs) + "$";
        edges += getBranch(id, n->gtBranch.node(), EdgeTag::GT);
        edges += getBranch(id, n->eqBranch.node(), EdgeTag::EQ);
        edges += getBranch(id, n->ngeBranch.node(), EdgeTag::INC);
        break;
      }
      case BranchTag::T_POLY: {
        style += "polynode,";
        bool first = true;
        label += "$";
        for (const auto& [var, coeff] : *n->varCoeffPairs) {
          label += coeff<0 ? "-" : (first ? "" : "+");
          first = false;
          auto a = std::abs(coeff);
          if (a==1) {
            label += to_tikz_term(TermList::var(var));
          } else {
            label += Int::toString(a) + "\\cdot " + to_tikz_term(TermList::var(var));
          }
        }
        if (n->w) {
          label += n->w<0 ? "-" : (first ? "" : "+");
          label += Int::toString((int)std::abs(n->w));
        }
        label += "$";
        edges += getBranch(id, n->gtBranch.node(), EdgeTag::GT);
        edges += getBranch(id, n->eqBranch.node(), EdgeTag::EQ);
        edges += getBranch(id, n->ngeBranch.node(), EdgeTag::INC);
        break;
      }
    }
    nodes += "n" + Int::toString(id) + " [\n  style = \"" + style + "\"\n  label = \"" + label + "\"\n];\n";
  }
  return "digraph {\nnodesep = 0;\nsep = 0;\nranksep = 0;\nesep = 0;\n" + nodes + "\n" + edges + "}\n";
}

OrderingComparator::OrderingComparator(const Ordering& ord)
: _ord(ord), _source(false), _fail(_source), _curr(&_source), _prev(nullptr)
{
}

OrderingComparator::~OrderingComparator() = default;

bool OrderingComparator::check(const SubstApplicator* applicator)
{
  _curr = &_source;
  _prev = nullptr;
  ASS(_curr);
  for (;;) {
    expand();
    auto node = _curr->node();
    ASS(node->ready);

    if (node->tag == BranchTag::T_RESULT) {
      return node->result;
    }

    Ordering::Result comp = Ordering::INCOMPARABLE;
    if (node->tag == BranchTag::T_TERM) {

      comp = _ord.isGreaterOrEq(AppliedTerm(node->lhs,applicator,true),AppliedTerm(node->rhs,applicator,true));
      // _trace.set({ node->lhs, node->rhs, comp });

    } else {
      ASS(node->tag == BranchTag::T_POLY);

      const auto& kbo = static_cast<const KBO&>(_ord);
      auto weight = node->w;
      ZIArray<int> varDiffs;
      for (const auto& [var, coeff] : *node->varCoeffPairs) {
        AppliedTerm tt(TermList::var(var), applicator, true);

        VariableIterator vit(tt.term);
        while (vit.hasNext()) {
          auto v = vit.next();
          varDiffs[v.var()] += coeff;
          // since the counts are sorted in descending order,
          // this can only mean we will fail
          if (varDiffs[v.var()]<0) {
            goto loop_end;
          }
        }
        int64_t w = kbo.computeWeight(tt);
        weight += coeff*w;
        // due to descending order of counts,
        // this also means failure
        if (coeff<0 && weight<0) {
          goto loop_end;
        }
      }

      if (weight > 0) {
        comp = Ordering::GREATER;
      } else if (weight == 0) {
        comp = Ordering::EQUAL;
      }
    }
loop_end:
    _prev = _curr;
    _curr = &node->getBranch(comp);
  }
  ASSERTION_VIOLATION;
}

void OrderingComparator::insert(const OrderingConstraints& ordCons)
{
  static Ordering::Result ordVals[] = { Ordering::EQUAL, Ordering::GREATER, Ordering::INCOMPARABLE };
  // we mutate current fail node and add a new one
  // cout << "before " << *this << endl;
  auto curr = &_fail;
  Branch newFail(false);

  curr->node()->~Node();
  curr->node()->ready = false;

  if (ordCons.isNonEmpty()) {
    curr->node()->tag = T_TERM;
    curr->node()->lhs = ordCons[0].lhs;
    curr->node()->rhs = ordCons[0].rhs;
    for (unsigned i = 0; i < 3; i++) {
      if (ordVals[i] != ordCons[0].rel) {
        curr->node()->getBranch(ordVals[i]) = newFail;
      }
    }
    curr = &curr->node()->getBranch(ordCons[0].rel);
    for (unsigned i = 1; i < ordCons.size(); i++) {
      auto [lhs,rhs,rel] = ordCons[i];
      *curr = OrderingComparator::Branch(lhs, rhs);
      for (unsigned i = 0; i < 3; i++) {
        if (ordVals[i] != rel) {
          curr->node()->getBranch(ordVals[i]) = newFail;
        }
      }
      curr = &curr->node()->getBranch(rel);
    }
    *curr = Branch(true);
  } else {
    curr->node()->tag = T_RESULT;
    curr->node()->result = true;
    curr->node()->ready = true;
  }

  _fail = newFail;
  ASS(_fail.node()->ready);
  ASS(!_fail.node()->result);
}

void OrderingComparator::expand()
{
  ASS(_curr->node());
  while (!_curr->node()->ready)
  {
    auto node = _curr->node();

    ASS_NEQ(node->tag, BranchTag::T_RESULT);

    if (node->tag == BranchTag::T_POLY) {
      // if refcnt > 1 we copy the node and
      // we can also safely use the original
      if (node->refcnt > 1) {
        auto varCoeffPairs = new Stack<VarCoeffPair>(*node->varCoeffPairs);
        *_curr = Branch(node->w, varCoeffPairs);
        _curr->node()->gtBranch = node->gtBranch;
        _curr->node()->eqBranch = node->eqBranch;
        _curr->node()->ngeBranch = node->ngeBranch;
      }

      processPolyCase();
      continue;
    }

    // Use compare here to filter out as many
    // precomputable comparisons as possible.
    auto comp = _ord.compare(node->lhs,node->rhs);
    if (comp != Ordering::INCOMPARABLE) {
      if (comp == Ordering::LESS) {
        *_curr = node->ngeBranch;
      } else if (comp == Ordering::GREATER) {
        *_curr = node->gtBranch;
      } else {
        *_curr = node->eqBranch;
      }
      continue;
    }
    // If we have a variable, we cannot expand further.
    if (node->lhs.isVar() || node->rhs.isVar()) {
      processVarCase();
      continue;
    }

    expandTermCase();
  }
}

void OrderingComparator::expandTermCase()
{
  ASS(_curr->node() && !_curr->node()->ready);
  _curr->node()->ready = true;
}

void OrderingComparator::processPolyCase()
{
  auto node = _curr->node();
  auto varCoeffPairs = node->varCoeffPairs;

  auto trace = getCurrentTrace();

  // If this happens this branch is invalid.
  // Since normal execution cannot reach it,
  // we can put a "success" here to make things
  // easier in subsumption/simplification.
  if (!trace) {
    *_curr = Branch(true);
    return;
  }

  unsigned pos = 0;
  unsigned neg = 0;
  for (unsigned i = 0; i < varCoeffPairs->size();) {
    auto& [var, coeff] = (*varCoeffPairs)[i];
    for (unsigned j = i+1; j < varCoeffPairs->size();) {
      auto& [var2, coeff2] = (*varCoeffPairs)[j];
      Ordering::Result res;
      if (trace->get(TermList::var(var), TermList::var(var2), res) && res == Ordering::EQUAL) {
        coeff += coeff2;
        swap((*varCoeffPairs)[j],varCoeffPairs->top());
        varCoeffPairs->pop();
        continue;
      }
      j++;
    }
    if (coeff == 0) {
      swap((*varCoeffPairs)[i],varCoeffPairs->top());
      varCoeffPairs->pop();
      continue;
    } else if (coeff > 0) {
      pos++;
    } else {
      neg++;
    }
    i++;
  }

  // TODO HERE

  if (node->w == 0 && pos == 0 && neg == 0) {
    *_curr = node->eqBranch;
    return;
  }
  if (node->w >= 0 && neg == 0) {
    *_curr = node->gtBranch;
    return;
  }
  if (node->w <= 0 && pos == 0) {
    *_curr = node->ngeBranch;
    return;
  }
  sort(varCoeffPairs->begin(),varCoeffPairs->end(),[](const auto& e1, const auto& e2) {
    return e1.second>e2.second;
  });
  _curr->node()->ready = true;
  _curr->node()->trace = trace;
}

void OrderingComparator::processVarCase()
{
  auto node = _curr->node();
  auto trace = getCurrentTrace();

  // If this happens this branch is invalid.
  // Since normal execution cannot reach it,
  // we can put a "success" here to make things
  // easier in subsumption/simplification.
  if (!trace) {
    *_curr = Branch(true);
    return;
  }
  Ordering::Result val;
  if (trace->get(node->lhs, node->rhs, val)) {
    if (val == Ordering::GREATER) {
      *_curr = node->gtBranch;
    } else if (val == Ordering::EQUAL) {
      *_curr = node->eqBranch;
    } else {
      *_curr = node->ngeBranch;
    }
    return;
  }
  // if refcnt > 1 we copy the node and
  // we can also safely use the original
  if (node->refcnt > 1) {
    *_curr = Branch(node->lhs, node->rhs);
    _curr->node()->eqBranch = node->eqBranch;
    _curr->node()->gtBranch = node->gtBranch;
    _curr->node()->ngeBranch = node->ngeBranch;
  }
  _curr->node()->ready = true;
  _curr->node()->trace = trace;
}

const OrderingComparator::Trace* OrderingComparator::getCurrentTrace()
{
  ASS(!_curr->node()->ready || _curr->node() == _fail.node());

  if (!_prev) {
    return Trace::getEmpty(_ord);
  }

  ASS(_prev->node()->ready);
  ASS(_prev->node()->trace);

  switch (_prev->node()->tag) {
    case BranchTag::T_TERM: {
      auto lhs = _prev->node()->lhs;
      auto rhs = _prev->node()->rhs;
      Ordering::Result res;
      if (_curr == &_prev->node()->eqBranch) {
        res = Ordering::EQUAL;
      } else if (_curr == &_prev->node()->gtBranch) {
        res = Ordering::GREATER;
      } else {
        res = Ordering::INCOMPARABLE;
      }
      return Trace::set(_prev->node()->trace, { lhs, rhs, res });
    }
    case BranchTag::T_RESULT:
    case BranchTag::T_POLY: {
      return _prev->node()->trace;
    }
  }
  ASSERTION_VIOLATION;
}

OrderingComparator::Branch::~Branch()
{
  setNode(nullptr);
}

OrderingComparator::Branch::Branch(const Branch& other)
{
  setNode(other._node);
}

OrderingComparator::Branch& OrderingComparator::Branch::operator=(const Branch& other)
{
  if (&other==this) {
    return *this;
  }
  setNode(other.node());
  return *this;
}

OrderingComparator::Branch::Branch(Branch&& other)
  : _node(other._node)
{
  other._node = nullptr;
}

OrderingComparator::Branch& OrderingComparator::Branch::operator=(Branch&& other)
{
  if (&other==this) {
    return *this;
  }
  swap(_node,other._node);
  return *this;
}

OrderingComparator::Node::~Node()
{
  if (tag == T_POLY) {
    delete varCoeffPairs;
  }
  ready = false;
  trace = nullptr;
}

void OrderingComparator::Node::incRefCnt()
{
  refcnt++;
}

void OrderingComparator::Node::decRefCnt()
{
  ASS(refcnt>=0);
  refcnt--;
  if (refcnt==0) {
    delete this;
  }
}

}
