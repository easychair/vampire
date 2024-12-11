#include "LinearConstraint.hpp"

#include <algorithm>

using namespace Lib;
using namespace std;

#define LC_LOG_LEVEL 5
#define LC_LOG

void LinearConstraint::reset()
{
  nPosVars = 0;
  nNegVars = 0;
  sumCoeffs = 0;

  capacities.clear();
  requirements.clear();

  flowMatrix.clear();
  greaterThanY.clear();
}

void Lib::LinearConstraint::setProblem(const Stack<std::pair<VarNum, Coeff>>& posVars,
                                       const Stack<std::pair<VarNum, Coeff>>& negVars)
{
  // In the affine function, the coefficients are partitioned into positive and negative ones
  reset();

  for (pair<VarNum, Coeff> var : posVars) {
    Coeff coeff = var.second;
    ASS_NEQ(coeff, 0);
    ASS( inverted || (coeff > 0));
    ASS(!inverted || (coeff < 0));
    Flow capa = abs(coeff);
    capacities.push_back(capa);
    sumCoeffs += capa;
  }

  for (pair<VarNum, Coeff> var : negVars) {
    Coeff coeff = var.second;
    ASS_NEQ(coeff, 0);
    ASS( inverted || (coeff < 0));
    ASS(!inverted || (coeff > 0));
    Flow req = abs(coeff);
    requirements.push_back(req);
    sumCoeffs -= req;
  }

  if (sumCoeffs < 0) {
    ASS(!inverted);
    inverted = true;
    swap(capacities, requirements);
  }

  nPosVars = capacities.size();
  nNegVars = requirements.size();
  greaterThanY.reshape(nNegVars, nPosVars);
  flowMatrix.reshape(nPosVars, nNegVars);
}

void Lib::LinearConstraint::setOrdering(const Stack<std::pair<VarNum, Coeff>>& posVars,
                                        const Stack<std::pair<VarNum, Coeff>>& negVars,
                                        const Kernel::TermPartialOrdering* partialOrdering)
{
  Result res;
  ASS(inverted || nPosVars == posVars.size());
  ASS(inverted || nNegVars == negVars.size());
  ASS(!inverted || nNegVars == posVars.size());
  ASS(!inverted || nPosVars == negVars.size());
  for (unsigned i = 0; i < nPosVars; i++) {
    for (unsigned j = 0; j < nNegVars; j++) {
      auto x = TermList::var(!inverted ? posVars[i].first : negVars[i].first);
      auto y = TermList::var(!inverted ? negVars[j].first : posVars[j].first);
      partialOrdering->get(x, y, res);
      if (res == Kernel::Result::GREATER)
        greaterThanY.set(j, i, true);
    }
  }
}

void Lib::LinearConstraint::setOrdering(const Stack<std::pair<VarNum, Coeff>>& posVars,
                                        const Stack<std::pair<VarNum, Coeff>>& negVars,
                                        const std::vector<std::vector<bool>>& partialOrdering)
{
  ASS(inverted || nPosVars == posVars.size());
  ASS(inverted || nNegVars == negVars.size());
  ASS(!inverted || nNegVars == posVars.size());
  ASS(!inverted || nPosVars == negVars.size());
  for (unsigned i = 0; i < nPosVars; i++) {
    for (unsigned j = 0; j < nNegVars; j++) {
      VarNum x = !inverted ? posVars[i].first : negVars[i].first;
      VarNum y = !inverted ? negVars[j].first : posVars[j].first;
      if (partialOrdering[x][y])
        greaterThanY.set(j, i, true);
    }
  }
}

bool LinearConstraint::pruneLevel0()
{
  if (!inverted && constant > sumCoeffs)
    return true;
  if (inverted && constant < sumCoeffs)
    return true;
  return false;
}

bool LinearConstraint::pruneLevel1()
{
  // if a variable y does not have any x that can simplify it, then prune
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++)
    if (greaterThanY.getSetOnRow(yAlias).size() == 0)
      return true;
  return false;
}

bool LinearConstraint::preProcess()
{
  bool progress = true;
  while (progress) {
    progress = false;
    for (VarAlias yAlias=0; yAlias < requirements.size(); yAlias++) {
      // if no variable x is greater than y, then nothing can be concluded
      vector<pair<unsigned, bool>> xs = greaterThanY.getSetOnRow(yAlias);
      if (xs.size() == 0)
        return false;
      // if there is only one variable x greater than y, transfer the flow to y
      if (xs.size() > 1)
        continue;

      VarAlias xAlias = xs[0].first;
      ASS_L(xAlias, nPosVars);
      if (capacities[xAlias] < requirements[yAlias])
        return false;

      Flow transfer = requirements[yAlias];
      capacities[xAlias] -= transfer;
      requirements[yAlias] = 0;

      // remove the variable y
      removeYVariable(yAlias);
      yAlias--;
      if (capacities[xAlias] == 0)
        removeXVariable(xAlias);
      progress = true;
    }
    // further, if a variable x is only greater than one variable y, then transfer the flow to y
    for (VarAlias xAlias=0; xAlias < capacities.size(); xAlias++) {
      vector<unsigned> smallerYs;
      for (VarAlias yAlias = 0; yAlias < requirements.size(); yAlias++)
        if (greaterThanY.get(yAlias, xAlias))
          smallerYs.push_back(yAlias);

      if (smallerYs.size() == 1) {
        VarAlias yAlias = smallerYs[0];
        ASS_L(yAlias, nNegVars);
        Flow transfer = min(capacities[xAlias], requirements[yAlias]);
        capacities[xAlias]   -= transfer;
        requirements[yAlias] -= transfer;
        // the flow should be irrelevant for the rest of the search. We should be able to remove it.
        flowMatrix.set(xAlias, yAlias, flowMatrix.get(xAlias, yAlias) + transfer);
        progress = true;

        if (requirements[yAlias] == 0)
          removeYVariable(yAlias);
        if (capacities[xAlias] == 0)
          removeXVariable(xAlias--);

        progress = true;
      }
    }
  }
  return true;
}

bool LinearConstraint::dfsX(VarAlias x, std::vector<VarAlias> &path)
{
  seenX[x] = true;
  auto transfers = flowMatrix.getSetOnRow(x);
  if (capacities[x] > 0) {
    path.push_back(x);
    return true;
  }

  if (transfers.size() == 0) {
    return false;
  }

  for (auto p : transfers) {
    VarAlias y = p.first;
    ASS(p.second != 0);
    if (!seenY[y] && dfsY(y, path)) {
      path.push_back(x);
      return true;
    }
  }
  return false;
}

bool LinearConstraint::dfsY(VarAlias y, std::vector<VarAlias> &path)
{
  seenY[y] = true;
  auto xs = greaterThanY.getSetOnRow(y);
  if (xs.size() == 0) {
    return false;
  }
  for (auto p : xs) {
    VarAlias x = p.first;
    if (!seenX[x] && dfsX(x, path)) {
      path.push_back(y);
      return true;
    }
  }
  return false;
}

LinearConstraint::Flow LinearConstraint::findPath(VarAlias sink, std::vector<VarAlias> &path)
{
  seenX.clear();
  seenX.resize(nPosVars, false);
  seenY.clear();
  seenY.resize(nNegVars, false);

  if (!dfsY(sink, path))
    return 0;
  // reverse the path
  reverse(path.begin(), path.end());
  ASS(path.size() % 2 == 0);
  // The path should look like [y0, x0, y1, ..., xn]
  Flow bottleneck = min(capacities[path.back()], requirements[path[0]]);
  for (unsigned i = 1; i < path.size() - 1; i += 2)
    bottleneck = min(bottleneck, flowMatrix.get(path[i], path[i + 1]));
  ASS(bottleneck != 0);

  return bottleneck;
}

bool LinearConstraint::search()
{
  // Uses the maximum flow problem to find maches between the ys and the xs
  // TODO here find heuristic to better choose the order of ys
  for (VarAlias y = 0; y < nNegVars; y++) {
    // requirement left for y
    auto xs = greaterThanY.getSetOnRow(y);
    ASS(xs.size() > 1);

    // first try to find some x with non zero capacity
    // TODO here find heuristic to better choose the order of xs
    for (auto p : xs) {
      Flow b = requirements[y];
      ASS(b > 0);
      VarAlias x = p.first;
      // capacity left for x
      Flow a = capacities[x];
      ASS(a >= 0);
      if (a == 0)
        continue;
      // transferred value
      Flow t = min(a, b);
      Flow prev = flowMatrix.get(x, y);
      flowMatrix.set(x, y, prev + t);
      capacities[x] -= t;
      requirements[y] -= t;
      if (b == t)
        break;
    }
    if (requirements[y] == 0)
      continue;

    // we failed to simplify y with the remaining x
    // we need to redirect some of the flow
    vector<VarAlias> path;
    while (requirements[y] > 0) {
      // We search for a path from which to redirect the flow.
      /**
       *     0/1    0/1    1/2
       *     x0     x1     x2
       *     |    / |    /
       *    1|  0/ 1|  1/
       *     |  /   |  /
       *     y0     y1
       *     1/2    2/2
       * Here, for example, y0 needs one more, but cannot draw if from either x0 nor x1
       * So we find the path [y0, x1, y1, x2] such that
       * - x2 now gives to y1,
       * - x1 now gives to y0
       *     0/1    0/1    0/2
       *     x0     x1     x2
       *     |    / |    /
       *    1|  1/ 0|  2/
       *     |  /   |  /
       *     y0     y1
       *     2/2    2/2
       */
      // transfer the flow from x to y
      Flow t = findPath(y, path);
      ASS(t >= 0);
      if (t == 0) {
        return false;
      }
      // for each upstream edge, increase the flow.
      // for each downstream edge, decrease the flow.
      for (unsigned i = 0; i < path.size() - 1; i++) {
        if (i % 2 == 0) {
          // upstream edge
          VarAlias y = path[i];
          VarAlias x = path[i+1];
          flowMatrix.set(x, y, flowMatrix.get(x, y) + t);
        } else {
          // downstream edge
          VarAlias x = path[i];
          VarAlias y = path[i+1];
          flowMatrix.set(x, y, flowMatrix.get(x, y) - t);
          ASS(flowMatrix.get(x, y) >= 0)
        }
      }
      capacities[path.back()] -= t;
      requirements[y] -= t;
    }
    if (requirements[y] > 0)
      return false;
  }
  return true;
}

Ordering::Result LinearConstraint::solve()
{
  cout << to_string() << endl;
  if (!preProcess())
    return Result::INCOMPARABLE;
  if (nNegVars == 0 || search()) {
    cout << "returns: " << (inverted ? Result::LESS : Result::GREATER) << endl;
    return inverted ? Result::LESS : Result::GREATER;
  }
  return Result::INCOMPARABLE;
}

void LinearConstraint::removeXVariable(VarAlias xAlias)
{
  ASS_L(xAlias, nPosVars);
  capacities[xAlias] = capacities.back();
  capacities.pop_back();

  greaterThanY.delCol(xAlias);
  flowMatrix.delRow(xAlias);

  nPosVars--;
}

void LinearConstraint::removeYVariable(VarAlias yAlias)
{
  ASS_L(yAlias, nNegVars);
  requirements[yAlias] = requirements.back();
  requirements.pop_back();

  greaterThanY.delRow(yAlias);
  flowMatrix.delCol(yAlias);

  nNegVars--;
}

LinearConstraint::LinearConstraint() :
  flowMatrix(0, 0, 0),
  greaterThanY(0, 0, false)
{
}

std::string LinearConstraint::to_string() const
{
  string s = "LinearConstraint\n";
  s += "nPosVars: ";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    s += "  x" + std::to_string(xAlias) + ": " + std::to_string(capacities[xAlias]) + "  ";
  }
  s += "\nnNegVars: ";
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
    s += "  y" + std::to_string(yAlias) + ": " + std::to_string(requirements[yAlias]) + "  ";
  }
  s += "\nflowMatrix: \n";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
      Flow flow = flowMatrix.get(xAlias, yAlias);
      if (flow != 0)
        s += "  x" + std::to_string(xAlias) + " -> y" + std::to_string(yAlias) + ": " + std::to_string(flowMatrix.get(xAlias, yAlias)) + "\n";
    }
  }
  s += "Partial Ordering :\n";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++)
      if (greaterThanY.get(yAlias, xAlias))
        s += "  x" + std::to_string(xAlias) + " > y" + std::to_string(yAlias) + "  ";
    s += "\n";
  }

  return s;
}


Result Lib::LinearConstraint::getSign(const Constant& c,
                                      const Lib::Stack<std::pair<VarNum, Coeff>>& posVars,
                                      const Lib::Stack<std::pair<VarNum, Coeff>>& negVars,
                                      const Kernel::TermPartialOrdering* partialOrdering)
{
  cout << ".";
  bool failed = false;
  inverted = false;
  constant = c;

  setProblem(posVars, negVars);

  if (nPosVars == 0 && nNegVars == 0) {
    if (constant == 0)
      return Result::EQUAL;
    if (constant > 0)
      return Result::LESS;
    return Result::GREATER;
  }

  if (pruneLevel0())
    failed = true;
  else {
    setOrdering(posVars, negVars, partialOrdering);
    if (pruneLevel1())
      failed = true;
  }

  Result result = Result::INCOMPARABLE;
  if (!failed)
    result = solve();

  if (result != Result::INCOMPARABLE)
    return result;

  // if we failed, we might need to try to solve the reverse problem if the sum of coefficients is zero
  if (sumCoeffs != 0)
    return Result::INCOMPARABLE;

  ASS(!inverted);
  inverted = true;
  setProblem(negVars, posVars);
  if (pruneLevel0())
    return Result::INCOMPARABLE;
  setOrdering(posVars, negVars, partialOrdering);
  if (pruneLevel1())
    return Result::INCOMPARABLE;
  return solve();
}

Result Lib::LinearConstraint::getSign(const Constant& c,
                                      const Lib::Stack<std::pair<VarNum, Coeff>>& posVars,
                                      const Lib::Stack<std::pair<VarNum, Coeff>>& negVars,
                                      const std::vector<std::vector<bool>> partialOrdering)
{
  bool failed = false;
  inverted = false;
  constant = c;

  setProblem(posVars, negVars);

  if (nPosVars == 0 && nNegVars == 0) {
    if (constant == 0)
      return Result::EQUAL;
    if (constant > 0)
      return Result::LESS;
    return Result::GREATER;
  }

  if (pruneLevel0())
    failed = true;
  else {
    setOrdering(posVars, negVars, partialOrdering);
    if (pruneLevel1())
      failed = true;
  }

  Result result = Result::INCOMPARABLE;
  if (!failed)
    result = solve();

  if (result != Result::INCOMPARABLE)
    return result;

  // if we failed, we might need to try to solve the reverse problem if the sum of coefficients is zero
  if (sumCoeffs != 0)
    return Result::INCOMPARABLE;

  ASS(!inverted);
  inverted = true;
  setProblem(negVars, posVars);
  if (pruneLevel0())
    return Result::INCOMPARABLE;
  setOrdering(posVars, negVars, partialOrdering);
  if (pruneLevel1())
    return Result::INCOMPARABLE;
  return solve();
}
