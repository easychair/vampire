#include "LinearConstraint.hpp"

#include <algorithm>

using namespace Ordering;
using namespace std;

#define LC_LOG_LEVEL 5
#define LC_LOG

void Ordering::LinearConstraint::reset()
{
  nPosVars = 0;
  nNegVars = 0;
  sumCoeffs = 0;
  inverted = false;

  capacities.clear();
  requirements.clear();

  flowMatrix.clear();
  greaterThanY.clear();
}
void LinearConstraint::setProblem(const vector<pair<VarNum, Coeff>>& affineFunc, bool invert)
{
  // In the affine function, the coefficients are partitioned into positive and negative ones
  reset();

  unsigned i = 0;
  int coeff;
  while (i < affineFunc.size() && (coeff = affineFunc[i].second) > 0) {
    sumCoeffs += coeff;
    capacities.push_back(coeff);
    i++;
  }
  nPosVars = i;
  while (i < affineFunc.size()) {
    coeff = affineFunc[i].second;
    sumCoeffs += coeff;
    ASS(coeff < 0);
    requirements.push_back(-coeff);
    i++;
  }
  nNegVars = i - nPosVars;

  // If the sum of coefficients is negative, flip the requirements and the capacities
  if (sumCoeffs < 0 || invert) {
    ASS(!inverted);
    inverted = true;
    swap(capacities, requirements);
    swap(nPosVars, nNegVars);
  }

  greaterThanY.reshape(nNegVars, nPosVars);
  flowMatrix.reshape(nPosVars, nNegVars);
}

void Ordering::LinearConstraint::setOrdering(const std::vector<std::pair<VarNum, Coeff>>& affineFunc,
                                             const Kernel::TermPartialOrdering partialOrdering)
{
  Result res;
  unsigned xStart = 0;
  unsigned xEnd   = nPosVars;
  unsigned yStart = nPosVars;
  unsigned yEnd   = nPosVars + nNegVars;
  if (inverted) {
    xStart = nNegVars;
    xEnd   = nPosVars + nNegVars;
    yStart = 0;
    yEnd   = nNegVars;
  }

  for (VarAlias xIndex = xStart; xIndex < xEnd; xIndex++) {
    for (VarAlias yIndex = yStart; yIndex < yEnd; yIndex++) {
      auto x = TermList::var(affineFunc[xIndex].first);
      auto y = TermList::var(affineFunc[yIndex].first);
      partialOrdering.get(x, y, res);
      if (res == Kernel::Result::GREATER)
        greaterThanY.set(yIndex - nPosVars, xIndex, true);
    }
  }
}

void Ordering::LinearConstraint::setOrdering(const std::vector<std::pair<VarNum, Coeff>>& affineFunc,
                                             const std::vector<std::vector<bool>> partialOrdering)
{
  unsigned xStart = 0;
  unsigned xEnd   = nPosVars;
  unsigned yStart = nPosVars;
  unsigned yEnd   = nPosVars + nNegVars;
  if (inverted) {
    xStart = nNegVars;
    xEnd   = nPosVars + nNegVars;
    yStart = 0;
    yEnd   = nNegVars;
  }

  for (VarAlias xIndex = xStart; xIndex < xEnd; xIndex++) {
    for (VarAlias yIndex = yStart; yIndex < yEnd; yIndex++) {
      VarNum x = affineFunc[xIndex].first;
      VarNum y = affineFunc[yIndex].first;
      if (partialOrdering[x][y])
        greaterThanY.set(yIndex - yStart, xIndex - xStart, true);
    }
  }
}

bool Ordering::LinearConstraint::pruneLevel0()
{
  if (!inverted && constant > sumCoeffs)
    return true;
  if (inverted && constant < sumCoeffs)
    return true;
  return false;
}

bool Ordering::LinearConstraint::pruneLevel1()
{
  // if a variable y does not have any x that can simplify it, then prune
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++)
    if (greaterThanY.getSetOnRow(yAlias).size() == 0)
      return true;
  return false;
}

bool Ordering::LinearConstraint::preProcess()
{
  bool progress = true;
  while (progress) {
    progress = false;
    for (VarAlias yAlias=0; yAlias < requirements.size(); yAlias++) {
      // if no variable x is greater than y, then nothing can be concluded
      vector<pair<unsigned, bool>> greaterXs = greaterThanY.getSetOnRow(yAlias);
      if (greaterXs.size() == 0)
        return false;
      // if there is only one variable x greater than y, transfer the flow to y
      if (greaterXs.size() > 1)
        continue;

      VarAlias xAlias = greaterXs[0].first;
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

bool Ordering::LinearConstraint::dfsX(VarAlias x, std::vector<VarAlias> &path)
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

bool Ordering::LinearConstraint::dfsY(VarAlias y, std::vector<VarAlias> &path)
{
  seenY[y] = true;
  auto greaterXs = greaterThanY.getSetOnRow(y);
  if (greaterXs.size() == 0) {
    return false;
  }
  for (auto p : greaterXs) {
    VarAlias x = p.first;
    if (!seenX[x] && dfsX(x, path)) {
      path.push_back(y);
      return true;
    }
  }
  return false;
}

Flow Ordering::LinearConstraint::findPath(VarAlias sink, std::vector<VarAlias> &path)
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

bool Ordering::LinearConstraint::search()
{
  // Uses the maximum flow problem to find maches between the ys and the xs
  // TODO here find heuristic to bette choose the order of ys
  for (VarAlias y = 0; y < nNegVars; y++) {
    // requirement left for y
    auto greaterXs = greaterThanY.getSetOnRow(y);
    ASS(greaterXs.size() > 1);

    // first try to find some x with non zero capacity
    // TODO here find heuristic to bette choose the order of xs
    for (auto p : greaterXs) {
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

Ordering::Comparison Ordering::LinearConstraint::solve()
{
  if (!preProcess())
    return Comparison::Incomparable;
  if (search())
    return inverted ? Comparison::Less : Comparison::Greater;
  return Comparison::Incomparable;
}

void Ordering::LinearConstraint::removeXVariable(VarAlias xAlias)
{
  ASS_L(xAlias, nPosVars);
  capacities[xAlias] = capacities.back();
  capacities.pop_back();

  greaterThanY.delCol(xAlias);
  flowMatrix.delRow(xAlias);

  nPosVars--;
}

void Ordering::LinearConstraint::removeYVariable(VarAlias yAlias)
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

std::string Ordering::LinearConstraint::to_string() const
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

Ordering::Comparison LinearConstraint::getSign(const vector<pair<VarNum, Coeff>>& affineFunc,
                                               const TermPartialOrdering partialOrdering,
                                               Constant c)
{
  ASS(false);
  return Incomparable;
}

Ordering::Comparison LinearConstraint::getSign(const vector<pair<VarNum, Coeff>>& affineFunc,
                                               const vector<vector<bool>> partialOrdering,
                                               Constant c)
{
  bool failed = false;
  constant = c;

  setProblem(affineFunc);
  if (pruneLevel0())
    failed = true;
  else {
    setOrdering(affineFunc, partialOrdering);
    if (pruneLevel1())
      failed = true;
  }
  Ordering::Comparison result = Comparison::Incomparable;
  if (!failed)
    result = solve();

  if (result != Comparison::Incomparable)
    return result;

  // if we failed, we might need to try to solve the reverse problem if the sum of coefficients is zero
  if (sumCoeffs != 0)
    return Comparison::Incomparable;

  setProblem(affineFunc, true);
  if (pruneLevel0())
    return Comparison::Incomparable;
  setOrdering(affineFunc, partialOrdering);
  if (pruneLevel1())
    return Comparison::Incomparable;
  return solve();
}
