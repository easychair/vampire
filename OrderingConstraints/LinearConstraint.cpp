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
  capacities.clear();
  requirements.clear();

  flowMatrix.reset();

  greaterThanY.reset();
  smallerThanX.reset();
}
bool LinearConstraint::setProblem(const vector<pair<VarNum, Coeff>>& affineFunc)
{
  // In the affine function, the coefficients are partitioned into positive and negative ones
  reset();

  sumCoeffs = 0;

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
    ASS(coeff < 0);
    requirements.push_back(-coeff);
    i++;
  }
  nNegVars = i - nPosVars;

  // If the sum of coefficients is negative, flip the requirements and the capacities
  if (sumCoeffs < 0) {
    vector<Flow>& tmp = capacities;
    capacities = requirements;
    requirements = tmp;

    nPosVars = nPosVars ^ nNegVars;
    nNegVars = nPosVars ^ nNegVars;
    nPosVars = nPosVars ^ nNegVars;
  }

  if (pruneLevel0())
    return false;
  greaterThanY.reshape(nNegVars, nPosVars);
  smallerThanX.reshape(nPosVars, nNegVars);
  flowMatrix.reshape(nPosVars, nNegVars);
  return true;
}

bool Ordering::LinearConstraint::setOrdering(const std::vector<std::pair<VarNum, Coeff>>& affineFunc, const Kernel::TermPartialOrdering partialOrdering)
{
  Result res;
  for (VarAlias xIndex = 0; xIndex < nPosVars; xIndex++) {
    auto x = TermList::var(affineFunc[xIndex].first);
    for (VarAlias yIndex = nPosVars; yIndex < nPosVars + nNegVars; yIndex++) {
      auto y = TermList::var(affineFunc[yIndex].first);
      partialOrdering.get(x, y, res);
      if (res == Kernel::Result::GREATER) {
        VarAlias xAlias = xIndex;
        VarAlias yAlias = yIndex - nPosVars;
        greaterThanY.set(xAlias, yAlias, true);
        smallerThanX.set(yAlias, xAlias, true);
      }
    }
  }

  if (pruneLevel1())
    return false;
  return true;
}

bool Ordering::LinearConstraint::setOrdering(const std::vector<std::pair<VarNum, Coeff>>& affineFunc,
                                             const std::vector<std::vector<bool>> partialOrdering)
{
  for (VarAlias xIndex = 0; xIndex < nPosVars; xIndex++) {
    for (VarAlias yIndex = nPosVars; yIndex < nPosVars + nNegVars; yIndex++) {
      if (partialOrdering[xIndex][yIndex]) {
        VarAlias xAlias = xIndex;
        VarAlias yAlias = yIndex - nPosVars;
        cout << "Setting " << xAlias << " > " << yAlias << endl;
        greaterThanY.set(yAlias, xAlias, true);
        smallerThanX.set(xAlias, yAlias, true);
      }
    }
  }

  if (pruneLevel1())
    return false;
  return true;
}

bool Ordering::LinearConstraint::pruneLevel0()
{
  if (sumCoeffs > 0 && constant > sumCoeffs)
    return true;
  if (sumCoeffs < 0 && constant < sumCoeffs)
    return true;
  return false;
}

bool Ordering::LinearConstraint::pruneLevel1()
{
  // if a variable y does not have any x that can simplify it, then prune
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
    if (greaterThanY.getSetOnRow(yAlias).size() == 0)
      return true;
  }
  return false;
}

bool Ordering::LinearConstraint::preProcess()
{
  bool progress = true;
  cout << "Preprocessing" << endl;
  while (progress) {
    cout << to_string() << endl;
    progress = false;
    for (VarAlias yAlias=0; yAlias < requirements.size(); yAlias++) {
      // if no variable x is greater than y, then nothing can be concluded
      vector<pair<unsigned, bool>> greaterXs = greaterThanY.getSetOnRow(yAlias);
      if (greaterXs.size() == 0)
        return false;
      // if there is only one variable x greater than y, transfer the flow to y
      if (greaterXs.size() > 1) {
        continue;
      }
      VarAlias xAlias = greaterXs[0].first;
      if (capacities[xAlias] < requirements[yAlias]) {
        return false;
      }
      cout << "y" << yAlias << " has to be simplified by x" << xAlias << endl;
      Flow transfer = requirements[yAlias];
      capacities[xAlias] -= transfer;
      requirements[yAlias] = 0;
      // the flow should be irrelevant for the rest of the search. We should be able to remove it.
      flowMatrix.set(xAlias, yAlias, transfer);

      // remove the variable y
      removeYVariable(yAlias);
      yAlias--;
      if (capacities[xAlias] == 0) {
        removeXVariable(xAlias);
      }
      cout << to_string() << endl;
      progress = true;
    }

    // further, if a variable x is only greater than one variable y, then transfer the flow to y
    for (VarAlias xAlias=0; xAlias < capacities.size(); xAlias++) {
      vector<pair<unsigned, bool>> smallerYs = smallerThanX.getSetOnRow(xAlias);
      if (smallerYs.size() == 1) {
        cout << "x" << xAlias << " can freely simplify y" << smallerYs[0].first << endl;
        VarAlias yAlias = smallerYs[0].first;
        Flow transfer = min(capacities[xAlias], requirements[yAlias]);
        capacities[xAlias] -= transfer;
        requirements[yAlias] -= transfer;
        // the flow should be irrelevant for the rest of the search. We should be able to remove it.
        cout << "Transferring " << transfer << " from x" << xAlias << " to y" << yAlias << endl;
        flowMatrix.set(xAlias, yAlias, flowMatrix.get(xAlias, yAlias, 0) + transfer);
        progress = true;

        // remove the variable y
        if (requirements[yAlias] == 0) {
          removeYVariable(yAlias);
        }
        if (capacities[xAlias] == 0) {
          removeXVariable(xAlias);
          xAlias--;
        }
        cout << to_string() << endl;
        progress = true;
      }
    }
  }
  return true;
}

void Ordering::LinearConstraint::removeXVariable(VarAlias xAlias)
{
  capacities[xAlias] = capacities.back();

  smallerThanX.swapRows(xAlias, nPosVars - 1);
  // we don't need to reshape the matrices. We will just ignore the last row
  // they will just be a bit too large

  // we however need to remove the xAliases from the partial ordering
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
    greaterThanY.del(xAlias, yAlias);
  }
  capacities.pop_back();
  nPosVars--;
}

void Ordering::LinearConstraint::removeYVariable(VarAlias yAlias)
{
  requirements[yAlias] = requirements.back();

  greaterThanY.swapRows(yAlias, nNegVars - 1);
  // as above, we don't need to reshape the matrix. We will just ignore the last column

  // we however need to remove the yAliases from the partial ordering
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    smallerThanX.del(xAlias, yAlias);
  }

  nNegVars--;
  requirements.pop_back();
}

LinearConstraint::LinearConstraint() :
  flowMatrix(0, 0),
  greaterThanY(0, 0),
  smallerThanX(0, 0)
{
}

std::string Ordering::LinearConstraint::to_string() const
{
  string s = "LinearConstraint\n";
  s += "nPosVars: \n";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    s += "  x" + std::to_string(xAlias) + ": " + std::to_string(capacities[xAlias]) + "\n";
  }
  s += "nNegVars: \n";
  for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
    s += "  y" + std::to_string(yAlias) + ": " + std::to_string(requirements[yAlias]) + "\n";
  }
  s += "flowMatrix: \n";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
      Flow flow = flowMatrix.get(xAlias, yAlias, 0);
      if (flow != 0)
        s += "  x" + std::to_string(xAlias) + " -> y" + std::to_string(yAlias) + ": " + std::to_string(flowMatrix.get(xAlias, yAlias, 0)) + "\n";
    }
  }
  s += "Partial Ordering :\n";
  for (VarAlias xAlias=0; xAlias < nPosVars; xAlias++) {
    for (VarAlias yAlias=0; yAlias < nNegVars; yAlias++) {
      if (greaterThanY.get(xAlias, yAlias, false)) {
        s += "  x" + std::to_string(xAlias) + " > y" + std::to_string(yAlias) + "\n";
      }
    }
  }
  return s;
}

Ordering::Comparison LinearConstraint::getSign(const vector<pair<VarNum, Coeff>>& affineFunc,
                                                         const TermPartialOrdering partialOrdering,
                                                         Constant c)
{
  constant = c;
  if (!setProblem(affineFunc))
    return Comparison::Incomparable;
  if (!setOrdering(affineFunc, partialOrdering))
    return Comparison::Incomparable;
  // todo deal with the case when sumCoeffs is 0
  if (!preProcess())
    return Comparison::Incomparable;
  if (nNegVars == 0)
    return sumCoeffs > 0 ? Comparison::Greater : Comparison::Less;


  return Incomparable;
}

Ordering::Comparison LinearConstraint::getSign(const vector<pair<VarNum, Coeff>>& affineFunc,
                                                         const vector<vector<bool>> partialOrdering,
                                                         Constant c)
{
  constant = c;
  if (!setProblem(affineFunc)) {
    cout << "setProblem failed" << endl;
    return Comparison::Incomparable;
  }
  cout << "setProblem done" << endl;
  cout << to_string() << endl;
  if (!setOrdering(affineFunc, partialOrdering)) {
    cout << "setOrdering failed" << endl;
    return Comparison::Incomparable;
  }
  cout << "setOrdering done" << endl;
  cout << to_string() << endl;
  if (!preProcess())
    return Comparison::Incomparable;
  if (nNegVars == 0) {
    cout << "Success!" << endl;
    cout << to_string() << endl;
    return sumCoeffs > 0 ? Comparison::Greater : Comparison::Less;
  }

  return Incomparable;
}
