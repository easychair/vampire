/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#include "Test/UnitTesting.hpp"
#include "Test/SyntaxSugar.hpp"
#include "Test/TestUtils.hpp"
#include "Test/GenerationTester.hpp"

#include "OrderingConstraints/SparseMatrix.hpp"
#include "OrderingConstraints/LinearConstraint.hpp"

#include <algorithm>
#include <random>

using namespace std;
using namespace Lib;
using namespace Test;

using VarNum   = LinearConstraint::VarNum;
using Coeff    = LinearConstraint::Coeff;
using Flow     = LinearConstraint::Flow;
using VarAlias = LinearConstraint::VarAlias;

TEST_FUN(MatrixRetrieval01) {
  SparseMatrix<int> matrix(3, 3, -1);
  matrix.set(0, 0, 1);

  ASS_EQ(matrix.get(0, 0), 1);
  for (unsigned i = 0; i < 3; i++)
    for (unsigned j  = 0; j < 3; j++)
      if (i != 0 || j != 0)
        ASS_EQ(matrix.get(i, j), -1);
}

TEST_FUN(MatrixRetrieval02) {
  SparseMatrix<int> matrix(1, 5, -1);
  for (unsigned i = 0; i < 5; i++) {
    matrix.set(0, i, i);
  }

  for (int i = 0; i < 5; i++)
    ASS_EQ(matrix.get(0, i), i);

}

TEST_FUN(MatrixNothingSet) {
  SparseMatrix<int> matrix(3, 3, -1);
  for (unsigned i = 0; i < 3; i++)
    for (unsigned j  = 0; j < 3; j++)
      ASS_EQ(matrix.get(i, j), -1);
}

TEST_FUN(MatrixCheckSliceOrder01) {
  SparseMatrix<int> matrix(1, 5, 0);
  matrix.set(0, 3, 8);
  vector<pair<unsigned, int>>& nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 1);
  ASS_EQ(nonZeros[0].first, 3);
  ASS_EQ(nonZeros[0].second, 8);

  matrix.set(0, 1, 6);
  nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 2);
  ASS_EQ(nonZeros[0].first, 1);
  ASS_EQ(nonZeros[0].second, 6);
  ASS_EQ(nonZeros[1].first, 3);
  ASS_EQ(nonZeros[1].second, 8);

  matrix.set(0, 4, 9);
  nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 3);
  ASS_EQ(nonZeros[0].first, 1);
  ASS_EQ(nonZeros[0].second, 6);
  ASS_EQ(nonZeros[1].first, 3);
  ASS_EQ(nonZeros[1].second, 8);
  ASS_EQ(nonZeros[2].first, 4);
  ASS_EQ(nonZeros[2].second, 9);

  matrix.set(0, 1, 11);
  nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 3);
  ASS_EQ(nonZeros[0].first, 1);
  ASS_EQ(nonZeros[0].second, 11);
  ASS_EQ(nonZeros[1].first, 3);
  ASS_EQ(nonZeros[1].second, 8);
  ASS_EQ(nonZeros[2].first, 4);
  ASS_EQ(nonZeros[2].second, 9);
}

TEST_FUN(MatrixCheckSliceOrder02) {
  SparseMatrix<int> matrix(1, 5, 0);
  matrix.set(0, 1, 6);
  matrix.set(0, 3, 8);
  matrix.set(0, 1, 11);
  matrix.set(0, 4, 9);
  vector<pair<unsigned, int>>& nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 3);
  ASS_EQ(nonZeros[0].first, 1);
  ASS_EQ(nonZeros[0].second, 11);
  ASS_EQ(nonZeros[1].first, 3);
  ASS_EQ(nonZeros[1].second, 8);
  ASS_EQ(nonZeros[2].first, 4);
  ASS_EQ(nonZeros[2].second, 9);
}

TEST_FUN(MatrixCheckDeletion) {
  SparseMatrix<int> matrix(1, 5, -1);
  for (int i = 0; i < 5; i++)
    matrix.set(0, i, 5 + i);
  matrix.set(0, 3, -1);
  auto nonZeros = matrix.getSetOnRow(0);
  ASS_EQ(nonZeros.size(), 4);
  for (int i = 0; i < 5; i++) {
    if (i == 3)
      continue;
    unsigned index = i - (i >= 3);
    ASS_EQ(nonZeros[index].first, (VarNum) i);
    ASS_EQ(nonZeros[index].second, 5 + i);
  }
}

TEST_FUN(MatrixCheckSetGetOptimization00) {
  SparseMatrix<int> matrix(1, 6, -1);
  matrix.set(0, 1, 6);
  matrix.set(0, 5, 5);
  matrix.set(0, 3, 8);
  ASS_EQ(matrix.get(0, 3), 8);
  matrix.set(0, 5, matrix.get(0, 5) + 1);
  ASS_EQ(matrix.get(0, 5), 6);
}

TEST_FUN(MatrixCheckSetGetOptimization01) {
  SparseMatrix<int> matrix(1, 2, -1);
  matrix.set(0, 0, 1);
  ASS_EQ(matrix.get(0, 0), 1);
  matrix.set(0, 0, -1);
  ASS_EQ(matrix.get(0, 0), -1);
}


/**************************************************************************************************
 *                                   LINEAR CONSTRAINTS
 *************************************************************************************************/

#define SETUP(nPosVars, nNegVars)             \
LinearConstraint lc;                          \
vector<vector<bool>> partialOrdering;         \
Stack<pair<VarNum, Coeff>> posVars;          \
Stack<pair<VarNum, Coeff>> negVars;          \
unsigned nPos = nPosVars;                     \
unsigned nNeg = nNegVars;                     \
unsigned total = nPos + nNeg;                 \
partialOrdering.resize(total);                \
for (unsigned i = 0; i < total; i++) {        \
  partialOrdering[i] = vector(total, false);  \
}

TEST_FUN(Trivial00) {
  /**
   * c >? 0
   */
  SETUP(0, 0);

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::EQUAL);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Trivial01) {
  /**
   * 0 + X0 - X1 >? 0
   * X0 > X1
   */
  SETUP(1, 1);
  partialOrdering[0][1] = true;

  posVars.push(make_pair(0, 1));
  negVars.push(make_pair(1, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Trivial02) {
  /**
   * 0 + X0 - X1 >? 0
   * X0 > X1
   */
  SETUP(1, 1);
  partialOrdering[1][0] = true;

  posVars.push(make_pair(0, 1));
  negVars.push(make_pair(1, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Trivial03) {
  /**
   * poly c + 2 * X0 - X1
   * X0 X1
   */
  SETUP(1, 1);
  partialOrdering[0][1] = true;

  posVars.push(make_pair(0, 2));
  negVars.push(make_pair(1, -1));

  ASS_EQ(lc.getSign( 2, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1 , posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Trivial04) {
  /**
   * c + X0 + X1 >? 0
   */
  SETUP(2, 0);

  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));

  ASS_EQ(lc.getSign( 3, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 2, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Trivial05) {
  /**
   * c - X0 - X1 >? 0
   */
  SETUP(1, 1);

  negVars.push(make_pair(0, -1));
  negVars.push(make_pair(1, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-2, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-3, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Imbalanced00) {
  /**
   * poly X0 + X1 - 2 * X2
   * X0 X2
   * X1 X2
   */
  SETUP(2, 1);
  partialOrdering[0][2] = true;
  partialOrdering[1][2] = true;

  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));
  negVars.push(make_pair(2, -2));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Imbalanced01) {
  /**
   * poly X0 + X1 - 2 * X2
   * X2 X0
   * X2 X1
   */
  SETUP(2, 1);
  partialOrdering[2][0] = true;
  partialOrdering[2][1] = true;

  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));
  negVars.push(make_pair(2, -2));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Imbalanced02) {
  /**
   * poly 2 * X0 - X1 - X2
   * X0 X1
   * X0 X2
   */
  SETUP(1, 2);
  partialOrdering[0][1] = true;
  partialOrdering[0][2] = true;

  posVars.push(make_pair(0, 2));
  negVars.push(make_pair(1, -1));
  negVars.push(make_pair(2, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Imbalanced03) {
  /**
   * poly 2 * X0 - X1 - X2
   * X1 X0
   * X2 X0
   */
  SETUP(1, 2);
  partialOrdering[1][0] = true;
  partialOrdering[2][0] = true;

  posVars.push(make_pair(0, 2));
  negVars.push(make_pair(1, -1));
  negVars.push(make_pair(2, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Imbalanced04) {
  /**
   * poly 0 + X0 - X1 - X2
   * X1 X0
   * X2 X0
   */
  SETUP(1, 2);
  partialOrdering[1][0] = true;
  partialOrdering[2][0] = true;


  posVars.push(make_pair(0, 1));
  negVars.push(make_pair(1, -1));
  negVars.push(make_pair(2, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-2, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Imbalanced05) {
  /**
   * poly 1 + X0 + X1 - X2
   * X0 X2
   */
  SETUP(2, 1);
  partialOrdering[0][2] = true;


  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));
  negVars.push(make_pair(2, -1));

  ASS_EQ(lc.getSign( 2, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Real00) {
  /**
   */
  /**
   * 0 + 2X0 + 2X1 - 3Y0 - Y1
   * X0 > Y0
   * X1 > Y0
   * X1 > Y1
   */
  SETUP(2, 2);

  partialOrdering[0][2] = true;
  partialOrdering[1][2] = true;
  partialOrdering[1][3] = true;

  posVars.push(make_pair(0, 2));
  posVars.push(make_pair(1, 2));
  negVars.push(make_pair(2, -3));
  negVars.push(make_pair(3, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Symmetric00) {
  /**
   * poly X0 + X1 - X2 - X3
   * X0 X2
   * X1 X3
   */
  SETUP(2, 2);
  partialOrdering[0][2] = true;
  partialOrdering[1][3] = true;

  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));
  negVars.push(make_pair(2, -1));
  negVars.push(make_pair(3, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Symmetric01) {
  /**
   * poly X0 + X1 - X2 - X3
   * X2 X0
   * X3 X1
   */
  SETUP(2, 2);
  partialOrdering[2][0] = true;
  partialOrdering[3][1] = true;

  posVars.push(make_pair(0, 1));
  posVars.push(make_pair(1, 1));
  negVars.push(make_pair(2, -1));
  negVars.push(make_pair(3, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::LESS);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
}

TEST_FUN(Tricky00) {
  /**
   * This problems should not be preprocessed
   * poly 0 + 2 * X0 + 2 * X1 + 2 * X2 + 2 * X3 - 2 * Y0 - 3 * Y1 - 1 * Y2 - 2 * Y3
   *
   * 2    2      2    2
   * x0   x1     x2   x3
   * | \ / |   / | \ / |
   * |  X  |  /  |  X  |
   * | / \ | /   | / \ |
   * y0   y1     y2   y3
   * 2    3      1    2
   */
  SETUP(4, 4);
  partialOrdering[0][0 + 4] = true; // X0 > Y0
  partialOrdering[0][1 + 4] = true; // X0 > Y1
  partialOrdering[1][0 + 4] = true; // X1 > Y0
  partialOrdering[1][1 + 4] = true; // X1 > Y1
  partialOrdering[2][1 + 4] = true; // X2 > Y1
  partialOrdering[2][2 + 4] = true; // X2 > Y2
  partialOrdering[2][3 + 4] = true; // X2 > Y3
  partialOrdering[3][2 + 4] = true; // X3 > Y2
  partialOrdering[3][3 + 4] = true; // X3 > Y3

  posVars.push(make_pair(0, 2));
  posVars.push(make_pair(1, 2));
  posVars.push(make_pair(2, 2));
  posVars.push(make_pair(3, 2));
  negVars.push(make_pair(4, -2));
  negVars.push(make_pair(5, -3));
  negVars.push(make_pair(6, -2));
  negVars.push(make_pair(7, -1));

  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}

TEST_FUN(Tricky01) {
  /**
   * poly 0 + 2 * X4 + 2 * X2 + 2 * X5 + 2 * X1 + 2 * X0 + 2 * X3
   *  - 2 * Y0 - 3 * Y1 - 2 * Y2 - 2 * Y3 - 2 * Y4 - 1 * Y5
   * 2    2      2    2
   * x0   x1     x2   x3     x4   x5
   * | \ / |   / | \ / |   / | \ / |
   * |  X  |  /  |  X  |  /  |  X  |
   * | / \ | /   | / \ | /   | / \ |
   * y0   y1     y2   y3     y4   y5
   * 2    3      2    2      1    2
   */
  SETUP(6, 6);
  partialOrdering[0][0 + 6] = true; // X0 > Y0
  partialOrdering[0][1 + 6] = true; // X0 > Y1
  partialOrdering[1][0 + 6] = true; // X1 > Y0
  partialOrdering[1][1 + 6] = true; // X1 > Y1
  partialOrdering[2][1 + 6] = true; // X2 > Y1
  partialOrdering[2][2 + 6] = true; // X2 > Y2
  partialOrdering[2][3 + 6] = true; // X2 > Y3
  partialOrdering[3][2 + 6] = true; // X3 > Y2
  partialOrdering[3][3 + 6] = true; // X3 > Y3
  partialOrdering[4][3 + 6] = true; // X4 > Y3
  partialOrdering[4][4 + 6] = true; // X4 > Y4
  partialOrdering[4][5 + 6] = true; // X4 > Y5
  partialOrdering[5][4 + 6] = true; // X5 > Y4
  partialOrdering[5][5 + 6] = true; // X5 > Y5


  posVars.push(make_pair(0, 2));
  posVars.push(make_pair(1, 2));
  posVars.push(make_pair(2, 2));
  posVars.push(make_pair(3, 2));
  posVars.push(make_pair(4, 2));
  posVars.push(make_pair(5, 2));
  negVars.push(make_pair(0 + 6, -3));
  negVars.push(make_pair(1 + 6, -2));
  negVars.push(make_pair(2 + 6, -2));
  negVars.push(make_pair(3 + 6, -2));
  negVars.push(make_pair(4 + 6, -2));
  negVars.push(make_pair(5 + 6, -1));


  ASS_EQ(lc.getSign( 1, posVars, negVars, partialOrdering), Result::INCOMPARABLE);
  ASS_EQ(lc.getSign( 0, posVars, negVars, partialOrdering), Result::GREATER);
  ASS_EQ(lc.getSign(-1, posVars, negVars, partialOrdering), Result::GREATER);
}
