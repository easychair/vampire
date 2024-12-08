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

using namespace std;
using namespace Lib;
using namespace Test;
using namespace Ordering;

TEST_FUN(MatrixRetrieval) {
  SparseMatrix<int> matrix(3, 3);
  matrix.set(0, 0, 1);

  ASS_EQ(matrix.get(0, 0, 0), 1);
  for (unsigned i = 0; i < 3; i++)
    for (unsigned j  = 0; j < 3; j++)
      if (i != 0 || j != 0)
        ASS_EQ(matrix.get(i, j, 0), 0);
}

TEST_FUN(MatrixCheckSliceOrder01) {
  SparseMatrix<int> matrix(1, 5);
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
  SparseMatrix<int> matrix(1, 5);
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

TEST_FUN(MatrixCheckSetGetOptimization) {
  SparseMatrix<int> matrix(1, 5);
  matrix.set(0, 1, 6);
  matrix.set(0, 5, 5);
  matrix.set(0, 3, 8);
  ASS_EQ(matrix.get(0, 3, -1), 8);
  matrix.set(0, 5, matrix.get(0, 5, -1) + 1);
  ASS_EQ(matrix.get(0, 5, -1), 6);
}


TEST_FUN(LinearConstraintPreprocessed) {
  /**
   *  0 + X0 - X1 >? 0
   *  X0 > X1
   */
  LinearConstraint lc;
  vector<vector<bool>> partialOrdering;
  partialOrdering.resize(2);
  partialOrdering[0].resize(2);
  partialOrdering[1].resize(2);
  partialOrdering[0][0] = false;
  partialOrdering[0][1] = true;
  partialOrdering[1][0] = false;
  partialOrdering[1][1] = false;


  vector<pair<VarNum, Coeff>> affineFunc;
  affineFunc.push_back(make_pair(0, 1));
  affineFunc.push_back(make_pair(1, -1));

  ASS(lc.getSign(affineFunc, partialOrdering, 0) == Greater);

  /**
   * 0 + 2X0 + 2X1 - 3Y0 - Y1
   * X0 > Y0
   * X1 > Y0
   * X1 > Y1
   */
  partialOrdering.clear();
  partialOrdering.resize(4);
  for (unsigned i = 0; i < 4; i++) {
    partialOrdering[i] = vector(4, false);
  }
  partialOrdering[0][2] = true;
  partialOrdering[1][2] = true;
  partialOrdering[1][3] = true;

  affineFunc.clear();
  affineFunc.push_back(make_pair(0, 2));
  affineFunc.push_back(make_pair(1, 2));
  affineFunc.push_back(make_pair(2, -3));
  affineFunc.push_back(make_pair(3, -1));

  ASS(lc.getSign(affineFunc, partialOrdering, 0) == Greater);

  /**
   * poly 0 + 2 * X0 + 2 * X1 + 2 * X2 + 2 * X3 - 2 * Y0 - 3 * Y1 - 1 * Y2 - 2 * Y3
   * X0 Y0
   * X0 Y1
   * X1 Y0
   * X1 Y1
   * X2 Y2
   * X2 Y3
   * X2 Y1
   * X3 Y2
   * X3 Y3
   */
  partialOrdering.clear();
  partialOrdering.resize(8);
  for (unsigned i = 0; i < 8; i++) {
    partialOrdering[i] = vector(8, false);
  }
  partialOrdering[0][4] = true;
  partialOrdering[0][5] = true;
  partialOrdering[1][4] = true;
  partialOrdering[1][5] = true;
  partialOrdering[2][5] = true;
  partialOrdering[2][6] = true;
  partialOrdering[2][7] = true;
  partialOrdering[3][6] = true;
  partialOrdering[3][7] = true;

  affineFunc.clear();
  affineFunc.push_back(make_pair(0, 2));
  affineFunc.push_back(make_pair(1, 2));
  affineFunc.push_back(make_pair(2, 2));
  affineFunc.push_back(make_pair(3, 2));
  affineFunc.push_back(make_pair(4, -2));
  affineFunc.push_back(make_pair(5, -3));
  affineFunc.push_back(make_pair(6, -2));
  affineFunc.push_back(make_pair(7, -1));

  ASS(lc.getSign(affineFunc, partialOrdering, 0) == Greater);
}
