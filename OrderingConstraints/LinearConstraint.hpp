#pragma once

#include <vector>

#include "SparseMatrix.hpp"
#include "Kernel/OrderingComparator.hpp"

namespace Lib
{
  class LinearConstraint
  {
    public:
    using Flow = unsigned;
    using VarNum = unsigned;
    using VarAlias = unsigned;
    using Coeff = int;
    using Constant = int;

    private:
      /**
       * @brief contains the alias of a variable. If posVars[i] = x, then i is the alias of x in this part of the code.
       * @details used to interface with the caller
       * @details the vector only contains the aliases for positively occuring variables in the affine constraint.
       */
      // std::vector<Variable> posVars;
      // std::vector<Variable> negVars;
      unsigned nPosVars = 0;
      unsigned nNegVars = 0;

      std::vector<Flow> capacities;

      std::vector<Flow> requirements;

      Constant constant;

      int sumCoeffs = 0;
      bool inverted = false;

      unsigned allocated_size = 0;


      /**
       * The element flowMatrix[i][j] contains the flow from posVars[i] to negVars[j]
       */
      Lib::SparseMatrix<Flow> flowMatrix;

      // The element greaterThanY[i][j] is true if i is greater than j
      Lib::SparseMatrix<bool> greaterThanY;

      inline unsigned index(VarAlias x, VarAlias y) { return x * nNegVars + y;}

      void reset();


      bool pruneLevel0();
      bool pruneLevel1();

      bool preProcess();

      std::vector<bool> seenX;
      std::vector<bool> seenY;

      bool dfsX(VarAlias x, std::vector<VarAlias>& path);

      bool dfsY(VarAlias y, std::vector<VarAlias>& path);

      Flow findPath(VarAlias sink, std::vector<VarAlias>& path);

      bool search();

      Result solve();

      void removeXVariable(VarAlias xAlias);

      void removeYVariable(VarAlias yAlias);

      void setProblem(const Stack<std::pair<VarNum, Coeff>>& posVars,
                      const Stack<std::pair<VarNum, Coeff>>& negVars);

      void setOrdering(const Stack<std::pair<VarNum, Coeff>>& posVars,
                       const Stack<std::pair<VarNum, Coeff>>& negVars,
                       const Kernel::TermPartialOrdering* partialOrdering);

      void setOrdering(const Stack<std::pair<VarNum, Coeff>>& posVars,
                       const Stack<std::pair<VarNum, Coeff>>& negVars,
                       const std::vector<std::vector<bool>>& partialOrdering);
    public:
      LinearConstraint();

      std::string to_string() const;

      /**
       * @brief returns a comparison between the affine function and the constant @p c
       */
      Result getSign(const Constant& c,
                     const Lib::Stack<std::pair<VarNum,Coeff>>& posVars,
                     const Lib::Stack<std::pair<VarNum,Coeff>>& negVars,
                     const Kernel::TermPartialOrdering* partialOrdering);

      /**
       * @brief returns a comparison between the affine function and the constant @p c
       * @details the points of this method is for unit testing
       */
      Result getSign(const Constant& c,
                     const Lib::Stack<std::pair<VarNum,Coeff>>& posVars,
                     const Lib::Stack<std::pair<VarNum,Coeff>>& negVars,
                     const std::vector<std::vector<bool>> partialOrdering);


  };
}
