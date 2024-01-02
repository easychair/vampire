/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */
#ifndef SUBSUMPTION_LOGGER_HPP
#define SUBSUMPTION_LOGGER_HPP

#include "Forwards.hpp"
#include "Lib/DHMap.hpp"
#include "Lib/STL.hpp"
#include "Shell/TPTPPrinter.hpp"

#include "./subsat/subsat.hpp"

#include "SATSubsumption/SATSubsumptionConfig.hpp"
#include "SATSubsumption/SubsumptionLogger.hpp"

#include <iostream>
#include <fstream>

namespace SATSubsumption {


class ForwardSubsumptionLogger {
public:
  using Clause = Kernel::Clause;
  using Literal = Kernel::Literal;
private:
  std::ofstream m_file_slog;  // *.slog
  std::ofstream m_file_tptp;  // *.p, write to this through m_tptp
  Shell::TPTPPrinter m_tptp;  // writes to m_file_tptp

  using LitIdx = uint32_t;
  using ClauseIdx = uint32_t;

  LitIdx m_next_lit = 1;      // sequence number of logged literals
  ClauseIdx m_next_clause = 1;   // sequence number of logged clauses

  Clause* m_main_premise = nullptr;

  DHMap<Literal*, LitIdx> m_logged_lits;

  struct ClauseInfo {
    ClauseIdx idx = 0;
    vvector<Literal*> lits;
    bool is_equal_to(Clause* cl) const;
  };
  DHMap<Clause*, ClauseInfo> m_logged_clauses;

  LitIdx logLiteral(Literal* lit);
  ClauseIdx logClause(Clause* cl);

public:
  /// @param filename Output file path without extension.
  ForwardSubsumptionLogger(vstring const& filename_base);

  void beginLoop(Clause* main_premise);
  void logSubsumption(Clause* side_premise, Clause* main_premise, int result);
  void logSubsumptionResolution(Clause* side_premise, Clause* main_premise, Clause* result);

  void flush();
};

struct SubsumptionInstance {
  Kernel::Clause* side_premise = nullptr;  // base
  // Kernel::Clause* main_premise = nullptr;  // instance
  int result = -999;
};

struct SubsumptionResolutionInstance {
  Kernel::Clause* side_premise = nullptr;  // base
  // Kernel::Clause* main_premise = nullptr;  // instance
  int result = -999;
};

struct ForwardSubsumptionLoop {
  // all of these should have the same main premise
  vvector<SubsumptionInstance> subsumptions;
  vvector<SubsumptionResolutionInstance> subsumption_resolutions;
  Kernel::Clause* main_premise = nullptr;  // instance

  bool empty() const {
    return subsumptions.empty() && subsumption_resolutions.empty();
  }

  void reset() {
    *this = ForwardSubsumptionLoop();
  }
};

struct SubsumptionBenchmark {
  vvector<ForwardSubsumptionLoop> fwd_loops;

  void load_from(Kernel::UnitList const* units, vstring const& slog_filename);
};


}  // namespace

#endif  // SUBSUMPTION_LOGGER_HPP