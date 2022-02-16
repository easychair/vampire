#ifndef SMTSUBSUMPTION_HPP
#define SMTSUBSUMPTION_HPP

#include "Kernel/Clause.hpp"
#include "Lib/STL.hpp"
#include "Inferences/InferenceEngine.hpp"

#include <memory>

namespace SMTSubsumption {


struct SubsumptionInstance
{
  Kernel::Clause* side_premise;  ///< also called "base clause"
  Kernel::Clause* main_premise;  ///< also called "instance clause"
  unsigned int number;
  int result;
};

struct SubsumptionResolutionInstance
{
  Kernel::Clause* side_premise;  ///< also called "base clause"
  Kernel::Clause* main_premise;  ///< also called "instance clause"
  /// Index into main_premise; the literal for which SR is tested
  /// - only relevant for the original SR... SMT-SR always tests all of them.
  /// - value UINT_MAX means 'all of them' (denoted as '*' in the logfile).
  unsigned res_lit;
  unsigned int number;
  int result;
};

// one round of subsumption/subsumption resolution
// (corresponds to one iteration of the saturation algorithm)
struct SubsumptionRound {
  /// index into vector of subsumption instances.
  /// one past the last instance for this round.
  unsigned s_end;
  /// index into vector of subsumption resolution instances.
  /// one past the last instance for this round.
  unsigned sr_end;
};

struct SubsumptionBenchmark {
  vvector<SubsumptionInstance> subsumptions;
  vvector<SubsumptionResolutionInstance> subsumptionResolutions;
  vvector<SubsumptionRound> rounds;
};


class SMTSubsumptionImpl2;

class SMTSubsumptionImpl3;

class ProofOfConcept {
  public:
    CLASS_NAME(ProofOfConcept);
    USE_ALLOCATOR(ProofOfConcept);

    ProofOfConcept();
    ~ProofOfConcept();

    void test(Kernel::Clause* side_premise, Kernel::Clause* main_premise);
    void benchmark_run(SubsumptionBenchmark b);
    void benchmark_micro(SubsumptionBenchmark b);

    bool checkSubsumption(Kernel::Clause* base, Kernel::Clause* instance);
    bool checkSubsumptionResolution(Kernel::Clause* base, Kernel::Clause* instance, Kernel::Clause* conclusion);

    void setupMainPremise(Kernel::Clause* instance);
    bool setupSubsumption(Kernel::Clause* base);
    bool solve();

  private:
    void add_common_benchmark_args(vvector<vstring>& args);
    void init_benchmark(vvector<vstring> the_args);

  private:
    std::unique_ptr<SMTSubsumptionImpl2> m_subsat_impl;
    std::unique_ptr<SMTSubsumptionImpl3> m_subsat_impl3;
};


}

#endif /* !SMTSUBSUMPTION_HPP */