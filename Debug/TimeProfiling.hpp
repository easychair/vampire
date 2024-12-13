/*
 * This file is part of the source code of the software program
 * Vampire. It is protected by applicable
 * copyright laws.
 *
 * This source code is distributed under the licence found here
 * https://vprover.github.io/license.html
 * and in the source directory
 */


#ifndef __TimeProfiling__
#define __TimeProfiling__

#include "Lib/Stack.hpp"
#include "Lib/Option.hpp"
#include "Debug/Assertion.hpp"
#include <chrono>
#include "Lib/MacroUtils.hpp"

namespace Shell {

#if VTIME_PROFILING

/**
 * Enables the runtime of the current block to be measured.
 * Time traces are stored in a tree starting with a root node that measures the runtime of the whole
 * probram execution, and subtrees that correspond to calls to TIME_TRACE.
 * The parameter \param name is a human-readable name to be presented when the time trace is being
 * outputed.
 * If you want two blocks to contribute to the same node in the tree, they have to occur in the same
 * subtree, and have to have the same const char* name. Note that just specifying the same literal
 * twice does not guarantee that const char* are the same, even though that is the case for some
 * compilers. Therefore one should use one static const char* for these kind of measurements.
 *
 * In addition to the time profiling tree a flattened version of the tree is being outputted at the
 * end of a program execution, that presents the sum of all the runtimes a node was encountered on
 * any branch of the tree.
 *
 * It is important to note that calls to this macro may preform allocations (i.e. the first time
 * when it is called on a certain level of the trace tree). Therefore it should be avoided in
 * recursive functions.
 * Further it should be noted that the macro introduces some overhead, hence it should also be
 * avoided to be used in parts of the codebase that are called very often and only perform short
 * tasks.
 * ```
 */
#define TIME_TRACE(name)                                                                            \
  Shell::TimeTrace::ScopedTimer CONCAT_IDENTS(__time_trace_, __LINE__) (name);

/**
 * same as TIME_TRACE, but measures the time it takes for an expression to be evaluated. the
 * arguments ... are the experssion to be time profiled.
 */
#define TIME_TRACE_EXPR(name, ...)                                                                  \
  [&](){ TIME_TRACE(name); return __VA_ARGS__; }()

/**
 * sets a new node as the root of the time trace. this is useful when launching child process.
 * The subtree will be set to the original root, when the call of TIME_TRACE_NEW_ROOT goes out of
 * scope.  Therefore the statistics must be outputted before that to see any effect of this call.
 */
#define TIME_TRACE_NEW_ROOT(name)                                                                   \
  TIME_TRACE(name)                                                                                  \
  Shell::TimeTrace::ScopedChangeRoot CONCAT_IDENTS(__change_root_, __LINE__);

#else // !VTIME_PROFILING

#define TIME_TRACE(name) {}
#define TIME_TRACE_EXPR(name, ...) __VA_ARGS__
#define TIME_TRACE_NEW_ROOT(name)

#endif // VTIME_PROFILING

#if VTIME_PROFILING

/**
 * A class to trace time for particular blocks. this class shall never be used directly,
 * as it is preprocessed away when the flag VTIME_PROFILING is set to false.
 * Use the macros TIME_TRACE, TIME_TRACE_EXPR, and TIME_TRACE_NEW_ROOT instead.
 */
class TimeTrace
{
public:
  // Let's fake a big enum like the one we used to have using a bunch of constexpr's
  // (NB: TimeTrace can only group TIME_TRACE calls with identical identifiers as pointers
  //  so always going through one place to declare a TIME_TRACE-able call site sounds like a nice routine)
  static inline constexpr const char* const CLAUSE_GENERATION = "clause generation";
  static inline constexpr const char* const CONSEQUENCE_FINDING = "consequence finding";
  static inline constexpr const char* const FMB_DEFINITION_INTRODUCTION = "fmb definition introduction";
  static inline constexpr const char* const HYPER_SUP = "hyper superposition";
  static inline constexpr const char* const LITERAL_ORDER_AFTERCHECK = "literal order aftercheck";
  static inline constexpr const char* const PARSING = "parsing";
  static inline constexpr const char* const PASSIVE_CONTAINER_MAINTENANCE = "passive container maintenance";
  static inline constexpr const char* const PREPROCESSING = "preprocessing";
  static inline constexpr const char* const PROPERTY_EVALUATION = "property evaluation";
  static inline constexpr const char* const AVATAR_SAT_SOLVER = "SAT solver";
  static inline constexpr const char* const SHUFFLING = "shuffling things";
  static inline constexpr const char* const SINE_SELECTION = "sine selection";
  static inline constexpr const char* const TERM_SHARING = "term sharing";
  
private:
  using Clock = std::chrono::steady_clock;
  using Duration = typename Clock::duration;
  using TimePoint = typename Clock::time_point;

  TimeTrace(TimeTrace     &&) = delete;
  TimeTrace(TimeTrace const&) = delete;

  class Measurements {
    Duration _sum;
    unsigned _cnt;

  public:
    void add(Duration d) {
      _cnt += 1;
      _sum += d;
    }
    void remove(Duration d) {
      _cnt -= 1;
      _sum -= d;
    }
    Duration sum() const { return _sum; }
    unsigned cnt() const { return _cnt; }
    Duration avg() const { return sum() / cnt(); }
    void extend(Measurements other) {
      _sum += other._sum;
      _cnt += other._cnt;
    }
  };


  struct Node {
    USE_ALLOCATOR(Node)
    const char* name;
    Lib::Stack<std::unique_ptr<Node>> children;
    Measurements measurements;
    Node(const char* name) : name(name), children(), measurements() {}
    struct NodeFormatOpts ;
    void printPrettyRec(std::ostream& out, NodeFormatOpts& opts);
    void printPrettySelf(std::ostream& out, NodeFormatOpts& opts);
    void serialize(std::ostream& out);
    Duration totalDuration() const;

    Node clone() const;

    Node flatten();
    void _focus(const char* name, Node& newRoot);
    Node focus(const char* name);
    void extendWith(Node const& n);
    struct FlattenState;
    void flatten_(FlattenState&);
  };

  friend std::ostream& operator<<(std::ostream& out, Duration const& self);
  TimeTrace();

public:
  static TimeTrace _instance;
  static TimeTrace& instance() { return _instance; }

  class ScopedTimer {
    TimeTrace& _trace;
#if VDEBUG
    TimePoint _start;
    const char* _name;
#endif
  public:
    ScopedTimer(const char* name);
    ScopedTimer(TimeTrace& trace, const char* name);
    ~ScopedTimer();
  };


  class ScopedChangeRoot {
    TimeTrace& _trace;
  public:
    ScopedChangeRoot();
    ScopedChangeRoot(TimeTrace& trace);
    ~ScopedChangeRoot();
  };

  void printPretty(std::ostream& out);
  void serialize(std::ostream& out);
  void setEnabled(bool);
private:

  Node _root;
  Lib::Stack<Node*> _tmpRoots;
  Lib::Stack<std::tuple<Node*, TimePoint>> _stack;
  bool _enabled;
};

#endif // VTIME_PROFILING

} // namespace Shell

#endif // __TimeProfiling__
