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
 * @file Timer.cpp
 * Implements class Timer.
 */

#include <iostream>
#include <mutex>
#include <thread>

#include "Lib/Environment.hpp"
#include "Shell/Statistics.hpp"
#include "Shell/UIHelper.hpp"
#include "System.hpp"

#include "Timer.hpp"

using namespace std::chrono_literals;
// tradeoff: faster ticks means more overhead but more accurate LRS/timeout
const auto TICK_INTERVAL = 1ms;

// can measure instructions
#if VAMPIRE_PERF_EXISTS
// includes for checking instruction count
#include <asm/unistd.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <unistd.h>

const long long MEGA = 1 << 20;

// static data for measuring instructions
static int PERF_FD = -1; // the file descriptor we later read the info from
static long long LAST_INSTRUCTION_COUNT_READ = -1;

// conveniece wrapper around a syscall (cf. https://linux.die.net/man/2/perf_event_open )
static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags)
{
  int ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,group_fd, flags);
  return ret;
}
#endif

enum LimitType {
  TIME_LIMIT,
  INSTRUCTION_LIMIT
};

// ensures that exactly one of the timer thread and the parent process tries to exit
// (the recursive bit allows us to be overprotective
// and call disableLimitEnforcement more than once (from the main thread)
// without the consequence of locking ourselves)
static std::recursive_mutex EXIT_LOCK;

// called by timer_thread to exit the entire process
// functions called here should be thread-safe
[[noreturn]] static void limitReached(LimitType whichLimit)
{
  using namespace Shell;

  // for debugging crashes of limitReached: it is good to know what was called by vampire proper just before the interrupt
  // Debug::Tracer::printStack(cout);

  const char* REACHED[2] = {"Time limit reached!\n","Instruction limit reached!\n"};
  const char* STATUS[2] = {"% SZS status Timeout for ","% SZS status InstrOut for "};

  // if we get this lock we can assume that the parent won't also try to exit
  EXIT_LOCK.lock();

  // NB unsynchronised output:
  // probably OK as we don't output anything in other parts of Vampire during search
  reportSpiderStatus('t');
  if (outputAllowed()) {
    addCommentSignForSZS(std::cout);
    std::cout << REACHED[whichLimit];

    if (UIHelper::portfolioParent) { // the boss
      addCommentSignForSZS(std::cout);
      std::cout << "Proof not found in time ";
      Timer::printMSString(std::cout,Timer::elapsedMilliseconds());

#if VAMPIRE_PERF_EXISTS
      if (LAST_INSTRUCTION_COUNT_READ > -1) {
        std::cout << " nor after " << LAST_INSTRUCTION_COUNT_READ << " (user) instruction executed.";
      }
#endif
      std::cout << std::endl;

      if (szsOutputMode()) {
        std::cout << STATUS[whichLimit] << (env.options ? env.options->problemName().c_str() : "unknown") << std::endl;
      }
    } else // the actual child
      if (env.statistics) {
        env.statistics->print(std::cout);
    }
  }

  System::terminateImmediately(1);
}

static std::chrono::time_point<std::chrono::steady_clock> START_TIME;

// TODO could maybe be more efficient if we special-case the no-instruction-limit case:
// then, we could simply sleep until the time limit
[[noreturn]] void timer_thread()
{
  unsigned limit = env.options->timeLimitInDeciseconds();
  while(true) {
    if(limit && Timer::elapsedDeciseconds() >= limit) {
      limitReached(TIME_LIMIT);
    }

#if VAMPIRE_PERF_EXISTS
    if(env.options->instructionLimit() || env.options->simulatedInstructionLimit()) {
      Timer::updateInstructionCount();
      if (env.options->instructionLimit() && LAST_INSTRUCTION_COUNT_READ >= MEGA*(long long)env.options->instructionLimit()) {
        // in principle could have a race on terminationReason, seems unlikely/harmless in practice
        env.statistics->terminationReason = Shell::Statistics::TIME_LIMIT;
        limitReached(INSTRUCTION_LIMIT);
      }
    }
#endif

    // doesn't matter if we 'sleep in' a bit, we get the real time from the system clock
    std::this_thread::sleep_for(TICK_INTERVAL);
  }
  ASSERTION_VIOLATION
}

namespace Lib {
namespace Timer {

void reinitialise() {
  // might (probably have) locked this in the parent process, release it for the child
  //
  // I am not sure of the semantics of placement-new for std::recursive_mutex,
  // but nobody else seems to be either - if you know, tell me! - Michael
  ::new (&EXIT_LOCK) std::recursive_mutex;

  START_TIME = std::chrono::steady_clock::now();

#if VAMPIRE_PERF_EXISTS // if available, (re)initialize the perf reading
  /*
   *
   * NOTE: we need to do this before initializing the actual timer
   * (otherwise tick() could start asking the uninitialized PERF_FD!)
   */
  LAST_INSTRUCTION_COUNT_READ = -1;
  PERF_FD = -1;
  struct perf_event_attr pe;

  memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.type = PERF_TYPE_HARDWARE;
    pe.size = sizeof(struct perf_event_attr);
    pe.config = PERF_COUNT_HW_INSTRUCTIONS;
    pe.disabled = 1;
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;

  PERF_FD = perf_event_open(&pe, 0, -1, -1, 0);
  if (PERF_FD == -1) {
    // we want this to be guarded by env.options->instructionLimit()
    // not to bother with the error people who don't even know about instruction limiting
    if(env.options->instructionLimit() || env.options->simulatedInstructionLimit()) {
      std::cerr
        << "perf_event_open failed (instruction limiting will be disabled): "
        << std::strerror(errno)
        << "\n(If you are seeing 'Permission denied' ask your admin to run 'sudo sysctl -w kernel.perf_event_paranoid=-1' for you.)"
        << std::endl;
    }
  } else {
    ioctl(PERF_FD, PERF_EVENT_IOC_RESET, 0);
    ioctl(PERF_FD, PERF_EVENT_IOC_ENABLE, 0);
  }
#endif

  std::thread(timer_thread).detach();
}

void disableLimitEnforcement() {
  EXIT_LOCK.lock();
}

// return elapsed time after `START_TIME`
// must be thread-safe as it is called by the main process and timer_thread
long elapsedMilliseconds() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::steady_clock::now() - START_TIME
  ).count();
}

/**
 * Print string representing @b ms of milliseconds to @b str
 */
void printMSString(std::ostream &str, int ms)
{
  if(ms<0) {
    str << '-';
    ms = -ms;
  }

  int sec=ms/1000;
  int msonly=ms%1000;
  if(sec) {
    str<<sec;
  }
  else {
    str<<'0';
  }
  str<<'.';
  if(msonly<100) {
    str<<'0';
    if(msonly<10) {
      str<<'0';
      if(!msonly) {
	str<<'0';
      }
    }
  }
  str<<msonly<<" s";
}

std::string msToSecondsString(int ms)
{
  return Int::toString(static_cast<float>(ms)/1000)+" s";
}

bool instructionLimitingInPlace()
{
#if VAMPIRE_PERF_EXISTS
  return (PERF_FD >= 0);
#else
  return false;
#endif
}

long elapsedMegaInstructions() {
#if VAMPIRE_PERF_EXISTS
  return (LAST_INSTRUCTION_COUNT_READ >= 0) ? LAST_INSTRUCTION_COUNT_READ / MEGA : 0;
#else
  return 0;
#endif
}

// called by the main process and timer_thread: should be thread-safe
void updateInstructionCount()
{
#if VAMPIRE_PERF_EXISTS
  if (PERF_FD >= 0) {
    // we could also decide not to guard this read by env.options->instructionLimit(),
    // to get info about instructions burned even when not instruction limiting
    read(PERF_FD, &LAST_INSTRUCTION_COUNT_READ, sizeof(long long));
  }
#endif
}

} //namespace Timer
} //namespace Lib
