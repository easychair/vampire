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
 * @file Assertion.hpp
 * Defines several functions that
 * could be used to make assertions.
 */

#ifndef __Assertion__
#define __Assertion__

#if VDEBUG
#include <iostream>
#include "Tracer.hpp"

namespace Shell {
bool outputAllowed(bool debug);
}

namespace Debug {
namespace Assertion {
[[noreturn]] void abortAfterViolation();
[[noreturn]] void violated(const char* file, int line, const char* condition);

template <typename T>
[[noreturn]] void violated(const char* file, int line, const char* condition,
                           const T& rep, const char* repStr);

template <typename T, typename U>
[[noreturn]] void violated(const char* file, int line, const char* condition,
                           const T& rep, const char* repStr, const U& rep2, const char* repStr2);

template <typename T, typename U>
[[noreturn]] void violatedEquality(const char* file, int line, const char* val1Str,
                                   const char* val2Str, const T& val1, const U& val2);

template <typename T, typename U>
[[noreturn]] void violatedNonequality(const char* file, int line, const char* val1Str,
                                      const char* val2Str, const T& val1, const U& val2);

template <typename T, typename U>
[[noreturn]] void violatedComparison(const char* file, int line, const char* val1Str,
                                     const char* val2Str, const T& val1, const U& val2, bool strict, bool greater);

template <typename T>
[[noreturn]] void violatedMethod(const char* file, int line, const T& obj,
                                 const char* objStr, const char* methodStr, const char* prefix);

[[noreturn]] void violatedStrEquality(const char* file, int line, const char* val1Str,
                                      const char* val2Str, const char* val1, const char* val2);

[[noreturn]] void reportAssertValidException(const char* file, int line, const char* obj);
} // namespace Assertion
} // namespace Debug

#define ASS(Cond)                                            \
  {                                                          \
    if (!(Cond)) {                                           \
      Debug::Assertion::violated(__FILE__, __LINE__, #Cond); \
    }                                                        \
  }

#define ASS_REP(Cond, ReportedVal)                                                      \
  {                                                                                     \
    if (!(Cond)) {                                                                      \
      Debug::Assertion::violated(__FILE__, __LINE__, #Cond, ReportedVal, #ReportedVal); \
    }                                                                                   \
  }

#define ASS_REP2(Cond, ReportedVal, ReportedVal2)                                                                    \
  {                                                                                                                  \
    if (!(Cond)) {                                                                                                   \
      Debug::Assertion::violated(__FILE__, __LINE__, #Cond, ReportedVal, #ReportedVal, ReportedVal2, #ReportedVal2); \
    }                                                                                                                \
  }

#define ASS_EQ(VAL1, VAL2)                                                              \
  {                                                                                     \
    if (!((VAL1) == (VAL2))) {                                                          \
      Debug::Assertion::violatedEquality(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2); \
    }                                                                                   \
  }

#define ASS_NEQ(VAL1, VAL2)                                                                \
  {                                                                                        \
    if (!((VAL1) != (VAL2))) {                                                             \
      Debug::Assertion::violatedNonequality(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2); \
    }                                                                                      \
  }

#define ASS_STR_EQ(VAL1, VAL2)                                                             \
  {                                                                                        \
    if (strcmp((VAL1), (VAL2))) {                                                          \
      Debug::Assertion::violatedStrEquality(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2); \
    }                                                                                      \
  }

#define ASS_G(VAL1, VAL2)                                                                             \
  {                                                                                                   \
    if (!((VAL1) > (VAL2))) {                                                                         \
      Debug::Assertion::violatedComparison(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2, true, true); \
    }                                                                                                 \
  }

#define ASS_L(VAL1, VAL2)                                                                              \
  {                                                                                                    \
    if (!((VAL1) < (VAL2))) {                                                                          \
      Debug::Assertion::violatedComparison(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2, true, false); \
    }                                                                                                  \
  }

#define ASS_GE(VAL1, VAL2)                                                                             \
  {                                                                                                    \
    if (!((VAL1) >= (VAL2))) {                                                                         \
      Debug::Assertion::violatedComparison(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2, false, true); \
    }                                                                                                  \
  }

#define ASS_LE(VAL1, VAL2)                                                                              \
  {                                                                                                     \
    if (!((VAL1) <= (VAL2))) {                                                                          \
      Debug::Assertion::violatedComparison(__FILE__, __LINE__, #VAL1, #VAL2, VAL1, VAL2, false, false); \
    }                                                                                                   \
  }


#define ASS_ALLOC_TYPE(PTR, TYPE)

#define ASS_METHOD(OBJ, METHOD)                                                       \
  {                                                                                   \
    if (!((OBJ).METHOD)) {                                                            \
      Debug::Assertion::violatedMethod(__FILE__, __LINE__, (OBJ), #OBJ, #METHOD, ""); \
    }                                                                                 \
  }

#define ASSERT_VALID(obj)                                                     \
  {                                                                           \
    try {                                                                     \
      (obj).assertValid();                                                    \
    }                                                                         \
    catch (...) {                                                             \
      Debug::Assertion::reportAssertValidException(__FILE__, __LINE__, #obj); \
    }                                                                         \
  }

#define ASSERTION_VIOLATION                                 \
  {                                                         \
    Debug::Assertion::violated(__FILE__, __LINE__, "true"); \
  }

#define ASSERTION_VIOLATION_REP(Val) ASS_REP(false, Val)

#define ASSERTION_VIOLATION_REP2(Val1, Val2) ASS_REP2(false, Val1, Val2)

#define ASS_NO_EXCEPT(...) \
  try { __VA_ARGS__ }\
  catch (Exception& e) { e.cry(std::cout); ASSERTION_VIOLATION } \
  catch (...)          {                   ASSERTION_VIOLATION } \

#define RELEASE_CODE(X) {}
#define DEBUG_CODE(X) X
#define ALWAYS(Cond) ASS(Cond)
#define NEVER(Cond) ASS(!(Cond))
#else // ! VDEBUG

/* __UNREACHABLE: this point in the code is statically unreachable */
#if defined(__GNUC__)
// __builtin_unreachable if the GNU dialect is availble: GCC, Clang, ICC
#define __UNREACHABLE        \
  {                          \
    __builtin_unreachable(); \
  }
#elif defined(_MSC_VER)
// __assume(0) if Microsoft-y
#define __UNREACHABLE \
  {                   \
    __assume(0);      \
  }
#endif
#ifndef __UNREACHABLE
// otherwise, infinite loop - UB and should be optimised out
#define __UNREACHABLE \
  {                   \
    while (true) {    \
    }                 \
  }
#endif

#define RELEASE_CODE(X) X
#define DEBUG_CODE(X) {}

#define ASS(Cond)  {}
#define ALWAYS(Cond) if (!(Cond)) __UNREACHABLE;
#define NEVER(Cond) if ((Cond)) __UNREACHABLE;

#define ASS_REP(Cond, ReportedVal) {}
#define ASS_REP2(Cond, ReportedVal, ReportedVal2) {}

#define ASS_EQ(VAL1,VAL2) {}
#define ASS_NEQ(VAL1,VAL2) {}
#define ASS_STR_EQ(VAL1,VAL2) {}

#define ASS_G(VAL1,VAL2) {}
#define ASS_L(VAL1,VAL2) {}
#define ASS_GE(VAL1,VAL2) {}
#define ASS_LE(VAL1,VAL2) {}

#define ASS_ALLOC_TYPE(PTR,TYPE) {}
#define ASS_METHOD(OBJ,METHOD) {}

#define ASSERTION_VIOLATION __UNREACHABLE
#define ASSERTION_VIOLATION_REP(Val) ASSERTION_VIOLATION
#define ASSERTION_VIOLATION_REP2(Val1,Val2)  ASSERTION_VIOLATION

#define ASSERT_VALID(obj) {}

#define ASS_NO_EXCEPT(...) __VA_ARGS__

#endif // VDEBUG


#if VDEBUG

template <typename T>
void Debug::Assertion::violated(const char* file, int line, const char* cond,
                                const T& rep, const char* repStr)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition in file " << file << ", line " << line
         << " violated:\n"
         << cond << "\n"
         << "Value of " << repStr << " is: " << rep
         << "\n----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violated

template <typename T, typename U>
void Debug::Assertion::violated(const char* file, int line, const char* cond,
                                const T& rep, const char* repStr, const U& rep2, const char* repStr2)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition in file " << file << ", line " << line
         << " violated:\n"
         << cond << "\n"
         << "Value of " << repStr << " is: " << rep << "\n"
         << "Value of " << repStr2 << " is: " << rep2
         << "\n----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violated

template <typename T, typename U>
void Debug::Assertion::violatedEquality(const char* file, int line, const char* val1Str,
                                        const char* val2Str, const T& val1, const U& val2)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition " << val1Str << " == " << val2Str << " in file " << file << ", line " << line
              << " was violated, as:\n"
              << val1Str << " == " << val1 << "\n"
              << val2Str << " == " << val2 << "\n"
              << "----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violatedEquality

template <typename T, typename U>
void Debug::Assertion::violatedNonequality(const char* file, int line, const char* val1Str,
                                           const char* val2Str, const T& val1, const U& val2)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition " << val1Str << " != " << val2Str << " in file " << file << ", line " << line
              << " was violated, as:\n"
              << val1Str << " == " << val1 << "\n"
              << val2Str << " == " << val2 << "\n"
              << "----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violatedNonequality

template <typename T, typename U>
void Debug::Assertion::violatedComparison(const char* file, int line, const char* val1Str,
                                          const char* val2Str, const T& val1, const U& val2, bool strict, bool greater)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition " << val1Str;
    if (strict) {
      if (greater) {
        std::cout << " > ";
      }
      else {
        std::cout << " < ";
      }
    }
    else if (greater) {
      std::cout << " >= ";
    }
    else {
      std::cout << " <= ";
    }

    std::cout << val2Str << " in file " << file << ", line " << line
              << " was violated, as:\n"
              << val1Str << " == " << val1 << "\n"
              << val2Str << " == " << val2 << "\n"
              << "----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violatedComparison

template <typename T>
void Debug::Assertion::violatedMethod(const char* file, int line, const T& obj,
                                      const char* objStr, const char* methodStr, const char* prefix)
{
  if (Shell::outputAllowed(true)) {
    std::cout << "Condition " << prefix << "(" << objStr << ")." << methodStr << " in file "
              << file << ", line " << line << " was violated for:\n"
              << objStr << " == " << obj << "\n"
              << "----- stack dump -----\n";
    Tracer::printStack();
    std::cout << "----- end of stack dump -----\n";
  }
  abortAfterViolation();
} // Assertion::violatedMethod

#endif // VDEBUG

/** expression version of ASSERTION_VIOLATION */
template<class T> T assertionViolation()
{ ASSERTION_VIOLATION }

#if defined(__clang__)
#  define __ALLOW_UNUSED(...)                                                             \
    _Pragma("GCC diagnostic push")                                                        \
    _Pragma("GCC diagnostic ignored \"-Wunused\"")                                        \
    __VA_ARGS__                                                                           \
    _Pragma("GCC diagnostic pop")                                                         \

#elif defined(__GNUC__) || defined(__GNUG__)

#  define __ALLOW_UNUSED(...)                                                             \
    _Pragma("GCC diagnostic push")                                                        \
    _Pragma("GCC diagnostic ignored \"-Wunused-but-set-variable\"")                       \
    __VA_ARGS__                                                                           \
    _Pragma("GCC diagnostic pop")                                                         \

#else
#  define __ALLOW_UNUSED(...) __VA_ARGS__
#endif
 
#endif // __Assertion__
