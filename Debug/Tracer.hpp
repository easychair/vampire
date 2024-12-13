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
 * @file Tracer.hpp
 * Implements call tracing.
 * @since 01/05/2002 Manchester
 * @since 24/10/2002 Manchester, changed after talking with Shura
 * @since 08/12/2005 Redmond, moved to the Debug namespace with the purpose
 *                   of making global to Vampire
 */


#ifndef __Tracer__
#define __Tracer__

#include <iostream>
#include <iomanip>
#include "Lib/Output.hpp"

namespace Debug {

namespace Tracer {
  // print the current stack to stdout
  void printStack();
};

template<class... As>
struct _printDbg {
  void operator()(const As&... msg);
};

template<> struct _printDbg<>{
  void operator()() { }
};

template<class A, class... As> struct _printDbg<A, As...>{
  void operator()(const A& a, const As&... as) {
    std::cout << a;
    _printDbg<As...>{}(as...);
  }
};

template<class... A> void printDbg(const char* file, int line, const A&... msg)
{
  unsigned width = 60;
  std::cout << "[ debug ] ";
  for (const char* c = file; *c != 0 && width > 0; c++, width--) {
    std::cout << *c;
  }
  for (unsigned i = 0; i < width; i++) {
    std::cout << ' ';
  }
  std::cout <<  "@" << std::setw(5) << line << ":";
  _printDbg<A...>{}(msg...);
  std::cout << std::endl; 
}

} // namespace Debug

#if VDEBUG

#ifdef ABSOLUTE_SOURCE_DIR
#define __REL_FILE__  (&__FILE__[sizeof(ABSOLUTE_SOURCE_DIR) / sizeof(ABSOLUTE_SOURCE_DIR[0])])
#else
#define __REL_FILE__  __FILE__
#endif

#  define DBG(...) { Debug::printDbg(__REL_FILE__, __LINE__, __VA_ARGS__); }
#  define WARN(...) { DBG("WARNING: ", __VA_ARGS__); }
#  define DBGE(x) DBG(#x, " = ", x)
#else // ! VDEBUG
#  define WARN(...) {}
#  define DBG(...) {}
#  define DBGE(x) {}
#endif

#endif // Tracer
