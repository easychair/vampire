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
 * @file StringUtils.hpp
 * Defines class StringUtils.
 */

#ifndef __StringUtils__
#define __StringUtils__

#include "DHMap.hpp"
#include "Lib/Output.hpp"
#include <cstdlib>

namespace Lib {


template<class A> struct StringParser;

class StringUtils {
public:
  static std::string replaceChar(std::string str, char src, char target);
  static std::string sanitizeSuffix(std::string str);
  static bool isPositiveInteger(std::string str);
  static bool isPositiveDecimal(std::string str);
  static void replaceAll(std::string& where, const std::string& from, const std::string& to);

  static bool starts_with(const std::string& str, const std::string& what) {
    return str.rfind(what,0) == 0;
  }
  static void splitStr(const char* str, char delimiter, Stack<std::string>& strings);
  static void dropEmpty(Stack<std::string>& strings);
  static bool readEquality(const char* str, char eqChar, std::string& lhs, std::string& rhs);
  static bool readEqualities(const char* str, char delimiter, char eqChar, DHMap<std::string,std::string>& pairs);
  template<class A>
  static A parse(std::string const& str)
  { return StringParser<A>{}(str); }

  static size_t distance(const std::string &s1, const std::string &s2);
};

template<> struct StringParser<int>
{
  int operator()(std::string const& str)
  { return atoi(str.c_str()); }
};


template<> struct StringParser<float>
{
  float operator()(std::string const& str)
  { return atof(str.c_str()); }
};

template<class... Cs>
std::string outputToString(Cs const&... xs) {
  std::stringstream out;
  out << Kernel::Output::cat(xs...);
  return out.str();
}

// TODO deprecate
template<class A>
std::string toString(A const& a) {
  std::stringstream out;
  out << a;
  return out.str();
}



}

#endif // __StringUtils__
