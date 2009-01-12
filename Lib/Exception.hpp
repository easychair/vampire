/**
 * Class Exception.hpp. Defines Vampire exceptions.
 *
 * @since 03/12/2003, Manchester
 */

#ifndef __Exception__
#define __Exception__

#include <string>
#include <iostream>

using namespace std;

namespace Lib {

/**
 * Abstract class Exception. It is the base class for
 * several concrete classes.
 */
class Exception
{
 public:
  /** Create an exception with a given error message */
  explicit Exception (const char* msg)
    : _message(msg)
  {}
  explicit Exception (const string msg)
    : _message(msg.c_str())
  {}
  virtual void cry (ostream&);
  virtual ~Exception() {}
 protected:
  /** Default constructor, required for some subclasses, made protected
   * so that it cannot be called directly */
  Exception () {}
  /** The error message */
  const char* _message;
}; // Exception


/**
 * Class UserErrorException. A UserErrorException is thrown
 * when a user error occurred, for example, a file name is
 * specified incorrectly; an invalid option to Vampire
 * was given, or there is a syntax error in the input file.
 */
class UserErrorException
  : public Exception
{
 public:
  UserErrorException (const char* msg)
    : Exception(msg)
  {}
  UserErrorException (const string& msg)
    : Exception(msg)
  {}
  void cry (ostream&);
}; // UserErrorException

/**
 * Class InvalidOperationException.
 */
class InvalidOperationException
  : public Exception
{
 public:
   InvalidOperationException (const char* msg)
    : Exception(msg)
  {}
   InvalidOperationException (const string& msg)
    : Exception(msg)
  {}
  void cry (ostream&);
}; // InvalidOperationException

/**
 * Class NotImplementedException.
 */
class NotImplementedException
  : public Exception
{
 public:
   NotImplementedException (const char* file,int line)
    : Exception(""), file(file), line(line)
  {}
   void cry (ostream&);
 private:
   const char* file;
   int line;
}; // InvalidOperationException


}

#define VAMPIRE_EXCEPTION \
  throw Lib::Exception(__FILE__,__LINE__)
#define USER_ERROR(msg) \
  throw Lib::UserErrorException(msg)
#define INVALID_OPERATION(msg) \
  throw Lib::InvalidOperationException(msg)
#define NOT_IMPLEMENTED \
  throw Lib::NotImplementedException(__FILE__, __LINE__)

#endif // __Exception__



