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
 * @file Shell/Lexer.cpp
 * Implements class Lexer
 *
 * @since 27/07/2004 Torrevieja
 */

#include "Debug/Assertion.hpp"

#include "Lib/Int.hpp"

#include "Lexer.hpp"

using namespace std;
using namespace Shell;
using namespace Lib;

/**
 * Initialise a lexer.
 * @since 27/07/2004 Torrevieja
 */
Lexer::Lexer (istream& in)
  : _charBuffer(512),
    _charCursor(0),
    _stream(in),
    _eof(false),
    _lineNumber(1),
    _lookAheadChar(0)
{
  readNextChar();
} // Lexer::Lexer


/**
 * Reads next character into _lastCharacter.
 *
 * @return true if such a character exists
 * @since 14/07/2004 Turku
 */
bool Lexer::readNextChar ()
{
  if (_lookAheadChar) {
    _lastCharacter = _lookAheadChar;
    _lookAheadChar = 0;
    if (_lastCharacter == -1) {
      _eof = true;
      return false;
    }
    return true;
  }

  if (_eof) {
    return false;
  }

  _lastCharacter = _stream.get();
  if (_lastCharacter == -1) {
    _eof = true;
    return false;
  }

  if (_lastCharacter == '\n') {
    _lineNumber++;
  }
  return true;
} // Lexer::readNextChar


/**
 * Read a number (either integer or floating point).
 * @since 15/07/2004 Turku
 */
void Lexer::readNumber (Token& token)
{
  if (_lastCharacter == '-') {
    saveLastChar();
    readNextChar();
  }
  readUnsignedInteger();
  if (_lastCharacter == '.') {
    saveLastChar();
    if (readNextChar()) {
      if (isDigit(_lastCharacter)) {
	readUnsignedInteger();
	token.tag = TT_REAL;
	saveTokenText(token);
	return;
      }
    }
    saveTokenText(token);
    throw LexerException((std::string)"incorrect number format in " + token.text,
			 *this);
  }
  token.tag = TT_INTEGER;
  saveTokenText(token);
} // Lexer::readNumber


/**
 * Save last character in the character buffer (if capacity permits).
 * @since 15/07/2004 Turku
 */
void Lexer::saveLastChar ()
{
  _charBuffer[_charCursor++]  = (char)_lastCharacter;
} // Lexer::saveLastChar


/**
 * Save character in the character buffer (if capacity permits).
 * @since 25/08/2004 Torrevieja
 */
void Lexer::saveChar (int character)
{
  _charBuffer[_charCursor++]  = (char)character;
} // Lexer::saveChar


/**
 * Save the content of the current character buffer as the text
 * of the given token.
 * @since 15/07/2004 Turku
 */
void Lexer::saveTokenText (Token& token)
{
  _charBuffer[_charCursor] = 0;
  token.text = _charBuffer.content();
} // Lexer::saveTokenText


/**
 * Read an unsigned integer.
 * @since 15/07/2004 Turku
 */
void Lexer::readUnsignedInteger ()
{
  saveLastChar();

  while (readNextChar() && isDigit(_lastCharacter)) {
    saveLastChar();
  }
} // Lexer::readUnsignedInteger

/**
 * Create a new lexer exception.
 * @since 15/07/2004 Turku
 */
LexerException::LexerException (std::string message,const Lexer& lexer)
  : _message (message)
{
  if (lexer.isAtEndOfFile()) {
    _message += " at end of input";
    return;
  }
  int line = lexer.lineNumber();
  if (lexer.lastCharacter() == '\n') {
    line--;
  }
  _message += " in line ";
  _message += Int::toString(line);
} // LexerException::LexerException


/**
 * Write itself to an ostream.
 * @since 15/07/2004 Turku
 */
void LexerException::cry (std::ostream& out) const
{
  out << "Lexer exception: " << _message << '\n';
} // LexerException::LexerException


/**
 * Read a sequence of characters cs without saving them.
 * @since 25/08/2004 Torrevieja
 */
void Lexer::readSequence (const char* cs)
{
  while (*cs) {
    readNextChar();
    if (lastCharacter() != *cs) {
      throw LexerException((std::string)cs + 
			   " expected",*this);
    }
    cs++;
  }
  readNextChar();
} // Lexer::readSequence


/**
 * Look ahead one character and return it.
 * @since 27/11/2006 Haifa
 */
int Lexer::lookAhead()
{
  ASS(! _lookAheadChar); // cannot look ahead by two characters!

  _lookAheadChar = _stream.get();
  return _lookAheadChar;
} // Lexer::lookAhead()

