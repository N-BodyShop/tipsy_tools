// simanexception.cpp - part of SimAn Simulation Analysis Library
//
//
// Copyright (c) Andrew Pontzen 2005, 2006
//
// SimAn is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// SimAn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public Licence for more details.
//
// You should have received a copy of the GNU General Public Licence
// along with SimAn; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA









#include "siman.hpp"
#include <boost/lexical_cast.hpp>

namespace siman {

ostream & operator <<(ostream &os, SimanException &e) {
  os << e.getText();
  return os;
}

SimanException::SimanException() : d("") {

}

SimanException::SimanException(const string & desc)  : d(desc) {

}

  OutOfRange::OutOfRange(int n) : n(n) { }

  string OutOfRange::getText() const throw() {
    return "Index "+boost::lexical_cast<string>(n)+" out of range";
  }

string SimanException::getText() const throw() {
  if(d=="")
    return "Unknown Siman exception";
  else
    return d;
}

const char* SimanException::what() const throw() {
  return getText().c_str();
}

  ConvergenceFailure::ConvergenceFailure() : SimanException("Convergence Failure") {
  }
  
  ConvergenceFailure::ConvergenceFailure(const string & c) : SimanException("Convergence Failure: "+c) {
  }

string UnknownCommand::getText() const throw() {
  return "Unknown command " + commandName;
}

UnknownCommand::UnknownCommand(string command) : commandName(command) {

}


UnknownArray::UnknownArray(string variable) : variableName(variable) {

}

string UnknownArray::getText() const throw() {
  return "Unknown array " + variableName;
}


ReadOnly::ReadOnly(string variable) : variableName(variable) {

}

string ReadOnly::getText() const throw() {
  return "Read only variable " + variableName;
}


TypeError::TypeError(string var): variableName(var) {

}

string TypeError::getText() const throw() {
  return "Variable " + variableName + " cannot be cast to correct type";
}


SyntaxError::SyntaxError(string si) {
  syntax=si;
}

string SyntaxError::getText() const throw() {
  if(syntax=="")
    return "Syntax error";
  else
    return "Syntax: " + syntax;
}

FileError::FileError(string fn) : fileName(fn) {

}

string FileError::getText() const throw() {
  return "Error in file " + fileName;
}

UnitsError::UnitsError(const Unit &u1i, const Unit &u2i) : u1(u1i), u2(u2i) {
  
}

string UnitsError::getText() const throw() {
  ostringstream oss;
  oss << "Incompatible units " << u1 << " cf " << u2;
  return oss.str();
}

  UnknownUnit::UnknownUnit(string sn) : n(sn) {

}

string UnknownUnit::getText() const throw() {
  if(n=="")
    return "Unknown internal unit specifier";
  else
    return "Unknown unit "+n;
}

MismatchedParenthesis::MismatchedParenthesis(char type) {
  t = type;
}

string MismatchedParenthesis::getText() const throw() {
  return "Mismatched "+t;
}

} // namespace siman
