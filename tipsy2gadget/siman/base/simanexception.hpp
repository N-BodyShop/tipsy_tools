// simanexception.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __SIMANEXCEPTION_H_INCLUDED

#define __SIMANEXCEPTION_H_INCLUDED

namespace siman {

  class SimanException : public std::exception {
  public:
    ~SimanException() throw() { }; 
    SimanException();
    SimanException(const std::string & desc);
    virtual std::string getText() const throw();
    virtual const char* what() const throw();
  private:
    std::string d;
  };

  class ConvergenceFailure : public SimanException {
  public:
    ~ConvergenceFailure() throw() { };
    ConvergenceFailure();
    ConvergenceFailure(const std::string &desc);
  };

  class OutOfRange : public SimanException {
  public:
    ~OutOfRange() throw() { } ;
    int n;
    OutOfRange(int n);
    std::string getText() const throw();
  };

  class UnknownCommand : public SimanException {
  public:
    ~UnknownCommand() throw() { }; 
    std::string commandName;
    UnknownCommand(std::string command);
    std::string getText() const throw();
  
  };

  class MismatchedArrayLength : public SimanException {
  public:
    ~MismatchedArrayLength() throw() { };
    MismatchedArrayLength() { };
    std::string getText() const throw() { return "Mismatched array lengths"; };

  };

  class UnknownArray : public SimanException {
  public:
    ~UnknownArray() throw() { }; 
    std::string variableName;
    UnknownArray(std::string variable);
    std::string getText() const throw();
  };


  class ReadOnly : public SimanException {
  public:
    ~ReadOnly() throw() { }; 
    std::string variableName;
    ReadOnly(std::string variable);
    std::string getText() const throw();
  };


  class TypeError : public SimanException {
  public:
    ~TypeError() throw() { }; 
    std::string variableName;
    TypeError(std::string variable);
    std::string getText() const throw();
  };

  class SyntaxError : public SimanException {
  public:
    ~SyntaxError() throw() { }; 
    SyntaxError(std::string si="");
    std::string getText() const throw();
    std::string syntax;
  };

  class FileError : public SimanException {
  public:
    ~FileError() throw() { }; 
    FileError(std::string fn);
    std::string getText() const throw();
    std::string fileName;
  };

  class UnitsError : public SimanException {
  public:
    ~UnitsError() throw() { }; 
    UnitsError(const Unit & u1i, const Unit &u2i);
    std::string getText() const throw();
    Unit u1;
    Unit u2;
  };

  class UnknownUnit : public SimanException {
  public:
    ~UnknownUnit() throw() { }; 
    UnknownUnit(std::string name="");
    std::string getText() const throw();
    std::string n;
  };

  class MismatchedParenthesis : public SimanException {
  public:
    ~MismatchedParenthesis() throw() { }; 
    MismatchedParenthesis(char type);
    std::string getText() const throw();
    char t;
  };

  extern std::ostream & operator <<(std::ostream &os, SimanException &e);

}

#endif
