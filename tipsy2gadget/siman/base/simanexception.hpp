//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __SIMANEXCEPTION_H_INCLUDED

#define __SIMANEXCEPTION_H_INCLUDED

class CSimanException {
public:
  CSimanException();
  virtual string getText();
};

class CUnknownCommand : public CSimanException {
public:
  std::string commandName;
  CUnknownCommand(std::string command);
  string getText();
  
};

class CUnknownVariable : public CSimanException {
public:
  std::string variableName;
  CUnknownVariable(std::string variable);
  string getText();
};


class CReadOnly : public CSimanException {
public:
  std::string variableName;
  CReadOnly(std::string variable);
  string getText();
};


class CTypeError : public CSimanException {
public:
  std::string variableName;
  CTypeError(std::string variable);
  string getText();
};

class CSyntaxError : public CSimanException {
public:
  CSyntaxError(string si="");
  string getText();
  string syntax;
};

class CFileError : public CSimanException {
public:
  CFileError(string fn);
  string getText();
  string fileName;
};

class CUnitsError : public CSimanException {
public:
  CUnitsError(const units::CUnit & u1i, const units::CUnit &u2i);
  string getText();
  units::CUnit u1;
  units::CUnit u2;
};

class CUnknownUnit : public CSimanException {
public:
  CUnknownUnit();
  string getText();
};

class CMismatchedParenthesis : public CSimanException {
public:
  CMismatchedParenthesis(char type);
  string getText();
  char t;
};

extern std::ostream & operator <<(std::ostream &os, CSimanException &e);

#endif
