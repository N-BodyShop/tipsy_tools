//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

ostream & operator <<(ostream &os, CSimanException &e) {
  os << e.getText();
  return os;
}

CSimanException::CSimanException() {

}

string CSimanException::getText() {
  return "Unknown exception";
}

string CUnknownCommand::getText() {
  return "Unknown command " + commandName;
}

CUnknownCommand::CUnknownCommand(string command) : commandName(command) {

}


CUnknownVariable::CUnknownVariable(string variable) : variableName(variable) {

}

string CUnknownVariable::getText() {
  return "Unknown variable " + variableName;
}


CReadOnly::CReadOnly(string variable) : variableName(variable) {

}

string CReadOnly::getText() {
  return "Read only variable " + variableName;
}


CTypeError::CTypeError(string var): variableName(var) {

}

string CTypeError::getText() {
  return "Variable " + variableName + " cannot be cast to correct type";
}


CSyntaxError::CSyntaxError(string si) {
  syntax=si;
}

string CSyntaxError::getText() {
  if(syntax=="")
    return "Syntax error";
  else
    return "Syntax: " + syntax;
}

CFileError::CFileError(string fn) : fileName(fn) {

}

string CFileError::getText() {
  return "Error in file " + fileName;
}

CUnitsError::CUnitsError(const units::CUnit &u1i, const units::CUnit &u2i) : u1(u1i), u2(u2i) {
  
}

string CUnitsError::getText() {
  ostringstream oss;
  oss << "Incompatible units " << u1 << " cf " << u2;
  return oss.str();
}

CUnknownUnit::CUnknownUnit() {

}

string CUnknownUnit::getText() {
  return "Unknown internal unit specifier";
}

CMismatchedParenthesis::CMismatchedParenthesis(char type) {
  t = type;
}

string CMismatchedParenthesis::getText() {
  return "Mismatched "+t;
}
