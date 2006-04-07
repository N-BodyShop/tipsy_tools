//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CScripted
//
// A basic macro language, e.g. to describe file subsets


#ifndef __SCRIPTED_H_INCLUDED

#define __SCRIPTED_H_INCLUDED


#include "simsnap.hpp"
#include "particle.hpp"
#include "subset.hpp"
#include "geometry.hpp"
#include "union.hpp"

#include <string>
#include <stack>
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;



class CScripted : public CSimanObject {
 public:

  /// Default constructor, loads a script with given filename and executes each line through dispatch().
  CScripted(string filename, bool noInteractivity=false);
  CScripted();

  ~CScripted();

  string className();

  /// run command given by string, with optional parameters to be streamed in
  ///
  /// @param command = full line command
  CSimanObject * doCommand(string command);  
 
  void infoVars();

  CSimanObject *getReturnValue();

  /// stores details of object so it may be garbage-collected later
  /// @param pPush = pointer to object to push
  void haveCreatedObj(CSimanObject* pPush);

  /// removes details of object so no attempt to delete it is made later
  /// @param pDel = pointer to object to delete
  void haveDeletedObj(CSimanObject *pDel);

  /// get a named variable
  /// @param name - name of variable
  /// @param require - set to require ability to cast (see @ref class_type_flags); if set to 0, returns object regardless of its type
  /// @throw CUnknownVariable if the named variable is not found
  /// @throw CTypeError if the variable is found but cannot be cast according to require
  CSimanObject* getNamedVar(string name, unsigned int require=0);

  /// read from a stream, like >>, but returns things between quotes in one go
  static string readQuoted(istream &is);

  /// set a named variable
  /// @param name - name of variable
  /// @param var - pointer to object to associate with name
  /// @throw CUnknownVariable if a suitable path to the named variable is not found
  /// @throw CTypeError if the variable should not be associated with the object given
  void setNamedVar(string name, CSimanObject *var);


  /// erase a named variable
  /// @param name - name of variable
  /// @throw CUnknownVariable if a suitable path to the named variable is not found
  void eraseNamedVar(string name);

  /// split a string according to some delimiter
  /// @param[in] input - input string
  /// @param[in] delim - delimiter string separating requried output
  /// @param[out] results - vector of resulting strings 
  void splitString(const string &input, const char &delim, vector<string> &results);

  /// garbage collector
  int gc();
  bool deleteSafe(CSimanObject *pOj);

  void pollStdIn();
  void mainLoop();
 private:

  

  bool load(bool noInteractivity);

  set<CSimanObject* > createdObjs;
 
  map<string,CSimanObject*> namedStates;

  CSimanObject *pRetVal;

  string filename;

};

#endif
