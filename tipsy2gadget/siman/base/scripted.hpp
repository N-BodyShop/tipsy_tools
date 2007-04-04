// scripted.hpp - part of SimAn Simulation Analysis Library
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








#ifndef __SCRIPTED_H_INCLUDED

#define __SCRIPTED_H_INCLUDED



#include "simsnap.hpp"
#include "particle.hpp"
#include "subset.hpp"
#include "geometry.hpp"
#include "union.hpp"

#include <boost/python.hpp>

#include <string>
#include <stack>
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <sstream>
#include <map>

namespace siman {


class Scripted : public SimanObject {
  friend class Config;
 public:

  /// Default constructor, loads a script with given filename
  Scripted(std::string filename, bool noInteractivity=false);
  Scripted();

  virtual ~Scripted();

  void pollStdIn();
  void mainLoop();
  
  void doCommand(std::string com);

  /// Take a SimanObject, load it into the interpreter,
  /// and make the interpreter do lifetime management (i.e.
  /// the interpreter will delete the object when it is
  /// no longer referenced by the interpreter.)
  void injectAndManage(std::string varname, SimanObject * pObj);

  /// Take a SimanObject, load it into the interpreter as
  /// variable name "varname" as a reference only. (i.e.
  /// the interperter is not allowed to delete the object
  /// if it becomes dereferenced).
  void injectAsReference(std::string varname, SimanObject & obj);

  SimanObject *getReturnValue();
  
 private:


  static void initPython();
  static void deInitPython();
  
  SimanObject *pRetVal;

  static boost::python::object main_namespace;
  static boost::python::object main_module;

  bool load(bool noInteractivity);

  std::string filename;
  
};

}
#endif
