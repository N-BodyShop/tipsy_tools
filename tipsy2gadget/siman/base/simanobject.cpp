// simanobject.cpp - part of SimAn Simulation Analysis Library
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

namespace siman {

string SimanObject::className() {
  return "SimanObject";
}

list<string> SimanObject::availCommands() {
  list<string> com;
  com.push_back(string("is"));
  return com;
}


list<string> SimanObject::availMembers() {
  list<string> com;
  return com;
}

SimanObject * SimanObject::dispatch(string command, istream *stream, Scripted *pS) {
  if(command=="is") {
    cout << className() << endl;
    return NULL;
  }
  throw(UnknownCommand(command));
}

SimanObject * SimanObject::getMember(const string & var) {
  
  map<string, void*>::iterator i = valMap.find(var);
  

  throw(UnknownVariable(var));
}


void SimanObject::setMember(const string & var, SimanObject *obj) {

  throw(UnknownVariable(var));
}

unsigned int SimanObject::supports() {
  return 0;
}

bool SimanObject::supports(unsigned int flag) {
  return ((flag & supports())==flag);
}


void SimanObject::registerVal(string name, void *ptr, unsigned int type, bool readOnly) {
  valMap[name] = ptr;
  supportsMap[name] = type;
  valReadOnlyMap[name] = readOnly;
}

bool SimanObject::references(SimanObject *obj) {
  return false;
}

} // namespace siman
