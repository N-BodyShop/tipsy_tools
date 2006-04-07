//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

string CSimanObject::className() {
  return "CSimanObject";
}

list<string> CSimanObject::availCommands() {
  list<string> com;
  com.push_back(string("is"));
  return com;
}


list<string> CSimanObject::availMembers() {
  list<string> com;
  return com;
}

CSimanObject * CSimanObject::dispatch(string command, istream *stream, CScripted *pS) {
  if(command=="is") {
    cout << className() << endl;
    return NULL;
  }
  throw(CUnknownCommand(command));
}

CSimanObject * CSimanObject::getMember(const string & var) {
  
  map<string, void*>::iterator i = valMap.find(var);
  
  if(i!=valMap.end()) {
    
   
    if((supportsMap[var] & valueTypeInt)==valueTypeInt) {
      CValue *v = new CValue(*(static_cast<int*>((*i).second)));
      return v;
    }


    if((supportsMap[var] & valueTypeFlt)==valueTypeFlt) {
      
      CValue *v = new CValue((double) (*(static_cast<float*>((*i).second))));
      return v;
    }


    if((supportsMap[var] & valueTypeDbl)==valueTypeDbl) {
            
      CValue *v = new CValue((*(static_cast<double*>((*i).second))));
      return v;
    }

  }


  throw(CUnknownVariable(var));
}


void CSimanObject::setMember(const string & var, CSimanObject *obj) {

  map<string, void*>::iterator i = valMap.find(var);
  
  if(i!=valMap.end()) {
    if(valReadOnlyMap[var])
      throw(CReadOnly(var));

    if(!obj->supports(supportsMap[var]))
      throw(CTypeError(var));
    
    

    if((supportsMap[var] & valueTypeInt)==valueTypeInt) {
      int *v = static_cast<int*>((*i).second);
      *v = (dynamic_cast<CValue*>(obj))->getInt();
    } else if((supportsMap[var] & valueTypeFlt)==valueTypeFlt) {
      float *v = static_cast<float*>((*i).second);
      *v = (float) ((dynamic_cast<CValue*>(obj))->getDouble());
    } else if((supportsMap[var] & valueTypeDbl)==valueTypeDbl) {
      double *v = static_cast<double*>((*i).second);
      *v = ((dynamic_cast<CValue*>(obj))->getDouble());
    } else {
      cerr << "No copy routine written yet for this object, sorry." << endl;
    }

    return;

  }

  // now we've reached the end of the line!

  throw(CUnknownVariable(var));
}

unsigned int CSimanObject::supports() {
  return 0;
}

bool CSimanObject::supports(unsigned int flag) {
  return ((flag & supports())==flag);
}


void CSimanObject::registerVal(string name, void *ptr, unsigned int type, bool readOnly) {
  valMap[name] = ptr;
  supportsMap[name] = type;
  valReadOnlyMap[name] = readOnly;
}

bool CSimanObject::references(CSimanObject *obj) {
  return false;
}
