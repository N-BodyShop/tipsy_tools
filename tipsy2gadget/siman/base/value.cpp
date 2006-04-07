//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//

#include "siman.hpp"

CValue::CValue() {

}

string CValue::className() {
  return "CValue";
}

unsigned int CValue::supports() {
  if((type & valueTypeDbl)!=0)
    return type|valueTypeFlt; // value never actually stored as float
  else
    return type;
}


CSimanObject * CValue::dispatch(string command, istream *stream, CScripted *pS) {
  if(command=="is") {
    cout.flush();
    if((type & valueTypeInt)!=0)
      cout << "(CValue integer) " << intData << endl;
    else if((type & valueTypeDbl)!=0)
      cout << "(CValue floating point) " << doubleData << endl;
    else
      cout << "(CValue string) " << stringData << endl;
    return NULL;
  }
  throw(CUnknownCommand(command));
}


CValue::CValue(string val, bool tryInterpret, bool mustInterpret) {
  stringData = val;
  type=valueTypeString;
  if(tryInterpret && isNumeric(val)) {
    
    if(val.find(".",0)!=string::npos) {
      // cannot be an integer
      doubleData = atof(val.c_str());
      type=valueTypeDbl;
    } else {

      intData = atoi(val.c_str());
      type=valueTypeInt;
    } 
  } else {
    if(mustInterpret)
      throw(CTypeError("unnamed"));
  }
}

CValue::CValue(int val) {
  intData = val;
  type=valueTypeInt;
  
  ostringstream ss(stringData);
  ss << val;
}

CValue::CValue(double val) {
  doubleData = val;
  type=valueTypeDbl;
 
  ostringstream ss(stringData);
  ss << val;
}

double CValue::getDouble(bool throwIfNotFP) {
  
  if((type & valueTypeDbl)!=0)
    return doubleData;
  
  if((type & valueTypeInt)!=0 && !throwIfNotFP)
    return (double) intData;

  // fallthrough

  throw CTypeError("unnamed");

}

int CValue::getInt(bool throwIfFP) {
  
  if((type & valueTypeInt)!=0)
    return intData;
  
  if((type & valueTypeDbl)!=0 && !throwIfFP)
    return (int) doubleData;

  // fallthrough

  throw CTypeError("unnamed");

}

string CValue::getString() {
  return stringData;
}

bool CValue::isNumeric(string det) {
  string::iterator i;
  for(i=det.begin();i!=det.end();i++) {
    if((*i<'0' || *i>'9') && *i!='.')
      return false;
  }
  return true;
}
