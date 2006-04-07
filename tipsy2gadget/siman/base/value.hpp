//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __VALUE_H_INCLUDED

#define __VALUE_H_INCLUDED

class CValue : public CSimanObject {
public:
  CValue();

  string className();
  unsigned int supports();
  CSimanObject * dispatch(string command, istream *stream, CScripted *pS);

  /// construct a CValue from the contents of val
  /// @param tryInterpret - if true (default), try to interpret this as input for a different type
  /// @param forceInterpret - if true (not default), throw CTypeError if this does not qualify as input for a different type
  CValue(string val, bool tryInterpret=true, bool forceIntepret=false);
  CValue(int intv);
  CValue(double doublev);

  /// get double value of this object. Always throws CTypeError if this cannot be interpreted as a double.
  /// @param throwIfNotFP - if true (not default) throw CTypeError if this was actually an integer.
  double getDouble(bool throwIfNotFP=false);

  
  /// get double value of this object. Always throws CTypeError if this cannot be interpreted as an integer.
  /// @param throwIfFP - if true (not default) throw CTypeError if this was actually floating point and has been truncated.
  int getInt(bool throwIfFP=false);

  /// get string value of this object. Always possible.
  string getString();
  

private:
  string stringData;
  double doubleData;
  int intData;

  
  unsigned int type;

  static bool isNumeric(string d);
};

#endif
