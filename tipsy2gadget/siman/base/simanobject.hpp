//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __SIMANOBJECT_H_INCLUDED

#define __SIMANOBJECT_H_INCLUDED

class CScripted;

class CSimanObject {
public:
  virtual string className();
  virtual std::list<std::string> availCommands();
  virtual std::list<std::string> availMembers();
  virtual CSimanObject * getMember(const std::string & var);
  virtual void setMember(const std::string & var, CSimanObject *obj);
  virtual CSimanObject * dispatch(string command, std::istream *stream, CScripted *script);
  virtual bool references(CSimanObject *obj);

  /// @defgroup class_type_flags class type flags
  /// CSimanObject contains a virtual function, supports(), which returns an unsigned int
  /// containing flags indicating to which classes the run-time object can be cast
  /// This is useful at runtime, e.g. in scripting support.
  /// @{
  virtual unsigned int supports(); ///< returns a bitmask representing types to which this object may be cast
  virtual bool supports(unsigned int flag); ///< returns ((supports | flag)==true)
  
  static const unsigned int SimSnap = 1; ///< if set, this object can be cast to a CSimSnap
  static const unsigned int ColourMap = 2; ///< if set, this object can be cast to a CColourMap
  static const unsigned int Visualise = 4; ///< if set, object can be cast to a CVisualise
  static const unsigned int valueTypeInt = 8; ///< special flag meaning a CValue containing an integer; the automatic variable system extracts this integer rather than using the object
  static const unsigned int valueTypeFlt = 16; ///< special flag meaning a CValue containing a float (@see typeInt)
  static const unsigned int valueTypeDbl = 32; ///< special flag meaning a CValue containing a double (@see typeInt)
  static const unsigned int valueTypeString = 64; ///< special flag meaning a CValue containing a uninterpreted string (@see typeInt)
  /// @}

protected:
  void registerVal(string name, void *ptr, unsigned int type, bool readOnly=false);
 
  map<string, void*> valMap;
  map<string, unsigned int> supportsMap;
  map<string, bool> valReadOnlyMap;
  
};

#endif
