//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// ColumnList


#ifndef __COLUMNLIST_H_INCLUDED

#define __COLUMNLIST_H_INCLUDED

class CColumnList {
  friend class CGrid;

public:
  CColumnList(CSimSnap *parentSim, float x1, float x2, int nxi, bool assign=true);

  float power(float wavenumber);

  ~CColumnList();

  virtual void realize();

  virtual CSimSnap * operator[](int index);

protected:
  CColumnList();

private:
  int nx;
  float dx;
  CSimSnap **pColumnContents;
};


template <typename CVirtualClass, typename CReturnClass>
class CTemp : public CVirtualClass {

  // workaround class - to allow a dereferenced CColumnList
  // to be returned from [] operations on CColumnGrid, but
  // prevent the object being destructed higher up the
  // stack...

  // One could, of course, make a copy constructor for
  // CColumnList, but this is actually quite expensive if
  // there are a lot of accesses going on, i.e.
  //
  //   pipeRow[x][y] 
  // 
  // would be an O(n) operation, so any 2D operation would
  // become O(n^3)!

  friend class CGrid; // shame, but used in constructor. sorry.

public:
  CTemp(CVirtualClass *refTo);
  ~CTemp();
  
  void realize();
  CReturnClass operator[](int index);
  
private:
  CVirtualClass *realV;
  
};

#define CTempColumnList CTemp<CColumnList, CSimSnap*>

template <typename CVirtualClass, typename CReturnClass> 
CTemp<CVirtualClass,CReturnClass>::CTemp(CVirtualClass *orig) {
  realV = orig;
}

template <typename CVirtualClass, typename CReturnClass> 
CReturnClass CTemp<CVirtualClass, CReturnClass>::operator[](int n) {
  return (*realV)[n];
}

template <typename CVirtualClass, typename CReturnClass> 
CTemp<CVirtualClass, CReturnClass>::~CTemp<CVirtualClass, CReturnClass>() {
  // do nothing - especially do NOT delete(realColumnList)
  realV = NULL;
}

template <typename CVirtualClass, typename CReturnClass> 
void CTemp<CVirtualClass, CReturnClass>::realize() {
  
  realV->realize();

}

#endif // __COLUMNLIST_H_INCLUDED
