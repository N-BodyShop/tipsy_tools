//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CUnion
//
// simple child class of CSubset to take two CSubsets and return their union,
// *WITH NO DUPLICATION CHECKING* (as yet)

#ifndef __UNION_H_INCLUDED

#define __UNION_H_INCLUDED

#include "subset.hpp"
#include <list>
#include <algorithm>

class CUnion : public CSubset {

 public:
  CUnion(CSimSnap *pSim);

  

  ~CUnion();

  virtual int getNumParticles();
  virtual CParticle* getParticle(int id);
  virtual void releaseParticle(CParticle *);
  
  virtual void pushParticle(int n);

  void add(CSimSnap *pSub);

  void cache();

protected:

  list<CSimSnap*> children;

};


#endif // UNION_H_INCLUDED
