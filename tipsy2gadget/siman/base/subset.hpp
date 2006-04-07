//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CSubset
//
// Selects a subset of a simulation, and is derived itself from CSimSnap so you 
// can use this subset as though it were its own file.

#ifndef __SUBSET_H_INCLUDED
#define __SUBSET_H_INCLUDED


#include <vector>

using namespace std;

class CSubset : public CSimSnap {

 public:

  CSubset(CSimSnap *simulation, CFilter &filter);
  CSubset(CSimSnap *simulation, int numPart, int *pParticleIDs);
  CSubset(CSimSnap *simulation);

  string className();
  CSimanObject * getMember(const string & member);

  ~CSubset();


  // Standard functions which probably do not require recoding for
  // child classes, if they store the data in the provided "protected"
  // variables

  virtual int getNumParticles();
  virtual CParticle* getParticle(int id);
  virtual void releaseParticle(CParticle *);
  virtual bool isLoaded();
  virtual void pushParticle(int n);
  virtual int deReference(int i, int n=1);

  CSimSnap* getParent();

  void cache();

protected:

  CSubset();

  int nParticles;
  vector<int> particleRefs;

private:
    
  int nHaloID;

};

#endif
