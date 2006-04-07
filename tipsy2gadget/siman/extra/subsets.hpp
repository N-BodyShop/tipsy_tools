//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __SUBSETS_H_INCLUDED

#define __SUBSETS_H_INCLUDED

class CSubsets {

public:
  virtual int getGroupParticleList(int groupID, int **particlearray);
  virtual void getGroupCentre(int groupID, float *cx, float *cy, float *cz);
  virtual int getGroupParticleLen(int groupID);
  auto_ptr<CSimSnap> getGroup(const CSimSnap & snap, int groupID);
  
  
};

#endif
