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

using namespace siman;

class Subsets {

protected:
  Subsets(SimSnap *snap) ;
public:
  virtual unsigned int getNumGroups() const;
  virtual void getGroupParticleList(int groupID, std::vector<unsigned int> &particles);
  virtual void getGroupCentre(int groupID, float *cx, float *cy, float *cz);
  virtual int getGroupParticleLen(int groupID);
  std::auto_ptr<SimSnap> getGroup( int groupID);
  virtual ~Subsets() { } ;

protected:
  SimSnap *snap;
  
};

#endif
