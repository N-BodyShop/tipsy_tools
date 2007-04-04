//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "../base.hpp"
#include "../extra.hpp"

using namespace siman;

class GrpFile: public Subsets {

public:

  GrpFile(SimSnap *snap, std::string fname);
  virtual ~GrpFile();

  virtual void getGroupParticleList(int groupID,std::vector<unsigned int> & particles);

  virtual int getGroupParticleLen(int groupID);

  virtual void getGroupCentre(int groupID, float *cx, float *cy, float *cz);
  virtual unsigned int getNumGroups() const;

private:

  

  unsigned int size;
  unsigned int max;
  int *arr;

  
  
};
