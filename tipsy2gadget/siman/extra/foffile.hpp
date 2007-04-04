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

class FoFFile: public Subsets {

public:

  FoFFile(SimSnap *snap, char *path, char *snapshot_name, int snap_id);
  virtual ~FoFFile();

  virtual void getGroupParticleList(int groupID,std::vector<unsigned int> &particles);

  virtual int getGroupParticleLen(int groupID);

  virtual void getGroupCentre(int groupID, float *cx, float *cy, float *cz);

private:

  void load();

  char fname_list[1024];
  char fname_cat[1024];
  
  bool loaded;
  int n_groups;
  int *group_offset;
  int *group_len;
  float *group_mass;
  float *group_com;
  
};
