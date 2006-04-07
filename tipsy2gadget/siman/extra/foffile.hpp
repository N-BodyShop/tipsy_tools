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

class CFoFFile: public CSubsets {

public:

  CFoFFile(char *path, char *snapshot_name, int snap_id);
  virtual ~CFoFFile();

  virtual int getGroupParticleList(int groupID,int **particlearray);

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
