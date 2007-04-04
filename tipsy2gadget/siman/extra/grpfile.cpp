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


GrpFile::GrpFile(SimSnap *snap_i, string fname) : Subsets(snap_i)
{  
  ifstream fin(fname.c_str());
  fin >> size;
  arr = new int[size];
  for(unsigned int i=0; i<size; i++) {
    fin >> arr[i];
    
  }

  max = 0;
  for(unsigned int i=0; i<size; i++) {
    if(arr[i]>max)
      max=arr[i];
  }

}

GrpFile::~GrpFile() {
  delete arr;
}

void GrpFile::getGroupParticleList(int group_id, vector<unsigned int> & particlearray) {
  if(group_id>max)
    throw(std::out_of_range(boost::lexical_cast<string>(group_id)));

  for(unsigned int i=0; i<size; i++) {
    if(arr[i]==group_id)
      particlearray.push_back(i);
  }
  if(getVerbose()>3)
    cerr << "GrpFile: " << particlearray.size() << " particles in gp " << group_id << " (of " << size << " tot)" <<  endl;
}

void GrpFile::getGroupCentre(int group_id, float *cx, float *cy, float *cz) {
  
}

unsigned int GrpFile::getNumGroups() const {

  return max;
  
}

int GrpFile::getGroupParticleLen(int group_id)
{
  int l=0;
  for(unsigned int i=0; i<size; i++) {
    if(arr[i]==group_id)
      l++;
  }
  return l;
}
