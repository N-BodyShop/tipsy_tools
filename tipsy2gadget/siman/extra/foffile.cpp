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


FoFFile::FoFFile(SimSnap *snap_i, char *path, char *snapshot_name, int snapshot) : Subsets(snap_i)
{  
  sprintf(fname_cat,  "%s/groups_catalogue/fof_special_catalogue_%03d", path, snapshot);
  sprintf(fname_list,  "%s/groups_indexlist/fof_special_indexlist_%03d", path, snapshot);

  load();
}

FoFFile::~FoFFile() {

}

void FoFFile::load() {
 
  FILE *fh;

  // open file

  if(!(fh=fopen(fname_cat,"r"))) {
    fprintf(stderr,"FoFFile: failed to open cat file '%s'\n",fname_cat);
    return;
  }

  // read number of groups

  fread(&n_groups,sizeof(int),1,fh);

  std::cerr << "FoFFile: allocating memory for " << n_groups << " groups...";

  // allocate memory
  group_offset =(int*) malloc(n_groups*sizeof(int));
  group_len = (int*) malloc(n_groups*sizeof(int));
  group_mass = (float*) malloc(n_groups*sizeof(float));
  group_com = (float*) malloc(n_groups*sizeof(float)*3);

  std::cerr << "done!\n";

  // read information on groups
  fread(group_len,sizeof(int),n_groups,fh); 
  fread(group_offset,sizeof(int),n_groups,fh); 
  fread(group_mass,sizeof(int),n_groups,fh);
  fread(group_com,sizeof(int),n_groups,fh); 

  fclose(fh); // that's all we need from that file...

}

void FoFFile::getGroupParticleList(int group_id, vector<unsigned int> & particles) {

  FILE *fh;
  

  if(!(fh=fopen(fname_list,"r"))) {
    fprintf(stderr,"read_fof_file: Error opening file\n");
    exit(0);
  }


  
 

  fseek(fh,sizeof(int) + group_offset[group_id]*4,SEEK_CUR);

  int *particlearray = new int[group_len[group_id]];

  fread(particlearray,4,group_len[group_id],fh);
 
  // particle array from Volker Springel's fof_special is
  // 1-based, not 0-based!!

  for(int n=0;n<group_len[group_id];n++) {
    particles.push_back(particlearray[n]-1);
   }
  
  
  delete[] particlearray;

  

  fclose(fh);

  
  

}

void FoFFile::getGroupCentre(int group_id, float *cx, float *cy, float *cz) {
  
  *cx = group_com[3*group_id];
  *cy = group_com[3*group_id+1];
  *cz = group_com[3*group_id+2];

}


int FoFFile::getGroupParticleLen(int group_id)
{
  return group_len[group_id];
}
