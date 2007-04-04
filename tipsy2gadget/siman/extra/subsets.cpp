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

#define VIRTUALONLY(s) std::cerr << "Subsets: Base class erroneously called for virtual method " << s

Subsets::Subsets(SimSnap *snap) : snap(snap) { }

void Subsets::getGroupParticleList(int groupID, vector<unsigned int> & particlearray) {
  VIRTUALONLY("getGroupParticleList");
}

void Subsets::getGroupCentre(int groupID, float* cx, float* cy, float *cz) {
  VIRTUALONLY("getGroupCentre");
}

unsigned int Subsets::getNumGroups() const {
  VIRTUALONLY("getNumGroups");
  return 0;
}

auto_ptr<SimSnap> Subsets::getGroup(int haloid) {

  vector<unsigned int> tempRefs;


  getGroupParticleList(haloid,tempRefs);

  // cheeky const_cast follows
  //
  // part of the larger problem of missing consts and 
  // pointers vs references throughout code

  if(getVerbose()>3) 
    cerr << "Subsets: construct group with " << tempRefs.size() << " ptcls" <<  endl;

  Subset *ss = new Subset(snap,tempRefs);
  

  return auto_ptr<SimSnap>(ss);


  /* OLD CODE FOLLOWS
     
  halos->getGroupCentre(nHaloID,&cx,&cy,&cz);

  cx*=simulation->getDistanceUnits();
  cy*=simulation->getDistanceUnits();
  cz*=simulation->getDistanceUnits();
  
  float BoxSize = simulation->getBoxSize();

  pParticles= new Particle[nParticles]();

  for(int n=0;n<nParticles;n++) {
    Particle *temp = simulation->getParticle(particle_list[n]);

    // copy data into our own array
    pParticles[n] = *temp;

    // give particle back to simulation
    simulation->releaseParticle(temp);

    // recentre...

    pParticles[n].x-=cx;
    pParticles[n].y-=cy;
    pParticles[n].z-=cz;

    // wrap around...

    if(pParticles[n].x<-BoxSize*0.5) pParticles[n].x+=BoxSize;
    if(pParticles[n].y<-BoxSize*0.5) pParticles[n].y+=BoxSize;
    if(pParticles[n].z<-BoxSize*0.5) pParticles[n].z+=BoxSize;
    
    if(pParticles[n].x>BoxSize*0.5)  pParticles[n].x-=BoxSize;
    if(pParticles[n].y>BoxSize*0.5)  pParticles[n].y-=BoxSize;
    if(pParticles[n].z>BoxSize*0.5)  pParticles[n].z-=BoxSize;


  }

  free(particle_list);

  */

  
}

int Subsets::getGroupParticleLen(int groupID)
{
  VIRTUALONLY("getGroupParticleLen");
  return 0;
}

