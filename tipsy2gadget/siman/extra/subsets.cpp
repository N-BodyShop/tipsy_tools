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

#define VIRTUALONLY(s) std::cerr << "CSubsets: Base class erroneously called for virtual method " << s

int CSubsets::getGroupParticleList(int groupID, int **particlearray) {
  VIRTUALONLY("getGroupParticleList");
  return 0;
}

void CSubsets::getGroupCentre(int groupID, float* cx, float* cy, float *cz) {
  VIRTUALONLY("getGroupCentre");
}


auto_ptr<CSimSnap> CSubsets::getGroup(const CSimSnap & simulation, int haloid) {

  int *pTempRefs;

  int nParticles = getGroupParticleList(haloid,&pTempRefs);

  // cheeky const_cast follows
  //
  // part of the larger problem of missing consts and 
  // pointers vs references throughout code

  CSubset *ss = new CSubset(&(const_cast<CSimSnap&>(simulation)),nParticles,pTempRefs);

  delete[] pTempRefs;

  return auto_ptr<CSimSnap>(ss);


  /* OLD CODE FOLLOWS
     
  halos->getGroupCentre(nHaloID,&cx,&cy,&cz);

  cx*=simulation->getDistanceUnits();
  cy*=simulation->getDistanceUnits();
  cz*=simulation->getDistanceUnits();
  
  float BoxSize = simulation->getBoxSize();

  pParticles= new CParticle[nParticles]();

  for(int n=0;n<nParticles;n++) {
    CParticle *temp = simulation->getParticle(particle_list[n]);

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

int CSubsets::getGroupParticleLen(int groupID)
{
  VIRTUALONLY("getGroupParticleLen");
  return 0;
}

