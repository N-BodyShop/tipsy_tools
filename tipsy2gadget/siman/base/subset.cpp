//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

// CSubset Implementation

using namespace std;

CSubset::~CSubset() {
 
}

CSubset::CSubset() {

  pAutoParent = NULL;
  nHaloID = -1;
} 

CSubset::CSubset(CSimSnap *simulation) {

  // start off with NO particles!
  pAutoParent = simulation;
}

CSubset::CSubset(CSimSnap *simulation,CFilter &filter) {
  pAutoParent = simulation;
  CParticle *consider;

  for(int n=0;n<pAutoParent->getNumParticles();n++) {
    CParticle *pParticle = simulation->getParticle(n);

    if(filter.includes(*pParticle))
      particleRefs.push_back(n);
    
    simulation->releaseParticle(pParticle);
  }
}

CSubset::CSubset(CSimSnap *simulation, int numParticles, int *pPartList) {
  pAutoParent = simulation;
  particleRefs.insert(particleRefs.end(),&pPartList[0],&pPartList[numParticles-1]);
}


string CSubset::className() {
  return "CSubset";
}

CSimSnap *CSubset::getParent() {
  return pAutoParent;
}

CSimanObject * CSubset::getMember(const string & member) {
  if(member=="parent")
    return pAutoParent;

  return CSimSnap::getMember(member);
}

int CSubset::getNumParticles() {
  return particleRefs.size();
}

void CSubset::cache() {
  
}

void CSubset::pushParticle(int n) {
  particleRefs.push_back(n);
}

bool CSubset::isLoaded() {
  return true;
}

CParticle* CSubset::getParticle(int id) {

  // could do checking of id, but a conditional may disturb pipelining
  // so have not for now...

  CParticle* particle = pAutoParent->getParticle(particleRefs[id]);
  
  return particle;

}

void CSubset::releaseParticle(CParticle *particle) {
  // probably do nothing; should streamline this...
  pAutoParent->releaseParticle(particle);

}

int CSubset::deReference(int i, int n) {

  if(n==1) 
    return particleRefs[i];
  else
    return pAutoParent->deReference(particleRefs[i], n-1);
}
