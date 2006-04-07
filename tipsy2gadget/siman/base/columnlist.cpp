//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

// 
// CColumnList
//

CColumnList::CColumnList(CSimSnap *parentSim, float x1, float x2, int nxi, bool autoAssign) {

  nx = nxi;

  pColumnContents = (CSimSnap **) malloc(sizeof(void*) * nx);
  
  dx = (x2-x1)/(float)nx;
    
  for(int n=0;n<nx;n++) {

    pColumnContents[n] = new CSubset(parentSim);
    
  }

  if(autoAssign) {
    unsigned int numParticles = parentSim->getNumParticles();
    for(int n=0;n<numParticles;n++) {
      int x_ref;
      CParticle *particle = parentSim->getParticle(n);
      
      x_ref = (int) ((particle->x - x1)/dx);
      
      
      CSimSnap *pSnap;
      
      if(x_ref>=0 && x_ref<nx) {
	
	pSnap = ((*this)[x_ref]);
	
	pSnap->pushParticle(n);
      }
      
      parentSim->releaseParticle(particle);
    }
  }
}

CColumnList::CColumnList() {
  // only available to child classes, which is CTempColumnList...
  // see grid.hpp for more information!
  pColumnContents = NULL;
}

float CColumnList::power(float wavenumber) {
  complex <float> sum;
  for(int n=0;n<nx;n++) {
    float cellMass = (*this)[n]->getTotalMass();
    complex <float> exponent(0,wavenumber*(float)n*dx);
    sum += cellMass * exp(exponent);
  }

  return pow(abs(sum),2);
}


void CColumnList::realize() {

  // dereference initially passed simulation

  for(int n=0;n<nx;n++) {
    CBaseSimSnap *pReal;
    pReal = new CBaseSimSnap(pColumnContents[n]);
    delete pColumnContents[n];
    pColumnContents[n] = pReal;
  }
}

CColumnList::~CColumnList() {
 
  if(pColumnContents!=NULL) {
    for(int n=0;n<nx;n++) {
      delete pColumnContents[n];
    }
    
    free(pColumnContents);
  }
}

CSimSnap * CColumnList::operator[](int n) {

  return pColumnContents[n];
}
