// columnlist.cpp - part of SimAn Simulation Analysis Library
//
//
// Copyright (c) Andrew Pontzen 2005, 2006
//
// SimAn is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// SimAn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public Licence for more details.
//
// You should have received a copy of the GNU General Public Licence
// along with SimAn; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA









#include "siman.hpp"

namespace siman {

ColumnList::ColumnList(SimSnap &parentSim, float x1, float x2, int nxi, bool autoAssign) {

  nx = nxi;

  pColumnContents = (SimSnap **) malloc(sizeof(void*) * nx);
  
  dx = (x2-x1)/(float)nx;
    
  for(int n=0;n<nx;n++) {

    pColumnContents[n] = new Subset(parentSim);
    
  }

  if(autoAssign) {
    unsigned int numParticles = parentSim.getNumParticles();
    for(unsigned int n=0;n<numParticles;n++) {
      int x_ref;
      Particle *particle = parentSim.getParticle(n);
      
      x_ref = (int) ((particle->x - x1)/dx);
      
      
      SimSnap *pSnap;
      
      if(x_ref>=0 && x_ref<nx) {
	
	pSnap = &((*this)[x_ref]);
	
	(static_cast<Subset*>(pSnap))->pushParticle(n);
      }
      

    }
  }
}

float ColumnList::power(float wavenumber) {
  complex <float> sum;
  for(int n=0;n<nx;n++) {
    float cellMass = (*this)[n].getTotalMass();
    complex <float> exponent(0,wavenumber*(float)n*dx);
    sum += cellMass * exp(exponent);
  }

  return std::pow(abs(sum),2);
}


void ColumnList::realize() {

  // dereference initially passed simulation

  for(int n=0;n<nx;n++) {
    BaseSimSnap *pReal;
    pReal = new BaseSimSnap(pColumnContents[n]);
    delete pColumnContents[n];
    pColumnContents[n] = pReal;
  }
}

ColumnList::~ColumnList() {
 
  if(pColumnContents!=NULL) {
    for(int n=0;n<nx;n++) {
      delete pColumnContents[n];
    }
    
    free(pColumnContents);
  }
}

SimSnap & ColumnList::operator[](int n) {

  return *pColumnContents[n];
}

} // namespace siman
