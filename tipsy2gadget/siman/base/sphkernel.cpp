// sphkernel.cpp - part of SimAn Simulation Analysis Library
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

float SPHKernel::operator()(float r, float h) const {

  // Default: spline kernel

  // Create child classes and override operator() to provide new kernels

  // See Hernquist & Katz 1989

  float fac1 = PI * h*h*h;
  float z = r/h;
  float fac2;
  if(z<1) {
    fac2 = 1-(3./2.)*z*z + (3./4.) * z*z*z;
  } else if(z<2) {
    fac2 = (1./4.)*(2.-z)*(2.-z)*(2.-z);
  } else {
    fac2 = 0.;
  }

  return fac2/fac1;
}

float SPHSquareKernel::operator()(float r, float h) const {
  if(r/h<1) return 1/((4.*PI/3.)*h*h*h);
  return 0;
}

SPHProjectionKernel::SPHProjectionKernel(const SPHKernel &from) {
  
  // These could probably do with generalisation:
  num_vals = 20;
  max_z = 2.2;
  
  pVals = new float[num_vals];

  for(int n=0; n<num_vals; n++) {
    float z = max_z*(float)n/(float)num_vals;
    pVals[n]=0;
    float delta_d = 0.01;
    for(float d=0;d<max_z;d+=delta_d) {
      float distance = sqrt(z*z+d*d);
      pVals[n]+=from(distance,1.) * delta_d;
    }
  }
  
}

SPHProjectionKernel::~SPHProjectionKernel() {
  delete[] pVals;
}

float SPHProjectionKernel::operator()(float r, float h) const {
  /*
  if(r<=2*h) 
    return 1./(4*PI*h*h);
  else
    return 0.;
  */

  float z = r/h;
  int z_i_low = (int)(num_vals*z/max_z);
  float z_i_resid = (z - max_z*z_i_low/(float)num_vals);
  if(z_i_low<num_vals) {
    float val_low = pVals[z_i_low];
    float val_high=0;
    if(z_i_low+1<num_vals) {
      val_high = pVals[z_i_low+1];
    }
    return 2*(val_low*(1.-z_i_resid) + val_high*z_i_resid)/(h*h);
  } else {
    return 0.;
  }
}

} // namespace siman
