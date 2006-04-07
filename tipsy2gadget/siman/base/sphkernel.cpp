//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// SPHKERNEL.CPP

#include "siman.hpp"


float CSPHKernel::operator()(float r, float h) const {

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

float CSPHSquareKernel::operator()(float r, float h) const {
  if(r/h<1) return 1/((4.*PI/3.)*h*h*h);
  return 0;
}

CSPHProjectionKernel::CSPHProjectionKernel(const CSPHKernel &from) {
  
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

CSPHProjectionKernel::~CSPHProjectionKernel() {
  delete[] pVals;
}

float CSPHProjectionKernel::operator()(float r, float h) const {
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
