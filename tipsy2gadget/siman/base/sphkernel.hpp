//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// SPH kernel classes

// used to represent different SPH kernels

#ifndef __SPHKERNEL_H_INCLUDED

#define __SPHKERNEL_H_INCLUDED
	
class CSPHKernel {
public:
  virtual float operator()(float r, float h) const;

};

class CSPHSquareKernel : public CSPHKernel {
public:
  virtual float operator()(float r, float h) const;
};

class CSPHProjectionKernel : public CSPHKernel {
public:
  CSPHProjectionKernel(const CSPHKernel &from);
  ~CSPHProjectionKernel();
  virtual float operator()(float r, float h) const;
  
private:
  float max_z;
  int num_vals;
  float *pVals;
};

#endif // SPHKERNEL_H_INCLUDED
