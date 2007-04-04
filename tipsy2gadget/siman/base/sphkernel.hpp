// sphkernel.hpp - part of SimAn Simulation Analysis Library
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







// used to represent different SPH kernels

#ifndef __SPHKERNEL_H_INCLUDED

#define __SPHKERNEL_H_INCLUDED

namespace siman {
	
class SPHKernel {
public:
  virtual float operator()(float r, float h) const;
  virtual ~SPHKernel() { }
};

class SPHSquareKernel : public SPHKernel {
public:
  virtual float operator()(float r, float h) const;
};

class SPHProjectionKernel : public SPHKernel {
public:
  SPHProjectionKernel(const SPHKernel &from);
  virtual ~SPHProjectionKernel();
  virtual float operator()(float r, float h) const;
  
private:
  float max_z;
  int num_vals;
  float *pVals;
};

} // namespace siman

#endif // SPHKERNEL_H_INCLUDED
