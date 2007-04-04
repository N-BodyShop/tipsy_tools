// geometry.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __GEOMETRY_H_INCLUDED

#define __GEOMETRY_H_INCLUDED


#include "simsnap.hpp"
#include "particle.hpp"

namespace siman {

class Geometry
{
public: 
  
  Geometry(SimSnap *parentSimIn);
  
  void reCentre(float cxi, float cyi, float czi);
  void reCentreVel(float cvxi, float cvyi, float cvzi);

  bool wrapping;
  bool transform;

  void setRotateX(float rxi);
  void setRotateY(float ryi);
  void setRotateZ(float rzi);

  void setRotate(float rxi, float ryi, float rzi);
  void setMatrix(float *mat);
  void applyToParticle(int p);
  void apply();

private:

  void updateTransformCoeffs();

  float mxx, mxy, mxz;
  float myx, myy, myz;
  float mzx, mzy, mzz;

  float cx, cy, cz;
  float rx, ry, rz;
  float cvx, cvy, cvz;

  SimSnap *pData;
  
};

}

#endif



