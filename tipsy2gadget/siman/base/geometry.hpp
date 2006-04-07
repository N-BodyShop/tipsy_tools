//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __GEOMETRY_H_INCLUDED

#define __GEOMETRY_H_INCLUDED


#include "simsnap.hpp"
#include "particle.hpp"

class CGeometry
{
public: 
  
  CGeometry(CSimSnap *parentSimIn);
  
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

  CSimSnap *pData;
  
};

#endif



