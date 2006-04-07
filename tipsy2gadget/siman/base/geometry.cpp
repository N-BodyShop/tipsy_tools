//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"


CGeometry::CGeometry(CSimSnap *pSim)
{ 
  rx = 0;
  ry = 0;
  rz = 0;
  cvx=cvy=cvz=cx=cy=cz=0;

  updateTransformCoeffs();
  transform = false;
  wrapping = false;
  pData = pSim;
}

void CGeometry::applyToParticle(int n) {
  CParticle *pParticle = pData->getParticle(n);
  

  pParticle->x -=cx;
  pParticle->y -=cy;
  pParticle->z -=cz;

  pParticle->vx-=cvx;
  pParticle->vy-=cvy;
  pParticle->vz-=cvz;

  if(wrapping) {
    float boxsize = pData->getBoxSize();

    if(pParticle->x > boxsize/2)
      pParticle->x -= boxsize;
    if(pParticle->y > boxsize/2)
      pParticle->y -= boxsize;
    if(pParticle->z > boxsize/2)
      pParticle->z -= boxsize;

    if(pParticle->x < -boxsize/2)
      pParticle->x += boxsize;
    if(pParticle->y < -boxsize/2)
      pParticle->y += boxsize;
    if(pParticle->z < -boxsize/2)
      pParticle->z += boxsize;
    
  }

  // rotation - Not done with matrix multiplication, as
  //            processor pipelining probably does better
  //            without loops!

  if(transform) {
    const float tx = pParticle->x, ty = pParticle->y, tz = pParticle->z;

    pParticle->x = mxx * tx + mxy * ty + mxz * tz;
    pParticle->y = myx * tx + myy * ty + myz * tz;
    pParticle->z = mzx * tx + mzy * ty + mzz * tz;

  }

}

void CGeometry::apply() {
  cerr << "CGeometry: applying transformations...";
  unsigned int np = pData->getNumParticles();
  for(unsigned int n=0;n<np;n++) {
    applyToParticle(n);
  }
  cerr << "done!" << endl;
}

void CGeometry::reCentre(float cxi, float cyi, float czi) {
  cx+= cxi;
  cy+= cyi;
  cz+= czi;
}

void CGeometry::reCentreVel(float cvxi, float cvyi, float cvzi) {
  cvx+= cvxi;
  cvy+= cvyi;
  cvz+= cvzi;
}


void CGeometry::setRotateX(float rxi) {
  rx = rxi;
  transform = true;
  updateTransformCoeffs();
}

void CGeometry::setRotateY(float ryi) {
  ry = ryi;
  
  transform = true;
  updateTransformCoeffs();
}

void CGeometry::setRotateZ(float rzi) {
  rz = rzi;
  
  transform = true;
  updateTransformCoeffs();
}

void CGeometry::setRotate(float rxi, float ryi, float rzi) {
  rx = rxi;
  ry = ryi;
  rz = rzi;
  
  transform = true;
  updateTransformCoeffs();
}

void CGeometry::setMatrix(float *mat) {
  transform=true;
  mxx = mat[0];
  mxy = mat[1];
  mxz = mat[2];
  myx = mat[3];
  myy = mat[4];
  myz = mat[5];
  mzx = mat[6];
  mzy = mat[7];
  mzz = mat[8];
}

void CGeometry::updateTransformCoeffs() {

  mxx = cos(rz)*cos(ry); // cos(rz)*cos(ry) - sin(rx)*sin(ry)*sin(rz);
  mxy = sin(rz)*cos(rx)-sin(rx)*cos(rz)*sin(ry); // cos(rx)*sin(rz);
  mxz = sin(rz)*sin(rx)+cos(rz)*sin(ry)*cos(rx); // cos(rz)*sin(ry) + sin(rx)*sin(rz)*cos(ry);
  
  myx = -sin(rz)*cos(ry); //-sin(rz)*cos(ry) - sin(ry)*sin(rx)*cos(rz);
  myy = cos(rz)*cos(rx)+sin(rx)*sin(rz)*sin(ry); // cos(rx)*cos(rz);
  myz = cos(rz)*sin(rx)-sin(rz)*sin(ry)*cos(rx); // -sin(ry)*sin(rz) + cos(ry)*sin(rx)*cos(rz);

  mzx = -sin(ry); // -sin(ry)*cos(rz);
  mzy = -sin(rx)*cos(ry); // -sin(rx);
  mzz = cos(rx)*cos(ry);

}
