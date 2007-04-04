// geometry.cpp - part of SimAn Simulation Analysis Library
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

Geometry::Geometry(SimSnap *pSim)
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

void Geometry::applyToParticle(int n) {
  Particle *pParticle = pData->getParticle(n);
  

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

void Geometry::apply() {
  cerr << "Geometry: applying transformations...";
  unsigned int np = pData->getNumParticles();
  for(unsigned int n=0;n<np;n++) {
    applyToParticle(n);
  }
  cerr << "done!" << endl;
}

void Geometry::reCentre(float cxi, float cyi, float czi) {
  cx+= cxi;
  cy+= cyi;
  cz+= czi;
}

void Geometry::reCentreVel(float cvxi, float cvyi, float cvzi) {
  cvx+= cvxi;
  cvy+= cvyi;
  cvz+= cvzi;
}


void Geometry::setRotateX(float rxi) {
  rx = rxi;
  transform = true;
  updateTransformCoeffs();
}

void Geometry::setRotateY(float ryi) {
  ry = ryi;
  
  transform = true;
  updateTransformCoeffs();
}

void Geometry::setRotateZ(float rzi) {
  rz = rzi;
  
  transform = true;
  updateTransformCoeffs();
}

void Geometry::setRotate(float rxi, float ryi, float rzi) {
  rx = rxi;
  ry = ryi;
  rz = rzi;
  
  transform = true;
  updateTransformCoeffs();
}

void Geometry::setMatrix(float *mat) {
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

void Geometry::updateTransformCoeffs() {

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

}
