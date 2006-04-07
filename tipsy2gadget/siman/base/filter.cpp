//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

using namespace std;

// CFilter

bool CFilter::includes(CParticle &particle) {
  cerr << "CFilter: base member virtual function called..." << endl;
  return true;
}

// PIPE filter class

CColumn::CColumn(float x1i, float y1i, float x2i, float y2i)
{
  x1 = (x1i<x2i)?x1i:x2i;
  x2 = (x1i<x2i)?x2i:x1i;
  
  y1 = (y1i<y2i)?y1i:y2i;
  y2 = (y1i<y2i)?y2i:y1i;
  
}

bool CColumn::includes(CParticle &particle) {
  // is particle within bounds?
  return (particle.x>x1 && particle.y>y1 && particle.y<y2 && particle.x<x2);
}



// SPHERE filter class

CSphere::CSphere(float xci, float yci, float zci, float ri) : 
  xc(xci), yc(yci), zc(zci), r(ri) {

}

bool CSphere::includes(CParticle &p) {

  float distance = sqrt((p.x-xc)*(p.x-xc) + 
			(p.y-yc)*(p.y-yc) +
			(p.z-zc)*(p.z-zc));
  
  return (distance<r);

}

// PARTICLE TYPE filter class

CParticleTypeFilter::CParticleTypeFilter(int typei) {
  type = typei;
}

bool CParticleTypeFilter::includes(CParticle &p) {
  return ((p.type & type) > 0);
}

// DENSITY CUT class

CDensityCutFilter::CDensityCutFilter(float cutAti) {
  cutAt = cutAti;
}

bool CDensityCutFilter::includes(CParticle &p) {
  return(p.rho > cutAt);
}

// RANDOM filter

CRandomFilter::CRandomFilter(float probi) {
  prob = probi;
}

bool CRandomFilter::includes(CParticle &p) {
  float s = (float)rand()/(float)RAND_MAX;
  if(s<prob) return true;
  return false;
}

// MODULO filter

CModuloFilter::CModuloFilter(int num, int offset) {
  cur=offset;
  mod=num;
}

bool CModuloFilter::includes(CParticle &p) {
  cur++;
  if(cur>mod) cur=0;
  return cur==0;
}
