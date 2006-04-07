//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CFilter
//
// Class which provides set of conditions for inclusion of
// a particle
//
// This file includes some simple examples, such as a pipe
// and a sphere.

#include "particle.hpp"

#ifndef __FILTER_H_INCLUDED

#define __FILTER_H_INCLUDED

class CFilter {

public:

  virtual bool includes(CParticle &particle);

};




class CColumn: public CFilter {

public:
  CColumn(float x1i, float y1i, float x2i, float y2i);

  bool includes(CParticle &particle);

private:

  float x1,y1,x2,y2;

};


class CSphere: public CFilter {

public:
  CSphere(float xci, float yci, float zci, float ri);

  bool includes(CParticle &particle);

private:

  float xc, yc, zc, r;

};

class CParticleTypeFilter: public CFilter {
public:
  CParticleTypeFilter(int type);
  bool includes(CParticle &particle);

private:

  int type;
};

class CDensityCutFilter: public CFilter {
public:
  CDensityCutFilter(float cutAti);
  bool includes(CParticle &particle);

private:

  float cutAt;
};

class CRandomFilter: public CFilter {
public:
  CRandomFilter(float probability);
  bool includes(CParticle &particle);

private:
  float prob;
};

class CModuloFilter : public CFilter {
public:
  CModuloFilter(int number, int offset);
  bool includes(CParticle &particle);
private:
  int mod, cur;
};

#endif // __FILTER_H_INCLUDED
