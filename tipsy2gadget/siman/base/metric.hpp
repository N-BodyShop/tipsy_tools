//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//

#ifndef __METRIC_H_INCLUDED

#define __METRIC_H_INCLUDED

class CMetric {
public:
  CMetric();
  virtual float operator()(const CParticle &p1, const CParticle &p2) const;
};

class CMetric2D : public CMetric {
public:
  CMetric2D();
  virtual float operator()(const CParticle &p1, const CParticle &p2) const;
};

#endif
