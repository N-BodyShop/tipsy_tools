//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//

#include "siman.hpp"

CMetric::CMetric() {
}

float CMetric::operator()(const CParticle &p1, const CParticle &p2) const {
  return sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z));
}


CMetric2D::CMetric2D() {
}

float CMetric2D::operator()(const CParticle &p1, const CParticle &p2)  const {
  return sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}
