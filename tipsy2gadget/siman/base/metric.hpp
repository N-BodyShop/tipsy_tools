// metric.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __METRIC_H_INCLUDED

#define __METRIC_H_INCLUDED

namespace siman {

/// Metric class
///
/// A metric class encodes just that, a way of measuring distances
/// between two points. This is experimental at the moment. 
/// The base class, Metric, provides a 3D Euclidean metric, i.e.
///
/// \f$d^2 = x^2+y^2+z^2\f$
///
class Metric {
public:
  Metric();
  virtual ~Metric() { } ;
  /// returns the distance between p1 and p2
  virtual float operator()(const Particle &p1, const Particle &p2) const;
};

/// Metric2D class
///
/// Returns a metric, discarding z information, so that
///
/// \f$d^2 = x^2+y^2\f$
///

class Metric2D : public Metric {
public:
  Metric2D();

  virtual float operator()(const Particle &p1, const Particle &p2) const;
};

} // namespace siman

#endif
