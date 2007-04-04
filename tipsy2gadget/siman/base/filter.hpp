// filter.hpp - part of SimAn Simulation Analysis Library
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







#include "particle.hpp"
#include <boost/shared_ptr.hpp>

#ifndef __FILTER_H_INCLUDED

#define __FILTER_H_INCLUDED

namespace siman {

  class NotFilter;
  class AndFilter;
  class OrFilter;

  class Filter {

  public:
    virtual ~Filter() { } 
    virtual bool includes(const Particle &particle) const;
    virtual Filter* copy() const; ///< Return a new identical copy of this filter
    NotFilter operator!() const; ///< Return a filter which returns false when (this) would return true, and vice-versa
    AndFilter operator&(const Filter &f2) const; ///< Combine filters, returning true only when both filters return true
    OrFilter operator|(const Filter &f2) const; ///< Combine filters, returning true if either or both filters return true
  };




  class Column: public Filter {

  public:
    Column(float x1i, float y1i, float x2i, float y2i);
    bool includes(const Particle &particle) const;
    Filter* copy() const;

  private:

    float x1,y1,x2,y2;

  };


  class Sphere: public Filter {

  public:
    Sphere(float xci, float yci, float zci, float ri); ///< Construct sphere centred on xci,yci,zci, radius ri
    Sphere(float ri); ///< Construct sphere centered on 0,0,0, radius ri
    bool includes(const Particle &particle) const;
    Filter* copy() const;

  private:

    float xc, yc, zc, r;

  };

  class VelocitySphereFilter: public Filter {
  public:
    VelocitySphereFilter(float xvi, float yvi, float zvi, float ri); ///< Construct velocity-space sphere
    VelocitySphereFilter(float ri); ///< Construct velocity-space sphere centre on 0,0,0
    bool includes(const Particle &particle) const;
    Filter *copy() const;
  private:
    float xc, yc, zc, r;
  };

  class ParticleTypeFilter: public Filter {
  public:
    ParticleTypeFilter(int type);
    bool includes(const Particle &particle) const;
    Filter* copy() const;

  private:

    int type;
  };

  class DensityCutFilter: public Filter {
  public:
    DensityCutFilter(float cutAti);
    bool includes(const Particle &particle) const;
    Filter* copy() const;

  private:

    float cutAt;
  };

  class RandomFilter: public Filter {
  public:
    RandomFilter(float probability);
    bool includes(const Particle &particle) const;
    Filter* copy() const;

  private:
    float prob;
  };

  class ModuloFilter : public Filter {
  public:
    ModuloFilter(int number, int offset);
    bool includes(const Particle &particle) const;
    Filter* copy() const;
  private:
    int mod;
    int mutable cur;
  };

  class NotFilter : public Filter {
  public:
    NotFilter(const Filter &fi);
    NotFilter(const boost::shared_ptr<const Filter> &f);
    bool includes(const Particle &particle) const;
    Filter* copy() const;
  private:
    boost::shared_ptr<const Filter> f;
  };

  class AndFilter : public Filter {
  public:
    AndFilter(const Filter &f1i, const Filter &f2i);
    AndFilter(const boost::shared_ptr<const Filter> &f1,const boost::shared_ptr<const Filter> &f2 );
    bool includes(const Particle &particle) const;
    Filter* copy() const;
  private:
    boost::shared_ptr<const Filter> f1;
    boost::shared_ptr<const Filter> f2;
  };

  class OrFilter : public Filter {
  public:
    OrFilter(const Filter &f1i, const Filter &f2i);
    OrFilter(const boost::shared_ptr<const Filter> &f1,const boost::shared_ptr<const Filter> &f2 );
    bool includes(const Particle &particle) const;
    Filter* copy() const;
  private:
    boost::shared_ptr<const Filter> f1;
    boost::shared_ptr<const Filter> f2;
  };

}

#endif // __FILTER_H_INCLUDED
