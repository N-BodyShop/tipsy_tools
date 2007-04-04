// visualiser.hpp - part of SimAn Simulation Analysis Library
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


#ifndef __ANNOTATE_H_INCLUDED

#define __ANNOTATE_H_INCLUDED

#include "../base.hpp"

namespace siman {
  class Annotate {
  protected:
    Annotate();
  public:
    virtual void plot(float lengthscale);
    virtual ~Annotate();
  };
  
  class AnnotateVector : public Annotate {
  public:
    AnnotateVector(const SimanVec &len);
    AnnotateVector(const SimanVec &start, const SimanVec &len);
    virtual void plot(float ls);
    virtual ~AnnotateVector();
  protected:
    SimanVec start;
    SimanVec len;
  };

  class AnnotatePlane : public Annotate {
  public:
    AnnotatePlane(const SimanVec &len);
    AnnotatePlane(const SimanVec &start, const SimanVec &len);
    virtual void plot(float ls);
    virtual ~AnnotatePlane();
  protected:
    SimanVec start;
    SimanVec len;
    SimanVec v1;
    SimanVec v2;
  };

} // namespace siman
#endif // ANNOTATE_H_INCLUDED
