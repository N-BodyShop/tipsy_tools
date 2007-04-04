// simanvec.hpp - part of SimAn Simulation Analysis Library
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




#ifndef __SIMANVEC_H_INCLUDED

#define __SIMANVEC_H_INCLUDED

namespace siman {

  class SimanVec : public std::vector<double> {
  public:
    SimanVec();
    SimanVec(const std::vector<double> &);
    SimanVec(double x, double y, double z);

    void operator/=(double r);
    void operator*=(double r);

    double abs() const;
    
    double dot(const SimanVec &) const throw(MismatchedArrayLength);
    SimanVec proj(const SimanVec &) const throw(MismatchedArrayLength);
 
  };
}

#endif
