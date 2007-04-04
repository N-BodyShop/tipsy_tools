// coordinate.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __COORDINATE_H_INCLUDED

#define __COORDINATE_H_INCLUDED

namespace siman {

template <typename T> 
class coordinate {
public:
  T x;
  T y;
  T z;

  coordinate(T xi,T yi,T zi);

  void operator=(const coordinate<T> &from);
};

template <typename T>
coordinate<T>::coordinate(T xi,T yi,T zi) {
  x=xi;
  y=yi;
  z=zi;
}

template <typename T>
void coordinate<T>::operator=(const coordinate<T> &from) {
  x=from.x;
  y=from.y;
  z=from.z;
}

}

#endif
