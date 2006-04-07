//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __COORDINATE_H_INCLUDED

#define __COORDINATE_H_INCLUDED

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
#endif
