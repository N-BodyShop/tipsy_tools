//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
//
// rational.hpp declares CRational, a very simple implementation
// for rational numbers, used for CUnits where floats are no good

#ifndef __RATIONAL_H_INCLUDED

#define __RATIONAL_H_INCLUDED

class CRational {

public:
  CRational(string s);
  CRational(int p, int q);
  CRational(int p);
  CRational();

  float flt() const;
  double dbl() const;
  int int_part() const;
  float flt_part() const;
  double dbl_part() const;

  void fromString(string s);
  static int gcd(int a, int b);
  void rationalize();

  CRational operator+(const CRational &addTo) const;
  CRational operator-(const CRational &sub) const;
  CRational operator/(const CRational &divBy) const;
  CRational operator*(const CRational &mulBy) const;
  void operator*=(const CRational &mul);
  void operator+=(const CRational &add);
  void operator-=(const CRational &sub);

  CRational operator-() const;
  
  void operator=(const CRational &copy);
  bool operator==(const CRational &compare) const;
  bool operator!=(const CRational &compare) const;
  int p;
  int q;

};

std::istream & operator>>(std::istream &is, CRational &in);
std::ostream & operator<<(std::ostream &os, const CRational &out);
CRational operator*(int a, CRational b);

#endif
