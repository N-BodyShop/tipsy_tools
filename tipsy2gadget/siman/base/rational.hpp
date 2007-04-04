// rational.hpp - part of SimAn Simulation Analysis Library
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


#ifndef __RATIONAL_H_INCLUDED

#define __RATIONAL_H_INCLUDED

namespace siman {

class Rational {

public:
  explicit Rational(std::string s);
  explicit Rational(int p, int q);
  Rational(int p);
  Rational();

  float flt() const;
  double dbl() const;
  int int_part() const;
  float flt_part() const;
  double dbl_part() const;

  void fromString(std::string s);
  static int gcd(int a, int b);
  void rationalize();

  std::string toString();

  Rational operator+(const Rational &addTo) const;
  Rational operator-(const Rational &sub) const;
  Rational operator/(const Rational &divBy) const;
  Rational operator*(const Rational &mulBy) const;
  void operator*=(const Rational &mul);
  void operator+=(const Rational &add);
  void operator-=(const Rational &sub);

  Rational operator-() const;
  
  void operator=(const Rational &copy);
  bool operator==(const Rational &compare) const;
  bool operator!=(const Rational &compare) const;
  int p;
  int q;

};

std::istream & operator>>(std::istream &is, Rational &in);
std::ostream & operator<<(std::ostream &os, const Rational &out);
Rational operator*(int a, Rational b);

}

#endif
