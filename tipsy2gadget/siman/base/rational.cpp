// rational.cpp - part of SimAn Simulation Analysis Library
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










#include "siman.hpp"

namespace siman {

Rational::Rational() : p(0), q(1) {

}


Rational::Rational(int n) : p(n), q(1) {

}

Rational::Rational(int n, int d) : p(n), q(d) {
  rationalize();
}

float Rational::flt() const {
  return ((float)p/(float)q);
}


double Rational::dbl() const {
  return ((double)p/(double)q);
}

int Rational::int_part() const {
  return (int) flt();
}

float Rational::flt_part() const {
  return flt() - (float) int_part();
}


double Rational::dbl_part() const {
  return dbl() - (double) int_part();
}


void Rational::rationalize() {
  if(p==0) {
    q=1;
    return;
  }

  if(q==0) {
    // could throw an exception here?
    return;
  }

  int gcdn = gcd(p,q);
  p/=gcdn;
  q/=gcdn;
  if(q<0) {
    p=-p;
    q=-q;
  }
}

Rational Rational::operator+(const Rational &a) const {
  Rational out;
  out.p=p*a.q+q*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}


Rational Rational::operator-(const Rational &a) const {
  Rational out;
  out.p=p*a.q-q*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}

Rational Rational::operator-() const {
  Rational out;
  out.p=-p;
  out.q=q;
  return out;
}

Rational Rational::operator*(const Rational &a) const {
  Rational out;
  out.p=p*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}


Rational Rational::operator/(const Rational &a) const {
  Rational out;
  out.p=p*a.q;
  out.q=q*a.p;
  out.rationalize();
  return out;
}

void Rational::operator*=(const Rational &a) {
  p*=a.p;
  q*=a.q;
  rationalize();
}


void Rational::operator+=(const Rational &a) {
  p=p*a.q + a.p*q;
  q*=a.q;
  rationalize();
}


void Rational::operator-=(const Rational &a) {
  p=p*a.q - a.p*q;
  q*=a.q;
  rationalize();
}

void Rational::operator=(const Rational &c) {
  p=c.p;
  q=c.q;
}

bool Rational::operator==(const Rational &e) const {
  return ((e.p==p) && (e.q==q));
}

bool Rational::operator!=(const Rational &e) const {
  return ((e.p!=p) || (e.q!=q));
}


int Rational::gcd(int a, int b) {
  if(a<0) a=-a;
  if(b<0) b=-b;

  int gcd=0;
  int min = a;
  if(b<a) min=b;
  
  for(int n=1; n<=min; n++) {
    if(a%n==0 && b%n==0)
      gcd=n;
  }

  return gcd;
}


Rational::Rational(string s) {
   
  fromString(s);

}

void Rational::fromString(string s) {
  
  string::size_type pos_div = s.find("/",0);
   
  if(pos_div!=string::npos) {
    q = atoi(s.substr(pos_div+1,s.length()-pos_div-1).c_str());
    p = atoi(s.substr(0,pos_div).c_str());
  } else {
    p = atoi(s.c_str());
    q = 1;
  }

  rationalize();
}

  string Rational::toString() {
    ostringstream s;
    s << (*this);
    return s.str();
  }

std::istream & operator>>(std::istream &is, Rational &in) {
  string s;
  is >> s;
 
  in.fromString(s);

  return is;
}

std::ostream & operator<<(std::ostream &os, const Rational &out) {
  if(out.q!=1) 
    os << out.p << "/" << out.q;
  else
    os << out.p;
  return os;
}

Rational operator*(int a, Rational b) {
  return b*a; // uses implicit constructor Rational(a);
}

} // namespace siman
