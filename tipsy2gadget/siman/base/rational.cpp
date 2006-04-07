//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
//

#include "siman.hpp"

CRational::CRational() : p(0), q(1) {

}


CRational::CRational(int n) : p(n), q(1) {

}

CRational::CRational(int n, int d) : p(n), q(d) {

}

float CRational::flt() const {
  return ((float)p/(float)q);
}


double CRational::dbl() const {
  return ((double)p/(double)q);
}

int CRational::int_part() const {
  return (int) flt();
}

float CRational::flt_part() const {
  return flt() - (float) int_part();
}


double CRational::dbl_part() const {
  return dbl() - (double) int_part();
}


void CRational::rationalize() {
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

CRational CRational::operator+(const CRational &a) const {
  CRational out;
  out.p=p*a.q+q*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}


CRational CRational::operator-(const CRational &a) const {
  CRational out;
  out.p=p*a.q-q*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}

CRational CRational::operator-() const {
  CRational out;
  out.p=-p;
  out.q=q;
  return out;
}

CRational CRational::operator*(const CRational &a) const {
  CRational out;
  out.p=p*a.p;
  out.q=q*a.q;
  out.rationalize();
  return out;
}


CRational CRational::operator/(const CRational &a) const {
  CRational out;
  out.p=p*a.q;
  out.q=q*a.p;
  out.rationalize();
  return out;
}

void CRational::operator*=(const CRational &a) {
  p*=a.p;
  q*=a.q;
  rationalize();
}


void CRational::operator+=(const CRational &a) {
  p=p*a.q + a.p*q;
  q*=a.q;
  rationalize();
}


void CRational::operator-=(const CRational &a) {
  p=p*a.q - a.p*q;
  q*=a.q;
  rationalize();
}

void CRational::operator=(const CRational &c) {
  p=c.p;
  q=c.q;
}

bool CRational::operator==(const CRational &e) const {
  return ((e.p==p) && (e.q==q));
}

bool CRational::operator!=(const CRational &e) const {
  return ((e.p!=p) || (e.q!=q));
}


int CRational::gcd(int a, int b) {
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


CRational::CRational(string s) {
   
  fromString(s);

}

void CRational::fromString(string s) {
  
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

std::istream & operator>>(std::istream &is, CRational &in) {
  string s;
  is >> s;
 
  in.fromString(s);

  return is;
}

std::ostream & operator<<(std::ostream &os, const CRational &out) {
  if(out.q!=1) 
    os << out.p << "/" << out.q;
  else
    os << out.p;
  return os;
}

CRational operator*(int a, CRational b) {
  return b*a; // uses implicit constructor CRational(a);
}
