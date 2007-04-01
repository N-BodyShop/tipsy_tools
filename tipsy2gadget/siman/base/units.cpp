//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

using namespace units;


std::istream & units::operator>>(std::istream &is, CUnit &unit) {
  unit = CUnit::streamIn(is);
  return is;
}

CUnit CUnit::streamIn(std::istream &is, int paren) {
  string u;
  is >> u;
  

  if(u[0]=='(') {
    paren++;
    u = u.substr(1,u.length()-1);
  }

  if(u[u.length()-1]==')') {
    paren--;
    u=u.substr(0,u.length()-1);
  }

  string::size_type pos_exp = u.find("^",0);
  CRational exp = 1;
  
  if(pos_exp!=string::npos) {
    istringstream is2(u.substr(pos_exp+1,u.length()-pos_exp-1));
    is2 >> exp;
    u=u.substr(0,pos_exp);
  }

  transform(u.begin(),u.end(),u.begin(),(int(*)(int))tolower);
  // preprocess u, remove parenthesis and split out exponent

  unsigned int meaning=0;
  if(u=="m_p")
    meaning = mass_protons;
  if(u=="g") 
    meaning = mass_g;
  if(u=="kg")
    meaning = mass_kg;
  if(u=="msol")
    meaning = mass_Msol;
  if(u=="cm")
    meaning = len_cm;
  if(u=="m")
    meaning = len_m;
  if(u=="km")
    meaning = len_km;
  if(u=="au")
    meaning = len_au;
  if(u=="pc")
    meaning = len_pc;
  if(u=="kpc")
    meaning = len_kpc;
  if(u=="mpc")
    meaning = len_mpc;
  if(u=="s")
    meaning = time_s;
  if(u=="yr")
    meaning = time_yr;
  if(u=="Myr")
    meaning = time_Myr;
  if(u=="Gyr") 
    meaning = time_Gyr;
  if(u=="h")
    meaning = nodim_h;
  if(u=="a")
    meaning = nodim_a;
  if(u=="1+z") {
    meaning = nodim_a;
    exp = -exp;
  }
  
  double num;
  if(meaning==0) {
    num = atof(u.c_str());
    if(num==0.0) {
      cerr << "CUnits: didn't understand " << u;
    } else {
      if(paren==0) 
	return CUnit()*num;
      else
	return streamIn(is,paren)*num;
    }
  }

  if(paren==0)
    return pow(CUnit(meaning),exp);
  else
    return pow(CUnit(meaning),exp)*streamIn(is,paren);

}

std::ostream & units::operator<<(std::ostream &os, const CUnit &unit) {
  
  bool previous=false;
  os << "(";
  if(unit.exponent!=0 || unit.multiplier!=1) {
    os << (unit.multiplier*pow((double)10,unit.exponent)) << " ";
  }
  if(unit.mass_pwr!=0) {
    switch(unit.mass_type) {
    case mass_protons:
      os << "m_p";
      break;
    case mass_g:
      os << "g";
      break;
    case mass_kg:
      os << "kg";
      break;
    case mass_Msol:
      os << "Msol";
      break;
    }
    if(unit.mass_pwr!=1) os << "^" << unit.mass_pwr;
    previous=true;
  }

  if(unit.len_pwr!=0) {
    if(previous) os << " ";
    switch(unit.len_type) {
    case len_cm:
      os << "cm";
      break;
    case len_m:
      os << "m";
      break;
    case len_km:
      os << "km";
      break;
    case len_au:
      os << "au";
      break;
    case len_pc:
      os << "pc";
      break;
    case len_kpc:
      os << "kpc";
      break;
    case len_mpc:
      os << "Mpc";
      break;
    default:
      os << "(len)";
      break;
    }

    if(unit.len_pwr!=1) os << "^" << unit.len_pwr;
    previous = true;

  }

  if(unit.time_pwr!=0) {
    if(previous) os << " ";
    switch(unit.time_type) {
    case time_s:
      os << "s";
      break;
    case time_yr:
      os << "yr";
      break;
    case time_Myr:
      os << "Myr";
      break;
    case time_Gyr:
      os << "Gyr";
      break;
    }
    if(unit.time_pwr!=1) os << "^" << unit.time_pwr;
    previous = true;
  }

  if(unit.scalefact_pwr!=0) {
    if(previous) os << " ";
    os << "a";
    if(unit.scalefact_pwr!=1) os << "^" << unit.scalefact_pwr;
    previous = true;
  }

  if(unit.h_pwr!=0) {
    if(previous) os << " ";
    os << "h";
    if(unit.h_pwr!=1) os << "^" << unit.h_pwr;
    previous = true;
  }
 
  os << ")";

  return os;


}


CUnit::CUnit(): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), multiplier(1.), scalefact_pwr(0), h_pwr(0), exponent(0) {
}

CUnit::CUnit(const string &st): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), multiplier(1.), scalefact_pwr(0), h_pwr(0), exponent(0) {
  set(st);
}

CUnit::CUnit(double num) : len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), multiplier(num), scalefact_pwr(0), h_pwr(0), exponent(0) {
 
}

CUnit::CUnit(unsigned int len_typei,int len_pwri,unsigned int mass_typei,int mass_pwri,unsigned int time_typei,int time_pwri):
  len_type(len_typei), len_pwr(len_pwri), mass_type(mass_typei), mass_pwr(mass_pwri), time_type(time_typei), time_pwr(time_pwri), scalefact_pwr(0), h_pwr(0), multiplier(1.), exponent(0)  {

}


CUnit::CUnit(unsigned int len_typei,int len_pwri,unsigned int mass_typei,int mass_pwri,unsigned int time_typei,int time_pwri, double multi, int expi):
  len_type(len_typei), len_pwr(len_pwri), mass_type(mass_typei), mass_pwr(mass_pwri), time_type(time_typei), time_pwr(time_pwri), multiplier(multi), exponent(expi), scalefact_pwr(0), h_pwr(0)  {

}

CUnit::CUnit(unsigned int identifier): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), multiplier(1.), exponent(0), scalefact_pwr(0), h_pwr(0) {
  set(identifier);
}

CUnit::CUnit(const CUnit& copy):  len_type(copy.len_type), len_pwr(copy.len_pwr), mass_type(copy.mass_type), mass_pwr(copy.mass_pwr), time_type(copy.time_type), time_pwr(copy.time_pwr), multiplier(copy.multiplier), exponent(copy.exponent), scalefact_pwr(copy.scalefact_pwr), h_pwr(copy.h_pwr) {

}

CUnit::CUnit(ifstream *file, int vernum) {
  file->read((char*) &len_type,sizeof(int));
  file->read((char*) &mass_type,sizeof(int));
  file->read((char*) &time_type,sizeof(int));
  if(vernum<3) {
    // older files just have integer powers
    int len_pwr_i, mass_pwr_i, time_pwr_i, scalefact_pwr_i, h_pwr_i;
    file->read((char*) &len_pwr_i,sizeof(int));
    file->read((char*) &mass_pwr_i,sizeof(int));
    file->read((char*) &time_pwr_i,sizeof(int));
    file->read((char*) &scalefact_pwr_i,sizeof(int));
    file->read((char*) &h_pwr_i,sizeof(int));
    len_pwr=len_pwr_i;
    mass_pwr=mass_pwr_i;
    time_pwr=time_pwr_i;
    scalefact_pwr=scalefact_pwr_i;
    h_pwr=h_pwr_i;
  } else {
    
    file->read((char*) &len_pwr.p,sizeof(int));
    file->read((char*) &mass_pwr.p,sizeof(int));
    file->read((char*) &time_pwr.p,sizeof(int));
    file->read((char*) &scalefact_pwr.p,sizeof(int));
    file->read((char*) &h_pwr.p,sizeof(int));
    
    file->read((char*) &len_pwr.q,sizeof(int));
    file->read((char*) &mass_pwr.q,sizeof(int));
    file->read((char*) &time_pwr.q,sizeof(int));
    file->read((char*) &scalefact_pwr.q,sizeof(int));
    file->read((char*) &h_pwr.q,sizeof(int));
  }

  file->read((char*) &multiplier,sizeof(double));
  file->read((char*) &exponent,sizeof(int));
}


void CUnit::nativeWrite(ofstream *file) {
  file->write((char*) &len_type,sizeof(int));
  file->write((char*) &mass_type,sizeof(int));
  file->write((char*) &time_type,sizeof(int));

  file->write((char*) &len_pwr.p,sizeof(int));
  file->write((char*) &mass_pwr.p,sizeof(int));
  file->write((char*) &time_pwr.p,sizeof(int));
  file->write((char*) &scalefact_pwr.p,sizeof(int));
  file->write((char*) &h_pwr.p,sizeof(int));


  file->write((char*) &len_pwr.q,sizeof(int));
  file->write((char*) &mass_pwr.q,sizeof(int));
  file->write((char*) &time_pwr.q,sizeof(int));
  file->write((char*) &scalefact_pwr.q,sizeof(int));
  file->write((char*) &h_pwr.q,sizeof(int));

  file->write((char*) &multiplier,sizeof(double));
  file->write((char*) &exponent,sizeof(int));
}


void CUnit::set(const string & st) {
  istringstream iss(st);
  (*this)=streamIn(iss);
}


void CUnit::set(unsigned int identifier) {

  scalefact_pwr =0 ;
  h_pwr =0;

  multiplier=1.;
  exponent=0;

   if((identifier & len_mask)!= 0) {
  
    len_pwr = 1;
    len_type = (identifier & len_mask);
    identifier-=len_type;
  }

  if((identifier & mass_mask)!=0) {
  
    mass_pwr = 1;
    mass_type = (identifier & mass_mask);
    identifier-=mass_type;
   
  }

  if((identifier & time_mask)!=0) {
  
    time_pwr = 1;
    time_type = (identifier & time_mask);
    identifier-=time_type;
   
  }

  if((identifier & nodim_h)!=0) {
    h_pwr = 1;
    identifier-=nodim_h;
  }

  if((identifier & nodim_a)!=0) {
    scalefact_pwr = 1;
    identifier-=nodim_a;
  }

  if(identifier!=0) {

    // derived units

    switch(identifier) {
    case den_protonsPerCm3:
      mass_type = mass_protons;
      len_type = len_cm;
      mass_pwr = 1;
      len_pwr = -3;
      break;
    case den_gramsPerCm3:
      mass_type = mass_g;
      len_type = len_cm;
      mass_pwr=1;
      len_pwr=-3;
      break;
    case den_MsolPerKpc3:
      mass_type = mass_Msol;
      len_type = len_kpc;
      mass_pwr = 1;
      len_pwr = -3;
      break;
    case colden_protonsPerCm2:
      mass_type = mass_protons;
      len_type = len_cm;
      mass_pwr = 1;
      len_pwr = -2;
      break;
    case vel_kmPerS:
      len_type = len_km;
      time_type = time_s;
      len_pwr = 1;
      time_pwr = -1;
      break;
    case vel_cmPerS:
      len_type = len_cm;
      time_type = time_s;
      len_pwr = 1;
      time_pwr = -1;
      break;
    case energy_J:
      mass_type = mass_kg;
      time_type = time_s;
      len_type = len_m;
      mass_pwr = 1;
      len_pwr = 2;
      time_pwr = -2;
      break;
    case energy_erg:
      mass_type = mass_g;
      time_type = time_s;
      len_type = len_cm;
      mass_pwr = 1;
      len_pwr = 2;
      time_pwr = -2;
      break;
    }
  }
}


bool CUnit::canConvert(const CUnit & to) const {
  if(to.len_pwr == len_pwr && to.mass_pwr == mass_pwr && to.time_pwr == time_pwr)
    return true;
  else
    return false;
}

void CUnit::assertDimensionless() {
  
  if(!canConvert(CUnit())) {
    throw(CUnitsError(*this,CUnit()));
  }
}

double CUnit::convertTo(const CUnit & to, CSimSnap *pSim) const {

  if(!canConvert(to)) {
    throw(CUnitsError(*this,to));
  }
  
  CUnit ratio_unit = (*this)/to;

  int c_exponent = ratio_unit.exponent;
  double mantissa_ratio = ratio_unit.multiplier;

  ratio_unit.assertDimensionless();

  if((h_pwr-to.h_pwr)!=0 || (scalefact_pwr-to.scalefact_pwr)!=0) {
    float a=0.5, h=0.7;
    if(pSim==NULL) {
      cerr << "CUnit: WARNING - no context given for cosmological conversion; assuming a=" << a << ", h=" << h << endl;
    } else {
      a = 1./(1.+pSim->getRedshift());
      h = pSim->getHubble();
    }
    mantissa_ratio*=pow((double)a,(scalefact_pwr-to.scalefact_pwr).dbl());
    mantissa_ratio*=pow((double)h,(h_pwr-to.h_pwr).dbl());
  }
  // cout << mantissa_ratio << " " << exponent << endl;
  if(abs(c_exponent)>60)
    cerr << "CUnit: warning - large exponent convertion ratio produced, may lead to inaccuracies" << endl;

  return mantissa_ratio*pow((double)10.,c_exponent);
}

double CUnit::convertTo(const unsigned int to, CSimSnap *pSim) const {
  return convertTo(CUnit(to), pSim);
}

void CUnit::operator*=(const double mul) {
  multiplier*=mul;
}

void CUnit::operator/=(const double div) {
  multiplier/=div;
}

double CUnit::getMantissa(int type) {
  if((type & len_mask)==type) {
    return unit_mantissa_len[(type>>len_shift)-1];
  }
  if((type & mass_mask)==type) {
    return unit_mantissa_mass[(type>>mass_shift)-1];
  }
  if((type & time_mask)==type) {
    return unit_mantissa_time[(type>>time_shift)-1];
  }
  CUnknownUnit e;
  throw(e);
}

int CUnit::getExponent(int type) {
  if((type & len_mask)==type) {
    return unit_exponent_len[(type>>len_shift)-1];
  }
  if((type & mass_mask)==type) {
    return unit_exponent_mass[(type>>mass_shift)-1];
  }
  if((type & time_mask)==type) {
    return unit_exponent_time[(type>>time_shift)-1];
  }
  CUnknownUnit e;
  throw(e);
}

void CUnit::getMultiplicationRule(unsigned int from_type,unsigned int to_type, CRational pwr, double &multiplier, int &exponent) const {
  
  double mul_pwr = pwr.dbl();

  multiplier = pow(getMantissa(from_type)/getMantissa(to_type),pwr.dbl());
  int exp_shift = getExponent(from_type) - getExponent(to_type);  
  
  exponent=exp_shift*pwr.int_part();
  multiplier*=pow(10.,exp_shift*pwr.dbl_part());

  //  cout << "GMR: " << getMantissa(from_type) << " " << getMantissa(to_type) << " " << pwr << " " << multiplier << " " << exponent << " " << multiplier << endl;
}



CUnit CUnit::operator*(const CUnit &mul) const {
  CUnit ret;
  
  ret.multiplier = multiplier*mul.multiplier;
  ret.exponent = exponent * mul.exponent;

  ret.len_pwr=len_pwr+mul.len_pwr;

  if(len_pwr!=0)
    ret.len_type=len_type;
  else if(mul.len_pwr!=0)
    ret.len_type=mul.len_type;
  
  ret.mass_pwr=mass_pwr+mul.mass_pwr;

  if(mass_pwr!=0)
    ret.mass_type=mass_type;
  else if(mul.mass_pwr!=0)
    ret.mass_type=mul.mass_type;
 

  
  ret.time_pwr=time_pwr+mul.time_pwr;

  if(time_pwr!=0)
    ret.time_type=time_type;
  else if(mul.time_pwr!=0 && ret.time_pwr!=0)
    ret.time_type=mul.time_type;
 
  
  // if return mass type is not the same as mul mass type,
  // do conversion
  
  if(ret.mass_type!=mul.mass_type && mul.mass_pwr!=0) {
    double im;
    int ie;
    getMultiplicationRule(mul.mass_type,ret.mass_type,mul.mass_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;
  }

  // if return mass type is not same as (*this) mass type,
  // do conversion

  if(ret.mass_type!=mass_type && mass_pwr!=0) {
    double im;
    int ie;
    getMultiplicationRule(mass_type,ret.mass_type,mass_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;
  }


  
  // if return time type is not the same as mul time type,
  // do conversion
  
  if(ret.time_type!=mul.time_type && mul.time_pwr!=0) {
    
    double im;
    int ie;
    getMultiplicationRule(mul.time_type,ret.time_type,mul.time_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;

   
  }

  // if return time type is not same as (*this) time type,
  // do conversion

  if(ret.time_type!=time_type && time_pwr!=0) {
    
    double im;
    int ie;
    getMultiplicationRule(time_type,ret.time_type,time_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;

  }

  
  // if return len type is not the same as mul len type,
  // do conversion
  
  if(ret.len_type!=mul.len_type && mul.len_pwr!=0) {

    double im;
    int ie;
    getMultiplicationRule(mul.len_type,ret.len_type,mul.len_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;

  }

  // if return len type is not same as (*this) len type,
  // do conversion

  if(ret.len_type!=len_type && len_pwr!=0) {
    
    double im;
    int ie;
    getMultiplicationRule(len_type,ret.len_type,len_pwr,im,ie);
    ret.multiplier*=im;
    ret.exponent+=ie;

  }

  if(ret.time_pwr==0) ret.time_type=0;
  if(ret.mass_pwr==0) ret.mass_type=0;
  if(ret.len_pwr==0) ret.len_type=0;

  ret.scalefact_pwr = scalefact_pwr + mul.scalefact_pwr;
  ret.h_pwr = h_pwr + mul.h_pwr;

  return ret;
}


CUnit CUnit::operator/(const CUnit &div) const {
  return (*this)*(1./div);
}

CUnit units::operator/(const double c, const CUnit &unit) {
  CUnit ret(unit);
  ret.multiplier=c/(unit.multiplier);
  ret.exponent=-unit.exponent;
  ret.mass_pwr=-unit.mass_pwr;
  ret.len_pwr=-unit.len_pwr;
  ret.time_pwr=-unit.time_pwr;
  ret.h_pwr=-unit.h_pwr;
  ret.scalefact_pwr = -unit.scalefact_pwr;

  return ret;
}

bool CUnit::operator==(const CUnit & comp) const {
  if(!canConvert(comp)) return false;

  if(convertTo(comp)==1.) return true;

  return false;

}

units::CUnit pow(const units::CUnit & unit, const CRational &power) {
  CUnit ret(unit);

  int integ_power = power.int_part();
  double remainder_power = power.dbl_part();

  

  ret.exponent*=integ_power;
  ret.multiplier*= pow(ret.multiplier,power.dbl())*pow(10.,((double)ret.exponent)*remainder_power);

  ret.mass_pwr*=power;
  ret.len_pwr*=power;
  ret.time_pwr*=power;
  ret.scalefact_pwr*=power;
  ret.h_pwr*=power;
  return ret;
}

