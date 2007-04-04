// units.cpp - part of SimAn Simulation Analysis Library
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

  std::istream & operator>>(std::istream &is, Unit &unit) {
    unit = Unit::streamIn(is);
    return is;
  }

  Unit Unit::streamIn(std::istream &is) {
    char line[1024];
    is.getline(line,1024);
    string line_s(line);
    istringstream input_isolated(line_s);
    return streamInIsolated(input_isolated);
  }

  Unit Unit::streamInIsolated(std::istream &is) {
    string u;
    
    if(!(is >> u)) return Unit(); // eof or similar
    
    string::size_type pos_exp = u.find("^",0);
    Rational exp = 1;
  
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
    if(u=="msol" || u=="m_sol")
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
    if(u=="myr")
      meaning = time_Myr;
    if(u=="Gyr") 
      meaning = time_Gyr;
    if(u=="erg" || u=="ergs")
      meaning = energy_erg;
    if(u=="ev" || u=="eV")
      meaning = energy_eV;
    if(u=="j")
      meaning = energy_J;
    if(u=="h")
      meaning = nodim_h;
    if(u=="a")
      meaning = nodim_a;
    if(u=="c")
      meaning = vel_light;
    if(u=="1+z") {
      meaning = nodim_a;
      exp = -exp;
    }
  
    double num;
    if(meaning==0) {
      num = atof(u.c_str());
      if(num==0.0) {
	throw UnknownUnit(u);
      } else {
	return streamIn(is)*num;
      }
    }


    return pow(Unit(meaning),exp)*streamIn(is);

  }

  std::ostream & operator<<(std::ostream &os, const Unit &unit) {
  
    bool previous=false;
  
    if(unit.exponent!=0 || unit.multiplier!=1) {
      os << (unit.multiplier*std::pow((double)10,unit.exponent)) << " ";
    }
    if(unit.mass_pwr!=0) {
      switch(unit.mass_type) {
      case Unit::mass_protons:
	os << "m_p";
	break;
      case Unit::mass_g:
	os << "g";
	break;
      case Unit::mass_kg:
	os << "kg";
	break;
      case Unit::mass_Msol:
	os << "Msol";
	break;
      }
      if(unit.mass_pwr!=1) os << "^" << unit.mass_pwr;
      previous=true;
    }

    if(unit.len_pwr!=0) {
      if(previous) os << " ";
      switch(unit.len_type) {
      case Unit::len_cm:
	os << "cm";
	break;
      case Unit::len_m:
	os << "m";
	break;
      case Unit::len_km:
	os << "km";
	break;
      case Unit::len_au:
	os << "au";
	break;
      case Unit::len_pc:
	os << "pc";
	break;
      case Unit::len_kpc:
	os << "kpc";
	break;
      case Unit::len_mpc:
	os << "Mpc";
	break;
      default:
	os << "(?len)";
	break;
      }

      if(unit.len_pwr!=1) os << "^" << unit.len_pwr;
      previous = true;

    }

    if(unit.time_pwr!=0) {
      if(previous) os << " ";
      switch(unit.time_type) {
      case Unit::time_s:
	os << "s";
	break;
      case Unit::time_yr:
	os << "yr";
	break;
      case Unit::time_Myr:
	os << "Myr";
	break;
      case Unit::time_Gyr:
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
 

    return os;


  }


  Unit::Unit(): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), scalefact_pwr(0),  h_pwr(0),  multiplier(1.),  exponent(0) {
  }

  Unit::Unit(const string &st): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0),  scalefact_pwr(0), h_pwr(0), multiplier(1.), exponent(0) {
    set(st);
  }

  Unit::Unit(double num) : len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0),  scalefact_pwr(0), h_pwr(0),multiplier(num), exponent(0) {
 
  }

  Unit::Unit(unsigned int len_typei,int len_pwri,unsigned int mass_typei,int mass_pwri,unsigned int time_typei,int time_pwri):
    len_type(len_typei), len_pwr(len_pwri), mass_type(mass_typei), mass_pwr(mass_pwri), time_type(time_typei), time_pwr(time_pwri), scalefact_pwr(0), h_pwr(0), multiplier(1.), exponent(0)  {

  }


  Unit::Unit(unsigned int len_typei,int len_pwri,unsigned int mass_typei,int mass_pwri,unsigned int time_typei,int time_pwri, double multi, int expi):
    len_type(len_typei), len_pwr(len_pwri), mass_type(mass_typei), mass_pwr(mass_pwri), time_type(time_typei), time_pwr(time_pwri), scalefact_pwr(0), h_pwr(0), multiplier(multi), exponent(expi)  {

  }

  Unit::Unit(unsigned int identifier): len_type(0), len_pwr(0), mass_type(0), mass_pwr(0), time_type(0), time_pwr(0), scalefact_pwr(0), h_pwr(0),  multiplier(1.), exponent(0) {
    set(identifier);
  }

  Unit::Unit(const Unit& copy):  len_type(copy.len_type), len_pwr(copy.len_pwr), mass_type(copy.mass_type), mass_pwr(copy.mass_pwr), time_type(copy.time_type), time_pwr(copy.time_pwr), scalefact_pwr(copy.scalefact_pwr), h_pwr(copy.h_pwr), multiplier(copy.multiplier), exponent(copy.exponent) {

  }

  Unit::Unit(ifstream *file, int vernum) {
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


  void Unit::nativeWrite(ofstream *file) {
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


  void Unit::set(const string & st) {
    istringstream iss(st);
    (*this)=streamInIsolated(iss);
  }

  string Unit::toString() const {
    ostringstream oss;
    oss << (*this);
    return oss.str();
  }

  void Unit::set(unsigned int identifier) {

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
      case vel_light:
	len_type = len_m;
	time_type = time_s;
	len_pwr = 1;
	time_pwr = -1;
	multiplier = 2.99792458;
	exponent = 8;
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
      case energy_eV:
	mass_type = mass_g;
	time_type = time_s;
	len_type = len_cm;
	mass_pwr = 1;
	len_pwr = 2;
	time_pwr = -2;
	multiplier = 1.60217646;
	exponent = -12;
	break;
      }
    }
  }


  bool Unit::canConvert(const Unit & to) const {
    if(to.len_pwr == len_pwr && to.mass_pwr == mass_pwr && to.time_pwr == time_pwr)
      return true;
    else
      return false;
  }

  void Unit::assertDimensionless() {
  
    if(!canConvert(Unit())) {
      throw(UnitsError(*this,Unit()));
    }
  }

  double Unit::convertTo(const Unit & to, const SimSnap *pSim) const {

    if(!canConvert(to)) {
      throw(UnitsError(*this,to));
    }
  
    Unit ratio_unit = (*this)/to;

    int c_exponent = ratio_unit.exponent;
    double mantissa_ratio = ratio_unit.multiplier;

    ratio_unit.assertDimensionless();

    if((h_pwr-to.h_pwr)!=0 || (scalefact_pwr-to.scalefact_pwr)!=0) {
      float a=0.5, h=0.7;
      if(pSim==NULL) {
	cerr << "Unit: WARNING - no context given for cosmological conversion; assuming a=" << a << ", h=" << h << endl;
      } else {
	a = 1./(1.+pSim->getRedshift());
	h = pSim->getHubble();
      }
      mantissa_ratio*=std::pow((double)a,(scalefact_pwr-to.scalefact_pwr).dbl());
      mantissa_ratio*=std::pow((double)h,(h_pwr-to.h_pwr).dbl());
    }
    // cout << mantissa_ratio << " " << exponent << endl;
    if(abs(c_exponent)>60)
      cerr << "Unit: warning - large exponent convertion ratio produced, may lead to inaccuracies" << endl;

    return mantissa_ratio*std::pow((double)10.,c_exponent);
  }

  double Unit::convertTo(const unsigned int to, const SimSnap *pSim) const {
    return convertTo(Unit(to), pSim);
  }

  void Unit::operator*=(const double mul) {
    multiplier*=mul;
  }

  void Unit::operator/=(const double div) {
    multiplier/=div;
  }

  void Unit::operator*=(const Unit &mul) {
    (*this) = (*this)*mul; // inefficient, revise!
  }


  void Unit::operator/=(const Unit &div) {
    (*this) = (*this)/div; // inefficient, revise!
  }

  double Unit::getMantissa(unsigned int type) {

    static const double unit_mantissa_len[7]={1,1,1,1.49598,3.08568025,3.08568025,3.08568025};
    static const double unit_mantissa_mass[4]={1.67262158,1,1,1.98892};
    static const double unit_mantissa_time[16]={1,3.1556926,3.1556926,3.1556926};

    if((type & len_mask)==type) {
      return unit_mantissa_len[(type>>len_shift)-1];
    }
    if((type & mass_mask)==type) {
      return unit_mantissa_mass[(type>>mass_shift)-1];
    }
    if((type & time_mask)==type) {
      return unit_mantissa_time[(type>>time_shift)-1];
    }
    UnknownUnit e;
    throw(e);
  }

  int Unit::getExponent(unsigned int type) {
    static const int unit_exponent_len[7]={-2,0,3,11,16,19,22};
    static const int unit_exponent_mass[4]={-27,-3,0,30};
    static const int unit_exponent_time[4]={0,7,13,16};

    if((type & len_mask)==type) {
      return unit_exponent_len[(type>>len_shift)-1];
    }
    if((type & mass_mask)==type) {
      return unit_exponent_mass[(type>>mass_shift)-1];
    }
    if((type & time_mask)==type) {
      return unit_exponent_time[(type>>time_shift)-1];
    }
    UnknownUnit e;
    throw(e);
  }

  void Unit::getMultiplicationRule(unsigned int from_type,unsigned int to_type, Rational pwr, double &multiplier, int &exponent) const {
  
    multiplier = std::pow(getMantissa(from_type)/getMantissa(to_type),pwr.dbl());
    int exp_shift = getExponent(from_type) - getExponent(to_type);  
  
    exponent=exp_shift*pwr.int_part();
    multiplier*=std::pow(10.,exp_shift*pwr.dbl_part());

  }



  Unit Unit::operator*(const Unit &mul) const {
    Unit ret;
  
    ret.multiplier = multiplier*mul.multiplier;
    ret.exponent = exponent + mul.exponent;

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


  Unit Unit::operator/(const Unit &div) const {
    return (*this)*(1./div);
  }

  Unit operator/(const double c, const Unit &unit) {
    Unit ret(unit);
    ret.multiplier=c/(unit.multiplier);
    ret.exponent=-unit.exponent;
    ret.mass_pwr=-unit.mass_pwr;
    ret.len_pwr=-unit.len_pwr;
    ret.time_pwr=-unit.time_pwr;
    ret.h_pwr=-unit.h_pwr;
    ret.scalefact_pwr = -unit.scalefact_pwr;

    return ret;
  }

  bool Unit::operator==(const Unit & comp) const {
    if(!canConvert(comp)) return false;

    if(convertTo(comp)==1.) return true;

    return false;

  }

  Unit pow(const Unit & unit, const Rational &power) {
    Unit ret(unit);

    int integ_power = power.int_part();
    double remainder_power = power.dbl_part();

  

    ret.exponent*=integ_power;
    ret.multiplier= std::pow(ret.multiplier,power.dbl())*std::pow(10.,((double)unit.exponent)*remainder_power);

    ret.mass_pwr*=power;
    ret.len_pwr*=power;
    ret.time_pwr*=power;
    ret.scalefact_pwr*=power;
    ret.h_pwr*=power;
    return ret;
  }

} // namespace siman
