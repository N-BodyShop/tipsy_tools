// units.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __UNITS_H_INCLUDED

#define __UNITS_H_INCLUDED

#include "siman.hpp"

namespace siman {

class SimSnap;

class Unit {
public:
  unsigned int len_type;
  Rational len_pwr;
  unsigned int mass_type;
  Rational mass_pwr;
  unsigned int time_type;
  Rational time_pwr;
  
  // dimensionless scalings:
  Rational scalefact_pwr;
  Rational h_pwr;
  
  double multiplier;
  int exponent;
  
  Unit(unsigned int len_type,int len_pwr,unsigned int mass_type,int mass_pwr,unsigned int time_type,int time_pwr);
  Unit(unsigned int len_type,int len_pwr,unsigned int mass_type,int mass_pwr,unsigned int time_type,int time_pwr,double multiplier, int exponent);
  
  /// constuct a Unit corresponding to an identifier (see namespace units). 
  /// Explicit to prevent confusing automatic casts (1!=Unit(1)...)
  explicit Unit(unsigned int identifier);

  Unit(const std::string &st);
  Unit(const Unit & copy);
  Unit(double num);
  Unit(std::ifstream *file, int vernum);
  Unit();
  
  /// convert to string - mainly for python wrapper (does not seem to work with << as expected - .def(str(self)) fails)
  std::string toString() const;
 
  static Unit streamIn(std::istream &is); ///< streams in a unit, ending at end of line
  static Unit streamInIsolated(std::istream &is); ///< streams in a unit, ending at EOF
  void set(unsigned int identifier); ///< sets a unit to an internal identifier
  void set(const std::string & st); ///< sets a unit to a given string
  
  void nativeWrite(std::ofstream *file);
  
  /// throws UnitError if this object has dimensions in length, time or mass
  void assertDimensionless();
  
  /// returns true if this object and \param to are compatible
  bool canConvert(const Unit &to) const;
  
  /// 
  double convertTo(const Unit &to, const SimSnap *pSim=NULL) const;
  double convertTo(const unsigned int to, const SimSnap *pSim=NULL) const;
  
  void operator*=(const double mul);
  void operator/=(const double div);
  void operator*=(const Unit &mul);
  void operator/=(const Unit &div);
  
  bool operator==(const Unit &check) const;
  
  Unit operator*(const Unit &mul) const;
  Unit operator/(const Unit &div) const;
  
  




  // on-the-fly conversion data:  

  static const unsigned int len_shift = 0;
  static const unsigned int len_mask = 15;
  static const unsigned int len_cm = 1;
  static const unsigned int len_m = 2;
  static const unsigned int len_km = 3;
  static const unsigned int len_au = 4;
  static const unsigned int len_pc = 5;
  static const unsigned int len_kpc = 6;
  static const unsigned int len_mpc = 7;

 

 
  static const unsigned int mass_shift = 4;
  static const unsigned int mass_mask = 15 << 4;
  static const unsigned int mass_protons = 1 << 4;
  static const unsigned int mass_g = 2 << 4;
  static const unsigned int mass_kg = 3 << 4;
  static const unsigned int mass_Msol = 4 << 4;

 
 

  static const unsigned int time_shift = 8;
  static const unsigned int time_mask = 15 << 8;
  static const unsigned int time_s = 1 << 8;
  static const unsigned int time_yr = 2 << 8;
  static const unsigned int time_Myr = 3 << 8;
  static const unsigned int time_Gyr = 4 << 8;

  
  

  // derived units

  static const unsigned int vel_kmPerS = 1 << 12;
  static const unsigned int vel_mPerS = 2 << 12;
  static const unsigned int vel_cmPerS = 3 << 12;
  static const unsigned int den_protonsPerCm3 = 4 << 12;
  static const unsigned int den_MsolPerKpc3 = 5 << 12;
  static const unsigned int energy_J = 6 << 12;
  static const unsigned int energy_erg=  7 << 12;
  static const unsigned int den_gramsPerCm3 = 8 << 12;
  static const unsigned int colden_protonsPerCm2 = 9 << 12;
  static const unsigned int vel_light = 10 << 12;
  static const unsigned int energy_eV = 11 << 12;

  // dimensionless units

  static const unsigned int nodim_h = 1 << 16;
  static const unsigned int nodim_a = 2 << 16;
  static const unsigned int nodim_mask = 15 << 16;


  private:
  /// Internal multiplication rule access
  /// @param[in] type_shift - the overall bitshift for the type in question
  /// @param[in] from_type - the subtype (i.e. unit) to convert from
  /// @param[in] to_type - the subtype (unit) to convert TO
  /// @param[out] multiplier - the multiplier to apply to the mantissa
  /// @param[out] exponent - the power of 10 to apply
  void getMultiplicationRule(unsigned int from_type,unsigned int to_type, Rational pwr, double &multiplier, int &exponent) const;
  
  /// Internal multiplication rule access
  /// @param[in] type - the unit identifier
  /// @returns the mantissa of the unit in SI (as currently defined in units.hpp)
  static double getMantissa(unsigned int type);
  
  /// Internal multiplication rule access
  /// @param[in] type - the unit identifier
  /// @returns the exponent of the unit in SI (as currently defined in units.hpp)
  static int getExponent(unsigned int type);
};



std::istream & operator>>(std::istream &is, Unit &unit);
Unit operator/(const double c, const Unit &unit);
std::ostream & operator<<(std::ostream &os, const Unit &unit);
Unit pow(const Unit &unit, const Rational & power);

} // namespace siman

#endif
