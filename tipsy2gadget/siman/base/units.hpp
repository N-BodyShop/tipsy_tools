//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __UNITS_H_INCLUDED

#define __UNITS_H_INCLUDED

#include "siman.hpp"

class CSimSnap;

namespace units {

  // Some constants

  const double protonMassInG = 1.6726e-24;
  const double boltzmannInErgPerK = 1.3806503e-16;
  const double planckInEvS = 4.136e-15;
  const double planckInErgS21 = 6.62e-6; // in erg seconds / 10^(-21)
  const double evInErg21 = 1.6e9; // 1ev in 10^(-21)ergs

  // LENGTHS:

  const double KpcInCm = 3.086e21;

  // ENERGIES:

  const double kmPerS2InJoules = 1.e6;
  const double kmPerS2InErgs = 1.e13;

  // VELOCITIES:

  const double cInKpcPerS = 9.716e-12;

  // COLUMN DENSITIES:

  // Msol/kpc^2 -> kg / m^2

  const double kgPerM2 = 2.089e-9;

  // Msol/kpc^2 -> protons/cm^2

  const double protonsPerCm2 = 1.249e14;


  // DENSITIES:

  // Msol/kpc^3 -> protons/cm^3

  const double protonsPerCm3 = 4.0473e-8;





  // on-the-fly conversion data:  

  const unsigned int len_shift = 0;
  const unsigned int len_mask = 15;
  const unsigned int len_cm = 1;
  const unsigned int len_m = 2;
  const unsigned int len_km = 3;
  const unsigned int len_au = 4;
  const unsigned int len_pc = 5;
  const unsigned int len_kpc = 6;
  const unsigned int len_mpc = 7;

  const double unit_mantissa_len[7]={1,1,1,1.49598,3.08568025,3.08568025,3.08568025};
  const int unit_exponent_len[7]={-2,0,3,11,16,19,22};

 
  const unsigned int mass_shift = 4;
  const unsigned int mass_mask = 15 << 4;
  const unsigned int mass_protons = 1 << 4;
  const unsigned int mass_g = 2 << 4;
  const unsigned int mass_kg = 3 << 4;
  const unsigned int mass_Msol = 4 << 4;

  const double unit_mantissa_mass[4]={1.67262158,1,1,1.98892};
  const int unit_exponent_mass[4]={-27,-3,0,30};

  const unsigned int time_shift = 8;
  const unsigned int time_mask = 15 << 8;
  const unsigned int time_s = 1 << 8;
  const unsigned int time_yr = 2 << 8;
  const unsigned int time_Myr = 3 << 8;
  const unsigned int time_Gyr = 4 << 8;

  const double unit_mantissa_time[16]={1,3.1556926,3.1556926,3.1556926};
  const int unit_exponent_time[4]={0,7,13,16};

  // derived units

  const unsigned int vel_kmPerS = 1 << 12;
  const unsigned int vel_mPerS = 2 << 12;
  const unsigned int vel_cmPerS = 3 << 12;
  const unsigned int den_protonsPerCm3 = 4 << 12;
  const unsigned int den_MsolPerKpc3 = 5 << 12;
  const unsigned int energy_J = 6 << 12;
  const unsigned int energy_erg=  7 << 12;
  const unsigned int den_gramsPerCm3 = 8 << 12;
  const unsigned int colden_protonsPerCm2 = 9 << 12;
  
  // dimensionless units

  const unsigned int nodim_h = 1 << 16;
  const unsigned int nodim_a = 2 << 16;
  const unsigned int nodim_mask = 15 << 16;

  class CUnit {
  public:
    unsigned int len_type;
    CRational len_pwr;
    unsigned int mass_type;
    CRational mass_pwr;
    unsigned int time_type;
    CRational time_pwr;

    // dimensionless scalings:
    CRational scalefact_pwr;
    CRational h_pwr;

    double multiplier;
    int exponent;

    CUnit(unsigned int len_type,int len_pwr,unsigned int mass_type,int mass_pwr,unsigned int time_type,int time_pwr);
    CUnit(unsigned int len_type,int len_pwr,unsigned int mass_type,int mass_pwr,unsigned int time_type,int time_pwr,double multiplier, int exponent);

    /// constuct a CUnit corresponding to an identifier (see namespace units). 
    /// Explicit to prevent confusing automatic casts (1!=CUnit(1)...)
    explicit CUnit(unsigned int identifier);

    CUnit(const string &st);
    CUnit(const CUnit & copy);
    CUnit(double num);
    CUnit(std::ifstream *file, int vernum);
    CUnit();

    static CUnit  streamIn(std::istream &is, int paren=0);

    void set(unsigned int identifier);
    void set(const string & st);

    void nativeWrite(std::ofstream *file);

    /// throws CUnitError if this object has dimensions in length, time or mass
    void assertDimensionless();

    /// returns true if this object and \param to are compatible
    bool canConvert(const CUnit &to) const;

    /// 
    double convertTo(const CUnit &to, CSimSnap *pSim=NULL) const;
    double convertTo(const unsigned int to, CSimSnap *pSim=NULL) const;

    void operator*=(const double mul);
    void operator/=(const double div);

    bool operator==(const CUnit &check) const;

    CUnit operator*(const CUnit &mul) const;
    CUnit operator/(const CUnit &div) const;

  private:
    /// Internal multiplication rule access
    /// @param[in] type_shift - the overall bitshift for the type in question
    /// @param[in] from_type - the subtype (i.e. unit) to convert from
    /// @param[in] to_type - the subtype (unit) to convert TO
    /// @param[out] multiplier - the multiplier to apply to the mantissa
    /// @param[out] exponent - the power of 10 to apply
    void getMultiplicationRule(unsigned int from_type,unsigned int to_type, CRational pwr, double &multiplier, int &exponent) const;
    
    /// Internal multiplication rule access
    /// @param[in] type - the unit identifier
    /// @returns the mantissa of the unit in SI (as currently defined in units.hpp)
    static double getMantissa(int type);

    /// Internal multiplication rule access
    /// @param[in] type - the unit identifier
    /// @returns the exponent of the unit in SI (as currently defined in units.hpp)
    static int getExponent(int type);
  };
  
  
  std::ostream & operator<<(std::ostream &os, const CUnit &unit);
  std::istream & operator>>(std::istream &is, CUnit &unit);
  CUnit operator/(const double c, const CUnit &unit);
  
};

units::CUnit pow(const units::CUnit &unit, const CRational & power);


#define PI 3.141592

#endif
