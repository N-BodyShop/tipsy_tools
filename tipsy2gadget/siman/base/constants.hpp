// constants.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __CONSTANTS_H_INCLUDED

#define __CONSTANTS_H_INCLUDED


namespace siman {
  namespace constants {
  
  // Some constants

  static const double protonMassInG = 1.6726e-24;
  static const double boltzmannInErgPerK = 1.3806503e-16;
  static const double planckInEvS = 4.136e-15;
  static const double planckInErgS21 = 6.62e-6; // in erg seconds / 10^(-21)
  static const double evInErg21 = 1.6e9; // 1ev in 10^(-21)ergs

  // LENGTHS:

  static const double KpcInCm = 3.086e21;

  // ENERGIES:

  static const double kmPerS2InJoules = 1.e6;
  static const double kmPerS2InErgs = 1.e13;

  // VELOCITIES:

  static const double cInKpcPerS = 9.716e-12;

  // COLUMN DENSITIES:

  // Msol/kpc^2 -> kg / m^2

  static const double kgPerM2 = 2.089e-9;

  // Msol/kpc^2 -> protons/cm^2

  static const double protonsPerCm2 = 1.249e14;


  // DENSITIES:

  // Msol/kpc^3 -> protons/cm^3

  static const double protonsPerCm3 = 4.0473e-8;

  
    static const float heliumY = 0.24;
   

  }

  
// fixme: should be in constants namespace
static const double PI = 3.1415192;

}


#endif
