// particle.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __PARTICLE_H_INCLUDED

#define __PARTICLE_H_INCLUDED

namespace siman {

class Particle {

public:

  // GENERIC PROPERTIES

  float x,y,z;
  float vx,vy,vz;
  float mass;
  
  unsigned char type;

  // GAS PROPERTIES

  float temp, rho, u, ne, metal;

  // CONSTRUCTORS / DESTRUCTORS

  Particle();
  Particle(float x1, float y1, float z1);
  Particle(float x1, float y1, float z1, float vx1, float vy1, float vz1, float mass1, unsigned int type1, float temp1=0, float rho1=0, float u1=0, float ne1=0, float metal=0);

  /// construct particle from data in file at current pointer
  /// @param vernum - version number of the file format being read
  Particle(std::ifstream *file, int vernum);
  ~Particle();

  // HANDY THINGS

  inline float distanceTo(const Particle &them) const  {
    float d2 = (them.x-x) * (them.x-x) + (them.y-y) * (them.y-y) + (them.z-z)*(them.z-z);
    return sqrt(d2);
  }



  inline float squaredDistanceTo(const Particle &them) const {
    float d2 = (them.x-x) * (them.x-x) + (them.y-y) * (them.y-y) + (them.z-z)*(them.z-z);
    return d2;
  }



  // OPERATOR METHODS

  void operator=(const Particle &copyfrom);
  bool operator==(const Particle &comp) const;
  bool operator!=(const Particle &comp) const;

  // FILE HANDLING

  /// dump particle to file
  void nativeWrite(std::ofstream *file) const;

  // TYPES
  // - no longer follows gadget convention

  static const unsigned char gas = 1;
  static const unsigned char dm = 2;
  static const unsigned char star = 4;
  static const unsigned char anyType = 7;


};

}

#endif
