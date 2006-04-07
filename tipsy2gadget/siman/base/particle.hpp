//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//

#ifndef __PARTICLE_H_INCLUDED

#define __PARTICLE_H_INCLUDED

class CParticle {

public:

  // GENERIC PROPERTIES

  float x,y,z;
  float vx,vy,vz;
  float mass;
  
  unsigned int type;

  // GAS PROPERTIES

  float temp, rho, u, ne, metal, nHp;

  // CONSTRUCTORS / DESTRUCTORS

  CParticle();
  CParticle(float x1, float y1, float z1);
  CParticle(float x1, float y1, float z1, float vx1, float vy1, float vz1, float mass1, unsigned int type1, float temp1=0, float rho1=0, float u1=0, float ne1=0, float metal=0);

  /// construct particle from data in file at current pointer
  /// @param vernum - version number of the file format being read
  CParticle(std::ifstream *file, int vernum);
  ~CParticle();

  // HANDY THINGS

  inline float distanceTo(const CParticle &them)  {
    float d2 = (them.x-x) * (them.x-x) + (them.y-y) * (them.y-y) + (them.z-z)*(them.z-z);
    return sqrt(d2);
  }



  inline float squaredDistanceTo(const CParticle &them) {
    float d2 = (them.x-x) * (them.x-x) + (them.y-y) * (them.y-y) + (them.z-z)*(them.z-z);
    return d2;
  }



  // OPERATOR METHODS

  void operator=(CParticle copyfrom);
  bool operator==(const CParticle &comp) const;
  bool operator!=(const CParticle &comp) const;

  // FILE HANDLING

  /// dump particle to file
  void nativeWrite(std::ofstream *file);

  // TYPES
  // - no longer follows gadget convention

  static const int gas = 1;
  static const int dm = 2;
  static const int star = 4;
  static const int anyType = 7;


};

#endif
