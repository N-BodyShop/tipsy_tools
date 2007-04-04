// particle.cpp - part of SimAn Simulation Analysis Library
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

  Particle::Particle(float x1, float y1, float z1, float vx1, float vy1, float vz1, float mass1, unsigned int type1, float temp1, float rho1, float u1, float ne1, float metal1) {
    x = x1;
    y = y1; 
    z = z1;
    vx = vx1;
    vy = vy1;
    vz = vz1;
    mass = mass1;
    temp = temp1;
    rho = rho1;
    u = u1;
    ne = ne1;
    type = type1;
    metal = metal1;
 
  }


  Particle::Particle(float x1, float y1, float z1) :
    x(x1), y(y1), z(z1), vx(0), vy(0), vz(0), mass(0), type(0), temp(0), rho(0), u(0), ne(0), metal(0) {
    
    
  }

  Particle::Particle() :
    x(0), y(0), z(0), vx(0), vy(0), vz(0), mass(0), type(0), temp(0), rho(0), u(0), ne(0), metal(0)
  {
 
  }


  Particle::~Particle() {

  }

  void Particle::operator=(const Particle &copyfrom) {
    this->x = copyfrom.x;
    this->y = copyfrom.y;
    this->z = copyfrom.z;
    this->vx = copyfrom.vx;
    this->vy = copyfrom.vy;
    this->vz = copyfrom.vz;
    this->mass = copyfrom.mass;
    this->temp = copyfrom.temp;
    this->rho = copyfrom.rho;
    this->u = copyfrom.u;
    this->ne = copyfrom.ne;
    this->type = copyfrom.type;
  }


  bool Particle::operator==(const Particle &comp) const 
  {
    if(x!=comp.x || y!=comp.y || z!=comp.z || vx!=comp.vx || vy!=comp.vy || vz!=comp.vz ||
       mass!=comp.mass || temp!=comp.temp || rho!=comp.rho || u!=comp.u || ne!=comp.ne ||
       type!=comp.type) return false;
  
    return true;
  }

  bool Particle::operator!=(const Particle &comp) const {
    return !((*this)==comp);
  }

  void Particle::nativeWrite(ofstream *file) const {
    file->write((char*) &x, sizeof(float));
    file->write((char*) &y, sizeof(float));
    file->write((char*) &z, sizeof(float));
    file->write((char*) &vx, sizeof(float));
    file->write((char*) &vy, sizeof(float));
    file->write((char*) &vz, sizeof(float));
    file->write((char*) &mass, sizeof(float));
    file->write((char*) &type, sizeof(int));
    file->write((char*) &temp, sizeof(float));
    file->write((char*) &rho, sizeof(float));
    file->write((char*) &u, sizeof(float));
    file->write((char*) &ne, sizeof(float));
    file->write((char*) &metal, sizeof(float));
  }

  Particle::Particle(ifstream *file, int vernum) {
  
    file->read((char*) &x, sizeof(float));
    file->read((char*) &y, sizeof(float));
    file->read((char*) &z, sizeof(float));
    file->read((char*) &vx, sizeof(float));
    file->read((char*) &vy, sizeof(float));
    file->read((char*) &vz, sizeof(float));
    file->read((char*) &mass, sizeof(float));
    file->read((char*) &type, sizeof(int));
    file->read((char*) &temp, sizeof(float));
    file->read((char*) &rho, sizeof(float));
    file->read((char*) &u, sizeof(float));
    file->read((char*) &ne, sizeof(float));
    if(vernum<4) {
      // no way to save certain data in these old files without mucking up structures
      float discard;
      file->read((char*) &discard, sizeof(float));
    }
    file->read((char*) &metal, sizeof(float));
  }

} // namespace siman
