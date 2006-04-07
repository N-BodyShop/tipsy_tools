//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

CParticle::CParticle(float x1, float y1, float z1, float vx1, float vy1, float vz1, float mass1, unsigned int type1, float temp1, float rho1, float u1, float ne1, float metal1) {
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


CParticle::CParticle(float x1, float y1, float z1) {
  x = x1;
  y = y1; 
  z = z1;
 
 }

CParticle::CParticle() {
 
}


CParticle::~CParticle() {

}

void CParticle::operator=(CParticle copyfrom) {
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


bool CParticle::operator==(const CParticle &comp) const 
{
  if(x!=comp.x || y!=comp.y || z!=comp.z || vx!=comp.vx || vy!=comp.vy || vz!=comp.vz ||
     mass!=comp.mass || temp!=comp.temp || rho!=comp.rho || u!=comp.u || ne!=comp.ne ||
     type!=comp.type) return false;
  
  return true;
}

bool CParticle::operator!=(const CParticle &comp) const {
  return !((*this)==comp);
}

void CParticle::nativeWrite(ofstream *file) {
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
  file->write((char*) &nHp, sizeof(float));
  file->write((char*) &metal, sizeof(float));
}

CParticle::CParticle(ifstream *file, int vernum) {
  
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
  file->read((char*) &nHp, sizeof(float));
  file->read((char*) &metal, sizeof(float));
}
