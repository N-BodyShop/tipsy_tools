// filter.cpp - part of SimAn Simulation Analysis Library
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


  // Filter

  bool Filter::includes(const Particle &particle) const {
    return true;
  }
  
  Filter* Filter::copy() const {
    return new Filter();
  }

  NotFilter Filter::operator!() const {
    return NotFilter(*this);
  }

  AndFilter Filter::operator&(const Filter &f2) const {
    return AndFilter(*this, f2);
  }

  OrFilter Filter::operator|(const Filter &f2) const {
    return OrFilter(*this, f2);
  }


  // PIPE filter class

  Column::Column(float x1i, float y1i, float x2i, float y2i)
  {
    x1 = (x1i<x2i)?x1i:x2i;
    x2 = (x1i<x2i)?x2i:x1i;
  
    y1 = (y1i<y2i)?y1i:y2i;
    y2 = (y1i<y2i)?y2i:y1i;
  
  }

  bool Column::includes(const Particle &particle) const {
    // is particle within bounds?
    return (particle.x>x1 && particle.y>y1 && particle.y<y2 && particle.x<x2);
  }

  Filter* Column::copy() const {
    return new Column(x1,y1,x2,y2);
  }



  // SPHERE filter class

  Sphere::Sphere(float xci, float yci, float zci, float ri) : 
    xc(xci), yc(yci), zc(zci), r(ri) {

  }

  Sphere::Sphere(float ri) :
    xc(0), yc(0), zc(0), r(ri) {

  }

  bool Sphere::includes(const Particle &p) const {

    float distance = sqrt((p.x-xc)*(p.x-xc) + 
			  (p.y-yc)*(p.y-yc) +
			  (p.z-zc)*(p.z-zc));
  
    return (distance<r);

  }

  Filter* Sphere::copy() const {
    return new Sphere(xc,yc,zc,r);
  }


  VelocitySphereFilter::VelocitySphereFilter(float xci, float yci, float zci, float ri) : 
    xc(xci), yc(yci), zc(zci), r(ri) {

  }

  VelocitySphereFilter::VelocitySphereFilter(float ri) :
    xc(0), yc(0), zc(0), r(ri) {

  }

  bool VelocitySphereFilter::includes(const Particle &p) const {

    float distance = sqrt((p.vx-xc)*(p.vx-xc) + 
			  (p.vy-yc)*(p.vy-yc) +
			  (p.vz-zc)*(p.vz-zc));
  
    return (distance<r);

  }

  Filter* VelocitySphereFilter::copy() const {
    return new VelocitySphereFilter(xc,yc,zc,r);
  }

  // PARTICLE TYPE filter class

  ParticleTypeFilter::ParticleTypeFilter(int typei) {
    type = typei;
  }

  bool ParticleTypeFilter::includes( const Particle &p) const {
    return ((p.type & type) > 0);
  }

  Filter* ParticleTypeFilter::copy() const {
    return new ParticleTypeFilter(type);
  }

  // DENSITY CUT class

  DensityCutFilter::DensityCutFilter(float cutAti) {
    cutAt = cutAti;
  }

  bool DensityCutFilter::includes(const Particle &p) const {
    return(p.rho > cutAt);
  }

  Filter* DensityCutFilter::copy() const {
    return new DensityCutFilter(cutAt);
  }

  // RANDOM filter

  RandomFilter::RandomFilter(float probi) {
    prob = probi;
  }

  bool RandomFilter::includes( const Particle &p) const {
    float s = (float)rand()/(float)RAND_MAX;
    if(s<prob) return true;
    return false;
  }

  Filter* RandomFilter::copy() const {
    return new RandomFilter(prob);
  }

  // MODULO filter

  ModuloFilter::ModuloFilter(int num, int offset) {
    cur=offset;
    mod=num;
  }

  bool ModuloFilter::includes( const Particle &p) const {
    bool inc = (cur==0);
    cur++;
    if(cur>=mod) cur=0;
    return inc;
  }

  Filter* ModuloFilter::copy() const {
    return new ModuloFilter(cur,mod);
  }
  // NOT filter (modifier)

  NotFilter::NotFilter(const Filter &fi) : f(fi.copy()) {
    
  }

  NotFilter::NotFilter(const boost::shared_ptr<const Filter> &fi) : f(fi) {

  }

  bool NotFilter::includes( const Particle &p) const {
    return !(f->includes(p));
  }

  Filter* NotFilter::copy() const {
    return new NotFilter(f);
  }


  AndFilter::AndFilter(const Filter &f1i, const Filter &f2i) : f1(f1i.copy()), f2(f2i.copy()) {
    
  }
  
  AndFilter::AndFilter(const boost::shared_ptr<const Filter> &f1i, const boost::shared_ptr<const Filter> &f2i) : f1(f1i), f2(f2i) {

  }

  Filter* AndFilter::copy() const {
    return new AndFilter(f1,f2);
  }


  bool AndFilter::includes(const Particle &p) const {
    return f1->includes(p) && f2->includes(p);
  }


  OrFilter::OrFilter(const Filter &f1i, const Filter &f2i) : f1(f1i.copy()), f2(f2i.copy()) {

  }

    
  OrFilter::OrFilter(const boost::shared_ptr<const Filter> &f1i, const boost::shared_ptr<const Filter> &f2i) : f1(f1i), f2(f2i) {

  }

  Filter* OrFilter::copy() const {
    return new OrFilter(f1,f2);
  }


  bool OrFilter::includes( const Particle &p) const {
    return f1->includes(p) || f2->includes(p);
  }


} // namespace siman
