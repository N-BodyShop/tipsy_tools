// simanarray.cpp - part of SimAn Simulation Analysis Library
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


#include "../base.hpp"
#include <boost/lexical_cast.hpp>

namespace siman {

  SimSnap * SimanArray::getOwner() {
    return pOwner;
  }

  float& SimanArray::operator[](unsigned int i) {
#ifndef SIMAN_UNSAFE_FASTER
    if(i>=len)
      throw(std::out_of_range(boost::lexical_cast<string>(i)));
#endif
    pOwner->bumpVersion();
    return array[i];
  }

  const float& SimanArray::operator[](unsigned int i) const {
#ifndef SIMAN_UNSAFE_FASTER
    if(i>=len)
      throw(std::out_of_range(boost::lexical_cast<string>(i)));
#endif
    return array[i];
  }

  void SimanArray::unset(unsigned int n) {
    (*this)[n] = std::numeric_limits<float>::quiet_NaN();
  }

  bool SimanArray::isSet(unsigned int n) {
    return((*this)[n]==(*this)[n]);
  }

  bool SimanArray::isAllSet() {
    for(unsigned int n=0; n!=size();n++) {
      if(!isSet(n))
	return false;
    }
    return true;
  }

  SimanArray::~SimanArray() {

#ifdef SIMAN_TRACE
    cerr << "~SimanArray " << this << "/" << array << endl;
#endif
    delete[] array;
  }

  unsigned int SimanArray::size() const {
    return len;
  }

  SimanArray::SimanArray(unsigned int leni, SimSnap *pOwner, const Unit & unitsi) : units(unitsi), pOwner(pOwner) {
    len=leni;
    array = new float[len];
    for(unsigned int n=0; n<len; n++)
      unset(n);
  }

  SimanArray::SimanArray() {
    len=0;
    array = NULL;
    pOwner=NULL;
  }

  void SimanArray::read(istream &fin, int file_ver_num) {
    fin.read((char*) array, sizeof(float)*len);
  }

  void SimanArray::readTipsy(string fname) {
    ifstream fin(fname.c_str());
    unsigned int n;
    fin >> n;
    if(n!=size())
      throw SimanException("Mismatched size");
    
    for(n=0;n<size();n++) {
      fin >> (*this)[n];
    }
  }

  void SimanArray::write(ostream &file) const {
    file.write((char*) array, sizeof(float)*len);
  }

  void SimanArray::operator=(const SimanArray &copyfrom) {
    //  if(copyfrom.pOwner!=pOwner)
    //  throw(SimanException("Unable to change owner of an existing array"));
    if(copyfrom.size()!=size())
      throw MismatchedArrayLength();
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]=copyfrom[n];
    }
    units=copyfrom.getUnits();
  }

  
  SimanArray::SimanArray (const SimanArray &obj) {
    len=obj.size();
    array = new float[len];
    for(unsigned int n=0; n<len; n++)
      array[n]=obj[n];
    units=obj.getUnits();
    pOwner=obj.pOwner;
  }
      

  Unit SimanArray::getUnits() const {
    return units;
  }
  
  void SimanArray::convertUnits(const Unit & to) throw(UnitsError) {
    double ratio = getUnits().convertTo(to,pOwner);
    for(unsigned int n=0; n<size(); n++)
      (*this)[n]*=ratio;
    units=to;
  }

  void SimanArray::operator*=(const SimanArray &o) throw(MismatchedArrayLength) {
    if(o.size()!=size())
      throw MismatchedArrayLength();
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]*=o[n];
    }
    units*=o.getUnits();
  }


  void SimanArray::operator/=(const SimanArray &o) throw(MismatchedArrayLength) {
    if(o.size()!=size())
      throw MismatchedArrayLength();
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]/=o[n];
    }
    units/=o.getUnits();
  }


  void SimanArray::operator+=(const SimanArray &o) throw(MismatchedArrayLength, UnitsError) {
    if(o.size()!=size())
      throw MismatchedArrayLength();

    double convrat = o.getUnits().convertTo(getUnits()); // will throw UnitsError if no good

    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]+=o[n]*convrat;
    }

  }


  void SimanArray::operator-=(const SimanArray &o) throw(MismatchedArrayLength, UnitsError) {
    if(o.size()!=size())
      throw MismatchedArrayLength();

    double convrat = o.getUnits().convertTo(getUnits()); // will throw UnitsError if no good

    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]-=o[n]*convrat;
    }

  }

  void SimanArray::operator+=(float o) {
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]+=o;
    }
  }


  void SimanArray::operator-=(float o) {
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]-=o;
    }
  }


  void SimanArray::operator*=(float o) {
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]*=o;
    }
  }


  void SimanArray::operator/=(float o) {
    for(unsigned int n=0; n<size(); n++) {
      (*this)[n]/=o;
    }
  }

  SimanArray SimanArray::operator+(const SimanArray &o) const throw(MismatchedArrayLength, UnitsError)  {
    SimanArray n(*this);
    n+=o;
    return n;
  }


  SimanArray SimanArray::operator-(const SimanArray &o) const throw(MismatchedArrayLength, UnitsError)  {
    SimanArray n(*this);
    n-=o;
    return n;
  }


  SimanArray SimanArray::operator*(const SimanArray &o) const throw(MismatchedArrayLength) {
    SimanArray n(*this);
    n*=o;
    return n;
  }


  SimanArray SimanArray::operator/(const SimanArray &o) const throw(MismatchedArrayLength) {
    SimanArray n(*this);
    n/=o;
    return n;
  }


  SimanArray SimanArray::operator+(float f) const {
    SimanArray n(*this);
    n+=f;
    return n;
  }


  SimanArray SimanArray::operator-(float f) const {
    SimanArray n(*this);
    n-=f;
    return n;
  }



  SimanArray SimanArray::operator*(float f) const {
    SimanArray n(*this);
    n*=f;
    return n;
  }


  SimanArray SimanArray::operator/(float f) const {
    SimanArray n(*this);
    n/=f;
    return n;
  }

  SimanArray SimanArray::operator*(const Unit &u) const {
    SimanArray n(*this);
    n.units*=u;
    return n;
  }


  SimanArray SimanArray::operator/(const Unit &u) const {
    SimanArray n(*this);
    n.units/=u;
    return n;
  }

  SimanArray SimanArray::power(const Rational &p) const {
    SimanArray ret(*this);
    for(unsigned int n=0; n<size(); n++) {
      ret[n]=std::pow(ret[n],p.flt());
    }
    ret.units = siman::pow(ret.units,p);
    return ret;
  }

  float SimanArray::max() {
    float max = numeric_limits<float>::min();
    for(unsigned int n=0; n<size(); n++) {
      if(max<(*this)[n] && isSet(n))
	max = (*this)[n];
    }
    return max;
  }


  float SimanArray::min() {
    float min = numeric_limits<float>::max();
    for(unsigned int n=0; n<size(); n++) {
      if(min>(*this)[n] && isSet(n))
	min = (*this)[n];
    }
    return min;
  }

  float SimanArray::mean() {
    double mean=0.;
    double ratio = (double)1./(double)size();
    for(unsigned int n=0; n<size(); n++) {
      if(isSet(n))
	mean+=((double)((*this)[n]))*ratio;
    }
    return (float) mean;
  }


  float SimanArray::total() {
    float tot=0.;
    for(unsigned int n=0; n<size(); n++) {
      if(isSet(n)) {
	tot+=(*this)[n];
	if(tot!=tot) {
	  static bool ferr = true;
	  if(ferr) cerr << n << " " << (*this)[n] << " " << (*this)[n-1] << " " << tot;
	}
      }
    }
   
    return tot;
  }
  

  /// SUBSCRIPTED OVERRIDES:


  float& SimanArraySubscripted::operator[](unsigned int i) {

    return (*origArray)[sub->deReference(i,origArray->getOwner())];
  }
  
   const float& SimanArraySubscripted::operator[](unsigned int i) const {

     return ((const SimanArray &)(*origArray))[sub->deReference(i,origArray->getOwner())];
  }

   
  SimanArraySubscripted::SimanArraySubscripted(SimanArray &origArray_in, Subset *pOwneri)
  {
    origArray = &origArray_in;

    sub=pOwneri;
    pOwner=pOwneri;
    
  }



   void SimanArraySubscripted::read(istream &fin, int file_ver_num) {
    throw SimanException("Read into subscripted array not currently supported.");
  }

   void SimanArraySubscripted::write(ostream &file) const {
    
     for(unsigned int i=0; i<size(); i++) {
      file.write((char*) &((*this)[i]), sizeof(float));
    }
  }

   unsigned int SimanArraySubscripted::size() const {
     return sub->getNumParticles();
  }

  void SimanArraySubscripted::convertUnits(const Unit & to) throw(UnitsError) {
    origArray->convertUnits(to);
  }

  Unit SimanArraySubscripted::getUnits() const {
    return origArray->getUnits();
  }

  // VIRTUAL overrides

  SimanArrayVirtual::SimanArrayVirtual(SimSnap *pOwneri, float Particle::* member) : member(member) { 
    pOwner=pOwneri;
    
  }

  float & SimanArrayVirtual::operator[](unsigned int i) {
    Particle *p=pOwner->getParticle(i);
    return p->*(member);
  }

  const float & SimanArrayVirtual::operator[](unsigned int i) const {
    const Particle *p = pOwner->getConstParticle(i);
    return p->*(member);
  }

  Unit SimanArrayVirtual::getUnits() const {
    if(member == &Particle::x || member == &Particle::y || member == &Particle::z) { 
      return pOwner->getDistanceUnits();
    } else if(member == &Particle::vx || member == &Particle::vy || member == &Particle::vz) {
      return pOwner->getVelocityUnits();
    } else if(member == &Particle::mass) {
      return pOwner->getMassUnits();
    } else if(member == &Particle::rho) {
      return pOwner->getDensityUnits();
    } else if(member== &Particle::u) {
      return pOwner->getEnergyUnits();
    } else if(member==&Particle::ne) {
      return Unit();
    }
    return Unit();
  }


  void SimanArrayVirtual::convertUnits(const Unit &to) throw(UnitsError) {
    if(member == &Particle::x || member == &Particle::y || member == &Particle::z) {
      
      pOwner->convertUnits(to,pOwner->getMassUnits(),pOwner->getVelocityUnits(),pOwner->getDensityUnits(),pOwner->getEnergyUnits());
    } else if(member == &Particle::vx || member == &Particle::vy || member == &Particle::vz) {
   
      pOwner->convertUnits(pOwner->getDistanceUnits(),pOwner->getMassUnits(),to,pOwner->getDensityUnits(),pOwner->getEnergyUnits());
    } else if(member == &Particle::mass) {

      pOwner->convertUnits(pOwner->getDistanceUnits(),to,pOwner->getVelocityUnits(),pOwner->getDensityUnits(),pOwner->getEnergyUnits());  
    } else if(member == &Particle::rho) {
      pOwner->convertUnits(pOwner->getDistanceUnits(),pOwner->getMassUnits(),pOwner->getVelocityUnits(),to,pOwner->getEnergyUnits());
    } else if(member== &Particle::u) {
      pOwner->convertUnits(pOwner->getDistanceUnits(),pOwner->getMassUnits(),pOwner->getVelocityUnits(),pOwner->getDensityUnits(),to);
    } else if(member==&Particle::ne) {
      // No in-built provision for this, we'll need to do it manually
      // (even if it's pretty nonsensical, someone might decide to do a
      // multiplication this way...)
      SimanArray::convertUnits(to);
      
    }
  }

  unsigned int SimanArrayVirtual::size() const {
    return pOwner->getNumParticles();
  }

  
}
