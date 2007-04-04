// subset.cpp - part of SimAn Simulation Analysis Library
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
#include <typeinfo>
#include <boost/lexical_cast.hpp>

namespace siman {

  // Subset Implementation


  Subset::~Subset() {
#ifdef SIMAN_TRACE
    cerr << "Subset::~Subset" << endl;
#endif
    map<SimanArray*,SimanArraySubscripted*>::iterator i;
    for(i=createdReferencers.begin();i!=createdReferencers.end();i++) {
      delete (*i).second;
    }
  }

  Subset::Subset() {

#ifdef SIMAN_TRACE
    cerr << "Subset::Subset" << endl; 
#endif
    pAutoParent = NULL;
  } 
 
  Subset::Subset(SimSnap *simulation) {

    // start off with NO particles!
    pAutoParent = simulation;
  }

  Subset::Subset(SimSnap &sim) {
    pAutoParent = &sim;
  }

  Subset::Subset(SimSnap *simulation,const Filter &filter) {
    pAutoParent = simulation;

    for(unsigned int n=0;n<pAutoParent->getNumParticles();n++) {
      const Particle *pParticle = simulation->getConstParticle(n);

      if(filter.includes(*pParticle))
	particleRefs.push_back(n);
    

    }
  }


  Subset::Subset(const SimSnap &simulation,const Filter &filter) {

    // const_cast - not particularly nice. We don't generally abuse this trust. The
    // constructor needs to be like this for boost.python. Would be nice
    // to sort this out, but it needs quite a bit of things doing.

    pAutoParent = const_cast<SimSnap*>(&simulation);

    for(unsigned int n=0;n<pAutoParent->getNumParticles();n++) {
      const Particle *pParticle = pAutoParent->getConstParticle(n);

      if(filter.includes(*pParticle))
	particleRefs.push_back(n);
    

    }
  }

  Subset::Subset(SimSnap *simulation, int numParticles, int *pPartList) {
    pAutoParent = simulation;
    particleRefs.insert(particleRefs.end(),&pPartList[0],&pPartList[numParticles-1]);
  }

  Subset::Subset(SimSnap *simulation, const vector<unsigned int> & particles) {
    pAutoParent = simulation;
    particleRefs = particles;
  }

  SimSnap *Subset::getParent() const {
    return pAutoParent;
  }

  unsigned int Subset::getNumParticles() const {
    return particleRefs.size();
  }

  void Subset::pushParticle(int n) {
    particleRefs.push_back(n);
  }

  Particle* Subset::getParticle(unsigned int id) {
#ifndef SIMAN_UNSAFE_FASTER
    if(id>=getNumParticles())
      throw(std::out_of_range(boost::lexical_cast<string>(id)));
#endif

    // could do checking of id, but a conditional may disturb pipelining
    // so have not for now...

    Particle* particle = pAutoParent->getParticle(particleRefs[id]);
  
    return particle;

  }


  const Particle* Subset::getConstParticle(unsigned int id) const {
#ifndef SIMAN_UNSAFE_FASTER
    if(id>=getNumParticles())
      throw(std::out_of_range(boost::lexical_cast<string>(id)));
#endif
    return pAutoParent->getConstParticle(particleRefs[id]);

  }

  void Subset::releaseParticle(unsigned int id) {
    
    pAutoParent->releaseParticle(particleRefs[id]);

  }

  int Subset::deReference(int i, int n) const {

    if(n==1) 
      return particleRefs[i];
    else
      return pAutoParent->deReference(particleRefs[i], n-1);
  }


  int Subset::deReference(int i, SimSnap *pObj) const {
    if(pObj==this)
      return i;
    else
      return pAutoParent->deReference(particleRefs[i],pObj);

  }

  SimanArray & Subset::getSubscriptedArrayFor(SimanArray &parent) const {
    if(createdReferencers.count(&parent)!=0)
      return *(createdReferencers[&parent]);
    else {
      SimanArraySubscripted * subscriptedArr = 
	new SimanArraySubscripted(parent,const_cast<Subset*>(this));
      createdReferencers[&parent]=subscriptedArr;
      return *subscriptedArr;
    }
  }


  SimanArray & Subset::createArray(string name, string fullname, Unit arrunits) {
    SimanArray & parentArray = pAutoParent->createArray(name,fullname,arrunits);

    return getSubscriptedArrayFor(parentArray);

  }

  void Subset::destroyArray(string name) {
    pAutoParent->destroyArray(name);
    // currently does not destroy cached SimanArraySubscripted
  }

  SimanArray & Subset::getArray(string name) {
    return getSubscriptedArrayFor(pAutoParent->getArray(name));
  }

  SimanArray & Subset::getArray(int n) {
    return getSubscriptedArrayFor(pAutoParent->getArray(n));
  }


  SimanArray const & Subset::getConstArray(string name) const {
    return getSubscriptedArrayFor(const_cast<SimanArray &>(pAutoParent->getConstArray(name)));
  }

  SimanArray const & Subset::getConstArray(int n) const {
    return getSubscriptedArrayFor(const_cast<SimanArray &>(pAutoParent->getConstArray(n)));
  }

} // namespace siman
