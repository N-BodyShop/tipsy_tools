// union.cpp - part of SimAn Simulation Analysis Library
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
#include <boost/lexical_cast.hpp>

namespace siman {


  Union::Union(SimSnap *parent) {

    // N.B. the notional "parent" of this union will be used for functions such
    // as getRedshift, etc...
    //
    // There is no strict need for unified particles to come from the same
    // simulation, however...

    pAutoParent = parent;

  }

  unsigned int Union::getNumParticles() const {
    unsigned int n = 0;
 
    for(list<SimSnap*>::const_iterator itr = children.begin(); itr!= children.end(); itr++) {
      if(*itr!=NULL) {
      
	n+=(*itr)->getNumParticles();
      }
    }

    return n;
  
  }

  Particle * Union::getParticle(unsigned int id) {
 
    for(list<SimSnap*>::iterator itr = children.begin(); itr!= children.end(); itr++) {
      int nPart = (*itr)->getNumParticles();

      if(nPart>id) 
	return (*itr)->getParticle(id);

      id-=nPart;
    }

    // fell off the end of the list

    assert(false);

    return NULL;
  }

  const Particle * Union::getConstParticle(unsigned int id) const {
    for(list<SimSnap*>::const_iterator itr = children.begin(); itr!= children.end(); itr++) {
      int nPart = (*itr)->getNumParticles();

      if(nPart>id) 
	return (*itr)->getConstParticle(id);

      id-=nPart;
    }

    // fell off the end of the list

    assert(false);

    return NULL;
 
  }

  pair<SimSnap*,unsigned int> Union::getSnapAndRef(unsigned int id) const {
    unsigned int idx=id;
    for(list<SimSnap*>::const_iterator itr = children.begin(); itr!= children.end(); itr++) {
      int nPart = (*itr)->getNumParticles();
      
      if(nPart>id) 
	return pair<SimSnap*,unsigned int>((*itr),id);
      
      id-=nPart;
    }
    
    throw(std::out_of_range(boost::lexical_cast<string>(idx)));

  }

  void Union::add(SimSnap *add) {
    //assert(add->getParent() == pSimulation);
    children.push_front(add);
  }

  void Union::add(SimSnap &add) {
    children.push_front(&add);
  }

  list<SimSnap*> & Union::getList() {
    return children;
  }

  
  int Union::deReference(int i, int n) const {
    
    if(n==0)
      return i;
    else {
      pair<SimSnap*,unsigned int> p=getSnapAndRef(i);
      return (p.first)->deReference(p.second,n-1);
    }

  }

  int Union::deReference(int i, SimSnap *pObj) const {
    if(pObj==this)
      return i;
    else {
      pair<SimSnap*,unsigned int> p=getSnapAndRef(i);
      return (p.first)->deReference(p.second,pObj);
    }
  }


  Union::~Union() {

  }

} // namespace siman
