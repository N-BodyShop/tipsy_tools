//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

using namespace std;

CUnion::CUnion(CSimSnap *parent) {

  // N.B. the notional "parent" of this union will be used for functions such
  // as getRedshift, etc...
  //
  // There is no strict need for unified particles to come from the same
  // simulation, however...

  pAutoParent = parent;

}

int CUnion::getNumParticles() {
  int n = 0;
 
  for(list<CSimSnap*>::iterator itr = children.begin(); itr!= children.end(); itr++) {
    if(*itr!=NULL) {
      
      n+=(*itr)->getNumParticles();
    }
  }

  return n;
  
}

CParticle * CUnion::getParticle(int id) {
 
  for(list<CSimSnap*>::iterator itr = children.begin(); itr!= children.end(); itr++) {
    int nPart = (*itr)->getNumParticles();

    if(nPart>id) 
      return (*itr)->getParticle(id);

    id-=nPart;
  }

  // fell off the end of the list

  assert(false);

  return NULL;
}

void CUnion::releaseParticle(CParticle *p) {
  pAutoParent->releaseParticle(p);
}

void CUnion::pushParticle(int n) {
  assert(false);
  // invalid operation

}

void CUnion::add(CSimSnap *add) {
  //assert(add->getParent() == pSimulation);
  children.push_front(add);
}

CUnion::~CUnion() {

}
