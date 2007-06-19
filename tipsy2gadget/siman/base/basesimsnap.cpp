//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CBaseSimSnap implementation

#include "siman.hpp"

using namespace std;


CBaseSimSnap::CBaseSimSnap(string filename) {
  
  cerr << "CBaseSimSnap: reading " << filename << "...";
  ifstream file(filename.c_str(),ios::binary);

  char id[4];

  file.read(id, 4);
  
  if(id[0]!='S' || id[1]!='I' || id[2]!='B' || id[3]!='I')
    return; // should throw something here


  int file_ver_num = 1;

  file.read((char*) &file_ver_num, sizeof(int)); // version number of filetype

  file.read((char*) &numParticles,sizeof(int));
  cerr << endl << numParticles << endl;

  lenUnits = units::CUnit(&file,file_ver_num);
  velUnits = units::CUnit(&file,file_ver_num);
  denUnits = units::CUnit(&file,file_ver_num);
  enUnits = units::CUnit(&file,file_ver_num);
  massUnits = units::CUnit(&file,file_ver_num);
  
  file.read((char*) &(boxSize), sizeof(float));
  file.read((char*) &(hubble),sizeof(float));
  file.read((char*) &(redshift), sizeof(float));
  file.read((char*) &(om_m0), sizeof(float));
  file.read((char*) &(om_lam0),sizeof(float));

  pParticles = (CParticle **) malloc(sizeof(void*)*numParticles);

  for(unsigned int n=0; n<numParticles; n++) {
    pParticles[n] = new CParticle(&file,file_ver_num);
  }

  if(file_ver_num>1) {
    char *temp_chararr= NULL;
    int nArr;
    file.read((char*) &nArr, sizeof(int));
    cerr << "done main section!" << endl << "CBaseSimSnap: " << nArr << " extra arrays to read" << endl;
    for(int n=0; n<nArr; n++) {
      int shortname_len;
      int longname_len;
      string shortname;
      string longname;
      file.read((char*) &shortname_len, sizeof(int));
      temp_chararr = new char[shortname_len+1];
      file.read(temp_chararr,shortname_len+1);
      shortname = temp_chararr;
      delete[] temp_chararr;

      file.read((char*) &longname_len, sizeof(int));
      temp_chararr = new char[longname_len+1];
      file.read(temp_chararr,longname_len+1);
      longname = temp_chararr;
      delete[] temp_chararr;

      cerr << "CBaseSimSnap: array " << n << " - " << longname << "(" << shortname << ")...";
      float *newArr = createArray(shortname,longname,units::CUnit(&file, file_ver_num));
      file.read((char*) newArr, sizeof(float)*numParticles);
      cerr << "done!" << endl;
    }
  } else {
    cerr << "done!" << endl;
  }

  cerr << lenUnits << " | " << velUnits << " | " << denUnits << " | " << enUnits << " | " << massUnits << endl;

}

CBaseSimSnap::CBaseSimSnap(CSimSnap *copyFrom) {
  numParticles = copyFrom->getNumParticles();

  pParticles = (CParticle **) malloc(sizeof(void*)*numParticles);
  
  for(unsigned int n=0; n<numParticles; n++) {
    CParticle *toCopy = copyFrom->getParticle(n);
    pParticles[n] = new CParticle(*toCopy);
    copyFrom->releaseParticle(toCopy);
  }

  lenUnits = copyFrom->getDistanceUnits();
  velUnits = copyFrom->getVelocityUnits();
  boxSize = copyFrom->getBoxSize();
  hubble = copyFrom->getHubble();
  redshift = copyFrom->getRedshift();


}

CBaseSimSnap::CBaseSimSnap(CParticle **pParticlesi, int numParticlesi, units::CUnit len, units::CUnit mass) {

  
  pParticles = pParticlesi;
  numParticles = numParticlesi;

  
  if(len==units::CUnit()) 
    lenUnits = units::CUnit(units::len_kpc);
  else
    lenUnits = len;

  if(mass==units::CUnit()) 
    massUnits = units::CUnit(units::mass_Msol);
  else
    massUnits = mass;


  velUnits = units::CUnit(units::vel_kmPerS);

  denUnits = massUnits/pow(lenUnits,3);
  enUnits = pow(velUnits,2);

 
}

CBaseSimSnap::CBaseSimSnap() {

  pParticles=NULL;
  numParticles = 0;

  lenUnits = units::CUnit(units::len_kpc);
  massUnits = units::CUnit(units::mass_Msol);
  velUnits = units::CUnit(units::vel_kmPerS);
  denUnits = massUnits/pow(lenUnits,3);
  enUnits = pow(velUnits,2);


}


CBaseSimSnap* CBaseSimSnap::makeSlab(unsigned int nPart, float depth, float extent, float density, float temp, int genType, units::CUnit length, units::CUnit mass) {

  // make a uniform density slab out of nPart particles

  CParticle **pParticles;

  if(genType==genRandom) {
    
    pParticles = (CParticle **) malloc(sizeof(void*)*nPart);
    float partMass = density*extent*extent*depth/((float)nPart);
    for(unsigned int n=0; n<nPart; n++) {
      pParticles[n] = new CParticle(((float)rand()/(float)RAND_MAX)*extent-extent/2,((float)rand()/(float)RAND_MAX)*extent-extent/2,((float)rand()/(float)RAND_MAX)*depth-depth/2,0,0,0,partMass,CParticle::gas,temp,density);
    }
  } else {
    float c = pow((double)(nPart/(depth*extent*extent)),1./3.);
    int nz= (int)(c*depth);
    int nx= (int)(c*extent);
    int ny = nx;
    nPart = nz*nx*ny;
    float partMass = density*extent*extent*depth/((float)nPart);
    pParticles = (CParticle **) malloc(sizeof(void*)*nPart);
    int i=0;
    for(int x=0;x<nx;x++) {
      for(int y=0;y<ny;y++) {
	for(int z=0; z<nz;z++) {
	  pParticles[i] = new CParticle(x/c-extent/2.,y/c-extent/2.,z/c-depth/2.,0,0,0,partMass,CParticle::gas,temp,density);
	  i++;
	}
      }
    }
    cerr << "CBaseSimSnap: generated regular slab with " << nPart << " particles (" << nx << " x " << ny << " x " << nz << ")" << endl;
  }
  return new CBaseSimSnap(pParticles,nPart,length,mass);
  
}

void CBaseSimSnap::allocateMemory() {
  cerr << "CBaseSimSnap: allocating memory...";
  pParticles = (CParticle **) malloc(sizeof(void*)*numParticles);
  for(unsigned int n=0; n<numParticles; n++) {
    pParticles[n] = new CParticle();
  }
  cerr << "done!" << endl;
}



void CBaseSimSnap::convertUnits(units::CUnit distance, units::CUnit mass, units::CUnit velocity, units::CUnit density, units::CUnit energy) {

  using namespace units;

  if(distance==CUnit()) distance=CUnit(len_kpc);
  if(mass==CUnit()) mass=CUnit(mass_Msol);
  if(velocity==CUnit()) velocity=CUnit(vel_kmPerS);
  if(density==CUnit()) density=mass/(distance*distance*distance);
  if(energy==CUnit()) energy=velocity*velocity;
  
  cerr << "CBaseSimSnap: convertUnits - dist: "<<distance<<" mass: "<<mass<<" vel: "<<velocity<<" den: "<<density<<" ener: "<<energy<<endl;

  double dist_ratio = lenUnits.convertTo(distance,this);
  double mass_ratio = massUnits.convertTo(mass,this);
  double vel_ratio = velUnits.convertTo(velocity,this);
  double den_ratio = denUnits.convertTo(density,this);
  double en_ratio = enUnits.convertTo(energy,this);

  cerr << "CBaseSimSnap: conversion ratios respectively: " << dist_ratio << " " << mass_ratio << " " << vel_ratio << " " << den_ratio << " " << en_ratio << endl;
  
  cerr << "CBaseSimSnap: performing conversion...";
  
  boxSize*=dist_ratio;

  for(int i=numParticles-1;i>=0;--i) {
    CParticle *p = pParticles[i];
    p->x*=dist_ratio;
    p->y*=dist_ratio;
    p->z*=dist_ratio;
    p->vx*=vel_ratio;
    p->vy*=vel_ratio;
    p->vz*=vel_ratio;
    p->rho*=den_ratio;
    p->u*=en_ratio;
    p->mass*=mass_ratio;
  }

  lenUnits = distance;
  massUnits = mass;
  velUnits = velocity;
  denUnits = density;
  enUnits = energy;

  cerr << "done!" << endl;

}

CBaseSimSnap::~CBaseSimSnap() {
  for(unsigned int n=0; n<numParticles; n++) {
    delete pParticles[n];
  }
  free(pParticles);
}

bool CBaseSimSnap::Load(bool fakeload) {
  return true;
}

CParticle* CBaseSimSnap::getParticle(int n) {
  return pParticles[n];
}

void CBaseSimSnap::releaseParticle(CParticle *) {
  // Do Nothing...
}

int CBaseSimSnap::getNumParticles() {
  return numParticles;
}

bool CBaseSimSnap::isLoaded() {
  return true;
}

float CBaseSimSnap::getBoxSize() {
  return boxSize;
}

float CBaseSimSnap::getHubble() {
  return hubble;
}

void CBaseSimSnap::setHubble(float h) {
  hubble=h;
}


void CBaseSimSnap::setRedshift(float z) {
  redshift= z;
}


void CBaseSimSnap::setBoxSize(float h) {
  boxSize=h;
}

float CBaseSimSnap::getRedshift() {
  return redshift;
}

float CBaseSimSnap::getOmegaM0() {
  return om_m0;
}

float CBaseSimSnap::getOmegaLambda0() {
  return om_lam0;
}

int CBaseSimSnap::deReference(int i, int n) {
  // if(n!=1)
  //  cerr << "CBaseSimSnap: Warning - tried to deReference beyond base level" << endl;
  return i;
}



void CBaseSimSnap::nativeWrite(CSimSnap *sim, string filename) {

  cerr << "CBaseSimSnap: writing " << filename << "...";
  ofstream file(filename.c_str(),ios::binary);

  char id[4]={'S','I','B','I'}; // identifier "SIman BInary" for autoloader

  file.write(id, 4);
  
  int file_ver_num = 3;

  file.write((char*) &file_ver_num, sizeof(int)); // version number of filetype

  unsigned int nPart = sim->getNumParticles();
  
  cerr << "(" << nPart << " particles)..." << endl;
  file.write((char*) &(nPart),sizeof(int));

  sim->getDistanceUnits().nativeWrite(&file);
  sim->getVelocityUnits().nativeWrite(&file);
  sim->getDensityUnits().nativeWrite(&file);
  sim->getEnergyUnits().nativeWrite(&file);
  sim->getMassUnits().nativeWrite(&file);
  
  float w = sim->getBoxSize();
  file.write((char*) &(w), sizeof(float));
  w=sim->getHubble();
  file.write((char*) &(w), sizeof(float));
  w=sim->getRedshift();
  file.write((char*) &(w), sizeof(float));
  w=sim->getOmegaM0();
  file.write((char*) &(w), sizeof(float));
  w=sim->getOmegaLambda0();
  file.write((char*) &(w), sizeof(float));
  

  for(unsigned int n=0; n<nPart; n++) {
    sim->getParticle(n)->nativeWrite(&file);
  }

  int numa = sim->getNumArrays();
  file.write((char*) &numa, sizeof(int));

  if(sim->getNumArrays()>0) {
    cerr << "done main section!" << endl;
    for(int n=0;n<sim->getNumArrays();n++) {
      string shortname = sim->getArrayName(n);
      int shortname_len = shortname.size();
      string longname = sim->getArrayLongName(n);
      int longname_len = longname.size();
      cerr << "CBaseSimSnap: writing " << longname << " (" << shortname << ")...";
      file.write((char*) &shortname_len, sizeof(int));
      file.write(shortname.c_str(),shortname.size()+1);
      file.write((char*) &longname_len, sizeof(int));
      file.write(longname.c_str(),longname.size()+1);
      sim->getArrayUnits(n).nativeWrite(&file);
      file.write((char*) sim->getArray(n), sizeof(float)*sim->getNumParticles());
      cerr << "done!" << endl;
    }
  } else {
    cerr << "done!" << endl;
  }
}


float * CBaseSimSnap::createArray(string name, string fullname, units::CUnit arrunits) {
  if(extraMap.count(name)>0) {
    cerr << "CBaseSimSnap::createArray: warning - array " << name << " already exists; returning existing array" << endl;
    return extraMap[name];
  } else {
    float *newArray = new float[getNumParticles()];
    extraMap[name] = newArray;
    extraMapName[name] = fullname;
    extraMapUnits[name] = arrunits;
    return newArray;
  }
}

void CBaseSimSnap::destroyArray(string name) {
  float *deleting = extraMap[name];
  delete[] deleting;
  extraMap.erase(name);
  extraMapName.erase(name);
}

units::CUnit CBaseSimSnap::getArrayUnits(string name) {
  return extraMapUnits[name];
}

units::CUnit CBaseSimSnap::getArrayUnits(int n) {
  map<string, units::CUnit>::iterator i;
  i=extraMapUnits.begin();
  for(int c=n-1;c>=0;c--)
    i++;
  return (*i).second;
}

string CBaseSimSnap::getArrayLongName(string shortname) {
  return extraMapName[shortname];
}

string CBaseSimSnap::getArrayLongName(int n) {
  map<string, string>::iterator i;
  i=extraMapName.begin();
  for(int c=n-1;c>=0;c--)
    i++;
  return (*i).second;
}


string CBaseSimSnap::getArrayName(int n) {
  map<string, float*>::iterator i;
  i=extraMap.begin();
  for(int c=n-1;c>=0;c--)
    i++;
  return (*i).first;
}


float * CBaseSimSnap::getArray(string name) {
    map<string, float*>::iterator i = extraMap.find(name);
    if(i != extraMap.end())
	return i->second;
    else
	return NULL;
}

float * CBaseSimSnap::getArray(int n) {
  map<string, float*>::iterator i;
  i=extraMap.begin();
  for(int c=n-1;c>=0;c--)
    i++;
  return (*i).second;
}

int CBaseSimSnap::getNumArrays() {
  return extraMap.size();
}
