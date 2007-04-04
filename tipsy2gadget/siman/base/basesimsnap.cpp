// basesimsnap.cpp - part of SimAn Simulation Analysis Library
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

  BaseSimSnap::BaseSimSnap(string filename) {
  
    if(getVerbose()>=1)
      cerr << "BaseSimSnap: reading " << filename << "...";
    ifstream file(filename.c_str(),ios::binary);

    char id[4];

    file.read(id, 4);
  
    if(id[0]!='S' || id[1]!='I' || id[2]!='B' || id[3]!='I')
      return; // should throw something here


    int file_ver_num = 1;

    file.read((char*) &file_ver_num, sizeof(int)); // version number of filetype

    file.read((char*) &numParticles,sizeof(int));

    lenUnits = Unit(&file,file_ver_num);
    velUnits = Unit(&file,file_ver_num);
    denUnits = Unit(&file,file_ver_num);
    enUnits = Unit(&file,file_ver_num);
    massUnits = Unit(&file,file_ver_num);
  
    file.read((char*) &(boxSize), sizeof(float));
    file.read((char*) &(hubble),sizeof(float));
    file.read((char*) &(redshift), sizeof(float));
    file.read((char*) &(om_m0), sizeof(float));
    file.read((char*) &(om_lam0),sizeof(float));

    pParticles = (Particle **) malloc(sizeof(void*)*numParticles);

    for(unsigned int n=0; n<numParticles; n++) {
      pParticles[n] = new Particle(&file,file_ver_num);
    }

    if(file_ver_num>1) {
      char *temp_chararr= NULL;
      int nArr;
      file.read((char*) &nArr, sizeof(int));
      if(getVerbose()>=2)
	cerr << "done main section!" << endl << "BaseSimSnap: " << nArr << " extra arrays to read" << endl;
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
	if(getVerbose()>=2)
	  cerr << "BaseSimSnap: array " << n << " - " << longname << "(" << shortname << ")...";
	SimanArray & newArr = createArray(shortname,longname,Unit(&file, file_ver_num));
	newArr.read(file,file_ver_num);
      
	if(getVerbose()>=2)
	  cerr << "done!" << endl;
      }
    } 

    if(getVerbose()>=2)
      cerr << "Units: " << lenUnits << " | " << velUnits << " | " << denUnits << " | " << enUnits << " | " << massUnits << endl;

    if(getVerbose()==1)
      cerr << "done!" << endl;

    setupReservedArrays();

  }

  BaseSimSnap::BaseSimSnap(const SimSnap *copyFrom) {
    pAutoParent=NULL;
    numParticles = copyFrom->getNumParticles();
    int numArrays = copyFrom->getNumArrays();
    for(int i=0; i<numArrays; i++) {
      SimanArray & newarr = createArray(copyFrom->getArrayName(i),copyFrom->getArrayLongName(i),
					     copyFrom->getArrayUnits(i));
      newarr = copyFrom->getConstArray(i);
    }
    

    pParticles = (Particle **) malloc(sizeof(void*)*numParticles);
  
    for(unsigned int n=0; n<numParticles; n++) {
      const Particle *toCopy = copyFrom->getConstParticle(n);
      pParticles[n] = new Particle(*toCopy);
    
    }

    lenUnits = copyFrom->getDistanceUnits();
    velUnits = copyFrom->getVelocityUnits();
    denUnits = copyFrom->getDensityUnits();
    massUnits = copyFrom->getMassUnits();
    enUnits = copyFrom->getEnergyUnits();

    boxSize = copyFrom->getBoxSize();
    hubble = copyFrom->getHubble();
    redshift = copyFrom->getRedshift();

    setupReservedArrays();


  }

  BaseSimSnap::BaseSimSnap(Particle **pParticlesi, int numParticlesi, Unit len, Unit mass) {

  
    pParticles = pParticlesi;
    numParticles = numParticlesi;

  
    if(len==Unit()) 
      lenUnits = Unit("kpc");
    else
      lenUnits = len;

    if(mass==Unit()) 
      massUnits = Unit("Msol");
    else
      massUnits = mass;


    velUnits = Unit("km s^-1");

    denUnits = massUnits/pow(lenUnits,3);
    enUnits = pow(velUnits,2);
    
    setupReservedArrays();
 
  }

  BaseSimSnap::BaseSimSnap() {

    pParticles=NULL;
    numParticles = 0;

    lenUnits = Unit("kpc");
    massUnits = Unit("msol");
    velUnits = Unit("km s^-1");
    denUnits = massUnits/pow(lenUnits,3);
    enUnits = pow(velUnits,2);

    setupReservedArrays();

  }


  BaseSimSnap* BaseSimSnap::makeSlab(unsigned int nPart, float depth, float extent, float density, float temp, int genType, Unit length, Unit mass) {

    // make a uniform density slab out of nPart particles

    Particle **pParticles;

    if(genType==genRandom) {
    
      pParticles = (Particle **) malloc(sizeof(void*)*nPart);
      float partMass = density*extent*extent*depth/((float)nPart);
      for(unsigned int n=0; n<nPart; n++) {
	pParticles[n] = new Particle(((float)rand()/(float)RAND_MAX)*extent-extent/2,((float)rand()/(float)RAND_MAX)*extent-extent/2,((float)rand()/(float)RAND_MAX)*depth-depth/2,0,0,0,partMass,Particle::gas,temp,density);
      }
    } else {
      float c = std::pow((double)(nPart/(depth*extent*extent)),1./3.);
      int nz= (int)(c*depth);
      int nx= (int)(c*extent);
      int ny = nx;
      nPart = nz*nx*ny;
      float partMass = density*extent*extent*depth/((float)nPart);
      pParticles = (Particle **) malloc(sizeof(void*)*nPart);
      int i=0;
      for(int x=0;x<nx;x++) {
	for(int y=0;y<ny;y++) {
	  for(int z=0; z<nz;z++) {
	    pParticles[i] = new Particle(x/c-extent/2.,y/c-extent/2.,z/c-depth/2.,0,0,0,partMass,Particle::gas,temp,density);
	    i++;
	  }
	}
      }
      if(getVerbose()>=1)
	cerr << "BaseSimSnap: generated regular slab with " << nPart << " particles (" << nx << " x " << ny << " x " << nz << ")" << endl;
    }
    return new BaseSimSnap(pParticles,nPart,length,mass);
  
  }

  void BaseSimSnap::setupReservedArrays() {
   
    extraMap["x"]=new SimanArrayVirtual(this,&Particle::x);
    extraMap["y"]=new SimanArrayVirtual(this,&Particle::y);
    extraMap["z"]=new SimanArrayVirtual(this,&Particle::z);
    extraMap["vx"]=new SimanArrayVirtual(this,&Particle::vx);
    extraMap["vy"]=new SimanArrayVirtual(this,&Particle::vy);
    extraMap["vz"]=new SimanArrayVirtual(this,&Particle::vz);
    extraMap["temp"]=new SimanArrayVirtual(this,&Particle::temp);
    extraMap["mass"]=new SimanArrayVirtual(this,&Particle::mass);
    extraMap["rho"]=new SimanArrayVirtual(this,&Particle::rho);
    extraMap["ne"]=new SimanArrayVirtual(this,&Particle::ne);
    extraMap["metal"]=new SimanArrayVirtual(this,&Particle::metal);
    extraMap["u"]=new SimanArrayVirtual(this,&Particle::u);
    extraMapName["x"]="auto.x";
    extraMapName["y"]="auto.y";
    extraMapName["z"]="auto.z";
    extraMapName["vx"]="auto.vx";
    extraMapName["vy"]="auto.vy";
    extraMapName["vz"]="auto.vz";
    extraMapName["temp"]="auto.temp";
    extraMapName["mass"]="auto.mass";
    extraMapName["rho"]="auto.rho";
    extraMapName["ne"]="auto.ne";
    extraMapName["metal"]="auto.metal";
    extraMapName["u"]="auto.u";
    extraMap_reserved_end=extraMap.size();
  }

  void BaseSimSnap::allocateMemory() {
    if(getVerbose()>=2)
      cerr << "BaseSimSnap: allocating memory...";
    pParticles = (Particle **) malloc(sizeof(void*)*numParticles);
    for(unsigned int n=0; n<numParticles; n++) {
      pParticles[n] = new Particle();
    }
    if(getVerbose()>=2)
      cerr << "done!" << endl;
  }



  void BaseSimSnap::convertUnits(Unit distance, Unit mass, Unit velocity, Unit density, Unit energy) {

    bumpVersion();

    if(distance==Unit()) distance=Unit("kpc");
    if(mass==Unit()) mass=Unit("Msol");
    if(velocity==Unit()) velocity=Unit("km s^-1");
    if(density==Unit()) density=mass/(distance*distance*distance);
    if(energy==Unit()) energy=velocity*velocity;
  
    if(getVerbose()>=1)
      cerr << "BaseSimSnap: convertUnits - dist: "<<distance<<" mass: "<<mass<<" vel: "<<velocity<<" den: "<<density<<" ener: "<<energy<<"...";

    double dist_ratio = lenUnits.convertTo(distance,this);
    double mass_ratio = massUnits.convertTo(mass,this);
    double vel_ratio = velUnits.convertTo(velocity,this);
    double den_ratio = denUnits.convertTo(density,this);
    double en_ratio = enUnits.convertTo(energy,this);


    if(getVerbose()>=2)
      cerr << endl << "BaseSimSnap: conversion ratios respectively: " << dist_ratio << " " << mass_ratio << " " << vel_ratio << " " << den_ratio << " " << en_ratio << "...";
  
 
  
    boxSize*=dist_ratio;

    for(int i=numParticles-1;i>=0;--i) {
      Particle *p = pParticles[i];
      p->x*=dist_ratio;
      p->y*=dist_ratio;
      p->z*=dist_ratio;
      if(pSmooth!=NULL) pSmooth[i]*=dist_ratio;
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
  
    if(getVerbose()>=1)
      cerr << "done!" << endl;

  }

  BaseSimSnap::~BaseSimSnap() {
    
    if(getVerbose()>2)
      cerr << "BaseSimSnap::~BaseSimSnap [" << this << "]" << endl;

    for(unsigned int n=0; n<numParticles; n++) {
      delete pParticles[n];
    }
    free(pParticles);
    
    map<string, SimanArray*>::iterator i;
    i=extraMap.begin();
    for(i=extraMap.begin() ; i!=extraMap.end(); i++) {
      delete (*i).second;
    }
  }

  bool BaseSimSnap::Load(bool fakeload) {
    return true;
  }

  Particle* BaseSimSnap::getParticle(unsigned int n) {
    bumpVersion();
#ifndef SIMAN_UNSAFE_FASTER
    if(n>=getNumParticles())
      throw(std::out_of_range(boost::lexical_cast<string>(n)));
#endif
    return pParticles[n];
  }


  const Particle * BaseSimSnap::getConstParticle(unsigned int n) const {
#ifndef SIMAN_UNSAFE_FASTER
    if(n>=getNumParticles())
      throw(std::out_of_range(boost::lexical_cast<string>(n)));
#endif
    return pParticles[n];
  }


  void BaseSimSnap::releaseParticle(unsigned int n) {
    // Do Nothing...
  }

  unsigned int BaseSimSnap::getNumParticles() const {
    return numParticles;
  }

  float BaseSimSnap::getBoxSize() const {
    return boxSize;
  }

  float BaseSimSnap::getHubble() const {
    return hubble;
  }

  void BaseSimSnap::setHubble(float h) {
    hubble=h;
  }


  void BaseSimSnap::setRedshift(float z) {
    redshift= z;
  }


  void BaseSimSnap::setBoxSize(float h) {
    boxSize=h;
  }
  
  void BaseSimSnap::setOmegaM0(float a) {
    om_m0=a;
  }

  void BaseSimSnap::setOmegaLambda0(float a) {
    om_lam0=a;
  }

  float BaseSimSnap::getRedshift() const {
    return redshift;
  }

  float BaseSimSnap::getOmegaM0() const {
    return om_m0;
  }

  float BaseSimSnap::getOmegaLambda0() const {
    return om_lam0;
  }

  int BaseSimSnap::deReference(int i, int n) const {
    // if(n!=1)
    //  cerr << "BaseSimSnap: Warning - tried to deReference beyond base level" << endl;
    return i;
  }


  int BaseSimSnap::deReference(int i, SimSnap *pRel) const {

    if(pRel!=static_cast<const SimSnap*>(this))
      throw(SimanException("Dereferencing to non-ancestorial object"));
    else
      return i;
  }


  void BaseSimSnap::nativeWrite(const SimSnap *sim, string filename) {

    if(getVerbose()>=1)
      cerr << "BaseSimSnap: writing " << filename << "...";
    ofstream file(filename.c_str(),ios::binary);

    char id[4]={'S','I','B','I'}; // identifier "SIman BInary" for autoloader

    file.write(id, 4);
  
    int file_ver_num = 4;

    file.write((char*) &file_ver_num, sizeof(int)); // version number of filetype

    unsigned int nPart = sim->getNumParticles();
  
    if(getVerbose()>1)
      cerr << "BaseSimSnap: writing " << nPart << " particles..." << endl;

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
      sim->getConstParticle(n)->nativeWrite(&file);
    }

    int numa = sim->getNumArrays();
    file.write((char*) &numa, sizeof(int));

    if(sim->getNumArrays()>0) {
      if(getVerbose()>=2)
	cerr << "done main section!" << endl;

      for(int n=0;n<sim->getNumArrays();n++) {
	string shortname = sim->getArrayName(n);
	int shortname_len = shortname.size();
	string longname = sim->getArrayLongName(n);
	int longname_len = longname.size();
	if(getVerbose()>=2)
	  cerr << "BaseSimSnap: writing " << longname << " (" << shortname << ")...";
	file.write((char*) &shortname_len, sizeof(int));
	file.write(shortname.c_str(),shortname.size()+1);
	file.write((char*) &longname_len, sizeof(int));
	file.write(longname.c_str(),longname.size()+1);
	sim->getArrayUnits(n).nativeWrite(&file);
	sim->getConstArray(n).write(file);

	if(getVerbose()>=2)
	  cerr << "done!" << endl;
      }

      if(getVerbose()==1)
	cerr << "done!" << endl;

    } else {
      if(getVerbose()>1)
	cerr << "done!" << endl;
    }
  }


  SimanArray & BaseSimSnap::createArray(string name, string fullname, Unit arrunits) {
    if(name=="x" || name=="y" || name=="z" || name=="vx" || name=="vy" || name=="vz" || name=="temp" || name=="rho" || name=="mass" || name=="ne" || name=="type")
      throw(SimanException("Reserved name `"+name+"'"));

    bumpVersion();
    if(extraMap.count(name)>0) {

      if(getVerbose()>3) // not really a problem, never report unless on silly-verbose mode
	cerr << "BaseSimSnap::createArray: warning - array " << name << " already exists; returning existing array" << endl;

      return *(extraMap[name]);
    } else {
      SimanArray* newArray = new SimanArray(getNumParticles(),this,arrunits);
      extraMap[name] = newArray;
      extraMapName[name] = fullname;
      extraList.push_back(extraMap.find(name));
      return *newArray;
    }
  }

  void BaseSimSnap::destroyArray(string name) {
     if(name=="x" || name=="y" || name=="z" || name=="vx" || name=="vy" || name=="vz" || name=="temp" || name=="rho" || name=="mass" || name=="ne" || name=="type")
      throw(SimanException("Reserved name `"+name+"'"));

     map<string,SimanArray*>::iterator i = extraMap.find(name);
     if(i==extraMap.end())
       throw UnknownArray(name);
     
     std::list<std::map<std::string, SimanArray*>::iterator>::iterator j = find(extraList.begin(),extraList.end(),i);
     extraList.erase(j);

     SimanArray *deleting = (*i).second;
     delete deleting;
     
     extraMap.erase(i);

     extraMapName.erase(name);
   
     bumpVersion();
  }

  int BaseSimSnap::getArrayIndex(string name) const {
    
    map<string,SimanArray*>::const_iterator i = extraMap.find(name);
    if(i==extraMap.end())
      throw UnknownArray(name);
    
    std::list<std::map<std::string, SimanArray*>::iterator>::const_iterator j;
    int index=0;
    for(j=extraList.begin();j!=extraList.end();j++, index++) {
      if((*j)==i) return index;
    }
    throw UnknownArray(name);
  }

  Unit BaseSimSnap::getArrayUnits(string name) const {
    return getConstArray(name).getUnits();
  }

  Unit BaseSimSnap::getArrayUnits(int n) const {
    return getConstArray(n).getUnits();
  }

  string BaseSimSnap::getArrayLongName(string shortname) const {
    map<string,string>::const_iterator i = 
      extraMapName.find(shortname);
    if(i==extraMapName.end())
      throw(UnknownArray(shortname));

    return (*i).second;
  }

  string BaseSimSnap::getArrayLongName(int n) const {
    std::list<std::map<std::string, SimanArray*>::iterator>::const_iterator i;
    int o=0;
    if((unsigned int)n>=extraList.size())
      throw UnknownArray("");

    for(i=extraList.begin(); o!=n; o++, i++);
    
    return getArrayLongName((*(*i)).first);
  }

  

  string BaseSimSnap::getArrayName(int n) const {
    
    std::list<std::map<std::string, SimanArray*>::iterator>::const_iterator i; // Tee hee
    int o=0;
    if((unsigned int)n>=extraList.size())
      throw UnknownArray("");

    for(i=extraList.begin(); o!=n; i++, o++);
    
    return (*(*i)).first;
   
  }


  SimanArray & BaseSimSnap::getArray(string name) {

    map<string,SimanArray*>::const_iterator i = extraMap.find(name);
    if(i==extraMap.end())
      throw(UnknownArray(name));
    else {
      return *((*i).second);
    }

  }

  SimanArray & BaseSimSnap::getArray(int n) {
    std::list<std::map<std::string, SimanArray*>::iterator>::iterator i; // Tee hee
    
    int o=0;
    if((unsigned int)n>=extraList.size())
      throw UnknownArray("");

    for(i=extraList.begin(); o!=n; i++, o++);
    
    return *((*(*i)).second);
   
    
  }


  const SimanArray & BaseSimSnap::getConstArray(string name) const {
    map<string,SimanArray*>::const_iterator i = extraMap.find(name);
    if(i==extraMap.end())
      throw(UnknownArray(name));
    else
      return *((*i).second);
  }

  const SimanArray & BaseSimSnap::getConstArray(int n) const {
    std::list<std::map<std::string, SimanArray*>::iterator>::const_iterator i; // Tee hee
    
    int o=0;
    if((unsigned int)n>=extraList.size())
      throw UnknownArray("");

    for(i=extraList.begin(); o!=n; i++, o++);
    
    return *((*(*i)).second);

  }

  int BaseSimSnap::getNumArrays() const {
    return extraList.size();
  }

}
