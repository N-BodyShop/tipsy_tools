// simsnap.cpp - part of SimAn Simulation Analysis Library
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
#include "endian.hpp"
#include <limits>

namespace siman {


  SimSnap::SimSnap() {
#ifdef SIMAN_TRACE
    cerr << "SimSnap::SimSnap" << endl;
#endif 
    pAutoParent=NULL;
    nPartSPH = 6; // number of particles to use in SPH-like calculations
    pSmooth = NULL;
    pKernel=NULL;
    version=0;
  }

  SimSnap::~SimSnap() {
#ifdef SIMAN_TRACE
    cerr << "SimSnap::~SimSnap" << endl;
#endif
    if(pSmooth!=NULL) free(pSmooth);
    if(pKernel!=NULL) delete(pKernel);
  }

  // AUTO-PARENTING FUNCTIONS
  //
  // These have no meaning in the base class, but derived classes can set pAutoParent
  // to save writing overrides to pass functions on to parent files (e.g. see Subset)

  SimSnap* SimSnap::getParent() const {
    return pAutoParent;
  }

  float SimSnap::getOmegaM0() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getOmegaM0()";
      return 0.3;
    } else {
      return pAutoParent->getOmegaM0();
    }
  }

  float SimSnap::getOmegaLambda0() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getOmegaLambda0()";
      return 0.7;
    } else {
      return pAutoParent->getOmegaLambda0();
    }
  }

  int SimSnap::deReference(int i, int n) const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method deReference();";
      return 0;
    } else {
      return pAutoParent->deReference(i,n);
    }
  }
  
  int SimSnap::deReference(int i, SimSnap *pRel) const {

    if(pRel!=this && pAutoParent!=NULL)
      return deReference(i,pRel);
    
    if(pRel!=this)
      throw(SimanException("Dereferencing to non-ancestorial object"));
    else
      return i;
  }

  void SimSnap::convertUnits(Unit distance, Unit mass, Unit velocity, Unit density, Unit energy) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method convertUnits" << endl;
      return;
    } else {
      pAutoParent->convertUnits(distance,mass,velocity,density,energy);
    }
  }

  void SimSnap::TempFromU() throw(UnitsError) {

    if(getVerbose()>1)
      cerr << "SimSnap: calculating temperatures from U and Ne..."; 

    Unit unit_spec_en_cgs= Unit("erg g^-1");
    Unit unit_spec_en_int = getEnergyUnits();
  
    double en_conv = unit_spec_en_int.convertTo(unit_spec_en_cgs);

 
    const float Y = 0.26;
    const int numParticles = getNumParticles();

    for(int i=0;i<numParticles;i++) {
    
      Particle *p = getParticle(i);


      if(p->type==Particle::gas)  {
      
	double MeanWeight = 1/(1-0.75*Y+p->ne*(1-Y)) * constants::protonMassInG;
      
	// internal energy -> cgs units
      
	float u  = p->u * en_conv;
	// = UnitEnergy_in_cgs/ UnitMass_in_g;
      
	float gamma= 5.0/3;
      
	// temperature in K:
      
	p->temp = MeanWeight/constants::boltzmannInErgPerK * (gamma-1) * u;
      
      }
    
    } // for i

    if(getVerbose()>1)
      cerr << "done!" << endl;
  }

  void SimSnap::UFromTemp() throw(UnitsError) {
   
    if(getVerbose()>1) 
      cerr << "SimSnap: calculating U from temperatures and Ne...";

    Unit unit_spec_en_cgs= Unit("erg g^-1");
    Unit unit_spec_en_int = getEnergyUnits();
  
    double en_conv = unit_spec_en_int.convertTo(unit_spec_en_cgs);
  
    const float Y = 0.26;
    const int numParticles = getNumParticles();

    for(int i=0;i<numParticles;i++) {
    
      Particle *p = getParticle(i);
      if(p->type==Particle::gas)  {
      
	double MeanWeight = 1/(1-0.75*Y+p->ne*(1-Y)) * constants::protonMassInG;
      
	float gamma= 5.0/3;
      
	// temperature in K:
      
      
	float u = (p->temp * constants::boltzmannInErgPerK)/(MeanWeight*(gamma-1.));
	p->u = u/ en_conv;
      
	// cerr << MeanWeight << " " << p->temp << " " << en_conv << "\t" << p->u << endl;
      }
    
    } // for i

    if(getVerbose()>1)
      cerr << "done!" << endl;

  }

  void SimSnap::makeMu() {
    SimanArray &mu = createArray("mu","mean atomic mass");
    float Y = getHeliumMassFrac();
    for(unsigned int i=0; i<getNumParticles(); i++) {
      const Particle *p=getConstParticle(i);
      mu[i] =  1/(1-0.75*Y+p->ne*(1-Y));
    }
  }

  void SimSnap::makeP() {
    Unit targetUnit("kg s^-2 m^-1");
    SimanArray &pr = createArray("p","pressure",targetUnit);
    Unit ratio("8254.41 m^2 s^-2");
    float Y = getHeliumMassFrac();

    double convRatio = ratio.convertTo(targetUnit/getDensityUnits());
    
    for(unsigned int i=0; i<getNumParticles(); i++) {
      const Particle *p=getConstParticle(i);
      float mu = 1/(1.-0.75*Y+p->ne*(1-Y));
      pr[i]=convRatio * p->rho * p->temp / mu;
    }
  }

  unsigned int SimSnap::getNumParticles() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getNumParticles()\n";
      return 0;
    } else {
      return pAutoParent->getNumParticles();
    }
  }


  Unit SimSnap::getDistanceUnits() const {
    if(pAutoParent==NULL) {
  
      return lenUnits;
    } else {
      return pAutoParent->getDistanceUnits();
    }
  }

  Unit SimSnap::getVelocityUnits() const {
    if(pAutoParent==NULL) {
 
      return velUnits;
    } else {
      return pAutoParent->getVelocityUnits();
    }
  }


  Unit SimSnap::getDensityUnits() const {
    if(pAutoParent==NULL) {
   
      return denUnits;
    } else {
      return pAutoParent->getDensityUnits();
    }
  }


  Unit SimSnap::getEnergyUnits() const {
    if(pAutoParent==NULL) {
    
      return enUnits;
    } else {
      return pAutoParent->getEnergyUnits();
    }
  }


  Unit SimSnap::getMassUnits() const {
    if(pAutoParent==NULL) {
  
      return massUnits;
    } else {
      return pAutoParent->getMassUnits();
    }
  }


  float SimSnap::getBoxSize() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getBoxSize()\n";
      return 0;
    } else {
      return pAutoParent->getBoxSize();
    }

  }

  float SimSnap::getHubble() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getHubble()\n";
      return 0;
    } else {
      return pAutoParent->getHubble();
    }

  }

  void SimSnap::setHubble(float h) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method setHubble()\n";
    
    } else {
      pAutoParent->setHubble(h);
    }

  }


  void SimSnap::setOmegaLambda0(float h) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method setOmegaLambda0()\n";
    
    } else {
      pAutoParent->setOmegaLambda0(h);
    }

  }


  void SimSnap::setOmegaM0(float h) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method setOmegaM0()\n";
    
    } else {
      pAutoParent->setOmegaM0(h);
    }

  }


  void SimSnap::setBoxSize(float h) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method setBoxSize()\n";
    
    } else {
      pAutoParent->setBoxSize(h);
    }

  }


  void SimSnap::setRedshift(float z) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method setRedshift()" << endl;
    
    } else {
      pAutoParent->setRedshift(z);
    }

  }

  float SimSnap::getRedshift() const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getRedshift()\n";
      return 0;
    } else {
      return pAutoParent->getRedshift();
    }
  }


  float SimSnap::scaleToAngDiDist(float a_em) const {
    return scaleToAngDiDist(1./(1.+getRedshift()),a_em);
  }

  float SimSnap::scaleToComovingDist(float a_em) const {
    return scaleToComovingDist(1./(1.+getRedshift()),a_em);
  }

  float SimSnap::scaleToComovingDist(float a_current, float a_em) const {
    float delta_a = 0.001;
    float r = 0.;
    for(float ai = a_em; ai<a_current; ai+=delta_a) {
      r+=delta_a/(ai*ai*scaleToHubble1e10yr(ai));
    }
    float c = Unit("c").convertTo(getDistanceUnits()/Unit("1e10 yr"),this);
    return r * c;
  }

  float SimSnap::scaleToAngDiDist(float a_current, float a_em) const {
    
    return a_em*scaleToComovingDist(a_current,a_em);
    
  }

  float SimSnap::scaleToHubble1e10yr(float a) const {
    if(a<0) a=1./(1.+getRedshift());
    return 1.02268944 * getHubble()*sqrt(getOmegaM0()/(a*a*a) + getOmegaLambda0() + (1.-getOmegaM0()-getOmegaLambda0())/(a*a));
        
  }

  
  float SimSnap::scaleToAge1e10yr(float a) const {
    if(a<0) a=1./(1.+getRedshift());
    float tage=0.;
    float delta_a = 0.001;
    for(float ai = delta_a; ai<a; ai+=delta_a) {
      tage+=delta_a/(ai*scaleToHubble1e10yr(ai));
    }
    return tage;
  }


  float SimSnap::scaleToHubble(float a) const {
    double hub = scaleToHubble1e10yr(a);
    return hub*Unit("1e10 yr^-1").convertTo(getVelocityUnits()/getDistanceUnits(),this);
    // hub is currently in units of 1.e10 years
    
  }


  float SimSnap::scaleToAge(float a) const {
    float tage = scaleToAge1e10yr(a);
    return tage*Unit("1e10 yr").convertTo(getDistanceUnits()/getVelocityUnits(),this);
  }
  


  Particle * SimSnap::getParticle(unsigned int id) {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getParticle()\n";
      return 0;
    } else {
      return pAutoParent->getParticle(id);
    }

  }


  Particle & SimSnap::getParticleRef(unsigned int id) {
    return *getParticle(id);
  }


  const Particle * SimSnap::getConstParticle(unsigned int id) const {
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method getConstParticle()\n";
      return 0;
    } else {
      return pAutoParent->getConstParticle(id);
    }

  }

  const Particle & SimSnap::getConstParticleRef(unsigned int id) const {
    return *getConstParticle(id);
  }

  void SimSnap::releaseParticle(unsigned int id) {
  
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method releaseParticle()\n";
    } else {
      pAutoParent->releaseParticle(id);
    }

  
  }

  void SimSnap::releaseAllParticles() {
    
    if(pAutoParent==NULL) {
      cerr << "SimSnap: Base class erroneously called for virtual method releaseAllParticles()\n";
    } else {
      pAutoParent->releaseAllParticles();
    }

  }

  float SimSnap::getTotalMass() const {

    float mass = 0;
    const Particle *particle;

    for(unsigned int n=0;n<getNumParticles();n++) {
      particle = getConstParticle(n);
      mass +=particle->mass;
   
    }

    return mass;
  }

  
  void SimSnap::centreOn(vector<double> coords) {
    if(coords.size()!=3)
      throw SimanException("SimSnap::centreOn - expected length 3 array");
    centreOn(coords[0],coords[1],coords[2]);
  }
  

  void SimSnap::centreOn(float cx, float cy, float cz) {
    
    if(getBoxSize()>0) 
      transform(Wrap(-getBoxSize()/2,+getBoxSize()/2) * Translation(-cx,-cy,-cz));
    else
      transform(Translation(-cx,-cy,-cz));
  }

  
  
  void SimSnap::centreOnVel(vector<double> coords) {
    if(coords.size()!=3)
      throw SimanException("SimSnap::centreOnVel - expected length 3 array");
    centreOnVel(coords[0],coords[1],coords[2]);
  }
  
  void SimSnap::centreOnVel(float cx, float cy, float cz) {
    
    transform(VelocityTranslation(-cx,-cy,-cz));

  }
  

   
  vector<double> SimSnap::centreOfMass() const {
    vector<double> res;
    float x,y,z;
    centreOfMass(x,y,z);
    res.push_back((double)x);
    res.push_back((double)y);
    res.push_back((double)z);
    return res;
  }
   

  void SimSnap::centreOfMass(float &cx, float &cy, float &cz) const {
 

    double cxd=0., cyd=0., czd=0.;
  
    double totMass = getTotalMass();
  
    int n=0;
    int numParticles = getNumParticles();



    const Particle *particle;

    // If the centre of mass is actually somewhere a long way off, but the
    // particles have little spread, numerical errors quickly mount. Thus
    // we subtract a "representative" position vector before averaging,
    // and add it back on at the end.

    particle = getConstParticle(0);

    double cxd_or = particle->x;
    double cyd_or = particle->y;
    double czd_or = particle->z;

    for(n=0;n<numParticles;n++) {
      particle = getConstParticle(n);

      cxd+=((double)particle->x-cxd_or)*((double)particle->mass/totMass);
      cyd+=((double)particle->y-cyd_or)*((double)particle->mass/totMass);
      czd+=((double)particle->z-czd_or)*((double)particle->mass/totMass);


    }

 
    cx = (float)(cxd + cxd_or);
    cy = (float)(cyd + cyd_or);
    cz = (float)(czd + czd_or);
 
  }
  
  void SimSnap::centreOfMassVel(float &cvx, float &cvy, float &cvz) const {

    cvx = 0.;
    cvy = 0.;
    cvz = 0.;

    double cxd=0., cyd=0., czd=0.;
  
    double totMass = getTotalMass();
  
    int n=0;
    int numParticles = getNumParticles();


    const Particle *particle;

    for(n=0;n<numParticles;n++) {
      particle = getConstParticle(n);

      cxd+=(double)particle->vx*(double)particle->mass/totMass;
      cyd+=(double)particle->vy*(double)particle->mass/totMass;
      czd+=(double)particle->vz*(double)particle->mass/totMass;


    }


    // cxd/=(double)numParticles;
    //cyd/=(double)numParticles;
    //czd/=(double)numParticles;

    cvx = (float)cxd;
    cvy = (float)cyd;
    cvz = (float)czd;
 
  }

  
  
  vector<double> SimSnap::centreOfMassVel() const {
    vector<double> res;
    float x,y,z;
    centreOfMassVel(x,y,z);
    res.push_back((double)x);
    res.push_back((double)y);
    res.push_back((double)z);
    return res;
  }
   


  vector<double> SimSnap::shrinkSphereCentre(float sf, int mp, float ir) const {
    vector<double> res;
    float x,y,z;
    shrinkSphereCentre(x,y,z,sf,mp,ir);
    res.push_back((double)x);
    res.push_back((double)y);
    res.push_back((double)z);
    return res;
  }


  vector<double> SimSnap::shrinkSphereCentreVel(float sf, int mp, float ir) const {
    vector<double> res;
    float x,y,z;
    shrinkSphereCentreVel(x,y,z,sf,mp,ir);
    res.push_back((double)x);
    res.push_back((double)y);
    res.push_back((double)z);
    return res;
  }

  void SimSnap::shrinkSphereCentre(float &cx, float &cy, float &cz, float shrinkFactor, int minParticles, float initial_radius) const {

    // Find centre using standard shrinking sphere method (Power et al 2003/4?)
    //
    // shrinkFactor - length scale factor to shrink sphere by at each stage
    // minParticles - assume centre found when minParticles remains
    // initial_radius - defaults to negative (=> determine box size), but in
    //                  recursive section is used to pass current radius
    //                  (thus no expensive calculation to determine it is
    //                   necessary!)


    // get COM, passing cx/cy/cz by reference

    centreOfMass(cx,cy,cz);

    //  cerr << "SS: " << getNumParticles() << " " << shrinkFactor << " " <<minParticles << endl;;

    if(getNumParticles()<=(unsigned int)minParticles)
      return;

 
    float r = 0.;

    if(initial_radius>0.)  
      r = initial_radius*shrinkFactor; 
    else 
      r = getApparentBoxSize()*shrinkFactor;


    Sphere sphereFilter(cx,cy,cz,r);  
    
    // not sure how to avoid following const_cast. Generating the
    // subset is effectively const on its argument, but since
    // a subset can later update the argument which it holds,
    // it's not sensible to take a const argument. Here we know
    // we will not make any changes so it's ok.
    Subset sphereParticles(const_cast<SimSnap*>(this),sphereFilter);


    sphereParticles.shrinkSphereCentre(cx,cy,cz,shrinkFactor,minParticles,r);

    // recursion from here will automatically return our required cx,cy,cz
  
  
  }



  void SimSnap::shrinkSphereCentreVel(float &cvx, float &cvy, float &cvz, float shrinkFactor, int minParticles, float initial_radius) const {

    centreOfMassVel(cvx,cvy,cvz);

    if(getNumParticles()<=(unsigned int)minParticles)
      return;
 
    float r = 0.;

    if(initial_radius>0.)  
      r = initial_radius*shrinkFactor; 
    else 
      r = getMaximalVelocity()*shrinkFactor;


    VelocitySphereFilter sphereFilter(cvx,cvy,cvz,r);  
    
    Subset sphereParticles(const_cast<SimSnap*>(this),sphereFilter);

    sphereParticles.shrinkSphereCentreVel(cvx,cvy,cvz,shrinkFactor,minParticles,r);  
  }


  float SimSnap::power(float lambda) const {
    complex <float> sum;

    unsigned int numParticles = this->getNumParticles();
 
    for(unsigned int n=0; n<numParticles; n++) {
      const Particle *p=getConstParticle(n);
      complex <float> exponent(0,lambda*p->x);
      sum += exp(exponent);
      
    }
    
    
    return std::pow(abs(sum),2);
  }


  void SimSnap::transform(const Transformation &t) {
    
    vector<SimanArray*> extraAddr;
    vector<float> extraVal;
    vector<string> extraData;
    
    t.extraBlocks(extraData, *this);

    // store references for each extra block requested
    if(extraData.size()!=0) {
      vector<string>::const_iterator i;
      for(i=extraData.begin();i!=extraData.end();i++) {
	extraAddr.push_back(&(getArray(*i)));
	extraVal.push_back(0);
      }
    }
    
    int num = getNumParticles();
    for(int n=0; n<num; n++) {
      
      for(unsigned int i=0; i<extraAddr.size(); i++) {

	// update previous values of extra arrays
	if(n!=0)
	  (*(extraAddr[i]))[n-1]=extraVal[i];

	// get values for this particle
	extraVal[i]=(*(extraAddr[i]))[n];
      }

      t.transform(*(getParticle(n)),extraVal);     
    }
   
          
    if(num!=0) {
      for(unsigned int i=0; i<extraAddr.size(); i++) {
	
	// update previous values of extra arrays
	(*(extraAddr[i]))[num-1]=extraVal[i];
	
      }
    }
    
  }

  SimSnap * SimSnap::copy() const {
    return new BaseSimSnap(this);
  }

  SimSnap * SimSnap::copyTransform(const Transformation &t) const {
    SimSnap *s = new BaseSimSnap(this);
    s->transform(t);
    return s;
  }

  Subset * SimSnap::subset(const Filter &f) {
    return new Subset(this,f);
  }

  SimSnap * SimSnap::copySubset(const Filter &f) const {
    Subset *s = new Subset(const_cast<SimSnap*>(this),f);
    SimSnap *r = s->copy();
    delete s;
    return r;
  }
  
 

  int SimSnap::determineType(std::string filename) {


  

    // currently very crude file type detector

    ifstream file(filename.c_str(),ios::binary);
  
    if(file.bad()) return unknown;

    int test_word;
    char* tw = (char*) &test_word;

    file.read((char*) &test_word,4);

    if(test_word==sizeof(gadget_header)) 
      return gadget;

    if(tw[0]=='S' && tw[1]=='I' && tw[2]=='B' && tw[3]=='I')
      return native;

    if(tw[0]=='#' && tw[1]=='!')
      return script;

    endian::flip4(&test_word);

    if(test_word==sizeof(gadget_header)) 
      return gadget + endianWrong;

  
    if(tw[0]=='#' && tw[1]=='!')
      return script;

    endian::flip4(&test_word);



    tipsy::header test_header;

    file.seekg(0,ios_base::beg);
  
    file.read((char*) &test_header,sizeof(test_header));

    if(test_header.nsph + test_header.ndark + test_header.nstar == test_header.nbodies)
      return tipsy;

    endian::flip4(&test_header.nsph);
    endian::flip4(&test_header.ndark);
    endian::flip4(&test_header.nstar);
    endian::flip4(&test_header.nbodies);

  
    if(test_header.nsph + test_header.ndark + test_header.nstar == test_header.nbodies)
      return tipsy + endianWrong;


    // don't know what this is

    return unknown;


  }

  SimSnap * SimSnap::loadFile(std::string filename) {
    int type = determineType(filename);

    if((type & gadget) > 0) {
      GadgetSnap *pGadg;
      if(getVerbose()>1)
	cerr << "SimSnap: trying to load " << filename << " as a gadget file" << endl;
      pGadg = new GadgetSnap(filename.c_str());
    
      return pGadg;
    }

    if((type & tipsy) > 0) {
      TipsySnap *pTipsy;
      if(getVerbose()>1)
	cerr << "SimSnap: trying to load " << filename << " as a tipsy-binary file" << endl;
      pTipsy = new TipsySnap(filename.c_str());
    
      return pTipsy;
    }

    if((type & native) > 0) {
      BaseSimSnap *pNative;
      if(getVerbose()>1)
	cerr << "SimSnap: trying to load " << filename << " as a SimAn file" << endl;
      pNative = new BaseSimSnap(filename);
      return pNative;
    }

    if((type & script) > 0) {
     
      SimSnap *pRet = NULL;
      Scripted script(filename);
      pRet = dynamic_cast<SimSnap*>(script.getReturnValue()); 
      // pRet=NULL if this is impossible
   
      if(pRet) return pRet; else {
	if(getVerbose()>0)
	  cerr << "SimSnap: returned type from script " << filename << " is not derived from SimSnap" << endl;
	return NULL;
      }
    }

    // fallthrough

    if(getVerbose()>0) {
      cerr << "SimSnap: don't know how to load file "<<filename << endl;
    }
    throw(FileError(filename));
    return NULL;
  }

  int SimSnap::getNearestNeighbour(const Particle &us) const {
  
    float bestDistance(0.);

    int nParts = getNumParticles();
    int nearest = -1;

    for(int n=0;n<nParts;n++) {
    
      const Particle *candidate = getConstParticle(n);

      if(&us!=candidate) {
    
	float distance = candidate->distanceTo(us);
      

	if(distance<bestDistance || nearest==-1) {
	  bestDistance = distance;
	  nearest = n;
	}

      }

    }

    return nearest;

  }

  void SimSnap::diagnostics() {
    cerr << "SimSnap at " << this << "." << endl;
    cerr << "Total mass " << getTotalMass() << endl;
    cerr << "Number of particles " << getNumParticles() << endl;
    cerr << "Reported box size " << getBoxSize() << " / Calculated box size = " << getApparentBoxSize() << endl;
    cerr << "Cosmo: z = " << getRedshift() << " | h = " << getHubble() << endl;
    cerr << "Unit system: " << getDistanceUnits() << " | " << getMassUnits() << " | " <<  getVelocityUnits() << " | " <<  getDensityUnits() << " | " << getEnergyUnits() << endl;
    cerr << "Internal revision " << getVersion() << endl;
    if(pAutoParent==NULL)
      cerr << "(No parents)" << endl;
    else
      cerr << "Parent at " << pAutoParent << endl;
    cerr << getNumArrays() << " arrays" << endl;
    for(int n=0; n<getNumArrays(); n++) {
      cerr << "("<<n<<") " << getArrayName(n) << endl;
    }
    
  }

  int SimSnap::getNearestNeighbour(int from) const {

    const Particle *us = getConstParticle(from);

    int nearest = getNearestNeighbour(*us);
    

    return nearest;
    
  }


  // below is for sorting when performing nearest neighbour searches

  /*bool operator < (pair<int,double> a,pair<int,double> b) {
    return a.second<b.second;
    }
  */

  // nearest neighbour searches:

  bool compareSecond::operator()(pair<int,double> p ,pair<int,double> q) const
  {
    return q.second < p.second;
  }

  float SimSnap::localPPscale(float x, float y, float z) const {
    Particle point(x,y,z);
    list<pair<int,double> > nn = getNearestNeighbours(point, 2);
    return (*nn.begin()).second;
  }


  double SimSnap::getSPHColDen(float x, float y, const Unit & inunits) const {
    double colden_tot = 0;

    if(pSmooth==NULL)
      initialiseSPH(32);

    float x1,x2,y1,y2,z1,z2;
    getExactBoundaries(x1,x2,y1,y2,z1,z2);
    double colden_conv = (getMassUnits()/pow(getDistanceUnits(),2)).convertTo(inunits,this);

    double delta  = 0.;
    double accum = 0.;


    for(float z=z1;z<z2;z+=delta) {
      Particle point(x,y,z);
      list<pair<int,double > > ptcls = getNearestNeighbours(point,nPartSPH);
      delta = ((*(ptcls.begin())).second)/5.; 
      // N.B. is this too large?
   
     
      float rho=0;
      float smooth = ((*(ptcls.begin())).second)/2.; 
     
      list<pair<int,double > >::iterator i;
      const Particle *p;

      
      for(i = ptcls.begin();i!=ptcls.end();i++) {
	int n = (*i).first;
	// float smooth_i = getSmooth(n);
	float k = (*pKernel)((*i).second,smooth);
	p = getConstParticle(n);

	rho += p->mass * k;
      }
      

      /*
      // cheating version
      int n = (*(--ptcls.end())).first;
      
      rho = p->rho;
      */

      colden_tot+=rho*delta*colden_conv;

    }  

    return colden_tot;
    

  }

  list<pair<int,double> > SimSnap::getNearestNeighbours(const Particle &from, int number, const Metric &metric) const {
  
    // returns linked list of n nearest neighbours, 
    // starting with the FURTHEST
  
    // Note an optimised version of this routine
    // is available by constructing a Grid around
    // the SimSnap object, and calling its getNearestNeighbour
    // override. 
    //
    // However, this routine should be pretty efficient also.

    int nParts = getNumParticles();
    list< pair<int,double> > candidates;

    int n_elements = 0; 

    // better to keep track of length - more 
    // efficient than traversing linked list
    // every time we need to know the length! 
    //
    // (STL does not define its length measuring complexity)

    compareSecond comp; // comparing object for insertion sort (STL)

#ifdef SIMAN_OMP
#pragma omp parallel for default(shared)
#endif
    for(int n=0;n<nParts;n++) {
      const Particle *p = getConstParticle(n);
    
      double dist = metric(*p,from);
      if(dist<(*candidates.begin()).second || (n_elements<number-1)) {

#ifdef SIMAN_OMP
	#pragma omp critical 
	{
#endif

	// i.e. this is within the boundaries of interest, or it's
	//      one of the first n particles to consider

	// now do insertion sort!
    
	// 1. setup insertion list

	pair<int,double> add(n,dist);
	list< pair<int,double> > add_list;
	add_list.push_front(add);
   

	// 2. perform insertion
      
	candidates.merge(add_list,comp);
	n_elements++;
      

	if(n_elements>number) {
	  candidates.pop_front();
	
	  // can only be one over, only adding one at a time!
	
	  n_elements--;
	}

#ifdef SIMAN_OMP
	}
#endif

      }
    }


    return candidates;
  }


  void SimSnap::initialiseSPH(int nSPH, SPHKernel *pK) const {

    if(pK!=NULL)
      pKernel = pK;
    else
      pKernel = new SPHKernel();

    nPartSPH = nSPH;

    if(pSmooth!=NULL) free(pSmooth);

    int numPart = getNumParticles();
  
  
    pSmooth = (float*) malloc(sizeof(float)*numPart);
    memset(pSmooth,0,sizeof(float)*numPart);


  
  }


  float SimSnap::getSmooth(int n, list< pair<int,double> > *nnList ) const {

    if(pAutoParent!=NULL) return getSmooth(deReference(n,1),nnList);

    if(pSmooth==NULL) 
      initialiseSPH();

  
    if(pSmooth[n]==0) {
     
      // obtain a first estimate for the smoothing length.
      // this may not be suitable for smoe purposes - consistentSmooth
      // is the current way to get a better set of smoothing lengths,
      // operating using (6) from Springel 2005

      if(nnList!=NULL) {
	// provided with nearest neighbour list
	pSmooth[n] = ((*(nnList->begin())).second)/2.;
      } else {
	// need to find nearest neighbour list
	const Particle *p=getConstParticle(n);
	list< pair<int,double> > nearestNe = getNearestNeighbours(*p,nPartSPH);
	pSmooth[n] = ((*(nearestNe.begin())).second)/2.; 
      }
    }

    return pSmooth[n];
  }

  float SimSnap::consistentSmooth(int n, float tol) const {
    const Particle *p=getConstParticle(n);
    list< pair<int,double> > nearestNe = getNearestNeighbours(*p,nPartSPH);
    float mass=0;
    float delta=1,xdelta=2;

    do {
      xdelta = delta;
      mass =0;
      for(list< pair<int,double> >::iterator  i=nearestNe.begin(); i!=nearestNe.end(); i++) {
	const Particle *p2 = getConstParticle((*i).first);
	mass+=p2->mass;
      }
    
      // note that if no initial estimate for the smoothing length has been performed,
      // estimateDensity will do it automatically through calling getSmooth
      float rho = estimateDensity(n,&nearestNe);
    
      float h = std::pow(3*mass/(4*PI*rho),1./3.)/2.;
      delta = abs(h-pSmooth[n])/h;
      pSmooth[n]=h;
    
    } while (delta>tol && delta<xdelta && tol!=0);
    if(delta>xdelta) cerr << "SimSnap: consistentSmooth divergence" << endl;
    return delta;
  }

  void SimSnap::consistentSmooth(float tol) const {
    int nPart = getNumParticles();
    for(int n=0;n<nPart;n++) {
      consistentSmooth(n,tol);
    }

  }


  float SimSnap::estimateDensity(const float x, const float y, const float z) const {


    Particle point(x,y,z);

    list< pair<int,double> > ptcls = getNearestNeighbours(point,nPartSPH);
    list< pair<int,double> >::iterator i;

    float rho = 0.;
    float smooth = ((*(ptcls.begin())).second)/2.;

    for(i = ptcls.begin();i!=ptcls.end();i++) {
      int n = (*i).first;
      const Particle *p = getConstParticle(n);
      float mass = p->mass;
  
      rho += mass * (*pKernel)((*i).second,smooth);
      
    }

    return rho;

  }


  float SimSnap::getHeliumMassFrac() const {
    return constants::heliumY;
  }

  float SimSnap::estimateDensity(int pid, list< pair<int,double> > *nnList) const {

    const Particle *point = getConstParticle(pid);
  

    list< pair<int,double> > ptcls;

    if(nnList==NULL) {
      ptcls = getNearestNeighbours(*point,nPartSPH);
      nnList = &ptcls;
    }
    
    float smooth = getSmooth(pid,nnList);

    list< pair<int,double> >::iterator i;
  
    float rho = 0.;
  
    for(i = nnList->begin();i!=nnList->end();i++) {
      int n = (*i).first;
      const Particle *p = getConstParticle(n);
    
      rho += p->mass * (*pKernel)((*i).second,smooth);
      
    }
  

    return rho;
  
  
  }

  SimanVec SimSnap::angMom() const {
    double L_x;
    double L_y;
    double L_z;

    for(unsigned int n=0; n<getNumParticles(); n++) {
      const Particle *p=getConstParticle(n);
      L_x+=p->mass*(p->y*p->vz - p->z*p->vy);
      L_y+=p->mass*(p->z*p->vx - p->x*p->vz);
      L_z+=p->mass*(p->x*p->vy - p->y*p->vx);
    }

    vector<double> ret;
    ret.push_back(L_x);
    ret.push_back(L_y);
    ret.push_back(L_z);

    return ret;
  }

  vector<double> SimSnap::rotCurve(vector<double> R_bins, double scale_height) const {
    
    SimanVec am = angMom();
    vector<double> ret(R_bins.size()-1,0);
    vector<double> mass(R_bins.size()-1,0);

    for(unsigned int n=0; n<getNumParticles(); n++) {
      const Particle *p=getConstParticle(n);
      SimanVec R(p->x,p->y,p->z);
      R.proj(am);
      double rotv = SimanVec(p->vx,p->vy,p->vz).proj(am).proj(R).abs();

      double modR = R.abs();

      for(unsigned int i=0; i<R_bins.size()-1; i++) {
	if(modR>R_bins[i] && modR<R_bins[i+1]) {
	  ret[i]+=rotv*p->mass;
	  mass[i]+=p->mass;
	  break;
	}
      }
    }
    for(unsigned int i=0; i<R_bins.size()-1; i++) {
      ret[i]/=mass[i];
    }

    return ret;
     


  }


  float SimSnap::getApparentBoxSize() const {
    float x1,x2,y1,y2,z1,z2;
    getExactBoundaries(x1,x2,y1,y2,z1,z2);
    float w = x2-x1;
    if(y2-y1>w) w=y2-y1;
    if(z2-z1>w) w=z2-z1;
    return w;
  }

  float SimSnap::getMaximalVelocity() const {
    float vmax =0;
    for(unsigned int n=0; n<getNumParticles(); n++) {
      const Particle *p = getConstParticle(n);
      float v = p->vx*p->vx + p->vy*p->vy + p->vz*p->vz;
      if(v>vmax)
	vmax =v;
    }
    return sqrt(vmax);
  }

  float SimSnap::criticalDensity() const {

    float rho_crit_0 = 277.619737 * getHubble()*getHubble(); // Msol/kpc^3
    
    // Now correct for evolution of hubble const:
    float rho_crit = rho_crit_0 * (getOmegaM0() * std::pow((float)1.+getRedshift(),(float)3) + getOmegaLambda0());
    
    // Finally convert to sim units
    return rho_crit * Unit("Msol kpc^-3").convertTo(getDensityUnits(),this);
  }

  float SimSnap::criticalDensity0() const {
    float rho_crit_0 = 277.619737 * getHubble()*getHubble(); // Msol/kpc^3
    return rho_crit_0 * Unit("Msol kpc^-3").convertTo(getDensityUnits(),this);
  }

  float SimSnap::getVirialRadius(float delta) const throw(ConvergenceFailure) {
    // uses binary chop algorithm
    float r_max = getApparentBoxSize();
    float r_min = 0;
    float density_aim = delta*criticalDensity0()*std::pow((float)(1.+getRedshift()),(float)3.);
    
    float r_try = (r_min+r_max)*0.5;
    
    float closeness = numeric_limits<float>::max();
    float x_closeness=closeness;


    int iter=0;

    do {



      Subset s(const_cast<SimSnap*>(this),Sphere(r_try));
      float den = s.getTotalMass()/(4.*PI/3. * std::pow(r_try,3));

      if(getVerbose()>4)
	cerr << r_try << "\t" << r_min << "\t" << r_max << "\t" << s.getTotalMass() << "\t" << den << "\t" << density_aim << endl;
      if (den<density_aim) r_max=r_try;
      if (den>density_aim) r_min=r_try;
      x_closeness = closeness;
      closeness = abs(den-density_aim)/density_aim;
      r_try = (r_min+r_max)*0.5;
      iter++;
    } while(closeness>0.0001 && iter<3000);
    if(closeness>0.0001) throw ConvergenceFailure();
    return r_try;
  }

  void SimSnap::getExactBoundaries(float &x1, float &x2, float &y1, float &y2, float &z1, float &z2) const {
    x2 = y2 = z2 = -FLT_MAX;
    x1 = y1 = z1 = FLT_MAX;

    int numParticles = getNumParticles();

    for(int n=0; n<numParticles; n++) {
      const Particle *p=getConstParticle(n);

      if(p->x < x1) x1 = p->x;
      if(p->y < y1) y1 = p->y;
      if(p->z < z1) z1 = p->z;

      if(p->x > x2) x2 = p->x;
      if(p->y > y2) y2 = p->y;
      if(p->z > z2) z2 = p->z;
    

    }
  
  }

  void SimSnap::write(string filename, int filetype) const {
    if(filetype==gadget) {
      GadgetSnap::nativeWrite(this, filename);
    } else if(filetype==native) {
      BaseSimSnap::nativeWrite(this,filename);
    } else {
      cerr << "SimSnap: No write code available for this filetype" << endl;
      return;
    }

  }


  SimanArray & SimSnap::createArray(string name, string fullname, Unit arrunits) {
    return pAutoParent->createArray(name,fullname,arrunits);
  }

  void SimSnap::destroyArray(string name) {
    pAutoParent->destroyArray(name);
  }

  Unit SimSnap::getArrayUnits(string name) const {
    return pAutoParent->getArrayUnits(name);
  }

  int SimSnap::getArrayIndex(string name) const {
    return pAutoParent->getArrayIndex(name);
  }

  Unit SimSnap::getArrayUnits(int n) const {
    return pAutoParent->getArrayUnits(n);
  }

  string SimSnap::getArrayLongName(string shortname) const {
    return pAutoParent->getArrayLongName(shortname);
  }

  string SimSnap::getArrayLongName(int n) const {
    return pAutoParent->getArrayLongName(n);
  }


  string SimSnap::getArrayName(int n) const {
    return pAutoParent->getArrayName(n);
  }

  SimanArray & SimSnap::getArray(string name) {
    return pAutoParent->getArray(name);
  }

  SimanArray & SimSnap::getArray(int n) {
    return pAutoParent->getArray(n);
  }


  const SimanArray & SimSnap::getConstArray(string name) const {
    return pAutoParent->getConstArray(name);
  }

  const SimanArray & SimSnap::getConstArray(int n) const {
    return pAutoParent->getConstArray(n);
  }

  int SimSnap::getNumArrays() const {
    return pAutoParent->getNumArrays();
  }

  /// VERSIONING
  
  long SimSnap::getVersion() const {
    if(pAutoParent!=NULL)
      return pAutoParent->getVersion();
    else
      return version;
  }

  void SimSnap::bumpVersion() {
    // assumed to only be called from "real" data-holding SimSnaps, to avoid
    // overhead of "if".
    version++;
  }

  void SimSnap::addHubbleFlow() {
    double dhub = Unit("100 h km s^-1 Mpc^-1").convertTo(getVelocityUnits()/getDistanceUnits(),this);
    for(unsigned int n=0; n<getNumParticles(); n++) {
      Particle &p(getParticleRef(n));
      p.vx+=p.x*dhub;
      p.vy+=p.y*dhub;
      p.vz+=p.z*dhub;
    }
  }

  void SimSnap::subtractHubbleFlow() {
    double dhub = Unit("100 h km s^-1 Mpc^-1").convertTo(getVelocityUnits()/getDistanceUnits(),this);
    for(unsigned int n=0; n<getNumParticles(); n++) {
      Particle &p(getParticleRef(n));
      p.vx-=p.x*dhub;
      p.vy-=p.y*dhub;
      p.vz-=p.z*dhub;
    }
  }

  vector<double> SimSnap::makeNaiveAbsProfile(int tau_size, double min_lambda, double max_lambda, double lambda_cen, double gamma, double osc_strength) const {
    float x1,x2,y1,y2,z1,z2;
    getExactBoundaries(x1,x2,y1,y2,z1,z2);
    float area = (x2-x1)*(y2-y1);
    double c = Unit("c").convertTo(getVelocityUnits());
    // double delta_lambda = (max_lambda-min_lambda)/tau_size;
    double colden_conv = (getMassUnits()/(getDistanceUnits()*getDistanceUnits())).convertTo(Unit("m_p cm^-2"));
    double lorentz_width = gamma * (lambda_cen/1.e5) * (lambda_cen/1.e5) / (4 * 3.141592 * 3.e8); // in angstroms
    vector<double> tau(tau_size);

    for(unsigned int n=0; n<getNumParticles(); n++) {
      const Particle &p(getConstParticleRef(n));
      if(p.type==Particle::gas) {
	double our_lambda = lambda_cen + lambda_cen * p.vz / c;
	double colden = p.mass/area * colden_conv;
	double gaussian_width = lambda_cen * sqrt(p.temp * constants::boltzmannInErgPerK / (constants::protonMassInG*1.e10) ) / 2.99e5;
	vector<double> add = siman::voigt(tau_size, min_lambda, max_lambda, our_lambda, gaussian_width, lorentz_width);
	for(unsigned int i=0; i<tau_size; i++) {
	  tau[i]+=double(8.87e-21)*double(lambda_cen*lambda_cen)*double(osc_strength*colden)*add[i];
	  // Constant: see 28/06 entry in blue book
	  // e^2 lambda^2 f / (4 epsilon_0 * m_e c^2)  =  8.87e-15 lambda^2 in SI   [m]^3
	  //                                           =  8.87e-5 once phi is in [angstroms]^-1
	  //                                           =  8.87e-25 once lambda_cen is in angstroms
	  //                                           =  8.87e-21 for result in cm^2
	}
      }
    }
    return tau;
  }

  
  vector<double> SimSnap::makeSPHAbsProfile(float x, float y, int tau_size, double min_lambda, double max_lambda, double lambda_cen, double gamma, double osc_strength) const {

    double colden_tot = 0;
    double colden_met = 0;
    double cloud_den = 0; // accumulate at least this much column density before sampling velocity field (probably should always be zero)
    if(pSmooth==NULL)
      initialiseSPH(10);

    float x1,x2,y1,y2,z1,z2;
    vector<double> tau(tau_size+2);
    getExactBoundaries(x1,x2,y1,y2,z1,z2);
    double colden_conv = (getMassUnits()/pow(getDistanceUnits(),2)).convertTo(Unit("m_p cm^-2"));
    double lorentz_width = gamma * (lambda_cen/1.e5) * (lambda_cen/1.e5) / (4 * 3.141592 * 3.e8); // in angstroms
    double c = Unit("c").convertTo(getVelocityUnits());
    double delta  = 0.;
    double accum = 0.;

    for(float z=z1;z<z2;z+=delta) {
      Particle point(x,y,z);
      list<pair<int,double > > ptcls = getNearestNeighbours(point,nPartSPH);
      delta = ((*(ptcls.begin())).second)/5.; 
      // N.B. is this too large?
   
      float vz_mean = 0.;
      float vz2_mean = 0.;
      float rho=0;
      float temp_mean=0;
      float metal_mean = 0;
      float mass_tot =0 ;

      float smooth = ((*(ptcls.begin())).second)/2.; 
      float norm = ((4.*PI/3.) * smooth*smooth*smooth)/float(nPartSPH);
      
      list<pair<int,double > >::iterator i;
      const Particle *p;

      /*
      for(i = ptcls.begin();i!=ptcls.end();i++) {
	int n = (*i).first;
	// float smooth_i = getSmooth(n);
	float k = (*pKernel)((*i).second,smooth);
	p = getConstParticle(n);

	rho += p->mass * k;
	mass_tot+=p->mass;
	vz_mean += p->vz * p->mass;
	vz2_mean += p->vz * p->vz * p->mass; 
	temp_mean  += p->temp * p->mass;
	metal_mean +=p->metal * p->mass; 
      }
      vz_mean/=mass_tot;
      vz2_mean/=mass_tot;
      temp_mean/=mass_tot;
      metal_mean/=mass_tot;
      */


      // cheating version
      int n = (*(--ptcls.end())).first;
      p = getConstParticle(n);
      rho = p->rho;
      vz_mean = p->vz;
      vz2_mean = p->vz*p->vz;
      temp_mean = p->temp;
      metal_mean = p->metal;

      float veldisp = lambda_cen * sqrt(vz2_mean - vz_mean*vz_mean) / c;
      
      double our_lambda = lambda_cen + lambda_cen * vz_mean / c;
      double colden = rho*delta * colden_conv;
      double gaussian_width = lambda_cen * sqrt(temp_mean * constants::boltzmannInErgPerK / (constants::protonMassInG*1.e10) ) / 2.99e5;
      // gaussian_width = sqrt(gaussian_width*gaussian_width + veldisp*veldisp); // turbulent and thermal components add in quadrature

      colden_tot+=colden;
      colden_met+=metal_mean*colden/0.76;
      /*
      accum+=colden;
      if(accum>cloud_den) {*/
	vector<double> add = siman::voigt(tau_size, min_lambda, max_lambda, our_lambda, gaussian_width, lorentz_width);
	for(unsigned int i=0; i<tau_size; i++) {
	  tau[i]+=double(8.87e-21)*double(lambda_cen*lambda_cen)*double(osc_strength*colden)*add[i];
	  // Constant: see naive, above
	  int idex = ((our_lambda-min_lambda)/(max_lambda-min_lambda))*tau_size;
	  //	if(idex==i)
	  // tau[i]+=delta*p->rho;
	}

	//	cerr << cloud_den << " " << accum << " " << z << " " << delta << " " << vz_mean << " " << rho << " " << colden << " " << veldisp << " " << gaussian_width << " " << ptcls.size() << endl;
	accum = 0.;
	/* } */
    }  
    tau[tau_size] = colden_tot;
    tau[tau_size+1] = colden_met;
    return tau;
    

  }



} // namespace siman
