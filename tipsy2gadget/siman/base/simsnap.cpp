//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"
#include "endian.hpp"

using namespace std;

CSimSnap::CSimSnap() {
  pAutoParent=NULL;
  nPartSPH = 6; // number of particles to use in SPH-like calculations
  pSmooth = NULL;
  pKernel=NULL;
}

CSimSnap::~CSimSnap() {
  if(pSmooth!=NULL) free(pSmooth);
  if(pKernel!=NULL) delete(pKernel);
}

// AUTO-PARENTING FUNCTIONS
//
// These have no meaning in the base class, but derived classes can set pAutoParent
// to save writing overrides to pass functions on to parent files (e.g. see CSubset)

string CSimSnap::className() {
  return "CSimSnap";
}

unsigned int CSimSnap::supports() {
  return SimSnap;
}

CSimSnap* CSimSnap::getParent() {
  return pAutoParent;
}

float CSimSnap::getOmegaM0() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getOmegaM0()";
    return 0.3;
  } else {
    return pAutoParent->getOmegaM0();
  }
}

float CSimSnap::getOmegaLambda0() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getOmegaLambda0()";
    return 0.7;
  } else {
    return pAutoParent->getOmegaLambda0();
  }
}

int CSimSnap::deReference(int i, int n) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method deReference();";
    return 0;
  } else {
    return pAutoParent->deReference(i,n);
  }
}


/*
void CSimSnap::makeUnitsPhysical() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method makeUnitsPhysical()\n";
  } else {
    pAutoParent->makeUnitsPhysical();
  }
}


bool CSimSnap::areUnitsPhysical() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method areUnitsPhysical()\n";
    return false;
  } else {
    return pAutoParent->areUnitsPhysical();
  }
}
*/

void CSimSnap::convertUnits(units::CUnit distance, units::CUnit mass, units::CUnit velocity, units::CUnit density, units::CUnit energy) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method convertUnits" << endl;
    return;
  } else {
    pAutoParent->convertUnits(distance,mass,velocity,density,energy);
  }
}

void CSimSnap::TempFromU() {

  
  cerr << "CSimSnap: calculating temperatures from U and Ne..."; 

  units::CUnit unit_spec_en_cgs= units::CUnit(units::energy_erg)/units::CUnit(units::mass_g);
  units::CUnit unit_spec_en_int = getEnergyUnits();
  
  double en_conv = unit_spec_en_int.convertTo(unit_spec_en_cgs);

  cerr << unit_spec_en_int << " - enToErgsPerG=" << en_conv << endl;
 
  const float Y = 0.26;
  const int numParticles = getNumParticles();

  for(int i=0;i<numParticles;i++) {
    
    CParticle *p = getParticle(i);


    if(p->type==CParticle::gas)  {
      
      double MeanWeight = 1/(1-0.75*Y+p->ne*(1-Y)) * units::protonMassInG;
      
      // internal energy -> cgs units
      
      float u  = p->u * en_conv;
      // = UnitEnergy_in_cgs/ UnitMass_in_g;
      
      float gamma= 5.0/3;
      
      // temperature in K:
      
      p->temp = MeanWeight/units::boltzmannInErgPerK * (gamma-1) * u;
      
    }
    
  } // for i

  cerr << "done!" << endl;
}

void CSimSnap::UFromTemp() {
   
  cerr << "CSimSnap: calculating U from temperatures and Ne...";

  units::CUnit unit_spec_en_cgs= units::CUnit(units::energy_erg)/units::CUnit(units::mass_g);
  units::CUnit unit_spec_en_int = getEnergyUnits();
  
  double en_conv = unit_spec_en_int.convertTo(unit_spec_en_cgs);
  
  const float Y = 0.26;
  const int numParticles = getNumParticles();

  for(int i=0;i<numParticles;i++) {
    
    CParticle *p = getParticle(i);
    if(p->type==CParticle::gas)  {
      
      double MeanWeight = 1/(1-0.75*Y+p->ne*(1-Y)) * units::protonMassInG;
      
      float gamma= 5.0/3;
      
      // temperature in K:
      
      
      float u = (p->temp * units::boltzmannInErgPerK)/(MeanWeight*(gamma-1.));
      p->u = u/ en_conv;
      
      // cerr << MeanWeight << " " << p->temp << " " << en_conv << "\t" << p->u << endl;
    }
    
  } // for i

  cerr << "done!" << endl;

}

int CSimSnap::getNumParticles() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getNumParticles()\n";
    return 0;
  } else {
    return pAutoParent->getNumParticles();
  }
}


units::CUnit CSimSnap::getDistanceUnits() {
  if(pAutoParent==NULL) {
  
    return lenUnits;
  } else {
    return pAutoParent->getDistanceUnits();
  }
}

units::CUnit CSimSnap::getVelocityUnits() {
  if(pAutoParent==NULL) {
 
    return velUnits;
  } else {
    return pAutoParent->getVelocityUnits();
  }
}


units::CUnit CSimSnap::getDensityUnits() {
  if(pAutoParent==NULL) {
   
    return denUnits;
  } else {
    return pAutoParent->getDensityUnits();
  }
}


units::CUnit CSimSnap::getEnergyUnits() {
  if(pAutoParent==NULL) {
    
    return enUnits;
  } else {
    return pAutoParent->getEnergyUnits();
  }
}


units::CUnit CSimSnap::getMassUnits() {
  if(pAutoParent==NULL) {
  
    return massUnits;
  } else {
    return pAutoParent->getMassUnits();
  }
}


float CSimSnap::getBoxSize() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getBoxSize()\n";
    return 0;
  } else {
    return pAutoParent->getBoxSize();
  }

}

float CSimSnap::getHubble() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getHubble()\n";
    return 0;
  } else {
    return pAutoParent->getHubble();
  }

}

void CSimSnap::setHubble(float h) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method setHubble()\n";
    
  } else {
    pAutoParent->setHubble(h);
  }

}


void CSimSnap::setBoxSize(float h) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method setBoxSize()\n";
    
  } else {
    pAutoParent->setBoxSize(h);
  }

}


void CSimSnap::setRedshift(float z) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method setRedshift()" << endl;
    
  } else {
    pAutoParent->setRedshift(z);
  }

}

float CSimSnap::getRedshift() {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getRedshift()\n";
    return 0;
  } else {
    return pAutoParent->getRedshift();
  }
}

void CSimSnap::pushParticle(int n) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method pushParticle()\n";
  } else {
    pAutoParent->pushParticle(n);
  }
}


CParticle * CSimSnap::getParticle(int id) {
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method getParticle()\n";
    return 0;
  } else {
    return pAutoParent->getParticle(id);
  }

}

void CSimSnap::releaseParticle(CParticle *p) {
  
  if(pAutoParent==NULL) {
    cerr << "CSimSnap: Base class erroneously called for virtual method releaseParticle()\n";
  } else {
    pAutoParent->releaseParticle(p);
  }

  
}


float CSimSnap::getTotalMass() {

  float mass = 0;
  CParticle *particle;

  for(int n=0;n<getNumParticles();n++) {
    particle = getParticle(n);
    mass +=particle->mass;
   
  }

  return mass;
}


float CSimSnap::getH0Mass() {

  float mass = 0;
  CParticle *particle;
  
  const float Yhe = getHeliumMassFrac();

  for(int n=0;n<getNumParticles();n++) {
    particle = getParticle(n);
    mass +=particle->mass*(1-Yhe)*(1-particle->nHp);
 
  }

  return mass;
}


bool CSimSnap::references(CSimanObject *p) {
  if(p==pAutoParent) return true;
  return CSimanObject::references(p);
}

void CSimSnap::centreOfMass(float &cx, float &cy, float &cz) {
  cx = 0.;
  cy = 0.;
  cz = 0.;

  double cxd=0., cyd=0., czd=0.;
  
  double totMass = getTotalMass();
  
  int n=0;
  int numParticles = getNumParticles();

  double meanMass = totMass/(double)numParticles;

  

  CParticle *particle;

  for(n=0;n<numParticles;n++) {
    particle = getParticle(n);

    cxd+=(double)particle->x*(double)particle->mass/totMass;
    cyd+=(double)particle->y*(double)particle->mass/totMass;
    czd+=(double)particle->z*(double)particle->mass/totMass;

    releaseParticle(particle);
  }

  // cxd/=(double)numParticles;
  //cyd/=(double)numParticles;
  //czd/=(double)numParticles;

  cx = (float)cxd;
  cy = (float)cyd;
  cz = (float)czd;
 
}


void CSimSnap::centreOfMassVel(float &cvx, float &cvy, float &cvz) {

  cvx = 0.;
  cvy = 0.;
  cvz = 0.;

  double cxd=0., cyd=0., czd=0.;
  
  double totMass = getTotalMass();
  
  int n=0;
  int numParticles = getNumParticles();

  double meanMass = totMass/(double)numParticles;

  

  CParticle *particle;

  for(n=0;n<numParticles;n++) {
    particle = getParticle(n);

    cxd+=(double)particle->vx*(double)particle->mass/totMass;
    cyd+=(double)particle->vy*(double)particle->mass/totMass;
    czd+=(double)particle->vz*(double)particle->mass/totMass;

    releaseParticle(particle);
  }

  // cxd/=(double)numParticles;
  //cyd/=(double)numParticles;
  //czd/=(double)numParticles;

  cvx = (float)cxd;
  cvy = (float)cyd;
  cvz = (float)czd;
 
}


void CSimSnap::shrinkSphereCentre(float &cx, float &cy, float &cz, float shrinkFactor, int minParticles, float initial_radius) {

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

  if(getNumParticles()<=minParticles)
    return;

 
  float r = 0.;

  if(initial_radius>0.)  
    r = initial_radius*shrinkFactor; 
  else 
    r = getBoxSize()*shrinkFactor;


  CSphere sphereFilter(cx,cy,cz,r);  
  CSubset sphereParticles(this,sphereFilter);


  sphereParticles.shrinkSphereCentre(cx,cy,cz,shrinkFactor,minParticles,r);

  // recursion from here will automatically return our required cx,cy,cz
  
  
}

float CSimSnap::power(float lambda) {
  complex <float> sum;

  unsigned int numParticles = this->getNumParticles();
 
  for(unsigned int n=0; n<numParticles; n++) {
    CParticle *p=getParticle(n);
    complex <float> exponent(0,lambda*p->x);
    sum += exp(exponent);
    releaseParticle(p);
  }
    
    
  return pow(abs(sum),2);
}



int CSimSnap::determineType(std::string filename) {


  ENDINIT;

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

  END4(&test_word);

  if(test_word==sizeof(gadget_header)) 
    return gadget + endianWrong;

  
  if(tw[0]=='#' && tw[1]=='!')
    return script;

  END4(&test_word);



  tipsy_header test_header;

  file.seekg(0,ios_base::beg);
  
  file.read((char*) &test_header,sizeof(test_header));

  if(test_header.nsph + test_header.ndark + test_header.nstar == test_header.nbodies)
    return tipsy;

  END4(&test_header.nsph);
  END4(&test_header.ndark);
  END4(&test_header.nstar);
  END4(&test_header.nbodies);

  
  if(test_header.nsph + test_header.ndark + test_header.nstar == test_header.nbodies)
    return tipsy + endianWrong;


  // don't know what this is

  return unknown;


}

CSimSnap * CSimSnap::loadFile(std::string filename) {
  int type = determineType(filename);

  if((type & gadget) > 0) {
    CGadgetFile *pGadg;
    cerr << "CSimSnap: trying to load " << filename << " as a gadget file" << endl;
    pGadg = new CGadgetFile(filename.c_str());
    
    return pGadg;
  }

  if((type & tipsy) > 0) {
    CTipsyFile *pTipsy;
    cerr << "CSimSnap: trying to load " << filename << " as a tipsy-binary file" << endl;
    pTipsy = new CTipsyFile(filename.c_str());
    
    return pTipsy;
  }

  if((type & native) > 0) {
    CBaseSimSnap *pNative;
    cerr << "CSimSnap: trying to load " << filename << " as a SimAn file" << endl;
    pNative = new CBaseSimSnap(filename);
    return pNative;
  }

  if((type & script) > 0) {
    CScripted *pScript;
    pScript = new CScripted(filename);
    if(pScript->getReturnValue()->supports(CSimanObject::SimSnap))
      return (CSimSnap*) pScript->getReturnValue();
  }

  // fallthrough

  cerr << "CSimSnap: don't know how to load file "<<filename<< " - Returning NULL pointer." << endl;

  return NULL;
}

int CSimSnap::getNearestNeighbour(const CParticle &us) {
  
  float bestDistance(0.), distance(0.);

  int nParts = getNumParticles();
  int nearest = -1;

  for(int n=0;n<nParts;n++) {
    
    CParticle *candidate = getParticle(n);

    if(&us!=candidate) {
    
      float distance = candidate->distanceTo(us);
      

      if(distance<bestDistance || nearest==-1) {
	bestDistance = distance;
	nearest = n;
      }

    }
    releaseParticle(candidate);
  }

  return nearest;

}

int CSimSnap::getNearestNeighbour(int from) {

  CParticle *us = getParticle(from);

  int nearest = getNearestNeighbour(*us);
  releaseParticle(us);

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

float CSimSnap::localPPscale(float x, float y, float z) {
  CParticle point(x,y,z);
  list<pair<int,double> > nn = getNearestNeighbours(point, 2);
  return (*nn.begin()).second;
}

float CSimSnap::getSPHColDen(float x, float y, float z1, float z2, double units, bool neutral) {
  
  
  float delta = 0;
  double colden=0;

  CMetric2D metric;

  list<pair<int,double> > ptcls = getNearestNeighbours(CParticle(x,y,z1),nPartSPH,metric);

  list<pair<int,double> >::iterator i;
  

  CSPHProjectionKernel Kernel2D(*pKernel);
  CParticle cenpix(x,y,0);
  const float Yhe = getHeliumMassFrac();
  float smooth = ((*(ptcls.begin())).second)/2.;
  
  for(i=ptcls.begin();i!=ptcls.end();i++) {
    CParticle *p = getParticle((*i).first);
    if(neutral) {
      colden+= p->mass*(1-Yhe)*(1-p->nHp) * Kernel2D((*i).second,smooth);
    } else {
      colden+=p->mass * Kernel2D((*i).second,smooth);
    }
    // cout << p->x << "\t" << p->y << "\t" << p->z << "\t" << (*i).second << "\t" << Kernel2D(metric(*p,cenpix),smooth) << endl;
  }

  
 
  /* A more obvious, but much slower, approach:

  if(neutral) {
    for(float z=z1;z<z2;z+=delta) {
      CParticle point(x,y,z);
      list<pair<int,double > > ptcls = getNearestNeighbours(point,nPartSPH,metric);
      delta = ((*(ptcls.begin())).second)/2.; 
      // N.B. is this too large?
      colden+=(delta*estimateDensityH0(x,y,z,&ptcls))*units;
      // cout << "-" << estimateDensityH0(x,y,z,&ptcls) << endl;
    }
  } else {
    for(float z=z1;z<z2;z+=delta) {
      colden+=(delta*estimateDensity(x,y,z))*units;
      delta = localPPscale(x,y,z+delta);
    }
  }
  */

  if(colden>1.) {
    return log(colden)/log(10.);
  } else {
    return 0.;
  }
}

list<pair<int,double> > CSimSnap::getNearestNeighbours(const CParticle &from, int number, const CMetric &metric) {
  
  // returns linked list of n nearest neighbours, 
  // starting with the FURTHEST
  
  // Note an optimised version of this routine
  // is available by constructing a CGrid around
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


  for(int n=0;n<nParts;n++) {
    CParticle *p = getParticle(n);
    
    double dist = metric(*p,from);
    
    if(dist<(*candidates.begin()).second || (n_elements<number-1)) {

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

      releaseParticle(p);
    }
  }

  return candidates;
}


void CSimSnap::initialiseSPH(int nSPH, CSPHKernel *pK) {

  if(pK!=NULL)
    pKernel = pK;
  else
    pKernel = new CSPHKernel();

  nPartSPH = nSPH;

  if(pSmooth!=NULL) free(pSmooth);

  int numPart = getNumParticles();
  
  cerr << "CSimSnap: allocating memory for SPH calculations...";
  
  pSmooth = (float*) malloc(sizeof(float)*numPart);
  memset(pSmooth,0,sizeof(float)*numPart);

  cerr << "done!" << endl;

  
}


float CSimSnap::getSmooth(int n, list< pair<int,double> > *nnList ) {

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
      CParticle *p=getParticle(n);
      list< pair<int,double> > nearestNe = getNearestNeighbours(*p,nPartSPH);
      pSmooth[n] = ((*(nearestNe.begin())).second)/2.; 
    }
  }

  return pSmooth[n];
}

float CSimSnap::consistentSmooth(int n, float tol) {
  CParticle *p=getParticle(n);
  list< pair<int,double> > nearestNe = getNearestNeighbours(*p,nPartSPH);
  float mass=0;
  float delta=1,xdelta=2;

  do {
    xdelta = delta;
    mass =0;
    for(list< pair<int,double> >::iterator  i=nearestNe.begin(); i!=nearestNe.end(); i++) {
      CParticle *p2 = getParticle((*i).first);
      mass+=p2->mass;
    }
    
    // note that if no initial estimate for the smoothing length has been performed,
    // estimateDensity will do it automatically through calling getSmooth
    float rho = estimateDensity(n,&nearestNe);
    
    float h = pow(3*mass/(4*PI*rho),1./3.)/2.;
    delta = abs(h-pSmooth[n])/h;
    pSmooth[n]=h;
    
  } while (delta>tol && delta<xdelta && tol!=0);
  if(delta>xdelta) cerr << "CSimSnap: consistentSmooth divergence" << endl;
   return delta;
}

void CSimSnap::consistentSmooth(float tol) {
  int nPart = getNumParticles();
  for(int n=0;n<nPart;n++) {
    float delta =  consistentSmooth(n,tol);

  }

}


float CSimSnap::estimateDensity(const float x, const float y, const float z) {


  CParticle point(x,y,z);

  list< pair<int,double> > ptcls = getNearestNeighbours(point,nPartSPH);
  list< pair<int,double> >::iterator i;

  float rho = 0.;
  float smooth = ((*(ptcls.begin())).second)/2.;

  for(i = ptcls.begin();i!=ptcls.end();i++) {
    int n = (*i).first;
    CParticle *p = getParticle(n);
    float mass = p->mass;
  
    rho += mass * (*pKernel)((*i).second,smooth);
    releaseParticle(p);
  }

  return rho;

}


float CSimSnap::getHeliumMassFrac() {
  return 0.24;
}

float CSimSnap::estimateDensityH0(const float x, const float y, const float z, list< pair<int,double> > *nnList) {

  const float Yhe = getHeliumMassFrac();

  CParticle point(x,y,z);
  
  list< pair<int,double> > ptcls;

  if(nnList==NULL) {
    ptcls = getNearestNeighbours(point,nPartSPH);
    nnList = &ptcls;
  }
  
  float smooth = ((*(nnList->begin())).second)/2.;

  list< pair<int,double> >::iterator i;

  float rho = 0.;

  for(i = nnList->begin();i!=nnList->end();i++) {
    int n = (*i).first;
    CParticle *p = getParticle(n);
    float mass = p->mass*(1-p->nHp)*(1-Yhe);
  
    rho += mass * (*pKernel)((*i).second,smooth);
    releaseParticle(p);
  }

  return rho;

}

float CSimSnap::estimateDensity(int pid, list< pair<int,double> > *nnList) {

  CParticle *point = getParticle(pid);
  

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
    CParticle *p = getParticle(n);
    
    rho += p->mass * (*pKernel)((*i).second,smooth);
    releaseParticle(p);
  }
  
  releaseParticle(point);
  return rho;
  
  
}

float CSimSnap::getApparentBoxSize() {
  float x1,x2,y1,y2,z1,z2;
  getExactBoundaries(x1,x2,y1,y2,z1,z2);
  float w = x2-x1;
  if(y2-y1>w) w=y2-y1;
  if(z2-z1>w) w=z2-z1;
  return w;
}

void CSimSnap::getExactBoundaries(float &x1, float &x2, float &y1, float &y2, float &z1, float &z2) {
  x2 = y2 = z2 = -FLT_MAX;
  x1 = y1 = z1 = FLT_MAX;

  int numParticles = getNumParticles();

  for(int n=0; n<numParticles; n++) {
    CParticle *p=getParticle(n);

    if(p->x < x1) x1 = p->x;
    if(p->y < y1) y1 = p->y;
    if(p->z < z1) z1 = p->z;

    if(p->x > x2) x2 = p->x;
    if(p->y > y2) y2 = p->y;
    if(p->z > z2) z2 = p->z;
    
    releaseParticle(p);
  }
  
}

void CSimSnap::write(string filename, int filetype) {
  if(filetype==gadget) {
    CGadgetFile::nativeWrite(this, filename);
  } else if(filetype==native) {
    CBaseSimSnap::nativeWrite(this,filename);
  } else {
    cerr << "CSimSnap: No write code available for this filetype" << endl;
    return;
  }

  string filename_units = filename+".units";
  
  ofstream file_units(filename_units.c_str());
  file_units << getDistanceUnits() << " " << getMassUnits() << " " << getVelocityUnits() << " " << getDensityUnits() << " " << getEnergyUnits() << endl;
}


float * CSimSnap::createArray(string name, string fullname, units::CUnit arrunits) {
  return pAutoParent->createArray(name,fullname,arrunits);
}

void CSimSnap::destroyArray(string name) {
  pAutoParent->destroyArray(name);
}

units::CUnit CSimSnap::getArrayUnits(string name) {
  return pAutoParent->getArrayUnits(name);
}

units::CUnit CSimSnap::getArrayUnits(int n) {
  return pAutoParent->getArrayUnits(n);
}

string CSimSnap::getArrayLongName(string shortname) {
  return pAutoParent->getArrayLongName(shortname);
}

string CSimSnap::getArrayLongName(int n) {
  return pAutoParent->getArrayLongName(n);
}


string CSimSnap::getArrayName(int n) {
  return pAutoParent->getArrayName(n);
}


float * CSimSnap::getArray(string name) {
  return pAutoParent->getArray(name);
}

float * CSimSnap::getArray(int n) {
  return pAutoParent->getArray(n);
}

int CSimSnap::getNumArrays() {
  return pAutoParent->getNumArrays();
}

CSimanObject * CSimSnap::dispatch(string command, istream *stream, CScripted *pS) {
 
  
  if(command=="simplify") {
    float prob;
   
    *stream >> prob;
    CRandomFilter filter(prob);
    CSubset *pNew = new CSubset(this,filter);
    
    return pNew;
  } else if(command=="denunits") {
   
    cerr << getDensityUnits();
    return NULL;
  } else if(command=="enunits") {
    cerr << getEnergyUnits();
    return NULL;
  } else if(command=="massunits") {
    cerr << getMassUnits();
    return NULL;
  } else if(command=="velunits") {
    cerr << getVelocityUnits(); 
    return NULL;
  } else if(command=="distunits") {
    cerr << getDistanceUnits();
    return NULL;
  } else if(command == "select") {
    int mask = 0;
    while(!stream->eof()) {
      string type;
      *stream >> type;
      transform(type.begin(),type.end(),type.begin(),(int(*)(int))tolower);
      if(type=="gas") {
	mask+=CParticle::gas;
      } else if(type=="star" || type=="stars") {
	mask+=CParticle::star;
      } else if(type=="dm" || type=="dark") {
	mask+=CParticle::dm;
      } else {
	cerr << "CScripted: Unknown particle type '" << type << "'" << endl;
      }
    }
    CParticleTypeFilter desc(mask);
    CSubset* pNew = new CSubset(this,desc);
    return pNew;
  } else if(command == "physical") {
    convertUnits();
    return NULL;
  } else if(command == "sphere") {
    float x,y,z,r;
    *stream >> x >> y >> z >> r;
    CSphere sphereDesc(x,y,z,r);
    CSubset *pNew = new CSubset(this,sphereDesc);
    return pNew;
  } else if(command=="transform") {
    CGeometry geom(this);
    float mat[9];
    for(int n=0; n<9; n++) {
      *stream >> mat[n];
    }
    geom.setMatrix(mat);
    geom.apply();
    return NULL;
  } else if (command=="recentre") {
    
    string method;
    float x,y,z;
    *stream >> method;
    if(method=="shrink") {
      shrinkSphereCentre(x,y,z);
    } else if(method=="on") {
      *stream >> x >> y >> z;
    } else {
      centreOfMass(x,y,z);
    }
    
    CGeometry geom(this);
    geom.reCentre(x,y,z);
    geom.apply();
    return NULL;
  } else if (command=="save") {
    string type;
    *stream >> type;
    unsigned int type_id;
    if(type=="gadget") {
      type_id=gadget;
    } 
    if(type=="tipsy") {
      type_id=tipsy;
    }
    if(type=="native") {
      type_id=native;
    }
    string fname;
    if(type_id!=0) {
      *stream >> fname;
    } else {
      type_id=native;
      fname=type;
    }
    if(fname!="") {
      write(fname,type_id);
    } else {
      CSyntaxError n;
      throw(n);
    }
    return NULL;
  } else if (command=="z") {
    if(!stream->eof()) {
      float z;
      *stream >> z;
      setRedshift(z);
    }
    cerr << "Redshift " << getRedshift() << endl;
    return NULL;
  }

  
  return CSimanObject::dispatch(command,stream,pS);
  
}
