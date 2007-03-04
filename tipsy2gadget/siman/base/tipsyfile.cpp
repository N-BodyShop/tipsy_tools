//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// TIPSYFILE.CPP

#include "siman.hpp"
#include "xdrfile.hpp"

using namespace std;

// CONSTRUCTOR / DESTRUCTORS

CTipsyFile::~CTipsyFile() {


}

CTipsyFile::CTipsyFile(const char *fname_in) {

  strcpy(fname,fname_in);
  load();
}



// LOAD TIPSY FILE

bool CTipsyFile::load() {
  
  bool xdrdata = false;

  // returns true for success

  ifstream file(fname,ios::binary);

  if(file.is_open()) {

    // READ HEADER
      // Note that this isn't perfectly portable.  E.g. if a system
      // happens to have 64 bit ints, this will fail.

    file.read((char*) &header,sizeof(header));
    
    if(header.nsph + header.ndark + header.nstar != header.nbodies || header.ndim>100) {
	XDR xdrs;
	xdrmem_create(&xdrs, reinterpret_cast<char *>(&header), sizeof(header),
		      XDR_DECODE);
	xdr_template(&xdrs, &header);
	xdr_destroy(&xdrs);
    
      if(header.nsph + header.ndark + header.nstar == header.nbodies && header.ndim<100) {
	cerr << "CTipsyFile: Using Standard tipsy (XDR) format, auto-flipping (performance warning)" << endl;
	xdrdata = true;
       	
      } else {
	cerr << "CTipsyFile: Don't recognise this file format" << endl;
	throw(CFileError(fname));
      }

    }

    // some file have a mysterious "zero" word, others don't
    // if this one doesn't, seek back to counteract reading it

    if(header.zero!=0)
      file.seekg(-4,ios_base::cur);

    if(header.ndim!=3) {
      cerr << "CTipsyFile: " << header.ndim << " dimensions not supported" << endl;
      return false;
    }

    cerr << "CTipsyFile: time " << header.time << "; G: " << header.nsph << " D: " << header.ndark << " S: " << header.nstar << endl;


    // ALLOCATE MEMORY

    numParticles = header.nsph+header.ndark+header.nstar;

    allocateMemory();


    // READ PARTICLE DATA

    cerr << "CTipsyFile: reading...";

    unsigned int i=0;

    tipsy_gas_particle t_gp;
    tipsy_dark_particle t_dp;
    tipsy_star_particle t_sp;

    float *tform=NULL;
 
    if(header.nstar>0)
      tform = createArray("tform","TIPSY tform");

    float *softening = NULL;
    softening = createArray("softening", "Gravitational softening");
    
    float *potential = NULL;
    potential = createArray("potential", "Gravitational potential");
    
    for(unsigned int n=header.nsph;n!=0;--n) {
      file.read((char*) &t_gp,sizeof(tipsy_gas_particle));
      if(xdrdata) {
	  XDR xdrs;
	  xdrmem_create(&xdrs, reinterpret_cast<char *>(&t_gp),
			sizeof(tipsy_gas_particle), XDR_DECODE);
	  xdr_template(&xdrs, &t_gp);
	  xdr_destroy(&xdrs);
	  }
      pParticles[i]->x=t_gp.pos[0];
      pParticles[i]->y=t_gp.pos[1];
      pParticles[i]->z=t_gp.pos[2];
      pParticles[i]->vx=t_gp.vel[0];
      pParticles[i]->vy=t_gp.vel[1];
      pParticles[i]->vz=t_gp.vel[2];
      pParticles[i]->mass=t_gp.mass;
      pParticles[i]->rho=t_gp.rho;
      pParticles[i]->temp=t_gp.temp;
      pParticles[i]->metal=t_gp.metals;
      //cerr << pParticles[i]->temp << endl;

      pParticles[i]->ne=-1;
     
      pParticles[i]->type=CParticle::gas;
      softening[i] = t_gp.hsmooth;
      potential[i] = t_gp.phi;
      i++;
    }

    for(unsigned int n=header.ndark;n!=0;--n) {
      file.read((char*) &t_dp,sizeof(tipsy_dark_particle));
      if(xdrdata) {
	  XDR xdrs;
	  xdrmem_create(&xdrs, reinterpret_cast<char *>(&t_dp),
			sizeof(tipsy_dark_particle), XDR_DECODE);
	  xdr_template(&xdrs, &t_dp);
	  xdr_destroy(&xdrs);
	  }
      pParticles[i]->x=t_dp.pos[0];
      pParticles[i]->y=t_dp.pos[1];
      pParticles[i]->z=t_dp.pos[2];
      pParticles[i]->vx=t_dp.vel[0];
      pParticles[i]->vy=t_dp.vel[1];
      pParticles[i]->vz=t_dp.vel[2];
      pParticles[i]->mass=t_dp.mass;
      pParticles[i]->type=CParticle::dm;
      softening[i] = t_dp.eps;
      potential[i] = t_dp.phi;
      i++;
    }

    for(unsigned int n=header.nstar;n!=0;--n) {
      file.read((char*) &t_sp,sizeof(tipsy_star_particle));
      if(xdrdata) {
	  XDR xdrs;
	  xdrmem_create(&xdrs, reinterpret_cast<char *>(&t_sp),
			sizeof(tipsy_star_particle), XDR_DECODE);
	  xdr_template(&xdrs, &t_sp);
	  xdr_destroy(&xdrs);
	  }
      pParticles[i]->x=t_sp.pos[0];
      pParticles[i]->y=t_sp.pos[1];
      pParticles[i]->z=t_sp.pos[2];
      pParticles[i]->vx=t_sp.vel[0];
      pParticles[i]->vy=t_sp.vel[1];
      pParticles[i]->vz=t_sp.vel[2];
      pParticles[i]->mass=t_sp.mass;
      pParticles[i]->metal=t_gp.metals;
      pParticles[i]->type=CParticle::star;
      softening[i] = t_sp.eps;
      potential[i] = t_sp.phi;
      tform[i]=t_sp.tform;
      i++;
    }
    
    cerr << "done!" << endl;

    
    
    string fname_units = (string)fname + ".units";
    
	if(!siman::fileExists(fname_units))
      fname_units = "tipsy.units";

    ifstream file_units(fname_units.c_str());


    boxSize = 1.;
    
    double distanceIntToPhys;
    double massIntToPhys;
    double velIntToPhys;

    if(file_units.is_open()) {
      cerr << "CTipsyFile: Loading conversion units from " << fname_units;

      file_units >> lenUnits;
      file_units >> massUnits;
      file_units >> velUnits;
      file_units >> denUnits;
      file_units >> enUnits;

      

      cerr << ": " << lenUnits << " | " <<  massUnits << " | " << velUnits << " | " << denUnits << " | " << enUnits << endl;
      
      
      // read in from FILE 
    } else {
      cerr << "CTipsyFile: >>> WARNING - physical unit conversion file " << fname_units << " was not found." << endl;
      
      cerr << "                Please store in requested file - MSolUnit KpcUnit Km/sUnit" << endl;
      cerr << "                (This doesn't matter if you're not doing any physics with this file at the moment)" << endl;
     
    }

    
    redshift =  1./header.time-1.;
    cerr << "CTipsyFile: redshift = " << redshift << endl;
    om_m0 = 0.3;
    om_lam0 = 0.7;


    string fname_HI = (string) fname + ".HI";
    ifstream file_HI(fname_HI.c_str());
    if(file_HI.is_open()) {
      cerr << "CTipsyFile: loading HI data from " << fname_HI << "...";
      int n_HI_dat;
      file_HI >> n_HI_dat;
      if(n_HI_dat!=numParticles) {
	cerr << "length mismatch ("<<numParticles << "/" << n_HI_dat<<"), aborted HI load" << endl;
      } else {
	float Xfrac = 1-getHeliumMassFrac();
	for(int n=0; n<numParticles;n++) {
	  float nH0;
	  file_HI >> nH0;
	  // nH0 is stored as the number of particles PER 1 proton mass of total density

	  // we want n_HII as a number proportion of N_Htot:
	  pParticles[n]->nHp=(Xfrac-nH0)/Xfrac;

	  // one electron for each ionised H:
	  pParticles[n]->ne=pParticles[n]->nHp;
	} 
	cerr << "done!" << endl;
	
	string fname_HeI = (string) fname + ".HeI";
	ifstream file_HeI(fname_HeI.c_str());
	string fname_HeII = (string) fname + ".HeII";
	ifstream file_HeII(fname_HeII.c_str());

	if(file_HeI.is_open() && file_HeII.is_open()) {
	  int n_HeI_dat, n_HeII_dat;
	  file_HeI >> n_HeI_dat;
	  file_HeII >> n_HeII_dat;
	  
	  if(n_HeI_dat == numParticles && n_HeII_dat == numParticles) {
	    cerr << "CTipsyFile: loading HeI & HeII data...";
	    float nY = getHeliumMassFrac()/(4*(1-getHeliumMassFrac()));  
	    for(int n=0; n<numParticles; n++) {
	      float HeI, HeII, HeIII;
	      file_HeI >> HeI;
	      file_HeII >> HeII;
	      HeIII = getHeliumMassFrac()/4-HeI-HeII;
	      pParticles[n]->ne+=(HeII+2*HeIII)*(1+4*nY); // hopefully - damn this is confusing!
	    }
	    cerr << "done!" << endl;
	  } else {
	    cerr << "length mismatch, failed to load" << endl;
	  }
	} else {
	  cerr << "CTipsyFile: warning - no HeI/HeII data, ne may be inaccurate" << endl;
	}

	UFromTemp();
	      
      } // if HI file is sane
      
    } else {
      cerr << "CTipsyFile: Did not find " << fname_HI << "; ionisation information unavailable." << endl;
    }

    return true;
  } else {
    
    cerr << "CTipsyFile: Failed to open file for reading" << endl;
    throw(CFileError(fname));
    return false;
  }

}

void CTipsyFile::nativeWrite(CSimSnap *s, string filename) {
  cerr << "CTipsyFile: writing " << filename << "..." << endl;

  // open file

  FILE *tipsyFp = fopen(filename.c_str(), "a");  // Create file
  if(tipsyFp == NULL) {
      assert(0);
      }
  fclose(tipsyFp);
  tipsyFp = fopen(filename.c_str(), "rb+");
  if(tipsyFp == NULL) {
      assert(0);
      }
  XDR xdrs;

  // Always use XDR format
  xdrstdio_create(&xdrs, tipsyFp, XDR_ENCODE);

  ofstream file(filename.c_str(),ios::binary);

  // count particles

  unsigned int numStar=0, numDM=0, numGas=0, numTot = s->getNumParticles();

  CParticle *p;
  unsigned int n;

  for(n=0; n<numTot; n++) {
    p=s->getParticle(n);
    switch(p->type) { 

    case CParticle::gas:
      numGas+=1;
      break;

    case CParticle::dm:
      numDM+=1;
      break;

    case CParticle::star:
      numStar+=1;
      break;

    }
    s->releaseParticle(p);
  }

  // create virtual simulations for ease of ordering
  // (tipsy files insist on specific ordering of particles)

  CParticleTypeFilter fDM(CParticle::dm);
  CParticleTypeFilter fStar(CParticle::star);
  CParticleTypeFilter fGas(CParticle::gas);

  CSubset sGas(s,fGas);
  CSubset sDM(s,fDM);
  CSubset sStar(s,fStar);

  // create union in correct order for writing!

  CUnion sOrdered(s);
  sOrdered.add(&sStar);
  sOrdered.add(&sDM);
  sOrdered.add(&sGas);
  
  assert(sOrdered.getNumParticles() == numTot);
  assert(numGas + numDM + numStar == numTot);

  float *tform = sOrdered.getArray("tform");
  if(tform == NULL && numStar > 0) {
      cerr << "WARNING: no star formation time" << endl;
      tform = sOrdered.createArray("tform", "TIPSY tform");
      }

  float *eps = sOrdered.getArray("softening");
  if(eps == NULL) {
      cerr << "WARNING: no softening" << endl;
      eps = sOrdered.createArray("softening", "Gravitational softening");
      }

  float *potential = sOrdered.getArray("potential");
  if(eps == NULL) {
      cerr << "WARNING: no softening" << endl;
      potential = sOrdered.createArray("potential", "Gravitational potential");
      }
  
  tipsy_header header;
  
  header.time = 1/(1+s->getRedshift());
  header.nbodies = numTot;
  header.ndim = 3;
  header.nsph = numGas;
  header.ndark = numDM;
  header.nstar = numStar;
  header.zero = 0;
  xdr_template(&xdrs, &header);

  for(n=0; n<numGas; n++) {
    p =  sOrdered.getParticle(n);
    
    tipsy_gas_particle gp;
    
    gp.pos[0] = p->x;
    gp.pos[1] = p->y;
    gp.pos[2] = p->z;
    gp.vel[0] = p->vx;
    gp.vel[1] = p->vy;
    gp.vel[2] = p->vz;
    gp.mass = p->mass;
    gp.rho = p->rho;
    gp.temp = p->temp;
    gp.hsmooth = eps[n]; // XXX The arrays may not be sorted as the particles.
    gp.metals = p->metal;
    gp.phi = potential[n];
    
    if(!xdr_template(&xdrs, &gp)) {
	assert(0);
	}
    sOrdered.releaseParticle(p);
  }

  for(n=numGas; n<numGas+numDM; n++) {
    p =  sOrdered.getParticle(n);
    
    tipsy_dark_particle dp;
    
    dp.pos[0] = p->x;
    dp.pos[1] = p->y;
    dp.pos[2] = p->z;
    dp.vel[0] = p->vx;
    dp.vel[1] = p->vy;
    dp.vel[2] = p->vz;
    dp.mass = p->mass;
    dp.eps = eps[n];		// XXX are these sorted?
    dp.phi = potential[n];
    
    if(!xdr_template(&xdrs, &dp)) {
	assert(0);
	}
    sOrdered.releaseParticle(p);
  }

  for(n=numGas+numDM; n<numTot; n++) {
    p =  sOrdered.getParticle(n);
    
    tipsy_star_particle sp;
    
    sp.pos[0] = p->x;
    sp.pos[1] = p->y;
    sp.pos[2] = p->z;
    sp.vel[0] = p->vx;
    sp.vel[1] = p->vy;
    sp.vel[2] = p->vz;
    sp.mass = p->mass;
    sp.metals = p->metal;
    sp.tform = tform[n-(numGas+numDM)];
    sp.eps = eps[n];
    sp.phi = potential[n];
    
    if(!xdr_template(&xdrs, &sp)) {
	assert(0);
	}
    sOrdered.releaseParticle(p);
  }
}


