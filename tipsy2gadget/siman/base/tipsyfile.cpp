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
#include "endian.hpp"

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
  
  bool flipendian = false;

  // returns true for success

  ifstream file(fname,ios::binary);

  ENDINIT;
  
  if(file.is_open()) {

    // READ HEADER

    file.read((char*) &header,sizeof(header));
    
    if(header.nsph + header.ndark + header.nstar != header.nbodies || header.ndim>100) {
    
      END4(&header.nsph);
      END4(&header.ndark);
      END4(&header.nstar);
      END4(&header.nbodies);
      END4(&header.ndim);
      END8(&header.time);
    
      if(header.nsph + header.ndark + header.nstar == header.nbodies && header.ndim<100) {
	cerr << "CTipsyFile: Wrong endian, auto-flipping (performance warning)" << endl;
	flipendian = true;
       	
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

    for(unsigned int n=header.nsph;n!=0;--n) {
      file.read((char*) &t_gp,sizeof(tipsy_gas_particle));
      if(flipendian) swapDataEndian(&t_gp,4,sizeof(tipsy_gas_particle));      
      pParticles[i]->x=t_gp.pos[0];
      pParticles[i]->y=t_gp.pos[1];
      pParticles[i]->z=t_gp.pos[2];
      pParticles[i]->vx=t_gp.vel[0];
      pParticles[i]->vy=t_gp.vel[1];
      pParticles[i]->vz=t_gp.vel[2];
      pParticles[i]->mass=t_gp.mass;
      pParticles[i]->rho=t_gp.rho;
      pParticles[i]->temp=t_gp.temp;
      //cerr << pParticles[i]->temp << endl;

      pParticles[i]->ne=-1;
     
      pParticles[i]->type=CParticle::gas;
      i++;
    }

    for(unsigned int n=header.ndark;n!=0;--n) {
      file.read((char*) &t_dp,sizeof(tipsy_dark_particle));
      if(flipendian) swapDataEndian(&t_dp,4,sizeof(tipsy_dark_particle));
      pParticles[i]->x=t_dp.pos[0];
      pParticles[i]->y=t_dp.pos[1];
      pParticles[i]->z=t_dp.pos[2];
      pParticles[i]->vx=t_dp.vel[0];
      pParticles[i]->vy=t_dp.vel[1];
      pParticles[i]->vz=t_dp.vel[2];
      pParticles[i]->mass=t_dp.mass;
      pParticles[i]->type=CParticle::dm;
      i++;
    }

    for(unsigned int n=header.nstar;n!=0;--n) {
      file.read((char*) &t_sp,sizeof(tipsy_star_particle));
      if(flipendian) swapDataEndian(&t_sp,4,sizeof(tipsy_star_particle));
      pParticles[i]->x=t_sp.pos[0];
      pParticles[i]->y=t_sp.pos[1];
      pParticles[i]->z=t_sp.pos[2];
      pParticles[i]->vx=t_sp.vel[0];
      pParticles[i]->vy=t_sp.vel[1];
      pParticles[i]->vz=t_sp.vel[2];
      pParticles[i]->mass=t_sp.mass;
      pParticles[i]->type=CParticle::star;
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

