// tipsysnap.cpp - part of SimAn Simulation Analysis Library
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

namespace siman {


  // CONSTRUCTOR / DESTRUCTORS

  TipsySnap::~TipsySnap() {
    if(getVerbose()>2)
      cerr << "~TipsySnap" << endl;
  }

  TipsySnap::TipsySnap(const char *fname_in) {

    strcpy(fname,fname_in);
    load();
  }



  // LOAD TIPSY FILE

  bool TipsySnap::load() {
  
    bool flipendian = false;

    // returns true for success

    ifstream file(fname,ios::binary);

  
  
    if(file.is_open()) {

      if(getVerbose()==1) 
	cerr << "TipsySnap: reading " << fname << "...";

      // READ HEADER

      file.read((char*) &header,sizeof(header));
    
      if(header.nsph + header.ndark + header.nstar != header.nbodies || header.ndim>100) {
    
	endian::flip4(&header.nsph);
	endian::flip4(&header.ndark);
	endian::flip4(&header.nstar);
	endian::flip4(&header.nbodies);
	endian::flip4(&header.ndim);
	endian::flip8(&header.time);
    
	if(header.nsph + header.ndark + header.nstar == header.nbodies && header.ndim<100) {
	  if(getVerbose()>1)
	    cerr << "TipsySnap: Wrong endian, auto-flipping (performance warning)" << endl;
	  flipendian = true;
       	
	} else {
	  if(getVerbose()>0)
	    cerr << "TipsySnap: Don't recognise this file format" << endl;
	  throw(FileError(fname));
	}

      }

      // some file have a mysterious "zero" word, others don't
      // if this one doesn't, seek back to counteract reading it

      if(header.zero!=0)
	file.seekg(-4,ios_base::cur);

      if(header.ndim!=3) {
	if(getVerbose()>0)
	  cerr << "TipsySnap: " << header.ndim << " dimensions not supported" << endl;
	throw(FileError(fname));
      }
    
      if(getVerbose()>1)
	cerr << "TipsySnap: time " << header.time << "; G: " << header.nsph << " D: " << header.ndark << " S: " << header.nstar << endl;


      // ALLOCATE MEMORY

      numParticles = header.nsph+header.ndark+header.nstar;

      allocateMemory();


      // READ PARTICLE DATA

      

      boxSize = 1.;
    
      string fname_units = (string)fname + ".units";
    
      if(!siman::fileExists(fname_units))
	fname_units = siman::withPathOf(fname,"tipsy.units");

      ifstream file_units(fname_units.c_str());

      if(file_units.is_open()) {
	if(getVerbose()>1)
	  cerr << "TipsySnap: Loading conversion units from " << fname_units;

	file_units >> lenUnits;
	file_units >> massUnits;
	file_units >> velUnits;
	file_units >> denUnits;
	file_units >> enUnits;
	file_units >> hubble;
      
	if(getVerbose()>1)
	  cerr << ": " << lenUnits << " | " <<  massUnits << " | " << velUnits << " | " << denUnits << " | " << enUnits << endl;
      
      
	// read in from FILE 
      } else {
	if(getVerbose()>0) {
	  cerr << "TipsySnap: >>> WARNING - physical unit conversion file " << fname_units << " was not found." << endl;
	
	  cerr << "                Please store in requested file, one per line, - MSolUnit KpcUnit Km/sUnit." << endl;
	  cerr << "                Append value of h = H_0/(100 km/s/Mpc) in final line if desired" << endl;
	}
	lenUnits=Unit("kpc");
	massUnits=Unit("Msol");
	velUnits=Unit("km s^-1");
	denUnits=Unit("Msol kpc^-3");
	enUnits=velUnits*velUnits;
	
      }

    
      redshift =  1./header.time-1.;
      if(getVerbose()>1)
	cerr << "TipsySnap: redshift = " << redshift << endl;
      om_m0 = 0.3;
      om_lam0 = 0.7;

      if(getVerbose()>1)
	cerr << "TipsySnap: reading...";

      
      unsigned int i=0;

      tipsy::gas_particle t_gp;
      tipsy::dark_particle t_dp;
      tipsy::star_particle t_sp;

      SimanArray &tform(createArray("tform","Time of formation",Unit("1.e8 yr")));
      // SimanArray &smooth(createArray("smooth","Smoothing length",lenUnits));
      double conv = (lenUnits/velUnits).convertTo(Unit("1.e8 yr"),this);

      initialiseSPH(32);

      for(unsigned int n=header.nsph;n!=0;--n) {
	file.read((char*) &t_gp,sizeof(tipsy::gas_particle));
	if(flipendian) endian::swapDataEndian(&t_gp,4,sizeof(tipsy::gas_particle));      
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
	pParticles[i]->metal=t_gp.metals;   
	pParticles[i]->type=Particle::gas;
	pSmooth[i]=t_gp.hsmooth;
	// smooth[i]=t_gp.hsmooth;
	i++;
      }

      for(unsigned int n=header.ndark;n!=0;--n) {
	file.read((char*) &t_dp,sizeof(tipsy::dark_particle));
	if(flipendian) endian::swapDataEndian(&t_dp,4,sizeof(tipsy::dark_particle));
	pParticles[i]->x=t_dp.pos[0];
	pParticles[i]->y=t_dp.pos[1];
	pParticles[i]->z=t_dp.pos[2];
	pParticles[i]->vx=t_dp.vel[0];
	pParticles[i]->vy=t_dp.vel[1];
	pParticles[i]->vz=t_dp.vel[2];
	pParticles[i]->mass=t_dp.mass;
	pParticles[i]->type=Particle::dm;
	i++;
      }

      for(unsigned int n=header.nstar;n!=0;--n) {
	file.read((char*) &t_sp,sizeof(tipsy::star_particle));
	if(flipendian) endian::swapDataEndian(&t_sp,4,sizeof(tipsy::star_particle));
	pParticles[i]->x=t_sp.pos[0];
	pParticles[i]->y=t_sp.pos[1];
	pParticles[i]->z=t_sp.pos[2];
	pParticles[i]->vx=t_sp.vel[0];
	pParticles[i]->vy=t_sp.vel[1];
	pParticles[i]->vz=t_sp.vel[2];
	pParticles[i]->mass=t_sp.mass;
	pParticles[i]->metal=t_sp.metals;
	pParticles[i]->type=Particle::star;
	tform[i]=t_sp.tform*conv;
	i++;
      }
    
      if(getVerbose()>1)
	cerr << "done!" << endl;

    
    
      SimanArray & nHpArray = createArray("nHII","HII per (HI+HII)");
      SimanArray & nHepArray = createArray("nHeII","HeII per (HI+HII)");
      SimanArray & nHeppArray = createArray("nHeIII","HeIII per (HI+HII)");

      string fname_HI = (string) fname + ".HI";
      ifstream file_HI(fname_HI.c_str());
      if(file_HI.is_open()) {
	if(getVerbose()>1)
	  cerr << "TipsySnap: loading HI data from " << fname_HI << "...";
	int n_HI_dat;
	file_HI >> n_HI_dat;
	if((unsigned int)n_HI_dat!=numParticles) {
	  if(getVerbose()>1)
	    cerr << "length mismatch ("<<numParticles << "/" << n_HI_dat<<"), aborted HI load" << endl;
	} else {
	  float Xfrac = 1-getHeliumMassFrac();
	  for(unsigned int n=0; n<numParticles;n++) {
	    float nH0;
	    file_HI >> nH0;
	    // nH0 is stored as the number of particles PER 1 proton mass of total density

	    // we want n_HII as a number proportion of N_Htot:
	    nHpArray[n]=(Xfrac-nH0)/Xfrac;
	   
	    // one electron for each ionised H:
	    pParticles[n]->ne=nHpArray[n];
	  } 
	  if(getVerbose()>1)
	    cerr << "done!" << endl;
	
	  string fname_HeI = (string) fname + ".HeI";
	  ifstream file_HeI(fname_HeI.c_str());
	  string fname_HeII = (string) fname + ".HeII";
	  ifstream file_HeII(fname_HeII.c_str());

	  if(file_HeI.is_open() && file_HeII.is_open()) {
	    unsigned int n_HeI_dat, n_HeII_dat;
	    file_HeI >> n_HeI_dat;
	    file_HeII >> n_HeII_dat;
	  
	    if(n_HeI_dat == numParticles && n_HeII_dat == numParticles) {
	      if(getVerbose()>1)
		cerr << "TipsySnap: loading HeI & HeII data...";
	      float nY = getHeliumMassFrac()/(4*(1-getHeliumMassFrac()));  
	      for(unsigned int n=0; n<numParticles; n++) {
		float HeI, HeII, HeIII;
		file_HeI >> HeI;
		file_HeII >> HeII;

		HeI = HeI/Xfrac;
		HeII = HeII/Xfrac;
		HeIII = (1-Xfrac)/(4*Xfrac) - HeI - HeII;
		pParticles[n]->ne+=(HeII+2*HeIII); // hopefully 
		nHepArray[n] = HeII;
		nHeppArray[n] = HeIII;

	      }
	      if(getVerbose()>1)
		cerr << "done!" << endl;
	    } else {
	      if(getVerbose()>1)
		cerr << "length mismatch, failed to load" << endl;
	    }
	  } else {
	    if(getVerbose()>1)
	      cerr << "TipsySnap: warning - no HeI/HeII data, ne may be inaccurate" << endl;
	    if(getVerbose()==1)
	      cerr << " [inaccurate-ion] ";
	  }

	  UFromTemp();
	      
	} // if HI file is sane
      
      } else {
	if(getVerbose()>1)
	  cerr << "TipsySnap: Did not find " << fname_HI << "; ionisation information unavailable." << endl;
	if(getVerbose()==1)
	  cerr << " [no-ion] ";
      }

      if(getVerbose()==1) 
	cerr << "done!" << endl;
	
      return true;
    } else {
      if(getVerbose()>=1)
	cerr << "TipsySnap: Failed to open file for reading" << endl;
      throw(FileError(fname));
      return false;
    }

  }
  
  void TipsySnap::nativeWrite(const SimSnap *s, string filename) {


    if(getVerbose()>1)
      cerr << "TipsySnap: writing " << filename << "..." << endl;
  
    // open file

    ofstream file(filename.c_str(),ios::binary);
 
    // create virtual simulations with only one particle type

    ParticleTypeFilter fDM(Particle::dm);
    ParticleTypeFilter fStar(Particle::star);
    ParticleTypeFilter fGas(Particle::gas);

    Subset sGas(const_cast<SimSnap*>(s),fGas);
    Subset sDM(const_cast<SimSnap*>(s),fDM);
    Subset sStar(const_cast<SimSnap*>(s),fStar);

    tipsy::header header;
  
    // fill in header

    memset(&header,0,sizeof(header));


    header.time = 1/(1+s->getRedshift());
    header.nbodies = s->getNumParticles();
    header.ndim = 3;
    header.nsph = sGas.getNumParticles();
    header.ndark = sDM.getNumParticles();
    header.nstar = sStar.getNumParticles();

    file.write((char*)&header,sizeof(header));

  
    tipsy::gas_particle t_gp;
    tipsy::dark_particle t_dp;
    tipsy::star_particle t_sp;

    if(getVerbose()>1)
      cerr << "TipsySnap: writing particle data... gas ";

      
    for(unsigned int n=0; n<sGas.getNumParticles(); n++) {

      const Particle *p = sGas.getConstParticle(n);
      t_gp.pos[0]=p->x;
      t_gp.pos[1]=p->y;
      t_gp.pos[2]=p->z;
      t_gp.vel[0]=p->vx;
      t_gp.vel[1]=p->vy;
      t_gp.vel[2]=p->vz;
      t_gp.mass=p->mass;
      t_gp.rho=p->rho;
      t_gp.temp=p->temp;
      t_gp.hsmooth=0;
      t_gp.metals=p->metal;
      t_gp.phi=0;
      file.write((char*) &t_gp,sizeof(tipsy::gas_particle));
    }

    if(getVerbose()>1)
      cerr << "dm ";
    for(unsigned int n=0; n<sDM.getNumParticles(); n++) {
	
      const Particle *p = sDM.getConstParticle(n);
      t_dp.pos[0]=p->x;
      t_dp.pos[1]=p->y;
      t_dp.pos[2]=p->z;
      t_dp.vel[0]=p->vx;
      t_dp.vel[1]=p->vy;
      t_dp.vel[2]=p->vz;
      t_dp.mass=p->mass;
      t_dp.eps=0;
      t_dp.phi=0;
      file.write((char*) &t_dp, sizeof(tipsy::dark_particle));
    }

    const SimanArray *pTformArr = NULL;
      
    double conv = 0.;

    try {
      pTformArr = &(sStar.getConstArray("tform"));
      conv = pTformArr->getUnits().convertTo(s->getDistanceUnits()/s->getVelocityUnits());

    } catch(UnknownArray &e) { }
    catch(UnitsError &e) { 
      if(getVerbose()>0)
	cerr << "TipsySnap: warning - unknown units for tform, writing without conversion" << endl;
      conv = 1.;
    }
 
    if(getVerbose()>1)
      cerr << "stars ";
    for(unsigned int n=0; n<sStar.getNumParticles(); n++) {
	
      const Particle *p = sStar.getConstParticle(n);
      t_sp.pos[0]=p->x;
      t_sp.pos[1]=p->y;
      t_sp.pos[2]=p->z;
      t_sp.vel[0]=p->vx;
      t_sp.vel[1]=p->vy;
      t_sp.vel[2]=p->vz;
      t_sp.mass=p->mass;
      t_sp.metals = p->metal;
      if(pTformArr!=NULL)
	t_sp.tform =  (*pTformArr)[n];
      else
	t_sp.tform = 0;
      t_sp.eps=0;
      t_sp.phi=0;
	
      file.write((char*) &t_sp, sizeof(tipsy::star_particle));
    }

    if(getVerbose()>1)
      cerr << "done!" << endl;

      
    string filename_units = filename+".units";
      
    ofstream file_units(filename_units.c_str());
    file_units << s->getDistanceUnits() << endl << s->getMassUnits() << endl  << s->getVelocityUnits() << endl << s->getDensityUnits() << endl << s->getEnergyUnits() << endl << s->getHubble() << endl;
      

  }

} // namespace siman
