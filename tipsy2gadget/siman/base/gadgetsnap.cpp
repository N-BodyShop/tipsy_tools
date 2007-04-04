// gadgetsnap.cpp - part of SimAn Simulation Analysis Library
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
#include <boost/lexical_cast.hpp>

namespace siman {

#define SKIP fread(&dummy, sizeof(dummy), 1, fd); 


// CONSTRUCTOR / DESTRUCTORS

GadgetSnap::~GadgetSnap() {
  
  if(Id!=NULL) free((void*)Id);
  
}

GadgetSnap::GadgetSnap(const char *fname_in, const bool swapEndian_in) {
  strcpy(fname,fname_in);
  
  Id=NULL;
  swapEndian = swapEndian_in;  // N.B. deprecated - this is autodetected in Load()
  load();
}

GadgetSnap::GadgetSnap(const char *path, const char *snapshot, const int snapid, const bool swapEndian_in) {

    sprintf(fname,  "%s/%s_%03d", path, snapshot,snapid);
    
    swapEndian = swapEndian_in; // N.B. deprecated - this is autodetected in Load()
    load();
}


// SIMPLE VALUE RETURNS


#define CHECKLENGTHFIELD(expect,section) fread(&checklen,sizeof(int),1,fd); if(swapEndian) { endian::flip4(&checklen); } if(checklen!=expect) { cerr << endl << "GadgetSnap: >>>Error in length field (" << section << ": says " << checklen << ", expected " << expect << ")" << endl; fileDebug(fd); }


// LOAD GADGET FILE

bool GadgetSnap::load() {
  
  // adapted from Volker Springel's readgadget.c
  // a little "hybrid" in C/C++ thinking...

  // returns true for success
  FILE *fd;
  char   buf[200];
  int    k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;
  
  unsigned int checklen;

  sprintf(buf,"%s",fname);
  
  if(!(fd=fopen(buf,"r")))
    {
      if(getVerbose()>0)
	cerr << "GadgetSnap: Error opening file " << buf << endl;
      FileError fe((string) fname);
      throw(fe);
      
      return false;
    }

  if(getVerbose()>1)
    cerr << "GadgetSnap: reading " << buf << endl;
 


  fread(&dummy, sizeof(dummy), 1, fd);
  if(dummy!=sizeof(header)) {
    endian::flip4(&dummy);
    if(dummy==sizeof(header)) {
      if(getVerbose()>1)
	cerr << "GadgetSnap: wrong endian, auto-flipping (performance warning)" << endl;
      swapEndian=true;
    } else {
      if(getVerbose()>0)
	cerr << "GadgetSnap: Don't recognise this file format" << endl;
      FileError fe((string) fname);
      throw(fe); 
      return false;
    }
  } else {
    swapEndian = false;
  }

  fread(&header, sizeof(header), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);

      
  if(swapEndian==1) {
    endian::swapDataEndian(&header, 4, 6*4);           // 6 integers
    endian::swapDataEndian(&header.mass, 8, 8*8);     // 8 doubles
    endian::swapDataEndian(&header.flag_sfr, 4, 40);  // 10 more integers
    endian::swapDataEndian(&header.BoxSize,8, 4*8);   // 4 more doubles
  }
     
  if(getVerbose()>1) {
    cerr << "GadgetSnap: z=" << header.redshift << "; t=" << header.time << "; h=" << header.HubbleParam << "; L=" << header.BoxSize << endl;
    cerr << "GadgetSnap: npart: G=" << header.npart[0] << "; DM=" << header.npart[1] << "; S=" << header.npart[4];
    cerr << "GadgetSnap: header masses   G: " << header.mass[0] << "   DM: "<< header.mass[1] << "  S: " << header.mass[4] << endl;
    cerr << "GadgetSnap: non-public flags - cooling " << string(header.flag_cooling?"YES":"NO") << " stellar age " << string(header.flag_stellarage?"YES":"NO") << " metals " << string(header.flag_metals?"YES":"NO");
  }

  for(k=0, numParticles=0, ntot_withmasses=0; k<5; k++)
    numParticles+= header.npart[k];
      

  for(k=0, ntot_withmasses=0; k<5; k++)
    {
      if(header.mass[k]==0)
	ntot_withmasses+= header.npart[k];
    }

  boxSize = header.BoxSize;
  hubble = header.HubbleParam;
  redshift = header.redshift;
      
  // UNITS:
  
      
  string fname_units = (string)fname + ".units";
  if(!siman::fileExists(fname_units))
    fname_units = withPathOf(fname,"gadget.units");
  
  ifstream file_units(fname_units.c_str());
  
  
  
  if(file_units.is_open()) {
    if(getVerbose()>1)
      cerr << "GadgetSnap: Loading conversion units from " << fname_units;
    
    file_units >> lenUnits;
    file_units >> massUnits;
    file_units >> velUnits;
    file_units >> denUnits;
    file_units >> enUnits;
  } else {

    if(getVerbose()>0) {
      cerr << "GadgetSnap: Warning - using default gadget units system which assumes GADGET was run in COSMOLOGICAL mode" << endl;
      cerr << "             Avoid this warning by storing units system in gadget.units or " + (string)fname + ".units" << endl;
    }
      // the following may be incorrect depending on configuration of GADGET
    
    
    lenUnits.set("kpc a h^-1");
    //lenUnits/=(header.HubbleParam*(header.redshift+1.));
    
    massUnits.set("1.e10 Msol h^-1");
    // massUnits*=1.e10/(header.HubbleParam);
    
    velUnits.set("km s^-1 a^1/2");
    
    /*
      velUnits.set(vel_kmPerS);
      velUnits*=sqrt(1./(1+header.redshift));
    */
    
    denUnits = massUnits/(lenUnits*lenUnits*lenUnits);
    
    enUnits.set("km^2 s^-2");
  }

  om_m0 = header.Omega0;
  om_lam0 = header.OmegaLambda;
  
  allocateMemory();

  pc = 0;

  if(getVerbose()>1)
    cerr << "GadgetSnap: Reading positions...";
  // fileDebug(fd);
  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"PosStart");
  for(k=0,pc_new=pc;k<6;k++)
    {
      if(getVerbose()>1)
	cerr << ".";
     
      for(n=0;n<header.npart[k];n++)
	{
	  fread(&(pParticles[pc_new]->x), sizeof(float), 1, fd);
	  fread(&(pParticles[pc_new]->y), sizeof(float), 1, fd);
	  fread(&(pParticles[pc_new]->z), sizeof(float), 1, fd);
	  if(swapEndian) {
	    endian::flip4(&(pParticles[pc_new]->x));
	    endian::flip4(&(pParticles[pc_new]->y));
	    endian::flip4(&(pParticles[pc_new]->z));
	  }
	  pc_new++;
	}
    }
  // fileDebug(fd);
  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"PosEnd");
  if(getVerbose()>1)
    cerr << "done!" << endl;

  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"VelStart");
  // fileDebug(fd);
  
  if(getVerbose()>1)
    cerr << "GadgetSnap: Reading velocities...";
  for(k=0,pc_new=pc;k<6;k++)
    {
      if(getVerbose()>1)
	cerr << ".";

      for(n=0;n<header.npart[k];n++)
	{
	  fread(&(pParticles[pc_new]->vx), sizeof(float), 3, fd);
	  if(swapEndian) {
	    endian::flip4(&(pParticles[pc_new]->vx));
	    endian::flip4(&(pParticles[pc_new]->vy));
	    endian::flip4(&(pParticles[pc_new]->vz));
	  }
	  pc_new++;
	}
    }
  if(getVerbose()>1)
    cerr << "done!" << endl;

  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"VelEnd");
  CHECKLENGTHFIELD(sizeof(int)*numParticles,"PID-Start");

  if(getVerbose()>1)
    cerr << "GadgetSnap: Skipping PID data (not currently supported)...";

  int ignore;
  for(k=0,pc_new=pc;k<6;k++)
    {
      if(getVerbose()>1)
	cerr << ".";
      for(n=0;n<header.npart[k];n++)
	{
	  fread(&ignore, sizeof(int), 1, fd);
	      
	  pc_new++;
	}
    }

  if(getVerbose()>1)
    cerr << "done!" << endl;

  CHECKLENGTHFIELD(sizeof(int)*numParticles,"PID-End");
           
  if(ntot_withmasses>0) {
    CHECKLENGTHFIELD(sizeof(float)*ntot_withmasses,"Mass-Start");
    if(getVerbose()>1)
      cerr << "GadgetSnap: Reading mass data...";
  }
      
      
  
  for(k=0, pc_new=pc; k<6; k++) {
    if(getVerbose()>1)
      cerr << ".";
    int writetype =0;
    switch(k) {
    case 0:
      writetype=Particle::gas;
      break;
    case 1:
      writetype=Particle::dm;
      break;
    case 4:
      writetype=Particle::star;
      break;
    }
    
    for(n=0;n<header.npart[k];n++)
      {
	pParticles[pc_new]->type=writetype;
	if(header.mass[k]==0) {
	  fread(&(pParticles[pc_new]->mass), sizeof(float), 1, fd);
	  if(swapEndian) { endian::flip4(&(pParticles[pc_new]->mass)); }
	}
	else {
	  pParticles[pc_new]->mass = header.mass[k];		
	}
	pc_new++;
	
      }
  }

  
  if(ntot_withmasses>0) {
    if(getVerbose()>1)
      cerr << "done!" << endl;
    CHECKLENGTHFIELD(sizeof(float)*ntot_withmasses,"Mass-End");
  }
      

  if(getVerbose()>1)
    cerr << "GadgetSnap: Reading gas data...";

      
      

  if(header.npart[0]>0)
    {
      if(getVerbose()>1)
	cerr << "energy ";
      
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"U-Start");
      for(n=0, pc_sph=pc; n<header.npart[0];n++)
	{
	  fread(&pParticles[pc_sph]->u, sizeof(float), 1, fd);
	  if(swapEndian) {endian::flip4(&(pParticles[pc_sph]->u)); }
	  pc_sph++;
	}
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"U-End");
      
      if(getVerbose()>1)
	cerr << "density ";

      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Rho-Start");
      for(n=0, pc_sph=pc; n<header.npart[0];n++)
	{
	  fread(&(pParticles[pc_sph]->rho), sizeof(float), 1, fd);
	  if(swapEndian) { endian::flip4(&(pParticles[pc_sph]->rho)); }
	  pc_sph++;
	}
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Rho-End");
    }

  // Any extra blocks?

      
  string fname_blocks = (string)fname + ".blocks";
  if(!siman::fileExists(fname_blocks))
    fname_blocks = "gadget.blocks";
  
  ifstream file_blocks(fname_blocks.c_str());
  
  int line = 1;
  if(file_blocks.is_open() && file_blocks.good()) {
    if(getVerbose()>1)
      cerr << endl << "GadgetSnap: Using block definition from " << fname_blocks << endl;
    

    char block_directive[1024];

    while(file_blocks.getline(block_directive,1024)) {
      
      vector<string> toks;
      tokenize(string(block_directive),toks,"|",true);
      
      if(toks.size()>=2) {
	
	string unit_s = "";
	if(toks.size()>2)
	  unit_s = toks[2];
	
	SimanArray *par=NULL;
	
	try {
	  par = &createArray(toks[0],toks[0],Unit(unit_s));
	} catch(SimanException &e) {
	  par = &getArray(toks[0]);
	}
	
	SimanArray &ar(*par);
	
	bool hasGas = false, hasDM = false, hasStars = false;
	if(toks[1].find("g")!=string::npos || toks[1].find("G")!=string::npos) hasGas=true;
	if(toks[1].find("d")!=string::npos || toks[1].find("D")!=string::npos) hasDM = true;
	if(toks[1].find("st")!=string::npos || toks[1].find("St")!=string::npos || toks[1].find("ST")!=string::npos) hasStars = true;
	
	int ntot = hasGas?header.npart[0]:0 + hasDM?header.npart[1]:0 + hasStars?header.npart[4]:0;
	
	CHECKLENGTHFIELD(sizeof(float)*ntot,"RB-Start");
	unsigned int offset = 0;
	for(int type=0; type<6; type++) {
	  if((type==0 && hasGas) || (type==1 && hasDM) || (type==4 & hasStars)) {
	    
	    for(n=0; n<header.npart[type];n++)
	      {
		fread(&(ar[offset]), sizeof(float), 1, fd);
		if(swapEndian) 
		{ endian::flip4(&(ar[offset])); 
		}
		offset++;
	      }
	  } else {
	    offset+=header.npart[type];
	  }
	}
	
	CHECKLENGTHFIELD(sizeof(float)*ntot,"RB-End");
      }
    }
  } else {

    // no block definition file
    
    
      if(header.flag_cooling) {
	if(getVerbose()>1)    
	  cerr << "ne ";
	
	CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Ne-Start");
	for(n=0, pc_sph=pc; n<header.npart[0];n++) {
	  
	  fread(&(pParticles[pc_sph]->ne), sizeof(float), 1, fd);
	  if(swapEndian) { endian::flip4(&(pParticles[pc_sph]->ne)); }
	  
	  
	  pc_sph++;
	}
	
	CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Ne-End");
	
      }

      if(header.flag_metals) {
	if(getVerbose()>1)    
	  cerr << "metals ";
	
	CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Met-Start");
	for(n=0, pc_sph=pc; n<header.npart[0];n++) {
	  
	  fread(&(pParticles[pc_sph]->metal), sizeof(float), 1, fd);
	  if(swapEndian) { endian::flip4(&(pParticles[pc_sph]->metal)); }
	  
	  
	  pc_sph++;
	}
	pc_sph+=header.npart[1];
	for(n=0; n<header.npart[4];n++) {
	  
	  fread(&(pParticles[pc_sph]->metal), sizeof(float), 1, fd);
	  if(swapEndian) { endian::flip4(&(pParticles[pc_sph]->metal)); }
	  
	  
	  pc_sph++;
	}
	
	//	int pc_star = header.npart[0]+header.npart[2];
	//int comp;
	
	CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Met-End");
	
      }
  }

  
  if(getVerbose()>1)
    cerr << "done!" << endl;

  /*  
  int ukblock = 0;
 
    
  while(!feof(fd)) {
    int lengthfield  = 0;
    int offset=0;
    fread(&lengthfield,sizeof(int),1,fd);
    fseek(fd,lengthfield,SEEK_CUR);
    CHECKLENGTHFIELD(lengthfield,"UKB-end");
    ukblock++;
  }

  if(ukblock>0 && getVerbose()>1) {
    cerr << "GadgetSnap: " << ukblock << " unidentified blocks at end of file";
  }
  */
  TempFromU();
  
  fclose(fd);
  return true;
  
}
      
void GadgetSnap::fileDebug(FILE *fd) {
  fpos_t init_pos;
  fgetpos(fd,&init_pos);
  
  printf("[%X] ",(unsigned int) ftell(fd));

  for(int n=0;n<20;n++) {
    unsigned char contents;
    fread(&contents,1,1,fd);
    printf("%02x ",(unsigned int)contents);
  }
  printf("\n");

  
  fsetpos(fd,&init_pos);
}



int GadgetSnap::Reorder(void)
{
  /*
  int i,j;
  int idsource, idsave, dest;
  gadget_particle psave, psource;


  printf("reordering....\n");

  for(i=0; i<NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  fprintf(stderr,"GadgetSnap: done.\n");

  free(Id);

  fprintf(stderr,"GadgetSnap: space for particle ID freed\n");
  */
	return 0;
}


void GadgetSnap::writeField(std::ofstream *file, char *buf, int size) {

  // FORTRAN-style field

  file->write((char*)&size,sizeof(size));
  file->write(buf,size);
  file->write((char*)&size,sizeof(size));
}


void GadgetSnap::nativeWrite(const SimSnap *s, string filename) {


  if(getVerbose()>1)
    cerr << "GadgetSnap: writing " << filename << "..." << endl;
  

  // not sure how this ought to be determined, as it
  // is a gadget-only flag:

  bool cooling=true;
  bool metals=true;

  // open file

  ofstream file(filename.c_str(),ios::binary);
  int sizefield = sizeof(gadget_header);

  // count particles

  unsigned int numStar=0, numDM=0, numGas=0, numTot = s->getNumParticles();
  float massStar=0, massDM=0, massGas=0;
  bool oneMassStar = true, oneMassDM = true, oneMassGas = true;

  const Particle *p;
  unsigned int n;

  for(n=0; n<numTot; n++) {
    p=s->getConstParticle(n);
    switch(p->type) { 

    case Particle::gas:
      numGas+=1;
      if(p->mass!=massGas && massGas!=0) oneMassGas=false;
      massGas = p->mass;
      break;

    case Particle::dm:
      numDM+=1;
      if(p->mass!=massDM && massDM!=0) oneMassDM=false;
      massDM = p->mass;
      break;

    case Particle::star:
      numStar+=1;
      if(p->mass!=massStar && massStar!=0) oneMassStar=false;
      massStar = p->mass;
      break;

    }

  }


  // create virtual simulations for ease of ordering
  // (gadget files insist on specific ordering of particles)

  ParticleTypeFilter fDM(Particle::dm);
  ParticleTypeFilter fStar(Particle::star);
  ParticleTypeFilter fGas(Particle::gas);

  Subset sGas(const_cast<SimSnap*>(s),fGas);
  Subset sDM(const_cast<SimSnap*>(s),fDM);
  Subset sStar(const_cast<SimSnap*>(s),fStar);

  // create union in correct order for writing!

  Union sOrdered(const_cast<SimSnap*>(s));
  sOrdered.add(&sStar);
  sOrdered.add(&sDM);
  sOrdered.add(&sGas);
  
 
  
  gadget_header header;
  
  // fill in header

  memset(&header,0,sizeof(header));

  header.HubbleParam = s->getHubble();
  header.redshift = s->getRedshift();
  header.time = 1/(1+s->getRedshift());
  header.BoxSize = s->getBoxSize();
  header.Omega0 = s->getOmegaM0();
  header.OmegaLambda = s->getOmegaLambda0();
  header.flag_cooling = cooling;
  header.flag_metals = metals;

  if(oneMassDM)
    header.mass[1] = massDM;

  if(oneMassGas)
    header.mass[0] = massGas;

  if(oneMassStar)
    header.mass[4] = massStar;

  header.npart[1] = header.npartTotal[1] =  numDM;
  header.npart[0] = header.npartTotal[0] = numGas;
  header.npart[4] = header.npartTotal[4] = numStar;

  writeField(&file,(char*)&header,sizeof(header));
  
  if(getVerbose()>1)
    cerr << "GadgetSnap: writing position data...";

  sizefield = sizeof(float) * 3 * numTot;
  
  file.write((char*)&sizefield,sizeof(int));

  for(n=0; n<numTot; n++) {
    p =  sOrdered.getConstParticle(n);
    file.write((char*)&(p->x),sizeof(float));
    file.write((char*)&(p->y),sizeof(float));
    file.write((char*)&(p->z),sizeof(float));
    
  }
  
  file.write((char*)&sizefield,sizeof(int));

  if(getVerbose()>1)
    cerr << "done!" << endl << "GadgetSnap: writing velocity data...";


  file.write((char*)&sizefield,sizeof(int));
  
  for(n=0; n<numTot; n++) {
    p =  sOrdered.getConstParticle(n);
    file.write((char*)&(p->vx),sizeof(float));
    file.write((char*)&(p->vy),sizeof(float));
    file.write((char*)&(p->vz),sizeof(float));
   
  }
  
  file.write((char*)&sizefield,sizeof(int));
  
  if(getVerbose()>1)
    cerr << "done!" << endl << "GadgetSnap: writing PID data...";
  
  sizefield = numTot * sizeof(int);

  file.write((char*)&sizefield,sizeof(int));
  
  for(n=0; n<numTot; n++) {
    file.write((char*)&n,sizeof(int));
  }
  
  file.write((char*)&sizefield,sizeof(int));

  if(getVerbose()>1)
    cerr << "done!" << endl;


  // Write masses, where required

  sizefield = ((oneMassDM?0:numDM) + (oneMassGas?0:numGas) +
	       (oneMassStar?0:numStar))*sizeof(float);

  if(sizefield>0) {
    
    if(getVerbose()>1)
      cerr << "GadgetSnap: writing mass data... (" << (oneMassDM?"":"D") << (oneMassStar?"":"S") << (oneMassGas?"":"G") << ")";

    int writeMassMask = (oneMassDM?0:Particle::dm) + (oneMassStar?0:Particle::star) + (oneMassGas?0:Particle::gas);
    
    file.write((char*)&sizefield,sizeof(int));
    
    
    for(n=0; n<numTot; n++) {
      p = s->getConstParticle(n);
      if((p->type & writeMassMask) > 0) {
	file.write((char*)&(p->mass),sizeof(float));
      }
     
    }

    file.write((char*)&sizefield,sizeof(int));

    if(getVerbose()>1)
      cerr << "done!" << endl;
  }

  if(numGas>0) {
    if(getVerbose()>1)
      cerr << "GadgetSnap: writing gas data...";
    sizefield = numGas * sizeof(float);

    file.write((char*)&sizefield,sizeof(int));

    for(n=0; n<numGas; n++) {
      p = sGas.getConstParticle(n);
      file.write((char*)&(p->u),sizeof(float));
      
    }

    file.write((char*)&sizefield,sizeof(int));

    file.write((char*)&sizefield,sizeof(int));

    for(n=0; n<numGas; n++) {
      p = sGas.getConstParticle(n);
      file.write((char*)&(p->rho),sizeof(float));
      
    }

    file.write((char*)&sizefield,sizeof(int));


    
    if(cooling) {
      
      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getConstParticle(n);
	file.write((char*)&(p->ne),sizeof(float));

      }
      file.write((char*)&sizefield,sizeof(int));
    }

     
    if(metals) {
      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getConstParticle(n);
	file.write((char*)&(p->metal),sizeof(float));
	
      }
      file.write((char*)&sizefield,sizeof(int));

      // also for stars
      sizefield=numStar * sizeof(float);

      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numStar; n++) {
	p = sStar.getConstParticle(n);
	file.write((char*)&(p->metal),sizeof(float));
	
      }
      file.write((char*)&sizefield,sizeof(int));
    } 
   
    if(getVerbose()>1)
      cerr << "done!" << endl;

  }

  if(getVerbose()>1)
    cerr << "GadgetSnap: finished with " << filename << endl;
   
  
  string filename_units = filename+".units";
  
  ofstream file_units(filename_units.c_str());
  file_units << s->getDistanceUnits() << endl << s->getMassUnits() << endl  << s->getVelocityUnits() << endl << s->getDensityUnits() << endl << s->getEnergyUnits() << endl;
  
}

} // namespace siman
