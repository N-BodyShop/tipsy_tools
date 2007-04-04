//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// GADGETFILE.CPP
//
// NOTE - this code could be optimized a lot more!
// e.g. moving endianness transforms into separate block
// when loading to avoid endless branches.

#include "siman.hpp"
#include "endian.hpp"

#define SKIP fread(&dummy, sizeof(dummy), 1, fd); 

using namespace std;

// CONSTRUCTOR / DESTRUCTORS

CGadgetFile::~CGadgetFile() {
  
  if(Id!=NULL) free((void*)Id);
  
}

CGadgetFile::CGadgetFile(const char *fname_in, const bool swapEndian_in) {
  strcpy(fname,fname_in);
  
  Id=NULL;
  swapEndian = swapEndian_in;  // N.B. deprecated - this is autodetected in Load()
  load();
}

CGadgetFile::CGadgetFile(const char *path, const char *snapshot, const int snapid, const bool swapEndian_in) {

    sprintf(fname,  "%s/%s_%03d", path, snapshot,snapid);
    
    swapEndian = swapEndian_in; // N.B. deprecated - this is autodetected in Load()
    load();
}


// SIMPLE VALUE RETURNS


#define CHECKLENGTHFIELD(expect,section) fread(&checklen,sizeof(int),1,fd); if(swapEndian) { END4(&checklen); } if(checklen!=expect) { cerr << endl << "CGadgetFile: >>>Error in length field (" << section << ": says " << checklen << ", expected " << expect << ")" << endl; fileDebug(fd); }


// LOAD GADGET FILE

bool CGadgetFile::load() {
  
  // if(Loaded) return true;

  ENDINIT;

  // returns true for success
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  
  unsigned int checklen;

  int files = 1; // for now can not handle spanned files

  sprintf(buf,"%s",fname);
  
  if(!(fd=fopen(buf,"r")))
    {
      
      fprintf(stderr,"CGadgetFile: Error opening file `%s`\n",buf);
      perror("CGadgetFile");
      
      return false;
    }

  fprintf(stderr,"CGadgetFile: reading `%s' ...\n",buf); fflush(stderr);


  fread(&dummy, sizeof(dummy), 1, fd);
  if(dummy!=sizeof(header)) {
    END4(&dummy);
    if(dummy==sizeof(header)) {
      fprintf(stderr,"CGadgetFile: -> looks like endianness is wrong; auto-flipping.\n");
      swapEndian=true;
    } else {
      fprintf(stderr,"CGadgetFile: -> sorry, I don't understand this file.\n");
      return false;
    }
  } else {
    swapEndian = false;
  }

  fread(&header, sizeof(header), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);

      
  if(swapEndian==1) {
    swapDataEndian(&header, 4, 6*4);           // 6 integers
    swapDataEndian(&header.mass, 8, 8*8);     // 8 doubles
    swapDataEndian(&header.flag_sfr, 4, 40);  // 10 more integers
    swapDataEndian(&header.BoxSize,8, 4*8);   // 4 more doubles
  }
     
  fprintf(stderr,"CGadgetFile: z=%.1f; t=%.1f; h=%.2f; L=%.1f\n",header.redshift,header.time,header.HubbleParam,header.BoxSize);
      
  fprintf(stderr,"CGadgetFile: n_part - G:%d DM:%d S:%d\n",header.npart[0],header.npart[1],header.npart[4]);
  cerr << "CGadgetFile: header masses   G: " << header.mass[0] << "   DM: "<< header.mass[1] << "  S: " << header.mass[4] << endl;

  fprintf(stderr,"CGadgetFile: non-public flags - %s cooling   %s stellar ages   %s metallcities\n",header.flag_cooling?"YES":"NO",header.flag_stellarage?"YES":"NO",header.flag_metals?"YES":"NO");

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
  if(siman::fileExists(fname_units))
	  fname_units = "gadget.units";
  
  ifstream file_units(fname_units.c_str());
  
  
  
  if(file_units.is_open()) {
    cerr << "CGadgetFile: Loading conversion units from " << fname_units;
    
    file_units >> lenUnits;
    file_units >> massUnits;
    file_units >> velUnits;
    file_units >> denUnits;
    file_units >> enUnits;
  } else {

    cerr << "CGadgetFile: Warning - using default gadget units system which assumes GADGET was run in COSMOLOGICAL mode" << endl;
    cerr << "             Avoid this warning by storing units system in gadget.units or " + (string)fname + ".units" << endl;
    
      // the following may be incorrect depending on configuration of GADGET
    
    
    lenUnits.set(units::len_kpc);
    lenUnits/=(header.HubbleParam*(header.redshift+1.));
    
    massUnits.set(units::mass_Msol);
    massUnits*=1.e10/(header.HubbleParam);
    
    velUnits.set("(km s^-1 a^1/2)");
    
    /*
      velUnits.set(units::vel_kmPerS);
      velUnits*=sqrt(1./(1+header.redshift));
    */
    
    denUnits = massUnits/(lenUnits*lenUnits*lenUnits);
    
    enUnits = units::CUnit(units::vel_kmPerS)*units::CUnit(units::vel_kmPerS);
  }

  om_m0 = header.Omega0;
  om_lam0 = header.OmegaLambda;
  
  allocateMemory();

  pc = 0;

  fprintf(stderr,"CGadgetFile: Reading position data");
  // fileDebug(fd);
  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"PosStart");
  for(k=0,pc_new=pc;k<6;k++)
    {
      printf(".");
      fflush(stdout);
      for(n=0;n<header.npart[k];n++)
	{
	  fread(&(pParticles[pc_new]->x), sizeof(float), 1, fd);
	  fread(&(pParticles[pc_new]->y), sizeof(float), 1, fd);
	  fread(&(pParticles[pc_new]->z), sizeof(float), 1, fd);
	  if(swapEndian) {
	    END4(&(pParticles[pc_new]->x));
	    END4(&(pParticles[pc_new]->y));
	    END4(&(pParticles[pc_new]->z));
	  }
	  pc_new++;
	}
    }
  // fileDebug(fd);
  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"PosEnd");
  printf("done!\n");

  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"VelStart");
  // fileDebug(fd);
  fprintf(stderr,"CGadgetFile: Reading velocity data");
  for(k=0,pc_new=pc;k<6;k++)
    {
      fprintf(stderr,".");
      fflush(stdout);
      for(n=0;n<header.npart[k];n++)
	{
	  fread(&(pParticles[pc_new]->vx), sizeof(float), 3, fd);
	  if(swapEndian) {
	    END4(&(pParticles[pc_new]->vx));
	    END4(&(pParticles[pc_new]->vy));
	    END4(&(pParticles[pc_new]->vz));
	  }
	  pc_new++;
	}
    }
  printf("done!\n");
  CHECKLENGTHFIELD(sizeof(float)*3*numParticles,"VelEnd");
    

  CHECKLENGTHFIELD(sizeof(int)*numParticles,"PID-Start");
  fprintf(stderr,"CGadgetFile: Skipping PID data (not supported)");
  int ignore;
  for(k=0,pc_new=pc;k<6;k++)
    {
      fprintf(stderr,".");
      for(n=0;n<header.npart[k];n++)
	{
	  fread(&ignore, sizeof(int), 1, fd);
	      
	  pc_new++;
	}
    }
  fprintf(stderr,"done!\n");
  CHECKLENGTHFIELD(sizeof(int)*numParticles,"PID-End");


  //  fileDebug(fd);
      
      
  if(ntot_withmasses>0) {
    CHECKLENGTHFIELD(sizeof(float)*ntot_withmasses,"Mass-Start");
    fprintf(stderr,"CGadgetFile: Reading mass data");
  }
      
      
  
  for(k=0, pc_new=pc; k<6; k++) {
    cerr << ".";
    int writetype =0;
    switch(k) {
    case 0:
      writetype=CParticle::gas;
      break;
    case 1:
      writetype=CParticle::dm;
      break;
    case 4:
      writetype=CParticle::star;
      break;
    }
    
    for(n=0;n<header.npart[k];n++)
      {
	pParticles[pc_new]->type=writetype;
	if(header.mass[k]==0) {
	  fread(&(pParticles[pc_new]->mass), sizeof(float), 1, fd);
	  if(swapEndian) { END4(&(pParticles[pc_new]->mass)); }
	}
	else {
	  pParticles[pc_new]->mass = header.mass[k];		
	}
	pc_new++;
	
      }
  }

  
  if(ntot_withmasses>0) {
    cerr << "done!" << endl;
    CHECKLENGTHFIELD(sizeof(float)*ntot_withmasses,"Mass-End");
  }
      

  cerr << "CGadgetFile: Reading gas data...";

      
      

  if(header.npart[0]>0)
    {
      cerr << "energy ";
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"U-Start");
      for(n=0, pc_sph=pc; n<header.npart[0];n++)
	{
	  fread(&pParticles[pc_sph]->u, sizeof(float), 1, fd);
	  if(swapEndian) {END4(&(pParticles[pc_sph]->u)); }
	  pc_sph++;
	}
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"U-End");
	  
      cerr << "density ";

      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Rho-Start");
      for(n=0, pc_sph=pc; n<header.npart[0];n++)
	{
	  fread(&(pParticles[pc_sph]->rho), sizeof(float), 1, fd);
	  if(swapEndian) { END4(&(pParticles[pc_sph]->rho)); }
	  pc_sph++;
	}
      CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Rho-End");
	  
      if(header.flag_cooling)
	{
	  cerr << "electrons ";
	      
	  CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Ne-Start");
	  for(n=0, pc_sph=pc; n<header.npart[0];n++)
	    {
	      fread(&(pParticles[pc_sph]->ne), sizeof(float), 1, fd);
	      if(swapEndian) 
		{ END4(&(pParticles[pc_sph]->ne)); 
		}
	      pc_sph++;
	    }
	  CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Ne-End");
	      
	}
      if(header.flag_metals) {
	    
	cerr << "metals ";
	CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Met-Start");
	for(n=0, pc_sph=pc; n<header.npart[0];n++) {
	      
	  fread(&(pParticles[pc_sph]->metal), sizeof(float), 1, fd);
	  if(swapEndian) { END4(&(pParticles[pc_sph]->metal)); }
		  
		  
	      pc_sph++;
	}
	  CHECKLENGTHFIELD(sizeof(float)*header.npart[0],"Met-End");
	      
      }
      cerr << "done!" << endl;
          
      TempFromU();
     

    } // if(gas in sim)
   fclose(fd);
   return true;
      
}
      
void CGadgetFile::fileDebug(FILE *fd) {
  fpos_t init_pos;
  fgetpos(fd,&init_pos);
  
  printf("[%X] ",ftell(fd));

  for(int n=0;n<20;n++) {
    unsigned char contents;
    fread(&contents,1,1,fd);
    printf("%02x ",(unsigned int)contents);
  }
  printf("\n");

  
  fsetpos(fd,&init_pos);
}



int CGadgetFile::Reorder(void)
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

  fprintf(stderr,"CGadgetFile: done.\n");

  free(Id);

  fprintf(stderr,"CGadgetFile: space for particle ID freed\n");
  */
	return 0;
}


void CGadgetFile::writeField(std::ofstream *file, char *buf, int size) {

  // FORTRAN-style field

  file->write((char*)&size,sizeof(size));
  file->write(buf,size);
  file->write((char*)&size,sizeof(size));
}


void CGadgetFile::nativeWrite(CSimSnap *s, string filename) {


  cerr << "CGadgetFile: writing " << filename << "..." << endl;
  

  // not sure how this ought to be determined, as it
  // is a gadget-only flag:

  bool cooling=true;
  bool metals=true;
  bool stellar_age=true;

  // open file

  ofstream file(filename.c_str(),ios::binary);
  int sizefield = sizeof(gadget_header);

  // count particles

  unsigned int numStar=0, numDM=0, numGas=0, numTot = s->getNumParticles();
  float massStar=0, massDM=0, massGas=0;
  bool oneMassStar = true, oneMassDM = true, oneMassGas = true;

  CParticle *p;
  unsigned int n;

  for(n=0; n<numTot; n++) {
    p=s->getParticle(n);
    switch(p->type) { 

    case CParticle::gas:
      numGas+=1;
      if(p->mass!=massGas && massGas!=0) oneMassGas=false;
      massGas = p->mass;
      break;

    case CParticle::dm:
      numDM+=1;
      if(p->mass!=massDM && massDM!=0) oneMassDM=false;
      massDM = p->mass;
      break;

    case CParticle::star:
      numStar+=1;
      if(p->mass!=massStar && massStar!=0) oneMassStar=false;
      massStar = p->mass;
      break;

    }
    s->releaseParticle(p);
  }


  // create virtual simulations for ease of ordering
  // (gadget files insist on specific ordering of particles)

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
  header.flag_stellarage = stellar_age;

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

  cerr << "CGadgetFile: writing position data...";

  sizefield = sizeof(float) * 3 * numTot;
  
  file.write((char*)&sizefield,sizeof(int));

  for(n=0; n<numTot; n++) {
    p =  sOrdered.getParticle(n);
    file.write((char*)&(p->x),sizeof(float));
    file.write((char*)&(p->y),sizeof(float));
    file.write((char*)&(p->z),sizeof(float));
    sOrdered.releaseParticle(p);
  }
  
  file.write((char*)&sizefield,sizeof(int));

  cerr << "done!" << endl << "CGadgetFile: writing velocity data...";


  file.write((char*)&sizefield,sizeof(int));
  
  for(n=0; n<numTot; n++) {
    p =  sOrdered.getParticle(n);
    file.write((char*)&(p->vx),sizeof(float));
    file.write((char*)&(p->vy),sizeof(float));
    file.write((char*)&(p->vz),sizeof(float));
    sOrdered.releaseParticle(p);
  }
  
  file.write((char*)&sizefield,sizeof(int));
  
  cerr << "done!" << endl << "CGadgetFile: writing PID data...";
  
  sizefield = numTot * sizeof(int);

  file.write((char*)&sizefield,sizeof(int));
  
  for(n=0; n<numTot; n++) {
    file.write((char*)&n,sizeof(int));
  }
  
  file.write((char*)&sizefield,sizeof(int));

  cerr << "done!" << endl;


  // Write masses, where required

  sizefield = ((oneMassDM?0:numDM) + (oneMassGas?0:numGas) +
	       (oneMassStar?0:numStar))*sizeof(float);

  if(sizefield>0) {
    
    cerr << "CGadgetFile: writing mass data... (" << (oneMassDM?"":"D") << (oneMassStar?"":"S") << (oneMassGas?"":"G") << ")";

    int writeMassMask = (oneMassDM?0:CParticle::dm) + (oneMassStar?0:CParticle::star) + (oneMassGas?0:CParticle::gas);
    
    file.write((char*)&sizefield,sizeof(int));
    
    
    for(n=0; n<numTot; n++) {
      p = s->getParticle(n);
      if((p->type & writeMassMask) > 0) {
	file.write((char*)&(p->mass),sizeof(float));
      }
      s->releaseParticle(p);
    }

    file.write((char*)&sizefield,sizeof(int));

    cerr << "done!" << endl;
  }

  if(numGas>0) {
    cerr << "CGadgetFile: writing gas data...";
    sizefield = numGas * sizeof(float);

    file.write((char*)&sizefield,sizeof(int));

    for(n=0; n<numGas; n++) {
      p = sGas.getParticle(n);
      file.write((char*)&(p->u),sizeof(float));
      sGas.releaseParticle(p);
    }

    file.write((char*)&sizefield,sizeof(int));

    file.write((char*)&sizefield,sizeof(int));

    for(n=0; n<numGas; n++) {
      p = sGas.getParticle(n);
      file.write((char*)&(p->rho),sizeof(float));
      sGas.releaseParticle(p);
    }

    file.write((char*)&sizefield,sizeof(int));


    
    if(cooling) {
      
      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getParticle(n);
	file.write((char*)&(p->ne),sizeof(float));
	sGas.releaseParticle(p);
      }
      file.write((char*)&sizefield,sizeof(int));

      // XXX I'm guessing what this next Gadget field is.  TD's IDL
      // script says "Nh"  --TRQ
      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getParticle(n);
	file.write((char*)&(p->nHp),sizeof(float));
	sGas.releaseParticle(p);
      }
      file.write((char*)&sizefield,sizeof(int));

      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getParticle(n);
	file.write((char*)&(p->nHp),sizeof(float));
	sGas.releaseParticle(p);
      }
      file.write((char*)&sizefield,sizeof(int));

      // Smoothing length
      file.write((char*)&sizefield,sizeof(int));
      float *hsml = sGas.getArray("smoothlength");
      if(hsml == NULL) {
	  cerr << "WARNING: no smoothinglength" << endl;
	  hsml = sGas.createArray("smoothlength", "SPH smoothing length");
	  }
      
      for(n=0; n<numGas; n++) {
	file.write((char*)&(hsml[n]),sizeof(float));
      }
      file.write((char*)&sizefield,sizeof(int));

      // SFR
      // XXX I assume this is star formation rate.  What is meant by
      // this and what are its units?  --TRQ

      file.write((char*)&sizefield,sizeof(int));
      float *sfr = sGas.getArray("SFR");
      if(sfr == NULL) {
	  cerr << "WARNING: no SFR" << endl;
	  sfr = sGas.createArray("SFR", "Star formation rate?");
	  }
      
      for(n=0; n<numGas; n++) {
	file.write((char*)&(sfr[n]),sizeof(float));
      }
      file.write((char*)&sizefield,sizeof(int));
    }

    // XXX How does gadget define ages? --TRQ
    if(stellar_age) {
	
	sizefield = numStar*sizeof(float);
	float *age = sStar.getArray("stellarage");
	if(age == NULL) {
	  cerr << "WARNING: no stellar age" << endl;
	  age = sStar.createArray("stellarage", "Stellar Age");
	  }
	file.write((char*)&sizefield,sizeof(int));
	for(n=0; n<numStar; n++) {
	    file.write((char*)&(age[n]),sizeof(float));
	    }
	file.write((char*)&sizefield,sizeof(int));
	}
	
    if(metals) {
	sizefield = (numGas+numStar)*sizeof(float);
	
      file.write((char*)&sizefield,sizeof(int));
      for(n=0; n<numGas; n++) {
	p = sGas.getParticle(n);
	file.write((char*)&(p->metal),sizeof(float));
	sGas.releaseParticle(p);
      }
      for(n=0; n<numStar; n++) {
	p = sStar.getParticle(n);
	file.write((char*)&(p->metal),sizeof(float));
	sStar.releaseParticle(p);
      }
      file.write((char*)&sizefield,sizeof(int));
    } 
   

    cerr << "done!" << endl;

  }

  cerr << "CGadgetFile: finished with " << filename << endl;
   
  
}

