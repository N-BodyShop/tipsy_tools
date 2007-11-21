/* 
 * Obtained by trq from Rubert Croft via Tiziana De Mateo.
 * Modified significantly by trq.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <endian.h>

/*
 * In place swap of data
 */
void swapEndian(void *data, int size, int count)
{
    char *cdata = (char *) data;
    int iCount;
    
    for(iCount = 0; iCount < count; iCount++) {
	int i;
	char temp;
	
	for(i = 0; i < size/2; i++) {
	    temp = cdata[size-i-1];
	    cdata[size-i-1] = cdata[i];
	    cdata[i] = temp;
	    }
	cdata += size;
	}
    }
	
double G,
    Hubble; 			/* H_0 in 100 km/s/Mpc */

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

double softenings[6];		/* Gravitational softenings for each
				   particle type */

#define TYPE_GAS 0
#define TYPE_HALO 1
#define TYPE_DISK 2
#define TYPE_BULGE 3
#define TYPE_STAR 4
#define TYPE_BNDRY 5

int     NumPart, NumPartFiltered;


struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass, Rho, Temp, Ne, Hsml;  /* temp */
  int    Flag;
} *P;

int   *Id;


double  Time, Redshift;

#include "tipsydefs.h"

struct gas_particle gasp;
struct dark_particle darkp;
struct star_particle starp;

void free_memory(void);
void allocate_memory(void);
void load_snapshot(char *fname, int files, int type);
void update_tipsy_header(char *output_fname);
void filter_darkmatter();
void filter_star();
void filter_gas();
int output_tipsy_gas(char *output_fname);
int output_tipsy_star(char *output_fname, int bBH);
int output_tipsy_dark(char *output_fname, int type);
void set_unit();


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char output_fname[200], input_fname[200], basename[200];
  int  type, files;
  
  if(argc != 10 && argc != 11) 
    {
      fprintf(stderr,
	      "usage: snap2tipsy h eps_gas eps_dark eps_disk eps_bulge eps_star eps_bndry(BH) <snapshotfile> <tipsyfile> [N_files]\n");
      fprintf(stderr, "	where h is H_0 in 100 km/s/Mpc and softenings are in kpc/h\n");
      exit(-1);
    }

  Hubble = atof(argv[1]);
  for(type = 0; type < 6; type++) {
      softenings[type] = atof(argv[2 + type]);
      }
  
  strcpy(basename, argv[8]);
  strcpy(output_fname, argv[9]);

  if(argc == 11 ) 
    {
      files = atoi(argv[10]);	/* # of files per snapshot */
    }
  else 
    files = 1;

  sprintf(input_fname, "%s", basename);

  printf("loading gas particles...\n");

  load_snapshot(input_fname, files, type=0);

  set_unit();

  filter_gas();
  output_tipsy_gas(output_fname);
  free_memory(); 

  printf("loading dark matter particles...\n");
  load_snapshot(input_fname, files, type=1);
  filter_darkmatter();
  output_tipsy_dark(output_fname, type);
  free_memory();

  /* Handle funny component particles from Gadget.  Here we assume
     that these types are collisionless only.  */
  printf("loading disk particles... ");
  load_snapshot(input_fname, files, type=2);
  filter_darkmatter();
  printf("and changing them to dark...\n");
  output_tipsy_dark(output_fname, type);
  free_memory();

  printf("loading bulge particles... ");
  load_snapshot(input_fname, files, type=3);
  filter_darkmatter();
  printf("and changing them to dark...\n");
  output_tipsy_dark(output_fname, type);
  free_memory();

  printf("loading star particles... ");
  load_snapshot(input_fname, files, type=4);
  filter_star();
  output_tipsy_star(output_fname, 0);
  free_memory();

  printf("loading black hole particles... ");
  load_snapshot(input_fname, files, type=5);
  filter_darkmatter();
  printf("and changing them to star...\n");
  output_tipsy_star(output_fname, 1);
  free_memory();

  update_tipsy_header(output_fname);

  exit(0);
}

struct units 
{
  double Length_in_cm; 
  double Mass_in_g;
  double Velocity_in_cm_per_s;
  double Time_in_s;
  double Time_in_Megayears;
  double Density_in_cgs;
  double Pressure_in_cgs;
  double CoolingRate_in_cgs;
  double Energy_in_cgs;
    double Natural_vel_in_cgs;		/* sqrt(M/r) with G=1 */
    } Unit;

void set_unit()
{
    /*
     * Assume length in kpc and mass in 1e10 Solarmasses
     */
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_KPC  3.085678e21
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7
#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define  H_MASSFRAC    0.76

  Unit.Length_in_cm = CM_PER_KPC/Hubble; 
  Unit.Mass_in_g    = SOLAR_MASS*1.0e10/Hubble;
  Unit.Velocity_in_cm_per_s = 1e5;
  Unit.Time_in_s= Unit.Length_in_cm / Unit.Velocity_in_cm_per_s;
  Unit.Time_in_Megayears= Unit.Time_in_s/SEC_PER_MEGAYEAR;
  Unit.Density_in_cgs=Unit.Mass_in_g/pow(Unit.Length_in_cm,3);
  Unit.Pressure_in_cgs=Unit.Mass_in_g/Unit.Length_in_cm/pow(Unit.Time_in_s,2);
  Unit.CoolingRate_in_cgs=Unit.Pressure_in_cgs/Unit.Time_in_s;
  Unit.Energy_in_cgs=Unit.Mass_in_g * pow(Unit.Length_in_cm,2) / pow(Unit.Time_in_s,2);
  Unit.Natural_vel_in_cgs = sqrt(GRAVITY * Unit.Mass_in_g/Unit.Length_in_cm);
  
  fprintf(stdout, "dMsolUnit is %g; dKpcUnit is %g;\n", 1e10/Hubble, 1.0/Hubble);
  fprintf(stdout, "velocity is in units of %g km/s; time is in units of %g Gyr\n",
	  Unit.Natural_vel_in_cgs/1e5,
	  Unit.Length_in_cm/Unit.Natural_vel_in_cgs/(1e9*SEC_PER_YEAR));
  
    }

void
filter_gas()
{
  int  i;

  double Xh=H_MASSFRAC;  /* mass fraction of hydrogen */
  double MeanWeight, rhomean;

  /* first, let's compute the temperature */

  Hubble = HUBBLE * Unit.Time_in_s;



  for(i=1; i<=NumPart; i++)
    {
      MeanWeight= 4.0/(3*Xh+1+4*Xh* P[i].Ne) * PROTONMASS;
      
      P[i].Temp *= MeanWeight/BOLTZMANN * GAMMA_MINUS1 
       	            * Unit.Energy_in_cgs/ Unit.Mass_in_g;
      
      /* temp now in kelvin */
    }

  if(header1.Omega0 != 0.0) {
      rhomean= header1.Omega0*3*Hubble*Hubble/(8*M_PI*G);
      }
  else {
      rhomean = 1.0;
      }

  NumPartFiltered= 0;

  for(i=1; i<= NumPart; i++)
    {

      /* these are just silly values for the moment, so that 
	 everything will be selected */

      if(P[i].Temp > -100000  && P[i].Rho/rhomean > -1000.)
	{
	  P[i].Flag= 1; /* take it */
	  NumPartFiltered++;
	}
      else
	P[i].Flag= 0; /* discard */
    }

}


void
filter_darkmatter()
{
  int  i;

  NumPartFiltered= 0;
  for(i=1; i<= NumPart; i++)
    {
      P[i].Flag= 1; /* keep- set this to 0 to discard */
      NumPartFiltered++; 
    }
}

void
filter_star()
{
  int  i;

  NumPartFiltered= 0;
  for(i=1; i<= NumPart; i++)
    {
      P[i].Flag= 1; /* keep- set this to 0 to discard */
      NumPartFiltered++; 
    }
}




int output_tipsy_gas(char *output_fname)
{
  int   i;
  FILE *outfile;
  double vscale = Unit.Velocity_in_cm_per_s/Unit.Natural_vel_in_cgs;

  if(!(outfile = fopen(output_fname,"w")))
    {
      printf("can't open file `%s'\n", output_fname);
      exit(0);
    }
  
  printf("velocity scaling: %g\n", vscale);

  /* Load tipsy header */
  header.time = header1.time;
  header.nbodies = NumPartFiltered; 
  header.ndim =  3;
  header.nsph =  NumPartFiltered;
  header.ndark = 0;
  header.nstar = 0;  /*  header1.npartTotal[4]; */

  fwrite(&header, sizeof(header), 1, outfile);
  
  for(i=1; i<=NumPart; i++) 
    {
      if(P[i].Flag)
	{
	  gasp.mass =   P[i].Mass;
	  /*      printf("gas mass %f ...\n",P[i].Mass); */
	  gasp.pos[0] = P[i].Pos[0];
	  gasp.pos[1] = P[i].Pos[1];
	  gasp.pos[2] = P[i].Pos[2];
	  gasp.vel[0] = P[i].Vel[0]*vscale;
	  gasp.vel[1] = P[i].Vel[1]*vscale;
	  gasp.vel[2] = P[i].Vel[2]*vscale;
	  gasp.temp =   P[i].Temp;
	  gasp.hsmooth = softenings[TYPE_GAS];
	  gasp.rho = P[i].Rho;
	  gasp.metals = 0.0;
	  gasp.phi = 0.0;

	  fwrite(&gasp, sizeof(struct gas_particle), 1, outfile) ;
	}
    }

  fclose(outfile);
  
  return 0;
}

int output_tipsy_star(char *output_fname, int bBH)
{
  int   i;
  FILE *outfile;
  double vscale = Unit.Velocity_in_cm_per_s/Unit.Natural_vel_in_cgs;

  if(!(outfile = fopen(output_fname,"a")))
    {
      printf("can't open file `%s'\n", output_fname);
      exit(0);
    }
  

  for(i=1; i<=NumPart; i++) 
    {
      if(P[i].Flag)
	{
	  starp.mass =   P[i].Mass;
	  starp.pos[0] = P[i].Pos[0];
	  starp.pos[1] = P[i].Pos[1];
	  starp.pos[2] = P[i].Pos[2];
	  starp.vel[0] = P[i].Vel[0]*vscale;
	  starp.vel[1] = P[i].Vel[1]*vscale;
	  starp.vel[2] = P[i].Vel[2]*vscale;
	  starp.metals =   0.0;		/* XXX this is incomplete */
	  if(bBH) {
	      /* Negative tForm signals black hole to GASOLINE */
	      starp.tform = -1;
	      starp.eps = softenings[TYPE_BNDRY];
	      }
	  else {
	      starp.tform = 0.0;
	      starp.eps = softenings[TYPE_STAR];
	      }
	  starp.phi = 0.0;

	  fwrite(&starp, sizeof(struct star_particle), 1, outfile) ;
	}
    }

  fclose(outfile);

  header.nbodies += NumPartFiltered; 
  header.nstar += NumPartFiltered;
  
  printf("done.\n");

  return 0;
}


void update_tipsy_header(char *output_fname)
{
  FILE *outfile;

  printf("updating tipsy header...\n");

  if(!(outfile = fopen(output_fname,"r+")))
    {
      printf("can't open file `%s'\n", output_fname);
      exit(0);
    }

  printf("header.nsph= %d...\n",header.nsph);
  printf("header.ndark= %d...\n",header.ndark);
  printf("header.nstar= %d...\n",header.nstar);

  fwrite(&header, sizeof(header), 1, outfile);
  fclose(outfile);

}


int output_tipsy_dark(char *output_fname, int type)
{
  int   i;
  FILE *outfile;
  double vscale = Unit.Velocity_in_cm_per_s/Unit.Natural_vel_in_cgs;

  if(!(outfile = fopen(output_fname,"a")))
    {
      printf("can't open file `%s'\n", output_fname);
      exit(0);
    }
  
  for(i=1; i<=NumPart; i++) 
    {
      if(P[i].Flag)
	{
	  darkp.mass =   P[i].Mass;
	  darkp.pos[0] = P[i].Pos[0];
	  darkp.pos[1] = P[i].Pos[1];
	  darkp.pos[2] = P[i].Pos[2];
	  darkp.vel[0] = P[i].Vel[0]*vscale;
	  darkp.vel[1] = P[i].Vel[1]*vscale;
	  darkp.vel[2] = P[i].Vel[2]*vscale;
	  darkp.eps = softenings[type];
	  darkp.phi = 0.;
	  
	  fwrite(&darkp, sizeof(struct dark_particle), 1, outfile) ;
	}
    }

  fclose(outfile);

  header.nbodies += NumPartFiltered; 
  header.ndark += NumPartFiltered;

  printf("done.\n");

  return 0;
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
void load_snapshot(char *fname, int files, int type)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;
  int nread;
  int swap = 0;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  fprintf(stderr,"can't open file `%s`\n",buf);
	  exit(0);
	}

      fread(&dummy, sizeof(dummy), 1, fd);
      if(dummy!=sizeof(header1)) {
	  swap = 1;
	  fprintf(stderr, "Trying endian swap\n");
	  swapEndian(&dummy, sizeof(dummy), 1);
	  assert(dummy == sizeof(header1));
	  }
      
      nread = fread(&header1, sizeof(header1), 1, fd);
      if(nread != 1) {
	fprintf(stderr, "Bad header read of %s\n", buf);
	exit(-1);
	}
      if(swap) {
	  swapEndian(&header1, 4, 6);           // 6 integers
	  swapEndian(&header1.mass, 8, 8);     // 8 doubles
	  swapEndian(&header1.flag_sfr, 4, 10);  // 10 more integers
	  swapEndian(&header1.BoxSize,8, 4);   // 4 more doubles
	  }
      
      fprintf(stderr, "BoxSize: %g\n", header1.BoxSize);
      fprintf(stderr, "Hubble parameter: %g\n", header1.HubbleParam);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  NumPart= header1.npart[type];
	}
      else
	{
	  NumPart= header1.npartTotal[type];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  if(k==type)
	    {
	      for(n=0;n<header1.npart[k];n++)
		{
		  nread = fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
		  assert(nread == 3);
		  if(swap)
		      swapEndian(&P[pc_new].Pos[0], sizeof(float), 3);
		  pc_new++;
		}
	    }
	  else
	    fseek(fd, sizeof(float)*3*header1.npart[k], SEEK_CUR);
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  if(k==type)
	    {
	      for(n=0;n<header1.npart[k];n++)
		{
		  nread = fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
		  assert(nread == 3);
		  if(swap)
		      swapEndian(&P[pc_new].Vel[0], sizeof(float), 3);
		  pc_new++;
		}
	    }
	  else
	    fseek(fd, sizeof(float)*3*header1.npart[k], SEEK_CUR);
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  if(k==type)
	    {
	      for(n=0;n<header1.npart[k];n++)
		{
		  nread = fread(&Id[pc_new], sizeof(int), 1, fd);
		  assert(nread == 1);
		  if(swap) 
		      swapEndian(&Id[pc_new], sizeof(int), 1);
		  pc_new++;
		}
	    }
	  else
	    fseek(fd, sizeof(int)*header1.npart[k], SEEK_CUR);
	}
      SKIP;

      
      /*----------------------------------------------------------*/

      if(ntot_withmasses>0)
	SKIP;

      /*      fseek(fd, sizeof(float)*ntot_withmasses, SEEK_CUR); */
      
	for(k=0, pc_new=pc; k<6; k++)
	{
	  if (k==type) {

	    for(n=0;n<header1.npart[k];n++) {
		if(header1.mass[k]==0) {
		    nread = fread(&P[pc_new].Mass, sizeof(float), 1, fd);
		    assert(nread == 1);
		    if(swap)
			swapEndian(&P[pc_new].Mass, sizeof(float), 1);
		    }
		else
		    P[pc_new].Mass= header1.mass[k];
		pc_new++;
		}
	      }
	  else
		if(header1.mass[k]==0)
		  {
		      for(n=0;n<header1.npart[k];n++) {
			  fseek(fd,sizeof(float),SEEK_CUR);
		      }
		  }
	    }
      

      if(ntot_withmasses>0)
	SKIP;

      /*----------------------------------------------------------*/
      /* Gas variables */

      if(header1.npart[0]>0 && type==0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      nread = fread(&P[pc_sph].Temp, sizeof(float), 1, fd);
	      if(nread != 1) {
		  if(pc_sph == pc) {
		      fprintf(stderr, "WARNING: No SPH temperatures\n");
		      return;
		      }
		  else {
		      fprintf(stderr, "Short read of temperatures\n");
		      exit(-1);
		      }
		  }
	      if(swap)
		  swapEndian(&P[pc_sph].Temp, sizeof(float), 1);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      nread = fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      if(nread != 1) {
		  if(pc_sph == pc || n == 1) {
		      fprintf(stderr, "WARNING: No SPH densities\n");
		      n = 0;
		      while(n < header1.npart[0]) {
			  P[pc_sph].Rho = 0.0; /* fill in with zeros */
			  n++;
			  pc_sph++;
			  }
		      return;
		      }
		  else {
		      fprintf(stderr, "Short read of densities at n = %d\n",
			      n);
		      exit(-1);
		      }
		  }
	      if(swap)
		  swapEndian(&P[pc_sph].Rho, sizeof(float), 1);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		    nread = fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		    assert(nread == 1);
		    if(swap)
			swapEndian(&P[pc_sph].Ne, sizeof(float), 1);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	

      /* now a dummy read of the neutral hydrogen into the hsm array */

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  nread = fread(&P[pc_sph].Hsml, sizeof(float), 1, fd);
		  assert(nread == 1);
		  if(swap)
		      swapEndian(&P[pc_sph].Hsml, sizeof(float), 1);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Hsml= 1.0;
		pc_sph++;
	      }
	


      /* now really read hsm into the hsm array */

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		    nread = fread(&P[pc_sph].Hsml, sizeof(float), 1, fd);
		    assert(nread == 1);
		    if(swap)
			swapEndian(&P[pc_sph].Hsml, sizeof(float), 1);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Hsml= 1.0;
		pc_sph++;
	      }
	
	}

      /*----------------------------------------------------------*/
      /* Extra variables */
      /* Unfortunately, we have no clue what could come here.
       * There could be POT, ACCEL, DTENTR, TSTP
       */

      fclose(fd);
    }

  Time= header1.time;
  Redshift= header1.time;
}




/* this routine allocates the memory for the 
 * particle data.
 */
void allocate_memory(void)
{
/*  fprintf(stderr,"allocating memory...\n");*/

  fprintf(stderr,"NumPart=%d \n",NumPart);

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

}


void free_memory(void)
{
  Id++;
  P++;
  free(Id);
  free(P);
}
