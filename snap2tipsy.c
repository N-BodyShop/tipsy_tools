#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7
#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define  H_MASSFRAC    0.76


double  UnitLength_in_cm,
        UnitMass_in_g,
               UnitVelocity_in_cm_per_s,
               UnitTime_in_s,
               UnitTime_in_Megayears,
               UnitDensity_in_cgs,
               UnitPressure_in_cgs,
               UnitCoolingRate_in_cgs,
               UnitEnergy_in_cgs,
               G,
               Hubble;




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


#define Real float
#define MAXDIM 3

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} gasp;


struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} darkp;


struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} starp;


struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;
struct dump header ;



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char output_fname[200], input_fname[200], basename[200], hsmlfile[200];
  int  type, files;
  
  if(argc != 2 && argc != 3 && argc !=4) 
    {
      fprintf(stderr,"usage: snap2tipsy <snapshotfile> [N_files]\n");
      exit(-1);
    }

  strcpy(basename, argv[1]);
  hsmlfile[0]= 0;

  if(argc >= 3 ) 
    {
      files = atoi(argv[2]);	/* # of files per snapshot */
    }
  else 
    files = 1;

  sprintf(input_fname, "%s", basename);
  sprintf(output_fname, "%s.bin", basename);

  printf("loading gas particles...\n");

  load_snapshot(input_fname, files, type=0);

  filter_gas();

  printf("filtered gas particles...\n");

  output_tipsy_gas(output_fname);

  free_memory(); 

  printf("loading dark matter particles...\n");

  load_snapshot(input_fname, files, type=1);

  filter_darkmatter();

  printf("filtered dark matter particles...\n");

  output_tipsy_dark(output_fname);

  free_memory();

  printf("loading star particles...\n");

  load_snapshot(input_fname, files, type=4);

  filter_star();

  output_tipsy_star(output_fname);

  free_memory();

  update_tipsy_header(output_fname);

  exit(0);
}


filter_gas()
{
  int  i;

  double Xh=H_MASSFRAC;  /* mass fraction of hydrogen */
  double HubbleParam= 0.7;
  double MeanWeight, rhomean;

  /* first, let's compute the temperature */

  UnitLength_in_cm = 3.085678e21; 
  UnitMass_in_g    = 1.989e43;
  UnitVelocity_in_cm_per_s = 1e5;

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears= UnitTime_in_s/SEC_PER_MEGAYEAR;
  G=GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
  UnitDensity_in_cgs=UnitMass_in_g/pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs=UnitMass_in_g/UnitLength_in_cm/pow(UnitTime_in_s,2);
  UnitCoolingRate_in_cgs=UnitPressure_in_cgs/UnitTime_in_s;
  UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
  Hubble = HUBBLE * UnitTime_in_s;



  for(i=1; i<=NumPart; i++)
    {
      MeanWeight= 4.0/(3*Xh+1+4*Xh* P[i].Ne) * PROTONMASS;
      
      P[i].Temp *= MeanWeight/BOLTZMANN * GAMMA_MINUS1 
       	            * UnitEnergy_in_cgs/ UnitMass_in_g;
      
      /* temp now in kelvin */
    }

  rhomean= header1.Omega0*3*Hubble*Hubble/(8*M_PI*G);


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

  if(!(outfile = fopen(output_fname,"w")))
    {
      printf("can't open file `%s'\n", outfile);
      exit(0);
    }
  
  printf("writing gas particles to tipsy file...\n");

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
	  gasp.vel[0] = P[i].Vel[0];
	  gasp.vel[1] = P[i].Vel[1];
	  gasp.vel[2] = P[i].Vel[2];
	  gasp.temp =   P[i].Temp;
	  gasp.hsmooth = P[i].Hsml;
	  gasp.rho = P[i].Rho;
	  gasp.metals = 0.0;
	  gasp.phi = 0.0;

	  fwrite(&gasp, sizeof(struct gas_particle), 1, outfile) ;
	}
    }

  fclose(outfile);
  
  printf("done.\n");

  return 0;
}

int output_tipsy_star(char *output_fname)
{
  int   i;
  FILE *outfile;

  if(!(outfile = fopen(output_fname,"a")))
    {
      printf("can't open file `%s'\n", outfile);
      exit(0);
    }
  
  printf("writing star particles to tipsy file...\n");

  for(i=1; i<=NumPart; i++) 
    {
      if(P[i].Flag)
	{
	  starp.mass =   P[i].Mass;
	  /*	      printf("star mass %f ...\n",P[i].Mass); */
	  starp.pos[0] = P[i].Pos[0];
	  starp.pos[1] = P[i].Pos[1];
	  starp.pos[2] = P[i].Pos[2];
	  starp.vel[0] = P[i].Vel[0];
	  starp.vel[1] = P[i].Vel[1];
	  starp.vel[2] = P[i].Vel[2];
	  starp.metals =   0.0;
	  starp.tform = 0.0;
	  starp.eps = 0.0;
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


int update_tipsy_header(char *output_fname)
{
  int   i;
  FILE *outfile;

  printf("updating tipsy header...\n");

  if(!(outfile = fopen(output_fname,"r+")))
    {
      printf("can't open file `%s'\n", outfile);
      exit(0);
    }

  printf("header.nsph= %d...\n",header.nsph);
  printf("header.ndark= %d...\n",header.ndark);
  printf("header.nstar= %d...\n",header.nstar);

  fwrite(&header, sizeof(header), 1, outfile);
  fclose(outfile);

}


int output_tipsy_dark(char *output_fname)
{
  int   i;
  FILE *outfile;

  printf("writing dark matter particles to tipsy file...\n");

  if(!(outfile = fopen(output_fname,"a")))
    {
      printf("can't open file `%s'\n", outfile);
      exit(0);
    }
  
  for(i=1; i<=NumPart; i++) 
    {
      if(P[i].Flag)
	{
	  darkp.mass =   header1.mass[1];
	  darkp.pos[0] = P[i].Pos[0];
	  darkp.pos[1] = P[i].Pos[1];
	  darkp.pos[2] = P[i].Pos[2];
	  darkp.vel[0] = P[i].Vel[0];
	  darkp.vel[1] = P[i].Vel[1];
	  darkp.vel[2] = P[i].Vel[2];
	  darkp.eps = 0.;
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
int load_snapshot(char *fname, int files, int type)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

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

      fprintf(stderr,"reading `%s' ...\n",buf);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

  fprintf(stderr,"header1.npart[k]=%d \n",header1.npart[type]);

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
		  fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
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
		  fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
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
		  fread(&Id[pc_new], sizeof(int), 1, fd);
		  pc_new++;
		}
	    }
	  else
	    fseek(fd, sizeof(int)*header1.npart[k], SEEK_CUR);
	}
      SKIP;

      
      /*----------------------------------------------------------*/

      printf("doing masses stuff \n");


      if(ntot_withmasses>0)
	SKIP;

      /*      fseek(fd, sizeof(float)*ntot_withmasses, SEEK_CUR); */
      
	for(k=0, pc_new=pc; k<6; k++)
	{
	  if (k==type) 

	    {

	for(n=0;n<header1.npart[k];n++)

	{
	
	if(header1.mass[k]==0) 

	fread(&P[pc_new].Mass, sizeof(float), 1, fd);

	else

	      P[pc_new].Mass= header1.mass[k];
	      pc_new++;

	      }
	    }
	else
	  
    if(header1.mass[k]==0)

      {
	  for(n=0;n<header1.npart[k];n++)
	    {
	  fseek(fd,sizeof(float),SEEK_CUR);
	  }
      }
	      }
      

      if(ntot_withmasses>0)
	SKIP;

      /*----------------------------------------------------------*/
      

      if(header1.npart[0]>0 && type==0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Temp, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
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
		  fread(&P[pc_sph].Hsml, sizeof(float), 1, fd);
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
		  fread(&P[pc_sph].Hsml, sizeof(float), 1, fd);
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


      fclose(fd);
    }
  fprintf(stderr,"done.\n");

  Time= header1.time;
  Redshift= header1.time;
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
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


int free_memory(void)
{
  Id++;
  P++;
  free(Id);
  free(P);
}









  











