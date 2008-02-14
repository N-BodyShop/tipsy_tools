/*
 * Unsure of provenance.  I (trq) believe this is from Dusan Keres.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>

#define int4byte int

struct io_header_1
{
  int4byte npart[6];           /*!< number of particles of each type in this fil
e */
  double mass[6];              /*!< mass of particles of each type. If 0, then t
he masses are explicitly
                                 stored in the mass-block of the snapshot file,
otherwise they are omitted */
  double time;                 /*!< time of snapshot file */
  double redshift;             /*!< redshift of snapshot file */
  int4byte flag_sfr;           /*!< flags whether the simulation was including s
tar formation */
  int4byte flag_feedback;      /*!< flags whether feedback was included (obsolet
e) */
  int4byte npartTotal[6];      /*!< total number of particles of each type in th
is snapshot. This can be
                                 different from npart if one is dealing with a m
ulti-file snapshot. */
  int4byte flag_cooling;       /*!< flags whether cooling was included  */
  int4byte num_files;          /*!< number of files in multi-file snapshot */
  double BoxSize;              /*!< box-size of simulation in case periodic boun
daries were used */
  double Omega0;               /*!< matter density in units of critical density
*/
  double OmegaLambda;          /*!< cosmological constant parameter */
  double HubbleParam;          /*!< Hubble parameter in units of 100 km/sec/Mpc
*/
  int4byte flag_stellarage;    /*!< flags whether the file contains formation ti
mes of star particles */
  int4byte flag_metals;        /*!< flags whether the file contains metallicity
values for gas and star
                                 particles */
  char fill[88];               /*!< fills to 256 Bytes */
}
header1;  /*!< holds header for snapshot files */

int     NumPart, Ngas, Nhot;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
} *P;

double  Time;

float masses[100];
int nmass,masscount[100];

double totMass=0.;
float Lambda,hubble,boxsize;
double unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity;
float etaold,t0,H0;
int startflag=1;
float mass_factor,length_factor,vel_factor;
float redshift,aex;

int main(int argc, char **argv)
{
	FILE *outp;

	if( argc != 4 ) {
		fprintf(stderr,"usage: treebi2snap infile BoxSize(Mpc/h) Hubble_Param(0.01*H0) > outfile\n");
		exit(-1);
	}
        if( (outp = fopen(argv[1],"r")) == NULL ) {
		fprintf(stderr,"Cannot open %s\n",argv[1]);
		exit(-1);
	}
	boxsize = atof(argv[2]);
	hubble = atof(argv[3]);

	load_header(outp);
	load_data(outp);

	write_snapshot();

	fprintf(stderr,"Ntot= %d  Ngas= %d  z= %g\n",NumPart,Ngas,header1.redshift);
	exit(0);
}

load_header(outp)
FILE *outp;
{
	int i;
	double oldmass,newmass;
	int NStar;
	float junk;
	char pname[80];

/* Read in TREEBI header info */
	Nhot = 0;
	fgets(pname,80,outp);
	sscanf(pname,"%d %d %d %d",&NumPart,&Ngas,&NStar,&Nhot);
	if( NStar != 0 ) {
		fprintf(stderr,"No stars allowed-- converts only m,r,v\n");
		exit(-1);
	}
	fgets(pname,80,outp); sscanf(pname," %g",&junk);
	fgets(pname,80,outp); sscanf(pname," %lg",&Time);

/* Read in mass info */
	oldmass = 0.;
	nmass = 0;
	if( Ngas == 0 ) nmass = 1;
	for( i=0; i<100; i++ ) masses[i] = 0.;
	for( i=0; i<NumPart; i++ ) {
		fscanf(outp,"%lg",&newmass);
		if( oldmass != newmass ) {
			masses[nmass] = newmass;
			masscount[nmass] = i;
			oldmass = newmass;
			nmass++;
		}
		totMass += newmass;
	}
	masscount[nmass] = NumPart;
	for(i=1; i<=nmass; i++ ) fprintf(stderr,"masses %d: %d %g\n",i,masscount[i],masses[i-1]);

	cosmopar(Time);

/* Unit conversion from TREEBI to gadget standard */
	mass_factor = unit_Mass/1.989e43*hubble; // Convert to 10^10 M_o/h
	length_factor = unit_Length/3.085678e21*hubble ; // Convert to kpc/h
	/* Convert to km/s, include sqrt(a) factor from Gadget */
	vel_factor = unit_Velocity/1.e5*sqrt(aex); 
	fprintf(stderr,"factors (a=%g): m=%g l=%g v=%g\n",aex,mass_factor,length_factor,vel_factor);

/* Load info into gadget header */
	header1.npart[0] = Ngas;
	header1.npart[1] = NumPart-Ngas-Nhot;
	header1.npart[2] = Nhot;
	for(i=3;i<6;i++) header1.npart[i] = 0;
	for(i=0;i<6;i++) header1.npartTotal[i] = header1.npart[i];
	for(i=0;i<6;i++) header1.mass[i] = mass_factor*masses[i];
	header1.time = aex;
	header1.redshift = redshift;
	header1.flag_sfr = 1;	/* do you want star formation? */
	header1.flag_feedback = 1;	/* what sort of feedback? */
	header1.flag_cooling = 1;	/* do you want cooling? */
	header1.num_files = 1;	        /* single file snapshots */
	header1.BoxSize = boxsize*1.e3;
	header1.Omega0 = totMass;
	header1.OmegaLambda = Lambda;
	header1.HubbleParam = hubble;
	header1.flag_stellarage = 1;
	header1.flag_metals = 1;

	return 0;
}

load_data(outp)
FILE *outp;
{
	int i;
	float m;
	char pname[80];

	if(!(P=malloc((NumPart+1)*sizeof(struct particle_data))))    {
		fprintf(stderr,"failed to allocate memory for %d particles.\n",NumPart);
		exit(-1);
	}

	rewind(outp);
	fgets(pname,80,outp); fgets(pname,80,outp); fgets(pname,80,outp);
	for( i=1; i<=NumPart; i++ ) fscanf(outp,"%g ",&m);
	for( i=1; i<=NumPart; i++ ) fscanf(outp,"%g %g %g",&P[i].Pos[0],&P[i].Pos[1],&P[i].Pos[2]);
	for( i=1; i<=NumPart; i++ ) fscanf(outp,"%g %g %g",&P[i].Vel[0],&P[i].Vel[1],&P[i].Vel[2]);
//	for( i=1; i<=128; i++ ) fprintf(stderr,"%g-%g ",(P[i].Pos[0]+0.5)-1.*(i-1)/128,P[i].Vel[0]/((P[i].Pos[0]+0.5)-1.*(i-1)/128));
	for( i=1,m=0.; i<=NumPart; i++ ) {
		P[i].Pos[0] = length_factor*(P[i].Pos[0]+0.5);
		P[i].Pos[1] = length_factor*(P[i].Pos[1]+0.5);
		P[i].Pos[2] = length_factor*(P[i].Pos[2]+0.5);
		P[i].Vel[0] *= vel_factor;
		P[i].Vel[1] *= vel_factor;
		P[i].Vel[2] *= vel_factor;
		if( P[i].Vel[1] > m ) m=P[i].Vel[1];
	}
//	fprintf(stderr,"\n\n");
//	for( i=1; i<=128; i++ ) fprintf(stderr,"%g-%g ",(P[i].Pos[0])-length_factor*(i-1)/128,P[i].Vel[0]/(P[i].Pos[0]-length_factor*(i-1)/128));

	return 0;
}


/* this routine loads particle data into Gadget's default
 * binary file format.
 */
int write_snapshot()
{
  FILE *fd;
  int    i,k,dummy;
  int    n,pc,pc_new,pc_sph;
	float zero=0.;
	int files=1;
	int idnum=0;

#define SKIP fwrite(&dummy, sizeof(dummy), 1, fd);

	fd = stdout;
	dummy = sizeof(header1);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      fprintf(stderr,"outputting...");

      SKIP;
      fwrite(&header1, sizeof(header1), 1, fd);
      SKIP;
      dummy=0;
      for(k = 0; k < 6; k++){
	dummy = dummy + 3*header1.npart[k]*sizeof(float);
      }      
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fwrite(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fwrite(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      dummy=dummy/3;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      idnum++;
	      fwrite(&idnum, sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fwrite(&zero, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fwrite(&zero, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fwrite(&zero, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;
	}

    }
	fprintf(stderr,"done.\n");

	return 0;
}


cosmounits()
{
    double Pi=3.14159265358979323846;
    double km=1.E5;
    double Mpc=3.085678e24;
    double m_p=1.6726231E-24;     /* proton mass */
    double k_B=1.380622E-16;      /* Boltzman constant */

	if( 1.-totMass > 1.e-6 ) {
		fprintf(stderr,"Setting Lambda = %g\n",1.-totMass);
		Lambda = 1.-totMass;
	}
	else Lambda = 0.0;

    H0=sqrt(8*Pi/3);
	t0 = 2./(3*H0);

    unit_Time=H0*Mpc/(100*hubble*km);
    unit_Density=1.8791E-29*hubble*hubble;
    unit_Length=boxsize*Mpc/hubble;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length;
    unit_Velocity=unit_Length/unit_Time;

	fprintf(stderr,"COSMO PARAMS:  L=%g h^-1Mpc, h=%g, Omega=%g, t=%g\n",boxsize,hubble,totMass,Time);
	fprintf(stderr,"UNITS: T=%g rho=%g L=%g M=%g v=%g\n",unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity);

	return 0;
}

cosmopar(t)
float t;
{
	float t1;
	float tol=2.e-6;
	float a0,astar,eta,etalast;
	float f,fprime;
	float hub,aexhub,aex3;
	int it;

	if( startflag ) {
		cosmounits();
		startflag = 0;
		etaold = 1.;
	}
	if( fabs(totMass - 1.0) < tol ) {
		t1 = t/t0;
		aex3 = t1*t1;
		aex = pow((aex3),1./3.);
		hub = 2.0/3.0/t;
		aexhub = aex*hub;
		redshift = 1./aex - 1.;
	}
	else if( totMass < 1.0 ) {
		if( Lambda > 0.0 ) {	/* Flat, low-density universe */
                eta = sqrt(1.-totMass)*1.5*H0*t;
                aex = pow(sqrt(totMass/(1.-totMass))*sinh(eta),2./3);
		aex3 = aex*aex*aex;
		hub = H0*sqrt(totMass/aex3+Lambda);
		aexhub = aex*hub;
		redshift = 1./aex - 1.;
/*	fprintf(stderr,"cosmopar: %g %g %g %g %g\n",t,aex,hub,totMass,Lambda);*/
		}
		else {			/* Open universe */
        	a0=1./H0/sqrt(1.-totMass);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(sinh(eta)-eta)-t;
        		fprime=astar*(cosh(eta)-1.);
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in cosmopar %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		aex = astar*(cosh(eta)-1.)/a0;
		aex3 = aex*aex*aex;
		etaold = eta;
		redshift = 1./aex - 1.;
		hub=H0*(1.+redshift)*sqrt(1.+totMass*redshift);
		aexhub = aex*hub;
		}
	}
	else if( totMass > 1.0 ) {
        	a0=1./H0/sqrt(totMass-1.);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(eta-sin(eta))-t;
        		fprime=astar*(1.-cos(eta));
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in cosmopar %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		aex = astar*(1.-cos(eta))/a0;
		aex3 = aex*aex*aex;
		etaold = eta;
		redshift = 1./aex - 1.;
		hub=H0*(1.+redshift)*sqrt(1.+totMass*redshift);
		aexhub = aex*hub;
	}

	return 0;
}

