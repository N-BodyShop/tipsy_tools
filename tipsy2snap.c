/*
 * Taken from treebi2snap and modified by trq.
 * Reads tipsy standard and writes snap.
 *
 * Unsure of provenance of treebi2snap.  I (trq) believe this is from
 * Dusan Keres.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>

#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <assert.h>

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

float masses[100];
int nmass,masscount[100];

double totMass=0.;
float Lambda,hubble,boxsize;
double unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity;
float etaold,t0,H0;
int startflag=1;
float mass_factor,length_factor,vel_factor;
float redshift,aex;

int load_header(FILE *outp);
int load_data(FILE *outp);
int write_snapshot();
int cosmounits();

int main(int argc, char **argv)
{
	if( argc != 3 ) {
	    fprintf(stderr,
		    "usage: tipsy2snap BoxSize(Mpc/h) Hubble_Param(0.01*H0) < infile > outfile\n");
	    exit(-1);
	}
	boxsize = atof(argv[1]);
	hubble = atof(argv[2]);

	load_header(stdin);
	load_data(stdin);

	write_snapshot();

	fprintf(stderr,"Ntot= %d  Ngas= %d  z= %g\n",NumPart,Ngas,header1.redshift);
	exit(0);
}

XDR xdrs;
struct dump header;

int
load_header(FILE *outp)
{
	int i;
	int NStar;

/* Read in tipsy header info */
	xdrstdio_create(&xdrs, outp, XDR_DECODE);
	if(xdr_header(&xdrs, &header) != 1) {
	    fprintf(stderr, "Bad header\n");
	    exit(-1);
	    }
	
	assert(header.ndim == 3);
	
	Nhot = 0;
	NumPart = header.nbodies;
	Ngas = header.nsph;
	NStar = header.nstar;
	if( NStar != 0 ) {
		fprintf(stderr,"No stars allowed-- converts only m,r,v\n");
		exit(-1);
	}
	aex = header.time;

	cosmounits();
	
/* Unit conversion from TREEBI to gadget standard */
	mass_factor = unit_Mass/1.989e43*hubble; // Convert to 10^10 M_o/h
	length_factor = unit_Length/3.085678e21*hubble ; // Convert to kpc/h
	/* Convert to km/s, include sqrt(a) factor from Gadget */
	vel_factor = unit_Velocity/1.e5*sqrt(aex); 
	fprintf(stderr,"factors (a=%g): m=%g l=%g v=%g\n",aex,mass_factor,length_factor,vel_factor);

/* Load info into gadget header */
	header1.npart[0] = Ngas;
	header1.npart[1] = header.ndark;
	header1.npart[2] = Nhot;
	for(i=3;i<6;i++) header1.npart[i] = 0;
	for(i=0;i<6;i++) header1.npartTotal[i] = header1.npart[i];
	for(i=0;i<6;i++) header1.mass[i] = 0.0;	/* masses will be
						   specifed on a per
						   particle basis */
	header1.time = aex;
	header1.redshift = 1.0/aex - 1.0;
	header1.flag_sfr = 1;	/* do you want star formation? */
	header1.flag_feedback = 1;	/* what sort of feedback? */
	header1.flag_cooling = 1;	/* do you want cooling? */
	header1.num_files = 1;	        /* single file snapshots */
	header1.BoxSize = boxsize*1.e3;
	header1.Omega0 = totMass; /* XXX needs updating below */
	header1.OmegaLambda = Lambda;
	header1.HubbleParam = hubble;
	header1.flag_stellarage = 1;
	header1.flag_metals = 1;

	return 0;
}

struct gas_particle *gas;
struct dark_particle *dark;
struct star_particle *star;

int
load_data(FILE *outp)
{
	int i;
	double newmass;

	if(!(gas=malloc((header.nsph)*sizeof(struct gas_particle))))    {
		fprintf(stderr,
			"failed to allocate memory for %d gas particles.\n",
			header.nsph);
		exit(-1);
	}
	if(!(dark=malloc((header.ndark)*sizeof(struct dark_particle))))    {
		fprintf(stderr,
			"failed to allocate memory for %d dark particles.\n",
			header.ndark);
		exit(-1);
	}
	if(!(star=malloc((header.nstar)*sizeof(struct star_particle))))    {
	    fprintf(stderr,
		    "failed to allocate memory for %d star particles.\n",
		    header.nstar);
		exit(-1);
	}

	for(i = 0; i < header.nsph; i++) {
	    xdr_gas(&xdrs, &(gas[i]));
	    gas[i].pos[0] = length_factor*(gas[i].pos[0]+0.5);
	    gas[i].pos[1] = length_factor*(gas[i].pos[1]+0.5);
	    gas[i].pos[2] = length_factor*(gas[i].pos[2]+0.5);
	    gas[i].vel[0] *= vel_factor;
	    gas[i].vel[1] *= vel_factor;
	    gas[i].vel[2] *= vel_factor;
	    gas[i].mass *= mass_factor;
	    newmass = gas[i].mass;
	    totMass += newmass;
	    }
	for(i = 0; i < header.ndark; i++) {
	    xdr_dark(&xdrs, &(dark[i]));
	    dark[i].pos[0] = length_factor*(dark[i].pos[0]+0.5);
	    dark[i].pos[1] = length_factor*(dark[i].pos[1]+0.5);
	    dark[i].pos[2] = length_factor*(dark[i].pos[2]+0.5);
	    dark[i].vel[0] *= vel_factor;
	    dark[i].vel[1] *= vel_factor;
	    dark[i].vel[2] *= vel_factor;
	    if(i == 0)
      fprintf(stderr, "Particle mass: %g\n", dark[0].mass);
	    dark[i].mass *= mass_factor;
	    newmass = dark[i].mass;
	    totMass += newmass;
	    }
	for(i = 0; i < header.nstar; i++) {
	    xdr_star(&xdrs, &(star[i]));
	    star[i].pos[0] = length_factor*(star[i].pos[0]+0.5);
	    star[i].pos[1] = length_factor*(star[i].pos[1]+0.5);
	    star[i].pos[2] = length_factor*(star[i].pos[2]+0.5);
	    star[i].vel[0] *= vel_factor;
	    star[i].vel[1] *= vel_factor;
	    star[i].vel[2] *= vel_factor;
	    star[i].mass *= mass_factor;
	    newmass = star[i].mass;
	    totMass += newmass;
	    }
	masscount[nmass] = NumPart;
	header1.Omega0 = totMass/mass_factor;
	header1.OmegaLambda = 1.0 - header1.Omega0;
	
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
      fprintf(stderr, "Particle mass: %g\n", dark[0].mass);

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
		if(k == 0) {
		    fwrite(&gas[n].pos[0], sizeof(float), 3, fd);
		    }
		else if(k == 1) {
		    fwrite(&dark[n].pos[0], sizeof(float), 3, fd);
		    }
		else {
		    assert(0);	/* incomplete implementation */
		    }
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
		if(k == 0) {
		    fwrite(&gas[n].vel[0], sizeof(float), 3, fd);
		    }
		else if(k == 1) {
		    fwrite(&dark[n].vel[0], sizeof(float), 3, fd);
		    }
		else {
		    assert(0);	/* incomplete implementation */
		    }
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

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
		if(k == 0) {
		    fwrite(&gas[n].mass, sizeof(float), 1, fd);
		    }
		else if(k == 1) {
		    fwrite(&dark[n].mass, sizeof(float), 1, fd);
		    }
		else {
		    assert(0);	/* incomplete implementation */
		    }
	      pc_new++;
	    }
	}
      SKIP;

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
		fwrite(&gas[n].temp, sizeof(float), 1, fd);
		pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
		fwrite(&gas[n].rho, sizeof(float), 1, fd);
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

int
cosmounits()
{
    double Pi=3.14159265358979323846;
    double km=1.E5;
    double Mpc=3.085678e24;
#if 0
    double m_p=1.6726231E-24;     /* proton mass */
    double k_B=1.380622E-16;      /* Boltzman constant */
#endif

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

        fprintf(stderr,"COSMO PARAMS:  L=%g h^-1Mpc, h=%g, Omega=%g\n",
		boxsize,hubble,totMass);
        fprintf(stderr,"UNITS: T=%g rho=%g L=%g M=%g v=%g\n",unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity);

        return 0;
}
