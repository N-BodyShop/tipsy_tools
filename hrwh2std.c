/*
 * Convert format of Heitmann, Ricker, Warren & Habib (2005 ApJS 160)
 * to tipsy format.
 */
#include <stdio.h>
#include <math.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

int
main(int argc, char *argv[])
{

    struct hrwh_part 
    {
	float x;
	float vx;
	float y;
	float vy;
	float z;
	float vz;
	float mass;
	int itag;
	} part;

    struct dark_particle dark;
    struct dump header;
    XDR xdrs;
    int npart;
    double z;
    double h0;
    double boxsize;
    double mconv;
    double xconv;
    double vconv;
    int i;
    
    if(argc != 5) {
	fprintf(stderr, "hrwh2std npart z boxsize(Mpc) H0(km/s/Mpc)\n");
	return 1;
	}

    npart = atoi(argv[1]);
    z = atof(argv[2]);
    boxsize = atof(argv[3]);
    h0 = atof(argv[4]);
    
    header.time = 1.0/(1.0+z);
    header.nbodies = npart;
    header.nsph = 0;
    header.ndark = npart;
    header.nstar = 0;
    header.ndim = 3;
    
    /*
      You have: 3/(8 pi G)
      You want: sunmass/(km/s/megaparsec)^2/megaparsec^3
        * 27749438
        / 3.6036765e-08
    */

    mconv = 3.6036765e-08/(h0*h0*boxsize*boxsize*boxsize);
    xconv = 1.0/boxsize;
    /* ratio of hubble velocity across the volume */
    vconv = sqrt(8*M_PI/3.0)/(h0*boxsize);
    
    xdrstdio_create(&xdrs, stdout, XDR_ENCODE);
    xdr_header(&xdrs, &header);
    for(i = 0; i < npart; i++) {
	fread((char *) &part, sizeof(part), 1, stdin);
	dark.mass = part.mass*mconv;
	dark.pos[0] = part.x*xconv - 0.5;
	dark.pos[1] = part.y*xconv - 0.5;
	dark.pos[2] = part.z*xconv - 0.5;
	dark.vel[0] = part.vx*vconv;
	dark.vel[1] = part.vy*vconv;
	dark.vel[2] = part.vz*vconv;
	dark.eps = 0.0;
	dark.phi = 0.0;
	xdr_dark(&xdrs, &dark);
	}
    xdr_destroy(&xdrs);
    return 0;
}
    
    
    
    
    
