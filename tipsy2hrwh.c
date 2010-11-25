/*
 * Convert format of Heitmann, Ricker, Warren & Habib (2005 ApJS 160)
 * to tipsy format.
 */
#include <stdio.h>
#include <math.h>
#include <assert.h>
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
    double h0;
    double boxsize;
    double mconv;
    double xconv;
    double vconv;
    int i;

    if(argc != 3) {
	fprintf(stderr, "hrwh2std boxsize(Mpc) H0(km/s/Mpc)\n");
	return 1;
	}

    boxsize = atof(argv[1]);
    h0 = atof(argv[2]);
    
    xdrstdio_create(&xdrs, stdin, XDR_DECODE);
    if(xdr_header(&xdrs, &header) != 1)
	assert(0);
    fprintf(stderr, "%g %d %d %d %d %d\n", header.time,
	    header.nbodies, header.ndim, header.nsph,
	    header.ndark, header.nstar);
    assert(header.ndim == 3);
    assert(header.nsph == 0);
    assert(header.nstar == 0);

    /*
      You have: 3/(8 pi G)
      You want: sunmass/(km/s/megaparsec)^2/megaparsec^3
        * 27749438
        / 3.6036765e-08
    */

    mconv = 1.0/(3.6036765e-08/(h0*h0*boxsize*boxsize*boxsize));
    xconv = boxsize;
    /* ratio of hubble velocity across the volume */
    vconv = (h0*boxsize)/sqrt(8*M_PI/3.0);
    
    for(i = 0; i < header.ndark; i++) {
	xdr_dark(&xdrs, &dark);
	part.mass = dark.mass*mconv;
	part.x = (dark.pos[0] + 0.5)*xconv;
	part.y = (dark.pos[1] + 0.5)*xconv;
	part.z = (dark.pos[2] + 0.5)*xconv;
	part.vx = dark.vel[0]*vconv;
	part.vy = dark.vel[1]*vconv;
	part.vz = dark.vel[2]*vconv;
	fwrite((char *) &part, sizeof(part), 1, stdout);
	}
    xdr_destroy(&xdrs);
    return 0;
    }
