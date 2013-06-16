#include <stdio.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

/* Sanitiy check to be sure all data is written */

int
main()
{
    struct gas_particle gas;
    struct dark_particle dark;
    struct star_particle star;
    int i;
    XDR xdrs;
    struct dump header;
    double mass_min, mass_max;
    double soft_min, soft_max;
    double mass_gas_t, mass_dark_t, mass_star_t;
	
    
    soft_min = mass_min = 1e38;
    soft_max = mass_max = -1e38;
    mass_gas_t = mass_dark_t = mass_star_t = 0.0;
    xdrstdio_create(&xdrs, stdin, XDR_DECODE);
	if(xdr_header(&xdrs, &header) != 1)
	  return -1;
	fprintf(stderr, "%g %d %d %d %d %d\n", header.time,
		header.nbodies, header.ndim, header.nsph,
		header.ndark, header.nstar);


	for(i = 0; i < header.nsph; i++) {
	    xdr_gas(&xdrs, &gas);
	    mass_gas_t += gas.mass;
	    if(gas.mass > mass_max)
		mass_max = gas.mass;
	    if(gas.hsmooth > soft_max)
		soft_max = gas.hsmooth;
	    if(gas.mass < mass_min)
		mass_min = gas.mass;
	    if(gas.hsmooth < soft_min)
		soft_min = gas.hsmooth;
	}
	fprintf(stderr, "Gas: mmin %g mmax %g smin %g smax %g\n", mass_min, mass_max,
		soft_min, soft_max);
    soft_min = mass_min = 1e38;
    soft_max = mass_max = -1e38;
	for(i = 0; i < header.ndark; i++) {
	    xdr_dark(&xdrs, &dark);
	    mass_dark_t += dark.mass;
	    if(dark.mass > mass_max)
		mass_max = dark.mass;
	    if(dark.eps > soft_max)
		soft_max = dark.eps;
	    if(dark.mass < mass_min)
		mass_min = dark.mass;
	    if(dark.eps < soft_min)
		soft_min = dark.eps;
	}
	fprintf(stderr, "Dark: mmin %g mmax %g smin %g smax %g\n", mass_min, mass_max,
		soft_min, soft_max);
    soft_min = mass_min = 1e38;
    soft_max = mass_max = -1e38;
	for(i = 0; i < header.nstar; i++) {
	    xdr_star(&xdrs, &star);
	    mass_star_t += star.mass;
	    if(star.mass > mass_max)
		mass_max = star.mass;
	    if(star.eps > soft_max)
		soft_max = star.eps;
	    if(star.mass < mass_min)
		mass_min = star.mass;
	    if(star.eps < soft_min)
		soft_min = star.eps;
	}
	fprintf(stderr, "Star: mmin %g mmax %g smin %g smax %g\n", mass_min, mass_max,
		soft_min, soft_max);
	fprintf(stderr, "read time %f\n",header.time) ;
	fprintf(stderr, "total gas, dark, star, all: %g %g %g %g\n",
		mass_gas_t, mass_dark_t, mass_star_t,
		mass_gas_t + mass_dark_t + mass_star_t) ;
  xdr_destroy(&xdrs);
  return 0;
}
