#include <stdio.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

static XDR xdrs;

int xdr_header()
{
  int pad = 0;
  
  if(xdr_double(&xdrs, &header.time) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nbodies) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.ndim) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nsph) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.ndark) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nstar) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &pad) != TRUE)
    return 0;
  return 1;
}

void xdr_gas()
{
  if(sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) gas_particles,
	     header.nsph*(sizeof(*gas_particles)/sizeof(Real)),
	     sizeof(Real), xdr_float);
    }
}  

void xdr_dark()
{
  if(sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) dark_particles,
	     header.ndark*(sizeof(*dark_particles)/sizeof(Real)),
	     sizeof(Real), xdr_float);
    }
}  

void xdr_star()
{
  if(sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) star_particles,
	     header.nstar*(sizeof(*star_particles)/sizeof(Real)),
	     sizeof(Real), xdr_float);
    }
}  

int
main()
{
  xdrstdio_create(&xdrs, stdout, XDR_ENCODE);
  forever {
	if(fread((char *)&header,sizeof(header),1,stdin) != 1)
	  break;
	if(gas_particles != NULL) free(gas_particles);
	if(header.nsph != 0) {
	    gas_particles = (struct gas_particle *)
				malloc(header.nsph*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		printf("<sorry, no memory for gas particles, master>\n") ;
		return -1;
	    }
	}
	else
	  gas_particles = NULL;

	if(dark_particles != NULL) free(dark_particles);
	if(header.ndark != 0) {
	    dark_particles = (struct dark_particle *)
				malloc(header.ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
	else
	  dark_particles = NULL;

	if(star_particles != NULL) free(star_particles);
	if(header.nstar != 0) {
	    star_particles = (struct star_particle *)
				malloc(header.nstar*sizeof(*star_particles));
	    if(star_particles == NULL) {
		printf("<sorry, no memory for star particles, master>\n") ;
		return -1;
	    }
	}
	else
	  star_particles = NULL;

	fread((char *)gas_particles,sizeof(struct gas_particle),
	      header.nsph, stdin) ;
	fread((char *)dark_particles,sizeof(struct dark_particle),
	       header.ndark,stdin) ;
	fread((char *)star_particles,sizeof(struct star_particle),
	      header.nstar, stdin) ;
  
	xdr_header();
	xdr_gas();
	xdr_dark();
	xdr_star();
	
	fprintf(stderr, "read time %f\n",header.time) ;
      }
  xdr_destroy(&xdrs);
  return 0;
}
