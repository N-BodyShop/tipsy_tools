#include <stdio.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

static XDR xdrs;

int xdr_header()
{
  int pad;
  
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

int xdr_gas(struct gas_particle *gas)
{
  if(sizeof(Real) == sizeof(float))
    {
      if(xdr_float(&xdrs, &gas->mass) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->pos[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->pos[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->pos[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->vel[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->vel[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->vel[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->rho) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->temp) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->hsmooth) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->metals) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &gas->phi) != TRUE)
	return 0;
      return 1;
    }
}  

int xdr_dark(struct dark_particle *dark)
{
  if(sizeof(Real) == sizeof(float))
    {
      if(xdr_float(&xdrs, &dark->mass) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->pos[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->pos[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->pos[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->vel[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->vel[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->vel[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->eps) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &dark->phi) != TRUE)
	return 0;
      return 1;
    }
}  

int xdr_star(struct star_particle *star)
{
  if(sizeof(Real) == sizeof(float))
    {
      if(xdr_float(&xdrs, &star->mass) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->pos[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->pos[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->pos[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->vel[0]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->vel[1]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->vel[2]) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->metals) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->tform) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->eps) != TRUE)
	return 0;
      if(xdr_float(&xdrs, &star->phi) != TRUE)
	return 0;
      return 1;
    }
}  

int
main()
{
    struct gas_particle gas;
    struct dark_particle dark;
    struct star_particle star;
    int i;

  xdrstdio_create(&xdrs, stdin, XDR_DECODE);
  forever {
	if(xdr_header() != 1)
	  break;
	fprintf(stderr, "%g %d %d %d %d %d\n", header.time,
		header.nbodies, header.ndim, header.nsph,
		header.ndark, header.nstar);


	fwrite((char *)&header,sizeof(header),1,stdout) ;
	for(i = 0; i < header.nsph; i++) {
	    xdr_gas(&gas);
	    fwrite((char *)&gas,sizeof(struct gas_particle), 1, stdout) ;
	}
	for(i = 0; i < header.ndark; i++) {
	    xdr_dark(&dark);
	    fwrite((char *)&dark,sizeof(struct dark_particle), 1, stdout) ;
	}
	for(i = 0; i < header.nstar; i++) {
	    xdr_star(&star);
	    fwrite((char *)&star,sizeof(struct star_particle), 1, stdout) ;
	}
	
	fprintf(stderr, "read time %f\n",header.time) ;
      }
  xdr_destroy(&xdrs);
  return 0;
}
