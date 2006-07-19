#include <stdio.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>


int
main()
{
    struct gas_particle gas;
    struct dark_particle dark;
    struct star_particle star;
    int i;
    XDR xdrs;
    struct dump header;

  xdrstdio_create(&xdrs, stdin, XDR_DECODE);
  forever {
	if(xdr_header(&xdrs, &header) != 1)
	  break;
	fprintf(stderr, "%g %d %d %d %d %d\n", header.time,
		header.nbodies, header.ndim, header.nsph,
		header.ndark, header.nstar);


	fwrite((char *)&header,sizeof(header),1,stdout) ;
	for(i = 0; i < header.nsph; i++) {
	    xdr_gas(&xdrs, &gas);
	    fwrite((char *)&gas,sizeof(struct gas_particle), 1, stdout) ;
	}
	for(i = 0; i < header.ndark; i++) {
	    xdr_dark(&xdrs, &dark);
	    fwrite((char *)&dark,sizeof(struct dark_particle), 1, stdout) ;
	}
	for(i = 0; i < header.nstar; i++) {
	    xdr_star(&xdrs, &star);
	    fwrite((char *)&star,sizeof(struct star_particle), 1, stdout) ;
	}
	
	fprintf(stderr, "read time %f\n",header.time) ;
      }
  xdr_destroy(&xdrs);
  return 0;
}
