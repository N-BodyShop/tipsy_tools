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
    struct dump header;
    int i;
    XDR xdrs;

  setvbuf(stdin, NULL, _IOFBF, 32*4096);
  setvbuf(stdout, NULL, _IOFBF, 32*4096);
  xdrstdio_create(&xdrs, stdout, XDR_ENCODE);
  forever {
	if(fread((char *)&header,sizeof(header),1,stdin) != 1)
	  break;
	fprintf(stderr, "reading time %f\n",header.time) ;
	xdr_header(&xdrs, &header);

	for(i = 0; i < header.nsph; i++) {
	    fread((char *)&gas,sizeof(struct gas_particle), 1, stdin) ;
	    xdr_gas(&xdrs, &gas);
	}
	for(i = 0; i < header.ndark; i++) {
	    fread((char *)&dark,sizeof(struct dark_particle), 1, stdin) ;
	    xdr_dark(&xdrs, &dark);
	}
	for(i = 0; i < header.nstar; i++) {
	    fread((char *)&star,sizeof(struct star_particle), 1, stdin) ;
	    xdr_star(&xdrs, &star);
	}
	
	fprintf(stderr, "read time %f\n",header.time) ;
      }
  xdr_destroy(&xdrs);
  return 0;
}
