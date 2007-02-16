/* $Header$
 * $Log$
 * Revision 1.1  2007/02/16 02:29:53  stinson
 * Created new, possibly useful program to convert xdr formated (standard)
 * files to tipsy ascii format.
 *
 * Revision 1.1.1.1  1995/01/20 20:02:45  trq
 * Initial revision
 *
 */
/*
 */

#include <stdio.h>
#include <malloc.h>
#include "tipsydefs.h"

int
main()
{
    struct gas_particle *gas_particles, *gp, *lastgp;
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct star_particle *star_particles, *sp, *lastsp;
    struct dump header;
    int i;
    XDR xdrs;
    
    xdrstdio_create(&xdrs, stdin, XDR_DECODE);
    while(1) 
      {
        if(xdr_header(&xdrs, &header) != 1)
          break;
	if(header.nsph != 0) {
	    gas_particles = (struct gas_particle *)
				malloc(header.nsph*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		printf("<sorry, no memory for gas particles, master>\n") ;
		return -1;
	    }
	}
	if(header.ndark != 0) {
	    dark_particles = (struct dark_particle *)
				malloc(header.ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
	if(header.nstar != 0) {
	    star_particles = (struct star_particle *)
				malloc(header.nstar*sizeof(*star_particles));
	    if(star_particles == NULL) {
		printf("<sorry, no memory for star particles, master>\n") ;
		return -1;
	    }
        }
    
/*	fread((char *)gas_particles,sizeof(struct gas_particle),
			 header.nsph,stdin) ;
	fread((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,stdin) ;
	fread((char *)star_particles,sizeof(struct star_particle),
			 header.nstar,stdin) ;*/
        for(i = 0; i < header.nsph; i++) xdr_gas(&xdrs, &gas_particles[i]);
        for(i = 0; i < header.ndark; i++) xdr_dark(&xdrs, &dark_particles[i]);
        for(i = 0; i < header.nstar; i++) xdr_star(&xdrs, &star_particles[i]);

	lastgp = gas_particles + header.nsph ;
	lastdp = dark_particles + header.ndark ;
	lastsp = star_particles + header.nstar ;
	fprintf(stdout, "%d %d %d\n" ,header.nbodies, header.nsph,
	       header.nstar) ;
	fprintf(stdout,"%d\n",header.ndim) ;
	fprintf(stdout,"%g\n",header.time) ;
	for(gp=gas_particles; gp < lastgp ; gp++){
	    fprintf(stdout,"%g\n",gp->mass);
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fprintf(stdout,"%g\n",dp->mass);
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    fprintf(stdout,"%g\n",sp->mass);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->pos[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(stdout,"%g\n",dp->pos[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->pos[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->pos[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(stdout,"%g\n",dp->pos[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->pos[1]);
	}
	if (header.ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fprintf(stdout,"%g\n",gp->pos[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		fprintf(stdout,"%g\n",dp->pos[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fprintf(stdout,"%g\n",sp->pos[2]);
	    }
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->vel[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
		fprintf(stdout,"%g\n",dp->vel[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->vel[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->vel[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(stdout,"%g\n",dp->vel[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->vel[1]);
	}
	if (header.ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fprintf(stdout,"%g\n",gp->vel[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		fprintf(stdout,"%g\n",dp->vel[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fprintf(stdout,"%g\n",sp->vel[2]);
	    }
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(stdout,"%g\n",dp->eps);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->eps);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->rho);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->temp);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->hsmooth);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(stdout,"%g\n",gp->metals);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n",sp->metals);
	}
	for(sp=star_particles ; sp < lastsp ; sp++) {
	    fprintf(stdout,"%g\n", sp->tform);
	}
	for(gp=gas_particles; gp < lastgp ; gp++){
	    fprintf(stdout,"%g\n",gp->phi);
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fprintf(stdout,"%g\n",dp->phi);
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    fprintf(stdout,"%g\n",sp->phi);
	}
	fprintf(stderr,"read time %f\n",header.time) ;
	if(header.nsph != 0) {
	  free(gas_particles);
	}
	if(header.ndark != 0) {
	  free(dark_particles);
	}
	if(header.nstar != 0) {
	  free(star_particles);
	}
    }
    return(0);
}
