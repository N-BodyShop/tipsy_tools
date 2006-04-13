/* $Header$
 * $Log$
 * Revision 1.1  2006/04/13 21:13:32  trq
 * Read simple files.
 *
 * Revision 1.1  2001/05/17 19:33:34  trq
 * Print out tipsy files in a "simple" format.
 *
 * Revision 1.1.1.1  1995/01/20  20:02:45  trq
 * Initial revision
 *
 */
/*
 */

#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include "tipsydefs.h"

int
main()
{
    struct gas_particle *gas_particles, *gp, *lastgp;
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct star_particle *star_particles, *sp, *lastsp;
    struct dump header;
    
    scanf("%d %d %d\n" ,&header.nbodies, &header.nsph, &header.nstar) ;
	header.ndark = header.nbodies - (header.nsph + header.nstar);
	header.ndim = 3;
	header.time = 0.0;

	fwrite((char *)&header,sizeof(header),1,stdout);
	
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
    
	lastgp = gas_particles + header.nsph ;
	for(gp=gas_particles; gp < lastgp;  gp++){
	    scanf("%g %g %g %g %g %g %g %g %g\n",&gp->mass,
		    &gp->pos[0], &gp->pos[1], &gp->pos[2],
		    &gp->vel[0], &gp->vel[1], &gp->vel[2], &gp->rho, &gp->temp);
	}

	fwrite((char *)gas_particles,sizeof(struct gas_particle),
			 header.nsph,stdout) ;

	lastdp = dark_particles + header.ndark ;
	for(dp=dark_particles; dp < lastdp;  dp++){
	    scanf("%g %g %g %g %g %g %g\n",&dp->mass,
		    &dp->pos[0], &dp->pos[1], &dp->pos[2],
		    &dp->vel[0], &dp->vel[1], &dp->vel[2]);
	    dp->eps = 0.0;
	}

	fwrite((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,stdout) ;

	lastsp = star_particles + header.nstar ;
	for(sp=star_particles; sp < lastsp;  sp++){
	    scanf("%g %g %g %g %g %g %g\n",&sp->mass,
		    &sp->pos[0], &sp->pos[1], &sp->pos[2],
		    &sp->vel[0], &sp->vel[1], &sp->vel[2]);
	}
	fwrite((char *)star_particles,sizeof(struct star_particle),
			 header.nstar,stdin) ;

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
    return(0);
}
