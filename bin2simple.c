/* $Header$
 * $Log$
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
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct dump header;
    
    while(1) 
      {
	if(fread((char *)&header,sizeof(header),1,stdin) == 0)
	  break;
	assert(header.nsph == 0);
	if(header.ndark != 0) {
	    dark_particles = (struct dark_particle *)
				malloc(header.ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
	assert(header.nstar == 0);
    
	fread((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,stdin) ;

	lastdp = dark_particles + header.ndark ;
	fprintf(stdout, "%d\n" ,header.nbodies) ;
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fprintf(stdout,"%g %g %g %g %g %g %g\n",dp->mass,
		    dp->pos[0], dp->pos[1], dp->pos[2],
		    dp->vel[0], dp->vel[1], dp->vel[2]);
	}
	fprintf(stderr,"read time %f\n",header.time) ;
	if(header.ndark != 0) {
	  free(dark_particles);
	}
    }
    return(0);
}
