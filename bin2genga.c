/*
 * Copied from bin2ascii.c
 * Convert TIPSY format into GENGA format, with
 * << x y z m vx vy vz r >>
 */
/*
 */

#include <stdio.h>
#include <malloc.h>
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
	if(header.ndark != 0) {
	    dark_particles = (struct dark_particle *)
				malloc(header.ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
    
	fread((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,stdin) ;

	lastdp = dark_particles + header.ndark ;

    for(dp=dark_particles; dp< lastdp; dp++) {
        fprintf(stdout,"%g %g %g %g %g %g %g %g\n", dp->pos[0], dp->pos[1], dp->pos[2],
                dp->mass, dp->vel[0], dp->vel[1], dp->vel[2], dp->eps*2);
    }

	if(header.ndark != 0) {
	  free(dark_particles);
	}
    }
    return(0);
}
