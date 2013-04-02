/* $Header$
 * $Log$
 * Revision 1.4  2013/04/02 22:45:29  trq
 * Initialize gas_particles, dark_particles, and star_particles.
 *
 * Revision 1.3  2006/07/19 00:10:06  trq
 * Fixed up conversion programs to match changes to tipsydefs.h
 *
 * Revision 1.2  1996/07/26 20:09:35  trq
 * Switched to accept output from genesis with the new star forming
 * algorithm.
 *
 * Revision 1.1.1.1  1995/01/20  20:02:46  trq
 * Initial revision
 *
 * Revision 1.1  93/09/14  09:51:22  trq
 * Initial revision
 * 
 * Revision 2.2  91/06/02  12:21:32  trq
 * Added star stuff from Neal
 * 
 * Revision 2.1  1991/04/18  18:11:35  trq
 * Added checks for malloc(0).
 *
 */
#include <stdio.h>
#include <malloc.h>
#include "tipsydefs.h"

int
main()
{
    int ndim ;
    int nbodies ;
    int ngas ;
    int ndark ;
    int nstar ;
    int count ;
    struct gas_particle *gp, *lastgp, *gas_particles = NULL;
    struct dark_particle *dp, *lastdp, *dark_particles = NULL;
    struct star_particle *sp, *lastsp, *star_particles = NULL;

    forever {
	count=fscanf(stdin, "%d%*[, \t\n]%d%*[, \t\n]%d"
		    ,&header.nbodies, &header.nsph, &header.nstar) ;
	if ( (count == EOF) || (count==0) ){
	    break ;
	}
	fscanf(stdin,"%d",&header.ndim) ;
	fscanf(stdin,"%lf",&header.time) ;
	ndim=header.ndim;
	nbodies=header.nbodies;
	ngas=header.nsph;
	nstar = header.nstar ;
	ndark = header.ndark = nbodies - nstar - ngas ;
	if(gas_particles != NULL) free(gas_particles);
	if(ngas != 0) {
	    gas_particles =
		(struct gas_particle *) malloc(ngas*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		fprintf(stderr,
			"<sorry, no memory for gas particles, master>\n") ;
		return -1;
	    }
	}
	else
	  gas_particles = NULL;
	if(dark_particles != NULL) free(dark_particles);
	if(ndark != 0) {
	    dark_particles =
		(struct dark_particle *) malloc(ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		fprintf(stderr,
			"<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
	else
	  dark_particles = NULL;
	if(star_particles != NULL) free(star_particles);
	if(nstar != 0) {
	    star_particles =
		 (struct star_particle *)malloc(nstar*sizeof(*star_particles));
	    if(star_particles == NULL) {
		fprintf(stderr,
		       "<sorry, no memory for star particles, master>\n") ;
		return -1;
	    }
	}
	else
	  star_particles = NULL;

	lastgp = gas_particles + ngas ;
	lastdp = dark_particles + ndark ;
	lastsp = star_particles + nstar ;

	for(gp=gas_particles; gp < lastgp ; gp++){
	    fscanf(stdin,"%f%*[, \t\n]",&gp->mass);
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fscanf(stdin,"%f%*[, \t\n]",&dp->mass);
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    fscanf(stdin,"%f%*[, \t\n]",&sp->mass);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->pos[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&dp->pos[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->pos[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->pos[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&dp->pos[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->pos[1]);
	}
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fscanf(stdin,"%f%*[, \t\n]",&gp->pos[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		fscanf(stdin,"%f%*[, \t\n]",&dp->pos[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fscanf(stdin,"%f%*[, \t\n]",&sp->pos[2]);
	    }
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->vel[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
		fscanf(stdin,"%f%*[, \t\n]",&dp->vel[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->vel[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->vel[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&dp->vel[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->vel[1]);
	}
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fscanf(stdin,"%f%*[, \t\n]",&gp->vel[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		count = fscanf(stdin,"%f%*[, \t\n]",&dp->vel[2]);
		if ( (count == EOF) || (count==0) ) 
		  {
		    fprintf(stderr, "error: short read on vz\n");
		    return -1;
		  }
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fscanf(stdin,"%f%*[, \t\n]",&sp->vel[2]);
	    }
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    count = fscanf(stdin,"%f%*[, \t\n]",&dp->eps);
	    if ( (count == EOF) || (count==0) ) 
	      {
		fprintf(stderr, "error: short read on epsilons\n");
		return -1;
	      }
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->eps);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->rho);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->temp);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->hsmooth);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&gp->metals);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->metals);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(stdin,"%f%*[, \t\n]",&sp->tform);
	}
	for(gp=gas_particles; gp < lastgp ; gp++){
	    count = fscanf(stdin,"%f%*[, \t\n]",&gp->phi);
	    if ( (count == EOF) || (count==0) )
	      return -1;
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    count = fscanf(stdin,"%f%*[, \t\n]",&dp->phi);
	    if ( (count == EOF) || (count==0) ) 
	      {
		fprintf(stderr, "error: short read on potentials\n");
		return -1;
	      }
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    count = fscanf(stdin,"%f%*[, \t\n]",&sp->phi);
	    if ( (count == EOF) || (count==0) )
	      return -1;
	}
	fwrite((char *)&header,sizeof(header),1,stdout) ;
	fwrite((char *)gas_particles,sizeof(struct gas_particle),ngas,stdout) ;
	fwrite((char *)dark_particles,sizeof(struct dark_particle),
	       ndark,stdout) ;
	fwrite((char *)star_particles,sizeof(struct star_particle),
	   nstar,stdout) ;
	fprintf(stderr, "read time %f\n",header.time) ;
    }
    return 0;
}
