#include <stdio.h>
#include "tipsydefs.h"

int free();
void *malloc();

int
main(argc,argv)
char *argv[];
int argc;
{
    struct file {
	char name[50];
	FILE *ptr ;
	int pipe;
    } ;
    struct file binaryfile ;
    double time ;
    static double currtime = 0.0;
    static long currpos = 0L ;
    static long lastpos = 0L ;


    sscanf(argv[1],"%s",binaryfile.name) ;
    sscanf(argv[2],"%lf",&time) ;
    if( (binaryfile.ptr = fopen(binaryfile.name,"r"))== NULL) {
	fprintf(stderr,"<sorry, binary file %s does not exist, %s>\n",
	    binaryfile.name) ;
    }
    else{
	if ((float)currtime > (float)time){
	    fseek(binaryfile.ptr,0L,0);
	    currtime=0.0;
	    currpos=0;
	}
	forever {
	    if(fread((char *)&header,sizeof(header),1,binaryfile.ptr) 
	       != 1) {
		fprintf(stderr,"<sorry time too large, using %f>\n",
				 (float)currtime) ;
		break ;
	    }
	    currtime = header.time ;
	    currpos = ftell(binaryfile.ptr) - sizeof(header);
	    if ( (float)header.time >= (float)time ) 
		break ;
	    fseek(binaryfile.ptr,
		sizeof(gas_particles[0])*header.nsph +
		sizeof(dark_particles[0])*header.ndark +
		sizeof(star_particles[0])*header.nstar,
		1) ;
	}	
	fseek(binaryfile.ptr,currpos,0) ;
	lastpos = currpos ;
	fread((char *)&header,sizeof(header),1,binaryfile.ptr) ;

	if(gas_particles != NULL) free(gas_particles);
	if(header.nsph != 0) {
	    gas_particles = (struct gas_particle *)
				malloc(header.nsph*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		fprintf(stderr,
		    "<sorry, no memory for gas particles>\n") ;
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
		fprintf(stderr,
		    "<sorry, no memory for dark particles>\n") ;
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
		fprintf(stderr,
		    "<sorry, no memory for star particles>\n") ;
		return -1;
	    }
	}
	else
	  star_particles = NULL;

	fread((char *)gas_particles,sizeof(struct gas_particle),
			 header.nsph,binaryfile.ptr) ;
	fread((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,binaryfile.ptr) ;
	fread((char *)star_particles,sizeof(struct star_particle),
			 header.nstar,binaryfile.ptr) ;
	currpos = lastpos ;
	fseek(binaryfile.ptr,currpos,0) ;
	currtime = header.time ;
	if ((float)time != (float)currtime){
	    fprintf(stderr,"<used time %f, hope you don't mind>\n",
		   (float)currtime);
	}
	fwrite((char *)&header,sizeof(header),1,stdout) ;
	fwrite((char *)gas_particles,sizeof(struct gas_particle),header.nsph,stdout) ;
	fwrite((char *)dark_particles,sizeof(struct dark_particle),
	       header.ndark,stdout) ;
	fwrite((char *)star_particles,sizeof(struct star_particle),
	   header.nstar,stdout) ;
    }
    return 0;
}
