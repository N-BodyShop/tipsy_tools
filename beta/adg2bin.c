#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsydefs.h"
#include "hughdefs.h"

#define DEBUG

void main(void)
{
	struct header1 h1;
	struct header2 h2;
	int iRecStart,iRecEnd;
	int i,j,nDark;
	struct dark_particle *pdp;
	struct dump h;

	fread(&h1,sizeof(struct header1),1,stdin);
	fread(&h2,sizeof(struct header2),1,stdin);
#ifdef DEBUG
	fprintf(stderr,"RECORD 1:%d %d\n",h1.iRecStart,h1.iRecEnd);
	fprintf(stderr,"nobj:\t\t\t%d\n",h1.nobj);
	fprintf(stderr,"ips:\t\t\t%d\n",h1.ips);
	fprintf(stderr,"aexp:\t\t\t%f\n",h1.aexp);
	fprintf(stderr,"p:\t\t\t%f\n",h1.p);
	fprintf(stderr,"dpnp1:\t\t\t%f\n",h1.dpnp1);
	fprintf(stderr,"dpn:\t\t\t%f\n",h1.dpn);

	fprintf(stderr,"RECORD 2:%d %d\n",h2.iRecStart,h2.iRecEnd);
	fprintf(stderr,"nobj:\t\t\t%d\n",h2.u.common.nobj);
	fprintf(stderr,"ips:\t\t\t%d\n",h2.u.common.ips);
	fprintf(stderr,"istart:\t\t\t%d\n",h2.u.common.istart);
	fprintf(stderr,"irun:\t\t\t%d\n",h2.u.common.irun);
	fprintf(stderr,"ipout:\t\t\t%d\n",h2.u.common.ipout);
	fprintf(stderr,"ipstop:\t\t\t%d\n",h2.u.common.ipstop);
	fprintf(stderr,"ipdump:\t\t\t%d\n",h2.u.common.ipdump);
	fprintf(stderr,"iseed1:\t\t\t%d\n",h2.u.common.iseed1);
	fprintf(stderr,"iseed2:\t\t\t%d\n",h2.u.common.iseed2);
	fprintf(stderr,"L:\t\t\t%d\n",h2.u.common.L);
	fprintf(stderr,"iru:\t\t\t%d\n",h2.u.common.iru);
	fprintf(stderr,"nlmx:\t\t\t%d\n",h2.u.common.nlmx);
	fprintf(stderr,"alpha:\t\t\t%f\n",h2.u.common.alpha);
	fprintf(stderr,"perr:\t\t\t%f\n",h2.u.common.perr);
	fprintf(stderr,"aexp:\t\t\t%f\n",h2.u.common.aexp);
	fprintf(stderr,"p:\t\t\t%f\n",h2.u.common.p);
	fprintf(stderr,"dpnp1:\t\t\t%f\n",h2.u.common.dpnp1);
	fprintf(stderr,"dpn:\t\t\t%f\n",h2.u.common.dpn);
	fprintf(stderr,"Ust:\t\t\t%f\n",h2.u.common.Ust);
	fprintf(stderr,"Tst:\t\t\t%f\n",h2.u.common.Tst);
	fprintf(stderr,"Uh:\t\t\t%f\n",h2.u.common.Uh);
	fprintf(stderr,"Usum:\t\t\t%f\n",h2.u.common.Usum);
	fprintf(stderr,"Usm:\t\t\t%f\n",h2.u.common.Usm);
	fprintf(stderr,"T:\t\t\t%f\n",h2.u.common.T);
	fprintf(stderr,"Tsum:\t\t\t%f\n",h2.u.common.Tsum);
	fprintf(stderr,"delI:\t\t\t%f\n",h2.u.common.delI);
	fprintf(stderr,"delC:\t\t\t%f\n",h2.u.common.delC);
	fprintf(stderr,"dfCU:\t\t\t%f\n",h2.u.common.dfCU);
	fprintf(stderr,"rLbox:\t\t\t%f\n",h2.u.common.rLbox);
	fprintf(stderr,"delta:\t\t\t%f\n",h2.u.common.delta);
	fprintf(stderr,"soft:\t\t\t%f\n",h2.u.common.soft);
	fprintf(stderr,"sftinit:\t\t%f\n",h2.u.common.sftinit);
	fprintf(stderr,"sftmin:\t\t\t%f\n",h2.u.common.sftmin);
#endif
	/*
	 ** Allocate tipsy particles.
	 */
	nDark = h1.nobj;
	pdp = malloc(nDark*sizeof(struct dark_particle));
	assert(pdp != NULL);
	/*
	 ** Start reading positions into tipsy dark particles.
	 */
	fread(&iRecStart,sizeof(int),1,stdin);
	for (i=0;i<nDark;++i) {
		fread(pdp[i].pos,sizeof(float),3,stdin);
		}
	fread(&iRecEnd,sizeof(int),1,stdin);
#ifdef DEBUG
	fprintf(stderr,"RECORD 3:%d %d POSITIONS\n",iRecStart,iRecEnd);
#endif
	/*
	 ** Start reading velocities into tipsy dark particles.
	 */
	fread(&iRecStart,sizeof(int),1,stdin);
	for (i=0;i<nDark;++i) {
		fread(pdp[i].vel,sizeof(float),3,stdin);
		}
	fread(&iRecEnd,sizeof(int),1,stdin);
#ifdef DEBUG
	fprintf(stderr,"RECORD 3:%d %d VELOCITIES\n",iRecStart,iRecEnd);
#endif
	/*
	 ** Finally set up the remaining fields in the tipsy dark particle.
	 */
	for (i=0;i<nDark;++i) {
		pdp[i].mass = 1.0;
		pdp[i].eps = 0.0;
		pdp[i].phi = 0.0;
		/*
		 ** Shift the box so that it has a period of 1 and is centered at 0.
		 */
		for (j=0;j<3;++j) {
			pdp[i].pos[j] = (pdp[i].pos[j]-1.0)/h2.u.common.L - 0.5;
			}
		}
	/*
	 ** Set up the tipsy header and output.
	 */
	h.time = h1.aexp;
	h.nbodies = nDark;
	h.ndim = 3;
	h.nsph = 0;
	h.ndark = nDark;
	h.nstar = 0;
	fwrite(&h,sizeof(struct dump),1,stdout);
	fwrite(pdp,sizeof(struct dark_particle),nDark,stdout);
	}



