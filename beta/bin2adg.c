#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <math.h>
#include "tipsydefs.h"
#include "hughdefs.h"

#define DEBUG

void main(int argc,char **argv)
{
	int i,j,nDark;
	struct dump h;
	struct dark_particle *pdp;
	struct header1 h1;
	struct header2 h2;
	int iRec;
	float r[3],v[3];
	/*
	 ** Parameters
	 */
	float alpha = 1.5;
	int L = 32;
	float perr = 6.0;
	float rLbox = 80.0;
	float soft = 0.6;
	
	fread(&h,sizeof(struct dump),1,stdin);
	assert(h.nbodies == h.ndark);
	nDark = h.ndark;
	pdp = malloc(nDark*sizeof(struct dark_particle));
	fread(pdp,sizeof(struct dark_particle),nDark,stdin);

	h1.iRecStart = REC_HEADER1;
	h1.nobj = nDark;
	h1.ips = 0;
	h1.aexp = h.time;
	h1.p = 3.0/2.0/alpha*pow(h1.aexp,alpha);
	h1.dpnp1 = 0.0;
	h1.dpn = 0.0;
	h1.iRecEnd = REC_HEADER1;
	fwrite(&h1,sizeof(struct header1),1,stdout);

	h2.iRecStart = REC_HEADER2;
	h2.u.common.nobj = h1.nobj;
	h2.u.common.ips = h1.ips;
	h2.u.common.istart = 0;
	h2.u.common.irun = 0;
	h2.u.common.ipout = 0;
	h2.u.common.ipstop = 0;
	h2.u.common.ipdump = 0;
	h2.u.common.iseed1 = 0;
	h2.u.common.iseed2 = 0;
	h2.u.common.L = L;
	h2.u.common.iru = 0;
	h2.u.common.nlmx = 0;
	h2.u.common.alpha = alpha;
	h2.u.common.perr = perr;
	h2.u.common.aexp = h1.aexp;
	h2.u.common.p = h1.p;
	h2.u.common.dpnp1 = 0.0;
	h2.u.common.dpn = 0.0;
	h2.u.common.Ust = 0.0;
	h2.u.common.Tst = 0.0;
	h2.u.common.Uh = 0.0;
	h2.u.common.Usum = 0.0;
	h2.u.common.Usm = 0.0;
	h2.u.common.T = 0.0;
	h2.u.common.Tsum = 0.0;
	h2.u.common.delI = 0.0;
	h2.u.common.delC = 0.0;
	h2.u.common.dfCU = 0.0;
	h2.u.common.rLbox = rLbox;
	h2.u.common.delta = 0.0;
	h2.u.common.soft = soft;
	h2.u.common.sftinit = soft;
	h2.u.common.sftmin = 0.0;	
	h2.iRecEnd = REC_HEADER2;
	fwrite(&h2,sizeof(struct header2),1,stdout);
	/*
	 ** Now write out positions, making transformation to silly coords.
	 ** Assume tipsy binary has a period of 1 and box centered at 0.
	 */
	iRec = 3*nDark*sizeof(float);
	fwrite(&iRec,sizeof(int),1,stdout);
	for (i=0;i<nDark;++i) {
		for (j=0;j<3;++j) {
			r[j] = (pdp[i].pos[j]+0.5)*L + 1.0;
			}
		fwrite(r,sizeof(float),3,stdout);
		}
	fwrite(&iRec,sizeof(int),1,stdout);
	/*
	 ** Now write out the velocities.
	 */
	fwrite(&iRec,sizeof(int),1,stdout);
	for (i=0;i<nDark;++i) {
		for (j=0;j<3;++j) {
			v[j] = pdp[i].vel[j];
			}
		fwrite(v,sizeof(float),3,stdout);
		}
	fwrite(&iRec,sizeof(int),1,stdout);	
	}





