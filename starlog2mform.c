#include <stdio.h>
#include <rpc/xdr.h>
#include <assert.h>

typedef struct sfEvent 		/* Holds statistics of the star
				   formation event */
{
    int iOrdStar;
    int iOrdGas;
    double timeForm;
    double rForm[3];
    double vForm[3];
    double massForm;
    double rhoForm;
    double TForm;
    } SFEVENT;

int main(int argc, char **argv)
{
    XDR xdrs;
    int magic;
    FILE *fpStarlog;
    FILE *fpIOrd;
    long nParts;
    long *iord;
    double *mform;
    long maxIOrd;
    int i, j;
    
    assert(argc == 3);
    fpStarlog = fopen(argv[1], "r");
    assert(fpStarlog != NULL);
    fpIOrd = fopen(argv[2], "r");
    assert(fpIOrd != NULL);

    xdrstdio_create(&xdrs,fpStarlog,XDR_DECODE);
    xdr_int(&xdrs, &magic);
    assert(magic == sizeof(SFEVENT));

    fscanf(fpIOrd, "%ld", &nParts);

    iord = malloc(nParts*sizeof(*iord));
    assert(iord != NULL);
    
    for(i = 0; i < nParts; i++) {
        int nread = fscanf(fpIOrd, "%ld", &iord[i]);
        assert(nread == 1);
        }
    maxIOrd = iord[nParts-1];
    mform = calloc(maxIOrd+1, sizeof(*mform));

    while(1) {
	SFEVENT SfEv;
	SFEVENT *pSfEv = &SfEv;

	if(xdr_int(&xdrs, &(pSfEv->iOrdStar)) == 0)
	    break;		/* End of file (hopefully) */
	xdr_int(&xdrs, &(pSfEv->iOrdGas));
	xdr_double(&xdrs, &(pSfEv->timeForm));
	xdr_double(&xdrs, &(pSfEv->rForm[0]));
	xdr_double(&xdrs, &(pSfEv->rForm[1]));
	xdr_double(&xdrs, &(pSfEv->rForm[2]));
	xdr_double(&xdrs, &(pSfEv->vForm[0]));
	xdr_double(&xdrs, &(pSfEv->vForm[1]));
	xdr_double(&xdrs, &(pSfEv->vForm[2]));
	xdr_double(&xdrs, &(pSfEv->massForm));
	xdr_double(&xdrs, &(pSfEv->rhoForm));
	xdr_double(&xdrs, &(pSfEv->TForm));
        if(SfEv.iOrdStar <= maxIOrd)
            mform[SfEv.iOrdStar] = SfEv.massForm;
	}
    xdr_destroy(&xdrs);
    printf("%ld\n", nParts);
    j = 0;
    for(i = 0; i < nParts; i++) {
        while(j < iord[i])
            j++;
        printf("%g\n", mform[j]);
        }
    return 0;
}

