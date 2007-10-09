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

    xdrstdio_create(&xdrs,stdin,XDR_DECODE);
    xdr_int(&xdrs, &magic);
    assert(magic == sizeof(SFEVENT));
    
    printf("# iOrder iOrdGas timeForm rFormX rFormY rFormZ vFormX vFormY vFromZ massForm rhoForm TempForm\n");
    
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
	printf("%d %d %g %g %g %g %g %g %g %g %g %g\n", SfEv.iOrdStar,
	       SfEv.iOrdGas, SfEv.timeForm, SfEv.rForm[0], SfEv.rForm[1],
	       SfEv.rForm[2], SfEv.vForm[0], SfEv.vForm[1], SfEv.vForm[2],
	       SfEv.massForm, SfEv.rhoForm, SfEv.TForm);
	}
    xdr_destroy(&xdrs);
    return 0;
    }
