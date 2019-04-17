/* Program to trim a "standard" tipsy file to the correct length.
 * This can be used to clean up a file that was overwritten, but not
 * truncated; the primary example is when changa or gasoline
 * terminates in the middle of writing an output, and is then rerun
 * with a slightly different star formation history, resulting in an
 * overwritten file with a smaller number of stars.
 *
 * Usage: trimstd filename
 */
#include <stdio.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <unistd.h>
#include <errno.h>


int
main(int argc, char **argv)
{
    FILE *pFD;
    XDR xdrs;
    off_t nFSize;
    off_t nDSize;
    
    if(argc != 2) {
        fprintf(stderr, "Usage: trimstd filename\n");
        exit(-1);
    }
    pFD = fopen(argv[1], "r+");
    if(pFD == NULL) {
        fprintf(stderr, "Can't open file: %s\n", argv[1]);
        exit(-1);
    }
    xdrstdio_create(&xdrs, pFD, XDR_DECODE);
    if(xdr_header(&xdrs, &header) != 1) {
        fprintf(stderr, "Can't read header from file: %s\n", argv[1]);
        exit(-1);
    }
    fprintf(stderr, "Header data: %g %d %d %d %d %d\n", header.time,
            header.nbodies, header.ndim, header.nsph, header.ndark,
            header.nstar);
    nDSize = 32 + 48L*header.nsph + 36L*header.ndark + 44L*header.nstar;
    fseek(pFD, 0L, SEEK_END);
    nFSize = ftell(pFD);
    if(nDSize == nFSize) {
        fprintf(stderr, "File is correct size, nothing to be done\n");
        xdr_destroy(&xdrs);
        fclose(pFD);
        return 0;
    }
    if(nDSize > nFSize) {
        fprintf(stderr, "File is truncated! Perhaps a quota overflow?\n");
        fprintf(stderr, "Desired size: %ld, actual size: %ld\n", nDSize, nFSize);
        xdr_destroy(&xdrs);
        fclose(pFD);
        return -1;
    }
    xdr_destroy(&xdrs);
    fclose(pFD);
    fprintf(stderr, "Truncating file to %ld bytes\n", nDSize);
    if(truncate(argv[1], nDSize) != 0) {
        fprintf(stderr, "Truncate failed, errno %d\n", errno);
        return -1;
    }
    return 0;
}
