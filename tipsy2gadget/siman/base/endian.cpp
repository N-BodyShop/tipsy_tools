//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include <stdio.h>

void swapDataEndian(void *data_vd,int unitsize,int datalen)
{

  char *data = (char*) data_vd;
  if(datalen % unitsize !=0) {
    printf("datalen = %d; unitsize = %d",datalen,unitsize);
    printf("swapEndian: length of data is invalid given unit size, returning\n");
    
    return;
  }

  if(unitsize % 2 !=0) {
    printf("swapEndian: cannot swap non-even number of fields, returning");
    return;
  }

  char temp = 0;
  int offset=0;
  int n;

  while(offset < datalen) {

    for(n=0;n<=unitsize/2-1;n++) {
      temp = data[offset+unitsize-n-1];
      data[offset+unitsize-n-1] = data[offset+n];
      data[offset+n] = temp;
    }
    

    offset+=unitsize;
  }

}
