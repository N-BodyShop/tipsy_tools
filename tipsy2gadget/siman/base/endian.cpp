// endian.cpp - part of SimAn Simulation Analysis Library
//
//
// Copyright (c) Andrew Pontzen 2005, 2006
//
// SimAn is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// SimAn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public Licence for more details.
//
// You should have received a copy of the GNU General Public Licence
// along with SimAn; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA









#include <stdio.h>

namespace siman {

namespace endian {
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



  char end_tp; 
  char* end4_ptr; 
  char* end8_ptr;


}

} // namespace siman
