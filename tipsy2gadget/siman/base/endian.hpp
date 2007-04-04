// endian.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __ENDIAN_HPP_INCLUDED

#define __ENDIAN_HPP_INCLUDED

namespace siman {
  namespace endian {

  extern char end_tp; 
  extern char* end4_ptr; 
  extern char* end8_ptr;
  
  // Generic swap endian routine
  
  extern void swapDataEndian(void *data_vd, int unitsize, int datalen);
  
  template <typename T> 
  inline void flip4(T* ENDPTR) {
    endian::end4_ptr = (char*) ENDPTR; 
    end_tp = endian::end4_ptr[0]; 
    endian::end4_ptr[0] = endian::end4_ptr[3]; 
    endian::end4_ptr[3] = end_tp; 
    end_tp = endian::end4_ptr[1]; 
    endian::end4_ptr[1]=endian::end4_ptr[2]; 
    endian::end4_ptr[2]=end_tp;

  }

  template <typename T>
  inline void flip8(T* ENDPTR) {
    endian::end8_ptr = (char*) ENDPTR; 
    end_tp = endian::end8_ptr[0]; 
    endian::end8_ptr[0] = endian::end8_ptr[7]; 
    endian::end8_ptr[7] = end_tp; 
    end_tp = endian::end8_ptr[1]; 
    endian::end8_ptr[1] = endian::end8_ptr[6]; 
    endian::end8_ptr[6] = end_tp;
    end_tp = endian::end8_ptr[2];
    endian::end8_ptr[2] = endian::end8_ptr[5]; 
    endian::end8_ptr[5] = end_tp; 
    end_tp = endian::end8_ptr[3]; 
    endian::end8_ptr[3] = endian::end8_ptr[4]; 
    endian::end8_ptr[4] = end_tp;
  }
}
}

#endif
