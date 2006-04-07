//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __ENDIAN_HPP_INCLUDED

#define __ENDIAN_HPP_INCLUDED

// Optimized swap endian routines.
// Note that for loops and if branches really slow things down, so these
// routines are much preferable to using a more generalized procedure.

#define ENDINIT register char end_tp; char* end4_ptr; char* end8_ptr

#define END4(ENDPTR) end4_ptr = (char*) ENDPTR; end_tp = end4_ptr[0]; end4_ptr[0] = end4_ptr[3]; end4_ptr[3] = end_tp; end_tp = end4_ptr[1]; end4_ptr[1]=end4_ptr[2]; end4_ptr[2]=end_tp

#define END8(ENDPTR)  end8_ptr = (char*) ENDPTR; end_tp = end8_ptr[0]; end8_ptr[0] = end8_ptr[7]; end8_ptr[7] = end_tp; end_tp = end8_ptr[1]; end8_ptr[1] = end8_ptr[6]; end8_ptr[6] = end_tp; end_tp = end8_ptr[2]; end8_ptr[2] = end8_ptr[5]; end8_ptr[5] = end_tp; end_tp = end8_ptr[3]; end8_ptr[3] = end8_ptr[4]; end8_ptr[4] = end_tp

// Generic swap endian routine

extern void swapDataEndian(void *data_vd, int unitsize, int datalen);

#endif
