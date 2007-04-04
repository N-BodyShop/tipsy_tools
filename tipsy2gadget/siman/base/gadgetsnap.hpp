// gadgetsnap.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __GADGETFILE_H_INCLUDED

#define __GADGETFILE_H_INCLUDED

namespace siman {

// STRUCTURE DEFINITIONS

typedef struct {
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_stellarage;
  int      flag_metals;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 2*4];  /* fills to 256 Bytes */
} gadget_header;

typedef struct {
  
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne, Metal;
} gadget_particle;

class GadgetSnap : public BaseSimSnap {
 public:
  GadgetSnap(const char *filename, bool swapendian = false);
  GadgetSnap(const char *path, const char *snapname, int snapshot, bool swapendian=false);
  virtual ~GadgetSnap();

  int Reorder();

  
  static void nativeWrite(const SimSnap *sim, std::string filename);

 private:
  
  bool load();
  static void writeField(std::ofstream *file, char *buf, int size);

  // Keep all data private - firstly to prevent accidental
  // overwriting, and secondly so that the class may be
  // extended later to access different formats or
  // read direct from disk for large files

  bool Loaded;

  gadget_header header;
  
  int *Id;
  char fname[1024];
  bool swapEndian;

  FILE *fileHandle;
  fpos_t startOfData;

    
 private:
  int allocate_memory();
  void fileDebug(FILE *fd);
};

}

#endif

