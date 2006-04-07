//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// GADGETFILE.HPP

#ifndef __GADGETFILE_H_INCLUDED

#define __GADGETFILE_H_INCLUDED

#include <string>

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

class CGadgetFile : public CBaseSimSnap {
 public:
  CGadgetFile(const char *filename, bool swapendian = false);
  CGadgetFile(const char *path, const char *snapname, int snapshot, bool swapendian=false);
  ~CGadgetFile();

  int Reorder();

  
  static void nativeWrite(CSimSnap *sim, std::string filename);

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

#endif

