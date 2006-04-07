//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// TIPSYFILE.HPP

#include "simsnap.hpp"
#include "particle.hpp"

#ifndef __TIPSYFILE_H_INCLUDED

#define __TIPSYFILE_H_INCLUDED

// tipsy structures:

#define TIPSY_MAXDIM 3

typedef struct _tipsy_gas_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float rho;
    float temp;
    float hsmooth;
    float metals ;
    float phi ;
} tipsy_gas_particle;


typedef struct _tipsy_dark_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float eps;
    float phi ;
} tipsy_dark_particle;


typedef struct _tipsy_star_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float metals ;
    float tform ;
    float eps;
    float phi ;
} tipsy_star_particle;


typedef struct _tipsy_header {
  double time ;
  int nbodies ;
  int ndim ;
  int nsph ;
  int ndark ;
  int nstar ;
  int zero;
} tipsy_header;



class CTipsyFile : public CBaseSimSnap {
 public:
  CTipsyFile(const char *filename);

  ~CTipsyFile();

  bool load();

  /// write a file in native TIPSY format
  ///
  /// @param sim - simulation data
  /// @param filename - filename
  /// This procedure also writes "filename".units which contains the necessary
  /// internal -> physical conversion units
  static void nativeWrite(CSimSnap *sim, std::string filename);
  
 private:

  bool Loaded;

  int NumPart;
  int Ngas;
  
  tipsy_header header;

  char fname[1024];


  FILE *fileHandle;
  fpos_t startOfData;

};

#endif
