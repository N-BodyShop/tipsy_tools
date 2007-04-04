// tipsysnap.hpp - part of SimAn Simulation Analysis Library
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







#include "simsnap.hpp"
#include "particle.hpp"

#ifndef __TIPSYFILE_H_INCLUDED

#define __TIPSYFILE_H_INCLUDED

namespace siman {

  namespace tipsy {
    // tipsy structures:
    
#define TIPSY_MAXDIM 3
    
    typedef struct {
      float mass;
      float pos[TIPSY_MAXDIM];
      float vel[TIPSY_MAXDIM];
      float rho;
      float temp;
      float hsmooth;
      float metals ;
      float phi ;
    } gas_particle;
    
    
    typedef struct {
      float mass;
      float pos[TIPSY_MAXDIM];
      float vel[TIPSY_MAXDIM];
      float eps;
      float phi ;
    } dark_particle;
    
    
    typedef struct {
      float mass;
      float pos[TIPSY_MAXDIM];
      float vel[TIPSY_MAXDIM];
      float metals ;
      float tform ;
      float eps;
      float phi ;
    } star_particle;
    
    
    typedef struct {
      double time ;
      int nbodies ;
      int ndim ;
      int nsph ;
      int ndark ;
      int nstar ;
      int zero;
    } header;
    
  } // namespace tipsy

class TipsySnap : public BaseSimSnap {
 public:
  TipsySnap(const char *filename);

  virtual ~TipsySnap();

  bool load();

  /// write a file in native TIPSY format
  ///
  /// @param sim - simulation data
  /// @param filename - filename
  /// This procedure also writes "filename".units which contains the necessary
  /// internal -> physical conversion units
  static void nativeWrite(const SimSnap *sim, std::string filename);
  
 private:

  bool Loaded;

  int NumPart;
  int Ngas;
  
  tipsy::header header;

  char fname[1024];


  FILE *fileHandle;
  fpos_t startOfData;

};

} // namespace siman

#endif
