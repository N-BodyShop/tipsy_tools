//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// Grid, ColumnGrid and supporting classes
//
// Basically for line-of-sight column density grids etc.

#ifndef __GRID_H_INCLUDED

#define __GRID_H_INCLUDED

class CGrid : public CSimSnap {
public:

  CGrid(CSimSnap *parentSim, float x1, float x2, int nxi, float y1, float y2, int nyi, float z1, float z2, int nzi, int maxRefine =0, int minPart = 0, bool displayStats=true);
  CGrid(CSimSnap *parentSim, int nxi, int nyi, int nzi, int maxRefine=0, int minPart=0, bool displayStats=true); // auto-ranging version 

  ~CGrid();

  void realize();

  virtual CTempColumnGrid operator[](int index);
  virtual CSimSnap* operator[](const coordinate<int> &coords);

  float getDx();
  float getDy();
  float getDz();
  float getX1();
  float getY1();
  float getZ1();
  int getNx();
  int getNy();
  int getNz();
  
  float getMeanRefineDepth();
  float getMeanLeafParticles();

  CSimSnap *getRegion(float x1, float x2, float y1, float y2, float z1, float z2);
  
  list<pair<int,double> > getNearestNeighbours(const CParticle &from,const int number, const CMetric &metric=CMetric());

  string className();

private:

  // helper functions for quick getNearestNeighbours code

  void doAdd(queue<coordinate<int> > &processQueue,bool* traversed,int cx, int cy, int cz);
  void addIfNecessary(queue<coordinate<int> > &processQueue,double furthest,const CParticle &from, bool *traversed,int cx, int cy, int cz);
  
  // information about grid

  int nx, ny, nz;
  CColumnGrid **pColumnGrids;
  float dx,dy,dz,x1,y1,z1; // pixel dx, dy, dz and bottom-left x,y,z respectively

  list<CSimSnap *> autoClean; // a list of CSimSnaps to delete on destruction
  // Note that ones explicitly referenced within the grid are cleaned up
  // elsewhere (in CColumnList). So this is basically for when we refine
  // and insert CGrids into the grid cells, and are left with "floating"
  // CSubsets.


  // called by constructors:
  void construct(int maxRefine, int minPart,bool displayStats);

};


#endif // __GRID_H_INCLUDED
