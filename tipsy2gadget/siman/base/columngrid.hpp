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

#ifndef __COLUMNGRID_H_INCLUDED

#define __COLUMNGRID_H_INCLUDED

class CColumnGrid {
  
public:
  CColumnGrid(CSimSnap *parentSim, float x1, float x2, int nxi, float y1, float y2, int ny, bool assign = true);
  // assign = false to just set up the structures, without assigning particles
  

  ~CColumnGrid();

  void realize();

  double getColDen(float x, float y, double units, bool neutral);
  double getColDen(int x, int y, double units, bool neutral);

  
#ifdef SIMAN_FITS

template <typename hdutype> 
void CColumnGrid::columnDensityImage(hdutype *fitsFloatHDU,double units, bool neutral) {

  // Renders a column density image into fitsFloatHDU
  // Multiplies by units first, to keep numbers nice.
  //
  // Note everything has internal DOUBLE precision, whilst the
  // rendered image is a FLOAT.

  
  long nelements = nx * ny; // number of pixels to render

  valarray<float> imgData(nelements);  // the array which will store the image
  
  

  for(int x=0;x<nx;x++) {

    for(int y=0;y<ny;y++) {

         
      imgData[x+y*ny] = (float) getColDen(x,y,units,neutral);
     
      
    } // for x 
  } // for y

  fitsFloatHDU->addKey("CTYPE1","x/kpc","X-axis");
  fitsFloatHDU->addKey("CTYPE2","y/kpc","Y-axis");
  fitsFloatHDU->addKey("CDELT1",dx,"dx per pixel");
  fitsFloatHDU->addKey("CDELT2",dy,"dy per pixel");
  fitsFloatHDU->addKey("CRPIX1",0,"Reference pixel");
  fitsFloatHDU->addKey("CRPIX2",0,"Reference pixel");
  fitsFloatHDU->addKey("CRVAL1",x1+dx/2.,"Reference pixel value");
  fitsFloatHDU->addKey("CRVAL2",y1+dy/2.,"Reference pixel value");
  fitsFloatHDU->write(1,nelements,imgData);

}

#endif // SIMAN_FITS

  virtual CTempColumnList operator[](int index);
  
  // get area type information
  
  float getDx();
  float getDy();  

protected:

  CColumnGrid(); // for compatability with CTemp

private:


  /// this is a helper function for getColDen
  /// it adds pixels in a square at distance >=dist to union u
  /// until at least minpix particles are contained within the union
  /// (poor man's sph)
  int addToUnion(int cx, int cy, int dist, int minpix, CUnion &u);

  int nx, ny;
  CColumnList **pColumnLists;
 
  float dx,dy,x1,y1; // pixel dx, dy and bottom-left x,y respectively


};

#define CTempColumnGrid CTemp<CColumnGrid, CTempColumnList > 


#endif // __COLUMNGRID_H_INCLUDED
