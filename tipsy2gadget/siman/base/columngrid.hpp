// columngrid.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __COLUMNGRID_H_INCLUDED

#define __COLUMNGRID_H_INCLUDED

namespace siman {

class ColumnGrid {
  
public:
  ColumnGrid(SimSnap &parentSim, float x1, float x2, int nxi, float y1, float y2, int ny, bool assign = true);
  // assign = false to just set up the structures, without assigning particles
  

  virtual ~ColumnGrid();

  void realize();

  double getColDen(float x, float y, std::string on_property);
  double getColDen(int x, int y, std::string on_property);

  
#ifdef SIMAN_FITS

template <typename hdutype> 
void columnDensityImage(hdutype *fitsFloatHDU,Unit out_unit,std::string on_property) {

  // Renders a column density image into fitsFloatHDU
 
  //
  // Note everything has internal DOUBLE precision, whilst the
  // rendered image is a FLOAT.

  Unit prop_unit = (*this)[0][0].getArray(on_property).getUnits();
  Unit dist_unit = (*this)[0][0].getDistanceUnits();

  double conv_ratio = (prop_unit/(dist_unit*dist_unit)).convertTo(out_unit);
  
 
  long nelements = nx * ny; // number of pixels to render

  valarray<float> imgData(nelements);  // the array which will store the image
  
  

  for(int x=0;x<nx;x++) {

    for(int y=0;y<ny;y++) {

      double result = getColDen(x,y,on_property)*conv_ratio;
     
      imgData[x+y*ny] = (float)(log(result)/2.3025851); // log10
     
      
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

  virtual ColumnList & operator[](int index);
  
  // get area type information
  
  float getDx();
  float getDy();  


private:


  /// this is a helper function for getColDen
  /// it adds pixels in a square at distance >=dist to union u
  /// until at least minpix particles are contained within the union
  /// (poor man's sph)
  int addToUnion(int cx, int cy, int dist, int minpix, Union &u);

  int nx, ny;
  ColumnList **pColumnLists;
 
  float dx,dy,x1,y1; // pixel dx, dy and bottom-left x,y respectively


};

}



#endif // __COLUMNGRID_H_INCLUDED
