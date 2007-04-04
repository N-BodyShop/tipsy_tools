// columngrid.cpp - part of SimAn Simulation Analysis Library
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









#include "siman.hpp"

namespace siman {

  ColumnGrid::ColumnGrid(SimSnap &parentSim, float x1i, float x2, int nxi, float y1i, float y2, int nyi, bool autoAssign) {
  
    x1 = x1i;
    y1 = y1i;

    nx = nxi;
    ny = nyi;

    pColumnLists = (ColumnList **) malloc(sizeof(void*) * nx);

    dx = (x2-x1)/(float)nx;
    dy = (y2-y1)/(float)ny;

    for(int n=0;n<nx;n++) {
    
      pColumnLists[n] = new ColumnList(parentSim, y1, y2, ny, false);
      // final parameter prevents assignment of particles, which we 
      // do ourselves later on...
    }

    int numParticles = parentSim.getNumParticles();

    if(autoAssign) {
      for(int n=0;n<numParticles;n++) {
	int x_ref,y_ref;
	Particle *particle = parentSim.getParticle(n);
      
	x_ref = (int) ((particle->x - x1)/dx); // always rounds down - good!
	y_ref = (int) ((particle->y - y1)/dy); // ditto...
      
	SimSnap *pSnap;
      
	if(x_ref>=0 && y_ref>=0 && x_ref<nx && y_ref<ny) {
	
	  pSnap = &((*this)[x_ref][y_ref]);
	
	  // this next bit is probably a bit naughty - actually, the dereferencing
	  // can return different child classes of SimSnap. But at this
	  // stage in the construction, it should still be a Subset.
	  static_cast<Subset*>(pSnap)->pushParticle(n);
	}
      

      }
    }
  
  }


  double ColumnGrid::getColDen(float x, float y, string on_prop) 
  {
    int cx = (int)((x-x1)/dx);
    int cy = (int)((y-y1)/dy);
    return getColDen(cx,cy,on_prop);

  }

  int ColumnGrid::addToUnion(int cx, int cy, int dist, int minpix, Union &u) {
    int numPixelsAdded = 0;
  
    int minx = cx-dist;
    int maxx=cx+dist;
    if(minx<0) minx=0;
    if(maxx>nx-1) maxx=nx-1;

    for(int x=minx; x<=maxx; x++) {
      if(cy-dist>=0) {
      
	u.add(&((*this)[x][cy-dist]));
	numPixelsAdded++;
     
      }
      if(cy+dist<ny) {
     
	u.add(&((*this)[x][cy+dist]));
	numPixelsAdded++;
     
      }
    }

    int miny = cy-dist+1;
    int maxy = cy+dist-1;
    if(miny<0) miny=0;
    if(maxy>ny-1) maxy=ny-1;

    for(int y=miny; y<=maxy; y++) {
      if(cx-dist>=0) {
	u.add(&((*this)[cx-dist][y]));
	numPixelsAdded++;
      }
      if(cx+dist<nx) {
	u.add(&((*this)[cx+dist][y]));
	numPixelsAdded++;
	
      }
    }

    if(u.getNumParticles()<minpix)
      return numPixelsAdded + addToUnion(cx,cy,dist+1,minpix,u);
    else
      return numPixelsAdded;

  }

  double ColumnGrid::getColDen(int x, int y, string on_prop) {
  
    SimSnap *s = &((*this)[x][y]);
    Union u(s->getParent());



    int totpixav = 1;

    if(s->getNumParticles()<10) {
      s= &u;
    
      u.add((*this)[x][y]);
      totpixav+=addToUnion(x,y,1,10,u);
    

      /*
	if(x>0) {
	totpixav++;
	u.add((*this)[x-1][y]);
	if(y>0) {
	u.add((*this)[x-1][y-1]);
	totpixav++;
	}
	if(y<ny-1) {
	u.add((*this)[x-1][y+1]);
	totpixav++;
	}
	}
	if(x<nx-1) {
	totpixav++;
	u.add((*this)[x+1][y]);
      
	if(y>0) {
	u.add((*this)[x+1][y-1]);
	totpixav++;
	}
	if(y<ny-1) {
	u.add((*this)[x+1][y+1]);
	totpixav++;
	}
	}
	if(y>0) {
	totpixav++;
	u.add((*this)[x][y-1]);
	}
	if(y<ny-1) {
	totpixav++;
	u.add((*this)[x][y+1]);
	}

      */
    }
    
    double area = dx * dy *totpixav;  // individual area of a pixel
    
    
    float mass;

  
    mass = s->getArray(on_prop).total();
  
    double colden = ((double) (mass / area));
  
   
    return colden;
    

  }


  float ColumnGrid::getDx() {
    return dx;
  }

  float ColumnGrid::getDy() {
    return dy;
  }

  void ColumnGrid::realize() {

    // dereference initially passed simulation

    for(int n=0;n<nx;n++) {
      pColumnLists[n]->realize();
    }
  }

  ColumnGrid::~ColumnGrid() {

    if(pColumnLists!=NULL) {
      for(int n=0;n<nx;n++) {
      
	delete pColumnLists[n];
      }
    

      free(pColumnLists);

      pColumnLists = NULL;
    }
  }

  ColumnList & ColumnGrid::operator[](int n) {
    return *(pColumnLists[n]);
  }

}
