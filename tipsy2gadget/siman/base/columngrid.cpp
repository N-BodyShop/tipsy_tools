//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"



//
// CColumnGrid
//



CColumnGrid::CColumnGrid(CSimSnap *parentSim, float x1i, float x2, int nxi, float y1i, float y2, int nyi, bool autoAssign) {
  
  x1 = x1i;
  y1 = y1i;

  nx = nxi;
  ny = nyi;

  pColumnLists = (CColumnList **) malloc(sizeof(void*) * nx);

  dx = (x2-x1)/(float)nx;
  dy = (y2-y1)/(float)ny;

  for(int n=0;n<nx;n++) {
    float cur_x = n * dx + x1;
    
    pColumnLists[n] = new CColumnList(parentSim, y1, y2, ny, false);
    // final parameter prevents assignment of particles, which we 
    // do ourselves later on...
  }

  int numParticles = parentSim->getNumParticles();

  if(autoAssign) {
    for(int n=0;n<numParticles;n++) {
      int x_ref,y_ref;
      CParticle *particle = parentSim->getParticle(n);
      
      x_ref = (int) ((particle->x - x1)/dx); // always rounds down - good!
      y_ref = (int) ((particle->y - y1)/dy); // ditto...
      
      CSimSnap *pSnap;
      
      if(x_ref>=0 && y_ref>=0 && x_ref<nx && y_ref<ny) {
	
	pSnap = ((*this)[x_ref][y_ref]);
	
	// this next bit is probably a bit naughty - actually, the dereferencing
	// can return different child classes of CSimSnap. But at this
	// stage in the construction, it should still be a CSubset.
	pSnap->pushParticle(n);
      }
      
      parentSim->releaseParticle(particle);
    }
  }
  
}

CColumnGrid::CColumnGrid() {
  // only for use during construction of CTemp<CColumnGrid>...
  pColumnLists = NULL;
}


double CColumnGrid::getColDen(float x, float y, double units, bool neutral) 
{
  int cx = (int)((x-x1)/dx);
  int cy = (int)((y-y1)/dy);
  return getColDen(cx,cy,units, neutral);

}

int CColumnGrid::addToUnion(int cx, int cy, int dist, int minpix, CUnion &u) {
  int numPixelsAdded = 0;
  
  int minx = cx-dist;
  int maxx=cx+dist;
  if(minx<0) minx=0;
  if(maxx>nx-1) maxx=nx-1;

  for(int x=minx; x<=maxx; x++) {
    if(cy-dist>=0) {
      
      u.add((*this)[x][cy-dist]);
      numPixelsAdded++;
     
    }
    if(cy+dist<ny) {
     
      u.add((*this)[x][cy+dist]);
      numPixelsAdded++;
     
    }
  }

  int miny = cy-dist+1;
  int maxy = cy+dist-1;
  if(miny<0) miny=0;
  if(maxy>ny-1) maxy=ny-1;

  for(int y=miny; y<=maxy; y++) {
    if(cx-dist>=0) {
      u.add((*this)[cx-dist][y]);
      numPixelsAdded++;
    }
    if(cx+dist<nx) {
      u.add((*this)[cx+dist][y]);
      numPixelsAdded++;
	
    }
  }

  if(u.getNumParticles()<minpix)
    return numPixelsAdded + addToUnion(cx,cy,dist+1,minpix,u);
  else
    return numPixelsAdded;

}

double CColumnGrid::getColDen(int x, int y, double units, bool neutral) {
  
  CSimSnap *s = (*this)[x][y];
  CUnion u((*this)[x][y]);

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

  if(neutral) mass= s->getH0Mass();
  else mass = s->getTotalMass();
  
  double colden = ((double) (mass / area)) * units;
  
  if(colden>1)
    return log(colden)/log(10.);
  else
    return 0.;

}


float CColumnGrid::getDx() {
  return dx;
}

float CColumnGrid::getDy() {
  return dy;
}

void CColumnGrid::realize() {

  // dereference initially passed simulation

  for(int n=0;n<nx;n++) {
    pColumnLists[n]->realize();
  }
}

CColumnGrid::~CColumnGrid() {

  if(pColumnLists!=NULL) {
    for(int n=0;n<nx;n++) {
      
      delete pColumnLists[n];
    }
    

    free(pColumnLists);

    pColumnLists = NULL;
  }
}

CTempColumnList CColumnGrid::operator[](int n) {
  CTempColumnList r(pColumnLists[n]);
  return r;
}
