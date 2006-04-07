//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "siman.hpp"

string CGrid::className() {
  return "CGrid";
}

CGrid::CGrid( CSimSnap *parentSim, float x1i, float x2, int nxi, float y1i, float y2, int nyi, float z1i, float z2, int nzi, int maxRefine, int minPart, bool displayStats) {
  
  // Constructs a grid with obvious parameters (covers x1->x2 in nx cells, etc)
  //
  // If any cells contain more than minPart, the grid autorefines up
  // to maxRefine levels (default: maxRefine=0, so no refinement)

  
  x1 = x1i;
  y1 = y1i;
  z1 = z1i;

  nx = nxi;
  ny = nyi;
  nz = nzi;

  dx = (x2-x1)/(float)nx; 
  dy = (y2-y1)/(float)ny; 
  dz = (z2-z1)/(float)nz;

  pAutoParent = parentSim;

  construct(maxRefine,minPart,displayStats);
}


CGrid::CGrid(CSimSnap *parentSim, int nxi, int nyi, int nzi, int maxRefine, int minPart, bool displayStats)
  : nx(nxi), ny(nyi), nz(nzi) {

  // constructs a grid with automatic ranges (x1->x2 etc) determined by simulation size
  //
  // otherwise the same as other constructor
  
  float x2,y2,z2;

  parentSim->getExactBoundaries(x1,x2,y1,y2,z1,z2);
  
  dx = (x2-x1)/(float)nx; 
  dy = (y2-y1)/(float)ny; 
  dz = (z2-z1)/(float)nz;
  
  pAutoParent=parentSim;

  construct(maxRefine,minPart,displayStats);
  
}


void CGrid::construct(int maxRefine, int minPart, bool displayStats) {

  pColumnGrids = (CColumnGrid **) malloc(sizeof(void*) * nx);


  for(int n=0;n<nx;n++) {
    float cur_x = n * dx + x1;
    
    pColumnGrids[n] = new CColumnGrid(pAutoParent, y1, y1+ny*dy, ny, z1, z1+nz*dz, nz, false);
    // final parameter here prevents ColumnGrid data from being built up, only
    // the structures are put in place
    
  }

  int numParticles = pAutoParent->getNumParticles();

  int x_ref,y_ref,z_ref;

  for(int n=0;n<numParticles;n++) {
    
    CParticle *particle = pAutoParent->getParticle(n);
    
    x_ref = (int) ((particle->x - x1)/dx); // always rounds down - good!
    y_ref = (int) ((particle->y - y1)/dy); // ditto...
    z_ref = (int) ((particle->z - z1)/dz); 

    CSimSnap *pSnap;
    
    if(x_ref>=0 && y_ref>=0 && z_ref>=0 && x_ref<nx && y_ref<ny && z_ref<nz) {
      
      pSnap = ((*this)[x_ref][y_ref][z_ref]);
      
      // this next bit is probably a bit naughty - actually, the dereferencing
      // can return different child classes of CSimSnap. But at this
      // stage in the construction, it should still be a CSubset.
      pSnap->pushParticle(n);
    }
    
    pAutoParent->releaseParticle(particle);
  }

  if(maxRefine>0 && minPart>3) {
    
    // Walk the grid and see if any cells could do with
    // refining...
    
    for(x_ref=0;x_ref<nx;x_ref++) {
      for(y_ref=0;y_ref<ny;y_ref++) {
	for(z_ref=0;z_ref<nz;z_ref++) {
	  
	  CSimSnap *parent = ((*this)[x_ref][y_ref][z_ref]);

	  if(parent->getNumParticles()>minPart) {

	    autoClean.push_back(parent);
	    
	    // will remove reference from grid, so need to clean up later
	    
	    // Create new grid, not forgetting to let it know it's a level 
	    // down from us. Otherwise refinement criterion (minPart)
	    // stays the same.
	    
	    CGrid *replace = new CGrid(parent,x1+x_ref*dx,x1+(x_ref+1)*(dx),nx,
				       y1+y_ref*dy,y1+(y_ref+1)*dy,ny,
				       z1+z_ref*dz,z1+(z_ref+1)*dz,nz,
				       maxRefine-1,
				       minPart,false);
	    
	    // The following is a bodge, requiring
	    //
	    //  (1) knowledge of inner workings (CTempColumnList, which is
	    //      an alias for a CTemp with filled in template parameters)
	    //  (2) the fact that we are a friend class of CColumnList, so can
	    //      access and overwrite its private member pColumnContents
	    //
	    // Sorry! There is probably a better way?
	    
	    CTempColumnList ref = (CTempColumnList) ((*this)[x_ref][y_ref]);
	    ref.realV->pColumnContents[z_ref] = replace;
	   
	  } // if parent needs refining

	} // for z

      } // for y

    } // for x
	} else {

  } // if we can refine anything

  if(displayStats) {
    cerr << "CGrid: mean refine depth = " << getMeanRefineDepth() << endl;
    cerr << "CGrid: mean leaf particles = " << getMeanLeafParticles() << endl;
  } 
  
}

float CGrid::getMeanRefineDepth() {
  int num_cells = nx*ny*nz;
  float rd=0;

  for(int x=0;x<nx;x++) {
    for(int y=0;y<ny;y++) {
      for(int z=0;z<nz;z++) {
	CSimSnap *cell = (*this)[x][y][z];
	if(cell->className()=="CGrid") {
	  rd+=(1.+((CGrid*)cell)->getMeanRefineDepth())/((float)num_cells);
	} else {
	  rd+=1./((float)num_cells);
	}
      }
    }
  }
  return rd;

}

float CGrid::getMeanLeafParticles() {
 int num_cells = nx*ny*nz;
  float rd=0;

  for(int x=0;x<nx;x++) {
    for(int y=0;y<ny;y++) {
      for(int z=0;z<nz;z++) {
	CSimSnap *cell = (*this)[x][y][z];
	if(cell->className()=="CGrid") {
	  rd+=((CGrid*)cell)->getMeanLeafParticles()/((float)num_cells);
	} else {
	  rd+=(cell->getNumParticles())/((float)num_cells);
	}
      }
    }
  }
  return rd;

}



void CGrid::doAdd(queue<coordinate<int> > &processQueue,bool* traversed,int cx, int cy, int cz) {

  // helper for getNearestNeighbours
  //
  // adds relevent CSimSnap to queue IF it has not already
  // been processed, and exists!

  if(cx<0 || cx>=nx || cy<0 || cy>=ny || cz<0 || cz>=nz) {
    // cerr << "CGrid: reject request for ("<<nx<<","<<cy<<","<<cz<<")" << endl;
    return;
    // out of range
    
  }

  if(traversed[cx+cy*nx+cz*nx*ny]==true) return;
  // already looked at (or in queue)

  traversed[cx+cy*nx+cz*nx*ny] = true;

  processQueue.push(coordinate<int>(cx,cy,cz));
  
}


void CGrid::addIfNecessary(queue<coordinate<int> > &processQueue,double furthest,const CParticle &from, bool *traversed,int cx, int cy, int cz) {
  
  // The cube of interest occupies the space (cx->cx+1)*(cy->cy+1)*(cz->cz+1)
  //
  // we want to know the distance from the particle "from" to the 
  // nearest point within the cube of interest. If this is closer
  // than "furthest", the cube needs to be queued since it may
  // contain particles of interest.
  
  // Here we calculate the index coordinates of the given point,
  // using proper rounding (hence +0.5) rather than index rounding
  // which is used elsewhere. Hence we get the closest grid point,
  // rather than the "box bottom left corner"

  int pcx = (int)((from.x-x1)/dx+0.5);
  int pcy = (int)((from.y-y1)/dy+0.5);
  int pcz = (int)((from.z-z1)/dz+0.5);
  
  int ccx, ccy, ccz; // index coordinates of closest corner
  
  if(pcx>cx) ccx=cx+1; else ccx=cx;
  if(pcy>cy) ccy=cy+1; else ccy=cy;
  if(pcz>cz) ccz=cz+1; else ccz=cz;
  
  // Now we recreate the coordinates of the closest point on the
  // cube, starting with the corner but adding back in 
  // information about dimensions where the point lies within the
  // boundaries of the planes defined by the continuation of the
  // faces

  float px,py,pz;
  px = x1+((float)ccx)*dx;
  py = y1+((float)ccy)*dy;
  pz = z1+((float)ccz)*dz;

  // at this point, we need the "cube" coordinates rather than
  // the nearest intersection coordinates
  pcx = (int)((from.x-x1)/dx);
  pcy = (int)((from.y-y1)/dy);
  pcz = (int)((from.z-z1)/dz);

  if(pcx==cx) px=from.x;
  if(pcy==cy) py=from.y;
  if(pcz==cz) pz=from.z;

  // calculate distance:
  
  float distance = sqrt((px-from.x)*(px-from.x)+(py-from.y)*(py-from.y)+(pz-from.z)*(pz-from.z));
  
  
  
  if(distance<furthest) {
    //    cerr << "Must add: " << px << " " << py << " " << pz << " " << distance << " " << furthest << endl;
    doAdd(processQueue,traversed,cx,cy,cz);
  }

}

list< pair<int,double> > CGrid::getNearestNeighbours(const CParticle &from, const int number, const CMetric &metric) {
  
  // N.B. - a more efficient version of this code might use a PriorityQueue to
  // keep track of the closest possible distance for particles in each box,
  // and terminate the loop when there are no more queued which could
  // change the list.
  //
  // Currently, it only checks at the point of queueing whether the box
  // could change the list, which is an OK first approximation

  if(metric(CParticle(0,0,0),CParticle(3,4,0))<5) {
    cerr << "Temporary error. Incompatible metric." << endl;
  }

  bool *traversed = new bool[nx*ny*nz];

  for(int n=0;n<nx*ny*nz;n++)
    traversed[n]=false;

  compareSecond comp; // comparing object for insertion sort (STL)

  int cx = (int)((from.x-x1)/dx);
  int cy = (int)((from.y-y1)/dy);
  int cz = (int)((from.z-z1)/dz);

  if(cx<0) cx=0;
  if(cy<0) cy=0;
  if(cz<0) cz=0;
  if(cx>=nx) cx=nx-1;
  if(cy>=ny) cy=ny-1;
  if(cz>=nz) cz=nz-1;

  list<pair<int,double> > candidateList; // will eventually become our final list


  // start with a queue of boxes to check: just the one which the given particle
  // is actually in (or closest to, if it's outside the box)

  queue<coordinate<int> > processQueue;
  processQueue.push(coordinate<int>(cx,cy,cz));
  traversed[cx+cy*nx+cz*nx*ny]=true; 
       // to prevent this box reentering the queue


  int n=0; // used within loop
  list< pair<int,double> >::iterator i; // ditto

  while(!processQueue.empty()) {

    coordinate <int> coords = processQueue.front();
    CSimSnap *consider = (*this)[coords];
    
    processQueue.pop();

    // 1. Get n nearest neighbours in cell at top of processQueue
    list< pair<int,double> > current;
    current = consider->getNearestNeighbours(from,number);
   

    // 1.5 Dereference identifiers as returned by getNearestNeighbours
    //     (they refer to child cell IDs, not parent IDs)

    for(i=current.begin();i!=current.end();i++) {
      (*i).first = consider->deReference((*i).first);
    }

    // 2. insertion sort as in CSimSnap, but with entire list from 
    // chosen cell

    // cerr << ": " << "(" << processQueue.size() << ") " << current.size() << "+" << candidateList.size() << endl;

    candidateList.merge(current,comp);
    
    // inspect our list to see if we need to truncate it
    n=candidateList.size();
    
    if(n>number) {
      
      for(int z=0;z<n-number;z++) {
	candidateList.pop_front();
      }
      // candidateList.erase(candidateList.begin(),i);
      
      
    }
    

    // Now add the neighbours to the queue for processing next, if necessary

    if(n<number) {

      // need to find more particles - MUST extend search
      // 
      // doAdd has error checking - just try to add all
      // surrounding blocks, even if they don't exist!
      //
      // doAdd will also not add a block that is already
      // looked at (flagged in traversed)
      // 
      // This also means doAdd will not look at the current
      // block.
      for(int ix=-1;ix<=1;ix++) {
	for(int iy=-1;iy<=1;iy++) {
	  for(int iz=-1;iz<=1;iz++) {
	    doAdd(processQueue,traversed,coords.x+ix,coords.y+iy,coords.z+iz);
	  }
	}
      }
   

    } else {

      // have enough particles, but could there be other blocks
      // which have particles closer to the chosen point?
      //
      // If it is necessary to check another block, addIfNecessary
      // calls doAdd - so all appropriate error checking is therefore
      // handled.

      double furthest = (*candidateList.begin()).second;
      // cerr << "-" << furthest << endl ;
       for(int ix=-1;ix<=1;ix++) {
	for(int iy=-1;iy<=1;iy++) {
	  for(int iz=-1;iz<=1;iz++) {
	    addIfNecessary(processQueue,furthest,from,traversed,coords.x+ix,coords.y+iy,coords.z+iz);
	  } // for ix
	} // for iy
       } // for iz
       
    } // else

    
  } // while(!processQueue.empty())
  
  // cout << endl << "--EXIT--" << endl << endl;
  /*
  cout << "TRAVERSED LIST: "<< endl;
  for(int x=0;x<nx;x++) {
    for(int y=0;y<ny;y++) {
      for(int z=0;z<nz;z++) {
	int n = x+y*nx+z*nx*ny;
	if(traversed[n]==true) {
	  cout << x <<"," << y << "," << z << endl;
	}
      }
    }
  }
  cout << endl;*/
  delete[] traversed;
  return candidateList;
  
}


float CGrid::getDx() {
  return dx;
}

float CGrid::getDy() {
  return dy;
}

float CGrid::getDz() {
  return dz;
}

float CGrid::getX1() {
  return x1;
}

float CGrid::getY1() {
  return y1;
}

float CGrid::getZ1() {
  return z1;
}

int CGrid::getNx() {
  return nx;
}

int CGrid::getNy() {
  return ny;
}

int CGrid::getNz() {
  return nz;
}

void CGrid::realize() {

  // dereference initially passed simulation

  for(int n=0;n<nx
;n++) {
    pColumnGrids[n]->realize();
  }
}

CGrid::~CGrid() {

  if(pColumnGrids!=NULL) {
    for(int n=0;n<nx;n++) {
      
      delete pColumnGrids[n];
    }
    

    free(pColumnGrids);

    pColumnGrids = NULL;
  }

  list<CSimSnap*>::iterator i;
  for(i=autoClean.begin();i!=autoClean.end();i++) {
    delete (*i);
  }

}

CTempColumnGrid CGrid::operator[](int n) {
  CTempColumnGrid r(pColumnGrids[n]);
  return r;
}

CSimSnap* CGrid::operator[](const coordinate<int> &coord) {
  return (*this)[coord.x][coord.y][coord.z];
}

CSimSnap* CGrid::getRegion(float bx1, float bx2, float by1, float by2, float bz1, float bz2) {

  int nx2 = (int) ((bx1-x1)/dx);
  int nx1 = (int) ((bx2-x1)/dx);
  int ny2 = (int) ((by1-y1)/dy);
  int ny1 = (int) ((by2-y1)/dy);
  int nz2 = (int) ((bz1-z1)/dz);
  int nz1 = (int) ((bz2-z1)/dz);

  // only access data in range!

  if(nx1<0) nx1=0;
  if(nx2<0) nx2=0;
  if(ny1<0) ny1=0;
  if(ny2<0) ny2=0;
  if(nz1<0) nz1=0;
  if(nz2<0) nz2=0;
  
  if(nx1>nx-1) nx1=nz-1;
  if(nx2>nx-1) nx2=nz-1;
  if(ny1>ny-1) ny1=ny-1;
  if(ny2>ny-1) ny2=ny-1;
  if(nz1>nz-1) nz1=nz-1;
  if(nz2>nz-1) nz2=nz-1;

  CUnion *region = new CUnion(pAutoParent);

  for(int x=nx1;x<=nx2;x++) {
    for(int y=ny1;y<=ny2;y++) {
      for(int z=nz1;z<=nz2;z++) {
	region->add((*this)[x][y][z]);
      }
    }
  }
  

  return region;

}

