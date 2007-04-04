// grid.cpp - part of SimAn Simulation Analysis Library
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
#include <typeinfo>
#include <boost/lexical_cast.hpp>

namespace siman {

  std::ostream & operator<<(std::ostream & lh, const GridWalkIterator & gw) {
    lh << gw.gridStack.size() << " ";
  }

  GridWalkIterator::GridWalkIterator( Grid & g ) {
    gridStack.push(&g);
    xStack.push(0);
    yStack.push(0);
    zStack.push(0);
    descStack.push(false);
  }

  GridWalkIterator::GridWalkIterator() {

  }
  
  coordinate<int> GridWalkIterator::location() const {
    coordinate<int> t(xStack.top(),yStack.top(),zStack.top());
    
    return t;
  }

  Grid & GridWalkIterator::parent() const {
    return *(gridStack.top());
  }

  GridWalkIterator & GridWalkIterator::operator++() {
    if(gridStack.size()==0)
      return (*this);
    Grid & g = *(gridStack.top());
    int x = xStack.top();
    int y = yStack.top();
    int z = zStack.top();
    bool descended = descStack.top();
    if(!descended) {
      if(typeid(g[x][y][z])==typeid(Grid)) {
	// descend: first remind us not to descend again
	descStack.pop();
	descStack.push(true);

	// now push on new frame
	gridStack.push(static_cast<Grid*>(&(g[x][y][z])));
	xStack.push(0);
	yStack.push(0);
	zStack.push(0);
	descStack.push(false);
	return (*this);
      }
    }

    xStack.pop();
    yStack.pop();
    zStack.pop();
    descStack.pop();
    x++;
    if(x>g.getNx()-1) {
      x=0;
      y++;
      if(y>g.getNy()-1) {
	y=0;
	z++;
	if(z>g.getNz()-1) {
	  gridStack.pop();
	  return ++(*this);
	}
      }
    }
    xStack.push(x);
    yStack.push(y);
    zStack.push(z);
    descStack.push(false);

    return (*this);
  }

  bool GridWalkIterator::operator==(const GridWalkIterator & o) {
    if(gridStack.size()!=o.gridStack.size())
      return false;
    if(gridStack.size()!=0) {
      if(gridStack.top()!=o.gridStack.top())
	return false;
      if(xStack.top()!=o.xStack.top())
	return false;
      if(yStack.top()!=o.yStack.top())
	return false;
      if(zStack.top()!=o.zStack.top())
	return false;
      if(descStack.top()!=o.descStack.top())
	return false;
    }
    return true;
  }

  SimSnap & GridWalkIterator::operator*() {
    Grid &g = *(gridStack.top());
    return g[xStack.top()][yStack.top()][zStack.top()];
  }

  bool GridWalkIterator::operator!=(const GridWalkIterator &o) {
    return !((*this)==o);
  }

  void GridWalkIterator::operator=(const GridWalkIterator &o) {
    gridStack = o.gridStack;
    xStack = o.xStack;
    yStack = o.yStack;
    zStack = o.zStack;
    descStack = o.descStack;
  }

  GridDirectionalIterator GridWalkIterator::directionalBegin(int dir) {
    return GridDirectionalIterator(*this,dir);
  }

  GridDirectionalIterator GridWalkIterator::directionalEnd() {
    return GridDirectionalIterator();
  }

  GridBaseWalkIterator::GridBaseWalkIterator() : GridWalkIterator() {
    
  }

  GridBaseWalkIterator::GridBaseWalkIterator(Grid &g) : GridWalkIterator(g) {
    findNextBase();
  }

  

  void GridBaseWalkIterator::findNextBase() {
	   
    if(descStack.size()==0)
      return; 

    while(typeid(*(*this))==typeid(Grid)) {
      GridWalkIterator::operator++();
      if(descStack.size()==0)
	return; 
    }
  }
  
  GridBaseWalkIterator & GridBaseWalkIterator::operator++() {
    GridWalkIterator::operator++();
    findNextBase();
    return(*this);
  }

  GridDirectionalIterator::GridDirectionalIterator() : GridWalkIterator() {

  }

  GridDirectionalIterator::GridDirectionalIterator(const GridWalkIterator & from, int dir) : GridWalkIterator(from), direction(dir) {
   
    require_nx = gridStack.top()->getNx();
    require_ny = gridStack.top()->getNy();
    require_nz = gridStack.top()->getNz();
  }

  
  bool GridDirectionalIterator::operator==(const GridWalkIterator & o) {
    if(gridStack.size()!=o.gridStack.size())
      return false;
    if(gridStack.size()!=0) {
      if(gridStack.top()!=o.gridStack.top())
	return false;
      if(xStack.top()!=o.xStack.top())
	return false;
      if(yStack.top()!=o.yStack.top())
	return false;
      if(zStack.top()!=o.zStack.top())
	return false;
   
    }
    return true;
  }

  
  void GridDirectionalIterator::operator=(const GridDirectionalIterator &o) {
    GridWalkIterator::operator=(o);
    direction = o.direction;
    downXstack = o.downXstack;
    downYstack = o.downYstack;
    downZstack = o.downZstack;
    require_nx = o.require_nx;
    require_ny = o.require_ny;
    require_nz = o.require_nz;
  }

  void GridDirectionalIterator::upFrame() {
  
    if(xStack.size()!=0) {
      int x = xStack.top();
      int y = yStack.top();
      int z = zStack.top();
      
      downXstack.push(x);
      downYstack.push(y);
      downZstack.push(z);
      
      gridStack.pop();
      xStack.pop();
      yStack.pop();
      zStack.pop();
    }
    
  }
  
  void GridDirectionalIterator::downFrame() {
    Grid *pG = gridStack.top();
    int x = xStack.top();
    int y = yStack.top();
    int z = zStack.top();

    while(typeid((*pG)[x][y][z])==typeid(Grid) && downXstack.size()!=0) {
      
      pG = static_cast<Grid*>(&((*pG)[x][y][z]));
      x = downXstack.top();
      y = downYstack.top();
      z = downZstack.top();
      
      downXstack.pop();
      downYstack.pop();
      downZstack.pop();
      

      switch(direction) {
      case xp:
	x=0;
	break;
      case yp:
	y=0;
	break;
      case zp:
	z=0;
	break;
      case xm:
	x=require_nx-1;
	break;
      case ym:
	y=require_ny-1;
	break;
      case zm:
	z=require_nz-1;
	break;
      }
	
  
      
      if(pG->getNx()!=require_nx || pG->getNy()!=require_ny || pG->getNz()!=require_nz)
	throw SimanException("Perverse grid structure "+boost::lexical_cast<string>(pG->getNx()));

      xStack.push(x);
      yStack.push(y);
      zStack.push(z);
      gridStack.push(pG);
    }
  }

  GridDirectionalIterator & GridDirectionalIterator::operator++() {
    if(gridStack.size()==0)
      return *this;
    int x = xStack.top();
    int y = yStack.top();
    int z = zStack.top();
    switch(direction) {
    case xp:
      x++;
      if(x>=require_nx) {
	upFrame();
	// have gone out of range in this grid, try jumping upframe, then...
 
	return operator++();
	// attempt transformation again
      }

      xStack.pop();
      xStack.push(x);
      break;
    case xm:
      x--;
      if(x<0) {
	upFrame();
       	return operator++();
      }
      xStack.pop();
      xStack.push(x);
      break;
    case yp:
      y++;
      if(y>=require_ny) {
	upFrame();
	return operator++();
      }
      
      yStack.pop();
      yStack.push(y);
      break;
    case ym:
      y--;
      if(y<0) {
	upFrame();
	return operator++();
      }
      
      yStack.pop();
      yStack.push(y);
      break;
    case zp:
      z++;
      if(z>=require_nz) {
	upFrame();
	return operator++();
      }
      
      zStack.pop();
      zStack.push(z);
      break;
    case zm:
      z--;
      if(z<0) {
	upFrame();
	return operator++();
      }
      
      zStack.pop();
      zStack.push(z);
      break;
    }

    
    
    // go down the tree as far as possible...
    downFrame();

    return (*this);
  }

  Grid::Grid( SimSnap &parentSim, float x1i, float x2, int nxi, float y1i, float y2, int nyi, float z1i, float z2, int nzi, int maxRefine, int minPart, bool displayStats) {
  
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

    pAutoParent = &parentSim;

    construct(maxRefine,minPart,displayStats);
  }


  Grid::Grid(SimSnap &parentSim, int nxi, int nyi, int nzi, int maxRefine, int minPart, bool displayStats)
    : nx(nxi), ny(nyi), nz(nzi) {

    // constructs a grid with automatic ranges (x1->x2 etc) determined by simulation size
    //
    // otherwise the same as other constructor
  
    float x2,y2,z2;

    parentSim.getExactBoundaries(x1,x2,y1,y2,z1,z2);
  
    // ensure all particles fall inside boundaries
    x1-=(x2-x1)*0.01;
    x2+=(x2-x1)*0.01;
    y1-=(y2-y1)*0.01;
    y2+=(y2-y1)*0.01;
    z1-=(z2-z1)*0.01;
    z2+=(z2-z1)*0.01;
    
    dx = (x2-x1)/(float)nx; 
    dy = (y2-y1)/(float)ny; 
    dz = (z2-z1)/(float)nz;
  
    pAutoParent=&parentSim;

    construct(maxRefine,minPart,displayStats);
  
  }


  void Grid::construct(int maxRefine, int minPart, bool displayStats) {

    pColumnGrids = (ColumnGrid **) malloc(sizeof(void*) * nx);


    for(int n=0;n<nx;n++) {
        
      pColumnGrids[n] = new ColumnGrid(*pAutoParent, y1, y1+ny*dy, ny, z1, z1+nz*dz, nz, false);
      // final parameter here prevents ColumnGrid data from being built up, only
      // the structures are put in place
    
    }

    int numParticles = pAutoParent->getNumParticles();

    int x_ref,y_ref,z_ref;

    for(int n=0;n<numParticles;n++) {
    
      const Particle *particle = pAutoParent->getConstParticle(n);
    
      x_ref = (int) ((particle->x - x1)/dx); // always rounds down - good!
      y_ref = (int) ((particle->y - y1)/dy); // ditto...
      z_ref = (int) ((particle->z - z1)/dz); 

      SimSnap *pSnap;
    
      if(x_ref>=0 && y_ref>=0 && z_ref>=0 && x_ref<nx && y_ref<ny && z_ref<nz) {
      
	pSnap = &((*this)[x_ref][y_ref][z_ref]);
      
	// this next bit is probably a bit naughty - actually, the dereferencing
	// can return different child classes of SimSnap. But at this
	// stage in the construction, it should still be a Subset.
	(static_cast<Subset*>(pSnap))->pushParticle(n);
      }
    

    }

    if(maxRefine>0 && minPart>3) {
    
      // Walk the grid and see if any cells could do with
      // refining...
    
      for(x_ref=0;x_ref<nx;x_ref++) {
	for(y_ref=0;y_ref<ny;y_ref++) {
	  for(z_ref=0;z_ref<nz;z_ref++) {
	  
	    SimSnap *parent = &((*this)[x_ref][y_ref][z_ref]);

	    if(parent->getNumParticles()>(unsigned int)minPart) {

	      autoClean.push_back(parent);
	    
	      // will remove reference from grid, so need to clean up later
	    
	      // Create new grid, not forgetting to let it know it's a level 
	      // down from us. Otherwise refinement criterion (minPart)
	      // stays the same.
	    
	      Grid *replace = new Grid(*parent,x1+x_ref*dx,x1+(x_ref+1)*(dx),nx,
				       y1+y_ref*dy,y1+(y_ref+1)*dy,ny,
				       z1+z_ref*dz,z1+(z_ref+1)*dz,nz,
				       maxRefine-1,
				       minPart,false);
	    
	    
	      // This requires us to be a friend class of ColumnList
	      // - not ideal. Don't think there is a more sensible way.

	      ((*this)[x_ref][y_ref]).pColumnContents[z_ref] = replace;
	   
	    } // if parent needs refining

	  } // for z

	} // for y

      } // for x
    } else {

    } // if we can refine anything

    if(displayStats && getVerbose()>=2) {
      cerr << "Grid: mean refine depth = " << getMeanRefineDepth() << endl;
      cerr << "Grid: mean leaf particles = " << getMeanLeafParticles() << endl;
    } 
  
  }

  float Grid::getMeanRefineDepth() {
    int num_cells = nx*ny*nz;
    float rd=0;

    for(int x=0;x<nx;x++) {
      for(int y=0;y<ny;y++) {
	for(int z=0;z<nz;z++) {
	  SimSnap *cell = &((*this)[x][y][z]);
	  if(typeid(*cell)==typeid(Grid)) {
	    rd+=(1.+((Grid*)cell)->getMeanRefineDepth())/((float)num_cells);
	  } else {
	    rd+=1./((float)num_cells);
	  }
	}
      }
    }
    return rd;

  }

  GridWalkIterator Grid::walkBegin() {
    return GridWalkIterator(*this);
  }

  GridWalkIterator Grid::walkEnd() {
    return GridWalkIterator();
  }

  GridBaseWalkIterator Grid::baseWalkBegin() {
    return GridBaseWalkIterator(*this);
  }

  GridBaseWalkIterator Grid::baseWalkEnd() {
    return GridBaseWalkIterator();
  }

  float Grid::getMeanLeafParticles() {
    int num_cells = nx*ny*nz;
    float rd=0;

    for(int x=0;x<nx;x++) {
      for(int y=0;y<ny;y++) {
	for(int z=0;z<nz;z++) {
	  SimSnap *cell = &((*this)[x][y][z]);
	  if(typeid(*cell)==typeid(Grid)) {
	    rd+=((Grid*)cell)->getMeanLeafParticles()/((float)num_cells);
	  } else {
	    rd+=(cell->getNumParticles())/((float)num_cells);
	  }
	}
      }
    }
    return rd;

  }



  void Grid::doAdd(queue<coordinate<int> > &processQueue,bool* traversed,int cx, int cy, int cz) {

    // helper for getNearestNeighbours
    //
    // adds relevent SimSnap to queue IF it has not already
    // been processed, and exists!

    if(cx<0 || cx>=nx || cy<0 || cy>=ny || cz<0 || cz>=nz) {
      // cerr << "Grid: reject request for ("<<nx<<","<<cy<<","<<cz<<")" << endl;
      return;
      // out of range
    
    }

    if(traversed[cx+cy*nx+cz*nx*ny]==true) return;
    // already looked at (or in queue)

    traversed[cx+cy*nx+cz*nx*ny] = true;

    processQueue.push(coordinate<int>(cx,cy,cz));
  
  }


  void Grid::addIfNecessary(queue<coordinate<int> > &processQueue,double furthest,const Particle &from, bool *traversed,int cx, int cy, int cz) {
  
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

  list< pair<int,double> > Grid::getNearestNeighbours(const Particle &from, int number, const Metric &metric) {
  
    // N.B. - a more efficient version of this code might use a PriorityQueue to
    // keep track of the closest possible distance for particles in each box,
    // and terminate the loop when there are no more queued which could
    // change the list.
    //
    // Currently, it only checks at the point of queueing whether the box
    // could change the list, which is an OK first approximation

    if(metric(Particle(0,0,0),Particle(3,4,0))<5) {
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
      SimSnap *consider = &((*this)[coords]);
    
      processQueue.pop();

      // 1. Get n nearest neighbours in cell at top of processQueue
      list< pair<int,double> > current;
      current = consider->getNearestNeighbours(from,number);
   

      // 1.5 Dereference identifiers as returned by getNearestNeighbours
      //     (they refer to child cell IDs, not parent IDs)

      for(i=current.begin();i!=current.end();i++) {
	(*i).first = consider->deReference((*i).first);
      }

      // 2. insertion sort as in SimSnap, but with entire list from 
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


  float Grid::getDx() {
    return dx;
  }

  float Grid::getDy() {
    return dy;
  }

  float Grid::getDz() {
    return dz;
  }

  float Grid::getX1() {
    return x1;
  }

  float Grid::getY1() {
    return y1;
  }

  float Grid::getZ1() {
    return z1;
  }

  int Grid::getNx() {
    return nx;
  }

  int Grid::getNy() {
    return ny;
  }

  int Grid::getNz() {
    return nz;
  }

  void Grid::realize() {

    // dereference initially passed simulation

    for(int n=0;n<nx
	  ;n++) {
      pColumnGrids[n]->realize();
    }
  }

  Grid::~Grid() {

    if(pColumnGrids!=NULL) {
      for(int n=0;n<nx;n++) {
      
	delete pColumnGrids[n];
      }
    

      free(pColumnGrids);

      pColumnGrids = NULL;
    }

    list<SimSnap*>::iterator i;
    for(i=autoClean.begin();i!=autoClean.end();i++) {
      delete (*i);
    }

  }

  ColumnGrid & Grid::operator[](int n) {
    return *(pColumnGrids[n]);
  }

  SimSnap& Grid::operator[](const coordinate<int> &coord) {
    return (*this)[coord.x][coord.y][coord.z];
  }

  void Grid::addRegion(float bx1, float bx2, float by1, float by2, float bz1, float bz2, Union *region, float mindelta) {
    
    int nx1 = (int) ((bx1-x1)/dx);
    int nx2 = (int) ((bx2-x1)/dx)+1;
    int ny1 = (int) ((by1-y1)/dy);
    int ny2 = (int) ((by2-y1)/dy)+1;
    int nz1 = (int) ((bz1-z1)/dz);
    int nz2 = (int) ((bz2-z1)/dz)+1;

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

    //  cout << nx1 << " " << nx2 << " " << ny1 << " " << ny2 << " " << nz1 << " " << nz2 << endl;
    if(nx2-nx1>1 && ny2-ny1>1 && nz2-nz1>1 || dx<mindelta) {
      // staying on this level is good enough
      for(int x=nx1;x<=nx2;x++) {
	for(int y=ny1;y<=ny2;y++) {
	  for(int z=nz1;z<=nz2;z++) {
	    region->add((*this)[x][y][z]);
	  }
	}
      }
    } else {
      // better to refine
      for(int x=nx1;x<=nx2;x++) {
	for(int y=ny1;y<=ny2;y++) {
	  for(int z=nz1;z<=nz2;z++) {
	    if(typeid((*this)[x][y][z])==typeid(Grid)) {
	    
	      Grid *pG = static_cast<Grid*>(&((*this)[x][y][z]));
	      pG->addRegion(bx1,bx2,by1,by2,bz1,bz2,region);
	    } else { 
	    
	      region->add((*this)[x][y][z]);
	    }
	  }
	}
      }
    }

  }

  auto_ptr<Union> Grid::getRegion(float bx1, float bx2, float by1, float by2, float bz1, float bz2, float mindelta) {
  
    Union *region = new Union(pAutoParent);
    addRegion(bx1,bx2,by1,by2,bz1,bz2,region);
    return auto_ptr<Union>(region);

  }

}
