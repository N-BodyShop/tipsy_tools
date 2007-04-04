// grid.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __GRID_H_INCLUDED

#define __GRID_H_INCLUDED

namespace siman {
  
  class Grid;
  class GridDirectionalIterator;
  class GridWalkIterator;

  std::ostream & operator<<(std::ostream & lh, const GridWalkIterator & gw);

  class GridWalkIterator {
    
    friend class GridDirectionalIterator;
  public:
    friend std::ostream & operator<<(std::ostream &lh, const GridWalkIterator & gw);
    GridWalkIterator(); ///< construct end() iterator for any Grid
    GridWalkIterator(Grid &g); ///< construct begin() iterator for this Grid
    virtual ~GridWalkIterator() { };
    virtual GridWalkIterator & operator++();
    bool operator!=(const GridWalkIterator &g);
    virtual bool operator==(const GridWalkIterator &g);
    void operator=(const GridWalkIterator &g);
    SimSnap & operator*();
    GridDirectionalIterator directionalBegin(int dir);
    GridDirectionalIterator directionalEnd();
    coordinate<int> location() const;
    Grid & parent() const ;
  protected:
    std::stack<Grid*> gridStack;
    std::stack<int> xStack, yStack, zStack;
    std::stack<bool> descStack;
    
  };

  class GridBaseWalkIterator : public GridWalkIterator {
  public:
    GridBaseWalkIterator();
    GridBaseWalkIterator(Grid &g);
    virtual ~GridBaseWalkIterator() { };
    GridBaseWalkIterator & operator++();
  protected:
    void findNextBase();
  };

  class GridDirectionalIterator : public GridWalkIterator {
  public:
    GridDirectionalIterator(const GridWalkIterator &from, int direction); 
    GridDirectionalIterator();
    static const int xp = 0;
    static const int xm = 1;
    static const int yp = 2;
    static const int ym = 3;
    static const int zp = 4;
    static const int zm = 5;
    GridDirectionalIterator & operator++();
    void operator=(const GridDirectionalIterator &o);
    virtual bool operator==(const GridWalkIterator &g);
  protected:
    void upFrame();
    void downFrame();
    std::stack<int> downXstack, downYstack, downZstack;
    int require_nx, require_ny, require_nz;
    int direction;
  };

  class Grid : public SimSnap {
  public:

    typedef GridWalkIterator walkIterator;
    typedef GridBaseWalkIterator baseWalkIterator;
    typedef GridDirectionalIterator directionalIterator;

    Grid(SimSnap &parentSim, float x1, float x2, int nxi, float y1, float y2, int nyi, float z1, float z2, int nzi, int maxRefine =0, int minPart = 0, bool displayStats=true);
    Grid(SimSnap &parentSim, int nxi, int nyi, int nzi, int maxRefine=0, int minPart=0, bool displayStats=true); // auto-ranging version 

    virtual ~Grid();

    void realize();

    virtual ColumnGrid & operator[](int index);
    virtual SimSnap& operator[](const coordinate<int> &coords);

    GridWalkIterator walkBegin();
    GridWalkIterator walkEnd();
    GridBaseWalkIterator baseWalkBegin();
    GridBaseWalkIterator baseWalkEnd();

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

   
    void addRegion(float x1, float x2, float y1, float y2, float z1, float z2, Union *reg, float mindelta=-1.); 


    std::auto_ptr<Union> getRegion(float x1, float x2, float y1, float y2, float z1, float z2, float mindelta=-1.);
  
    std::list<std::pair<int,double> > getNearestNeighbours(const Particle &from,int number, const Metric &metric=Metric());

  private:

    // helper functions for quick getNearestNeighbours code

    void doAdd(std::queue<coordinate<int> > &processQueue,bool* traversed,int cx, int cy, int cz);
    void addIfNecessary(std::queue<coordinate<int> > &processQueue,double furthest,const Particle &from, bool *traversed,int cx, int cy, int cz);
  
    // information about grid

    int nx, ny, nz;
    ColumnGrid **pColumnGrids;
    float dx,dy,dz,x1,y1,z1; // pixel dx, dy, dz and bottom-left x,y,z respectively

    std::list<SimSnap *> autoClean; // a list of SimSnaps to delete on destruction
    // Note that ones explicitly referenced within the grid are cleaned up
    // elsewhere (in ColumnList). So this is basically for when we refine
    // and insert Grids into the grid cells, and are left with "floating"
    // Subsets.


    // called by constructors:
    void construct(int maxRefine, int minPart,bool displayStats);

  };

} // namespace siman 

#endif // __GRID_H_INCLUDED
