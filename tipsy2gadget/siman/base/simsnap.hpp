// simsnap.hpp - part of SimAn Simulation Analysis Library
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








#ifndef __SIMSNAP_H_INCLUDED

#define __SIMSNAP_H_INCLUDED

#include "siman.hpp"

namespace siman {

  class Subset;
  class Filter;


  /// SimSnap is effectively an abstract class defining an interface for
  /// extracting information from a simulation. It is not actually a pure 
  /// virtual class, however, since it provides internal support for "derived"
  /// classes (e.g. if you create a \ref Subset of a file) - and also some
  /// higher level features like \ref getTotalMass(). 
  ///
  /// Note that, ultimately, SimSnap might be renamed SimSnapSPH to allow for
  /// a more generalised base class which can represent more general forms
  /// of simulation, e.g. Eulerian grid-based simulations.

  class SimSnap : public SimanObject {
    friend class SimanArray;
  public:

    SimSnap();
    virtual ~SimSnap();

  
    /// getParticle returns a pointer to the particle identified by id.
    /// Note that a pointer has explicitly been chosen over a reference here, since
    /// getParticle is likely to be used in loops, where reseating (i.e. reassignment of
    /// the pointer) is desirable. 
    ///
    /// Some child classes of SimSnap will create objects for each call to getParticle,
    /// if they are performing transformations on the data. However, they will *never*
    /// create duplicates (i.e. if you get particle i then get particle i again, no
    /// extra memory will be used.) You can guarantee minimum memory usage by calling
    /// releaseParticle(id) within your loop. However, there is a speed performance 
    /// hit associated with this.
    ///
    /// You can also optionally call releaseAllParticles() when you are done, to free
    /// up any memory that can be freed. However, this will be called in object destructors
    /// where necessary.
    ///
    /// @param id - ID of particle, in range 0..getNumParticles()-1

    virtual Particle * getParticle(unsigned int id);
    virtual Particle & getParticleRef(unsigned int id);

    /// @see getParticle
    /// returns a const Particle * for read-only routines. This saves CPU for Transformer classes
    /// which need to de-transform to write-back to original underlying SimSnaps.
    virtual const Particle *getConstParticle(unsigned int id) const;
    virtual const Particle &getConstParticleRef(unsigned int id) const;

    /// Releases a particle returned by getParticle(). Optional - @see getParticle()
    virtual void releaseParticle(unsigned int id);

    /// Releases all particles previously returned by getParticle(). Optional - @see getParticle()
    virtual void releaseAllParticles();

    /// This is incremented each time a change is made to the SimSnap object, so that you can
    /// tell if anything has happened between accesses. Most useful for Visualiser classes.
    /// Note that it is not totally reliable; for instance, if you getParticle() and hold the
    /// pointer, then make a change to it later, the version is only incremented on the getParticle
    /// call. There is no obvious way around this without severe CPU overheads.
    long getVersion() const;

    /// Returns the number of particles in this object. The IDs are zero-based, so run to getNumParticles()-1.
    virtual unsigned int getNumParticles() const;


    virtual Unit getDistanceUnits() const;   ///< get the units of distances (i.e. particles' positions)
    virtual Unit getVelocityUnits() const;   ///< get the units of velocities 

    virtual Unit getDensityUnits() const;
    virtual Unit getEnergyUnits() const;
    virtual Unit getMassUnits() const;

    /// returns this simulation's value of Y, the helium mass fraction
    float getHeliumMassFrac() const;

    /// convertUnits converts all standard units from whatever they are currently
    /// to those specified in the passed units. If not specified (or a dimensionless
    /// unit is passed), they default to the units shown in brackets.
    /// @param distance - [kpc]
    /// @param mass - [msol]
    /// @param density - [msol kpc^-3]
    /// @param energy - [km^2 s^-2]
    virtual void convertUnits(Unit distance = Unit(), 
			      Unit mass = Unit(),
			      Unit velocity = Unit(),
			      Unit density = Unit(),
			      Unit energy = Unit());

    virtual void UFromTemp() throw(UnitsError); ///< use Temp and NE to determine U for gas particles
    virtual void TempFromU() throw(UnitsError); ///< use U and NE to determine Temp for gas particles
    virtual void makeMu();  ///< create "mu" array from ne data (mean atomic mass)
    virtual void makeP(); ///< create "P" (pressure) array from other thermodynamic data 
   
    /// returns the box size of this simulation according to the original file
    virtual float getBoxSize() const;
    virtual float getApparentBoxSize() const; ///< calculates from maximal particle positions the apparent box size of this sim
    virtual float getMaximalVelocity() const; ///< calculates maximal velocity magnitude
    virtual float getHubble() const; ///< returns h=H_0/(100 km/s/Mpc) according to the original file
    virtual float getRedshift() const; ///< returns the redshift z
    virtual float getOmegaM0() const; ///< returns omega_M0 
    virtual float getOmegaLambda0() const; ///< returns omega_Lambda0 

    virtual void setHubble(float h); ///< sets h=H_0/(100 km/s/Mpc)
    virtual void setBoxSize(float b); ///< sets the box size in units given by getDistanceUnits()
    virtual void setRedshift(float z); ///< sets the redshift z
    virtual void setOmegaM0(float om0); ///< sets OmegaM0
    virtual void setOmegaLambda0(float ol0); /// <sets OmegaLambda0

    virtual float scaleToHubble1e10yr(float a=-1.) const; ///< returns hubble constant in units of [1/10^10 year] for given scalefactor and current cosmology. Negative a=>evaluate for current redshift.
    virtual float scaleToHubble(float a=-1.) const; ///< returns hubble constant in [vel/distance] units for given scalefactor and current cosmology. Negative a=>evaluate for current redshift.
    virtual float scaleToAge(float a=-1.) const; ///< calculates age in [distance/vel] units, in current cosmology [integrates numerically]; Negative a=>evaluate for current redshift. N.B. NOT CURRENTLY VERY ACCURATE
    virtual float scaleToAge1e10yr(float a=-1.) const; ///< calculates age in [10^10 year] units, in current cosmology [integrates numerically]; Negative a=>evaluate for current redshift. N.B. NOT CURRENTLY VERY ACCURATE
    virtual float scaleToAngDiDist(float a_em) const; ///< calculates current angular diameter distance of an object placed at scalefactor a_em in [distance] units
    virtual float scaleToAngDiDist(float a_now, float a_em) const; ///< calculates angular diameter distance from scalefactor a_now of an object placed at scalefactor a_em in [distance] units
    
    virtual float scaleToComovingDist(float a_em) const; ///< calculates comoving distance of an object placed at scalefactor a_em in [distance] units
    virtual float scaleToComovingDist(float a_now, float a_em) const; ///< calculates comoving distance from scalefactor a_now of an object placed at scalefactor a_em in [distance] units
    
    
    
    virtual void addHubbleFlow(); ///< add the hubble flow according to getHubble()
    virtual void subtractHubbleFlow(); ///< subtract the hubble flow according to getHubble()
    virtual std::vector<double> makeNaiveAbsProfile(int n_elements, double min_lambda, double max_lambda, double lambda_cen, double gamma, double osc_strength) const;
    virtual std::vector<double> makeSPHAbsProfile(float x, float y, int n_elements, double min_lambda, double max_lambda, double lambda_cen, double gamma, double osc_strength) const;

    ///\defgroup ExtraData Extra data handlers:
    /// Methods for associating extra data with a SimSnap
    ///@{
    virtual SimanArray & createArray(std::string name, std::string fullname, Unit inunits = Unit());
    virtual void destroyArray(std::string name);
    virtual SimanArray & getArray(std::string name);
    virtual SimanArray & getArray(int n);

    virtual const SimanArray & getConstArray(std::string name) const ;
    virtual const SimanArray & getConstArray(int n) const;

    virtual std::string getArrayLongName(std::string name) const;
    virtual std::string getArrayLongName(int n)const ;
    virtual std::string getArrayName(int n) const;
    virtual Unit getArrayUnits(std::string name) const;
    virtual Unit getArrayUnits(int n) const;
    virtual int getNumArrays() const;
    virtual int getArrayIndex(std::string name) const;
    ///@}

  
    /// copy a simsnap, disinheriting it from any parents. Caller is responsible for
    /// deleting returned object.
    virtual SimSnap * copy() const;

    /// perform a transformation on this simsnap 
    virtual void transform(const Transformation &t);

    /// perform a transformation on a copy of this object, and return the transformed
    /// copy. Caller is responsible for deleting returned object.
    virtual SimSnap * copyTransform(const Transformation &t) const;

    /// Create a subset from a filter. Equivalent to using constructor Subset(..)
    virtual Subset * subset(const Filter &f);
  
    /// Create a new simsnap from a subset corresponding to a filter. Equivalent to using constructor BaseSimSnap(Subset(..))
    virtual SimSnap * copySubset(const Filter &f) const;

    /// deReference returns the particle which i refers to in the nth parent
    /// If it reaches the base of the inheritance tree before getting to nth parent,
    /// it returns the value in the base. 
    virtual int deReference(int i, int n=1) const;
    
    /// Alternative version of deReference which returns the
    /// index of particle i in this SimSnap is in the ancestor
    /// pointed to by pRel.
    virtual int deReference(int i, SimSnap *pRel) const;

    // Here are some functions which have a default implementation
    // which should work fine, but you might like to override them
    // anyway if you have a custom implementation which works
    // better. (e.g. getNearestNeighbours->grid implementation in grid.cpp)

    /// Get the index of the particle nearest the particle with index from
    virtual int getNearestNeighbour(int from) const;
    
    /// Get the index of the particle nearest to the particle from
    virtual int getNearestNeighbour(const Particle &from) const;

    /// Get a specified number of nearest neighbours according to the
    /// user-specified metric (which defaults to standard Euclidean 3D)
    ///
    /// @param from - the particle from which the nearest neighbours need to be found
    /// @param number - the number of particles to get
    /// @param metric - the metric object to use
    /// @returns list of pairs: [ptcl index, distance] in descending order of distance
    virtual std::list<std::pair<int,double> > getNearestNeighbours(const Particle &from, int number, const Metric &metric = Metric()) const;

    // SPH-like procedures

    /// Allocate extra memory to store SPH smoothing lengths, and optionally
    /// specify the number of particles and type of kernel to use for all
    /// subsequent SPH operations.
    virtual void initialiseSPH(int nSPH=40, SPHKernel *kern=NULL) const;

    /// Estimate the density at a given point using SPH
    virtual float estimateDensity(float x, float y, float z) const;
 
    /// Estimate density of particle (index n). Optionally you can
    /// pass a list of the nearest neighbours, if this is already known,
    /// to save recalculation
    virtual float estimateDensity(int n, std::list< std::pair<int,double> > *nnList=NULL) const;
    
    /// Get the smoothing length of particle index n. Optionally you can
    /// pass a list of the nearest neighbours, if this is already known,
    /// to save recalculation. If @ref consistentSmooth() has been called,
    /// this returns the consistent smoothing length. Otherwise it
    /// returns an estiamte of the smoothing length.
    virtual float getSmooth(int n,std::list< std::pair<int,double> > *nnList=NULL) const;
   
    /// Find smoothing lengths for all particles in simulation using
    /// self-consistent method (such that integral of density = mass)
    virtual void consistentSmooth(float tol=0) const;

    /// Find smoothing lengths for particle index n in simulation
    /// using self-consistent method.
    virtual float consistentSmooth(int n, float tol=0) const;

    /// Return the characteristic scale length between particles aroudn
    /// point specified by x,y,z.
    virtual float localPPscale(float x, float y, float z) const;

    /// Use SPH to calculate a hair-line column density 
    /// at coordinates x,y in given untis.
    virtual double getSPHColDen(float x, float y, const Unit &inunits) const;

    /// quick sanity check output
    virtual void diagnostics();
  
#ifdef SIMAN_FITS

    template <typename hdutype> 
    void SPHColumnDensityImage(hdutype *fitsFloatHDU, float x1, float x2, int nx, float y1, float y2, int ny, double units) const
    {
      // Renders a column density image into fitsFloatHDU
      // Multiplies by units first, to keep numbers nice.
    
      // Note everything has internal DOUBLE precision, whilst the
      // rendered image is a FLOAT.
    
      float discard,z1,z2;
      getExactBoundaries(discard,discard,discard,discard,z1,z2);
  
      long nelements = nx * ny; // number of pixels to render
    
      valarray<float> imgData(nelements);  // the array which will store the image
    
      float dx = (x2-x1)/(float)nx;
      float dy = (y2-y1)/(float)ny;

      int totel = nx*ny;
      int curel = 0;

      for(int x=0;x<nx;x++) {
      
	for(int y=0;y<ny;y++) {
	
	  //  imgData[x+y*ny] = (float) getSPHColDen((float)x*dx+x1,(float)y*dy+y1,z1,z2,units);

	  curel++;
	
	  std::cerr << "SimSnap: SPH column density calculation " << (curel*100)/totel << "%..." << "\r";
	  std::cerr.flush();
	
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

#endif //SIMAN_FITS
  


    virtual float getTotalMass() const; ///< Get the total mass in this simulation, in units of \ref getMassUnits()
 
    virtual void centreOn(float cx, float cy, float cz); ///< Recentre this simulation on the given coordinates
    virtual void centreOn(std::vector<double> coords); ///< Recentre simulation on given coordinates (x,y,z)
    virtual void centreOnVel(float cvx, float cvy, float cvz); ///< Recentre the velocity of this simulation on the given coordinates
    virtual void centreOnVel(std::vector<double> coords); ///< Recentre simulation on given velocity coordinates (vx,vy,vz);

    /// Calculate the centre of mass of this simulation and return it in parameters passed by reference
    virtual void centreOfMass(float &cx, float &cy, float &cz) const;
    virtual std::vector<double> centreOfMass() const; ///< Get centre of mass as vector (x,y,z)

    /// Calculate the centre of momentum of this simulation and return it as a velocity in parameters passed by reference
    /// NOTE - this is by no means equivalent to the velocity of the centre of mass
    virtual void centreOfMassVel(float &cvx, float &cvy, float &cvz) const;
    virtual std::vector<double> centreOfMassVel() const;

    /// Calculate the angular momentum of this simulation. Should be centred and velocity centred first!
    SimanVec angMom() const;

    /// Calculate a rotation curve at given points of R taking into account particles within scale_height
    /// of plane defined by angular momentum and halo centering.
    std::vector<double> rotCurve(std::vector<double> R, double scale_height) const;

    /// Calculate the shrinking sphere centre of particles in this simulation, see e.g. Power (2002)
    /// \param[out] cx, cy, cz - returns centre 
    /// \param sizechange - the ratio by which to change the size of the shrinking sphere at each stage
    /// \param minParticles - the number of particles at which the refinement ends and the final centre is calculated
    /// \param initialRadius - intial radius of the sphere in units of \ref getDistanceUnits() - negative to enclose entire simulation
    virtual void shrinkSphereCentre(float &cx, float &cy, float &cz, float sizechange=0.8, int minParticles=20 , float initialRadius=-1.) const;
    virtual std::vector<double> shrinkSphereCentre(float sizechange=0.8, int minParticles=50, float initialRadius=-1.) const;

    /// Calculate shrink sphere centre of particles in momentum space (and return in velocity coordinates)
    virtual void shrinkSphereCentreVel(float &cvx, float &cvy, float &cvz, float sizechange=0.8, int minParticles=20, float initialRadius=-1.) const;
    virtual std::vector<double> shrinkSphereCentreVel(float sizechange=0.8, int minParticles=50, float initialRadius=-1.) const;
    
    /// Get the exact boundaries of the simulation
    virtual void getExactBoundaries(float &x1, float &x2, float &y1, float &y2, float &z1, float &z2) const;

    /// Get the virial radius of a halo (which must already be centred, e.g. use centreOn(shrinkSphereCentre9)).
    /// @param delta - the mean overdensity inside the virial radius, relative to the universe critical
    ///                density (*not* the mean density; so multiply by omega_M if this is your intention)
    virtual float getVirialRadius(float delta=178) const throw(ConvergenceFailure);

    /// Return the critical density of the Universe with snapshots' omega_M0 and omega_Lambda0 
    /// at the redshift of this snapshot, in units of getDensityUnits().
    float criticalDensity() const;

    /// Return the critical density of the simulated Universe at z=0
    float criticalDensity0() const;
    
    virtual float power(float lambda) const;

    virtual SimSnap* getParent() const;


    // The following are not designed for overriding


    void write(std::string filename, int filetype) const;

    static int determineType(std::string filename);
    static SimSnap* loadFile(std::string filename);

    static const int unknown = 1;

    static const int tipsy = 2;

    static const int gadget = 4;

    static const int script = 8;

    static const int native = 16; 

    static const int endianWrong = 256;

    float *tau; /// To be replaced by general "property" system later on
  protected:
    SimSnap* pAutoParent; ///< set this to not NULL in child classes to automatically pass on calls to (e.g.) getBoxsize, getRedshift etc

    mutable int nPartSPH; 
    mutable float *pSmooth;
    mutable SPHKernel *pKernel; ///< kernel information for SPH

    Unit denUnits;
    Unit lenUnits;
    Unit enUnits;
    Unit velUnits;
    Unit massUnits;

    /// call this from derived classes when the data might change
    void bumpVersion();

  private:
  
    long version;

  
  };

  /// used for sorting when constructing nearest neighbour lists
  class compareSecond : public std::binary_function<std::pair<int,double>,std::pair<int,double>,bool> {
  public:   
    bool operator()(std::pair<int,double> p ,std::pair<int,double> q) const;

  };

} // namespace siman

#endif
