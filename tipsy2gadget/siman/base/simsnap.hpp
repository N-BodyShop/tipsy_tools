//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// SIMSNAP
//
// Simsnap holds a set of particles and associated properties.
//
// I am aiming to make it a random access iterator
// acting like a vector<CParticle>, ultimately.
//
// At the moment, this is not the case however.


#ifndef __SIMSNAP_H_INCLUDED

#define __SIMSNAP_H_INCLUDED

#include "siman.hpp"


class CSimSnap : public CSimanObject {

public:

  CSimSnap();
  virtual ~CSimSnap();

  virtual unsigned int supports();
  virtual string className();  ///< returns class name
  virtual bool references(CSimanObject *p);
  virtual CSimanObject * dispatch(string command, std::istream *stream, CScripted *script);
   


  virtual CParticle * getParticle(int id);
  virtual void releaseParticle(CParticle *); // by definition does nothing
  
  virtual int getNumParticles();

  virtual units::CUnit getDistanceUnits();
  virtual units::CUnit getVelocityUnits();
  virtual units::CUnit getDensityUnits();
  virtual units::CUnit getEnergyUnits();
  virtual units::CUnit getMassUnits();
  float getHeliumMassFrac();


  virtual void convertUnits(units::CUnit distance = units::CUnit(), 
			    units::CUnit mass = units::CUnit(),
			    units::CUnit velocity = units::CUnit(),
			    units::CUnit density = units::CUnit(),
			    units::CUnit energy = units::CUnit());

  virtual void UFromTemp(); ///< use Temp and NE to determine U for gas particles
  virtual void TempFromU(); ///< use U and NE to determine Temp for gas particles

  virtual float getBoxSize();
  virtual float getApparentBoxSize(); ///< calculates from maximal particle positions the apparent box size of this sim
  virtual float getHubble();
  virtual float getRedshift();
  virtual float getOmegaM0();
  virtual float getOmegaLambda0();

  virtual void setHubble(float h);
  virtual void setBoxSize(float b);
  virtual void setRedshift(float z);

  ///\defgroup ExtraData Extra data handlers:
  /// Methods for associating extra data with a SimSnap
  /// @todo Currently these do not propogate up and down a stack of SimSnaps, which is unfortunate.
  ///        Fixing this will be a bit involved, especially for classes like CSubset
  ///@{
  virtual float *createArray(string name, string fullname, units::CUnit inunits = units::CUnit());
  virtual void destroyArray(string name);
  virtual float *getArray(string name);
  virtual float *getArray(int n);
  virtual string getArrayLongName(string name);
  virtual string getArrayLongName(int n);
  virtual string getArrayName(int n);
  virtual units::CUnit getArrayUnits(string name);
  virtual units::CUnit getArrayUnits(int n);
  virtual int getNumArrays();
  ///@}

  

  // Problematic: this is only implemented in a select few
  // child classes, so should only be called when you know
  // it's fine! But casts (CSubset*) bla-> thing() could be
  // equally problematic.

  virtual void pushParticle(int n);

  // Another slightly odd one: deReference returns the particle
  // which i refers to in the nth parent, but may well
  // crash and burn for obvious reasons!

  virtual int deReference(int i, int n=1);
  
  // Here are some functions which have a default implementation
  // which should work fine, but you might like to override them
  // anyway if you have a custom implementation which works
  // better. (e.g. getNearestNeighbours->grid implementation in grid.cpp)

  virtual int getNearestNeighbour(int from);
  virtual int getNearestNeighbour(const CParticle &from);
  virtual list<pair<int,double> > getNearestNeighbours(const CParticle &from, int number, const CMetric &metric = CMetric());

  // SPH-like procedures

  virtual void initialiseSPH(int nSPH=40, CSPHKernel *kern=NULL);
  virtual float estimateDensity(float x, float y, float z);
  virtual float estimateDensityH0(float x, float y, float z, list< pair<int,double > > *nnList=NULL);
  virtual float estimateDensity(int n, list< pair<int,double> > *nnList=NULL);
  virtual float getSmooth(int n,list< pair<int,double> > *nnList=NULL);
  virtual void consistentSmooth(float tol=0);
  virtual float consistentSmooth(int n, float tol=0);

  virtual float localPPscale(float x, float y, float z);
  virtual float getSPHColDen(float x, float y, float z1, float z2,double units, bool neutral);
 
  
  #ifdef SIMAN_FITS

  template <typename hdutype> 
  void SPHColumnDensityImage(hdutype *fitsFloatHDU, float x1, float x2, int nx, float y1, float y2, int ny, double units, bool neutral) 
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
	
	imgData[x+y*ny] = (float) getSPHColDen((float)x*dx+x1,(float)y*dy+y1,z1,z2,units,neutral);

	curel++;
	
	cout << "CSimSnap: SPH column density calculation " << (curel*100)/totel << "%..." << "\r";
	cout.flush();
	
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
  
  // Simple procedures

  virtual float getTotalMass();
  virtual float getH0Mass();

  virtual void centreOfMass(float &cx, float &cy, float &cz);
  virtual void centreOfMassVel(float &cvx, float &cvy, float &cvz);

  virtual void shrinkSphereCentre(float &cx, float &cy, float &cz, float sizechange=0.5, int minParticles=20 , float initialRadius=-1.);

  virtual void getExactBoundaries(float &x1, float &x2, float &y1, float &y2, float &z1, float &z2);

  virtual float power(float lambda);

  virtual CSimSnap* getParent();


  // The following are not designed for overriding


  void write(std::string filename, int filetype);

  static int determineType(std::string filename);
  static CSimSnap* loadFile(std::string filename);

  static const int unknown = 1;

  static const int tipsy = 2;

  static const int gadget = 4;

  static const int script = 8;

  static const int native = 16; 

  static const int endianWrong = 256;

protected:
  CSimSnap* pAutoParent; // set this to not NULL in child classes to automatically pass on calls to (e.g.) getBoxsize, getRedshift etc



  int nPartSPH; 
  float *pSmooth;
  CSPHKernel *pKernel; // kernel information for SPH

  units::CUnit denUnits;
  units::CUnit lenUnits;
  units::CUnit enUnits;
  units::CUnit velUnits;
  units::CUnit massUnits;

private:

  
};

// for nearest neighbour searches

class compareSecond : public binary_function<pair<int,double>,pair<int,double>,bool> {
public:
  // used for sorting when constructing nearest neighbour lists
  bool operator()(pair<int,double> p ,pair<int,double> q) const;

};



#endif
