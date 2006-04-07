//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CBaseSimSnap Class
//
// Takes a copy of (possibly highly nested!) CSimSnaps, but stores
// the data as actual CParticles which uses extra memory but may
// lead to a huge speed increase for certain operations.

#ifndef __REALIZE_H_INCLUDED

#define __REALIZE_H_INCLUDED


class CBaseSimSnap : public CSimSnap {

public:
  CBaseSimSnap(CParticle **pParticleArray, int nParticles, units::CUnit length = units::CUnit(), units::CUnit mass = units::CUnit() );
  CBaseSimSnap(CSimSnap *copySim);
  CBaseSimSnap(std::string filename);
  CBaseSimSnap();

  ~CBaseSimSnap();

  static CBaseSimSnap* makeSlab(unsigned int nParticles, float depth, float extent, float density, float temp, int gentype=genRandom, units::CUnit length = units::CUnit(), units::CUnit mass = units::CUnit() );

  virtual CParticle * getParticle(int id);
  virtual void releaseParticle(CParticle *);
  
  virtual bool Load(bool fakeload=false);
  virtual bool isLoaded();
  
  virtual int getNumParticles();

  virtual float getBoxSize();
  virtual float getHubble();
  virtual float getRedshift();
  virtual float getOmegaM0();
  virtual float getOmegaLambda0();

  virtual void setHubble(float hubble);
  virtual void setBoxSize(float boxsize);
  virtual void setRedshift(float z);

  virtual void convertUnits(  units::CUnit distance = units::CUnit(), 
			      units::CUnit mass = units::CUnit(), 
			      units::CUnit velocity = units::CUnit(), 
			      units::CUnit density = units::CUnit(), 
			      units::CUnit energy = units::CUnit());


  virtual int deReference(int i, int n=1);
  
  static const int genRandom = 1;
  static const int genGrid = 2;
  

  static void nativeWrite(CSimSnap *sim, std::string filename);

  
  ///\ingroup ExtraData
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

protected:

  void allocateMemory();
  CParticle **pParticles;
  unsigned int numParticles;
    
  float boxSize;
  float hubble;
  float redshift;
  float om_m0;
  float om_lam0;

      
  /// \ingroup ExtraData
  ///@{
  map<string, float* > extraMap;
  map<string, string> extraMapName;
  map<string, units::CUnit> extraMapUnits;
  ///@}


};

#endif // _REALIZE_H_INCLUDED
