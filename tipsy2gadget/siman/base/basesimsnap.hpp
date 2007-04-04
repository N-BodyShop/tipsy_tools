// basesimsnap.hpp - part of SimAn Simulation Analysis Library
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








#ifndef __BASESIMSNAP_H_INCLUDED

#define __BASESIMSNAP_H_INCLUDED

namespace siman {

class BaseSimSnap : public SimSnap {

public:
  BaseSimSnap(Particle **pParticleArray, int nParticles, Unit length = Unit(), Unit mass = Unit() );
  BaseSimSnap(const SimSnap *copySim);
  BaseSimSnap(std::string filename);
  BaseSimSnap();

  virtual ~BaseSimSnap();

  static BaseSimSnap* makeSlab(unsigned int nParticles, float depth, float extent, float density, float temp, int gentype=genRandom, Unit length = Unit(), Unit mass = Unit() );

  virtual Particle * getParticle(unsigned int id);
  virtual const Particle * getConstParticle(unsigned int id) const;
  virtual void releaseParticle(unsigned int id);

  virtual void setupReservedArrays();
  
  virtual bool Load(bool fakeload=false);
  
  virtual unsigned int getNumParticles() const;

  virtual float getBoxSize() const;
  virtual float getHubble() const;
  virtual float getRedshift() const;
  virtual float getOmegaM0() const;
  virtual float getOmegaLambda0() const;

  virtual void setHubble(float hubble);
  virtual void setBoxSize(float boxsize);
  virtual void setRedshift(float z);
  virtual void setOmegaM0(float om0);
  virtual void setOmegaLambda0(float ol0);

  virtual void convertUnits(  Unit distance = Unit(), 
			      Unit mass = Unit(), 
			      Unit velocity = Unit(), 
			      Unit density = Unit(), 
			      Unit energy = Unit());


  virtual int deReference(int i, int n=1) const;
  virtual int deReference(int i, SimSnap * pObj) const;
  static const int genRandom = 1;
  static const int genGrid = 2;
  

  static void nativeWrite(const SimSnap *sim, std::string filename);

  
  ///\ingroup ExtraData
  ///@{
  virtual SimanArray & createArray(std::string name, std::string fullname, Unit inunits = Unit());
  virtual void destroyArray(std::string name);
  virtual SimanArray & getArray(std::string name);
  virtual SimanArray & getArray(int n);
  virtual const SimanArray & getConstArray(std::string name) const;
  virtual const SimanArray & getConstArray(int n) const;
  virtual std::string getArrayLongName(std::string name) const;
  virtual std::string getArrayLongName(int n) const;
  virtual std::string getArrayName(int n) const;
  virtual Unit getArrayUnits(std::string name) const;
  virtual Unit getArrayUnits(int n) const;
  virtual int getNumArrays() const;
  virtual int getArrayIndex(std::string n) const;
  ///@}

protected:

  void allocateMemory();
  Particle **pParticles;
  unsigned int numParticles;
    
  float boxSize;
  float hubble;
  float redshift;
  float om_m0;
  float om_lam0;

      
  /// \ingroup ExtraData
  ///@{
  std::map<std::string, SimanArray* > extraMap;
  std::map<std::string, std::string> extraMapName;
  std::list<std::map<std::string, SimanArray*>::iterator> extraList;
  int extraMap_reserved_end;
  ///@}


};

} // namespace siman

#endif // __BASESIMSNAP_H_INCLUDED
