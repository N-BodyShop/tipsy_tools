// subset.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __SUBSET_H_INCLUDED
#define __SUBSET_H_INCLUDED


#include <vector>

namespace siman {

class Subset : public SimSnap {

 public:

  /// construct subset containing particles which satisfy filter
  Subset(const SimSnap &simulation, const Filter &filter);

  /// for compatibility you can pass a pointer to a SimSnap instead of a reference, if you like
  Subset(SimSnap *simulation, const Filter &filter);

  /// construct subset based on array of particle IDs
  Subset(SimSnap *simulation, int numPart, int *pParticleIDs);

  /// construct subset based on vector of particle IDs
  Subset(SimSnap *simulation, const std::vector<unsigned int> & particles);

  /// construct subset with no members, derived from given parent (pointer version)
  Subset(SimSnap *simulation);

  /// construct subset with no members, derived from given parent
  Subset(SimSnap &simulation);

  virtual ~Subset();


  // Standard functions which probably do not require recoding for
  // child classes, if they store the data in the provided "protected"
  // variables

  virtual unsigned int getNumParticles() const;
  virtual Particle* getParticle(unsigned int id);
  virtual const Particle* getConstParticle(unsigned int id) const;
  virtual void releaseParticle(unsigned int id);
 
  virtual void pushParticle(int n);
  virtual int deReference(int i, int n=1) const;
  virtual int deReference(int i, SimSnap *pObj) const;

  SimSnap* getParent() const;

  
  ///\ingroup ExtraData
  ///@{
  virtual SimanArray & createArray(std::string name, std::string fullname, Unit inunits = Unit());
  virtual void destroyArray(std::string name);
  virtual SimanArray & getArray(std::string name);
  virtual SimanArray & getArray(int n);
  virtual const SimanArray & getConstArray(std::string name) const;
  virtual const SimanArray & getConstArray(int n) const;
  // others are handled by auto-parenting of SimSnap
  ///@}


protected:

  Subset();

  /// returns the array corresponding to the parent's array
  /// - creates it if necessary, but reuses previously created
  /// instances where possible.
  SimanArray & getSubscriptedArrayFor(SimanArray &parnet) const;
  std::vector<unsigned int> particleRefs;

private:
  
  /// Holds details of referencer objects to be deleted when subset is destroyed
  mutable std::map<SimanArray*,SimanArraySubscripted*> createdReferencers;

};

}
#endif
