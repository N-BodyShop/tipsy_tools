// union.hpp - part of SimAn Simulation Analysis Library
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







#ifndef __UNION_H_INCLUDED

#define __UNION_H_INCLUDED


namespace siman {

class Union : public Subset {

 public:
  Union(SimSnap *pSim);

  

  virtual ~Union();

  virtual unsigned int getNumParticles() const;
  virtual Particle* getParticle(unsigned int id);
  virtual const Particle* getConstParticle(unsigned int id) const;
  virtual std::pair<SimSnap *,unsigned int> getSnapAndRef(unsigned int id) const;
 
  void add(SimSnap *pSub);
  void add(SimSnap &sub);

  virtual int deReference(int i, int n=1) const;
  virtual int deReference(int i, SimSnap *pObj) const;


  std::list<SimSnap*> & getList();

protected:

  std::list<SimSnap*> children;

};

}

#endif // UNION_H_INCLUDED
