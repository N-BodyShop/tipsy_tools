// simanarray.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __SIMANARRAY_H_INCLUDED

#define __SIMANARRAY_H_INCLUDED

namespace siman {

  class SimSnap;
  class Subset;

  class SimanArray {
  public:
    friend SimanArray operator*(const Unit &u, const SimanArray &c);
    friend SimanArray operator/(const Unit &u, const SimanArray &c);

    virtual float& operator[](unsigned int i) ;
    virtual const float& operator[](unsigned int i) const ;
    virtual void unset(unsigned int n) ;
    virtual bool isSet(unsigned int n) ;

    virtual bool isAllSet() ;
 
    virtual ~SimanArray() ;
    virtual unsigned int size() const ;
 

    SimanArray(unsigned int leni, SimSnap *pOwner, const Unit& unit) ;


    /// Return the owner SimSnap of this array
    SimSnap * getOwner();

    virtual void read(std::istream &fin, int file_ver_num) ;
    void readTipsy(std::string filename);


    virtual void write(std::ostream &file) const;


    SimanArray(const SimanArray& copyfrom);
    void operator=(const SimanArray &copyfrom) ;

    /// Multiply this array by the contents of another array. Modfies the units of this
    /// array as appropriate. Arrays must be the same length, otherwise throws
    /// MismatchedArrayLength.
    void operator*=(const SimanArray &mul) throw(MismatchedArrayLength);


    /// Divide this array by the contents of another array. Modfies the units of this
    /// array as appropriate. Arrays must be the same length, otherwise throws
    /// MismatchedArrayLength.
    void operator/=(const SimanArray &div) throw(MismatchedArrayLength);


    /// Add to this array the contents of another array. 
    /// Units must be dimensionally compatible (else throws UnitsError).
    /// Array length must be the same (else throws MismatchedArrayLength).
    void operator+=(const SimanArray &add) throw(UnitsError, MismatchedArrayLength);


    /// Subtract from this array the contents of another array. 
    /// Units must be dimensionally compatible (else throws UnitsError).
    /// Array length must be the same (else throws MismatchedArrayLength).
    void operator-=(const SimanArray &sub) throw(UnitsError, MismatchedArrayLength);

    /// Multiply this array by a fixed value
    void operator*=(float mul);

    /// Divide this array by a fixed value
    void operator/=(float div);

    /// Add to this array a fixed value
    /// It is the caller's responsibility to ensure the value is in the 
    /// appropriate units, determined by calling getUnits()
    void operator+=(float add);

    /// Subtract from this array a fixed value
    /// It is the caller's responsibility to ensure the value is in the 
    /// appropriate units, determined by calling getUnits()
    void operator-=(float sub);


    /// Multiply this array with another, returning the result in a new
    /// SimanArray (with appropriate units). Lengths must match 
    /// (otherwise throws MismatchedArrayLength).
    SimanArray operator*(const SimanArray &mul2) const throw(MismatchedArrayLength);

  
    /// Multiply this array by another, returning the result in a new
    /// SimanArray (with appropriate units). Lengths must match 
    /// (otherwise throws MismatchedArrayLength).
    SimanArray operator/(const SimanArray &div2) const throw(MismatchedArrayLength);


    /// Add another array to this array, returning the result in a new
    /// SimanArray (with same units as this, i.e. RH side). Lengths must match 
    /// (otherwise throws MismatchedArrayLength). Units must match dimensionally
    /// (else throws UnitsError).
    SimanArray operator+(const SimanArray &add2) const throw(MismatchedArrayLength,UnitsError);


    /// Subtratct another array from this array, returning the result in a new
    /// SimanArray (with same units as this, i.e. LH side). Lengths must match 
    /// (otherwise throws MismatchedArrayLength). Units must match dimensionally
    /// (else throws UnitsError).
    SimanArray operator-(const SimanArray &sub2) const throw(MismatchedArrayLength,UnitsError);
 
    SimanArray operator*(float mul) const;
    SimanArray operator/(float div) const;
    SimanArray operator+(float add) const;
    SimanArray operator-(float sub) const;

    SimanArray operator*(const Unit& u) const;
    SimanArray operator/(const Unit& u) const;

    /// Take power of this array, returning the result in a new array.
    SimanArray power(const Rational &p) const;

    /// Return maximum element
    float max();
    
    /// Return minimum element
    float min();

    /// Return mean
    float mean();

    /// Return total
    float total();

    virtual Unit getUnits() const;
    virtual void convertUnits(const Unit & to) throw(UnitsError);

  protected:
    SimanArray() ;
    unsigned int len;
    Unit units;
    SimSnap *pOwner;
  private:

  
 
      
    float* array;
  };

  inline SimanArray operator*(float m, const SimanArray &c) {
    return c*m;
  }


  inline SimanArray operator/(float m, const SimanArray &c) {
    SimanArray n(c);
    for(unsigned int i=0; i<n.size(); i++) {
      n[i]=m/n[i];
    }
    return n;
  }


  inline SimanArray operator+(float m, const SimanArray &c) {
    return c+m;
  }


  inline SimanArray operator-(float m, const SimanArray &c) {
    SimanArray n(c);
    for(unsigned int i=0; i<n.size(); i++) {
      n[i]=m-n[i];
    }
    return n;
  }

  inline SimanArray operator*(const Unit &u, const SimanArray &c) {
    SimanArray n(c);
    n.units*=u;
    return n;
  }


  inline SimanArray operator/(const Unit &u, const SimanArray &c) {
    SimanArray n(c);
    n.units=u/n.units;
    return n;
  }

  class SimanArraySubscripted : public SimanArray {
  public:
    virtual float& operator[](unsigned int i) ;
    virtual const float& operator[](unsigned int i) const ;

    virtual ~SimanArraySubscripted() { }
  
    SimanArraySubscripted(SimanArray &origArray_in, Subset *pForSubset);
   

    virtual void read(std::istream &fin, int file_ver_num);
    virtual void write(std::ostream &file) const;
    
    virtual Unit getUnits() const;
    virtual void convertUnits(const Unit & to) throw(UnitsError);
    virtual unsigned int size() const;
  
    unsigned int len;
  
  private:
    SimanArray *origArray;
    const Subset *sub;
 
  };


  /// A SimanArrayVirtual holds a pointer to member float of an abstract Particle,
  /// then uses it to dereference against a SimSnap set of particles to produce
  /// the effect of having an "array" of individual particle properties with little
  /// CPU overhead. The easiest way to instantiate a SimanArrayVirtual is to use
  /// SimSnap::getArray("..") where .. is one of the reserved names:
  /// x, y, z, vx, vy, vz, mass, rho, temp, u, ne
  class SimanArrayVirtual : public SimanArray {

  public:
    virtual float& operator[](unsigned int i) ;
    virtual const float& operator[](unsigned int i) const ;

    virtual ~SimanArrayVirtual() { }
  
    /// Constructs a virtual array of particle members of pOwner, referring
    /// to the element "member" of each particle.
    SimanArrayVirtual(SimSnap *pOwner, float Particle::* member);
       
    virtual Unit getUnits() const;
    virtual void convertUnits(const Unit & to) throw(UnitsError);
    virtual unsigned int size() const;
  
  
  private:
    float Particle::* member;
  
 
  };


}
#endif
