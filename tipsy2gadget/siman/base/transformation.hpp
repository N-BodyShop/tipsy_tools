// transformation.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __TRANSFORMATION_H_INCLUDED

#define __TRANSFORMATION_H_INCLUDED

namespace siman {
  class Transformation {
  public:
    Transformation() { };
    virtual ~Transformation() { } ;
    /// Perform the transformation - does not write back to SimSnap, but outputs data new Particle
    /// @param p - the particle data to update
    /// @param extra_in - any extra requested data from the array blocks (updatable)
   
    virtual void transform(Particle &p, std::vector<float> &extra ) const;

    /// Called to find out what (if any) extra data should be supplied in extra_in to @ref transform()
    /// Note - names may not be empty when this is called if the overarching Transformer object is dealing
    /// with a list of transformations. Thus, it is important that the Transformation child class takes note
    /// of the index within the vector of the data it will receive (which corresponds exactly with the index
    /// of the insertion of the name, so is accessible with size() before calling push_back()). See @ref MassTransformation
    /// implementation for an example of how this is done.
    ///
    /// @param[out] names - a list of the array names of extra data required for the transformation
    /// @param forSim - the simulation on which the transformation is about to be performed. Note
    /// that taking references/pointers to anything within forSim is not recommended, as this can
    /// all change between now and when "transform" gets called.
    virtual void extraBlocks(std::vector<std::string> & names, SimSnap & forSim) const;


    /// Because we want to be able to pass temporaries to Transformer (syntax like:
    ///
    ///   Transformer t(sim,Translation(1,1,5));
    ///
    ///  ) we need to be able to copy. Copy constructors are no use in this respect, since
    /// we'd end up with a base class, not the derived class. Thus we must be able to make a
    /// new copy of transformation objects, for which reason you must override copy() in any
    /// derived classes.
    virtual Transformation *copy() const;


  };

  /// Translation objects do what they say on the tin!
  class Translation : public Transformation {
  public:
    Translation(float x, float y, float z);
    virtual ~Translation() { } ;

    virtual void transform(Particle &p, std::vector<float> &extra ) const;
  
    virtual Transformation * copy() const;

  protected:
    float tx, ty, tz;
  };


  /// Wrap objects perform wrap-arounds so that the final value
  /// of all coordinates is clamped between min_coord and max_coord
  ///
  /// Useful after translating a periodic box, i.e.:
  ///
  ///  wrap(-boxlen/2,+boxlen/2) * translation(dx,dy,dz)
  class Wrap : public Transformation {
  public:
    Wrap(float min_coord, float max_coord);
    virtual ~Wrap() { } ;
    
    virtual void transform(Particle &p, std::vector<float> &extra) const;
    virtual Transformation * copy() const;

  protected:
    float min, max, delta;
  };


  /// A VelocityTranslation object when applied adds a constant 
  /// vector to the velocities of each particle

  class VelocityTranslation : public Transformation {
  public:
    VelocityTranslation(float vx, float vy, float vz);
    virtual ~VelocityTranslation() { } ;

    virtual void transform(Particle &p, std::vector<float> &extra ) const;
  
    virtual Transformation * copy() const;

  protected:
    float tx, ty, tz;
  };

  

  /// OrthoTrans encapsulates orthogonal matrix transformations, e.g. rotations
  class OrthoTrans : public Transformation {
  public:
    /// This constructor initialises the matrix from 9 floats. Note the constructor does not check explicitly
    /// whether the matrix is orthogonal and undefined results could occur if it is not. (Especially on inverseTransform).
    OrthoTrans(float mxx, float mxy, float mxz, float myx, float myy, float myz, float mzx, float mzy, float mzz);
    OrthoTrans(const std::vector<float> &matrix);
    /// Construct with identity matrix
    OrthoTrans(); 
    virtual ~OrthoTrans() { } ;

    virtual void transform(Particle &p,  std::vector<float> &extra ) const;
  
    virtual Transformation * copy() const;

  protected:
    float mxx,mxy,mxz,myx,myy,myz,mzx,mzy,mzz;
  };


  /// Type of orthogonal transform which rotates about the X axis
  class RotationX : public OrthoTrans {
  public:
    RotationX(float r);
  };


  /// Type of orthogonal transform which rotates about the Y axis
  class RotationY : public OrthoTrans {
  public:
    RotationY(float r);
  };

  /// Type of orthogonal transform which rotates about the Z axis
  class RotationZ : public OrthoTrans {
  public:
    RotationZ(float r);
  };


  class MassTransformation : public Transformation {
  public:
    /// constructs a transformation object which changes mass and rho attributes of each particle
    /// to reflect only the types of matter indicated in the flags. You can combine these flags
    /// to get mass of more than one species, e.g. MassTransformation(HI+HII+DM) sets mass of all particles
    /// so that dark matter particles are unchanged, gas particles have ionised and neutral H but no He,
    /// and star particles have zero mass.
    ///
    /// The output masses are written to the array specified by name write_[mass/rho]_as; these must already
    /// exist with suitable target units. Or if you want to write directly back to the particle
    /// property, leave the string blank ("")
    
    MassTransformation(int type, std::string write_mass_as="", std::string write_rho_as="");
    virtual void transform(Particle &p, std::vector<float> &extra ) const;

    virtual void extraBlocks(std::vector<std::string> & names, SimSnap &forsim) const;
    virtual Transformation * copy() const;


    
    const static int HI=1; ///< if set, include HI mass
    const static int HII=2; ///< if set, include HII mass
    const static int HeI=4; ///< if set, include HeI mass
    const static int HeII=8; ///< if set, include HeII mass
    const static int HeIII=16; ///< if set, include HeIII mass
    const static int dm=32; ///< if set, include DM mass
    const static int star=64; ///< if set, include star mass
  protected:
    int type;
    SimanArray *HIIarr;
    SimanArray *HeIIarr;
    SimanArray *HeIIIarr;
    std::string write_mass_as;
    std::string write_rho_as;

    /// the following refer to the index within the extra array of the required data
    mutable int nHII_ref, nHeII_ref, nHeIII_ref, write_mass_as_ref, write_rho_as_ref;

  };

  class LuminosityTransformation : public Transformation {
  public:
    LuminosityTransformation(std::string band);
    virtual void transform(Particle &p, std::vector<float> &extra ) const;

    virtual void extraBlocks(std::vector<std::string> & names, SimSnap &forsim) const;
    virtual Transformation * copy() const;
  protected:
    std::vector<std::vector<double> > tm_vals;
    std::vector<double> t_arr;
    std::vector<double> m_arr;
    mutable int l_ref;
    mutable int age_ref;
    mutable double t0, time_conv, lum_conv;
    std::string filter;
  };


  class TwoTransformation : public Transformation {
    
  public:
    /// a transformation object which applies t2, then t1
    TwoTransformation(const Transformation &t1, const Transformation &t2);
    TwoTransformation(const boost::shared_ptr<Transformation> & apT1, const boost::shared_ptr<Transformation> & apT2);
    virtual void transform(Particle &p, std::vector<float> &extra_in ) const;
  
    virtual void extraBlocks(std::vector<std::string> & names, SimSnap &forSim) const;
    virtual Transformation * copy() const;

    virtual ~TwoTransformation();
    
  protected:
    boost::shared_ptr<Transformation> apT1; ///< encapsulated pointer to component transformation. Using shared_ptr guarantees safety when copying etc.
    boost::shared_ptr<Transformation> apT2;
    //std::list<Transformation*> transList;
  };
  
  inline TwoTransformation operator*(const Transformation &t1, const Transformation &t2) {
    return TwoTransformation(t1,t2);
  }

}


#endif
