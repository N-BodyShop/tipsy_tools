// transformation.cpp - part of SimAn Simulation Analysis Library
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









#include "../base.hpp"

namespace siman {

  void Transformation::transform(Particle &p,  vector<float> &extra ) const 
  {

  }


  void Transformation::extraBlocks(vector<string> & names, SimSnap & forsim) const {

  }

  Transformation * Transformation::copy() const {
    return new Transformation();
  }

  Translation::Translation(float x, float y, float z) : tx(x), ty(y), tz(z) {
  }


  void Translation::transform(Particle &p, vector<float> &extra ) const {
    p.x+=tx;
    p.y+=ty;
    p.z+=tz;
  }

  
  Transformation * Translation::copy() const {
    return new Translation(tx,ty,tz);
  }




  VelocityTranslation::VelocityTranslation(float x, float y, float z) : tx(x), ty(y), tz(z) {
  }


  void VelocityTranslation::transform(Particle &p, vector<float> &extra ) const {
    p.vx+=tx;
    p.vy+=ty;
    p.vz+=tz;
  }

  
  Transformation * VelocityTranslation::copy() const {
    return new VelocityTranslation(tx,ty,tz);
  }



  OrthoTrans::OrthoTrans(float mxx, float mxy, float mxz, float myx, float myy, float myz, float mzx, float mzy, float mzz) :
    mxx(mxx), mxy(mxy), mxz(mxz), myx(myx), myy(myy), myz(myz), mzx(mzx), mzy(mzy), mzz(mzz) {

  }

  OrthoTrans::OrthoTrans(const vector<float> &matrix) : 
    mxx(matrix[0]), mxy(matrix[3]), mxz(matrix[6]), myx(matrix[1]), myy(matrix[4]), myz(matrix[7]), mzx(matrix[2]), mzy(matrix[5]), mzz(matrix[8]) {
      
  }

  OrthoTrans::OrthoTrans() : mxx(1), mxy(0), mxz(0), myx(0), myy(1), myz(0), mzx(0), mzy(0), mzz(1) {

  }

  void OrthoTrans::transform(Particle &p, vector<float> &extra ) const {
    Particle in = p;
    p.x = mxx * in.x + mxy * in.y + mxz * in.z;
    p.y = myx * in.x + myy * in.y + myz * in.z;
    p.z = mzx * in.x + mzy * in.y + mzz * in.z;
   
    p.vx = mxx * in.vx + mxy * in.vy + mxz * in.vz;
    p.vy = myx * in.vx + myy * in.vy + myz * in.vz;
    p.vz = mzx * in.vx + mzy * in.vy + mzz * in.vz;
   
  }

  
  Transformation * OrthoTrans::copy() const {
    return new OrthoTrans(mxx,mxy,mxz,myx,myy,myz,mzx,mzy,mzz);
  }


  RotationX::RotationX(float r) : OrthoTrans() {
    myy = cos(r);
    mzz = myy;
    myz = sin(r);
    mzy = -myz;
  }



  RotationY::RotationY(float r) : OrthoTrans() {
    mxx = cos(r);
    mzz = mxx;
    mzx = sin(r);
    mxz = -mzx;
  }


  RotationZ::RotationZ(float r) : OrthoTrans() {
    mxx = cos(r);
    myy = mxx;
    mxy = sin(r);
    myx = -mxy;
  }

  MassTransformation::MassTransformation(int type, string write_mass_as, string write_rho_as) : type(type), write_mass_as(write_mass_as), write_rho_as(write_rho_as) {
  }


  void MassTransformation::extraBlocks(vector<string> & names, SimSnap & forSim) const {
 
    if((type&HI)!=0 || (type&HII)!=0) {
      nHII_ref = names.size();
      names.push_back("nHII");
    }

    if((type&HeI)!=0 || (type&HeII)!=0) {
      nHeII_ref =names.size();
      names.push_back("nHeII");
    }

    if((type&HeI)!=0 || (type&HeIII)!=0) {
      nHeIII_ref =names.size();
      names.push_back("nHeIII");
    }

    if(write_mass_as!="") {
      write_mass_as_ref = names.size();
      names.push_back(write_mass_as);
    } else {
      write_mass_as_ref = -1;
    }


    if(write_rho_as!="") {
      write_rho_as_ref = names.size();
      names.push_back(write_rho_as);
    } else {
      write_rho_as_ref = -1;
    }
      
    
  }

  void MassTransformation::transform(Particle &p, vector<float> &extra) const {
    
    Particle in = p;

    float *pMassOut = &(p.mass);
    if(write_mass_as_ref!=-1)
      pMassOut = &(extra[write_mass_as_ref]);

    float *pRhoOut = &(p.rho);
    if(write_rho_as_ref!=-1)
      pRhoOut = &(extra[write_rho_as_ref]);

    
    if(in.type==Particle::dm) {
      if((type&dm)==0) 
	(*pMassOut)=0;
	
    }

    else if(in.type==Particle::star) {
      if((type&star)==0)
	(*pMassOut)=0;
    }

    else if(in.type==Particle::gas) {
      (*pMassOut)=0;
      (*pRhoOut)=0;

      if((type&HI)!=0) {
	(*pMassOut)+=in.mass*(1-constants::heliumY)*(1-extra[nHII_ref]);
	(*pRhoOut)+=in.rho*(1-constants::heliumY)*(1-extra[nHII_ref]);
      }
      if((type&HII)!=0) {
	(*pMassOut)+=in.mass*(1-constants::heliumY)*extra[nHII_ref];
	(*pRhoOut)+=in.rho*(1-constants::heliumY)*extra[nHII_ref];
      }


      float num_he_per_H = constants::heliumY/(4.*(1.-constants::heliumY));

      if((type&HeI)!=0) {
	(*pRhoOut)+=in.rho*(1-constants::heliumY)*(num_he_per_H-extra[nHeII_ref]-extra[nHeIII_ref]);
	(*pMassOut)+=in.mass*(1-constants::heliumY)*(num_he_per_H-extra[nHeII_ref]-extra[nHeIII_ref]);
      }
      if((type&HeII)!=0) {
	
	(*pRhoOut)+=in.rho*(1-constants::heliumY)*extra[nHeII_ref];
	(*pMassOut)+=in.mass*(1-constants::heliumY)*extra[nHeII_ref];
      }

      if((type&HeIII)!=0) {
	
	(*pRhoOut)+=in.rho*(1-constants::heliumY)*extra[nHeIII_ref];
	(*pMassOut)+=in.mass*(1-constants::heliumY)*extra[nHeIII_ref];
      }
      
    }
    
  }

  LuminosityTransformation::LuminosityTransformation(string filter) : filter(filter){
    
    string pathname = (string) getenv("SIMAN_DATA");
    pathname+="/LUM_"+filter;
    
    ifstream file_lumdata(pathname.c_str());

    if(file_lumdata.is_open()) {

      char p[1024];
      file_lumdata.getline(p,1024);
      if(p[0]!='#') throw(SimanException("Format error in "+pathname));
      vector<string> mlabels;
      tokenize(string(p),mlabels," ",true);
      int num_mt = mlabels.size()-1;
      for(int i=0;i<num_mt;i++) {
	m_arr.push_back(boost::lexical_cast<double>(mlabels.at(i+1)));
      }
      while(file_lumdata.good()) {
	double time, lum;
	vector<double> row;
	file_lumdata >> time;
	for(int i=0; i<num_mt; i++) {
	  file_lumdata >> lum;
	  lum = magToLum(lum,filter)/1.e6;
	  row.push_back(lum); 
	}
	t_arr.push_back(time);
	tm_vals.push_back(row);
      }
    } else throw(SimanException("Unable to find "+pathname+" for luminosity transform"));
    
  }

  void LuminosityTransformation::extraBlocks(vector<string> &names, SimSnap & forsim) const {
    l_ref = names.size();
    names.push_back("L_"+filter);
    forsim.createArray("L_"+filter,filter+" luminosity",Unit(""));
    forsim.createArray("sage","stellar age",Unit("yr"));
    age_ref = names.size();
    names.push_back("tform");
    names.push_back("sage");
    t0 = forsim.scaleToAge() * (forsim.getDistanceUnits()/forsim.getVelocityUnits()).convertTo(forsim.getConstArray("tform").getUnits(),&forsim);
    time_conv = forsim.getConstArray("tform").getUnits().convertTo(Unit("yr"),&forsim);
    lum_conv = forsim.getMassUnits().convertTo(Unit("Msol"),&forsim);
  }

  void LuminosityTransformation::transform(Particle &p, vector<float> &extra) const {
    if(p.type==Particle::star) {
      static int dcheck  = 0;
      double time = time_conv * (t0-extra[age_ref]);
      extra[age_ref+1] = time;
      extra[l_ref] = p.mass * lum_conv * ((float) siman::interpolate2d(t_arr,m_arr,time,p.metal,tm_vals));
      dcheck++;
    } else {
      extra[l_ref] = 0.;
    }
  }
  
  Transformation * LuminosityTransformation::copy() const {
    return new LuminosityTransformation(filter);
  }

  Transformation * MassTransformation::copy() const {
    return new MassTransformation(type);
  }

  
  TwoTransformation::TwoTransformation(const Transformation &t1, const Transformation &t2) {
    apT1=boost::shared_ptr<Transformation>(t1.copy());
    apT2=boost::shared_ptr<Transformation>(t2.copy());
  }


  
  TwoTransformation::TwoTransformation(const boost::shared_ptr<Transformation> &t1, const boost::shared_ptr<Transformation> &t2) : 
    apT1(t1), apT2(t2) {

  }

  void TwoTransformation::transform( Particle &p,  vector<float> &extra ) const {
    Particle temp;
    apT2->transform(p,extra);
    apT1->transform(p,extra);
  }

  void TwoTransformation::extraBlocks(vector<string> & names, SimSnap &forsim) const {
    apT2->extraBlocks(names,forsim);
    apT1->extraBlocks(names,forsim);
  }

  Transformation * TwoTransformation::copy() const {
    return new TwoTransformation(apT1,apT2);
  }
  

  TwoTransformation::~TwoTransformation() {
 
  }

  
  Wrap::Wrap(float min, float max) : min(min), max(max), delta(max-min) {

  }


  void Wrap::transform(Particle &p, vector<float> &extra ) const {
    while(p.x<min) p.x+=delta;
    while(p.x>max) p.x-=delta;

    while(p.y<min) p.y+=delta;
    while(p.y>max) p.y-=delta;

    while(p.z<min) p.z+=delta;
    while(p.z>max) p.z-=delta;

  }

  
  Transformation * Wrap::copy() const {
    return new Wrap(min,max);
  }

  

}
