// annotate.cpp - part of SimAn Simulation Analysis Library
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









#ifdef SIMAN_VIS

#include "../base.hpp"
#include <typeinfo>
#include <boost/lexical_cast.hpp>

#define PI 3.141592


namespace siman {

  Annotate::Annotate() {

  }

  void Annotate::plot(float ls) {

  }

  Annotate::~Annotate() {
    
  }



  AnnotateVector::AnnotateVector(const SimanVec &start, const SimanVec &len) : start(start), len(len) {
    
  }

  AnnotateVector::~AnnotateVector() {

  }

  void AnnotateVector::plot(float ls) {
    glColor4f(1.0,1.0,1.0,1.0);
    glBegin(GL_LINES);
    glVertex3f(start[0],start[1],start[2]);
    glVertex3f(start[0]+len[0],start[1]+len[1],start[2]+len[2]);
    glEnd();
  }




  AnnotatePlane::AnnotatePlane(const SimanVec &start, const SimanVec &len) : start(start), len(len) {
    v1 = SimanVec(1,0,0).proj(len);
    if(v1.abs()<0.0001) {
      v1=SimanVec(0,0,1).proj(len);
    }
    v2 = SimanVec(0,1,0).proj(len).proj(v1);
    if(v2.abs()<0.0001) {
      v2=SimanVec(0,0,1).proj(len).proj(v1);
    }

    
    //   cerr << v1[0] << " " << v1[1] << " " << v1[2] << " " << v1.abs() << endl;
    // cerr << v2[0] << " " << v2[1] << " " << v2[2] << " " << v2.abs() << endl;

    v1/=5.*v1.abs();
    v2/=5.*v2.abs();
    
	
  }

  AnnotatePlane::~AnnotatePlane() {

  }

  void AnnotatePlane::plot(float ls) {

    glColor4f(1.0,1.0,0.7,0.4);
    glBegin(GL_QUADS);
    glVertex3f(start[0]+(-v1[0]-v2[0])*ls,start[1]+(-v1[1]-v2[1])*ls,start[2]+(-v1[2]-v2[2])*ls);
    glVertex3f(start[0]+(-v1[0]+v2[0])*ls,start[1]+(-v1[1]+v2[1])*ls,start[2]+(-v1[2]+v2[2])*ls);
    glVertex3f(start[0]+(v1[0]+v2[0])*ls ,start[1]+( v1[1]+v2[1])*ls,start[2]+(v1[2]+v2[2])*ls);
    glVertex3f(start[0]+(v1[0]-v2[0])*ls ,start[1]+( v1[1]-v2[1])*ls,start[2]+(v1[2]-v2[2])*ls);
    glEnd();
   
  }


} // namespace siman

#endif // SIMAN_VIS
 
