// glutinterfaceautomated.cpp - part of SimAn Simulation Analysis Library
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

#ifdef SIMAN_VIS

namespace siman {

  GlutInterfaceAutomated::GlutInterfaceAutomated(SimSnap *pSimI, unsigned int flags, Scripted *pCLI_i) :
    GlutInterface(pSimI, flags, pCLI_i) { cerr << "Init" << endl; };
  
  void GlutInterfaceAutomated::fullReset() { 
    
    GlutInterface::fullReset();
    
    cerr << "fR" << endl;

    flags |= autoRotate + axes + promo + writeTiff;
    flags &= ~dispColourBar; 
    auto_rx=0.004;
    auto_ry=0.01;
    pointScale=2.;
    setTrans(getTrans()/4.);

    numPartTarget=800000;
    numPartUpdatedTarget=numPartTarget*2;
    
  }

  bool GlutInterfaceAutomated::tick() {
    //  frame++;
    auto_rx=0.007*sin((float)frame/500.);
    auto_ry=0.007*cos((float)frame/323.);
    scale_aim*=1.005;
   
    //    if(frame>500 && frame<1000) scale_aim*=1.02;

    if(frame>1500) exit(0);
    return GlutInterface::tick();
  }

  
  
}

#endif
