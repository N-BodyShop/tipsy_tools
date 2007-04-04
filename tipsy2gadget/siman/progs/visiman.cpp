// visiman.cpp - part of SimAn Simulation Analysis Library
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


using namespace siman;

int main(int argc, char **argv)
{

#ifdef SIMAN_VIS
  siman::interactive_banner();

    
  if(argc!=2) {
    std::cerr << "Syntax: visiman sim_file" << std::endl;
    std::cerr << std::endl << "sim_file = path to halo file" << std::endl;
    exit(0);
  }


  SimSnap *pF = SimSnap::loadFile(argv[1]);

  GlutInterface vis(pF,GlutInterface::dispColourBar);
   
  float col1[4] = {0,0,1,1};
  float col2[4] = {1,0,0,1};
  float col3[4] = {1,1,0,1};
  float col4[4] = {0,1,1,1};
  
  ColourMapGradient gascm(&vis,col1, col2);
  gascm.setColourBy(ColourMapGradient::useTemp);
  gascm.setLogScale(true);
  gascm.setHideOutOfRange(false);
  std::vector<double> r;
  r.push_back(3);
  r.push_back(6);
  gascm.setRange(r);

  ColourMapGradient starcm(&vis,col3,col4);
  starcm.setColourBy("tform");

  ColourMap dmcm(&vis);
  dmcm.setReferenceColour(0,1,0,0.1);
  
  ColourMapByType cm(&vis,&dmcm,&gascm,&starcm);

  cm.autoRange();

  vis.setColourMap(cm);
  vis.run(); 


  delete pF;

#else
  std::cerr << "Sorry, SimAn was not compiled with OpenGL support" << std::endl;
#endif
}


