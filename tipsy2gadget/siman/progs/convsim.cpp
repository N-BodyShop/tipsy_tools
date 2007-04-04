// convsim.cpp - part of SimAn Simulation Analysis Library
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


#include <siman/base.hpp>

using namespace siman;

int main(int argc, char **argv)
{

    
  std::string path;
 
  unsigned int type=0;
  
  if(argc<3) {
    std::cerr << "Syntax: convsim [type] [file_list]" << std::endl;
    exit(0);
  }

  if(strcmp(argv[1],"gadget")==0)
    type = SimSnap::gadget;

  if(strcmp(argv[1],"tipsy")==0)
    type = SimSnap::tipsy;
  
  if(strcmp(argv[1],"siman")==0)
    type = SimSnap::native;

  if(type==0) {
    std::cerr << "No match for type " << argv[1] << std::endl;
    exit(0);
  }


  
  for(int n=2;n<argc;n++) {
    path = argv[n];
    
    SimSnap *pF = SimSnap::loadFile(path);

    // see if convsim.units exists; if so, read in the units and convert.
    // If an error occurs, ignore file and continue writing without conversion
    if(fileExists("convsim.units")) {
      try {
	Unit lenUnits("(kpc h^-1 a)");
	Unit massUnits("(1.e10 msol h^-1)");
	Unit velUnits("(km s^-1 a^1/2)");
	
	Unit denUnits = massUnits/(lenUnits*lenUnits*lenUnits);
	
	Unit enUnits("(km^2 s^-2)");
	
	std::ifstream file_units("convsim.units");
	
	if(file_units.is_open()) {
	  if(getVerbose()>1)
	    std::cerr << "Loading output units from convsim.units" << std::endl;
	  file_units >> lenUnits;
	  file_units >> massUnits;
	  file_units >> velUnits;
	  file_units >> denUnits;
	  file_units >> enUnits;
	  pF->convertUnits(lenUnits,massUnits,velUnits,denUnits,enUnits);
	}
      } catch (UnitsError &e) {
	std::cerr << "Unable to convert units" << std::endl << e;
      }
    }

    pF->write(path+"."+argv[1],type);
    delete pF;
  }
}


