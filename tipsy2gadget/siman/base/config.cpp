// config.cpp - part of SimAn Simulation Analysis Library
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
  Config::Config() {

#ifdef SIMAN_VIS
    pCurVis=NULL;
    glutInitialised=false;
#endif

    std::string filename = "config";
    
    if(!fileExists(filename)) {
      filename = std::string(getenv("SIMAN_DATA"))+std::string("/")+filename;
      
    }
    
    std::ifstream infile(filename.c_str());
    
    char line[2048];
    
    while(infile.getline(line,2048)) {
      if(line[0]!='#') {
	istringstream ss(line);
	string p,val;
	ss >> p >> val;
	(*this)[p]=val;
      }
    }

    
    try {
      setVerbose(boost::lexical_cast<int>(config["verbosity"]));
   
    } catch(boost::bad_lexical_cast &) {
   
    }


  }

  Config::~Config() {
    Scripted::deInitPython();
  }

}
