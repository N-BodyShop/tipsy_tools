// siman.cpp - part of SimAn Simulation Analysis Library
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

  siman::interactive_banner();
  std::cerr << "(RL)" << std::endl;
  Scripted scr;

  for(int n=1; n<argc; n++) {
    SimSnap *pF = SimSnap::loadFile(argv[n]);
    char vname = 'a'+(char)(n-1);
    std::string vname_s = " ";
    vname_s[0]=vname;
    scr.injectAndManage(vname_s,pF);

    std::cout << "Loaded " << argv[n] << " into var " << vname_s << std::endl;
  }


  scr.mainLoop();

}


