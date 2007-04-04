// simanobject.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __SIMANOBJECT_H_INCLUDED

#define __SIMANOBJECT_H_INCLUDED

namespace siman {

/// In earlier versions of SimAn, SimanObject performed more useful functions.
/// These are now superceded by better use of RTTI. SimanObject remains as a
/// way of passing generic objects around (particularly between python and C++).

class SimanObject {
public:
  virtual ~SimanObject() { } ;
protected:
  SimanObject() { } ;
};

}
#endif
