// columnlist.hpp - part of SimAn Simulation Analysis Library
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








#ifndef __COLUMNLIST_H_INCLUDED

#define __COLUMNLIST_H_INCLUDED

namespace siman {

class ColumnList {
  friend class Grid;

public:
  ColumnList(SimSnap &parentSim, float x1, float x2, int nxi, bool assign=true);

  float power(float wavenumber);

  virtual ~ColumnList();

  virtual void realize();

  virtual SimSnap & operator[](int index);

private:
  int nx;
  float dx;
  SimSnap **pColumnContents;
};

}

#endif // __COLUMNLIST_H_INCLUDED
