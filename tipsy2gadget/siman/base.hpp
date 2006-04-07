//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __SIMLIB_H_INCLUDED

#define __SIMLIB_H_INCLUDED

// STL & other standard library requirements

#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <utility>
#include <functional>
#include <list>
#include <valarray>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <complex>
#include <queue>
#include <map>
#include <set>

using namespace std;
// CCFits
#ifdef SIMAN_FITS
#include <CCfits/CCfits>
using namespace CCfits;
#endif

// OMP
#ifdef SIMAN_OMP
#include <omp.h>
#endif

// siman
#include "base/rational.hpp"
#include "base/units.hpp"
#include "base/simanexception.hpp"
#include "base/simanobject.hpp"
#include "base/value.hpp"
#include "base/coordinate.hpp"
#include "base/particle.hpp"
#include "base/metric.hpp"
#include "base/sphkernel.hpp"
#include "base/simsnap.hpp"
#include "base/basesimsnap.hpp"
#include "base/filter.hpp"
#include "base/subset.hpp"
#include "base/gadgetfile.hpp"
#include "base/tipsyfile.hpp"
#include "base/scripted.hpp"
#include "base/columnlist.hpp"
#include "base/columngrid.hpp"
#include "base/grid.hpp"
#include "base/units.hpp"
#include "base/geometry.hpp"
#include "base/sphkernel.hpp"

namespace siman {
  extern CScripted config;
  bool fileExists(string filename);
}

#endif
