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
#include <limits>

#include <boost/shared_ptr.hpp>

// CCFits
#ifdef SIMAN_FITS
#include <CCfits/CCfits>
#endif

// OMP
#ifdef SIMAN_OMP
#include <omp.h>
#endif

// siman

#include "base/constants.hpp"
#include "base/rational.hpp"
#include "base/units.hpp"
#include "base/simanexception.hpp"
#include "base/simanobject.hpp"
#include "base/simanvec.hpp"
#include "base/particle.hpp"
#include "base/simanarray.hpp"
#include "base/coordinate.hpp"
#include "base/metric.hpp"
#include "base/sphkernel.hpp"
#include "base/transformation.hpp"
#include "base/simsnap.hpp"
#include "base/basesimsnap.hpp"
#include "base/filter.hpp"
#include "base/subset.hpp"
#include "base/gadgetsnap.hpp"
#include "base/tipsysnap.hpp"
#include "base/scripted.hpp"
#include "base/columnlist.hpp"
#include "base/columngrid.hpp"
#include "base/grid.hpp"
#include "base/sphkernel.hpp"
#include "base/config.hpp"


#ifdef SIMAN_VIS
#include "base/glew.h"
#include <GL/freeglut.h>
#include "base/visualiser.hpp"
#include "base/annotate.hpp"
#include "base/colourmap.hpp"
#include "base/glutinterface.hpp"
#include "base/glutinterfaceautomated.hpp"
#endif

#ifdef SIMAN_EXTRA
#include "extra.hpp"
#endif

namespace siman {
  //  extern Scripted config;
  
  /// determines whether a given file exists and is readable
  bool fileExists(std::string filename);

  /// extracts a path from extract_path and appends filename, e.g.
  ///
  ///   withPathOf("/foo/bar/pling","file") --> /foo/bar/file
  ///
  std::string withPathOf(std::string extract_path, std::string filename);
  
  double interpolate(const std::vector<std::pair<double, double> > & source, double val_first);
  double interpolate2d(const std::vector<double> & xvals, const std::vector<double> & yvals, double xv, double yv, const std::vector<std::vector<double> > &fvals);

  /// setVerbose sets the amount of output for SimAn classes to
  /// produce. When implementing classes, use getVerbose() to
  /// determine the level and produce appropriate output.
  ///
  /// @param level - 0 = no output; 1 = output on errors; 2 = verbose; 3 = debug only verbose
  void setVerbose(int level);

  /// Get the current verbosity level - @see setVerbose()
  int getVerbose();

  /// Returns a string describing the current version of the library
  std::string version_string();

  /// Prints a copyright message- as required by GPL
  void interactive_banner();

  /// The configuration settings can be accessed via the one and only Config object.
  extern Config config;
  
  double lumToMag(double luminosity,std::string band);
  double magToLum(double magnitude,std::string band);
  
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = "|",  bool strip_whitespace = true);

  /// Returns a voigt profile, normalized such that the integral is 1
  /// 
  /// @param nbins - number of bins in wavelength
  /// @param min_lambda - value of lambda in bin 0
  /// @param max_lambda - value of lambda in bin nbins-1
  /// @param lambda_centre - value of centre of voigt profile
  /// @param gaussian_width - FWHM of gaussian, in angstroms
  /// @param lorentz_width - characteristic width of lorentz, in angstroms ( phi_lor  = lorentz_width / (lambda^2 + lorentz_width^2) )
  std::vector<double> voigt(int nbins, double min_lambda, double max_lambda, double lambda_centre, double gaussian_width, double lorentz_width);

  /// Returns the optical depth of an absorption line as a function of lambda
  ///
  /// @param nbins - number of bins in wavelength
  /// @param min_lambda - value of lambda in bin 0
  /// @param max_lambda - value of lambda in bin nbins-1
  /// @param lambda_centre - value of centre of voigt profile
  /// @param gamma - transition damping gamma (s^-1)
  /// @param v - velocity width of line (km s^-1)
  /// @param fN - product of oscillator strength with column density (cm^-2)
  std::vector<double> absorption_line( int nbins, double min_lambda, double max_lambda, double lambda_centre, double gamma, double v, double fN);
}

#ifdef SIMAN_COMPILING_LIBRARY
using namespace std; // temporary measure

#endif
#endif
