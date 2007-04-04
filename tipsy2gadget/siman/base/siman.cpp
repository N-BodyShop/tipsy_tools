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









#include "siman.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>


namespace siman {

  Config config;

  int __verbosity_level=1;

  void setVerbose(int level) {
    __verbosity_level = level;
  }

  int getVerbose() {
    return __verbosity_level;
  }

  /// A quick function to determine if a file exists
  bool fileExists(string filename) {
    ifstream fs(filename.c_str());
    return(fs.is_open());
  }

  /// Given a set of pairs (a,b) and a value of a, this function
  /// performs linear interpolation between the nearest neighbours
  /// of a. The set is passed as a vector which does not need to be
  /// sorted (hence the routine may be a little slow for repeated
  /// use).
  ///
  /// @param vector - represents set of pairs (a,b)
  /// @param val_first - value of a 
  /// @returns interpolated value of b
  ///
  /// @todo Implement alternative version of this routine which takes
  /// a sorted array and is therefore much more efficient
  ///
  /// @todo Templatise - no need for a, b to be doubles
  double interpolate(const vector<pair<double, double> > & source, double val_first) {
    vector<pair<double, double> >::const_iterator i = source.begin();
    
    pair<double, double> v, close_low, close_high;
    
    close_low.first=std::numeric_limits<int>::min();
    close_high.first=std::numeric_limits<int>::max();
    
    for(i=source.begin();i!=source.end();i++) {
      v=*i;
      if(v.first>close_low.first && v.first<val_first) {
	close_low=v;
      }
      if(v.first>val_first && v.first<close_high.first) {
	close_high=v;
      }
    }
    
    double distance = (val_first-close_low.first)/(close_high.first-close_low.first);
    
    return close_low.second*(1.-distance) + close_high.second*(distance);
  }

  /// Given a vector of values (vals) assorted in ascending order, and a target value v,
  /// this function returns the bracketing indexing in the vector for that value. It also
  /// returns prophi, which gives the fractional closeness to the higher value.
  void bracketval(const vector<double> &vals, double v, int &hi, int &lo, double &prophi) {
    hi=vals.size()-1;
    lo=0;
    while(abs(hi-lo)>1) {
      int md = (hi+lo)/2;
      if(vals[md]>v) hi=md; else lo=md;
    }
    if(vals[hi]!=vals[lo])
      prophi = (v-vals[lo])/(vals[hi]-vals[lo]);
    else
      prophi = 0.5;
  
  }

  double interpolate2d(const vector<double> & xvals, const vector<double> & yvals, double xv, double yv, const vector<vector<double> > &fvals) {
    int xhi, xlo, yhi, ylo;
    double xphi, yphi;
    bracketval(xvals,xv,xhi,xlo,xphi);
    bracketval(yvals,yv,yhi,ylo,yphi);
    return (fvals[xhi][yhi]*xphi+fvals[xlo][yhi]*(1.-xphi))*yphi + (fvals[xhi][ylo]*xphi+fvals[xlo][ylo]*(1.-xphi))*(1.-yphi);
  }
  
  string version_string() {
    return "SimAn Preview Release 0.1";
  }

  string withPathOf(string original, string filename) {
    boost::filesystem::path path(original);
    return (path.branch_path()/filename).string();
  }

  /// Prints a copyright message- as required by GPL
  void interactive_banner() {
    cerr << version_string() << endl;
    cerr << " Copyright (c) Andrew Pontzen 2005-7" << endl;
    cerr << " SimAn comes with ABSOLUTELY NO WARRANTY; for details" << endl;
    cerr << " see the LICENCE file or source code headers. This is free " << endl;
    cerr << " software; you are welcome to redistribute it under certain" << endl;
    cerr << " conditions; please refer to the same location for details." << endl;
    
  }


  double lumToMag(double lum, string band) {
    double solar = 5;
    try {
      solar = boost::lexical_cast<double>(siman::config["solarMag_"+band]);
    } catch(boost::bad_lexical_cast &) {
      if(getVerbose()>0)
	cerr << "lumToMag: No solar data for band, using solar abs. mag. = 5" << endl;
    }
    return solar-2.5*log10(lum);
  }

  double magToLum(double mag, string band) {
      double solar = 5;
    try {
      solar = boost::lexical_cast<double>(siman::config["solarMag_"+band]);
    } catch(boost::bad_lexical_cast &) {
      if(getVerbose()>0)
	cerr << "magToLum: No solar data for band, using solar abs. mag. = 5" << endl;
    }
    return std::pow(10.,(solar-mag)/2.5);
  }


  void tokenize(const string& str,
		vector<string>& tokens,
		const string& delimiters, 
		bool strip_whitespace)
  {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
      {
        // Found a token, add it to the vector.
	string consider = str.substr(lastPos, pos - lastPos);
	if(strip_whitespace) boost::algorithm::trim(consider);
        if(consider.size()>0)
	  tokens.push_back(consider);
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
      }
  }

/* These are apparently not defined in MSVC */
#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif
#if !defined(M_LN2)
#define M_LN2 0.69314718055994530942
#endif



#ifdef SIMAN_HUMLICECK_VOIGT

static double humlicek_v12(double x, double y)
    /** Approximation of Voigt profile by Humlicek's 12-point formula.
     *
     * J. Humlicek, J. Quant. Spectrosc. Radiat. Transfer, 21(1978), 309.
     *
     * Voigt-Profil:
     * V(x, y) = 2/pi^(1.5) * y^2/FWHM_L * \int[-inf,+inf](e^(-y^2)/(x+y)^2+...)
     */
{
    static const double T_v12[6] = {
	0.314240376254359,    0.947788391240164,    1.597682635152605,
	2.27950708050106,     3.020637025120890,    3.889724897869782
    };
    static const double alpha_v12[6] = {
       -1.393236997981977,   -0.231152406188676,    0.155351465642094,
       -6.21836623696556e-3, -9.190829861057113e-5, 6.275259577497896e-7
    };
    static const double beta_v12[6] = {
	1.011728045548831,   -0.751971469674635,    1.255772699323164e-2,
	1.0022008145159e-2,  -2.42068134815573e-4,  5.008480613664573e-7
    };
    static const double  y0_v12 = 1.50;
    double yp, xp, xm, sum, yp2, xp2, xm2;
    int k;

    sum = 0.;
    yp = y + y0_v12;
    yp2 = yp * yp;
    if((y > 0.85) || (fabs(x) < (18.1 * y + 1.65))) {
	/* Bereich I */
	for(k=0; k<6; k++) {
	    xp = x + T_v12[k];
	    xm = x - T_v12[k];
	    sum += ((alpha_v12[k] * xm + beta_v12[k] * yp) / (xm * xm + yp2)
                    + (beta_v12[k] * yp - alpha_v12[k] * xp) / (xp * xp + yp2));
	}
    } else {
	/* Bereich II */
	for(k=0; k<6; k++) {
	    xp = x + T_v12[k];
	    xp2 = xp * xp;
	    xm = x - T_v12[k];
	    xm2 = xm * xm;
	    sum += (((beta_v12[k] * (xm2 - y0_v12 * yp) - alpha_v12[k] * xm * (yp + y0_v12))
                     / ((xm2 + yp2) * (xm2 + y0_v12 * y0_v12)))
                    + ((beta_v12[k] * (xp2 - y0_v12 * yp) + alpha_v12[k] * xp * (yp + y0_v12))
                       / ((xp2 + yp2) * (xp2 + y0_v12 * y0_v12))));
	}
	if(fabs(x) < 100.)
	    sum = y * sum + exp(-x*x);
	else
	    sum *= y;
    }
    return sum;
}



  vector<double> voigt(int nbins, double min_lambda, double max_lambda, double lambda_centre, double gaussian_width, double lorentz_width)
  {
    /* Transform into reduced coordinates and call Humlicek's 12 point
     * formula:
     *     x = 2 \sqrt{\ln2} \frac{\nu-\nu_0}{\Delta\nu_G}
     *     y = \sqrt{\ln2} \frac{\Delta\nu_L}{\Delta\nu_G}
     */
    vector<double> res;
    double delta_lambda = (max_lambda-min_lambda)/nbins;
    double yh = sqrt(M_LN2) * lorentz_width / gaussian_width;
    for(int i=0; i<nbins; i++) {
      double xh = 1.665109 * (i*delta_lambda +min_lambda-lambda_centre) / gaussian_width;
      res.push_back(0.9394372/gaussian_width * humlicek_v12(xh, yh));
    }

    return res;
  }
  

#else // Harris Voigt

  

  vector<double> voigt(int nbins, double min_lambda, double max_lambda, double lambda_centre, double gaussian_width, double lorentz_width) {
    
    // data file stolen from vpguess
#include "h1h3.h"
    
    vector<double> result;
    
    double delta_lambda = (max_lambda-min_lambda)/nbins;
    double a = lorentz_width/gaussian_width;

    double factor = 0.564191/gaussian_width; // factor to ensure area stays normalized 

    for(int i=0; i<nbins; i++) {
      double lambda = delta_lambda*i+min_lambda;
      double u = (lambda-lambda_centre)/gaussian_width;
     
      if (u < 0.0)
	u = -u;
      
      double usq = u * u;
      
      if (u > 19.99) {
	/* Use asymptotic approximation */
	result.push_back(factor * (a / (1.772453851 * usq) * 
	  (1.0 + 1.5 / usq + 3.75 / (usq*usq)
	   - a*a / usq * (1.0 + 5.0 / usq + 26.25 / (usq*usq)))));
      } else {
	
	const double DU = 0.01;
	
	double emusq = exp(-usq);
	double ud = u / DU;
	int n1 = (int)ud;
	int n2 = n1 + 1;
	double dif1 = (double)n2 - ud;
	double dif2 = ud - (double)n1;
	
	result.push_back (factor * (emusq + a * (h1[n1] * dif1 + h1[n2] * dif2
				       + a * ((1.0 - 2.0 * usq) * emusq 
					      + a * (h3[n1] * dif1 + h3[n2] * dif2
						     + a * (0.5 - 2.0 * usq + 2.0/3.0 * usq*usq)
						     * emusq)))));
      }
    }
    return result;
  }

#endif

  vector<double> absorption_line(int nbins, double min_lambda, double max_lambda, double lambda_centre, double gamma, double v, double fN) {
    double l_width = gamma * (lambda_centre/1.e5) * (lambda_centre/1.e5) / (4 * 3.141592 * 3.e8);
    vector<double> res = voigt(nbins,min_lambda,max_lambda,lambda_centre,lambda_centre*v/3.e5,l_width);
    for(int i=0; i<nbins; i++) {
      res[i]*=fN*double(8.87e-21)*lambda_centre*lambda_centre;
    }
    return res; 
  }


} // namespace siman
