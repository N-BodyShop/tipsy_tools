
#include "../base.hpp"
#include <boost/lexical_cast.hpp>

namespace siman {

  SimanVec::SimanVec() : vector<double>() {

  }

  SimanVec::SimanVec(const std::vector<double> &in) : vector<double>(in) {

  }

  SimanVec::SimanVec(double x, double y, double z) : vector<double>() {
    push_back(x);
    push_back(y);
    push_back(z);
  }

  double SimanVec::dot(const SimanVec & v2) const throw(MismatchedArrayLength) {
    if(v2.size()!=size())
      throw MismatchedArrayLength();
    double r = 0;
    for(unsigned int i=0; i<size(); i++) {
     
      r+=(*this)[i]*v2[i];
      
    }
    return r;
  }

  SimanVec SimanVec::proj(const SimanVec &v2) const throw(MismatchedArrayLength) {

    // Projects using tensor I - (v2) [x] (v2)
    double v2n = (v2).dot(v2);


    if(v2.size()!=size())
      throw MismatchedArrayLength();
    
    SimanVec re;
    
    for(unsigned int i=0; i<size(); i++) {
      double t=0;
      for(unsigned int j=0; j<size(); j++) {
	if(i==j) {
	 
	  t+=(*this)[j]*(1.-v2[j]*v2[i]/v2n);
	} else {
	 
	  t+=(*this)[j]*(-v2[j]*v2[i]/v2n);
	}
      }
     
      re.push_back(t);
    }
    
    return re;

  }

  double SimanVec::abs() const {
    return sqrt((*this).dot(*this));
  }

  void SimanVec::operator/=(double r) {
    for(unsigned int i=0; i<size(); i++) {
      (*this)[i]/=r;
    }
  }


  void SimanVec::operator*=(double r) {
    for(unsigned int i=0; i<size(); i++) {
      ((*this)[i])*=r;
    }
  }


}
