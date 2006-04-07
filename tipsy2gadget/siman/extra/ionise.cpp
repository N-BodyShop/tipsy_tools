//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#include "../base.hpp"
#include "../extra.hpp"

using namespace std;

CIonise::CIonise(double j21i, double alphai, CSimSnap *pSimi) {

  j21 = j21i;
  alpha = alphai;
  
  nu_min = 13.6;
  nu_max = 100;

  nu_quant = 30;

  pJ=NULL;

  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(units::CUnit(units::den_protonsPerCm3));
  enToErgsPerG = pSim->getEnergyUnits().convertTo(units::CUnit(units::energy_erg)/units::CUnit(units::mass_g));

  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv = 0.001;

  flags = initJAndAlpha;
  
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
  
}


CIonise::CIonise(CSimSnap *pSimi, float z) {

  if(z<0)
    z=pSimi->getRedshift();


  string pathname = (string) getenv("SIMAN_DATA");
  pathname+="/UVNORM";
  
  ifstream file_uv_norm(pathname.c_str());

  
  j21 = -1;
  alpha = -1.6;

  if(file_uv_norm.good()) {

    float z_f=0, j21_f=0, z_xf=100, j21_xf=0;
    while(!file_uv_norm.eof()) {
      
      file_uv_norm >> z_f >> j21_f;
      if(z_f<z) {
	float interp = (z_xf-z)/(z_xf-z_f);
	j21 = j21_xf*(1.-interp)+j21_f*interp;
	cerr << "CIonise: Setting j21 = " << j21 << endl;
	break;
      }
      j21_xf = j21_f;
      z_xf = z_f;
    }
    if(j21<0) { 
      j21=j21_f;
      cerr << "CIonise: Setting j21 = " << j21 << " (redshift out of range)" <<  endl;
    }
  } else {
    j21=0.1;
    cerr << "CIonise: warning - couldn't open " << pathname<< ", using j21=" << j21 << endl;
  }
  
  nu_min = 13.6;
  nu_max = 100;

  nu_quant = 30;

  pJ=NULL;

  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(units::CUnit(units::den_protonsPerCm3));
  enToErgsPerG = pSim->getEnergyUnits().convertTo(units::CUnit(units::energy_erg)/units::CUnit(units::mass_g));

  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv = 0.001;

  flags = initJAndAlpha;
  
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
  
}

CIonise::CIonise(double gHg, double gHeg, double gHepg, CSimSnap *pSimi) {
  gammaHg = gHg;
  gammaHeg = gHeg;
  gammaHePg = gHepg;
  j21 = 0;
  alpha = 1;
  flags = initGammas;
  ne_thin_conv = 0.001;
  
  
  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(units::CUnit(units::den_protonsPerCm3));
  enToErgsPerG = pSim->getEnergyUnits().convertTo(units::CUnit(units::energy_erg)/units::CUnit(units::mass_g));

}

CIonise::CIonise(float *pJi, float nu_mini, float nu_maxi, int n_j, CSimSnap *pSimi) : pJ(pJi), nu_min(nu_mini), nu_max(nu_maxi), nu_quant(n_j) {
  
  
  
  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(units::CUnit(units::den_protonsPerCm3));
  enToErgsPerG = pSim->getEnergyUnits().convertTo(units::CUnit(units::energy_erg)/units::CUnit(units::mass_g));
  
  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv=0.001;
  flags = initJArray;
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
}



double CIonise::J(double nu) {
  // return J appropriate to how we have been setup.
  const float nu_0 = 13.6;

  if((flags&initJArray)!=0) {
    // interpolate
    int nu_index = nuToIndex(nu);
    float nu_low = nuFromIndex(nu_index);
    float nu_high = nu_low*nu_ratio;
    float frac = (nu-nu_low)/(nu_high-nu_low);

    if(nu_index<0 || nu_index>nu_quant-1)
      return 0;

    if(nu_index<nu_quant-1)
      return pJ[nu_index]*frac+pJ[nu_index+1]*(1.-frac);

    return pJ[nu_index];
  } 
  if((flags&initJAndAlpha)!=0) { 
    return  j21 * pow(nu/nu_0,alpha);
  }
  return 0.;
}

double CIonise::J(int nu_index) {
  const float nu_0=13.6;
  if((flags&initJArray)!=0) {
    return pJ[nu_index];
  }
  if((flags&initJAndAlpha)!=0) {
    return j21 * pow(nuFromIndex(nu_index)/nu_0,alpha);
  }
  return 0;
}


double CIonise::nuFromIndex(int index) {
  return nu_min*pow(nu_ratio,index);
}

int CIonise::nuToIndex(double nu ) {
  return (int)(log(nu/nu_min)/log(nu_ratio)); // always rounds down -> base index
}


// Cross-sections from Osterbrock (1989)

double CIonise::xsec_H(double nu) {
  const double A0 = 6.30e-18; 
  const double nu_thresh = 13.6; // ev
  
  if(nu<nu_thresh) return 0;
  else {
    
    double epsilon = sqrt(nu/nu_thresh - 0.99999); // almost 1, to prevent NaNs
    
    double a = A0 * pow((nu_thresh/nu),4)*exp(4-4*atan(epsilon)/epsilon);
    a/=(1-exp(-2*PI/epsilon));
    
    return a;
  }

}

double CIonise::xsec_HeP(double nu) {
  const double A0 = 6.30e-18 / 4.;  // cm^2
  const double nu_thresh = 54.4; // ev
  
  if(nu<nu_thresh) return 0; 
  else {

    double epsilon = sqrt(nu/nu_thresh - 0.999999);
    
    double a = A0 * pow((nu_thresh/nu),4)*exp(4-4*atan(epsilon)/epsilon);
    a/=(1-exp(-2*PI/epsilon));
    
    return a;
  }

}

double CIonise::xsec_He(double nu) {

  // const double aT = 1.58e-18;  // cm^2 - OLD
  const double aT = 7.42e-18;

  const double nu_thresh = 24.6; // ev
  if(nu<nu_thresh) return 0;
  else {
    // double a = aT * ( 1.34 * pow(nu/nu_thresh,-2.99) - 0.34 * pow(nu/nu_thresh,-3.99)); - OLD
    double a = aT * (1.66*pow(nu/nu_thresh,-2.05) - 0.66 * pow(nu/nu_thresh,-3.05));
    return a;
  }
}

// integrands to obtain photoionisation gammas:

double CIonise::gamma_H_integrand(double nu) {
  // could optimize if necessary
  
  const double nu_thresh = 13.6; // in ev

  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(units::planckInEvS * units::evInErg21*nu) * xsec_H(nu);
 
}

double CIonise::gamma_HeP_integrand(double nu) {
  
  const double nu_thresh = 54.4; // in ev
  const double nu_0 = 13.6; // in ev

  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(units::planckInEvS * units::evInErg21*nu) * xsec_HeP(nu);

}


double CIonise::gamma_He_integrand(double nu) {
 
  const double nu_thresh = 54.4; // in ev
  const double nu_0 = 13.6; // in ev

  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(units::planckInEvS * units::evInErg21*nu) * xsec_He(nu);

}

void CIonise:: calculateGammas() {
  gammaHg = gammaHeg = gammaHePg = 0;
  
  double nu = 13.6;
  
  //  cerr << "CIonise: calculating gammas...";

  for(nu=nu_min; nu<nu_max; nu*=nu_ratio) {
    double delta_nu = (nu_ratio - 1.) * nu;
    gammaHg += gamma_H_integrand(nu+delta_nu/2) * delta_nu;
    gammaHeg += gamma_He_integrand(nu+delta_nu/2) * delta_nu;
    gammaHePg += gamma_HeP_integrand(nu+delta_nu/2) * delta_nu;
  }
  // gammaHg = gammaHeg = gammaHePg = 0;
  // cerr << "done!" << endl;
  // 
}

double CIonise::electronFraction(CParticle *p) {
  double nH0, nHe0, nHp, nHep, nHepp, ne;
  calculateFractions(p, nH0, nHe0, nHp, nHep, nHepp, ne);
  return ne/(nH0+nHp);
}

double CIonise::electronDensity(CParticle *p) {
  double nH0, nHe0, nHp, nHep, nHepp, ne;
  calculateFractions(p, nH0, nHe0, nHp, nHep, nHepp, ne);
  return ne;
}

void CIonise::calculateFractions(CParticle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg, double local_gammaHeg, double local_gammaHePg) {

  
  
  if(local_gammaHeg == 0) local_gammaHeg = gammaHeg;
  if(local_gammaHg ==0 ) local_gammaHg= gammaHg;
  if(local_gammaHePg==0) local_gammaHePg=gammaHePg;


  if((flags & useTemp)!=0) {
    // assume temperature is correct
    calculateFractionsAtT(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);

  } else {
    // assumes internal energy, bot not necessarily temperature, is
    // correct
    
    const double PROTONMASS_OVER_BOLTZMANN = 1.2114e-8; // K s^2 / cm^2
    
    int maxIterations = 20;
    ne = 0;
    
    nH0 = 1;
    nHe0 = Y/(4*(1-Y));
    nHp = nHep = nHepp = 0;

    int n=0;
    float xne;
    do {
      
      xne=ne;
      double MeanWeight = 1/(1-0.75*Y+ne/(nH0+nHp)*(1-Y))*units::protonMassInG;
      // temperature in K:
      //cerr << Y << "\t" << MeanWeight << "\t" << p->u << "\t -  " << p->temp;
      p->temp = (MeanWeight/units::boltzmannInErgPerK) * (2.0/3) * p->u * enToErgsPerG;
      //cerr << "\t" << p->temp << "\t" << nH0 << "\t" << nHe0 << "\t" << nHp << "\t" << nHep << "\t" << nHepp << "\t" << ne << endl;
      calculateFractionsAtT(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);
      n++;
    } while(n<maxIterations && abs(ne-xne)>0.01);
    
  }
    
}


void CIonise::calculateFractionsAtT(CParticle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg, double local_gammaHeg, double local_gammaHePg) {
  
  // optically thin photoionisation equilibrium calculation
  // based on Katz, Weinberg, Hernquist 1996
  //
  // performed for particle n
  //
  // see also cooling.c in non-public version of Gadget
 
  // constants
  const double logTmin = 0.0;
  const double logTmax = 9.0;
  
  const int max_iterations = 100;

  double temp = p->temp;
 
  //temp*=1.5;
  double logT = log(temp)/log(10.);

  // rho is in units of Msol/kpc^3.
  double nH = ((double)p->rho*denToProtonsPerCm3)*(1-Y);
 

  ne=nH/100000000;

  if(logT<=logTmin) {
    // totally neutral

    nH0 = nH;
    nHe0 = nH*Y/(4-4*Y);
    nHp = 0;
    nHep = 0;
    nHepp = 0;
    ne = 0;
    return;
  }

  if(logT>=logTmax) {
    nH0 = 0;
    nHe0 = 0;
    nHp = nH;
    nHep = 0;
    nHepp = nH*Y/(4-4*Y);
    ne=nHp+2*nHepp;
    return;
  }


  // Recombination Rates

  double aHp = 3.34e-10 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7));

  double aHep = 1.5e-10 * pow(temp,-0.6353);
  double aHepp = 1.34e-9 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7)); 
  double ad = 1.9e-3 * pow(temp,-1.5) *exp(-4.7e5/temp) * (1+0.3*exp(-9.4e4/temp));


  // Collisonal Ionisation Rates

  double gammaH0 =  5.85e-11 * pow(temp,0.5) * exp(-157809.1/temp) / (1+pow(temp/1.e5,.5));
  double gammaHe0 = 2.38e-11 * pow(temp,0.5) * exp(-285335.4/temp) / (1+pow(temp/1.e5,.5));
  double gammaHeP =  5.68e-12 * pow(temp,0.5) * exp(-631515.0/temp) / (1+pow(temp/1.e5,.5));
  
  double ne_x=ne;
  // gammaH0 = 0;
  // gammaHe0 = gammaHeP = 0;
  
  

  // Physical constants for photoionisation calculations:


  int n_iterations = 0;
  
  do {
  
    ne_x = ne;
    nH0 = nH*aHp/(aHp+gammaH0+local_gammaHg/ne); // (33)
    nHp=nH-nH0; // (34)
    nHep = (Y/(4-4*Y))*nH/(1+(aHep+ad)/(gammaHe0+local_gammaHeg/ne)+(gammaHeP+local_gammaHePg/ne)/aHepp); // (35)
    nHe0=nHep*(aHep + ad)/(gammaHe0+local_gammaHeg/ne); // (36)
    nHepp=nHep*(gammaHeP+local_gammaHePg/ne)/aHepp; // (37)
    ne=nHp+nHep+2*nHepp; // (38)

    
    n_iterations++;

    
  } while(n_iterations < max_iterations && (abs(ne_x-ne)/nH)>ne_thin_conv);
  
}

void CIonise::setFlag(const int flag) {
  if((flags & flag)==0)
    flags+=flag;
}

void CIonise::unsetFlag(const int flag) {
  if((flags & flag)!=0)
    flags-=flag;
}

void CIonise::thickPPRadiative(CSimSnap *pSim) {
  
  if((flags & initGammas)!=0) {
    cerr << "CIonise: cannot call thickRadiative; spectrum unconstrained";
    return;
  }

  int numPart = pSim->getNumParticles();

  float *nPhot = (float*) malloc(sizeof(float) * numPart * nu_quant );
  float *delta_nPhot = (float*) malloc(sizeof(float) * numPart * nu_quant );

  cerr << "CIonise::thickPPRadiative: initialising...";

  double nu = nu_min;
  for(int nu_index=0; nu_index<nu_quant; nu_index++) {
    double delta_nu = (nu_ratio-1.)*nu;
    double N_Nu = 4*PI*J(nu+delta_nu/2)/(units::planckInErgS21*nu);	    
    for(int pn = 0; pn<numPart; pn++) {
      nPhot[pn*nu_quant+nu_index] = N_Nu;
    }
    nu*=nu_ratio;
  }

  cerr << "done!" << endl;

  int iter = 1;

  do {
    cerr << "CIonise::thickPPRadiative: iteration " << iter << endl;
    cerr << "CIonise::thickPPRadiative:   calculating local gammas & electron fractions...";
    for(int pn=0; pn<numPart; pn++) {
      
      double nu=nu_min;
      for(int nu_index = 0; nu_index<nu_quant; nu_index++) { 
	double delta_nu = (nu_ratio-1.)*nu;
	gammaHg += nPhot[nu_index] * xsec_H(nu+delta_nu/2) * delta_nu;
	gammaHeg += nPhot[nu_index] * xsec_He(nu+delta_nu/2) * delta_nu;
	gammaHePg += nPhot[nu_index] * xsec_HeP(nu+delta_nu/2) * delta_nu;
	nu*=nu_ratio;
      }


      CParticle *p = pSim->getParticle(pn);

      double nH0, nHe0, nHp, nHep, nHepp, ne;
     
      // TODO: this should calculate correct fractions for given U,
      // and return sensible temperature

      calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne);
      double temp = p->temp;

      
      nu=nu_min;
      for(int nu_index=0; nu_index<nu_quant; nu_index++) {
	double delta_nu = (nu_ratio-1.)*nu;
	delta_nPhot[nu_index+pn*nu_quant] = - nPhot[nu_index] * (nH0 * xsec_H(nu+delta_nu*0.5) + nHe0 * xsec_He(nu+delta_nu*0.5) +
								 nHep * xsec_HeP(nu+delta_nu*0.5));
	nu*=nu_ratio;
      }

      //           = number of absorptions per unit volume in this region
      //             due to this particle
      
      //  N.B. self-absorption is not taken into account, i.e. each particle must
      // separately be optically thin!!

    }
    
    cerr << "done!" << endl << "CIonise::thickPPRadiative:    propogating delta_nPhot...";

    double convRatio = ((pSim->getMassUnits()/pSim->getDensityUnits())/pow(pSim->getDistanceUnits(),2)).convertTo(units::len_cm);

    double delta_conv = 0;

    nu=nu_min;
    for(int nu_index=0; nu_index<nu_quant; nu_index++) {
      cerr << nu_index;
      
      
      double N_nu = 4*PI*J(nu)/(units::planckInErgS21*nu);	    
      
      for(int pn1=numPart-1; pn1!=0; --pn1) { // efficient for!
	
	CParticle *p1 = pSim->getParticle(pn1);
	
	nPhot[nu_index+pn1*nu_quant] = N_nu;
	
	for(int pn2=numPart-1; pn2!=0; --pn2) {
	
	  if(pn1!=pn2) {
	    
	    CParticle *p2 = pSim->getParticle(pn2);
	    float vol = p2->mass/p2->rho; // probably this estimation is not ideal
	    double dSquared = p2->squaredDistanceTo(*p1);
	    
	    // delta_nPhot[] is in units of s^-1 cm^-3
	    // vol is in units of (simlen)^3
	    // d is in units of simlen
	    // vol/d^2 is in units of kpc, 
	    nPhot[nu_index+pn1*nu_quant]+= (delta_nPhot[nu_index+pn2*nu_quant]*vol/(4*PI*dSquared))*convRatio;

	  } // if pn1!=pn2
	} // for pn2 (=source particle)
	if(nPhot[nu_index+pn1*nu_quant]<0) {
	  // cerr << "CIonise::thickPPRadiative:   nPhot is negative (" << nPhot[nu_index+pn1*nu_quant] << " vs " << N_nu << " for particle " << pn1 << endl;
	  nPhot[nu_index+pn1*nu_quant]=0;
	}

	delta_conv += ((N_nu-nPhot[nu_index+pn1*nu_quant])/N_nu)/nu_quant;
      } // for pn1 (=particle to update)
      
      nu*=nu_ratio;
    } // for nu_index
    cerr << "CIonise::thickPPRadiative: end of iteration " << iter << "; delta_conv = " << delta_conv << endl;
  } while(0==0);

}

void CIonise::thinRadiative(CSimSnap *pSim) {
  int numPart=  pSim->getNumParticles();
  for(int n=0; n<numPart; n++) {
    CParticle *p = pSim->getParticle(n);
    double nH0, nHe0, nHp, nHep, nHepp, ne;
    calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne);
    p->ne = ne/(nH0+nHp);
    p->nHp = nHp/(nH0+nHp);
  }
}

float CIonise::quickNegExp(float ex) {

  const double res[240] = {1,0.951229,0.904837,0.860708,0.818731,0.778801,0.740818,0.704688,0.67032,0.637628,0.606531,0.57695,0.548812,0.522046,0.496585,0.472366,0.449329,0.427415,0.40657,0.386741,0.367879,0.349938,0.332871,0.316637,0.301194,0.286505,0.272532,0.25924,0.246597,0.23457,0.22313,0.212248,0.201897,0.19205,0.182684,0.173774,0.165299,0.157237,0.149569,0.142274,0.135335,0.128735,0.122457,0.116484,0.110803,0.105399,0.100259,0.0953693,0.0907181,0.0862937,0.0820851,0.0780818,0.0742737,0.0706513,0.0672056,0.063928,0.0608102,0.0578444,0.0550233,0.0523398,0.0497872,0.047359,0.0450493,0.0428522,0.0407623,0.0387743,0.0368832,0.0350844,0.0333733,0.0317457,0.0301974,0.0287247,0.0273238,0.0259912,0.0247236,0.0235178,0.0223708,0.0212798,0.020242,0.0192548,0.0183157,0.0174224,0.0165727,0.0157645,0.0149956,0.0142643,0.0135686,0.0129068,0.0122774,0.0116786,0.011109,0.0105672,0.0100518,0.00956161,0.00909528,0.0086517,0.00822975,0.00782837,0.00744658,0.0070834,0.00673794,0.00640933,0.00609674,0.0057994,0.00551656,0.00524751,0.00499158,0.00474814,0.00451657,0.00429629,0.00408676,0.00388745,0.00369785,0.0035175,0.00334595,0.00318277,0.00302754,0.00287989,0.00273943,0.00260583,0.00247874,0.00235785,0.00224286,0.00213347,0.00202942,0.00193044,0.00183629,0.00174674,0.00166155,0.00158051,0.00150343,0.00143011,0.00136036,0.00129401,0.0012309,0.00117087,0.00111377,0.00105945,0.00100778,0.000958627,0.000911874,0.000867401,0.000825098,0.000784857,0.000746579,0.000710168,0.000675532,0.000642586,0.000611247,0.000581436,0.000553079,0.000526105,0.000500446,0.000476039,0.000452822,0.000430738,0.00040973,0.000389747,0.000370739,0.000352658,0.000335458,0.000319098,0.000303535,0.000288732,0.00027465,0.000261255,0.000248513,0.000236393,0.000224864,0.000213897,0.000203465,0.000193542,0.000184103,0.000175124,0.000166583,0.000158459,0.000150731,0.000143379,0.000136387,0.000129735,0.000123408,0.000117389,0.000111664,0.000106218,0.000101038,9.611e-05,9.14226e-05,8.69639e-05,8.27226e-05,7.86882e-05,7.48505e-05,7.12e-05,6.77275e-05,6.44244e-05,6.12823e-05,5.82936e-05,5.54505e-05,5.27462e-05,5.01737e-05,4.77267e-05,4.5399e-05,4.31849e-05,4.10787e-05,3.90753e-05,3.71695e-05,3.53568e-05,3.36324e-05,3.19921e-05,3.04318e-05,2.89476e-05,2.75358e-05,2.61929e-05,2.49155e-05,2.37003e-05,2.25444e-05,2.14449e-05,2.0399e-05,1.94042e-05,1.84578e-05,1.75576e-05,1.67013e-05,1.58868e-05,1.5112e-05,1.43749e-05,1.36739e-05,1.3007e-05,1.23726e-05,1.17692e-05,1.11952e-05,1.06492e-05,1.01298e-05,9.63579e-06,9.16585e-06,8.71882e-06,8.2936e-06,7.88911e-06,7.50436e-06,7.13836e-06,6.79022e-06,6.45906e-06};

  int n=(int) (-ex/0.05);
  if(n>239 || n<0) return 0.;
  
  return (float) res[n];
  
}

#define INDEX(Ax,Ay,Az,ADir,AQuant) Ax+nx*(Ay+ny*(Az+nz*(AQuant+nu_quant*ADir)))

void CIonise::thickRadiative(CGrid *pGrid) {

  if((flags & initGammas)!=0) {
    cerr << "CIonise: cannot call thickRadiative; spectrum unconstrained";
    return;
  }

  // flags are set once a warning has been given once, to prevent it happening again

  bool emiss_too_big_warning = false;
  bool abs_too_big_warning = false;

  int nx = pGrid->getNx();
  int ny = pGrid->getNy();
  int nz = pGrid->getNz();

  double dx = pGrid->getDx();
  double dy = pGrid->getDy();
  double dz = pGrid->getDz();

  double x1 = pGrid->getX1();
  double y1 = pGrid->getY1();
  double z1 = pGrid->getZ1();

  const int xp = 0;
  const int xm = 1;
  const int yp = 2;
  const int ym = 3;
  const int zp = 4;
  const int zm = 5;
  const int ndir = 6;

  const float nu_0 = 13.6;

  long est_mem = (nx*ny*nz*ndir*nu_quant*4*3)/(1024*1024);

  cout << "Estimated memory requirements: " << est_mem << "M" << endl;

  float *n_phot = new float[nx*ny*nz*ndir*nu_quant];
  float *n_phot_out = new float[nx*ny*nz*ndir*nu_quant];

  float *gammaArray=NULL; // for output information only, not used in procedure
  float *tauArray, *n136Array, *n246Array, *n544Array;

  if((flags & writeGam)!=0) {
    units::CUnit gamma_units = 1.e-15/units::CUnit(units::time_s);
    gammaArray = pGrid->createArray("gamH0g","Gamma_(H0+photon)",gamma_units);
    tauArray = pGrid->createArray("tau","tau");
    n136Array = pGrid->createArray("n13.6","n13.6");
    n246Array = pGrid->createArray("n24.6","n24.6");
    n544Array = pGrid->createArray("n54.4","n54.4");
  } 

    
  // set up array with optically thin approximation
  
  for(int dir=0; dir<ndir; dir++) {
    bool setRad = true;
    switch(dir) { 
    case xp:
      if((flags&noIlluminationX0)!=0) setRad=false;	
      break;
    case xm:
      if((flags&noIlluminationXN)!=0) setRad=false;
      break;
    case yp:
      if((flags&noIlluminationY0)!=0) setRad=false;
      break;
    case ym:
      if((flags&noIlluminationYN)!=0) setRad=false;
      break;
    case zp:
      if((flags&noIlluminationZ0)!=0) setRad=false;
      break;
    case zm:
      if((flags&noIlluminationZN)!=0) setRad=false;
      break;
    }
    double nu = nu_min;
    for(int nu_index = 0; nu_index<nu_quant; nu_index++) {
      
      float N_nu;
      double delta_nu = (nu_ratio - 1.) * nu;

      if(setRad) N_nu = 4*PI*J(nu+delta_nu/2)/(units::planckInErgS21*nu) /(float)ndir;
      else N_nu = 0;

      for(int x=nx-1;x>=0;--x) {
	for(int y=ny-1; y>=0; --y) {
	  for(int z=nz-1; z>=0; --z) {
			    
	    n_phot[INDEX(x,y,z,dir,nu_index)] = N_nu ; 
	 
	  } // for z
	} // for y
      } // for x
      nu*=nu_ratio;
    } // for nu
  } // for direction
  
  float *cur_n_phot = new float[nu_quant];

  int iter =0;
  int maxiter  =200;
  int miniter = (nx>ny)?nx:ny;
  if(miniter<nz) miniter=nz;

  double maxdelta=0;

  float max_grid_space = (dx>dy)?dx:dy;
  if(max_grid_space<dz) max_grid_space=dz;

  
  // zero optical depth array, which is used for interpolating photon numbers
  // within cells
  float *tau = new float[nx*ny*nz*nu_quant];
  for(int x=0;x<nx;x++) {
    for(int y=0;y<ny;y++) {
      for(int z=0;z<nz;z++) {
	for(int nu_index=0; nu_index<nu_quant; nu_index++) {
	  tau[INDEX(x,y,z,0,nu_index)]=0.;
	}
      }
    }
  }


  cerr << "CIonise: caching element cross-sections...";
  double *xsec_H_dnu = new double[nu_quant];
  double *xsec_He_dnu = new double[nu_quant];
  double *xsec_HeP_dnu = new double[nu_quant];

  double nu = nu_min;
  for(int nu_index = 0; nu_index<nu_quant; nu_index++) {
    
    double delta_nu = (nu_ratio-1.)*nu;
   
    xsec_H_dnu[nu_index] = xsec_H(nu+delta_nu/2) * delta_nu;
    xsec_He_dnu[nu_index] = xsec_He(nu+delta_nu/2) * delta_nu;
    xsec_HeP_dnu[nu_index] = xsec_HeP(nu+delta_nu/2) * delta_nu;
    
    nu*=nu_ratio;
  }
  
  cerr << "done!" << endl;
    
  do {
    if((flags & verbose)!=0)
      cerr << "CIonise::thickRadiative: Iteration " << iter << endl;
    
    // write-out backup

     
    std::ostringstream filename;
    filename << "ionise.backup." << iter;
    
    pGrid->write(filename.str(),CSimSnap::native);

    #ifdef SIMAN_OMP
    #pragma omp parallel for private(emiss_too_big_warning,abs_too_big_warning)
    #endif 
    for(int x=0;x<nx;x++) {
      for(int y=0;y<ny;y++) {
	for(int z=0; z<nz; z++) {
	
	  // now do equilibrium calculation for each particle in this cell


	  CSimSnap *thisCell = (*pGrid)[x][y][z];
	  int num_parts = thisCell->getNumParticles();

	  double ne_mean=0, nh0_mean=0, nhe0_mean=0, nhep_mean=0;
	  double Hrecomb=0, Herecomb=0, Heprecomb=0;
	  double local_gammaHg=0, local_gammaHeg=0, local_gammaHePg=0;


	  
	  for(int n=0; n<num_parts; n++) {
	    local_gammaHg=local_gammaHeg=local_gammaHePg=0;

	    CParticle *p = thisCell->getParticle(n);
	    
	    // STEP 1.
	    // Calculate distances through this cell that the various radiation paths
	    // have taken. This will require generalisation for ndir>6.

	    float dir_path_len[ndir];
	    dir_path_len[xp] = p->x-(x*dx+x1);
	    dir_path_len[xm] = dx-dir_path_len[xp];
	    dir_path_len[yp] = p->y-(y*dy+y1);
	    dir_path_len[ym] = dy-dir_path_len[ym];
	    dir_path_len[zp] = p->z-(z*dz+z1);
	    dir_path_len[zm] = dz-dir_path_len[zm];

	    


	    // STEP 2.
	    // Based on tau being constant locally, calculate local photon number
	    // in all wavebands, and use that to calculate the photoionisation gammas
	    
	    double nu = nu_min;
	    
	    for(int nu_index = 0; nu_index<nu_quant; nu_index++) {
	      
	      double delta_nu = (nu_ratio-1.)*nu;
	      // again, the following will require generalisation if n_dir>6
	      
	      cur_n_phot[nu_index]=0;
	      for(int dir=0; dir<ndir; dir++) {
		// proper exp is too slow
		float t= -dir_path_len[dir]*tau[INDEX(x,y,z,0,nu_index)];
		if(t>0) cerr << "AWOOGA!" << t << endl;
		float result = quickNegExp(t);
		
		cur_n_phot[nu_index]+=n_phot[INDEX(x,y,z,dir,nu_index)]*result;
	
	      } // for dir
	      
	      local_gammaHg += cur_n_phot[nu_index] * xsec_H_dnu[nu_index];
	      local_gammaHeg += cur_n_phot[nu_index] * xsec_He_dnu[nu_index];
	      local_gammaHePg += cur_n_phot[nu_index] * xsec_HeP_dnu[nu_index];
	      
	      nu*=nu_ratio;
	    } // for nu_index
	    

	    // STEP 3.
	    // Calculate fractions
    
	    double nH0, nHe0, nHp, nHep, nHepp, ne;

	    // TODO: this should calculate correct fractions for given U,
	    // and return sensible temperature
	    calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);
	    double temp = p->temp;

	    ne_mean+=ne/(double)num_parts;
	    nh0_mean+=nH0/(double)num_parts;
	    nhe0_mean+=nHp/(double)num_parts;
	    nhep_mean+=nHep/(double)num_parts;

	  
	    // Recombination Rates... 
	    // TODO: these are currently being calculated (at least)
	    // twice in one cycle (duplication in calculateFractions)

	    double aHp = 3.34e-10 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7));	  
	    double aHep = 1.5e-10 * pow(temp,-0.6353);
	    double aHepp = 1.34e-9 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7)); 
	    double ad = 1.9e-3 * pow(temp,-1.5) *exp(-4.7e5/temp) * (1+0.3*exp(-9.4e4/temp));


	    // Recombination rates directly to ground state (emit Ly-alpha)
	    double aHp1 = 1.5e-11 / sqrt(temp); // see Padmanabhan, vol1, p322
	    
	    Hrecomb+=(nHp* ne*aHp1)/(double)num_parts;
	    
	    //  Herecomb+=(nHep*ne*(aHep+ad))/(double)num_parts;
	    // Heprecomb+=(nHepp*ne*aHepp)/(double)num_parts;
	    
	    p->ne = ne/(nH0+nHp); 
	    p->nHp = nHp/(nH0+nHp);
	    
	    
	    if(gammaArray!=NULL) {
	      // write out gamma H0 for information
	      gammaArray[thisCell->deReference(n,1)]=(float)(log(local_gammaHg*(double)1.e15));
	      tauArray[thisCell->deReference(n,1)] = tau[INDEX(x,y,z,0,2)];
	      n136Array[thisCell->deReference(n,1)] = cur_n_phot[2];
	      n246Array[thisCell->deReference(n,1)] = cur_n_phot[nuToIndex(24.6)];
	      n544Array[thisCell->deReference(n,1)] = cur_n_phot[nuToIndex(54.4)];
	    }
	  }  // for n (particles in this cell)

	  
	  if(x==nx/2 && y==ny/2 && (flags&verbose)!=0) {
	    // midplane diagnosis
	    cerr << "(" << local_gammaHg << "," << ne_mean << "," << num_parts << ")";
	    if(z==nz-1) cerr << endl; else cerr << "\t";
	  }

	  // now calculate output fluxes using inferred optical depth

	  float j_emiss_tot = 0;

	  double nu = nu_min;
	  for(int nu_index = 0; nu_index<nu_quant; nu_index++) { 
	    
	    double delta_nu = (1-nu_ratio)*nu;

	    double local_tau = (nh0_mean * xsec_H(nu+delta_nu*0.5) + nhe0_mean * xsec_He(nu+delta_nu*0.5) +
				nhep_mean * xsec_HeP(nu+delta_nu*0.5)) * units::KpcInCm;
	    
	     
	    // local_tau = number of absorptions in this waveband per kpc per photon
	    tau[INDEX(x,y,z,0,nu_index)] = local_tau;
	    

	    double j_emiss_tot =0; // density of reemitted photons per length of box beam traversed (per cm^3 per kpc)

	    /*
	    if(nu+delta_nu/2>=13.6 && nu-delta_nu/2<13.6) {
	      // hydrogen II->I re-combination emits here
	      j_emiss_tot = Hrecomb*units::KpcInCm/((float)ndir);
	      if(j_emiss_tot*max_grid_space>cur_n_phot[nu_index]/5. && !(emiss_too_big_warning) && (flags&verbose)!=0) {
		cerr << "CIonise: warning - emission is comparable or greater than incident radiation; convergence may be effected" << endl;
		cerr << "CIonise:           Suggest finer grid resolution? ("<<j_emiss_tot*dx<<" / " << cur_n_phot[nu_index] << ")" << endl;
		cerr << "CIonise:           (This warning will be supressed for further iterations)" << endl;
		emiss_too_big_warning = true;
	      }
	    

	    }
	  
	    if(nu<=24.6 && nu+delta_nu>24.6) {
	      // helium II->I re-combination emits here
	      j_emiss_tot = Herecomb/(units::cInKpcPerS*ndir);
	    }
	  
	    if(nu<=54.4 && nu+delta_nu>54.4) {
	      // helium III->II re-combination emits here
	      j_emiss_tot = Heprecomb/(units::cInKpcPerS*ndir);
	    }
	    */

	    //	    cout << tau[INDEX(x,y,z,0,2)] << endl;
	    j_emiss_tot = 0; 
	    	  
	    if(local_tau*max_grid_space>0.2 && !(abs_too_big_warning) && (flags&verbose)!=0) {
	      cerr << "CIonise: warning - absorption is comparable or greater than incident radiation; convergence may be effected" << endl;
	      cerr << "CIonise:           Suggest finer grid resolution? (local_tau = " << local_tau*max_grid_space << ")" << endl;
	      cerr << "CIonise:           (This warning will be supressed for further iterations)" << endl;
	      abs_too_big_warning = true;
	    }
	    
	
	    // TODO: this will require generalisation if more than 6 directions
	    // are to be allowed
	    
	    n_phot_out[INDEX(x,y,z,xm,nu_index)]= n_phot[INDEX(x,y,z,xm,nu_index)]*exp(-dx*local_tau) + j_emiss_tot*dx;
	    n_phot_out[INDEX(x,y,z,xp,nu_index)]= n_phot[INDEX(x,y,z,xp,nu_index)]*exp(-dx*local_tau) + j_emiss_tot*dx;
	    n_phot_out[INDEX(x,y,z,yp,nu_index)]= n_phot[INDEX(x,y,z,yp,nu_index)]*exp(-dy*local_tau) + j_emiss_tot*dy;
	    n_phot_out[INDEX(x,y,z,ym,nu_index)]= n_phot[INDEX(x,y,z,ym,nu_index)]*exp(-dy*local_tau) + j_emiss_tot*dy;
	    n_phot_out[INDEX(x,y,z,zp,nu_index)]= n_phot[INDEX(x,y,z,zp,nu_index)]*exp(-dz*local_tau) + j_emiss_tot*dz;
	    n_phot_out[INDEX(x,y,z,zm,nu_index)]= n_phot[INDEX(x,y,z,zm,nu_index)]*exp(-dz*local_tau) + j_emiss_tot*dz;
	  
	    nu = nu * nu_ratio;

	  } // for nu_index (updating all directions of output flux)
 

	} // for z
      } // for y
    } // for x
  	    
  
    // Propogate radiation (copy output fluxes into input fluxes)

    
    int x_max =0, y_max=0, z_max=0;

    double delta,tot_nphot=0,tot_nphot_at_maxdelta=0;

#ifdef SIMAN_OMP
#pragma omp parallel for private(delta,tot_nphot) default(shared)
#endif
    for(int x=0;x<nx;x++) {
      for(int y=0;y<ny;y++) {
	for(int z=0; z<nz; z++) {

	  delta = 0;
	  tot_nphot=0;

	  for(int nu_index=0; nu_index<nu_quant; nu_index++) {
	   

	    // TODO: this will require generalisation for ndir>6

	  
	    if(x<nx-1) {
	      delta+= abs(n_phot[INDEX(x,y,z,xm,nu_index)] - n_phot_out[INDEX(x+1,y,z,xm,nu_index)]);
	      n_phot[INDEX(x,y,z,xm,nu_index)] = n_phot_out[INDEX(x+1,y,z,xm,nu_index)] ;
	    } 
	  
	    if(x>0) {
	      delta+=abs(n_phot[INDEX(x,y,z,xp,nu_index)] - n_phot_out[INDEX(x-1,y,z,xp,nu_index)]);
	      n_phot[INDEX(x,y,z,xp,nu_index)] = n_phot_out[INDEX(x-1,y,z,xp,nu_index)] ;
	    }

	    
	    if(y<ny-1) {
	      delta+=abs(n_phot[INDEX(x,y,z,ym,nu_index)] - n_phot_out[INDEX(x,y+1,z,ym,nu_index)]);
	      n_phot[INDEX(x,y,z,ym,nu_index)] = n_phot_out[INDEX(x,y+1,z,ym,nu_index)] ;
	    }
	    if(y>0) {
	      delta+=abs(n_phot[INDEX(x,y,z,yp,nu_index)] - n_phot_out[INDEX(x,y-1,z,yp,nu_index)]);
	      n_phot[INDEX(x,y,z,yp,nu_index)] = n_phot_out[INDEX(x,y-1,z,yp,nu_index)] ;
	    }
	    
	    if(z<nz-1) {
	      delta+=abs( n_phot[INDEX(x,y,z,zm,nu_index)] - n_phot_out[INDEX(x,y,z+1,zm,nu_index)]);
	      n_phot[INDEX(x,y,z,zm,nu_index)] = n_phot_out[INDEX(x,y,z+1,zm,nu_index)] ;
	    }
	    if(z>0) {
	      delta+=abs(n_phot[INDEX(x,y,z,zp,nu_index)] - n_phot_out[INDEX(x,y,z-1,zp,nu_index)]);
	      n_phot[INDEX(x,y,z,zp,nu_index)] = n_phot_out[INDEX(x,y,z-1,zp,nu_index)] ;
	    }
	    
	    
	    if((flags & periodicBoundsX)!=0) {
	     
	      if(x==nx-1) {
		delta+= abs(n_phot[INDEX(x,y,z,xm,nu_index)] - n_phot_out[INDEX(0,y,z,xm,nu_index)]);
		n_phot[INDEX(x,y,z,xm,nu_index)] = n_phot_out[INDEX(0,y,z,xm,nu_index)] ;
	      } 
	      
	      if(x==0) {
		delta+=abs(n_phot[INDEX(x,y,z,xp,nu_index)] - n_phot_out[INDEX(nx-1,y,z,xp,nu_index)]);
		n_phot[INDEX(x,y,z,xp,nu_index)] = n_phot_out[INDEX(nx-1,y,z,xp,nu_index)] ;
	      }
	      
	    }

	    if((flags & periodicBoundsY)!=0) {
	      
	      if(y==ny-1) {
		delta+=abs(n_phot[INDEX(x,y,z,ym,nu_index)] - n_phot_out[INDEX(x,0,z,ym,nu_index)]);
		n_phot[INDEX(x,y,z,ym,nu_index)] = n_phot_out[INDEX(x,0,z,ym,nu_index)] ;
	      }
	      if(y==0) {
		delta+=abs(n_phot[INDEX(x,y,z,yp,nu_index)] - n_phot_out[INDEX(x,ny-1,z,yp,nu_index)]);
		n_phot[INDEX(x,y,z,yp,nu_index)] = n_phot_out[INDEX(x,ny-1,z,yp,nu_index)] ;
	      }
	    }

	    if((flags & periodicBoundsZ)!=0) {

	      if(z==nz-1) {
		delta+=abs( n_phot[INDEX(x,y,z,zm,nu_index)] - n_phot_out[INDEX(x,y,0,zm,nu_index)]);
		n_phot[INDEX(x,y,z,zm,nu_index)] = n_phot_out[INDEX(x,y,0,zm,nu_index)] ;
	      }
	      if(z==0) {
		delta+=abs(n_phot[INDEX(x,y,z,zp,nu_index)] - n_phot_out[INDEX(x,y,nz-1,zp,nu_index)]);
		n_phot[INDEX(x,y,z,zp,nu_index)] = n_phot_out[INDEX(x,y,nz-1,zp,nu_index)] ;
	      }
	   
	    }
	    

	    // more efficient to do summation explicitly than in loop:
	    tot_nphot+=n_phot[INDEX(x,y,z,xm,nu_index)]+n_phot[INDEX(x,y,z,ym,nu_index)]+n_phot[INDEX(x,y,z,zm,nu_index)]
	      + n_phot[INDEX(x,y,z,xp,nu_index)]+n_phot[INDEX(x,y,z,yp,nu_index)]+n_phot[INDEX(x,y,z,zp,nu_index)];

	

	  } // for nu_index
	  
#ifdef SIMAN_OMP
#pragma omp critical (deltacomp)
#endif
	  if(delta>maxdelta) {
	    maxdelta=delta/tot_nphot;
	    tot_nphot_at_maxdelta = tot_nphot;
	    x_max = x;
	    y_max = y;
	    z_max = z;
	  }
	} // for z
      } // for y
    } // for x
    // (implicit end of parallel region)

    if((flags&verbose)!=0) cerr << endl << "Maxdelta = " << maxdelta << " ("<<x_max<<","<<y_max<<","<<z_max<<": "<<tot_nphot_at_maxdelta<<")"<< endl;

    iter++;
  } while((iter<miniter || maxdelta>0.001) && iter<maxiter);

  delete[] tau;
  delete[] n_phot;
  delete[] n_phot_out;
  delete[] cur_n_phot;
  delete[] xsec_H_dnu;
  delete[] xsec_He_dnu;
  delete[] xsec_HeP_dnu;

}

