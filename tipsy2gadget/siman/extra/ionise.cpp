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
#include <limits>

using namespace std;

namespace gasoline {

  const double  CL_eV_per_K = 8.6173423e-5;
  const double  CL_MAX_NEG_EXP_ARG = -500; 
  /*-----------------------------------------------------------------
   *     Collisional Ionization rates
   *-----------------------------------------------------------------*/
  /*     H + e- -> H+ + 2e-  Janev et al. 1987 (Abel 1996) */
  double clRateCollHI( double T ) {
    double TL,arg;

    TL = log(T*CL_eV_per_K);
    arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 +
						  TL*(1.56315498 +
						      TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 +
												   TL*(1.11954395e-4  + TL*(-2.03914985e-6))))))));
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return exp ( arg );
  }

  /*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
  double clRateCollHeI( double T ) {
    double TL,arg;

    TL = log(T*CL_eV_per_K);
    arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 +
						  TL*(3.05803875 +
						      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 +
												   TL*(2.06723616e-4  + TL*(-3.64916141e-6))))))));
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return exp( arg );
  }

  /*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
  double clRateCollHeII( double T ) {
    double TL,arg;

    TL = log(T*CL_eV_per_K);
    arg = -68.71040990  + TL*(43.93347633   + TL*(-18.4806699 +
						  TL*(4.70162649 +
						      TL*(-0.76924663    + TL*(8.113042e-2   + TL*(-5.32402063e-3 +
												   TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))));
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return exp( arg );
  }

  /*-----------------------------------------------------------------
   *     Radiative Recombination rates
   *-----------------------------------------------------------------*/
  /*     H+ + e- -> H + gam  Verner & Ferland 1996 */
  double clRateRadrHII( double T ) {
    double Tsq = sqrt(T);

    return 7.982e-11/( Tsq*0.563615 *
		       pow(1+Tsq*0.563615,0.252) * pow(1+Tsq*1.192167e-3,1.748));
  }

  /*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
  double clRateRadrHeII( double T ) {
    /*
     * Note that these functions do not meet perfectly at 1e6 -- 2% difference
     * The derivatives are different there also: So the apparent error is large
     */
    double Tsq = sqrt(T);
    if (T < 1e6)
      return  3.294e-11/( Tsq*0.253673 *
			  pow(1+Tsq*0.253673,0.309) * pow(1+Tsq*1.649348e-4,1.691));
    else
      return  9.356e-10/( Tsq*4.841607 *
			  pow(1+Tsq*4.841607,0.2108) * pow(1+Tsq*4.628935e-4,1.7892));
  }

  /*     He+ + e- -> He + gam  dielectronic  Aldovandi&Pequignot 1973 (Black 1981) */
  double clRateDielHeII( double T ) {
    double T_inv = 1.0/T,arg;

    arg = -4.7e5*T_inv;
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return 1.9e-3*pow(T,-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv));
  }

  /*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
  double clRateRadrHeIII( double T ) {
    double Tsq = sqrt(T);

    return 1.891e-10/( Tsq*0.326686 *
		       pow(1+Tsq*0.326686,0.2476) * pow(1+Tsq*6.004084e-4,1.7524));
  }

  // --------------------------------------------------------------
  // Lots of pretty constants collated
  // --------------------------------------------------------------

  const double CL_Cbremss1 (1.426e-27);
  const double CL_al       (0.79464);
  const double CL_bl       (0.1243);
  const double CL_ar       (2.13164);
  const double CL_br       (-0.1240);
  const double CL_B_gm         (6.022e23*(938.7830/931.494));
  const double CL_k_Boltzmann  (1.38066e-16);
  const double CL_eV_erg      ( 1.60219e-12);
  const double CL_eHI     (13.60*CL_eV_erg);
  const double  CL_eHeI    (24.59*CL_eV_erg);
  const double CL_eHeII   (54.42*CL_eV_erg);
  const double CL_E2HeII  (3.0*13.6*CL_eV_erg);


  /*-----------------------------------------------------------------
   *     Collisional cooling
   *-----------------------------------------------------------------*/

  const double clCoolCollHI(CL_eHI*CL_B_gm) ;
  const double clCoolCollHeI(CL_eHeI*CL_B_gm) ;
  const double clCoolCollHeII(CL_eHeII*CL_B_gm);
  const double clCoolDielHeII((CL_E2HeII+CL_eHeI)*CL_B_gm);

  /*-----------------------------------------------------------------
   *     Bremsstrahlung   
   *-----------------------------------------------------------------*/



  double clCoolBrem1( double T ) {
    double Tlog10, Tsq;

    Tlog10 = log10(T);
    Tsq = sqrt(T);
    if (T < 3.2e5) 
      return Tsq*CL_Cbremss1*(CL_al+CL_bl*Tlog10)*CL_B_gm;
    else   
      return Tsq*CL_Cbremss1*(CL_ar+CL_br*Tlog10)*CL_B_gm;
  }

  const double CL_alog4   (0.602059991);
  const double CL_alII    (4.0*(CL_al-CL_bl*CL_alog4));
  const double CL_blII    (4.0*CL_bl);
  const double CL_arII    (4.0*(CL_ar-CL_br*CL_alog4));
  const double CL_brII    (4.0*CL_br);

  double clCoolBrem2( double T ) {
    double Tlog10, Tsq;

    Tlog10 = log10(T);
    Tsq = sqrt(T);

    if (T<12.8e5) 
      return Tsq*CL_Cbremss1*(CL_alII+CL_blII*Tlog10)*CL_B_gm;
    else
      return Tsq*CL_Cbremss1*(CL_arII+CL_brII*Tlog10)*CL_B_gm;
  }

  /*-----------------------------------------------------------------
   *     Cooling multiplier for radiative recombination
   *-----------------------------------------------------------------*/
  const double CL_aHII  (0.0215964);
  const double CL_b     (0.270251);

  double clCoolRadrHII( double T ) {

    double Tpow;
    
    Tpow=pow(T,CL_b);
    /* return CL_B_gm*(CL_eHI+exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T); */
    /* Though 13.6eV is lost to the Gas as radiation, calculating the
     * Energy using u = 3/2 k T requires we don't subtract it here.
     */
    return CL_B_gm*(exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T);
  }
 
  double clCoolRadrHeII( double T ) {

    double Tpow;
    
    Tpow=pow(T,CL_b);
    /* return CL_B_gm*(CL_eHeI+exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T); */
    return CL_B_gm*(exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T);
  }

  double clCoolRadrHeIII( double T ) {
    double Tpow;
    
    Tpow=pow(T,CL_b);
    /* return CL_B_gm*(CL_eHeII+exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T); */
    return CL_B_gm*(exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T);
  }

  /*-----------------------------------------------------------------
   *     Line Cooling
   *-----------------------------------------------------------------*/
  /*      CEN (1992, Ap.J.Suppl 78,341) ADVOCATES MULTIPLYING EACH OF 
   *      THESE RATES BY Cen_correctn - HE CLAIMS THIS GIVES THE RIGHT
   *      HIGH T LIMIT FOR PROCESSES INVOLVING A FREE EL INTERACTING 
   *      WITH AN ORBITAL ELECTRON ?? */
  const double CL_aHI  (7.5e-19);
  const double CL_bHI  (1.18348e05);

  double clCoolLineHI( double T ) {
    double T_inv, arg;
    double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

    T_inv=1.0/T;
    arg = -CL_bHI*T_inv;
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return CL_B_gm*CL_aHI*exp( arg )*Cen_correctn;
  }

  const double CL_aHeI   (9.10e-27);
  const double CL_bHeI   (1.3179e04);
  const double CL_p_HeI  (0.1687);

  double clCoolLineHeI( double T ) {
    double T_inv,arg;
    double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

    T_inv=1.0/T;
    arg = -CL_bHeI*T_inv;
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return CL_B_gm*CL_aHeI*exp(-CL_bHeI*T_inv)*pow(T_inv,CL_p_HeI)*Cen_correctn;
  }

  const double CL_aHeII   (5.54e-17);
  const double CL_bHeII   (4.73638e05);
  const double CL_p_HeII  (0.397);

  double clCoolLineHeII( double T ) {

    double T_inv,arg;
    double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

    T_inv=1.0/T;
    arg = -CL_bHeII*T_inv;
    if (arg < CL_MAX_NEG_EXP_ARG) return 0;
    return CL_B_gm*CL_aHeII*exp(-CL_bHeII*T_inv)*pow(T_inv,CL_p_HeII)*Cen_correctn;
  }



  const double CL_Ccomp0 (0.565e-9);
  const double CL_Tcmb0  (2.735);
  const double CL_Ccomp  (CL_Ccomp0*CL_Tcmb0);

}


Ionise::Ionise(double j21i, double alphai, SimSnap *pSimi) {

  j21 = j21i;
  alpha = alphai;
  
  nu_min = 13.6;
  nu_max = 100;

  nu_quant = 30;

  pJ=NULL;

  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(Unit("m_p cm^-3"),pSim);
  enToErgsPerG = pSim->getEnergyUnits().convertTo(Unit("ergs g^-1"),pSim);

  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv = 0.001;

  flags = initJAndAlpha;
  
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
  
}

double Ionise::coolingRate(double nden, double T, double ne, double nHI, double nHII, double nHeI, double nHeII, double nHeIII) 
{
  // output is ergs/g/second

  double rate = ne * (
		      gasoline::clCoolBrem1(T) * (nHII+nHeII) +
		      gasoline::clCoolBrem2(T) * (nHeIII) +
		      gasoline::clCoolRadrHII(T) * (nHII) * gasoline::clRateRadrHII(T) +
		      gasoline::clCoolRadrHeII(T) * (nHeII) * gasoline::clRateRadrHeII(T) +
		      gasoline::clCoolRadrHeIII(T) * (nHeIII) * gasoline::clRateRadrHeIII(T) +
		      gasoline::clCoolCollHI * nHI * gasoline::clRateCollHI(T) +
		      gasoline::clCoolCollHeI * nHeI * gasoline::clRateCollHeI(T) +
		      gasoline::clCoolCollHeII * nHeII * gasoline::clRateCollHeII(T) +
		      gasoline::clCoolDielHeII * nHeII * gasoline::clRateDielHeII(T) +
		      gasoline::clCoolLineHI(T) * nHI +
		      gasoline::clCoolLineHeI(T) * nHeI +
		      gasoline::clCoolLineHeII(T) * nHeII) * nden;
 
  return rate;

}



double Ionise::heatingRate(double nden, double nHI, double nHeI, double nHeII, double epsHg, double epsHeg, double epsHePg) {

  // nden in units cm^-3
  // epsHg is in units of eV
  // output is in units of ergs/g/second
  // (hence conversion factor)

  return (nHI   * epsHg +
         nHeI  * epsHeg +
	  nHeII * epsHePg) * 9.5788341e11;

}

void Ionise::cooling() {
  SimanArray &ca(pSim->createArray("netcool","Net Cooling Rate",Unit("ergs g^-1 s^-1")));
  SimanArray &co(pSim->createArray("cool","Cooling Rate",Unit("ergs g^-1 s^-1")));
  SimanArray &ho(pSim->createArray("heat","Photoheating rate",Unit("ergs g^-1 s^-1")));
  SimanArray &ct(pSim->createArray("cooltime","Cooling Time",Unit("1e9 yr")));

  
  SimanArray &a_nHII(pSim->getArray("nHII"));
  SimanArray &a_nHeII(pSim->getArray("nHeII"));
  SimanArray &a_nHeIII(pSim->getArray("nHeIII"));
  SimanArray &a_ne(pSim->getArray("ne"));
  SimanArray &rho(pSim->getArray("rho"));
  SimanArray &temp(pSim->getArray("temp"));
  SimanArray &u(pSim->getArray("u"));

  double fac_rho = rho.getUnits().convertTo(Unit("m_p cm^-3"),pSim);
  double fac_u = u.getUnits().convertTo(Unit("erg g^-1"),pSim);

  SimanArray *epsHIg=NULL, *epsHeIg=NULL, *epsHeIIg=NULL;
  double eHIg, eHeIg, eHeIIg;

  
  double z = pSim->getRedshift(); 

  double compton = pow((1+z)*(gasoline::CL_Ccomp),4.0)*gasoline::CL_B_gm;

  try {
    epsHIg = & (pSim->getArray("epsHIg"));
    epsHeIg = &(pSim->getArray("epsHeIg"));
    epsHeIIg = & (pSim->getArray("epsHeIIg"));
  } catch(UnknownArray &e) {
    cerr << "Note: no heating rates in arrays; using optically thin values instead";
    eHIg = epsilonHg;
    eHeIg = epsilonHeg;
    eHeIIg = epsilonHePg;
    
  }


  for(unsigned int i=0; i<pSim->getNumParticles(); i++) {
    double nHI =0, nHII=0,nHeI=0,nHeII=0,nHeIII=0, ne=0, T=0, tot_n;
    double nHe, nH;
    
    nH = 0.76*fac_rho*rho[i];
    nHe = 0.24*fac_rho*rho[i]/4.;
   
    nHII = a_nHII[i]*nH;
    nHeII = a_nHeII[i]*nH;
    nHeIII = a_nHeIII[i]*nH;
    
    nHI = nH-nHII;
    nHeI = nHe-nHeII-nHeIII;

    T = temp[i];
    ne = a_ne[i]*nH;

    double massfac = 1/(nH+nHe*4); 

    // cooling fns require number density per baryon density

    ne*=massfac;
    nH*=massfac;
    nHe*=massfac;
    nHI*=massfac;
    nHII*=massfac;
    nHeI*=massfac;
    nHeII*=massfac;
    nHeIII*=massfac;
    

    if(epsHeIIg!=NULL) {
      eHIg = double((*epsHIg)[i])*1.e-15;
      eHeIg = double((*epsHeIg)[i])*1.e-15;
      eHeIIg = double((*epsHeIIg)[i])*1.e-15;
    } 
    
    co[i] = ne * (T-2.735)*(1+z) * compton + coolingRate(fac_rho*rho[i],T,ne,nHI,nHII,nHeI,nHeII,nHeIII) ;
    ho[i] = heatingRate(fac_rho*rho[i],nHI,nHeI,nHeII, eHIg, eHeIg, eHeIIg);
    
    ca[i]=  co[i] - ho[i];
    
    ct[i] = (double(u[i])*fac_u/double(ca[i]))/double(3.15556926e16); // factor is s->1e9 yr
  }

}

Ionise::Ionise(SimSnap *pSimi, float z) {

  if(z<0)
    z=pSimi->getRedshift();


  string pathname = (string) getenv("SIMAN_DATA");
  pathname+="/UVNORM";
  
  ifstream file_uv_norm(pathname.c_str());

  
  j21 = -1;
  alpha = -1.6;

  nu_min = 13.6;
  nu_max = 100;
  nu_quant = 30;
  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  pJ=NULL;
  flags = initJAndAlpha;

  if(file_uv_norm.good()) {
    string p;
    file_uv_norm >> p;
    if(p!="#") throw(SimanException("Format error in UVNORM"));
  
    int use_col = -1;
    int use_block = 0;
    while(!file_uv_norm.eof() && use_col==-1) {
      file_uv_norm >> p;
      use_block++;
      if(p!="lambd/z") throw(SimanException("Format error in UVNORM"));
      int z_index;
      for(z_index=0;z_index<10;z_index++) {
	float cur_z;
	file_uv_norm >> cur_z;
	if(cur_z>=z-0.002 && use_col==-1) {
	  use_col = z_index;
	  
	}
      }

      if(use_col==-1) {
	// skip to beginning of next block
	while(!file_uv_norm.eof() && p!="#")
	  file_uv_norm >> p;


      } else {
	// read in data from this block

	cerr << "Ionise: using UVNORM block " << use_block << ", column " << use_col << endl;
	double lambda,  nu, j_nu; 

	vector< pair<double,double> > nu_j_vec;

	for(int n=0; n<432; n++) {
	 
	  file_uv_norm >> lambda;
	 
	  for(int i=0; i<=use_col; i++) {
	    file_uv_norm >> j_nu;
	  }

	  for(int i=use_col+1; i<10; i++) {
	    file_uv_norm >> p;	  
	  }
	  nu = 1.23984 / (0.0001*lambda); // hc/1eV = 1.23984 microns

	  if(nu>13.5 && nu <110) {
	    nu_j_vec.push_back(pair<double,double>(nu,j_nu*1.e21));
	  }
	}

	pJ = new float[nu_quant];
	int nu_index=0;
	for(nu=nu_min; nu<nu_max; nu*=nu_ratio) {
	  double delta_nu = (nu_ratio - 1.) * nu;
	  pJ[nu_index]=(float) siman::interpolate(nu_j_vec,nu);
	  //  cout << nu << "\t" << pJ[nu_index] << endl;
	  nu_index++;
	}
	flags = initJArray;
      }
    }
    /*
      Old code for j / alpha normalization

    float z_f=0, j21_f=0, z_xf=100, j21_xf=0;
    while(!file_uv_norm.eof()) {
      
      file_uv_norm >> z_f >> j21_f;
      if(z_f<z) {
	float interp = (z_xf-z)/(z_xf-z_f);
	j21 = j21_xf*(1.-interp)+j21_f*interp;
	cerr << "Ionise: Setting j21 = " << j21 << endl;
	break;
      }
      j21_xf = j21_f;
      z_xf = z_f;
    }
    if(j21<0) { 
      j21=j21_f;
      cerr << "Ionise: Setting j21 = " << j21 << " (redshift out of range)" <<  endl;
    }
    */
  } else {
    j21=0.1;
    cerr << "Ionise: warning - couldn't open " << pathname<< ", using j21=" << j21 << endl;
  }
  
  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(Unit("m_p cm^-3"),pSim);
  enToErgsPerG = pSim->getEnergyUnits().convertTo(Unit("erg g^-1"),pSim);

  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv = 0.001;

 
  
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
  
}

Ionise::Ionise(double gHg, double gHeg, double gHepg, SimSnap *pSimi) {
  gammaHg = gHg;
  gammaHeg = gHeg;
  gammaHePg = gHepg;
  j21 = 0;
  alpha = 1;
  flags = initGammas;
  ne_thin_conv = 0.001;
  
  
  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(Unit("m_p cm^-3"),pSim);
  enToErgsPerG = pSim->getEnergyUnits().convertTo(Unit("erg g^-1"),pSim);

}

Ionise::Ionise(float *pJi, float nu_mini, float nu_maxi, int n_j, SimSnap *pSimi) : pJ(pJi), nu_min(nu_mini), nu_max(nu_maxi), nu_quant(n_j) {
  
  
  
  pSim=pSimi;
  Y=pSim->getHeliumMassFrac();
  denToProtonsPerCm3 = pSim->getDensityUnits().convertTo(Unit("m_p cm^-3"),pSim);
  enToErgsPerG = pSim->getEnergyUnits().convertTo(Unit("erg g^-1"),pSim);
  
  nu_ratio = pow((double)nu_max/nu_min,1./(double)nu_quant);
  ne_thin_conv=0.001;
  flags = initJArray;
  calculateGammas();
  cerr << "gamma hg = " << gammaHg << "\t" << "gamma Heg = " << gammaHeg << "\t gammaHePg = " << gammaHePg << endl;
}



double Ionise::J(double nu) {
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

double Ionise::J(int nu_index) {
  const float nu_0=13.6;
  if((flags&initJArray)!=0) {
    return pJ[nu_index];
  }
  if((flags&initJAndAlpha)!=0) {
    return j21 * pow(nuFromIndex(nu_index)/nu_0,alpha);
  }
  return 0;
}


double Ionise::nuFromIndex(int index) {
  return nu_min*pow(nu_ratio,index);
}

int Ionise::nuToIndex(double nu ) {
  return (int)(log(nu/nu_min)/log(nu_ratio)); // always rounds down -> base index
}


// Cross-sections from Osterbrock (1989)

double Ionise::xsec_H(double nu) {
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

double Ionise::xsec_HeP(double nu) {
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

double Ionise::xsec_He(double nu) {
  const double Eth = 24.59;
  const double E0 = 13.61;
  const double sig0 = 9.492e-16; // cm^2
  const double ya = 1.469;
  const double P = 3.188;
  const double yw= 2.039;
  const double y0=4.434e-1;
  const double y1=2.136;

  if(nu<Eth) return 0;
  else {
    double x = nu/E0-y0;
    double y = sqrt(x*x+y1*y1);
    return sig0 * ((x-1.)*(x-1.)+yw*yw)*pow(y,0.5*P-5.5)*pow((1+sqrt(y/ya)),-P);
  }

  /*

  // const double aT = 1.58e-18;  // cm^2 - OLD
  const double aT = 7.42e-18;

  const double nu_thresh = 24.6; // ev
  if(nu<nu_thresh) return 0;
  else {
    // double a = aT * ( 1.34 * pow(nu/nu_thresh,-2.99) - 0.34 * pow(nu/nu_thresh,-3.99)); - OLD
    double a = aT * (1.66*pow(nu/nu_thresh,-2.05) - 0.66 * pow(nu/nu_thresh,-3.05));
    return a;
  }
  */

}

// integrands to obtain photoionisation gammas:

double Ionise::gamma_H_integrand(double nu) {
  // could optimize if necessary
  
  const double nu_thresh = 13.6; // in ev

  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(constants::planckInEvS * constants::evInErg21*nu) * xsec_H(nu); 
  // N.B. extra "planck constant" is there because integral will be performed over nu in eV, not nu in Hz
 
}

double Ionise::gamma_HeP_integrand(double nu) {
  
  const double nu_thresh = 54.4; // in ev


  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(constants::planckInEvS * constants::evInErg21*nu) * xsec_HeP(nu);

}


double Ionise::gamma_He_integrand(double nu) {
 
  const double nu_thresh = 24.6; // in ev


  if(nu<nu_thresh)
    return 0;
  else
    return 4*PI*J(nu)/(constants::planckInEvS * constants::evInErg21*nu) * xsec_He(nu);

}

void Ionise:: calculateGammas() {
  gammaHg = gammaHeg = gammaHePg = epsilonHg = epsilonHeg = epsilonHePg = 0;
  
  double nu = 13.6;
  
  //  cerr << "Ionise: calculating gammas...";

  for(nu=nu_min; nu<nu_max; nu*=nu_ratio) {
    double delta_nu = (nu_ratio - 1.) * nu;

    //    cerr << nu << "\t" << delta_nu << "\t" << xsec_H(nu)*delta_nu << endl;

    gammaHg += gamma_H_integrand(nu+delta_nu/2) * delta_nu;
    gammaHeg += gamma_He_integrand(nu+delta_nu/2) * delta_nu;
    gammaHePg += gamma_HeP_integrand(nu+delta_nu/2) * delta_nu;
    epsilonHg += gamma_H_integrand(nu+delta_nu/2) * delta_nu * (-13.6+nu);
    epsilonHeg += gamma_He_integrand(nu+delta_nu/2) * delta_nu *(-24.6+nu);
    epsilonHePg += gamma_HeP_integrand(nu+delta_nu/2) * delta_nu *(-54.4+nu);
    
  }
  // gammaHg = gammaHeg = gammaHePg = 0;
  // cerr << "done!" << endl;
  // 
}

double Ionise::electronFraction(Particle *p) {
  double nH0, nHe0, nHp, nHep, nHepp, ne;
  calculateFractions(p, nH0, nHe0, nHp, nHep, nHepp, ne);
  return ne/(nH0+nHp);
}

double Ionise::electronDensity(Particle *p) {
  double nH0, nHe0, nHp, nHep, nHepp, ne;
  calculateFractions(p, nH0, nHe0, nHp, nHep, nHepp, ne);
  return ne;
}

void Ionise::calculateFractions(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne) {
  return calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne,gammaHg,gammaHeg,gammaHePg);
}

void Ionise::calculateFractions(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg, double local_gammaHeg, double local_gammaHePg) {

  
  
  if((flags & useTemp)!=0) {
    
    // assume temperature is correct
    calculateFractionsAtT(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);

  } else {
    // assumes internal energy, bot not necessarily temperature, is
    // correct
    
    
    int maxIterations = 20;
    ne = p->ne;
    
    nH0 = 1;
    nHe0 = Y/(4*(1-Y));
    nHp = nHep = nHepp = 0;

    int n=0;
    float xtemp;
  
    do {
      xtemp = p->temp;
      double MeanWeight = 1/(1-0.75*Y+(ne/(nH0+nHp))*(1-Y))*constants::protonMassInG;
      
      p->temp = (MeanWeight/constants::boltzmannInErgPerK) * (2.0/3) * p->u * enToErgsPerG;

      calculateFractionsAtT(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);
      n++;
    } while(n<maxIterations && abs((p->temp-xtemp)/(p->temp))>0.01);
    
    

  }
    
}


void Ionise::calculateFractionsAtT(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg, double local_gammaHeg, double local_gammaHePg) {
  
  // optically thin photoionisation equilibrium calculation
  // based on Katz, Weinberg, Hernquist 1996
  //
  // performed for particle n
  //
  // see also cooling.c in non-public version of Gadget
 
  // constants
  const double logTmin = 0.0;
  const double logTmax = 9.0;
  
  const int max_iterations = 1000;

  double temp = p->temp;
 
  //temp*=1.5;
  double logT = log(temp)/log(10.);

  // rho is in units of Msol/kpc^3.
  double nH = ((double)p->rho*denToProtonsPerCm3)*(1-Y);
  

  ne=nH/100000000;

  if(logT<=logTmin) {
    // expect totally neutral

    nH0 = nH;
    nHe0 = nH*Y/(4-4*Y);
    nHp = 0;
    nHep = 0;
    nHepp = 0;
    ne = 0;
 
  }

  if(logT>=logTmax) {
    // must be fully ionised
    nH0 = 0;
    nHe0 = 0;
    nHp = nH;
    nHep = 0;
    nHepp = nH*Y/(4-4*Y);
    ne=nHp+2*nHepp;
    return;
  }


  // Recombination Rates

  
  // without OTS approximation:
  // double aHp = 3.34e-10 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7));

  // with OTS:
  double lam = 2.*157807./temp;
  double aHp = 2.753e-14 * pow(lam,1.5)/pow(1.+pow(lam/2.740,0.407),2.242);

  double aHep = 1.5e-10 * pow(temp,-0.6353);
  double aHepp = 1.34e-9 * pow(temp,-0.7)/(1.+pow(temp/1.e6,0.7)); 
  double ad = 1.9e-3 * pow(temp,-1.5) *exp(-4.7e5/temp) * (1+0.3*exp(-9.4e4/temp));
  

  /*
  double aHp = gasoline::clRateRadrHII(temp);
  double aHep = gasoline::clRateRadrHeII(temp);
  double aHepp = gasoline::clRateRadrHeIII(temp);
  double ad = gasoline::clRateDielHeII(temp);
  */
 

  // Collisonal Ionisation Rates

  
  double gammaH0 =  5.85e-11 * pow(temp,0.5) * exp(-157809.1/temp) / (1+pow(temp/1.e5,.5));
  double gammaHe0 = 2.38e-11 * pow(temp,0.5) * exp(-285335.4/temp) / (1+pow(temp/1.e5,.5));
  double gammaHeP =  5.68e-12 * pow(temp,0.5) * exp(-631515.0/temp) / (1+pow(temp/1.e5,.5));


  /*
  double gammaH0 = gasoline::clRateCollHI(temp);
  double gammaHe0 = gasoline::clRateCollHeI(temp);
  double gammaHeP = gasoline::clRateCollHeII(temp);
  */

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
  // cerr << n_iterations << gammaH0*ne*nH0 + local_gammaHg*nH0 << " " << aHp*ne*ne << endl;
}

void Ionise::setFlag(const int flag) {
  if((flags & flag)==0)
    flags+=flag;
}

void Ionise::unsetFlag(const int flag) {
  if((flags & flag)!=0)
    flags-=flag;
}

void Ionise::thickPPRadiative(SimSnap *pSim) {
  
  if((flags & initGammas)!=0) {
    cerr << "Ionise: cannot call thickRadiative; spectrum unconstrained";
    return;
  }

  int numPart = pSim->getNumParticles();

  float *nPhot = (float*) malloc(sizeof(float) * numPart * nu_quant );
  float *delta_nPhot = (float*) malloc(sizeof(float) * numPart * nu_quant );

  cerr << "Ionise::thickPPRadiative: initialising...";

  double nu = nu_min;
  for(int nu_index=0; nu_index<nu_quant; nu_index++) {
    double delta_nu = (nu_ratio-1.)*nu;
    double N_Nu = 4*PI*J(nu+delta_nu/2)/(constants::planckInErgS21*nu);	    
    for(int pn = 0; pn<numPart; pn++) {
      nPhot[pn*nu_quant+nu_index] = N_Nu;
    }
    nu*=nu_ratio;
  }

  cerr << "done!" << endl;

  int iter = 1;

  do {
    cerr << "Ionise::thickPPRadiative: iteration " << iter << endl;
    cerr << "Ionise::thickPPRadiative:   calculating local gammas & electron fractions...";
    for(int pn=0; pn<numPart; pn++) {
      
      double nu=nu_min;
      for(int nu_index = 0; nu_index<nu_quant; nu_index++) { 
	double delta_nu = (nu_ratio-1.)*nu;
	gammaHg += nPhot[nu_index] * xsec_H(nu+delta_nu/2) * delta_nu;
	gammaHeg += nPhot[nu_index] * xsec_He(nu+delta_nu/2) * delta_nu;
	gammaHePg += nPhot[nu_index] * xsec_HeP(nu+delta_nu/2) * delta_nu;
	nu*=nu_ratio;
      }


      Particle *p = pSim->getParticle(pn);

      double nH0, nHe0, nHp, nHep, nHepp, ne;
     
      // TODO: this should calculate correct fractions for given U,
      // and return sensible temperature

      calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne);
      

      
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
    
    cerr << "done!" << endl << "Ionise::thickPPRadiative:    propogating delta_nPhot...";

    double convRatio = ((pSim->getMassUnits()/pSim->getDensityUnits())/pow(pSim->getDistanceUnits(),2)).convertTo(Unit("cm"),pSim);

    double delta_conv = 0;

    nu=nu_min;
    for(int nu_index=0; nu_index<nu_quant; nu_index++) {
      cerr << nu_index;
      
      
      double N_nu = 4*PI*J(nu)/(constants::planckInErgS21*nu);	    
      
      for(int pn1=numPart-1; pn1!=0; --pn1) { // efficient for!
	
	Particle *p1 = pSim->getParticle(pn1);
	
	nPhot[nu_index+pn1*nu_quant] = N_nu;
	
	for(int pn2=numPart-1; pn2!=0; --pn2) {
	
	  if(pn1!=pn2) {
	    
	    Particle *p2 = pSim->getParticle(pn2);
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
	  // cerr << "Ionise::thickPPRadiative:   nPhot is negative (" << nPhot[nu_index+pn1*nu_quant] << " vs " << N_nu << " for particle " << pn1 << endl;
	  nPhot[nu_index+pn1*nu_quant]=0;
	}

	delta_conv += ((N_nu-nPhot[nu_index+pn1*nu_quant])/N_nu)/nu_quant;
      } // for pn1 (=particle to update)
      
      nu*=nu_ratio;
    } // for nu_index
    cerr << "Ionise::thickPPRadiative: end of iteration " << iter << "; delta_conv = " << delta_conv << endl;
  } while(0==0);

}

void Ionise::thinRadiative() {
  int numPart=  pSim->getNumParticles();
  SimanArray & nHII = pSim->createArray("nHII","nHII");
  for(int n=0; n<numPart; n++) {
    Particle *p = pSim->getParticle(n);
    if(p->type==Particle::gas) {
      double nH0, nHe0, nHp, nHep, nHepp, ne;
      calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne);
      p->ne = ne/(nH0+nHp);
      nHII[n] = nHp/(nH0+nHp);
    }
  }
}

void Ionise::haehneltRadiative() {
  int numPart=  pSim->getNumParticles();
  SimanArray & nHII = pSim->createArray("nHII","nHII");
  double crit = 1.e-2 * Unit("m_p cm^-3").convertTo(pSim->getDensityUnits());
  cerr << "limit is " << crit;
  for(int n=0; n<numPart; n++) {
    Particle *p = pSim->getParticle(n);
    if(p->type==Particle::gas) {
      double nH0, nHe0, nHp, nHep, nHepp, ne;
      if(p->rho>crit) 
	calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne,0,0,0);
      else
	calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne);

      p->ne = ne/(nH0+nHp);
      nHII[n] = nHp/(nH0+nHp);
    }
  }
}

float Ionise::quickNegExp(float ex) {

  const double res[480] = {1.0,0.951229424501,0.904837418036,0.860707976425,0.818730753078,0.778800783071,0.740818220682,0.704688089719,0.670320046036,0.637628151622,0.606530659713,0.57694981038,0.548811636094,0.522045776761,0.496585303791,0.472366552741,0.449328964117,0.427414931949,0.406569659741,0.386741023455,0.367879441171,0.349937749111,0.332871083698,0.316636769379,0.301194211912,0.28650479686,0.272531793034,0.259240260646,0.246596963942,0.234570288094,0.223130160148,0.212247973827,0.201896517995,0.192049908621,0.182683524053,0.17377394345,0.165298888222,0.157237166314,0.149568619223,0.142274071587,0.135335283237,0.128734903588,0.122456428253,0.116484157773,0.110803158362,0.105399224562,0.100258843723,0.0953691622155,0.0907179532894,0.0862935864994,0.0820849986239,0.0780816660012,0.0742735782143,0.0706512130604,0.0672055127397,0.0639278612067,0.0608100626252,0.0578443208748,0.0550232200564,0.0523397059484,0.0497870683679,0.0473589243911,0.0450492023936,0.042852126867,0.0407622039784,0.0387742078317,0.0368831674012,0.0350843541008,0.0333732699603,0.0317456363781,0.0301973834223,0.0287246396542,0.0273237224473,0.0259911287788,0.0247235264703,0.023517745856,0.0223707718562,0.0212797364384,0.0202419114458,0.0192547017754,0.0183156388887,0.0174223746395,0.0165726754018,0.0157644164849,0.0149955768205,0.014264233909,0.0135685590122,0.0129068125805,0.0122773399031,0.0116785669704,0.0111089965382,0.0105672043839,0.0100518357446,0.00956160193054,0.0090952771017,0.00865169520312,0.00822974704902,0.00782837754923,0.00744658307092,0.00708340892905,0.00673794699909,0.00640933344626,0.00609674656552,0.00579940472684,0.00551656442076,0.00524751839918,0.00499159390691,0.00474815099941,0.00451658094261,0.00429630469075,0.00408677143846,0.00388745724348,0.00369786371648,0.00351751677491,0.00334596545747,0.00318278079651,0.00302755474538,0.00287989915809,0.00273944481877,0.00260584051841,0.00247875217667,0.00235786200649,0.00224286771949,0.00213348177004,0.0020294306363,0.00193045413623,0.00183630477703,0.00174674713626,0.00166155727317,0.00158052216874,0.00150343919298,0.00143011559831,0.00136036803755,0.00129402210547,0.00123091190267,0.00117087962079,0.00111377514784,0.00105945569291,0.00100778542905,0.000958635153694,0.000911881965555,0.000867408957307,0.000825104923266,0.000784864081311,0.000746585808377,0.000710174388843,0.000675538775194,0.000642592360356,0.00061125276113,0.000581441612194,0.000553084370148,0.000526110127116,0.000500451433441,0.000476044129022,0.000452827182887,0.000430742540576,0.00040973497898,0.000389751968253,0.000370743540459,0.000352662164628,0.000335462627903,0.000319101922481,0.000303539138079,0.000288735359628,0.000274653569972,0.000261258557302,0.000248516827108,0.000236396518429,0.000224867324179,0.000213900415368,0.000203468369011,0.000193545099558,0.000184105793668,0.000175126848158,0.000166585810988,0.000158461325116,0.000150733075095,0.000143381736276,0.000136388926482,0.000129737160046,0.000123409804087,0.000117391036919,0.00011166580849,0.000106219802746,0.000101039401837,9.61116520614e-05,9.14242314782e-05,8.69654190944e-05,8.27240655566e-05,7.86895652718e-05,7.48518298877e-05,7.12012630669e-05,6.77287364909e-05,6.44255670344e-05,6.12834950532e-05,5.82946637309e-05,5.54515994322e-05,5.27471930155e-05,5.01746820562e-05,4.77276339368e-05,4.53999297625e-05,4.31857490603e-05,4.10795552253e-05,3.90760816757e-05,3.71703186841e-05,3.53575008504e-05,3.36330951857e-05,3.19927897777e-05,3.04324830084e-05,2.89482732982e-05,2.75364493497e-05,2.61934808678e-05,2.49160097315e-05,2.37008415978e-05,2.25449379132e-05,2.14454083166e-05,2.03995034112e-05,1.94046078899e-05,1.84582339958e-05,1.75580153011e-05,1.67017007902e-05,1.58871492309e-05,1.51123238199e-05,1.437528709e-05,1.36741960657e-05,1.30072976541e-05,1.23729242618e-05,1.17694896249e-05,1.11954848426e-05,1.06494746038e-05,1.01300935986e-05,9.63604310396e-06,9.16608773625e-06,8.71905236227e-06,8.29381916076e-06,7.8893248272e-06,7.50455791508e-06,7.13855630669e-06,6.79040480738e-06,6.45923285705e-06,6.14421235333e-06,5.84455558087e-06,5.55951324165e-06,5.28837258136e-06,5.03045560711e-06,4.78511739213e-06,4.55174446308e-06,4.32975326609e-06,4.11858870754e-06,3.91772276602e-06,3.72665317208e-06,3.54490215219e-06,3.37201523414e-06,3.20756011058e-06,3.05112555804e-06,2.90232040865e-06,2.76077257204e-06,2.62612810488e-06,2.49805032587e-06,2.37621897385e-06,2.26032940698e-06,2.15009184098e-06,2.04523062452e-06,1.94548354994e-06,1.85060119758e-06,1.76034631216e-06,1.67449320943e-06,1.59282721194e-06,1.51514411214e-06,1.44124966183e-06,1.37095908638e-06,1.30409662276e-06,1.24049507996e-06,1.179995421e-06,1.12244636523e-06,1.06770401003e-06,1.015631471e-06,9.66098539667e-07,9.18981357898e-07,8.741621082e-07,8.31528719104e-07,7.90974584929e-07,7.52398299216e-07,7.15703401159e-07,6.80798134398e-07,6.47595217584e-07,6.16011626132e-07,5.85968384611e-07,5.57390369269e-07,5.30206120182e-07,5.04347662568e-07,4.79750336813e-07,4.5635263679e-07,4.34096056064e-07,4.12924941587e-07,3.92786354548e-07,3.73629937989e-07,3.55407790889e-07,3.3807434839e-07,3.21586267858e-07,3.05902320502e-07,2.90983288284e-07,2.76791865854e-07,2.63292567263e-07,2.50451637233e-07,2.3823696675e-07,2.26618012777e-07,2.15565721875e-07,2.05052457561e-07,1.95051931198e-07,1.85539136262e-07,1.76490285808e-07,1.67882753e-07,1.59695014519e-07,1.51906596757e-07,1.44498024611e-07,1.37450772792e-07,1.307472195e-07,1.2437060236e-07,1.18304976508e-07,1.12535174719e-07,1.07046769484e-07,1.01826036931e-07,9.68599225093e-08,9.21360083457e-08,8.76424821944e-08,8.33681078996e-08,7.93021972991e-08,7.54345834984e-08,7.17555954487e-08,6.82560337633e-08,6.49271477154e-08,6.17606133558e-08,5.87485126993e-08,5.58833139252e-08,5.31578525442e-08,5.05653134834e-08,4.80992140445e-08,4.57533876945e-08,4.35219686456e-08,4.13993771879e-08,3.93803057371e-08,3.7459705563e-08,3.56327741646e-08,3.3894943262e-08,3.22418673726e-08,3.06694129456e-08,2.91736480261e-08,2.77508324224e-08,2.63974083546e-08,2.51099915574e-08,2.38853628184e-08,2.27204599277e-08,2.16123700215e-08,2.05583222976e-08,1.95556810879e-08,1.86019392669e-08,1.76947119835e-08,1.68317306967e-08,1.6010837504e-08,1.52299797447e-08,1.44872048677e-08,1.37806555489e-08,1.31085650471e-08,1.24692527858e-08,1.18611201513e-08,1.12826464955e-08,1.07323853328e-08,1.02089607236e-08,9.71106383386e-09,9.23744966197e-09,8.78693392581e-09,8.35839010137e-09,7.95074660588e-09,7.56298411827e-09,7.19413303033e-09,6.84327102222e-09,6.50952075617e-09,6.19204768266e-09,5.89005795366e-09,5.60279643754e-09,5.32954483087e-09,5.06961986232e-09,4.82237158407e-09,4.58718174665e-09,4.36346225294e-09,4.1506536877e-09,3.94822391865e-09,3.75566676594e-09,3.57250073638e-09,3.3982678195e-09,3.23253234224e-09,3.07487987959e-09,2.92491621827e-09,2.78226637102e-09,2.64657363891e-09,2.51749871944e-09,2.39471885807e-09,2.27792704121e-09,2.16683122846e-09,2.06115362244e-09,1.96062997408e-09,1.8650089219e-09,1.77405136347e-09,1.68752985751e-09,1.60522805519e-09,1.52694015913e-09,1.45247040881e-09,1.38163259108e-09,1.31424957448e-09,1.25015286639e-09,1.18918219163e-09,1.13118509177e-09,1.07601654385e-09,1.02353859776e-09,9.73620031301e-10,9.26136022057e-10,8.8096783527e-10,8.38002526948e-10,7.97132661439e-10,7.58256042791e-10,7.21275459208e-10,6.86098439969e-10,6.52637024203e-10,6.2080754094e-10,5.90530399894e-10,5.61729892442e-10,5.34334002312e-10,5.08274225511e-10,4.83485399021e-10,4.59905537865e-10,4.37475680108e-10,4.16139739422e-10,3.95844364843e-10,3.76538807361e-10,3.58174793028e-10,3.40706402243e-10,3.24089954929e-10,3.08283901314e-10,2.9324871803e-10,2.78946809287e-10,2.65342412864e-10,2.52401510685e-10,2.40091743752e-10,2.28382331236e-10,2.17243993508e-10,2.06648878921e-10,1.9657049417e-10,1.86983638043e-10,1.77864338406e-10,1.69189792262e-10,1.60938308724e-10,1.53089254788e-10,1.45623003729e-10,1.38520886031e-10,1.31765142701e-10,1.25338880861e-10,1.19226031509e-10,1.13411309337e-10,1.07880174513e-10,1.02618796317e-10,9.76140185636e-11,9.28533267014e-11,8.83248165212e-11,8.40171643886e-11,7.99195989295e-11,7.60218740961e-11,7.23142435459e-11,6.87874362713e-11,6.54326334173e-11,6.22414462291e-11,5.92058950766e-11,5.63183895007e-11,5.35717092336e-11,5.09589861438e-11,4.84736870627e-11,4.61095974481e-11,4.38608058445e-11,4.17216891016e-11,3.96868983133e-11};

  int n=(int) (-ex/0.05);
  if(n>479 || n<0) return 0.;
  
  return (float) res[n];
  
}

#define INDEX(Ax,Ay,Az,ADir,AQuant) Ax+nx*(Ay+ny*(Az+nz*(AQuant+nu_quant*ADir)))

void Ionise::propogateMeanTaus(Grid *pGrid) {
  for(int nu=0; nu<nu_quant; nu++) {
    pGrid->tau[nu]=0;
  }

  int num_cells = pGrid->getNx()*pGrid->getNy()*pGrid->getNz();

  for(int x=0; x<pGrid->getNx(); x++) {
    for(int y=0; y<pGrid->getNy(); y++) {
      for(int z=0; z<pGrid->getNz(); z++) {
	SimSnap *consider = &((*pGrid)[x][y][z]);
	if(typeid(*consider)==typeid(Grid))
	  propogateMeanTaus(static_cast<Grid*>(consider));
	for(int nu=0; nu<nu_quant; nu++) {
	  pGrid->tau[nu]+=(consider->tau[nu])/num_cells;

	}

      }
    }
  }
	
}

void Ionise::thickRadiative(Grid *pGrid) {

  if((flags & initGammas)!=0) {
    cerr << "Ionise: cannot call thickRadiative; spectrum unconstrained";
    return;
  }

  // flags are set once a warning has been given once, to prevent it happening again


  const int xp = 0;
  const int xm = 1;
  const int yp = 2;
  const int ym = 3;
  const int zp = 4;
  const int zm = 5;
  const int ndir = 6;

  float *external_tau = new float[ndir*nu_quant+ndir];

 

  SimanArray* gammaArray, *epsilonArray, *tauArray, *gammaHeIArray, *gammaHeIIArray, *epsilonHeIArray, *epsilonHeIIArray,  *nHIIArray, *numDenArray, *nHeIIArray, *nHeIIIArray;

  nHIIArray = &(pGrid->createArray("nHII","num HII per (HI+HII)"));
  nHeIIArray = &(pGrid->createArray("nHeII","num HeII per (HI+HII)"));
  nHeIIIArray = &(pGrid->createArray("nHeIII","num HeIII per (HI+HII)"));

  if((flags & writeGam)!=0) {
    Unit gamma_units = 1.e-15/Unit("s");
    Unit epsilon_units = Unit("eV") * gamma_units;
    gammaArray = &(pGrid->createArray("gamHIg","Gamma_(HI+photon)",gamma_units));
    gammaHeIArray = &(pGrid->createArray("gamHeIg","Gamma_(HeI+photon)",gamma_units));
    gammaHeIIArray = &(pGrid->createArray("gamHeIIg","Gamma_(HeII+photon)",gamma_units));
    epsilonArray = &(pGrid->createArray("epsHIg","Epsilon_(HI+photon)",epsilon_units));
    epsilonHeIArray = &(pGrid->createArray("epsHeIg","Epsilon_(HeI+photon)",epsilon_units));
    epsilonHeIIArray = &(pGrid->createArray("epsHeIIg","Epsilon_(HeII+photon)",epsilon_units));
    tauArray = &(pGrid->createArray("tau","tau",Unit("kpc^-1")));
    numDenArray = &(pGrid->createArray("numDen","total number density"));
  } 

  float *cur_n_phot = new float[nu_quant];

  int iter =0;
  int maxiter  =40;
  int miniter = 20; // just a guess for now

  double maxdelta=0;

  pGrid->tau = new float[nu_quant];

  Grid::walkIterator i;
  
  for(i=pGrid->walkBegin(); i!=pGrid->walkEnd(); ++i) {
    
    (*i).tau = new float[nu_quant];
    for(int n=0; n<nu_quant; n++)
      (*i).tau[n]=0;
  }

  cerr << "Ionise: caching element cross-sections...";
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
  
  double max_local_tau = 0;
  double max_change_tau = 0;

  do {
    
    max_change_tau = 0;
    max_local_tau =  0;

    if((flags & verbose)!=0)
      cerr << "Ionise::thickRadiative: Iteration " << iter << endl;
    
    // write-out backup

     
    std::ostringstream filename;
    filename << "ionise.backup"; // << iter;
    
    pGrid->write(filename.str(),SimSnap::native);

    Grid::baseWalkIterator wi;

    
    #ifdef SIMAN_OMP
    #pragma omp parallel for 
    #endif 
    for(wi=pGrid->baseWalkBegin(); wi!=pGrid->baseWalkEnd(); ++wi) {
    	
      // now do equilibrium calculation for each particle in this cell


      SimSnap *thisCell = &(*wi);

      Grid *immediateParent = &(wi.parent());
      coordinate<int> location = wi.location();

      int x = location.x, y = location.y, z=location.z;

      int num_parts = thisCell->getNumParticles();

      

      
      double dx = immediateParent->getDx();
      double dy = immediateParent->getDy();
      double dz = immediateParent->getDz();
      
      double x1 = immediateParent->getX1();
      double y1 = immediateParent->getY1();
      double z1 = immediateParent->getZ1();
      
      
      double ne_mean=0, nh0_mean=0, nhe0_mean=0, nhep_mean=0;
      double local_gammaHg=0, local_gammaHeg=0, local_gammaHePg=0;
      double local_epsilonHg=0, local_epsilonHeg=0, local_epsilonHePg=0;
      
      
      float* our_tau = thisCell->tau;
      
      float total_tau[ndir];


      for(int dir=0; dir<ndir; dir++) {
	for(int nu_index=0; nu_index<nu_quant; nu_index++) {
	  external_tau[dir+ndir*nu_index]=0;
	  
	  if((flags & (64<<dir))!=0)
	    external_tau[dir+ndir*nu_index]=numeric_limits<float>::infinity();
	}
	
	// Now note that we need to iterate in the OPPOSITE direction
	// to our own "direction" index which gives us the direction of the RADIATION flux we're interested in

	int iterate_dir = 0;
	switch(dir) {
	case xm:
	  iterate_dir=xp;
	  break;
	case xp:
	  iterate_dir=xm;
	  break;
	case ym:
	  iterate_dir=yp;
	  break;
	case yp:
	  iterate_dir=ym;
	  break;
	case zm:
	  iterate_dir=zp;
	  break;
	case zp:
	  iterate_dir=zm;
	  break;
	}

	for(Grid::directionalIterator di = ++(wi.directionalBegin(iterate_dir)); di!=wi.directionalEnd(); ++di) {

	  for(int nu_index=0; nu_index<nu_quant; nu_index++) {
	    switch(dir) {
	    case xm: case xp:
	      external_tau[dir+ndir*nu_index]+=(*di).tau[nu_index]*di.parent().getDx();
	      break;
	    case ym: case yp:
	      external_tau[dir+ndir*nu_index]+=(*di).tau[nu_index]*di.parent().getDy();
	      break;
	    case zm: case zp:
	      external_tau[dir+ndir*nu_index]+=(*di).tau[nu_index]*di.parent().getDz();
	      break;
	    } // switch dir
	  } // for nu_index
	
	} // for di

      } // for dir

      for(int n=0; n<num_parts; n++) {

	int n_de =  thisCell->deReference(n,pGrid->getParent()); // just for writing out feedback information, not used by physical calculations

	local_gammaHg=local_gammaHeg=local_gammaHePg=0;
	
	Particle *p = thisCell->getParticle(n);
	
	
	
	
	// Based on tau being constant locally, calculate local photon number
	// in all wavebands, and use that to calculate the photoionisation gammas
	
	double nu = nu_min;
	
	for(int nu_index = 0; nu_index<nu_quant; nu_index++) {
	  
	  
	  	  
	  // contribution within this cell
	  
	  total_tau[xp] = external_tau[xp+ndir*nu_index] + (p->x-(x*dx+x1))*our_tau[nu_index];
	  total_tau[xm] = external_tau[xm+ndir*nu_index] + (((x+1)*dx+x1)-p->x)*our_tau[nu_index];
	  total_tau[yp] = external_tau[yp+ndir*nu_index] + (p->y-(y*dy+y1))*our_tau[nu_index];
	  total_tau[ym] = external_tau[ym+ndir*nu_index] + (((y+1)*dy+y1)-p->y)*our_tau[nu_index];
	  total_tau[zp] = external_tau[zp+ndir*nu_index] + (p->z-(z*dz+z1))*our_tau[nu_index];
	  total_tau[zm] = external_tau[zm+ndir*nu_index] + (((z+1)*dz+z1)-p->z)*our_tau[nu_index];
	  
		  
	  // again, the following will require generalisation if n_dir>6
	  
	  float mul = 0;
	  
	  for(int dir=0; dir<ndir; dir++) {
	    // proper exp is too slow	
	    mul+=quickNegExp(-total_tau[dir]);
	  } // for dir
	  

	  cur_n_phot[nu_index]=0.66666*PI*J(nu_index)/(constants::planckInErgS21*nu)*mul; // constant 0.66 = 4/ndir
	  
	  local_gammaHg += cur_n_phot[nu_index] * xsec_H_dnu[nu_index];
	  local_gammaHeg += cur_n_phot[nu_index] * xsec_He_dnu[nu_index];
	  local_gammaHePg += cur_n_phot[nu_index] * xsec_HeP_dnu[nu_index];

	  local_epsilonHg += (nu - 13.6) * cur_n_phot[nu_index] * xsec_H_dnu[nu_index];
	  local_epsilonHeg += (nu - 24.6) * cur_n_phot[nu_index] * xsec_He_dnu[nu_index];
	  local_epsilonHePg += (nu - 54.4) * cur_n_phot[nu_index] * xsec_HeP_dnu[nu_index];
	  
	  
	  nu*=nu_ratio;
	  
	  if(tauArray!=NULL && nu_index==0) {
	    (*tauArray)[n_de] = mul; // total_tau[xp];
	  }

	} // for nu_index
	
	
	// STEP 3.
	// Calculate fractions
	
	double nH0, nHe0, nHp, nHep, nHepp, ne;
	
	// TODO: this should calculate correct fractions for given U,
	// and return sensible temperature
	calculateFractions(p,nH0,nHe0,nHp,nHep,nHepp,ne,local_gammaHg,local_gammaHeg,local_gammaHePg);

	
	ne_mean+=ne;
	nh0_mean+=nH0;
	nhe0_mean+=nHe0;
	nhep_mean+=nHep;
	

	p->ne = ne/(nH0+nHp); 
	//p->nHp = nHp/(nH0+nHp);
	

	if(gammaArray!=NULL) {
	  // write out gamma H0 for information
	 
	  (*gammaArray)[n_de]=(float)(local_gammaHg*(double)1.e15);
	  (*gammaHeIArray)[n_de] = (float)(local_gammaHeg*(double)1.e15);
	  (*gammaHeIIArray)[n_de] = (float)(local_gammaHePg*(double)1.e15);

	 
	  (*epsilonArray)[n_de]=(float)(local_epsilonHg*(double)1.e15);
	  (*epsilonHeIArray)[n_de] = (float)(local_epsilonHeg*(double)1.e15);
	  (*epsilonHeIIArray)[n_de] = (float)(local_epsilonHePg*(double)1.e15);



	  (*numDenArray)[n_de] = nH0+nHe0+nHp+nHep+nHepp;
	}
   
	(*nHIIArray)[n_de] = nHp/(nH0+nHp);
	(*nHeIIArray)[n_de] = nHep/(nH0+nHp);
	(*nHeIIIArray)[n_de] = nHepp/(nH0+nHp);
      }  // for n (particles in this cell)
      
      if(num_parts>0) {
	ne_mean/=num_parts;
	nh0_mean/=num_parts;
	nhe0_mean/=num_parts;
	nhep_mean/=num_parts;
      }

      // now calculate output fluxes using inferred optical depth
      
      
      double nu = nu_min;
      // cerr << dx << " " << xsec_H(nu_min+(nu_ratio-1.)*nu_min*0.5) << " " <<  max_local_tau << " " << nh0_mean <<" " << nhe0_mean <<" " <<  nhep_mean << endl;
      for(int nu_index = 0; nu_index<nu_quant; nu_index++) { 
	
	double delta_nu = (nu_ratio-1.)*nu;
	
	double local_tau = (nh0_mean * xsec_H(nu+delta_nu*0.5) + nhe0_mean * xsec_He(nu+delta_nu*0.5) +
			    nhep_mean * xsec_HeP(nu+delta_nu*0.5)) * constants::KpcInCm;
	
	


	double tau_change = abs(thisCell->tau[nu_index]-local_tau);
	
	// local_tau = number of absorptions in this waveband per kpc per photon
	
	//tau[INDEX(x,y,z,0,nu_index)] = local_tau;
	thisCell->tau[nu_index] = local_tau;

	if(local_tau*dx  > max_local_tau)
	   max_local_tau = local_tau*dx;
	
	if(tau_change*dx > max_change_tau)
	  max_change_tau = tau_change*dx;
	
	nu = nu * nu_ratio;
	
      } // for nu_index (updating all directions of output flux)
      
    } // for wi (walk base iterator)


    cerr << "propogateMeanTaus...";
    propogateMeanTaus(pGrid);
    cerr << "done!" << endl;

    cout << "max_local_tau = " << max_local_tau << " ( maximum optical depth of any cell in any dir. ) " << endl;
    cout << "max_change_tau = " << max_change_tau << " ( maximum *change* in optical depth of any cell in any dir. ) "<< endl;
    iter++;
  } while(((max_change_tau/max_local_tau)>0.001 || iter<miniter) && iter<maxiter);


  delete[] external_tau;
  delete[] cur_n_phot;
  delete[] xsec_H_dnu;
  delete[] xsec_He_dnu;
  delete[] xsec_HeP_dnu;

  
  for(i=pGrid->walkBegin(); i!=pGrid->walkEnd(); ++i) {
    delete[] (*i).tau;
  }

}

