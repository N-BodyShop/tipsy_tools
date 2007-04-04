//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __IONISE_H_INCLUDED

#define __IONISE_H_INCLUDED


using namespace siman;

class Ionise {

public:

  
  /** \brief Constructor with background J passed in array

  \param pJ - pointer to J array in units of 10^-21 ergs/cm^2/s/Hz with n_j elements, equally spaced bins in log space
  \param nu_min - frequency of first bin in eV
  \param nu_max - frequency of last bin in eV
  \param n_j - number of bins, i.e. length of array passed in pJ
  \param pSim - pointer to simulation on which to act    
  */
  Ionise(float *pj,float nu_min, float nu_max, int n_j, SimSnap *pSim);

  /** \brief Constructor with background J passed as normalisation and log slope
      
  \param j21 - value of j at the lyman limit (13.6eV)
  \param alpha - log slope
  \param pSim - pointer to simulation on which to act
  */
  Ionise(double j21, double alpha, SimSnap *pSim);


  /** \brief Constructor with background J read in from data file
      
  \param pSim - pointer to simulation on which to act
  \param z - redshift to load background J from; otherwise z will be determined from pSim
  */
  Ionise(SimSnap *pSim, float z=-1);


  /** \brief Constructor with background J passed in form of optically thin photo-ionisation gammas
      
  \param gam_Hg - photoionisation gamma for neutral hydrogen in s^-1
  \param gam_Heg - photoionisation gamma for neutral helium in s^-1
  \param gam_Hepg - photionisation gamma for singly ionised helium in s^-1
  \param pSim - pointer to simulation on which to act

  Note that if a Ionise is constructed in this way, optically thick methods are not available, since
  the spectrum is not constrained
  */
  Ionise(double gam_Hg, double gam_Heg, double gam_Hepg, SimSnap *pSim);

  /** \brief set a flag to control behaviour of Ionise
   */
  void setFlag(const int flag);

  /** \brief unset a flag to control behaviour of Ionise
   */
  void unsetFlag(const int flag);

  
  inline double nuFromIndex(int nu_index); //!< convert to indexed nu
  inline int nuToIndex(double nu); //!< convert to index from value of nu

  double J(double nu); //!< calculate J at frequency nu according to initial specifications
  double J(int nu_index); //!< calculate J at frequency indexed by nu_index

  //@{
  //! x-section/cm^2 to ionisation at frequency nu/eV :
  //! from Osterbrock 1989
  double xsec_H(double nu); 
  double xsec_HeP(double nu);
  double xsec_He(double nu);
  //@}

  double gamma_H_integrand(double nu);
  double gamma_HeP_integrand(double nu);
  double gamma_He_integrand(double nu);
  
  double coolingRate(double nden, double T, double ne, double nHI, double nHII, double nHeI, double nHeII, double nHeIII);
  double heatingRate(double nden, double nHI, double nHeI, double nHeII, double epsHg, double epsHeg, double epsHePg);

  void cooling();
  
  double j21;
  double alpha;
  float nu_min; ///< minimum nu in eV
  float nu_max; ///< maximum nu in eV
  int nu_quant; ///< number of log-spaced bins between nu_min and nu_max
  float nu_ratio; ///< ratio of each bin base level to the next

  float Y; // default helium fraction by mass

  float ne_thin_conv; // convergence parameter for thin, const temp calculations
  
  double denToProtonsPerCm3;
  double enToErgsPerG;

  // photoionisation constants

  double gammaHg;
  double gammaHeg;
  double gammaHePg;
  double epsilonHg;
  double epsilonHeg;
  double epsilonHePg;

  
  // flags:

  int flags;

  void calculateGammas(); //!< called internally to calculate gammas from supplied radiation field
private:
  //@{
  //! these flags are set at initialisation and cannot be changed
  static const int initJAndAlpha = 1;
  static const int initGammas = 2;
  static const int initJArray = 4;
  //@}

  void propogateMeanTaus(Grid* pGrid);

public:
  /// \defgroup Flags
  /// Flags which effect operation of Ionise
  //@{
  /// \defgroup PubFlags Public Available Flags:
  /// public flags which can be changed
  /// @see setFlag()
  /// @see unsetFlag()
  //@{
  static const int periodicBoundsX = 8; //!< radiation field is periodic in x-direction
  static const int periodicBoundsY = 16; //!< radiation field is periodic in y-direction
  static const int periodicBoundsZ = 32; //!< radiation field is periodic in z-direction
  static const int periodicBoundsAll = 8+16+32;

  static const int noIlluminationX0 = 64; //!< diffuse radiation field is not set at low x boundary
  static const int noIlluminationXN = 128; //!< diffuse radiation fields is not set at high x boundary
  static const int noIlluminationY0 = 256;  //!< diffuse radiation field is not set at low y boundary
  static const int noIlluminationYN = 512; //!< diffuse radiation fields is not set at high y boundary
  static const int noIlluminationZ0 = 1024;  //!< diffuse radiation field is not set at low z boundary
  static const int noIlluminationZN = 2048; //!< diffuse radiation fields is not set at high z boundary 
  static const int noIlluminationAll = 64+128+256+512+1024+2048; //!< easier to set J=0!

  static const int verbose = 16; //!< set for more diagnostics

  static const int writeGam = 4096; //!< if set, thickRadiative creates an array called gamH0g and writes H0 photoionisation rate for each particle
  static const int useTemp = 8192; //!< if set, assume temperatures rather than internal energies are accurate

  //@}
  //@}

  double electronFraction(Particle *p);
  double electronDensity(Particle *p);

  void calculateFractions(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne);
  void calculateFractions(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg, double local_gammaHeg , double local_gammaHePg) ;
  

  void thinRadiative();
  void haehneltRadiative();
  void thickRadiative(Grid *pGrid);
  void thickPPRadiative(SimSnap *pSimSnap);

  static float quickNegExp(float ex); //!< returns very rough approximation to exp(ex) for ex<0  
private:
  
  //! stores pointer to J array
  /// @see Ionise()
  float *pJ;


  SimSnap *pSim; ///< stores pointer to simulation on which we are to act


  void calculateFractionsAtT(Particle *p, double &nH0, double &nHe0, double &nHp, double &nHep, double &nHepp, double &ne, double local_gammaHg=0, double local_gammaHeg =0, double local_gammaHePg=0); //!< called internally - iteratively if finding correct temp for given U, or just once if flag useTemp is set

};

#endif // __IONISE_H_INCLUDED
