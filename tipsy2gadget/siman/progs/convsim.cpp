#include <siman/base.hpp>

int main(int argc, char **argv)
{

    
  string path;
 
  unsigned int type=0;
  
  if(argc<3) {
    cerr << "Syntax: convsim [type] [file_list]" << endl
	 << "[type] can be one of gadget, siman, tipsy, or tipsycosmo."
	 << endl << "tipsycosmo converts the units to standard PKDGRAV "
	 << "cosmological units." << endl;
    exit(1);
  }

  if(strcmp(argv[1],"gadget")==0)
    type = CSimSnap::gadget;

  if(strcmp(argv[1],"tipsy")==0)
    type = CSimSnap::tipsy;
  
  if(strcmp(argv[1],"tipsycosmo")==0)
    type = CSimSnap::tipsy;
  
  if(strcmp(argv[1],"siman")==0)
    type = CSimSnap::native;

  if(type==0) {
    cerr << "No match for type " << argv[1] << endl;
    exit(0);
  }

  using namespace units;

  cout << "Enter h=";
  float h;
  cin >> h;


  
  for(int n=2;n<argc;n++) {
    path = argv[n];
    
    CSimSnap *pF = CSimSnap::loadFile(path);
    cout << "Calling setHubble...";
    cout.flush();
    pF->setHubble(h);
    cout << "done!" << endl;
    cout.flush();
    float a =1./(1.+ pF->getRedshift());
    cout << "a=" << a;
    
    float x1,x2,d;

    pF->getExactBoundaries(x1,x2,d,d,d,d);
    //    pF->setBoxSize(x2-x1);

    cerr << "box="<< x2-x1 << endl;

    CUnit lenUnits;
    CUnit massUnits;
    CUnit velUnits;
    CUnit denUnits;
    CUnit enUnits;

    if(strcmp(argv[1], "tipsycosmo") == 0) {
	/*
	 * Scale to standard PKDGRAV Units
	 */
	lenUnits = pF->getDistanceUnits();
	lenUnits /= pF->getBoxSize();
	/*
	 * Start with critical mass in 1 kiloparsec^3 for H_0=100
	 * km/s/Mpc Cosmology
	 */
	massUnits = CUnit("(277.49438 msol)");
	massUnits *= h*h;      	// Scale by h^2
	// Multiply by volume
	massUnits *= pow(lenUnits.convertTo(CUnit("(kpc)"), pF), 3.0);
	velUnits = CUnit("(100 km s^-1 a)");
	/*
	 * Scale by h and make the Hubble velocity across the volume
	 * be sqrt(8 PI/3)
	 */
	velUnits *= h*lenUnits.convertTo(CUnit("(Mpc)"), pF)
	    /sqrt(8.0*M_PI/3.0);
	denUnits = massUnits/(lenUnits*lenUnits*lenUnits);
	}
    else {
	lenUnits = CUnit("(kpc h^-1 a)");
	massUnits = CUnit("(1.e10 msol h^-1)");
	velUnits = CUnit("(km s^-1 a^1/2)");

	denUnits = massUnits/(lenUnits*lenUnits*lenUnits);

	enUnits = CUnit("(km^2 s^-2)");
	}

    try {
      pF->convertUnits(lenUnits,massUnits,velUnits,denUnits,enUnits);
      pF->write(path+"."+argv[1],type);
    } catch (CUnitsError &e) {
      cerr << e;
    }
    delete pF;
  }
}


