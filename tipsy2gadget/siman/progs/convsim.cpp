#include <siman/base.hpp>

int main(int argc, char **argv)
{

    
  string path;
 
  unsigned int type=0;
  
  if(argc<3) {
    cerr << "Syntax: convsim [type] [file_list]" << endl;
    exit(0);
  }

  if(strcmp(argv[1],"gadget")==0)
    type = CSimSnap::gadget;

  if(strcmp(argv[1],"tipsy")==0)
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
    pF->setBoxSize(x2-x1);

    cerr << "box="<< x2-x1;

    CUnit lenUnits("(kpc h^-1 a)");
    CUnit massUnits("(1.e10 msol h^-1)");
    CUnit velUnits("(km s^-1 a^1/2)");
    
    CUnit denUnits = massUnits/(lenUnits*lenUnits*lenUnits);
    
    CUnit enUnits("(km^2 s^-2)");
   

    try {
      pF->convertUnits(lenUnits,massUnits,velUnits,denUnits,enUnits);
      pF->write(path+"."+argv[1],type);
    } catch (CUnitsError &e) {
      cerr << e;
    }
    delete pF;
  }
}


