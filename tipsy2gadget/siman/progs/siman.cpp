
#include <siman/base.hpp>

CScripted settings("config",true);

int main(int argc, char **argv)
{

  CScripted scr;

  for(int n=1; n<argc; n++) {
    CSimSnap *pF = CSimSnap::loadFile(argv[n]);
    char vname = 'a'+(char)n;
    string vname_s = " ";
    vname_s[0]=vname;
    scr.setNamedVar(vname_s,pF);
    cout << "Loaded " << argv[n] << " into var " << vname_s << endl;
  }


  scr.mainLoop();

}


