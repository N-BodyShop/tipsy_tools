#include <siman/base.hpp>

using namespace std;

int main(int argc, char **argv)
{
  
  char path[1024];
  
  if(argc < 2) {
    cerr << "Syntax: glancesim sim_file" << endl;
    cerr << endl << "sim_file = path to simulation file to inspect." << endl;
    exit(0);
  }

  strcpy(path,argv[1]);

  CSimSnap *pF = CSimSnap::loadFile(path);
  pF->convertUnits();
  
  cerr << "omm0=" << pF->getOmegaM0() << "\toml0 = " << pF->getOmegaLambda0() << "\tbox = " << pF->getBoxSize() << " / apparent " << pF->getApparentBoxSize() << "\tz = " << pF->getRedshift() << endl;
  cerr << endl;
  cerr << "total HI mass: " << pF->getH0Mass() << endl;

  CParticleTypeFilter ptf(CParticle::gas);
  CSubset gas(pF,ptf);
  
  for(int n=0;n<20;n++) {
    CParticle *p = pF->getParticle(n);
    cout << p->vx << " " << p->vy << " " << p->vz << endl;
  }

  delete pF;
}
