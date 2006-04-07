#include <siman/base.hpp>

int main(int argc, char **argv)
{
  CSimSnap *pF=CSimSnap::loadFile(argv[1]);
  float *arrNlim1=pF->getArray("n13.6");
  float *arrNlim2=pF->getArray("n24.6");
  float *arrNlim3=pF->getArray("n54.4");
  double rhoconv=1;
  try {
    rhoconv=pF->getDensityUnits().convertTo((string)"(m_p cm^-3)");
  } catch (CUnitsError &e1) {
    cerr << e1;
  }
  for(int n=0;n<pF->getNumParticles();n++) {
    
    cout << pF->getParticle(n)->rho*rhoconv << "\t" << arrNlim1[n] << "\t" << arrNlim2[n] << "\t" << arrNlim3[n] <<  endl;
  }
}








  



