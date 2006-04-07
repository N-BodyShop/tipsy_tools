
#ifdef SIMAN_VIS
#include <siman/base.hpp>
#include <siman/visual.hpp>
#endif

int main(int argc, char **argv)
{

#ifdef SIMAN_VIS
  char path[1024];
  
  if(argc < 2) {
    cerr << "Syntax: vishalo sim_file" << endl;
    cerr << endl << "sim_file = path to halo file" << endl;
    exit(0);
  }

  strcpy(path,argv[1]);

  CSimSnap *pF = CSimSnap::loadFile(path);

  pF->convertUnits();

  cout << pF->getParticle(4)->x << " " << pF->getParticle(8)->y << " " << pF->getParticle(16)->z << endl;

  cerr << pF->getTotalMass() << pF->getMassUnits() << endl;

  if(argc>2) {
    CSimSnap *pF2 = CSimSnap::loadFile(argv[2]);
    if(pF2->getNumParticles()==pF->getNumParticles()) {
      float* dNH0 = pF->createArray("dNH0","delta_NH0");
      for(int n=0; n<pF2->getNumParticles(); n++) {
	dNH0[n] = pow(pF->getParticle(n)->nHp-pF2->getParticle(n)->nHp,2);
      }
    }
    delete pF2;
  }

  float x,y,z;
  float ax,ay,az;
  
  CParticleTypeFilter ptf(CParticle::gas);
  CSubset gas(pF,ptf);
 
  
  for(int n=0; n<gas.getNumParticles(); n++) {
    CParticle *p = gas.getParticle(n);
    
    gas.releaseParticle(p);
    
  }



  CGrid *g = new CGrid(pF,20,20,20);
  CVisualise vis(g);
  
  float *colours = pF->createArray("old_ne","input ne");
  
  for(int n=0;n<pF->getNumParticles();n++) {
    colours[n]=pF->getParticle(n)->ne;
  }
  
  //CIonise ionise(0.1,-1,pF->getDensityUnits(),pF->getEnergyUnits()/pF->getMassUnits(),0.);
  //ionise.setFlag(CIonise::useTemp);
  //  ionise.thinRadiative(pF); // get ne fraction assuming given ionising background
  // pF->UFromTemp();
  // this should have roughly recreated the U, ne values found internally in gasoline
  
  // now see what the output temps/ne would have looked like without a UV background (or shielded)

  pF->UFromTemp();
  //  CIonise ionise(pF);
  //  ionise.setFlag(CIonise::useTemp);
  // ionise.thinRadiative(pF);
  
  float col1[4] = {0,0,1,1};
  float col2[4] = {1,0,0,1};
  float col3[4] = {1,1,0,1};
  float col4[4] = {0,1,1,1};
  
  CColourMapGradient gascm(&vis,col1, col2);
 
  CColourMapGradient starcm(&vis,col3,col4);
  
  //CColourMap starcm(pF);
  //starcm.setReferenceColour(1,0,0,0.1);

  CColourMap dmcm(&vis);
  dmcm.setReferenceColour(0,1,0,0.1);

  //CColourMap gascm(&bla);
  //gascm.setReferenceColour(0,0,1,0.1);

  CColourMapByType cm(&vis,&dmcm,&gascm,&starcm);

  cm.autoRange();

  double mean_temp = 0.;
  int n_t = 0;
  for(int n=0;n<pF->getNumParticles();n++) {
    CParticle *p = pF->getParticle(n);
    if(p->type == CParticle::gas);
    mean_temp+=p->temp/1.e4;
    
    pF->releaseParticle(p);
    n_t++;
  }

  vis.colourise(&cm);
  vis.run(); 

  cout << "Finished!" << endl;

  delete pF;

#endif
}


