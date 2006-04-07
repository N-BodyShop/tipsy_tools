// MAKEABS.CPP
//
// Absorption profiles stuff


#ifdef SIMAN_FITS

#include <siman/base.hpp>
#include <siman/extra.hpp>

#include <time.h>

using namespace std;
using namespace CCfits;




clock_t _clock_start = clock();

void startclock(void) {
  _clock_start = clock();
}

void finishclock(void) {
  
  cerr << "(Time: " << (double)(clock()-_clock_start)/((double)CLOCKS_PER_SEC) << ")" << endl;
}

#endif

int main(int argc, char **argv)
{

#ifdef SIMAN_FITS

  /*
  while(true) {
    units::CUnit r1,r2;
    cin >> r1;
    cout << r1 << endl;
    cin >> r2;
    cout << r2 << endl;
    cout << "RAT: " << r1/r2 << endl;
    try {
      double cr = r1.convertTo(r2);
      cout << "CONVRAT: " << cr << endl;
    } catch(CUnitsError &e1) {
      cout << e1 << endl;
    }
  }
  */

  char path[1024]="/home/app26/fabio/MW1.512g11.00512";
  float spsize;
  float imsize;

  if(argc < 4) {
    cerr << "Syntax: makeabs sim_file sphere_size image_size [zeroflux]" << endl;
    exit(0);
    }
  
  spsize=atof(argv[2]);
  imsize=atof(argv[3]);
  

  CSimSnap *pFL = CSimSnap::loadFile(argv[1]);
  
  CGeometry geom(pFL);
  geom.setRotate(0.1,0,0);
  geom.apply();

  CSphere sp_desc(0,0,0,spsize);

  CSubset sphere(pFL,sp_desc);

  CParticleTypeFilter f(CParticle::gas);
  CSubset *pF = new CSubset(&sphere,f);
  
  cout << pF->getNumParticles() << "TC " << endl;
  
  if(argc>4) {
    // zero-flux recalculation of ionisation state
    cerr << endl << ">>> Assuming j=0 <<<" << endl << endl;
    CIonise ionise(0,-1,pF);
    ionise.thinRadiative(pF);
  }
  
  float x1 = -imsize;
  float x2=  imsize;

  int size = 300;

  

  cerr << "COLUMN-GRID FITS OUTPUT" << endl;
  
  ostringstream oss;
  oss << argv[1];

  if(argc>4) oss << ".j0";

  oss << ".HImap." << spsize << "." << imsize << ".fits";

  long dimensions[2] = {size,size};

  FITS fitsFile(oss.str(), FLOAT_IMG, 2, dimensions);

  
  ExtHDU *imageExt;
  vector<long> dim2(2,size);
  
  pF->initialiseSPH();
  for(int n=0; n<31; n++) {
    CGeometry geom(pFL);
    geom.setRotate(3.1415*((float)rand()/(float)RAND_MAX),3.1415*((float)rand()/(float)RAND_MAX),3.1415*((float)rand()/(float)RAND_MAX));
    geom.apply();
    
    CColumnGrid cg(pF,x1,x2,size,x1,x2,size);
    ostringstream oss;
    oss << "rot = " << (float)n/10.;

    /*
    CGrid griddedRef(pF,x1,x2,20,x1,x2,20,x1,x2,20,0,100);
    griddedRef.initialiseSPH();
    */

    imageExt = fitsFile.addImage(oss.str(),FLOAT_IMG,dim2);
    cg.columnDensityImage(imageExt,units::protonsPerCm2,true);

    /*    
    oss << " SPH";
    
    imageExt = fitsFile.addImage(oss.str(),FLOAT_IMG,dim2);
    pF->SPHColumnDensityImage(imageExt,x1,x2,size,x1,x2,size,units::protonsPerCm2,true);
    
    */
  }


  cerr << "ADAPTIVE GRID SPH FITS OUTPUT:" << endl;
  startclock();
  //griddedRef.SPHColumnDensityImage(imageExt,x1,x2,size,y1,y2,size,units::protonsPerCm2/(double)1.e20,true);
  finishclock();  


  /*

  for(float rho=0;rho<100000;rho+=10000) {
    
    CDensityCutFilter denCut(rho);
    CSubset ionGuessSim(&virtualSim,denCut);

    ostringstream title_stream;

    title_stream << "Cut: rho = " << rho << ". Npart = " << ionGuessSim.getNumParticles();


    imageExt = fitsFile.addImage(title_stream.str(),FLOAT_IMG,dim2);
    
    pFiltColGrid = new CColumnGrid(&ionGuessSim,-10,10,size,-10,10,size);
    
    pFiltColGrid->columnDensityImage(imageExt,units::protonsPerCm2/(double)1.e20);
    delete pFiltColGrid;

    // Now render a rotated version

    title_stream << " (rotated ry = 1 rad)";
    CGeometry rotateSim(&ionGuessSim);
    
    rotateSim.setRotateY(1);
    
    pFiltColGrid = new CColumnGrid(&rotateSim,-10,10,size,-10,10,size);
    
    imageExt = fitsFile.addImage(title_stream.str(),FLOAT_IMG,dim2);

    pFiltColGrid->columnDensityImage(imageExt,units::protonsPerCm2/(double)1.e20);
    delete pFiltColGrid;


  }
    



 
  CParticle centre(0,0,0);
  for(int n=0;n<20;n++) {
    
    CParticle *p = virtualSim.getParticle(n);

    float es = pow((double) p->mass/p->rho, (double) 1./3.);
    cerr << "E: " << es << " ";
    
    cerr << "C: " << p->distanceTo(centre) << " ";

    CParticle *neighbour = virtualSim.getParticle(virtualSim.getNearestNeighbour(n));

    float ne = p->distanceTo(*neighbour);
    cerr << "N: " << ne << " R: " << es/ne << endl;

    virtualSim.releaseParticle(neighbour);
    virtualSim.releaseParticle(p);
  }
  */
#endif
}


