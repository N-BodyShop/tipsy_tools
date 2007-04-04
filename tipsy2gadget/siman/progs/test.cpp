// test.cpp - part of SimAn Simulation Analysis Library
//
//
// Copyright (c) Andrew Pontzen 2005, 2006
//
// SimAn is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// SimAn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public Licence for more details.
//
// You should have received a copy of the GNU General Public Licence
// along with SimAn; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA


#define SIMAN_TRACE

#include <siman/base.hpp>

using namespace siman;

int numCellsRecursive(Grid &g, bool count_particles=false, bool base_only = false) {
  int result = 0;
  for(int x=0;x<g.getNx();x++) {
    for(int y=0;y<g.getNy();y++) {
      for(int z=0;z<g.getNz();z++) {
	if(typeid(g[x][y][z])==typeid(Grid))
	  result+=numCellsRecursive(static_cast<Grid&>(g[x][y][z]),count_particles,base_only);
	
	if(typeid(g[x][y][z])!=typeid(Grid) || base_only==false) {
	  if(count_particles)
	    result+=g[x][y][z].getNumParticles();
	  else
	    result++;
	}
      }
    }
  }
  return result;
}

void testDirectionalGridIterators() {
  std::cerr << "Test directional grid iterators...";
  SimSnap * testDat = BaseSimSnap::makeSlab(100000,30,30,30,0);
  SimSnap * testDatF = testDat->subset(Sphere(-10,0,0,5) | Sphere(10,0,0,5) | RandomFilter(0.01));
  Grid *g = new Grid(*testDatF,3,3,3,5,5);
  
  int numCellsIterator =0;
 
  Grid::walkIterator i;
  for(i=g->walkBegin();i!=g->walkEnd() && numCellsIterator<500;++i) {
    numCellsIterator++;
   
    for(unsigned int n=0; n<(*i).getNumParticles(); n++)
      (*i).getParticle(n)->temp=1;
  }
 
  GridDirectionalIterator id; 
  

  for( id = i.directionalBegin(GridDirectionalIterator::xp); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
 
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=5;
    }
  }




  for( id = i.directionalBegin(GridDirectionalIterator::yp); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
    
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=5;
    }
  }


  for( id = i.directionalBegin(GridDirectionalIterator::zp); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=5;
    }
  }



  for( id = i.directionalBegin(GridDirectionalIterator::xm); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
  
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=-5;
    }
  }


  for( id = i.directionalBegin(GridDirectionalIterator::ym); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=-5;
    }
  }


  for( id = i.directionalBegin(GridDirectionalIterator::zm); id!=i.directionalEnd(); ++id) {
    SimSnap & s = *id;
    for(unsigned int n=0; n<s.getNumParticles(); n++) {
      s.getParticle(n)->temp=-5;
    }
  }

 
  delete g;
  delete testDatF;
  delete testDat;

  std::cerr << "OK" << std::endl;
}

void testGridIterators() {
  std::cerr << "Test grid iterators...";
  SimSnap * testDat = BaseSimSnap::makeSlab(10000,10,10,10,10);
  Grid *g = new Grid(*testDat,2,2,2,10,10);
  
  int numCellsIterator =0;
  int numPtclsIterator =0;
  for(Grid::walkIterator i=g->walkBegin();i!=g->walkEnd();++i) {
    numCellsIterator++;
    numPtclsIterator+=(*i).getNumParticles();
  }
 
  if(numCellsIterator!=numCellsRecursive(*g) || numPtclsIterator!=numCellsRecursive(*g,true,false)) {
    std::cerr << "FAILED - " << numCellsIterator << " " << numCellsRecursive(*g) << " - number of traversed cells should be independent of method" << std::endl;
    std::cerr << "         " << numPtclsIterator << " " << numCellsRecursive(*g,true,false) << " - number of exposed particles should be indep of method" << std::endl;
    delete g;
    delete testDat;
    return;
  }

  // now try base walk
  numCellsIterator = numPtclsIterator = 0;

  for(Grid::baseWalkIterator i=g->baseWalkBegin();i!=g->baseWalkEnd();++i) {
    numCellsIterator++;
    numPtclsIterator+=(*i).getNumParticles();
  }
  if(numCellsIterator!=numCellsRecursive(*g,false,true) || numPtclsIterator!=numCellsRecursive(*g,true,true)) {
    std::cerr << "FAILED (base) - " << numCellsIterator << " " << numCellsRecursive(*g) << " - number of base traversed cells should be independent of method" << std::endl;
    std::cerr << "                " << numPtclsIterator << " " << numCellsRecursive(*g,true,true) << " - number of exposed particles should be indep of method" << std::endl;
    delete g;
    delete testDat;
    return;
  }

  delete g;
  delete testDat;

  std::cerr << "OK" << std::endl;
}
void testFilteringArrays() {
  std::cerr << "Test filtering arrays...";
  SimSnap * testDat = BaseSimSnap::makeSlab(100,10,10,10,10);
  
  SimanArray & ar = testDat->createArray("testar","testar");
  for(unsigned int n=0; n<testDat->getNumParticles(); n++) 
    ar[n]=(float) n;

  Subset sub(testDat,ModuloFilter(5,0));
  // sub should contain every 5th particle

  SimanArray &subar(sub.getArray("testar"));
  for(unsigned int n=0; n<sub.getNumParticles(); n++) {
    if(subar[n]!=ar[n*5]) {
      std::cerr << n << "not " << subar[n] << "\t - FAIL " << std::endl;
      delete testDat;
      return;
    }
    subar[n]=0;
    if(ar[n*5]!=0) {
      std::cerr << n << "not zero\t" << "set FAIL" << std::endl;
      delete testDat;
      return;
    }

    subar[n] = (float) n;
    
    
  }

  sub.write("simantest_dat",SimSnap::native);
  
  SimSnap *f = SimSnap::loadFile("simantest_dat");
  SimanArray &reloaded_subar = f->getArray("testar");
  
  for(unsigned int n=0; n<f->getNumParticles(); n++) {
    if(reloaded_subar[n]!=(float) n) {
      std::cerr << n << " " << reloaded_subar[n] << " ; reload array FAIL" << std::endl;
    }
  }
  
  std::cerr << "OK" << std::endl;


  
  delete testDat;
}

void testDensityEstimation(SimSnap *pF) {
  std::cerr << "Test density estimator..." << std::endl;
  Subset sub(pF,ParticleTypeFilter(Particle::gas));
  sub.convertUnits();
  if(sub.getNumParticles()<1000) {
    std::cerr << "FAILED - not enough gas particles in test file" << std::endl;
    return;
  }
  Grid gridGas(sub,10,10,10);
  
  
  float cum_ratio=0;
  for(int n=0; n<1000; n++) {
    std::cerr << n << "/1000\r";
    float ratio = gridGas.estimateDensity(n)/gridGas.getParticle(n)->rho;
    cum_ratio+=ratio/((float)1000);
  }

  std::cerr << "  ratio (estimated/specified) av. 1000 particles: " << cum_ratio << std::endl;
  if(cum_ratio>1.1 || cum_ratio<0.9)
    std::cerr << "ODDITY - ratio is not close to unity" << std::endl;
  std::cerr << "Test finished OK " << std::endl;
}

bool uConv(const Unit &a, const Unit &b, double expect, SimSnap *context=NULL) {
  try {
    double rat = a.convertTo(b,context);
    if(std::abs(rat-expect)>expect/1.e6) {
      std::cerr << "FAILED - " << a << " / " << b << ": " << rat << ", expected " << expect << std::endl;
      return false;
    }
    return true;
  } catch(UnitsError &e) {
    std::cerr << "FAILED - " << e <<std::endl;
    return false;
  }
}

void testUnitConversion() {
  std::cerr << "Test unit conversion...";
  
  Unit a("km s^-1");
  Unit b("kpc Msol^-1");
  
  try {

    double rat = a.convertTo(b);
    std::cerr << "FAILED - allowed conversion between incompatible units: " << rat <<  std::endl;
    return;
  } catch(UnitsError &e) {
    
  }
  
  // comparisons derived from google calculator :)

  if (!uConv(Unit("s"),Unit("yr"),3.16887646e-8)) return;
  if (!uConv(Unit("m"),Unit("kpc"),3.24077649e-20)) return;
  if (!uConv(Unit("kpc^3 Msol^-1 yr^-5"),Unit("m^3 kg^-1 s^-5"),4.72022e-10)) return;
  

  // cosmological corrections test

  std::auto_ptr<SimSnap> s(BaseSimSnap::makeSlab(100,10,10,10,10));
  s->setRedshift(3.);
  s->setHubble(0.7);

  if(!uConv(Unit("kpc a"),Unit("kpc"),0.25,s.get())) return;
  if(!uConv(Unit("kpc a h^-1"),Unit("kpc"),0.25/0.7,s.get())) return;
  if(!uConv(Unit("Msol h^-1")/Unit("kpc^3 a^3 h^-3"),Unit("Msol kpc^-3"),0.7*0.7/(0.25*0.25*0.25),s.get())) return;
  
  std::cerr << "OK" << std::endl;
  
}

void testGasFractionExtraction() {
  std::cerr << "Test MassTransformation...";
  SimSnap *testDat = BaseSimSnap::makeSlab(100,10,10,10,10);
  SimanArray & nHII(testDat->createArray("nHII","HII per HI+HII"));
  SimanArray & nHeII(testDat->createArray("nHeII","HeII per HI+HII"));
  SimanArray & nHeIII(testDat->createArray("nHeIII","HeIII per HI+HII"));
  
  for(unsigned int n=0;n<testDat->getNumParticles();n++) {
    nHII[n]=0.5;
    nHeII[n]=0.1;
    nHeIII[n]=0.01;
  }
  
  SimSnap *t = testDat->copyTransform(MassTransformation(MassTransformation::HI));
  
  if(std::abs(t->getTotalMass()-testDat->getTotalMass()*(1.-testDat->getHeliumMassFrac())/2.)>t->getTotalMass()*1.e-7) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << t->getTotalMass() << " in MassTransformation, should be " << testDat->getTotalMass()*(1.-testDat->getHeliumMassFrac())/2. << std::endl;
    delete t;
    delete testDat;
    return;
  }

  delete t;
  


  delete testDat;

  std::cerr << "OK" << std::endl;
}

void testTransformationMultiplication() {
  std::cerr << "Test transformation multiplication...";
  
  Translation t1(1,2,3);
  Translation t2(3,4,5);
  
  SimSnap * testDat = BaseSimSnap::makeSlab(100,10,10,10,10);

  SimSnap * td = testDat->copyTransform(t1*t2);

  for(unsigned int n=0; n<testDat->getNumParticles(); n++) {
    const Particle *p1 = td->getConstParticle(n);
    Particle p2(*(testDat->getConstParticle(n)));
    
    p2.x+=4;
    p2.y+=6;
    p2.z+=8;

    if(p1->distanceTo(p2)>1.e-6) {
      std::cerr << p1->x << " " << p2.x << " ";
      std::cerr << "FAILED on translation";
      delete td;
      delete testDat;
      
      return;
    }
  }
  
  delete td;

  RotationX rx(PI/2);
  RotationY ry(PI/2);
  RotationZ rz(PI/2);

  td = testDat->copyTransform(rx*rx*rx*rx);



  
  for(unsigned int n=0; n<testDat->getNumParticles(); n++) {
    const Particle *p1 = td->getConstParticle(n);
    const Particle *p2 = testDat->getConstParticle(n);
    
    if(p1->distanceTo(*p2)>1.e-2) {
      std::cerr << "FAILED on rotation" << std::endl;
	
      std::cerr << p1->x << " " << p2->x << std::endl;
      std::cerr << p1->y << " " << p2->y << std::endl;
      std::cerr << p1->z << " " << p2->z << std::endl;
      delete td;
    
      delete testDat;
      
      return;
    }
  }
  
  
  delete td;
  

  delete testDat;
  
  std::cerr << "OK" << std::endl;
}


int main(int argc, char **argv)
{
  std::cerr << "testsiman - - - - - - - - - - - " << std::endl << "Please note this is a very incomplete set of regression tests at the moment" << std::endl << std::endl;
  siman::setVerbose(0); // get Siman to be really quiet so we can see the wood for the trees
  if(argc<2) {
    std::cerr << "**** Provide a filename to perform file-based tests ****" << std::endl;
  } else {
    std::cerr << "**** TESTS USING FILE ****" << std::endl;
    std::cerr << "**** NOTE: failure in this section can reflect file errors OR SimAn errors ****" << std::endl;
    SimSnap *pF = SimSnap::loadFile(argv[1]);
    testDensityEstimation(pF);
    delete pF;
  }
  
  std::cerr << "**** TESTS NOT REQUIRING A FILE ****" << std::endl;

  testUnitConversion();
  testFilteringArrays();
  testTransformationMultiplication();
  testGasFractionExtraction();
  testGridIterators();
  testDirectionalGridIterators();

}








  



