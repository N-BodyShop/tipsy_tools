// Absorption profiles stuff

// simlib includes

#include <siman/base.hpp>
#include <siman/extra.hpp>


ostream & operator <<(ostream &os, list<pair <int,double> > &o) {
  os << "{ ";
  list< pair<int,double> >::iterator i;
 
  for(i=o.begin(); i!=o.end(); i++) {
    if(i!=o.begin()) os << ", ";
    os << "(" << (*i).first << "," << (*i).second << ")";
  }
  os << "} ";
  return os;
}

bool operator ==(list<pair <int,double> > &one, list<pair <int,double> > &two) {
  if(one.size()!=two.size())
    return false;

  list< pair<int,double> >::iterator i1, i2;
 
  for(i1=one.begin(), i2=two.begin(); i1!=one.end(); i1++, i2++) {
    if((*i1).first!=(*i2).first)
      return false;
 
  }

  return true;

}

bool operator !=(list<pair <int,double> > &one, list<pair <int,double> > &two) {
  return !(one==two);
}

void densityCheck(CSimSnap *pGas) {
  
  clock_t start;
  float delta1,delta2;
  
  float x1,x2,y1,y2,z1,z2;
 
  cerr << pGas->getTotalMass() << " " << pGas->getMassUnits() << " in gas" << endl;
  pGas->getExactBoundaries(x1,x2,y1,y2,z1,z2);
  

  int numPart = pGas->getNumParticles();
  
  
  cout << "smooth\trho_p\trho_measured\tratio\tdelta_time" << endl;

  CParticle *p;
  float rho_estimate_g;

  for(int n=numPart-1;n!=0;--n) {
    
    
    p = pGas->getParticle(n);
    if(p->x>x1+1 && p->x<x2-1 && p->y>y1+1 && p->y<y2-1 && p->z>z1+1 && p->z<z2-1) {
    
      start = clock();

      //griddedGas.consistentSmooth(n,1.e-6);
      rho_estimate_g = pGas->estimateDensity(n);
      delta2=(float)(clock()-start)/(float)CLOCKS_PER_SEC;
    
      cout << pGas->getSmooth(n) << "\t" << p->rho << "  \t" << rho_estimate_g << "\t" << p->rho/rho_estimate_g << "\t"  << delta2 <<  endl;
    
    } 
  }

}

void thinIoniseCheck(CSimSnap *pGas) {

    
  // CIonise ionise(8.55e-15,7.24e-15,1.71e-16,pGas->getDensityUnits(),(float)0.);
  CIonise ionise(0.1,1,pGas);
  
  ionise.setFlag(CIonise::periodicBoundsX);
  ionise.setFlag(CIonise::periodicBoundsY);

  int i=0;
  CParticle p;
  p.rho = 1 * units::CUnit(units::den_protonsPerCm3).convertTo(pGas->getDensityUnits());
  cout << "temp dependence test using rho = " << p.rho << pGas->getDensityUnits() << endl;

  for(float t=0;t<30000;t+=1000) {
    p.temp = t;

    double nH0, nHe0, nHp, nHep, nHepp, ne;
    ionise.calculateFractions(&p, nH0, nHe0, nHp, nHep, nHepp, ne);
    
    
    double my_ne = ne/(nH0+nHp);
    
    
    cout <<  p.temp << "\t" << my_ne <<endl;
    
    //cout << "(" << nH0 << "\t" << nHe0 << "\t" << nHp << "\t" << nHep << "\t" << nHepp << ")" << endl;
  
  }

  /*
  cout << endl << "with sim data:" << endl;

  for(int n=0;n<pGas->getNumParticles() && i<10;n+=100) {
    CParticle *p = pGas->getParticle(n);
    
    
    
    cout << p->ne << "\t" << p->temp;
    
    // p->temp*=1.4993 ;
    double nH0, nHe0, nHp, nHep, nHepp, ne;
    ionise.calculateFractions(p, nH0, nHe0, nHp, nHep, nHepp, ne);
    
    
    double my_ne = ne/(nH0+nHp);
    
    p->ne = my_ne;
    
    cout << "\t" << my_ne << "\t" << p->temp << endl;
    
    cout << "(" << nH0 << "\t" << nHe0 << "\t" << nHp << "\t" << nHep << "\t" << nHepp << ")" << endl;
    i++;
    
    pGas->releaseParticle(p);
    
  }
  */

}


void thickIoniseCheck(CGrid *pGas) {

  float j_min=13;
  float j_max=14;
  int n_div=100;
  
  float j[n_div];
  
  
  for(int n=0; n<n_div; n++) {
    j[n]=0.1;
  }
  

 
  CIonise ionise(0.01,-1.8,pGas);

  ionise.setFlag(CIonise::writeGam);
  ionise.setFlag(CIonise::verbose);
  ionise.thickRadiative(pGas);
  


}

void reemissionCheck() {

  
  for(float temp=10000; temp<30000; temp+=5000) {
    for(float rho=0.000001; rho<10; rho*=10) {
      CSimSnap *pGas = CBaseSimSnap::makeSlab(5000,1.,1.,rho,temp,CBaseSimSnap::genRandom,units::CUnit(units::len_cm),units::CUnit(units::mass_g));
      
      CGrid griddedGas(pGas,3,3,3,0,0,false);
      
      CIonise ionise(0,-1,pGas);
      ionise.setFlag(CIonise::periodicBoundsAll);
      
      ionise.thickRadiative(&griddedGas);
      
      cout << "\t" << pGas->getParticle(0)->ne;
      
      delete pGas;
    }
    cout << endl;
  }

 
}

void ratesCheck() {
  
  CSimSnap *pF = CBaseSimSnap::makeSlab(20,0.02,0.02,0.1,10); // of little importance, just for instantiation of CIonise

  string pathname = (string) getenv("SIMAN_DATA");
  pathname+="/UVNORM";
  
  ifstream file_uv_norm(pathname.c_str());


  if(file_uv_norm.good()) {

    float z_f=0, j21_f=0, z_xf=100, j21_xf=0;
    while(!file_uv_norm.eof()) {
      
      file_uv_norm >> z_f >> j21_f;

      CIonise i(j21_f,-1.6,pF);

      cout << z_f << "\t" <<  i.gammaHg << "\t" <<  i.gammaHeg << "\t" <<  i.gammaHePg << endl;
      
    }
  } else cerr << "file error" << endl;

  delete pF;
}

void workBack() {
  CSimSnap *pF = CBaseSimSnap::makeSlab(20,0.02,0.02,0.1,10); // of little importance, just for instantiation of CIonise

  double z,skip;
  float j=0.1,alpha=-0.1;
  CIonise first_guess(j,alpha,pF);
 
  
  while(!cin.eof()) {

    double gammaHg, gammaHeg, gammaHePg;
    double out_gammaHg=1, out_gammaHeg=1, out_gammaHePg=1;

    cin >> z;
    cin >> skip;

    //cout << "target gammaHg="; 
    cin >> gammaHg;
    cin >> skip >> gammaHeg >> skip;
    //cout << "target gammaHePg="; 
    cin >> gammaHePg;
    cin >> skip;

  

    float a_jump=0.01;
    j=first_guess.j21=0.1;
    alpha=first_guess.alpha=-0.1;


    while(out_gammaHePg>gammaHePg) {
      alpha-=a_jump;
      first_guess.alpha=alpha;
      first_guess.calculateGammas();

      j *= gammaHg/first_guess.gammaHg;
      first_guess.j21=j;
      first_guess.calculateGammas();


      out_gammaHg = first_guess.gammaHg;
      out_gammaHeg = first_guess.gammaHeg;
      out_gammaHePg = first_guess.gammaHePg;
    }

    cout << j << "\t" << alpha << "\t" << out_gammaHg << "\t" << out_gammaHeg << "\t" << out_gammaHePg << "\t" << out_gammaHeg/gammaHeg << endl;

  }
  delete pF;
}


int main(int argc, char **argv)
{

  char path[1024];

  CSimSnap *pGas=NULL;
  CSimSnap *pF = NULL;
  int ngrid = 10;

  if(argc < 2) {
    float den = 1 * units::CUnit(units::den_protonsPerCm3).convertTo(units::CUnit(units::den_MsolPerKpc3));

    //cout << "den = " << den << endl;

    pGas = CBaseSimSnap::makeSlab(50000,0.0002,0.0002,den,10000,CBaseSimSnap::genRandom);
    cout << pGas->getDensityUnits() << endl;

  } else {
    
    strcpy(path,argv[1]);
    
    pF = CSimSnap::loadFile(path);
    
    pF->convertUnits();
    
    
    CSphere sp_desc(0,0,0,200);
    CSubset *sphere = new CSubset(pF,sp_desc);

    cerr << sphere->getNumParticles() << endl;

    CParticleTypeFilter gasfilter(CParticle::gas);

    pGas = new CSubset(sphere,gasfilter);

    cerr << pGas->getNumParticles() << "/" << pGas->getTotalMass() << " " << pGas->getMassUnits() << " tot " << endl;

    if(argc>2)
      ngrid=atoi(argv[2]);

  }


  CGrid griddedGas(pGas,ngrid,ngrid,ngrid);

  griddedGas.initialiseSPH(40);

  //densityCheck(&griddedGas);
  // thinIoniseCheck(&griddedGas);

  thickIoniseCheck(&griddedGas);

  // reemissionCheck();


  ostringstream oss;
  oss << argv[1] << ".ionise-" << ngrid;

 
  pGas->write(oss.str(),CSimSnap::native);
  
  
  delete pGas;
}


