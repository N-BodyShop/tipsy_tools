//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifdef SIMAN_VIS

#include "../base.hpp"
#include "../visual.hpp"

using namespace std;

map <int, pair<CColourMap*,int> > colourmap_menu_map;

void colourMapMenuCallback(int v) {
  pair<CColourMap*, int> obj = colourmap_menu_map[v];
  (obj.first)->menuCallback(obj.second);
}

int allocMenuID(CColourMap *call, int reason) {
  static int mid = 100;
  pair<CColourMap*, int> obj(call,reason);
  colourmap_menu_map[mid]=obj;
  mid++;
  return mid-1;
}

CColourMap::CColourMap(CVisualise *pV) {
  refA=1.;
  refG=1.;
  refB=1.;
  refR=1.;
  pVis = pV;
 
  pParent=NULL;
  display = true;
}

unsigned int CColourMap::supports() {
  return ColourMap;
}

string CColourMap::className() {
  return "CColourMap";
}

CSimanObject * CColourMap::dispatch(string command, std::istream *stream, CScripted *script) {
  if(command=="display") {
    if(!stream->eof()) {
      string s;
      *stream >> s;
      if(s=="on" || s=="1" || s=="yes")
	display = true;
      if(s=="off" || s=="0" || s=="no")
	display=false;
    } else {
      display=!display;
    }
    cout << "Display " << (display?"on":"off") << endl;
    pVis->invalidate();
    return NULL;
  }
  if(command=="hide") {
    display=false;
    pVis->invalidate();
    return NULL;
  }
  if(command=="show") {
    display=true;
    pVis->invalidate();
    return NULL;
  }
  return CSimanObject::dispatch(command,stream,script);
}

int CColourMap::buildMenu() {
  int oldm = glutGetMenu();
  int menu = glutCreateMenu(colourMapMenuCallback);
  
  glutAddMenuEntry("one",1);
  glutAddMenuEntry("two",2);
  glutSetMenu(oldm);
  return menu;


}


void CColourMap::menuCallback(int id) {
  
}


void CColourMap::autoRange() {
  if(pParent==NULL) {
    for(int n=0;n<pVis->pSim->getNumParticles();n++) {
      autoRange_Consider(n);
    }
  } else {
    pParent->autoRange();
  }
}

void CColourMap::autoRange_Consider(int n) {
  // nothing appropriate to do in base class
}


void CColourMap::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1) {
  
  glBegin(GL_QUADS);
  
  glColor4f(refR,refG,refB,1.);
  glVertex2f(x0,y0);
  glVertex2f(x1,y0);
  glVertex2f(x1,y1);
  glVertex2f(x0,y1);

  glEnd();
  
  glColor4f(1,1,1,1);
  CVisualise::writeText(x0,y1,label);
  
  
  glFlush();
    
}



bool CColourMap::operator()(int index, float *gl4f, CParticle *p)
{
  gl4f[0] = refR;
  gl4f[1] = refG;
  gl4f[2] = refB;
  gl4f[3] = refA;
  return display;
}

void CColourMap::setReferenceAlpha(float al) {
  refA=al;
}

float CColourMap::getReferenceAlpha() {
  return refA;
}

void CColourMap::setReferenceColour(float r, float g, float b, float al) {
  refR=r;
  refG=g;
  refB=b;
  refA=al;
}

bool CColourMap::click(int button, int state, float x, float y) {

  if(state == GLUT_UP && button==0) {
    display = !(display);
    return true;
  }

  return false;
}

// Continuuous colour maps


CContinuuousColourMap::CContinuuousColourMap(CVisualise *pV, float mini, float maxi) : CColourMap(pV) {

  if(mini!=maxi) {
    min=mini;
    max=maxi;
  } else {

    // autorange required

    min = 0;
    max = 1;
   
  }

  hideOutOfRange=false;
  colourBy=0;

}


bool CContinuuousColourMap::click(int button, int state, float x, float y) {
  if(state==GLUT_DOWN && button==0) {
    x_button_down =x;
    return false;
  }

  if(state==GLUT_UP && button==0) {
    if(abs(x_button_down-x)>0.01) {
      float max_x = max, min_x = min;
      float range_bot = (x_button_down<x)?x_button_down:x;
      float range_top = (x_button_down>x)?x_button_down:x;

      max = min_x + range_top * (max_x-min_x);
      min = min_x + range_bot * (max_x-min_x);
      return true;
    } else {
      display=!display;
      return true;
    }
  }

  
  if(button==1 && state==GLUT_UP) {
    this->autoRange();
    return true;
  }

  if(button==2 && state==GLUT_UP) {
    colourBy++;
    if(colourBy>use_numOptions+pVis->pSim->getNumArrays()) colourBy=0;
    setColourBy(colourBy);
    this->autoRange();
    return true;
  }

  return false;
}

void CContinuuousColourMap::setColourBy(int by) {
  colourBy = by;
  if(colourBy>use_numOptions) {
    map = pVis->pSim->getArray(colourBy-use_numOptions-1);
  }
  label = getColourLabel(colourBy);
}


void CContinuuousColourMap::autoRange() {
  min = FLT_MAX;
  max = -FLT_MAX;

  CColourMap::autoRange();
}

void CContinuuousColourMap::autoRange_Consider(int n) {
  float val = particleToColourIndex(n,pVis->pSim->getParticle(n));
  if(val>max) max=val;
  if(val<min) min=val;
}


CSimanObject * CContinuuousColourMap::dispatch(string command, std::istream *stream, CScripted *script) {
  if(command=="clip") {
    if(!stream->eof()) {
      string clipmode;
      *stream >> clipmode;
      if(clipmode[0]=='h' && hideOutOfRange!=true) {
	pVis->invalidate();
	hideOutOfRange = true;
      }
      if(clipmode[0]=='s' && hideOutOfRange!=false) {
	pVis->invalidate();
	hideOutOfRange = false;
      }
    }
    cerr << "Clipping mode: " << (hideOutOfRange?"hide":"saturate") << endl;
    return NULL;
  }
  if(command=="range") {
    if(!stream->eof()) {
      *stream >> min;
      *stream >> max;
      pVis->invalidate();
    } else {
      cerr << "[" << min << "," << max << "]" << endl;
    }
    return NULL;
  }
  
  if(command=="colourby") {
    if(stream->eof())
      throw(CSyntaxError("colourby [var]"));
    string t;
    *stream >> t;
    transform(t.begin(),t.end(),t.begin(),(int(*)(int))tolower);
    if(t=="ne")
      setColourBy(useNe);
    if(t=="nhp")
      setColourBy(useNHp);
    if(t=="temp")
      setColourBy(useTemp);
    if(t=="rho" || t=="den")
      setColourBy(useRho);
    if(t=="u" || t=="en")
      setColourBy(useU);
    if(t=="z" || t=="met" || t=="metal") 
      setColourBy(useMetal);
    return NULL;
  }
  return CColourMap::dispatch(command,stream,script);
}

float CContinuuousColourMap::particleToColourIndex(int index, CParticle *p) {
  if(colourBy>use_numOptions)
    return map[index];

  switch(colourBy) {
 
  case useNe:
    return p->ne;
    break;
  case useNHp:
    return p->nHp;
    break;
  case useNH0:
    return 1-p->nHp;
    break;
  case useTemp:
    return log(p->temp)/log(10.);
    break;
  case useRho:
    return log(p->rho)/log(10.);
    break;
  case useU:
    return log(p->u)/log(10.);
    break;
  case useMetal:
    return p->metal;
    break;
  }
}

bool CContinuuousColourMap::operator()(int index, float *gl4f, CParticle *p) {
  return (*this)(particleToColourIndex(index,p),gl4f);
}

bool CContinuuousColourMap::operator()(float n, float *gl4f) {
  return false;
}


void CContinuuousColourMap::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1 ) {


  // default rendering for continuous colour bars

  float col[4]={0,0,0,0};
  
  glBegin(GL_QUADS);
  float deltax = 0.01;
  int ixx=10000 ,ix=-1;
  
  for(float x = x0; x<x1; x+=deltax) {

    (*this)((x-x0)*(this->max-this->min)/(x1-x0)+this->min,col);
    glColor4f(col[0],col[1],col[2],1.);
    
    glVertex2f(x+deltax,y1);
    glVertex2f(x,y1);
    glVertex2f(x,y0);
    glVertex2f(x+deltax,y0);
  
  }
  
  glEnd();
  
  if(max!=min) {
    glColor4f(0,0,0,1);
    
    for(float x = x0; x<x1; x+=deltax) {
      
      ix = (int)((x-x0)*extx);
      
      if((ix%100)<ixx) {
	ostringstream ss;
	ss << (x-x0) * (this->max-this->min)/(x1-x0)+this->min;
	CVisualise::writeText(x,y0+(y1-y0)/2,ss.str());
      }
      
      ixx = ix%100;
      
    }
  }
  
  glColor4f(1,1,1,1);
  CVisualise::writeText(x0,y1,label);

 
  glFlush();
}

string CContinuuousColourMap::getColourLabel(int n) {
   if(n>use_numOptions) {
     return pVis->pSim->getArrayLongName(n-use_numOptions-1);
  } else {
     ostringstream ss;
     switch(n) {
     case useNe:
       return "ne";
       break;
     case useNHp:
       return "nHp";
       break;
     case useNH0:
       return "nH0";
       break;
     case useTemp:
       return "log10(T/K)";
       break;
     case useRho:
       ss << "log10(rho/" << pVis->pSim->getDensityUnits() << ")";
       return ss.str();
       break;
     case useU:
       ss << "log10(u/" << pVis->pSim->getEnergyUnits()  << ")";
       return ss.str();
       break;
     case useMetal:
       ss << "metal";
       return ss.str();
       break;
     }
   }
}




int CContinuuousColourMap::buildMenu() {
  int oldm = glutGetMenu();
  int mainmenu = glutCreateMenu(colourMapMenuCallback);
  int colourbymenu = glutCreateMenu(colourMapMenuCallback);
  
  int maxCol = use_numOptions + pVis->pSim->getNumArrays();
  for(int n=0; n<=maxCol;n++) {
    glutAddMenuEntry(getColourLabel(n).c_str(),allocMenuID(this,n));
  }

  glutSetMenu(mainmenu);
  glutAddSubMenu("Colour by",colourbymenu);
  glutAddMenuEntry("c.one",allocMenuID(this,100));
  glutAddMenuEntry("c.two",allocMenuID(this,101));
  glutSetMenu(oldm);
  return mainmenu;


}

void CContinuuousColourMap::menuCallback(int id) {
  if(id<100) {
    colourBy = id;
    if(colourBy>use_numOptions) {
      map = pVis->pSim->getArray(colourBy-use_numOptions-1);
    }
    label = getColourLabel(colourBy);
    autoRange();
    pVis->invalidate();
  }
}


// "Spectrum" continuuous colour map

CColourMapSpectrum::CColourMapSpectrum(CVisualise *pV,  float mini, float maxi) : CContinuuousColourMap(pV, mini, maxi) {

}

bool CColourMapSpectrum::operator()(float val, float *gl4f)
{
  if(hideOutOfRange) {
    if(val<min || val>max) return false;
  } else {
    if(val<min) val=min;
    if(val>max) val=max;
  }

  float ref = (val-min)/(max-min);
  gl4f[0] = refR * (1-ref);
  
  gl4f[1] = refG * sin(3.1415*ref);
  gl4f[2] = refB * ref;

  gl4f[3] = refA;

  return display;
}


// "Gradient" continuuous colour map

CColourMapGradient::CColourMapGradient(CVisualise *pV, float *ref1, float *ref2, float mini, float maxi) : CContinuuousColourMap(pV,mini,maxi) {

  ref1R = ref1[0];
  ref1G = ref1[1];
  ref1B = ref1[2];
  ref1A = ref1[3];

  ref2R = ref2[0];
  ref2G = ref2[1];
  ref2B = ref2[2];
  ref2A = ref2[3];

}


CSimanObject * CColourMapGradient::dispatch(string command, std::istream *stream, CScripted *script) {

  if(command=="col1") {
    if(!stream->eof()) {
      *stream >> ref1R >> ref1G >> ref1B >> ref1A;
    }
    cout << "col1 (R,G,B,A) = (" << ref1R << "," << ref1G << "," << ref1B << "," << ref1A << ")" << endl;
    pVis->invalidate();
    return NULL;
  }

  if(command=="col2") {
    if(!stream->eof()) {
      *stream >> ref2R >> ref2G >> ref2B >> ref2A;
    }
    cout << "col2 (R,G,B,A) = (" << ref2R << "," << ref2G << "," << ref2B << "," << ref2A << ")" << endl;
    pVis->invalidate();
    return NULL;
  }


  if(command=="alphagrade") {
    if(!stream->eof()) {
      *stream >> ref2R >> ref2G >> ref2B >> ref1A >> ref2A;
      ref1R=ref2R;
      ref1G=ref2G;
      ref1B=ref2B;
    }
    pVis->invalidate();
    return NULL;
  }
  

  return CContinuuousColourMap::dispatch(command,stream,script);
}


bool CColourMapGradient::operator()(float val, float *gl4f)
{
  
  if(hideOutOfRange) {
    if(val<min || val>max) return false;
  } else {
    if(val<min) val=min;
    if(val>max) val=max;
  }

  float ref = (val-min)/(max-min);
  gl4f[0] = ref1R * (1-ref) + ref2R * ref;
  
  gl4f[1] = ref1G * (1-ref) + ref2G * ref;
  gl4f[2] = ref1B * (1-ref) + ref2B * ref;

  gl4f[3] = (ref1A * (1-ref) + ref2A * ref) * refA; // allow reference alpha to change

  return display;
}


// Colour Map by Type (discrete)


CColourMapByType::CColourMapByType(CVisualise *pV, CColourMap *dm, CColourMap *gas, CColourMap *stars) : CColourMap(pV) {
  pDM = dm;
  pGas = gas;
  pStar = stars;
  dmAl = pDM->getReferenceAlpha();
  gasAl = pGas->getReferenceAlpha();
  starAl = pStar->getReferenceAlpha();
  pDM->pParent = this;
  pGas->pParent = this;
  pStar->pParent = this;
  
  
}


int CColourMapByType::buildMenu() {
  
  int oldm = glutGetMenu();
  int menu = glutCreateMenu(colourMapMenuCallback);
  
  glutAddSubMenu("DM",pDM->buildMenu());
  glutAddSubMenu("Gas",pGas->buildMenu());
  glutAddSubMenu("Stars",pStar->buildMenu());
  glutSetMenu(oldm);
  return menu;


}

bool CColourMapByType::operator()(int index, float *gl4f, CParticle *p)
{
  bool freeAtEnd = false;
  if(p==NULL) {
    p = pVis->pSim->getParticle(index);
    freeAtEnd=true;
  }

  bool r=false;

  if(p->type==CParticle::dm)
    r=(*pDM)(index, gl4f, p);
  
  if(p->type==CParticle::gas)
    r=(*pGas)(index, gl4f, p);
  
  if(p->type==CParticle::star)
    r=(*pStar)(index, gl4f, p);
  
  if(freeAtEnd) {
    pVis->pSim->releaseParticle(p);
  }
  return r;
}


CSimanObject * CColourMapByType::getMember(const string & var) {
  if(var=="dm")
    return pDM;
  if(var=="gas")
    return pGas;
  if(var=="stars")
    return pStar;

  return CSimanObject::getMember(var);
}

void CColourMapByType::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1) {
 
  if(pDM->display) pDM->renderColourBar(extx,exty,0,0.33,y0,y1*0.66);
  if(pGas->display) pGas->renderColourBar(extx,exty,0.33,0.66,y0,y1*0.66);
  if(pStar->display) pStar->renderColourBar(extx,exty,0.66,1.0,y0,y1*0.66);

  string title = "DM";
  if(!(pDM->display)) title+=" not displayed";
  CVisualise::writeText(0.01,y1,title);

  title = "Gas";
  if(!(pGas->display)) title+=" not displayed";
  CVisualise::writeText(0.33,y1,title);

  title = "Stars";
  if(!(pStar->display)) title+=" not displayed";
  CVisualise::writeText(0.66,y1,title);
  

}

bool CColourMapByType::click(int button, int state, float x, float y) {
 
  if(button==0 || button==2) {
    if(x<0.33) return pDM->click(button,state,x/0.33,y);
    if(x>0.33 && x<0.66) return pGas->click(button,state,(x-0.33)/0.33,y);
    if(x>0.66) return pStar->click(button,state,(x-0.66)/0.33,y);
  }

  if(button==1)
    this->autoRange();
}
void CColourMapByType::setReferenceAlpha(float a) {
  
  pDM->setReferenceAlpha(a*dmAl);
  pGas->setReferenceAlpha(a*gasAl);
  pStar->setReferenceAlpha(a*starAl);

  
}

void CColourMapByType::autoRange_Consider(int n) {
  CParticle *p = pVis->pSim->getParticle(n);

  if(p->type==CParticle::gas)
    pGas->autoRange_Consider(n);

  if(p->type==CParticle::star)
    pStar->autoRange_Consider(n);

  if(p->type==CParticle::dm) 
    pDM->autoRange_Consider(n);

  pVis->pSim->releaseParticle(p);
}

#endif
