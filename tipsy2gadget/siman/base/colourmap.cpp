// colourmap.cpp - part of SimAn Simulation Analysis Library
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









#ifdef SIMAN_VIS

#include "../base.hpp"
#include <boost/lexical_cast.hpp>

namespace siman {

  
  extern void colourMapMenuCallback(int v);

  extern void resetMenuIDs();

  extern int allocMenuID(ColourMap *call, int reason);

  ColourMap::ColourMap(Visualiser *pV) {
    refA=1.;
    refG=1.;
    refB=1.;
    refR=1.;
    pVis = pV;
 
    menuID=-1;
    pParent=NULL;
    display = true;
  }

  
  void ColourMap::setCol(const vector<double> &arr) {
    if(arr.size()!=4 && arr.size()!=3)
      throw(SyntaxError("ColourMap::setCol: requires 3 or 4-element array"));
    vector<double>::const_iterator i = arr.begin();
    refR=*i;
    refG=*(++i);
    refB=*(++i);
    i++;
    if(i!=arr.end())
      refA=*i;
    pVis->invalidate();
  }

  vector<double> ColourMap::getCol() {
    vector<double> r;
    r.push_back(refR);
    r.push_back(refG);
    r.push_back(refB);
    r.push_back(refA);
    return r;
  }

  int ColourMap::buildMenu() {
    return 0;
  }


  void ColourMap::menuCallback(int id) {
  
  }


  void ColourMap::autoRange() {
    if(pParent==NULL) {
      for(unsigned int n=0;n<pVis->pSim->getNumParticles();n++) {
	autoRange_Consider(n);
      }
    } else {
      pParent->autoRange();
    }
  }

  string ColourMap::getCurrentLabel() {
    return label;
  }
  void ColourMap::autoRange_Consider(unsigned int n) {
    // nothing appropriate to do in base class
  }


  void ColourMap::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1, bool renderLab) {
  
    glBegin(GL_QUADS);
  
    glColor4f(refR,refG,refB,1.);
    glVertex2f(x0,y0);
    glVertex2f(x1,y0);
    glVertex2f(x1,y1);
    glVertex2f(x0,y1);

    glEnd();
  
    if(renderLab) {
      glColor4f(1,1,1,1);
      Visualiser::rasterText(x0,y1,label);
  
    }
  
    glFlush();
    
  }



  bool ColourMap::operator()(int index, float *gl4f, const Particle *p)
  {
    gl4f[0] = refR;
    gl4f[1] = refG;
    gl4f[2] = refB;
    gl4f[3] = refA;
    return display;
  }

  void ColourMap::setReferenceAlpha(float al) {
    refA=al;
  }

  float ColourMap::getReferenceAlpha() {
    return refA;
  }

  void ColourMap::setDisplay(bool visible) {
    display = visible;
    pVis->invalidate();
  }
  
  bool ColourMap::getDisplay() {
    return display;
  }
  

  void ColourMap::setReferenceColour(float r, float g, float b, float al) {
    refR=r;
    refG=g;
    refB=b;
    refA=al;
  }

  bool ColourMap::click(int button, int state, float x, float y) {

    if(state == GLUT_UP && button==0) {
      display = !(display);
      return true;
    }

    return false;
  }

  // Continuuous colour maps


  ContinuuousColourMap::ContinuuousColourMap(Visualiser *pV, float mini, float maxi) : ColourMap(pV) {

    if(mini!=maxi) {
      min=mini;
      max=maxi;
    } else {

      // autorange required

      min = 0;
      max = 1;
   
    }

    hideOutOfRange=true;
    colourBy=0;
    logScale=false;

  }


  bool ContinuuousColourMap::click(int button, int state, float x, float y) {
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
      
      return true;
    }

    return false;
  }

  void ContinuuousColourMap::setColourBy(int by) {
    colourBy = by;
    if(colourBy>use_numOptions) {
      map = &(pVis->pSim->getConstArray(colourBy-use_numOptions-1));
    }
    label = getColourLabel(colourBy);
    this->autoRange();
    pVis->invalidate();
  }

  void ContinuuousColourMap::setColourBy(string col_b) {
    if(col_b=="ne") colourBy = useNe;
    else if(col_b=="rho") colourBy = useRho;
    else if(col_b=="temp") colourBy = useTemp;
    else if(col_b=="u") colourBy = useU;
    else if(col_b=="metal") colourBy = useMetal;
    else {
      try {
	colourBy = use_numOptions+1+pVis->pSim->getArrayIndex(col_b);
      } catch (UnknownArray &e) {
	cerr << "ContinuuousColourMap: Can't find a colour index " << col_b << endl;
	return;
      }
      
    }
    setColourBy(colourBy); 
  
  }

  void ContinuuousColourMap::setRange(const vector<double> &arr) {
    if(arr.size()!=2)
      throw(SyntaxError("ContinuuousColourMap::setRange: requires 2-element array"));
    min=*(arr.begin());
    max=*(++arr.begin());
    pVis->invalidate();
  }

  vector<double> ContinuuousColourMap::getRange() {
    vector<double> r;
    r.push_back(min);
    r.push_back(max);
    return r;
  }

  void ContinuuousColourMap::autoRange() {
    min = FLT_MAX;
    max = -FLT_MAX;

    ColourMap::autoRange();
  }

  void ContinuuousColourMap::autoRange_Consider(unsigned int n) {
    float val = particleToColourIndex(n,pVis->pSim->getConstParticle(n));
    if(val>max) max=val;
    if(val<min) min=val;
    if(colourBy==useZero) {
      max=1.;
      min=-1.;
    }
  }


  float ContinuuousColourMap::particleToColourIndex(int index, const Particle *p) {
    float r=0;
    
    if(colourBy>use_numOptions)
      r=((const SimanArray &)(*map))[index];

    switch(colourBy) {
 
    case useNe:
      r=p->ne;
      break;
  
    case useTemp:
      r=p->temp;
      break;
    case useRho:
      r=p->rho;
      break;
    case useU:
      r=p->u;
      break;
    case useZero:
      return 0.;
      break;
    case useMetal:
      r=p->metal;
      break;
    }

    if(logScale)
      return log(r)/2.302585;
    else
      return r;
  }

  bool ContinuuousColourMap::operator()(int index, float *gl4f, const Particle *p) {
    return (*this)(particleToColourIndex(index,p),gl4f);
  }

  bool ContinuuousColourMap::operator()(float n, float *gl4f) {
    return false;
  }


  void ContinuuousColourMap::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1, bool renderLab ) {

    float text_height = 12./(float)exty;

    // default rendering for continuous colour bars

    glColor4f(0,0,0,1);
    glBegin(GL_QUADS);
    glVertex2f(x0,y0);
    glVertex2f(x1,y0);
    glVertex2f(x1,y1);
    glVertex2f(x0,y1);
    glEnd();

    float col[4]={0,0,0,0};
  
    glBegin(GL_QUADS);
    float deltax = 0.01;
  
    for(float x = x0; x<x1; x+=deltax) {

      (*this)((x-x0)*(this->max-this->min)/(x1-x0)+this->min,col);
      glColor4f(col[0],col[1],col[2],1.);
    
      glVertex2f(x+deltax,y0+text_height);
      glVertex2f(x,y0+text_height);
      glVertex2f(x,y1);
      glVertex2f(x+deltax,y1);
  
    }
  
    glEnd();
  
    if(max!=min) {
    

      int nticks_targ = (int)((x1-x0)*(float)extx/30.);
      float pertick_targ= abs(max-min)/(float)nticks_targ;
      float pertick_use = std::pow(10.,(int)(log(pertick_targ)/log(10.)+0.5));
    
      if(pertick_use/5>pertick_targ)
	pertick_use/=5;
      else if(pertick_use/2>pertick_targ)
	pertick_use/=2;

      

      float rpos_get[4];
      glColor4f(1,1,1,1);
      float begin = pertick_use*((int)(min/pertick_use)+1);
      if(pertick_use>0) {
	for(float v = begin; v<max; v+=pertick_use) {
	  glBegin(GL_LINES);
	  float x = x0+((v-min)/(max-min))*(x1-x0);
	  glVertex2f(x,y0);
	  glVertex2f(x,y1);

	  /*
	    // uncomment for minor ticks

	  for(float sub=v+pertick_use/5; sub<v+pertick_use; sub+=pertick_use/5) {
	    float x = x0+((sub-min)/(max-min))*(x1-x0);    	
	    glVertex2f(x,y0);
	    glVertex2f(x,y0+(y1-y0)/4);
	
	  }
	  */
	  glEnd();
	
	  ostringstream ss;
	  ss << v;
	  Visualiser::rasterText(x+3./(float)extx,y0+1./exty,ss.str());
	  glGetFloatv(GL_CURRENT_RASTER_POSITION,rpos_get);
	}
      }

    }
  
    if(renderLab) {
      glColor4f(1,1,1,1);
      Visualiser::rasterText(x0+(x1-x0)/2.,y1,label);
    }
 
    glFlush();
  }

  string ContinuuousColourMap::getColourLabel(int n) {
    
    ostringstream ss;
    
    if(logScale)
      ss << "log10(";

    if(n>use_numOptions) {
      ss << pVis->pSim->getArrayLongName(n-use_numOptions-1) << "/" << pVis->pSim->getArrayUnits(n-use_numOptions-1);
    } else {
      switch(n) {
      case useNe:
	ss << "ne";
	break;
      case useTemp:
	ss << "T/K";
	break;
      case useRho:
	ss << "rho/" << pVis->pSim->getDensityUnits();
	return ss.str();

      case useU:
	ss << "u/" << pVis->pSim->getEnergyUnits();

	break;
      case useMetal:
	ss << "metal/solar";

	break;
      case useZero:
	return "zero";
	break;
      }
      
      
    }
    
    if(logScale)
      ss << ")";

    return ss.str();

  }




  int ContinuuousColourMap::buildMenu() {
    
    int oldm = glutGetMenu();

    // destroy previous version of menu, if it exists
    if(menuID!=-1)
      glutDestroyMenu(menuID);

    menuID = glutCreateMenu(colourMapMenuCallback);
    int colourbymenu = glutCreateMenu(colourMapMenuCallback);
  
    int maxCol = use_numOptions + pVis->pSim->getNumArrays();
    for(int n=0; n<=maxCol;n++) {
      glutAddMenuEntry(string(((n==colourBy)?"(*)":"( )") + getColourLabel(n)).c_str(),allocMenuID(this,n));
    }

    glutSetMenu(menuID);
    glutAddSubMenu("Colour by",colourbymenu);
    glutAddMenuEntry((string(logScale?"(*)":"( )") + "Log scale").c_str(),allocMenuID(this,100));
    glutAddMenuEntry((string(hideOutOfRange?"(*)":"( )") + "Hide particles when out of range").c_str(),allocMenuID(this,101));
    glutSetMenu(oldm);
    return menuID;


  }

  
  bool ContinuuousColourMap::getHideOutOfRange() {
    return hideOutOfRange;
  }

  void ContinuuousColourMap::setHideOutOfRange(bool to) {
    hideOutOfRange=to;
    pVis->invalidate();
  }
  
  bool ContinuuousColourMap::getLogScale() {
    return logScale;
  }

  void ContinuuousColourMap::setLogScale(bool to) {
    logScale=to;
    label=getColourLabel(colourBy);
    autoRange();
    pVis->invalidate();
  }

  void ContinuuousColourMap::menuCallback(int id) {
    if(id<100) {
      colourBy = id;
      if(colourBy>use_numOptions) {
	map = &(pVis->pSim->getConstArray(colourBy-use_numOptions-1));
      }
      label = getColourLabel(colourBy);
      autoRange();
      pVis->invalidate();
    }
    if(id==100) {
      logScale=!logScale;
      label=getColourLabel(colourBy);
      autoRange();
      pVis->invalidate();
      
    }
    if(id==101) {
      hideOutOfRange=!hideOutOfRange;
      pVis->invalidate();
    }
  }


  // "Spectrum" continuuous colour map

  ColourMapSpectrum::ColourMapSpectrum(Visualiser *pV,  float mini, float maxi) : ContinuuousColourMap(pV, mini, maxi) {

  }


  ColourMapSpectrum::ColourMapSpectrum(Visualiser *pV) : ContinuuousColourMap(pV, 0, 1) {

  }

  bool ColourMapSpectrum::operator()(float val, float *gl4f)
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

  ColourMapGradient::ColourMapGradient(Visualiser *pV, float *ref1, float *ref2, float mini, float maxi) : ContinuuousColourMap(pV,mini,maxi) {

    ref1R = ref1[0];
    ref1G = ref1[1];
    ref1B = ref1[2];
    ref1A = ref1[3];

    ref2R = ref2[0];
    ref2G = ref2[1];
    ref2B = ref2[2];
    ref2A = ref2[3];

  }

  ColourMapGradient::ColourMapGradient(Visualiser *pV, const vector<double> &one, const vector<double> &two) : ContinuuousColourMap(pV,0,1) {
    setCol1(one);
    setCol2(two);
  }



  void ColourMapGradient::setCol1(const vector<double> &arr) {
    if(arr.size()!=4 && arr.size()!=3)
      throw(SyntaxError("ColourMapGradient::setCol1: requires 3 or 4-element array"));
    vector<double>::const_iterator i = arr.begin();
    ref1R=*i;
    ref1G=*(++i);
    ref1B=*(++i);
    i++;
    if(i!=arr.end())
      ref1A=*i;
    pVis->invalidate();
  }

  vector<double> ColourMapGradient::getCol1() {
    vector<double> r;
    r.push_back(ref1R);
    r.push_back(ref1G);
    r.push_back(ref1B);
    r.push_back(ref1A);
    return r;
  }



  void ColourMapGradient::setCol2(const vector<double> &arr) {
    if(arr.size()!=4 && arr.size()!=3)
      throw(SyntaxError("ColourMapGradient::setCol2: requires 3 or 4-element array"));
    vector<double>::const_iterator i = arr.begin();
    ref2R=*i;
    ref2G=*(++i);
    ref2B=*(++i);
    i++;
    if(i!=arr.end())
      ref2A=*i;
    pVis->invalidate();
  }

  vector<double> ColourMapGradient::getCol2() {
    vector<double> r;
    r.push_back(ref2R);
    r.push_back(ref2G);
    r.push_back(ref2B);
    r.push_back(ref2A);
    return r;
  }


  bool ColourMapGradient::operator()(float val, float *gl4f)
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


  ColourMapByType::ColourMapByType(Visualiser *pV, ColourMap *dm, ColourMap *gas, ColourMap *stars) : ColourMap(pV) {
    pDM = dm;
    pGas = gas;
    pStar = stars;
    refA=0.2;
    dmAl = pDM->getReferenceAlpha();
    gasAl = pGas->getReferenceAlpha();
    starAl = pStar->getReferenceAlpha();
    pDM->setReferenceAlpha(getReferenceAlpha()*dmAl);
    pGas->setReferenceAlpha(getReferenceAlpha()*gasAl);
    pStar->setReferenceAlpha(getReferenceAlpha()*starAl);
    pDM->pParent = this;
    pGas->pParent = this;
    pStar->pParent = this;
  
  
  }

  ColourMap & ColourMapByType::getGasMap() {
    return *pGas;
  }

  ColourMap & ColourMapByType::getDMMap() {
    return *pDM;
  }

  ColourMap & ColourMapByType::getStarMap() {
    return *pStar;
  }


  void ColourMapByType::setStarMap(ColourMap &cm) {
    pStar=&cm;
    pVis->invalidate();
  }


  int ColourMapByType::buildMenu() {
  
    
    // destroy previous version of menu, if it exists
    if(menuID!=-1)
      glutDestroyMenu(menuID);

    int oldm = glutGetMenu();
    menuID = glutCreateMenu(colourMapMenuCallback);
  
    int men_temp = pDM->buildMenu();
    if(men_temp!=0) glutAddSubMenu("DM",men_temp);

    men_temp = pGas->buildMenu();
    if(men_temp!=0) glutAddSubMenu("Gas",men_temp);

    men_temp = pStar->buildMenu();
    if(men_temp!=0) glutAddSubMenu("Star",men_temp);

    glutSetMenu(oldm);
    return menuID;


  }

  bool ColourMapByType::operator()(int index, float *gl4f, const Particle *p)
  {

    if(p==NULL) {
      p = pVis->pSim->getConstParticle(index);
    }

    bool r=false;

    if(p->type==Particle::dm)
      r=(*pDM)(index, gl4f, p);
  
    if(p->type==Particle::gas)
      r=(*pGas)(index, gl4f, p);
  
    if(p->type==Particle::star)
      r=(*pStar)(index, gl4f, p);

    return r;
  }


  void ColourMapByType::renderColourBar(int extx, int exty, float x0, float x1, float y0,float y1, bool renderLab) {
    float text_height = 16./(float)exty;

    glColor4f(0,0,0,0.5);
    glBegin(GL_QUADS);
    glVertex2f(x0,y0);
    glVertex2f(x1,y0);
    glVertex2f(x1,y1);
    glVertex2f(x0,y1);
    glEnd();


    if(pDM->display) pDM->renderColourBar(extx,exty,0,0.33,y0,y1-text_height,false);
    if(pGas->display) pGas->renderColourBar(extx,exty,0.33,0.66,y0,y1-text_height,false);
    if(pStar->display) pStar->renderColourBar(extx,exty,0.66,1.0,y0,y1-text_height,false);

    if(renderLab) {
       float text_height = 12./(float)exty;
       
      string title = "DM: ";
      if(!(pDM->display)) title+=" not displayed";
      else title+=pDM->getCurrentLabel();
      title+=" alpha="+boost::lexical_cast<string>(refA*dmAl);
      
      Visualiser::rasterText(0.01,y1-text_height,title);
      
      title = "Gas: ";
      if(!(pGas->display)) title+=" not displayed";
      else title+=pGas->getCurrentLabel();
      title+=" alpha="+boost::lexical_cast<string>(refA*gasAl);

      Visualiser::rasterText(0.33,y1-text_height,title);
   
      
      title = "Stars: ";
      if(!(pStar->display)) title+=" not displayed";
      else title+=pStar->getCurrentLabel();
      title+=" alpha="+boost::lexical_cast<string>(refA*starAl);
      Visualiser::rasterText(0.66,y1-text_height,title);
    }

    glColor4f(1,1,1,1);
    glBegin(GL_LINES);
    glVertex2f(x0,y1);
    glVertex2f(x1,y1);
    glEnd();

  }

  bool ColourMapByType::click(int button, int state, float x, float y) {
 
    if(button==0 || button==2) {
    
      if(x<0.33) return pDM->click(button,state,x/0.33,y);
      if(x>0.33 && x<0.66) return pGas->click(button,state,(x-0.33)/0.33,y);
      if(x>0.66) return pStar->click(button,state,(x-0.66)/0.33,y);
      
    }

    if(button==1) {
      this->autoRange();
      pVis->invalidate();
    }

    if(button==3 || button==5) {
      // pVis->setTrans(pVis->getTrans()*1.1);
     
      if(x<0.33) dmAl*=1.1;
      if(x>0.33 && x<0.66) gasAl*=1.1;
      if(x>0.66) starAl*=1.1;
      setReferenceAlpha(refA);
      pVis->setTrans(pVis->getTrans()); // force update - either invalidate or update vertex shader
      
      
    }

    if(button==4) {
      // pVis->setTrans(pVis->getTrans()/1.1);

      if(x<0.33) dmAl/=1.1;
      if(x>0.33 && x<0.66) gasAl/=1.1;
      if(x>0.66) starAl/=1.1;
      setReferenceAlpha(refA);
      pVis->setTrans(pVis->getTrans()); // force update - either invalidate or update vertex shader
      
     }

    return false;
  }

  void ColourMapByType::setReferenceAlpha(float a) {
    refA=a;
    pDM->setReferenceAlpha(a*dmAl);
    pGas->setReferenceAlpha(a*gasAl);
    pStar->setReferenceAlpha(a*starAl);

  
  }

  void ColourMapByType::autoRange_Consider(unsigned int n) {
    const Particle *p = pVis->pSim->getConstParticle(n);

    if(p->type==Particle::gas)
      pGas->autoRange_Consider(n);

    if(p->type==Particle::star)
      pStar->autoRange_Consider(n);

    if(p->type==Particle::dm) 
      pDM->autoRange_Consider(n);

  }



}

#endif
