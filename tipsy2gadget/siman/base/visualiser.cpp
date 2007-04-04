// visualiser.cpp - part of SimAn Simulation Analysis Library
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
#include <typeinfo>
#include <boost/lexical_cast.hpp>

#define PI 3.141592


namespace siman {

  Visualiser::Visualiser(SimSnap *pSimi, unsigned int flagsi) {

    pSim=pSimi;
    flags = flagsi;
    
    pCM = NULL;
    
    noUpdatesLastDisplay=false;
    latestUpdate=-1;
    trans=0.4;

    numberLists = 10;

    shader_prog = 0;

    use_fragshader = false;

    pGrid=NULL;


  }

  void Visualiser::redisplay() {
    cerr << "Visualiser: don't know how to redisplay from base class" << endl;
  }

  void Visualiser::printInfoLog(GLuint b)
  {
    int infologLength = 0;
    int charsWritten  = 0;
    char *infoLog;
    
    glGetObjectParameterivARB(b, GL_OBJECT_INFO_LOG_LENGTH_ARB,
			      &infologLength);
    
    if (infologLength > 0)
      {
	infoLog = (char *)malloc(infologLength);
	glGetInfoLogARB(b, infologLength, &charsWritten, infoLog);
	fprintf(stderr,"%s\n",infoLog);
	free(infoLog);
      }
  }

  void Visualiser::setUpVertexShade() {

    if(!use_fragshader)
      return;

    if(shader_prog!=0) {
      disableVertexShader();
      glDeleteObjectARB(shader_prog);
      glDeleteObjectARB(shader_obj_vert);
      glDeleteObjectARB(shader_obj_frag);
    }
   
    string shader_name = "pixels";
    
    if((flags&splatting)!=0)
      shader_name = "splatting";

    // vectors override splatting

    if((flags&vectors)!=0)
      shader_name = "vectors";
        
    char *vs = (char*)malloc(20480);
    char *fs = (char*)malloc(20480);

    memset(vs,0,20480);
    memset(fs,0,20480);


    string filename = std::string(getenv("SIMAN_DATA"))+std::string("/shaders/vshade_")+shader_name;
    FILE* fh = fopen(filename.c_str(),"r");
    if(fh==NULL) {
      cerr << "Can't read shader " << filename;
      free(vs);
      free(fs);
      return;
    }
    fread(vs,1,20480,fh);
    fclose(fh);
   
    filename = std::string(getenv("SIMAN_DATA"))+std::string("/shaders/fshade_")+shader_name;
    fh = fopen(filename.c_str(),"r");
    if(fh==NULL) {
      cerr << "Can't read shader " << filename;
      free(vs);
      free(fs);
      return;
    }
    fread(fs,1,20480,fh);
    fclose(fh);
        
    // create shaders and compile source
    shader_obj_frag = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);	
    shader_obj_vert = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);	    
    glShaderSourceARB(shader_obj_frag, 1, (const GLcharARB**) &fs,NULL);
    glShaderSourceARB(shader_obj_vert, 1, (const GLcharARB**) &vs,NULL);   
    glCompileShaderARB(shader_obj_frag);
    glCompileShaderARB(shader_obj_vert);
    
    
    // create program and link shaders
    shader_prog = glCreateProgramObjectARB();
    glAttachObjectARB(shader_prog,shader_obj_frag);
    glAttachObjectARB(shader_prog,shader_obj_vert);
    glLinkProgramARB(shader_prog);

    free(vs);
    free(fs);

    // setup prorgam for use
    


    GLint result;
    
    glGetObjectParameterivARB(shader_prog,GL_OBJECT_LINK_STATUS_ARB,&result);

    if(result!=1) {
      cerr << "Visualiser: ARB vertex shader failed to compile/link." << endl;
      printInfoLog(shader_obj_frag);
      printInfoLog(shader_obj_vert);
      printInfoLog(shader_prog);
    } else {
      cerr << "Visualiser: using vertex shader '" << shader_name << "'" <<  endl;
      glUseProgramObjectARB(shader_prog);
    }

  }

  void Visualiser::disableVertexShader() {
    
    glUseProgramObjectARB(0);

    
  }

  void Visualiser::enableVertexShader() {
    glUseProgramObjectARB(shader_prog);
  }

  void Visualiser::deleteShaders() {

    if(use_fragshader) {
      glDeleteObjectARB(shader_prog);
    }
  }

  Visualiser::~Visualiser() {
    
    if(pGrid!=NULL && static_cast<SimSnap*>(pGrid)!=pSim)
      delete pGrid;

    list<Annotate*>::iterator i;
    for(i=annotationList.begin();i!=annotationList.end();i++) {
      delete (*i);
    }

 }

  bool Visualiser::hasFlag(unsigned int flag) {
    return((flags&flag)==flag);
  }

  void Visualiser::setFlag(unsigned int flag, bool val) {
    if(val) {
      // set
      flags=flags|flag;
    } else {
      flags=flags & ~flag;
    }
    if((flag&noAutoSimp)!=0)
      fullReset();
    if((flag&vectors)!=0 || (flag&splatting)!=0) {
      invalidate();
      setUpVertexShade();
    }
    if((flag&ortho)!=0)
      initProjection();
    updateWidgets();
    
  }

  void Visualiser::toggleFlag(unsigned int flag) {
    setFlag(flag,!hasFlag(flag));
  }

  void Visualiser::fullReset() {
  
    static bool firstInit = true;

    maxmass = 0.;
    
    for(unsigned int n=0; n<pSim->getNumParticles(); n++) {
      const Particle *p = pSim->getConstParticle(n);
      if(p->mass>maxmass) maxmass=p->mass;
    }

    if(firstInit) {
      // Required initialisation once GL context is available, but only wants to be done
      // once. Slightly messy to have it here, might relocate at some point.
      GLenum res = glewInit();
      if(res!=GLEW_OK) {
	cerr << glewGetErrorString(res);
	exit(0);
      }
      
      
      if(GL_TRUE!=glewGetExtension("GL_ARB_fragment_shader")) {
	use_fragshader=false;
	cerr << "Visualiser: Vertex shader unavailable; some visualisation functions may be slow or not fully functional" << endl;
      } else {
	use_fragshader=true;
	setUpVertexShade();
      }
      firstInit=false;
    }

    try {
      numPartTarget = boost::lexical_cast<int>(config["vis_numPartTarget"]);
      numPartUpdatedTarget = boost::lexical_cast<int>(config["vis_numPartUpdatedTarget"]);
      numberSplitsPerCell = boost::lexical_cast<int>(config["vis_numberSplitsPerCell"]);
      displayWhileUpdating = boost::lexical_cast<bool>(config["vis_displayWhileUpdating"]);
      numberLists = 10;

      if(numberSplitsPerCell>59) numberSplitsPerCell=59;
      if(numberSplitsPerCell<2) numberSplitsPerCell=2;
    } catch(boost::bad_lexical_cast &) {
      cerr << "Visualiser: Configuration file bad, using defaults" << endl;
      numPartTarget = 400000;
      numPartUpdatedTarget = 50000;
      numberSplitsPerCell = 10;
      numberLists = 10;
      displayWhileUpdating = false;
    }

    
    trans  = 0.5;
    pointScale = 3;
    vectorScale = 0.001;
    units = "pc";
    unit_scaler =Unit(units).convertTo(pSim->getDistanceUnits(),pSim);
    // initModelView();

    if(pGrid!=NULL) {
      if(static_cast<SimSnap*>(pGrid)!=pSim) delete pGrid;
      pGrid=NULL;
    }

    if((flags&noAutoSimp)==0) {
      int verb=getVerbose();
      if(verb>1)
	cerr << "Visualiser: Building grid structure...(" << pSim->getVersion() << ")";
      setVerbose(0);
      
      pGrid = dynamic_cast<Grid*>(pSim);
      if(pGrid==NULL)
	pGrid = new Grid(*pSim,3,3,3,10,numPartTarget/5);
      setVerbose(verb);
      if(verb>1)
	cerr << "done!"  <<  endl;

      latestUpdate = 1;
      initialGridWalk();
    } else {
      // not using a grid section, just naively splitting into lists
      updateList = 0;
    }

    oldVersion=pSim->getVersion();

    invalidate();
    
  }


  void Visualiser::drawMetaData(SimSnap *sim, float mindelta) {
  
    if(sim==NULL) return;
    if(typeid(*sim)==typeid(Grid)) {
      Grid * ourGrid = static_cast<Grid*>(sim);
      float dx = ourGrid->getDx(), dy = ourGrid->getDy(), dz = ourGrid->getDz();
      float x1 = ourGrid->getX1(), y1 = ourGrid->getY1(), z1 = ourGrid->getZ1();
      int   nx = ourGrid->getNx(), ny = ourGrid->getNy(), nz = ourGrid->getNz();
      float x2 = x1+nx*dx, y2 = y1+ny*dy, z2 = z1+nz*dz;

      int x,y,z;

      glColor4f(1,1,1,0.7);
    
      for(x=0;x<=nx;x++) {
	for(y=0;y<=ny;y++) {
	  glBegin(GL_LINES);
	  glVertex3f(x1+x*dx,y1+y*dy,z1);
	  glVertex3f(x1+x*dx,y1+y*dy,z2);
	  glEnd();    
	}
      }
    
      for(x=0;x<=nx;x++) {
	for(z=0;z<=nz;z++) {
	  glBegin(GL_LINES);
	  glVertex3f(x1+x*dx,y1,z1+z*dz);
	  glVertex3f(x1+x*dx,y2,z1+z*dz);
	  glEnd();    
	}
      }
    
      for(z=0;z<=nz;z++) {
	for(y=0;y<=ny;y++) {
	  glBegin(GL_LINES);
	  glVertex3f(x1,y1+y*dy,z1+z*dz);
	  glVertex3f(x2,y1+y*dy,z1+z*dz);
	  glEnd();    
	}
      }

      
          
      for(x=0;x<nx;x++) {
	for(y=0;y<ny;y++) {
	  for(z=0;z<nz;z++) {
	    SimSnap *consider = &((*ourGrid)[x][y][z]);
	    
	    if(typeid(*consider)==typeid(Grid)) {
	      if(static_cast<Grid*>(consider)->getDx()>mindelta)
		drawMetaData(consider,mindelta);
	    } /*else {
	      float opacity = consider->getNumParticles()/1000.;
	      opacity*=trans/(dx);
	      glColor4f(1.,1.,1.,opacity);
	      glTranslatef(x1+x*dx+dx/2.,y1+y*dy+dy/2.,z1+z*dz+dz/2.);
	      glutSolidCube(dx);
	      glTranslatef(-x1-x*dx-dx/2.,-y1-y*dy-dy/2.,-z1-z*dz-dz/2.);
	      
	      }*/
	  }
	}
      }
      
      
    }
  
  }

  bool Visualiser::displayIncomplete() {
    return (!noUpdatesLastDisplay);
  }

  void Visualiser::invalidate() {
    if(pSim->getVersion()!=oldVersion)
      fullReset();
    else {
      if(pGrid==NULL) {
	
	if(updateList>0) restartUpdateAtEnd = true;
	else updateList=1;
      } else {
	
	latestUpdate++;
      }
    }
  }

  void Visualiser::doScale(float ratio) {
      
    if((flags & lockZoom)>0) {
    	pointScale *=ratio;
    }

    glMatrixMode(GL_MODELVIEW);
  
    GLdouble tranMatrix[16];
  
    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
  
    GLdouble transx=tranMatrix[12],transy=tranMatrix[13],transz=tranMatrix[14];
  
    tranMatrix[3] = 0;
    tranMatrix[7] = 0;
    tranMatrix[11] = 0;
    tranMatrix[12] = 0;
    tranMatrix[13] = 0;
    tranMatrix[14] = 0;
    tranMatrix[15] = 1; // probably is always 1 - best to be sure...
  
    // Now do post-multiplying stuff:
    glLoadIdentity();
  
    // perform rotation about point axesOffset away
    glTranslated( 0,0,axesOffset);
  
    glScalef(ratio,ratio,ratio);
    // move back to correct location
    glTranslated( transx, transy, transz-axesOffset);
  
    glMultMatrixd(tranMatrix);

  
  }

  void Visualiser::screenMove(float mx, float my, float mz) {
    GLdouble tranMatrix[16];
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
    glLoadIdentity();
    glTranslatef(mx,my,mz);
    glMultMatrixd(tranMatrix);
    glutPostRedisplay();
  }

  void Visualiser::screenRotate(float rx, float ry) {
    glMatrixMode(GL_MODELVIEW);

    // we want to pre-multiply by the following matricies:
  
    GLdouble multMatrixX[16] = {cos(rx),0,-sin(rx),0,0,1,0,0,sin(rx),0,cos(rx),0,0,0,0,1};
 
    GLdouble multMatrixY[16] = {1,0,0,0,0,cos(ry),sin(ry),0,0,-sin(ry),cos(ry),0,0,0,0,1};

    // but glMultMatrix post-multiplies, so firstly get current matrix:

    GLdouble tranMatrix[16];

    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);

    // Remove translation information, which we do not want to rotate

    GLdouble transx=tranMatrix[12],transy=tranMatrix[13],transz=tranMatrix[14];
 
    tranMatrix[3] = 0;
    tranMatrix[7] = 0;
    tranMatrix[11] = 0;
    tranMatrix[12] = 0;
    tranMatrix[13] = 0;
    tranMatrix[14] = 0;
    tranMatrix[15] = 1; // probably is always 1 - best to be sure...
  
    // Now do post-multiplying stuff:
    glLoadIdentity();

    // perform rotation about point axesOffset away
    glTranslated( 0,0,axesOffset);

    glMultMatrixd(multMatrixX); 

    glMultMatrixd(multMatrixY);
  
    // move back to correct location
    glTranslated( transx, transy, transz-axesOffset);

    glMultMatrixd(tranMatrix);
  
  }

  void Visualiser::setColourMap(ColourMap & cm) {
    pCM = &cm;
    pCM->setReferenceAlpha(trans);
    invalidate();
  }

  ColourMap & Visualiser::getColourMap() {
    return *pCM;
  }


  void Visualiser::setSimSnap(SimSnap & ss) {
    pSim = &ss;
    fullReset();
  }

  SimSnap & Visualiser::getSimSnap() {
    return *pSim;
  }

  void Visualiser::rasterText(float x, float y, string text, void* fonttype) {

    glRasterPos3f(x,y,0);
    
    int len = text.length();
    for(int i=0; i<len; i++) {
      glutBitmapCharacter(fonttype, text[i]);
    }
  }


  void Visualiser::rasterText(string text) {
    int len = text.length();
    for(int n=0; n<len; n++) {
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,(int) text[n]);
    }
  }


  float Visualiser::getApproxLengthScale() {
    // returns an approximate length scale in internal simulation units, for the current view
    
    float tranMatrix[16];

    glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);

    float sc = sqrt(tranMatrix[0]*tranMatrix[0] + tranMatrix[4]*tranMatrix[4] + tranMatrix[8]*tranMatrix[8]);

    //  float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
    // det = pow((double)det,(double)1./3.);
  
    if((flags&ortho)==0)
      return 120./sc;
    else
      return 20./sc;
  }

  void Visualiser::setLengthScale(float ls) {
    float ratio = ls/getApproxLengthScale();
    doScale(ratio);
  }

  SimanVec Visualiser::getCameraPos() {
    vector<float> r;
    float x,y,z;
    getCameraPos(x,y,z);
    
    return SimanVec(x,y,z);

  }

  void Visualiser::setCameraPos(const SimanVec &pos) {
    
    double tranMatrix[16];

    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);

    // to get camera coordinates, need to invert model view matrix
    // (3x3 part) and multiply by camera coords stored in final column

    // note we assume matrix is orthogonal (apart from magnification part) to calculate this inverse...

    double det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
    det = std::pow((double)det,(double)1./3.);

    tranMatrix[14]-=axesOffset;

    
    tranMatrix[12]=0-(tranMatrix[0]*pos[0]+tranMatrix[4]*pos[1]+tranMatrix[8]*pos[2]);
    tranMatrix[13]=0-(tranMatrix[1]*pos[0]+tranMatrix[5]*pos[1]+tranMatrix[9]*pos[2]);
    tranMatrix[14]=axesOffset-(tranMatrix[2]*pos[0]+tranMatrix[6]*pos[1]+tranMatrix[10]*pos[2]);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(tranMatrix);

    
  
  }

  void Visualiser::getCameraPos(float &camx, float &camy, float &camz)
  {
  
    float tranMatrix[16];

    glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);

    // to get camera coordinates, need to invert model view matrix
    // (3x3 part) and multiply by camera coords stored in final column

    // note we assume matrix is orthogonal (apart from magnification part) to calculate this inverse...

    float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
    det = std::pow((double)det,(double)2./3.);

    tranMatrix[14]-=axesOffset;
    camx=-(tranMatrix[0]*tranMatrix[12]+tranMatrix[1]*tranMatrix[13]+tranMatrix[2]*tranMatrix[14])/(det);
    camy=-(tranMatrix[4]*tranMatrix[12]+tranMatrix[5]*tranMatrix[13]+tranMatrix[6]*tranMatrix[14])/(det);
    camz=-(tranMatrix[8]*tranMatrix[12]+tranMatrix[9]*tranMatrix[13]+tranMatrix[10]*tranMatrix[14])/(det);
  
  }

  vector<float> Visualiser::getRotationMatrix() {
  
    float tranMatrix[16];

    glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
  
    float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
    det = std::pow((double)det,(double)1./3.);
  
    vector<float> rotMat;
    for(int c=0; c<3; c++)
      for(int r=0; r<3; r++) 
	rotMat.push_back(tranMatrix[c*4+r]/det);

    return rotMat;

  }



  void Visualiser::vertexOnAxis(char which, float along, float below) {
    switch(which) {
    case 'x':
      glVertex3f(along,-below,0);
      break;
    case 'y':
      glVertex3f(below,along,0);
      break;
    case 'z':
      glVertex3f(0,-below,along);
      break;
    }
  }


  void Visualiser::rasterPosOnAxis(char which, float along, float below) {
    switch(which) {
    case 'x':
      glRasterPos3f(along,-below,0);
      break;
    case 'y':
      glRasterPos3f(below,along,0);
      break;
    case 'z':
      glRasterPos3f(0,-below,along);
      break;
    }
  }

  void Visualiser::drawAxis(char which, float max, int nticks) {

    // calculating overall scaling into "kilo"units, "Mega"units, etc...

    float logscale = log(max)/log(10.);
    int mult=0;

    if(logscale>=2 && logscale <5)
      mult = 3;
    if(logscale>=5)
      mult = 6;

    float rel_scale = std::pow(10.,-mult);

  
    glBegin(GL_LINES);
    vertexOnAxis(which,-max*unit_scaler,0);
    vertexOnAxis(which,max*unit_scaler,0);
    for(int n=-nticks; n<=nticks; n++) {
      if(n!=0) {
	vertexOnAxis(which,(float)n*max*unit_scaler/(float)nticks,0);
	vertexOnAxis(which,(float)n*max*unit_scaler/(float)nticks,max*unit_scaler/20);
      }
    }
    glEnd();
  
  
    stringstream labelstream(stringstream::out);

    for(int n=-nticks; n<=nticks; n++) {
      if(n!=0) {
	rasterPosOnAxis(which, (float)n*max*unit_scaler/(float)nticks,max*unit_scaler/20);
	labelstream << rel_scale*(float)n*max/(float)nticks;
	rasterText(labelstream.str());
	labelstream.str("");
      }
    }


    // overall axis label:

    labelstream << which;

    labelstream << "/";

    switch(mult) { 
    case 3:
      labelstream << 'k';
      break;
    case 6:
      labelstream << "M";
      break;
    }

    labelstream << units;
 
    rasterPosOnAxis(which,(max+max/(float)(nticks*2))*unit_scaler,0);
    rasterText(labelstream.str());
  
  }


  void Visualiser::drawAxes() {
  
    float range = getApproxLengthScale()/(unit_scaler*3); // approximate range of axes in visualiser "unit" units.
    float logrange = log(range)/log(10.);
  
    float snap_range = std::pow(10.,(int)logrange);

    if(snap_range>range*5) snap_range/=5;
    if(snap_range>range*2) snap_range/=2;
    if(snap_range<range/2) snap_range*=2;
    if(snap_range<range/5) snap_range*=5;

    glColor4f(1,1,1,0.7);

    // we always want the axes to be axesOffset in front of the camera

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    GLdouble tranMatrix[16];

    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
  
    // zero translation information
    tranMatrix[3] = 0;
    tranMatrix[7] = 0;
    tranMatrix[11] = 0;
    tranMatrix[12] = 0;
    tranMatrix[13] = 0;
    tranMatrix[14] = 0;
    tranMatrix[15] = 1; // probably is always 1 - best to be sure...

    glLoadIdentity();
    glTranslatef(0.,0.,axesOffset);
    glMultMatrixd(tranMatrix);

    drawAxis('x',snap_range,5);
    drawAxis('y',snap_range,5);
    drawAxis('z',snap_range,5);

    glPopMatrix();

  }

  void Visualiser::setTrans(float t)
  {
    trans=t;
    
    if (pCM != NULL) {
      pCM->setReferenceAlpha(trans);
    }
    if(!use_fragshader) {
      invalidate();
      // no fragshader - transparency has to be updated in individual particles
    } else {
      redisplay();
    }    
  }

  float Visualiser::getTrans() {
    return trans;
  }

  void Visualiser::setVectorScale(float v) {

    vectorScale = v;
    if(!use_fragshader) {
      invalidate();
    } 
  }

  
  float Visualiser::getVectorScale() {

    return vectorScale;
  }

  void Visualiser::renderSection(SimSnap *pRender, int start, int end, float lensc) {
    noUpdatesLastDisplay=false;

    float col_dm[4]   = {0., 1., 0., trans};
    float col_gas[4]  = {0., 0., 1., trans};
    float col_star[4] = {1., 1., 0., trans};
    float col_unknown[4] = {1., 1., 1., trans};
    float *col=col_unknown;
  
    if((flags&vectors)!=0) {
      glBegin(GL_LINES);
    } else if((flags&splatting)!=0 && use_fragshader) {
      glBegin(GL_QUADS);
    } else
      glBegin(GL_POINTS);

    GLint typeInShade=0;
    GLint velInShade=0;
    GLint deltaInShade=0;
    GLint lenscInShade=0;
    GLint massInShade = 0;
    if(use_fragshader) {
      typeInShade = glGetAttribLocationARB(shader_prog,"type");
      deltaInShade = glGetAttribLocationARB(shader_prog,"delta");
      velInShade = glGetAttribLocationARB(shader_prog,"vel");
      massInShade = glGetAttribLocationARB(shader_prog,"mass");
      lenscInShade = glGetAttribLocationARB(shader_prog,"lensc");
    }
    for(int n=start;n<end;n++) {
    
      bool plot = true;
    
      const Particle *p = pRender->getConstParticle(n);
    
      if(pCM!=NULL) {
      
	// colour map present
      
	if (!(*pCM)(pRender->deReference(n,pSim),col_unknown,p)) plot=false;
      
      } else {
      
	// simple default colour scheme, colours by
	// particle type
      
	switch (p->type)
	  {
	  case Particle::dm:
	    col = col_dm;
	    break;
	  case Particle::gas:
	    col = col_gas;
	    break;
	  case Particle::star:
	    col = col_star;
	    break;
	  default:
	    col = col_unknown;
	    break;
	  }
      } // pCM==NULL
    
      if(use_fragshader) {
	col[3]=1.; // transparency information will be added in GPU hardware
	glVertexAttrib1fARB(typeInShade,(p->type==Particle::dm)?1.2:0+(p->type==Particle::gas)?2.2:0+(p->type==Particle::star)?3.2:0);
	glVertexAttrib1fARB(lenscInShade,lensc);
      }

      if(plot) {
	     
	if((flags&vectors)==0) {
	  
	  glColor4fv(col);
	  if((flags&splatting)!=0 && use_fragshader) {
	    glVertexAttrib1fARB(massInShade,p->mass/maxmass);
	    
	    glVertexAttrib2fARB(deltaInShade,-1,-1);
	    glVertex3f(p->x,p->y,p->z);
	    glVertexAttrib2fARB(deltaInShade,-1,+1);
	    glVertex3f(p->x,p->y,p->z);
	    glVertexAttrib2fARB(deltaInShade,+1,+1);
	    glVertex3f(p->x,p->y,p->z);
	    glVertexAttrib2fARB(deltaInShade,+1,-1);
	    glVertex3f(p->x,p->y,p->z);
	  } else glVertex3f(p->x,p->y,p->z);

	} else {
	  glColor4fv(col);
	  glVertex3f(p->x,p->y,p->z);
	  if(use_fragshader) {
	    glVertexAttrib3fARB(velInShade,p->vx,p->vy,p->vz);
	    glVertex3f(p->x,p->y,p->z);
	    glVertexAttrib3fARB(velInShade,0,0,0);
	  } else {
	    glVertex3f(p->x+p->vx*vectorScale,p->y+p->vy*vectorScale,p->z+p->vz*vectorScale);
	  }
	}
      }
    }

    glEnd();
  
  }

  int Visualiser::getListFor(SimSnap *s) {
    map<SimSnap*,int>::iterator i;
    i=gridCellListMap.find(s);
    if(i!=gridCellListMap.end())
      return (*i).second;
    
    // not already assigned
    lastGLList+=numberSplitsPerCell;
    gridCellListMap[s] = lastGLList;
    return lastGLList;
  }

  void Visualiser::createSubListCallsForGrid(Grid & s) {
    int listBase = getListFor(&s);
    int nx=s.getNx(), ny=s.getNy(), nz=s.getNz();
    for(int refine=0; refine<numberSplitsPerCell; refine++) {
      glNewList(listBase+refine,GL_COMPILE);
      for(int x=0; x<nx; x++)
	for(int y=0; y<ny; y++)
	  for(int z=0;z<nz;z++)
	    glCallList(getListFor(&(s[x][y][z]))+refine);
      glEndList();
    }
  }

  int Visualiser::updateBaseListForGrid(SimSnap &s, int detail_0, int detail_1, float lensc) {
    
    int listBase = getListFor(&s);
    int num_up=0;

    for(int r=detail_0; r<=detail_1; r++) {
      int uid = r + listBase;
      
      ModuloFilter f(numberSplitsPerCell,r);
      Subset *pPartRender = new Subset(s,f);
	
      glNewList(uid, GL_COMPILE);
      
      renderSection(pPartRender,0,pPartRender->getNumParticles(),lensc);
      glEndList();
      num_up+=pPartRender->getNumParticles();

      delete pPartRender;    
    }
    //      cout << "UBL: " << &s << " " << listBase << " " << s.getNumParticles() << endl;
    return num_up;
 
  }

  int Visualiser::updateListForGrid(SimSnap &s, int max_update, int detail, float lensc) {
   
 
    pair<int,int> & lastUpdate = requiresUpdateMap[&s];

    if(lastUpdate.first<latestUpdate) lastUpdate.second=0; // redraw ALL detail levels 0..detail!
    if(lastUpdate.second<detail) {

      noUpdatesLastDisplay=false;
      int updates=0;
      
      if(typeid(s)==typeid(Grid)) {
	
	//       	cout << "--- Enter Grid " << &s << " ---" << endl;
	Grid *g = static_cast<Grid*>(&s);
	int nx=g->getNx(), ny=g->getNy(), nz=g->getNz();
	updates=0;

	bool exit_loop_early=false;

	for(int x=0; x<nx && !exit_loop_early; x++)
	  for(int y=0; y<ny && !exit_loop_early; y++) 
	    for(int z=0;z<nz && !exit_loop_early; z++) {
	      if(updates>max_update)
		exit_loop_early=true;
	      else {
		if(typeid((*g)[x][y][z])==typeid(Grid)) {
		  updates+=updateListForGrid((*g)[x][y][z],max_update-updates,detail,1.);
		} else {
		  // update and mark base cell as being done to given detail level
		  pair<int,int> & localLastUpdate = requiresUpdateMap[&((*g)[x][y][z])];
		  if(localLastUpdate.first<latestUpdate) localLastUpdate.second=0;
		  
		  if(localLastUpdate.second<detail) {
		    updates+=updateBaseListForGrid((*g)[x][y][z],lastUpdate.second,detail,
						 (*g).getDx()/std::pow((float)g->getNumParticles(),(float)0.33333));
		    
		    localLastUpdate.first=latestUpdate;
		    localLastUpdate.second=detail;
		  }
		}
	      }

	    }

	
	if(!exit_loop_early) {
	  // mark this entire cell as being done to this detail level
	  lastUpdate.first=latestUpdate;
	  lastUpdate.second=detail;
	}
	// cout << "--- End Grid " << &s << " / " << updates << " ---"  << endl;
	return updates;
	
      } else {

	int num=updateBaseListForGrid(s,lastUpdate.second,detail,lensc);
	lastUpdate.first=latestUpdate;
	lastUpdate.second=detail;
	return num;
      }
      
    }

    return 0;
  }
  
  
 

  void Visualiser::initialGridWalk(SimSnap &s) {
    

    pair<int,int> a;
    a.first=-1;
    a.second=0;

    requiresUpdateMap[&s]=a;
    
    if(typeid(s)==typeid(Grid)) {
      
      Grid *g = static_cast<Grid*>(&s);
      createSubListCallsForGrid(*g);
      int nx=g->getNx(), ny=g->getNy(), nz=g->getNz();
      for(int x=0; x<nx; x++)
	for(int y=0; y<ny; y++)
	  for(int z=0;z<nz;z++) 
	    initialGridWalk((*g)[x][y][z]);
	 
    } else {
      getListFor(&s); // allocate a display list
    }
    
    
  }

  void Visualiser::initialGridWalk() {
    lastGLList=0;
    requiresUpdateMap.clear();
    gridCellListMap.clear();
   
    initialGridWalk(*pGrid);
  }

  void Visualiser::transformSim() {
    float camx, camy, camz;
    getCameraPos(camx,camy,camz);
    Translation translate(-camx,-camy,-camz);
    OrthoTrans transform(getRotationMatrix());
    pSim->transform(transform*translate);
   

    // reset everything, but keep scale of visual the same

    float tranMatrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
    float sc = sqrt(tranMatrix[0]*tranMatrix[0] + tranMatrix[4]*tranMatrix[4] + tranMatrix[8]*tranMatrix[8]);
    
    fullReset();
    initModelView();

    glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
    float sc2 = sqrt(tranMatrix[0]*tranMatrix[0] + tranMatrix[4]*tranMatrix[4] + tranMatrix[8]*tranMatrix[8]);
    
    doScale(sc/sc2);

  }

  void Visualiser::paint(string extraText) {

    if(use_fragshader) {
      enableVertexShader();
      if(pCM!=NULL) if(typeid(*pCM)==typeid(ColourMapByType)) {
	ColourMapByType *pCMbT = static_cast<ColourMapByType*>(pCM);
	GLint trans_p =glGetUniformLocationARB(shader_prog,"transvec");
	// the "type vector" - see vertex code for explanation
	glUniform3fARB(trans_p,pCMbT->getDMMap().getReferenceAlpha(),pCMbT->getGasMap().getReferenceAlpha(),pCMbT->getStarMap().getReferenceAlpha());

      } else {
	GLint trans_p =glGetUniformLocationARB(shader_prog,"transvec");
	glUniform3fARB(trans_p,trans,trans,trans);
      }
      glUniform1fARB(glGetUniformLocationARB(shader_prog,"vecscale"),vectorScale);
      if((flags & splatting)!=0) {
	float tranMatrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
	float sc2 = std::pow((float)(tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) - tranMatrix[4]*(tranMatrix[1]*tranMatrix[10]-tranMatrix[9]*tranMatrix[2]) +tranMatrix[8] * (tranMatrix[1]*tranMatrix[6] - tranMatrix[5]*tranMatrix[2])),(float)0.3333);
	glUniform1fARB(glGetUniformLocationARB(shader_prog,"detmodelview"),sc2);
      }
    }
    noUpdatesLastDisplay=true;

    int nParticles;

    static int temp = 0;
    
    temp++;

    // get some information about the modelview transformation:

    float camx, camy, camz;
    getCameraPos(camx,camy,camz);



  
    float v_ls = getApproxLengthScale()/4.; // approximate range of view in simulation internal units

    // clear frame
    glClearColor(0, 0, 0, 0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
    // setup pixel size
    if((flags & vectors)!=0) 
      glLineWidth(pointScale);
    else
      glPointSize(pointScale);
     
    nParticles = pSim->getNumParticles();  

    int numPartPlotted =0;
    int numPartUpdated = 0; 
    int numPartInView = 0;
      
    if((flags & noParticles)==0) 
      {

	// Perform particle rendering


	// the following are used when a grid is being rendered:
	

    
	if(pGrid==NULL) {
	  for(int n=1;n<=numberLists;n++) {
	    if(n!=updateList) glCallList(n);
	  }
      
	  if (updateList>0 )
	    {
	 	  
	      int start,end;
	  
	      start = (nParticles*(updateList-1))/(numberLists);
	      end =   (nParticles*updateList)/(numberLists);
	  
	      glNewList(updateList, GL_COMPILE);
	      renderSection(pSim,start,end,1.);
	      glEndList();
	      updateList++;
	      numPartUpdated+=end-start;
	  
	      if(updateList>numberLists) {
		if(restartUpdateAtEnd)
		  updateList=1; 
		else
		  updateList=0;
		restartUpdateAtEnd=false;
	    
	      }
	  
	    }
      
	} else {

	  
	  
	  auto_ptr<Union> un = pGrid->getRegion(camx-v_ls,camx+v_ls,camy-v_ls,camy+v_ls,camz-v_ls,camz+v_ls,v_ls/10.);

	  list<SimSnap*> &plotList(un.get()->getList());

	  
	  list<SimSnap*>::iterator i;

	
	  for(i=plotList.begin();i!=plotList.end();i++) {
	    numPartInView += (*i)->getNumParticles();
	  }

      
	  float tot = ((float)numPartTarget)/((float)numPartInView);
	  if(tot>1.) tot = 1;

	  
	  int usenum = (int)(numberSplitsPerCell*tot);
	  if(usenum<1) usenum=1;
	  
	  for(i=plotList.begin();i!=plotList.end();i++) {
	    int listBase = getListFor(*i);
	    int numPartUpdatedX = numPartUpdated;

	    numPartPlotted+=(int)((*i)->getNumParticles()*usenum/numberSplitsPerCell);
	  
	    
	    if(numPartUpdated<numPartUpdatedTarget) {
	      numPartUpdated+=updateListForGrid(**i,numPartUpdatedTarget-numPartUpdated,usenum,1.);

	    }

	    if(numPartUpdatedX<numPartUpdatedTarget || displayWhileUpdating) {
	      static const int lists[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
	      glListBase(listBase);
	      glCallLists(usenum,GL_INT,lists);
	    }
	 	   
	  }


	  if((flags & plotMeta)>0) {
	   
	    if(pGrid!=NULL) drawMetaData(pGrid,v_ls/10.);
	  }



	  

	  //cout << tot << endl;
      
	  // always have centre cell up to date (would be nice to order updating
	  // more user-oriented in general, but this is a quick fix)
      
	  //	  int base = getListFor(*pGrid[g_cx][g_cy][g_cz]);
	  //updateListForGrid(*pGrid[g_cx][g_cy][g_cz]);

	  /*
	  for(int offset=0;offset<numberSplitsPerCell*tot;offset++) {
	    int uid = offset*numberLists+base;
	    if(requiresUpdate[uid]<latestUpdate) {
	      updateListForGrid(*pGrid[g_cx][g_cy][g_cz]);
	      requiresUpdate[uid]=latestUpdate;
	      numPartUpdated+=(int)((*(static_cast<Grid*>(pSim)))[g_cx][g_cy][g_cz].getNumParticles()*tot);
	    }
	  }
	  */
	  
	} // if auto-simplifying using grid

	  
      } // if plotting particles


    // other optional things to plot

    
    if(use_fragshader) {
      disableVertexShader();
      
    }

    glLineWidth(1.);
    glPointSize(1.);
  
  
    if((flags & axes) > 0) 
      drawAxes();


    list<Annotate*>::iterator i;
    for(i=annotationList.begin();i!=annotationList.end();i++) {
      (*i)->plot(v_ls);
    }

    // Overplot informative things:

    if(pCM!=NULL && (flags & dispColourBar)!=0) {
      
      TwoDViewport();
      
      pCM->renderColourBar(extx, exty, 0, 1, 0, 40./exty);
      
      glColor4f(1,1,1,1);
      if(numPartUpdated>0) {
	noUpdatesLastDisplay=false;
	rasterText(0.01,0.90,"Updating...",GLUT_BITMAP_HELVETICA_18);
      }
      
      glColor4f(1,1,1,1);
      ostringstream ss;
      ss << "SimAn/Vis Cen: x=" << setprecision(3 ) << camx/unit_scaler << " y=" << camy/unit_scaler <<" z=" << camz/unit_scaler << " [" << units << "]";
      rasterText(0.01,1.0-14./(float)exty,ss.str());
      
      ss.str("");
      ss << "numPartTarget=" << numPartTarget << "; numPartInView=" << numPartInView <<"; numPartPlotted=" << numPartPlotted<< "; numPartUpdated=" << numPartUpdated << " ";
      rasterText(0.01,1.0-28./(float)exty,ss.str()+extraText);
      
      
      UnTwoDViewport();
    } else if((flags & promo)!=0) {
      TwoDViewport();
      glColor4f(1,1,1,1);
      rasterText(0.01,1.0-14./(float)exty,"Real Time Rendering with SimAn");
      rasterText(0.01,1.0-28./(float)exty,"www.ast.cam.ac.uk/~app26/siman");
      UnTwoDViewport();
    }
    
    
    glutSwapBuffers();
    
  }


  void Visualiser::TwoDViewport() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glViewport(0,0,extx,exty);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.,1.0,0,1.0);
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //  glDisable(GL_BLEND);
  }

  void Visualiser::UnTwoDViewport() {
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    //  glEnable(GL_BLEND);
    initBlendingAndFog();
  }

  void Visualiser::initBlendingAndFog() {
  
    if(glBlendEquation!=NULL) {
      if((flags & blendMax)!=0) {
	glBlendEquation(GL_MAX);
      } else {
	glBlendEquation(GL_FUNC_ADD);
      }
    }


    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);

    /*
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogf(GL_FOG_DENSITY, 0.35f);
    glFogf(GL_FOG_START, 30.0f);
    glFogf(GL_FOG_END, 80.0f);
    glEnable(GL_FOG);
    */
  }

  void Visualiser::initModelView() {

   
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );
    glTranslatef( 0.0, 0.0, axesOffset );


    float boxsize = pSim->getApparentBoxSize();
    doScale(80./boxsize); 

  
  }

  void Visualiser::initProjection() { 
    
   
    glViewport( 0, 0, extx, exty );

    float aspectRatio = (float)exty/(float)extx;
  
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
  
    if((flags & ortho)>0) {
      glOrtho(-3.0/aspectRatio,3.0/aspectRatio,-3.0,3.0,1.0,80.0);
    } else {
      glFrustum( -1.0, 1.0, -aspectRatio, aspectRatio, 1.0, 80.0 );
    }
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();


  

  }

  void Visualiser::annotateVector(SimanVec len) {
    annotationList.push_back(new AnnotateVector(SimanVec(0,0,0),len));
    
  }


  void Visualiser::annotateVector(SimanVec start,  SimanVec len) {
    annotationList.push_back(new AnnotateVector(start,len));
  }


  void Visualiser::annotatePlane(SimanVec normal) {
    annotationList.push_back(new AnnotatePlane(SimanVec(0,0,0),normal));
  }

  void Visualiser::annotatePlane(SimanVec centre, SimanVec normal) {
    annotationList.push_back(new AnnotatePlane(centre,normal));
  }

  void Visualiser::popAnnotate() {
    if(annotationList.size()>0) {
      list<Annotate*>::iterator i = annotationList.end();
      i--;
      delete *i;
      annotationList.erase(i);
    } else {
      cerr << "(No annotations)" << endl;
    }
  }

} // namespace siman

#endif // SIMAN_VIS
