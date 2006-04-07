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

#ifdef SIMAN_TIF
#include <tiffio.h>
#endif

#define PI 3.141592

using namespace std;


// ALERT, ALERT...
// GLOBAL VARIABLE APPROACHING
//
// not nice... for now we can only have one CVisualise
// object working at once!!  This is because of C callbacks...
// 
// could extend this later to have a list of CVisualises?

CVisualise *pCurVis;


// CALLBACKS:
void timerCallback(int value) {
  glutTimerFunc(30,&(timerCallback),0);
  pCurVis->tick();
}

void dragCallback(int x, int y) {
  pCurVis->drag(x,y);
}

void motionCallback(int x, int y) {
  pCurVis->motion(x,y);
}

void clickCallback(int button, int state, int x, int y) {
  pCurVis->click(button,state,x,y);
}

void keyCallback(unsigned char key, int x, int y) {
  pCurVis->key(key,x,y);
}

void specialCallback(int key, int x, int y) {
  pCurVis->specialKey(key,x,y);
}

void reshapeCallback(int extx, int exty) {
  pCurVis->reshape(extx, exty);
}

void updateCallback() {
  pCurVis->update();
}

void idleCallback() {
  pCurVis->idle();
}

void menuCallback(int v) {
  pCurVis->menuItem(v);
}

string CVisualise::className() {
  return "CVisualise";
}

unsigned int CVisualise::supports() {
  return Visualise;
}

bool CVisualise::references(CSimanObject *p) {
  if(p==pCM || p==pSim) return true;
  return CSimanObject::references(p);
}

CVisualise::CVisualise(CSimSnap *pSimi, unsigned int flagsi, CScripted *pCLI_i) {

  pSim=pSimi;
  flags = flagsi;
  requiresUpdate=NULL;
  pCM = NULL;
  fullReset();

  if(pCLI_i==NULL) {
    pCLI = new CScripted();
    pCLI->setNamedVar("vis",this);
  } else {
    pCLI = pCLI_i;
  }


  registerVal("trans",reinterpret_cast<void*>(&trans),valueTypeFlt);
  registerVal("scale_aim",reinterpret_cast<void*>(&scale_aim),valueTypeFlt);
  registerVal("scale_limiter",reinterpret_cast<void*>(&scale_aim),valueTypeFlt);
  registerVal("scale",reinterpret_cast<void*>(&scale_aim),valueTypeFlt);

}

void CVisualise::fullReset() {
  
  numberOfLists = 10;
  numberOfSplitsPerCell = 10;

  if(pSim->className()=="CGrid") {
    CGrid *pG = (CGrid*) pSim;
    numberOfLists = pG->getNx()*pG->getNy()*pG->getNz();
    if(requiresUpdate!=NULL)
      delete[] requiresUpdate;

    requiresUpdate = new int[numberOfLists*numberOfSplitsPerCell];
    for(int n=0; n<numberOfLists * numberOfSplitsPerCell; n++) {
      requiresUpdate[n]=0;
    }
    latestUpdate=1;
  }
  pCurVis = this;
  scale = 0.0005;
  scale_aim = 0.0005;
  auto_rx = 0;
  auto_ry = 0;
  velocity = 0;
  mouseDown = false;
  running = false;
  frame=0;
  trans  = 0.5;
  pointScale = 3;
  flagMenu = -1;
  pLocalScale = NULL;
  units = "pc";
  unit_scaler = 1.e-3;
  scale_limiter = 2;
  initialisedModelView = false;
  display=true;

}



void CVisualise::run(bool blocking) {

  // need to update lists
  updateList =1;

  // initialise GLUT

  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  //glutInitDisplayMode( GLUT_RGB );

  glutInitWindowPosition(100,100);
  glutInitWindowSize(400,400);

  char **argv = NULL;
  int arg = 0;

  glutInit(&arg,argv);

  glutCreateWindow("Siman: Visualisation");

  // setup GLUT callbacks.
  
  glutDisplayFunc(&(updateCallback));
  glutReshapeFunc(&(reshapeCallback));
  glutTimerFunc(30,&(timerCallback),0);
  glutIdleFunc(&(idleCallback));
  glutMotionFunc(&(dragCallback));
  glutPassiveMotionFunc(&(motionCallback));
  glutMouseFunc(&(clickCallback));
  glutKeyboardFunc(&(keyCallback));
  glutSpecialFunc(&(specialCallback));

  
  
  buildMenu();
  

  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
  
  running=true;
  glutMainLoop();

  running=false;
  initialisedModelView=false;
}

void CVisualise::idle() {

  
  pCLI->pollStdIn();

}

void CVisualise::addFlagMenu(string text, int forflag) {

  if((flags & forflag)>0)
    text = "(*)" + text;
  else
    text = "( )" + text;

  glutAddMenuEntry(text.c_str(),forflag);

}


void CVisualise::buildMenu() {
  
  if(flagMenu!=-1) glutDestroyMenu(flagMenu);

  flagMenu = glutCreateMenu(menuCallback);
  
  
  addFlagMenu("Fly-through",flyThrough);
  addFlagMenu("Auto-Rotate",autoRotate);
  addFlagMenu("Pixels",pixels);
  if((flags & pixels)>0) {
    addFlagMenu("Lock Pixel size with Zoom",lockZoom);
  } else {
    addFlagMenu("Scale Poly size with Local Scalelength",sphScaling);
  }
  addFlagMenu("No particles",noParticles);
  addFlagMenu("Show meta data",plotMeta);
  addFlagMenu("Ortho mode",ortho);
  addFlagMenu("Draw axes",axes);
  addFlagMenu("Vectors",vectors);
  addFlagMenu("Blend by Max",blendMax);
  addFlagMenu("Auto-Simplify",autoSimp);
  if(pCM!=NULL) {
    colMapMenu = pCM->buildMenu();
    glutAddSubMenu("Colourmap",colMapMenu);
  }
  mainMenu = flagMenu;
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  
}

void CVisualise::menuItem(int v) {
  if((flags & v)>0) {
    flags-=v;
  } else {
    flags+=v;
  }

  if(v==ortho) {
    reshape(extx,exty);
  }
  
  if( ((flags & pixels)==0) && ((flags & lockZoom)>0) ) {
    cerr << "CVisualise: disabling inappropriate flag lockZoom" <<endl;
    flags-=lockZoom;
  }
  
  
  if( ((flags & pixels)>0) && ((flags & sphScaling)>0) ) {
    cerr << "CVisualise: disabling inappropriate flag sphScaling" <<endl;
    flags-=sphScaling;
  }
  
  initForRun();
  buildMenu();
 
  glutPostRedisplay();
}

void CVisualise::putVisSphere(float x, float y, float z, float r, float *col, float camx, float camy, float camz) {

   // for efficiency, glBegin and glEnd lie outside
  // this routine and must have already been called!
  // (GL_POINTS for pixels, or GL_QUADS otherwise)
    
  
  if(col!=NULL) {
    glColor4fv(col);
  }


  if((flags & pixels)!=0) {
    
    // just put a pixel
    // 
   
    
    glVertex3f(x,y,z);
    
  } else {
    
    // put a polygon with normal towards the viewer

    
    float normalx=camx - x, normaly = camy - y, normalz = camz-z;
    
    
    float normal = sqrt(normalx*normalx+normaly*normaly+normalz*normalz);
    
    
    normalx/=normal;
    normaly/=normal;
    normalz/=normal;

    float ax=-(normalz+normaly)/normalx,ay=1,az=1;
    
    float normala = sqrt(ax*ax+2);
    ax/=normala;
    ay/=normala;
    az/=normala;
    
    float bx=normaly*az - normalz*ay, by = normalz*ax - normalx * az, bz = normalx * ay - normaly * ax;

    // renormalise

    
    ax*= r * pointScale;
    ay*= r * pointScale;
    az*= r * pointScale;
    bx*= r * pointScale;
    by*= r * pointScale;
    bz*= r * pointScale;
    


    glVertex3f(x-ax-bx,y-ay-by,z-az-bz);
    glVertex3f(x+ax-bx,y+ay-by,z+az-bz);
    glVertex3f(x+ax+bx,y+ay+by,z+az+bz);
    glVertex3f(x-ax+bx,y-ay+by,z-az+bz);    

  }


}


void CVisualise::drawMetaData(CSimSnap *sim) {
  
  if(sim==NULL) return;
  if(sim->className()=="CGrid") {
    CGrid * ourGrid = (CGrid*) sim;
    float dx = ourGrid->getDx(), dy = ourGrid->getDy(), dz = ourGrid->getDz();
    float x1 = ourGrid->getX1(), y1 = ourGrid->getY1(), z1 = ourGrid->getZ1();
    int   nx = ourGrid->getNx(), ny = ourGrid->getNy(), nz = ourGrid->getNz();
    float x2 = x1+nx*dx, y2 = y1+ny*dy, z2 = z1+nz*dz;

    int x,y,z;
    
    
    


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

    /*
    
    bool refined = false;
    
    for(x=0;x<nx;x++) {
      for(y=0;y<ny;y++) {
	for(z=0;z<nz;z++) {
	  CSimSnap *consider = (*ourGrid)[x][y][z];
	 
	  if(consider->className()=="CGrid") {
	    drawMetaData(consider);
	  } else {
	    float opacity = consider->getNumParticles()/1000.;
	    opacity*=trans/(dx);
	    glColor4f(1.,1.,1.,opacity);
	    glTranslatef(x1+x*dx+dx/2.,y1+y*dy+dy/2.,z1+z*dz+dz/2.);
	    glutSolidCube(dx);
	    glTranslatef(-x1-x*dx-dx/2.,-y1-y*dy-dy/2.,-z1-z*dz-dz/2.);
	    	    
	  }
	}
      }
    }
    
    */
  }
  
}


void CVisualise::invalidate() {
  if(pSim->className()!="CGrid") {
    if(updateList>0) restartUpdateAtEnd = true;
    else updateList=1;
  } else {

    latestUpdate++;
  }
}

void CVisualise::doScale(float ratio) {
      
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


void CVisualise::tick() {
 
  bool needsRender = false;
  
  if(scale_aim!=scale) {

    float xscale = scale;

    scale = (9.*scale + scale_aim)/10.;

    if(fabs((scale_aim-scale)/scale)<0.005) scale = scale_aim;

    if(scale/xscale>scale_limiter) {
     
      scale = xscale * scale_limiter;
    }

    if(scale/xscale<1./scale_limiter) {
     
      scale = xscale/scale_limiter;
    }
    doScale(scale/xscale);
   
    needsRender = true;
  }

  if((flags & flyThrough)!=0) {
    flyThroughAdvance();
  } else {
    if(auto_rx!=0 || auto_ry!=0) {
      dragRotate(auto_rx,auto_ry);
      needsRender = true;
    }
  }

  if((flags & writeTiff)!=0) {
    
    glutPostRedisplay();
    frame++;
    scale = 18.;
    scale_aim = 18.;

    ostringstream ss;
    
    ss << setw(4) << setfill('0') << frame << ".tif";
    
#ifdef TIFF_LIBRARY_PRESENT
    writeTiffFile(ss.str().c_str(),ss.str().c_str(),0,0,extx,exty,COMPRESSION_NONE);
#endif

  }

  if(needsRender || updateList>0) glutPostRedisplay();
    
}

void CVisualise::flyThroughAdvance() {

  float rx = auto_rx, ry = auto_ry;


  if(mouseDown && velocity<0.3) velocity+=0.01;
  if(!mouseDown && velocity>0) velocity-=0.005;


  glMatrixMode(GL_MODELVIEW);
    // we want to pre-multiply by the following matricies:
  
  GLdouble multMatrixX[16] = {cos(rx),0,-sin(rx),0,0,1,0,0,sin(rx),0,cos(rx),0,0,0,0,1};
 
  GLdouble multMatrixY[16] = {1,0,0,0,0,cos(ry),sin(ry),0,0,-sin(ry),cos(ry),0,0,0,0,1};

  // but glMultMatrix post-multiplies, so firstly get current matrix:

  GLdouble tranMatrix[16];

  glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);

  // N.B. will also be rotating translation information!
  
  // Now do post-multiplying stuff:
  
  glLoadIdentity();

  glTranslatef(0.,0.,velocity);

  glMultMatrixd(multMatrixX); 

  glMultMatrixd(multMatrixY);

  glMultMatrixd(tranMatrix);
  
  
  
  //if(velocity>0) velocity-=0.005;

  glutPostRedisplay();

}

void CVisualise::motion(int x, int y) {
  static int xy=0;
  if((flags & flyThrough)!=0) {
    auto_rx = 0.3 * (float) (x-(extx/2))/(float)(extx);
    auto_ry = 0.3 * (float) (y-(exty/2))/(float)(exty);
  }
  if(y>0.9*exty && xy<=0.9*exty) {
    glutSetMenu(colMapMenu);   
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }
  if(y<=0.9*exty && xy>=0.9*exty) {
    glutSetMenu(mainMenu);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }
  xy = y;
  
}

void CVisualise::drag(int x, int y) {

  if(y<=0.9*exty) {
    if(dragModif==0) {
      if((flags & flyThrough)!=0) {
	auto_rx = 0.3 * (float) (x-(extx/2))/(float)(extx);
	auto_ry = 0.3 * (float) (y-(exty/2))/(float)(exty);
      } else {
	float rx = 6. * (float) (x-xi)/(float)(extx);
	float ry = 6. * (float) (y-yi)/(float)(exty);
	
	if((flags & autoRotate) != 0) {
	  auto_rx = 3. * (float) (x-xi)/(float)(extx);
	  auto_ry = 3. * (float) (y-yi)/(float)(exty);
	} else {
	  dragRotate(rx,ry);
	  glutPostRedisplay();
	}
	
	
      
      }
    } else if(dragModif==GLUT_ACTIVE_CTRL) {
      
      	float mx = 30. * (float) (x-xi)/(float)(extx);
	float my = 30. * (float) (yi-y)/(float)(exty);
	
	dragMove(mx,my);
	glutPostRedisplay();
    }
    
    
    xi = x;
    yi = y;
    
  }
}

void CVisualise::dragMove(float mx, float my) {
  GLdouble tranMatrix[16];
  glMatrixMode(GL_MODELVIEW);
  glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
  glLoadIdentity();
  glTranslatef(mx,my,0.);
  glMultMatrixd(tranMatrix);
  glutPostRedisplay();
}

void CVisualise::dragRotate(float rx, float ry) {
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

void CVisualise::click(int button, int state, int x, int y) {
	
	// Windows oddities:

	if((button & 8)==8) button-=8;
	if((button & 16)==16) button-=16;

  if(y>0.9*exty) {
    if(pCM->click(button,state,(float)x/(float)extx,(float)y/(float)exty))
      invalidate();
    glutPostRedisplay();
    return;
  } else {

    int modif = glutGetModifiers();

	
    if(state==GLUT_DOWN && button==0) {
      auto_rx = 0.;
      auto_ry = 0.;
      xi = x;
      yi = y;
      mouseDown = true;
      dragModif = modif;
    }
    
    if(state==GLUT_UP && button==0) {
      mouseDown = false;
      dragModif = 0;
    }
    

    
    if(modif==GLUT_ACTIVE_CTRL) {
      if(button == 3 || button==5) {
	trans*=1.1;
	
	// if(trans>1) trans = 1;
	invalidate();
	glutPostRedisplay();
      }
      
      if(button ==4) {
	trans/=1.1;
	if(trans<0.000005) trans = 0.000005;
	invalidate();
	glutPostRedisplay();
      }
      
    } else if (modif==GLUT_ACTIVE_SHIFT) {
      if(button == 3 || button==5) {
	pointScale*=1.1;
	
	glutPostRedisplay();
      }
      
      if(button ==4) {
	pointScale/=1.1;
	
	glutPostRedisplay();
    }
      
    } else {
      if(button == 3 || button==5) {
	
	scale_aim*=1.3;
	glutPostRedisplay();
      }
      
      if(button == 4) {
	scale_aim/=1.3;
	glutPostRedisplay();
      }
    }
  }

}

#ifdef SIMAN_OPENDAY
extern CSimSnap *pOpenDay1;
extern CSimSnap *pOpenDay2;
int curmode = 0;
#endif

void CVisualise::key(unsigned char key, int x, int y) {
	
  GLdouble tranMatrix[16];

  if(keyBind.find(key)!=keyBind.end()) {
    string com = keyBind[key];
    cerr << com << endl << "Siman>";
    pCLI->doCommand(keyBind[key]);
    return;
  }

  switch (key)
  {
  case 'i':
    scale_aim*=1.3;
    glutPostRedisplay();
    break;
  case 'o':
    scale_aim/=1.3;
    glutPostRedisplay();
    break;
  case 'a':
    menuItem(axes);
    break;
    
  case '+':
    pointScale *= 1.1;
    glutPostRedisplay();
    break;
  case '-':
    pointScale /= 1.1;
    glutPostRedisplay();
    break;
    
  case 'T':
    trans*=1.1;
    if(trans>1.0) trans = 1.0;
    changeTrans(trans);
    glutPostRedisplay();
    break;
    
  case 't':
    trans/=1.1;
    if(trans<0.000005) trans = 0.000005;
    changeTrans(trans);
    glutPostRedisplay();
    break;
  case ';':
    glMatrixMode(GL_MODELVIEW);
   
    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
    
    glLoadIdentity();

    glTranslatef(0.,0.,1.);
    
    glMultMatrixd(tranMatrix);
    
    glutPostRedisplay();
    break;
  case ' ':
    interruptWait =true;
    break;
  case '/':
    
    glMatrixMode(GL_MODELVIEW);
   
    glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
    
    glLoadIdentity();

    glTranslatef(0.,0.,-1.);
    
    glMultMatrixd(tranMatrix);
    
    glutPostRedisplay();
    break;
  }

}


void CVisualise::specialKey(int key, int x, int y) {
	
  GLdouble tranMatrix[16];
  
  glMatrixMode(GL_MODELVIEW);
  
  glGetDoublev(GL_MODELVIEW_MATRIX,tranMatrix);
  
  
  
  switch (key)
  {
 
  case GLUT_KEY_UP:
    glLoadIdentity();
    glTranslatef(0.,-1.,0.);
    glMultMatrixd(tranMatrix);
    glutPostRedisplay();
    break;
 
  case GLUT_KEY_DOWN:
    glLoadIdentity();
    glTranslatef(0.,1.,0.);
    glMultMatrixd(tranMatrix);
    glutPostRedisplay();
    break;
 
  case GLUT_KEY_LEFT:
    glLoadIdentity();
    glTranslatef(1.,0.,0.);
    glMultMatrixd(tranMatrix);
    glutPostRedisplay();
    break;
 
  case GLUT_KEY_RIGHT:
    glLoadIdentity();
    glTranslatef(-1.,0.,0.);
    glMultMatrixd(tranMatrix);
    glutPostRedisplay();
    break;
 
  }

}

void CVisualise::colourise(CColourMap* cm) {

  pCM = cm;

}

void CVisualise::writeText(float x, float y, string text, void* fonttype) {
  glRasterPos2f(x,y);
  for(int i=0; i<text.length(); i++) {
    glutBitmapCharacter(fonttype, text[i]);
  }
}

float CVisualise::getApproxLengthScale() {
  // returns an approximate length scale for the current view
    
  float tranMatrix[16];

  glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);

  float sc = sqrt(tranMatrix[0]*tranMatrix[0] + tranMatrix[4]*tranMatrix[4] + tranMatrix[8]*tranMatrix[8]);

  //  float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
  // det = pow((double)det,(double)1./3.);
  
  return 80./sc;
}


void CVisualise::getCameraPos(float &camx, float &camy, float &camz)
{
  
  float tranMatrix[16];

  glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);

  // to get camera coordinates, need to invert model view matrix
  // (3x3 part) and multiply by camera coords stored in final column

  // note we assume matrix is orthogonal (apart from magnification part) to calculate this inverse...

  float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
  det = pow((double)det,(double)2./3.);

  tranMatrix[14]-=axesOffset;
  camx=-(tranMatrix[0]*tranMatrix[12]+tranMatrix[1]*tranMatrix[13]+tranMatrix[2]*tranMatrix[14])/(det);
  camy=-(tranMatrix[4]*tranMatrix[12]+tranMatrix[5]*tranMatrix[13]+tranMatrix[6]*tranMatrix[14])/(det);
  camz=-(tranMatrix[8]*tranMatrix[12]+tranMatrix[9]*tranMatrix[13]+tranMatrix[10]*tranMatrix[14])/(det);
  
}

void CVisualise::getModelRot(float &rotx, float &roty, float &rotz) {

  float tranMatrix[16];

  glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
  
  float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
  det = pow((double)det,(double)1./3.);
  
  roty = asin(tranMatrix[2]/det);
  rotx = asin(-tranMatrix[6]/(det*cos(roty)));
  rotz = acos(tranMatrix[0]/(det*cos(roty)));
  
}

void CVisualise::dumpTranMatrix() {
  
  float tranMatrix[16];

  glGetFloatv(GL_MODELVIEW_MATRIX,tranMatrix);
  
  float det = tranMatrix[0]*(tranMatrix[5]*tranMatrix[10]-tranMatrix[9]*tranMatrix[6]) + tranMatrix[4]*(tranMatrix[9]*tranMatrix[2]-tranMatrix[1]*tranMatrix[10]) + tranMatrix[8]*(tranMatrix[1]*tranMatrix[6]-tranMatrix[5]*tranMatrix[2]);
    
  det = pow((double)det,(double)1./3.);
  
  for(int c=0; c<3; c++)
    for(int r=0; r<3; r++) 
      cout << tranMatrix[c*4+r]/det << " ";

}



void CVisualise::vertexOnAxis(char which, float along, float below) {
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


void CVisualise::rasterPosOnAxis(char which, float along, float below) {
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

void CVisualise::rasterText(string text) {
  int len = text.length();
  for(int n=0; n<len; n++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,(int) text[n]);
  }
}

void CVisualise::drawAxis(char which, float max, int nticks) {

  // calculating overall scaling into "kilo"units, "Mega"units, etc...

  float logscale = log(max/unit_scaler)/log(10.);
  int mult=0;

  if(logscale>=2 && logscale <5)
    mult = 3;
  if(logscale>=5)
    mult = 6;

  float rel_scale = pow(10.,-mult);

  
  glBegin(GL_LINES);
  vertexOnAxis(which,-max,0);
  vertexOnAxis(which,max,0);
  for(int n=-nticks; n<=nticks; n++) {
    if(n!=0) {
      vertexOnAxis(which,(float)n*max/(float)nticks,0);
      vertexOnAxis(which,(float)n*max/(float)nticks,max/20);
    }
  }
  glEnd();
  
  
  stringstream labelstream(stringstream::out);

  for(int n=-nticks; n<=nticks; n++) {
    if(n!=0) {
      rasterPosOnAxis(which, (float)n*max/(float)nticks,max/20);
      labelstream << rel_scale*(float)n*max/unit_scaler/(float)nticks;
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
 
  rasterPosOnAxis(which,max+max/(float)(nticks*2),0);
  rasterText(labelstream.str());
  
}


void CVisualise::drawAxes() {
  
  float range = getApproxLengthScale()/4;
  float logrange = log(range)/log(10.);
  
  float snap_range = pow(10.,(int)logrange);

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

void CVisualise::changeTrans(float t)
{

}

void CVisualise::renderSection(CSimSnap *pRender, int start, int end) {
  

  float col_dm[4]   = {0., 1., 0., trans};
  float col_gas[4]  = {0., 0., 1., trans};
  float col_star[4] = {1., 1., 0., trans};
  float col_unknown[4] = {1., 1., 1., trans};
  float *col=col_unknown;
  
  if((flags&vectors)!=0) {
    glLineWidth(pointScale);
    glBegin(GL_LINES);
  } else if((flags&pixels)!=0) 
    glBegin(GL_POINTS);
  else
    glBegin(GL_QUADS);
  
  if (pCM != NULL) {
    pCM->setReferenceAlpha(trans);
  }	  
  
  for(int n=start;n<end;n++) {
    
    float r=0.2;
    
    bool plot = true;
    
    CParticle *p = pRender->getParticle(n);
    
    if(pCM!=NULL) {
      
      // colour map present
      
      if (!(*pCM)(pRender->deReference(n,1000),col_unknown,p)) plot=false;
      
    } else {
      
      // simple default colour scheme, colours by
      // particle type
      
      switch (p->type)
	{
	case CParticle::dm:
	  col = col_dm;
	  break;
	case CParticle::gas:
	  col = col_gas;
	  break;
	case CParticle::star:
	  col = col_star;
	  break;
	default:
	  col = col_unknown;
	  break;
	}
    } // pCM==NULL
    
    if(plot && !(flags&vectors) ) putVisSphere(p->x,p->y,p->z,r,col,0,0,0);
    if(plot && (flags&vectors)) {
      glColor4fv(col);
      glVertex3f(p->x,p->y,p->z);
      glVertex3f(p->x+p->vx*0.1*pointScale,p->y+p->vy*0.1*pointScale,p->z+p->vz*0.1*pointScale);
      
    }
    
  }

  
  glEnd();
  
}

void CVisualise::updateListForGrid(int uid, int x, int y, int z, int r) {

  CSimSnap * pRender=(*(static_cast<CGrid*>(pSim)))[x][y][z];
  glDeleteLists(uid,1);
  glNewList(uid, GL_COMPILE);
  CModuloFilter f(numberOfSplitsPerCell,r);
  CSubset *pPartRender = new CSubset(pRender,f);
  // cout << pPartRender->getNumParticles() << "\t" << pRender->getNumParticles() << endl;
  renderSection(pPartRender,0,pPartRender->getNumParticles());
  delete pPartRender;    
  glEndList();
}


void CVisualise::update() {

  if(!display) {
    // clear frame
    glClearColor(0.1, 0, 0.1, 1.);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TwoDViewport();

    writeText(0.01,0.90,"Please wait, loading data...",GLUT_BITMAP_HELVETICA_18);
        
    UnTwoDViewport();
 
    glutSwapBuffers();

    return;
  }

  int nParticles;
  CSimSnap *pRender=pSim;
  
  
  float camx, camy, camz;
  getCameraPos(camx,camy,camz);
  
  bool deleteAfterRender = false; // delete pRender after finishing
  

  flags |= pixels; // no such thing as non-pixel mode now

  float v_ls; // visual length-scale

  // clear frame
  glClearColor(0, 0, 0, 0);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // setup pixel size in case we are rendering with pixels
  glPointSize(pointScale);
  
  //cout << camx << " " << camy << " " << camz << endl;
  
  
  nParticles = pRender->getNumParticles();  

  int plotted_ptcls =0;
  int num_updates_this_cycle = 0;

  if((flags & noParticles)==0) 
  {


    // the following are used when a grid is being rendered:
    CGrid *pGrid = NULL;
    int g_nx, g_ny, g_nz;
    float g_dx, g_dy, g_dz;
    float g_x1, g_y1, g_z1;
    int g_cx, g_cy, g_cz;



    bool gridded = pRender->className()=="CGrid";
    

    if(!gridded) { 
      for(int n=1;n<=numberOfLists;n++) {
	if(n!=updateList) glCallList(n);
      }
      
      if (updateList>0)
	{
	  // find camera position in data coords, which is useful in a number of
	  // calculations
	  
	  int start,end;
	  
	  start = (nParticles*(updateList-1))/(numberOfLists);
	  end =   (nParticles*updateList)/(numberOfLists);
	  
	  glNewList(updateList, GL_COMPILE);
	  renderSection(pRender,start,end);
	  glEndList();
	  updateList++;
	  num_updates_this_cycle++;
	  
	  if(updateList>numberOfLists) {
	    if(restartUpdateAtEnd)
	      updateList=1; 
	    else
	      updateList=0;
	    restartUpdateAtEnd=false;
	    
	  }
	  
	}
      
    } else {
      flags|=autoSimp;

      pGrid = (CGrid*) pRender;
      
      // We now calculate the following:
      //
      // g_n* = number of cells on * axis
      // g_d* = cell spacing on * axis
      // g_*1 = position of lowest cell on * axis
      // g_c* = cell coordinate of our current visual centre

      g_nx = pGrid->getNx();
      g_ny = pGrid->getNy();
      g_nz = pGrid->getNz();
      
      g_dx = pGrid->getDx();
      g_dy = pGrid->getDy();
      g_dz = pGrid->getDz();
      
      g_x1 = pGrid->getX1();
      g_y1 = pGrid->getY1();
      g_z1 = pGrid->getZ1();


      g_cx = (int) ((camx-g_x1)/g_dx);
      g_cy = (int) ((camy-g_y1)/g_dy);
      g_cz = (int) ((camz-g_z1)/g_dz);

      // clip:

      if(g_cx<0) g_cx=0;
      if(g_cy<0) g_cy=0;
      if(g_cz<0) g_cz=0;
      if(g_cx>=g_nx) g_cx=g_nx-1;
      if(g_cy>=g_ny) g_cy=g_ny-1;
      if(g_cz>=g_nz) g_cz=g_nz-1;

      v_ls = getApproxLengthScale();

      if((flags & autoSimp)!=0) { 
	int x_ls = (int) (v_ls/g_dx) + 1;
	int y_ls = (int) (v_ls/g_dy) + 1;
	int z_ls = (int) (v_ls/g_dz) + 1;

	int min_x = g_cx-x_ls;
	int max_x = g_cx+x_ls;

	
	if(min_x<0) min_x=0;
	if(max_x>g_nx-1) max_x=g_nx-1;

	int min_y = g_cy-y_ls;
	int max_y = g_cy+y_ls;

	if(min_y<0) min_y=0;
	if(max_y>g_ny-1) max_y=g_ny-1;

	int min_z = g_cz-z_ls;
	int max_z = g_cz+z_ls;

	if(min_z<0) min_z=0;
	if(max_z>g_nz-1) max_z=g_nz-1;

	// old way to estimate simplification:
	// float tot = (g_dx*g_nx)/(10.*v_ls);
	// if(tot>1) tot=1;

	// new way  - is there a more efficient way of doing this?
	
	int npart_est = 0;
	int npart_targ = 400000;
	
	
	for(int plotx=min_x;plotx<=max_x;plotx++) {
	  for(int ploty=min_y;ploty<=max_y;ploty++) {
	    for(int plotz=min_z;plotz<=max_z;plotz++) {
	      npart_est +=( (*(static_cast<CGrid*>(pSim)))[plotx][ploty][plotz]->getNumParticles());

	    }
	  }
	}
	

	float tot = ((float)npart_targ)/((float)npart_est);
	if(tot>1.) tot = 1;
	if(tot<0.1) tot=0.1;
	//cout << tot << endl;

	// always have centre cell up to date (would be nice to order updating
	// more user-oriented in general, but this is a quick fix)

	for(int offset=0;offset<numberOfSplitsPerCell*tot;offset++) {
	  int uid = offset*numberOfLists+g_cx+g_nx*(g_cy+g_ny*g_cz);
	  if(requiresUpdate[uid]<latestUpdate) {
	    
	    updateListForGrid(uid,g_cx,g_cy,g_cz,offset);
	    requiresUpdate[uid]=latestUpdate;
	    num_updates_this_cycle+=(int)((*(static_cast<CGrid*>(pSim)))[g_cx][g_cy][g_cz]->getNumParticles()*tot);
	  }
	}


	for(int plotx=min_x;plotx<=max_x;plotx++) {
	  for(int ploty=min_y;ploty<=max_y;ploty++) {
	    for(int plotz=min_z;plotz<=max_z;plotz++) {
	      for(int offset=0;offset<numberOfSplitsPerCell*tot;offset++) {
		int uid = offset*numberOfLists+plotx+g_nx*(ploty+g_ny*plotz);
		int npart =(int)( (*(static_cast<CGrid*>(pSim)))[plotx][ploty][plotz]->getNumParticles()*tot);
		if(requiresUpdate[uid]<latestUpdate && num_updates_this_cycle<50000) {
		 
		  updateListForGrid(uid,plotx,ploty,plotz,offset);
		  requiresUpdate[uid]=latestUpdate;
		  num_updates_this_cycle+=npart;
		  
		}
		glCallList(uid);
		plotted_ptcls+=npart;
	      }
	    }
	  }
	}
	
      } else {
	
	for(int n=1;n<=numberOfLists;n++) {
	  if(n!=updateList) glCallList(n);
	}
      }
    }

	  
  }

  
  //  glDisable(GL_DEPTH_TEST);

  if((flags & plotMeta)>0) {
    glColor4f(1,1,1,trans);
    drawMetaData(pRender);
  }

  //  glEnable(GL_DEPTH_TEST);

  
  float rotx, roty, rotz;
  getModelRot(rotx,roty,rotz);

  TwoDViewport();
  if(pCM!=NULL) {
    pCM->renderColourBar(extx, exty);
  }
  if(num_updates_this_cycle>0) {
    writeText(0.01,0.90,"Updating...",GLUT_BITMAP_HELVETICA_18);
  }
  
  
  ostringstream ss;

  
  ss << "SIMAN: x=" << camx << " y=" << camy<<" z=" << camz << " rot=[" << rotx << ", " << roty << ", " << rotz << "] sl = " << v_ls << " / p.p. = " << plotted_ptcls << " / u.p. = " << num_updates_this_cycle;
  writeText(0.00,0.95,ss.str());
  ss.clear();
  
  UnTwoDViewport();

  if((flags & axes) > 0) 
    drawAxes();

  glutSwapBuffers();

}


void CVisualise::TwoDViewport() {
  glPushMatrix();
  glLoadIdentity();
  glViewport(0,0,extx,exty);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.,1.0,0,1.0);
  glDisable(GL_BLEND);
}

void CVisualise::UnTwoDViewport() {
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_BLEND);
}

CSimanObject * CVisualise::getMember(const string & member) {
  if(member=="sim")
    return pSim;
  if(member=="cm")
    return pCM;

  return CSimanObject::getMember(member);
}

void CVisualise::setMember(const string & member,CSimanObject * var) {
  if(member=="sim") {
    if(!var->supports(SimSnap))
      throw(CTypeError(member));
    pSim = (CSimSnap *) var;
    invalidate();
    
  } else if(member=="cm") {
    if(!var->supports(ColourMap))
      throw(CTypeError(member));
    pCM = (CColourMap*) var;
    invalidate();
  } else CSimanObject::setMember(member,var);
}

CSimanObject * CVisualise::dispatch(string command, istream *stream, CScripted *pS) {
  
  if(command=="open") {
    if(!running) { 
      run();
      cerr << "Siman> ";
    }
    return NULL;
  }
  if(command=="close") {
    if(running) {
      glutLeaveMainLoop();
      cerr << "\r" << endl;
    }
    return NULL;
  }
  
  if(command=="rotmat") {
    dumpTranMatrix();
    return NULL;
  }

  if(command=="bind") {
    string l;
    *stream >> l;
    if(l.size()!=1 || stream->eof()) {
      throw(CSyntaxError("bind [key] [command]"));
    }
    string r;
    getline(*stream,r);
    keyBind[l[0]]=r;
    return NULL;
  }
  if(command=="autort") {
    if(!stream->eof()) {
      *stream >> auto_rx >> auto_ry;
    }
    cerr << "auto_rx = " << auto_rx << "\tauto_ry = " << auto_ry << endl;
    return NULL;
  }
  if(command=="wait") {
    return NULL;
  }
      
  return CSimanObject::dispatch(command,stream,pS);
  
}

void CVisualise::initForRun() {
  
  if(pLocalScale==NULL && ((flags & sphScaling)>0)) {
    cerr << "CVisualise: calculating local particle sizes...";
    pLocalScale = (float*) malloc(sizeof(float)*pSim->getNumParticles());
    for(int n=0; n<pSim->getNumParticles(); n++) {
      CParticle *p = pSim->getParticle(n);
      pLocalScale[n] = pSim->localPPscale(p->x, p->y, p->z);
      pSim->releaseParticle(p);
    }
    cerr << "done!" << endl;
  }
  
  glEnable( GL_CULL_FACE );
  
  /*
  if((flags & blendMax)!=0) {
    glBlendEquation(GL_MAX);
  } else {
    glBlendEquation(GL_FUNC_ADD);
  }
  */

  glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  glEnable(GL_BLEND);

  
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogf(GL_FOG_DENSITY, 0.35f);
  glFogf(GL_FOG_START, 30.0f);
  glFogf(GL_FOG_END, 80.0f);
  glEnable(GL_FOG);
  
}

void CVisualise::initModelView() {
  
  cerr << "CVisualise::initModelView()" << endl;


  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );
  glTranslatef( 0.0, 0.0, axesOffset );
  
  doScale(scale);
}


void CVisualise::reshape( int nWidth, int nHeight )
{

  initForRun();

  axesOffset = -40.0;

  extx = nWidth;
  exty = nHeight;

  GLfloat fAspect = (GLfloat) nHeight / (GLfloat) nWidth;
  GLfloat fPos[ 4 ] = { 0.0f, 0.0f, 10.0f, 0.0f };
  
  
  glViewport( 0, 0, nWidth, nHeight );
  
    
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();

  if((flags & ortho)>0) {
    glOrtho(-3.0/fAspect,3.0/fAspect,-3.0,3.0,1.0,80.0);
  } else {
    glFrustum( -1.0, 1.0, -fAspect, fAspect, 1.0, 80.0 );
  }
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  
 
  if(!initialisedModelView) {
    // first resize call since window opened (or we are explicitly
    // asked to reset the modelview by user)
    initModelView();
    initialisedModelView=true;
  }
}

int CVisualise::writeTiffFile(const char *filename, const char *description,
  int x, int y, int width, int height, int compression)
{
#ifdef TIFF_LIBRARY_PRESENT
  TIFF *file;
  GLubyte *image, *p;
  int i;

  file = TIFFOpen(filename, "w");
  if (file == NULL) {
    return 1;
  }
  image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = image;
  for (i = height - 1; i >= 0; i--) {
    if (TIFFWriteScanline(file, p, i, 0) < 0) {
      free(image);
      TIFFClose(file);
      return 1;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  TIFFClose(file);
#endif
  return 0;
}


#endif // SIMAN_VIS
