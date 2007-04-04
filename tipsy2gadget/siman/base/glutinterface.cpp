// glutinterface.cpp - part of SimAn Simulation Analysis Library
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


#ifdef SIMAN_TIF
#include <tiffio.h>
#endif

#ifdef SIMAN_PNG
#include <png.h>
#endif

#ifdef SIMAN_READLINE
#include <readline/readline.h>
#endif

#include <typeinfo>



namespace siman {

  // ALERT, ALERT...
  // GLOBAL VARIABLE APPROACHING
  //
  // not nice... for now we can only have one Visualiser
  // object working at once!!  This is because of C callbacks...
  // 
  // could extend this later to have a list of Visualisers?

  GlutInterface *pCurVis=NULL;
  bool glutInitialised=false;
  
  int dynamic_glut_menu_alloc = 100;
  
  map <int, pair<ColourMap*,int> > colourmap_menu_map;

#ifdef SIMAN_READLINE
  int readlineGlutHook() {
    if(pCurVis->pSim->getVersion()!=pCurVis->oldVersion) {
      cerr << "\r";
      pCurVis->fullReset();
      cerr << "[SimAn]>>>";
    }
    glutMainLoopEvent();
    return 0;
  }
#endif

  void colourMapMenuCallback(int v) {
    pair<ColourMap*, int> obj = colourmap_menu_map[v];
    (obj.first)->menuCallback(obj.second);
  }

  void resetMenuIDs() {
    colourmap_menu_map.clear();
    dynamic_glut_menu_alloc=100;
  }

  int allocMenuID(ColourMap *call, int reason) {  
    pair<ColourMap*, int> obj(call,reason);
    colourmap_menu_map[dynamic_glut_menu_alloc]=obj;
    dynamic_glut_menu_alloc++;
    return dynamic_glut_menu_alloc-1;
  }



  // CALLBACKS:
  void timerCallback(int value) {
  
    if(pCurVis->tick()) glutTimerFunc(30,&(timerCallback),0);
  
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

  GlutInterface::GlutInterface(SimSnap *pSimi, unsigned int flagsi, Scripted *pCLI_i)
    : Visualiser(pSimi,flagsi), pCLI(pCLI_i) {
    
    running=false;
    deleteCLIonExit=false;
    pos_aim=0;
  }

  GlutInterface::~GlutInterface() {
    if(deleteCLIonExit)
      delete pCLI;

    
    //~Visualiser();
  }
  
  void GlutInterface::reinitiateTimer() {
    
    glutTimerFunc(30,&(timerCallback),0);
  }

  void GlutInterface::fullReset() {
  

    pCurVis = this;
    mouseDown = false;
    
    frame=0;
    flagMenu = -1;
   
  
    scale = 1;
    scale_aim = 1;
    scale_limiter = 2;

    forceFrame=false;

    auto_rx = 0;
    auto_ry =0;
    
    reinitiateTimer();

    ;
    Visualiser::fullReset();
    ;

  }

  void GlutInterface::redisplay() {
    glutPostRedisplay();
  }


  void GlutInterface::poll() {
    glutMainLoopEvent();
  }

  void GlutInterface::setup() {
    
    if(!glutInitialised) {
    
      // initialise GLUT
      glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
      glutInitWindowPosition(100,100);
      glutInitWindowSize(600,600);
  
      char **argv = NULL;
      int arg = 0;
      glutInit(&arg,argv);
      glutInitialised=true;

    }

    pCurVis=this;

    glutCreateWindow("Siman: Visualisation");
  

    // setup GLUT callbacks
    glutDisplayFunc(&(updateCallback));
    glutReshapeFunc(&(reshapeCallback));
    glutTimerFunc(30,&(timerCallback),0);
#ifndef SIMAN_READLINE
    glutIdleFunc(&(idleCallback));
#endif
    glutMotionFunc(&(dragCallback));
    glutPassiveMotionFunc(&(motionCallback));
    glutMouseFunc(&(clickCallback));
    glutKeyboardFunc(&(keyCallback));
    glutSpecialFunc(&(specialCallback));

  
    // build menu
    buildMenu();
  
    // want control back when the window is closed - only freeglut (not glut) does this
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    ;
   
    fullReset();
    initialisedModelView=false;

  }

  void GlutInterface::run(bool withCLI) {

    if(pCurVis!=NULL) {
      cerr << "GlutInterface: This version is unable to concurrently visualise two simulations" << endl;
      return;
    }
    setup();

    // initialize scripting interactivity if asked for
    if(pCLI==NULL && withCLI) {
      pCLI = new Scripted();
      deleteCLIonExit=true;
      pCLI->injectAsReference("vis",*this);
      cerr << "GlutInterface: Object `vis' created in shell" << endl;
    }  

#ifdef SIMAN_READLINE

    // implement hook
    if(pCLI!=NULL) {
      rl_event_hook = &readlineGlutHook;
      rl_set_keyboard_input_timeout(10);
      pCLI->mainLoop();
    } else {
      fullReset();
      glutMainLoop();
    }
#else
     // do main loop
    
    glutMainLoop();
#endif
    // deinitialize
    running=false;
    pCurVis=NULL;
    glutInitialised=false;
  }

  void GlutInterface::idle() {
    
    if(pCLI!=NULL)
      pCLI->pollStdIn();

    if(pSim->getVersion()!=oldVersion) {
      cerr << "\r";
      fullReset();
      cerr << "[SimAn]>>> ";
    }
  }

  void GlutInterface::addFlagMenu(string text, int forflag) {

    if((flags & forflag)>0)
      text = "(*)" + text;
    else
      text = "( )" + text;

    glutAddMenuEntry(text.c_str(),forflag);

  }


  void GlutInterface::invalidate() {
    if(!running) return;
    Visualiser::invalidate();
    updateWidgets();
    forceFrame=true;
    if(running) reinitiateTimer();
  }

  void GlutInterface::updateWidgets() {
    // cannot do this if we're in a menu callback
    reinitiateTimer();
    buildMenuOnNextFrame = true;

  }

  void GlutInterface::setCameraPos(const SimanVec & vec) {
    pos_aim_x=vec[0];
    pos_aim_y=vec[1];
    pos_aim_z=vec[2];
    pos_aim=1;
    reinitiateTimer();
  }
  
  void GlutInterface::setLengthScale(float ls) {
    scale_aim=scale*getApproxLengthScale()/ls;
    reinitiateTimer();
  }


  void GlutInterface::buildMenu() {

    resetMenuIDs();

    if(flagMenu!=-1) glutDestroyMenu(flagMenu);
    flagMenu = glutCreateMenu(menuCallback);
  
    addFlagMenu("Auto-Rotate",autoRotate);
    addFlagMenu("Lock Pixel size with Zoom",lockZoom);
    addFlagMenu("No particles",noParticles);
    addFlagMenu("Show meta data",plotMeta);
    addFlagMenu("Ortho mode",ortho);
    addFlagMenu("Draw axes",axes);
    addFlagMenu("Vectors (experimental)",vectors);
    addFlagMenu("Splatting (experimental",splatting);

    if(glBlendEquation!=NULL) {
      addFlagMenu("Blend by Max",blendMax);
    }

    addFlagMenu("Never Auto-Simplify data",noAutoSimp);
    addFlagMenu("Display colourbar",dispColourBar);

    if(pCM!=NULL) {
      colMapMenu = pCM->buildMenu();     
      glutAddSubMenu("Colourmap",colMapMenu);
    }

    mainMenu = flagMenu;

    if((flags&noContextMenu)==0)
      glutAttachMenu(GLUT_RIGHT_BUTTON);
  
  }

  void GlutInterface::setFlag(unsigned int flag, bool v) {
    Visualiser::setFlag(flag,v);
    if((flag&ortho)==ortho)
      reshape(extx,exty);
    if((flag&noContextMenu)==noContextMenu)
      buildMenu();
    glutPostRedisplay();
  }

  void GlutInterface::menuItem(int v) {
    toggleFlag(v);
  }



  bool GlutInterface::tick() {
 
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

    if(pos_aim>0) {
      float x, y, z;
      getCameraPos(x,y,z);
      x=(9.*x + pos_aim_x)/10.;
      y=(9.*y + pos_aim_y)/10.;
      z=(9.*z + pos_aim_z)/10.;
     
      
      if(pos_aim>40) {
	x = pos_aim_x;
	y = pos_aim_y;
	z = pos_aim_z;
	pos_aim=-1;
      }

      pos_aim++;


      Visualiser::setCameraPos(SimanVec(x,y,z));
      needsRender=true;
    }

    if((auto_rx!=0 || auto_ry!=0) && ((flags&autoRotate)!=0)) {
      screenRotate(auto_rx*20./fps,auto_ry*20./fps);
      needsRender = true;
    }

   

    
    if(needsRender || displayIncomplete() || forceFrame) {
      glutPostRedisplay();

      frame++;

      forceFrame=false;
      return true;
    } else return false;

   
  }

  void GlutInterface::motion(int x, int y) {
    static int xy=0;

    if((flags&noContextMenu)==0) {
      if(y>0.9*exty && xy<=0.9*exty) {
	glutSetMenu(colMapMenu);   
	glutAttachMenu(GLUT_RIGHT_BUTTON);
      }
      if(y<=0.9*exty && xy>=0.9*exty) {
	glutSetMenu(mainMenu);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
      }
    }

    xy = y;
  
  }

  void GlutInterface::drag(int x, int y) {

    if(y<=0.9*exty) {
      if(dragModif==0) {
	
	float rx = 6. * (float) (x-xi)/(float)(extx);
	float ry = 6. * (float) (y-yi)/(float)(exty);
	
	if((flags & autoRotate) != 0) {
	  reinitiateTimer();
	  auto_rx = rx/2.;
	  auto_ry = ry/2.;
	} else {
	  screenRotate(rx,ry);
	  glutPostRedisplay();
	}
	
      } else if(dragModif==GLUT_ACTIVE_CTRL) {
      
      	float mx = 30. * (float) (x-xi)/(float)(extx);
	float my = 30. * (float) (yi-y)/(float)(exty);
	
	screenMove(mx,my);
	glutPostRedisplay();
      }
    
    
      xi = x;
      yi = y;
    
    }
  }

  void GlutInterface::click(int button, int state, int x, int y) {
	
    // Windows oddities:

    if((button & 8)==8) button-=8;
    if((button & 16)==16) button-=16;

    if(y>0.9*exty) {
      if(pCM->click(button,state,(float)x/(float)extx,(float)y/(float)exty)) {
	invalidate();
	return;
      }
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
	  setTrans(getTrans()*1.1);
	  glutPostRedisplay();
	}
      
	if(button ==4) {
	  setTrans(getTrans()/1.1);
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
      
      } else if(modif==GLUT_ACTIVE_ALT) {
	if(button == 3 || button==5) {
	  setVectorScale(getVectorScale()*1.1);
	  glutPostRedisplay();

	}
      
	if(button ==4) {
	  setVectorScale(getVectorScale()/1.1);
	  glutPostRedisplay();

	}
      } else {
	if(button == 3 || button==5) {
	  scale_aim*=1.3;
	  reinitiateTimer();
	}
      
	if(button == 4) {
	  scale_aim/=1.3;
	  reinitiateTimer();
	}
      }
    }

  }

  void GlutInterface::bind(unsigned char key, string pycom) {
    keyBind[key]=pycom;
  }
  
  void GlutInterface::bind(string key, string pycom) {
    keyBind[key[0]]=pycom;
  }

  void GlutInterface::unBind(unsigned char key) {
    keyBind.erase(key);
  }

  void GlutInterface::unBind(string key) {
    keyBind.erase(key[0]);
  }

  void GlutInterface::key(unsigned char key, int x, int y) {
	
    
    if(keyBind.find(key)!=keyBind.end()) {
      string com = keyBind[key];
      pCLI->doCommand(keyBind[key]);
      return;
    }
    

    if((flags&noContextMenu)==0) {
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
	  setTrans(getTrans()*1.1);
	  glutPostRedisplay();
	  break;
    
	case 't':
	  setTrans(getTrans()/1.1);
	  glutPostRedisplay();
	  break;
	case ';':
	  screenMove(0,0,1);
	  glutPostRedisplay();
	  break;
     
	case '/':
	  screenMove(0,0,-1);
	  glutPostRedisplay();
	  break;
	}
    }

  }


  void GlutInterface::specialKey(int key, int x, int y) {
  
    switch (key)
      {
 
      case GLUT_KEY_UP:
	screenMove(0,-1,0);
	glutPostRedisplay();
	break;
 
      case GLUT_KEY_DOWN:
	screenMove(0,1,0);
	glutPostRedisplay();
	break;
 
      case GLUT_KEY_LEFT:
	screenMove(1,0,0);
	break;
 
      case GLUT_KEY_RIGHT:
	screenMove(-1,0,0);
	break;
 
      }

  }

  void GlutInterface::update() {
    if(buildMenuOnNextFrame) {
      buildMenu();
      buildMenuOnNextFrame=false;
    }

    static int xtime=0;
    if(!running) {
      initModelView();
      fullReset();
      running=true;
    }

    if((flags & writeTiff)==0) {
      int time = glutGet(GLUT_ELAPSED_TIME);
      fps = 1000.0/((float)(time-xtime));
      xtime = time;
    } else {
      fps = 25.0;
    }
    
    Visualiser::paint("FPS: "+boost::lexical_cast<string>(fps));
    glutSwapBuffers();

    
    if((flags & writeTiff)!=0) {
      
      
#if defined(SIMAN_PNG) || defined(SIMAN_TIF )
      
      ostringstream ss;
      
      ss << setw(4) << setfill('0') << frame;
      
#ifdef SIMAN_TIF
      ss << ".tif";
#else
      ss << ".png";
#endif
      
      cout << "frame = " << ss.str().c_str() << "                  \r"; 
      writeTiffFile(ss.str().c_str(),ss.str().c_str(),0,0,extx,exty,0);
#endif
      
      }
      
  }



  void GlutInterface::reshape( int width, int height )
  {

    initBlendingAndFog();

    axesOffset = -40.0;

    extx = width;
    exty = height;



 
    if(!initialisedModelView) {
      // first resize call since window opened (or we are explicitly
      // asked to reset the modelview by user)
      initModelView();
      initialisedModelView=true;
    }

    initProjection();
  }

  int GlutInterface::writeTiffFile(const char *filename, const char *description,
				   int x, int y, int width, int height, int compression)
  {
#ifdef SIMAN_TIF
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
    TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
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

#ifdef SIMAN_PNG
    FILE *fp = fopen(filename,"wb");
    if(fp==NULL) {
      cerr << "fopen fail" << endl;
      return 1;
    }
    png_structp png_ptr;
    png_infop info_ptr;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    info_ptr = png_create_info_struct(png_ptr);
    png_set_compression_level(png_ptr,Z_BEST_SPEED);

    if (info_ptr == NULL) {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
      cerr << "png_create_info_struct fail" << endl;
      return 1;
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
      cerr << "png_jmpbuf fail" << endl;
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return 1;
    }
    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr,info_ptr,width,height,8,PNG_COLOR_TYPE_RGB,
		 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
		 PNG_FILTER_TYPE_BASE);
    
    png_write_info(png_ptr,info_ptr);

    GLubyte *image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);

    png_bytep row_pointers[height];
    for (int k = 0; k < height; k++)
      row_pointers[k] = image + (height-k-1)*width*3;
 
    png_write_image(png_ptr, row_pointers);
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    free(image);

    fclose(fp); 
    
#endif
    return 0;
  }


} // namespace siman

#endif // SIMAN_VIS
