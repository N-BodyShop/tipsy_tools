//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// VISUALISE class - uses GLUT to visualise a CSimSnap



#ifndef __VISUALISE_H_INCLUDED

#define __VISUALISE_H_INCLUDED

#include "../base.hpp"


class CVisualise : public CSimanObject {
  friend class CColourMap;
  friend class CContinuuousColourMap;
  friend class CColourMapByType;
public:
  CVisualise(CSimSnap *pSimI, unsigned int flags=0, CScripted *pCLI_i=NULL);

  string className();
  virtual unsigned int supports();
  bool references(CSimanObject *p);
  
  void fullReset();

  void run(bool blocking=true);

  void drawMetaData(CSimSnap *sim);

  void update();
  void tick();
  void reshape(int x, int y);
  void drag(int x, int y);
  void motion(int x, int y);
  void click(int button, int state, int x, int y);
  void idle();
  void key(unsigned char key, int x, int y);
  void specialKey(int key, int x, int y);

  void doScale(float ratio);
  void dragRotate(float rx, float ry);
  void dragMove(float mx, float my);
  void menuItem(int v);
  void buildMenu();

  void colourise(CColourMap *pCMi);

  void flyThroughAdvance();
  static void writeText(float x, float y, string text, void* fonttype=GLUT_BITMAP_HELVETICA_12);

  void drawAxes();

  void invalidate();

  float getApproxLengthScale();
  void getCameraPos(float &x, float &y, float &z);
  void getModelRot(float &rx, float &ry, float &rz);

  void putVisSphere(float x, float y, float z, float r, float *col=NULL, float camx = 0, float camy = 0, float camz=0);
  
  void addFlagMenu(string text, int forflag);

  void changeTrans(float t);

  void renderSection(CSimSnap *render, int start, int end);

  void TwoDViewport();
  void UnTwoDViewport();
  void initModelView();

  void updateListForGrid(int uid, int x, int y, int z, int r);
  
  CSimanObject * dispatch(string command, istream *stream, CScripted *pS);
  CSimanObject * getMember(const string & member);
  void setMember(const string & member, CSimanObject *pOb);

  void dumpTranMatrix(); ///< temporary command. Send rotation matrix to cout.

protected:

  void rasterText(string text);
  void rasterPosOnAxis(char which, float along, float below);
  void vertexOnAxis(char which, float along, float below);
  void drawAxis(char which, float max, int nticks);
  
public:

  static const int flyThrough = 1; ///< fly-through or static rotate about (0,0,0)/
  static const int autoRotate = 2; ///< if static mode, does rotation continue after mouse activity ceases?
  static const int pixels = 4;     ///< if on, display GL points rather than spheres for each particle.
  static const int writeTiff = 8;  
  static const int lockZoom = 16;  // on = lock pixel & zoom sizes
  static const int noParticles = 32; // on = don't plot particles
  static const int plotMeta = 64; // on = plot metadata (e.g. grids)
  static const int ortho = 128; // on = orthogonal projection (no perspective)
  static const int sphScaling = 256; // on = scale by local particle-particle scale lengths
  
  static const int axes = 512; // on = draw axes
  static const int vectors = 1024; // on = draw vector field (currently velocity only)
  static const int blendMax = 2048; ///< on = use maximum value of two pixels, rather than additive blending
  static const int autoSimp = 4096; ///< on = automatically simplify regions by removing particles
private:

  
  int updateList; ///< =0 when no GL list needs be changed, or >0 to index currently updating list
  bool restartUpdateAtEnd; ///<set to TRUE if the start of the update is now invalidated

  int numberOfLists; ///< number of GL lists to be used if no CGrid is passed
  int numberOfSplitsPerCell; ///< if a CGrid is passed, total number of lists = numberOfLists*numberOfSplitsPerCell, the former determined automatically (=nx*ny*nz)
  int *requiresUpdate; ///< pointer to array storing update requirements for individual lists in Grid mode
  int latestUpdate; ///< if requiresUpdate[id]<latestUpdate, it needs updating!

  void initForRun();

  static int writeTiffFile(const char *fn, const char *desc, int x, int y, int width, int height, int compression);

  CSimSnap *pSim;
  CScripted *pCLI;

  float auto_rx, auto_ry;
  int xi, yi;

  int extx, exty;

  float scale; ///< the current scaling of the display
  float scale_aim; ///< the current scaling we are aiming for
  float scale_limiter; ///< the maximum ratio by which scale can change in one frame

  float velocity; ///< for fly-through mode
  CGrid *pGrid;   ///M also for fly-through, when only a fast-selecting subset of data will be displayed...

  bool mouseDown; ///< true if left button is down
  int dragModif; ///< set to modif at start of drag (i.e. stores keyboard modifiers)

  float axesOffset; ///< distance from view plane to axes (about which we rotate, zoom, etc.)

  bool running; ///< true if currently in main GLUT loop, otherwise false
  bool initialisedModelView; ///< false if initModelView() needs to be called on window resize

  bool display; ///< false if the data should not currently be accessed (e.g. during updates)

  unsigned int flags;

  unsigned int frame;

  string units;       /// string representing base units, which are obtained by multiplying the position vector by unit_scaler
  float unit_scaler;  /// multiply by this quantity to obtain coordinates in units defined above
  
  bool interruptWait; ///< set to false by dispatch() if a wait is called; set to true to interrupt that wait

  // changeable rendering constants

  float trans;
  float pointScale;


  // From GLUT:

  int mainMenu;
  int colMapMenu;
  int flagMenu;


  // Scaling data:

  float *pLocalScale;

  // Colourising data:
 
  CColourMap *pCM;

  /// Key-binding map
  map<char,string> keyBind;

};


#endif // VISUALISE_H_INCLUDED
