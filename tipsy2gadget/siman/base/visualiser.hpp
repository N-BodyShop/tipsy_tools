// visualiser.hpp - part of SimAn Simulation Analysis Library
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









#ifndef __VISUALISE_H_INCLUDED

#define __VISUALISE_H_INCLUDED

#include "../base.hpp"

namespace siman {

  class ColourMap;
  class ContinuuousColourMap;
  class ColourMapByType;
  class Annotate;

/// OpenGL renderer for SimSnap
///
/// Requires a child class to set up windows, provide interactivity etc.
/// Currently this would be GlutInterface, but it is likely to be superceded by
/// a better interface soon. However, the base class Visualiser should remain largely
/// the same, so you can rely on member functions herein.

class Visualiser : public SimanObject {
  friend class ColourMap;
  friend class ContinuuousColourMap;
  friend class ColourMapByType;
public:
  Visualiser(SimSnap *pSimI, unsigned int flags=0);
  virtual ~Visualiser();
  
  virtual void fullReset(); ///< called if total reinitialization is required (e.g. because the simulation has been changed)

  void paint(std::string extraText=""); ///< render a frame

  /// rotate model by rx, ry radians where rx and ry are angles about x,y screen axes respectively
  void screenRotate(float rx, float ry);

  /// move model by mx, my, mz screen units
  void screenMove(float mx, float my, float mz=0);

  /// get visual length scale in simulation units
  float getApproxLengthScale();

  /// set visual length scale in simulation units
  virtual void setLengthScale(float);

  /// translate modelview
  /// @param x,y - translation in plane of view
  /// @param z - translation perpendicular to plane of view
  void translate(float x, float y, float z);

  /// annotate diagram with vector
  void annotateVector(SimanVec len);
  void annotateVector(SimanVec start, SimanVec len);
 
  /// annotate diagram with plane
  void annotatePlane(SimanVec normal);
  void annotatePlane(SimanVec centre, SimanVec normal);

  /// annotation manipulation
  void popAnnotate();

  /// get the rotation matrix (no scale factor), column major form (i.e. id=c*3+r)
  std::vector<float> getRotationMatrix(); 

    
  void setTrans(float t);
  float getTrans();
  
  void setVectorScale(float v);
  float getVectorScale();

  virtual void updateWidgets()=0; ///< pure virtual, in child classes this should be implemented to re-evaluate all widgets with correct data

  

  /// Get top-level colour map currently in use
  ColourMap & getColourMap();

  /// Replace top-level colour map. 
  /// @param[in] CM - must reference a ColourMap object which remains accessible until the Visualise is deleted or setColourMap is called again
  void setColourMap(ColourMap & CM);

  /// Get the data currently in use
  SimSnap & getSimSnap();

  /// Replace the data in use. Triggers a full reset of the Visualiser object. (Scale information, rotation, camera position etc are lost.)
  ///@param[in] snap - must reference a SimSnap object which remains accessible until the Visualise is deleted or setSimSnap is called again
  void setSimSnap(SimSnap & snap);

  /// Check whether a specific flag is set
  virtual bool hasFlag(unsigned int flag);

  /// Set a particular flag to given state
  virtual void setFlag(unsigned int flag, bool state);

  /// Toggle a flag
  virtual void toggleFlag(unsigned int flag);

  /// Called to indicate that all display lists must be updated (e.g. because simulation data has changed)
  virtual void invalidate();

  /// Called to indicate a repaint is necessary, but only from existing display lists
  virtual void redisplay();
 
  static const int autoRotate = 1; ///< if static mode, does rotation continue after mouse activity ceases?
  static const int writeTiff = 2;  ///< use tifflib to write a TIFF of each frame
  static const int lockZoom = 4;  ///< set = lock pixel & zoom sizes
  static const int noParticles = 8; ///< set = don't plot particles
  static const int plotMeta = 16; ///< set = plot metadata (e.g. grids)
  static const int ortho = 32; ///< set = orthogonal projection (no perspective)
  static const int axes = 64; ///< set = draw axes
  static const int vectors = 128; ///< set = draw vector field (currently velocity only)
  static const int blendMax = 256; ///< set = use maximum value of two pixels, rather than additive blending
  static const int noAutoSimp = 1024; ///< set = don't use grid structure to perform auto-simplification (not recommended)
  static const int dispColourBar = 2048; ///< set = show colour bar at bottom of screen
  static const int promo = 4096; ///< set = display a promotional message for SimAn!
  static const int splatting = 8192; ///< set = use splatting vertex shader
  static const int noContextMenu = 1<<14; ///< set = disable context menu

  /// Perform a transformation on the simulation so that the centering and rotation matches the current view
  /// View is then centred on new origin and rotation-reset so that position does not appear to change 
  void transformSim();

  void getCameraPos(float &x, float &y, float &z);  ///< Get camera position in units of current simulation
  SimanVec getCameraPos(); ///< Get camera position in units of current simulation
  virtual void setCameraPos(const SimanVec & pos); ///<Set camera position in units of current simulation
  
  bool displayIncomplete(); ///< returns true if the last render did not reflect the internal state well
  
protected:

  
  /// Draw meta data at resolution in simulation units given by mindelta (set negative to plot max resolution always)
  void drawMetaData(SimSnap *sim, float mindelta=-1.); 

  void doScale(float ratio);
 
  void drawAxes();

  void disableVertexShader();
  void enableVertexShader();
  
  void printInfoLog(GLuint obj);

  void renderSection(SimSnap *render, int start, int end, float lensc);

  /// update all GL display lists belonging to s or its subgrids. Returns number of particles updated. 
  /// @param s - simsnap which is likely to be a Grid
  /// @param maxParts - attempt not to process more than this number of particles (so that rendering does not take too long)
  /// @param detail - maximum detail level of current view (inclusive). More detailed lists will not be rendered.
  int updateListForGrid(SimSnap &s, int maxParts, int detail, float lensc);

  /// update the GL display list, plotting actual data (not references to subcells)
  /// @param s - simsnap of part to plot
  /// @param detail_0 - lower detail level to update from (inclusive)
  /// @param detail_1 - upper detail level to update to (inclusive)
  int updateBaseListForGrid(SimSnap &s, int detail_0, int detail_1, float lensc); 

  void createSubListCallsForGrid(Grid &g); ///< update the subgrid, using calls to subgrids


  /// temporarily set up a 2D viewport with x=[0,1] and y=[0,1] over the whole display plane
  void TwoDViewport();

  /// revert to 3D transformations
  void UnTwoDViewport();

  /// initialise MODELVIEW matrix
  void initModelView();

  /// draw some text at current GL raster position
  static void rasterText(std::string text);
  /// draw some text at given position with given fonttype (optional)
  /// @todo uses GLUT, whereas Visualise class should only rely on GL
  static void rasterText(float x, float y, std::string text, void* fonttype=GLUT_BITMAP_HELVETICA_12);

  /// set the raster position to correspond to a position along an axis
  /// @param which - 'x', 'y', or 'z' depending on axis of interest
  /// @param along - distance along the axis to position, in pre-GL-transformation (i.e. "physical") units (pSim->getDistanceUnits())
  /// @param below - distance below teh axis to position, also in pre-GL-transformation "physical" units
  void rasterPosOnAxis(char which, float along, float below);
  
  /// Execute glVertex for a position on an axis.
  /// @see rasterPosOnAxis() for parameters
  void vertexOnAxis(char which, float along, float below); 

  /// Draw a complete axis, and label it.
  /// @param which - 'x', 'y', or 'z' depending on axis of interest
  /// @param max - length of axis, in pre-GL-transformation "physical" (pSim->getDistanceUnits()) units
  /// @param nticks - number of ticks/numerical labels to render
  void drawAxis(char which, float max, int nticks);
  
  /// initialise GL_BLEND mode and GL_FOG
  void initBlendingAndFog();

  /// initialise GL_PROJECTION matrix
  /// @param sizex, sizey - size of viewport in pixels
  void initProjection();
  
  /// get the GL list base corresponding to the given SimSnap.
  /// above this number, there are numberSplitsPerCell lists available for distributing the
  /// plotting data for the given cell
  int getListFor(SimSnap *s);
  

  SimSnap *pSim;
  float maxmass; ///< maximum particle mass in PSim
  Grid *pGrid; ///< gridded version of pSim
  long oldVersion; ///< last known version of the simulation, to allow auto-updating on data changes

  float axesOffset; ///< distance from view plane to axes (about which we rotate, zoom, etc.)
  

  unsigned int flags; ///< stores flags

  unsigned int frame; ///< stores a frame number for writing frame files

  std::string units;       ///< string representing base units, which are obtained by multiplying the position vector by unit_scaler
  float unit_scaler;  ///< multiply by this quantity to obtain coordinates in units defined above
  
 

  float pointScale; ///< the global point size for rendering
  float vectorScale; ///< multiplier for vector display
  ColourMap *pCM; ///< the colourmap in use


  int extx, exty; ///< the size of the viewport

  bool noUpdatesLastDisplay; ///< true if the last frame paint did not result in any list updates


  int numPartTarget; ///< the number of particles to render in each frame
  int numPartUpdatedTarget; ///< the number of particles to update in each frame, when updating is necessary
  bool displayWhileUpdating; ///< set to true to display cells which are currently updating. This looks smoother but is bad for performance under certain conditions.
private:

  
  void initialGridWalk(); ///< ensure grid structures are reflected in GL display lists
  void initialGridWalk(SimSnap &s); ///< ensure structures are reflected in GL display lists for part of grid corresponding to s


  float trans; ///< the global transparency multiplier (gets passed to colourmaps when selecting RGBA values)
  
  bool use_fragshader; ///< true if fragment shader is available
  
  int updateList; ///< =0 when no GL list needs be changed, or >0 to index currently updating list
  bool restartUpdateAtEnd; ///<set to TRUE if the start of the update is now invalidated

  int numberLists; ///< number of GL lists to be used if no Grid is passed
  int numberSplitsPerCell; ///< if a Grid is passed, total number of lists = numberLists*numberSplitsPerCell, the former determined automatically (=nx*ny*nz)
  
  std::map<SimSnap*,std::pair<int,int> > requiresUpdateMap; ///< last update map. First of pair refers to value of latestUpdate at last update, second refers to maximum "layer" of detail (0..numberSplitsPerCell)
  int latestUpdate; ///< if requiresUpdateMap[sim]<latestUpdate, the related display list needs updating!

  std::map<SimSnap*,int> gridCellListMap; ///< map accessed by getListFor()
  
  std::list<Annotate*> annotationList;
  
  int lastGLList; ///< most recent GL list base to be assigned; slots above are free


  void setUpVertexShade();
  void deleteShaders();
  GLhandleARB shader_prog;
  GLuint shader_obj_vert;
  GLuint shader_obj_frag;

  



};

} // namespace siman
#endif // VISUALISE_H_INCLUDED
