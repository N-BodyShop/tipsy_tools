// glutinterface.hpp - part of SimAn Simulation Analysis Library
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


#ifndef __GLUTINTERFACE_H_INCLUDED

#define __GLUTINTERFACE_H_INCLUDED

#include "../base.hpp"

namespace siman {

#ifdef SIMAN_READLINE
  extern int readlineGlutHook();
#endif

class GlutInterface : public Visualiser {

#ifdef SIMAN_READLINE
  friend int readlineGlutHook();
#endif
  friend class ColourMap;
  friend class ContinuuousColourMap;
  friend class ColourMapByType;
public:
  GlutInterface(SimSnap *pSimI, unsigned int flags=0, Scripted *pCLI_i=NULL);
  virtual ~GlutInterface();
  
  void run(bool withCLI=true);

  virtual void fullReset();
  virtual void setFlag(unsigned int flag, bool v);

  virtual void updateWidgets();
  virtual void invalidate();
  /// @defgroup glut_callbacks
  /// @{
  void update();
  virtual bool tick(); ///< returns true if it wants to be called again (i.e. animation ongoing)
  void reshape(int x, int y);
  void drag(int x, int y);
  void motion(int x, int y);
  void click(int button, int state, int x, int y);
  void idle();
  void key(unsigned char key, int x, int y);
  void specialKey(int key, int x, int y);
  void menuItem(int v);
  /// @}

  void bind(std::string key, std::string py);
  void bind(unsigned char key, std::string py);
  void unBind(std::string key);
  void unBind(unsigned char key);

  void reinitiateTimer(); ///< called to ensure the timer is running - e.g. if a zoom operation starts (it may turn off if nothing is happening)
 
  void redisplay();

  void buildMenu();  

  void addFlagMenu(std::string text, int forflag);

  void setup(); ///< called to prepare for entry to main glut loop
  void poll(); ///< called to deal with one event in the glut queue
  void end(); ///< called after final poll

  virtual void setCameraPos(const SimanVec &pos); ///< override to perform smooth fly to new point
  virtual void setLengthScale(float ls); ///< override to perform smooth scale transition

protected:


  static int writeTiffFile(const char *fn, const char *desc, int x, int y, int width, int height, int compression);

  Scripted *pCLI; ///< pointer to command line interpreter (always python at the moment)
  bool deleteCLIonExit; ///< if true, issues delete pCLI in destructor

  int xi, yi; ///< initial positions at start of drag
  bool mouseDown; ///< true if left button is down
  int dragModif; ///< set to modif at start of drag (i.e. stores keyboard modifiers)

  bool running; ///< true if currently in main GLUT loop, otherwise false
  unsigned int frame;

  
  float scale; ///< the current scaling of the display
  float scale_aim; ///< the current scaling we are aiming for
  float scale_limiter; ///< the maximum ratio by which scale can change in one frame

  int pos_aim; ///< exceeds zero if currently flying to new position
  float pos_aim_x, pos_aim_y, pos_aim_z; ///< the current position we are aiming for

  float fps;

  bool initialisedModelView; ///< false if initModelView() needs to be called on window resize

  bool forceFrame; ///< true if the next tick must redraw the frame, even if there's no apparent reason to do so (i.e. pending update)

  // From GLUT:

  int mainMenu;
  int colMapMenu;
  int flagMenu;

  bool buildMenuOnNextFrame; ///< menu needs rebuilding, but were not in a state to do it when we realised
  float auto_rx, auto_ry; ///< automatic rotation about x,y per frame

  /// Key-binding map
  std::map<char,std::string> keyBind;

};


extern GlutInterface *pCurVis;
extern bool glutInitialised;
extern int dynamic_glut_menu_alloc;
  
} // namespace siman
#endif // GLUTINTERFACE_H_INCLUDED
