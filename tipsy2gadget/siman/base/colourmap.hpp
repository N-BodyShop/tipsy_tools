// colourmap.hpp - part of SimAn Simulation Analysis Library
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






#ifndef __COLOURMAP_H_INCLUDED

#define __COLOURMAP_H_INCLUDED

namespace siman {

class Visualiser;
	
class ColourMap : public SimanObject {
public:
  friend class ColourMapByType;

  ColourMap(Visualiser *pV);
  virtual ~ColourMap() { };

  virtual bool operator()(int index, float *gl4f, const Particle *p=NULL);

  virtual void autoRange();  // pick min and max values

  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05, bool renderColLabel=true);

  virtual void setReferenceColour(float r, float g, float b, float a);
  virtual void setReferenceAlpha(float a);
  float getReferenceAlpha();

  virtual bool click(int button, int state, float x, float y); 
       // returns true if changes have been made that require a rewrite of list

  virtual int buildMenu();

  virtual void menuCallback(int id);
  
  virtual std::string getCurrentLabel();

  virtual std::vector<double> getCol();
  virtual void setCol(const std::vector<double> &);
  
  virtual bool getDisplay();
  virtual void setDisplay(bool);
protected:
 
  bool display;
  std::string label;

  // simsnap
  Visualiser *pVis;
  ColourMap *pParent;

  // for interface

  int menuID;


  
  virtual void autoRange_Consider(unsigned int n);

  float refR, refG, refB, refA;

  bool autoranging;           
};

// 
// ContinuuousColourMap
//

class ContinuuousColourMap : public ColourMap {
public:
  ContinuuousColourMap(Visualiser *pV, float mini=0, float maxi=0);
  
  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05, bool renderColLabel=true);
  virtual bool operator()(int index, float *gl4f, const Particle *p=NULL);
  virtual bool operator()(float val, float *gl4f);
  
  virtual bool click(int button, int state, float x, float y); 
  virtual int buildMenu();

  virtual float particleToColourIndex(int n, const Particle *p=NULL);

  virtual void autoRange();
  virtual void autoRange_Consider(unsigned int n);

  virtual bool getHideOutOfRange();
  virtual void setHideOutOfRange(bool to);
  
  virtual bool getLogScale();
  virtual void setLogScale(bool to);


  virtual void menuCallback(int id);

  void setRange(std::vector<double> const &array);
  std::vector<double> getRange();
     
  void setColourBy(int colourBy);
  void setColourBy(std::string colourBy);

  std::string getColourLabel(int n);


  /// minimum and maximum values to map into colours
  float min;
  float max;

  /// if true, colours outside range [min,max] are not displayed; otherwise colour simply saturates
  bool hideOutOfRange;
  

  // for UI
  float x_button_down;

 
  /// \defgroup colmode Colouring Mode
  /// Definitions for determining what attribute to colour by
  /// @{
  int colourBy;  ///< defines what to colour by
  
  /// store details of extra data array if colouring by non-standard property
  const SimanArray *map;

  static const int useNe = 1;
  static const int useTemp = 2;
  static const int useRho = 0;
  static const int useU = 3;
  static const int useMetal = 4;
  static const int useZero = 5;
  static const int use_numOptions = 5;
  ///@}

  bool logScale;

};

class ColourMapSpectrum : public ContinuuousColourMap {
public:

  ColourMapSpectrum(Visualiser *pV, float, float);
  ColourMapSpectrum(Visualiser *pV);
  virtual bool operator()(float val, float *gl4f);

  
};

class ColourMapGradient: public ContinuuousColourMap {
public:
  
  void setCol1(std::vector<double> const &array);
  std::vector<double> getCol1();


  void setCol2(std::vector<double> const &array);
  std::vector<double> getCol2();

  ColourMapGradient(Visualiser *pV, float *ref1, float *ref2, float mini=0, float maxi=0);
  ColourMapGradient(Visualiser *pV, const std::vector<double> &, const std::vector<double> &);
  virtual bool operator()(float val, float *gl4f);
  
  float ref1R, ref1G, ref1B, ref1A;
  float ref2R, ref2G, ref2B, ref2A;
};



class ColourMapByType : public ColourMap {
public:
  ColourMapByType(Visualiser *pV, ColourMap *dm, ColourMap *gas, ColourMap *stars);
  bool operator()(int index, float *gl4f, const Particle *p=NULL);
  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05, bool renderColLabel=true);
  virtual int buildMenu();

  virtual void setReferenceAlpha(float alpha);
  virtual void autoRange_Consider(unsigned int n);

  virtual bool click(int button, int state, float x, float y);
 


  ColourMap & getGasMap();
  ColourMap & getDMMap();
  ColourMap & getStarMap();

  void setStarMap(ColourMap &cm);
private:
  ColourMap *pDM;
  ColourMap *pGas;
  ColourMap *pStar;

  float dmAl, gasAl, starAl;

  int *pData;
};

}

#endif // COLOURMAP_H_INCLUDED
