//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
#ifndef __COLOURMAP_H_INCLUDED

#define __COLOURMAP_H_INCLUDED

class CVisualise;
	
class CColourMap : public CSimanObject {
public:
  friend class CColourMapByType;

  CColourMap(CVisualise *pV);
  virtual bool operator()(int index, float *gl4f, CParticle *p=NULL);

  virtual void autoRange();  // pick min and max values

  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05);

  virtual void setReferenceColour(float r, float g, float b, float a);
  virtual void setReferenceAlpha(float a);
  float getReferenceAlpha();

  virtual bool click(int button, int state, float x, float y); 
       // returns true if changes have been made that require a rewrite of list

  virtual int buildMenu();

  virtual void menuCallback(int id);

  virtual unsigned int supports();
  virtual string className();
  virtual CSimanObject * dispatch(string command, std::istream *stream, CScripted *script);


  // reference values, will be used differently by different
  // child classes.
  

  string label;

  // simsnap
  CVisualise *pVis;
  CColourMap *pParent;

  bool display;

protected:
  virtual void autoRange_Consider(int n);

  float refR, refG, refB, refA;

  bool autoranging;           
};

// 
// CContinuuousColourMap
//

class CContinuuousColourMap : public CColourMap {
public:
  CContinuuousColourMap(CVisualise *pV, float mini=0, float maxi=0);
  
  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05);
  virtual bool operator()(int index, float *gl4f, CParticle *p=NULL);
  virtual bool operator()(float val, float *gl4f);
  
  virtual bool click(int button, int state, float x, float y); 
  virtual int buildMenu();

  virtual float particleToColourIndex(int n, CParticle *p=NULL);

  virtual void autoRange();
  virtual void autoRange_Consider(int n);

  virtual void menuCallback(int id);
     
  void setColourBy(int colourBy);

  string getColourLabel(int n);

  virtual CSimanObject * dispatch(string command, std::istream *stream, CScripted *script);

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
  float *map;

  static const int useNe = 1;
  static const int useNHp = 2;
  static const int useNH0 = 3;
  static const int useTemp = 4;
  static const int useRho = 0;
  static const int useU = 5;
  static const int useMetal = 6;
  static const int use_numOptions = 6;
  ///@}

};

class CColourMapSpectrum : public CContinuuousColourMap {
public:

  CColourMapSpectrum(CVisualise *pV, float mini=0, float maxi=0);
  virtual bool operator()(float val, float *gl4f);

  
};

class CColourMapGradient: public CContinuuousColourMap {
public:

  CColourMapGradient(CVisualise *pV, float *ref1, float *ref2, float mini=0, float maxi=0);
  virtual bool operator()(float val, float *gl4f);
  virtual CSimanObject * dispatch(string command, std::istream *stream, CScripted *script);

  float ref1R, ref1G, ref1B, ref1A;
  float ref2R, ref2G, ref2B, ref2A;
};



class CColourMapByType : public CColourMap {
public:
  CColourMapByType(CVisualise *pV, CColourMap *dm, CColourMap *gas, CColourMap *stars);
  bool operator()(int index, float *gl4f, CParticle *p=NULL);
  virtual void renderColourBar(int ext, int exty, float x0=0, float x1=1., float y0=0.,float y1=0.05);
  virtual int buildMenu();

  virtual void setReferenceAlpha(float alpha);
  virtual void autoRange_Consider(int n);

  virtual bool click(int button, int state, float x, float y);
 
  virtual CSimanObject* getMember(const string &var);

  CColourMap *pDM;
  CColourMap *pGas;
  CColourMap *pStar;

  float dmAl, gasAl, starAl;

  int *pData;
};

#endif // SPHKERNEL_H_INCLUDED
