//
// This file is part of SimAn
//
// Copyright (c) 2005-6 Andrew Pontzen
// SimAn may not (currently) be used in any form without
// prior permission. Please contact app26 (at) ast (dot) cam...
// with all enquiries
//
// CScripted implementation

#include "siman.hpp"


CScripted::CScripted(string filenamei, bool noInteractivity) {
 
  filename = filenamei;
  pRetVal = NULL;
  load(noInteractivity);
}

CScripted::CScripted() {
  pRetVal=NULL;
}

void CScripted::haveCreatedObj(CSimanObject *pushObject) {
  if(pushObject!=NULL) 
    createdObjs.insert(pushObject);
  
}


void CScripted::haveDeletedObj(CSimanObject *delObject) {
  createdObjs.erase(delObject);
  
}


string CScripted::className() {
  return "CScripted";
}

CSimanObject *CScripted::getReturnValue() {
  return pRetVal;
}

bool CScripted::deleteSafe(CSimanObject *pCandidate) {

  
  bool allowDelete = true;
  
  
  for(set<CSimanObject*>::iterator m = createdObjs.begin(); m!=createdObjs.end() && allowDelete; m++) {
    if((*m)->references(pCandidate)) {
      allowDelete=false;
    }
  }
  
  for(map<string,CSimanObject*>::iterator m = namedStates.begin(); m!=namedStates.end() && allowDelete; m++) {
    if(((*m).second)->references(pCandidate))
      allowDelete=false;
    if((*m).second==pCandidate)
      allowDelete=false;
  }
  
  return allowDelete;
}



int CScripted::gc() {
  // garbage collector
  int num_deletions=0;
  for(set<CSimanObject*>::iterator n = createdObjs.begin(); n!=createdObjs.end();) { // N.B. increment of n is dependent on result of loop
    CSimanObject *pCandidate=*n;
    bool allowDelete = deleteSafe(pCandidate);
    if(allowDelete) {

      haveDeletedObj(pCandidate);
      delete pCandidate;

      // N.B. Need to restart loop firstly because our iterator is now invalid, and secondly
      //      because dependencies may have changed!

      n=createdObjs.begin();
      num_deletions++;
    } else n++;
  }

  return num_deletions;

}
CScripted::~CScripted() {

  for(map<string,CSimanObject*>::iterator n = namedStates.begin(); n!=namedStates.end(); n++) {
     pair<string, CSimanObject*> keypair = *n;
     delete keypair.second;
     haveDeletedObj(keypair.second);
  }

  gc();
  pRetVal = NULL;
}

string CScripted::readQuoted(istream &from) {
  /*
    // NOT WORKING

  string out;
  string temp;
  from >> temp;
  bool quoted=false;
  bool escaped=false;
  bool first=true;
  
  string::iterator i_q1;
  string::iterator i_q2;

  out=temp;
  do {
    if(!first) out+=" ";
    for(string::iterator i = temp.begin(); i!=temp.end(); i++) {
      if(*i=='\\') {
	escaped=!escaped;
      } else {
	if(!escaped) {
	  if(*i=='"')
	    quoted=!quoted;
	}
	escaped=false;
      }
    }
    out+=temp;
    first=false;
  } while(quoted && !from.eof());
  
  if(quoted)
    throw(CMismatchedParenthesis('"'));

  return out;
  */
	return "";
}


bool CScripted::load(bool noInteractivity) {
 
	if(!siman::fileExists(filename))
    {
      const char* datadir_cstr = getenv("SIMAN_DATA");
      if(datadir_cstr!=NULL)
	filename = ((string)datadir_cstr)+"/"+filename;
      else
	cerr << "CScripted: warning - couldn't find " << filename << ", and environment variable SIMAN_DATA is not set" << endl;
    }
  

  ifstream file(filename.c_str());

  while(file.good() && pRetVal==NULL) {
    string line,command;
    getline(file,line);
    
    if(!noInteractivity && line!="" && line[0]!='#')
	cerr << "Siman> " << line << endl;

    doCommand(line);
  }

  if(pRetVal==NULL && !noInteractivity)
    mainLoop();

  return true;
}

void CScripted::setNamedVar(string name, CSimanObject *var) {
  
  vector<string> individuals;
  splitString(name,'.',individuals);
  if(individuals.size()==1) {
    namedStates[name] = var;
  } else {
    vector<string>::iterator i = individuals.begin();
    
    CSimanObject *obj=NULL;
    
    while(i!=individuals.end()-1) {
      
      string curname = *i;
      
      map<string,CSimanObject*>::iterator mapi = namedStates.find(curname);
      
      if(mapi==namedStates.end()) {
	if(obj!=NULL)
	  obj=obj->getMember(curname); 	// will throw(CUnknownVariable) if the variable does not exist there either
	else
	  throw(CUnknownVariable(curname));

      } else {
	obj = (*mapi).second;
      }
      
      i++;
      
    }
    
    obj->setMember(*i,var);
    
  }
		     
}


void CScripted::eraseNamedVar(string name) {
  
  vector<string> individuals;
  splitString(name,'.',individuals);
  if(individuals.size()==1) {
    namedStates.erase(name);
  } else {
    
    throw(CReadOnly(name));
    
  }
		     
}


CSimanObject* CScripted::getNamedVar(string name, unsigned int require) {
      
  vector<string> individuals;
  splitString(name,'.',individuals);
  
  vector<string>::iterator i = individuals.begin();
  
  CSimanObject *obj=NULL;

  while(i!=individuals.end()) {
  
    string curname = *i;
    map<string,CSimanObject*>::iterator mapi = namedStates.find(curname);
    
    if(mapi==namedStates.end()) {

      if(obj!=NULL)
	obj=obj->getMember(curname); 	// will throw(CUnknownVariable) if the variable does not exist there either
      else
	throw(CUnknownVariable(curname));

    } else {

      obj = (*mapi).second;
    }

    i++;
    
  }
  
  if(!obj->supports(require))
    throw(CTypeError(name));

  return obj;
}

void CScripted::infoVars() {
   
  map<string,CSimanObject*>::iterator i = namedStates.begin();
  cerr << namedStates.size() << " named objects; " << createdObjs.size() << " known distinct objects total" << endl;
  cerr << endl;
  for(i=namedStates.begin();i!=namedStates.end();i++) {
    cerr << (*i).first << "\t" << (*i).second->className() << endl;
  }

}


void CScripted::splitString(const string& input, 
			   const char &delimiter, vector<string>& results)
{
  istringstream istream(input);
  string res;
  while(getline(istream,res,delimiter)) {
    results.push_back(res);
  }
}

void CScripted::mainLoop() {
  cerr << "Siman> ";
  while(pRetVal==NULL) {
    string command,line;
    getline(cin,line);
    if(line!="" && line[0 ]!='#') {
      
      CSimanObject *discard = doCommand(line);
      if(discard!=NULL) {
	cerr << "Warning: returned object " << discard->className() << " discarded" << endl;
	delete discard;
      }
    }
    gc();
    cerr << "Siman> ";
  }
}



void CScripted::pollStdIn() {

  static bool firstCall = true;
  if(firstCall) {
    cerr << "Siman> ";
    firstCall = false;
  }

  // NOTE:
  //
  // Need to read only when we know the call won't block
  // 
  // the best solution here would be to use
  // if(cin.rdbuf()->in_avail()>0), however this does not
  // work as you would like to expect, at least for the
  // GNU libraries

#ifdef _WIN32
  bool retval=(cin.rdbuf()->in_avail()>0);
#else
  fd_set rfds;
  struct timeval tv;
  int retval;
  FD_ZERO(&rfds);
  FD_SET(0, &rfds);
  tv.tv_sec = 0;
  tv.tv_usec = 0;
  retval = select(1, &rfds, NULL, NULL, &tv);
#endif

  if (retval) {
    // input waiting
    string command,line;
    getline(cin,line);
    if(line!="" && line[0 ]!='#') {
      
      CSimanObject *discard = doCommand(line);
      if(discard!=NULL) {
	cerr << "Warning: returned object " << discard->className() << " discarded" << endl;
	delete discard;
      }
    }
    gc();
    cerr << "Siman> ";
    
  }
}

CSimanObject * CScripted::doCommand(string line) {
  
  CSimanObject *createdObj = NULL;
  CSimanObject *callObj = NULL;

  vector<string> individual_commands;
  splitString(line,';',individual_commands);
  if(individual_commands.size()>1) {
    for(vector<string>::iterator i = individual_commands.begin(); i!=individual_commands.end(); i++) {
      haveCreatedObj(doCommand(*i));
    }
    return NULL;
  }

  try {

    
    vector<string> com_eq;
    splitString(line,'=',com_eq);



    if(com_eq.size()>2) {
      // some versions of g++ seem to object to:
      // 
      // throw(CSyntaxError());
      // 
      // - no idea why

      CSyntaxError n;
      throw(n);
    }

    if(com_eq.size()==2) {
      vector<string>::iterator assign = com_eq.begin();
      vector<string>::iterator read = assign+1;
      CSimanObject *object = NULL;
      try {
	try {
	  object = new CValue(*read,true,true);
	  haveCreatedObj(object);
	} catch(CTypeError &te) {
      	  object = getNamedVar(*read);
	}
      } catch(CSimanException &ukc) {

	object = doCommand(*read);
      }
      if(object!=NULL) setNamedVar(*assign,object);
      return NULL;
    }

    // not an assignment, continue
    
    stringstream stream(line,stringstream::in);
    string command;
    stream >> command;
    transform(command.begin(),command.end(),command.begin(),(int(*)(int))tolower);

    // split up command into "." sections

    vector<string> com_sep;
    splitString(command,'.',com_sep);
    
    if(com_sep.size()>1) {
      vector<string>::iterator i;
      for(i=com_sep.begin();i!=com_sep.end()-1;i++) {
	if(callObj==NULL) callObj = (getNamedVar(*i));
	else callObj = callObj->getMember(*i);
      }
      command = (*i);
    }
  
    
    if(command[0] == '#') {
      // comment
      return NULL;  
    } 

    if(command=="return") { 
      string varname;
      stream >> varname;
      pRetVal = (CSimSnap*) getNamedVar(varname,CSimanObject::SimSnap);
    } else if(command=="new") {
      string type;
      stream >> type;
      if(type=="")
	throw(CSyntaxError("new [type] (...)"));

#ifdef SIMAN_VIS_PLUS

      if(type=="visualise") {
	string sim;
	stream >> sim;
	if(sim=="") 
	  throw(CSyntaxError("new visualise [simulation]"));
	CSimSnap* pSim = (CSimSnap*) getNamedVar(sim,CSimanObject::SimSnap);
	CVisualise *pVis;
	pVis=new CVisualise(pSim,0,this);
	return pVis;
      }
      if(type=="colourmap") {
	string viso;
	stream >> viso;
	if(viso=="")
	  throw(CSyntaxError("new colourmap [visualisation object] [colourmap type] (...)"));

	CVisualise *pV = (CVisualise*) getNamedVar(viso,CSimanObject::Visualise);
	string type;
	stream >> type;
	if(type=="")
	  throw(CSyntaxError("new colourmap [visualisation object] [colourmap type] (...)"));
	if(type[0]=='g') {
	  float ref1[4]={1,0,0,1};
	  float ref2[4]={0,0,1,1};
	  CColourMapGradient *pCM = new CColourMapGradient(pV,ref1,ref2);
	  return pCM;
	}
	if(type[0]=='s') {
	  CColourMap *pCM = new CColourMap(pV);
	  pCM->setReferenceColour(0,0,1,0.3);
	  return pCM;
	}
	if(type[0]=='t') {
	  string dm,gas,stars;
	  stream >> dm >> gas >> stars;
	  CColourMap *pDmCM = (CColourMap*) getNamedVar(dm,CSimanObject::ColourMap);
	  CColourMap *pGasCM = (CColourMap*) getNamedVar(gas,CSimanObject::ColourMap);
	  CColourMap *pStarCM = (CColourMap*) getNamedVar(stars,CSimanObject::ColourMap);
	  if(stars=="")
	    throw(CSyntaxError("new colourmap [visualisation object] type [dm colourmap] [gas colourmap] [star colourmap]"));
	  CColourMapByType *pCM = new CColourMapByType(pV,pDmCM,pGasCM,pStarCM);
	  return pCM;
	}
      }
#endif // SIMAN_VIS
    } else 
      if(command == "load") {
      string filename;
      stream >> filename;
      createdObj = CSimSnap::loadFile(filename);
    } else if(command == "vars") {
      infoVars();
    } else if(command == "union") {
      string var1,var2;
      stream >> var1 >> var2;
      CSimSnap* pV1 = (CSimSnap*) getNamedVar(var1,CSimanObject::SimSnap);
      CSimSnap* pV2 = (CSimSnap*) getNamedVar(var2,CSimanObject::SimSnap);
      CUnion *pU = new CUnion(pV1);
      pU->add(pV1);
      pU->add(pV2);
      while(!stream.eof()) {
	string varn;
	stream >> varn;
	CSimSnap *pVn = (CSimSnap*) getNamedVar(varn,CSimanObject::SimSnap);
	pU->add(pVn);
      }
      createdObj = pU;
      } else if(command=="delete" || command=="del") {
	string obj;
	stream >> obj;
	CSimanObject *pObj = getNamedVar(obj);
	eraseNamedVar(obj);
      } else if(command!="") {
      if(callObj!=NULL) 
	createdObj = callObj->dispatch(command,&stream,this);
      else
	throw(CUnknownCommand(command));
    } 
  } catch(CSimanException &uc) {
    cerr << uc.getText() << endl;
  }


  if(createdObj!=NULL) haveCreatedObj(createdObj);

  cerr.flush();
  cout.flush();

  return createdObj;
}
