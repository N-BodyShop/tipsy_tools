// scripted.cpp - part of SimAn Simulation Analysis Library
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





#include "siman.hpp"
#include "../python/python_definitions.hpp"
#include <Python.h> 
#include <node.h> 
#include <errcode.h> 
#include <grammar.h> 
#include <parsetok.h> 
#include <compile.h>
#include <eval.h>

#ifdef SIMAN_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

namespace siman {

  using namespace python;


  int interpreter_use_count(0);
 

  boost::python::object Scripted::main_namespace;
  boost::python::object Scripted::main_module;

  Scripted::Scripted(string filenamei, bool noInteractivity) {
    pRetVal=NULL;
    filename = filenamei;
    initPython();
    load(noInteractivity);
  }

  Scripted::Scripted() {
    pRetVal=NULL;
    initPython();
  }

  Scripted::~Scripted() {
    deInitPython();
  }


  SimanObject *Scripted::getReturnValue() {
    return pRetVal; // TODO
  }

  void Scripted::deInitPython() {
    interpreter_use_count--;
    if(interpreter_use_count==0) Py_Finalize();
  }

  void Scripted::doCommand(string st) {
    try {
      
      object result = object(handle<>(PyRun_String(
						   const_cast<char*>(st.c_str())
						   , Py_single_input
						   , main_namespace.ptr()
						   , main_namespace.ptr())
				      ));
      
      
    } catch(error_already_set) {
      PyErr_Print();
    }
 
  }

  void Scripted::initPython() {

    if(interpreter_use_count==0 && Py_IsInitialized())
      interpreter_use_count = 1;

    if(!Py_IsInitialized()) {
      char *path = getenv("PYTHONPATH");
      char *sdata = getenv("SIMAN_DATA");
      if(sdata!=NULL && path!=NULL) {
	string path_plus = (string)path+":"+(string)sdata;
     
	setenv("PYTHONPATH",path_plus.c_str(),1);
      } else if(sdata!=NULL) {
	setenv("PYTHONPATH",sdata,1);
      }

      if(interpreter_use_count==0) {
	// artificially only allow interpreter_use_count to reach
	// zero on extra deinit call in ~Config.
	interpreter_use_count++;
	if(PyImport_AppendInittab("siman",initsiman) == -1)
	  throw SimanException("Failed to import siman classes to python");
	Py_SetProgramName("siman");
	Py_Initialize();
	char *arg_none = "\0";
	PySys_SetArgv(1,&arg_none);
    
    
	PyObject *v;
    
	v = PySys_GetObject("ps1");
	if (v==NULL) {
	  PySys_SetObject("ps1", v = PyString_FromString(""));
	  Py_XDECREF(v);
	}
	v = PySys_GetObject("ps2");
	if (v == NULL) {
	  PySys_SetObject("ps2", v = PyString_FromString("[SimAn]... "));
	  Py_XDECREF(v);
	}
      }
    }

    try {
     
      main_module = object(handle<>(borrowed(PyImport_AddModule("__main__"))));
      main_namespace = main_module.attr("__dict__");

      static bool inserted=false;
      if(!inserted)
	object result = object(handle<>(PyRun_String(
						     "import siman"
						     , Py_single_input
						     , main_namespace.ptr()
						     , main_namespace.ptr())
					));
      inserted = true;
      //      cout << extract<float>(result) << endl;
    } catch(error_already_set) {
      PyErr_Print();
    }
   
    interpreter_use_count++;
  }

  void Scripted::injectAndManage(string varname, SimanObject * obj) {
    python::pWaiting = obj;
    string command = varname+"=siman._siman_internal_fetch_manage()";
    try {
    
      object result = object(handle<>(PyRun_String(
						   const_cast<char*>(command.c_str())
						   , Py_single_input
						   , main_namespace.ptr()
						   , main_namespace.ptr())
				      ));
      //      cout << extract<float>(result) << endl;
    } catch(error_already_set) {
      cerr << "Scripted: Python error whilst injecting variable..." << endl << "*********************" ;
      PyErr_Print();
      cerr << "************************" << endl;
    
    }
  

  }

  void Scripted::injectAsReference(string varname, SimanObject & obj) {
    python::pWaiting = &obj;
    string command = varname+"=siman._siman_internal_fetch_reference()";
    try {
    
      object result = object(handle<>(PyRun_String(
						   const_cast<char*>(command.c_str())
						   , Py_single_input
						   , main_namespace.ptr()
						   , main_namespace.ptr())
				      ));
    
    } catch(error_already_set) {
      cerr << "Scripted: Python error whilst injecting variable..." << endl << "*********************" ;
      PyErr_Print();
      cerr << "************************" << endl;
    
    }
  

  }


  bool Scripted::load(bool noInteractivity) {
 
    if(!siman::fileExists(filename))
      {
	const char* datadir_cstr = getenv("SIMAN_DATA");
	if(datadir_cstr!=NULL)
	  filename = ((string)datadir_cstr)+"/"+filename;
	else {
	  cerr << "Scripted: warning - couldn't find " << filename << ", and environment variable SIMAN_DATA is not set" << endl;
       
	}
      }
  
    // pre-process script
  
    ifstream preprocess_file(filename.c_str());
  
    string initobj;

    if(preprocess_file.is_open()) {
      while(!preprocess_file.eof()) {
	string temp;
	preprocess_file >> temp;
	if(temp=="SIMAN_INIT_OBJECT")
	  preprocess_file >> initobj;
      }
    } else {
      throw(FileError(filename));
    }

    if(initobj=="")
      throw(SyntaxError("Can't use file as initializer; no SIMAN_INIT_OBJECT found"));

    preprocess_file.close();

    // Re-open for PYTHON use


    FILE *fhand = fopen(filename.c_str(),"r");

    if(fhand==0)
      throw FileError(filename);
  
    try {
      PyRun_AnyFile(fhand,const_cast<char*>(filename.c_str()));
    
      fclose(fhand);
    
      PyObject* obj = PyRun_String(
				   const_cast<char*>((initobj+(string)"()").c_str())
				   , Py_eval_input
				   , main_namespace.ptr()
				   , main_namespace.ptr());
      if(obj!=NULL) {
	extract<SimanObject*> extractor(obj);
	if(extractor.check()) {
	  pRetVal = extractor();
	}
	else
	  throw(SyntaxError("Can't use file as initializer; Python script returned incompatible object"));
      } else {
	PyErr_Print();
	throw(SyntaxError("Can't use file as initializer; error in Python script"));
      }
    } catch(error_already_set) {
      PyErr_Print();
      throw(SyntaxError("Can't use file as initializer; error in Python script"));
    }
    return true;
  }




  void Scripted::mainLoop() {
 
    
#ifdef SIMAN_READLINE
    
    rl_readline_name="siman";
    rl_initialize();
    
 

    int i, j, done = 0;                          /* lengths of line, code */
    char ps1[] = "[SimAn]>>> ";
    char ps2[] = "[SimAn]... ";
    char *prompt = ps1;
    char *msg, *line, *code = NULL;
    PyObject *src;
    PyObject *exc, *val, *trb, *obj, *dum;
    
       
    
    while (!done)
      {
	line = readline (prompt);
	
	if(NULL==line)
	  {
	    done = 1;
	  }
	else
	  {
	    	
	    i = strlen (line);
	    
	    if(i>0 || code!=NULL)
	    {
	      	
	      
	      if (i > 0)
		add_history (line);                    /* save non-empty lines */
	      	
	      if (NULL == code)                        /* nothing in code yet */
		j = 0;
	      else
		j = strlen (code);
	    
	     
	      code = (char*) realloc (code, i + j + 2);
	      if (NULL == code)                        /* out of memory */
		exit (1);
	      
	      if (0 == j)                              /* code was empty, so */
		code[0] = '\0';                        /* keep strncat happy */
	      
	      strncat (code, line, i);                 /* append line to code */
	      code[i + j] = '\n';                      /* append '\n' to code */
	      code[i + j + 1] = '\0';
	   
	      src = Py_CompileString (code, "<stdin>", Py_single_input);       
	      
	      if (NULL != src)                         /* compiled just fine - */
		{
		  if (ps1  == prompt ||                  /* ">>> " or */
		      '\n' == code[i + j - 1])           /* "... " and double '\n' */
		    {                                               /* so execute it */
		    
		      dum = PyEval_EvalCode ((PyCodeObject *)src, main_namespace.ptr(), main_namespace.ptr());
		      Py_XDECREF (dum);
		      Py_XDECREF (src);
		      free (code);
		      code = NULL;
		      if (PyErr_Occurred ())
			PyErr_Print ();
		      prompt = ps1;
		    }
		}                                        /* syntax error or E_EOF? */
	      else if (PyErr_ExceptionMatches (PyExc_SyntaxError))           
		{
		 
		  PyErr_Fetch (&exc, &val, &trb);        /* clears exception! */
		  
		  if (PyArg_ParseTuple (val, "sO", &msg, &obj) &&
		      !strcmp (msg, "unexpected EOF while parsing")) /* E_EOF */
		    {
		      Py_XDECREF (exc);
		      Py_XDECREF (val);
		      Py_XDECREF (trb);
		      prompt = ps2;
		    }
		  else                                   /* some other syntax error */
		    {
		      PyErr_Restore (exc, val, trb);
		      PyErr_Print ();
		      free (code);
		      code = NULL;
		      prompt = ps1;
		    }
		}
	      else                                     /* some non-syntax error */
		{
		  PyErr_Print ();
		  free (code);
		  code = NULL;
		  prompt = ps1;
		}

	    } 
	  free (line);
	  }
      }
    

#else 
    // no readline support
    
    while(!cin.eof()) {
      cout << "[SimAn]>>> ";
      char line[1024];
      cin.getline(line,1024);
    
      if(!cin.eof()) {
    static bool firstCall = true;
	try {
	
	  object result = object(handle<>( PyRun_String(
							line
							, Py_single_input
							, main_namespace.ptr()
							, main_namespace.ptr())
					   ));

	} catch(error_already_set) {
	  PyErr_Print();
	}
      }
    
    }

#endif

  }



  void Scripted::pollStdIn() {

    static bool firstCall = true;
    if(firstCall) {
      cout << "[SimAn]>>> ";
      cout.flush();
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
      try {
      
	/*
	  object result = object(handle<>(PyRun_String(
	  line
	  , Py_single_input
	  , main_namespace.ptr()
	  , main_namespace.ptr())
	  ));
	*/
	PyRun_InteractiveOne(stdin,"stdin");
	//      cout << extract<float>(result) << endl;
      } catch(error_already_set) {
	PyErr_Print();
      }
      cout << "[SimAn]>>> ";
      cout.flush();
      cerr.flush();
   
    }
  }

} // namespace siman
