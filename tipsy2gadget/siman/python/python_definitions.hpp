#include <boost/python.hpp>
#include "scitbx/container_conversions.h"

#ifdef SIMAN_EXTRA
#include "../extra.hpp"
#endif

namespace siman {
namespace python {
  using namespace boost::python;

  // AUTO - CONVERTER FOR PYTHON STRING -> UNIT
  struct unit_from_python_str
  {
    unit_from_python_str() {
      boost::python::converter::registry::push_back(
						    &convertible,
						    &construct,
						    boost::python::type_id<Unit>());
	}
    
    static void* convertible(PyObject * obj_ptr) {
      if(!PyString_Check(obj_ptr)) return 0;
      try {
	Unit un(PyString_AsString(obj_ptr));
      } catch (UnknownUnit &u) {
	return 0;
      }
      return obj_ptr;
    }
    
    static void construct(PyObject *obj_ptr, 
			  boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      const char* value = PyString_AsString(obj_ptr);
      if(value==NULL) boost::python::throw_error_already_set();
      void *storage = (
		       (boost::python::converter::rvalue_from_python_storage<Unit>*)data)->storage.bytes;
      try {
	new (storage) Unit(value);
      } catch (UnknownUnit &u) {
	// FIXME - 
	boost::python::throw_error_already_set();
      }
      data->convertible = storage;
    }
  };
    
    
    void wrap_SimSnap_convertUnits_0(SimSnap* a) {
      a->convertUnits();
    }
  
  
  void wrap_SimSnap_convertUnits_1(SimSnap* a, const Unit &u1) {
    a->convertUnits(u1);
  }
  
  
  void wrap_SimSnap_convertUnits_2(SimSnap* a, const Unit &u1, const Unit &u2) {
    a->convertUnits(u1,u2);
  }
  
  
  void wrap_SimSnap_convertUnits_3(SimSnap* a, const Unit &u1, const Unit &u2, const Unit &u3) {
    a->convertUnits(u1,u2,u3);
  }


  void wrap_SimSnap_convertUnits_4(SimSnap* a, const Unit &u1, const Unit &u2, const Unit &u3, const Unit &u4) {
    a->convertUnits(u1,u2,u3,u4);
  }


  double wrap_Unit_convertTo_1(const Unit *u1, const Unit *u2) {
    return u1->convertTo(*u2);
  }


  double wrap_Unit_convertTo_2(const Unit *u1, const Unit *u2, const SimSnap *s) {
    return u1->convertTo(*u2,s);
  }

  void wrap_SimanArray_set(SimanArray *s, int n, float v) {
    (*s)[n]=v;
  }

  float wrap_SimanArray_get(const SimanArray *s, int n) {
    return (*s)[n];
  }

  ColumnGrid & wrap_Grid_get(Grid *g, int n) {
    return (*g)[n];
  }


  ColumnList & wrap_ColumnGrid_get(ColumnGrid *g, int n) {
    return (*g)[n];
  }

  SimSnap & wrap_ColumnList_get(ColumnList *s, int n) {
    return (*s)[n];
  }

  SimanArray & wrap_SimSnap_getArray(SimSnap *s, string t) {
    return s->getArray(t);
  }

  const SimanArray & wrap_SimSnap_getConstArray(const SimSnap *s, string t) {
    return s->getConstArray(t);
  }

  
  vector<double> wrap_SimSnap_shrinkSphereCentre_0(const SimSnap *s) {
    return s->shrinkSphereCentre();
  }
  

  vector<double> wrap_SimSnap_shrinkSphereCentre_1(const SimSnap *s,float sf) {
    return s->shrinkSphereCentre(sf);
  }
  

  vector<double> wrap_SimSnap_shrinkSphereCentre_2(const SimSnap *s, float sf, int mp) {
    return s->shrinkSphereCentre(sf,mp);
  }


  vector<double> wrap_SimSnap_shrinkSphereCentre_3(const SimSnap *s, float sf, int mp, float va) {
    return s->shrinkSphereCentre(sf,mp,va);
  }

 
  vector<double> wrap_SimSnap_shrinkSphereCentreVel_0(const SimSnap *s) {
    return s->shrinkSphereCentreVel();
  }
  

  vector<double> wrap_SimSnap_shrinkSphereCentreVel_1(const SimSnap *s,float sf) {
    return s->shrinkSphereCentreVel(sf);
  }
  

  vector<double> wrap_SimSnap_shrinkSphereCentreVel_2(const SimSnap *s, float sf, int mp) {
    return s->shrinkSphereCentreVel(sf,mp);
  }


  vector<double> wrap_SimSnap_shrinkSphereCentreVel_3(const SimSnap *s, float sf, int mp, float va) {
    return s->shrinkSphereCentreVel(sf,mp,va);
  }

  void wrap_SimSnap_initialiseSPH_0(const SimSnap *s) {
    s->initialiseSPH();
  }

  void wrap_SimSnap_initialiseSPH_1(const SimSnap *s, int nSPH) {
    s->initialiseSPH(nSPH);
  }

  float wrap_SimSnap_getVirialRadius_0(const SimSnap *s) {
    return s->getVirialRadius();
  }

  vector<double> wrap_SimSnap_centreOfMass(const SimSnap *s) {
    return s->centreOfMass();
  }

  vector<double> wrap_SimSnap_centreOfMassVel(const SimSnap *s) {
    return s->centreOfMassVel();
  }

  float wrap_SimSnap_getSmooth_1(const SimSnap *s , unsigned int id) {
    return s->getSmooth(id);
  }
 
  void wrap_Grid_ctor_4(Grid*& g, SimSnap *s, int nx, int ny, int nz) {
    cerr << g << endl;
    g = new Grid(*s,nx,ny,nz);
  }

  Grid * wrap_Grid_ctor_6(SimSnap *s, int nx, int ny, int nz, int maxRef, int minPart) {
    return new Grid(*s,nx,ny,nz,maxRef,minPart);
  }

  std::list<std::pair<int, double> > wrap_SimSnap_getNearestNeighbours_2(SimSnap *s, const Particle & a, int i) {
    return s->getNearestNeighbours(a,i);
  }

   
  // very hacky way of importing C++ managed data into python interpreter environment:
  SimanObject * pWaiting;
  
  SimanObject * fetch() {
    return pWaiting;
  }


  BOOST_PYTHON_MODULE(siman)
  {
    
    class_<SimanObject>("SimanObject",no_init);

    /*
    def("pm", &printMessage);
    def("pm3", &printMessageWith3);

    class_<test1>("test1",init<>());
    // class_<test2, bases<test1> >("test2",init<>());
    class_<test3 >("test3",init<test1 &>());
    */

    
    class_<Rational>("Rational",init<int>())
      .def(init<int,int>())
      .def(init<string>())
      .def("__str__",&Rational::toString)
      .def("__flt__",&Rational::dbl)
      ;
    

    class_<SimanArray >("SimanArray",init<const SimanArray &>())
      .def("__getitem__",&wrap_SimanArray_get)
      .def("__setitem__",&wrap_SimanArray_set)
      .def("__len__",&SimanArray::size)
      .def("copyIn",&SimanArray::operator=)
      .def("isSet",&SimanArray::isSet)
      .def("unset",&SimanArray::unset)
      .def("size",&SimanArray::size)
      .def("getUnits",&SimanArray::getUnits)
      .def("convertUnits",&SimanArray::convertUnits)
      .def("max",&SimanArray::max)
      .def("min",&SimanArray::min)
      .def("mean",&SimanArray::mean)
      .def("total",&SimanArray::total)
      .def("readTipsy",&SimanArray::readTipsy)
      .def("__pow__",&SimanArray::power)
      .def(self + self)
      .def(self - self)
      .def(self / self)
      .def(self * self)
      .def(self += float())
      .def(self -= float())
      .def(self *= float())
      .def(self /= float())
      .def(self + float())
      .def(float() + self)
      .def(self - float())
      .def(float() - self)
      .def(self * float())
      .def(float() * self)
      .def(self / float())
      .def(float() / self)
      .def(self += self)
      .def(self -= self)
      .def(self *= self)
      .def(self /= self)
      .def(Unit() * self)
      .def(self * Unit())
      .def(Unit() / self)
      .def(self / Unit())

      ;
    
    class_<SimanArrayVirtual, bases<SimanArray> >("SimanArrayVirtual",no_init);

    class_<SimanArraySubscripted, bases<SimanArray> >("SimanArraySubscripted",no_init);

    void (SimSnap::*SimSnap_wrap_centreOn)(vector<double>)  = &SimSnap::centreOn;
    void (SimSnap::*SimSnap_wrap_centreOn_3)(float,float,float) = &SimSnap::centreOn;

    void (SimSnap::*SimSnap_wrap_centreOnVel)(vector<double>) = &SimSnap::centreOnVel;
    void (SimSnap::*SimSnap_wrap_centreOnVel_3)(float,float,float) = &SimSnap::centreOnVel;

    float (SimSnap::*SimSnap_wrap_scaleToAngDiDist_1)(float) const = &SimSnap::scaleToAngDiDist;
    float (SimSnap::*SimSnap_wrap_scaleToAngDiDist_2)(float,float) const = &SimSnap::scaleToAngDiDist;

    float (SimSnap::*SimSnap_wrap_scaleToComovingDist_1)(float) const = &SimSnap::scaleToComovingDist;
    float (SimSnap::*SimSnap_wrap_scaleToComovingDist_2)(float,float) const = &SimSnap::scaleToComovingDist;



	
    int (SimSnap::*SimSnap_wrap_deReference)(int i, SimSnap *pRel) const = &SimSnap::deReference;

    object simsnap_class = 
    class_<SimSnap, boost::noncopyable, bases<SimanObject> , auto_ptr<SimSnap> >("SimSnap",no_init)
      .def("__getitem__",&SimSnap::getParticleRef,
	   with_custodian_and_ward_postcall<0,1,
	   return_value_policy<copy_non_const_reference> >())
      .def("getParticle",&SimSnap::getParticleRef,
	   with_custodian_and_ward_postcall<0,1,
	   return_value_policy<copy_non_const_reference> >())
      .def("getConstParticle",&SimSnap::getConstParticleRef,
	   with_custodian_and_ward_postcall<0,1,
	   return_value_policy<copy_const_reference> >())
      .def("getNumParticles",&SimSnap::getNumParticles)
      .def("getNearestNeighbours",&SimSnap::getNearestNeighbours)
      .def("getNearestNeighbours",&wrap_SimSnap_getNearestNeighbours_2)
      .def("__len__",&SimSnap::getNumParticles)
      .def("deReference",SimSnap_wrap_deReference)
      .def("convertUnits",&wrap_SimSnap_convertUnits_0)
      .def("convertUnits",&wrap_SimSnap_convertUnits_1)
      .def("convertUnits",&wrap_SimSnap_convertUnits_2)
      .def("convertUnits",&wrap_SimSnap_convertUnits_3)
      .def("convertUnits",&wrap_SimSnap_convertUnits_4)
      .def("convertUnits",&SimSnap::convertUnits)
      .def("getTotalMass",&SimSnap::getTotalMass)
      .def("transform",&SimSnap::transform)
      .def("copy",&SimSnap::copy, return_value_policy<manage_new_object>())
      .def("copyTransform", &SimSnap::copyTransform, return_value_policy<manage_new_object>())
      .def("getArray", &wrap_SimSnap_getArray, return_internal_reference<1> () )
      .def("getConstArray",&wrap_SimSnap_getConstArray, return_internal_reference<1>() )
      .def("getNumArrays",&SimSnap::getNumArrays)
      .def("getDistanceUnits",&SimSnap::getDistanceUnits)
      .def("getMassUnits",&SimSnap::getMassUnits)
      .def("getVelocityUnits",&SimSnap::getVelocityUnits)
      .def("getEnergyUnits",&SimSnap::getEnergyUnits)
      .def("getDensityUnits",&SimSnap::getDensityUnits)
      .def("createArray", &SimSnap::createArray, return_internal_reference<1>() )
      .def("subset",&SimSnap::subset, with_custodian_and_ward_postcall<0,1,
	   return_value_policy<manage_new_object> >())
      .def("copySubset",&SimSnap::copySubset, return_value_policy<manage_new_object>())
      .def("diagnostics",&SimSnap::diagnostics)
      .def("shrinkSphereCentre",&wrap_SimSnap_shrinkSphereCentre_0)
      .def("shrinkSphereCentre",&wrap_SimSnap_shrinkSphereCentre_1)
      .def("shrinkSphereCentre",&wrap_SimSnap_shrinkSphereCentre_2)
      .def("shrinkSphereCentre",&wrap_SimSnap_shrinkSphereCentre_3)
      .def("shrinkSphereCentreVel",&wrap_SimSnap_shrinkSphereCentreVel_0)
      .def("shrinkSphereCentreVel",&wrap_SimSnap_shrinkSphereCentreVel_1)
      .def("shrinkSphereCentreVel",&wrap_SimSnap_shrinkSphereCentreVel_2)
      .def("shrinkSphereCentreVel",&wrap_SimSnap_shrinkSphereCentreVel_3)
      .def("angMom",&SimSnap::angMom)
      .def("rotCurve",&SimSnap::rotCurve)
      .def("getSPHColDen",&SimSnap::getSPHColDen)
      .def("centreOn",SimSnap_wrap_centreOn)
      .def("centreOn",SimSnap_wrap_centreOn_3)
      .def("centreOfMass",&wrap_SimSnap_centreOfMass)
      .def("centreOnVel",SimSnap_wrap_centreOnVel)
      .def("centreOnVel",SimSnap_wrap_centreOnVel_3)
      .def("centreOfMassVel",&wrap_SimSnap_centreOfMassVel)
      .def("getVirialRadius",&SimSnap::getVirialRadius)
      .def("getVirialRadius",&wrap_SimSnap_getVirialRadius_0)
      .def("getSmooth",&wrap_SimSnap_getSmooth_1)
      .def("makeP",&SimSnap::makeP)
      .def("makeMu",&SimSnap::makeMu)
      .def("TempFromU",&SimSnap::TempFromU)
      .def("UFromTemp",&SimSnap::UFromTemp)
      .def("scaleToHubble",&SimSnap::scaleToHubble)
      .def("scaleToAge",&SimSnap::scaleToAge)
      .def("scaleToHubble1e10yr",&SimSnap::scaleToHubble1e10yr)
      .def("scaleToAge1e10yr",&SimSnap::scaleToAge1e10yr)
      .def("scaleToAngDiDist",SimSnap_wrap_scaleToAngDiDist_1)
      .def("scaleToAngDiDist",SimSnap_wrap_scaleToAngDiDist_2)
      .def("scaleToComovingDist",SimSnap_wrap_scaleToComovingDist_1)
      .def("scaleToComovingDist",SimSnap_wrap_scaleToComovingDist_2)
      .def("addHubbleFlow",&SimSnap::addHubbleFlow)
      .def("subtractHubbleFlow",&SimSnap::subtractHubbleFlow)
      .def("makeNaiveAbsProfile",&SimSnap::makeNaiveAbsProfile)
      .def("initialiseSPH",&SimSnap::initialiseSPH)
      .def("initialiseSPH",&wrap_SimSnap_initialiseSPH_1)
      .def("initialiseSPH",&wrap_SimSnap_initialiseSPH_0)
      .def("makeSPHAbsProfile",&SimSnap::makeSPHAbsProfile)
      .def("write",&SimSnap::write)
      .add_property("criticalDensity",&SimSnap::criticalDensity)
      .add_property("criticalDensity0",&SimSnap::criticalDensity0)
      .add_property("omegaM0",&SimSnap::getOmegaM0,&SimSnap::setOmegaM0)
      .add_property("omegaLambda0",&SimSnap::getOmegaLambda0,&SimSnap::setOmegaLambda0)
      .add_property("h",&SimSnap::getHubble,&SimSnap::setHubble)
      .add_property("hubble",&SimSnap::getHubble,&SimSnap::setHubble)
      .add_property("z",&SimSnap::getRedshift,&SimSnap::setRedshift)
      .add_property("redshift",&SimSnap::getRedshift,&SimSnap::setRedshift)
      ;

    simsnap_class.attr("tipsy") = 2;
    simsnap_class.attr("gadget") = 4;
    simsnap_class.attr("native") = 16;


    class_<BaseSimSnap, boost::noncopyable, bases<SimSnap> >("BaseSimSnap",init<>());
    class_<GadgetSnap, boost::noncopyable, bases<BaseSimSnap, SimSnap> >("GadgetSnap",no_init);
    class_<TipsySnap, boost::noncopyable, bases<BaseSimSnap, SimSnap> >("TipsySnap",no_init);
    class_<Grid, boost::noncopyable, bases<SimSnap> >("Grid",init<SimSnap&,int,int,int,optional<int,int,bool> >()[with_custodian_and_ward<1,2>()])
      .def("__getitem__", &wrap_Grid_get, return_internal_reference<1>())    
      ;
    

    double (ColumnGrid::*ColumnGrid_wrap_getColDen_int)(int, int, string) = &ColumnGrid::getColDen;
    double (ColumnGrid::*ColumnGrid_wrap_getColDen_float)(float, float, string) = &ColumnGrid::getColDen;

    class_<ColumnGrid, boost::noncopyable, bases<SimSnap> >("ColumnGrid",init<SimSnap&,float,float,int,float,float,int,optional<bool> >()
							    [with_custodian_and_ward<1,2>()])
      .def("__getitem__", &wrap_ColumnGrid_get, return_internal_reference<1>())
      .def("getColDen", ColumnGrid_wrap_getColDen_int)
      .def("getColDen", ColumnGrid_wrap_getColDen_float)
      ;

    class_<ColumnList, boost::noncopyable, bases<SimSnap> >("ColumnList",no_init)
      .def("__getitem__", &wrap_ColumnList_get, return_internal_reference<1>())
      ;

    def("loadFile", &SimSnap::loadFile, return_value_policy<manage_new_object>());
    
    def("lumToMag",&siman::lumToMag);
    def("magToLum",&siman::magToLum);
    def("voigt",&siman::voigt);
    def("absorption_line",&siman::absorption_line);

    class_<Filter>("Filter",no_init)
      .def(self & self)
      .def(self | self)
      .def(!self)
      ;

    class_<Sphere, bases<Filter> >("SphereFilter", init<float,float,float,float>() )
      .def(init<float>()) 
      ;

    
    class_<VelocitySphereFilter, bases<Filter> >("VelocitySphereFilter", init<float,float,float,float>() )
      .def(init<float>())
      ;

    class_<ParticleTypeFilter, bases<Filter> >("ParticleTypeFilter", init<int>());

    class_<RandomFilter, bases<Filter> >("RandomFilter", init<float>()) ;

    class_<AndFilter, bases<Filter> >("AndFilter",init<Filter &, Filter &>());
    class_<OrFilter, bases<Filter> >("OrFilter",init<Filter &, Filter &>());
    class_<NotFilter, bases<Filter> >("NotFilter",init<Filter &>());
    
    class_<Column, bases<Filter> >("ColumnFilter",init<float,float,float,float>());
    
    class_<Subset, boost::noncopyable, bases<SimSnap> >("Subset", init<SimSnap &, Filter &>()
							[with_custodian_and_ward<1,2>()]
							// keep parent SimSnap alive while Subset is
							);

    class_<Transformation>("Transformation", no_init)
      .def(self * self)
      ;

    class_<TwoTransformation, bases<Transformation> >("TwoTransformation", init<const Transformation &, const Transformation &>());

    class_<OrthoTrans, bases<Transformation> >("OrthoTrans", init<>())
      .def(init<float,float,float,float,float,float,float,float,float>());

    class_<RotationX, bases<OrthoTrans> >("RotationX", init<float>());
    class_<RotationY, bases<OrthoTrans> >("RotationY", init<float>());
    class_<RotationZ, bases<OrthoTrans> >("RotationZ", init<float>());

    class_<Translation, bases<Transformation> >("Translation", init<float,float,float>());
    class_<Wrap, bases<Transformation> >("Wrap", init<float,float>());
      
    class_<MassTransformation, bases<Transformation> >("MassTransformation", init<int>())
    .def(init<int,string,string>()) ;
    class_<LuminosityTransformation, bases<Transformation> >("LuminosityTransformation", init<string>());
    
    
    object particle_class = 
      class_<Particle>("Particle")
      .def(init<float,float,float>())
      .def_readwrite("x",&Particle::x)
      .def_readwrite("y",&Particle::y)
      .def_readwrite("z",&Particle::z)
      .def_readwrite("vx",&Particle::vx)
      .def_readwrite("vy",&Particle::vy)
      .def_readwrite("vz",&Particle::vz)
      .def_readwrite("mass",&Particle::mass)
      .def_readwrite("type",&Particle::type)
      .def_readwrite("rho",&Particle::rho)
      .def_readwrite("ne",&Particle::ne)
      .def_readwrite("u",&Particle::u)
      .def_readwrite("temp",&Particle::temp)
      .def_readwrite("metal",&Particle::metal)
      ;

    
    particle_class.attr("gas") = 1;
    particle_class.attr("dm") = 2;
    particle_class.attr("star") = 4;
    
    Unit (*pow_Unit_Rational)(const Unit &u, const Rational &r) = &siman::pow;
  
    class_<Unit>("Unit",init<string>())
      .def("__str__", &Unit::toString)
      .def("__pow__", pow_Unit_Rational)
      .def("convertTo", &wrap_Unit_convertTo_2)
      .def("convertTo", &wrap_Unit_convertTo_1)
      .def(self * self)
      .def(self / self)
      .def(self * float())
      .def(self / float())
     
      ;
    
#ifdef SIMAN_VIS


    void (Visualiser::*Visualiser_wrap_annotateVector_1)(SimanVec)  = &Visualiser::annotateVector;
    void (Visualiser::*Visualiser_wrap_annotateVector_2)(SimanVec, SimanVec)  = &Visualiser::annotateVector;
    void (Visualiser::*Visualiser_wrap_annotatePlane_1)(SimanVec)  = &Visualiser::annotatePlane;
    void (Visualiser::*Visualiser_wrap_annotatePlane_2)(SimanVec, SimanVec )  = &Visualiser::annotatePlane;

    SimanVec (Visualiser::*Visualiser_wrap_getCameraPos)(void )  = &Visualiser::getCameraPos;


    object visualiser_class = 
      class_<Visualiser, boost::noncopyable, bases<SimanObject> >("Visualiser",no_init)
      .add_property("colourMap", 
		    make_function(&Visualiser::getColourMap, return_internal_reference<1>()), 
		    make_function(&Visualiser::setColourMap, with_custodian_and_ward<1,2>()))
      .add_property("simSnap",
		    make_function(&Visualiser::getSimSnap, return_internal_reference<1>()),
		    make_function(&Visualiser::setSimSnap, with_custodian_and_ward<1,2>()))
      .add_property("cm", 
		    make_function(&Visualiser::getColourMap, return_internal_reference<1>()), 
		    make_function(&Visualiser::setColourMap, with_custodian_and_ward<1,2>()))
      .add_property("sim",
		    make_function(&Visualiser::getSimSnap, return_internal_reference<1>()),
		    make_function(&Visualiser::setSimSnap, with_custodian_and_ward<1,2>()))
      .def("transformSim", &Visualiser::transformSim)
      .def("setFlag", &Visualiser::setFlag)
      .def("hasFlag", &Visualiser::hasFlag)
      .def("toggleFlag", &Visualiser::toggleFlag)
      .def("av",Visualiser_wrap_annotateVector_1)
      .def("av",Visualiser_wrap_annotateVector_2)
      .def("ap",Visualiser_wrap_annotatePlane_1)
      .def("ap",Visualiser_wrap_annotatePlane_2)
      .def("pa",&Visualiser::popAnnotate)
      .add_property("camera", Visualiser_wrap_getCameraPos, &Visualiser::setCameraPos)
      ;
    
    visualiser_class.attr("autoRotate") = 1;
    visualiser_class.attr("lockZoom") = 4;
    visualiser_class.attr("noParticles") = 8;
    visualiser_class.attr("plotMeta") = 16;
    visualiser_class.attr("ortho") =  32;
    visualiser_class.attr("axes") = 64;
    visualiser_class.attr("vectors") = 128;
    visualiser_class.attr("blendMax") = 256;
    visualiser_class.attr("noAutoSimp") = 1024;
    visualiser_class.attr("dispColourBar") = 2048;
    visualiser_class.attr("promo") = 4096; 
    visualiser_class.attr("splatting") = 8192;
    visualiser_class.attr("noContextMenu") = 1<<14;
    
    void (GlutInterface::*GlutInterface_wrap_bind)(string,string) = &GlutInterface::bind;

    class_<GlutInterface, boost::noncopyable, bases<Visualiser> >("GlutInterface",init<SimSnap *>()
								  [with_custodian_and_ward<1,2>()]) 
      .def("bind", GlutInterface_wrap_bind)
      .def("run",&GlutInterface::run)
      .def("setup",&GlutInterface::setup)
      .def("poll",&GlutInterface::poll)
      .add_property("camera", Visualiser_wrap_getCameraPos, &GlutInterface::setCameraPos)
      .add_property("lengthScale", &Visualiser::getApproxLengthScale, &GlutInterface::setLengthScale)
      ;

    class_<ColourMap, boost::noncopyable, bases<SimanObject> >("ColourMap",init<Visualiser *>())
      .add_property("col",&ColourMap::getCol,&ColourMap::setCol)
      .add_property("display",&ColourMap::getDisplay,&ColourMap::setDisplay)
      ;

    void (ContinuuousColourMap::*ContinuuousColourMap_wrap_setColourBy)(string) = &ContinuuousColourMap::setColourBy;
     
    class_<ContinuuousColourMap, boost::noncopyable, bases<ColourMap> >("ContinuuousColourMap",no_init)
      .add_property("range",&ContinuuousColourMap::getRange,&ContinuuousColourMap::setRange)
      .add_property("logScale",&ContinuuousColourMap::getLogScale,&ContinuuousColourMap::setLogScale)
      .add_property("hideOutOfRange",&ContinuuousColourMap::getHideOutOfRange, &ContinuuousColourMap::setHideOutOfRange)
      .def("setColourBy",ContinuuousColourMap_wrap_setColourBy)
      ;
    
    class_<ColourMapGradient, boost::noncopyable, bases<ContinuuousColourMap> >("ColourMapGradient",
										init<Visualiser *, vector<double>, vector<double> >()
										[with_custodian_and_ward<1,2>()])
      .add_property("col1",&ColourMapGradient::getCol1,&ColourMapGradient::setCol1)
      .add_property("col2",&ColourMapGradient::getCol2,&ColourMapGradient::setCol2)
      ;
    
    class_<ColourMapSpectrum, boost::noncopyable, bases<ContinuuousColourMap> >("ColourMapSpectrum",init<Visualiser *>()
										[with_custodian_and_ward<1,2>()])
      ;
    
    class_<ColourMapByType, boost::noncopyable, bases<ColourMap> >("ColourMapByType",
								   init<Visualiser *,ColourMap *, ColourMap *, ColourMap *>()
								   [with_custodian_and_ward<1,2,
								    with_custodian_and_ward<1,3,
								    with_custodian_and_ward<1,4,
								    with_custodian_and_ward<1,5> > > >()])
      .add_property("gas",
		    make_function(&ColourMapByType::getGasMap,return_internal_reference<1>())
		    )
      .add_property("dm",
		    make_function(&ColourMapByType::getDMMap,return_internal_reference<1>())
		    )
      .add_property("star",
		    make_function(&ColourMapByType::getStarMap,return_internal_reference<1>()),
		    make_function(&ColourMapByType::setStarMap,with_custodian_and_ward<1,2>())
		    )
      ;
#endif 

#ifdef SIMAN_EXTRA
    
    object ionise_class = class_<Ionise>("Ionise",init<SimSnap*>()[with_custodian_and_ward<1,2>()])
      .def(init<float,float,SimSnap*>())
      .def(init<float,float,float,SimSnap*>())
      .def("thinRadiative",&Ionise::thinRadiative)
      .def("haehneltRadiative",&Ionise::haehneltRadiative)
      .def("thickRadiative",&Ionise::thickRadiative)
      .def("setFlag",&Ionise::setFlag)
      .def("unsetFlag",&Ionise::unsetFlag)
      .def("cooling",&Ionise::cooling)
      ;

    ionise_class.attr("verbose")=16;
    ionise_class.attr("writeGam")=4096;
    ionise_class.attr("useTemp")=8192;
    
    class_<Subsets>("Subsets",no_init)
      .def("getGroup",&Subsets::getGroup,with_custodian_and_ward_postcall<0,1>())
      .def("getNumGroups",&Subsets::getNumGroups)
      .def("__len__",&Subsets::getNumGroups)
      .def("__getitem__",&Subsets::getGroup,with_custodian_and_ward_postcall<0,1>())
      ;

    class_<GrpFile, bases<Subsets> >("GrpFile",init<SimSnap*,std::string>()[with_custodian_and_ward<1,2>()])
      ;

    class_<FoFFile, bases<Subsets> >("FoFFile",init<SimSnap*, char*, char*, int>()[with_custodian_and_ward<1,2>()])
      ;
#endif

    def("setVerbose",&setVerbose);
    def("getVerbose",&getVerbose);

    def("_siman_internal_fetch_manage" , &fetch, return_value_policy<manage_new_object>());
 
    def("_siman_internal_fetch_reference" , &fetch, return_value_policy<reference_existing_object>()); // dangerous but useful e.g. Visualise takes responsibility for its own var


    unit_from_python_str();

    boost::python::to_python_converter<
      std::vector<double>,
      scitbx::boost_python::container_conversions::to_tuple<
      std::vector<double> > > ();

    boost::python::to_python_converter<
      SimanVec,
      scitbx::boost_python::container_conversions::to_tuple<
      SimanVec > > ();
    
    boost::python::to_python_converter<
      std::vector<float>,
      scitbx::boost_python::container_conversions::to_tuple<
      std::vector<float> > > ();

    scitbx::boost_python::container_conversions::from_python_sequence<
    std::vector<double>,
      scitbx::boost_python::container_conversions::variable_capacity_policy>();


    scitbx::boost_python::container_conversions::from_python_sequence<
    std::vector<float>,
      scitbx::boost_python::container_conversions::variable_capacity_policy>();


    scitbx::boost_python::container_conversions::from_python_sequence<
    SimanVec,
      scitbx::boost_python::container_conversions::variable_capacity_policy>();

    
    
  }
    
  
}


}
