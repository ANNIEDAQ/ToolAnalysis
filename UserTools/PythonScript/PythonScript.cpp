#include "PythonScript.h"
PyMODINIT_FUNC test(void){

  return PyModule_Create(&StoreModule);

  }

PythonScript::PythonScript():Tool(){}


bool PythonScript::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("PythonScript",pythonscript);
  m_variables.Get("InitialiseFunction",initialisefunction);
  m_variables.Get("ExecuteFunction",executefunction);
  m_variables.Get("FinaliseFunction",finalisefunction);
  std::string scriptconfigfile="";
  m_variables.Get("ConfigurationsFile",scriptconfigfile);

  // Initialize a Store with the python script's variables
  thisscriptsconfigstore.Initialise(scriptconfigfile);

  PyImport_AppendInittab("Store", test);
  
  gstore=m_data;
  // each time we're running this python script, set the global gconfig to point to it's configurations store
  // so that calls to the PythonAPI for config values will return values from the correct store.
  gconfig=&thisscriptsconfigstore;

  // Initialising Python
  pyinit=0;
  if(!(m_data->CStore.Get("PythonInit",pyinit))){
    Py_Initialize();
  }

  //// if multiple python scipts this is to keep track of which is last for finalisation
  pyinit++;
  m_data->CStore.Set("PythonInit",pyinit);
  
  /// Starting python thread for this tool
  //pythread=Py_NewInterpreter();  
  //PyThreadState_Swap(pythread);

  // Loading store API into python env
  //Py_InitModule("Store", StoreMethods); ///python 2 version
  //PyModule_Create(&StoreMethods); ///python 3 version
  //PyImport_AppendInittab("Store", test);
  PyImport_ImportModule("Store");


  // Loading python script/module
  //  std::cout<<pythonscript.c_str()<<std::endl;
  //pName = PyString_FromString(pythonscript.c_str()); //python 2
  pName = PyUnicode_FromString(pythonscript.c_str());

  // Error checking of pName to be added later
  pModule = PyImport_Import(pName);
  //pModule = PyImport_ReloadModule(pName); 
  
  // deleting pname
  Py_DECREF(pName);

  //loading python fuctions
  if (pModule != NULL) {

    pFuncI = PyObject_GetAttrString(pModule, initialisefunction.c_str());
    pFuncE = PyObject_GetAttrString(pModule, executefunction.c_str());
    pFuncF = PyObject_GetAttrString(pModule, finalisefunction.c_str());

    if (pFuncI && pFuncE && pFuncF && PyCallable_Check(pFuncI) && PyCallable_Check(pFuncE) && PyCallable_Check(pFuncF)) {
      //pArgs = PyTuple_New(0);
       pArgs = Py_BuildValue("(i)",pyinit);
      
      /* we dont really want to pass any arguments but have the option here
	 for (int i = 0; i < 0; ++i) {
	 pValue = PyLong_FromLong(i);
	 if (!pValue) {
	 Py_DECREF(pArgs);
	 Py_DECREF(pModule);
	 fprintf(stderr, "Cannot convert argument\n");
	 return 1;
	 }
	 // pValue reference stolen here: 
	 PyTuple_SetItem(pArgs, i, pValue);
	 }*/
      
      
      // run python Initialise function
      pValue = PyObject_CallObject(pFuncI, pArgs);
      Py_DECREF(pArgs); //delete args after running
      
      
      if(pValue != NULL){
	if (!(PyLong_AsLong(pValue))){
	  std::cout<<"Python script returned internal error in initialise "<<std::endl;	
	  return false;
	}
	Py_DECREF(pValue);
      }
      else { //something went wrong with function call
	Py_DECREF(pFuncI);
	Py_DECREF(pModule);
	PyErr_Print();
	fprintf(stderr,"Python Initialise call failed\n");
	return false;
      }
    }
    else {
      if (PyErr_Occurred())
	PyErr_Print();
      fprintf(stderr, "Cannot find python functions");
      Py_XDECREF(pFuncE);
      Py_XDECREF(pFuncF);
      Py_DECREF(pModule);
    }
    Py_XDECREF(pFuncI);
    
  }
  else {
    
    PyErr_Print();
    fprintf(stderr, "Failed to load python script/module");
    return false;
  }
  
  
  
  
  return true;
}


bool PythonScript::Execute(){

  // make the config Store for this script accessible to the tool, should it need it in Execute
  gconfig=&thisscriptsconfigstore;

  //PyThreadState_Swap(pythread);

  if (pModule != NULL) {

    if (pFuncE && PyCallable_Check(pFuncE)) {
      pArgs = PyTuple_New(0);//no arguments
      
      pValue = PyObject_CallObject(pFuncE, pArgs); // call the function
      Py_DECREF(pArgs); // delete arguments
      
      if(pValue != NULL){
	if (!(PyLong_AsLong(pValue))){
	  std::cout<<"Python script returned internal error in execute "<<std::endl;
	  return false;
	}
	Py_DECREF(pValue);
      }
      else {
        Py_DECREF(pFuncE);
        Py_DECREF(pModule);
        PyErr_Print();
        fprintf(stderr,"Call to Python Execute failed\n");
        return false;
      }
    }
    else {
      if (PyErr_Occurred())
        PyErr_Print();
      fprintf(stderr, "Cannot python execute function \n");
      Py_XDECREF(pFuncE);
      Py_DECREF(pModule);
      return false;
    }
    
  }
 
  
  return true;
}


bool PythonScript::Finalise(){

  // make the config Store for this script accessible to the tool, should it need it in Finalise
  gconfig=&thisscriptsconfigstore;
  
  //PyThreadState_Swap(pythread);  
  
  if (pModule != NULL) {
    
    if (pFuncF && PyCallable_Check(pFuncF)) {
      pArgs = PyTuple_New(0);
      
      pValue = PyObject_CallObject(pFuncF, pArgs);
      Py_DECREF(pArgs);
      
      if (pValue != NULL) {
	if (!(PyLong_AsLong(pValue))){
	  std::cout<<"Python script returned internal error in finalise "<<std::endl;
          return false;
        }
        Py_DECREF(pValue);
      }
      else {
        Py_DECREF(pFuncF);
        Py_DECREF(pModule);
        PyErr_Print();
        fprintf(stderr,"Call to Python Finalise failed\n");
        return false;
      }
    }
    else {
      if (PyErr_Occurred())
        PyErr_Print();
      fprintf(stderr, "Cannot python finalise function \n");
      Py_XDECREF(pFuncF);
      Py_DECREF(pModule);
      return false;
    }

  }
  
  int tmpinit=0;
  m_data->CStore.Get("PythonInit",tmpinit);
  if(tmpinit==pyinit)  Py_Finalize();
  
  return true;
}
