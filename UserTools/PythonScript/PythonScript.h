#ifndef PythonScript_H
#define PythonScript_H

#include <string>
#include <iostream>
#include <Python.h>
#include <PythonAPI.h>

#include "Tool.h"

class PythonScript: public Tool {


 public:

  PythonScript();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  std::string pythonscript;
  std::string initialisefunction;
  std::string executefunction;
  std::string finalisefunction;

  PyObject *pName, *pModule, *pFuncI, *pFuncE, *pFuncF;
  PyObject *pArgs, *pValue;
  PyThreadState* pythread;

  Store thisscriptsconfigstore;

  int pyinit;

};


#endif
