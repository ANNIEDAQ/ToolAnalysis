#ifndef LoadRunInfo_H
#define LoadRunInfo_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class LoadRunInfo
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M. Nieslony $
* $Date: 2021/02/21 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/
class LoadRunInfo: public Tool {


 public:

  LoadRunInfo(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:

  //Configuration variables
  int verbosity = 2;
  std::string runinfofile;

  //Verbosity variables
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;



};


#endif
