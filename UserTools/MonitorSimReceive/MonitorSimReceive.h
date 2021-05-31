#ifndef MonitorSimReceive_H
#define MonitorSimReceive_H

#include <string>
#include <iostream>
#include <stdlib.h>     
#include <time.h>   

#include "MRDOut.h"
#include "Tool.h"

#include "boost/date_time/gregorian/gregorian.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/chrono.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <boost/filesystem.hpp>           // Get file sizes


class MonitorSimReceive: public Tool {


 public:

  MonitorSimReceive();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  std::string file_list;
  std::string mode;
  std::string outpath;
  int verbosity;
  
  BoostStore *indata = nullptr;
  BoostStore* MRDData = nullptr;
  BoostStore* PMTData = nullptr;
  BoostStore* TrigData = nullptr;

  std::vector<std::string> vec_filename;
  int i_loop;


};


#endif
