#ifndef MonitorReceive_H
#define MonitorReceive_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "boost/date_time/gregorian/gregorian.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/chrono.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class                                     
#include <boost/uuid/uuid_generators.hpp> // generators                                 
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.        
#include <boost/filesystem.hpp>		  // Get file sizes

class MonitorReceive: public Tool {


 public:

  MonitorReceive();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  int UpdateMonitorSources();

 private:

  int sources;
  zmq::socket_t* MonitorReceiver;
  boost::posix_time::ptime last;
  boost::posix_time::time_duration period;
  zmq::pollitem_t items[1];
  std::map<std::string,Store*> connections; 
  BoostStore* indata;
  BoostStore* MRDData;
  BoostStore* PMTData;
  BoostStore* TrigData;
  std::vector<std::string> loaded_files;
};


#endif
