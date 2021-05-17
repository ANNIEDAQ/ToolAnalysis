#ifndef BeamDecoder_H
#define BeamDecoder_H

#include <string>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <sstream>

#include "Tool.h"
#include "BeamStatus.h"
#include "MinibufferLabel.h"
#include "ANNIEconstants.h"
#include "BeamDataPoint.h"
#include "TimeClass.h"



/**
 * \class BeamDecoder
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M. Nieslony $
* $Date: 2021/02/17 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/
class BeamDecoder: public Tool {


 public:

  BeamDecoder(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool initialise_beam_db();
  BeamStatus get_beam_status(uint64_t ns_since_epoch, MinibufferLabel mb_label);
  

 private:

  int verbosity_;
  bool first_entry;

  bool NewCTCDataAvailable;
  std::map<uint64_t,std::vector<uint32_t>>* TimeToTriggerWordMap;	//Trigger information (needed to check the beam DB)
  std::map<uint64_t,BeamStatus> *BeamStatusMap;				//Map containing the beam status information

  std::string horn_current_device;
  std::string first_toroid;
  std::string second_toroid;

  BoostStore beam_db_store_;
  std::map<int, std::pair<uint64_t, uint64_t>> beam_db_index_;
  uint64_t start_ms_since_epoch_;
  uint64_t end_ms_since_epoch_;

};


#endif
