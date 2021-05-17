#ifndef TriggerDataDecoder_H
#define TriggerDataDecoder_H

#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "Tool.h"
#include "TriggerData.h"

/**
 * \class TriggerDataDecoder
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class TriggerDataDecoder: public Tool {

 public:

  TriggerDataDecoder(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool AddWord(uint32_t word);
  std::vector<int> LoadTriggerMask(std::string triggermask_file);
  std::map<int,std::string> LoadTriggerWords(std::string triggerwords_file);
  void CheckForRunChange();
 private:

  //std::vector<TriggerData> *Tdata = nullptr;
  TriggerData *Tdata = nullptr;
  std::map<uint64_t,std::vector<uint32_t>>* TimeToTriggerWordMap;
  std::map<uint64_t,std::vector<uint32_t>>* TimeToTriggerWordMapComplete;	//Info about all triggerwords
  bool UseTrigMask = false;
  std::string TriggerWordFile;
  std::string mode;
  bool storetrigoverlap;
  bool readtrigoverlap;
  int loop_nr;

  std::string TriggerMaskFile;
  std::vector<int> TriggerMask;

  std::string TriggerWordsFile;
  std::map<int,std::string> TriggerWords;

  std::vector<int> fiforesets;
  std::vector<uint32_t> processed_sources;
  std::vector<uint64_t> processed_ns;
  bool have_c1 = false;
  bool have_c2 = false;
  uint64_t c1 = 0;
  uint64_t c2 = 0;
  int CurrentRunNum;
  int CurrentSubrunNum;
  int CurrentPartNum;
  
  bool usecstore;

  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  int vv_debug=4;
  std::string logmessage;
};


#endif
