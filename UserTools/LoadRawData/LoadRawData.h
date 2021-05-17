#ifndef LoadRawData_H
#define LoadRawData_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "CardData.h"
#include "TriggerData.h"
#include "BoostStore.h"
#include "Store.h"

/**
 * \class LoadRawData
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class LoadRawData: public Tool {


 public:

  LoadRawData(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  void LoadPMTMRDData(); 
  void LoadTriggerData(); 
  void LoadRunInformation();
  void GetNextDataEntries();
  bool InitializeNewFile(); 
  int GetRunFromFilename();
  int GetSubRunFromFilename();
  int GetPartFromFilename();

 private:

  bool DummyRunInfo;
  std::vector<std::string> OrganizedFileList;
  std::string CurrentFile = "NONE";
  std::string BuildType;
  std::string Mode;
  std::string InputFile;
  std::vector<std::string> OrganizeRunParts(std::string InputFile); //Parses all run files in InputFile and returns a vector of file paths organized by part
  bool readtrigoverlap;
  bool storetrigoverlap;
  bool storerawdata;

  int FileNum = 0;
  int tanktotalentries;
  int trigtotalentries;
  int mrdtotalentries;
  bool TankEntriesCompleted;
  bool MRDEntriesCompleted;
  bool TrigEntriesCompleted;
  bool FileCompleted;
  int TankEntryNum = 0;
  int MRDEntryNum = 0;
  int TrigEntryNum = 0;

  //Run / Subrun / part info
  int extract_run;
  int extract_subrun;
  int extract_part;

  //Bools modified to determine which Decoder tools are paused downstream
  bool TankPaused;
  bool MRDPaused;
  bool CTCPaused;

  BoostStore *RawData = nullptr;
  BoostStore *PMTData = nullptr;
  BoostStore *MRDData = nullptr;
  BoostStore *TrigData = nullptr;
  std::vector<CardData>* Cdata = nullptr;
  //std::vector<TriggerData>* Tdata = nullptr;
  TriggerData* Tdata = nullptr;
  MRDOut* Mdata = nullptr;

  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  int vv_debug=4;
  std::string logmessage;

};


#endif
