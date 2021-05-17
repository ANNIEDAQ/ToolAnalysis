#include "MRDDataDecoder.h"

MRDDataDecoder::MRDDataDecoder():Tool(){}


bool MRDDataDecoder::Initialise(std::string configfile, DataModel &data){
  std::cout << "Initializing MRDDataDecoder tool" << std::endl;
  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  verbosity = 0;
  DaylightSavings = 0;	//Make the default no correction for daylight savings (mismatch between UTC & Fermi time is only a few weeks in the year)

  m_variables.Get("verbosity",verbosity);
  // was this run taken at a daylight savings time of year?
  m_variables.Get("DaylightSavingsSpring",DaylightSavings);
  // work out the conversion required to put MRD timestamps into UTC
  TimeZoneShift = 21600000;
  if(DaylightSavings) TimeZoneShift = 18000000;

  if (verbosity > 0) std::cout <<"MRDDataDecoder: TimeZoneShift: "<<TimeZoneShift<<std::endl;

  m_data->CStore.Get("MRDCrateSpaceToChannelNumMap",MRDCrateSpaceToChannelNumMap);
  m_data->CStore.Set("NewMRDDataAvailable",false);

  m_data->CStore.Set("PauseMRDDecoding",false);
  Log("MRDDataDecoder Tool: Initialized successfully",v_message,verbosity);
  return true;
}


bool MRDDataDecoder::Execute(){
  m_data->CStore.Set("NewMRDDataAvailable",false);

  bool NewEntryAvailable;
  m_data->CStore.Get("NewRawDataEntryAccessed",NewEntryAvailable);
  if(!NewEntryAvailable){ //Something went wrong processing raw data.  Stop and save what's left
    Log("MRDDataDecoder Tool: There's no new MRD data.  stop at next loop to save what data has been built.",v_warning,verbosity); 
    m_data->vars.Set("StopLoop",1);
    return true;
  }
  
  bool PauseMRDDecoding = false;
  m_data->CStore.Get("PauseMRDDecoding",PauseMRDDecoding);
  if (PauseMRDDecoding){
    std::cout << "MRDDataDecoder tool: Pausing MRD decoding to let Tank data catch up..." << std::endl;
    return true;
  }


  /////////////////// getting MRD Data ////////////////////
  Log("MRDDataDecoder Tool: Accessing MRDData from CStore",v_message,verbosity); 
  m_data->CStore.Get("MRDData",mrddata);
  std::string mrdTriggertype = "No Loopback";
  std::vector<unsigned long> chankeys;
  uint64_t timestamp = static_cast<uint64_t>(mrddata->TimeStamp);    //in ms since 1970/1/1
  // before anything else convert it to UTC ns
  timestamp = (timestamp+TimeZoneShift)*1E6;
  std::vector<std::pair<unsigned long, int>> ChankeyTimePairs;
  MRDEvents.emplace(timestamp,ChankeyTimePairs);
  
  bool cosmic_loopback = false;
  bool beam_loopback = false;
  int cosmic_tdc = -1;
  int beam_tdc = -1;
  std::vector<int> CrateSlotChannel_Beam{7,11,15};
  std::vector<int> CrateSlotChannel_Cosmic{7,11,14};
    
  //For each entry, loop over all crates and get data
  for (unsigned int i_data = 0; i_data < mrddata->Crate.size(); i_data++){
    int crate = mrddata->Crate.at(i_data);
    int slot = mrddata->Slot.at(i_data);
    int channel = mrddata->Channel.at(i_data);
    int hittimevalue = mrddata->Value.at(i_data);
    std::vector<int> CrateSlotChannel{crate,slot,channel};
    unsigned long chankey = 999;
    if (MRDCrateSpaceToChannelNumMap.find(CrateSlotChannel) != MRDCrateSpaceToChannelNumMap.end()){
      chankey = MRDCrateSpaceToChannelNumMap[CrateSlotChannel];
    }
    //std::cout <<"crate: "<<crate<<", slot: "<<slot<<", channel: "<<channel<<", chankey: "<<chankey<<std::endl;
    if (CrateSlotChannel != CrateSlotChannel_Beam && CrateSlotChannel != CrateSlotChannel_Cosmic && chankey != 999){
      std::pair <unsigned long,int> keytimepair(chankey,hittimevalue);  //chankey will be 0 when looking at loopback channels that don't have an entry in the mapping-->skip
      MRDEvents[timestamp].push_back(keytimepair);
    }
    if (crate == 7 && slot == 11 && channel == 14) {cosmic_loopback=true; cosmic_tdc = hittimevalue;}   //FIXME: don't hard-code the trigger channels?
    if (crate == 7 && slot == 11 && channel == 15) {beam_loopback=true; beam_tdc = hittimevalue;}     //FIXME: don't hard-code the trigger channels?
  }
  
  if (beam_loopback) mrdTriggertype = "Beam";
  if (cosmic_loopback) mrdTriggertype = "Cosmic";      //prefer cosmic loopback over beam loopback (cosmic event will always also have a beam loopback entry)

  CosmicLoopbackMap.emplace(timestamp,cosmic_tdc);
  BeamLoopbackMap.emplace(timestamp,beam_tdc);

  //Entry processing done.  Label the trigger type and increment index
  TriggerTypeMap.emplace(timestamp,mrdTriggertype);

  //MRD Data file fully processed.   
  //Push the map of TriggerTypeMap and FinishedMRDHits 
  //to the CStore for ANNIEEvent to start Building ANNIEEvents. 
  Log("MRDDataDecoder Tool: Saving Finished MRD Data into CStore.",v_debug, verbosity);
  //FIXME: add a check for if there is or is not an entry in the CStore
  m_data->CStore.Get("MRDEvents",CStoreMRDEvents);
  CStoreMRDEvents.insert(MRDEvents.begin(),MRDEvents.end());
  m_data->CStore.Set("MRDEvents",CStoreMRDEvents);

  m_data->CStore.Get("MRDEventTriggerTypes",CStoreTriggerTypeMap);
  CStoreTriggerTypeMap.insert(TriggerTypeMap.begin(),TriggerTypeMap.end());
  m_data->CStore.Set("MRDEventTriggerTypes",CStoreTriggerTypeMap);
 
  m_data->CStore.Get("MRDBeamLoopback",CStoreBeamLoopbackMap);
  CStoreBeamLoopbackMap.insert(BeamLoopbackMap.begin(),BeamLoopbackMap.end());
  m_data->CStore.Set("MRDBeamLoopback",CStoreBeamLoopbackMap);
 
  m_data->CStore.Get("MRDCosmicLoopback",CStoreCosmicLoopbackMap);
  CStoreCosmicLoopbackMap.insert(CosmicLoopbackMap.begin(),CosmicLoopbackMap.end());
  m_data->CStore.Set("MRDCosmicLoopback",CStoreCosmicLoopbackMap);
  
  m_data->CStore.Set("NewMRDDataAvailable",true);

  //Check the size of the WaveBank to see if things are bloating
  Log("MRDDataDecoder Tool: Size of MRDEvents in CStore (# MRD Triggers processed):" + 
          to_string(CStoreMRDEvents.size()),v_debug, verbosity);
  Log("MRDDataDecoder Tool: Size of TriggerTypeMap in CStore:" + 
          to_string(CStoreTriggerTypeMap.size()),v_debug, verbosity);

  std::cout << "MRD EVENT CSTORE ENTRIES SET SUCCESSFULLY.  Clearing MRDEvent vector in MRDDataDecoder tool." << std::endl;
  MRDEvents.clear();
  TriggerTypeMap.clear();
  BeamLoopbackMap.clear();
  CosmicLoopbackMap.clear();

  ////////////// END EXECUTE LOOP ///////////////
  return true;
}


bool MRDDataDecoder::Finalise(){
  Log("MRDDataDecoder tool exitting",v_message,verbosity);
  return true;
}
