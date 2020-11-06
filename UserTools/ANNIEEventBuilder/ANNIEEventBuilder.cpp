#include "ANNIEEventBuilder.h"

ANNIEEventBuilder::ANNIEEventBuilder():Tool(){}


bool ANNIEEventBuilder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer

  SavePath = "./";
  ProcessedFilesBasename = "ProcessedRawData";
  BuildType = "TankAndMRD";
  EventsPerPairing = 200;
  NumWavesInCompleteSet = 140;
  OrphanOldTankTimestamps = true;
  OldTimestampThreshold = 120; //seconds
  OrphanWarningValue = 20;
  ExecutesPerBuild = 50;
  MRDTankTimeTolerance = 10;     //ms
  CTCTankTimeTolerance = 100;    //ns
  CTCMRDTimeTolerance = 2000000; //ns
  DriftWarningValue = 5000000;   //ns
  pause_threshold = 5*60;        //s

  /////////////////////////////////////////////////////////////////
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("SavePath",SavePath);
  m_variables.Get("ProcessedFilesBasename",ProcessedFilesBasename);
  m_variables.Get("BuildType",BuildType);
  m_variables.Get("NumEventsPerPairing",EventsPerPairing);
  m_variables.Get("MinNumWavesInSet",NumWavesInCompleteSet);
  m_variables.Get("OrphanOldTankTimestamps",OrphanOldTankTimestamps);
  m_variables.Get("OldTimestampThreshold",OldTimestampThreshold);
  m_variables.Get("ExecutesPerBuild",ExecutesPerBuild);
  m_variables.Get("MRDTankTimeTolerance",MRDTankTimeTolerance);
  m_variables.Get("CTCTankTimeTolerance",CTCTankTimeTolerance);
  m_variables.Get("CTCMRDTimeTolerance",CTCMRDTimeTolerance);
  m_variables.Get("OrphanFileBase",OrphanFileBase);
  m_variables.Get("MaxStreamMatchingTimeSeparation",pause_threshold);
  pause_threshold*=1E9;

  if(BuildType == "TankAndMRD" || BuildType == "TankAndMRDAndCTC"){
    std::cout << "BuildANNIEEvent Building Tank and MRD-merged ANNIE events. " <<
        std::endl;
  }
  else if(BuildType == "Tank" || BuildType == "MRD"){
    std::cout << "BuildANNIEEvent Building " << BuildType << "ANNIEEvents only."
         << std::endl;
  }
  else{
    std::cout << "BuildANNIEEvent ERROR: BuildType not recognized! " <<
        "Please select Tank, MRD, or TankAndMRD" << std::endl;
    return false;
  }

  m_data->CStore.Get("TankPMTCrateSpaceToChannelNumMap",TankPMTCrateSpaceToChannelNumMap);
  m_data->CStore.Get("AuxCrateSpaceToChannelNumMap",AuxCrateSpaceToChannelNumMap);
  m_data->CStore.Get("MRDCrateSpaceToChannelNumMap",MRDCrateSpaceToChannelNumMap);

  if(BuildType == "Tank" || BuildType == "TankAndMRD" || BuildType == "TankAndMRDAndCTC"){
    int NumTankPMTChannels = TankPMTCrateSpaceToChannelNumMap.size();
    int NumAuxChannels = AuxCrateSpaceToChannelNumMap.size();
    if(verbosity>4) std::cout << "TOTAL TANK + AUX CHANNELS: " << NumTankPMTChannels + NumAuxChannels << std::endl;
    if(verbosity>4) std::cout << "CURRENT SET THRESHOLD FOR BUILDING PMT EVENTS: " << NumWavesInCompleteSet << std::endl;
  }

  //////////////////////initialize subrun index//////////////
  ProcessedStore = new BoostStore(false,2);
  ANNIEEvent = new BoostStore(false,2);
  ANNIEEventNum = 0;
  CurrentRunNum = -1;
  CurrentSubRunNum = -1;
  CurrentDriftMean = 0;
  
  OrphanStore = new BoostStore(false,2);

  return true;
}


bool ANNIEEventBuilder::Execute(){
  bool NewEntryAvailable;
  m_data->CStore.Get("NewRawDataEntryAccessed",NewEntryAvailable);
  if(!NewEntryAvailable){ //Something went wrong processing raw data.  Stop and save what's left
    Log("ANNIEEventBuilder Tool: There's no new PMT/MRD data.  Stopping loop, ANNIEEvent BoostStore will save.",v_warning,verbosity); 
    m_data->vars.Set("StopLoop",1);
  }
  
  // ensure we always build if the toolchain is stopping
  bool toolchain_stopping=false;
  m_data->vars.Get("StopLoop",toolchain_stopping);
  // likewise ensure we try to do all possible matching if we've hit the end of the file
  bool file_completed=false;
  m_data->CStore.Set("FileCompleted",file_completed);
  toolchain_stopping |= file_completed;
  if(toolchain_stopping){
    Log("ANNIEEventBuilder: StopLoop or FileCompleted detected, forcing building of any remaining events in the timestream",v_warning,verbosity);
  }
    
  ExecuteCount+=1;
  if((ExecuteCount<ExecutesPerBuild)&&(!toolchain_stopping)) return true;

  //See if the MRD and Tank are at the same run/subrun for building
  m_data->CStore.Get("RunInfoPostgress",RunInfoPostgress);
  int RunNumber;
  int SubRunNumber;
  uint64_t StarTime;
  int RunType;
  RunInfoPostgress.Get("RunNumber",RunNumber);
  RunInfoPostgress.Get("SubRunNumber",SubRunNumber);
  RunInfoPostgress.Get("RunType",RunType);
  RunInfoPostgress.Get("StarTime",StarTime);
  //If our first execute loop, Initialize Current run information with Tank decoding progress
  if(CurrentRunNum == -1 || CurrentSubRunNum == -1){
    CurrentRunNum = RunNumber;
    CurrentSubRunNum = SubRunNumber;
    CurrentStarTime = StarTime;
    CurrentRunType = RunType;
  }
  //If we're in a new run or subrun, make a new ANNIEEvent file. 
  if((CurrentRunNum != RunNumber) || (CurrentSubRunNum != SubRunNumber)){
    this->OpenNewANNIEEvent(RunNumber,SubRunNumber,StarTime,RunType);
  }

  //Built ANNIE events with only PMT data
  if (BuildType == "Tank"){
    //Check to see if there's new PMT data
    m_data->CStore.Get("NewTankPMTDataAvailable",IsNewTankData);
    if((!IsNewTankData)&&(!toolchain_stopping)){
      Log("ANNIEEventBuilder:: No new Tank Data.  Not building ANNIEEvent. ",v_message, verbosity);
      return true;
    }
    else if(IsNewTankData) this->ProcessNewTankPMTData();
    this->ManageOrphanage();

    std::vector<uint64_t> PMTEventsToDelete;
    for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : FinishedTankEvents){
      uint64_t PMTCounterTime = apair.first;
      std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
      this->BuildANNIEEventRunInfo(RunNumber,SubRunNumber,RunType,StarTime);
      this->BuildANNIEEventTank(PMTCounterTime, aWaveMap);
      this->SaveEntryToFile(CurrentRunNum,CurrentSubRunNum);
      //Erase this entry from the InProgressTankEventsMap
      if(verbosity>4) std::cout << "Counter time will be erased from InProgressTankEvents: " << PMTCounterTime << std::endl;
      PMTEventsToDelete.push_back(PMTCounterTime);
    }
    //Update the current InProgressTankEvents map
    for(unsigned int i=0; i< PMTEventsToDelete.size(); i++) FinishedTankEvents.erase(PMTEventsToDelete.at(i));
  }
  /*
    //Get the current InProgressTankEvents map
    if(verbosity>4) std::cout << "ANNIEEventBuilder: Getting waves and run info from CStore" << std::endl;
    bool got_inprogresstanks = m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);
    if(!got_inprogresstanks){
      std::cout << "ANNIEEventBuilder ERROR: No InProgressTankEvents pointer!" << std::endl;
      return false;
    }
 
    //Assume a whole processed file will have all it's PMT data finished
    std::vector<uint64_t> PMTEventsToDelete;
    for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : *InProgressTankEvents){
      uint64_t PMTCounterTime = apair.first;
      if(verbosity>4) std::cout << "Finished waveset has clock counter: " << PMTCounterTime << std::endl;
      std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
      if(verbosity>4) std::cout << "Number of waves for this counter: " << aWaveMap.size() << std::endl;
      //For this counter, need to have the number of TankPMT channels plus number of aux channels
      if(aWaveMap.size() >= (NumWavesInCompleteSet)){
        this->BuildANNIEEventRunInfo(RunNumber,SubRunNumber,RunType,StarTime);
        this->BuildANNIEEventTank(PMTCounterTime, aWaveMap);
        this->SaveEntryToFile(CurrentRunNum,CurrentSubRunNum);
        //Erase this entry from the InProgressTankEventsMap
        if(verbosity>4) std::cout << "Counter time will be erased from InProgressTankEvents: " << PMTCounterTime << std::endl;
        PMTEventsToDelete.push_back(PMTCounterTime);
      }
    }

    //Update the current InProgressTankEvents map
    for(unsigned int i=0; i< PMTEventsToDelete.size(); i++) InProgressTankEvents->erase(PMTEventsToDelete.at(i));
    if(verbosity>3) std::cout << "Current number of unfinished PMT Waveforms: " << InProgressTankEvents->size() << std::endl;

  }
  */
  
  //Built ANNIE events with only MRD data
  else if (BuildType == "MRD"){
    m_data->CStore.Get("NewMRDDataAvailable",IsNewMRDData);
    std::vector<uint64_t> MRDEventsToDelete;
    if((!IsNewMRDData)&&(!toolchain_stopping)){
      Log("ANNIEEventBuilder:: No new MRD Data.  Not building ANNIEEvent: ",v_message, verbosity);
      return true;
    }
    m_data->CStore.Get("MRDEvents",myMRDMaps.MRDEvents);
    m_data->CStore.Get("MRDEventTriggerTypes",myMRDMaps.MRDTriggerTypeMap);
    m_data->CStore.Get("MRDBeamLoopback",myMRDMaps.MRDBeamLoopbackMap);
    m_data->CStore.Get("MRDCosmicLoopback",myMRDMaps.MRDCosmicLoopbackMap);
    int NumMRDTimestamps = myTimeStream.BeamMRDTimestamps.size();
    m_data->CStore.Set("NumMRDTimestamps",NumMRDTimestamps);
    
   //Loop through MRDEvents and process each into ANNIEEvent.
    for(std::pair<uint64_t,std::vector<std::pair<unsigned long,int>>> apair : myMRDMaps.MRDEvents){
      uint64_t MRDTimeStamp = apair.first;
      std::vector<std::pair<unsigned long,int>> MRDHits = apair.second;
      std::string MRDTriggerType = myMRDMaps.MRDTriggerTypeMap[MRDTimeStamp];
      int beam_tdc = myMRDMaps.MRDBeamLoopbackMap[MRDTimeStamp];
      int cosmic_tdc = myMRDMaps.MRDCosmicLoopbackMap[MRDTimeStamp];
      this->BuildANNIEEventRunInfo(RunNumber,SubRunNumber,RunType,StarTime);
      this->BuildANNIEEventMRD(MRDHits, MRDTimeStamp, MRDTriggerType, beam_tdc, cosmic_tdc);
      this->SaveEntryToFile(RunNumber,SubRunNumber);
      //Erase this entry from the InProgressTankEventsMap
      MRDEventsToDelete.push_back(MRDTimeStamp);
    }
    for (unsigned int i=0; i< MRDEventsToDelete.size(); i++){
      myMRDMaps.MRDEvents.erase(MRDEventsToDelete.at(i));
      myMRDMaps.MRDTriggerTypeMap.erase(MRDEventsToDelete.at(i));
      myMRDMaps.MRDBeamLoopbackMap.erase(MRDEventsToDelete.at(i));
      myMRDMaps.MRDCosmicLoopbackMap.erase(MRDEventsToDelete.at(i));
    }
    m_data->CStore.Set("MRDEvents",myMRDMaps.MRDEvents);
    m_data->CStore.Set("MRDEventTriggerTypes",myMRDMaps.MRDTriggerTypeMap);
    m_data->CStore.Set("MRDBeamLoopback",myMRDMaps.MRDBeamLoopbackMap);
    m_data->CStore.Set("MRDCosmicLoopback",myMRDMaps.MRDCosmicLoopbackMap);
  }

  //Built ANNIE events with Tank and MRD data
  else if (BuildType == "TankAndMRD"){
    //Check if any In-progress tank events now have all waveforms
    m_data->CStore.Get("NewTankPMTDataAvailable",IsNewTankData);
    if(IsNewTankData) this->ProcessNewTankPMTData();
    this->ManageOrphanage();

    //Look through our MRD data for any new timestamps
    m_data->CStore.Get("MRDEventTriggerTypes",myMRDMaps.MRDTriggerTypeMap);
    m_data->CStore.Get("NewMRDDataAvailable",IsNewMRDData);
    m_data->CStore.Get("MRDEvents",myMRDMaps.MRDEvents);
    m_data->CStore.Get("MRDBeamLoopback",myMRDMaps.MRDBeamLoopbackMap);
    m_data->CStore.Get("MRDCosmicLoopback",myMRDMaps.MRDCosmicLoopbackMap);
    if(IsNewMRDData) this->ProcessNewMRDData();
    

    
    //Now, pair up PMT and MRD events...
    int NumTankTimestamps = myTimeStream.BeamTankTimestamps.size();
    int NumMRDTimestamps = myTimeStream.BeamMRDTimestamps.size();
    int MinStamps = std::min(NumTankTimestamps,NumMRDTimestamps);


    std::map<uint64_t,uint64_t> BeamTankMRDPairs; //Pairs of beam-triggered Tank PMT/MRD counters ready to be built if all PMT waveforms are ready (TankAndMRD mode only)
    //*10 provides a buffer between paired events and events at the tail of the timestamp streams
    //Without it, events get moved to the orphanage that probably shouldn't
    if((MinStamps > (EventsPerPairing*10))||toolchain_stopping){
      this->RemoveCosmics();
      BeamTankMRDPairs = this->PairTankPMTAndMRDTriggers();
      this->ManageOrphanage();

      std::vector<uint64_t> BuiltTankTimes;
      for(std::pair<uint64_t,uint64_t> cpair : BeamTankMRDPairs){
        uint64_t TankCounterTime = cpair.first;
        if(verbosity>4) std::cout << "TANK EVENT WITH TIMESTAMP " << TankCounterTime << " HAS REQUIRED MINIMUM NUMBER OF WAVES TO BUILD" << std::endl;
        std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = FinishedTankEvents.at(TankCounterTime);
        uint64_t MRDTimeStamp = cpair.second;
        if(verbosity>4) std::cout << "MRD TIMESTAMP: " << MRDTimeStamp << std::endl;
        std::vector<std::pair<unsigned long,int>> MRDHits = myMRDMaps.MRDEvents.at(MRDTimeStamp);
        std::string MRDTriggerType = myMRDMaps.MRDTriggerTypeMap[MRDTimeStamp];
        int beam_tdc = myMRDMaps.MRDBeamLoopbackMap[MRDTimeStamp];
        int cosmic_tdc = myMRDMaps.MRDCosmicLoopbackMap[MRDTimeStamp];
        if(verbosity>4) std::cout << "BUILDING AN ANNIE EVENT" << std::endl;
        this->BuildANNIEEventRunInfo(CurrentRunNum, CurrentSubRunNum, CurrentRunType,CurrentStarTime);
        this->BuildANNIEEventTank(TankCounterTime, aWaveMap);
        this->BuildANNIEEventMRD(MRDHits, MRDTimeStamp, MRDTriggerType, beam_tdc, cosmic_tdc);
        this->SaveEntryToFile(CurrentRunNum,CurrentSubRunNum);
        if(verbosity>4) std::cout << "BUILT AN ANNIE EVENT SUCCESSFULLY" << std::endl;
        //Erase this entry from maps/vectors used when pairing completed events 
        BuiltTankTimes.push_back(TankCounterTime);
        FinishedTankEvents.erase(TankCounterTime);
        myMRDMaps.MRDEvents.erase(MRDTimeStamp);
        myMRDMaps.MRDTriggerTypeMap.erase(MRDTimeStamp);
        myMRDMaps.MRDBeamLoopbackMap.erase(MRDTimeStamp);
        myMRDMaps.MRDCosmicLoopbackMap.erase(MRDTimeStamp);
      }
      for(int i=0; i<(int)BuiltTankTimes.size(); i++){
        BeamTankMRDPairs.erase(BuiltTankTimes.at(i));
      }
    }
  }

  //Build ANNIE events based on matching Tank,MRD, and CTC timestamps
  else if (BuildType == "TankAndMRDAndCTC"){
    m_data->CStore.Get("NewTankPMTDataAvailable",IsNewTankData);
    if(IsNewTankData) this->ProcessNewTankPMTData();
    this->ManageOrphanage();

    
    //Look through our MRD data for any new timestamps
    m_data->CStore.Get("MRDEventTriggerTypes",myMRDMaps.MRDTriggerTypeMap);
    m_data->CStore.Get("MRDEvents",myMRDMaps.MRDEvents);
    m_data->CStore.Get("MRDBeamLoopback",myMRDMaps.MRDBeamLoopbackMap);
    m_data->CStore.Get("MRDCosmicLoopback",myMRDMaps.MRDCosmicLoopbackMap);
    m_data->CStore.Get("NewMRDDataAvailable",IsNewMRDData);
    if(IsNewMRDData) this->ProcessNewMRDData();

    //Look through CTC data for any new timestamps
    m_data->CStore.Get("TimeToTriggerWordMap",TimeToTriggerWordMap);
    m_data->CStore.Get("NewCTCDataAvailable",IsNewCTCData);
    if(IsNewCTCData) this->ProcessNewCTCData();

    //Now, pair up PMT/MRD/Triggers...
    int NumTankTimestamps = myTimeStream.BeamTankTimestamps.size();
    int NumMRDTimestamps = myTimeStream.BeamMRDTimestamps.size();
    int NumTrigs = myTimeStream.CTCTimestamps.size();
    
    // get get the most recent timestamp from each TimeStream
    uint64_t most_recent_beam = *std::max_element(myTimeStream.BeamTankTimestamps.begin(),
                                                  myTimeStream.BeamTankTimestamps.end());
    uint64_t most_recent_mrd = *std::max_element(myTimeStream.BeamMRDTimestamps.begin(),
                                                  myTimeStream.BeamMRDTimestamps.end());
    uint64_t most_recent_ctc = *std::max_element(myTimeStream.CTCTimestamps.begin(),
                                                  myTimeStream.CTCTimestamps.end());
    
    // find which TimeStream is lagging the most, and what time it's currently read up to.
    std::vector<uint64_t> newest_timestamps{most_recent_beam,most_recent_mrd,most_recent_ctc};
    uint64_t slowest_stream_timestamp = *std::min_element(newest_timestamps.begin(),newest_timestamps.end());
    
    if(verbosity>4){
      std::cout<<"ANNIEEventBuilder: most recent timestamps from each stream are: "
           <<most_recent_beam<<", "<<most_recent_mrd<<", "<<most_recent_ctc
           <<", with differences from the slowest stream being: "
           <<(static_cast<int64_t>(most_recent_beam)-static_cast<int64_t>(slowest_stream_timestamp))
           <<", "
           <<(static_cast<int64_t>(most_recent_mrd)-static_cast<int64_t>(slowest_stream_timestamp))
           <<", "
           <<(static_cast<int64_t>(most_recent_ctc)-static_cast<int64_t>(slowest_stream_timestamp))
           <<std::endl;
    }
    
    // always ensure the slowest stream is unpaused
           if(slowest_stream_timestamp==most_recent_beam){
      m_data->CStore.Set("PauseTankDecoding",false);
    } else if(slowest_stream_timestamp==most_recent_mrd){
      m_data->CStore.Set("PauseMRDDecoding",false);
    } else if(slowest_stream_timestamp==most_recent_ctc){
      m_data->CStore.Set("PauseCTCDecoding",false);
    }
    
    // if any other streams are more than X seconds ahead of the slowest one, pause them...
    if((static_cast<int64_t>(most_recent_beam)-static_cast<int64_t>(slowest_stream_timestamp))>pause_threshold){
      // secondly, only pause a stream once we have at least EventsPerPairing timestamps in it,
      // otherwise doing so will prevent building attempts
      if(NumTankTimestamps>EventsPerPairing){
        m_data->CStore.Set("PauseTankDecoding",true);
        Log("ANNIEEventBuilder: Pausing tank stream",v_debug,verbosity);
      }
    }
    if((static_cast<int64_t>(most_recent_mrd)-static_cast<int64_t>(slowest_stream_timestamp))>pause_threshold){
      if(NumMRDTimestamps>EventsPerPairing){
        m_data->CStore.Set("PauseMRDDecoding",true);
        Log("ANNIEEventBuilder: Pausing mrd stream",v_debug,verbosity);
      }
    }
    if((static_cast<int64_t>(most_recent_ctc)-static_cast<int64_t>(slowest_stream_timestamp))>pause_threshold){
      if(NumTrigs>EventsPerPairing){
        m_data->CStore.Set("PauseCTCDecoding",true);
        Log("ANNIEEventBuilder: Pausing ctc stream",v_debug,verbosity);
      }
    }
    
    if(verbosity>3){
        std::cout << "Number of CTCTimes, MRDTimes, PMTTimes: " << 
            NumTrigs << "," << NumMRDTimestamps << "," << NumTankTimestamps << std::endl;
    }

    std::map<uint64_t,std::map<std::string,uint64_t>> ThisBuildMap; //key: CTC timestamp, value: vector of maps.  Each map has the keys: "Tank", "MRD", "CTC", or "LAPPD" with values: 
                                                                    //timestamp of that stream (exception: "CTC" key has the trigger word as the value).
    
    int MinStamps = (NumTankTimestamps<NumMRDTimestamps) ? NumTankTimestamps : NumMRDTimestamps;
    MinStamps = (NumTrigs<MinStamps) ? NumTrigs : MinStamps;
    
    // only attempt timestamp matching if we have enough timestamps to try to match,
    // OR if the toolchain is being stopped (reached and of file, for example)
    if((MinStamps>EventsPerPairing)||toolchain_stopping){
      if(verbosity>4) std::cout << "MERGING COSMIC/MRD PAIRS " << std::endl;
      ThisBuildMap = this->PairCTCCosmicPairs(ThisBuildMap,std::max(most_recent_mrd,most_recent_ctc),toolchain_stopping);
      this->ManageOrphanage();

      if(verbosity>4) std::cout << "BEGINNING STREAM MERGING " << std::endl;
      ThisBuildMap = this->MergeStreams(ThisBuildMap,slowest_stream_timestamp,toolchain_stopping);
      Log("ANNIEEventBuilder: Calling ManageOrphanage post MergeStreams",v_debug,verbosity);
      this->ManageOrphanage();
      Log("ANNIEEventBuilder: Done managing orphanage",v_debug,verbosity);

      for(std::pair<uint64_t, std::map<std::string,uint64_t>> buildmap_entries : ThisBuildMap){
        uint64_t CTCtimestamp = buildmap_entries.first;
        std::map<std::string,uint64_t> aBuildSet = buildmap_entries.second; 
        for(std::pair<std::string,uint64_t> buildset_entries : aBuildSet){
          this->BuildANNIEEventRunInfo(CurrentRunNum, CurrentSubRunNum, CurrentRunType,CurrentStarTime);
          //If we have PMT and MRD infor for this trigger, build it 
          std::string label = buildset_entries.first;
          if(label == "CTC"){
            uint64_t CTCWord = buildset_entries.second;
            this->BuildANNIEEventCTC(CTCtimestamp,CTCWord);
            TimeToTriggerWordMap->erase(CTCtimestamp);
          }
          if(label == "TankPMT"){
            uint64_t TankPMTTime = buildset_entries.second;
            if(verbosity>4) std::cout << "TANK EVENT WITH TIMESTAMP " << TankPMTTime << "HAS REQUIRED MINIMUM NUMBER OF WAVES TO BUILD" << std::endl;
            std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = FinishedTankEvents.at(TankPMTTime);
            this->BuildANNIEEventTank(TankPMTTime, aWaveMap);
            FinishedTankEvents.erase(TankPMTTime);
          }
          if(label == "MRD"){
            uint64_t MRDTimeStamp = buildset_entries.second;
            if(verbosity>4) std::cout << "MRD TIMESTAMP: " << MRDTimeStamp << std::endl;
            std::vector<std::pair<unsigned long,int>> MRDHits = myMRDMaps.MRDEvents.at(MRDTimeStamp);
            std::string MRDTriggerType = myMRDMaps.MRDTriggerTypeMap[MRDTimeStamp];
            int beam_tdc = myMRDMaps.MRDBeamLoopbackMap[MRDTimeStamp];
            int cosmic_tdc = myMRDMaps.MRDCosmicLoopbackMap[MRDTimeStamp];
            this->BuildANNIEEventMRD(MRDHits, MRDTimeStamp, MRDTriggerType, beam_tdc, cosmic_tdc);
            myMRDMaps.MRDEvents.erase(MRDTimeStamp);
            myMRDMaps.MRDTriggerTypeMap.erase(MRDTimeStamp);
            myMRDMaps.MRDBeamLoopbackMap.erase(MRDTimeStamp);
            myMRDMaps.MRDCosmicLoopbackMap.erase(MRDTimeStamp);
          }
        }
        this->SaveEntryToFile(CurrentRunNum,CurrentSubRunNum);
        if(verbosity>4) std::cout << "BUILT AN ANNIE EVENT SUCCESSFULLY" << std::endl;
      }
    }
    ThisBuildMap.clear();
  }

  m_data->CStore.Set("MRDEvents",myMRDMaps.MRDEvents);
  m_data->CStore.Set("MRDEventTriggerTypes",myMRDMaps.MRDTriggerTypeMap);
  m_data->CStore.Set("MRDBeamLoopback",myMRDMaps.MRDBeamLoopbackMap);
  m_data->CStore.Set("MRDCosmicLoopback",myMRDMaps.MRDCosmicLoopbackMap);
  myMRDMaps.MRDEvents.clear();
  myMRDMaps.MRDTriggerTypeMap.clear();
  myMRDMaps.MRDBeamLoopbackMap.clear();
  myMRDMaps.MRDCosmicLoopbackMap.clear();

  ExecuteCount = 0;
  return true;
}


bool ANNIEEventBuilder::Finalise(){
  if(verbosity>4) std::cout << "ANNIEEvent Finalising.  Closing any open ANNIEEvent Boostore" << std::endl;
  if(verbosity>2) std::cout << "ANNIEEventBuilder: Saving and closing file." << std::endl;
  ANNIEEvent->Close();
  ANNIEEvent->Delete();
  delete ANNIEEvent;
  if(verbosity>2) std::cout << "PMT/MRD Orphan number at finalise: " << myOrphanage.OrphanTankTimestamps.size() <<
          "," << myOrphanage.OrphanMRDTimestamps.size() << std::endl;
  //Save the current subrun and delete ANNIEEvent
  OrphanStore->Close();
  OrphanStore->Delete();
  delete OrphanStore;
  std::cout << "ANNIEEventBuilder Exitting" << std::endl;
  return true;
}

void ANNIEEventBuilder::ProcessNewCTCData(){
  for(std::pair<uint64_t,uint32_t> apair : *TimeToTriggerWordMap){
    uint64_t CTCTimeStamp = apair.first;
    uint32_t CTCWord = apair.second;
    myTimeStream.CTCTimestamps.push_back(CTCTimeStamp);
    if(verbosity>5)std::cout << "CTCTIMESTAMP,WORD" << CTCTimeStamp << "," << CTCWord << std::endl;
  }
  RemoveDuplicates(myTimeStream.CTCTimestamps);
  return;
}

void ANNIEEventBuilder::ProcessNewMRDData(){
  for(std::pair<uint64_t,std::vector<std::pair<unsigned long,int>>> apair : myMRDMaps.MRDEvents){
    uint64_t MRDTimeStamp = apair.first;
    myTimeStream.BeamMRDTimestamps.push_back(MRDTimeStamp);
    if(verbosity>5)std::cout << "MRDTIMESTAMPTRIGTYPE," << MRDTimeStamp << "," << myMRDMaps.MRDTriggerTypeMap.at(MRDTimeStamp) << std::endl;
  }
  RemoveDuplicates(myTimeStream.BeamMRDTimestamps);
  return;
}

void ANNIEEventBuilder::ProcessNewTankPMTData(){
  //Check if any In-progress tank events now have all waveforms
  m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);
  std::vector<uint64_t> InProgressTankEventsToDelete;
  std::map<uint64_t,std::string> TankOrphans;
  std::map<uint64_t,std::string> MRDOrphans;
  std::map<uint64_t,std::string> CTCOrphans;
  
  if(verbosity>5) std::cout << "ANNIEEventBuilder Tool: Processing new tank data " << std::endl;
  for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : *InProgressTankEvents){
    uint64_t PMTCounterTimeNs = apair.first;
    std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
    if(verbosity>4) std::cout << "Number of waves for this counter: " << aWaveMap.size() << std::endl;
    
    //Push back any new timestamps, then remove duplicates in the end
    if(PMTCounterTimeNs>NewestTankTimestamp){
      NewestTankTimestamp = PMTCounterTimeNs;
      if(verbosity>3)std::cout << "TANKTIMESTAMP," << PMTCounterTimeNs << std::endl;
    }
    //Events and delete it from the in-progress events
    int NumTankPMTChannels = TankPMTCrateSpaceToChannelNumMap.size();
    int NumAuxChannels = AuxCrateSpaceToChannelNumMap.size();
    if(aWaveMap.size() >= (NumWavesInCompleteSet)){
      FinishedTankEvents.emplace(PMTCounterTimeNs,aWaveMap);
      myTimeStream.BeamTankTimestamps.push_back(PMTCounterTimeNs);
      //Put PMT timestamp into the timestamp set for this run.
      if(verbosity>4) std::cout << "Finished waveset has clock counter: " << PMTCounterTimeNs << std::endl;
      InProgressTankEventsToDelete.push_back(PMTCounterTimeNs);
    }

    //If this InProgressTankEvent is too old, clear it
    //out from all TankTimestamp maps
    if(OrphanOldTankTimestamps && ((NewestTankTimestamp - PMTCounterTimeNs) > OldTimestampThreshold*1E9)){
      InProgressTankEventsToDelete.push_back(PMTCounterTimeNs);
      TankOrphans.emplace(PMTCounterTimeNs,"incomplete_tank_event");
    }
  }
  RemoveDuplicates(myTimeStream.BeamTankTimestamps);
  
  // move abandoned in-progress events to the orphanage
  this->MoveToOrphanage(TankOrphans, MRDOrphans, CTCOrphans);
 
  //Since timestamp pairing has been done for finished Tank Events,
  //Erase the finished Tank Events from the InProgressTankEventsMap
  if(verbosity>3)std::cout << "Now erasing finished timestamps from InProgressTankEvents" << std::endl;
  for (unsigned int j=0; j< InProgressTankEventsToDelete.size(); j++){
    InProgressTankEvents->erase(InProgressTankEventsToDelete.at(j));
  }
  if(verbosity>3) std::cout << "Current number of unfinished PMT Waveforms: " << InProgressTankEvents->size() << std::endl;

  return;
}

std::map<uint64_t,std::map<std::string,uint64_t>> ANNIEEventBuilder::PairCTCCosmicPairs(std::map<uint64_t,std::map<std::string,uint64_t>> BuildMap, uint64_t max_timestamp, bool force_matching){
  uint32_t CosmicWord = 36;  //TS_in(35) + 1, where 35 is the MRD_CR_Trigger
  std::vector<uint64_t> MRDCosmicTimes;
  std::vector<uint64_t> MRDStampsToDelete;
  std::map<uint64_t,uint64_t> PairedCTCMRDTimes;
  std::map<uint64_t,std::string> MRDOrphans;
  std::map<uint64_t,std::string> TankOrphans;
  std::map<uint64_t,std::string> CTCOrphans;

  for(auto&& aTS : myTimeStream.BeamMRDTimestamps) {
    if((aTS>max_timestamp)&&(!force_matching)) break; // don't try to match yet
    std::string MRDTriggerType = myMRDMaps.MRDTriggerTypeMap[aTS];
    if(MRDTriggerType == "Cosmic"){
      MRDCosmicTimes.push_back(aTS);
    }
  }
 
  //Now, loop through CTC timestamps and pair any cosmic triggers
  //With Cosmic candidates in time tolerance

  //First, pair up myTimeStream.CTCTimestamps and MRD/PMT timestamps
  std::map<uint64_t,uint64_t> PairedCTCCosmicTimes;

  //Form CTC-MRD pairs
  if(verbosity>3) std::cout << "Finding CTC-MRD pairs..." << std::endl;
  for(auto&& aMrdTS : MRDCosmicTimes){
    for(auto&& aCtcTS : myTimeStream.CTCTimestamps){
      if((aCtcTS>max_timestamp)&&(!force_matching)) continue;        // do not try to match just yet
      if(TimeToTriggerWordMap->at(aCtcTS) != CosmicWord) continue;   // not a cosmic CTC
      double TSDiff =  static_cast<double>(aMrdTS) - 
                       static_cast<double>(aCtcTS);
      if(verbosity>4) std::cout << "CosmicTS - CTCTS in nanoseconds is " << TSDiff << std::endl;
      if(TSDiff>CTCMRDTimeTolerance){ // MRD timestamp is later than this CTC one. Check next CTC timestamp.
        continue;
      } else if (TSDiff<(-1.*static_cast<double>(CTCMRDTimeTolerance))){ // We've crossed past where a pair would be found
        if(verbosity>4) std::cout << "NO CTC STAMP FOUND MATCHING COSMIC STAMP... MRD TO ORPHANAGE" << std::endl;
        MRDOrphans.emplace(aMrdTS,"mrd_cosmic_no_ctc");
        break;
      } else { //We've found a valid CTC-MRD pair!
        if(verbosity>4) std::cout << "FOUND A MATCHING CTC TIME FOR THIS COSMIC MRD TIMESTAMP. NICE, PAIR EM." << std::endl;
        PairedCTCMRDTimes.emplace(aCtcTS, aMrdTS);
        break;
      }
    }
  }

  //Neat.  Now, add timestamps to the buildmap.
  std::vector<uint64_t> BuiltCTCs;
  for(std::pair<uint64_t, uint64_t> aPair : PairedCTCMRDTimes){
    uint64_t CTCTimestamp = aPair.first;
    uint64_t CosmicTimestamp = aPair.second;
    std::map<std::string,uint64_t> aBuildSet;
    aBuildSet.emplace("CTC",TimeToTriggerWordMap->at(CTCTimestamp));
    aBuildSet.emplace("MRD",CosmicTimestamp);
    myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
         myTimeStream.BeamMRDTimestamps.end(),CosmicTimestamp), 
         myTimeStream.BeamMRDTimestamps.end());
    if(verbosity>4) std::cout << "BUILDING A CTC/COSMIC BUILD MAP ENTRY. CTC IS " << CTCTimestamp << std::endl;
    BuildMap.emplace(CTCTimestamp,aBuildSet);
    BuiltCTCs.push_back(CTCTimestamp);
  }
  

  //Delete CTCTimestamps that have a PMT or MRD pair from CTC timestamp tracker
  for(int j=0;j<(int)BuiltCTCs.size();j++){
    if(verbosity>4) std::cout << "REMOVING CTC TIME OF " << BuiltCTCs.at(j) << "FROM CTCTIMESTAMPS VECTOR" << std::endl;
    myTimeStream.CTCTimestamps.erase(std::remove(myTimeStream.CTCTimestamps.begin(),
        myTimeStream.CTCTimestamps.end(),BuiltCTCs.at(j)), 
        myTimeStream.CTCTimestamps.end());
  }

  //Move MRD timestamps with no pairs to the orphanage.  Just empty Tank and CTC vectors for 
  //Input to function.  TODO; could overload function
  this->MoveToOrphanage(TankOrphans, MRDOrphans, CTCOrphans);
  return BuildMap;
}

void ANNIEEventBuilder::RemoveCosmics(){
  std::vector<uint64_t> MRDStampsToDelete;
  for (int i=0;i<(int)(myTimeStream.BeamMRDTimestamps.size()); i++) {
    std::string MRDTriggerType = myMRDMaps.MRDTriggerTypeMap[myTimeStream.BeamMRDTimestamps.at(i)];
    if(verbosity>4) std::cout << "THIS MRD TRIGGER TYPE IS: " << MRDTriggerType << std::endl;
    if(MRDTriggerType == "Cosmic"){
      MRDStampsToDelete.push_back(myTimeStream.BeamMRDTimestamps.at(i));
      if(verbosity>3){
        std::cout << "NO BEAM COSMIC COSMIC IS: " << myTimeStream.BeamMRDTimestamps.at(i) << std::endl;
        std::cout << "DELETING FROM TIMESTAMPS TO PAIR/BUILD" << std::endl;
      }
    }
  }
  for (int j=0; j<(int)MRDStampsToDelete.size(); j++){
    myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
        myTimeStream.BeamMRDTimestamps.end(),MRDStampsToDelete.at(j)), 
        myTimeStream.BeamMRDTimestamps.end());
    myMRDMaps.MRDEvents.erase(MRDStampsToDelete.at(j));
    myMRDMaps.MRDTriggerTypeMap.erase(MRDStampsToDelete.at(j));
    myMRDMaps.MRDBeamLoopbackMap.erase(MRDStampsToDelete.at(j));
    myMRDMaps.MRDCosmicLoopbackMap.erase(MRDStampsToDelete.at(j));
  }
  return;
}


std::map<uint64_t,std::map<std::string,uint64_t>> ANNIEEventBuilder::MergeStreams(std::map<uint64_t,std::map<std::string,uint64_t>> BuildMap, uint64_t max_timestamp, bool force_matching){
  //This method takes timestamps from the BeamMRDTimestamps, BeamTankTimestamps, and
  //CTCTimestamps vectors (acquired as the building continues) and builds maps 
  //stored in the BuildMap and used to build ANNIEEvents.
 
  //First, pair up CTCTimestamps and MRD/PMT timestamps
  std::map<uint64_t,uint64_t> PairedCTCTankTimes; 
  std::map<uint64_t,uint64_t> PairedCTCMRDTimes;
  std::map<uint64_t,std::string> MRDOrphans;
  std::map<uint64_t,std::string> TankOrphans;
  std::map<uint64_t,std::string> CTCOrphans;

  // only attempt matching of any kind on timestamps older than the newest timestamp we have
  // from ALL streams - i.e. if the slowest stream has only read up to 4pm, do not try to do
  // matching on any timestamps newer than this
  //Form CTC-MRD pairs
  uint64_t max_ctc=0;
  if(verbosity>3) std::cout << "Finding CTC-MRD pairs..." << std::endl;
  for(auto&& aMrdTS : myTimeStream.BeamMRDTimestamps){
    if((aMrdTS>max_timestamp)&&(!force_matching)) continue;   // don't try to match just yet
    for(auto&& aCtcTS : myTimeStream.CTCTimestamps){
      if((aCtcTS>max_timestamp)&&(!force_matching)) continue; // don't try to match this just yet
      double TSDiff =  static_cast<double>(aMrdTS) - 
                       static_cast<double>(aCtcTS);
      if(verbosity>4) std::cout << "MRDTS - CTCTS in nanoseconds is " << TSDiff << std::endl;
      if(TSDiff>CTCMRDTimeTolerance){ // MRD timestamp is later than CTC timestamp; try next CTC timestamp
        continue;
      } else if (TSDiff<(-1.*static_cast<double>(CTCMRDTimeTolerance))){ //We've crossed past where a pair would be found
        if(verbosity>4) std::cout << "NO CTC STAMP FOUND MATCHING MRD STAMP... MRD TO ORPHANAGE" << std::endl;
        MRDOrphans.emplace(aMrdTS,"mrd_beam_no_ctc");
        break;
      } else { //We've found a valid CTC-MRD pair!
        if(verbosity>4) std::cout << "FOUND A MATCHING CTC TIME FOR THIS MRD TIMESTAMP. NICE, PAIR EM." << std::endl;
        PairedCTCMRDTimes.emplace(aCtcTS, aMrdTS);
        max_ctc=aCtcTS;
        break;
      }
    }
  }

  //Now form CTC-PMT pairs
  if(verbosity>3) std::cout << "Finding CTC-Tank pairs..." << std::endl;
  for(auto&& aTankTS : myTimeStream.BeamTankTimestamps){
    if((aTankTS>max_timestamp)&&(!force_matching)) continue; // don't try to match just yet
    for(auto&& aCtcTS : myTimeStream.CTCTimestamps){
      if((aCtcTS>max_timestamp)&&(!force_matching)) continue; // don't try to match just yet
      double TSDiff =  static_cast<double>(aTankTS) - 
                       static_cast<double>(aCtcTS);
      if(verbosity>4) std::cout << "TankTS - CTCTS in nanoseconds is " << TSDiff << std::endl;
      if(TSDiff>CTCTankTimeTolerance){ // Need to move forward in indices
        continue;
      } else if (TSDiff<(-1.*static_cast<double>(CTCTankTimeTolerance))){ //We've crossed past where a pair would be found
        if(verbosity>4) std::cout << "NO CTC STAMP FOUND MATCHING TANK STAMP... TANK TO ORPHANAGE" << std::endl;
        TankOrphans.emplace(aTankTS,"tank_no_ctc");
        break;
      } else { //We've found a valid CTC-Tank pair!
        if(verbosity>4) std::cout << "FOUND A MATCHING CTC TIME FOR THIS TANK TIMESTAMP. NICE, PAIR EM." << std::endl;
        PairedCTCTankTimes.emplace(aCtcTS, aTankTS);
        if(aCtcTS>max_ctc) max_ctc=aCtcTS;
        break;
      }
    }
  }
  int LargestCTCIndex = std::distance(myTimeStream.CTCTimestamps.begin(),std::find(myTimeStream.CTCTimestamps.begin(),myTimeStream.CTCTimestamps.end(),max_ctc));
  if(max_ctc==0) LargestCTCIndex=0;
  if(verbosity>4) std::cout << "LARGEST CTC INDEX WAS " << LargestCTCIndex << std::endl;
  if(verbosity>4) std::cout << "AND CTCTIMESTAMPS SIZE IS " << myTimeStream.CTCTimestamps.size() << std::endl;
  
  //Neat.  Now, add timestamps to the buildmap.
  //FIXME: Right now, this will only build PMT/MRD/CTC or PMT/CTC pairs.  
  // Any logic needed for pairing CTC/MRD pairs?  I don't think so, since cosmics 
  // handled in a different method..
  std::vector<uint64_t> BuiltCTCs;
  for(int i=0; i<(LargestCTCIndex-1); i++){
    uint64_t CTCKey = myTimeStream.CTCTimestamps.at(i);
    if(verbosity>4) std::cout << "TRYING TO BUILD A SET WITH CTCTIMESTAMP INDEX " << i << std::endl;

    bool have_tankmatch = false;
    bool have_mrdmatch = false;    
    std::map<uint64_t, uint64_t>::iterator it_tank = PairedCTCTankTimes.find(CTCKey);
    if(it_tank != PairedCTCTankTimes.end()) have_tankmatch = true;
    std::map<uint64_t, uint64_t>::iterator it_mrd = PairedCTCMRDTimes.find(CTCKey);
    if(it_mrd != PairedCTCMRDTimes.end()) have_mrdmatch = true;
    if(have_tankmatch && have_mrdmatch){
      std::map<std::string,uint64_t> aBuildSet;
      aBuildSet.emplace("CTC",TimeToTriggerWordMap->at(CTCKey));
      aBuildSet.emplace("TankPMT",it_tank->second);
      aBuildSet.emplace("MRD",it_mrd->second);
      myTimeStream.BeamTankTimestamps.erase(std::remove(myTimeStream.BeamTankTimestamps.begin(),
           myTimeStream.BeamTankTimestamps.end(),it_tank->second), 
           myTimeStream.BeamTankTimestamps.end());
      myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
           myTimeStream.BeamMRDTimestamps.end(),it_mrd->second), 
           myTimeStream.BeamMRDTimestamps.end());
      if(verbosity>4) std::cout << "BUILDING A BUILD MAP ENTRY. CTC IS " << CTCKey << std::endl;
      BuildMap.emplace(CTCKey,aBuildSet);
      BuiltCTCs.push_back(CTCKey);
    } else if (have_tankmatch && !have_mrdmatch){
      std::map<std::string,uint64_t> aBuildSet;
      aBuildSet.emplace("CTC",TimeToTriggerWordMap->at(CTCKey));
      aBuildSet.emplace("TankPMT",it_tank->second);
      myTimeStream.BeamTankTimestamps.erase(std::remove(myTimeStream.BeamTankTimestamps.begin(),
           myTimeStream.BeamTankTimestamps.end(),it_tank->second), 
           myTimeStream.BeamTankTimestamps.end());
      if(verbosity>4) std::cout << "BUILDING A PMT ONLY BUILD MAP ENTRY. CTC IS " << CTCKey << std::endl;
      BuildMap.emplace(CTCKey,aBuildSet);
      BuiltCTCs.push_back(CTCKey);
    } else if (!have_tankmatch && !have_mrdmatch){
      if(verbosity>4) std::cout << "NO MRD OR TANK TIMESTAMP FOR THIS CTC TIME... ORPHAN THE CTC" << std::endl;
      //uint64_t latest_tank_orphan = TankOrphans.at(TankOrphans.size()-1);
      //uint64_t latest_MRD_orphan = MRDOrphans.at(MRDOrphans.size()-1);
      //if(CTCKey < latest_tank_orphan && CTCKey < latest_MRD_orphan){
      CTCOrphans.emplace(CTCKey,"ctc_no_mrd_or_tank");
      //}
    } else {
      // beam ctc and mrd match, but no tank
      // at the least do we need to orphan the CTC and MRD events?
      CTCOrphans.emplace(CTCKey,"ctc_mrd_no_tank");
      MRDOrphans.emplace(it_mrd->second,"ctc_mrd_no_tank");
      if(verbosity>4) std::cout << "CTC TIMESTAMP " << CTCKey << " PAIRS WITH EITHER A PMT OR MRD..."  << std::endl;
    }
  }

  //Delete myTimeStream.CTCTimestamps that have a PMT or MRD pair from CTC timestamp tracker
  for(int j=0;j<(int)BuiltCTCs.size();j++){
    if(verbosity>4) std::cout << "REMOVING CTC TIME OF " << BuiltCTCs.at(j) << "FROM CTCTIMESTAMPS VECTOR" << std::endl;
    myTimeStream.CTCTimestamps.erase(std::remove(myTimeStream.CTCTimestamps.begin(),
        myTimeStream.CTCTimestamps.end(),BuiltCTCs.at(j)), 
        myTimeStream.CTCTimestamps.end());
  }

  //Move timestamps with no pairs to the orphanage
  this->MoveToOrphanage(TankOrphans, MRDOrphans, CTCOrphans);
  Log("ANNIEEventBuilder: Returning from Merging the Streams",v_debug,verbosity);

  return BuildMap;
}

void ANNIEEventBuilder::MoveToOrphanage(std::map<uint64_t, std::string> TankOrphans,
                                        std::map<uint64_t, std::string> MRDOrphans,
                                        std::map<uint64_t, std::string> CTCOrphans){
  if(verbosity>3) std::cout << "MOVING TIMESTAMPS WITH NO FAMILY TO ORPHANAGE" << std::endl;
  //Finally, we need to move data associated with our orphaned timestamps to the orphange
  for(auto&& nextorphan : CTCOrphans){
    uint64_t CTCOrphanStamp = nextorphan.first;
    
    // build map of orphan info to save to OrphanStore
    std::map<std::string,std::string> orphaninfo;
    orphaninfo.emplace("reason",nextorphan.second);
    orphaninfo.emplace("ctc_word",to_string(TimeToTriggerWordMap->at(CTCOrphanStamp)));
    
    // move to orphanage
    myOrphanage.OrphanCTCTimestamps.emplace(CTCOrphanStamp,orphaninfo);
    
    // remove the orphan from the timestream
    myTimeStream.CTCTimestamps.erase(std::remove(myTimeStream.CTCTimestamps.begin(),
         myTimeStream.CTCTimestamps.end(),CTCOrphanStamp), 
         myTimeStream.CTCTimestamps.end());
  }
  
  for(auto&& nextorphan : TankOrphans){
    uint64_t TankOrphanStamp = nextorphan.first;
    
    // build map of orphan info to save to OrphanStore
    std::map<std::string,std::string> orphaninfo;
    orphaninfo.emplace("reason",nextorphan.second);
    
    // move to orphanage
    myOrphanage.OrphanTankTimestamps.emplace(TankOrphanStamp,orphaninfo);
    
    // remove the orphan from the timestream
    myTimeStream.BeamTankTimestamps.erase(std::remove(myTimeStream.BeamTankTimestamps.begin(),
               myTimeStream.BeamTankTimestamps.end(),TankOrphanStamp), 
               myTimeStream.BeamTankTimestamps.end());
  }
  
  for(auto&& nextorphan : MRDOrphans){
    uint64_t MrdOrphanStamp = nextorphan.first;
    
    // build map of orphan info to save to OrphanStore
    std::map<std::string,std::string> orphaninfo;
    orphaninfo.emplace("reason",nextorphan.second);
    
    // move to orphanage
    myOrphanage.OrphanMRDTimestamps.emplace(MrdOrphanStamp,orphaninfo);
    
    // remove the orphan from the timestream
    myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
               myTimeStream.BeamMRDTimestamps.end(),MrdOrphanStamp), 
               myTimeStream.BeamMRDTimestamps.end());
  }
  if(verbosity>3) std::cout << "ORPHAN MOVEMENT COMPLETE" << std::endl;
  return;
}

void ANNIEEventBuilder::ManageOrphanage(){
  
  // TODO
  // try to find incomplete tank orphans to match ctc/mrd orphans
  // build a fixed type that says what data was present:
  // MRD? CTC? Tank? Right now without merging it'll be just one,
  // but later if we're able to merge orphans it could be more than one
//  std::cout<<"Going to manage the orphanage"<<std::endl;
  
  std::string OrphanFile = SavePath + OrphanFileBase + "R" + to_string(CurrentRunNum) + 
      "S" + to_string(CurrentSubRunNum);
//  std::cout<<"orphans will be saved to "<<OrphanFile<<std::endl;
//  std::cout<<"OrphanStore is "<<OrphanStore<<std::endl;
  
//  std::cout<<"managing tank orphans"<<std::endl;
  //For now, just save orphans then delete them
  for(auto&& nextorphan : myOrphanage.OrphanTankTimestamps){
    // copy the orphan info into OrphanStore
//    std::cout<<"set tank orphan info into OrphanStore"<<std::endl;
    OrphanStore->Set("EventType",std::string("Tank"));
    OrphanStore->Set("Timestamp",nextorphan.first);
    OrphanStore->Set("Reason",nextorphan.second.at("reason"));
    OrphanStore->Set("Info",nextorphan.second); // redundant for now, but maybe we'll add more info later
//    std::cout<<"Variables set, calling Save to write Orphans to file "<<OrphanFile<<std::endl;
    OrphanStore->Save(OrphanFile);
//    std::cout<<"Orphans have been saved, deleting the orphanage"<<std::endl;
    OrphanStore->Delete();
//    std::cout<<"orphanage has been deleted, erasing the orphan from the FinishedTankEvents"<<std::endl;
    
    // cleanup from events to process
    FinishedTankEvents.erase(nextorphan.first);
//    std::cout<<"done, move to next orphan"<<std::endl;
  }
//  std::cout<<"managing mrd orphans"<<std::endl;
  for(auto&& nextorphan : myOrphanage.OrphanMRDTimestamps){
    // copy the orphan info into OrphanStore
    OrphanStore->Set("EventType",std::string("MRD"));
    OrphanStore->Set("Timestamp",nextorphan.first);
    OrphanStore->Set("Reason",nextorphan.second.at("reason"));
    OrphanStore->Set("Info",nextorphan.second); // redundant for now, but maybe we'll add more info later
    OrphanStore->Save(OrphanFile);
    OrphanStore->Delete();
    
    // cleanup from events to process
    myMRDMaps.MRDEvents.erase(nextorphan.first);
    myMRDMaps.MRDTriggerTypeMap.erase(nextorphan.first);
    myMRDMaps.MRDBeamLoopbackMap.erase(nextorphan.first);
    myMRDMaps.MRDCosmicLoopbackMap.erase(nextorphan.first);
  }
//  std::cout<<"Managing CTC orphans"<<std::endl;
  for(auto&& nextorphan : myOrphanage.OrphanCTCTimestamps){
    // copy the orphan info into OrphanStore
    OrphanStore->Set("EventType",std::string("CTC"));
    OrphanStore->Set("Timestamp",nextorphan.first);
    OrphanStore->Set("Reason",nextorphan.second.at("reason"));
    OrphanStore->Set("Info",nextorphan.second); // CTC word
    OrphanStore->Save(OrphanFile);
    OrphanStore->Delete();
    
    // cleanup from events to process
    TimeToTriggerWordMap->erase(nextorphan.first);
  }
//  std::cout<<"all CTC orphans handled, clearing the orphan timestamps"<<std::endl;
  
  myOrphanage.OrphanTankTimestamps.clear();
  myOrphanage.OrphanMRDTimestamps.clear();
  myOrphanage.OrphanCTCTimestamps.clear();
//  std::cout<<"timestamps cleared, returning"<<std::endl;
  return;
}

std::map<uint64_t,uint64_t> ANNIEEventBuilder::PairTankPMTAndMRDTriggers(){
  if(verbosity>4) std::cout << "ANNIEEventBuilder Tool: Beginning to pair events" << std::endl;

  std::map<uint64_t,uint64_t> TankMRDTimePairs; //Pairs of beam-triggered Tank PMT/MRD counters ready to be built if all PMT waveforms are ready (TankAndMRD mode only)
  std::vector<uint64_t> TankStampsToDelete;
  std::vector<uint64_t> MRDStampsToDelete;
  int NumPairsToMake = EventsPerPairing;
  std::vector<double> ThisPairingTSDiffs;
  
  std::map<uint64_t,std::string> MRDOrphans;
  std::map<uint64_t,std::string> TankOrphans;
  std::map<uint64_t,std::string> CTCOrphans;

  int NumOrphans = 0;

  //Pair PMT and MRD timestamps and calculate how much PMTTime - MRDTime has drifted 
  if(verbosity>4) std::cout << "MEAN OF PMT-MRD TIME DIFFERENCE LAST LOOP: " << CurrentDriftMean << std::endl;
  for (int i=0;i<NumPairsToMake; i++) {
    double TSDiff = static_cast<double>(myTimeStream.BeamTankTimestamps.at(i)) - 
                    static_cast<double>(myTimeStream.BeamMRDTimestamps.at(i));
    if(verbosity>4){
      std::cout << "PAIRED TANK TIMESTAMP: " << myTimeStream.BeamTankTimestamps.at(i) << std::endl;
      std::cout << "PAIRED MRD TIMESTAMP: " << myTimeStream.BeamMRDTimestamps.at(i) << std::endl;
      std::cout << "DIFFERENCE BETWEEN PMT AND MRD TIMESTAMP (ns): " << 
      (static_cast<double>(myTimeStream.BeamTankTimestamps.at(i)) - static_cast<double>(myTimeStream.BeamMRDTimestamps.at(i))) << std::endl;
    }
    if(std::abs(TSDiff-CurrentDriftMean) > MRDTankTimeTolerance){ // PMT/MRD timestamps farther apart than set tolerance
      if(verbosity>3) std::cout << "DEVIATION OF " << MRDTankTimeTolerance << " ms DETECTED IN STREAMS" << std::endl;
      if(TSDiff > 0) {
        if(verbosity>3) std::cout << "MOVING MRD TIMESTAMP TO ORPHANAGE" << std::endl;
        MRDOrphans.emplace(myTimeStream.BeamMRDTimestamps.at(i),"mrd_beam_no_tank");
        NumOrphans +=1;
        myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
            myTimeStream.BeamMRDTimestamps.end(),myTimeStream.BeamMRDTimestamps.at(i)), 
            myTimeStream.BeamMRDTimestamps.end());
      } else {
        if(verbosity>3) std::cout << "MOVING TANK TIMESTAMP TO ORPHANAGE" << std::endl;
        TankOrphans.emplace(myTimeStream.BeamTankTimestamps.at(i),"tank_no_mrd");
        NumOrphans +=1;
        myTimeStream.BeamTankTimestamps.erase(std::remove(myTimeStream.BeamTankTimestamps.begin(),
            myTimeStream.BeamTankTimestamps.end(),myTimeStream.BeamTankTimestamps.at(i)), 
            myTimeStream.BeamTankTimestamps.end());
      }
      i-=1;
      NumPairsToMake-=1;
    } else {  //PMT/MRD timestamp times within tolerance; Leave stamps in stream to be paired
      ThisPairingTSDiffs.push_back(TSDiff);
    }
  }
  this->MoveToOrphanage(TankOrphans, MRDOrphans, CTCOrphans);

  //With the last set of pairs calculate what the mean drift is
  //TODO: Error handling for high drift and high orphan rates?
  double ThisPairingMean, ThisPairingVariance;
  if(ThisPairingTSDiffs.size()>2){
    ComputeMeanAndVariance(ThisPairingTSDiffs,ThisPairingMean,ThisPairingVariance);
    if((std::abs(CurrentDriftMean-ThisPairingMean)>DriftWarningValue) && (verbosity>=v_warning)){
      std::cout << "ANNIEEventBuilder tool: WARNING! Shift in drift greater than " << DriftWarningValue << " since last pairings." << std::endl;
    }
    if((NumOrphans > OrphanWarningValue) && (verbosity>=v_warning)){
      std::cout << "ANNIEEventBuilder tool: WARNING! High orphan rate detected.  More than " << OrphanWarningValue << " this pairing sequence." << std::endl;
    }
    CurrentDriftMean = ThisPairingMean;
    CurrentDriftVariance = ThisPairingVariance;
  }
  for (int i=0;i<(NumPairsToMake); i++) {
    TankMRDTimePairs.emplace(myTimeStream.BeamTankTimestamps.at(i),myTimeStream.BeamMRDTimestamps.at(i));
    TankStampsToDelete.push_back(myTimeStream.BeamTankTimestamps.at(i));
    MRDStampsToDelete.push_back(myTimeStream.BeamMRDTimestamps.at(i));
  }

  if(verbosity>4) std::cout << "DELETE PAIRED TIMESTAMPS FROM THE TIMESTAMP STREAMS " << std::endl;
  //Delete all paired timestamps still in the unpaired unbuilt timestamps vectors 
  for (int i=0; i<(NumPairsToMake); i++){
    myTimeStream.BeamTankTimestamps.erase(std::remove(myTimeStream.BeamTankTimestamps.begin(),
        myTimeStream.BeamTankTimestamps.end(),TankStampsToDelete.at(i)), 
        myTimeStream.BeamTankTimestamps.end());
    myTimeStream.BeamMRDTimestamps.erase(std::remove(myTimeStream.BeamMRDTimestamps.begin(),
               myTimeStream.BeamMRDTimestamps.end(),MRDStampsToDelete.at(i)), 
               myTimeStream.BeamMRDTimestamps.end());
  }
  return TankMRDTimePairs;
}


void ANNIEEventBuilder::BuildANNIEEventMRD(std::vector<std::pair<unsigned long,int>> MRDHits, 
        uint64_t MRDTimeStamp, std::string MRDTriggerType, int beam_tdc, int cosmic_tdc)
{
  std::cout << "Building an ANNIE Event (MRD), ANNIEEventNum = "<<ANNIEEventNum << std::endl;
  TDCData = new std::map<unsigned long, std::vector<Hit>>;

  //TODO: Loop through MRDHits at this timestamp and form the Hit vector.
  for (unsigned int i_value=0; i_value< MRDHits.size(); i_value++){
    unsigned long channelkey = MRDHits.at(i_value).first;
    int hitTimeADC = MRDHits.at(i_value).second;
    double hitTime = 4000. - 4.*(double)hitTimeADC;
    if (verbosity > 4) std::cout <<"creating hit with ADC value "<<hitTimeADC<<"and chankey "<<channelkey<<std::endl;
    if (verbosity > 3) std::cout <<"creating hit with time vlaue "<<hitTime<<"and chankey "<<channelkey<<std::endl;
    if (TDCData->count(channelkey)==0){
      std::vector<Hit> newhitvector;
      newhitvector.push_back(Hit(0,hitTime,1.));    //Hit(tubeid, time, charge). 1 TDC tick corresponds to 4ns, no charge information (set to 1)
      TDCData->emplace(channelkey,newhitvector);
    } else {
      TDCData->at(channelkey).push_back(Hit(0,hitTime,1.));
    }
  }

  std::map<std::string, int> mrd_loopback_tdc;
  mrd_loopback_tdc.emplace("BeamLoopbackTDC",beam_tdc);
  mrd_loopback_tdc.emplace("CosmicLoopbackTDC",cosmic_tdc);

  Log("ANNIEEventBuilder: TDCData size: "+std::to_string(TDCData->size()),v_debug,verbosity);

  ANNIEEvent->Set("TDCData",TDCData,true);
  TimeClass timeclass_timestamp(MRDTimeStamp);
  ANNIEEvent->Set("EventTimeMRD",timeclass_timestamp);
  ANNIEEvent->Set("MRDTriggerType",MRDTriggerType);
  ANNIEEvent->Set("MRDLoopbackTDC",mrd_loopback_tdc);
  return;
}

void ANNIEEventBuilder::BuildANNIEEventRunInfo(int RunNumber, int SubRunNumber,
        int RunType, uint64_t StartTime)
{
  if(verbosity>v_message)std::cout << "Building an ANNIE Event Run Info" << std::endl;
  ANNIEEvent->Set("EventNumber",ANNIEEventNum);
  ANNIEEvent->Set("RunNumber",RunNumber);
  ANNIEEvent->Set("SubrunNumber",SubRunNumber);
  ANNIEEvent->Set("RunType",RunType);
  ANNIEEvent->Set("RunStartTime",StartTime);
  //TODO: Things missing from ANNIEEvent that should be in before this tool finishes:
  //  - BeamStatus?  
  return;
}

void ANNIEEventBuilder::BuildANNIEEventCTC(uint64_t CTCTime, uint32_t CTCWord)
{
  if(verbosity>v_message)std::cout << "Building an ANNIE Event CTC Info" << std::endl;
  ANNIEEvent->Set("CTCTimestamp",CTCTime);
  ANNIEEvent->Set("TriggerWord",CTCWord);
  return;
}

void ANNIEEventBuilder::BuildANNIEEventTank(uint64_t ClockTime, 
        std::map<std::vector<int>, std::vector<uint16_t>> WaveMap)
{
  if(verbosity>v_message)std::cout << "Building an ANNIE Event" << std::endl;

  ///////////////LOAD RAW PMT DATA INTO ANNIEEVENT///////////////
  std::map<unsigned long, std::vector<Waveform<uint16_t>> > RawADCData;
  std::map<unsigned long, std::vector<Waveform<uint16_t>> > RawADCAuxData;
  for(std::pair<std::vector<int>, std::vector<uint16_t>> apair : WaveMap){
    int CardID = apair.first.at(0);
    int ChannelID = apair.first.at(1);
    int CrateNum=-1;
    int SlotNum=-1;
    this->CardIDToElectronicsSpace(CardID, CrateNum, SlotNum);
    std::vector<uint16_t> TheWaveform = apair.second;
    Waveform<uint16_t> TheWave(ClockTime, TheWaveform);
    //Placing waveform in a vector in case we want a hefty-mode minibuffer storage eventually
    std::vector<Waveform<uint16_t>> WaveVec{TheWave};
    
    std::vector<int> CrateSpace{CrateNum,SlotNum,ChannelID};
    unsigned long ChannelKey;
    if(TankPMTCrateSpaceToChannelNumMap.count(CrateSpace)>0){
      ChannelKey = TankPMTCrateSpaceToChannelNumMap.at(CrateSpace);
      RawADCData.emplace(ChannelKey,WaveVec);
    }
    else if (AuxCrateSpaceToChannelNumMap.count(CrateSpace)>0){
      ChannelKey = AuxCrateSpaceToChannelNumMap.at(CrateSpace);
      RawADCAuxData.emplace(ChannelKey,WaveVec);
    } else{
      Log("ANNIEEventBuilder:: Cannot find channel key for crate space entry: ",v_error, verbosity);
      Log("ANNIEEventBuilder::CrateNum "+to_string(CrateNum),v_error, verbosity);
      Log("ANNIEEventBuilder::SlotNum "+to_string(SlotNum),v_error, verbosity);
      Log("ANNIEEventBuilder::ChannelID "+to_string(ChannelID),v_error, verbosity);
      Log("ANNIEEventBuilder:: Passing over the wave; PMT DATA LOST",v_error, verbosity);
      continue;
    }
  }
  if(RawADCData.size() == 0){
    std::cout << "No Raw ADC Data in entry.  Not putting to ANNIEEvent." << std::endl;
  }
  std::cout << "Setting ANNIE Event information" << std::endl;
  ANNIEEvent->Set("RawADCData",RawADCData);
  ANNIEEvent->Set("RawADCAuxData",RawADCAuxData);
  ANNIEEvent->Set("EventTimeTank",ClockTime);
  if(verbosity>v_debug) std::cout << "ANNIEEventBuilder: ANNIE Event "+
      to_string(ANNIEEventNum)+" built." << std::endl;
  return;
}


void ANNIEEventBuilder::SaveEntryToFile(int RunNum, int SubRunNum)
{
  /*if(verbosity>4)*/ std::cout << "ANNIEEvent: Saving ANNIEEvent entry"+to_string(ANNIEEventNum) << std::endl;
  std::string Filename = SavePath + ProcessedFilesBasename + "R" + to_string(RunNum) + 
      "S" + to_string(SubRunNum);
  ANNIEEvent->Save(Filename);
  //std::cout <<"ANNIEEvent saved, now delete"<<std::endl;
  ANNIEEvent->Delete();		//Delete() will delete the last entry in the store from memory and enable us to set a new pointer (won't erase the entry from saved file)
  //std::cout <<"ANNIEEvent deleted"<<std::endl;
  ANNIEEventNum+=1;
  return;
}

void ANNIEEventBuilder::CardIDToElectronicsSpace(int CardID, 
        int &CrateNum, int &SlotNum)
{
  //CardID = CrateNum * 1000 + SlotNum.  This logic works if we have less than
  // 10 crates and less than 100 Slots (which we do).
  SlotNum = CardID % 100;
  CrateNum = CardID / 1000;
  return;
}

void ANNIEEventBuilder::OpenNewANNIEEvent(int RunNum, int SubRunNum,uint64_t StarT, int RunT){
  uint64_t StarTime;
  int RunType;
  if(verbosity>v_warning) std::cout << "ANNIEEventBuilder: New run or subrun encountered. Opening new BoostStore" << std::endl;
  if(verbosity>v_debug) std::cout << "ANNIEEventBuilder: Current run,subrun:" << CurrentRunNum << "," << CurrentSubRunNum << std::endl;
  if(verbosity>v_debug) std::cout << "ANNIEEventBuilder: Encountered run,subrun:" << RunNum << "," << SubRunNum << std::endl;
  ANNIEEvent->Close();
  ANNIEEvent->Delete();
  delete ANNIEEvent; ANNIEEvent = new BoostStore(false,2);
  CurrentRunNum = RunNum;
  CurrentSubRunNum = SubRunNum;
  CurrentRunType = RunT;
  CurrentStarTime = StarT;
  if((CurrentRunNum != RunNum)) CurrentDriftMean = 0;
}
