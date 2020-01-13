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

  /////////////////////////////////////////////////////////////////
  //FIXME: Eventually, RunNumber should be loaded from A run database
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("SavePath",SavePath);
  m_variables.Get("ProcessedFilesBasename",ProcessedFilesBasename);
  m_variables.Get("BuildType",BuildType);

  if(BuildType == "PMTAndMRD"){
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

  //////////////////////initialize subrun index//////////////
  ANNIEEvent = new BoostStore(false,2);
  ANNIEEventNum = 0;
  CurrentRunNum = -1;
  CurrentSubrunNum = -1;

  return true;
}


bool ANNIEEventBuilder::Execute(){

  if (BuildType == "Tank"){
    //Check to see if there's new PMT data
    m_data->CStore.Get("NewTankPMTDataAvailable",IsNewTankData);
    if(!IsNewTankData){
      Log("ANNIEEventBuilder:: No new Tank Data.  Not building ANNIEEvent. ",v_message, verbosity);
      return true;
    }
    //Get the current InProgressTankEvents map
    if(verbosity>4) std::cout << "ANNIEEventBuilder: Getting waves and run info from CStore" << std::endl;
    m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);
    m_data->CStore.Get("TankRunInfoPostgress",RunInfoPostgress);
    int RunNumber;
    int SubRunNumber;
    uint64_t StarTime;
    int RunType;
    RunInfoPostgress.Get("RunNumber",RunNumber);
    RunInfoPostgress.Get("SubRunNumber",SubRunNumber);
    RunInfoPostgress.Get("RunType",RunType);
    RunInfoPostgress.Get("StarTime",StarTime);

    //Initialize Current RunNum and SubrunNum. New ANNIEEvent for any New Run or Subrun
    if(CurrentRunNum == -1) CurrentRunNum = RunNumber;
    if(CurrentSubrunNum == -1) CurrentSubrunNum = SubRunNumber;
    //If we're in a new run or subrun, make a new ANNIEEvent file. 
    if((CurrentRunNum != RunNumber) || (CurrentSubrunNum != SubRunNumber)){
      if(verbosity>v_warning) std::cout << "New run or subrun encountered. Opening new BoostStore" << std::endl;
      if(verbosity>v_warning) std::cout << "ANNIEEventBuilder: Saving and closing file." << std::endl;
      ANNIEEvent->Close();
      ANNIEEvent->Delete();
      delete ANNIEEvent; ANNIEEvent = new BoostStore(false,2);
      CurrentRunNum = RunNumber;
      CurrentSubrunNum = SubRunNumber;
    }

    //Assume a whole processed file will have all it's PMT data finished
    std::vector<uint64_t> PMTEventsToDelete;
    for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : InProgressTankEvents){
      if(verbosity>4) std::cout << "Accessing next PMT counter" << std::endl;
      uint64_t PMTCounterTime = apair.first;
      if(verbosity>4) std::cout << "Finished waveset has clock counter: " << PMTCounterTime << std::endl;
      std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
      if(verbosity>4) std::cout << "Number of waves for this counter: " << aWaveMap.size() << std::endl;
      //For this counter, need to have the number of TankPMT channels plus number of aux channels
      int NumTankPMTChannels = TankPMTCrateSpaceToChannelNumMap.size();
      int NumAuxChannels = AuxCrateSpaceToChannelNumMap.size();
      if(aWaveMap.size() >= (NumTankPMTChannels + NumAuxChannels)){
        this->BuildANNIEEventTank(PMTCounterTime, aWaveMap,RunNumber,SubRunNumber,RunType,StarTime);
        this->SaveEntryToFile(CurrentRunNum,CurrentSubrunNum);
        //Erase this entry from the InProgressTankEventsMap
        if(verbosity>4) std::cout << "Counter time will be erased from InProgressTankEvents: " << PMTCounterTime << std::endl;
        PMTEventsToDelete.push_back(PMTCounterTime);
      }
    }

    for(unsigned int i=0; i< PMTEventsToDelete.size(); i++) InProgressTankEvents.erase(PMTEventsToDelete.at(i));
    //Update the current InProgressTankEvents map
    m_data->CStore.Set("InProgressTankEvents",InProgressTankEvents);

  } else if (BuildType == "MRD"){
    m_data->CStore.Get("NewMRDDataAvailable",IsNewMRDData);
    std::vector<unsigned long> MRDEventsToDelete;
    if(!IsNewMRDData){
      Log("ANNIEEventBuilder:: No new MRD Data.  Not building ANNIEEvent: ",v_message, verbosity);
      return true;
    }
    m_data->CStore.Get("MRDEvents",MRDEvents);
    m_data->CStore.Get("MRDEventTriggerTypes",TriggerTypeMap);
    
    m_data->CStore.Get("MRDRunInfoPostgress",RunInfoPostgress);
    int RunNumber;
    int SubRunNumber;
    uint64_t StarTime;
    int RunType;
    RunInfoPostgress.Get("RunNumber",RunNumber);
    RunInfoPostgress.Get("SubRunNumber",SubRunNumber);
    RunInfoPostgress.Get("RunType",RunType);
    RunInfoPostgress.Get("StarTime",StarTime);

   //Loop through MRDEvents and process each into ANNIEEvent.
    for(std::pair<unsigned long,std::vector<std::pair<unsigned long,int>>> apair : MRDEvents){
      unsigned long MRDTimeStamp = apair.first;
      std::vector<std::pair<unsigned long,int>> MRDHits = apair.second;
      std::string MRDTriggerType = TriggerTypeMap[MRDTimeStamp];
      this->BuildANNIEEventMRD(MRDHits, MRDTimeStamp, MRDTriggerType, RunNumber, SubRunNumber, RunType);
      this->SaveEntryToFile(RunNumber,SubRunNumber);
      //Erase this entry from the InProgressTankEventsMap
      MRDEventsToDelete.push_back(MRDTimeStamp);
    }
    for (unsigned int i=0; i< MRDEventsToDelete.size(); i++){
      MRDEvents.erase(MRDEventsToDelete.at(i));
      TriggerTypeMap.erase(MRDEventsToDelete.at(i));
    }
  }

  else if (BuildType == "TankAndMRD"){
    //See if the MRD and Tank are at the same run/subrun for building
    m_data->CStore.Get("TankRunInfoPostgress",RunInfoPostgress);
    int TankRunNumber;
    int TankSubRunNumber;
    uint64_t TankStarTime;
    int TankRunType;
    RunInfoPostgress.Get("RunNumber",TankRunNumber);
    RunInfoPostgress.Get("SubRunNumber",TankSubRunNumber);
    RunInfoPostgress.Get("RunType",TankRunType);
    RunInfoPostgress.Get("StarTime",TankStarTime);
    
    m_data->CStore.Get("MRDRunInfoPostgress",RunInfoPostgress);
    int MRDRunNumber;
    int MRDSubRunNumber;
    uint64_t MRDStarTime;
    int MRDRunType;
    RunInfoPostgress.Get("RunNumber",MRDRunNumber);
    RunInfoPostgress.Get("SubRunNumber",MRDSubRunNumber);
    RunInfoPostgress.Get("RunType",TankRunType);
    RunInfoPostgress.Get("StarTime",TankStarTime);
    
    //Initialize Current RunNum and SubrunNum. New ANNIEEvent for any New Run or Subrun
    if(CurrentRunNum == -1 || CurrentSubrunNum == -1){
      CurrentRunNum = TankRunNumber;
      CurrentSubrunNum = TankSubRunNumber;
      CurrentStarTime = TankStarTime;
      CurrentRunType = TankRunType;
      LowestRunNumber = TankRunNumber;
      LowestSubRunNumber = TankSubRunNumber;
    }
  
    // Check that Tank and MRD decoding are on the same run/subrun.  If not, take the
    // lowest Run Number and Subrun number to keep building 
    if (TankRunNumber != MRDRunNumber || TankSubRunNumber != MRDSubRunNumber){
      if((TankRunNumber*10000 + TankSubRunNumber) > (MRDRunNumber*10000 + MRDSubRunNumber)){
        LowestRunNumber = TankRunNumber;
        LowestSubRunNumber = TankSubRunNumber;
        m_data->CStore.Set("PauseMRDDecoding",true);
      } else {
        LowestRunNumber = MRDRunNumber;
        LowestSubRunNumber = MRDSubRunNumber;
        m_data->CStore.Set("PauseTankDecoding",true);
      }
    } else {
      m_data->CStore.Set("PauseMRDDecoding",false);
      m_data->CStore.Set("PauseTankDecoding",false);
    }

    //If we're in a new run or subrun, make a new ANNIEEvent file. 
    if((CurrentRunNum != LowestRunNumber) || (CurrentSubrunNum != LowestSubRunNumber)){
      if(verbosity>v_warning) std::cout << "PMT and MRD data have both finished run " <<
         CurrentRunNum << "," << CurrentSubRunNum << ".  Opening new BoostStore" << std::endl;
      if(verbosity>v_warning) std::cout << "ANNIEEventBuilder: Saving and closing file." << std::endl;
      ANNIEEvent->Close();
      ANNIEEvent->Delete();
      delete ANNIEEvent; ANNIEEvent = new BoostStore(false,2);
      CurrentRunNum = LowestRunNumber;
      CurrentSubrunNum = LowestSubRunNumber;
    }

    //First, see if any In-progress tank events now have all waveforms
    m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);
    std::vector<uint64_t> InProgressEventsToDelete;
    m_data->CStore.Get("NewTankPMTDataAvailable",IsNewTankData);
    if(!IsNewTankData){
      for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : InProgressTankEvents){
        if(verbosity>4) std::cout << "Accessing next PMT counter" << std::endl;
        uint64_t PMTCounterTime = apair.first;
        //See if this is a new timestamp; let's add it to our timestamp tracker if so
        if(std::find(RunTankTimestamps.begin(),RunTankTimestamps.end(), PMTCounterTime) == 
                     RunTankTimestamps.end()){
          RunTankTimestamps.push_back(PMTCounterTime);
        }
        std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
        //For this counter, need to have the number of TankPMT channels plus number of aux channels
        int NumTankPMTChannels = TankPMTCrateSpaceToChannelNumMap.size();
        int NumAuxChannels = AuxCrateSpaceToChannelNumMap.size();
        if(aWaveMap.size() >= (NumTankPMTChannels + NumAuxChannels)){
          FinishedTankEvents.emplace(PMTCounterTime,aWaveMap);
          FinishedTankTimestamps.push_back(PMTCounterTime);
          //Put PMT timestamp into the timestamp set for this run.
          if(verbosity>4) std::cout << "Finished waveset has clock counter: " << PMTCounterTime << std::endl;
          if(verbosity>4) std::cout << "Number of waves for this counter: " << aWaveMap.size() << std::endl;
          InProgressEventsToDelete.push_back(PMTCounterTime);
        }
        //Erase this entry from the InProgressTankEventsMap
        if(verbosity>4) std::cout << "Counter time will be erased from InProgressTankEvents: " << PMTCounterTime << std::endl;
        InProgressEventsToDelete.push_back(PMTCounterTime);
      }
    }

    //Now, look through our MRD data to see if anything is new
    m_data->CStore.Get("NewMRDDataAvailable",IsNewMRDData);
    m_data->CStore.Get("MRDEvents",MRDEvents);
    std::vector<unsigned long> MRDEventsToDelete;
    if(IsNewMRDData){
      for(std::pair<unsigned long,std::vector<std::pair<unsigned long,int>>> apair : MRDEvents){
        unsigned long MRDTimeStamp = apair.first;
        std::vector<std::pair<unsigned long,int>> MRDHits = apair.second;
        //Check if any timestamps are new
        if(std::find(RunMRDTimestamps.begin(),RunMRDTimestamps.end(), PMTCounterTime) == 
                     RunMRDTimestamps.end()){
          RunMRDTimestamps.push_back(PMTCounterTime);
        }
      }
    }

    //Now, come up with an algorithm to pair up PMT and MRD events...

    //Finally, Build the ANNIEEvent of any PMT/MRD data that is done/has been paired
    m_data->CStore.Get("MRDEventTriggerTypes",TriggerTypeMap);
    //Loop through MRDEvents and process each into ANNIEEvent.
    for(std::pair<unsigned long,std::vector<std::pair<unsigned long,int>>> apair : MRDEvents){
      unsigned long MRDTimeStamp = apair.first;
      std::vector<std::pair<unsigned long,int>> MRDHits = apair.second;
      std::string MRDTriggerType = TriggerTypeMap[MRDTimeStamp];
      this->BuildANNIEEvent(PMTCounterTime, aWaveMap, MRDHits, MRDTimeStamp, MRDTriggerType, CurrentRunNumber, CurrentSubRunNumber, CurrentRunType,CurrentStarTime);
      this->SaveEntryToFile(RunNumber,SubRunNumber);
      //Erase this entry from the InProgressTankEventsMap
      MRDEventsToDelete.push_back(MRDTimeStamp);
    }
    for (unsigned int i=0; i< MRDEventsToDelete.size(); i++){
      MRDEvents.erase(MRDEventsToDelete.at(i));
      TriggerTypeMap.erase(MRDEventsToDelete.at(i));
    }
  }
    for(unsigned int i=0; i< InProgressEventsToDelete.size(); i++) InProgressTankEvents.erase(InProgressEventsToDelete.at(i));
    //Update the current InProgressTankEvents map
    m_data->CStore.Set("InProgressTankEvents",InProgressTankEvents);



  return true;
}


bool ANNIEEventBuilder::Finalise(){
  if(verbosity>4) std::cout << "ANNIEEvent Finalising.  Closing any open ANNIEEvent Boostore" << std::endl;
  if(verbosity>2) std::cout << "ANNIEEventBuilder: Saving and closing file." << std::endl;
  ANNIEEvent->Close();
  ANNIEEvent->Delete();
  delete ANNIEEvent;
  //Save the current subrun and delete ANNIEEvent
  std::cout << "ANNIEEventBuilder Exitting" << std::endl;
  return true;
}


void ANNIEEventBuilder::BuildANNIEEventMRD(std::vector<std::pair<unsigned long,int>> MRDHits, 
        unsigned long MRDTimeStamp, std::string MRDTriggerType, int RunNum, int SubrunNum, int
        RunType)
{
  std::cout << "Building an ANNIE Event (MRD), ANNIEEventNum = "<<ANNIEEventNum << std::endl;

  TDCData = new std::map<unsigned long, std::vector<Hit>>;

  //TODO: Loop through MRDHits at this timestamp and form the Hit vector.
  for (unsigned int i_value=0; i_value< MRDHits.size(); i_value++){
    unsigned long channelkey = MRDHits.at(i_value).first;
    int hitTimeADC = MRDHits.at(i_value).second;
    if (TDCData->count(channelkey)==0){
      std::vector<Hit> newhitvector;
      if (verbosity > 3) std::cout <<"creating hit with time value "<<hitTimeADC*4<<"and chankey "<<channelkey<<std::endl;
      newhitvector.push_back(Hit(0,hitTimeADC*4.,1.));    //Hit(tubeid, time, charge). 1 TDC tick corresponds to 4ns, no charge information (set to 1)
      TDCData->emplace(channelkey,newhitvector);
    } else {
      if (verbosity > 3) std::cout <<"creating hit with time value "<<hitTimeADC*4<<"and chankey "<<channelkey<<std::endl;
      TDCData->at(channelkey).push_back(Hit(0,hitTimeADC*4.,1.));
    }
  }

  Log("ANNIEEventBuilder: TDCData size: "+std::to_string(TDCData->size()),v_debug,verbosity);

  ANNIEEvent->Set("TDCData",TDCData,true);
  ANNIEEvent->Set("RunNumber",RunNum);
  ANNIEEvent->Set("SubrunNumber",SubrunNum);
  ANNIEEvent->Set("RunType",RunType);
  ANNIEEvent->Set("EventNumber",ANNIEEventNum);
  TimeClass timeclass_timestamp((uint64_t)MRDTimeStamp*1000);  //in ns
  ANNIEEvent->Set("EventTime",timeclass_timestamp); //not sure if EventTime is also in UTC or defined differently
  ANNIEEvent->Set("MRDTriggerType",MRDTriggerType);
  return;
}

void ANNIEEventBuilder::BuildANNIEEventTank(uint64_t ClockTime, 
        std::map<std::vector<int>, std::vector<uint16_t>> WaveMap, int RunNum, int SubrunNum,
        int RunType, uint64_t StartTime)
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
    if(verbosity>v_debug) std::cout << "Converting card ID " << CardID << ", channel ID " <<
          ChannelID << " to electronics space" << std::endl;
    this->CardIDToElectronicsSpace(CardID, CrateNum, SlotNum);
    std::vector<uint16_t> TheWaveform = apair.second;
    //FIXME: We're feeding Waveform class expects a double, not a uint64_t (?)
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
  ANNIEEvent->Set("RunNumber",RunNum);
  ANNIEEvent->Set("SubrunNumber",SubrunNum);
  ANNIEEvent->Set("RunType",RunType);
  ANNIEEvent->Set("RunStartTime",StartTime);
  //TODO: Things missing from ANNIEEvent that should be in before this tool finishes:
  //  - EventTime
  //  - TriggerData
  //  - BeamStatus?  
  //  - RawLAPPDData
  if(verbosity>v_debug) std::cout << "ANNIEEventBuilder: ANNIE Event "+
      to_string(ANNIEEventNum)+" built." << std::endl;
  return;
}


void ANNIEEventBuilder::BuildANNIEEvent( uint64_t TankClockTime, 
        std::map<std::vector<int>, std::vector<uint16_t>> WaveMap, std::vector<std::pair<unsigned long,int>> MRDHits, 
        unsigned long MRDTimeStamp, std::string MRDTriggerType, int RunNum, int SubrunNum, int
        RunType, uint64_t RunStartTime)
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
    if(verbosity>v_debug) std::cout << "Converting card ID " << CardID << ", channel ID " <<
          ChannelID << " to electronics space" << std::endl;
    this->CardIDToElectronicsSpace(CardID, CrateNum, SlotNum);
    std::vector<uint16_t> TheWaveform = apair.second;
    //FIXME: We're feeding Waveform class expects a double, not a uint64_t (?)
    Waveform<uint16_t> TheWave(TankClockTime, TheWaveform);
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
  ANNIEEvent->Set("RawADCData",RawADCData);
  ANNIEEvent->Set("RawADCAuxData",RawADCAuxData);
  ANNIEEvent->Set("RunNumber",RunNum);
  ANNIEEvent->Set("SubrunNumber",SubrunNum);
  ANNIEEvent->Set("RunType",RunType);
  ANNIEEvent->Set("RunStartTime",StartTime);

  ///////////////LOAD RAW MRD DATA INTO ANNIEEVENT///////////////
  TDCData = new std::map<unsigned long, std::vector<Hit>>;
  //TODO: Loop through MRDHits at this timestamp and form the Hit vector.
  for (unsigned int i_value=0; i_value< MRDHits.size(); i_value++){
    unsigned long channelkey = MRDHits.at(i_value).first;
    int hitTimeADC = MRDHits.at(i_value).second;
    if (TDCData->count(channelkey)==0){
      std::vector<Hit> newhitvector;
      if (verbosity > 3) std::cout <<"creating hit with time value "<<hitTimeADC*4<<"and chankey "<<channelkey<<std::endl;
      newhitvector.push_back(Hit(0,hitTimeADC*4.,1.));    //Hit(tubeid, time, charge). 1 TDC tick corresponds to 4ns, no charge information (set to 1)
      TDCData->emplace(channelkey,newhitvector);
    } else {
      if (verbosity > 3) std::cout <<"creating hit with time value "<<hitTimeADC*4<<"and chankey "<<channelkey<<std::endl;
      TDCData->at(channelkey).push_back(Hit(0,hitTimeADC*4.,1.));
    }
  }

  Log("ANNIEEventBuilder: TDCData size: "+std::to_string(TDCData->size()),v_debug,verbosity);

  ANNIEEvent->Set("TDCData",TDCData,true);
  TimeClass timeclass_timestamp((uint64_t)MRDTimeStamp*1000);  //in ns
  ANNIEEvent->Set("MRDEventTime",timeclass_timestamp); //not sure if EventTime is also in UTC or defined differently
  ANNIEEvent->Set("MRDTriggerType",MRDTriggerType);

  //TODO: Things missing from ANNIEEvent that should be in before this tool finishes:
  //  - EventTime
  //  - TriggerData
  //  - BeamStatus?  
  //  - RawLAPPDData
  return;
}

void ANNIEEventBuilder::BuildANNIEEventTank( int RunNum, int SubrunNum,
        int RunType, uint64_t StartTime)
{
  return;
}



void ANNIEEventBuilder::SaveEntryToFile(int RunNum, int SubrunNum)
{
  /*if(verbosity>4)*/ std::cout << "ANNIEEvent: Saving ANNIEEvent entry"+to_string(ANNIEEventNum) << std::endl;
  std::string Filename = SavePath + ProcessedFilesBasename + "R" + to_string(RunNum) + 
      "S" + to_string(SubrunNum);
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
