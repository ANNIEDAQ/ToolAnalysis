#include "TriggerDataDecoder.h"

TriggerDataDecoder::TriggerDataDecoder():Tool(){}


bool TriggerDataDecoder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  verbosity = 0;
  TriggerMaskFile = "none";
  TriggerWordFile = "none";
  mode = "EventBuilding";

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("TriggerMaskFile",TriggerMaskFile);
  m_variables.Get("TriggerWordFile",TriggerWordFile);
  m_variables.Get("Mode",mode);

  if (mode != "EventBuilding" && mode != "Monitoring"){
    Log("TriggerDataDecoder tool: Specified mode of operation >> "+mode+" unknown. Use standard EventBuilding mode.",v_error,verbosity);
    mode = "EventBuilding";
  }

  CurrentRunNum = -1;
  CurrentSubrunNum = -1;

  TimeToTriggerWordMap = new std::map<uint64_t,uint32_t>;
  m_data->CStore.Set("PauseCTCDecoding",false);

  if(TriggerMaskFile!="none"){
    TriggerMask = LoadTriggerMask(TriggerMaskFile);
    if(TriggerMask.size()>0) UseTrigMask = true;
  }

  if(TriggerWordFile!="none"){
    TriggerWords = LoadTriggerWords(TriggerWordFile);    //maps trigger words to human-readable labels
    m_data->CStore.Set("TriggerWordMap",TriggerWords);
  }

  return true;
}


bool TriggerDataDecoder::Execute(){

  if (mode == "EventBuilding"){
    m_data->CStore.Set("NewCTCDataAvailable",false);
    bool PauseCTCDecoding = false;
    m_data->CStore.Get("PauseCTCDecoding",PauseCTCDecoding);
    if (PauseCTCDecoding){
      std::cout << "TriggerDataDecoder tool: Pausing trigger decoding to let Tank and MRD data catch up..." << std::endl;
      return true;
    }
    //Clear decoding maps if a new run/subrun is encountered
    this->CheckForRunChange();
  
    //Get the TriggerData vector pointer from the CStore
    Log("TriggerDataDecoder Tool: Accessing TrigData vector in CStore",v_debug, verbosity);
    bool got_tdata = m_data->CStore.Get("TrigData",Tdata);
    if(!got_tdata){
      if(verbosity>0) std::cout << "TriggerDataDecoder error: No TriggerData in CStore!" << std::endl;
      return false;
    }
    bool new_ts_available = false;
    std::vector<uint32_t> aTimeStampData = Tdata->TimeStampData;
    for(int i = 0; i < (int) aTimeStampData.size(); i++){
      if(verbosity>v_debug) std::cout<<"TriggerDataDecoder Tool: Loading next TrigData from entry's index " << i <<std::endl;
      new_ts_available = this->AddWord(aTimeStampData.at(i));
      if(new_ts_available){
        if(verbosity>4){
          std::cout << "PARSED TRIGGER TIME: " << processed_ns.back() << std::endl;
          std::cout << "PARSED TRIGGER WORD: " << processed_sources.back() << std::endl;
        }
        if(UseTrigMask){
          uint32_t recent_trigger_word = processed_sources.back();
          for(int j = 0; j<(int) TriggerMask.size(); j++){
            if(TriggerMask.at(j) == (int) recent_trigger_word){
              m_data->CStore.Set("NewCTCDataAvailable",true);
              if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
              TimeToTriggerWordMap->emplace(processed_ns.back(),processed_sources.back());
            }
          }
        } else {
          m_data->CStore.Set("NewCTCDataAvailable",true);
          if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
          TimeToTriggerWordMap->emplace(processed_ns.back(),processed_sources.back());
        }
      }
    }
  } 
  else if (mode == "Monitoring"){
    TimeToTriggerWordMap->clear();
    processed_sources.clear();
    processed_ns.clear();
    std::map<int,TriggerData> TrigData_Map;
    m_data->Stores["TrigData"]->Get("TrigDataMap",TrigData_Map);
    bool new_ts_available = false;
    for (int i_entry=0; i_entry < int(TrigData_Map.size()); i_entry++){
      TriggerData TData = TrigData_Map.at(i_entry);
      std::vector<uint32_t> aTimeStampData = TData.TimeStampData;
      for (int i=0; i < (int) aTimeStampData.size(); i++){
        new_ts_available = this->AddWord(aTimeStampData.at(i));
        if (new_ts_available){
          if(verbosity>4){
            std::cout << "PARSED TRIGGER TIME: " << processed_ns.back() << std::endl;
            std::cout << "PARSED TRIGGER WORD: " << processed_sources.back() << std::endl;
          }
          if(UseTrigMask){
            uint32_t recent_trigger_word = processed_sources.back();
            for(int j = 0; j<(int) TriggerMask.size(); j++){
              if(TriggerMask.at(j) == (int) recent_trigger_word){
                m_data->CStore.Set("NewCTCDataAvailable",true);
                if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
                TimeToTriggerWordMap->emplace(processed_ns.back(),processed_sources.back());
              }
            }
          } else {
            m_data->CStore.Set("NewCTCDataAvailable",true);
            if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
            TimeToTriggerWordMap->emplace(processed_ns.back(),processed_sources.back());
          }
        }     
      }
    }
  }

  if(verbosity>3) Log("TriggerDataDecoder Tool: size of TimeToTriggerWordMap: "+to_string(TimeToTriggerWordMap->size()),v_message,verbosity); 
 
  m_data->CStore.Set("TimeToTriggerWordMap",TimeToTriggerWordMap);

  return true;
}


bool TriggerDataDecoder::Finalise(){
  //delete TimeToTriggerWordMap;	//DONT delete TimeToTriggerWordMap since it wil be deleted by the CStore automatically
  std::cout << "TriggerDataDecoder tool exitting" << std::endl;
  return true;
}

bool TriggerDataDecoder::AddWord(uint32_t word){
  bool new_timestamp_available = false;
  uint64_t payload = 0;
  uint32_t wordid = 0;
  if(word == 0xF1F0E5E7){
    have_c1 = false;
    c1 = 0;
    have_c2 = false;
    c2 = 0;
    fiforesets.push_back(processed_sources.size()); 
    return new_timestamp_available;
  }

  wordid = word>>24;
  payload = word & 0x00FFFFFF;

  if(wordid == 0xC0){
    return new_timestamp_available;
  }
  else if (wordid == 0xC1 || wordid == 0xC3){
    c1 = (payload&0x0000FFFF)<<24;
    have_c1 = true;
  } else if (wordid == 0xC2 || wordid == 0xC4){
    c2 = payload << 40;
    have_c2 = true;
  } else {
    if(have_c1 && have_c2){
      uint64_t ns = (c1 + c2 + payload)*8;
      processed_sources.push_back(wordid);
      processed_ns.push_back(ns);
      new_timestamp_available = true;
    }
  }
  return new_timestamp_available;
}

void TriggerDataDecoder::CheckForRunChange()
{
  // Load RawData BoostStore to use in execute loop
  Store RunInfoPostgress;
  m_data->CStore.Get("RunInfoPostgress",RunInfoPostgress);
  int RunNumber;
  int SubRunNumber;
  uint64_t StarTime;
  int RunType;
  RunInfoPostgress.Get("RunNumber",RunNumber);
  RunInfoPostgress.Get("SubRunNumber",SubRunNumber);
  RunInfoPostgress.Get("RunType",RunType);
  RunInfoPostgress.Get("StarTime",StarTime);

  //If we have moved onto a new run number, we should clear the event building maps 
  if (CurrentRunNum == -1){
    CurrentRunNum = RunNumber;
    CurrentSubrunNum = SubRunNumber;
  }
  else if (RunNumber != CurrentRunNum){ //New run has been encountered
    Log("TriggerDataDecoder Tool: New run encountered.  Clearing event building maps",v_message,verbosity); 
    TimeToTriggerWordMap->clear();
  }
  else if (SubRunNumber != CurrentSubrunNum){ //New run has been encountered
    Log("TriggerDataDecoder Tool: New subrun encountered.",v_message,verbosity); 
    TimeToTriggerWordMap->clear();
  }
  return;
}

std::vector<int> TriggerDataDecoder::LoadTriggerMask(std::string triggermask_file){
  std::vector<int> trigger_mask;
  std::string fileline;
  ifstream myfile(triggermask_file.c_str());
  if (myfile.is_open()){
    while(getline(myfile,fileline)){
      if(fileline.find("#")!=std::string::npos) continue;
      std::cout << fileline << std::endl; //has our stuff;
      std::vector<std::string> dataline;
      boost::split(dataline,fileline, boost::is_any_of(","), boost::token_compress_on);
      uint32_t triggernum = std::stoul(dataline.at(0));
      if(verbosity>4) std::cout << "Trigger mask will have trigger number " << triggernum << std::endl;
      trigger_mask.push_back(triggernum);
    }
  } else {
    Log("TriggerDataDecoder Tool: Input trigger mask file not found. "
        " all triggers from CTC will attempt to be paired with PMT/MRD data. ",
        v_warning, verbosity);
  }
  return trigger_mask;
}

std::map<int,std::string> TriggerDataDecoder::LoadTriggerWords(std::string triggerwords_file){
  std::map<int,std::string> triggerwordmap;
  int triggerword;
  std::string triggerlabel;
  ifstream myfile(triggerwords_file.c_str());
  if (myfile.is_open()){
    while (!myfile.eof()){
      myfile >> triggerword >> triggerlabel;
      triggerwordmap.emplace(triggerword,triggerlabel);
      if (myfile.eof()) break;
    }
  } else {
    Log("TriggerDataDecoder Tool: Input trigger words file not found. "
	" Please check why "+triggerwords_file+" cannot be found.",
	v_warning,verbosity);
  }

  return triggerwordmap;
}
