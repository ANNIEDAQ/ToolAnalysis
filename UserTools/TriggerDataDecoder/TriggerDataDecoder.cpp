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
  readtrigoverlap = 0;
  storetrigoverlap = 0;
  usecstore = 1;

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("TriggerMaskFile",TriggerMaskFile);
  m_variables.Get("TriggerWordFile",TriggerWordFile);
  m_variables.Get("Mode",mode);
  m_variables.Get("ReadTrigOverlap",readtrigoverlap);
  m_variables.Get("StoreTrigOverlap",storetrigoverlap);
  m_variables.Get("UseCStore",usecstore);

  if (mode != "EventBuilding" && mode != "Monitoring"){
    Log("TriggerDataDecoder tool: Specified mode of operation >> "+mode+" unknown. Use standard EventBuilding mode.",v_error,verbosity);
    mode = "EventBuilding";
  }

  CurrentRunNum = -1;
  CurrentSubrunNum = -1;
  CurrentPartNum = -1;

  TimeToTriggerWordMap = new std::map<uint64_t,std::vector<uint32_t>>;
  TimeToTriggerWordMapComplete = new std::map<uint64_t,std::vector<uint32_t>>;
  m_data->CStore.Set("PauseCTCDecoding",false);

  if(TriggerMaskFile!="none"){
    TriggerMask = LoadTriggerMask(TriggerMaskFile);
    if(TriggerMask.size()>0) UseTrigMask = true;
  }

  if(TriggerWordFile!="none"){
    TriggerWords = LoadTriggerWords(TriggerWordFile);    //maps trigger words to human-readable labels
    m_data->CStore.Set("TriggerWordMap",TriggerWords);
  }

  loop_nr = 0;


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
  
    //Read in trigger overlap information from previous file (if desired)
    if (loop_nr == 0){
      if(readtrigoverlap){
        BoostStore ReadTrigOverlap;
        std::stringstream ss_trigoverlap;
        ss_trigoverlap << "TrigOverlap_R"<<CurrentRunNum<<"S"<<CurrentSubrunNum<<"p"<<CurrentPartNum;
        bool got_trig = ReadTrigOverlap.Initialise(ss_trigoverlap.str().c_str());
        if (got_trig && (CurrentPartNum!=0)){
          ReadTrigOverlap.Get("c1",c1);
          ReadTrigOverlap.Get("c2",c2);
          Log("TriggerDataDecoder tool: Got c1 = "+std::to_string(c1)+" and c2 = "+std::to_string(c2)+" from C1C2 store located at "+ss_trigoverlap.str(),v_message,verbosity);
          have_c1 = true;
          have_c2 = true;
        } else {
          Log("TriggerDataDecoder tool: Did not find trigger overlap file: "+ss_trigoverlap.str(),v_error,verbosity);
        }
      }
    }

    //Get the TriggerData vector pointer from the CStore
    Log("TriggerDataDecoder Tool: Accessing TrigData vector in CStore",v_debug, verbosity);
    bool got_tdata = m_data->CStore.Get("TrigData",Tdata);
    if(!got_tdata){
      if(verbosity>0) std::cout << "TriggerDataDecoder error: No TriggerData in CStore!" << std::endl;
      return false;
    }
    bool new_ts_available = false;
    std::vector<uint32_t> aTimeStampData = Tdata->TimeStampData;
    std::cout <<"aTimeStampData.size(): "<<aTimeStampData.size()<<std::endl;
    for(int i = 0; i < (int) aTimeStampData.size(); i++){
      if(verbosity>v_debug) std::cout<<"TriggerDataDecoder Tool: Loading next TrigData from entry's index " << i <<std::endl;
      new_ts_available = this->AddWord(aTimeStampData.at(i));
      if(new_ts_available){
        if(verbosity>4){
          std::cout << "PARSED TRIGGER TIME: " << processed_ns.back() << std::endl;
          std::cout << "PARSED TRIGGER WORD: " << processed_sources.back() << std::endl;
        }
        if (TimeToTriggerWordMapComplete->find(processed_ns.back()) != TimeToTriggerWordMapComplete->end()) TimeToTriggerWordMapComplete->at(processed_ns.back()).push_back(processed_sources.back());
        else {
          std::vector<uint32_t> timestamp_ns{processed_sources.back()};
          TimeToTriggerWordMapComplete->emplace(processed_ns.back(),timestamp_ns);
        }
        if(UseTrigMask){
          uint32_t recent_trigger_word = processed_sources.back();
          for(int j = 0; j<(int) TriggerMask.size(); j++){
            if(TriggerMask.at(j) == recent_trigger_word){
              m_data->CStore.Set("NewCTCDataAvailable",true);
              if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
              if (TimeToTriggerWordMap->find(processed_ns.back()) != TimeToTriggerWordMap->end()) TimeToTriggerWordMap->at(processed_ns.back()).push_back(processed_sources.back());
              else {
                std::vector<uint32_t> timestamp_ns{processed_sources.back()};
                TimeToTriggerWordMap->emplace(processed_ns.back(),timestamp_ns);
              }
            }
          }
        } else {
          m_data->CStore.Set("NewCTCDataAvailable",true);
          if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
          if (TimeToTriggerWordMap->find(processed_ns.back()) != TimeToTriggerWordMap->end()) TimeToTriggerWordMap->at(processed_ns.back()).push_back(processed_sources.back());
          else {
            std::vector<uint32_t> timestamp_ns{processed_sources.back()};
            TimeToTriggerWordMap->emplace(processed_ns.back(),timestamp_ns);
          }
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
              if(TriggerMask.at(j) == recent_trigger_word){
                m_data->CStore.Set("NewCTCDataAvailable",true);
                if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
                if (TimeToTriggerWordMap->find(processed_ns.back()) != TimeToTriggerWordMap->end()) TimeToTriggerWordMap->at(processed_ns.back()).push_back(processed_sources.back());
                else {
                  std::vector<uint32_t> timestamp_ns{processed_sources.back()};
                  TimeToTriggerWordMap->emplace(processed_ns.back(),timestamp_ns);
                }
              }
            }
          } else {
            m_data->CStore.Set("NewCTCDataAvailable",true);
            if(verbosity>4) std::cout << "TRIGGER WORD BEING ADDED TO TRIGWORDMAP" << std::endl;
            if (TimeToTriggerWordMap->find(processed_ns.back()) != TimeToTriggerWordMap->end()) TimeToTriggerWordMap->at(processed_ns.back()).push_back(processed_sources.back());
            else {
              std::vector<uint32_t> timestamp_ns{processed_sources.back()};
              TimeToTriggerWordMap->emplace(processed_ns.back(),timestamp_ns);
            }
          }
        }     
      }
    }
  }

  if(verbosity>3) Log("TriggerDataDecoder Tool: size of TimeToTriggerWordMap: "+to_string(TimeToTriggerWordMap->size()),v_message,verbosity); 
 
  if (usecstore){
    m_data->CStore.Set("TimeToTriggerWordMap",TimeToTriggerWordMap);
    m_data->CStore.Set("TimeToTriggerWordMapComplete",TimeToTriggerWordMapComplete);
  }

  loop_nr++;

  return true;
}


bool TriggerDataDecoder::Finalise(){
  //delete TimeToTriggerWordMap;	//DONT delete TimeToTriggerWordMap since it wil be deleted by the CStore automatically
  std::cout << "TriggerDataDecoder tool exitting" << std::endl;
/*
  if (storec1c2 != "none"){
    ofstream StoreC1C2(storec1c2.c_str());
    StoreC1C2 << "c1    "<<c1 << std::endl;
    StoreC1C2 << "c2    "<<c2 << std::endl;
    StoreC1C2.close();
    
    StoreC1C2.Set("c1",c1);
    StoreC1C2.Set("c2",c2);
    StoreC1C2.Save(storec1c2.c_str());
  }*/

  return true;
}

bool TriggerDataDecoder::AddWord(uint32_t word){
  bool new_timestamp_available = false;
  uint64_t payload = 0;
  uint32_t wordid = 0;
  //std::cout <<"word: "<<word<<std::endl;
  if(word == 0xF1F0E5E7){
    //std::cout <<"word = 0xF1F0E5E7"<<std::endl;
    have_c1 = false;
    c1 = 0;
    have_c2 = false;
    c2 = 0;
    fiforesets.push_back(processed_sources.size()); 
    return new_timestamp_available;
  }

  wordid = word>>24;
  payload = word & 0x00FFFFFF;
  //std::cout <<"wordid: "<<wordid<<", payload: "<<payload<<std::endl;

  if(wordid == 0xC0){
    //std::cout <<"wordid = 0xC0"<<std::endl;
    return new_timestamp_available;
  }
  else if (wordid == 0xC1 || wordid == 0xC3){
    //std::cout <<"wordid is 0xC1 or 0xC3"<<std::endl;
    c1 = (payload&0x0000FFFF)<<24;
    have_c1 = true;
  } else if (wordid == 0xC2 || wordid == 0xC4){
    //std::cout <<"wordid is 0xC2 or 0xC4"<<std::endl;
    c2 = payload << 40;
    have_c2 = true;
  } else {
    //std::cout <<"have_c1: "<<have_c1<<std::endl;
    //std::cout <<"have_c2: "<<have_c2<<std::endl;
    //std::cout <<"c1: "<<c1<<", c2: "<<c2<<std::endl;
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
  int PartNumber;
  uint64_t StarTime;
  int RunType;
  RunInfoPostgress.Get("RunNumber",RunNumber);
  RunInfoPostgress.Get("SubRunNumber",SubRunNumber);
  RunInfoPostgress.Get("PartNumber",PartNumber);
  RunInfoPostgress.Get("RunType",RunType);
  RunInfoPostgress.Get("StarTime",StarTime);

  //If we have moved onto a new run number, we should clear the event building maps 
  if (CurrentRunNum == -1){
    CurrentRunNum = RunNumber;
    CurrentSubrunNum = SubRunNumber;
    CurrentPartNum = PartNumber;
  }
  else if (RunNumber != CurrentRunNum){ //New run has been encountered
    Log("TriggerDataDecoder Tool: New run encountered.  Clearing event building maps",v_message,verbosity); 
    TimeToTriggerWordMap->clear();
  }
  else if (SubRunNumber != CurrentSubrunNum){ //New run has been encountered
    Log("TriggerDataDecoder Tool: New subrun encountered.",v_message,verbosity); 
    TimeToTriggerWordMap->clear();
  }
  else if (PartNumber != CurrentPartNum){ // New part file has been encountered
   if (storetrigoverlap){
     BoostStore StoreTrigOverlap;
     std::stringstream ss_trig_overlap;
     ss_trig_overlap << "TrigOverlap_R"<<RunNumber<<"S"<<SubRunNumber<<"p"<<PartNumber;
     bool store_exist = StoreTrigOverlap.Initialise(ss_trig_overlap.str().c_str());
     StoreTrigOverlap.Set("c1",c1);
     StoreTrigOverlap.Set("c2",c2);
     StoreTrigOverlap.Save(ss_trig_overlap.str().c_str());
   }
   CurrentPartNum = PartNumber;
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
