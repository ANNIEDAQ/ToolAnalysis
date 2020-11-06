#include "LoadRawData.h"

LoadRawData::LoadRawData():Tool(){}


bool LoadRawData::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  verbosity = 0;
  DummyRunInfo = false;

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("BuildType",BuildType);
  m_variables.Get("Mode",Mode);
  m_variables.Get("InputFile",InputFile);
  m_variables.Get("DummyRunInfo",DummyRunInfo);

  m_data= &data; //assigning transient data pointer
  
  if(Mode=="FileList"){
    if(verbosity>v_warning){
      std::cout << "LoadRawData tool: Running in file list mode. " <<
          "Files must be from same run and have sequential parts." << std::endl;
      std::cout << "LoadRawData tool: Organizing file list by part." << std::endl;
    }
    OrganizedFileList = this->OrganizeRunParts(InputFile);
    Log("LoadRawData tool: files to load have been organized.",v_message,verbosity);
  }

  //RawDataObjects
  RawData = new BoostStore(false,0);
  PMTData = new BoostStore(false,2);
  MRDData = new BoostStore(false,2);
  TrigData = new BoostStore(false,2);

  //RawDataEntryObjects
  Cdata = new std::vector<CardData>;
  Tdata = new TriggerData;
  Mdata = new MRDOut;


  FileNum = 0;
  TankEntryNum = 0;
  MRDEntryNum = 0;
  TrigEntryNum = 0;
  FileCompleted = false;
  TankEntriesCompleted = false;
  MRDEntriesCompleted = false;
  TrigEntriesCompleted = false;

  m_data->CStore.Set("FileProcessingComplete",false);
  return true;
}


bool LoadRawData::Execute(){
  m_data->CStore.Set("NewRawDataEntryAccessed",false);
  m_data->CStore.Set("NewRawDataFileAccessed",false);

  //Check if we've reached the end of our file list or single file 
  bool ProcessingComplete = false;
  if(FileCompleted) ProcessingComplete = this->InitializeNewFile();
  if(ProcessingComplete) {
    Log("LoadRawData Tool: All files have been processed.",v_message,verbosity);
    m_data->CStore.Set("FileProcessingComplete",true);
    return true;
  }

  m_data->CStore.Get("PauseTankDecoding",TankPaused);
  m_data->CStore.Get("PauseMRDDecoding",MRDPaused);
  m_data->CStore.Get("PauseCTCDecoding",CTCPaused);

  // Load RawData BoostStore to use in execute loop
  // Points RawDataObjects to current file
  if(Mode=="SingleFile"){
    if(FileNum>0){
      if(verbosity>v_message) std::cout << "LoadRawData tool: Single file has already been loaded" << std::endl;
      return true;
    } else if (TankEntryNum==0 && MRDEntryNum == 0 && TrigEntryNum == 0){
      Log("LoadRawData Tool: Loading Raw Data file as BoostStore",v_message,verbosity); 
      RawData->Initialise(InputFile.c_str());
      if(verbosity>4) RawData->Print(false);
      this->LoadPMTMRDData();
      this->LoadTriggerData();
      this->LoadRunInformation();
    } else {
      Log("LoadRawData Tool: Continuing Raw Data file processing",v_message,verbosity); 
    }
  } 
  
  else if (Mode=="FileList"){
    if(OrganizedFileList.size()==0){
      std::cout << "LoadRawData tool ERROR: no files in file list to parse!" << std::endl;
      return false;
    }
    if(FileCompleted || CurrentFile=="NONE"){
      Log("LoadRawData tool:   Moving to next file.",v_message,verbosity);
      if(verbosity>v_warning) std::cout << "LoadRawData tool: Next file to load: "+OrganizedFileList.at(FileNum) << std::endl;
      CurrentFile = OrganizedFileList.at(FileNum);
      Log("LoadRawData Tool: LoadingRaw Data file as BoostStore",v_debug,verbosity); 
      RawData->Initialise(CurrentFile.c_str());
      m_data->CStore.Set("NewRawDataFileAccessed",true);
      if(verbosity>4) RawData->Print(false);
      this->LoadPMTMRDData();
      this->LoadTriggerData();
      this->LoadRunInformation();
    } else {
     if(verbosity>v_message) std::cout << "LoadRawDataTool: continuing file " << OrganizedFileList.at(FileNum) << std::endl;
    }
    FileCompleted = false;
  } 
  
  else if (Mode == "Processing"){
    std::string State;
    m_data->CStore.Get("State",State);
    Log("LoadRawData tool: checking CStore for status of data stream",v_debug,verbosity);
    if (State == "PMTSingle" || State == "Wait"){
      //Single event file available for monitoring; not relevant for this tool
      if (verbosity > v_message) std::cout <<"LoadRawData: State is "<<State<< ". No new full data file available" << std::endl;
      return true; 
    } 
    else if (State == "DataFile"){
      if (verbosity > v_warning) std::cout<<"LoadRawData: New raw data file available."<<std::endl;
      m_data->Stores["PMTData"]->Get("FileData",RawData);
      if(verbosity>4) RawData->Print(false);
      this->LoadPMTMRDData();
      this->LoadTriggerData();
      this->LoadRunInformation();
    } 
    else {
      Log("LoadRawData Tool: The State >>> "+State+" <<< was not recognized. Please make sure you execute the MonitorReceive tool before the LoadRawData tool when operating in continuous mode",v_error,verbosity);
      return true;   
    }
  }

  else if (Mode == "Monitoring"){
    Log("LoadRawData tool: Monitoring mode not implemented!",v_warning,verbosity);
    return false;
  }

  //Print the current progress of event building in this file
  if(verbosity > v_message){
    std::cout << "LoadRawData Tool: Current progress in file processing: TankEntryNum = "<< TankEntryNum <<", fraction = "<<((double)TankEntryNum/(double)tanktotalentries)*100 << std::endl;
    std::cout << "LoadRawData Tool: Current progress in file processing: MRDEntryNum = "<< MRDEntryNum <<", fraction = "<<((double)MRDEntryNum/(double)mrdtotalentries)*100 << std::endl;
    std::cout << "LoadRawData Tool: Current progress in file processing: TrigEntryNum = "<< TrigEntryNum <<", fraction = "<<((double)TrigEntryNum/(double)trigtotalentries)*100 << std::endl;
    if(verbosity>10){
      // print human readable timestamps from the 3 most recently read entries
      // note that the timestamps aren't easily accessible from the raw data, so we'll lag one loop behind
      // and get the most recent timestamp from the respective decoding tools
      std::cout<<"Getting most recent timestamps"<<std::endl;
      bool get_ok;
      // PMT:
      std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > >* InProgressTankEvents=nullptr;
      get_ok = m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);
      TimeClass tt;
      if(get_ok && InProgressTankEvents){
        if(InProgressTankEvents->size()){
          uint64_t PMTCounterTimeNs = InProgressTankEvents->rbegin()->first;
          tt = TimeClass(PMTCounterTimeNs);
        }
      }
      // TrigData:
      std::map<uint64_t,uint32_t>* TimeToTriggerWordMap=nullptr;
      get_ok = m_data->CStore.Get("TimeToTriggerWordMap",TimeToTriggerWordMap);
      TimeClass cc;
      if(get_ok && TimeToTriggerWordMap){
        if(TimeToTriggerWordMap->size()){
          uint64_t CTCTimeStamp = TimeToTriggerWordMap->rbegin()->first;
          cc = TimeClass(CTCTimeStamp);
        }
      }
      // MRDData:
      std::map<uint64_t, std::string>  MRDTriggerTypeMap;
      get_ok = m_data->CStore.Get("MRDEventTriggerTypes",MRDTriggerTypeMap);
      TimeClass mm;
      if(get_ok){
        if(MRDTriggerTypeMap.size()){
          uint64_t MRDTimeStamp = MRDTriggerTypeMap.rbegin()->first;
          mm = TimeClass (MRDTimeStamp);
        }
      }
      std::cout<<"Lates times are:\nTANK: "<<tt.AsString()<<"\nMRD:  "<<mm.AsString()
               <<"\nCTC:  "<<cc.AsString()<<std::endl;
    }
  }

  // Unpause all other streams when any stream has read all of its events
  if(TankEntryNum == tanktotalentries){
    Log("LoadRawData Tool: ALL PMT ENTRIES COLLECTED.",v_debug, verbosity);
    TankEntriesCompleted = true;
    TankPaused = true;
    if(MRDEntryNum < mrdtotalentries) MRDPaused = false;
    if(TrigEntryNum < trigtotalentries) CTCPaused = false;
  }
  if(MRDEntryNum == mrdtotalentries){
    Log("LoadRawData Tool: ALL MRD ENTRIES COLLECTED.",v_debug, verbosity);
    MRDEntriesCompleted = true;
    MRDPaused = true;
    if(TankEntryNum < tanktotalentries) TankPaused = false;
    if(TrigEntryNum < trigtotalentries) CTCPaused = false;
  }
  if(TrigEntryNum == trigtotalentries){
    Log("LoadRawData Tool: ALL TRIG ENTRIES COLLECTED.",v_debug, verbosity);
    TrigEntriesCompleted = true;
    CTCPaused = true;
    if(TankEntryNum < tanktotalentries) TankPaused = false;
    if(MRDEntryNum < mrdtotalentries) MRDPaused = false;
  }

  m_data->CStore.Set("PauseTankDecoding",TankPaused);
  m_data->CStore.Set("PauseMRDDecoding",MRDPaused);
  m_data->CStore.Set("PauseCTCDecoding",CTCPaused);

  //Get next data entries; are saved to CStore for tools downstream
  this->GetNextDataEntries();

  //Set if the raw data file has been completed
  if (TankEntriesCompleted && BuildType == "Tank") FileCompleted = true;
  if (MRDEntriesCompleted && BuildType == "MRD") FileCompleted = true;
  if ((TankEntriesCompleted && MRDEntriesCompleted) && (BuildType == "TankAndMRD")) FileCompleted = true;
  if ((TrigEntriesCompleted && TankEntriesCompleted && MRDEntriesCompleted) && (BuildType == "TankAndMRDAndCTC")) FileCompleted = true;
    
  m_data->CStore.Set("NewRawDataEntryAccessed",true);
  m_data->CStore.Set("FileCompleted",FileCompleted);
  Log("LoadRawData tool: execution loop complete.",v_debug,verbosity);
  return true;
}


bool LoadRawData::Finalise(){
  RawData->Close();
  RawData->Delete();
  delete RawData;
  PMTData->Close();
  PMTData->Delete();
  delete PMTData;
  MRDData->Close();
  MRDData->Delete();
  delete MRDData;
  std::cout << "LoadRawData Tool Exitting" << std::endl;
  return true;
}

void LoadRawData::LoadRunInformation(){
  Log("LoadRawData Tool: Accessing run information data",v_message,verbosity); 
  Store Postgress;
  if(DummyRunInfo){
    Postgress.Set("RunNumber",-1);
    Postgress.Set("SubRunNumber",-1);
    Postgress.Set("RunType",-1);
    Postgress.Set("StarTime",-1);
  } else{
    BoostStore RunInfo(false,0);
    RawData->Get("RunInformation",RunInfo);
    if(verbosity>3) RunInfo.Print(false);
    RunInfo.Get("Postgress",Postgress);
    if(verbosity>3) Postgress.Print();
  }
  m_data->CStore.Set("RunInfoPostgress",Postgress);
}

void LoadRawData::LoadPMTMRDData(){
  if((BuildType == "TankAndMRD") || (BuildType == "Tank") || (BuildType == "TankAndMRDAndCTC")){
    Log("LoadRawData Tool: Accessing PMT Data in raw data",v_message,verbosity);
    RawData->Get("PMTData",*PMTData);
    PMTData->Header->Get("TotalEntries",tanktotalentries);
    Log("LoadRawData Tool: PMTData has "+std::to_string(tanktotalentries)+" entries",v_debug,verbosity);
    if(verbosity>3) PMTData->Print(false);
    if(verbosity>3) PMTData->Header->Print(false);
    Log("LoadRawData Tool: Setting PMTData into CStore",v_debug, verbosity);
  }
  if((BuildType == "TankAndMRD") || (BuildType == "MRD") || (BuildType == "TankAndMRDAndCTC")){
    Log("LoadRawData Tool: Accessing MRD Data in raw data",v_message,verbosity);
    RawData->Get("CCData",*MRDData);
    MRDData->Header->Get("TotalEntries",mrdtotalentries);
    Log("LoadRawData Tool: MRDData has "+std::to_string(mrdtotalentries)+" entries",v_debug,verbosity);
    if(verbosity>3) MRDData->Print(false);
  }
  return;
}

void LoadRawData::LoadTriggerData(){
  Log("LoadRawData Tool: Accessing Trigger Data in raw data",v_message,verbosity);
  RawData->Get("TrigData",*TrigData);
  if(verbosity>3) TrigData->Print(false);
  TrigData->Header->Get("TotalEntries",trigtotalentries);
  Log("LoadRawData Tool: TrigData has "+to_string(trigtotalentries)+" entries",v_message,verbosity);
  return;
}

std::vector<std::string> LoadRawData::OrganizeRunParts(std::string FileList)
{
  std::vector<std::string> OrganizedFiles;
  std::vector<std::string> UnorganizedFileList;
  std::vector<int> RunCodes;
  int ThisRunCode;
  //First, parse the lines and get all files.
  std::string line;
  ifstream myfile(FileList.c_str());
  if (myfile.is_open()){
    if(verbosity>v_debug) std::cout << "Lines in FileList being printed"  << std::endl; //has our stuff;
    while(getline(myfile,line)){
      if(verbosity>vv_debug) std::cout << "Parsing next line " << std::endl;
      if(line.find("#")!=std::string::npos) continue;
      if(verbosity>v_debug) std::cout << line << std::endl; //has our stuff;
      //std::vector<std::string> splitline;
      //boost::split(splitline,line, boost::is_any_of(","), boost::token_compress_on);
      //for(int i = 0;i<splitline.size(); i ++){
        //std::string filename = splitline.at(i);
        std::string filename = line;
	int rawfilerun=-1;
        int rawfilesubrun, rawfilepart;
        int numargs = 0;
	try{
	    std::reverse(line.begin(),line.end());
            //Part number is the first character now.
            size_t pplace = line.find("p");
            std::string part = line.substr(0,pplace);
            std::reverse(part.begin(),part.end());
            rawfilepart = std::stoi(part);
	    	// Get the
	    	char * p = std::strtok(const_cast<char*>(line.c_str()),"RSp");
	    	while(numargs <2){
	    		p = std::strtok(NULL,"RSp");
	    		std::string capturedstring(p);
                if(verbosity>vv_debug) std::cout << "Caotured string is " << capturedstring << std::endl;
	    		std::reverse(capturedstring.begin(),capturedstring.end());
	    		if (numargs == 0) rawfilesubrun = std::stoi(capturedstring.c_str());
	    		else if (numargs == 1) rawfilerun = std::stoi(capturedstring.c_str());
	    		numargs++;
	    	}
	    } catch (int e){
            std::cout <<"unrecognised input file pattern!: "<<filename
	    		<<"File will be ignored,"<<endl;
	    	//return;
	    	rawfilerun=-1;
	    	rawfilesubrun=-1;
	    	rawfilepart = -1;
	    }
        ThisRunCode = rawfilerun*1000000 + ((rawfilesubrun+1)*1000) + rawfilepart;
        if (rawfilerun!=-1){
          if(verbosity>v_warning) std::cout << "MRDDataDecoder Tool: will parse file : " << filename << std::endl;
          UnorganizedFileList.push_back(filename);
          RunCodes.push_back(ThisRunCode);
        } else {
          if(verbosity>v_error) std::cout << "MRDDataDecoder Tool: Problem parsing this file name: " << filename << std::endl;
        }
      //} // End parse filename in data line 
    } // End parsing each line in file
    if(verbosity>v_debug) std::cout << "MRDDataDecoder Tool: Organizing filenames: " << std::endl;
    //Now, organize files based on the part number array
    std::vector < std::pair<int,std::string>> SortingVector;
    for (int i = 0; i < (int) UnorganizedFileList.size(); i++){
      SortingVector.push_back( std::make_pair(RunCodes.at(i),UnorganizedFileList.at(i)));
    }
    std::sort(SortingVector.begin(),SortingVector.end());
    for (int j = 0; j < (int) SortingVector.size();j ++) {
      OrganizedFiles.push_back(SortingVector.at(j).second);
    }
  } else {
    Log("MRDDataDecoder Tool: FIle not found at nput file list path. "
        "no Organized Run List will be returned ",
        v_error, verbosity);
  }
  return OrganizedFiles;
}

bool LoadRawData::InitializeNewFile(){
  bool EndOfProcessing = false;
  FileNum += 1;
  RawData->Close(); RawData->Delete(); delete RawData; RawData = new BoostStore(false,0);
  MRDData->Close(); MRDData->Delete(); delete MRDData; MRDData = new BoostStore(false,2);
  TrigData->Close(); TrigData->Delete(); delete TrigData; TrigData = new BoostStore(false,2);
  PMTData->Close(); PMTData->Delete(); delete PMTData; PMTData = new BoostStore(false,2);

  TankEntryNum = 0;
  MRDEntryNum = 0;
  TrigEntryNum = 0;
  TankEntriesCompleted = false;
  MRDEntriesCompleted = false;
  TrigEntriesCompleted = false;
  m_data->CStore.Set("PauseTankDecoding",false);
  m_data->CStore.Set("PauseMRDDecoding",false);
  m_data->CStore.Set("PauseCTCDecoding",false);

  if(Mode == "SingleFile"){
    Log("LoadRawData Tool: Single file parsed.  Ending toolchain after this loop.",v_message, verbosity);
    m_data->vars.Set("StopLoop",1);
    EndOfProcessing = true;
  }
  if(Mode == "FileList" && FileNum == int(OrganizedFileList.size())){
    Log("LoadRawData Tool: Full file list parsed.  Ending toolchain after this loop.",v_message, verbosity);
    m_data->vars.Set("StopLoop",1);
    EndOfProcessing = true;
  }
  //No need to stop the loop in continous mode
  if (Mode == "Continous"){
    Log("MRDDataDecoder Tool: Full raw file parsed. Waiting until next raw file is available.",v_message,verbosity);
  }

  return EndOfProcessing;
}

void LoadRawData::GetNextDataEntries(){
  //Get next PMTData Entry
  if(BuildType == "Tank" || BuildType == "TankAndMRD" || BuildType == "TankAndMRDAndCTC"){
    if(!TankPaused && !TankEntriesCompleted){
      Log("LoadRawData Tool: Procesing PMTData Entry "+to_string(TankEntryNum)+"/"+to_string(tanktotalentries),v_debug, verbosity);
      PMTData->GetEntry(TankEntryNum);
      Log("LoadRawData Tool: Getting the PMT card data entry",v_debug, verbosity);
      PMTData->Get("CardData",*Cdata);
      Log("LoadRawData Tool: Setting PMT card data entry into CStore",v_debug, verbosity);
      m_data->CStore.Set("CardData",Cdata);
      Log("LoadRawData Tool: Setting Tank Entry Num CStore",v_debug, verbosity);
      m_data->CStore.Set("TankEntryNum",TankEntryNum);
      TankEntryNum+=1;
    }
  }

  //Get next MRDData Entry
  if(BuildType == "MRD" || BuildType == "TankAndMRD" || BuildType == "TankAndMRDAndCTC"){
    if(!MRDPaused && !MRDEntriesCompleted){
      Log("LoadRawData Tool: Procesing CCData Entry "+to_string(MRDEntryNum)+"/"+to_string(mrdtotalentries),v_debug, verbosity);
      MRDData->GetEntry(MRDEntryNum);
      MRDData->Get("Data",*Mdata);
      m_data->CStore.Set("MRDData",Mdata,true);
      MRDEntryNum+=1;
    }
  }

  //Get next TrigData Entry
  if(BuildType == "TankAndMRDAndCTC" && !TrigEntriesCompleted && !CTCPaused){
    Log("LoadRawData Tool: Procesing TrigData Entry "+to_string(TrigEntryNum)+"/"+to_string(trigtotalentries),v_debug, verbosity);
    TrigData->GetEntry(TrigEntryNum);
    TrigData->Get("TrigData",*Tdata);
    m_data->CStore.Set("TrigData",Tdata);
    TrigEntryNum+=1;
  }
  return;
}
