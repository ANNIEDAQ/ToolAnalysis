#include "MonitorReceive.h"

MonitorReceive::MonitorReceive():Tool(){}


bool MonitorReceive::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  std::string outpath="";
  m_variables.Get("OutPath",outpath);
  m_data->CStore.Set("OutPath",outpath);

  MonitorReceiver= new zmq::socket_t(*m_data->context, ZMQ_SUB);
  MonitorReceiver->setsockopt(ZMQ_SUBSCRIBE, "", 0);
 
  items[0].socket = *MonitorReceiver;
  items[0].fd = 0;
  items[0].events = ZMQ_POLLIN;
  items[0].revents =0;
  

  sources=UpdateMonitorSources();
 
  last= boost::posix_time::ptime(boost::posix_time::second_clock::local_time());
  period =boost::posix_time::time_duration(0,0,1,0);

  m_data->Stores["CCData"]=new BoostStore(false,0);
  m_data->Stores["PMTData"]=new BoostStore(false,0);
  m_data->Stores["TrigData"]=new BoostStore(false,0);
	
  indata=0;
  MRDData=0; 
  PMTData=0;
  TrigData=0;

  return true;
}


bool MonitorReceive::Execute(){


  boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
  boost::posix_time::time_duration duration(current - last);
  if(duration>period){
    last=current;
    sources=UpdateMonitorSources();
  }

  m_data->CStore.Set("HasCCData",false);
  m_data->CStore.Set("HasPMTData",false);
  m_data->CStore.Set("HasTrigData",false);
  m_data->CStore.Set("HasNewFile",false);
  m_data->CStore.Set("Above100",false);

  std::string State="Wait";
  m_data->CStore.Set("State",State);
 
  if(sources>0){
    zmq::poll(&items[0], 1, 100);
    
    if ((items [0].revents & ZMQ_POLLIN)) {
      //   std::cout<<"in poll in"<<std::endl;
      zmq::message_t command;
      MonitorReceiver->recv(&command);
      
      std::istringstream tmp(static_cast<char*>(command.data()));
      m_data->CStore.Set("State",tmp.str());
      
      if (tmp.str()=="MRDSingle"){
	//	std::cout<<"received single"<<std::endl;	
	MRDOut data;
	data.Receive(MonitorReceiver);
	m_data->Stores["CCData"]->Set("Single",data);      
	m_data->CStore.Set("HasCCData",true);
      }
      
      else if(tmp.str()=="DataFile"){
	
	zmq::message_t filepath;
	MonitorReceiver->recv(&filepath);

	std::istringstream iss(static_cast<char*>(filepath.data()));
	//Check if we already processed this file
	if (std::find(loaded_files.begin(),loaded_files.end(),iss.str()) != loaded_files.end()){
		//std::cout <<"MonitorReceive: File "<<iss.str()<<" already loaded before. Abort execution of MonitorReceive tool."<<std::endl;
		return true;
	}
	//Add current file to the list of already processed filenames
	loaded_files.push_back(iss.str());
	//Delete any filenames that are older than the last 100
	while (loaded_files.size()>100){
	  loaded_files.erase(loaded_files.begin());	
	}


	//std::cout<<"received data file="<<iss.str()<<std::endl;

	// Don't delete MRDData & PMTData in case they are stored in a BoostStore, otherwise there will be attempts to delete them twice!
	      
	      
	if (MRDData!=0){
	  m_data->Stores["CCData"]->Delete();
	}
	      
	if (PMTData!=0){
	  m_data->Stores["PMTData"]->Delete();
	  PMTData->Close();
	  PMTData->Delete();
	  delete PMTData;
	  PMTData=0;
	}

	if (TrigData!=0){
	  m_data->Stores["TrigData"]->Delete();
	  TrigData->Close();
	  TrigData->Delete();
	  delete TrigData;
	  TrigData=0;
	}
	      
	if(indata!=0){
	  indata->Close();
	  indata->Delete();
	  delete indata;
	  indata=0;
	}
	
	//Check the size of the current file
	uintmax_t current_filesize = boost::filesystem::file_size(iss.str().c_str());
	m_data->CStore.Set("CurrentFileSize",current_filesize);
	      
	double file_size_mb = current_filesize/1048576.;
	if (file_size_mb > 2000.) {
		std::cout <<"MonitorReceive: Received new file: "<<iss.str()<<std::endl;
		m_data->CStore.Set("HasNewFile",true);
		m_data->CStore.Set("CurrentFileName",iss.str());

		std::time_t current_filetime = boost::filesystem::last_write_time(iss.str().c_str());
		m_data->CStore.Set("CurrentFileTime",current_filetime);
		
		m_data->CStore.Set("Above100",true);

		return true;
        }
	      
	m_data->CStore.Set("Above100",false);
	      
	indata=new BoostStore(false,0); 
	indata->Initialise(iss.str());
	      
	std::cout <<"MonitorReceive: Received new file: "<<iss.str()<<std::endl;
	m_data->CStore.Set("HasNewFile",true);
	m_data->CStore.Set("CurrentFileName",iss.str());



	std::time_t current_filetime = boost::filesystem::last_write_time(iss.str().c_str());
	m_data->CStore.Set("CurrentFileTime",current_filetime);

	MRDData= new BoostStore(false,2);
	PMTData= new BoostStore(false,2);
	TrigData = new BoostStore(false,2);

	if (indata->Has("CCData")){
		try{
		  m_data->CStore.Set("HasCCData",true);	
		  indata->Get("CCData",*MRDData);
		  MRDData->Save("tmp");
		  m_data->Stores["CCData"]->Set("FileData",MRDData,false);
		} catch (...) {
		  Log("MonitorReceive: Did not find CCData in file! (Maybe corrupted!!!) Don't process CCData",0,0);
                  m_data->CStore.Set("HasCCData",false);
		}
	} else {
		m_data->CStore.Set("HasCCData",false);
	}
	if (indata->Has("TrigData")){
		try{
		  m_data->CStore.Set("HasTrigData",true);
		  indata->Get("TrigData",*TrigData);
		  long totalentries_trig;
		  TrigData->Header->Get("TotalEntries",totalentries_trig);
		  std::map<int,TriggerData> TrigData_Map;
		  for (int i_trig=0; i_trig < totalentries_trig; i_trig++){
			  TriggerData TData;
			  TrigData->GetEntry(i_trig);
			  TrigData->Get("TrigData",TData);
			  TrigData_Map.emplace(i_trig,TData);
		  }
		  m_data->Stores["TrigData"]->Set("TrigDataMap",TrigData_Map);
		} catch (...) {
		  Log("MonitorReceive: Did not find TrigData in file! (Maybe corrupted!!!) Don't process TrigData",0,0);
                  m_data->CStore.Set("HasTrigData",false);
		}
	} else {
		m_data->CStore.Set("HasTrigData",false);
	}
	if (indata->Has("PMTData")){
		try{
		  m_data->CStore.Set("HasPMTData",true);
		  indata->Get("PMTData",*PMTData);
		  long totalentries;
        	  PMTData->Header->Get("TotalEntries",totalentries);
        	  std::cout <<"MonitorReceive: Total entries: "<<totalentries<<std::endl;
        	  int ExecuteEntryNum=0;
        	  int EntriesToDo,CDEntryNum;
        	  if (totalentries < 14000) EntriesToDo = 70;      //don't process as many waveforms for AmBe runs (typically ~ 1000 entries)
        	  else EntriesToDo = 1000;               //otherwise do ~1000 entries out of ~15000 (or more)
		  CDEntryNum = 0;
        	  std::map<int,std::vector<CardData>> CardData_Map;
        	  while ((ExecuteEntryNum < EntriesToDo) && (CDEntryNum < totalentries)){
            		  std::vector<CardData> vector_CardData;
            		  PMTData->GetEntry(CDEntryNum);
            		  PMTData->Get("CardData",vector_CardData);
            		  CardData_Map.emplace(CDEntryNum,vector_CardData);
            		  ExecuteEntryNum++;
            		  CDEntryNum++;
        	  }
        	  m_data->Stores["PMTData"]->Set("CardDataMap",CardData_Map);  
		} catch (...) {
		  Log("MonitorReceive: Did not find PMTData in file! (Maybe corrupted!!!) Don't process PMTData",0,0);
                  m_data->CStore.Set("HasPMTData",false);
		}
	} else {
		m_data->CStore.Set("HasPMTData",false);
	}
      }
     }     
  }
  else usleep(100000);

  return true;
}


bool MonitorReceive::Finalise(){

  delete MonitorReceiver;
  MonitorReceiver=0;

  for ( std::map<std::string,Store*>::iterator it=connections.begin(); it!=connections.end(); ++it){
    delete it->second;
    it->second=0;
  }

  connections.clear();
	
  if (indata != 0) {indata->Close(); indata->Delete(); delete indata; indata = 0;}
  if (PMTData != 0) {PMTData->Close(); PMTData->Delete(); delete PMTData; PMTData = 0;}
  if (TrigData != 0) {TrigData->Close(); TrigData->Delete(); delete TrigData; TrigData = 0;}

  m_data->CStore.Remove("State");

  m_data->Stores["CCData"]->Remove("FileData");
  m_data->Stores["CCData"]->Remove("Single");
  m_data->Stores["PMTData"]->Remove("CardDataMap");
  m_data->Stores["TrigData"]->Remove("TrigDataMap");
  // m_data->Stores["PMTData"]->Remove("FileData");
	
  m_data->Stores["CCData"]->Close(); m_data->Stores["CCData"]->Delete();
  m_data->Stores["PMTData"]->Close(); m_data->Stores["PMTData"]->Delete();
  m_data->Stores["TrigData"]->Close(); m_data->Stores["TrigData"]->Delete();
  m_data->Stores.clear();

  return true;
}


int MonitorReceive::UpdateMonitorSources(){

  //std::cout<<"updating monitor sources"<<std::endl;
  boost::uuids::uuid m_UUID=boost::uuids::random_generator()();
  long msg_id=0;
  
  zmq::socket_t Ireceive (*m_data->context, ZMQ_DEALER);
  Ireceive.connect("inproc://ServiceDiscovery");
  
  
  zmq::message_t send(4);
  snprintf ((char *) send.data(), 4 , "%s" ,"All") ;
  
  
  Ireceive.send(send);
  
  
  zmq::message_t receive;
  Ireceive.recv(&receive);
  std::istringstream iss(static_cast<char*>(receive.data()));
  
  int size;
  iss>>size;
  
  for(int i=0;i<size;i++){
    
    Store *service = new Store;
    
    zmq::message_t servicem;
    Ireceive.recv(&servicem);
    
    std::istringstream ss(static_cast<char*>(servicem.data()));
    service->JsonParser(ss.str());
    
    std::string type;
    std::string uuid;
    service->Get("msg_value",type);
    service->Get("uuid",uuid);    
    
    if(type == "MonitorData" && connections.count(uuid)==0){
      connections[uuid]=service;
      std::string ip;
      std::string port;
      service->Get("ip",ip);
      service->Get("remote_port",port);
      std::string tmp="tcp://"+ ip +":"+port;
      MonitorReceiver->connect(tmp.c_str());
      //      MonitorReceiver->setsockopt(ZMQ_SUBSCRIBE, "", 0);
      //std::cout<<type<<" = "<<tmp<<std::endl;
    }  
    else{

      delete service;
      service=0;
    }
    
    
  }
  
  return connections.size();
}
