#include "LoadGeometry.h"

LoadGeometry::LoadGeometry():Tool(),adet(nullptr),AnnieGeometry(nullptr),LAPPD_channel_count(0){}


bool LoadGeometry::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Make the ANNIEEvent Store if it doesn't exist
  int recoeventexists = m_data->Stores.count("ANNIEEvent");
  if(recoeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);

  m_variables.Get("verbosity", verbosity);
  m_variables.Get("FACCMRDGeoFile", fFACCMRDGeoFile);
  m_variables.Get("TankPMTGeoFile", fTankPMTGeoFile);
  m_variables.Get("TankPMTGainFile", fTankPMTGainFile);
  m_variables.Get("AuxiliaryChannelFile", fAuxChannelFile);
  m_variables.Get("LAPPDGeoFile", fLAPPDGeoFile);
  m_variables.Get("DetectorGeoFile", fDetectorGeoFile);
  m_variables.Get("LAPPDChannelCount", LAPPD_channel_count);

  //Check files exist
  if(!this->FileExists(fDetectorGeoFile)){
     Log("LoadGeometry Tool: File for Detector Geometry does not exist!",v_error,verbosity);
     if (verbosity > 0) std::cout << "Filepath was... " << fDetectorGeoFile << std::endl;
     return false;
  }
  if(!this->FileExists(fFACCMRDGeoFile)){
    Log("LoadGeometry Tool: File for FACC/MRD Geometry does not exist!",v_error,verbosity);
    if (verbosity > 0) std::cout << "Filepath was... " << fFACCMRDGeoFile << std::endl;
    return false;
  }

  if(!this->FileExists(fLAPPDGeoFile)){
    Log("LoadGeometry Tool: File for the LAPPDs does not exist!",v_error,verbosity);
    if (verbosity > 0) std::cout << "Filepath was... " << fDetectorGeoFile << std::endl;
    return false;
  }

  if(!this->FileExists(fTankPMTGeoFile)){
    Log("LoadGeometry Tool: File for Tank PMT Geometry does not exist!",v_error,verbosity);
    if (verbosity > 0) std::cout << "Filepath was... " << fTankPMTGeoFile << std::endl;
    return false;
  }
  if(!this->FileExists(fTankPMTGainFile)){
    Log("LoadGeometry Tool: File for Tank PMT Gains does not exist!",v_error,verbosity);
    if (verbosity > 0) std::cout << "Filepath was... " << fTankPMTGainFile << std::endl;
    return false;
  }
  
  if(!this->FileExists(fAuxChannelFile)){
    Log("LoadGeometry Tool: File for Auxiliary Channels does not exist!",v_error,verbosity);
    if (verbosity > 0) std::cout << "Filepath was... " << fAuxChannelFile << std::endl;
    return false;
  }

  //Make the map of channel key to crate space info
  MRDCrateSpaceToChannelNumMap = new std::map<std::vector<int>,int>;
  MRDChannelNumToCrateSpaceMap = new std::map<int,std::vector<int>>;
  TankPMTCrateSpaceToChannelNumMap = new std::map<std::vector<int>,int>;
  ChannelNumToTankPMTSPEChargeMap = new std::map<int,double>;
  ChannelNumToTankPMTCrateSpaceMap = new std::map<int,std::vector<int>>;
  AuxCrateSpaceToChannelNumMap = new std::map<std::vector<int>,int>;
  AuxChannelNumToTypeMap = new std::map<int,std::string>;
  AuxChannelNumToCrateSpaceMap = new std::map<int,std::vector<int>>;
  LAPPDCrateSpaceToChannelNumMap = new std::map<std::vector<unsigned int>,int>;

  //Initialize the geometry using the geometry CSV file entries
  this->InitializeGeometry();

  //Load MRD Geometry Detector/Channel Information
  this->LoadFACCMRDDetectors();

  //Load TankPMT Geometry Detector/Channel Information
  this->LoadTankPMTDetectors();

  //Load TankPMT charge to PE conversion 
  this->LoadTankPMTGains();

  //Load auxiliary and spare channels
  this->LoadAuxiliaryChannels();

  //Load LAPPD Geometry Information
  this->LoadLAPPDs();

  m_data->Stores.at("ANNIEEvent")->Header->Set("AnnieGeometry",AnnieGeometry,true);

  m_data->CStore.Set("MRDCrateSpaceToChannelNumMap",MRDCrateSpaceToChannelNumMap);
  m_data->CStore.Set("MRDChannelNumToCrateSpaceMap",MRDChannelNumToCrateSpaceMap);
  m_data->CStore.Set("TankPMTCrateSpaceToChannelNumMap",TankPMTCrateSpaceToChannelNumMap);
  m_data->CStore.Set("ChannelNumToTankPMTCrateSpaceMap",ChannelNumToTankPMTCrateSpaceMap);
  m_data->CStore.Set("ChannelNumToTankPMTSPEChargeMap",ChannelNumToTankPMTSPEChargeMap);
  m_data->CStore.Set("AuxCrateSpaceToChannelNumMap",AuxCrateSpaceToChannelNumMap);
  m_data->CStore.Set("AuxChannelNumToCrateSpaceMap",AuxChannelNumToCrateSpaceMap);
  m_data->CStore.Set("AuxChannelNumToTypeMap",AuxChannelNumToTypeMap);
  m_data->CStore.Set("LAPPDCrateSpaceToChannelNumMap",LAPPDCrateSpaceToChannelNumMap);
   //AnnieGeometry->GetChannel(0); // trigger InitChannelMap

  return true;
}


bool LoadGeometry::Execute(){
  return true;
}


bool LoadGeometry::Finalise(){
  std::cout << "LoadGeometry tool exitting" << std::endl;
  return true;
}


void LoadGeometry::InitializeGeometry(){
  Log("LoadGeometry tool: Now loading DetectorGeoFile and initializing geometry",v_message,verbosity);
  //Get the Detector file data key
  std::string DetectorLegend = this->GetLegendLine(fDetectorGeoFile);
  std::vector<std::string> DetectorLegendEntries;
  boost::split(DetectorLegendEntries,DetectorLegend, boost::is_any_of(","), boost::token_compress_on);

  //Initialize at zero; will be set later after channels are loaded
  int numtankpmts = 0;
  int numlappds = 0;
  int nummrdpmts = 0;
  int numvetopmts = 0;

  //Initialize data that will be fed to Geometry (units in meters)
  int geometry_version = 0;
  double tank_xcenter = 0.0,tank_ycenter = 0.0,tank_zcenter = 0.0;
  double tank_radius = 0.0,tank_halfheight = 0.0, pmt_enclosed_radius = 0.0, pmt_enclosed_halfheight = 0.0;
  double mrd_width = 0.0,mrd_height = 0.0,mrd_depth = 0.0,mrd_start = 0.0;

  std::string line = "default";
  ifstream myfile(fDetectorGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#") != std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector data (should only be one line here)
    while(getline(myfile,line)){
      if(verbosity>3) std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> DataEntries;
      boost::split(DataEntries,line, boost::is_any_of(","), boost::token_compress_on);
      for (unsigned int i=0; i<DataEntries.size(); i++){
        //Check Legend at i, load correct data type
        int ivalue = 0;
        double dvalue = 0.0;
        if(DetectorLegendEntries.at(i) == "geometry_version") ivalue = std::stoi(DataEntries.at(i));
        else dvalue = std::stod(DataEntries.at(i));
        if (DetectorLegendEntries.at(i) == "geometry_version") geometry_version = ivalue;
        if (DetectorLegendEntries.at(i) == "tank_xcenter") tank_xcenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_ycenter") tank_ycenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_zcenter") tank_zcenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_radius") tank_radius = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_halfheight") tank_halfheight = dvalue;
        if (DetectorLegendEntries.at(i) == "pmt_enclosed_radius") pmt_enclosed_radius = dvalue;
        if (DetectorLegendEntries.at(i) == "pmt_enclosed_halfheight") pmt_enclosed_halfheight = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_width") mrd_width = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_height") mrd_height = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_depth") mrd_depth = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_start") mrd_start = dvalue;
      }
    }
    Position tank_center(tank_xcenter, tank_ycenter, tank_zcenter);
    // Initialize the Geometry
    AnnieGeometry = new Geometry(geometry_version,
                                 tank_center,
                                 tank_radius,
                                 tank_halfheight,
                                 pmt_enclosed_radius,
                                 pmt_enclosed_halfheight,
                                 mrd_width,
                                 mrd_height,
                                 mrd_depth,
                                 mrd_start,
                                 numtankpmts,
                                 nummrdpmts,
                                 numvetopmts,
                                 numlappds,
                                 geostatus::FULLY_OPERATIONAL);
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  myfile.close();
}

void LoadGeometry::LoadFACCMRDDetectors(){
  //First, get the MRD file data key
  Log("LoadGeometry tool: Now loading FACC/MRD detectors",v_message,verbosity);
  std::string MRDLegend = this->GetLegendLine(fFACCMRDGeoFile);
  std::vector<std::string> MRDLegendEntries;
  boost::split(MRDLegendEntries,MRDLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line;
  ifstream myfile(fFACCMRDGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    while(getline(myfile,line)){
      if(verbosity > 4) std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseMRDDataEntry(SpecLine,MRDLegendEntries);
      if(not add_ok){
        std::cerr<<"Faild to add Detector to Geometry!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: FACC/MRD Detector/Channel loading complete",v_message,verbosity);
}

bool LoadGeometry::ParseMRDDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> MRDLegendEntries){
  //Parse the line for information needed to fill the detector & channel classes
  int detector_num = 0,channel_num = 0,detector_system = 0,orientation = 0,layer = 0,side = 0,num = 0,
      rack = 0,TDC_slot = 0,TDC_channel = 0,discrim_slot = 0,discrim_ch = 0,
      patch_panel_row = 0,patch_panel_col = 0,amp_slot = 0,amp_channel = 0,
      hv_crate = 0,hv_slot = 0,hv_channel = 0,nominal_HV,polarity = 0;
  double x_center = 0.0,y_center = 0.0,z_center = 0.0,x_width = 0.0,y_width = 0.0,z_width = 0.0;
  std::string PMT_type = "default",cable_label = "default",paddle_label = "default";

  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing data line into variables",v_debug,verbosity);
  for (unsigned int i=0; i<SpecLine.size(); i++){
    int ivalue = 0;
    double dvalue = 0.0;
    std::string svalue = "default";
    for (unsigned int j=0; j<MRDIntegerValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (unsigned int j=0; j<MRDStringValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }
    for (unsigned int j=0; j<MRDDoubleValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDDoubleValues.at(j)){
        dvalue = std::stod(SpecLine.at(i));
        break;
      }
    }
    //Integers
    if (MRDLegendEntries.at(i) == "detector_num") detector_num = ivalue;
    if (MRDLegendEntries.at(i) == "channel_num") channel_num = ivalue;
    if (MRDLegendEntries.at(i) == "detector_system") detector_system = ivalue;
    if (MRDLegendEntries.at(i) == "orientation") orientation = ivalue;
    if (MRDLegendEntries.at(i) == "layer") layer = ivalue;
    if (MRDLegendEntries.at(i) == "side") side = ivalue;
    if (MRDLegendEntries.at(i) == "num") num = ivalue;
    if (MRDLegendEntries.at(i) == "rack") rack = ivalue;
    if (MRDLegendEntries.at(i) == "TDC_slot") TDC_slot = ivalue;
    if (MRDLegendEntries.at(i) == "TDC_channel") TDC_channel = ivalue;
    if (MRDLegendEntries.at(i) == "discrim_slot") discrim_slot = ivalue;
    if (MRDLegendEntries.at(i) == "discrim_ch") discrim_ch = ivalue;
    if (MRDLegendEntries.at(i) == "patch_panel_row") patch_panel_row = ivalue;
    if (MRDLegendEntries.at(i) == "patch_panel_col") patch_panel_col = ivalue;
    if (MRDLegendEntries.at(i) == "amp_slot") amp_slot = ivalue;
    if (MRDLegendEntries.at(i) == "amp_channel") amp_channel = ivalue;
    if (MRDLegendEntries.at(i) == "hv_crate") hv_crate = ivalue;
    if (MRDLegendEntries.at(i) == "hv_slot") hv_slot = ivalue;
    if (MRDLegendEntries.at(i) == "hv_channel") hv_channel = ivalue;
    if (MRDLegendEntries.at(i) == "nominal_HV") nominal_HV = ivalue;
    if (MRDLegendEntries.at(i) == "polarity") polarity = ivalue;
    //Doubles
    if (MRDLegendEntries.at(i) == "x_center") x_center = dvalue;
    if (MRDLegendEntries.at(i) == "y_center") y_center = dvalue;
    if (MRDLegendEntries.at(i) == "z_center") z_center = dvalue;
    if (MRDLegendEntries.at(i) == "x_width") x_width = dvalue;
    if (MRDLegendEntries.at(i) == "y_width") y_width = dvalue;
    if (MRDLegendEntries.at(i) == "z_width") z_width = dvalue;
    //Strings
    if (MRDLegendEntries.at(i) == "PMT_type") PMT_type = svalue;
    if (MRDLegendEntries.at(i) == "paddle_label") paddle_label = svalue;
    if (MRDLegendEntries.at(i) == "cable_label") cable_label = svalue;
  }

  // Parse whether this is an MRD or Veto Paddle
  std::string dettype = "unknown";
  if (detector_system == 0) dettype = "Veto";
  else if (detector_system == 1) dettype = "MRD";
  //FIXME Need the direction of the MRD PMT
  //FIXME: things that are not loaded in with the default det/channel format:
  //  - discrim_slot, discrim_ch
  //  - patch_panel_row, patch_panel_col, amp_slot, amp_channel
  //  - nominal_HV, polarity
  //  - cable_label, paddle_label
  //
  if(verbosity>4) std::cout << "Filling a FACC/MRD data line into Detector/Channel classes" << std::endl;
  Detector adet(detector_num,
                dettype,
                "MRD", //Change to orientation for PaddleDetector class?
                Position( x_center/100.,
                          y_center/100.,
                          z_center/100.),
                Direction(0.,
                          0.,
                          0.),
                PMT_type,
                detectorstatus::ON,
                0.);

  int MRD_x, MRD_y, MRD_z;
  // orientation 0=horizontal, 1=vertical
  MRD_x = (orientation) ? num  : side;
  MRD_y = (orientation) ? side : num;
  // veto layers are both cabled as z=0, with the layers differentiated by x=0, x=1
  // in practice of course, both span the same x, but are offset in z.
  if(layer>0) MRD_z = layer;
  else        MRD_z = side;

  Paddle apad( detector_num,
               MRD_x,
               MRD_y,
               MRD_z,
               orientation,
               Position( x_center/100.,
                         y_center/100.,
                         z_center/100.),
               std::pair<double,double>{x_center/100.-(x_width/200.), x_center/100.+(x_width/200.)},
               std::pair<double,double>{y_center/100.-(y_width/200.), y_center/100.+(y_width/200.)},
               std::pair<double,double>{z_center/100.-(z_width/200.), z_center/100.+(z_width/200.)});
  
  Channel pmtchannel( channel_num,
                      Position(0,0,0.),
                      -1, // stripside
                      -1, // stripnum
                      rack,
                      TDC_slot,
                      TDC_channel,
                      -1,                 // TDC has no level 2 signal handling
                      -1,
                      -1,
                      hv_crate,
                      hv_slot,
                      hv_channel,
                      channelstatus::ON);

  // Add this channel to the geometry
  if(verbosity>4) cout<<"Adding channel "<<channel_num<<" to detector "<<detector_num<<endl;
  adet.AddChannel(pmtchannel);

  // Also add this channel to the electronics map
  std::vector<int> crate_map{rack,TDC_slot,TDC_channel};
  if(MRDCrateSpaceToChannelNumMap->count(crate_map)==0){
    MRDCrateSpaceToChannelNumMap->emplace(crate_map, channel_num);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning an MRD channel_num to a crate space already defined!!! ",v_error, verbosity);
  }
  if(MRDChannelNumToCrateSpaceMap->count(channel_num)==0){
    MRDChannelNumToCrateSpaceMap->emplace(channel_num, crate_map);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning an MRD crate space to a channel number already defined!!! ",v_error, verbosity);
  }

  if(verbosity>5) cout<<"Adding detector to Geometry"<<endl;
  AnnieGeometry->AddDetector(adet);
  if(verbosity>4) cout<<"Adding paddle to Geometry"<<endl;
  AnnieGeometry->SetDetectorPaddle(detector_num, apad);
  return true;
}

void LoadGeometry::LoadAuxiliaryChannels(){
  //First, get the Tank PMT file legend key
  Log("LoadGeometry tool: Now loading Auxiliary channels",v_message,verbosity);
  std::string AuxChannelLegend = this->GetLegendLine(fAuxChannelFile);
  std::vector<std::string> AuxChannelLegendEntries;
  boost::split(AuxChannelLegendEntries,AuxChannelLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line = "default";
  ifstream myfile(fAuxChannelFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    while(getline(myfile,line)){
      if(verbosity > 3)std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseAuxChannelDataEntry(SpecLine,AuxChannelLegendEntries);
      if(not add_ok){
        std::cerr<<"Failed to add Aux Channel to Crate Space/Channel Key Map!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening the Auxiliary Channel File!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: Auxiliary Channel loading complete",v_message,verbosity);
}

bool LoadGeometry::ParseAuxChannelDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> AuxChannelLegendEntries){

  //Parse the line for information needed to fill the Tdetector & channel classes
  int channel_num = 0, signal_crate = 0,signal_slot = 0,signal_channel = 0;
  std::string channel_type = "NA";

  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing Auxiliary channel line into variables",v_debug,verbosity);
  for (unsigned int i=0; i<SpecLine.size(); i++){
    int ivalue = 0;
    std::string svalue = "default";
    for (unsigned int j=0; j<AuxChannelIntegerValues.size(); j++){
      if(AuxChannelLegendEntries.at(i) == AuxChannelIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (unsigned int j=0; j<AuxChannelStringValues.size(); j++){
      if(AuxChannelLegendEntries.at(i) == AuxChannelStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }

    //Integers
    if (AuxChannelLegendEntries.at(i) == "channel_num") channel_num = ivalue;
    if (AuxChannelLegendEntries.at(i) == "signal_crate") signal_crate = ivalue;
    if (AuxChannelLegendEntries.at(i) == "signal_slot") signal_slot = ivalue;
    if (AuxChannelLegendEntries.at(i) == "signal_channel") signal_channel = ivalue;
    //Strings
    if (AuxChannelLegendEntries.at(i) == "channel_type") channel_type = svalue;
  }

  // Also add this channel to the Tank PMT crate space electronics map
  std::vector<int> crate_map{signal_crate,signal_slot,signal_channel};
  if(AuxCrateSpaceToChannelNumMap->count(crate_map)==0){
    AuxCrateSpaceToChannelNumMap->emplace(crate_map, channel_num);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning an Auxiliary Channel channel_num to a crate space already defined!!! ",v_error, verbosity);
    Log("LoadGeometry Tool: ERROR DETAILS: Signal Crate = "+std::to_string(signal_crate)+", Signal Slot = "+std::to_string(signal_slot)+", Signal Channel = "+std::to_string(signal_channel),v_error,verbosity);
  }
  if(AuxChannelNumToCrateSpaceMap->count(channel_num)==0){
    AuxChannelNumToCrateSpaceMap->emplace(channel_num,crate_map);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning a crate space to an Auxiliary Channel already defined!!! ",v_error, verbosity);
    Log("LoadGeometry Tool: ERROR DETAILS: Signal Crate = "+std::to_string(signal_crate)+", Signal Slot = "+std::to_string(signal_slot)+", Signal Channel = "+std::to_string(signal_channel),v_error,verbosity);
  }
  if(AuxChannelNumToTypeMap->count(channel_num)==0){
    AuxChannelNumToTypeMap->emplace(channel_num, channel_type);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning an Auxiliary Channel Type to a channel already defined!!! ",v_error, verbosity);
    Log("LoadGeometry Tool: ERROR DETAILS: Signal Type = "+ channel_type +", channel key = "+std::to_string(channel_num),v_error,verbosity);
  }

  return true;
}


void LoadGeometry::LoadTankPMTDetectors(){
  //First, get the Tank PMT file legend key
  Log("LoadGeometry tool: Now loading TankPMT detectors",v_message,verbosity);
  std::string TankPMTLegend = this->GetLegendLine(fTankPMTGeoFile);
  std::vector<std::string> TankPMTLegendEntries;
  boost::split(TankPMTLegendEntries,TankPMTLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line = "default";
  ifstream myfile(fTankPMTGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    while(getline(myfile,line)){
      if(verbosity > 3)std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseTankPMTDataEntry(SpecLine,TankPMTLegendEntries);
      if(not add_ok){
        std::cerr<<"Failed to add Tank PMT Detector to Geometry!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening the Tank PMT file!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: Tank PMT Detector/Channel loading complete",v_message,verbosity);
}

bool LoadGeometry::ParseTankPMTDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> TankPMTLegendEntries){

  //Parse the line for information needed to fill the Tdetector & channel classes
  int detector_num = 0,channel_num = 0,panel_number = 0,signal_crate = 0,signal_slot = 0,signal_channel = 0,
      mt_crate = 0, mt_slot = 0, mt_channel = 0, hv_crate = 0,hv_slot = 0,hv_channel = 0,nominal_HV = 0,
      sb_num = 0, sb_channel = 0;
  double x_pos = 0.0,y_pos = 0.0,z_pos = 0.0,x_dir = 0.0,y_dir = 0.0,z_dir = 0.0;
  std::string detector_tank_location = "default",PMT_type = "default",cable_label = "default",detector_status = "default";

  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing data line into variables",v_debug,verbosity);
  for (unsigned int i=0; i<SpecLine.size(); i++){
    int ivalue = 0;
    double dvalue = 0.0;
    std::string svalue = "default";
    for (unsigned int j=0; j<TankPMTIntegerValues.size(); j++){
      if(TankPMTLegendEntries.at(i) == TankPMTIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (unsigned int j=0; j<TankPMTStringValues.size(); j++){
      if(TankPMTLegendEntries.at(i) == TankPMTStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }
    for (unsigned int j=0; j<TankPMTDoubleValues.size(); j++){
      if(TankPMTLegendEntries.at(i) == TankPMTDoubleValues.at(j)){
        dvalue = std::stod(SpecLine.at(i));
        break;
      }
    }

    //Integers
    if (TankPMTLegendEntries.at(i) == "detector_num") detector_num = ivalue;
    if (TankPMTLegendEntries.at(i) == "channel_num") channel_num = ivalue;
    if (TankPMTLegendEntries.at(i) == "panel_number") panel_number = ivalue;
    if (TankPMTLegendEntries.at(i) == "sb_num") sb_num = ivalue;
    if (TankPMTLegendEntries.at(i) == "sb_channel") sb_channel = ivalue;
    if (TankPMTLegendEntries.at(i) == "signal_crate") signal_crate = ivalue;
    if (TankPMTLegendEntries.at(i) == "signal_slot") signal_slot = ivalue;
    if (TankPMTLegendEntries.at(i) == "signal_channel") signal_channel = ivalue;
    if (TankPMTLegendEntries.at(i) == "mt_crate") mt_crate = ivalue;
    if (TankPMTLegendEntries.at(i) == "mt_slot") mt_slot = ivalue;
    if (TankPMTLegendEntries.at(i) == "mt_channel") mt_channel = ivalue;
    if (TankPMTLegendEntries.at(i) == "hv_crate") hv_crate = ivalue;
    if (TankPMTLegendEntries.at(i) == "hv_slot") hv_slot = ivalue;
    if (TankPMTLegendEntries.at(i) == "hv_channel") hv_channel = ivalue;
    if (TankPMTLegendEntries.at(i) == "nominal_HV") nominal_HV = ivalue;
    //Doubles
    if (TankPMTLegendEntries.at(i) == "x_pos") x_pos = dvalue;
    if (TankPMTLegendEntries.at(i) == "y_pos") y_pos = dvalue;
    if (TankPMTLegendEntries.at(i) == "z_pos") z_pos = dvalue;
    if (TankPMTLegendEntries.at(i) == "x_dir") x_dir = dvalue;
    if (TankPMTLegendEntries.at(i) == "y_dir") y_dir = dvalue;
    if (TankPMTLegendEntries.at(i) == "z_dir") z_dir = dvalue;
    //Strings
    if (TankPMTLegendEntries.at(i) == "detector_tank_location") detector_tank_location = svalue;
    if (TankPMTLegendEntries.at(i) == "PMT_type") PMT_type = svalue;
    if (TankPMTLegendEntries.at(i) == "cable_label") cable_label = svalue;
    if (TankPMTLegendEntries.at(i) == "detector_status") detector_status = svalue;
  }

  //Parse out the Detector Status for filling into Detector class
  detectorstatus detstatus = detectorstatus::OFF;
  channelstatus chanstatus = channelstatus::OFF;
  if(detector_status == "ON"){
    detstatus = detectorstatus::ON;
    chanstatus = channelstatus::ON;
  }
  else if(detector_status == "OFF"){
    detstatus = detectorstatus::OFF;
    chanstatus = channelstatus::OFF;
  }
  else if(detector_status == "UNSTABLE"){
    detstatus = detectorstatus::UNSTABLE;
    chanstatus = channelstatus::UNSTABLE;
  }
  else {
    Log("LoadGeometry Tool: Undefined status of Tank PMT detector",v_error,verbosity);
    if (verbosity > v_error) std::cout << "channel_num is " << channel_num << std::endl;
  }

  //FIXME: things that are not loaded in with the default det/channel format:
  //      - panel_number

  if(verbosity>4) std::cout << "Filling a Tank PMT data line into Detector/Channel classes" << std::endl;
  Detector adet(detector_num,
                "Tank",
                detector_tank_location,
                Position( x_pos,
                          y_pos,
                          z_pos),
                Direction(x_dir,
                          y_dir,
                          z_dir),
                PMT_type,
                detstatus,
                0.);

  Channel pmtchannel( channel_num,
                      Position(0,0,0.),
                      -1, // stripside
                      -1, // stripnum
                      signal_crate,
                      signal_slot,
                      signal_channel,
                      mt_crate,
                      mt_slot,
                      mt_channel,
                      hv_crate,
                      hv_slot,
                      hv_channel,
                      chanstatus); //channel status same as detector status here

  // Also add this channel to the Tank PMT electronics map
  std::vector<int> crate_map{signal_crate,signal_slot,signal_channel};
  if(TankPMTCrateSpaceToChannelNumMap->count(crate_map)==0){
    TankPMTCrateSpaceToChannelNumMap->emplace(crate_map, channel_num);
    ChannelNumToTankPMTCrateSpaceMap->emplace(channel_num,crate_map);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning a Tank PMT channel_num to a crate space already defined!!! ",v_error, verbosity);
    Log("LoadGeometry Tool: ERROR DETAILS: Signal Crate = "+std::to_string(signal_crate)+", Signal Slot = "+std::to_string(signal_slot)+", Signal Channel = "+std::to_string(signal_channel),v_error,verbosity);
  }

  // Add this channel to the geometry
  if(verbosity>4) cout<<"Adding channel "<<channel_num<<" to detector "<<detector_num<<endl;
  adet.AddChannel(pmtchannel);
  if(verbosity>5) cout<<"Adding detector to Geometry"<<endl;
  AnnieGeometry->AddDetector(adet);
  return true;
}



void LoadGeometry::LoadLAPPDs(){
  //First, get the LAPPD file data key
  Log("LoadGeometry tool: Now loading LAPPDs",v_message,verbosity);
  std::string LAPPDLegend = this->GetLegendLine(fLAPPDGeoFile);
  std::vector<std::string> LAPPDLegendEntries;
  boost::split(LAPPDLegendEntries,LAPPDLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line;
  ifstream myfile(fLAPPDGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    detector_num_store = 100000;
    counter = 0;
    while(getline(myfile,line)){
      if(verbosity>4) std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseLAPPDDataEntry(SpecLine,LAPPDLegendEntries);
      if(not add_ok){
        std::cerr<<"Faild to add Detector to Geometry!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: LAPPD Detector/Channel loading complete",v_message,verbosity);
}


bool LoadGeometry::ParseLAPPDDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> LAPPDLegendEntries){
  //Parse the line for information needed to fill the detector & channel classes
   int detector_num = 0,channel_strip_side = 0,channel_strip_num = 0;
   unsigned int channel_signal_crate = 0,channel_signal_card = 0,channel_signal_channel = 0,channel_level2_crate = 0,channel_level2_card = 0,channel_level2_channel = 0,channel_hv_crate = 0,
   channel_hv_card = 0,channel_hv_channel = 0,channel_num = 0;
   double detector_position_x = 0.0,detector_position_y = 0.0,detector_position_z = 0.0,detector_direction_x = 0.0,detector_direction_y = 0.0,detector_direction_z = 0.0,
   channel_position_x = 0.0,channel_position_y = 0.0,channel_position_z = 0.0;
   std::string detector_type = "default",detector_status = "default",channel_status = "default";
  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing data line into variables",v_debug,verbosity);
  for (unsigned int i=0; i<SpecLine.size(); i++){
    int ivalue = 0;
    unsigned int uivalue =0;
    double dvalue = 0.0;
    std::string svalue = "default";
    for (unsigned int j=0; j<LAPPDIntegerValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (unsigned int j=0; j<LAPPDStringValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }
    for (unsigned int j=0; j<LAPPDDoubleValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDDoubleValues.at(j)){
        dvalue = std::stod(SpecLine.at(i));
        break;
      }
    }
    for (unsigned int j=0; j<LAPPDUnIntValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDUnIntValues.at(j)){
        uivalue = std::stoul(SpecLine.at(i));
        break;
      }
    }
    //Integers
    if (LAPPDLegendEntries.at(i) == "detector_num") detector_num = ivalue;
    if (LAPPDLegendEntries.at(i) == "channel_strip_side") channel_strip_side = ivalue;
    if (LAPPDLegendEntries.at(i) == "channel_strip_num") channel_strip_num = ivalue;

    //Unsigned Integers
    if (LAPPDLegendEntries.at(i) == "channel_signal_crate") channel_signal_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_signal_card") channel_signal_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_signal_channel") channel_signal_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_crate") channel_level2_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_card") channel_level2_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_channel") channel_level2_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_crate") channel_hv_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_card") channel_hv_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_channel") channel_hv_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_num") channel_num = uivalue;

    //Doubles
    if (LAPPDLegendEntries.at(i) == "detector_position_x") detector_position_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_position_y") detector_position_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_position_z") detector_position_z = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_x") detector_direction_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_y") detector_direction_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_z") detector_direction_z = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_x") channel_position_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_y") channel_position_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_z") channel_position_z = dvalue;

    //Strings
    if (LAPPDLegendEntries.at(i) == "detector_type") detector_type = svalue;
    if (LAPPDLegendEntries.at(i) == "detector_status") detector_status = svalue;
    if (LAPPDLegendEntries.at(i) == "channel_status") channel_status = svalue;
  }

  if(verbosity>4) std::cout << "Filling a LAPPD data line into Detector/Channel classes" << std::endl;
  if(detector_num != detector_num_store){
  detectorstatus detstat = detectorstatus::OFF;
  if(detector_status == "OFF"){
    detstat = detectorstatus::OFF;
    }
    else if(detector_status == "ON"){
      detstat = detectorstatus::ON;
    }
    else if(detector_status == "UNSTABLE"){
      detstat = detectorstatus::UNSTABLE;
    }
    else{
      std::cerr << "The chosen detector status isn't available!!!" << std::endl;
    }
  //TODO Somewhere it has to be stated that the units are in [m] for LAPPDs for now
  adet = new Detector(464+detector_num,
                "LAPPD",
                "Barrel",
                Position(detector_position_x,
                        detector_position_y,
                        detector_position_z),
                Direction(detector_direction_x,
                          detector_direction_y,
                          detector_direction_z),
                detector_type,
                detstat,
                0.);
  detector_num_store = detector_num;
  }

  channelstatus channelstat = channelstatus::OFF;
  if(channel_status == "OFF"){
      channelstat = channelstatus::OFF;
      }
  else if(channel_status == "ON"){
      channelstat = channelstatus::ON;
        }
  else if(channel_status == "UNSTABLE"){
      channelstat = channelstatus::UNSTABLE;
      }
  else{
  std::cerr << "The chosen channel status isn't available!!!" << std::endl;
      }
  Channel lappdchannel(464+channel_num,
                      Position(channel_position_x,
                               channel_position_y,
                               channel_position_z),
                      channel_strip_side,
                      channel_strip_num,
                      channel_signal_crate,
                      channel_signal_card,
                      channel_signal_channel,
                      channel_level2_crate,
                      channel_level2_card,
                      channel_level2_channel,
                      channel_hv_crate,
                      channel_hv_card,
                      channel_hv_channel,
                      channelstat);

  // Also add this channel to the Tank PMT electronics map
  std::vector<unsigned int> crate_map{channel_signal_crate,channel_signal_card,channel_signal_channel};
  if(LAPPDCrateSpaceToChannelNumMap->count(crate_map)==0){
    LAPPDCrateSpaceToChannelNumMap->emplace(crate_map, channel_num);
  } else {
    Log("LoadGeometry Tool: ERROR: Tried assigning a Tank PMT channel_num to a crate space already defined!!! ",v_error, verbosity);
  }

  // Add this channel to the detector
  if(adet != nullptr){
  if(verbosity>4) cout<<"Adding channel "<<channel_num<<" to LAPPD "<<detector_num<<endl;
  adet->AddChannel(lappdchannel);
  }
  counter++;
  if(adet != nullptr && counter == LAPPD_channel_count){
  if(verbosity>5) cout<<"Adding LAPPD to Geometry"<<endl;
  AnnieGeometry->AddDetector(*adet);
  counter = 0;
  }
  return true;
}


bool LoadGeometry::FileExists(std::string name) {
  ifstream myfile(name.c_str());
  return myfile.good();
}


std::string LoadGeometry::GetLegendLine(std::string name) {
  if(verbosity>4) std::cout << "Getting legend of file: " << name << std::endl;
  std::string line;
  std::string legendline = "null";
  ifstream myfile(name.c_str());
  if (myfile.is_open()){
    while(std::getline(myfile,line)){
      if(verbosity>4) std::cout << line << std::endl;
      if(line.find("#") != std::string::npos) continue;
      if(line.find(LegendLineLabel) != std::string::npos){
        //Next line is the title line
        getline(myfile,line);
        legendline = line;
        if(verbosity>4) std::cout<<"Legend line loaded. Legend is: " << legendline << std::endl;
        break;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(legendline=="null"){
    Log("LoadGeometry tool: Legend line label not found!!!",v_error,verbosity);
  }
  myfile.close();
  return legendline;
}

void LoadGeometry::LoadTankPMTGains(){
  ifstream myfile(fTankPMTGainFile.c_str());
  std::string line;
  if (myfile.is_open()){
    //Loop over lines, collect all detector data (should only be one line here)
    while(getline(myfile,line)){
      if(verbosity>3) std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      std::vector<std::string> DataEntries;
      boost::split(DataEntries,line, boost::is_any_of(","), boost::token_compress_on);
      int channelkey = -9999;
      double SPECharge = -9999.;
      channelkey = std::stoi(DataEntries.at(0));
      SPECharge= std::stod(DataEntries.at(1));
      ChannelNumToTankPMTSPEChargeMap->emplace(channelkey,SPECharge);
    }
  }
  return;
}
