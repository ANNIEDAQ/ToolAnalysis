#include "CNNImage.h"

CNNImage::CNNImage():Tool(){}


bool CNNImage::Initialise(std::string configfile, DataModel &data){


  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();
  m_data= &data; //assigning transient data pointer

  //---------------------------------------------------------------
  //------------------Read in configuration------------------------
  //---------------------------------------------------------------
 
  //Default values
  includeTopBottom = false;
  dimensionLAPPD = 1;
  isData = 0;

  //User-specified values
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("DetectorConf",detector_config);
  m_variables.Get("DataMode",data_mode);
  m_variables.Get("SaveMode",save_mode);
  m_variables.Get("DimensionX",dimensionX);
  m_variables.Get("DimensionY",dimensionY);
  m_variables.Get("IncludeTopBottom",includeTopBottom);
  m_variables.Get("OutputFile",cnn_outpath);
  m_variables.Get("DimensionLAPPD",dimensionLAPPD);
  m_variables.Get("IsData",isData); 
 
  Log("CNNImage tool: Initialising...",v_message,verbosity);

  //Check for user options
  if (data_mode != "Normal" && data_mode != "Charge-Weighted" && data_mode != "TimeEvolution") data_mode = "Normal";
  if (save_mode != "Geometric" && save_mode != "PMT-wise") save_mode = "Geometric";
  nlappdX = dimensionLAPPD;
  nlappdY = dimensionLAPPD;
  npmtsX = 16;    //in every row, we have 16 PMTs (8 sides, 2 PMTs per side per row)
  if (!includeTopBottom) npmtsY = 6;     //we have 6 rows of PMTs (excluding top and bottom PMTs)
  else npmtsY = 10;		//for top & bottom PMTs, include 2 extra rows of PMTs at top and at the bottom
  
  //Output configuration
  Log("CNNImage tool: Data mode: "+data_mode+" [Normal/Charge-Weighted]",v_message,verbosity);
  Log("CNNImage tool: Save mode: "+save_mode+" [Geometric/PMT-wise]",v_message,verbosity);
  Log("CNNImage tool: LAPPD dimension: "+dimensionLAPPD,v_message,verbosity);

  //---------------------------------------------------------------
  //----------------------Get geometry-----------------------------
  //---------------------------------------------------------------

  //Get general detector properties
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  tank_radius = geom->GetTankRadius();
  tank_height = geom->GetTankHalfheight();
  double barrel_compression = 0.82;
  if (detector_config == "ANNIEp2v6") tank_height*=barrel_compression;
  if (tank_radius < 1.) tank_radius = 1.37504;
  Position detector_center = geom->GetTankCentre();
  tank_center_x = detector_center.X();
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();
  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  if (isData){
    m_data->CStore.Get("pmt_tubeid_to_channelkey",pmtid_to_channelkey);
  }

  //Read in PMT positions
  max_y = -100.;
  min_y = 100.;

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Tank").begin();
                                                    it != Detectors->at("Tank").end();
                                                  ++it){
    Detector* apmt = it->second;
    unsigned long detkey = it->first;
    pmt_detkeys.push_back(detkey);
    unsigned long chankey = apmt->GetChannels()->begin()->first;
    pmt_chankeys.push_back(chankey);
    Position position_PMT = apmt->GetDetectorPosition();
    x_pmt.insert(std::pair<int,double>(detkey,position_PMT.X()-tank_center_x));
    y_pmt.insert(std::pair<int,double>(detkey,position_PMT.Y()-tank_center_y));
    z_pmt.insert(std::pair<int,double>(detkey,position_PMT.Z()-tank_center_z));
    if (y_pmt[detkey]>max_y && apmt->GetTankLocation()!="OD") max_y = y_pmt.at(detkey);
    if (y_pmt[detkey]<min_y && apmt->GetTankLocation()!="OD") min_y = y_pmt.at(detkey);
    
    //Some debug information
    Log("CNNImage tool: Reading in detkey: "+std::to_string(detkey)+", chankey: "+std::to_string(chankey)+", position: ("+std::to_string(position_PMT.X())+","+std::to_string(position_PMT.Y())+","+std::to_string(position_PMT.Z())+")",v_debug,verbosity);
    Log("CNNImage tool: Rho PMT: "+std::to_string(sqrt(x_pmt.at(detkey)*x_pmt.at(detkey)+z_pmt.at(detkey)*z_pmt.at(detkey)))+", y PMT: "+std::to_string(y_pmt.at(detkey)),v_debug,verbosity);
  }

  Log("CNNImage tool: Max/min y-dimensions of PMT volume: Max y = "+std::to_string(max_y)+", min y = "+std::to_string(min_y),v_message,verbosity);

  //---------------------------------------------------------------
  //-------------------Order PMT positions-------------------------
  //---------------------------------------------------------------

  for (unsigned int i_pmt = 0; i_pmt < y_pmt.size(); i_pmt++){

    double x,y;
    unsigned long detkey = pmt_detkeys[i_pmt];
    Position pmt_pos(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]);
    unsigned long chankey = pmt_chankeys[i_pmt];
    Detector *apmt = geom->ChannelToDetector(chankey);
    if (apmt->GetTankLocation()=="OD") continue;  //don't include OD PMTs
    if ((y_pmt[detkey] >= max_y-0.001 || y_pmt[detkey] <= min_y+0.001) && !includeTopBottom) continue;     //don't include top/bottom PMTs if specified
 
    if (y_pmt[detkey] >= max_y-0.001) {
      ConvertPositionTo2D_Top(pmt_pos, x, y);
    }
    else if (y_pmt[detkey] <= min_y+0.001) {
	ConvertPositionTo2D_Bottom(pmt_pos, x, y);
    }
    else {
       ConvertPositionTo2D(pmt_pos, x, y);
    }
    x = (round(100*x)/100.);
    y = (round(100*y)/100.);
    if (fabs(x-0.52)<0.001) x = 0.53;
    if (fabs(x-0.47)<0.001) x = 0.48;
    if (y_pmt[detkey] <= min_y + 0.001){
      if (fabs(x-0.23)<0.001) x = 0.22;
      if (fabs(x-0.77)<0.001) x = 0.78;
    }
    if (y_pmt[detkey] >= max_y - 0.001){
      if (fabs(x-0.19)<0.001) x = 0.21;
      if (fabs(x-0.34)<0.001) x = 0.32;
      if (fabs(x-0.48)<0.001) x = 0.50;
      if (fabs(x-0.66)<0.001) x = 0.68;
    }
    if (y_pmt[detkey] >= max_y-0.001) vec_pmt2D_x_Top.push_back(x);
    else if (y_pmt[detkey] <= min_y+0.001) vec_pmt2D_x_Bottom.push_back(x);
    else vec_pmt2D_x.push_back(x);
    vec_pmt2D_y.push_back(y);
    Log("CNNImage tool: Detkey "+std::to_string(detkey)+", x (2D): "+std::to_string(x)+", y (2D): "+std::to_string(y),v_message,verbosity);
  }

  //Sort the positions (in x & y separately)
  Log("CNNImage tool: vec_pmt2D_* size: "+std::to_string(vec_pmt2D_x.size()),v_message,verbosity);
  std::sort(vec_pmt2D_x.begin(),vec_pmt2D_x.end());
  std::sort(vec_pmt2D_y.begin(),vec_pmt2D_y.end());
  std::sort(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end());
  std::sort(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end());
  vec_pmt2D_x.erase(std::unique(vec_pmt2D_x.begin(),vec_pmt2D_x.end()),vec_pmt2D_x.end());
  vec_pmt2D_y.erase(std::unique(vec_pmt2D_y.begin(),vec_pmt2D_y.end()),vec_pmt2D_y.end());
  vec_pmt2D_x_Top.erase(std::unique(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end()),vec_pmt2D_x_Top.end());
  vec_pmt2D_x_Bottom.erase(std::unique(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end()),vec_pmt2D_x_Bottom.end());
  
  //Output sorted positions
  if (verbosity >= v_message){
    std::cout <<"Sorted 2D position vectors: Barrel x "<<std::endl;
    for (unsigned int i_x=0;i_x<vec_pmt2D_x.size();i_x++){
      std::cout <<i_x<<": "<<vec_pmt2D_x.at(i_x)<<std::endl;
    }
    std::cout <<"Sorted 2D position vectors: Top x: "<<std::endl;
    for (unsigned int i_x=0;i_x<vec_pmt2D_x_Top.size();i_x++){
      std::cout <<i_x<<": "<<vec_pmt2D_x_Top.at(i_x)<<std::endl;
    }
    std::cout <<"Sorted 2D position vectors: Bottom x: "<<std::endl;
    for (unsigned int i_x=0;i_x<vec_pmt2D_x_Bottom.size();i_x++){
      std::cout <<i_x<<": "<<vec_pmt2D_x_Bottom.at(i_x)<<std::endl;
    }
    std::cout <<"Sorted 2D position vectors: y: "<<std::endl;
    for (unsigned int i_y=0;i_y<vec_pmt2D_y.size();i_y++){
      std::cout <<i_y<<": "<<vec_pmt2D_y.at(i_y)<<std::endl;
    }
  }

  //---------------------------------------------------------------
  //-------------------Get LAPPD positions-------------------------
  //---------------------------------------------------------------
 
  if (!isData){
  max_y_lappd = -100.;
  min_y_lappd = 100.;

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("LAPPD").begin();
                                                    it != Detectors->at("LAPPD").end();
                                                  ++it){
    Detector* anlappd = it->second;
    unsigned long detkey = it->first;
    lappd_detkeys.push_back(detkey);
    unsigned long chankey = anlappd->GetChannels()->begin()->first;
    Position position_LAPPD = anlappd->GetDetectorPosition();
    x_lappd.insert(std::pair<int,double>(detkey,position_LAPPD.X()-tank_center_x));
    y_lappd.insert(std::pair<int,double>(detkey,position_LAPPD.Y()-tank_center_y));
    z_lappd.insert(std::pair<int,double>(detkey,position_LAPPD.Z()-tank_center_z));
    if (y_lappd[detkey]>max_y_lappd) max_y_lappd = y_lappd.at(detkey);
    if (y_lappd[detkey]<min_y_lappd) min_y_lappd = y_lappd.at(detkey);

    //Print some debug information about LAPPD positions
    Log("CNNImage tool: Detkey: "+std::to_string(detkey)+", Chankey: "+std::to_string(chankey)+", position: ("+std::to_string(position_LAPPD.X())+","+std::to_string(position_LAPPD.Y())+","+std::to_string(position_LAPPD.Z())+")",v_debug,verbosity);
    Log("CNNImage tool: Rho LAPPD "+std::to_string(sqrt(x_lappd.at(detkey)*x_lappd.at(detkey)+z_lappd.at(detkey)*z_lappd.at(detkey)))+", y LAPPD: "+std::to_string(y_lappd.at(detkey)),v_debug,verbosity);
  }

  //---------------------------------------------------------------
  //-------------------Order LAPPD positions-----------------------
  //---------------------------------------------------------------

  for (unsigned int i_lappd = 0; i_lappd < y_lappd.size(); i_lappd++){

    double x,y;
    unsigned long detkey = lappd_detkeys[i_lappd];
    Position lappd_pos(x_lappd[detkey],y_lappd[detkey],z_lappd[detkey]);
    if (y_lappd[detkey] >= max_y || y_lappd[detkey] <= min_y) continue;     //don't include top/bottom/OD PMTs for now
    ConvertPositionTo2D(lappd_pos, x, y);
    vec_lappd2D_x.push_back(x);
    vec_lappd2D_y.push_back(y);
    Log("CNNImage tool: Detkey "+std::to_string(detkey)+", x (2D): "+std::to_string(x)+", y (2D): "+std::to_string(y),v_message,verbosity);
  }

  //Sort LAPPD positions
  Log("CNNImage tool: vec_lappd2D_* size: "+std::to_string(vec_lappd2D_x.size()),v_message,verbosity);
  std::sort(vec_lappd2D_x.begin(),vec_lappd2D_x.end());
  std::sort(vec_lappd2D_y.begin(),vec_lappd2D_y.end());
  vec_lappd2D_x.erase(std::unique(vec_lappd2D_x.begin(),vec_lappd2D_x.end()),vec_lappd2D_x.end());
  vec_lappd2D_y.erase(std::unique(vec_lappd2D_y.begin(),vec_lappd2D_y.end()),vec_lappd2D_y.end());

  //Output the sorted positions
  if (verbosity >=2){
    std::cout <<"CMMImage tool: Sorted LAPPD 2D position vectors: x-values: "<<std::endl;
    for (unsigned int i_x=0;i_x<vec_lappd2D_x.size();i_x++){
      std::cout <<i_x<<": "<<vec_lappd2D_x.at(i_x)<<std::endl;
    }
    std::cout <<"CMMImage tool: Sorted LAPPD 2D position vectors: y-values: "<<std::endl;
    for (unsigned int i_y=0;i_y<vec_lappd2D_y.size();i_y++){
      std::cout <<i_y<<": "<<vec_lappd2D_y.at(i_y)<<std::endl;
    }
  }
  }

  //---------------------------------------------------------------
  //-------------------Define root+csv output files----------------
  //---------------------------------------------------------------

  std::string str_root = ".root";
  std::string str_csv = ".csv";
  std::string str_time = "_time";
  std::string str_charge = "_charge";
  std::string str_lappd = "_lappd";
  std::string str_Rings = "_Rings";
  std::string str_MRD = "_MRD";
  std::string str_abs = "_abs";
  std::string str_first = "_first";
  std::string rootfile_name = cnn_outpath + str_root;
  std::string csvfile_name = cnn_outpath + str_charge + str_csv;
  std::string csvfile_abs_name = cnn_outpath + str_charge + str_abs + str_csv;
  std::string csvfile_time_name = cnn_outpath + str_time + str_csv;
  std::string csvfile_time_first_name = cnn_outpath + str_time + str_first + str_csv;
  std::string csvfile_time_abs_name = cnn_outpath + str_time + str_abs + str_csv;
  std::string csvfile_time_first_abs_name = cnn_outpath + str_time + str_first + str_abs + str_csv;
  std::string csvfile_lappd_name = cnn_outpath + str_lappd + str_charge + str_csv;
  std::string csvfile_lappd_abs_name = cnn_outpath + str_lappd + str_charge + str_abs + str_csv;
  std::string csvfile_lappd_time_name = cnn_outpath + str_lappd + str_time + str_csv;
  std::string csvfile_lappd_time_first_name = cnn_outpath + str_lappd + str_time + str_first + str_csv;
  std::string csvfile_lappd_time_abs_name = cnn_outpath + str_lappd + str_time + str_abs + str_csv;
  std::string csvfile_lappd_time_first_abs_name = cnn_outpath + str_lappd + str_time + str_first + str_abs + str_csv;
  std::string csvfile_MRD = cnn_outpath + str_MRD + str_csv;
  std::string csvfile_Rings = cnn_outpath + str_Rings + str_csv;

  file = new TFile(rootfile_name.c_str(),"RECREATE");
  outfile.open(csvfile_name.c_str());
  outfile_abs.open(csvfile_abs_name.c_str());
  outfile_time.open(csvfile_time_name.c_str());
  outfile_time_first.open(csvfile_time_first_name.c_str());
  outfile_time_abs.open(csvfile_time_abs_name.c_str());
  outfile_time_first_abs.open(csvfile_time_first_abs_name.c_str());
  outfile_lappd.open(csvfile_lappd_name.c_str());
  outfile_lappd_abs.open(csvfile_lappd_abs_name.c_str());
  outfile_lappd_time.open(csvfile_lappd_time_name.c_str());
  outfile_lappd_time_first.open(csvfile_lappd_time_first_name.c_str());
  outfile_lappd_time_abs.open(csvfile_lappd_time_abs_name.c_str());
  outfile_lappd_time_first_abs.open(csvfile_lappd_time_first_abs_name.c_str());
  outfile_Rings.open(csvfile_Rings.c_str());
  outfile_MRD.open(csvfile_MRD.c_str());

  return true;
}


bool CNNImage::Execute(){

  Log("CNNImage tool: Executing ...",v_message,verbosity);

  //Get ANNIEEvent store
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){
    std::cerr<<"CNNImage tool: No ANNIEEvent store!"<<std::endl;
  }
  
  //---------------------------------------------------------------
  //-------------------Get ANNIEEvent objects----------------------
  //---------------------------------------------------------------

  bool get_ok=false;
  if (!isData){
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles",mcparticles); //needed to retrieve true vertex and direction
    if (!get_ok) {Log("CNNImage tool: Error retrieving MCParticles object",v_error,verbosity); return false;}
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
    if (!get_ok) {Log("CNNImage tool: Error retrieving MCHits object",v_error,verbosity); return false;}
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits", MCLAPPDHits);
    if (!get_ok) {Log("CNNImage tool: Error retrieving MCLAPPDHits object",v_error,verbosity); return false;}
    get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);
    if (!get_ok) {Log("CNNImage tool: Error retrieving TDCData object (MC)",v_error,verbosity); return false;}
    get_ok = m_data->Stores.at("RecoEvent")->Get("NRings",nrings);
    if (!get_ok) {Log("CNNImage tool: Error retrieving NRings (MC)",v_error,verbosity); return false;}
  } else {
    get_ok = m_data->Stores["ANNIEEvent"]->Get("Hits", Hits);
    if (!get_ok) {Log("CNNImage tool: Error retrieving Hits object",v_error,verbosity); return false;}
    get_ok = m_data->Stores.at("RecoEvent")->Get("RecoDigit",RecoDigits);
    if(not get_ok){ Log("CNNImage Tool: Error retrieving RecoDigit,true from RecoEvent!",v_error,verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData_Data);
    if (!get_ok) {Log("CNNImage tool: Error retrieving TDCData object (Data)",v_error,verbosity); return false;}
  }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  if (!get_ok) {Log("CNNImage tool: Error retrieving EventNumber",v_error,verbosity); return false;}
  /*get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  if (!get_ok) {Log("CNNImage tool: Error retrieving RunNumber",v_error,verbosity); return false;}
  get_ok = m_data->Stores["ANNIEEvent"]->Get("SubRunNumber",subrunnumber);
  if (!get_ok) {Log("CNNImage tool: Error retrieving SubRunNumber",v_error,verbosity); return false;}
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);
  if (!get_ok) {Log("CNNImage tool: Error retrieving EventTime",v_error,verbosity); return false;}*/
  get_ok = m_data->CStore.Get("NumMrdTimeClusters",mrdeventcounter);
  if (!get_ok) {Log("CNNImage tool: Error retrieving NumMrdTimeClusters",v_error,verbosity); return false;}
  if (mrdeventcounter > 0 ){
    get_ok = m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
    if (!get_ok) {Log("CNNImage tool: Error retrieving MrdDigitChankeys",v_error,verbosity); return false;}
  }
  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (!get_ok) {Log("CNNImage tool: Error retrieving MrdTimeClusters",v_error,verbosity); return false;}
  
  if (!isData) Log("CNNImage tool: # Rings: "+std::to_string(nrings)+", # of MRD clusters: "+std::to_string(mrdeventcounter),v_message,verbosity);
  else Log("CNNImage tool: # of MRD clusters: "+std::to_string(mrdeventcounter),v_message,verbosity);

  //Clear variables & containers
  charge.clear();
  time.clear();
  time_first.clear();
  hitpmt_detkeys.clear();
  charge_lappd.clear();
  time_lappd.clear();
  time_first_lappd.clear();
  hits_lappd.clear();
  total_charge_lappd.clear();

  for (unsigned int i_pmt=0; i_pmt<pmt_detkeys.size();i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    charge.emplace(detkey,0.);
    time.emplace(detkey,0.);
    time_first.emplace(detkey,0.);
  }

  for (unsigned int i_lappd=0; i_lappd<lappd_detkeys.size(); i_lappd++){
    unsigned long detkey = lappd_detkeys[i_lappd];
    total_charge_lappd.emplace(detkey,0);
    std::vector<std::vector<double>> temp_lappdXY;
    std::vector<std::vector<int>> temp_int_lappdXY;;
    for (int lappdX=0; lappdX<dimensionLAPPD; lappdX++){
      std::vector<double> temp_lappdY;
      std::vector<int> temp_int_lappdY;
      temp_lappdY.assign(dimensionLAPPD,0.);
      temp_int_lappdY.assign(dimensionLAPPD,0);
      temp_lappdXY.push_back(temp_lappdY);
      temp_int_lappdXY.push_back(temp_int_lappdY);
    }
    charge_lappd.emplace(detkey,temp_lappdXY);
    time_lappd.emplace(detkey,temp_lappdXY);
    time_first_lappd.emplace(detkey,temp_lappdXY);
    hits_lappd.emplace(detkey,temp_int_lappdXY);
  }

  //---------------------------------------------------------------
  //-------------------Iterate over MCHits ------------------------
  //---------------------------------------------------------------

  int vectsize;
  if (!isData){
    vectsize = MCHits->size();
    Log("CNNImage tool: MCHits size: "+std::to_string(vectsize),v_message,verbosity);
    total_hits_pmts=0;
    for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits){
      unsigned long chankey = apair.first;
      Detector* thistube = geom->ChannelToDetector(chankey);
      unsigned long detkey = thistube->GetDetectorID();
      Log("CNNImage tool: Read in MCHits for chankey: "+std::to_string(chankey)+", detkey: "+std::to_string(detkey),v_debug,verbosity);
      if (thistube->GetDetectorElement()=="Tank"){
        if (thistube->GetTankLocation()=="OD") continue;
        hitpmt_detkeys.push_back(detkey);
        std::vector<MCHit>& Hits = apair.second;
        int hits_pmt = 0;
        for (MCHit &ahit : Hits){
          if (ahit.GetTime()>-10. && ahit.GetTime()<40.){
            charge[detkey] += ahit.GetCharge();
            if (data_mode == "Normal") time[detkey] += ahit.GetTime();
            else if (data_mode == "Charge-Weighted") time[detkey] += (ahit.GetTime()*ahit.GetCharge());
            if (hits_pmt == 0) time_first[detkey] = ahit.GetTime();
            hits_pmt++;
          }
        }
        if (data_mode == "Normal" && hits_pmt>0) time[detkey]/=hits_pmt;         //use mean time of all hits on one PMT
        else if (data_mode == "Charge-Weighted" && charge[detkey]>0.) time[detkey] /= charge[detkey];
        total_hits_pmts++;
      }
    }
    Log("CNNImage tool: Done with MCHits loop",v_message,verbosity);
  } else {
    vectsize = RecoDigits->size();
    Log("CNNImage tool: RecoDigits size: "+std::to_string(vectsize),v_message,verbosity);
    total_hits_pmts = 0;
    for (unsigned int i_digit = 0; i_digit < RecoDigits->size(); i_digit++){
      RecoDigit thisdigit = RecoDigits->at(i_digit);
      int digitID = thisdigit.GetDetectorID();
      unsigned long chankey = pmtid_to_channelkey[digitID];
      double digitQ = thisdigit.GetCalCharge();
      double digitT = thisdigit.GetCalTime();
      hitpmt_detkeys.push_back(chankey);
      time[chankey] = digitT;
      time_first[chankey] = digitT;
      charge[chankey] = digitQ;
    }
    Log("CNNImage tool: Done with RecoDigits loop",v_message,verbosity);
  }

  //---------------------------------------------------------------
  //-------------------Iterate over MCLAPPDHits -------------------
  //---------------------------------------------------------------

  total_hits_lappds=0;
  total_charge_lappds=0;
  min_time_lappds=9999;
  max_time_lappds=-9999;
  if (!isData){
    vectsize = MCLAPPDHits->size();
    Log("CNNImage tool: MCLAPPDHits size: "+std::to_string(vectsize),v_message,verbosity);

    for(std::pair<unsigned long, std::vector<MCLAPPDHit>>&& apair : *MCLAPPDHits){
      unsigned long chankey = apair.first;
      Detector* det = geom->ChannelToDetector(chankey);
      unsigned long detkey = det->GetDetectorID();
      std::vector<MCLAPPDHit>& hits = apair.second;
      Log("CNNImage tool: Chankey: "+std::to_string(chankey)+", detkey: "+std::to_string(detkey),v_debug,verbosity);

      for (MCLAPPDHit& ahit : hits){
        std::vector<double> temp_pos = ahit.GetPosition();
        double x_lappd = temp_pos.at(0)-tank_center_x;    //global x-position (in detector coordinates)
        double y_lappd = temp_pos.at(1)-tank_center_y;    //global y-position (in detector coordinates)
        double z_lappd = temp_pos.at(2)-tank_center_z;    //global z-position (in detector coordinates)
        std::vector<double> local_pos = ahit.GetLocalPosition();

        double x_local_lappd = local_pos.at(1);           //local parallel position within LAPPD
        double y_local_lappd = local_pos.at(0);           //local transverse position within LAPPD
        
	//Get corresponding bin
        double lappd_charge = 1.0;
        double t_lappd = ahit.GetTime();
        int binx_lappd = round((x_local_lappd+0.1)/0.2*dimensionLAPPD);       //local positions can be positive and negative, 10cm > x_local_lappd > -10cm 
        int biny_lappd = round((y_local_lappd+0.1)/0.2*dimensionLAPPD);
        
        //Discard unphysical bins
        if (binx_lappd < 0) binx_lappd=0;         //FIXME: hits outside of -10cm...10cm should probably be just discarded
        if (biny_lappd < 0) biny_lappd=0;         //FIXME: hits outside of -10cm...10cm should probably be just discarded
        if (binx_lappd > dimensionLAPPD-1) binx_lappd = dimensionLAPPD-1;     //FIXME: hits outside of -10cm...10cm should probably be just discarded
        if (biny_lappd > dimensionLAPPD-1) biny_lappd = dimensionLAPPD-1;     //FIXME: hits outside of -10cm...10cm should probably be just discarded
        
	//Write charge & time hits into the histogram
        if (t_lappd>-10. && t_lappd<40){
          charge_lappd[detkey].at(binx_lappd).at(biny_lappd) += lappd_charge;
          time_lappd[detkey].at(binx_lappd).at(biny_lappd) += t_lappd;
          hits_lappd[detkey].at(binx_lappd).at(biny_lappd)++;
          total_charge_lappd[detkey] += lappd_charge;
          total_charge_lappds += lappd_charge;
          if (max_time_lappds < t_lappd) max_time_lappds = t_lappd;
          if (min_time_lappds > t_lappd) min_time_lappds = t_lappd;
          if (hits_lappd[detkey].at(binx_lappd).at(biny_lappd) == 1) time_first_lappd[detkey].at(binx_lappd).at(biny_lappd) = t_lappd;
        }
      total_hits_lappds++;
      }
    }
  }

  for (unsigned int i_lappd=0; i_lappd<lappd_detkeys.size(); i_lappd++){
    unsigned long detkey = lappd_detkeys.at(i_lappd);
    for (int i_lappdX=0; i_lappdX<dimensionLAPPD; i_lappdX++){
      for (int i_lappdY=0; i_lappdY<dimensionLAPPD; i_lappdY++){
        if (hits_lappd[detkey].at(i_lappdX).at(i_lappdY)>0) time_lappd[detkey].at(i_lappdX).at(i_lappdY)/=hits_lappd[detkey].at(i_lappdX).at(i_lappdY); 
      }
    }
  }

  if (!isData && min_time_lappds>0.) min_time_lappds=0.;

  //---------------------------------------------------------------
  //------------- Determine max+min values ------------------------
  //---------------------------------------------------------------


  maximum_pmts = 0;
  max_time_pmts = 0;
  min_time_pmts = 999.;
  double max_firsttime_pmts = 0.;
  double min_firsttime_pmts = 999.;
  total_charge_pmts = 0;
  for (unsigned int i_pmt=0;i_pmt<hitpmt_detkeys.size();i_pmt++){
    unsigned long detkey = hitpmt_detkeys[i_pmt];
    if (charge[detkey]>maximum_pmts) maximum_pmts = charge[detkey];
    total_charge_pmts+=charge[detkey];
    if (time[detkey]>max_time_pmts) max_time_pmts = time[detkey];
    if (time[detkey]<min_time_pmts) min_time_pmts = time[detkey];
    if (time_first[detkey]>max_firsttime_pmts) max_firsttime_pmts = time_first[detkey];
    if (time_first[detkey]<min_firsttime_pmts) min_firsttime_pmts = time_first[detkey];
  }
  maximum_lappds = 0;
  double max_firsttime_lappds = 0.;
  double min_firsttime_lappds = 999.;
  for (unsigned int i_lappd = 0; i_lappd < lappd_detkeys.size();i_lappd++){
    unsigned long detkey = lappd_detkeys[i_lappd];
    for (int iX=0; iX<dimensionLAPPD; iX++){
      for (int iY=0; iY<dimensionLAPPD; iY++){
        if (charge_lappd[detkey].at(iX).at(iY) > maximum_lappds) maximum_lappds = charge_lappd[detkey].at(iX).at(iY);
        if (time_first_lappd[detkey].at(iX).at(iY) > max_firsttime_lappds) max_firsttime_lappds = time_first_lappd[detkey].at(iX).at(iY);
        if (time_first_lappd[detkey].at(iX).at(iY) < min_firsttime_lappds) min_firsttime_lappds = time_first_lappd[detkey].at(iX).at(iY);
      }
    }
  }
  double global_min_time = (min_time_pmts < min_time_lappds)? min_time_pmts : min_time_lappds;
  double global_max_time = (max_time_pmts > max_time_lappds)? max_time_pmts : max_time_lappds;
  double global_max_charge = (maximum_pmts > maximum_lappds)? maximum_pmts : maximum_lappds;
  double global_min_charge = 0.;
  double global_min_time_first = (min_firsttime_pmts < min_firsttime_lappds)? min_firsttime_pmts : min_firsttime_lappds;
  double global_max_time_first = (max_firsttime_pmts > max_firsttime_lappds)? max_firsttime_pmts : max_firsttime_lappds;

  Log("CNNImage tool: Min_time_pmts: "+std::to_string(min_time_pmts)+", min_time_lappds: "+std::to_string(min_time_lappds)+", global_min_time: "+std::to_string(global_min_time),v_message,verbosity);

  if (fabs(global_max_time-global_min_time)<0.01) global_max_time = global_min_time+1;
  if (fabs(global_max_time_first-global_min_time_first)<0.01) global_max_time_first = global_min_time_first+1;
  if (global_max_charge<0.001) global_max_charge=1;

  //---------------------------------------------------------------
  //-------------- Readout MRD ------------------------------------
  //---------------------------------------------------------------

  int num_mrd_paddles_cluster=0;
  int num_mrd_layers_cluster=0;
  int num_mrd_conslayers_cluster=0;
  int num_mrd_adjacent_cluster=0;
  double mrd_padperlayer_cluster=0.;
  bool layer_occupied_cluster[11] = {0};
  double mrd_paddlesize_cluster[11];

  if(MrdTimeClusters.size()==0){
    //No entries in MrdTimeClusters object, don't read out anything
    Log("CNNImage tool: No MRDClusters entries.",v_message,verbosity);
  } else {
    std::vector<std::vector<double>> mrd_hits;
    for (int i_layer=0; i_layer<11; i_layer++){
      std::vector<double> empty_hits;
      mrd_hits.push_back(empty_hits);
    }
    std::vector<int> temp_cons_layers;
    for (int i_cluster=0; i_cluster < (int) MrdTimeClusters.size(); i_cluster++){
      if (i_cluster > 0) continue; 		//only evaluate first cluster, if there are multiple
      std::vector<int> single_cluster = MrdTimeClusters.at(i_cluster);
      for (int i_pmt = 0; i_pmt < (int) single_cluster.size(); i_pmt++){
        int idigit = single_cluster.at(i_pmt);
        unsigned long chankey = mrddigitchankeysthisevent.at(idigit);
        Detector *thedetector = geom->ChannelToDetector(chankey);
        if(thedetector->GetDetectorElement()!="MRD") {
          continue;          // this is a veto hit, not an MRD hit.
        }
        num_mrd_paddles_cluster++;
	int detkey = thedetector->GetDetectorID();
        Paddle *apaddle = geom->GetDetectorPaddle(detkey);
        int layer = apaddle->GetLayer();
        layer_occupied_cluster[layer-1]=true;
        if (apaddle->GetOrientation()==1) {
          mrd_hits.at(layer-2).push_back(0.5*(apaddle->GetXmin()+apaddle->GetXmax()));
          mrd_paddlesize_cluster[layer-2]=apaddle->GetPaddleWidth();
        }
        else if (apaddle->GetOrientation()==0) {
          mrd_hits.at(layer-2).push_back(0.5*(apaddle->GetYmin()+apaddle->GetYmax()));
          mrd_paddlesize_cluster[layer-2]=apaddle->GetPaddleWidth();
        }
      }
      if (num_mrd_paddles_cluster > 0) {
        for (int i_layer=0;i_layer<11;i_layer++){
          if (layer_occupied_cluster[i_layer]==true) {
            num_mrd_layers_cluster++;
            if (num_mrd_conslayers_cluster==0) num_mrd_conslayers_cluster++;
            else {
              if (layer_occupied_cluster[i_layer-1]==true) num_mrd_conslayers_cluster++;
              else {
                temp_cons_layers.push_back(num_mrd_conslayers_cluster);
                num_mrd_conslayers_cluster=0;
              }
            }
          }
          for (unsigned int i_hitpaddle=0; i_hitpaddle<mrd_hits.at(i_layer).size(); i_hitpaddle++){
            for (unsigned int j_hitpaddle= i_hitpaddle+1; j_hitpaddle < mrd_hits.at(i_layer).size(); j_hitpaddle++){
              if (fabs(mrd_hits.at(i_layer).at(i_hitpaddle)-mrd_hits.at(i_layer).at(j_hitpaddle))<=mrd_paddlesize_cluster[i_layer]+0.001) num_mrd_adjacent_cluster++;
            }
          }	 		
        }
      }
      temp_cons_layers.push_back(num_mrd_conslayers_cluster);
      std::vector<int>::iterator it = std::max_element(temp_cons_layers.begin(),temp_cons_layers.end());
      if (it != temp_cons_layers.end()) num_mrd_conslayers_cluster = *it;
      else num_mrd_conslayers_cluster=0;
      if (num_mrd_layers_cluster > 0.) mrd_padperlayer_cluster = double(num_mrd_paddles_cluster)/num_mrd_layers_cluster;
    }
  }

  //---------------------------------------------------------------
  //-------------- Create CNN images ------------------------------
  //---------------------------------------------------------------

  // Define histogram as an intermediate step to the CNN

  std::stringstream ss_cnn, ss_title_cnn, ss_cnn_time, ss_title_cnn_time, ss_cnn_pmtwise, ss_title_cnn_pmtwise, ss_cnn_time_pmtwise, ss_title_cnn_time_pmtwise;
  std::stringstream ss_cnn_abs, ss_title_cnn_abs, ss_cnn_time_first, ss_title_cnn_time_first, ss_cnn_time_abs, ss_title_cnn_time_abs, ss_cnn_time_first_abs, ss_title_cnn_time_first_abs, ss_cnn_abs_pmtwise, ss_title_cnn_abs_pmtwise, ss_cnn_time_first_pmtwise, ss_title_cnn_time_first_pmtwise, ss_cnn_time_abs_pmtwise, ss_title_cnn_time_abs_pmtwise, ss_cnn_time_first_abs_pmtwise, ss_title_cnn_time_first_abs_pmtwise;

  ss_cnn<<"hist_cnn"<<evnum;
  ss_title_cnn<<"EventDisplay (CNN), Event "<<evnum;
  ss_cnn_abs<<"hist_cnn_abs"<<evnum;
  ss_title_cnn_abs<<"EventDisplay Charge (CNN), Event "<<evnum;
  ss_cnn_time<<"hist_cnn_time"<<evnum;
  ss_title_cnn_time<<"EventDisplay Time (CNN), Event "<<evnum;
  ss_cnn_time_first<<"hist_cnn_time_first"<<evnum;
  ss_title_cnn_time_first<<"EventDisplay First Times (CNN), Event "<<evnum;
  ss_cnn_time_abs<<"hist_cnn_time_abs"<<evnum;
  ss_title_cnn_time_abs<<"EventDisplay Absolute Times (CNN), Event "<<evnum;
  ss_cnn_time_first_abs<<"hist_cnn_time_first_abs"<<evnum;
  ss_title_cnn_time_first_abs<<"EventDisplay Absolute First Times(CNN), Event "<<evnum;
  ss_cnn_pmtwise<<"hist_cnn_pmtwise"<<evnum;
  ss_title_cnn_pmtwise<<"EventDisplay (CNN, pmt wise), Event "<<evnum;
  ss_cnn_abs_pmtwise<<"hist_cnn_abs_pmtwise"<<evnum;
  ss_title_cnn_abs_pmtwise<<"EventDisplay Charge (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_pmtwise << "hist_cnn_time_pmtwise"<<evnum;
  ss_title_cnn_time_pmtwise <<"EventDisplay Time (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_first_pmtwise << "hist_cnn_time_first_pmtwise"<<evnum;
  ss_title_cnn_time_first_pmtwise <<"EventDisplay First Times (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_abs_pmtwise << "hist_cnn_time_abs_pmtwise"<<evnum;
  ss_title_cnn_time_abs_pmtwise <<"EventDisplay Absolute Times (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_first_abs_pmtwise << "hist_cnn_time_first_abs_pmtwise"<<evnum;
  ss_title_cnn_time_first_abs_pmtwise <<"EventDisplay Absolute First Times (CNN, pmt wise), Event "<<evnum;
  
  TH2F *hist_cnn = new TH2F(ss_cnn.str().c_str(),ss_title_cnn.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_abs = new TH2F(ss_cnn_abs.str().c_str(),ss_title_cnn_abs.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time = new TH2F(ss_cnn_time.str().c_str(),ss_title_cnn_time.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time_first = new TH2F(ss_cnn_time_first.str().c_str(),ss_title_cnn_time_first.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time_abs = new TH2F(ss_cnn_time_abs.str().c_str(),ss_title_cnn_time_abs.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time_first_abs = new TH2F(ss_cnn_time_first_abs.str().c_str(),ss_title_cnn_time_first_abs.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_pmtwise = new TH2F(ss_cnn_pmtwise.str().c_str(),ss_title_cnn_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_abs_pmtwise = new TH2F(ss_cnn_abs_pmtwise.str().c_str(),ss_title_cnn_abs_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_pmtwise = new TH2F(ss_cnn_time_pmtwise.str().c_str(),ss_title_cnn_time_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_first_pmtwise = new TH2F(ss_cnn_time_first_pmtwise.str().c_str(),ss_title_cnn_time_first_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_abs_pmtwise = new TH2F(ss_cnn_time_abs_pmtwise.str().c_str(),ss_title_cnn_time_abs_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_first_abs_pmtwise = new TH2F(ss_cnn_time_first_abs_pmtwise.str().c_str(),ss_title_cnn_time_first_abs_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  
  // Define LAPPD histograms (1 per LAPPD)
  std::vector<TH2F*> hists_lappd, hists_abs_lappd, hists_time_lappd, hists_time_first_lappd, hists_time_abs_lappd, hists_time_first_abs_lappd;
  for (unsigned int i_lappd = 0; i_lappd < lappd_detkeys.size(); i_lappd++){
    std::stringstream ss_hist_lappd, ss_title_lappd, ss_hist_abs_lappd, ss_title_abs_lappd, ss_hist_time_lappd, ss_title_time_lappd, ss_hist_time_first_lappd, ss_title_time_first_lappd, ss_hist_time_abs_lappd, ss_title_time_abs_lappd, ss_hist_time_first_abs_lappd, ss_title_time_first_abs_lappd;
    ss_hist_lappd<<"hist_lappd"<<i_lappd<<"_ev"<<evnum;
    ss_title_lappd<<"EventDisplay (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    ss_hist_abs_lappd<<"hist_lappd_abs"<<i_lappd<<"_ev"<<evnum;
    ss_title_abs_lappd<<"EventDisplay Charge (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    ss_hist_time_lappd<<"hist_time_lappd"<<i_lappd<<"_ev"<<evnum;
    ss_title_time_lappd<<"EventDisplay Time (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    ss_hist_time_first_lappd<<"hist_time_first_lappd"<<i_lappd<<"_ev"<<evnum;
    ss_title_time_first_lappd<<"EventDisplay First Times (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    ss_hist_time_abs_lappd<<"hist_time_abs_lappd"<<i_lappd<<"_ev"<<evnum;
    ss_title_time_abs_lappd<<"EventDisplay Absolute Times (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    ss_hist_time_first_abs_lappd<<"hist_time_first_abs_lappd"<<i_lappd<<"_ev"<<evnum;
    ss_title_time_first_abs_lappd<<"EventDisplay Absolute First Times (CNN), LAPPD "<<i_lappd<<", Event "<<evnum;
    TH2F *hist_lappd = new TH2F(ss_hist_lappd.str().c_str(),ss_title_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_lappd.push_back(hist_lappd);
    TH2F *hist_abs_lappd = new TH2F(ss_hist_abs_lappd.str().c_str(),ss_title_abs_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_abs_lappd.push_back(hist_abs_lappd);
    TH2F *hist_time_lappd = new TH2F(ss_hist_time_lappd.str().c_str(),ss_title_time_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_time_lappd.push_back(hist_time_lappd);
    TH2F *hist_time_first_lappd = new TH2F(ss_hist_time_first_lappd.str().c_str(),ss_title_time_first_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_time_first_lappd.push_back(hist_time_first_lappd);
    TH2F *hist_time_abs_lappd = new TH2F(ss_hist_time_abs_lappd.str().c_str(),ss_title_time_abs_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_time_abs_lappd.push_back(hist_time_abs_lappd);
    TH2F *hist_time_first_abs_lappd = new TH2F(ss_hist_time_first_abs_lappd.str().c_str(),ss_title_time_first_abs_lappd.str().c_str(),dimensionLAPPD,0,dimensionLAPPD,dimensionLAPPD,0,dimensionLAPPD);    //resolution in x/y direction should be 1cm, over a total length of 20cm each
    hists_time_first_abs_lappd.push_back(hist_time_first_abs_lappd);
  }

  //---------------------------------------------------------------
  //---------------- Fill PMT images ------------------------------
  //---------------------------------------------------------------

  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];

    //Convert PMT position to 2D (x,y)-position
    double x,y;
    Position pmt_pos(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]);
    ConvertPositionTo2D(pmt_pos, x, y);
    
    //Get corresponding bins for 2D (x,y)-position
    int binx = hist_cnn->GetXaxis()->FindBin(x);
    int biny = hist_cnn->GetYaxis()->FindBin(y);
    Log("CNNImage tool: Filling binx: "+std::to_string(binx)+", biny: "+std::to_string(biny)+" with charge: "+std::to_string(charge[detkey])+", time: "+std::to_string(time[detkey]),v_message,verbosity);

    if (maximum_pmts < 0.001) maximum_pmts = 1.;
    if (fabs(max_time_pmts) < 0.001) max_time_pmts = 1.;

    //Fill geometric histograms
    double charge_fill = charge[detkey]/global_max_charge;
    double time_fill = 0.;
    double time_first_fill = 0.;
    if (charge_fill > 1e-10) {
     time_fill = (time[detkey]-global_min_time)/(global_max_time-global_min_time);
     time_first_fill = (time_first[detkey]-global_min_time_first)/(global_max_time_first-global_min_time_first);
    }

    hist_cnn->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+charge_fill);
    hist_cnn_abs->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+charge[detkey]);
    hist_cnn_time->SetBinContent(binx,biny,time_fill);
    hist_cnn_time_abs->SetBinContent(binx,biny,time[detkey]);
    hist_cnn_time_first->SetBinContent(binx,biny,time_first_fill);
    hist_cnn_time_first_abs->SetBinContent(binx,biny,time_first[detkey]);

    //Fill the pmt-wise histograms
    if ((y_pmt[detkey]>=max_y || y_pmt[detkey]<=min_y) && !includeTopBottom) continue;       //don't include endcaps in the pmt-wise histogram for now
    if (y_pmt[detkey]>=max_y) ConvertPositionTo2D_Top(pmt_pos,x,y);
    if (y_pmt[detkey]<=min_y) ConvertPositionTo2D_Bottom(pmt_pos,x,y);
    double xCorr, yCorr;
    xCorr = (round(100*x)/100.);
    yCorr = (round(100*y)/100.);
    if (fabs(xCorr-0.52)<0.001) xCorr = 0.53;
    if (fabs(xCorr-0.47)<0.001) xCorr = 0.48;
    std::vector<double>::iterator it_x, it_y;
    if (y_pmt[detkey]>=max_y){
      if (fabs(xCorr-0.19)<0.001) xCorr = 0.21;
      if (fabs(xCorr-0.34)<0.001) xCorr = 0.32;
      if (fabs(xCorr-0.48)<0.001) xCorr = 0.50;
      if (fabs(xCorr-0.66)<0.001) xCorr = 0.68;
      it_x = std::find(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end(),xCorr);
    }
    else if (y_pmt[detkey]<=min_y){
      if (fabs(xCorr-0.23)<0.001) xCorr = 0.22;
      if (fabs(xCorr-0.77)<0.001) xCorr = 0.78;
      it_x = std::find(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end(),xCorr);
    }
    else {
      it_x = std::find(vec_pmt2D_x.begin(),vec_pmt2D_x.end(),xCorr);
    }
    it_y = std::find(vec_pmt2D_y.begin(),vec_pmt2D_y.end(),yCorr);
    int index_x, index_y;
    if (y_pmt[detkey]>=max_y) index_x = std::distance(vec_pmt2D_x_Top.begin(),it_x);
    else if (y_pmt[detkey]<=min_y) index_x = std::distance(vec_pmt2D_x_Bottom.begin(),it_x);
    else index_x = std::distance(vec_pmt2D_x.begin(),it_x);
    index_y = std::distance(vec_pmt2D_y.begin(),it_y);
    hist_cnn_pmtwise->SetBinContent(index_x+1,index_y+1,charge_fill);
    hist_cnn_abs_pmtwise->SetBinContent(index_x+1,index_y+1,charge[detkey]);
    hist_cnn_time_pmtwise->SetBinContent(index_x+1,index_y+1,time_fill);
    hist_cnn_time_abs_pmtwise->SetBinContent(index_x+1,index_y+1,time[detkey]);
    hist_cnn_time_first_pmtwise->SetBinContent(index_x+1,index_y+1,time_first_fill);
    hist_cnn_time_first_abs_pmtwise->SetBinContent(index_x+1,index_y+1,time_first[detkey]);
  }
  
  //---------------------------------------------------------------
  //---------------- Fill LAPPD images ----------------------------
  //---------------------------------------------------------------

  for (unsigned int i_lappd=0; i_lappd<lappd_detkeys.size(); i_lappd++){
    unsigned long detkey = lappd_detkeys.at(i_lappd);
    for (int iX=0; iX < dimensionLAPPD; iX++){
      for (int iY=0; iY < dimensionLAPPD; iY++){
        double lappd_charge_fill = charge_lappd[detkey].at(iX).at(iY)/global_max_charge;
        double lappd_time_fill = 0.;
        double lappd_time_first_fill = 0.;
        if (lappd_charge_fill > 1e-10) {
          lappd_time_fill = (time_lappd[detkey].at(iX).at(iY)-global_min_time)/(global_max_time-global_min_time);
          lappd_time_first_fill = (time_first_lappd[detkey].at(iX).at(iY)-global_min_time)/(global_max_time-global_min_time);
        }

        if (lappd_time_fill < 0)  std::cout <<"Min LAPPD time: "<<min_time_lappds<<", Max LAPPD time: "<<max_time_lappds<<", time_lappd: "<<time_lappd[detkey].at(iX).at(iY)<<", fill time: "<<lappd_time_fill<<std::endl;
        hists_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,lappd_charge_fill);
        hists_abs_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,charge_lappd[detkey].at(iX).at(iY));
        hists_time_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,lappd_time_fill);
        hists_time_first_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,lappd_time_first_fill);
        hists_time_abs_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,time_lappd[detkey].at(iX).at(iY));
        hists_time_first_abs_lappd.at(i_lappd)->SetBinContent(iX+1,iY+1,time_first_lappd[detkey].at(iX).at(iY));
      }
    }
  }

  //---------------------------------------------------------------
  //---------------- Write information to csv----------------------
  //---------------------------------------------------------------

  //
  // Save information from histogram to csv file
  // (1 line corresponds to 1 event, histogram entries flattened out to a 1D array)
  //

  bool passed_eventselection;
  m_data->Stores["RecoEvent"]->Get("EventCutStatus",passed_eventselection);
  Log("CNNImage tool: Passed Event Selection: "+std::to_string(passed_eventselection),v_message,verbosity);

  if (passed_eventselection) {
    // Save Rings and MRD information
    if (!isData) outfile_Rings << nrings << endl;
    outfile_MRD << mrdeventcounter<< ","<< num_mrd_paddles_cluster <<","<< num_mrd_layers_cluster<<","<<num_mrd_conslayers_cluster<<","<<num_mrd_adjacent_cluster<<","<<mrd_padperlayer_cluster<< endl;

    // Save root histograms
    hist_cnn->Write();
    hist_cnn_abs->Write();
    hist_cnn_time->Write();
    hist_cnn_time_first->Write();
    hist_cnn_time_abs->Write();
    hist_cnn_time_first_abs->Write();
    hist_cnn_pmtwise->Write();
    hist_cnn_abs_pmtwise->Write();
    hist_cnn_time_pmtwise->Write();
    hist_cnn_time_first_pmtwise->Write();
    hist_cnn_time_abs_pmtwise->Write();
    hist_cnn_time_first_abs_pmtwise->Write();
    for (unsigned int i_lappd=0; i_lappd<lappd_detkeys.size();i_lappd++){
      hists_lappd.at(i_lappd)->Write();
      hists_abs_lappd.at(i_lappd)->Write();
      hists_time_lappd.at(i_lappd)->Write();
      hists_time_first_lappd.at(i_lappd)->Write();
      hists_time_abs_lappd.at(i_lappd)->Write();
      hists_time_first_abs_lappd.at(i_lappd)->Write();
    }

    if (save_mode == "Geometric"){
      for (int i_binY=0; i_binY < hist_cnn->GetNbinsY();i_binY++){
        for (int i_binX=0; i_binX < hist_cnn->GetNbinsX();i_binX++){
          outfile << hist_cnn->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn->GetNbinsX()-1 || i_binY!=hist_cnn->GetNbinsY()-1) outfile<<",";
          outfile_abs << hist_cnn_abs->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_abs->GetNbinsX()-1 || i_binY!=hist_cnn_abs->GetNbinsY()-1) outfile_abs<<",";
          outfile_time << hist_cnn_time->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time->GetNbinsX()-1 || i_binY!=hist_cnn_time->GetNbinsY()-1) outfile_time<<",";    
          outfile_time_first << hist_cnn_time_first->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_first->GetNbinsX()-1 || i_binY!=hist_cnn_time_first->GetNbinsY()-1) outfile_time_first<<",";    
          outfile_time_abs << hist_cnn_time_abs->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_abs->GetNbinsX()-1 || i_binY!=hist_cnn_time_abs->GetNbinsY()-1) outfile_time_abs<<",";    
          outfile_time_first_abs << hist_cnn_time_first_abs->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_first_abs->GetNbinsX()-1 || i_binY!=hist_cnn_time_first_abs->GetNbinsY()-1) outfile_time_first_abs<<",";    
        }
      }
    } else if (save_mode == "PMT-wise"){
      for (int i_binY=0; i_binY < hist_cnn_pmtwise->GetNbinsY();i_binY++){
        for (int i_binX=0; i_binX < hist_cnn_pmtwise->GetNbinsX();i_binX++){
          outfile << hist_cnn_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_pmtwise->GetNbinsY()-1) outfile<<",";
          outfile_abs << hist_cnn_abs_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_abs_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_abs_pmtwise->GetNbinsY()-1) outfile_abs<<",";
          outfile_time << hist_cnn_time_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_pmtwise->GetNbinsY()-1) outfile_time<<",";    
          outfile_time_first << hist_cnn_time_first_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_first_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_first_pmtwise->GetNbinsY()-1) outfile_time_first<<",";    
          outfile_time_abs << hist_cnn_time_abs_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_abs_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_abs_pmtwise->GetNbinsY()-1) outfile_time_abs<<",";    
          outfile_time_first_abs << hist_cnn_time_first_abs_pmtwise->GetBinContent(i_binX+1,i_binY+1);
          if (i_binX != hist_cnn_time_first_abs_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_first_abs_pmtwise->GetNbinsY()-1) outfile_time_first_abs<<",";    
        }
      }
      for (unsigned int i_lappd=0;i_lappd<lappd_detkeys.size();i_lappd++){
        for (int i_binY=0; i_binY < hists_lappd.at(i_lappd)->GetNbinsY();i_binY++){
          for (int i_binX=0; i_binX < hists_lappd.at(i_lappd)->GetNbinsX();i_binX++){
            outfile_lappd << hists_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            outfile_lappd_abs << hists_abs_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            if (i_binX != hists_lappd.at(i_lappd)->GetNbinsX()-1 || i_binY!=hists_lappd.at(i_lappd)->GetNbinsY()-1 || i_lappd != lappd_detkeys.size()-1) {
              outfile_lappd<<",";
              outfile_lappd_abs<<",";
            }
            outfile_lappd_time << hists_time_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            outfile_lappd_time_first << hists_time_first_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            outfile_lappd_time_abs << hists_time_abs_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            outfile_lappd_time_first_abs << hists_time_first_abs_lappd.at(i_lappd)->GetBinContent(i_binX+1,i_binY+1);
            if (i_binX != hists_time_lappd.at(i_lappd)->GetNbinsX()-1 || i_binY!=hists_time_lappd.at(i_lappd)->GetNbinsY()-1 || i_lappd != lappd_detkeys.size()-1) {
              outfile_lappd_time<<",";
              outfile_lappd_time_first<<",";
              outfile_lappd_time_abs<<",";
              outfile_lappd_time_first_abs<<",";
            }
          }
        }
      }
      outfile_lappd << std::endl;
      outfile_lappd_abs << std::endl;
      outfile_lappd_time << std::endl;
      outfile_lappd_time_first << std::endl;
      outfile_lappd_time_abs << std::endl;
      outfile_lappd_time_first_abs << std::endl;
    }
    outfile << std::endl;
    outfile_abs << std::endl;
    outfile_time << std::endl;
    outfile_time_first << std::endl;
    outfile_time_abs << std::endl;
    outfile_time_first_abs << std::endl;
  }

  //Delete histograms after accessing them
  delete hist_cnn;
  delete hist_cnn_abs;
  delete hist_cnn_time;
  delete hist_cnn_time_first;
  delete hist_cnn_time_abs;
  delete hist_cnn_time_first_abs;
  delete hist_cnn_pmtwise;
  delete hist_cnn_abs_pmtwise;
  delete hist_cnn_time_pmtwise;
  delete hist_cnn_time_first_pmtwise;
  delete hist_cnn_time_abs_pmtwise;
  delete hist_cnn_time_first_abs_pmtwise;
  for (unsigned int i_lappd=0; i_lappd<lappd_detkeys.size();i_lappd++){
    delete hists_lappd.at(i_lappd);
    delete hists_abs_lappd.at(i_lappd);
    delete hists_time_lappd.at(i_lappd);
    delete hists_time_first_lappd.at(i_lappd);
    delete hists_time_abs_lappd.at(i_lappd);
    delete hists_time_first_abs_lappd.at(i_lappd);
  }

  return true;
}


bool CNNImage::Finalise(){

  Log("CNNImage tool: Finalising ...",v_message,verbosity);
  file->Close();
  outfile.close();
  outfile_abs.close();
  outfile_time.close();
  outfile_time_first.close();
  outfile_time_abs.close();
  outfile_time_first_abs.close();
  outfile_lappd.close();
  outfile_lappd_abs.close();
  outfile_lappd_time.close();
  outfile_lappd_time_first.close();
  outfile_lappd_time_abs.close();
  outfile_lappd_time_first_abs.close();
  outfile_Rings.close();
  outfile_MRD.close();
  return true;
}

void CNNImage::ConvertPositionTo2D(Position xyz_pos, double &x, double &y){

    if (fabs(xyz_pos.Y()-max_y)<0.01){
      //top PMTs
      x=0.5-size_top_drawing*xyz_pos.X()/tank_radius;
      y=0.5+((0.45*tank_height)/tank_radius+1)*size_top_drawing-size_top_drawing*xyz_pos.Z()/tank_radius;
    } else if (fabs(xyz_pos.Y()-min_y)<0.01 || fabs(xyz_pos.Y()+1.30912)<0.01){
      //bottom PMTs
      x=0.5-size_top_drawing*xyz_pos.X()/tank_radius;
      y=0.5-(0.45*tank_height/tank_radius+1)*size_top_drawing+size_top_drawing*xyz_pos.Z()/tank_radius;
    } else {
      //barrel PMTs
      double phi=0.;
      if (xyz_pos.X()>0 && xyz_pos.Z()>0) phi = atan(xyz_pos.Z()/xyz_pos.X())+TMath::Pi()/2;
      else if (xyz_pos.X()>0 && xyz_pos.Z()<0) phi = atan(xyz_pos.X()/-xyz_pos.Z());
      else if (xyz_pos.X()<0 && xyz_pos.Z()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.Z()/xyz_pos.X());
      else if (xyz_pos.X()<0 && xyz_pos.Z()>0) phi = TMath::Pi()+atan(-xyz_pos.X()/xyz_pos.Z());
      else if (fabs(xyz_pos.X())<0.0001){
        if (xyz_pos.Z()>0) phi = TMath::Pi();
        else if (xyz_pos.Z()<0) phi = 2*TMath::Pi();
      }
      else if (fabs(xyz_pos.Z())<0.0001){
        if (xyz_pos.X()>0) phi = 0.5*TMath::Pi();
        else if (xyz_pos.X()<0) phi = 3*TMath::Pi()/2;
      }
      else phi = 0.;
      if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
      phi-=TMath::Pi();
      if (phi < - TMath::Pi()) phi = -TMath::Pi();
      if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
      x=0.5+phi*size_top_drawing;
      y=0.5+xyz_pos.Y()/tank_height*tank_height/tank_radius*size_top_drawing;
      }
}

void CNNImage::ConvertPositionTo2D_Top(Position xyz_pos, double &x, double &y){

	double rho = sqrt(xyz_pos.X()*xyz_pos.X()+xyz_pos.Z()*xyz_pos.Z());
	if (rho < 0.8) y = 0.9;
	else y = 0.8;

	double phi = 0.;
	if (xyz_pos.X()>0 && xyz_pos.Z()>0) phi = atan(xyz_pos.Z()/xyz_pos.X())+TMath::Pi()/2;
        else if (xyz_pos.X()>0 && xyz_pos.Z()<0) phi = atan(xyz_pos.X()/-xyz_pos.Z());
        else if (xyz_pos.X()<0 && xyz_pos.Z()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.Z()/xyz_pos.X());
        else if (xyz_pos.X()<0 && xyz_pos.Z()>0) phi = TMath::Pi()+atan(-xyz_pos.X()/xyz_pos.Z());
        else if (fabs(xyz_pos.X())<0.0001){
          if (xyz_pos.Z()>0) phi = TMath::Pi();
          else if (xyz_pos.Z()<0) phi = 2*TMath::Pi();
        }
        else if (fabs(xyz_pos.Z())<0.0001){
          if (xyz_pos.X()>0) phi = 0.5*TMath::Pi();
          else if (xyz_pos.X()<0) phi = 3*TMath::Pi()/2;
        }
        else phi = 0.;
        if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
        phi-=TMath::Pi();
        if (phi < - TMath::Pi()) phi = -TMath::Pi();
        if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
        x=0.5+phi*size_top_drawing; 

}

void CNNImage::ConvertPositionTo2D_Bottom(Position xyz_pos, double &x, double &y){

        double rho = sqrt(xyz_pos.X()*xyz_pos.X()+xyz_pos.Z()*xyz_pos.Z());
        if (rho < 0.6) y = -0.9;
        else y = -0.8;
 
	double phi = 0.;
        if (xyz_pos.X()>0 && xyz_pos.Z()>0) phi = atan(xyz_pos.Z()/xyz_pos.X())+TMath::Pi()/2;
        else if (xyz_pos.X()>0 && xyz_pos.Z()<0) phi = atan(xyz_pos.X()/-xyz_pos.Z());
        else if (xyz_pos.X()<0 && xyz_pos.Z()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.Z()/xyz_pos.X());
        else if (xyz_pos.X()<0 && xyz_pos.Z()>0) phi = TMath::Pi()+atan(-xyz_pos.X()/xyz_pos.Z());
        else if (fabs(xyz_pos.X())<0.0001){
          if (xyz_pos.Z()>0) phi = TMath::Pi();
          else if (xyz_pos.Z()<0) phi = 2*TMath::Pi();
        }
        else if (fabs(xyz_pos.Z())<0.0001){
          if (xyz_pos.X()>0) phi = 0.5*TMath::Pi();
          else if (xyz_pos.X()<0) phi = 3*TMath::Pi()/2;
        }
        else phi = 0.;
        if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
        phi-=TMath::Pi();
        if (phi < - TMath::Pi()) phi = -TMath::Pi();
        if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
        x=0.5+phi*size_top_drawing;

}
