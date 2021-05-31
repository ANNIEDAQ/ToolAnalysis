#include "MonitorTankTime.h"

MonitorTankTime::MonitorTankTime():Tool(){}


bool MonitorTankTime::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();
  m_data= &data; //assigning transient data pointer

  Log("Tool MonitorTankTime: Initialising....",v_message,verbosity);

  //-------------------------------------------------------
  //---------------Initialise config file------------------
  //-------------------------------------------------------

  m_data= &data; //assigning transient data pointer

  //only for debugging memory leaks, otherwise comment out
  //std::cout <<"List of Objects (beginning of initialise): "<<std::endl;
  //gObjectTable->Print();

  //-------------------------------------------------------
  //-----------------Get Configuration---------------------
  //-------------------------------------------------------

  m_variables.Get("OutputPath",outpath_temp);
  m_variables.Get("ActiveSlots",active_slots);
  m_variables.Get("SignalChannels",signal_channels);
  m_variables.Get("DisabledChannels",disabled_channels);
  m_variables.Get("StartTime", StartTime);
  m_variables.Get("UpdateFrequency",update_frequency);
  m_variables.Get("PlotConfiguration",plot_configuration);
  m_variables.Get("PathMonitoring",path_monitoring);
  m_variables.Get("ImageFormat",img_extension);
  m_variables.Get("ForceUpdate",force_update);
  m_variables.Get("DrawMarker",draw_marker);
  m_variables.Get("DrawSingle",draw_single);
  m_variables.Get("verbose",verbosity);

  if (verbosity > 2) std::cout <<"MonitorTankTime: Outpath (temporary): "<<outpath_temp<<std::endl;
  if (outpath_temp == "fromStore") m_data->CStore.Get("OutPath",outpath);
  else outpath = outpath_temp;
  if (verbosity > 2) std::cout <<"MonitorTankTime: Output path for plots is "<<outpath<<std::endl;
  if (update_frequency < 0.1) {
  if (verbosity > 0) std::cout <<"MonitorTankTime: Update Frequency of every "<<update_frequency<<" mins is too high. Setting default value of 5 mins."<<std::endl;
    update_frequency = 5.;
  }
  //Don't enable the drawing of single channel histograms/graphs by default --> too many plots!
  if (draw_single != 0 && draw_single !=1) {
    draw_single = 0;
  }
  //default should be no forced update of the monitoring plots every execute step
  if (force_update !=0 && force_update !=1) {
    force_update = 0;
  }
  //check if the image format is jpg or png
  if (!(img_extension == "png" || img_extension == "jpg" || img_extension == "jpeg")){
    img_extension = "jpg";
  }

  //-------------------------------------------------------
  //-----------------Get active slots----------------------
  //-------------------------------------------------------

  num_active_slots=0;
  num_active_slots_cr1=0;
  num_active_slots_cr2=0;
  num_active_slots_cr3=0;

  ifstream file(active_slots.c_str());
  unsigned int temp_crate, temp_slot;
  while (!file.eof()){

    file>>temp_crate>>temp_slot;
    Log("MonitorTankTime: Reading in active Slots: Crate "+std::to_string(temp_crate)+", Card "+std::to_string(temp_slot),v_message,verbosity);
    if (file.eof()) {
      if (verbosity >= 1) std::cout << std::endl;
      break;
    }
    if (int(temp_slot) < 2 || int(temp_slot) > num_slots_tank){
      Log("ERROR (MonitorTankTime): Specified slot "+std::to_string(temp_slot)+" out of range for VME crates [2...21]. Continue with next entry.",v_error,verbosity);
      continue;
    }
    if (!(std::find(crate_numbers.begin(),crate_numbers.end(),temp_crate)!=crate_numbers.end())) crate_numbers.push_back(temp_crate);
    std::vector<int>::iterator it = std::find(crate_numbers.begin(),crate_numbers.end(),temp_crate);
    int index = std::distance(crate_numbers.begin(), it);
    Log("MonitorTankTime: index = "+std::to_string(index),v_message,verbosity);
    switch (index) {
      case 0: {
        active_channel_cr1[temp_slot-1]=1;        //slot numbering starts at 1
        num_active_slots_cr1++;
        active_slots_cr1.push_back(temp_slot);
        break;
      }
      case 1: {
        active_channel_cr2[temp_slot-1]=1;
        num_active_slots_cr2++;
        active_slots_cr2.push_back(temp_slot);
        break;
      }
      case 2: {
        active_channel_cr3[temp_slot-1]=1;
        num_active_slots_cr3++;
        active_slots_cr3.push_back(temp_slot);
        break;
      }
    }
    std::vector<unsigned int> crateslot{temp_crate,temp_slot};
    map_slot_to_crateslot.emplace(num_active_slots,crateslot);
    map_crateslot_to_slot.emplace(crateslot,num_active_slots);
    for (unsigned int i_channel = 0; i_channel < u_num_channels_tank; i_channel++){
      std::vector<unsigned int> crateslotch{temp_crate,temp_slot,i_channel+1};
      map_ch_to_crateslotch.emplace(num_active_slots*num_channels_tank+i_channel,crateslotch);
      map_crateslotch_to_ch.emplace(crateslotch,num_active_slots*num_channels_tank+i_channel);
    }
    num_active_slots++;
  }
  file.close();
  num_active_slots = num_active_slots_cr1+num_active_slots_cr2+num_active_slots_cr3;

  Log("MonitorTankTime: Number of active Slots (Crate 1/2/3): "+std::to_string(num_active_slots_cr1)+" / "+std::to_string(num_active_slots_cr2)+" / "+std::to_string(num_active_slots_cr3),v_message,verbosity);
  Log("MonitorTankTime: Vector crate_numbers has size: "+std::to_string(crate_numbers.size()),v_debug,verbosity);
  
  //-------------------------------------------------------
  //-----------------Get RWM / BRF channels----------------
  //-------------------------------------------------------

  ifstream file_signal(signal_channels);
  std::string signal_name;
  unsigned int signal_crate, signal_slot, signal_channel;
  while(!file_signal.eof()){
	file_signal >> signal_name >> signal_crate >> signal_slot >> signal_channel;
  	if (file_signal.eof()) break;
	if (signal_name == "BRF") {
		Crate_BRF = signal_crate;
		Slot_BRF = signal_slot;
		Channel_BRF = signal_channel;
	} else if (signal_name == "RWM"){
		Crate_RWM = signal_crate;
		Slot_RWM = signal_slot;
		Channel_RWM = signal_channel;
	} else {
		Log("MonitorTankTime: Error trying to read in from signal txt file. Signal "+signal_name+"unknown. Please only specify RWM/BRF channels.",v_error,verbosity);
	}
  }
  file_signal.close();
  
  //-------------------------------------------------------
  //-----------------Get Disabled channels ----------------
  //-------------------------------------------------------

  ifstream file_disabled(disabled_channels);
  unsigned int dis_crate, dis_slot, dis_channel;
  while (!file_disabled.eof()){
  	file_disabled >> dis_crate >> dis_slot >> dis_channel;
	if (file_disabled.eof()) break;
	std::vector<unsigned int> temp_disabled{dis_crate,dis_slot,dis_channel};
	vec_disabled_channels.push_back(temp_disabled);
  }
  file_disabled.close();

  for (int i_ch=0; i_ch < num_active_slots*num_channels_tank; i_ch++){
    std::vector<unsigned int> crateslotch = map_ch_to_crateslotch[(unsigned int) i_ch];
    if (std::find(vec_disabled_channels.begin(),vec_disabled_channels.end(),crateslotch)!=vec_disabled_channels.end()) vec_disabled_global.push_back(i_ch);
  }

  //-------------------------------------------------------
  //-------------------Get geometry------------------------
  //-------------------------------------------------------

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  m_data->CStore.Get("PMTCrateSpaceToChannelNumMap",PMTCrateSpaceToChannelNumMap);
  Position center_position = geom->GetTankCentre();
  tank_center_x = center_position.X();
  tank_center_y = center_position.Y();
  tank_center_z = center_position.Z();


  //-------------------------------------------------------
  //------Setup time variables for periodic updates--------
  //-------------------------------------------------------

  Epoch = new boost::posix_time::ptime(boost::gregorian::from_string(StartTime));
  period_update = boost::posix_time::time_duration(0,update_frequency,0,0);
  last = boost::posix_time::ptime(boost::posix_time::second_clock::local_time());

  //-------------------------------------------------------
  //----------Initialize configuration/hists---------------
  //-------------------------------------------------------

  ReadInConfiguration();
  InitializeHists();
  //omit warning messages from ROOT: 1001 - info messages, 2001 - warnings, 3001 - errors
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  
  return true;
}


bool MonitorTankTime::Execute(){

  Log("Tool MonitorTankTime: Executing ....",v_message,verbosity);

  current = (boost::posix_time::second_clock::local_time());
  duration = boost::posix_time::time_duration(current - last); 
  current_stamp_duration = boost::posix_time::time_duration(current - *Epoch);
  current_stamp = current_stamp_duration.total_milliseconds();
  utc = (boost::posix_time::second_clock::universal_time());
  current_utc_duration = boost::posix_time::time_duration(utc-current);
  current_utc = current_utc_duration.total_milliseconds();
  utc_to_t = (ULong64_t) current_utc;
  utc_to_fermi = (ULong64_t) utc_to_t*MSEC_to_SEC*MSEC_to_SEC;

  //Log("MonitorTankTime: "+std::to_string(duration.total_milliseconds()/1000./60.)+" mins since last time plot",v_message,verbosity);

  //for testing purposes only: execute in every step / set State to DataFile
  //m_data->CStore.Set("State","DataFile");

  bool has_pmt;
  m_data->CStore.Get("HasPMTData",has_pmt);

  std::string State;
  m_data->CStore.Get("State",State);

  if (has_pmt){
  if (State == "Wait"){

    //-------------------------------------------------------
    //--------------No tool is executed----------------------
    //-------------------------------------------------------

    Log("MonitorTankTime: State is "+State,v_debug,verbosity);

  } else if (State == "DataFile"){

    Log("MonitorTankTime: New data file available.",v_message,verbosity);

    //-------------------------------------------------------
    //--------------MonitorTankTime executed-----------------
    //-------------------------------------------------------

    //get FinishedPMTWaves from DataDecoder tools
    m_data->CStore.Get("FinishedPMTWaves",FinishedPMTWaves);
    LoopThroughDecodedEvents(FinishedPMTWaves);

    //Write the event information to a file
    //TODO: change this to a database later on!
    //Check if data has already been written included in WriteToFile function
    WriteToFile();

    //draw last file plots
    DrawLastFilePlots();

    //Draw customly defined plots
    UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

  } else {
   	Log("MonitorTankTime: State not recognized: "+State,v_debug,verbosity);
  }
  }

  //-------------------------------------------------------------
  //---------------Draw customly defined plots-------------------
  //-------------------------------------------------------------
  
  // if force_update is specified, the plots will be updated no matter whether there has been a new file or not
 
  if (force_update) UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

  //-------------------------------------------------------
  //-----------Has enough time passed for update?----------
  //-------------------------------------------------------

  if(duration >= period_update){
    Log("MonitorTankTime: "+std::to_string(update_frequency)+" mins passed... Updating file history plot.",v_message,verbosity);

    last=current;
    DrawFileHistory(current_stamp,24.,"current_24h",1);     //show 24h history of Tank files
    PrintFileTimeStamp(current_stamp,24.,"current_24h");
    DrawFileHistory(current_stamp,2.,"current_2h",3);

  }
  
  //only for debugging memory leaks, otherwise comment out
  //std::cout <<"List of Objects (after execute step): "<<std::endl;
  //gObjectTable->Print();

  return true;

}


bool MonitorTankTime::Finalise(){

  Log("Tool MonitorTankTime: Finalising ....",v_message,verbosity);

  //delete all histograms/canvases/other objects that were created

  //help objects
  delete label_ped;
  delete label_sigma;
  delete label_rate;
  delete label_sigmadiff;
  delete label_peddiff;
  delete label_ratediff;
  delete label_fifo;
  delete label_cr1;
  delete label_cr2;
  delete label_cr3;
  delete line1;
  delete line2;
  delete text_crate1;
  delete text_crate2;
  delete text_crate3;
  delete leg_ped;
  delete leg_sigma;
  delete leg_rate;
  delete leg_freq;
  delete leg_temp;
  delete separate_crates;
  delete separate_crates2;

  for (unsigned int i_box = 0; i_box < vector_box_inactive.size(); i_box++){
    delete vector_box_inactive.at(i_box);
  }

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    for (unsigned int i_box = 0; i_box < vector_box_inactive_hitmap[i_slot].size(); i_box++){
      delete vector_box_inactive_hitmap[i_slot].at(i_box);
    }
  }

  //graphs
  for (unsigned int i_channel=0; i_channel < gr_ped.size(); i_channel++){
    delete gr_ped.at(i_channel);
    delete gr_sigma.at(i_channel);
    delete gr_rate.at(i_channel);
  }

  //multigraphs
  delete multi_ch_ped;
  delete multi_ch_sigma;
  delete multi_ch_rate;

  //histograms
  delete h2D_ped;
  delete h2D_sigma;
  delete h2D_rate; 
  delete h2D_peddiff;
  delete h2D_sigmadiff;
  delete h2D_ratediff;
  delete h2D_fifo1;
  delete h2D_fifo2;
  delete log_files;

  for (unsigned int i_channel = 0; i_channel < hChannels_temp.size(); i_channel++){
    delete hChannels_temp.at(i_channel);
    delete hChannels_freq.at(i_channel);
  }
  delete hChannels_RWM;
  delete hChannels_BRF;
  delete hChannels_temp_RWM;
  delete hChannels_temp_BRF;
  delete hist_vme;
  delete hist_vme_cluster;
  delete hist_vme_cluster_20;

  for (int i_crate = 0; i_crate < (int) hist_hitmap.size(); i_crate++){
    delete hist_hitmap.at(i_crate);
  }
  for (int i_slot=0; i_slot < (int) hist_hitmap_slot.size(); i_slot++){
    delete hist_hitmap_slot.at(i_slot);
  }

  //canvases
  delete canvas_ped;
  delete canvas_sigma;
  delete canvas_rate;
  delete canvas_peddiff;
  delete canvas_sigmadiff;
  delete canvas_ratediff;
  delete canvas_fifo;
  delete canvas_logfile_tank;
  delete canvas_ch_ped;
  delete canvas_ch_sigma;
  delete canvas_ch_rate_tank;
  delete canvas_ch_single_tank;
  delete canvas_file_timestamp_tank;
  delete canvas_vme;
  delete canvas_hitmap_tank;
  delete canvas_hitmap_tank_slot;
 

  for (unsigned int i_channel = 0; i_channel < canvas_Channels_temp.size(); i_channel++){
    delete canvas_Channels_temp.at(i_channel);
    delete canvas_Channels_freq.at(i_channel);
  }
  

  return true;

}

void MonitorTankTime::ReadInConfiguration(){

  Log("MonitorTankTime: ReadInConfiguration.",v_message,verbosity);

  //-------------------------------------------------------
  //----------------ReadInConfiguration -------------------
  //-------------------------------------------------------

  ifstream file(plot_configuration.c_str());

  std::string line;
  if (file.is_open()){
    while(std::getline(file,line)){
      if (line.find("#") != std::string::npos) continue;
      std::vector<std::string> values;
      std::stringstream ss;
        ss.str(line);
        std::string item;
        while (std::getline(ss, item, '\t')) {
            values.push_back(item);
        }
      if (values.size() < 4 ) {
        Log("ERROR (MonitorTankTime): ReadInConfiguration: Need at least 4 arguments in one line: TimeFrame - TimeEnd - FileLabel - PlotType1. Please look over the configuration file and adjust it accordingly.",v_error,verbosity);
        continue;
      }
      double double_value = std::stod(values.at(0));
      config_timeframes.push_back(double_value);
      config_endtime.push_back(values.at(1));
      config_label.push_back(values.at(2));
      std::vector<std::string> plottypes;
      for (unsigned int i=3; i < values.size(); i++){
        plottypes.push_back(values.at(i));
      }
      config_plottypes.push_back(plottypes);
    }
  } else {
    Log("ERROR (MonitorTankTime): ReadInConfiguration: Could not open file "+plot_configuration+"! Check if path is valid.",v_error,verbosity);
  }
  file.close();


  if (verbosity > 2){
    std::cout <<"---------------------------------------------------------------------"<<std::endl;
    std::cout << "MonitorTankTime: ReadInConfiguration: Read in the following data into configuration variables: "<<std::endl;
    for (unsigned int i_t=0; i_t < config_timeframes.size(); i_t++){
      std::cout <<config_timeframes.at(i_t)<<", "<<config_endtime.at(i_t)<<", "<<config_label.at(i_t)<<", ";
      for (unsigned int i_plot = 0; i_plot < config_plottypes.at(i_t).size(); i_plot++){
        std::cout <<config_plottypes.at(i_t).at(i_plot)<<", ";
      }
      std::cout<<std::endl;
    }
    std::cout <<"-----------------------------------------------------------------------"<<std::endl;
  }

  Log("MonitorTankTime: ReadInConfiguration: Parsing dates: ",v_message,verbosity);
  for (unsigned int i_date = 0; i_date < config_endtime.size(); i_date++){
    if (config_endtime.at(i_date) == "TEND_LASTFILE") {
      Log("MonitorTankTime: TEND_LASTFILE: Starting from end of last read-in file",v_message,verbosity);
      ULong64_t zero = 0; 
      config_endtime_long.push_back(zero);
    } else if (config_endtime.at(i_date).size()==15){
        boost::posix_time::ptime spec_endtime(boost::posix_time::from_iso_string(config_endtime.at(i_date)));
	boost::posix_time::time_duration spec_endtime_duration = boost::posix_time::time_duration(spec_endtime - *Epoch);
	ULong64_t spec_endtime_long = spec_endtime_duration.total_milliseconds();
	config_endtime_long.push_back(spec_endtime_long);
    } else {
      Log("MonitorTankTime: Specified end date "+config_endtime.at(i_date)+" does not have the desired format yyyymmddThhmmss. Please change the format in the config file in order to use this tool. Starting from end of last file",v_message,verbosity);
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    }
  }

}

void MonitorTankTime::InitializeHists(){

  Log("MonitorTankTime: Initialize Hists",v_message,verbosity);

  //-------------------------------------------------------
  //----------------InitializeHists -----------------------
  //-------------------------------------------------------

  std::vector<double> empty_vec;
  readfromfile_tend = 0;
  readfromfile_timeframe = 0.;


  str_ped = " Pedestal Mean (VME)";
  str_sigma = " Pedestal Sigma (VME)";
  str_rate = " Rate (VME)";
  str_peddiff = " Pedestal Difference";
  str_sigmadiff = " Sigma Difference";
  str_ratediff = " Rate Difference";
  str_fifo1 = " FIFO Overflow Error I";
  str_fifo2 = " FIFO Overflow Error II";
  crate_str = "cr";
  slot_str = "_slot";
  ch_str = "_ch";

  ss_title_ped << title_time.str() << str_ped;
  ss_title_sigma << title_time.str() << str_sigma;
  ss_title_rate << title_time.str() << str_rate;
  ss_title_peddiff << title_time.str() << str_peddiff;
  ss_title_sigmadiff << title_time.str() << str_sigmadiff;
  ss_title_ratediff << title_time.str() << str_ratediff;
  ss_title_fifo1 << title_time.str() << str_fifo1;
  ss_title_fifo2 << title_time.str() << str_fifo2;

  gROOT->cd();
  
  h2D_ped = new TH2F("h2D_ped",ss_title_ped.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);                     //Fitted gauss ADC distribution mean in 2D representation of channels, slots
  h2D_sigma = new TH2F("h2D_sigma",ss_title_sigma.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);               //Fitted gauss ADC distribution sigma in 2D representation of channels, slots
  h2D_rate = new TH2F("h2D_rate",ss_title_rate.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);                  //Rate in 2D representation of channels, slots
  h2D_peddiff = new TH2F("h2D_peddiff",ss_title_peddiff.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);         //Time difference of fitted pedestal mean values of all PMT channels
  h2D_sigmadiff = new TH2F("h2D_sigmadiff",ss_title_sigmadiff.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);   //Time difference of fitted sigma values of all PMT channels
  h2D_ratediff = new TH2F("h2D_ratediff",ss_title_ratediff.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);      //Time difference of rate values of all PMT channels
  h2D_fifo1 = new TH2F("h2D_fifo1",ss_title_fifo1.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);      //Number of FIFO Type I errors for all cards
  h2D_fifo2 = new TH2F("h2D_fifo2",ss_title_fifo2.str().c_str(),num_slots_tank,0,num_slots_tank,num_crates_tank*num_channels_tank,0,num_crates_tank*num_channels_tank);      //Number of FIFO Type II errors for all cards

  canvas_ped = new TCanvas("canvas_ped","Pedestal Mean (VME)",900,600);
  canvas_sigma = new TCanvas("canvas_sigma","Pedestal Sigma (VME)",900,600);
  canvas_rate = new TCanvas("canvas_rate","Signal Counts (VME)",900,600);
  canvas_peddiff = new TCanvas("canvas_peddiff","Pedestal Time Difference (VME)",900,600);
  canvas_sigmadiff = new TCanvas("canvas_sigmadiff","Sigma Time Difference (VME)",900,600);
  canvas_ratediff = new TCanvas("canvas_ratediff","Rate Time Difference (VME)",900,600);
  canvas_fifo = new TCanvas("canvas_fifo","FIFO Overflow Errors (VME)",900,600);
  canvas_ch_ped = new TCanvas("canvas_ch_ped","Channel Ped Canvas",900,600);
  canvas_ch_sigma = new TCanvas("canvas_ch_sigma","Channel Sigma Canvas",900,600);
  canvas_ch_rate_tank = new TCanvas("canvas_ch_rate_tank","Channel Rate Canvas",900,600);
  canvas_ch_single_tank = new TCanvas("canvas_ch_single_tank","Channel Canvas Single",900,600);
  canvas_logfile_tank = new TCanvas("canvas_logfile_tank","PMT File History",900,600); 
  canvas_file_timestamp_tank = new TCanvas("canvas_file_timestamp_tank","Timestamp Last file",900,600);
  canvas_vme = new TCanvas("canvas_vme","VME Canvas",900,600);
  canvas_hitmap_tank = new TCanvas("canvas_hitmap_tank","Hitmap VME",900,600);
  canvas_hitmap_tank_slot = new TCanvas("canvas_hitmap_tank_slot","Hitmap VME Slot",900,600);

  canvas_hitmap_tank->SetGridy();
  canvas_hitmap_tank->SetLogy();
  //canvas_hitmap_tank_slot->SetGridx();
  canvas_hitmap_tank_slot->SetLogy();
 
  int color_scheme[3] = {8,9,2};
  //Initialize hitmap histograms
  for (int i_crate = 0; i_crate < 3; i_crate++){
    //std::cout <<"crate: "<<i_crate<<std::endl;
    std::stringstream ss_hitmap_name, ss_hitmap_title;
    ss_hitmap_name <<"hist_hitmap_tank_"<<i_crate;
    ss_hitmap_title<<"Hitmap VME Crate "<<i_crate;
    TH1F *hist_hitmap_single = new TH1F(ss_hitmap_name.str().c_str(),ss_hitmap_title.str().c_str(),num_channels_tank*num_active_slots,0,num_channels_tank*num_active_slots);
    hist_hitmap_single->SetLineColor(color_scheme[i_crate]);
    hist_hitmap_single->SetFillColor(color_scheme[i_crate]);
    hist_hitmap_single->SetStats(0);
    hist_hitmap_single->GetYaxis()->SetTitle("# of entries");
    hist_hitmap_single->GetXaxis()->SetNdivisions(int(num_active_slots));
    hist_hitmap.push_back(hist_hitmap_single);
   }

   //std::cout <<"num_active_slots: "<<num_active_slots<<std::endl;
   for (int i_label=0; i_label < int(num_active_slots); i_label++){
    //std::cout <<"i_label: "<<i_label<<std::endl;
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_label];
    unsigned int i_crate = crateslot.at(0)-1;
    unsigned int slot = crateslot.at(1);
    //std::cout <<"crate: "<<i_crate<<std::endl;
    std::stringstream ss_slot;
    ss_slot << "slot "<< slot;
    for (int i_cr=0; i_cr < 3; i_cr++){
      hist_hitmap[i_cr]->GetXaxis()->SetBinLabel(num_channels_tank*(i_label+0.5),ss_slot.str().c_str());
    }
  }
  //std::cout <<"blub"<<std::endl;
  hist_hitmap[0]->SetTitleSize(0.3,"t");
  for (int i_crate = 0; i_crate < 3; i_crate++){
    //std::cout <<"crate: "<<i_crate<<std::endl;
    hist_hitmap[i_crate]->LabelsOption("v");
    hist_hitmap[i_crate]->SetTickLength(0,"X");
    hist_hitmap[i_crate]->SetLineWidth(2);
  }
  label_cr1 = new TPaveText(0.12,0.83,0.27,0.88,"NDC");
  label_cr1->SetFillColor(0);
  label_cr1->SetTextColor(color_scheme[0]);
  label_cr1->AddText("VME01");
  label_cr2 = new TPaveText(0.43,0.83,0.57,0.88,"NDC");
  label_cr2->SetFillColor(0);
  label_cr2->SetTextColor(color_scheme[1]);
  label_cr2->AddText("VME02");
  label_cr3 = new TPaveText(0.73,0.83,0.88,0.88,"NDC");
  label_cr3->SetFillColor(0);
  label_cr3->SetTextColor(color_scheme[2]);
  label_cr3->AddText("VME03");

  //std::cout <<"num_active_slots_cr1: "<<num_active_slots_cr1<<std::endl;
  separate_crates = new TLine(num_channels_tank*num_active_slots_cr1,0.8,num_channels_tank*num_active_slots_cr1,100);
  separate_crates->SetLineStyle(2);
  separate_crates->SetLineWidth(2);
  separate_crates2 = new TLine(num_channels_tank*(num_active_slots_cr1+num_active_slots_cr2),0.8,num_channels_tank*(num_active_slots_cr1+num_active_slots_cr2),100);
  separate_crates2->SetLineStyle(2);
  separate_crates2->SetLineWidth(2);
  f1 = new TF1("f1","x",0,num_active_slots);      //workaround to only have labels for every slot

  //-------------------------------------------------------
  //------------Initialize single slot hitmaps-------------
  //-------------------------------------------------------
  
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    //std::cout <<"i_slot: "<<i_slot<<std::endl;
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int crate_num = crateslot.at(0);
    unsigned int slot_num = crateslot.at(1);

    std::stringstream ss_title_single, ss_name_single;
    ss_name_single << "hist_hitmap_tank_cr"<<crate_num<<"_sl"<<slot_num;
    ss_title_single << "Hitmap VME Crate "<<crate_num<<" Slot "<<slot_num;

    int hist_color = color_scheme[crate_num-1];

    TH1F *hist_hitmap_slot_single = new TH1F(ss_name_single.str().c_str(),ss_title_single.str().c_str(),num_channels_tank,0,num_channels_tank);
    hist_hitmap_slot_single->SetFillColor(hist_color);
    hist_hitmap_slot_single->SetLineColor(hist_color);
    hist_hitmap_slot_single->SetLineWidth(2);
    hist_hitmap_slot_single->SetStats(0);
    hist_hitmap_slot_single->GetYaxis()->SetTitle("# of entries");
    hist_hitmap_slot_single->GetXaxis()->SetTitle("Channel #");
    hist_hitmap_slot.push_back(hist_hitmap_slot_single);

  }

  for (int i_active = 0; i_active<num_active_slots; i_active++){
    
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_active];
    unsigned int crate_num = crateslot.at(0);
    unsigned int slot_num = crateslot.at(1);

    std::stringstream ss_name_hist;
    std::stringstream ss_name_hist_temp;
    std::stringstream ss_title_hist;
    std::stringstream ss_name_canvas;
    std::stringstream ss_name_canvas_temp;
    std::stringstream ss_title_canvas;
    std::string crate_str="cr";
    std::string slot_str = "_slot";
    std::string ch_str = "_ch";
    std::string canvas_str = "canvas_";

    ss_name_canvas <<canvas_str<<crate_str<<crate_num<<slot_str<<slot_num;
    ss_name_canvas_temp << ss_name_canvas.str() << "_temp";
    ss_title_canvas << "Crate" <<crate_str<<"_Slot"<<slot_str;
    TCanvas *canvas_Channel_temp = new TCanvas(ss_name_canvas_temp.str().c_str(),ss_title_canvas.str().c_str(),900,600);
    canvas_Channels_temp.push_back(canvas_Channel_temp);
    TCanvas *canvas_Channel = new TCanvas(ss_name_canvas.str().c_str(),ss_title_canvas.str().c_str(),900,600);
    canvas_Channels_freq.push_back(canvas_Channel);

    //create frequency + buffer histograms
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      ss_name_hist.str("");
      ss_title_hist.str("");
      ss_name_hist<<crate_str<<crate_num<<slot_str<<slot_num<<ch_str<<i_channel;
      ss_title_hist << title_time.str() << " ADC Freq (VME Crate " << crate_num << " Slot " << slot_num <<")";
      
      TH1I *hChannel_freq = new TH1I(ss_name_hist.str().c_str(),ss_title_hist.str().c_str(),maximum_adc-minimum_adc,minimum_adc,maximum_adc);   //ADC pulse shapes for current event, 1 canvas per slot
      hChannel_freq->GetXaxis()->SetTitle("ADC");
      hChannel_freq->GetYaxis()->SetTitle("Counts");
      hChannel_freq->SetLineWidth(2);
      hChannel_freq->SetLineColor(i_channel+1);
      hChannel_freq->SetStats(0);
      hChannel_freq->GetYaxis()->SetTitleOffset(1.35);
      hChannels_freq.push_back(hChannel_freq);

      ss_name_hist_temp.str("");
      ss_title_hist_temp.str("");
      ss_name_hist_temp<<crate_str<<crate_num<<slot_str<<slot_num<<ch_str<<i_channel<<"_temp";
      ss_title_hist_temp << title_time.str() << " Buffer Temp (VME Crate " << crate_num << " Slot " << slot_num <<")";
      TH1F* hChannel_temp = new TH1F(ss_name_hist_temp.str().c_str(),ss_title_hist_temp.str().c_str(),4000,0,4000);   //default buffer size for initializing, will be readjusted according to read-out size
      hChannel_temp->GetXaxis()->SetTitle("Buffer Position");
      hChannel_temp->GetYaxis()->SetTitle("Volts");
      hChannel_temp->SetLineWidth(2);
      hChannel_temp->SetLineColor(i_channel+1);
      hChannel_temp->GetYaxis()->SetTitleOffset(1.35);
      hChannel_temp->SetStats(0);
      hChannels_temp.push_back(hChannel_temp);

      channels_mean.push_back(0);                             //starting mean/Sigma values for all channels set to 0, will be filled with the fit values from Gaussian fit
      channels_sigma.push_back(0);

    }
  }

  canvas_Channels_freq.at(num_active_slots-1)->Clear();       //this canvas gets otherwise drawn with an additional frequency histogram (last channel)

  //initialize BRF/RWM histograms

  hChannels_BRF = new TH1I("hChannels_BRF","ADC Freq BRF",4000,0,4000);
  hChannels_BRF->GetXaxis()->SetTitle("ADC");
  hChannels_BRF->GetYaxis()->SetTitle("Counts");
  hChannels_BRF->SetLineWidth(2);
  hChannels_BRF->SetStats(0);
  hChannels_BRF->GetYaxis()->SetTitleOffset(1.35);
   
  hChannels_RWM = new TH1I("hChannels_RWM","ADC Freq RWM",4000,0,4000);
  hChannels_RWM->GetXaxis()->SetTitle("ADC");
  hChannels_RWM->GetYaxis()->SetTitle("Counts");
  hChannels_RWM->SetLineWidth(2);
  hChannels_RWM->SetStats(0);
  hChannels_RWM->GetYaxis()->SetTitleOffset(1.35);

  hChannels_temp_BRF = new TH1F("hChannels_temp_BRF","Buffer Temp BRF",2000,0,2000);
  hChannels_temp_BRF->GetXaxis()->SetTitle("Buffer Position");
  hChannels_temp_BRF->GetYaxis()->SetTitle("Volts");
  hChannels_temp_BRF->SetLineWidth(2);
  hChannels_temp_BRF->GetYaxis()->SetTitleOffset(1.35);
  hChannels_temp_BRF->SetStats(0);

  hChannels_temp_RWM = new TH1F("hChannels_temp_RWM","Buffer Temp RWM",2000,0,2000);
  hChannels_temp_RWM->GetXaxis()->SetTitle("Buffer Position");
  hChannels_temp_RWM->GetYaxis()->SetTitle("Volts");
  hChannels_temp_RWM->SetLineWidth(2);
  hChannels_temp_RWM->GetYaxis()->SetTitleOffset(1.35);
  hChannels_temp_RWM->SetStats(0);

  //initialize file history histogram and canvas
  num_files_history = 10;
  log_files = new TH1F("log_files","PMT Files History",num_files_history,0,num_files_history);
  log_files->GetXaxis()->SetTimeDisplay(1);
  log_files->GetXaxis()->SetLabelSize(0.03);
  log_files->GetXaxis()->SetLabelOffset(0.03);
  log_files->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  log_files->GetYaxis()->SetTickLength(0.);
  log_files->GetYaxis()->SetLabelOffset(999);
  log_files->SetStats(0);

  //Clustered time histograms
  hist_vme = new TH1F("hist_vme","VME (last file)",100,0,2000);
  hist_vme->GetXaxis()->SetTitle("Time [ns]");
  hist_vme->GetYaxis()->SetTitle("#");
  hist_vme_cluster = new TH1F("hist_vme_cluster","VME Cluster (last file)",20,0,2000);
  hist_vme_cluster->GetXaxis()->SetTitle("Time [ns]");
  hist_vme_cluster->GetYaxis()->SetTitle("#");
  hist_vme_cluster_20 = new TH1F("hist_vme_cluster_20","VME Cluster (last 20 files)",20,0,2000);
  hist_vme_cluster_20->GetXaxis()->SetTitle("Time [ns]");
  hist_vme_cluster_20->GetYaxis()->SetTitle("#");

  std::vector<TBox*> empty_boxvector;
  for (int i_slot=0; i_slot < num_active_slots; i_slot++){
    vector_box_inactive_hitmap.emplace(i_slot,empty_boxvector);
  }

  //initialize vector inactive boxes
  for (int i_slot = 0; i_slot < num_slots_tank*num_crates_tank; i_slot++){
    int slot_in_crate;
    int crate_nr;
    if (i_slot < num_slots_tank) {slot_in_crate = i_slot; crate_nr=0;}
    else if (i_slot < 2*num_slots_tank) {slot_in_crate = i_slot -num_slots_tank; crate_nr=1;}
    else {slot_in_crate = i_slot - 2*num_slots_tank; crate_nr=2;}

    if (active_channel_cr1[slot_in_crate] == 0 && crate_nr==0) {
      TBox *box_inactive = new TBox(slot_in_crate,2*num_channels_tank,slot_in_crate+1,3*num_channels_tank);
      box_inactive->SetFillStyle(3004);
      box_inactive->SetFillColor(1);
      vector_box_inactive.push_back(box_inactive);
      for (int i_ch = 0; i_ch < num_channels_tank; i_ch++){
        std::vector<int> inactive_slotch{slot_in_crate+1,2*num_channels_tank+i_ch+1};
        inactive_xy.push_back(inactive_slotch);
      }
    }
    if (active_channel_cr2[slot_in_crate] == 0 && crate_nr==1) {
      TBox *box_inactive = new TBox(slot_in_crate,num_channels_tank,slot_in_crate+1,2*num_channels_tank);
      box_inactive->SetFillStyle(3004);
      box_inactive->SetFillColor(1);
      vector_box_inactive.push_back(box_inactive);
      for (int i_ch = 0; i_ch < num_channels_tank; i_ch++){
        std::vector<int> inactive_slotch{slot_in_crate+1,num_channels_tank+i_ch+1};
        inactive_xy.push_back(inactive_slotch);
      }
    }
    if (active_channel_cr3[slot_in_crate] == 0 && crate_nr==2) {
      TBox *box_inactive = new TBox(slot_in_crate,0,slot_in_crate+1,num_channels_tank);
      box_inactive->SetFillStyle(3004);
      box_inactive->SetFillColor(1);
      vector_box_inactive.push_back(box_inactive);
      for (int i_ch = 0; i_ch < num_channels_tank; i_ch++){
        std::vector<int> inactive_slotch{slot_in_crate+1,i_ch+1};
        inactive_xy.push_back(inactive_slotch);
      }
    }
  }
  for (unsigned int i_disabled = 0; i_disabled < vec_disabled_channels.size(); i_disabled++){
    unsigned int crate = vec_disabled_channels.at(i_disabled).at(0);
    unsigned int slot = vec_disabled_channels.at(i_disabled).at(1);
    unsigned int ch = vec_disabled_channels.at(i_disabled).at(2);
    unsigned int global_ch = 4 - ch;
    if (crate == 1) global_ch = 2*num_channels_tank + 4 - ch;
    else if (crate == 2) global_ch = num_channels_tank + 4 - ch;
    else if (crate == 3) global_ch = 4 - ch;
    TBox *box_inactive = new TBox(slot-1,global_ch,slot,global_ch+1);
    box_inactive->SetFillStyle(3004);
    box_inactive->SetFillColor(2);
    vector_box_inactive.push_back(box_inactive);
    TBox *box_inactive_hitmap = new TBox(ch-1,0.8,ch,1.);
    box_inactive_hitmap->SetFillStyle(3004);
    box_inactive_hitmap->SetFillColor(1);
    std::vector<unsigned int> crateslot{crate,slot};
    int active_slot = map_crateslot_to_slot[crateslot];
    vector_box_inactive_hitmap[active_slot].push_back(box_inactive_hitmap);
  }

  //initialize labels/text boxes/etc
  label_rate = new TLatex(0.905,0.92,"Rate [kHz]");
  label_rate->SetNDC(1);
  label_rate->SetTextSize(0.030);
  label_ped = new TLatex(0.905,0.92,"#mu_{Ped}");
  label_ped->SetNDC(1);
  label_ped->SetTextSize(0.030);
  label_sigma = new TLatex(0.905,0.92,"#sigma_{Ped}");
  label_sigma->SetNDC(1);
  label_sigma->SetTextSize(0.030);
  label_peddiff = new TLatex(0.905,0.92,"#Delta #mu_{Ped}");
  label_peddiff->SetNDC(1);
  label_peddiff->SetTextSize(0.030);
  label_sigmadiff = new TLatex(0.905,0.92,"#Delta #sigma_{Ped}");
  label_sigmadiff->SetNDC(1);
  label_sigmadiff->SetTextSize(0.030);
  label_ratediff = new TLatex(0.905,0.92,"#Delta R [kHz]");
  label_ratediff->SetNDC(1);
  label_ratediff->SetTextSize(0.030);
  label_fifo = new TLatex(0.905,0.92,"#");
  label_fifo->SetNDC(1);
  label_fifo->SetTextSize(0.030);
  std::stringstream ss_crate1, ss_crate2, ss_crate3;
  ss_crate1 << "ANNIEVME0"<<crate_numbers.at(0);
  ss_crate2 << "ANNIEVME0"<<crate_numbers.at(1);
  ss_crate3 << "ANNIEVME0"<<crate_numbers.at(2);
  text_crate1 = new TText(0.04,0.68,ss_crate1.str().c_str());
  text_crate1->SetNDC(1);
  text_crate1->SetTextSize(0.030);
  text_crate1->SetTextAngle(90.);
  text_crate2 = new TText(0.04,0.41,ss_crate2.str().c_str());
  text_crate2->SetNDC(1);
  text_crate2->SetTextSize(0.030);
  text_crate2->SetTextAngle(90.);
  text_crate3 = new TText(0.04,0.14,ss_crate3.str().c_str());
  text_crate3->SetNDC(1);
  text_crate3->SetTextSize(0.030);
  text_crate3->SetTextAngle(90.);
  line1 = new TLine(-1,4,num_slots_tank,4.);
  line1->SetLineWidth(2);
  line2 = new TLine(-1,8.,num_slots_tank,8.);
  line2->SetLineWidth(2);
  leg_freq = new TLegend(0.75,0.7,0.9,0.9);
  leg_temp = new TLegend(0.75,0.7,0.9,0.9);


  // Initialize time evolution graphs

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){

    std::stringstream ss_crate, ss_slot;

    unsigned int crate_num = map_slot_to_crateslot[i_slot].at(0);
    unsigned int slot_num = map_slot_to_crateslot[i_slot].at(1);
    ss_crate << crate_num;
    ss_slot << slot_num;

    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){

      std::stringstream ss_ch;
      ss_ch << i_channel+1;
      int line_color_ch = i_channel+1;

      //ped_file.push_back(empty_vec);

      //channel-wise graphs
      std::stringstream name_graph, name_graph_ped, name_graph_sigma, name_graph_rate;
      std::stringstream title_graph, title_graph_ped, title_graph_sigma, title_graph_rate;

      name_graph << "cr"<<ss_crate.str()<<"_sl"<<ss_slot.str()<<"_ch"<<ss_ch.str();
      name_graph_ped << name_graph.str()<<"_ped";
      name_graph_sigma << name_graph.str()<<"_sigma";
      name_graph_rate << name_graph.str()<<"_rate";
      title_graph << "Crate "<<ss_crate.str()<<", Slot "<<ss_slot.str()<<", Channel "<<ss_ch.str();
      title_graph_ped << title_graph.str() << " (Pedestal)";
      title_graph_sigma << title_graph.str() <<" (Ped Sigma)";
      title_graph_rate << title_graph.str() <<" (Rate)";

      TGraph *graph_ch_ped = new TGraph();
      TGraph *graph_ch_sigma = new TGraph();
      TGraph *graph_ch_rate = new TGraph();

      graph_ch_ped->SetName(name_graph_ped.str().c_str());
      graph_ch_ped->SetTitle(title_graph_ped.str().c_str());
      graph_ch_sigma->SetName(name_graph_sigma.str().c_str());
      graph_ch_sigma->SetTitle(title_graph_sigma.str().c_str());
      graph_ch_rate->SetName(name_graph_rate.str().c_str());
      graph_ch_rate->SetTitle(title_graph_rate.str().c_str());
      
      if (draw_marker) {
        graph_ch_ped->SetMarkerStyle(20);
        graph_ch_sigma->SetMarkerStyle(20);
        graph_ch_rate->SetMarkerStyle(20);
      }

      graph_ch_ped->SetMarkerColor(line_color_ch);
      graph_ch_sigma->SetMarkerColor(line_color_ch);
      graph_ch_rate->SetMarkerColor(line_color_ch);
      graph_ch_sigma->SetMarkerColor(line_color_ch);
      graph_ch_rate->SetMarkerColor(line_color_ch);
      graph_ch_ped->SetLineColor(line_color_ch);
      graph_ch_ped->SetLineWidth(2);
      graph_ch_ped->SetFillColor(0);
      graph_ch_ped->GetYaxis()->SetTitle("Pedestal");
      graph_ch_ped->GetXaxis()->SetTimeDisplay(1);
      graph_ch_ped->GetXaxis()->SetLabelSize(0.03);
      graph_ch_ped->GetXaxis()->SetLabelOffset(0.03);
      graph_ch_ped->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      graph_ch_sigma->SetLineColor(line_color_ch);
      graph_ch_sigma->SetLineWidth(2);
      graph_ch_sigma->SetFillColor(0);
      graph_ch_sigma->GetYaxis()->SetTitle("Sigma (Ped)");
      graph_ch_sigma->GetXaxis()->SetTimeDisplay(1);
      graph_ch_sigma->GetXaxis()->SetLabelSize(0.03);
      graph_ch_sigma->GetXaxis()->SetLabelOffset(0.03);
      graph_ch_sigma->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      graph_ch_rate->SetLineColor(line_color_ch);
      graph_ch_rate->SetLineWidth(2);
      graph_ch_rate->SetFillColor(0);
      graph_ch_rate->GetYaxis()->SetTitle("Rate [kHz]");
      graph_ch_rate->GetXaxis()->SetTimeDisplay(1);
      graph_ch_rate->GetXaxis()->SetLabelSize(0.03);
      graph_ch_rate->GetXaxis()->SetLabelOffset(0.03);
      graph_ch_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");

      gr_ped.push_back(graph_ch_ped);
      gr_sigma.push_back(graph_ch_sigma);
      gr_rate.push_back(graph_ch_rate);
    }
  }

  multi_ch_ped = new TMultiGraph();
  multi_ch_sigma = new TMultiGraph();
  multi_ch_rate = new TMultiGraph();
  leg_ped = new TLegend(0.7,0.7,0.88,0.88);
  leg_sigma = new TLegend(0.7,0.7,0.88,0.88);
  leg_rate = new TLegend(0.7,0.7,0.88,0.88);
  leg_ped->SetLineColor(0);
  leg_sigma->SetLineColor(0);
  leg_rate->SetLineColor(0);

}

void MonitorTankTime::LoopThroughDecodedEvents(std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t>>> finishedPMTWaves){

  Log("MonitorTankTime: LoopThroughDecodedEvents",v_message,verbosity);

  //-------------------------------------------------------
  //-------------LoopThroughDecodedEvents -----------------
  //-------------------------------------------------------

  ped_file.clear();
  sigma_file.clear();
  rate_file.clear();
  samples_file.clear();
  timestamp_file.clear();
  channels_times.clear();

  bool first_rwm_entry=true;
  bool first_brf_entry=true;
  int num_samples_first_rwm=0;
  int num_samples_first_brf=0;

  int i_timestamp = 0;
  for (std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t>>>::iterator it = finishedPMTWaves.begin(); it != finishedPMTWaves.end(); it++){

    uint64_t timestamp = it->first;
    uint64_t timestamp_temp = timestamp - utc_to_fermi;			//conversion from UTC time to Fermilab US time
    timestamp_file.push_back(timestamp_temp);
    std::vector<double> ped_file_temp, sigma_file_temp, rate_file_temp;
    std::vector<int> samples_file_temp;
    std::vector<std::vector<int>> channels_pmt_time;
    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      std::vector<int> pmt_hittimes;
      channels_pmt_time.push_back(pmt_hittimes);
    }
    ped_file_temp.assign(num_active_slots*num_channels_tank,0.);
    sigma_file_temp.assign(num_active_slots*num_channels_tank,0.);
    rate_file_temp.assign(num_active_slots*num_channels_tank,0.);
    samples_file_temp.assign(num_active_slots*num_channels_tank,0);
    channels_rate.assign(num_active_slots*num_channels_tank,0.);
    channels_mean.assign(num_active_slots*num_channels_tank,0.);
    channels_sigma.assign(num_active_slots*num_channels_tank,0.);

    std::map<std::vector<int>, std::vector<uint16_t>> afinishedPMTWaves = finishedPMTWaves.at(timestamp);
    for(std::pair<std::vector<int>, std::vector<uint16_t>> apair : afinishedPMTWaves){

      int CardID = apair.first.at(0);
      int ChannelID = apair.first.at(1);
      std::vector<uint16_t> awaveform = apair.second;
      int num_samples = int(awaveform.size()) - 50;
      int CrateNum, SlotNum;
      this->CardIDToElectronicsSpace(CardID, CrateNum, SlotNum);
      unsigned int uCrateNum = (unsigned int) CrateNum;
      unsigned int uSlotNum = (unsigned int) SlotNum;
      unsigned int uChannelID = (unsigned int) ChannelID;
      std::vector<unsigned int> CrateSlot{uCrateNum,uSlotNum};
      std::vector<unsigned int> CrateSlotCh{uCrateNum,uSlotNum,uChannelID};
      //check if read out crate/slot configuration is active according to configuration file...
      if (map_crateslot_to_slot.find(CrateSlot) == map_crateslot_to_slot.end()){
        if (uCrateNum == Crate_BRF && uSlotNum == Slot_BRF){    //Assume that BRF and RWM are in same VME slot
          if (ChannelID == int(Channel_BRF)) {
		hChannels_BRF->Reset();
		for (int i_buffer =0; i_buffer < num_samples; i_buffer++){
			hChannels_BRF->Fill(awaveform.at(i_buffer));
		}
		double BRF_mean = hChannels_BRF->GetMean();
		double BRF_sigma = hChannels_BRF->GetRMS();
		//hChannels_temp_BRF->SetBins(num_samples,0,num_samples);
		if (first_brf_entry){
			hChannels_temp_BRF->SetBins(num_samples,0,num_samples);
			num_samples_first_brf=num_samples;
			first_brf_entry=false;
		}
		int samples = (num_samples > num_samples_first_brf)? num_samples_first_brf : num_samples;
		for (int i_buffer=0; i_buffer < samples; i_buffer++){
			hChannels_temp_BRF->SetBinContent(i_buffer,hChannels_temp_BRF->GetBinContent(i_buffer)+(awaveform.at(i_buffer)-BRF_mean)*conversion_ADC_Volt);
			//hChannels_temp_BRF->SetBinContent(i_buffer,(awaveform.at(i_buffer)-BRF_mean)*conversion_ADC_Volt);
		}
	  }
          else if (ChannelID == int(Channel_RWM)) {
		hChannels_RWM->Reset();
		for (int i_buffer =0; i_buffer < num_samples; i_buffer++){
			hChannels_RWM->Fill(awaveform.at(i_buffer));
		}
		double RWM_mean = hChannels_RWM->GetMean();
		double RWM_sigma = hChannels_RWM->GetRMS();
                //hChannels_temp_RWM->SetBins(num_samples,0,num_samples);
		if (first_rwm_entry){
			hChannels_temp_RWM->SetBins(num_samples,0,num_samples);
			num_samples_first_rwm = num_samples;
			first_rwm_entry=false;
		}
		int samples = (num_samples > num_samples_first_rwm)? num_samples_first_rwm : num_samples;
		for (int i_buffer=0; i_buffer < samples; i_buffer++){
			hChannels_temp_RWM->SetBinContent(i_buffer,hChannels_temp_RWM->GetBinContent(i_buffer)+(awaveform.at(i_buffer)-RWM_mean)*conversion_ADC_Volt);
			//hChannels_temp_RWM->SetBinContent(i_buffer,(awaveform.at(i_buffer)-RWM_mean)*conversion_ADC_Volt);
		}
	}
        }
	else {
        	Log("MonitorTankTime ERROR: Slot read out from data (Cr"+std::to_string(CrateNum)+"/Sl"+std::to_string(SlotNum)+") should not be active according to config file. Check config file...",v_error,verbosity);
        }
      } else if (std::find(vec_disabled_channels.begin(),vec_disabled_channels.end(),CrateSlotCh)!=vec_disabled_channels.end()){
        
        int i_slot = map_crateslot_to_slot[CrateSlot];
        int i_channel = i_slot*num_channels_tank+ChannelID-1;
        hChannels_freq.at(i_channel)->Reset();
        hChannels_temp.at(i_channel)->Reset();
        hChannels_temp.at(i_channel)->SetBins(num_samples,0,num_samples);
        ped_file_temp.at(i_channel) = 0.;
        sigma_file_temp.at(i_channel) = 0.;
        rate_file_temp.at(i_channel) = 0.;
        
        //Don't plot the information for disabled channels
	continue;

      } else {
      int i_slot = map_crateslot_to_slot[CrateSlot];
      int i_channel = i_slot*num_channels_tank+ChannelID-1;
      hChannels_freq.at(i_channel)->Reset();            //only show the most recent plot for each PMT
      hChannels_temp.at(i_channel)->Reset();            //only show the most recent plot for each PMT
      hChannels_temp.at(i_channel)->SetBins(num_samples,0,num_samples);

      //Fill frequency histograms
      for (int i_buffer = 0; i_buffer < num_samples; i_buffer++){
        hChannels_freq.at(i_channel)->Fill(awaveform.at(i_buffer));
      }

      //fit pedestal values with Gaussian
      TF1 *fgaus = new TF1("fgaus","gaus",minimum_adc,maximum_adc);
      fgaus->SetParameter(1,hChannels_freq.at(i_channel)->GetMean());
      fgaus->SetParameter(2,hChannels_freq.at(i_channel)->GetRMS());
      TFitResultPtr gaussFitResult = hChannels_freq.at(i_channel)->Fit("fgaus","Q");
      Int_t gaussFitResultInt = gaussFitResult;
      if (gaussFitResultInt == 0){            //status variable 0 means the fit was ok
        //TF1 *gaus = (TF1*) hChannels_freq.at(i_channel)->GetFunction("gaus");
        //std::stringstream ss_gaus;
        //ss_gaus<<"gaus_"<<i_timestamp<<"_"<<i_channel;
        //gaus->SetName(ss_gaus.str().c_str());
        bool out_of_bounds = ((fgaus->GetParameter(1) < 300.) || (fgaus->GetParameter(1) > 400.) ||  (channels_sigma.at(i_channel) > 5.) || (channels_sigma.at(i_channel) < 0.5));
        bool sudden_change = (fabs(fgaus->GetParameter(1) - channels_mean.at(i_channel)) > 10 || fabs(fgaus->GetParameter(2)-channels_sigma.at(i_channel)) > 0.2);
        if (out_of_bounds || sudden_change) {	//if fit results are unphysical OR indicate a bad fit, use RMS & Mean instead
	  channels_mean.at(i_channel) = hChannels_freq.at(i_channel)->GetMean();
	  channels_sigma.at(i_channel) = hChannels_freq.at(i_channel)->GetRMS();
	} else {
          channels_mean.at(i_channel) = fgaus->GetParameter(1);
          channels_sigma.at(i_channel) = fgaus->GetParameter(2);
        }
      }else {     //if fit failed, use RMS & Mean instead
        channels_mean.at(i_channel) = hChannels_freq.at(i_channel)->GetMean();
	channels_sigma.at(i_channel) = hChannels_freq.at(i_channel)->GetRMS();
      }

      delete fgaus;

      //fill buffer plots
      for (int i_buffer = 0; i_buffer < num_samples; i_buffer++){
        hChannels_temp.at(i_channel)->SetBinContent(i_buffer,(awaveform.at(i_buffer)-channels_mean.at(i_channel))*conversion_ADC_Volt);
      }
      
      //Evaluate rates
      long sum = 0;
      std::vector<int> pmt_time;
      for (int i_buffer = 0; i_buffer < num_samples; i_buffer++){
        if (awaveform.at(i_buffer) > channels_mean[i_channel]+5*channels_sigma[i_channel]) {
		Log("MonitorTankTime tool: Found waveform entry > sigma: waveform = "+std::to_string(awaveform.at(i_buffer))+", mean+5sigma = "+std::to_string(channels_mean[i_channel]+5*channels_sigma[i_channel]),v_debug,verbosity);
		sum++;
		channels_pmt_time.at(i_channel).push_back(i_buffer*2);
		hist_vme->Fill(i_buffer*2);
      	}
      }
      channels_rate.at(i_channel) = sum;	//actually this is just the number of signal counts, convert to a rate later on

      ped_file_temp.at(i_channel) = channels_mean.at(i_channel);
      sigma_file_temp.at(i_channel) = channels_sigma.at(i_channel);
      rate_file_temp.at(i_channel) = channels_rate.at(i_channel);
      samples_file_temp.at(i_channel) = num_samples;
      }
      
    }

    channels_times.push_back(channels_pmt_time); 
    ped_file.push_back(ped_file_temp);
    sigma_file.push_back(sigma_file_temp);
    rate_file.push_back(rate_file_temp);
    samples_file.push_back(samples_file_temp);
    i_timestamp++;   
 
  }

  fifo1.clear();
  fifo2.clear();
  m_data->CStore.Get("FIFOError1",fifo1);
  m_data->CStore.Get("FIFOError2",fifo2);


}

void MonitorTankTime::WriteToFile(){

  Log("MonitorTankTime: WriteToFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------WriteToFile -------------------------
  //-------------------------------------------------------

 // t_file_start = timestamp_file.at(0)*8/1000000.;	//conversion from clock ticks to UTC in msec
 // t_file_end = timestamp_file.at(timestamp_file.size()-1)*8/1000000.;       //conversion from clock ticks to UTC in msec
  t_file_start = timestamp_file.at(0)/1000000.;	//conversion from ns to UTC in msec
  t_file_end = timestamp_file.at(timestamp_file.size()-1)/1000000.;       //conversion from ns to UTC in msec
  std::string file_start_date = convertTimeStamp_to_Date(t_file_start);
  std::stringstream root_filename;
  root_filename << path_monitoring << "PMT_" << file_start_date <<".root";
  
  Log("MonitorTankTime: ROOT filename: "+root_filename.str(),v_message,verbosity);

  std::string root_option = "RECREATE";
  if (does_file_exist(root_filename.str())) root_option = "UPDATE";
  TFile *f = new TFile(root_filename.str().c_str(),root_option.c_str());

  ULong64_t t_start, t_end, t_frame;
  std::vector<unsigned int> *crate = new std::vector<unsigned int>;
  std::vector<unsigned int> *slot = new std::vector<unsigned int>;
  std::vector<unsigned int> *channel = new std::vector<unsigned int>;
  std::vector<double> *ped = new std::vector<double>;
  std::vector<double> *sigma = new std::vector<double>;
  std::vector<double> *rate = new std::vector<double>;
  std::vector<int> *channelcount = new std::vector<int>;

  TTree *t;
  if (f->GetListOfKeys()->Contains("tankmonitor_tree")) {
    Log("MonitorTankTime: WriteToFile: Tree already exists",v_message,verbosity);
    t = (TTree*) f->Get("tankmonitor_tree");
    t->SetBranchAddress("t_start",&t_start);
    t->SetBranchAddress("t_end",&t_end);
    t->SetBranchAddress("crate",&crate);
    t->SetBranchAddress("slot",&slot);
    t->SetBranchAddress("channel",&channel);
    t->SetBranchAddress("ped",&ped);
    t->SetBranchAddress("sigma",&sigma);
    t->SetBranchAddress("rate",&rate);
    t->SetBranchAddress("channelcount",&channelcount);
  } else {
    t = new TTree("tankmonitor_tree","Tank Monitoring tree");
    Log("MonitorTankTime: WriteToFile: Tree is created from scratch",v_message,verbosity);
    t->Branch("t_start",&t_start);
    t->Branch("t_end",&t_end);
    t->Branch("crate",&crate);
    t->Branch("slot",&slot);
    t->Branch("channel",&channel);
    t->Branch("ped",&ped);
    t->Branch("sigma",&sigma);
    t->Branch("rate",&rate);
    t->Branch("channelcount",&channelcount);
  }

  int n_entries = t->GetEntries();
  bool omit_entries = false;
  for (int i_entry = 0; i_entry < n_entries; i_entry++){
    t->GetEntry(i_entry);
    if (t_start == t_file_start) {
      Log("WARNING (MonitorTankTime): WriteToFile: Wanted to write data from file that is already written to DB. Omit entries",v_warning,verbosity);
      omit_entries = true;
    }
  }
  //if data is already written to DB/File, do not write it again
  if (omit_entries) {

  //don't write file again, but still delete TFile and TTree object!!!
  f->Close();
  delete crate;
  delete slot;
  delete channel;
  delete ped;
  delete sigma;
  delete rate;
  delete channelcount;
  delete f;

  gROOT->cd();

  return;

  } 

  crate->clear();
  slot->clear();
  channel->clear();
  ped->clear();
  sigma->clear();
  rate->clear();
  channelcount->clear();
  
  t_start = t_file_start;
  t_end = t_file_end;
  t_frame = t_end - t_start;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(t_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_start/MSEC_to_SEC/SEC_to_MIN)%60,int(t_start/MSEC_to_SEC/1000.)%60,t_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  Log("MonitorTankTime: WriteToFile: Writing data to file: "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
    +"..."+std::to_string(endtime_tm.tm_year+1900)+"/"+std::to_string(endtime_tm.tm_mon+1)+"/"+std::to_string(endtime_tm.tm_mday)+"-"+std::to_string(endtime_tm.tm_hour)+":"+std::to_string(endtime_tm.tm_min)+":"+std::to_string(endtime_tm.tm_sec),v_message,verbosity);

  for (int i_channel = 0; i_channel < num_active_slots*num_channels_tank; i_channel++){
    std::vector<unsigned int> crateslotch_temp = map_ch_to_crateslotch[i_channel];
    unsigned int crate_temp = crateslotch_temp.at(0);
    unsigned int slot_temp = crateslotch_temp.at(1);
    unsigned int channel_temp = crateslotch_temp.at(2);    //channel numbering goes from 1-4
    double rate_temp = 0;
    double mean_temp = 0;
    double sigma_temp = 0.;
    double channelcount_temp = 0.;
    double samples_temp=0.;
    int num_zero_buffers=0;
    for (unsigned int i_t = 0; i_t < timestamp_file.size(); i_t++){
      //if (i_channel == 0) std::cout <<"ped_file.at i_t = "<<i_t<<": "<<ped_file.at(i_t).at(i_channel)<<", sigma_file = "<<sigma_file.at(i_t).at(i_channel)<<", rate_temp: "<<rate_file.at(i_t).at(i_channel)<<std::endl; 
      if (rate_file.at(i_t).at(i_channel) < 0.001 && sigma_file.at(i_t).at(i_channel) < 0.001 && ped_file.at(i_t).at(i_channel) < 0.001){
	num_zero_buffers++;
	continue;
}
      rate_temp += rate_file.at(i_t).at(i_channel);
      mean_temp += ped_file.at(i_t).at(i_channel);
      sigma_temp += sigma_file.at(i_t).at(i_channel);
      samples_temp += samples_file.at(i_t).at(i_channel);
    }
    int num_buffers = int(ped_file.size()) - num_zero_buffers;
    if (num_buffers>0) {
    mean_temp/=num_buffers;
    sigma_temp/=num_buffers;
    samples_temp/=num_buffers;
    }
    channelcount_temp = rate_temp;
    double t_acquisition = num_buffers*samples_temp*ADC_TO_NS;
    //if (t_frame>0.) rate_temp /= (t_frame/1000.);  //convert into units of 1/s
    if (t_acquisition > 0.) rate_temp /= (t_acquisition/1000000.);	//convert from us to seconds

    crate->push_back(crate_temp);
    slot->push_back(slot_temp);
    channel->push_back(channel_temp);
    ped->push_back(mean_temp);
    sigma->push_back(sigma_temp);
    rate->push_back(rate_temp);
    channelcount->push_back(channelcount_temp);
  }

  t->Fill();
  t->Write("",TObject::kOverwrite);           //prevent ROOT from making endless keys for the same tree when updating the tree
  f->Close();

  delete crate;
  delete slot;
  delete channel;
  delete ped;
  delete sigma;
  delete rate;
  delete channelcount;
  delete f;     //tree should get deleted automatically by closing file

  gROOT->cd();

}

void MonitorTankTime::ReadFromFile(ULong64_t timestamp_end, double time_frame){
  
  Log("MonitorTankTime: ReadFromFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------ReadFromFile ------------------------
  //-------------------------------------------------------

  ped_plot.clear();
  sigma_plot.clear();
  rate_plot.clear();
  channelcount_plot.clear();
  tstart_plot.clear();
  tend_plot.clear();
  labels_timeaxis.clear();

  //take the end time and calculate the start time with the given time_frame
  ULong64_t timestamp_start = timestamp_end - time_frame*MIN_to_HOUR*SEC_to_MIN*MSEC_to_SEC;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(timestamp_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_start/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_start/MSEC_to_SEC/1000.)%60,timestamp_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);

   Log("MonitorTankTime: ReadFromFile: Reading in data for time frame "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
    +" ... "+std::to_string(endtime_tm.tm_year+1900)+"/"+std::to_string(endtime_tm.tm_mon+1)+"/"+std::to_string(endtime_tm.tm_mday)+"-"+std::to_string(endtime_tm.tm_hour)+":"+std::to_string(endtime_tm.tm_min)+":"+std::to_string(endtime_tm.tm_sec),v_message,verbosity);

  std::stringstream ss_startdate, ss_enddate;
  ss_startdate << starttime_tm.tm_year+1900 << "-" << starttime_tm.tm_mon+1 <<"-"<< starttime_tm.tm_mday;
  ss_enddate << endtime_tm.tm_year+1900 << "-" << endtime_tm.tm_mon+1 <<"-"<< endtime_tm.tm_mday;

  int number_of_days;
  if (endtime_tm.tm_mon == starttime_tm.tm_mon) number_of_days = endtime_tm.tm_mday - starttime_tm.tm_mday;
  else {
    boost::gregorian::date endtime_greg(boost::gregorian::from_simple_string(ss_enddate.str()));
    boost::gregorian::date starttime_greg(boost::gregorian::from_simple_string(ss_startdate.str()));
    boost::gregorian::days days_dataframe = endtime_greg - starttime_greg;
    number_of_days = days_dataframe.days();
  }

  for (int i_day = 0; i_day <= number_of_days; i_day++){

    ULong64_t timestamp_i = timestamp_start+i_day*HOUR_to_DAY*MIN_to_HOUR*SEC_to_MIN*MSEC_to_SEC;
    std::string string_date_i = convertTimeStamp_to_Date(timestamp_i);
    std::stringstream root_filename_i;
    root_filename_i << path_monitoring << "PMT_" << string_date_i <<".root";
    bool tree_exists = true;

    if (does_file_exist(root_filename_i.str())) {
      TFile *f = new TFile(root_filename_i.str().c_str(),"READ");
      TTree *t;
      if (f->GetListOfKeys()->Contains("tankmonitor_tree")) t = (TTree*) f->Get("tankmonitor_tree");
      else { 
        Log("WARNING (MonitorTankTime): File "+root_filename_i.str()+" does not contain tankmonitor_tree. Omit file.",v_warning,verbosity);
        tree_exists = false;
      }

      if (tree_exists){

        Log("MonitorTankTime: Tree exists, start reading in data",v_message,verbosity);

        ULong64_t t_start, t_end;
        std::vector<unsigned int> *crate = new std::vector<unsigned int>;
        std::vector<unsigned int> *slot = new std::vector<unsigned int>;
        std::vector<unsigned int> *channel = new std::vector<unsigned int>;
        std::vector<double> *ped = new std::vector<double>;
        std::vector<double> *sigma = new std::vector<double>;
        std::vector<double> *rate = new std::vector<double>;
        std::vector<int> *channelcount = new std::vector<int>;

        int nevents;
        int nentries_tree;

        t->SetBranchAddress("t_start",&t_start);
        t->SetBranchAddress("t_end",&t_end);
        t->SetBranchAddress("crate",&crate);
        t->SetBranchAddress("slot",&slot);
        t->SetBranchAddress("channel",&channel);
        t->SetBranchAddress("ped",&ped);
        t->SetBranchAddress("sigma",&sigma);
        t->SetBranchAddress("rate",&rate);
        t->SetBranchAddress("channelcount",&channelcount);

        nentries_tree = t->GetEntries();
	      
	//Sort timestamps for the case that they are not in order

	std::vector<ULong64_t> vector_timestamps;
        std::map<ULong64_t,int> map_timestamp_entry;
	for (int i_entry = 0; i_entry < nentries_tree; i_entry++){
	  t->GetEntry(i_entry);
	  if (t_start >= timestamp_start && t_end <= timestamp_end){
	    vector_timestamps.push_back(t_start);
	    map_timestamp_entry.emplace(t_start,i_entry);	    
	  }
	}

	std::sort(vector_timestamps.begin(), vector_timestamps.end());
	std::vector<int> vector_sorted_entry;

	for (int i_entry = 0; i_entry < (int) vector_timestamps.size(); i_entry++){
	  vector_sorted_entry.push_back(map_timestamp_entry.at(vector_timestamps.at(i_entry)));
	}

        for (int i_entry = 0; i_entry < (int) vector_sorted_entry.size(); i_entry++){
		
	  int next_entry = vector_sorted_entry.at(i_entry);

          t->GetEntry(next_entry);
          if (t_start >= timestamp_start && t_end <= timestamp_end){
            ped_plot.push_back(*ped);
            sigma_plot.push_back(*sigma);
            rate_plot.push_back(*rate);
            channelcount_plot.push_back(*channelcount);
            tstart_plot.push_back(t_start);
            tend_plot.push_back(t_end);
            boost::posix_time::ptime boost_tend = *Epoch+boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
            struct tm label_timestamp = boost::posix_time::to_tm(boost_tend);
            
            TDatime datime_timestamp(1900+label_timestamp.tm_year,label_timestamp.tm_mon+1,label_timestamp.tm_mday,label_timestamp.tm_hour,label_timestamp.tm_min,label_timestamp.tm_sec);
            labels_timeaxis.push_back(datime_timestamp);
          }

        }

        delete crate;
        delete slot;
        delete channel;
        delete ped;
        delete sigma;
        delete rate;
        delete channelcount;

      }

      f->Close();
      delete f;
      gROOT->cd();

    } else {
      Log("MonitorTankTime: ReadFromFile: File "+root_filename_i.str()+" does not exist. Omit file.",v_warning,verbosity);
    }

  }

  //Set the readfromfile time variables to make sure data is not read twice for the same time window
  readfromfile_tend = timestamp_end;
  readfromfile_timeframe = time_frame;

}

void MonitorTankTime::DrawLastFilePlots(){

  Log("MonitorTankTime: DrawLastFilePlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawLastFilePlots -------------------
  //-------------------------------------------------------

  //draw time buffer plots
  DrawBufferPlots();

  //draw ADC frequency plots
  DrawADCFreqPlots();

  //draw FIFO error plots
  DrawFIFOPlots();

  //draw VME histogram plots
  DrawVMEHistogram();

  //Draw ped plots plots
  DrawPedPlotElectronics(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");
  //DrawPedPlotPhysical(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw ped Sigma plots
  DrawSigmaPlotElectronics(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");
  //DrawSigmaPlotPhysical(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw rate plots
  DrawRatePlotElectronics(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");
  //DrawRatePlotPhysical(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw hitmap plots
  DrawHitMap(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

}

void MonitorTankTime::UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes){

  Log("MonitorTankTime: UpdateMonitorPlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------UpdateMonitorPlots ------------------
  //-------------------------------------------------------

  //Draw the monitoring plots according to the specifications in the configfiles

  
  for (unsigned int i_time = 0; i_time < timeFrames.size(); i_time++){
    //std::cout <<"endTimes: "<<endTimes.at(i_time)<<", timeFrames: "<<timeFrames.at(i_time)<<", label: "<<fileLabels.at(i_time)<<std::endl;
  }

  for (unsigned int i_time = 0; i_time < timeFrames.size(); i_time++){

    ULong64_t zero = 0;
    if (endTimes.at(i_time) == zero) endTimes.at(i_time) = t_file_end;        //set 0 for t_file_end since we did not know what that was at the beginning of initialise


    for (unsigned int i_plot = 0; i_plot < plotTypes.at(i_time).size(); i_plot++){
      //std::cout <<"i_plot: "<<i_plot<<std::endl;
      if (plotTypes.at(i_time).at(i_plot) == "RateElectronics") DrawRatePlotElectronics(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "RatePhysical") DrawRatePlotPhysical(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "PedElectronics") DrawPedPlotElectronics(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "PedPhysical") DrawPedPlotPhysical(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "SigmaElectronics") DrawSigmaPlotElectronics(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "SigmaPhysical") DrawSigmaPlotPhysical(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "TimeEvolution") DrawTimeEvolution(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "TimeDifference") DrawTimeDifference(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "HitMap") DrawHitMap(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "FileHistory") DrawFileHistory(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time),1);
      else {
        if (verbosity > 0) std::cout <<"ERROR (MonitorTankTime): UpdateMonitorPlots: Specified plot type -"<<plotTypes.at(i_time).at(i_plot)<<"- does not exist! Omit entry."<<std::endl;
      }
    }
  }

}


void MonitorTankTime::DrawRatePlotElectronics(ULong64_t timestamp_end, double time_frame, std::string file_ending){


  Log("MonitorTankTime: DrawRatePlotElectronics",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;
 
  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //-------------------------------------------------------
  //-------------DrawRatePlotElectronics ------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  std::vector<double> overall_rates;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_rates.assign(num_active_slots*num_channels_tank,0);

  for (unsigned int i_file = 0; i_file < rate_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

  for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) continue;
      overall_rates.at(i_ch) += (rate_plot.at(i_file).at(i_ch)*current_timeframe);
    }
  }

  if (overall_timeframe > 0.){
    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      overall_rates.at(i_ch)/=overall_timeframe;
    }
  }

  long max_rate = 0;
  long min_rate = 999999;
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_rate->SetBinContent(x,y,overall_rates.at(i_active));
      if (overall_rates.at(i_active) > max_rate) max_rate = overall_rates.at(i_active);
      else if (overall_rates.at(i_active) < min_rate) min_rate = overall_rates.at(i_active);
    }
  }

  ss_title_rate.str("");
  ss_title_rate << "Rate VME (last "<<ss_timeframe.str()<<"h) "<<end_time.str()<<std::endl;
  h2D_rate->SetTitle(ss_title_rate.str().c_str());

  TPad *p_rate = (TPad*) canvas_rate->cd();
  h2D_rate->SetStats(0);
  h2D_rate->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_rate->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_rate->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_rate->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_rate->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_rate->LabelsOption("v");
  h2D_rate->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  line1->Draw("same");
  line2->Draw("same");
  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  label_rate->Draw();
  p_rate ->Update();

  if (h2D_rate->GetMaximum()>0.){
  if (abs(max_rate-min_rate)==0) h2D_sigma->GetZaxis()->SetRangeUser(min_rate-1,max_rate+1);
  else h2D_rate->GetZaxis()->SetRangeUser(1e-6,max_rate+0.5);
  TPaletteAxis *palette = 
  (TPaletteAxis*)h2D_rate->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.9);
  palette->SetX2NDC(0.92);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  }
  p_rate->Update();

  std::stringstream ss_rate;
  ss_rate<<outpath<<"PMT_Rate_Electronics_"<<file_ending<<".jpg";
  canvas_rate->SaveAs(ss_rate.str().c_str());
  Log("MonitorTankTime: Output path Rate plot (Electronics space): "+ss_rate.str(),v_message,verbosity);

  h2D_rate->Reset();
  canvas_rate->Clear();

}

void MonitorTankTime::DrawRatePlotPhysical(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawRatePlotPhysical",v_message,verbosity);

  //-------------------------------------------------------
  //-------------DrawRatePlotPhysical ---------------------
  //-------------------------------------------------------

  //TODO


}

void MonitorTankTime::DrawPedPlotElectronics(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawPedPlotElectronics",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //-------------------------------------------------------
  //-------------DrawPedPlotElectronics -------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  std::vector<double> overall_peds;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_peds.assign(num_active_slots*num_channels_tank,0);

  for (unsigned int i_file = 0; i_file < ped_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) continue;
      overall_peds.at(i_ch) += (ped_plot.at(i_file).at(i_ch)*current_timeframe);
    }
  }

  if (overall_timeframe > 0.){
    for (int i_ch = 0.; i_ch < num_active_slots*num_channels_tank; i_ch++){
      overall_peds.at(i_ch)/=overall_timeframe;
    }
  }

  double max_ped = 0;
  double min_ped = 999999;
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_ped->SetBinContent(x,y,overall_peds.at(i_active));
      if (overall_peds.at(i_active) > max_ped) max_ped = overall_peds.at(i_active);
      else if (overall_peds.at(i_active) < min_ped) min_ped = overall_peds.at(i_active);
    }
  }

  ss_title_ped.str("");
  ss_title_ped << "Pedestal VME (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  h2D_ped->SetTitle(ss_title_ped.str().c_str());

  TPad *p = (TPad*) canvas_ped->cd();
  h2D_ped->SetStats(0);
  h2D_ped->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_ped->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_ped->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_ped->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_ped->LabelsOption("v");
  h2D_ped->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  p->Update();
  if (h2D_ped->GetMaximum()>0.){
  if (abs(max_ped-min_ped)==0) h2D_ped->GetZaxis()->SetRangeUser(min_ped-1,max_ped+1);
  else h2D_ped->GetZaxis()->SetRangeUser(200.,400.);
//  else h2D_ped->GetZaxis()->SetRangeUser(1e-6,max_ped+0.5);
  TPaletteAxis *palette = 
  (TPaletteAxis*)h2D_ped->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.9);
  palette->SetX2NDC(0.92);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  }
  p->Update();

  line1->Draw("same");
  line2->Draw("same");
  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  label_ped->Draw();

  std::stringstream ss_ped;
  ss_ped<<outpath<<"PMT_Ped_Electronics_"<<file_ending<<".jpg";
  Log("MonitorTankTime: Output path Pedestal mean plot (Electronics space): "+ss_ped.str(),v_message,verbosity);
  canvas_ped->SaveAs(ss_ped.str().c_str());

  h2D_ped->Reset();
  canvas_ped->Clear();

}

void MonitorTankTime::DrawPedPlotPhysical(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawPedPlotPhysical",v_message,verbosity);

  //-------------------------------------------------------
  //-------------DrawPedPlotPhysical ----------------------
  //-------------------------------------------------------

  //TODO

}

void MonitorTankTime::DrawSigmaPlotElectronics(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawSigmaPlotElectronics",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //-------------------------------------------------------
  //-------------DrawSigmaPlotElectronics -----------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);
    
  std::vector<double> overall_sigmas;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_sigmas.assign(num_active_slots*num_channels_tank,0);

  for (unsigned int i_file = 0; i_file < sigma_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) continue;
      overall_sigmas.at(i_ch) += (sigma_plot.at(i_file).at(i_ch)*current_timeframe);
    }
  }

  if (overall_timeframe > 0.){
    for (int i_ch = 0.; i_ch < num_active_slots*num_channels_tank; i_ch++){
      overall_sigmas.at(i_ch)/=overall_timeframe;
    }
  }

  double max_sigma = 0;
  double min_sigma = 999999;
  
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_sigma->SetBinContent(x,y,overall_sigmas.at(i_active));
      if (overall_sigmas.at(i_active) > max_sigma) max_sigma = overall_sigmas.at(i_active);
      else if (overall_sigmas.at(i_active) < min_sigma) min_sigma = overall_sigmas.at(i_active);
    }
  }

  ss_title_sigma.str("");
  ss_title_sigma << "Sigma VME (last "<< ss_timeframe.str() << "h) "<<end_time.str();
  h2D_sigma->SetTitle(ss_title_sigma.str().c_str());

  TPad *p_sigma = (TPad*) canvas_sigma->cd();
  h2D_sigma->SetStats(0);
  h2D_sigma->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_sigma->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_sigma->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_sigma->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_sigma->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_sigma->LabelsOption("v");
  h2D_sigma->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  line1->Draw("same");
  line2->Draw("same");
  label_sigma->Draw();
  p_sigma ->Update();

  if (h2D_sigma->GetMaximum()>0.){
    if (abs(max_sigma-min_sigma)==0) h2D_sigma->GetZaxis()->SetRangeUser(min_sigma-1,max_sigma+1);
    else h2D_sigma->GetZaxis()->SetRangeUser(1e-6,max_sigma+0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_sigma->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }
  p_sigma->Update();
  std::stringstream ss_sigma;
  ss_sigma<<outpath<<"PMT_Sigma_Electronics_"<<file_ending<<".jpg";
  canvas_sigma->SaveAs(ss_sigma.str().c_str());
  Log("MonitorTankTime: Output path Pedestal Sigma plot (Electronics Space): "+ss_sigma.str(),v_message,verbosity);

  h2D_sigma->Reset();
  canvas_sigma->Clear();

}

void MonitorTankTime::DrawSigmaPlotPhysical(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawSigmaPlotPhysical",v_message,verbosity);


  //-------------------------------------------------------
  //-------------DrawSigmaPlotPhysical ----------------------
  //-------------------------------------------------------

  //TODO

}

void MonitorTankTime::DrawHitMap(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawHitMap",v_message,verbosity);

  //-------------------------------------------------------
  //---------------------DrawHitMap -----------------------
  //-------------------------------------------------------

   boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  std::stringstream ss_timeframe, ss_hitmap_title;
  ss_timeframe << round(time_frame*100.)/100.;

  //variables for scaling hitmaps
  long max_hitmap=0;
  std::vector<long> max_hitmap_slot;
  max_hitmap_slot.assign(num_active_slots,0); 

  for (int i_channel = 0; i_channel < num_active_slots*num_channels_tank; i_channel++){

    long sum_channel = 0;
    
    for (unsigned int i_file = 0; i_file < rate_plot.size(); i_file++){
      sum_channel+=channelcount_plot.at(i_file).at(i_channel);
    }

    int active_slot = i_channel / num_channels_tank;
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[active_slot];
    unsigned int i_crate = crateslot.at(0)-1;
    unsigned int slot = crateslot.at(1);
    unsigned int channel = i_channel % num_channels_tank;

    for (int i_cr = 0; i_cr < 3; i_cr++){
      if (int(i_crate) != i_cr) hist_hitmap[i_cr]->SetBinContent(i_channel+1,0.001);
    }
    if (sum_channel!=0){
      if (max_hitmap < sum_channel) max_hitmap = sum_channel;
      if (max_hitmap_slot.at(active_slot) < sum_channel) max_hitmap_slot.at(active_slot) = sum_channel;
      hist_hitmap[i_crate]->SetBinContent(i_channel+1,sum_channel);
      hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1,sum_channel);
    } else {
      hist_hitmap[i_crate]->SetBinContent(i_channel+1,0.001);
      hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1,0.001);
    }
  }

  //Plot the overall hitmap histogram

  canvas_hitmap_tank->cd();
  hist_hitmap[0]->GetYaxis()->SetRangeUser(0.8,max_hitmap+10);
  hist_hitmap[1]->GetYaxis()->SetRangeUser(0.8,max_hitmap+10);
  hist_hitmap[2]->GetYaxis()->SetRangeUser(0.8,max_hitmap+10);
  ss_hitmap_title << "Hitmap VME "<<end_time.str()<<" (last "<<ss_timeframe.str()<<"h)";
  hist_hitmap[0]->SetTitle(ss_hitmap_title.str().c_str());
  hist_hitmap[0]->Draw();
  hist_hitmap[1]->Draw("same");
  hist_hitmap[2]->Draw("same");
  separate_crates->SetY2(max_hitmap+10);
  separate_crates->SetLineStyle(2);
  separate_crates->SetLineWidth(2);
  separate_crates->Draw("same");
  separate_crates2->SetY2(max_hitmap+10);
  separate_crates2->SetLineStyle(2);
  separate_crates2->SetLineWidth(2);
  separate_crates2->Draw("same");
  label_cr1->Draw();
  label_cr2->Draw();
  label_cr3->Draw();
  canvas_hitmap_tank->Update();
  TGaxis *labels_grid = new TGaxis(0,canvas_hitmap_tank->GetUymin(),num_active_slots*num_channels_tank,canvas_hitmap_tank->GetUymin(),"f1",num_active_slots,"w");
  labels_grid->SetLabelSize(0);
  labels_grid->Draw("w");
  std::stringstream save_path_hitmap;
  save_path_hitmap << outpath <<"PMTHitmap_"<<file_ending<<"."<<img_extension;
  canvas_hitmap_tank->SaveAs(save_path_hitmap.str().c_str());

  delete labels_grid;

  //Plot the more detailed slot-wise hitmap distributions

  //std::cout <<"num_active_slots: "<<num_active_slots<<std::endl;
  for (int i_slot=0; i_slot < num_active_slots; i_slot++){
    //std::cout <<"hitmap single slot "<<i_slot<<std::endl;
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);

    canvas_hitmap_tank_slot->cd();
    std::stringstream ss_slot;
    ss_slot << i_slot;
    std::stringstream save_path_singlehitmap, ss_hitmap_title_slot;
    save_path_singlehitmap << outpath <<"PMTHitmap_Cr"<<crate<<"_Sl"<<slot<<"_"<<file_ending<<"."<<img_extension;
    ss_hitmap_title_slot << "Hitmap VME "<<end_time.str()<<" Cr "<<crate <<" Sl "<<slot<<" (last "<<ss_timeframe.str()<<"h)";
    //std::cout <<"savepath single slot: "<<save_path_singlehitmap.str()<<std::endl;
    hist_hitmap_slot.at(i_slot)->GetYaxis()->SetRangeUser(0.8,max_hitmap_slot.at(i_slot)+10);
    hist_hitmap_slot.at(i_slot)->SetTitle(ss_hitmap_title_slot.str().c_str());
    hist_hitmap_slot.at(i_slot)->Draw();

    for (unsigned int i_box = 0; i_box < vector_box_inactive_hitmap[i_slot].size(); i_box++){
      vector_box_inactive_hitmap[i_slot].at(i_box)->SetY2(max_hitmap_slot.at(i_slot)+10);
      vector_box_inactive_hitmap[i_slot].at(i_box)->Draw("same");
    }
    canvas_hitmap_tank_slot->Update();
    canvas_hitmap_tank_slot->SaveAs(save_path_singlehitmap.str().c_str());
    canvas_hitmap_tank_slot->Clear();
  }  


}

void MonitorTankTime::DrawTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawTimeEvolution",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  //-------------------------------------------------------
  //-------------DrawTimeEvolution ------------------------
  //-------------------------------------------------------

   if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  //looping over all files that are in the time interval, each file will be one data point

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  for (int i_channel = 0; i_channel < num_active_slots*num_channels_tank; i_channel++){
    gr_ped.at(i_channel)->Set(0);
    gr_sigma.at(i_channel)->Set(0);
    gr_rate.at(i_channel)->Set(0);
  }

  for (unsigned int i_file=0; i_file<ped_plot.size(); i_file++){

    //Updating channel graphs

    Log("MonitorTankTime: Stored data (file #"+std::to_string(i_file+1)+"): ",v_message,verbosity);
    for (int i_channel = 0; i_channel < num_active_slots*num_channels_tank; i_channel++){
      gr_ped.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),ped_plot.at(i_file).at(i_channel));
      gr_sigma.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),sigma_plot.at(i_file).at(i_channel));
      gr_rate.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),rate_plot.at(i_file).at(i_channel));
    }

  }

  // Drawing time evolution plots

  double max_canvas = 0;
  double min_canvas = 9999999.;
  double max_canvas_sigma = 0;
  double min_canvas_sigma = 99999999.;
  double max_canvas_rate = 0;
  double min_canvas_rate = 999999999.;

  int CH_per_CANVAS = 4;   //channels per canvas

  for (int i_channel = 0; i_channel<num_active_slots*num_channels_tank; i_channel++){
     

    std::stringstream ss_ch_ped, ss_ch_sigma, ss_ch_rate, ss_leg_time;
  
     if (i_channel%CH_per_CANVAS == 0 || i_channel == num_active_slots*num_channels_tank-1) {
      if (i_channel != 0){
        
        std::vector<unsigned int> crateslotchannel = map_ch_to_crateslotch[i_channel-1];
        unsigned int crate = crateslotchannel.at(0);
        unsigned int slot = crateslotchannel.at(1);
        unsigned int channel = crateslotchannel.at(2);

        ss_ch_ped.str("");
        ss_ch_sigma.str("");
        ss_ch_rate.str("");

        ss_ch_ped<<"Crate "<<crate<<" Slot "<<slot<<" ("<<ss_timeframe.str()<<"h) "<<end_time.str();
        ss_ch_sigma<<"Crate "<<crate<<" Slot "<<slot<<" ("<<ss_timeframe.str()<<"h) "<<end_time.str();
        ss_ch_rate<<"Crate "<<crate<<" Slot "<<slot<<" ("<<ss_timeframe.str()<<"h) "<<end_time.str();

        if (i_channel == num_active_slots*num_channels_tank - 1 && std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_channel)==vec_disabled_global.end()){
          ss_leg_time.str("");
          ss_leg_time<<"ch "<<channel + 1;
          multi_ch_ped->Add(gr_ped.at(i_channel));
          leg_ped->AddEntry(gr_ped.at(i_channel),ss_leg_time.str().c_str(),"l");
          multi_ch_sigma->Add(gr_sigma.at(i_channel));
          leg_sigma->AddEntry(gr_sigma.at(i_channel),ss_leg_time.str().c_str(),"l");
          multi_ch_rate->Add(gr_rate.at(i_channel));
          leg_rate->AddEntry(gr_rate.at(i_channel),ss_leg_time.str().c_str(),"l");
        }

        canvas_ch_ped->cd();
        multi_ch_ped->Draw("apl");
        multi_ch_ped->SetTitle(ss_ch_ped.str().c_str());
        multi_ch_ped->GetYaxis()->SetTitle("Pedestal");
        multi_ch_ped->GetXaxis()->SetTimeDisplay(1);
        multi_ch_ped->GetXaxis()->SetLabelSize(0.03);
        multi_ch_ped->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_ped->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_ped->GetXaxis()->SetTimeOffset(0.);
        leg_ped->Draw();
        canvas_ch_sigma->cd();
        multi_ch_sigma->Draw("apl");
        multi_ch_sigma->SetTitle(ss_ch_sigma.str().c_str());
        multi_ch_sigma->GetYaxis()->SetTitle("Sigma");
        multi_ch_sigma->GetYaxis()->SetTitleSize(0.035);
        multi_ch_sigma->GetYaxis()->SetTitleOffset(1.3);
        multi_ch_sigma->GetXaxis()->SetTimeDisplay(1);
        multi_ch_sigma->GetXaxis()->SetLabelSize(0.03);
        multi_ch_sigma->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_sigma->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_sigma->GetXaxis()->SetTimeOffset(0.);
        leg_sigma->Draw();
        canvas_ch_rate_tank->cd();
        multi_ch_rate->Draw("apl");
        multi_ch_rate->SetTitle(ss_ch_rate.str().c_str());
        multi_ch_rate->GetYaxis()->SetTitle("Rate [kHz]");
        multi_ch_rate->GetYaxis()->SetTitleSize(0.035);
        multi_ch_rate->GetYaxis()->SetTitleOffset(1.3);
        multi_ch_rate->GetXaxis()->SetTimeDisplay(1);
        multi_ch_rate->GetXaxis()->SetLabelSize(0.03);
        multi_ch_rate->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_rate->GetXaxis()->SetTimeOffset(0.);
        leg_rate->Draw();

        ss_ch_ped.str("");
        ss_ch_sigma.str("");
        ss_ch_rate.str("");
        ss_ch_ped<<outpath<<"PMTTimeEvolutionPed_Cr"<<crate<<"_Sl"<<slot<<"_"<<file_ending<<"."<<img_extension;
        ss_ch_sigma<<outpath<<"PMTTimeEvolutionSigma_Cr"<<crate<<"_Sl"<<slot<<"_"<<file_ending<<"."<<img_extension;
        ss_ch_rate<<outpath<<"PMTTimeEvolutionRate_Cr"<<crate<<"_Sl"<<slot<<"_"<<file_ending<<"."<<img_extension;

        canvas_ch_ped->SaveAs(ss_ch_ped.str().c_str());
        canvas_ch_sigma->SaveAs(ss_ch_sigma.str().c_str());
        canvas_ch_rate_tank->SaveAs(ss_ch_rate.str().c_str()); 

        for (int i_gr=0; i_gr < CH_per_CANVAS; i_gr++){
          int i_balance = (i_channel == num_active_slots*num_channels_tank-1)? 1 : 0;
          if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_channel-CH_per_CANVAS+i_gr+i_balance)!=vec_disabled_global.end()) continue;
          multi_ch_ped->RecursiveRemove(gr_ped.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
          multi_ch_sigma->RecursiveRemove(gr_sigma.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
          multi_ch_rate->RecursiveRemove(gr_rate.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
        }
      }

      leg_ped->Clear();
      leg_sigma->Clear();
      leg_rate->Clear();
      
      canvas_ch_ped->Clear();
      canvas_ch_sigma->Clear();
      canvas_ch_rate_tank->Clear();

     } 

    if (i_channel != num_active_slots*num_channels_tank-1){

      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_channel)!=vec_disabled_global.end()) continue;
      ss_leg_time.str("");
      ss_leg_time<<"ch "<<i_channel%(num_channels_tank) + 1;
      multi_ch_ped->Add(gr_ped.at(i_channel));
      if (gr_ped.at(i_channel)->GetMaximum()>max_canvas) max_canvas = gr_ped.at(i_channel)->GetMaximum();
      if (gr_ped.at(i_channel)->GetMinimum()<min_canvas) min_canvas = gr_ped.at(i_channel)->GetMinimum();
      leg_ped->AddEntry(gr_ped.at(i_channel),ss_leg_time.str().c_str(),"l");
      multi_ch_sigma->Add(gr_sigma.at(i_channel));
      leg_sigma->AddEntry(gr_sigma.at(i_channel),ss_leg_time.str().c_str(),"l");
      if (gr_sigma.at(i_channel)->GetMaximum()>max_canvas_sigma) max_canvas_sigma = gr_sigma.at(i_channel)->GetMaximum();
      if (gr_sigma.at(i_channel)->GetMinimum()<min_canvas_sigma) min_canvas_sigma = gr_sigma.at(i_channel)->GetMinimum();
      multi_ch_rate->Add(gr_rate.at(i_channel));
      leg_rate->AddEntry(gr_rate.at(i_channel),ss_leg_time.str().c_str(),"l");
      if (gr_rate.at(i_channel)->GetMaximum()>max_canvas_rate) max_canvas_rate = gr_rate.at(i_channel)->GetMaximum();
      if (gr_rate.at(i_channel)->GetMinimum()<min_canvas_rate) min_canvas_rate = gr_rate.at(i_channel)->GetMinimum();

    }

    //single channel time evolution plots

    if (draw_single){
     
      std::vector<unsigned int> crateslotchannel = map_ch_to_crateslotch[i_channel];
      unsigned int crate = crateslotchannel.at(0);
      unsigned int slot = crateslotchannel.at(1);
      unsigned int channel = crateslotchannel.at(2);

      canvas_ch_single_tank->cd();
      canvas_ch_single_tank->Clear();
      gr_ped.at(i_channel)->GetYaxis()->SetTitle("Pedestal");
      gr_ped.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_ped.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_ped.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_ped.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_ped.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_ped.at(i_channel)->Draw("apl");
      std::stringstream ss_ch_single;
      ss_ch_single<<outpath<<"PMTTimeEvolutionPed_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single_tank->SaveAs(ss_ch_single.str().c_str());

      canvas_ch_single_tank->Clear();
      gr_sigma.at(i_channel)->GetYaxis()->SetTitle("Sigma");
      gr_sigma.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_sigma.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_sigma.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_sigma.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_sigma.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_sigma.at(i_channel)->Draw("apl");
      ss_ch_single.str("");
      ss_ch_single<<outpath<<"PMTTimeEvolutionSigma_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single_tank->SaveAs(ss_ch_single.str().c_str());

      canvas_ch_single_tank->Clear();
      gr_rate.at(i_channel)->GetYaxis()->SetTitle("Rate [kHz]");
      gr_rate.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_rate.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_rate.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_rate.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_rate.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_rate.at(i_channel)->Draw("apl");
      ss_ch_single.str("");
      ss_ch_single<<outpath<<"PMTTimeEvolutionRate_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single_tank->SaveAs(ss_ch_single.str().c_str());

    }

  }


}

void MonitorTankTime::DrawTimeDifference(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTankTime: DrawTimeDifference",v_message,verbosity);
  
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //-------------------------------------------------------
  //----------------DrawTimeDifference --------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);
    
  double max_peddiff = 0.;
  double min_peddiff = 999999.;
  double max_sigmadiff = 0.;
  double min_sigmadiff = 9999999.;
  double max_ratediff = 0.;
  double min_ratediff = 9999999.; 

  //pedestal time difference plot

  std::vector<double> overall_peddiffs;
  overall_peddiffs.assign(num_active_slots*num_channels_tank,0.);

  if (ped_plot.size() >= 1){
    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) overall_peddiffs.at(i_ch) = -9999;
      else overall_peddiffs.at(i_ch) = ped_plot.at(ped_plot.size()-1).at(i_ch) - ped_plot.at(0).at(i_ch);
    }
    for (unsigned int i_inactive=0; i_inactive < inactive_xy.size(); i_inactive++){
      h2D_peddiff->SetBinContent(inactive_xy.at(i_inactive).at(0),inactive_xy.at(i_inactive).at(1),-9999);
    } 
  }

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_peddiff->SetBinContent(x,y,overall_peddiffs.at(i_active));
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_active)!=vec_disabled_global.end()) continue;
      if (overall_peddiffs.at(i_active) > max_peddiff) max_peddiff = overall_peddiffs.at(i_active);
      else if (overall_peddiffs.at(i_active) < min_peddiff) min_peddiff = overall_peddiffs.at(i_active);
    }
  }

  ss_title_peddiff.str("");
  ss_title_peddiff << "Pedestal Difference (last " << ss_timeframe.str() <<"h) "<<end_time.str();
  h2D_peddiff->SetTitle(ss_title_peddiff.str().c_str());

  TPad *p_peddiff = (TPad*) canvas_peddiff->cd();
  h2D_peddiff->SetStats(0);
  h2D_peddiff->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_peddiff->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_peddiff->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_peddiff->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_peddiff->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_peddiff->LabelsOption("v");
  h2D_peddiff->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  line1->Draw("same");
  line2->Draw("same");
  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  label_peddiff->Draw();
  p_peddiff->Update();

  if (h2D_peddiff->GetMaximum()>0.){
    double global_max = (fabs(max_peddiff)>fabs(min_peddiff))? fabs(max_peddiff) : fabs(min_peddiff);
    if (abs(max_peddiff-min_peddiff)==0) h2D_peddiff->GetZaxis()->SetRangeUser(min_peddiff-1,max_peddiff+1);
    else h2D_peddiff->GetZaxis()->SetRangeUser(-global_max-0.5,global_max+0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_peddiff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  } else if (fabs(h2D_peddiff->GetMaximum())<0.001){
    h2D_peddiff->GetZaxis()->SetRangeUser(-0.5,0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_peddiff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }
  p_peddiff->Update();

  std::stringstream ss_peddiff;
  ss_peddiff<<outpath<<"PMT_PedDifference_"<<file_ending<<".jpg";
  canvas_peddiff->SaveAs(ss_peddiff.str().c_str());
  Log("MonitorTankTime: Output path Pedestal time difference plot: "+ss_peddiff.str(),v_message,verbosity);

  h2D_peddiff->Reset();
  canvas_peddiff->Clear();

  //sigma time difference plot

  std::vector<double> overall_sigmadiffs;
  overall_sigmadiffs.assign(num_active_slots*num_channels_tank,0.);

  if (sigma_plot.size() >= 1){
    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) overall_sigmadiffs.at(i_ch) = -9999;
      else overall_sigmadiffs.at(i_ch) = sigma_plot.at(sigma_plot.size()-1).at(i_ch) - sigma_plot.at(0).at(i_ch);
    }
    for (unsigned int i_inactive=0; i_inactive < inactive_xy.size(); i_inactive++){
      h2D_sigmadiff->SetBinContent(inactive_xy.at(i_inactive).at(0),inactive_xy.at(i_inactive).at(1),-9999);
    }
  }

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_sigmadiff->SetBinContent(x,y,overall_sigmadiffs.at(i_active));
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_active)!=vec_disabled_global.end()) continue;
      if (overall_sigmadiffs.at(i_active) > max_sigmadiff) max_sigmadiff = overall_sigmadiffs.at(i_active);
      else if (overall_sigmadiffs.at(i_active) < min_sigmadiff) min_sigmadiff = overall_sigmadiffs.at(i_active);
    }
  }

  ss_title_sigmadiff.str("");
  ss_title_sigmadiff << "Sigma Difference (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  h2D_sigmadiff->SetTitle(ss_title_sigmadiff.str().c_str());

  TPad *p_sigmadiff = (TPad*) canvas_sigmadiff->cd();
  h2D_sigmadiff->SetStats(0);
  h2D_sigmadiff->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_sigmadiff->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_sigmadiff->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_sigmadiff->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_sigmadiff->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_sigmadiff->LabelsOption("v");
  h2D_sigmadiff->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  line1->Draw("same");
  line2->Draw("same");
  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  label_sigmadiff->Draw();

  p_sigmadiff->Update();

  if (h2D_sigmadiff->GetMaximum()>=0.){
    double global_max = (fabs(max_sigmadiff)>fabs(min_sigmadiff))? fabs(max_sigmadiff) : fabs(min_sigmadiff);
    if (abs(max_sigmadiff-min_sigmadiff)==0) h2D_sigmadiff->GetZaxis()->SetRangeUser(min_sigmadiff-1,max_sigmadiff+1);
    else h2D_sigmadiff->GetZaxis()->SetRangeUser(-global_max-0.5,global_max+0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_sigmadiff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
    } else if (fabs(h2D_sigmadiff->GetMaximum())<0.001){
    h2D_sigmadiff->GetZaxis()->SetRangeUser(-0.5,0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_sigmadiff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }

  std::stringstream ss_sigmadiff;
  ss_sigmadiff<<outpath<<"PMT_SigmaDifference_"<<file_ending<<".jpg";
  canvas_sigmadiff->SaveAs(ss_sigmadiff.str().c_str());
  Log("MonitorTankTime: Output path Pedestal Sigma Time Difference plot: "+ss_sigmadiff.str(),v_message,verbosity);

  h2D_sigmadiff->Reset();
  canvas_sigmadiff->Clear();

  //rate time difference plot

  std::vector<double> overall_ratediffs;
  overall_ratediffs.assign(num_active_slots*num_channels_tank,0.);

  if (rate_plot.size() >= 1){
    for (int i_ch = 0; i_ch < num_active_slots*num_channels_tank; i_ch++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_ch)!=vec_disabled_global.end()) overall_ratediffs.at(i_ch) = -9999; 
      else overall_ratediffs.at(i_ch) = rate_plot.at(rate_plot.size()-1).at(i_ch) - rate_plot.at(0).at(i_ch);
    }
    for (unsigned int i_inactive=0; i_inactive < inactive_xy.size(); i_inactive++){
      h2D_ratediff->SetBinContent(inactive_xy.at(i_inactive).at(0),inactive_xy.at(i_inactive).at(1),-9999);
    }
  }

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    unsigned int i_crate = crateslot.at(0);
    unsigned int slot = crateslot.at(1);
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      int i_active = i_slot*num_channels_tank+i_channel;
      int x = slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank - i_channel;      //top most crate is displayed at the top, then crate #2 in the middle and crate #3 at the bottom
      h2D_ratediff->SetBinContent(x,y,overall_ratediffs.at(i_active));
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_active)!=vec_disabled_global.end()) continue; 
      if (overall_ratediffs.at(i_active) > max_ratediff) max_ratediff = overall_ratediffs.at(i_active);
      else if (overall_ratediffs.at(i_active) < min_ratediff) min_ratediff = overall_ratediffs.at(i_active);
    }
  }

  ss_title_ratediff.str("");
  ss_title_ratediff << "Rate Difference (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  h2D_ratediff->SetTitle(ss_title_ratediff.str().c_str());

  TPad *p_ratediff = (TPad*) canvas_ratediff->cd();
  h2D_ratediff->SetStats(0);
  h2D_ratediff->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_ratediff->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_ratediff->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_ratediff->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_ratediff->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_ratediff->LabelsOption("v");
  h2D_ratediff->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  line1->Draw("same");
  line2->Draw("same");
  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  label_ratediff->Draw();

  p_ratediff->Update();

  if (h2D_ratediff->GetMaximum()>0.){
    double global_max = (fabs(max_ratediff)>fabs(min_ratediff))? fabs(max_ratediff) : fabs(min_ratediff);
    if (abs(max_ratediff-min_ratediff)==0) h2D_ratediff->GetZaxis()->SetRangeUser(min_ratediff-1,max_ratediff+1);
    else h2D_ratediff->GetZaxis()->SetRangeUser(-global_max-0.5,global_max+0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_ratediff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  } else if (fabs(h2D_ratediff->GetMaximum())<0.001){
    h2D_ratediff->GetZaxis()->SetRangeUser(-0.5,0.5);
    TPaletteAxis *palette = 
    (TPaletteAxis*)h2D_ratediff->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }

  std::stringstream ss_ratediff;
  ss_ratediff<<outpath<<"PMT_RateDifference_"<<file_ending<<".jpg";
  canvas_ratediff->SaveAs(ss_ratediff.str().c_str());
  Log("MonitorTankTime: Output path Rate Time Difference plot: "+ss_ratediff.str(),v_message,verbosity);

  h2D_ratediff->Reset();
  canvas_ratediff->Clear();

}

void MonitorTankTime::DrawBufferPlots(){

  Log("MonitorTankTime: DrawBufferPlots",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  //-------------------------------------------------------
  //------------------DrawBufferPlots ---------------------
  //-------------------------------------------------------

  std::stringstream ss_canvas_temp;

  for (int i_slot = 0; i_slot<num_active_slots; i_slot++){

    int max_freq = 0;
    double max_temp = -999999.;
    double min_temp = 999999.;
    leg_temp->Clear();

    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    int crate_num = crateslot.at(0);
    int slot_num = crateslot.at(1);

    canvas_Channels_temp.at(i_slot)->cd();
    canvas_Channels_temp.at(i_slot)->Clear();
    ss_title_hist_temp.str("");
    ss_title_hist_temp << title_time.str() << " Temp Crate " << crate_num << " Slot " << slot_num <<" (last File) "<<end_time.str();
    hChannels_temp.at(i_slot*num_channels_tank)->SetTitle(ss_title_hist_temp.str().c_str());

    bool first_ch = true;
    int i_first_ch=0;
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_slot*num_channels_tank+i_channel)!=vec_disabled_global.end()) continue;
      std::stringstream ss_channel;
      ss_channel << "Cr "<<crate_num<<"/Sl"<<slot_num<<"/Ch"<<i_channel+1;
      if (first_ch) {
        hChannels_temp.at(i_slot*num_channels_tank+i_channel)->Draw();
        first_ch=false;
        i_first_ch=i_channel;
      }
      else hChannels_temp.at(i_slot*num_channels_tank+i_channel)->Draw("same");
      if (hChannels_temp.at(i_slot*num_channels_tank+i_channel)->GetMaximum() > max_temp) max_temp = hChannels_temp.at(i_slot*num_channels_tank+i_channel)->GetMaximum();
      if (hChannels_temp.at(i_slot*num_channels_tank+i_channel)->GetMinimum() < min_temp) min_temp = hChannels_temp.at(i_slot*num_channels_tank+i_channel)->GetMinimum();
      leg_temp->AddEntry(hChannels_temp.at(i_slot*num_channels_tank+i_channel),ss_channel.str().c_str(),"l");
    }
    hChannels_temp.at(i_slot*num_channels_tank+i_first_ch)->GetYaxis()->SetRangeUser(min_temp,max_temp);
    leg_temp->Draw();
    ss_canvas_temp.str("");
    ss_canvas_temp << outpath << "PMT_Temp_Cr"<<crate_num<<"_Sl"<<slot_num<<".jpg";
    canvas_Channels_temp.at(i_slot)->SaveAs(ss_canvas_temp.str().c_str());

  }
  
  ss_title_hist_temp.str("");
  ss_title_hist_temp << "Temp RWM (last File) "<<end_time.str();
  hChannels_temp_RWM->SetTitle(ss_title_hist_temp.str().c_str());
  canvas_Channels_temp.at(0)->cd();
  canvas_Channels_temp.at(0)->Clear();
  hChannels_temp_RWM->Draw();
  ss_canvas_temp.str("");
  ss_canvas_temp << outpath << "PMT_Temp_RWM.jpg";
  canvas_Channels_temp.at(0)->SaveAs(ss_canvas_temp.str().c_str());

  ss_title_hist_temp.str("");
  ss_title_hist_temp << "Temp BRF (last File) "<<end_time.str();
  hChannels_temp_BRF->SetTitle(ss_title_hist_temp.str().c_str());
  canvas_Channels_temp.at(0)->Clear();
  hChannels_temp_BRF->Draw();
  ss_canvas_temp.str("");
  ss_canvas_temp << outpath << "PMT_Temp_BRF.jpg";
  canvas_Channels_temp.at(0)->SaveAs(ss_canvas_temp.str().c_str());

  //Resetting & Clearing Channel histograms, canvasses
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_slot*num_channels_tank+i_channel)!=vec_disabled_global.end()) continue;
      hChannels_temp.at(i_slot*num_channels_tank+i_channel)->Reset();
      hChannels_temp.at(i_slot*num_channels_tank+i_channel)->SetMaximum(-1111); //needed to reset max/minimum computation
      hChannels_temp.at(i_slot*num_channels_tank+i_channel)->SetMinimum(-1111); //needed to reset max/minimum computation
    }
  }

  hChannels_temp_RWM->Reset();
  hChannels_temp_BRF->Reset();

}

void MonitorTankTime::DrawADCFreqPlots(){

  Log("MonitorTankTime: DrawADCFreqPlots",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  //-------------------------------------------------------
  //------------------DrawADCFreqPlots --------------------
  //-------------------------------------------------------

  std::stringstream ss_title_hist_freq;
  std::stringstream ss_canvas_freq;

  for (int i_slot = 0; i_slot<num_active_slots; i_slot++){

    int max_freq = 0;
    double max_temp = -999999.;
    double min_temp = 999999.;
    leg_freq->Clear();

    std::vector<unsigned int> crateslot = map_slot_to_crateslot[i_slot];
    int crate_num = crateslot.at(0);
    int slot_num = crateslot.at(1);

    canvas_Channels_freq.at(i_slot)->cd();
    canvas_Channels_freq.at(i_slot)->Clear();
    ss_title_hist_freq.str("");
    ss_title_hist_freq << title_time.str() << " Freq Crate " << crate_num << " Slot " << slot_num <<" (last File) "<<end_time.str();
    hChannels_freq.at(i_slot*num_channels_tank)->SetTitle(ss_title_hist_freq.str().c_str());

    bool first_ch=true;
    int i_first_ch=0;
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_slot*num_channels_tank+i_channel)!=vec_disabled_global.end()) continue;
      std::stringstream ss_channel;
      ss_channel << "Cr "<<crate_num<<"/Sl"<<slot_num<<"/Ch"<<i_channel+1;
      if (first_ch) {
        hChannels_freq.at(i_slot*num_channels_tank+i_channel)->Draw();
        first_ch=false;
	i_first_ch = i_channel;
      }
      else hChannels_freq.at(i_slot*num_channels_tank+i_channel)->Draw("same");
      if (hChannels_freq.at(i_slot*num_channels_tank+i_channel)->GetMaximum() > max_freq) max_freq = hChannels_freq.at(i_slot*num_channels_tank+i_channel)->GetMaximum();
      leg_freq->AddEntry(hChannels_freq.at(i_slot*num_channels_tank+i_channel),ss_channel.str().c_str());
    }
    hChannels_freq.at(i_slot*num_channels_tank+i_first_ch)->GetYaxis()->SetRangeUser(0.,max_freq);
    leg_freq->Draw();
    ss_canvas_freq.str("");
    ss_canvas_freq << outpath << "PMT_Freq_Cr"<<crate_num<<"_Sl"<<slot_num<<".jpg";
    canvas_Channels_freq.at(i_slot)->SaveAs(ss_canvas_freq.str().c_str());
  }

  ss_title_hist_freq.str("");
  ss_title_hist_freq << "Freq RWM (last File) "<<end_time.str();
  hChannels_RWM->SetTitle(ss_title_hist_freq.str().c_str());
  canvas_Channels_freq.at(0)->cd();
  canvas_Channels_freq.at(0)->Clear();
  hChannels_RWM->Draw();
  ss_canvas_freq.str("");
  ss_canvas_freq << outpath << "PMT_Freq_RWM.jpg";
  canvas_Channels_freq.at(0)->SaveAs(ss_canvas_freq.str().c_str());

  ss_title_hist_freq.str("");
  ss_title_hist_freq << "Freq BRF (last File) "<<end_time.str();
  hChannels_RWM->SetTitle(ss_title_hist_freq.str().c_str());
  canvas_Channels_freq.at(0)->Clear();
  hChannels_BRF->Draw();
  ss_canvas_freq.str("");
  ss_canvas_freq << outpath << "PMT_Freq_BRF.jpg";
  canvas_Channels_freq.at(0)->SaveAs(ss_canvas_freq.str().c_str());

  
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    for (int i_channel = 0; i_channel < num_channels_tank; i_channel++){
      if (std::find(vec_disabled_global.begin(),vec_disabled_global.end(),i_slot*num_channels_tank+i_channel)!=vec_disabled_global.end()) continue;
      hChannels_freq.at(i_slot*num_channels_tank+i_channel)->Reset();
      hChannels_freq.at(i_slot*num_channels_tank+i_channel)->SetMaximum(-1111); //needed to reset max/minimum value computation for histogram
    }
  }

  hChannels_RWM->Reset();
  hChannels_BRF->Reset();

}

void MonitorTankTime::DrawFIFOPlots(){

  Log("MonitorTankTime: DrawFIFOPlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawFIFOPlots -----------------------
  //-------------------------------------------------------

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  int max_fifo1 = 0;
  int max_fifo2 = 0;
  int min_fifo1 = 99999;
  int min_fifo2 = 99999;

  for (unsigned int i_card = 0; i_card < fifo1.size(); i_card++){
    int i_crate,i_slot;
    CardIDToElectronicsSpace(fifo1.at(i_card),i_crate,i_slot);
    for (int i_ch=0; i_ch < num_channels_tank; i_ch++){
      int x = i_slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank-i_ch;
      h2D_fifo1->Fill(x-1,y-1);
    }
  }

  for (unsigned int i_card = 0; i_card < fifo2.size(); i_card++){
    int i_crate,i_slot;
    CardIDToElectronicsSpace(fifo2.at(i_card),i_crate,i_slot);
    for (int i_ch=0; i_ch < num_channels_tank; i_ch++){
      int x = i_slot;
      int y = num_channels_tank + (3-i_crate)*num_channels_tank-i_ch;
      h2D_fifo2->Fill(x-1,y-1);
    } 
  } 

  for (int i_crate = 0; i_crate < num_crates_tank; i_crate++){
    for (int i_slot = 0; i_slot < num_slots_tank; i_slot++){
	if (h2D_fifo1->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank) > max_fifo1) max_fifo1 = h2D_fifo1->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank);
	if (h2D_fifo1->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank) < min_fifo1) min_fifo1 = h2D_fifo1->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank);
	if (h2D_fifo2->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank) > max_fifo2) max_fifo2 = h2D_fifo2->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank);
	if (h2D_fifo2->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank) < min_fifo2) min_fifo2 = h2D_fifo2->GetBinContent(i_slot+1,num_channels_tank+(2-i_crate)*num_channels_tank);
    }
  }

  ss_title_fifo1.str("");
  ss_title_fifo1 << "FIFO Overflow Error I (last File) "<<end_time.str();
  h2D_fifo1->SetTitle(ss_title_fifo1.str().c_str());

  TPad *p_fifo1 = (TPad*) canvas_fifo->cd();
  h2D_fifo1->SetStats(0);
  h2D_fifo1->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_fifo1->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_fifo1->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_fifo1->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_fifo1->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_fifo1->LabelsOption("v");
  h2D_fifo1->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  line1->Draw("same");
  line2->Draw("same");
  label_fifo->Draw();
  p_fifo1 ->Update();

   if (h2D_fifo1->GetEntries()>0.){
    if (abs(max_fifo1-min_fifo1)==0) h2D_fifo1->GetZaxis()->SetRangeUser(min_fifo1-1,max_fifo1+1);
    else h2D_fifo1->GetZaxis()->SetRangeUser(1e-6,max_fifo1+0.5);
    TPaletteAxis *palette =
    (TPaletteAxis*)h2D_fifo1->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }
  p_fifo1->Update();
  std::stringstream ss_fifo1;
  ss_fifo1<<outpath<<"PMT_FIFO1_Electronics_lastFile.jpg";
  canvas_fifo->SaveAs(ss_fifo1.str().c_str());
  Log("MonitorTankTime: Output path FIFO I plot (Electronics Space): "+ss_fifo1.str(),v_message,verbosity);

  h2D_fifo1->Reset();
  canvas_fifo->Clear();

  ss_title_fifo2.str("");
  ss_title_fifo2 << "FIFO Overflow Error II (last File) "<<end_time.str();
  h2D_fifo2->SetTitle(ss_title_fifo2.str().c_str());

  TPad *p_fifo2 = (TPad*) canvas_fifo->cd();
  h2D_fifo2->SetStats(0);
  h2D_fifo2->GetXaxis()->SetNdivisions(num_slots_tank);
  h2D_fifo2->GetYaxis()->SetNdivisions(num_crates_tank*num_channels_tank);
  p_fifo2->SetGrid();
  for (int i_label = 0; i_label < int(num_slots_tank); i_label++){
    std::stringstream ss_slot;
    ss_slot<<(i_label+1);
    std::string str_slot = "slot "+ss_slot.str();
    h2D_fifo2->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
  }
  for (int i_label=0; i_label < int(num_channels_tank*num_crates_tank); i_label++){
    std::stringstream ss_ch;
    if (i_label < 4) ss_ch<<((3-i_label)%4)+1;
    else if (i_label < 8) ss_ch<<((7-i_label)%4)+1;
    else ss_ch<<((11-i_label)%4)+1;
    std::string str_ch = "ch "+ss_ch.str();
    h2D_fifo2->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
  }
  h2D_fifo2->LabelsOption("v");
  h2D_fifo2->Draw("colz");

  for (unsigned int i_box =0; i_box < vector_box_inactive.size(); i_box++){
    vector_box_inactive.at(i_box)->Draw("same");
  }

  text_crate1->Draw();
  text_crate2->Draw();
  text_crate3->Draw();
  line1->Draw("same");
  line2->Draw("same");
  label_fifo->Draw();
  p_fifo2 ->Update();
   if (h2D_fifo2->GetEntries()>0.){
    if (abs(max_fifo2-min_fifo2)==0) h2D_fifo2->GetZaxis()->SetRangeUser(min_fifo2-1,max_fifo2+1);
    else h2D_fifo2->GetZaxis()->SetRangeUser(1e-6,max_fifo2+0.5);
    TPaletteAxis *palette =
    (TPaletteAxis*)h2D_fifo2->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }
  p_fifo2->Update();
  std::stringstream ss_fifo2;
  ss_fifo2<<outpath<<"PMT_FIFO2_Electronics_lastFile.jpg";
  canvas_fifo->SaveAs(ss_fifo2.str().c_str());
  Log("MonitorTankTime: Output path FIFO I plot (Electronics Space): "+ss_fifo2.str(),v_message,verbosity);

  h2D_fifo2->Reset();
  canvas_fifo->Clear();

}

void MonitorTankTime::DrawVMEHistogram(){

  Log("MonitorTankTime: DrawVMEHistogram",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawVMEHistogram---------------------
  //-------------------------------------------------------

    boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  canvas_vme->cd();
  canvas_vme->Clear();
  canvas_vme->SetLogy();

  // Find clusters of hits

  // Plot histograms
  std::stringstream ss_vme_hist_title;
  ss_vme_hist_title <<"VME "<<end_time.str()<<" (last file) ";
  hist_vme->SetTitle(ss_vme_hist_title.str().c_str());
  hist_vme->Draw();
  std::stringstream ss_vme_hist;
  ss_vme_hist << outpath << "VMEHist_lastFile."<<img_extension;
  canvas_vme->SaveAs(ss_vme_hist.str().c_str());
  hist_vme->Reset();
  canvas_vme->Clear();
  std::vector<int> coinc_times_all;

  for (int i_ev = 0; i_ev < (int) channels_times.size(); i_ev++){

    int num_pmts = 0;
    std::vector<std::vector<int>> temp_times;
    //std::cout <<"i_ev: "<<i_ev<<std::endl;
    std::vector<std::vector<int>> channels_pmt_times = channels_times.at(i_ev);
    for (int i_ch = 0; i_ch < (int) channels_pmt_times.size(); i_ch++){
      std::vector<int> temp_times_single;
      std::vector<int> ch_hits = channels_pmt_times.at(i_ch);
      if (ch_hits.size()>0) {
        //std::cout <<"Channel "<<i_ch<<", hits: ";
        num_pmts++;
      }
      for (int i_hit = 0; i_hit < (int) ch_hits.size(); i_hit++){
        int hit = ch_hits.at(i_hit);
        //std::cout <<hit<<",";
        if (i_hit == 0) temp_times_single.push_back(ch_hits.at(i_hit));
        if (i_hit != 0 && fabs(ch_hits.at(i_hit)-ch_hits.at(i_hit-1)) > 100) temp_times_single.push_back(ch_hits.at(i_hit));
      }     
      if (ch_hits.size() > 0) {
        //std::cout << std::endl;
        std::sort(temp_times_single.begin(),temp_times_single.end());
        temp_times.push_back(temp_times_single);
    }
  }

  int vme_start=0;
  int min_cluster=3;
  int coincident_time = 50;
  std::vector<int> coinc_ch;

  if (int(temp_times.size()) >=min_cluster+1){

    //std::cout <<"Number of PMTs: "<<temp_times.size()<<std::endl;

    for (int i_temp=0; i_temp < (int) temp_times.size(); i_temp++){
      //std::cout <<"i_temp: "<<i_temp<<std::endl;
      
      for (int i=0; i< (int) temp_times.at(i_temp).size(); i++){
        std::vector<int> coinc_times;
        int num_coincident_pmts=0;
        //std::cout <<"i: "<<i<<", time: "<<temp_times.at(i_temp).at(i)<<std::endl;
        vme_start = temp_times.at(i_temp).at(i);
        for (int i_other=0; i_other < (int) temp_times.size() - i_temp-1; i_other++){
          int current = i_temp+i_other+1;
          for (int i_hit=0; i_hit < (int) temp_times.at(current).size(); i_hit++){
            if (fabs(vme_start - temp_times.at(current).at(i_hit)) < coincident_time){
              if (std::find(coinc_ch.begin(),coinc_ch.end(),current)==coinc_ch.end()){
                num_coincident_pmts++;
                coinc_ch.push_back(current);
                coinc_times.push_back(temp_times.at(current).at(i_hit));
              }
            }
          }
        }
        if (num_coincident_pmts>=min_cluster){
        //std::cout <<"num_coincident_pmts: "<<num_coincident_pmts<<std::endl;
          if (vme_start < 2000) {
            hist_vme_cluster->Fill(vme_start);
            coinc_times_all.push_back(vme_start);
          }
          for (int i_cl=0; i_cl < (int) coinc_times.size(); i_cl++){
            if (coinc_times.at(i_cl) < 2000){
		//std::cout <<"Fill "<<coinc_times.at(i_cl)<<std::endl;
		hist_vme_cluster->Fill(coinc_times.at(i_cl));  
                coinc_times_all.push_back(coinc_times.at(i_cl));
	    }
          }
        }
      }

    }
  }

  }

  ss_vme_hist_title.str("");
  ss_vme_hist_title <<"VME Cluster "<<end_time.str()<<" (last file) ";
  hist_vme_cluster->SetTitle(ss_vme_hist_title.str().c_str());
  hist_vme_cluster->Draw();
  ss_vme_hist.str("");
  ss_vme_hist << outpath << "VMEHist_Cluster_lastFile."<<img_extension;
  canvas_vme->SaveAs(ss_vme_hist.str().c_str());
  hist_vme_cluster->Reset();
  canvas_vme->Clear();

  overall_coinc_times.push_back(coinc_times_all);
  
  if (overall_coinc_times.size() > 20){
    overall_coinc_times.erase(overall_coinc_times.begin());
  }

  for (int i_coinc = 0; i_coinc < (int) overall_coinc_times.size(); i_coinc++){
    std::vector<int> single_coinc_times = overall_coinc_times.at(i_coinc);
    for (int i_hit=0; i_hit < (int) single_coinc_times.size(); i_hit++){
      hist_vme_cluster_20->Fill(single_coinc_times.at(i_hit));
    }
  }

  ss_vme_hist_title.str("");
  ss_vme_hist_title <<"VME Cluster "<<end_time.str()<<" (last 20 files) ";
  hist_vme_cluster_20->SetTitle(ss_vme_hist_title.str().c_str());
  hist_vme_cluster_20->Draw();
  ss_vme_hist.str("");
  ss_vme_hist << outpath << "VMEHist_Cluster_last20Files."<<img_extension;
  canvas_vme->SaveAs(ss_vme_hist.str().c_str());
  hist_vme_cluster_20->Reset();
  canvas_vme->Clear();

}

void MonitorTankTime::DrawFileHistory(ULong64_t timestamp_end, double time_frame, std::string file_ending, int _linewidth){

  Log("MonitorTankTime: DrawFileHistory",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawFileHistory ---------------------
  //-------------------------------------------------------

  //Creates a plot showing the time stamps for all the files within the last time_frame mins
  //The plot is updated with the update_frequency specified in the configuration file (default: 5 mins)

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  timestamp_end += utc_to_t;

  ULong64_t timestamp_start = timestamp_end - time_frame*MSEC_to_SEC*SEC_to_MIN*MIN_to_HOUR;
  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  canvas_logfile_tank->cd();
  log_files->SetBins(num_files_history,timestamp_start/MSEC_to_SEC,timestamp_end/MSEC_to_SEC);
  log_files->GetXaxis()->SetTimeOffset(0.);
  log_files->Draw();

  std::stringstream ss_title_filehistory;
  ss_title_filehistory << "PMT Files History (last "<<ss_timeframe.str()<<"h)";
    
  log_files->SetTitle(ss_title_filehistory.str().c_str());

  std::vector<TLine*> file_markers;
  for (unsigned int i_file = 0; i_file < tend_plot.size(); i_file++){
    TLine *line_file = new TLine((tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,0.,(tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,1.);
    line_file->SetLineColor(1);
    line_file->SetLineStyle(1);
    line_file->SetLineWidth(_linewidth);
    line_file->Draw("same"); 
    file_markers.push_back(line_file);
  }

  std::stringstream ss_logfiles;
  ss_logfiles << outpath << "PMT_FileHistory_" << file_ending << "." << img_extension;
  canvas_logfile_tank->SaveAs(ss_logfiles.str().c_str());

  for (unsigned int i_line = 0; i_line < file_markers.size(); i_line++){
    delete file_markers.at(i_line);
  }

  log_files->Reset();
  canvas_logfile_tank->Clear();


}

void MonitorTankTime::PrintFileTimeStamp(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  if (verbosity > 2) std::cout <<"MonitorTankTime: PrintFileTimeStamp"<<std::endl;

  //-------------------------------------------------------
  //-----------------PrintFileTimeStamp--------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << "Current time: "<<endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  TText *label_lastfile = nullptr;
  label_lastfile = new TText(0.04,0.66,end_time.str().c_str());
  label_lastfile->SetNDC(1);
  label_lastfile->SetTextSize(0.1);

  TLatex *label_timediff = nullptr;

  if (tend_plot.size() == 0){
    std::stringstream time_diff;
    time_diff << "#Delta t Last PMT File: >"<<time_frame<<"h";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  } else {
    ULong64_t timestamp_lastfile = tend_plot.at(tend_plot.size()-1);
    boost::posix_time::ptime filetime = *Epoch + boost::posix_time::time_duration(int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_lastfile/MSEC_to_SEC/1000.)%60,timestamp_lastfile%1000);
    boost::posix_time::time_duration t_since_file= boost::posix_time::time_duration(endtime - filetime);
    int t_since_file_min = int(t_since_file.total_milliseconds()/MSEC_to_SEC/SEC_to_MIN);
    std::stringstream time_diff;
    time_diff << "#Delta t Last PMT File: "<<t_since_file_min/60<<"h:"<<t_since_file_min%60<<"min";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  }

  canvas_file_timestamp_tank->cd();
  label_lastfile->Draw();
  label_timediff->Draw();
  std::stringstream ss_file_timestamp;
  ss_file_timestamp << outpath << "PMT_FileTimeStamp_" << file_ending << "." << img_extension;
  canvas_file_timestamp_tank->SaveAs(ss_file_timestamp.str().c_str());

  delete label_lastfile;
  delete label_timediff;

  canvas_file_timestamp_tank->Clear();

}

void MonitorTankTime::CardIDToElectronicsSpace(int CardID, 
        int &CrateNum, int &SlotNum)
{
  //CardID = CrateNum * 1000 + SlotNum.  This logic works if we have less than
  // 10 crates and less than 100 Slots (which we do).
  SlotNum = CardID % 100;
  CrateNum = CardID / 1000;
  return;
}

std::string MonitorTankTime::convertTimeStamp_to_Date(ULong64_t timestamp){

  //format of date is YYYY_MM-DD

  boost::posix_time::ptime convertedtime = *Epoch + boost::posix_time::time_duration(int(timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp/MSEC_to_SEC/1000.)%60,timestamp%1000);
  struct tm converted_tm = boost::posix_time::to_tm(convertedtime);
  std::stringstream ss_date;
  ss_date << converted_tm.tm_year+1900 << "_" << converted_tm.tm_mon+1 << "-" <<converted_tm.tm_mday;
  return ss_date.str();

}


bool MonitorTankTime::does_file_exist(std::string filename){

  std::ifstream infile(filename.c_str());
  bool file_good = infile.good();
  infile.close();
  return file_good;

}


