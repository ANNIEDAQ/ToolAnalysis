#include "MonitorMRDTime.h"

MonitorMRDTime::MonitorMRDTime():Tool(){}

bool MonitorMRDTime::Initialise(std::string configfile, DataModel &data){


  //-------------------------------------------------------
  //---------------Initialise config file------------------
  //-------------------------------------------------------

  if(configfile!="")  m_variables.Initialise(configfile);
  m_data= &data;

  //gObjectTable only for debugging memory leaks, otherwise comment out
  //std::cout <<"MonitorMRDTime: List of Objects (beginning of Initialise): "<<std::endl;
  //gObjectTable->Print();

  //-------------------------------------------------------
  //-----------------Get Configuration---------------------
  //-------------------------------------------------------

  update_frequency = 0.;

  m_variables.Get("OutputPath",outpath_temp);
  m_variables.Get("ActiveSlots",active_slots);
  m_variables.Get("InActiveChannels",inactive_channels);
  m_variables.Get("LoopbackChannels",loopback_channels);
  m_variables.Get("StartTime",StartTime);
  m_variables.Get("PlotConfiguration",plot_configuration);
  m_variables.Get("UpdateFrequency",update_frequency);
  m_variables.Get("PathMonitoring",path_monitoring);
  m_variables.Get("ImageFormat",img_extension);
  m_variables.Get("ForceUpdate",force_update);
  m_variables.Get("DrawMarker",draw_marker);
  m_variables.Get("DrawSingle",draw_single);
  m_variables.Get("verbose",verbosity);

  if (verbosity > 1) std::cout <<"Tool MonitorMRDTime: Initialising...."<<std::endl;

  //Update frequency specifies the frequency at which the File Log Histogram is updated
  //All other monitor plots are updated as soon as a new file is available for readout
  if (update_frequency < 0.1) {
    if (verbosity > 0) std::cout <<"MonitorMRDTime: Update Frequency of "<<update_frequency<<" mins is too low. Setting default value of 5 mins."<<std::endl;
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
  //check if the path for the monitoring files exists
  /*if (!boost::filesystem::exists(boost::filesystem::path(path_monitoring.c_str()))){
    if (verbosity > 0) std::cout <<"ERROR (MonitorMRDTime): Specified path for the monitoring files -"<<path_monitoring<<"- does not seem to exist. Using default directory -/monitoringfiles/-"<<std::endl;
    path_monitoring = "/monitoringfiles/";
  }*/
  std::cout <<"PathMonitoring: "<<path_monitoring<<std::endl;

  //-------------------------------------------------------
  //-----------------Load ANNIE geometry-------------------
  //-------------------------------------------------------

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  m_data->CStore.Get("MRDCrateSpaceToChannelNumMap",CrateSpaceToChannelNumMap);
  Position center_position = geom->GetTankCentre();
  tank_center_x = center_position.X();
  tank_center_y = center_position.Y();
  tank_center_z = center_position.Z();

  //Set up Epoch for converting the timestamps to actual times
  Epoch = new boost::posix_time::ptime(boost::gregorian::from_string(StartTime));

  //Evaluating output path for monitoring plots
  if (outpath_temp == "fromStore") m_data->CStore.Get("OutPath",outpath);
  else outpath = outpath_temp;
  if (verbosity > 2) std::cout <<"MonitorMRDTime: Output path for plots is "<<outpath<<std::endl;

  //-------------------------------------------------------
  //-----------------Get active slots----------------------
  //-------------------------------------------------------

  num_active_slots=0;
  num_active_slots_cr1=0;
  num_active_slots_cr2=0;

  ifstream file(active_slots.c_str());
  unsigned int temp_crate, temp_slot, temp_channel;
  int loopnum=0;
  while (!file.eof()){
    file>>temp_crate>>temp_slot;
    loopnum++;
    if (file.eof()) break;
      if (temp_crate-min_crate<0 || temp_crate-min_crate>=num_crates) {
        std::cout <<"ERROR (MonitorMRDTime): Specified crate out of range [7...8]. Continue with next entry."<<std::endl;
        continue;
      }
      if (temp_slot<1 || temp_slot>num_slots){
        std::cout <<"ERROR (MonitorMRDTime): Specified slot out of range [1...24]. Continue with next entry."<<std::endl;
        continue;
      }
      active_channel[temp_crate-min_crate][temp_slot-1]=1;		//slot start with number 1 instead of 0, crates with 7
      nr_slot.push_back((temp_crate-min_crate)*100+(temp_slot)); 
      std::vector<unsigned int> temp_crate_slot{temp_crate,temp_slot};
      CrateSlot_to_ActiveSlot.emplace(temp_crate_slot,num_active_slots);

      if (temp_crate-min_crate==0) {
        num_active_slots_cr1++;
        active_slots_cr1.push_back(temp_slot);
      }
      if (temp_crate-min_crate==1) {
        num_active_slots_cr2++;
        active_slots_cr2.push_back(temp_slot);
      }

      ActiveSlot_to_Crate.emplace(num_active_slots,temp_crate);
      ActiveSlot_to_Slot.emplace(num_active_slots,temp_slot);

      for (int i_channel = 0; i_channel < num_channels; i_channel++){
        unsigned int temp_channel = i_channel;
        int total_channel = num_active_slots*num_channels + i_channel;
        TotalChannel_to_Crate.emplace(total_channel,temp_crate);
        TotalChannel_to_Slot.emplace(total_channel,temp_slot);
        TotalChannel_to_Channel.emplace(total_channel,i_channel);
        std::vector<unsigned int> CrateSlotChannel{temp_crate,temp_slot,temp_channel};
        CrateSlotChannel_to_TotalChannel.emplace(CrateSlotChannel,total_channel);
      }

      num_active_slots++;
  }
  file.close();

  //-------------------------------------------------------
  //-------------Get inactive channels---------------------
  //-------------------------------------------------------

  ifstream file_inactive(inactive_channels.c_str());

  while (!file_inactive.eof()){
    file_inactive>>temp_crate>>temp_slot>>temp_channel;
    if (verbosity > 2) std::cout<<temp_crate<<" , "<<temp_slot<<" , "<<temp_channel<<std::endl;
    if (file_inactive.eof()) break;

    if (temp_crate-min_crate<0 || temp_crate-min_crate>=num_crates) {
      std::cout <<"ERROR (MonitorMRDTime): Specified crate out of range [7...8]. Continue with next entry."<<std::endl;
      continue;
    }
    if (temp_slot<1 || temp_slot>num_slots){
      std::cout <<"ERROR (MonitorMRDTime): Specified slot out of range [1...24]. Continue with next entry."<<std::endl;
      continue;
    }
    if (temp_channel<0 || temp_channel > num_channels){
      std::cout <<"ERROR (MonitorMRDTime): Specified channel out of range [0...31]. Continue with next entry."<<std::endl;
      continue;
    }

    if (temp_crate == min_crate){
      inactive_ch_crate1.push_back(temp_channel);
      inactive_slot_crate1.push_back(temp_slot);
    } else if (temp_crate == min_crate+1){
      inactive_ch_crate2.push_back(temp_channel);
      inactive_slot_crate2.push_back(temp_slot);
    } else {
      std::cout <<"ERROR (MonitorMRDTime): Crate # out of range, entry ("<<temp_crate<<"/"<<temp_slot<<"/"<<temp_channel<<") not added to inactive channel configuration." <<std::endl;
    }
  }
  file_inactive.close();

  //-------------------------------------------------------
  //-------------Get loopback channels---------------------
  //-------------------------------------------------------

  ifstream file_loopback(loopback_channels.c_str());
  std::string temp_loopback_name;

  while (!file_loopback.eof()){
    file_loopback>>temp_loopback_name>>temp_crate>>temp_slot>>temp_channel;
    if (verbosity > 2) std::cout<<temp_loopback_name<<": "<<temp_crate<<" , "<<temp_slot<<" , "<<temp_channel<<std::endl;
    if (file_loopback.eof()) break;

    if (temp_crate-min_crate<0 || temp_crate-min_crate>=num_crates) {
      std::cout <<"ERROR (MonitorMRDTime): Specified loopback crate out of range [7...8]. Continue with next entry."<<std::endl;
      continue;
    }
    if (temp_slot<1 || temp_slot>num_slots){
      std::cout <<"ERROR (MonitorMRDTime): Specified loopback slot out of range [1...24]. Continue with next entry."<<std::endl;
      continue;
    }
    if (temp_channel<0 || temp_channel > num_channels){
      std::cout <<"ERROR (MonitorMRDTime): Specified loopback channel out of range [0...31]. Continue with next entry."<<std::endl;
      continue;
    }
    loopback_name.push_back(temp_loopback_name);
    loopback_crate.push_back(temp_crate);
    loopback_slot.push_back(temp_slot);
    loopback_channel.push_back(temp_channel);
    std::vector<unsigned int> CrateSlotChannel{temp_crate,temp_slot,temp_channel};
    int i_active = CrateSlotChannel_to_TotalChannel[CrateSlotChannel];
    if (temp_loopback_name == "Beam") {
	beam_ch = i_active;
    }
    else if (temp_loopback_name == "Cosmic") {  
      cosmic_ch = i_active;
    }
  }
  file_loopback.close();

  //-------------------------------------------------------
  //----------Initialize storing containers----------------
  //-------------------------------------------------------

  InitializeVectors();

  //-------------------------------------------------------
  //----------Read in configuration option for plots-------
  //-------------------------------------------------------

  ReadInConfiguration();

  //-------------------------------------------------------
  //------Setup time variables for periodic updates--------
  //-------------------------------------------------------

  period_update = boost::posix_time::time_duration(0,int(update_frequency),0,0);
  last = boost::posix_time::ptime(boost::posix_time::second_clock::local_time());

  //Keep track whether there has been a MRDdata file or not
  bool_mrddata = false;

  // Omit warning messages from ROOT: 1001 - info messages, 2001 - warnings, 3001 - errors
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  
  return true;
}


bool MonitorMRDTime::Execute(){

  if (verbosity > 10) std::cout <<"MonitorMRDTime: Executing ...."<<std::endl;

  current = (boost::posix_time::second_clock::local_time());
  duration = boost::posix_time::time_duration(current - last);
  current_stamp_duration = boost::posix_time::time_duration(current - *Epoch);
  current_stamp = current_stamp_duration.total_milliseconds();
  utc = (boost::posix_time::second_clock::universal_time());
  current_utc_duration = boost::posix_time::time_duration(utc-current);
  current_utc = current_utc_duration.total_milliseconds();
  utc_to_t = (ULong64_t) current_utc;

  if (verbosity > 10) std::cout <<"MRDMonitorTime: "<<duration.total_milliseconds()/MSEC_to_SEC/SEC_to_MIN<<" mins since last time plot"<<std::endl;

  //-------------------------------------------------------
  //---------Checking the state of MRD data stream---------
  //-------------------------------------------------------

  bool has_cc;
  m_data->CStore.Get("HasCCData",has_cc);

  std::string State;
  m_data->CStore.Get("State",State);

  if (has_cc){
   if (State == "MRDSingle" || State == "Wait"){

    //MRDMonitorLive is executed
    if (verbosity > 2) std::cout <<"MRDMonitorTime: State is "<<State<<std::endl;

   } else if (State == "DataFile"){

    //MRDMonitorTime executed
   	if (verbosity > 1) std::cout<<"MRDMonitorTime: New data file available."<<std::endl;

    //Setting print to false is necessary in order for it to work properly
    m_data->Stores["CCData"]->Get("FileData",MRDdata);
    MRDdata->Print(false);
    bool_mrddata = true;

    //Read in information from MRD file store, fill into the storing containers (vectors)
    ReadInData();

    //Write the event information to a file
    //TODO: change this to a database later on!
    //Check if data has already been written included in WriteToFile function
    WriteToFile();

    //Plot plots only associated to current file
    DrawLastFilePlots();

    //Draw customly defined plots
    UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

   }else {

   	if (verbosity > 1) std::cout <<"MRDMonitorTime: State not recognized: "<<State<<std::endl;

   }
  }
  //-------------------------------------------------------------
  //---------------Draw customly defined plots-------------------
  //-------------------------------------------------------------

   // if force_update is specified, the plots will be updated no matter whether there has been a new file or not

   if (force_update) UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

  //-------------------------------------------------------------
  //---Has enough time passed for updating File history plot?----
  //-------------------------------------------------------------

  if(duration>=period_update){
    last=current;
    DrawFileHistory(current_stamp,24.,"current_24h",1);     //show 24h history of MRD files
    PrintFileTimeStamp(current_stamp,24.,"current_24h");
    DrawFileHistory(current_stamp,2.,"current_2h",3);     //show 2h history of MRD files
  }

  
  //gObjectTable only for debugging memory leaks, otherwise comment out
  //std::cout <<"List of Objects (after execute step): "<<std::endl;
  //gObjectTable->Print();

  return true;

}


bool MonitorMRDTime::Finalise(){

  if (verbosity > 1) std::cout <<"Tool MonitorMRDTime: Finalising ...."<<std::endl;

  //if (bool_mrddata) MRDdata->Delete();

  //delete all the pointer to objects that are still active

  //timing pointers
  delete Epoch;

  //other objects
  delete label_cr1;
  delete label_cr2;
  delete label_rate_cr1;
  delete label_rate_cr2;
  delete label_rate_facc;
  delete separate_crates;
  for (unsigned int i_box = 0; i_box < vector_box_inactive.size(); i_box++){
    delete vector_box_inactive.at(i_box);
  }
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    for (unsigned int i_box = 0; i_box < vector_box_inactive_hitmap[i_slot].size(); i_box++){
      delete vector_box_inactive_hitmap[i_slot].at(i_box);
    }
  }

  //legends
  delete leg_trigger;
  delete leg_noloopback;
  delete leg_eventtypes;
  delete leg_scatter;
  delete leg_scatter_trigger;
  delete leg_tdc;
  delete leg_rms;
  delete leg_rate;

  //delete pie charts
  delete pie_triggertype;
  delete pie_weirdevents;


  //histograms
  delete hist_hitmap_cr1;
  delete hist_hitmap_cr2;
  for (unsigned int i=0; i<hist_hitmap_slot.size();i++){
    delete hist_hitmap_slot.at(i);
  }
  for (unsigned int i_scatter = 0; i_scatter < hist_scatter.size(); i_scatter++){
    delete hist_scatter.at(i_scatter);
  }
  delete hist_tdc;
  delete hist_tdc_cluster;
  delete hist_tdc_cluster_20;
  delete hist_tdc_beam;
  delete hist_tdc_cluster_beam;
  delete hist_tdc_cluster_beam_20;
  delete hist_tdc_cosmic;
  delete hist_tdc_cluster_cosmic;
  delete hist_tdc_cluster_cosmic_20;
  delete log_files_mrd;
  delete rate_crate1;
  delete rate_crate2;
  delete rate_top;
  delete rate_side;
  delete rate_facc;

  //graphs
  for (unsigned int i_ch=0; i_ch<gr_tdc.size();i_ch++){
    delete gr_tdc.at(i_ch);
    delete gr_rms.at(i_ch);
    delete gr_rate.at(i_ch);
  }
  for (unsigned int i_trigger = 0; i_trigger < gr_trigger.size(); i_trigger++){
    delete gr_trigger.at(i_trigger);
  }
  delete gr_noloopback;
  delete gr_zerohits;
  delete gr_doublehits;

  //multigraphs
  delete multi_ch_tdc;
  delete multi_ch_rms;
  delete multi_ch_rate;
  delete multi_trigger;
  delete multi_eventtypes;

  //canvases
  delete canvas_hitmap;
  delete canvas_hitmap_slot;
  delete canvas_logfile_mrd;
  delete canvas_ch_tdc;
  delete canvas_ch_rms;
  delete canvas_ch_rate;
  delete canvas_ch_single;
  delete canvas_scatter;
  delete canvas_scatter_single;
  delete canvas_tdc;
  delete canvas_trigger;
  delete canvas_trigger_time;
  delete canvas_rate_electronics;
  delete canvas_rate_physical;
  delete canvas_rate_physical_facc;
  delete canvas_pie;
  delete canvas_file_timestamp;

  return true;
}

void MonitorMRDTime::ReadInConfiguration(){

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
        if (verbosity > 0) std::cout <<"ERROR (MonitorMRDTime): ReadInConfiguration: Need at least 4 arguments in one line: TimeFrame - TimeEnd - FileLabel - PlotType1. Please look over the configuration file and adjust it accordingly."<<std::endl;
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
    if (verbosity > 0) std::cout <<"ERROR (MonitorMRDTime): ReadInConfiguration: Could not open file "<<plot_configuration<<"! Check if path is valid..."<<std::endl;
  }
  file.close();


  if (verbosity > 2){
    std::cout <<"---------------------------------------------------------------------"<<std::endl;
    std::cout <<"MonitorMRDTime: ReadInConfiguration: Read in the following data into configuration variables: "<<std::endl;
    for (unsigned int i_t=0; i_t < config_timeframes.size(); i_t++){
      std::cout <<config_timeframes.at(i_t)<<", "<<config_endtime.at(i_t)<<", "<<config_label.at(i_t)<<", ";
      for (unsigned int i_plot = 0; i_plot < config_plottypes.at(i_t).size(); i_plot++){
        std::cout <<config_plottypes.at(i_t).at(i_plot)<<", ";
      }
      std::cout<<std::endl;
    }
    std::cout <<"-----------------------------------------------------------------------"<<std::endl;
  }

  if (verbosity > 2) std::cout <<"MonitorMRDTime: ReadInConfiguration: Parsing dates: "<<std::endl;
  for (unsigned int i_date = 0; i_date < config_endtime.size(); i_date++){
    if (config_endtime.at(i_date) == "TEND_LASTFILE") {
      if (verbosity > 2) std::cout <<"TEND_LASTFILE: Starting from end of last read-in file"<<std::endl;
      ULong64_t zero = 0; 
      config_endtime_long.push_back(zero);
    } else if (config_endtime.at(i_date).size()==15){
        boost::posix_time::ptime spec_endtime(boost::posix_time::from_iso_string(config_endtime.at(i_date)));
        boost::posix_time::time_duration spec_endtime_duration = boost::posix_time::time_duration(spec_endtime - *Epoch);
        ULong64_t spec_endtime_long = spec_endtime_duration.total_milliseconds();
        config_endtime_long.push_back(spec_endtime_long);
    } else {
      if (verbosity > 2) std::cout <<"Specified end date "<<config_endtime.at(i_date)<<" does not have the desired format YYYYMMDDTHHMMSS. Please change the format in the config file in order to use this tool. Starting from end of last file"<<std::endl;
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    }
  }


}

void MonitorMRDTime::ReadInData(){

  //-------------------------------------------------------
  //----------------ReadInData ----------------------------
  //-------------------------------------------------------

  n_doublehits = 0;
  n_zerohits = 0;
  n_noloopback = 0;
  std::vector<int> vector_nhits;

  //std::cout <<"Reading in Data"<<std::endl;

  for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){
    tdc_file.at(i_channel).clear();
    timestamp_file.at(i_channel).clear();
  }
  tdc_file_times.clear();
  tdc_file_ch.clear();

  if (verbosity > 1) std::cout <<"MonitorMRDTime: ReadInData..."<<std::endl;
  t_file_end = 0;
  int total_number_entries;  
  MRDdata->Header->Get("TotalEntries",total_number_entries);
  if (verbosity > 1) std::cout <<"MonitorMRDTime: MRDdata file total entries: "<<total_number_entries<<std::endl;
  int count=0;

  for (int i_event = 0; i_event <  total_number_entries; i_event++){

    //initialize vector for event times
    std::vector<int> tdc_file_times_single;
    std::vector<int> tdc_file_ch_single;

    //get MRDout data
    MRDdata->GetEntry(i_event);
    MRDdata->Get("Data",MRDout);

    //get event time
    boost::posix_time::ptime eventtime;
    ULong64_t timestamp = MRDout.TimeStamp;
    eventtime = *Epoch + boost::posix_time::time_duration(int(timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp/MSEC_to_SEC)%60,timestamp%1000);   //default last entry in miliseconds, need special compilation option for nanosec hh:mm:ss:msms

    //giving information about the time frame of data in current MRD file
    if (verbosity > 1){
      if (i_event == 0) std::cout <<"MonitorMRDTime: File time range: "<<eventtime.date()<<","<<eventtime.time_of_day();
      else if (i_event == total_number_entries-1) std::cout <<" ... "<<eventtime.date()<<","<<eventtime.time_of_day()<<std::endl;
    }

    //get end and start times for processed data file
    if (i_event == 0) t_file_start = MRDout.TimeStamp;
    else if (i_event == total_number_entries-1) t_file_end = MRDout.TimeStamp;

    //printing intormation about the event
    if (verbosity > 3){
      std::cout <<"------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      std::cout <<"Entry: "<<i_event<<", TimeStamp: "<<timestamp<<std::endl;
      std::cout <<"Slot size: "<<MRDout.Slot.size()<<", Crate size: "<<MRDout.Crate.size()<<", Channel size: "<<MRDout.Channel.size()<<std::endl;
      std::cout <<"OutN: "<<MRDout.OutN<<", Trigger: "<<MRDout.Trigger<<", Type size: "<<MRDout.Type.size()<<std::endl;
    }

    if (MRDout.Slot.size() == 0) n_zerohits++;

    //looping over all channels in event
    bool no_loopback = true;
    vector_nhits.assign(num_active_slots*num_channels,0.);

    for (unsigned int i_entry = 0; i_entry < MRDout.Slot.size(); i_entry++){

      //print out information if needed for debugging hardware
      if (verbosity > 3) std::cout <<"MonitorMRDTime: Channel entry: Crate "<<MRDout.Crate.at(i_entry)<<", Slot "<<MRDout.Slot.at(i_entry)<<", Channel "<<MRDout.Channel.at(i_entry)<<", TDC value: "<<MRDout.Value.at(i_entry)<<std::endl;

      int active_slot_nr;
      std::vector<int>::iterator it = std::find(nr_slot.begin(), nr_slot.end(), (MRDout.Slot.at(i_entry))+(MRDout.Crate.at(i_entry)-min_crate)*100);
      if (it == nr_slot.end()){
        if (verbosity > 0){
          std::cout <<"WARNING (MonitorMRDTime): Read-out Crate/Slot/Channel number not active according to configuration file. Check the configfile to process the data..."<<std::endl;
          std::cout <<"WARNING (MonitorMRDTime): Crate: "<<MRDout.Crate.at(i_entry)<<", Slot: "<<MRDout.Slot.at(i_entry)<<std::endl;
        }
        continue;
      }
      count++;
      active_slot_nr = std::distance(nr_slot.begin(),it);
      std::vector<unsigned int> crate_slot_temp{MRDout.Crate.at(i_entry),MRDout.Slot.at(i_entry)};
      int check_active_slot = CrateSlot_to_ActiveSlot[crate_slot_temp];

      //fill data in live vectors
      int ch = active_slot_nr*num_channels+MRDout.Channel.at(i_entry);
      tdc_file.at(ch).push_back(MRDout.Value.at(i_entry));
      timestamp_file.at(ch).push_back(MRDout.TimeStamp);
      vector_nhits[ch]++;
      tdc_file_times_single.push_back(MRDout.Value.at(i_entry));
      tdc_file_ch_single.push_back(ch);

      for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
        if (MRDout.Crate.at(i_entry) == loopback_crate.at(i_trigger) && MRDout.Slot.at(i_entry) == loopback_slot.at(i_trigger) && MRDout.Channel.at(i_entry) == loopback_channel.at(i_trigger)) no_loopback = false;
      }
    }

    if (no_loopback) n_noloopback++;
    for (int i_ch = 0; i_ch < num_active_slots*num_channels; i_ch++){
      if (vector_nhits[i_ch] > 1) n_doublehits++;
    }

    tdc_file_times.push_back(tdc_file_times_single);
    tdc_file_ch.push_back(tdc_file_ch_single);

    //clear MRDout member vectors afterwards
    MRDout.Value.clear();
    MRDout.Slot.clear();
    MRDout.Channel.clear();
    MRDout.Crate.clear();
    MRDout.Type.clear();
  }

  n_normalhits = total_number_entries - n_zerohits - n_doublehits;

}

void MonitorMRDTime::WriteToFile(){

  if (verbosity > 2) std::cout <<"MonitorMRDTime: WriteToFile..."<<std::endl;

  std::string file_start_date = convertTimeStamp_to_Date(t_file_start);
  std::stringstream root_filename;
  root_filename << path_monitoring << "MRD_" << file_start_date <<".root";

  if (verbosity > 2) std::cout <<"root filename: "<<root_filename.str() <<std::endl;

  std::string root_option = "RECREATE";
  if (does_file_exist(root_filename.str())) root_option = "UPDATE";
  TFile *f = new TFile(root_filename.str().c_str(),root_option.c_str());

  ULong64_t t_start, t_end, t_frame;
  std::vector<unsigned int> *crate = new std::vector<unsigned int>;
  std::vector<unsigned int> *slot = new std::vector<unsigned int>;
  std::vector<unsigned int> *channel = new std::vector<unsigned int>;
  std::vector<double> *tdc = new std::vector<double>;
  std::vector<double> *rms = new std::vector<double>;
  std::vector<double> *rate = new std::vector<double>;
  std::vector<int> *channelcount = new std::vector<int>;
  double rate_cosmic;
  double rate_beam;
  double rate_noloopback;
  double rate_normalhit;
  double rate_doublehit;
  double rate_zerohits;
  int nevents;

  TTree *t;
  if (f->GetListOfKeys()->Contains("mrdmonitor_tree")) {
    if (verbosity > 2) std::cout <<"MonitorMRDTime: WriteToFile: Tree already exists"<<std::endl;
    t = (TTree*) f->Get("mrdmonitor_tree");
    t->SetBranchAddress("t_start",&t_start);
    t->SetBranchAddress("t_end",&t_end);
    t->SetBranchAddress("crate",&crate);
    t->SetBranchAddress("slot",&slot);
    t->SetBranchAddress("channel",&channel);
    t->SetBranchAddress("tdc",&tdc);
    t->SetBranchAddress("rms",&rms);
    t->SetBranchAddress("rate",&rate);
    t->SetBranchAddress("channelcount",&channelcount);
    t->SetBranchAddress("rate_cosmic",&rate_cosmic);
    t->SetBranchAddress("rate_beam",&rate_beam);
    t->SetBranchAddress("rate_noloopback",&rate_noloopback);
    t->SetBranchAddress("rate_normalhit",&rate_normalhit);
    t->SetBranchAddress("rate_doublehit",&rate_doublehit);
    t->SetBranchAddress("rate_zerohits",&rate_zerohits);
    t->SetBranchAddress("nevents",&nevents);
  } else {
    t = new TTree("mrdmonitor_tree","MRD Monitoring tree");
    if (verbosity > 2) std::cout <<"MonitorMRDTime: WriteToFile: Tree is created from scratch"<<std::endl;
    t->Branch("t_start",&t_start);
    t->Branch("t_end",&t_end);
    t->Branch("crate",&crate);
    t->Branch("slot",&slot);
    t->Branch("channel",&channel);
    t->Branch("tdc",&tdc);
    t->Branch("rms",&rms);
    t->Branch("rate",&rate);
    t->Branch("channelcount",&channelcount);
    t->Branch("rate_cosmic",&rate_cosmic);
    t->Branch("rate_beam",&rate_beam);
    t->Branch("rate_noloopback",&rate_noloopback);
    t->Branch("rate_normalhit",&rate_normalhit);
    t->Branch("rate_doublehit",&rate_doublehit);
    t->Branch("rate_zerohits",&rate_zerohits);
    t->Branch("nevents",&nevents);
  }

  int n_entries = t->GetEntries();
  bool omit_entries = false;
  for (int i_entry = 0; i_entry < n_entries; i_entry++){
    t->GetEntry(i_entry);
    if (t_start == t_file_start) {
      if (verbosity > 0) std::cout <<"WARNING (MonitorMRDTime): WriteToFile: Wanted to write data from file that is already written to DB. Omit entries"<<std::endl;
      omit_entries = true;
    }
  }
  //if data is already written to DB/File, do not write it again
  if (omit_entries) {

    //don't write file again, but still delete TFile and TTree object!!! (and vectors)
    f->Close();
    delete crate;
    delete slot;
    delete channel;
    delete tdc;
    delete rms;
    delete rate;
    delete channelcount;
    delete f;

    gROOT->cd();

    return;

  } 

  crate->clear();
  slot->clear();
  channel->clear();
  tdc->clear();
  rms->clear();
  rate->clear();
  channelcount->clear();

  t_start = t_file_start;
  t_end = t_file_end;
  t_frame = t_end - t_start;

  rate_doublehit = n_doublehits/(t_frame/MSEC_to_SEC);
  rate_zerohits = n_zerohits/(t_frame/MSEC_to_SEC);
  rate_normalhit = n_normalhits/(t_frame/MSEC_to_SEC);
  rate_noloopback = n_noloopback/(t_frame/MSEC_to_SEC);

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(t_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_start/MSEC_to_SEC/SEC_to_MIN)%60,int(t_start/MSEC_to_SEC/1000.)%60,t_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  if (verbosity > 2) std::cout <<"MonitorMRDTime: WriteToFile: Writing data to file: "<<starttime_tm.tm_year+1900<<"/"<<starttime_tm.tm_mon+1<<"/"<<starttime_tm.tm_mday<<"-"<<starttime_tm.tm_hour<<":"<<starttime_tm.tm_min<<":"<<starttime_tm.tm_sec;
  if (verbosity > 2) std::cout <<"..."<<endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec<<std::endl;
  nevents=0;
  long n_beam=0;
  long n_cosmic=0;

  for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){
    unsigned int crate_temp = TotalChannel_to_Crate[i_channel];
    unsigned int slot_temp = TotalChannel_to_Slot[i_channel];
    unsigned int channel_temp = TotalChannel_to_Channel[i_channel];
    double rate_temp = tdc_file.at(i_channel).size() / (t_frame/MSEC_to_SEC);
    double mean_tdc = 0.;
    double rms_temp = 0.;
    for (unsigned int i_tdc = 0; i_tdc < tdc_file.at(i_channel).size(); i_tdc++){
      mean_tdc+=tdc_file.at(i_channel).at(i_tdc);
    }
    if (tdc_file.at(i_channel).size() > 0) {
      mean_tdc/=tdc_file.at(i_channel).size();
      for (unsigned int i_tdc = 0.; i_tdc < tdc_file.at(i_channel).size(); i_tdc++){
        rms_temp += pow((tdc_file.at(i_channel).at(i_tdc)-mean_tdc),2);
      }
      rms_temp=sqrt(rms_temp);
      rms_temp/=tdc_file.at(i_channel).size();
    } else {
      rms_temp = 0.;
    }

    for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
      if (crate_temp == loopback_crate.at(i_trigger) && slot_temp == loopback_slot.at(i_trigger) && channel_temp == loopback_channel.at(i_trigger)){
        if (loopback_name.at(i_trigger) == "Cosmic") n_cosmic+=tdc_file.at(i_channel).size();
        else if (loopback_name.at(i_trigger) == "Beam") n_beam+=tdc_file.at(i_channel).size();
      }
    }

    crate->push_back(crate_temp);
    slot->push_back(slot_temp);
    channel->push_back(channel_temp);
    tdc->push_back(mean_tdc);
    rms->push_back(rms_temp);
    rate->push_back(rate_temp);
    channelcount->push_back(tdc_file.at(i_channel).size());
    nevents+=tdc_file.at(i_channel).size();
  }

  if (fabs(t_frame) > 0.1) {
    rate_beam = n_beam/(t_frame/MSEC_to_SEC);
    rate_cosmic = n_cosmic/(t_frame/MSEC_to_SEC);
  } else {
    rate_beam = 0.;
    rate_cosmic = 0.;
  }

  t->Fill();
  t->Write("",TObject::kOverwrite);           //prevent ROOT from making endless keys for the same tree when updating the tree
  f->Close();

  delete crate;
  delete slot;
  delete channel;
  delete tdc;
  delete rms;
  delete rate;
  delete channelcount;
  delete f;     //tree should get deleted automatically by closing file

  gROOT->cd();

}

void MonitorMRDTime::ReadFromFile(ULong64_t timestamp_end, double time_frame){

  //-------------------------------------------------------
  //------------------ReadFromFile-------------------------
  //-------------------------------------------------------

  tdc_plot.clear();
  rms_plot.clear();
  rate_plot.clear();
  channelcount_plot.clear();
  tstart_plot.clear();
  tend_plot.clear();
  cosmicrate_plot.clear();
  beamrate_plot.clear();
  noloopbackrate_plot.clear();
  normalhitrate_plot.clear();
  doublehitrate_plot.clear();
  zerohitsrate_plot.clear();
  nevents_plot.clear();
  labels_timeaxis.clear();

  //take the end time and calculate the start time with the given time_frame
  ULong64_t timestamp_start = timestamp_end - time_frame*MIN_to_HOUR*SEC_to_MIN*MSEC_to_SEC;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(timestamp_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_start/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_start/MSEC_to_SEC/1000.)%60,timestamp_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);

  if (verbosity > 2) {
    std::cout <<"MonitorMRDTime: ReadFromFile: Reading in data for time frame "<<starttime_tm.tm_year+1900<<"/"<<starttime_tm.tm_mon+1<<"/"<<starttime_tm.tm_mday<<"-"<<starttime_tm.tm_hour<<":"<<starttime_tm.tm_min<<":"<<starttime_tm.tm_sec;
    std::cout <<" ... "<<endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec<<std::endl;
  }

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
    root_filename_i << path_monitoring << "MRD_" << string_date_i <<".root";
    bool tree_exists = true;

    if (does_file_exist(root_filename_i.str())) {
      TFile *f = new TFile(root_filename_i.str().c_str(),"READ");
      TTree *t;
      if (f->GetListOfKeys()->Contains("mrdmonitor_tree")) t = (TTree*) f->Get("mrdmonitor_tree");
      else { 
        if (verbosity > 1) std::cout <<"WARNING (MonitorMRDTime): File "<<root_filename_i.str()<<" does not contain mrdmonitor_tree. Omit file."<<std::endl;
        tree_exists = false;
      }

      if (tree_exists){

        if (verbosity > 2) std::cout <<"Tree exists, start reading in data"<<std::endl;

        ULong64_t t_start, t_end;
        std::vector<unsigned int> *crate = new std::vector<unsigned int>;
        std::vector<unsigned int> *slot = new std::vector<unsigned int>;
        std::vector<unsigned int> *channel = new std::vector<unsigned int>;
        std::vector<double> *tdc = new std::vector<double>;
        std::vector<double> *rms = new std::vector<double>;
        std::vector<double> *rate = new std::vector<double>;
        std::vector<int> *channelcount = new std::vector<int>;
        double rate_cosmic=0.;
        double rate_beam=0.;
        double rate_noloopback=0.;
        double rate_normalhit=0.;
        double rate_doublehit=0.;
        double rate_zerohits=0.;
        int nevents;
        int nentries_tree;

        t->SetBranchAddress("t_start",&t_start);
        t->SetBranchAddress("t_end",&t_end);
        t->SetBranchAddress("crate",&crate);
        t->SetBranchAddress("slot",&slot);
        t->SetBranchAddress("channel",&channel);
        t->SetBranchAddress("tdc",&tdc);
        t->SetBranchAddress("rms",&rms);
        t->SetBranchAddress("rate",&rate);
        t->SetBranchAddress("channelcount",&channelcount);
        t->SetBranchAddress("rate_cosmic",&rate_cosmic);
        t->SetBranchAddress("rate_beam",&rate_beam);
        t->SetBranchAddress("rate_noloopback",&rate_noloopback);
        t->SetBranchAddress("rate_normalhit",&rate_normalhit);
        t->SetBranchAddress("rate_doublehit",&rate_doublehit);
        t->SetBranchAddress("rate_zerohits",&rate_zerohits);
        t->SetBranchAddress("nevents",&nevents);

        nentries_tree = t->GetEntries();

        //Sort timestamps for the case that they are not in order
        //
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
            tdc_plot.push_back(*tdc);
            rms_plot.push_back(*rms);
            rate_plot.push_back(*rate);
            channelcount_plot.push_back(*channelcount);
            tstart_plot.push_back(t_start);
            tend_plot.push_back(t_end);
            cosmicrate_plot.push_back(rate_cosmic);
            beamrate_plot.push_back(rate_beam);
            noloopbackrate_plot.push_back(rate_noloopback);
            normalhitrate_plot.push_back(rate_normalhit);
            doublehitrate_plot.push_back(rate_doublehit);
            zerohitsrate_plot.push_back(rate_zerohits);
            nevents_plot.push_back(nevents);
            boost::posix_time::ptime boost_tend = *Epoch+boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
            struct tm label_timestamp = boost::posix_time::to_tm(boost_tend);
            TDatime datime_timestamp(1900+label_timestamp.tm_year,label_timestamp.tm_mon+1,label_timestamp.tm_mday,label_timestamp.tm_hour,label_timestamp.tm_min,label_timestamp.tm_sec);
            labels_timeaxis.push_back(datime_timestamp);
          }

        }

        delete crate;
        delete slot;
        delete channel;
        delete tdc;
        delete rms;
        delete rate;
        delete channelcount;

      }
      f->Close();
      delete f;
      gROOT->cd();

    } else {
      if (verbosity > 0) std::cout <<"MonitorMRDTime: ReadFromFile: File "<<root_filename_i.str()<<" does not exist. Omit file."<<std::endl;
    }

  }

  //set the readfromfile time variables to make sure data is not read twice for the same time window
  readfromfile_tend = timestamp_end;
  readfromfile_timeframe = time_frame;

}

std::string MonitorMRDTime::convertTimeStamp_to_Date(ULong64_t timestamp){

    //format of date is YYYY_MM-DD

    boost::posix_time::ptime convertedtime = *Epoch + boost::posix_time::time_duration(int(timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp/MSEC_to_SEC/1000.)%60,timestamp%1000);
    struct tm converted_tm = boost::posix_time::to_tm(convertedtime);
    std::stringstream ss_date;
    ss_date << converted_tm.tm_year+1900 << "_" << converted_tm.tm_mon+1 << "-" <<converted_tm.tm_mday;
    return ss_date.str();

}


bool MonitorMRDTime::does_file_exist(std::string filename){

  std::ifstream infile(filename.c_str());
  bool file_good = infile.good();
  infile.close();
  return file_good;

}

void MonitorMRDTime::InitializeVectors(){

  if (verbosity > 2) std::cout <<"MonitorMRDTime: Initialize vectors, histograms and canvases ..."<<std::endl;

  std::vector<int> empty_vec;
  std::vector<ULong64_t> empty_longvector;
  n_bins_scatter = 200;

  readfromfile_tend = 0;
  readfromfile_timeframe = 0.;

  inactive_crate1=0;
  inactive_crate2=0;

  //-------------------------------------------------------
  //-----------------Initialize canvases-------------------
  //-------------------------------------------------------

  gROOT->cd();
  canvas_ch_tdc = new TCanvas("canvas_ch_tdc","Channel TDC Canvas",900,600);
  canvas_ch_rms = new TCanvas("canvas_ch_rms","Channel RMS Canvas",900,600);
  canvas_ch_rate = new TCanvas("canvas_ch_rate","Channel Freq Canvas",900,600);
  canvas_ch_single = new TCanvas("canvas_ch_single","Channel Canvas Single",900,600);
  canvas_hitmap = new TCanvas("canvas_hitmap","Hitmap MRD",900,600);
  canvas_hitmap_slot = new TCanvas("canvas_hitmap_slot","Hitmap MRD Slot",900,600);
  canvas_logfile_mrd = new TCanvas("canvas_logfile_mrd","MRD File History",900,600);
  canvas_scatter = new TCanvas("canvas_scatter","Scatter Plot Canvas",900,600);
  canvas_scatter_single = new TCanvas("canvas_scatter_single","Scatter Plot Canvas Single",900,600);
  canvas_tdc = new TCanvas("canvas_tdc","TDC Canvas",900,600);
  canvas_rate_electronics = new TCanvas("canvas_rate_electronics","Rate Electronics Space",900,600);
  canvas_rate_physical = new TCanvas("canvas_rate_physical","Rate Physical Space",900,600);
  canvas_rate_physical_facc = new TCanvas("canvas_rate_physical_facc","Rate Physical Space FACC",900,600);
  canvas_trigger = new TCanvas("canvas_trigger","MRD Trigger Rates",900,600);
  canvas_trigger_time = new TCanvas("canvas_trigger_time","MRD Trigger Rates (Time)",900,600);
  canvas_pie = new TCanvas("canvas_pie","MRD Pie Chart Canvas",700,700);
  canvas_file_timestamp = new TCanvas("canvas_file_timestamp","Timestamp Last File",900,600);

  canvas_hitmap->SetGridy();
  canvas_hitmap->SetLogy();
  canvas_hitmap_slot->SetGridx();
  canvas_hitmap_slot->SetLogy();


  //-------------------------------------------------------
  //-----------------Initialize hitmap histograms----------
  //-------------------------------------------------------

  hist_hitmap_cr1 = new TH1F("hist_hitmap_fivemin_cr1","Hitmap Rack 7",num_channels*num_active_slots,0,num_channels*num_active_slots);
  hist_hitmap_cr2 = new TH1F("hist_hitmap_fivemin_cr2","Hitmap Rack 8",num_channels*num_active_slots,0,num_channels*num_active_slots);

  hist_hitmap_cr1->SetLineColor(8);
  hist_hitmap_cr2->SetLineColor(9);
  hist_hitmap_cr1->SetFillColor(8);
  hist_hitmap_cr2->SetFillColor(9);
  hist_hitmap_cr1->SetStats(0);
  hist_hitmap_cr2->SetStats(0);
  hist_hitmap_cr1->GetYaxis()->SetTitle("# of entries");
  hist_hitmap_cr2->GetYaxis()->SetTitle("# of entries");
  hist_hitmap_cr1->GetXaxis()->SetNdivisions(int(num_active_slots));
  hist_hitmap_cr2->GetXaxis()->SetNdivisions(int(num_active_slots));

  for (int i_label=0;i_label<int(num_active_slots);i_label++){
    if (i_label<num_active_slots_cr1){
      std::stringstream ss_slot;
      ss_slot<<(active_slots_cr1.at(i_label));
      std::string str_slot = "slot "+ss_slot.str();
      hist_hitmap_cr1->GetXaxis()->SetBinLabel(num_channels*(i_label+0.5),str_slot.c_str());
    }
    else {
      std::stringstream ss_slot;
      ss_slot<<(active_slots_cr2.at(i_label-num_active_slots_cr1));
      std::string str_slot = "slot "+ss_slot.str();
      hist_hitmap_cr1->GetXaxis()->SetBinLabel(num_channels*(i_label+0.5),str_slot.c_str());          
    }
  }
  hist_hitmap_cr1->LabelsOption("v");
  hist_hitmap_cr1->SetTickLength(0,"X");          //workaround to only have labels for every slot
  hist_hitmap_cr2->SetTickLength(0,"X");
  hist_hitmap_cr1->SetLineWidth(2);
  hist_hitmap_cr2->SetLineWidth(2);
  hist_hitmap_cr1->SetTitleSize(0.3,"t");

  label_cr1 = new TPaveText(0.02,0.93,0.17,0.98,"NDC");
  label_cr1->SetFillColor(0);
  label_cr1->SetTextColor(8);
  label_cr1->AddText("Rack 7");
  label_cr2 = new TPaveText(0.83,0.93,0.98,0.98,"NDC");
  label_cr2->SetFillColor(0);
  label_cr2->SetTextColor(9);
  label_cr2->AddText("Rack 8");
  separate_crates = new TLine(num_channels*num_active_slots_cr1,0.8,num_channels*num_active_slots_cr1,100);
  separate_crates->SetLineStyle(2);
  separate_crates->SetLineWidth(2);
  f1 = new TF1("f1","x",0,num_active_slots);       //workaround to only have labels for every slot

  //-------------------------------------------------------
  //------------Initialize single slot hitmaps-------------
  //-------------------------------------------------------

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){

    unsigned int crate = ActiveSlot_to_Crate[i_slot];
    unsigned int slot = ActiveSlot_to_Slot[i_slot];
    std::stringstream ss_title_single, ss_name_single;
    ss_name_single << "hist_hitmap_cr"<<crate<<"_sl"<<slot;
    ss_title_single << "Hitmap Crate "<<crate<<" Slot "<<slot;

    int hist_color;
    if (crate == min_crate) hist_color = 8;
    else hist_color = 9;

    TH1F *hist_hitmap_single = new TH1F(ss_name_single.str().c_str(),ss_title_single.str().c_str(),num_channels,0,num_channels);
    hist_hitmap_single->SetFillColor(hist_color);
    hist_hitmap_single->SetLineColor(hist_color);
    hist_hitmap_single->SetLineWidth(2);
    hist_hitmap_single->SetStats(0);
    hist_hitmap_single->GetYaxis()->SetTitle("# of entries");
    hist_hitmap_single->GetXaxis()->SetTitle("Channel #");
    hist_hitmap_slot.push_back(hist_hitmap_single);

  }

  //-------------------------------------------------------
  //------------Initialize scatter legends-----------------
  //-------------------------------------------------------

  leg_scatter = new TLegend(0.7,0.7,0.88,0.88);
  leg_scatter->SetNColumns(4);
  leg_scatter->SetLineColor(0);

  leg_scatter_trigger = new TLegend(0.7,0.7,0.85,0.78);
  leg_scatter_trigger->SetLineColor(0);

  //-------------------------------------------------------
  //------------Initialize file history hist---------------
  //-------------------------------------------------------

  num_files_history = 10;
  log_files_mrd = new TH1F("log_files_mrd","MRD Files History",num_files_history,0,num_files_history);
  log_files_mrd->GetXaxis()->SetTimeDisplay(1);
  log_files_mrd->GetXaxis()->SetLabelSize(0.03);
  log_files_mrd->GetXaxis()->SetLabelOffset(0.03);
  log_files_mrd->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  log_files_mrd->GetYaxis()->SetTickLength(0.);
  log_files_mrd->GetYaxis()->SetLabelOffset(999);
  log_files_mrd->SetStats(0);

  //-------------------------------------------------------
  //------------TDC histogram------------------------------
  //-------------------------------------------------------

  hist_tdc = new TH1F("hist_tdc","TDC (last file)",500,0,1000);
  hist_tdc->GetXaxis()->SetTitle("TDC");
  hist_tdc->GetYaxis()->SetTitle("#");
  hist_tdc_cluster = new TH1F("hist_tdc_cluster","TDC Cluster (last file)",20,0,1000);
  hist_tdc_cluster->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster->GetYaxis()->SetTitle("#");
  hist_tdc_cluster_20 = new TH1F("hist_tdc_cluster_20","TDC Cluster (last 20 files)",20,0,1000);
  hist_tdc_cluster_20->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster_20->GetYaxis()->SetTitle("#");
  hist_tdc_beam = new TH1F("hist_tdc_beam","TDC Beam (last file)",500,0,1000);
  hist_tdc_beam->GetXaxis()->SetTitle("TDC");
  hist_tdc_beam->GetYaxis()->SetTitle("#");
  hist_tdc_cluster_beam = new TH1F("hist_tdc_cluster_beam","TDC Cluster Beam (last file)",20,0,1000);
  hist_tdc_cluster_beam->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster_beam->GetYaxis()->SetTitle("#");
  hist_tdc_cluster_beam_20 = new TH1F("hist_tdc_cluster_beam_20","TDC Cluster Beam (last 20 files)",20,0,1000);
  hist_tdc_cluster_beam_20->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster_beam_20->GetYaxis()->SetTitle("#");
  hist_tdc_cosmic = new TH1F("hist_tdc_cosmic","TDC Cosmic (last file)",500,0,1000);
  hist_tdc_cosmic->GetXaxis()->SetTitle("TDC");
  hist_tdc_cosmic->GetYaxis()->SetTitle("#");
  hist_tdc_cluster_cosmic = new TH1F("hist_tdc_cluster_cosmic","TDC Cluster Cosmic (last file)",20,0,1000);
  hist_tdc_cluster_cosmic->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster_cosmic->GetYaxis()->SetTitle("#");
  hist_tdc_cluster_cosmic_20 = new TH1F("hist_tdc_cluster_cosmic_20","TDC Cluster Cosmic (last 20 files)",20,0,1000);
  hist_tdc_cluster_cosmic_20->GetXaxis()->SetTitle("TDC");
  hist_tdc_cluster_cosmic_20->GetYaxis()->SetTitle("#");

  //-------------------------------------------------------
  //------------Pie charts---------------------------------
  //-------------------------------------------------------

  pie_triggertype = new TPie("pie_triggertype","MRD Triggers",3);   //types are Cosmic, Beam, No loopback
  pie_weirdevents = new TPie("pie_weirdevents","MRD Event Types",3);  //types are Normal, Zero Hits, Multiple Clusters
  
  pie_triggertype->GetSlice(0)->SetFillColor(4);
  pie_triggertype->GetSlice(0)->SetTitle("beam");
  pie_triggertype->GetSlice(1)->SetFillColor(1);
  pie_triggertype->GetSlice(1)->SetTitle("cosmic");
  pie_triggertype->GetSlice(2)->SetFillColor(2);
  pie_triggertype->GetSlice(2)->SetTitle("no loopback");
  pie_triggertype->SetCircle(.5,.45,.3);
  pie_triggertype->GetSlice(1)->SetRadiusOffset(0.05);
  pie_triggertype->GetSlice(2)->SetRadiusOffset(0.1);
  pie_triggertype->SetLabelFormat("%val");
  pie_triggertype->SetLabelsOffset(0.02);
  leg_triggertype = (TLegend*) pie_triggertype->MakeLegend();
  leg_triggertype->SetY1(0.78);
  leg_triggertype->SetY2(0.91);

  pie_weirdevents->GetSlice(0)->SetFillColor(8);
  pie_weirdevents->GetSlice(0)->SetTitle("normal event");
  pie_weirdevents->GetSlice(1)->SetFillColor(2);
  pie_weirdevents->GetSlice(1)->SetTitle("0-size event");
  pie_weirdevents->GetSlice(2)->SetFillColor(9);
  pie_weirdevents->GetSlice(2)->SetTitle("multiple-hit event");
  pie_weirdevents->SetCircle(.5,.45,.3);
  pie_weirdevents->GetSlice(1)->SetRadiusOffset(0.05);
  pie_weirdevents->GetSlice(2)->SetRadiusOffset(0.1);
  pie_weirdevents->SetLabelFormat("%val");
  pie_weirdevents->SetLabelsOffset(0.02);
  leg_weirdevents = (TLegend*) pie_weirdevents->MakeLegend();
  leg_weirdevents->SetY1(0.78);
  leg_weirdevents->SetY2(0.91);

  //-------------------------------------------------------
  //------------Rate histograms----------------------------
  //-------------------------------------------------------

  rate_crate1 = new TH2F("rate_crate1","MRD Rates Rack 7",num_slots,0,num_slots,num_channels,0,num_channels);
  rate_crate2 = new TH2F("rate_crate2","MRD Rates Rack 8",num_slots,0,num_slots,num_channels,0,num_channels);

  rate_crate1->SetStats(0);
  rate_crate1->GetXaxis()->SetNdivisions(num_slots);
  rate_crate1->GetYaxis()->SetNdivisions(num_channels);
  for (int i_label=0;i_label<int(num_channels);i_label++){
    std::stringstream ss_slot, ss_ch;
    ss_slot<<(i_label+1);
    ss_ch<<(i_label);
    std::string str_ch = "ch "+ss_ch.str();
    rate_crate1->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
    if (i_label<num_slots){
      std::string str_slot = "slot "+ss_slot.str();
      rate_crate1->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
    }
  }
  rate_crate1->LabelsOption("v");
  rate_crate2->SetStats(0);
  rate_crate2->GetXaxis()->SetNdivisions(num_slots);
  rate_crate2->GetYaxis()->SetNdivisions(num_channels);
  for (int i_label=0;i_label<int(num_channels);i_label++){
    std::stringstream ss_slot, ss_ch;
    ss_slot<<(i_label+1);
    ss_ch<<(i_label);
    std::string str_ch = "ch "+ss_ch.str();
    rate_crate2->GetYaxis()->SetBinLabel(i_label+1,str_ch.c_str());
    if (i_label<num_slots){
      std::string str_slot = "slot "+ss_slot.str();
      rate_crate2->GetXaxis()->SetBinLabel(i_label+1,str_slot.c_str());
    }
  }
  rate_crate2->LabelsOption("v");

  label_rate_cr1 = new TLatex(0.88,0.93, "R [Hz]");
  label_rate_cr1->SetNDC();
  label_rate_cr1->SetTextSize(0.03);
  label_rate_cr2 = new TLatex(0.88,0.93, "R [Hz]");
  label_rate_cr2->SetNDC();
  label_rate_cr2->SetTextSize(0.03);
  label_rate_facc = new TLatex(0.88,0.93,"R [Hz]");
  label_rate_facc->SetNDC();
  label_rate_facc->SetTextSize(0.03);

  //set up histogram for rate plots in physical space

  rate_top = new TH2Poly("rate_top","MRD Rates Top View",1.6,3.,-2.,2.);
  rate_side = new TH2Poly("rate_side","MRD Rates Side View",1.6,3.,-2.,2.);
  rate_facc = new TH2Poly("rate_facc","FMV Rates",-1.66,-1.58,-2.5,2.5);
  rate_top->GetXaxis()->SetTitle("z [m]");
  rate_top->GetYaxis()->SetTitle("x [m]");
  rate_side->GetXaxis()->SetTitle("z [m]");
  rate_side->GetYaxis()->SetTitle("y [m]");
  rate_facc->GetXaxis()->SetTitle("z [m]");
  rate_facc->GetYaxis()->SetTitle("y [m]");
  rate_top->SetStats(0);
  rate_side->SetStats(0);
  rate_facc->SetStats(0);

  // Set custom bin shapes for the histograms

  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
                                                    it != Detectors->at("MRD").end();
                                                  ++it){
    Detector* amrdpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
    Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);

    double xmin = mrdpaddle->GetXmin();
    double xmax = mrdpaddle->GetXmax();
    double ymin = mrdpaddle->GetYmin();
    double ymax = mrdpaddle->GetYmax();
    double zmin = mrdpaddle->GetZmin();
    double zmax = mrdpaddle->GetZmax();
    int orientation = mrdpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
    int half = mrdpaddle->GetHalf();                  //0 or 1
    int side = mrdpaddle->GetSide();

    if (orientation == 0){
      //horizontal layers --> side representation
      if (half==0) rate_side->AddBin(zmin-tank_center_z-enlargeBoxes,ymin-tank_center_y,zmax-tank_center_z+enlargeBoxes,ymax-tank_center_y);
      else rate_side->AddBin(zmin+shiftSecRow-tank_center_z-enlargeBoxes,ymin-tank_center_y,zmax+shiftSecRow-tank_center_z+enlargeBoxes,ymax-tank_center_y);
    } else {
      //vertical layers --> top representation
      if (half==0) rate_top->AddBin(zmin-tank_center_z-enlargeBoxes,xmin-tank_center_x,zmax-tank_center_z+enlargeBoxes,xmax-tank_center_x);
      else  rate_top->AddBin(zmin+shiftSecRow-tank_center_z-enlargeBoxes,xmin-tank_center_x,zmax+shiftSecRow-tank_center_z+enlargeBoxes,xmax-tank_center_x);
    }
  }

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Veto").begin();
                                                  it != Detectors->at("Veto").end();
                                                ++it){
    Detector* afaccpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = afaccpmt->GetChannels()->begin()->first;
    Paddle *faccpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);

    double xmin = faccpaddle->GetXmin();
    double xmax = faccpaddle->GetXmax();
    double ymin = faccpaddle->GetYmin();
    double ymax = faccpaddle->GetYmax();
    double zmin = faccpaddle->GetZmin();
    double zmax = faccpaddle->GetZmax();
    int orientation = faccpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
    int half = faccpaddle->GetHalf();                  //0 or 1
    int side = faccpaddle->GetSide();

    //std::cout <<"facc chankey "<<chankey<<", xmin = "<<xmin<<", xmax = "<<xmax<<", ymin = "<<ymin<<", ymax = "<<ymax<<", zmin = "<<zmin-tank_center_z<<", zmax = "<<zmax-tank_center_z<<", half = "<<half<<", side = "<<side<<std::endl;

    rate_facc->AddBin(zmin-tank_center_z,ymin-tank_center_y,zmax-tank_center_z,ymax-tank_center_y);

  }

  //-------------------------------------------------------
  //----------Color inactive channels in grey--------------
  //-------------------------------------------------------

  std::vector<TBox*> empty_boxvector;
  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){
    vector_box_inactive_hitmap.emplace(i_slot,empty_boxvector);
  }

  inactive_crate1 = 0;
  inactive_crate2 = 0;
  for (int i_slot=0;i_slot<num_slots;i_slot++){
    if (active_channel[0][i_slot]==0){
      TBox *box_inactive = new TBox(i_slot,0,i_slot+1,num_channels);
      vector_box_inactive.push_back(box_inactive);
      box_inactive->SetFillStyle(3004);
      box_inactive->SetFillColor(1);
      inactive_crate1++;
    }
  }
  for (unsigned int i_ch = 0; i_ch < inactive_ch_crate1.size(); i_ch++){
    if (verbosity > 2) std::cout <<"inactive ch crate1, entry "<<i_ch<<std::endl;
    TBox *box_inactive = new TBox(inactive_slot_crate1.at(i_ch)-1,inactive_ch_crate1.at(i_ch),inactive_slot_crate1.at(i_ch),inactive_ch_crate1.at(i_ch)+1);
    vector_box_inactive.push_back(box_inactive);
    box_inactive->SetFillStyle(3004);
    box_inactive->SetFillColor(1);
    inactive_crate1++;
    TBox *box_inactive_hitmap = new TBox(inactive_ch_crate1.at(i_ch),0.8,inactive_ch_crate1.at(i_ch)+1,1.);
    box_inactive_hitmap->SetFillStyle(3004);
    box_inactive_hitmap->SetFillColor(1);
    std::vector<unsigned int> CrateSlot{min_crate,inactive_slot_crate1.at(i_ch)};
    int active_slot = CrateSlot_to_ActiveSlot[CrateSlot];
    vector_box_inactive_hitmap[active_slot].push_back(box_inactive_hitmap);
  }

  for (int i_slot=0;i_slot<num_slots;i_slot++){
    if (active_channel[1][i_slot]==0){
      TBox *box_inactive = new TBox(i_slot,0,i_slot+1,num_channels);
      vector_box_inactive.push_back(box_inactive);
      box_inactive->SetFillColor(1);
      box_inactive->SetFillStyle(3004);
      inactive_crate2++;
    }
  }
  for (unsigned int i_ch = 0; i_ch < inactive_ch_crate2.size(); i_ch++){
    if (verbosity > 2) std::cout <<"inactive ch crate2, entry "<<i_ch<<std::endl;
    TBox *box_inactive = new TBox(inactive_slot_crate2.at(i_ch)-1,inactive_ch_crate2.at(i_ch),inactive_slot_crate2.at(i_ch),inactive_ch_crate2.at(i_ch)+1);
    vector_box_inactive.push_back(box_inactive);
    box_inactive->SetFillStyle(3004);
    box_inactive->SetFillColor(1);
    inactive_crate2++;
    TBox *box_inactive_hitmap = new TBox(inactive_ch_crate2.at(i_ch),0.,inactive_ch_crate2.at(i_ch)+1,1.);
    box_inactive_hitmap->SetFillStyle(3004);
    box_inactive_hitmap->SetFillColor(1);
    std::vector<unsigned int> CrateSlot{min_crate+1,inactive_slot_crate2.at(i_ch)};
    int active_slot = CrateSlot_to_ActiveSlot[CrateSlot];
    vector_box_inactive_hitmap[active_slot].push_back(box_inactive_hitmap);
  }


  //-------------------------------------------------------
  //-------Initialize channel-wise graphs + hists ---------
  //-------------------------------------------------------


  for (int i_crate = 0; i_crate < num_crates; i_crate++){

    std::stringstream ss_crate;
    ss_crate << i_crate+min_crate;

    for (int i_slot = 0; i_slot < num_slots; i_slot++){

      std::stringstream ss_slot;
      ss_slot << i_slot+1;

      if (active_channel[i_crate][i_slot]==0) continue;           //don't plot values for slots that are not active

      for (int i_channel = 0; i_channel < num_channels; i_channel++){

        std::stringstream ss_ch;
        ss_ch << i_channel;
        int line_color_ch = color_scheme[i_channel%(num_channels/2)];

        tdc_file.push_back(empty_vec);
        timestamp_file.push_back(empty_longvector);

        //channel-wise graphs
        std::stringstream name_graph,name_graph_tdc, name_graph_rms, name_graph_rate, name_scatter;
        std::stringstream title_graph, title_graph_tdc, title_graph_rms, title_graph_rate, title_scatter;

        name_graph << "cr"<<ss_crate.str()<<"_sl"<<ss_slot.str()<<"_ch"<<ss_ch.str();
        name_graph_tdc << name_graph.str()<<"_tdc";
        name_graph_rms << name_graph.str()<<"_rms";
        name_graph_rate << name_graph.str()<<"_rate";
        name_scatter << name_graph.str()<<"_scatter";
        title_graph << "Crate "<<ss_crate.str()<<", Slot "<<ss_slot.str()<<", Channel "<<ss_ch.str();
        title_graph_tdc << title_graph.str() << " (TDC)";
        title_graph_rms << title_graph.str() <<" (RMS)";
        title_graph_rate << title_graph.str() <<" (Rate)";
        title_scatter << title_graph.str() << " (TDC Scatter)";


        TGraph *graph_ch_tdc = new TGraph();
        TGraph *graph_ch_rms = new TGraph();
        TGraph *graph_ch_rate = new TGraph();
        TH2F *hist_scatter_tdc = new TH2F(name_scatter.str().c_str(),title_scatter.str().c_str(),200,0,3600,200,0,200);

        graph_ch_tdc->SetName(name_graph_tdc.str().c_str());
        graph_ch_tdc->SetTitle(title_graph_tdc.str().c_str());
        graph_ch_rms->SetName(name_graph_rms.str().c_str());
        graph_ch_rms->SetTitle(title_graph_rms.str().c_str());
        graph_ch_rate->SetName(name_graph_rate.str().c_str());
        graph_ch_rate->SetTitle(title_graph_rate.str().c_str());
        hist_scatter_tdc->SetName(name_scatter.str().c_str());
        hist_scatter_tdc->SetTitle(title_scatter.str().c_str());
        
        if (draw_marker) {
          graph_ch_tdc->SetMarkerStyle(20);
          graph_ch_rms->SetMarkerStyle(20);
          graph_ch_rate->SetMarkerStyle(20);
        }
        graph_ch_tdc->SetMarkerColor(line_color_ch);
        graph_ch_rms->SetMarkerColor(line_color_ch);
        graph_ch_rate->SetMarkerColor(line_color_ch);
        graph_ch_tdc->SetLineColor(line_color_ch);
        graph_ch_tdc->SetLineWidth(2);
        graph_ch_tdc->SetFillColor(0);
        graph_ch_tdc->GetYaxis()->SetTitle("TDC");
        graph_ch_tdc->GetXaxis()->SetTimeDisplay(1);
        graph_ch_tdc->GetXaxis()->SetLabelSize(0.03);
        graph_ch_tdc->GetXaxis()->SetLabelOffset(0.03);
        graph_ch_tdc->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        graph_ch_rms->SetLineColor(line_color_ch);
        graph_ch_rms->SetLineWidth(2);
        graph_ch_rms->SetFillColor(0);
        graph_ch_rms->GetYaxis()->SetTitle("RMS (TDC)");
        graph_ch_rms->GetXaxis()->SetTimeDisplay(1);
        graph_ch_rms->GetXaxis()->SetLabelSize(0.03);
        graph_ch_rms->GetXaxis()->SetLabelOffset(0.03);
        graph_ch_rms->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        graph_ch_rate->SetLineColor(line_color_ch);
        graph_ch_rate->SetLineWidth(2);
        graph_ch_rate->SetFillColor(0);
        graph_ch_rate->GetYaxis()->SetTitle("Rate [Hz]");
        graph_ch_rate->GetXaxis()->SetTimeDisplay(1);
        graph_ch_rate->GetXaxis()->SetLabelSize(0.03);
        graph_ch_rate->GetXaxis()->SetLabelOffset(0.03);
        graph_ch_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        hist_scatter_tdc->SetMarkerColor(line_color_ch);
        hist_scatter_tdc->SetFillColor(0);
        hist_scatter_tdc->GetYaxis()->SetTitle("TDC");
        hist_scatter_tdc->SetStats(0);
        hist_scatter_tdc->SetMarkerStyle(20);
        hist_scatter_tdc->SetMarkerSize(0.3);
        hist_scatter_tdc->GetXaxis()->SetLabelSize(0.03);
        hist_scatter_tdc->GetXaxis()->SetLabelOffset(0.03);
        hist_scatter_tdc->GetXaxis()->SetTimeDisplay(1);
        hist_scatter_tdc->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");

        gr_tdc.push_back(graph_ch_tdc);
        gr_rms.push_back(graph_ch_rms);
        gr_rate.push_back(graph_ch_rate);
        hist_scatter.push_back(hist_scatter_tdc);
      }
    }
  }

  multi_ch_tdc = new TMultiGraph();
  multi_ch_rms = new TMultiGraph();
  multi_ch_rate = new TMultiGraph();
  leg_tdc = new TLegend(0.7,0.7,0.88,0.88);
  leg_rms = new TLegend(0.7,0.7,0.88,0.88);
  leg_rate = new TLegend(0.7,0.7,0.88,0.88);
  leg_tdc->SetNColumns(4);
  leg_rms->SetNColumns(4);
  leg_rate->SetNColumns(4);
  leg_tdc->SetLineColor(0);
  leg_rms->SetLineColor(0);
  leg_rate->SetLineColor(0);

  //define graphs for trigger time evolution plots
  leg_trigger = new TLegend(0.7,0.7,0.82,0.8);
  leg_noloopback = new TLegend(0.7,0.7,0.88,0.76);
  leg_eventtypes = new TLegend(0.7,0.7,0.88,0.8);
  for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
    std::stringstream name_graph_trigger, title_graph_trigger;
    name_graph_trigger << "trigger_"<<loopback_name.at(i_trigger);
    title_graph_trigger << "Rate "<<loopback_name.at(i_trigger)<<" trigger";
    TGraph *graph_trigger = new TGraph();
    graph_trigger->SetName(name_graph_trigger.str().c_str());
    graph_trigger->SetTitle(title_graph_trigger.str().c_str());
    int line_color;
    if (i_trigger == 0) line_color = 4;
    else if (i_trigger == 1) line_color = 1;
    else line_color = 2;
    if (draw_marker){
      graph_trigger->SetMarkerStyle(20); 
    }
    graph_trigger->SetMarkerColor(line_color);
    graph_trigger->SetLineColor(line_color);
    graph_trigger->SetLineWidth(2);
    graph_trigger->SetFillColor(0);
    graph_trigger->GetYaxis()->SetTitle("Rate [Hz]");
    graph_trigger->GetXaxis()->SetTimeDisplay(1);
    graph_trigger->GetXaxis()->SetLabelSize(0.03);
    graph_trigger->GetXaxis()->SetLabelOffset(0.03);
    graph_trigger->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    leg_trigger->AddEntry(graph_trigger,loopback_name.at(i_trigger).c_str(),"l");
    gr_trigger.push_back(graph_trigger);
  }

  gr_noloopback = new TGraph();
  gr_noloopback->SetName("trigger_noloopback");
  gr_noloopback->SetTitle("Rate No Loopback");
  gr_noloopback->SetLineColor(2);
  if (draw_marker){
    gr_noloopback->SetMarkerStyle(20);
  }
  gr_noloopback->SetMarkerColor(2);
  gr_noloopback->SetLineWidth(2);
  gr_noloopback->SetFillColor(0);
  gr_noloopback->GetYaxis()->SetTitle("Rate [Hz]");
  gr_noloopback->GetXaxis()->SetTimeDisplay(1);
  gr_noloopback->GetXaxis()->SetLabelSize(0.03);
  gr_noloopback->GetXaxis()->SetLabelOffset(0.03);
  gr_noloopback->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  leg_noloopback->AddEntry(gr_noloopback,"no loopback","l");

  gr_zerohits = new TGraph();
  gr_zerohits->SetName("trigger_zerohits");
  gr_zerohits->SetTitle("Rate 0-size events");
  gr_zerohits->SetLineColor(2);
  if (draw_marker){
    gr_zerohits->SetMarkerStyle(20);
  }
  gr_zerohits->SetMarkerColor(2);
  gr_zerohits->SetLineWidth(2);
  gr_zerohits->SetFillColor(0);
  gr_zerohits->GetYaxis()->SetTitle("Rate [Hz]");
  leg_eventtypes->AddEntry(gr_zerohits,"0-hit events","l");

  gr_doublehits = new TGraph();
  gr_doublehits->SetName("trigger_doublehits");
  gr_doublehits->SetTitle("Rate multiple-hit events");
  gr_doublehits->SetLineColor(9);
  if (draw_marker){
    gr_doublehits->SetMarkerStyle(20);
  }
  gr_doublehits->SetMarkerColor(9);
  gr_doublehits->SetLineWidth(2);
  gr_doublehits->SetFillColor(0);
  gr_doublehits->GetYaxis()->SetTitle("Rate [Hz]");
  leg_eventtypes->AddEntry(gr_doublehits,"multiple-hit events","l");

  multi_trigger = new TMultiGraph();
  multi_eventtypes = new TMultiGraph();

}


void MonitorMRDTime::DrawFileHistory(ULong64_t timestamp_end, double time_frame, std::string file_ending, int _linewidth){

  if (verbosity > 2) std::cout <<"MonitorMRDTime: Drawing File History plot"<<std::endl;

  //-------------------------------------------------------
  //------------------DrawFileHistory----------------------
  //-------------------------------------------------------

  //Creates a plot showing the time stamps for all the files within the last time_frame mins
  //The plot is updated with the update_frequency specified in the configuration file (default: 5 mins)

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  timestamp_end += utc_to_t;
  ULong64_t timestamp_start = timestamp_end - time_frame*MSEC_to_SEC*SEC_to_MIN*MIN_to_HOUR;

  canvas_logfile_mrd->cd();
  log_files_mrd->SetBins(num_files_history,timestamp_start/MSEC_to_SEC,timestamp_end/MSEC_to_SEC);
  log_files_mrd->GetXaxis()->SetTimeOffset(0.);
  log_files_mrd->Draw();

  std::vector<TLine*> file_markers;
  for (unsigned int i_file = 0; i_file < tend_plot.size(); i_file++){
   if ((tend_plot.at(i_file)+utc_to_t)>=timestamp_start && (tend_plot.at(i_file)+utc_to_t)<=timestamp_end){
      TLine *line_file = new TLine((tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,0.,(tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,1.);
      line_file->SetLineColor(1);
      line_file->SetLineStyle(1);
      line_file->SetLineWidth(_linewidth);
      line_file->Draw("same"); 
      file_markers.push_back(line_file);
    }
  }

  std::stringstream ss_logfiles;
  ss_logfiles << outpath << "MRD_FileHistory_" << file_ending << "." << img_extension;
  canvas_logfile_mrd->SaveAs(ss_logfiles.str().c_str());

  for (unsigned int i_line = 0; i_line < file_markers.size(); i_line++){
    delete file_markers.at(i_line);
  }

  log_files_mrd->Reset();
  canvas_logfile_mrd->Clear();

}

void MonitorMRDTime::PrintFileTimeStamp(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  if (verbosity > 2) std::cout <<"MonitorMRDTime: PrintFileTimeStamp"<<std::endl;

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
    time_diff << "#Delta t Last MRD File: >"<<time_frame<<"h";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  } else {
    ULong64_t timestamp_lastfile = tend_plot.at(tend_plot.size()-1);
    boost::posix_time::ptime filetime = *Epoch + boost::posix_time::time_duration(int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_lastfile/MSEC_to_SEC/1000.)%60,timestamp_lastfile%1000);
    boost::posix_time::time_duration t_since_file= boost::posix_time::time_duration(endtime - filetime);
    int t_since_file_min = int(t_since_file.total_milliseconds()/MSEC_to_SEC/SEC_to_MIN);
    std::stringstream time_diff;
    time_diff << "#Delta t Last MRD File: "<<t_since_file_min/60<<"h:"<<t_since_file_min%60<<"min";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  }

  canvas_file_timestamp->cd();
  label_lastfile->Draw();
  label_timediff->Draw();
  std::stringstream ss_file_timestamp;
  ss_file_timestamp << outpath << "MRD_FileTimeStamp_" << file_ending << "." << img_extension;
  canvas_file_timestamp->SaveAs(ss_file_timestamp.str().c_str());
  
  delete label_lastfile;
  delete label_timediff;

  canvas_file_timestamp->Clear();
}


void MonitorMRDTime::DrawLastFilePlots(){

  //-------------------------------------------------------
  //------------------DrawLastFilePlots--------------------
  //-------------------------------------------------------

  //Plots drawn using the information from the last file
  //Some of the plots are only shown here since they require a lot of detailed information that is not saved later on (e.g. the scatter plots)
  //This includes scatter TDC plots for all channels and more time-resolved time evolution plots for the trigger rate time evolution (beam, cosmic)
  //It also includes a plot of the multiple-hit-per-channel rate time evolution

  //Draw scatter plots
  //DrawScatterPlots();
  DrawScatterPlotsTrigger();

  //Draw hitmap plots
  DrawHitMap(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw TDC histogram plot
  DrawTDCHistogram();

  //Draw rate plots in 2D (complementary to hitmap plots), both in electronics and in physical space
  DrawRatePlotElectronics(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");
  DrawRatePlotPhysical(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw pie charts showing the event/trigger type distribution
  DrawPieChart(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

}

void MonitorMRDTime::UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes){

  //-------------------------------------------------------
  //--------------UpdateMonitorPlots-----------------------
  //-------------------------------------------------------

  //Draw the monitoring plots according to the specifications in the configfiles

  for (unsigned int i_time = 0; i_time < timeFrames.size(); i_time++){

    ULong64_t zero = 0;
    if (endTimes.at(i_time) == zero) endTimes.at(i_time) = t_file_end;        //set 0 for t_file_end since we did not know what that was at the beginning of initialise
    /*std::cout << (endTimes.at(i_time) == zero) << std::endl;
    std::cout << (endTimes.at(i_time) == 0) << std::endl;
    std::cout <<t_file_end<<std::endl;*/

    for (unsigned int i_plot = 0; i_plot < plotTypes.at(i_time).size(); i_plot++){

      if (plotTypes.at(i_time).at(i_plot) == "Hitmap") DrawHitMap(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "RateElectronics") DrawRatePlotElectronics(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "RatePhysical") DrawRatePlotPhysical(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "TimeEvolution") DrawTimeEvolution(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "PieChartTrigger") DrawPieChart(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "TriggerEvolution") DrawTriggerEvolution(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "FileHistory") DrawFileHistory(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time),1);
      else {
        if (verbosity > 0) std::cout <<"ERROR (MonitorMRDTime): UpdateMonitorPlots: Specified plot type -"<<plotTypes.at(i_time).at(i_plot)<<"- does not exist! Omit entry."<<std::endl;
      }
    }
  }

}


void MonitorMRDTime::DrawScatterPlots(){

  //-------------------------------------------------------
  //--------------DrawScatterPlots-------------------------
  //-------------------------------------------------------

  //Since the scatter plots are only shown for the last file, there is no read in stage for the data
  //Data is simply taken from the vectors already storing the information from the current file

  for (int i_channel=0; i_channel<num_active_slots*num_channels;i_channel++){
    if (hist_scatter.at(i_channel)->GetEntries()>0) hist_scatter.at(i_channel)->Reset();
    hist_scatter.at(i_channel)->SetBins(n_bins_scatter,0,t_file_end/MSEC_to_SEC-t_file_start/MSEC_to_SEC,n_bins_scatter,0,200);
    for (unsigned int i_entry=0; i_entry<tdc_file.at(i_channel).size(); i_entry++){
      hist_scatter.at(i_channel)->Fill((timestamp_file.at(i_channel).at(i_entry)-t_file_start)/MSEC_to_SEC,tdc_file.at(i_channel).at(i_entry));
    }
  }
  
  double max_canvas_scatter = 0.;
  double min_canvas_scatter = 999999999.;
  int CH_per_CANVAS = 16;   //channels per canvas
  int CANVAS_NR=0;

  for (int i_channel=0; i_channel<num_active_slots*num_channels;i_channel++){

    std::stringstream ss_ch_scatter, ss_leg_scatter;
    unsigned int crate = TotalChannel_to_Crate[i_channel];
    unsigned int slot = TotalChannel_to_Slot[i_channel];
    unsigned int channel = TotalChannel_to_Channel[i_channel];
    ss_leg_scatter.str("");
    ss_leg_scatter<<"ch "<<channel;

    if ((i_channel%CH_per_CANVAS == 0 || i_channel == num_active_slots*num_channels-1) && i_channel != 0){

      //save the canvases every CH_per_CANVAS channels

      std::string channel_range;
      if (CANVAS_NR == 0) channel_range = " Ch 0-15";
      else {
	//we already switched to the new slot with channel 32, so we need to consider the slot for the channel before
        slot = TotalChannel_to_Slot[i_channel-1];
	channel_range = " Ch 16-31";
      }

      ss_ch_scatter << "Crate "<<crate<<" Slot "<<slot<<channel_range<<" (last File)";

      if ( i_channel == num_active_slots*num_channels-1){
        canvas_scatter->cd();
        hist_scatter.at(i_channel)->Draw("same");
        leg_scatter->AddEntry(hist_scatter.at(i_channel),ss_leg_scatter.str().c_str(),"p");
        hist_scatter.at(i_channel-CH_per_CANVAS+1)->SetTitle(ss_ch_scatter.str().c_str());
        leg_scatter->Draw();
      } else {
        canvas_scatter->cd();
        hist_scatter.at(i_channel-CH_per_CANVAS)->SetTitle(ss_ch_scatter.str().c_str());
        leg_scatter->Draw();
      }

      std::string channel_range_name;
      if (CANVAS_NR == 0) channel_range_name = "Ch0-15";
      else channel_range_name = "Ch16-31";
      ss_ch_scatter.str("");
      ss_ch_scatter<<outpath<<"MRDScatter_lastFile_Cr"<<crate<<"_Sl"<<slot<<"_"<<channel_range_name<<"."<<img_extension;
      canvas_scatter->SaveAs(ss_ch_scatter.str().c_str());
      CANVAS_NR=(CANVAS_NR+1)%2;

      canvas_scatter->Clear();
      //if (i_channel != num_active_slots*num_channels-1) leg_scatter->Clear();
      leg_scatter->Clear();
    }

    //add histograms to the canvas
    canvas_scatter->cd();
    if (i_channel%CH_per_CANVAS == 0) {
      hist_scatter.at(i_channel)->GetXaxis()->SetTimeOffset(t_file_start+utc_to_t/MSEC_to_SEC);
      hist_scatter.at(i_channel)->Draw();
    } else {
      hist_scatter.at(i_channel)->Draw("same");
    }

    //add legend entry for histograms
    if (hist_scatter.at(i_channel)->GetMaximum()>max_canvas_scatter) max_canvas_scatter = hist_scatter.at(i_channel)->GetMaximum();
    if (hist_scatter.at(i_channel)->GetMinimum()<min_canvas_scatter) min_canvas_scatter = hist_scatter.at(i_channel)->GetMinimum();
    leg_scatter->AddEntry(hist_scatter.at(i_channel),ss_leg_scatter.str().c_str(),"p");

    if (draw_single){
      //save single channel scatter plots as well
      canvas_scatter_single->Clear();
      canvas_scatter_single->cd();
      hist_scatter.at(i_channel)->GetXaxis()->SetTimeOffset(t_file_start+utc_to_t/MSEC_to_SEC);
      hist_scatter.at(i_channel)->Draw();
      std::stringstream ss_ch_scatter_single;
      ss_ch_scatter_single<<outpath<<"MRDScatter_lastFile_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"."<<img_extension;
      canvas_scatter_single->SaveAs(ss_ch_scatter_single.str().c_str());
      }

  }
}

void MonitorMRDTime::DrawScatterPlotsTrigger(){

  //-------------------------------------------------------
  //--------------DrawScatterPlotsTrigger------------------
  //-------------------------------------------------------

  //Since the scatter plots are only shown for the last file, there is no read in stage for the data
  //Data is simply taken from the vectors already storing the information from the current file
  //This version of the function only plots the TDC scatter plot of the trigger loopback channel

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){

    unsigned int crate = loopback_crate.at(i_trigger);
    unsigned int slot = loopback_slot.at(i_trigger);
    unsigned int channel = loopback_channel.at(i_trigger);
    std::vector<unsigned int> CrateSlotChannel{crate,slot,channel};
    int total_ch = CrateSlotChannel_to_TotalChannel[CrateSlotChannel];
    if (hist_scatter.at(total_ch)->GetEntries()>0) hist_scatter.at(total_ch)->Reset();
    hist_scatter.at(total_ch)->SetBins(n_bins_scatter,0,t_file_end/MSEC_to_SEC-t_file_start/MSEC_to_SEC,n_bins_scatter,0,1000);    //show whole TDC acq window from 0 ... 1000 TDC units
    for (unsigned int i_entry=0; i_entry<tdc_file.at(total_ch).size(); i_entry++){
      hist_scatter.at(total_ch)->Fill((timestamp_file.at(total_ch).at(i_entry)-t_file_start)/MSEC_to_SEC,tdc_file.at(total_ch).at(i_entry));
    }
  }
  
  canvas_scatter->cd();
  canvas_scatter->Clear();

  for (unsigned int i_trigger=0; i_trigger < loopback_name.size(); i_trigger++){

    std::stringstream ss_ch_scatter;
    unsigned int crate = loopback_crate.at(i_trigger);
    unsigned int slot = loopback_slot.at(i_trigger);
    unsigned int channel = loopback_channel.at(i_trigger);
    std::vector<unsigned int> CrateSlotChannel{crate,slot,channel};
    int total_ch = CrateSlotChannel_to_TotalChannel[CrateSlotChannel];

    if (i_trigger == 0){
      ss_ch_scatter << "Trigger TDC "<<end_time.str()<<" (last File)";
      hist_scatter.at(total_ch)->GetXaxis()->SetTimeOffset((t_file_start+utc_to_t)/MSEC_to_SEC);
      hist_scatter.at(total_ch)->SetTitle(ss_ch_scatter.str().c_str());
      hist_scatter.at(total_ch)->Draw();
    }
    else hist_scatter.at(total_ch)->Draw("same");
    leg_scatter_trigger->AddEntry(hist_scatter.at(total_ch),loopback_name.at(i_trigger).c_str(),"lp");

  }
  leg_scatter_trigger->Draw();
  std::stringstream ss_scatter_trigger;
  ss_scatter_trigger<<outpath<<"MRDScatter_Triggers_lastFile."<<img_extension;
  canvas_scatter->SaveAs(ss_scatter_trigger.str().c_str());

  leg_scatter_trigger->Clear();

}

void MonitorMRDTime::DrawTDCHistogram(){

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::vector<int> coinc_times_mrd;	//Introduce a vector to store found coincidence times
  std::vector<int> coinc_times_mrd_beam;	//Introduce a vector to store found coincidence times (Beam events)
  std::vector<int> coinc_times_mrd_cosmic;	//Introduce a vector to store found coincidence times (Cosmic events)

  canvas_tdc->cd();
  canvas_tdc->Clear();
  canvas_tdc->SetLogy();

  for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){
      if (TotalChannel_to_Crate[i_channel] == loopback_crate.at(0) && TotalChannel_to_Slot[i_channel] == loopback_slot.at(0) && TotalChannel_to_Channel[i_channel] == loopback_channel.at(0)) continue; //Omit cosmic loopback signal
      if (TotalChannel_to_Crate[i_channel] == loopback_crate.at(1) && TotalChannel_to_Slot[i_channel] == loopback_slot.at(1) && TotalChannel_to_Channel[i_channel] == loopback_channel.at(1)) continue; //Omit beam loopback signal
    for (unsigned int i_tdc = 0; i_tdc < tdc_file.at(i_channel).size(); i_tdc++){
      hist_tdc->Fill(tdc_file.at(i_channel).at(i_tdc));
    }
  }

  std::vector<bool> is_beam;
  std::vector<bool> is_cosmic;  

  for (unsigned int i_entry=0; i_entry < tdc_file_times.size(); i_entry++){
    std::sort(tdc_file_times.at(i_entry).begin(),tdc_file_times.at(i_entry).end());
    bool cosmic = false;
    bool beam = false;
    for (int i_t=0; i_t < (int) tdc_file_times.at(i_entry).size(); i_t++){
      int channel = tdc_file_ch.at(i_entry).at(i_t);
      if (TotalChannel_to_Crate[channel] == loopback_crate.at(1) && TotalChannel_to_Slot[channel] == loopback_slot.at(1) && TotalChannel_to_Channel[channel] == loopback_channel.at(1)) beam = true;
      if (TotalChannel_to_Crate[channel] == loopback_crate.at(0) && TotalChannel_to_Slot[channel] == loopback_slot.at(0) && TotalChannel_to_Channel[channel] == loopback_channel.at(0)) cosmic = true;
    }
    if (beam && cosmic) {
      is_cosmic.push_back(1);
      is_beam.push_back(0);
    }
    else if (beam) {
      is_cosmic.push_back(0);
      is_beam.push_back(1);
    } else {
      is_cosmic.push_back(0);
      is_beam.push_back(0);
    }
  }

  int tdc_start=0;
  int min_cluster = 4;
  int min_tdc_separation = 8;
  int n_channels = 0;
  int tdc_current = 0;

  for (unsigned int i_entry = 0; i_entry < tdc_file_times.size(); i_entry++){
   std::vector<int> tdc_times = tdc_file_times.at(i_entry);
   bool beam = is_beam.at(i_entry);
   bool cosmic = is_cosmic.at(i_entry);
   for (unsigned int i_tdc=0; i_tdc < tdc_times.size(); i_tdc++){
     if (beam) hist_tdc_beam->Fill(tdc_times.at(i_tdc));
     else if (cosmic) hist_tdc_cosmic->Fill(tdc_times.at(i_tdc));
     if (i_tdc == 0) {
	tdc_start = tdc_times.at(0);
        tdc_current = tdc_times.at(0);
	n_channels=1;
        continue;
      } 
    if ((tdc_times.at(i_tdc)-tdc_current) <= min_tdc_separation){
        n_channels++;
	tdc_current = tdc_times.at(i_tdc);
	if (i_tdc==tdc_times.size()-1) {
	  if (n_channels >= min_cluster){
	    for (int i_ch = 0; i_ch < n_channels; i_ch++){
	      hist_tdc_cluster->Fill(tdc_times.at(i_tdc-i_ch));
              coinc_times_mrd.push_back(tdc_times.at(i_tdc-i_ch));
              if (cosmic){
	        hist_tdc_cluster_cosmic->Fill(tdc_times.at(i_tdc-i_ch));
                coinc_times_mrd_cosmic.push_back(tdc_times.at(i_tdc-i_ch));
              } else if (beam){
	        hist_tdc_cluster_beam->Fill(tdc_times.at(i_tdc-i_ch));
                coinc_times_mrd_beam.push_back(tdc_times.at(i_tdc-i_ch));
              }
            }
          }
	}
      } else {
        if (n_channels >= min_cluster){
          for (int i_ch = 1; i_ch <= n_channels; i_ch++){
            hist_tdc_cluster->Fill(tdc_times.at(i_tdc-i_ch));
            coinc_times_mrd.push_back(tdc_times.at(i_tdc-i_ch));
            if (cosmic){
              hist_tdc_cluster_cosmic->Fill(tdc_times.at(i_tdc-i_ch));
              coinc_times_mrd_cosmic.push_back(tdc_times.at(i_tdc-i_ch));
            } else if (beam){
              hist_tdc_cluster_beam->Fill(tdc_times.at(i_tdc-i_ch));
              coinc_times_mrd_beam.push_back(tdc_times.at(i_tdc-i_ch));
            }
	  }
        }
        n_channels = 1;
        tdc_current = tdc_times.at(i_tdc);
      }
    }  
  }
  
  std::stringstream ss_tdc_hist_title;
  ss_tdc_hist_title << "TDC "<<end_time.str()<<" (last file) ";
  hist_tdc->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc->Draw();
  std::stringstream ss_tdc_hist;
  ss_tdc_hist << outpath << "MRDTDCHist_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster "<<end_time.str()<<" (last file) ";
  hist_tdc_cluster->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_Cluster_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Beam "<<end_time.str()<<" (last file) ";
  hist_tdc_beam->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_beam->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHistBeam_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_beam->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster Beam "<<end_time.str()<<" (last file) ";
  hist_tdc_cluster_beam->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster_beam->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_ClusterBeam_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster_beam->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cosmic "<<end_time.str()<<" (last file) ";
  hist_tdc_cosmic->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cosmic->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHistCosmic_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cosmic->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster Cosmic "<<end_time.str()<<" (last file) ";
  hist_tdc_cluster_cosmic->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster_cosmic->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_ClusterCosmic_lastFile."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster_cosmic->Reset();
  canvas_tdc->Clear();
  
  overall_mrd_coinc_times.push_back(coinc_times_mrd);
  overall_mrd_coinc_times_beam.push_back(coinc_times_mrd_beam);
  overall_mrd_coinc_times_cosmic.push_back(coinc_times_mrd_cosmic);
  while (overall_mrd_coinc_times.size() > 20){
    overall_mrd_coinc_times.erase(overall_mrd_coinc_times.begin());
  }
  while (overall_mrd_coinc_times_beam.size() > 20){
    overall_mrd_coinc_times_beam.erase(overall_mrd_coinc_times_beam.begin());
  }
  while (overall_mrd_coinc_times_cosmic.size() > 20){
    overall_mrd_coinc_times_cosmic.erase(overall_mrd_coinc_times_cosmic.begin());
  }
  for (int i_coinc = 0; i_coinc < (int) overall_mrd_coinc_times.size(); i_coinc++){
    std::vector<int> single_coinc_times = overall_mrd_coinc_times.at(i_coinc);
    std::vector<int> single_coinc_times_beam = overall_mrd_coinc_times_beam.at(i_coinc);
    std::vector<int> single_coinc_times_cosmic = overall_mrd_coinc_times_cosmic.at(i_coinc);
    for (int i_hit=0; i_hit < (int) single_coinc_times.size(); i_hit++){
      hist_tdc_cluster_20->Fill(single_coinc_times.at(i_hit));
    }
    for (int i_hit=0; i_hit < (int) single_coinc_times_beam.size(); i_hit++){
      hist_tdc_cluster_beam_20->Fill(single_coinc_times_beam.at(i_hit));
    }
    for (int i_hit=0; i_hit < (int) single_coinc_times_cosmic.size(); i_hit++){
      hist_tdc_cluster_cosmic_20->Fill(single_coinc_times_cosmic.at(i_hit));
    }
  }
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster "<<end_time.str()<<" (last 20 files) ";
  hist_tdc_cluster_20->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster_20->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_Cluster_last20Files."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster_20->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster Beam "<<end_time.str()<<" (last 20 files) ";
  hist_tdc_cluster_beam_20->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster_beam_20->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_ClusterBeam_last20Files."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster_beam_20->Reset();
  canvas_tdc->Clear();
  ss_tdc_hist_title.str("");
  ss_tdc_hist_title << "TDC Cluster Cosmic "<<end_time.str()<<" (last 20 files) ";
  hist_tdc_cluster_cosmic_20->SetTitle(ss_tdc_hist_title.str().c_str());
  hist_tdc_cluster_cosmic_20->Draw();
  ss_tdc_hist.str("");
  ss_tdc_hist << outpath << "MRDTDCHist_ClusterCosmic_last20Files."<<img_extension;
  canvas_tdc->SaveAs(ss_tdc_hist.str().c_str());
  hist_tdc_cluster_cosmic_20->Reset();
  canvas_tdc->Clear();

}

  
void MonitorMRDTime::DrawHitMap(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------------DrawHitMap---------------------------
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

  //fill the hitmap histograms with the counts for each channel
  for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){

    long sum_channel=0;

    for (unsigned int i_file=0; i_file < tdc_plot.size(); i_file++){
      sum_channel+=channelcount_plot.at(i_file).at(i_channel);
    }

    unsigned int crate = TotalChannel_to_Crate[i_channel];
    unsigned int slot = TotalChannel_to_Slot[i_channel];
    unsigned int channel = TotalChannel_to_Channel[i_channel];
    std::vector<unsigned int> crate_slot{crate,slot};
    int active_slot = CrateSlot_to_ActiveSlot[crate_slot];

    if (crate == min_crate){

      hist_hitmap_cr2->SetBinContent(i_channel+1, 0.001);
      if (sum_channel!=0){

        if (max_hitmap < sum_channel) max_hitmap = sum_channel;
        if (max_hitmap_slot.at(active_slot) < sum_channel)  max_hitmap_slot.at(active_slot) = sum_channel;
        hist_hitmap_cr1->SetBinContent(i_channel+1,sum_channel);
        hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1,sum_channel);
      } else {
        hist_hitmap_cr1->SetBinContent(i_channel+1, 0.001);
        hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1, 0.001);
      }

    } else {

      hist_hitmap_cr1->SetBinContent(i_channel+1, 0.001);
      if (sum_channel!=0){
        if (max_hitmap < sum_channel) max_hitmap = sum_channel;
        if (max_hitmap_slot.at(active_slot) < sum_channel)  max_hitmap_slot.at(active_slot) = sum_channel;
        hist_hitmap_cr2->SetBinContent(i_channel+1,sum_channel);
        hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1,sum_channel);
      } else {
        hist_hitmap_cr2->SetBinContent(i_channel+1, 0.001);
        hist_hitmap_slot.at(active_slot)->SetBinContent(channel+1, 0.001);
      }
    }
  }

  //Plot the overall hitmap histogram

  canvas_hitmap->cd();
  hist_hitmap_cr1->GetYaxis()->SetRangeUser(0.8,max_hitmap+10);
  hist_hitmap_cr2->GetYaxis()->SetRangeUser(0.8,max_hitmap+10);
  ss_hitmap_title<<"Hitmap "<<end_time.str()<<" (last "<<ss_timeframe.str()<<"h)";
  hist_hitmap_cr1->SetTitle(ss_hitmap_title.str().c_str());
  hist_hitmap_cr1->Draw();
  hist_hitmap_cr2->Draw("same");
  separate_crates->SetY2(max_hitmap+10);     //adjust the range of the hitmap histogram
  separate_crates->SetLineStyle(2);
  separate_crates->SetLineWidth(2);
  separate_crates->Draw("same");
  label_cr1->Draw();
  label_cr2->Draw();
  canvas_hitmap->Update();
  TGaxis *labels_grid = new TGaxis(0,canvas_hitmap->GetUymin(),num_active_slots*num_channels,canvas_hitmap->GetUymin(),"f1",num_active_slots,"w");
  labels_grid->SetLabelSize(0);
  labels_grid->Draw("w");
  std::stringstream save_path_hitmap;
  save_path_hitmap << outpath <<"MRDHitmap_"<<file_ending<<"."<<img_extension;
  canvas_hitmap->SaveAs(save_path_hitmap.str().c_str());

  delete labels_grid;

  //Plot the more detailed slot-wise hitmap distributions

  for (int i_slot = 0; i_slot < num_active_slots; i_slot++){

    unsigned int crate = ActiveSlot_to_Crate[i_slot];
    unsigned int slot = ActiveSlot_to_Slot[i_slot];

    canvas_hitmap_slot->cd();
    std::stringstream ss_slot;
    ss_slot << i_slot;
    std::stringstream save_path_singlehitmap, ss_hitmap_title_slot;
    save_path_singlehitmap << outpath <<"MRDHitmap_Cr"<<crate<<"_Sl"<<slot<<"_"<<file_ending<<"."<<img_extension;
    ss_hitmap_title_slot << "Hitmap "<<end_time.str()<<" Cr "<<crate <<" Sl "<<slot<<" (last "<<ss_timeframe.str()<<"h)";
    hist_hitmap_slot.at(i_slot)->GetYaxis()->SetRangeUser(0.8,max_hitmap_slot.at(i_slot)+10);
    hist_hitmap_slot.at(i_slot)->SetTitle(ss_hitmap_title_slot.str().c_str());
    hist_hitmap_slot.at(i_slot)->Draw();

    for (unsigned int i_box = 0; i_box < vector_box_inactive_hitmap[i_slot].size(); i_box++){
      vector_box_inactive_hitmap[i_slot].at(i_box)->SetY2(max_hitmap_slot.at(i_slot)+10);
      vector_box_inactive_hitmap[i_slot].at(i_box)->Draw("same");
    }
    canvas_hitmap_slot->Update();
    canvas_hitmap_slot->SaveAs(save_path_singlehitmap.str().c_str());
    canvas_hitmap_slot->Clear();

  }

}

void MonitorMRDTime::DrawRatePlotElectronics(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------DrawRatePlotElectronics--------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::vector<double> overall_rates;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_rates.assign(num_active_slots*num_channels,0);
  double max_ch = 0.;
  double min_ch = 999999999.;

  for (unsigned int i_file = 0; i_file < rate_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

    for (int i_ch = 0; i_ch < num_active_slots*num_channels; i_ch++){
      overall_rates.at(i_ch) += (rate_plot.at(i_file).at(i_ch)*current_timeframe);
    }
  }

  if (overall_timeframe > 0.){
    for (int i_ch = 0.; i_ch < num_active_slots*num_channels; i_ch++){
      overall_rates.at(i_ch)/=overall_timeframe;
    }
  }

  for (int i_ch = 0; i_ch < num_active_slots*num_channels; i_ch++){
    
    unsigned int crate = TotalChannel_to_Crate[i_ch];
    unsigned int slot = TotalChannel_to_Slot[i_ch];
    unsigned int channel = TotalChannel_to_Channel[i_ch];
    if (crate == min_crate) rate_crate1->SetBinContent(slot,channel+1,overall_rates.at(i_ch));
    else if (crate == min_crate + 1) rate_crate2->SetBinContent(slot,channel+1,overall_rates.at(i_ch));
    if (overall_rates.at(i_ch) > max_ch) max_ch = overall_rates.at(i_ch); 
    if (overall_rates.at(i_ch) < min_ch) min_ch = overall_rates.at(i_ch); 
  }

  canvas_rate_electronics->cd();
  canvas_rate_electronics->Divide(2,1);

  TPad *p1 = (TPad*) canvas_rate_electronics->cd(1);
  p1->SetGrid();
  std::stringstream ss_rack7;
  ss_rack7 << "Rates Rack 7 "<<end_time.str()<<" ("<<file_ending<<")";
  rate_crate1->SetTitle(ss_rack7.str().c_str());
  rate_crate1->Draw("colz");
  for (int i_ch = 0; i_ch < inactive_crate1; i_ch++){
    vector_box_inactive.at(i_ch)->Draw("same");
  } 
  p1->Update();
  if (rate_crate1->GetMaximum()>0.){
    if (min_ch == max_ch) rate_crate1->GetZaxis()->SetRangeUser(min_ch-0.5,max_ch+0.5);
    else rate_crate1->GetZaxis()->SetRangeUser(1e-6,1.);
    TPaletteAxis *palette = 
    (TPaletteAxis*) rate_crate1->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }

  label_rate_cr1->Draw();
  std::stringstream ss_rack8;
  ss_rack8 << "Rates Rack 8 "<<end_time.str()<<" ("<<file_ending<<")";
  rate_crate2->SetTitle(ss_rack8.str().c_str());
  TPad *p2 = (TPad*) canvas_rate_electronics->cd(2);
  p2->SetGrid();
  rate_crate2->Draw("colz");
  for (int i_ch = 0; i_ch < inactive_crate2; i_ch++){
    vector_box_inactive.at(inactive_crate1+i_ch)->Draw("same");
  }
  p2->Update();

  if (rate_crate2->GetMaximum()>0.){
    if (min_ch == max_ch) rate_crate2->GetZaxis()->SetRangeUser(min_ch-0.5,max_ch+0.5);
    else rate_crate2->GetZaxis()->SetRangeUser(1e-6,1.);
    TPaletteAxis *palette = 
    (TPaletteAxis*)rate_crate2->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
  }
  label_rate_cr2->Draw();

  //p2->Modified();
  std::stringstream ss_rate_electronics;
  ss_rate_electronics<<outpath<<"MRDRates_Electronics_"<<file_ending<<"."<<img_extension;
  canvas_rate_electronics->SaveAs(ss_rate_electronics.str().c_str());

  rate_crate1->Reset();
  rate_crate2->Reset();
  canvas_rate_electronics->Clear();

}

void MonitorMRDTime::DrawRatePlotPhysical(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------DrawRatePlotPhysical-----------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::vector<double> overall_rates;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_rates.assign(num_active_slots*num_channels,0);
  double max_ch = 0.;
  double min_ch = 999999999.;

  for (unsigned int i_file = 0; i_file < rate_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

    for (int i_ch = 0; i_ch < num_active_slots*num_channels; i_ch++){
      overall_rates.at(i_ch) += (rate_plot.at(i_file).at(i_ch)*current_timeframe);
    }
  }

  if (overall_timeframe > 0.){
    for (int i_ch = 0.; i_ch < num_active_slots*num_channels; i_ch++){
      overall_rates.at(i_ch)/=overall_timeframe;
    }
  }


  std::stringstream ss_topTitle, ss_sideTitle, ss_faccTitle;
  ss_topTitle << "MRD Rates - Top "<<end_time.str()<<" ("<<file_ending<<")";
  ss_sideTitle << "MRD Rates - Side "<<end_time.str()<<" ("<<file_ending<<")";
  ss_faccTitle << "FMV Rates "<<end_time.str()<<" ("<<file_ending<<")";
  rate_side->SetTitle(ss_sideTitle.str().c_str());
  rate_top->SetTitle(ss_topTitle.str().c_str());
  rate_facc->SetTitle(ss_faccTitle.str().c_str());

  //fill the rate plot in detector space (Event Display - like)

  for (int i_ch = 0; i_ch < num_active_slots*num_channels; i_ch++){
  
    unsigned int crate_id = TotalChannel_to_Crate[i_ch];
    unsigned int slot_id = TotalChannel_to_Slot[i_ch];
    unsigned int channel_id = TotalChannel_to_Channel[i_ch];

    if (overall_rates.at(i_ch) > max_ch) max_ch = overall_rates.at(i_ch);
    if (overall_rates.at(i_ch) < min_ch) min_ch = overall_rates.at(i_ch);

    //fill EventDisplay plots
    std::vector<unsigned int> crate_slot_channel{crate_id,slot_id,channel_id};
    std::vector<int> crate_slot_channel_int(crate_slot_channel.begin(),crate_slot_channel.end());
    if (CrateSpaceToChannelNumMap->find(crate_slot_channel_int)!=CrateSpaceToChannelNumMap->end()){
      int chankey = CrateSpaceToChannelNumMap->at(crate_slot_channel_int);

      Detector *det = (Detector*) geom->GetDetector(chankey);
      int detkey = det->GetDetectorID();
      Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
      
      int orientation = mrdpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
      int half = mrdpaddle->GetHalf();                  //0 or 1
      Position paddle_pos = mrdpaddle->GetOrigin();
      double x_value = paddle_pos.X()-tank_center_x;
      double y_value = paddle_pos.Y()-tank_center_y;
      double z_value = paddle_pos.Z()-tank_center_z;

      std::string detector_element = det->GetDetectorElement();

      if (detector_element == "MRD"){
        if (half == 1) z_value+=shiftSecRow;

        if (orientation == 0) {
          rate_side->Fill(z_value,y_value,overall_rates.at(i_ch));
        }
        else{
          rate_top->Fill(z_value,x_value,overall_rates.at(i_ch));
        } 
      } else if (detector_element == "Veto"){
          rate_facc->Fill(z_value,y_value,overall_rates.at(i_ch));
      }
    }
  }

  canvas_rate_physical->Divide(2,1);
  canvas_rate_physical->cd(1);
  //rate_side->GetZaxis()->SetRangeUser(1e-6,max_ch);
  rate_side->GetZaxis()->SetRangeUser(1e-6,1.);
  rate_side->GetZaxis()->SetTitleOffset(1.3);
  rate_side->GetZaxis()->SetTitleSize(0.03);
  rate_side->Draw("colz L");                           //option L to show contours around bins (indicating where MRD paddles are)
  canvas_rate_physical->Update();
  TPaletteAxis *palette = 
  (TPaletteAxis*) rate_side->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.9);
  palette->SetX2NDC(0.92);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  label_rate_cr1->Draw();

  canvas_rate_physical->cd(2);
  //rate_top->GetZaxis()->SetRangeUser(1e-6,max_ch);
  rate_top->GetZaxis()->SetRangeUser(1e-6,1.);
  rate_top->GetZaxis()->SetTitleOffset(1.3);
  rate_top->GetZaxis()->SetTitleSize(0.03);
  rate_top->Draw("colz L");
  canvas_rate_physical->Update();
  TPaletteAxis *palette2 = 
  (TPaletteAxis*) rate_top->GetListOfFunctions()->FindObject("palette");
  palette2->SetX1NDC(0.9);
  palette2->SetX2NDC(0.92);
  palette2->SetY1NDC(0.1);
  palette2->SetY2NDC(0.9);
  label_rate_cr2->Draw();


  std::stringstream ss_ratephysical;
  ss_ratephysical<<outpath<<"MRDRates_Detector_"<<file_ending<<"."<<img_extension;
  canvas_rate_physical->SaveAs(ss_ratephysical.str().c_str());

  canvas_rate_physical_facc->cd();
  //rate_facc->GetZaxis()->SetRangeUser(1e-6,max_ch);
  rate_facc->GetZaxis()->SetRangeUser(1e-6,0.4);
  rate_facc->GetZaxis()->SetTitleOffset(1.3);
  rate_facc->GetZaxis()->SetTitleSize(0.03);
  rate_facc->Draw("colz L");                           //option L to show contours around bins (indicating where MRD paddles are)
  canvas_rate_physical_facc->Update();
  TPaletteAxis *palette_facc = 
  (TPaletteAxis*) rate_facc->GetListOfFunctions()->FindObject("palette");
  palette_facc->SetX1NDC(0.9);
  palette_facc->SetX2NDC(0.92);
  palette_facc->SetY1NDC(0.1);
  palette_facc->SetY2NDC(0.9);
  label_rate_facc->Draw();

  std::stringstream ss_ratephysical_facc;
  ss_ratephysical_facc<<outpath<<"MRDRates_DetectorFACC_"<<file_ending<<"."<<img_extension;
  canvas_rate_physical_facc->SaveAs(ss_ratephysical_facc.str().c_str());

  rate_top->Reset("M");
  rate_side->Reset("M");
  rate_facc->Reset("M");
  canvas_rate_physical->Clear();
  canvas_rate_physical_facc->Clear();

}


void MonitorMRDTime::DrawTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------DrawTimeEvolution--------------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  //looping over all files that are in the time interval, each file will be one data point

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){
    gr_tdc.at(i_channel)->Set(0);
    gr_rms.at(i_channel)->Set(0);
    gr_rate.at(i_channel)->Set(0);
  }

  for (unsigned int i_file=0; i_file<tdc_plot.size(); i_file++){

    //Updating channel graphs

    if (verbosity > 2) std::cout <<"MonitorMRDTime: Stored data (file #"<<i_file+1<<"): "<<std::endl;
    for (int i_channel = 0; i_channel < num_active_slots*num_channels; i_channel++){
      gr_tdc.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),tdc_plot.at(i_file).at(i_channel));
      gr_rms.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),rms_plot.at(i_file).at(i_channel));
      gr_rate.at(i_channel)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),rate_plot.at(i_file).at(i_channel));
    }

  }

  // Drawing time evolution plots

  double max_canvas = 0;
  double min_canvas = 9999999.;
  double max_canvas_rms = 0;
  double min_canvas_rms = 99999999.;
  double max_canvas_rate = 0;
  double min_canvas_rate = 999999999.;

  int CH_per_CANVAS = 16;   //channels per canvas
  int CANVAS_NR=0;

  for (int i_channel = 0; i_channel<num_active_slots*num_channels; i_channel++){

    unsigned int crate = TotalChannel_to_Crate[i_channel];
    unsigned int slot = TotalChannel_to_Slot[i_channel];
    unsigned int channel = TotalChannel_to_Channel[i_channel];

    std::stringstream ss_ch_tdc, ss_ch_rms, ss_ch_rate, ss_leg_time;
  
     if (i_channel%CH_per_CANVAS == 0 || i_channel == num_active_slots*num_channels-1) {
      if (i_channel != 0){

        ss_ch_tdc.str("");
        ss_ch_rms.str("");
        ss_ch_rate.str("");
        std::string channel_range;
        if (CANVAS_NR == 0) channel_range = " Channel 0-15";
        else {
        //we already switched to the new slot with channel 32, so we need to consider the slot for the channel before
          crate = TotalChannel_to_Crate[i_channel-1];
          slot = TotalChannel_to_Slot[i_channel-1];
          channel_range = " Channel 16-31";
        }
        ss_ch_tdc<<"Crate "<<crate<<" Slot "<<slot<<channel_range<<" ("<<ss_timeframe.str()<<"h)";
        ss_ch_rms<<"Crate "<<crate<<" Slot "<<slot<<channel_range<<" ("<<ss_timeframe.str()<<"h)";
        ss_ch_rate<<"Crate "<<crate<<" Slot "<<slot<<channel_range<<" ("<<ss_timeframe.str()<<"h)";

        if ( i_channel == num_active_slots*num_channels - 1){
          ss_leg_time.str("");
          ss_leg_time<<"ch "<<channel;
          multi_ch_tdc->Add(gr_tdc.at(i_channel));
          leg_tdc->AddEntry(gr_tdc.at(i_channel),ss_leg_time.str().c_str(),"l");
          multi_ch_rms->Add(gr_rms.at(i_channel));
          leg_rms->AddEntry(gr_rms.at(i_channel),ss_leg_time.str().c_str(),"l");
          multi_ch_rate->Add(gr_rate.at(i_channel));
          leg_rate->AddEntry(gr_rate.at(i_channel),ss_leg_time.str().c_str(),"l");
        }

        canvas_ch_tdc->cd();
        multi_ch_tdc->Draw("apl");
        multi_ch_tdc->SetTitle(ss_ch_tdc.str().c_str());
        multi_ch_tdc->GetYaxis()->SetTitle("TDC");
        multi_ch_tdc->GetXaxis()->SetTimeDisplay(1);
        multi_ch_tdc->GetXaxis()->SetLabelSize(0.03);
        multi_ch_tdc->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_tdc->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_tdc->GetXaxis()->SetTimeOffset(0.);
        leg_tdc->Draw();
        canvas_ch_rms->cd();
        multi_ch_rms->Draw("apl");
        multi_ch_rms->SetTitle(ss_ch_rms.str().c_str());
        multi_ch_rms->GetYaxis()->SetTitle("RMS (TDC)");
        multi_ch_rms->GetYaxis()->SetTitleSize(0.035);
        multi_ch_rms->GetYaxis()->SetTitleOffset(1.3);
        multi_ch_rms->GetXaxis()->SetTimeDisplay(1);
        multi_ch_rms->GetXaxis()->SetLabelSize(0.03);
        multi_ch_rms->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_rms->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_rms->GetXaxis()->SetTimeOffset(0.);
        leg_rms->Draw();
        canvas_ch_rate->cd();
        multi_ch_rate->Draw("apl");
        multi_ch_rate->SetTitle(ss_ch_rate.str().c_str());
        multi_ch_rate->GetYaxis()->SetTitle("Rate [Hz]");
        multi_ch_rate->GetYaxis()->SetTitleSize(0.035);
        multi_ch_rate->GetYaxis()->SetTitleOffset(1.3);
        multi_ch_rate->GetXaxis()->SetTimeDisplay(1);
        multi_ch_rate->GetXaxis()->SetLabelSize(0.03);
        multi_ch_rate->GetXaxis()->SetLabelOffset(0.03);
        multi_ch_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_ch_rate->GetXaxis()->SetTimeOffset(0.);
        leg_rate->Draw();

        std::string channel_range_name;
        if (CANVAS_NR == 0) channel_range_name = "Ch0-15";
        else channel_range_name = "Ch16-31";
        ss_ch_tdc.str("");
        ss_ch_rms.str("");
        ss_ch_rate.str("");
        ss_ch_tdc<<outpath<<"MRDTimeEvolutionTDC_Cr"<<crate<<"_Sl"<<slot<<"_"<<channel_range_name<<"_"<<file_ending<<"."<<img_extension;
        ss_ch_rms<<outpath<<"MRDTimeEvolutionRMS_Cr"<<crate<<"_Sl"<<slot<<"_"<<channel_range_name<<"_"<<file_ending<<"."<<img_extension;
        ss_ch_rate<<outpath<<"MRDTimeEvolutionRate_Cr"<<crate<<"_Sl"<<slot<<"_"<<channel_range_name<<"_"<<file_ending<<"."<<img_extension;

        canvas_ch_tdc->SaveAs(ss_ch_tdc.str().c_str());
        canvas_ch_rms->SaveAs(ss_ch_rms.str().c_str());
        canvas_ch_rate->SaveAs(ss_ch_rate.str().c_str()); 

        CANVAS_NR=(CANVAS_NR+1)%2;

        for (int i_gr=0; i_gr < CH_per_CANVAS; i_gr++){
          int i_balance = (i_channel == num_active_slots*num_channels-1)? 1 : 0;
          multi_ch_tdc->RecursiveRemove(gr_tdc.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
          multi_ch_rms->RecursiveRemove(gr_rms.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
          multi_ch_rate->RecursiveRemove(gr_rate.at(i_channel-CH_per_CANVAS+i_gr+i_balance));
        }
      }

      leg_tdc->Clear();
      leg_rms->Clear();
      leg_rate->Clear();
      
      canvas_ch_tdc->Clear();
      canvas_ch_rms->Clear();
      canvas_ch_rate->Clear();

     } 

    if (i_channel != num_active_slots*num_channels-1){

      ss_leg_time.str("");
      ss_leg_time<<"ch "<<i_channel%(2*CH_per_CANVAS);
      multi_ch_tdc->Add(gr_tdc.at(i_channel));
      if (gr_tdc.at(i_channel)->GetMaximum()>max_canvas) max_canvas = gr_tdc.at(i_channel)->GetMaximum();
      if (gr_tdc.at(i_channel)->GetMinimum()<min_canvas) min_canvas = gr_tdc.at(i_channel)->GetMinimum();
      leg_tdc->AddEntry(gr_tdc.at(i_channel),ss_leg_time.str().c_str(),"l");
      multi_ch_rms->Add(gr_rms.at(i_channel));
      leg_rms->AddEntry(gr_rms.at(i_channel),ss_leg_time.str().c_str(),"l");
      if (gr_rms.at(i_channel)->GetMaximum()>max_canvas_rms) max_canvas_rms = gr_rms.at(i_channel)->GetMaximum();
      if (gr_rms.at(i_channel)->GetMinimum()<min_canvas_rms) min_canvas_rms = gr_rms.at(i_channel)->GetMinimum();
      multi_ch_rate->Add(gr_rate.at(i_channel));
      leg_rate->AddEntry(gr_rate.at(i_channel),ss_leg_time.str().c_str(),"l");
      if (gr_rate.at(i_channel)->GetMaximum()>max_canvas_rate) max_canvas_rate = gr_rate.at(i_channel)->GetMaximum();
      if (gr_rate.at(i_channel)->GetMinimum()<min_canvas_rate) min_canvas_rate = gr_rate.at(i_channel)->GetMinimum();

    }

    //single channel time evolution plots

    if (draw_single){

      canvas_ch_single->cd();
      canvas_ch_single->Clear();
      gr_tdc.at(i_channel)->GetYaxis()->SetTitle("TDC");
      gr_tdc.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_tdc.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_tdc.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_tdc.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_tdc.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_tdc.at(i_channel)->Draw("apl");
      std::stringstream ss_ch_single;
      ss_ch_single<<outpath<<"MRDTimeEvolutionTDC_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single->SaveAs(ss_ch_single.str().c_str());

      canvas_ch_single->Clear();
      gr_rms.at(i_channel)->GetYaxis()->SetTitle("RMS (TDC)");
      gr_rms.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_rms.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_rms.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_rms.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_rms.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_rms.at(i_channel)->Draw("apl");
      ss_ch_single.str("");
      ss_ch_single<<outpath<<"MRDTimeEvolutionRMS_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single->SaveAs(ss_ch_single.str().c_str());

      canvas_ch_single->Clear();
      gr_rate.at(i_channel)->GetYaxis()->SetTitle("Rate [Hz]");
      gr_rate.at(i_channel)->GetXaxis()->SetTimeDisplay(1);
      gr_rate.at(i_channel)->GetXaxis()->SetLabelSize(0.03);
      gr_rate.at(i_channel)->GetXaxis()->SetLabelOffset(0.03);
      gr_rate.at(i_channel)->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gr_rate.at(i_channel)->GetXaxis()->SetTimeOffset(0.);
      gr_rate.at(i_channel)->Draw("apl");
      ss_ch_single.str("");
      ss_ch_single<<outpath<<"MRDTimeEvolutionRate_Cr"<<crate<<"_Sl"<<slot<<"_Ch"<<channel<<"_"<<file_ending<<"."<<img_extension;
      canvas_ch_single->SaveAs(ss_ch_single.str().c_str());
    }

  }

} 

void MonitorMRDTime::DrawTriggerEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------DrawTriggerEvolution-----------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  std::stringstream ss_trigger, ss_filename_trigger, ss_noloopback, ss_filename_noloopback, ss_eventtypes, ss_filename_eventtypes;
  ss_trigger << "Trigger rates (last "<<ss_timeframe.str()<<" h)";
  ss_filename_trigger << outpath << "MRDTriggertypes_timeevolution_"<<file_ending<<"."<<img_extension;
  ss_noloopback << "No Loopback Trigger Rate (last "<<ss_timeframe.str()<<" h)";
  ss_filename_noloopback << outpath << "MRDTrigger_noloopback_"<<file_ending<<"."<<img_extension;
  ss_eventtypes << "Event Types Trigger Rate (last "<<ss_timeframe.str()<<" h)";
  ss_filename_eventtypes << outpath << "MRDTrigger_eventtypes_"<<file_ending<<"."<<img_extension;

  for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
    gr_trigger.at(i_trigger)->Set(0);
  }

  for (unsigned int i_file=0; i_file<cosmicrate_plot.size(); i_file++){

    for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
      double rate_temp;
      if (loopback_name.at(i_trigger) == "Cosmic") rate_temp = cosmicrate_plot.at(i_file);
      else if (loopback_name.at(i_trigger) == "Beam") rate_temp = beamrate_plot.at(i_file);
      else rate_temp = noloopbackrate_plot.at(i_file);
      gr_trigger.at(i_trigger)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),rate_temp);
    }
  }

  for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
    multi_trigger->Add(gr_trigger.at(i_trigger));
  }

  canvas_trigger_time->cd();
  canvas_trigger_time->Clear();
        
  multi_trigger->Draw("apl");
  multi_trigger->SetTitle(ss_trigger.str().c_str());
  multi_trigger->GetYaxis()->SetTitle("Rate [Hz]");
  multi_trigger->GetXaxis()->SetTimeDisplay(1);
  multi_trigger->GetXaxis()->SetLabelSize(0.03);
  multi_trigger->GetXaxis()->SetLabelOffset(0.03);
  multi_trigger->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  multi_trigger->GetXaxis()->SetTimeOffset(0.);
  leg_trigger->Draw();
  canvas_trigger_time->SaveAs(ss_filename_trigger.str().c_str());

  for (unsigned int i_trigger = 0; i_trigger < loopback_name.size(); i_trigger++){
    multi_trigger->RecursiveRemove(gr_trigger.at(i_trigger));
  }

  canvas_trigger_time->Clear();

  gr_noloopback->Set(0);
  gr_zerohits->Set(0);
  gr_doublehits->Set(0);

  for (unsigned int i_file=0; i_file<noloopbackrate_plot.size(); i_file++){
    gr_noloopback->SetPoint(i_file,labels_timeaxis[i_file].Convert(),noloopbackrate_plot.at(i_file));
    if (zerohitsrate_plot.at(i_file)<0.0000001) gr_zerohits->SetPoint(i_file,labels_timeaxis[i_file].Convert(),0.);
    else gr_zerohits->SetPoint(i_file,labels_timeaxis[i_file].Convert(),zerohitsrate_plot.at(i_file));
    if (doublehitrate_plot.at(i_file)<0.0000001) gr_doublehits->SetPoint(i_file,labels_timeaxis[i_file].Convert(),0.);  
    else gr_doublehits->SetPoint(i_file,labels_timeaxis[i_file].Convert(),doublehitrate_plot.at(i_file));
  }
  gr_noloopback->SetTitle(ss_noloopback.str().c_str());
  gr_noloopback->Draw("apl");
  gr_noloopback->GetYaxis()->SetTitle("Rate [Hz]");
  gr_noloopback->GetXaxis()->SetTimeDisplay(1);
  gr_noloopback->GetXaxis()->SetLabelSize(0.03);
  gr_noloopback->GetXaxis()->SetLabelOffset(0.03);
  gr_noloopback->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  gr_noloopback->GetXaxis()->SetTimeOffset(0.);
  leg_noloopback->Draw();
  canvas_trigger_time->SaveAs(ss_filename_noloopback.str().c_str());

  canvas_trigger_time->Clear();
  multi_eventtypes->Add(gr_zerohits);
  multi_eventtypes->Add(gr_doublehits);
  multi_eventtypes->Draw("apl");
  multi_eventtypes->SetTitle(ss_eventtypes.str().c_str());
  multi_eventtypes->GetYaxis()->SetTitle("Rate [Hz]");
  multi_eventtypes->GetXaxis()->SetTimeDisplay(1);
  multi_eventtypes->GetXaxis()->SetLabelSize(0.03);
  multi_eventtypes->GetXaxis()->SetLabelOffset(0.03);
  multi_eventtypes->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  multi_eventtypes->GetXaxis()->SetTimeOffset(0.);
  leg_eventtypes->Draw();
  canvas_trigger_time->SaveAs(ss_filename_eventtypes.str().c_str());

  multi_eventtypes->RecursiveRemove(gr_zerohits);
  multi_eventtypes->RecursiveRemove(gr_doublehits);

}

void MonitorMRDTime::DrawPieChart(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  //-------------------------------------------------------
  //------------DrawPieChartTrigger------------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  int nevents_cosmic = 0;
  int nevents_beam = 0;
  int nevents_noloopback = 0;
  int nevents_normal= 0;
  int nevents_zerohits = 0;
  int nevents_doublehits = 0;

  for (unsigned int i_file = 0; i_file < rate_plot.size(); i_file++){

    ULong64_t current_timeframe = (tend_plot.at(i_file)-tstart_plot.at(i_file))/MSEC_to_SEC;
    nevents_cosmic += cosmicrate_plot.at(i_file)*current_timeframe;
    nevents_beam += beamrate_plot.at(i_file)*current_timeframe;
    nevents_noloopback += noloopbackrate_plot.at(i_file)*current_timeframe;
    nevents_normal += normalhitrate_plot.at(i_file)*current_timeframe;
    nevents_zerohits += zerohitsrate_plot.at(i_file)*current_timeframe;
    nevents_doublehits += doublehitrate_plot.at(i_file)*current_timeframe;

  }

  canvas_pie->cd();
  pie_triggertype->GetSlice(0)->SetValue(nevents_beam);
  pie_triggertype->GetSlice(1)->SetValue(nevents_cosmic);
  pie_triggertype->GetSlice(2)->SetValue(nevents_noloopback);
  std::stringstream ss_pie, ss_pie_title;
  ss_pie_title << "MRD Triggers "<<end_time.str()<<" (last "<<round(time_frame*100)/100.<<"h)";
  ss_pie << outpath <<"MRDTriggertypes_"<<file_ending<<"."<<img_extension;
  pie_triggertype->SetTitle(ss_pie_title.str().c_str());
  pie_triggertype->Draw("tsc");
  leg_triggertype->Draw();
  canvas_pie->SaveAs(ss_pie.str().c_str());
  canvas_pie->Clear();

  pie_weirdevents->GetSlice(0)->SetValue(nevents_normal);
  pie_weirdevents->GetSlice(1)->SetValue(nevents_zerohits);
  pie_weirdevents->GetSlice(2)->SetValue(nevents_doublehits);
  std::stringstream ss_pie2_title;
  ss_pie2_title << "MRD Event types "<<end_time.str()<<" (last "<<round(time_frame*100)/100.<<"h)";
  pie_weirdevents->SetTitle(ss_pie2_title.str().c_str());
  pie_weirdevents->Draw("tsc");
  leg_weirdevents->Draw();
  ss_pie.str("");
  ss_pie << outpath <<"MRDEventtypes_"<<file_ending<<"."<<img_extension;
  canvas_pie->SaveAs(ss_pie.str().c_str());
  canvas_pie->Clear();

}
