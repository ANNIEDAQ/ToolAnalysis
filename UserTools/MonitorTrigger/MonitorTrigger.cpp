#include "MonitorTrigger.h"

MonitorTrigger::MonitorTrigger():Tool(){}


bool MonitorTrigger::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  //-------------------------------------------------------
  //---------------Initialise config file------------------
  //-------------------------------------------------------

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  //Only for debugging memory leaks (ROOT-related), otherwise comment out
  //std::cout <<"List of Objects (Beginning of Initialise): "<<std::endl;
  //gObjectTable->Print();

  //-------------------------------------------------------
  //-----------------Get Configuration---------------------
  //-------------------------------------------------------

  m_variables.Get("OutputPath",outpath_temp);
  m_variables.Get("StartTime",StartTime);
  m_variables.Get("UpdateFrequency",update_frequency);
  m_variables.Get("PlotConfiguration",plot_configuration);
  m_variables.Get("PathMonitoring",path_monitoring);
  m_variables.Get("ImageFormat",img_extension);
  m_variables.Get("ForceUpdate",force_update);
  m_variables.Get("DrawMarker",draw_marker);
  m_variables.Get("TriggerMaskFile",triggermaskfile);
  m_variables.Get("TriggerWordFile",triggerwordfile);
  m_variables.Get("TriggerAlignFile",triggeralignfile);
  m_variables.Get("verbose",verbosity);

  if (verbosity > 2) std::cout <<"MonitorTrigger: Outpath (temporary): "<<outpath_temp<<std::endl;
  if (outpath_temp == "fromStore") m_data->CStore.Get("OutPath",outpath);
  else outpath = outpath_temp;
  if (verbosity > 2) std::cout <<"MonitorTrigger: Output path for plots is "<<outpath<<std::endl;
  if (update_frequency < 0.1) {
  if (verbosity > 0) std::cout <<"MonitorTrigger: Update Frequency of every "<<update_frequency<<" mins is too high. Setting default value of 5 mins."<<std::endl;
    update_frequency = 5.;
  }
  //default should be no forced update of the monitoring plots every execute step
  if (force_update !=0 && force_update !=1) {
    force_update = 0;
  }
  //check if the image format is jpg or png
  if (!(img_extension == "png" || img_extension == "jpg" || img_extension == "jpeg")){
    img_extension = "jpg";
  }
  if (triggermaskfile != "none") {
    TriggerMask = this->LoadTriggerMask(triggermaskfile);
  }
  if (triggerwordfile != "none") {
    TriggerWord = this->LoadTriggerWord(triggerwordfile);
  }
  if (triggeralignfile != "none"){
    TriggerAlign = this->LoadTriggerAlign(triggeralignfile);
  }
  
  num_triggerwords_selected = TriggerMask.size();

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


bool MonitorTrigger::Execute(){

  //Get current time, time since last timestamp, etc.
  current = (boost::posix_time::second_clock::local_time());
  duration = boost::posix_time::time_duration(current - last);
  current_stamp_duration = boost::posix_time::time_duration(current - *Epoch);
  current_stamp = current_stamp_duration.total_milliseconds();
  utc = (boost::posix_time::second_clock::universal_time());
  current_utc_duration = boost::posix_time::time_duration(utc-current);
  current_utc = current_utc_duration.total_milliseconds();
  utc_to_t = (ULong64_t) current_utc;
  utc_to_fermi = (ULong64_t) utc_to_t*MSEC_to_SEC*MSEC_to_SEC;

  bool has_trig;
  m_data->CStore.Get("HasTrigData",has_trig);

  std::string State;
  m_data->CStore.Get("State",State);

  if (has_trig){
  if (State == "Wait"){
    //-------------------------------------------------------
    //--------------No tool is executed----------------------
    //-------------------------------------------------------
    
    Log("MonitorTrigger: State is "+State,v_debug,verbosity);
  } else if (State == "DataFile"){

    Log("MonitorTrigger: New data file available.",v_message,verbosity);
    
    //-------------------------------------------------------
    //--------------MonitorTrigger executed------------------
    //-------------------------------------------------------

    //Get parsed trigger information
    m_data->CStore.Get("TimeToTriggerWordMap",TimeToTriggerWordMap);
    this->LoopThroughDecodedEvents(TimeToTriggerWordMap);

    //Write the event information to a file
    //TODO: change this to a database later on!
    //Check if data has already been written included in WriteToFile function
    this->WriteToFile();

    //Draw last file plots
    this->DrawLastFilePlots();

    //Draw customly defined plots
    this->UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

  } else {
    Log("MonitorTrigger: State not recognized: "+State,v_debug,verbosity);
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
    Log("MonitorTrigger: "+std::to_string(update_frequency)+" mins passed... Updating file history plot.",v_message,verbosity);

    last=current;
    //std::cout <<"DrawFileHistory (period update)"<<std::endl;
    this->DrawFileHistoryTrig(current_stamp,24.,"current_24h",1);     //show 24h history of Tank files
    this->PrintFileTimeStamp(current_stamp,24.,"current_24h");
    this->DrawFileHistoryTrig(current_stamp,2.,"current_2h",3);

  }

  //Only for debugging memory leaks, otherwise comment out
  //std::cout <<"List of objects (after Execute step): "<<std::endl;
  //gObjectTable->Print();

  return true;
}


bool MonitorTrigger::Finalise(){

  Log("Tool MonitorTrigger: Finalising ....",v_message,verbosity);

  //Deleting things
  
  //Histograms
  delete h_triggerword;
  delete h_triggerword_selected;
  delete h_rate_triggerword;
  delete h_rate_triggerword_selected;

  for (int i_trig=0; i_trig < num_triggerwords; i_trig++){
    delete h_timestamp.at(i_trig);
  }

  for (int i_align=0; i_align < (int) h_triggeralign.size(); i_align++){
    delete h_triggeralign.at(i_align);
  }

  //Legends
  delete leg_rate;

  //Graphs
  for (int i_trig=0; i_trig < num_triggerwords; i_trig++){
    delete gr_rate.at(i_trig);
  }

  //Multigraphs
  delete multi_trig_rate;

  //Canvas
  delete canvas_triggerword;
  delete canvas_triggerword_selected;
  delete canvas_rate_triggerword;
  delete canvas_rate_triggerword_selected;
  delete canvas_logfile_trig;
  delete canvas_file_timestamp_trig;
  delete canvas_timestamp;
  delete canvas_triggeralign;

  return true;
}

void MonitorTrigger::ReadInConfiguration(){

  Log("MonitorTrigger: ReadInConfiguration.",v_message,verbosity);

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
        Log("ERROR (MonitorTrigger): ReadInConfiguration: Need at least 4 arguments in one line: TimeFrame - TimeEnd - FileLabel - PlotType1. Please look over the configuration file and adjust it accordingly.",v_error,verbosity);
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
    Log("ERROR (MonitorTrigger): ReadInConfiguration: Could not open file "+plot_configuration+"! Check if path is valid.",v_error,verbosity);
  }
  file.close();

  if (verbosity > 2){
    std::cout <<"---------------------------------------------------------------------"<<std::endl;
    std::cout << "MonitorTrigger: ReadInConfiguration: Read in the following data into configuration variables: "<<std::endl;
    for (unsigned int i_t=0; i_t < config_timeframes.size(); i_t++){
      std::cout <<config_timeframes.at(i_t)<<", "<<config_endtime.at(i_t)<<", "<<config_label.at(i_t)<<", ";
      for (unsigned int i_plot = 0; i_plot < config_plottypes.at(i_t).size(); i_plot++){
        std::cout <<config_plottypes.at(i_t).at(i_plot)<<", ";
      }
      std::cout<<std::endl;
    }
    std::cout <<"-----------------------------------------------------------------------"<<std::endl;
  }

  Log("MonitorTrigger: ReadInConfiguration: Parsing dates: ",v_message,verbosity);
  for (unsigned int i_date = 0; i_date < config_endtime.size(); i_date++){
    if (config_endtime.at(i_date) == "TEND_LASTFILE") {
      Log("MonitorTrigger: TEND_LASTFILE: Starting from end of last read-in file",v_message,verbosity);
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    } else if (config_endtime.at(i_date).size()==15){
        boost::posix_time::ptime spec_endtime(boost::posix_time::from_iso_string(config_endtime.at(i_date)));
        boost::posix_time::time_duration spec_endtime_duration = boost::posix_time::time_duration(spec_endtime - *Epoch);
        ULong64_t spec_endtime_long = spec_endtime_duration.total_milliseconds();
        config_endtime_long.push_back(spec_endtime_long);
    } else {
      Log("MonitorTrigger: Specified end date "+config_endtime.at(i_date)+" does not have the desired format yyyymmddThhmmss. Please change the format in the config file in order to use this tool. Starting from end of last file",v_message,verbosity);
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    }
  }

}

void MonitorTrigger::InitializeHists(){

  Log("MonitorTrigger: Initialize Hists",v_message,verbosity);
  
  //-------------------------------------------------------
  //----------------InitializeHists -----------------------
  //-------------------------------------------------------

  readfromfile_tend = 0;
  readfromfile_timeframe = 0.;

  std::string str_triggerword, str_triggerword_selected, str_rate_triggerword, str_rate_triggerword_selected;
  std::stringstream ss_title_triggerword, ss_title_triggerword_selected, ss_title_rate_triggerword, ss_title_rate_triggerword_selected;
  str_triggerword = "Triggerword frequency";
  str_triggerword_selected = "Selected triggerword frequency";
  str_rate_triggerword = "Triggerword rate";
  str_rate_triggerword_selected = "Selected triggerword rate";

  ss_title_triggerword << title_time.str() << str_triggerword;
  ss_title_triggerword_selected << title_time.str() << str_triggerword_selected;
  ss_title_rate_triggerword << title_time.str() << str_rate_triggerword;
  ss_title_rate_triggerword_selected << title_time.str() << str_rate_triggerword_selected;

  gROOT->cd();

  h_triggerword = new TH1F("h_triggerword",ss_title_triggerword.str().c_str(),num_triggerwords,0,num_triggerwords);
  h_triggerword_selected = new TH1F("h_triggerword_selected",ss_title_triggerword_selected.str().c_str(),num_triggerwords_selected,0,num_triggerwords_selected);
  h_rate_triggerword = new TH1F("h_rate_triggerword",ss_title_rate_triggerword.str().c_str(),num_triggerwords,0,num_triggerwords);
  h_rate_triggerword_selected = new TH1F("h_rate_triggerword_selected",ss_title_rate_triggerword_selected.str().c_str(),num_triggerwords_selected,0,num_triggerwords_selected);
  num_files_history = 10;
  log_files_trig = new TH1F("log_files_trig","Trigger Files History",num_files_history,0,num_files_history);

  h_triggerword->SetStats(0);
  h_triggerword_selected->SetStats(0);
  h_rate_triggerword->SetStats(0);
  h_rate_triggerword_selected->SetStats(0);

  h_triggerword->GetXaxis()->SetTitle("Triggerword");
  h_triggerword->GetYaxis()->SetTitle("#");
  //h_triggerword_selected->GetXaxis()->SetTitle("Triggerword");
  h_triggerword_selected->GetYaxis()->SetTitle("#");
  h_rate_triggerword->GetXaxis()->SetTitle("Triggerword");
  h_rate_triggerword->GetYaxis()->SetTitle("Rate [Hz]");
  //h_rate_triggerword_selected->GetXaxis()->SetTitle("Triggerword");
  h_rate_triggerword_selected->GetYaxis()->SetTitle("Rate [Hz]");

  h_triggerword->SetLineWidth(2);
  h_triggerword_selected->SetLineWidth(2);
  h_rate_triggerword->SetLineWidth(2);
  h_rate_triggerword_selected->SetLineWidth(2);

  log_files_trig->GetXaxis()->SetTimeDisplay(1);
  log_files_trig->GetXaxis()->SetLabelSize(0.03);
  log_files_trig->GetXaxis()->SetLabelOffset(0.03);
  log_files_trig->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  log_files_trig->GetYaxis()->SetTickLength(0.);
  log_files_trig->GetYaxis()->SetLabelOffset(999);
  log_files_trig->SetStats(0);

  canvas_triggerword = new TCanvas("canvas_triggerword","Triggerwords",900,600);
  canvas_triggerword_selected = new TCanvas("canvas_triggerword_selected","Selected triggerwords",900,600);
  canvas_rate_triggerword = new TCanvas("canvas_rate_triggerword","Rate Triggerwords",900,600);
  canvas_rate_triggerword_selected = new TCanvas("canvas_rate_triggerword_selected","Rate Triggerwords selected",900,600);
  canvas_logfile_trig = new TCanvas("canvas_logfile_trig","Trigger File History",900,600);
  canvas_file_timestamp_trig = new TCanvas("canvas_file_timestamp_trig","Timestamp Last Trig File",900,600);
  canvas_timestamp = new TCanvas("canvas_timestamp","Timestamps",900,600);
  canvas_timeevolution = new TCanvas("canvas_timeevolution","Time Evolution",900,600);
  canvas_triggeralign = new TCanvas("canvas_triggeralign","Triggerword time alignment",900,600);

  for (int i_trig=0; i_trig < num_triggerwords; i_trig++){
    std::stringstream ss_timestamp, ss_title_timestamp;
    ss_timestamp << "h_timestamp_trig"<<i_trig;
    ss_title_timestamp << "Timestamps Triggerword "<<i_trig;
    TH1F *htemp_timestamp = new TH1F(ss_timestamp.str().c_str(),ss_title_timestamp.str().c_str(),num_time_bins,0,num_time_bins);
    htemp_timestamp->SetStats(0);
    htemp_timestamp->SetLineWidth(2);
    htemp_timestamp->GetYaxis()->SetTitle("#");
    htemp_timestamp->GetXaxis()->SetTimeDisplay(1);
    htemp_timestamp->GetXaxis()->SetLabelSize(0.03);
    htemp_timestamp->GetXaxis()->SetLabelOffset(0.03);
    htemp_timestamp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    h_timestamp.push_back(htemp_timestamp);
  }

  for (int i_align = 0; i_align < (int) TriggerAlign.size(); i_align++){
    std::stringstream ss_align, ss_title_align;
    std::vector<int> single_align = TriggerAlign.at(i_align);
    ss_align << "h_trig_align"<<i_align;
    ss_title_align << "Trigger alignment ("<<single_align.at(0)<<"/"<<single_align.at(1)<<")";
    TH1F *htemp_align = new TH1F(ss_align.str().c_str(),ss_title_align.str().c_str(),num_align_bins,single_align.at(2),single_align.at(3));
    htemp_align->SetStats(0);
    htemp_align->SetLineWidth(2);
    htemp_align->GetYaxis()->SetTitle("#");
    htemp_align->GetXaxis()->SetTitle("#Delta t [ns]");
    h_triggeralign.push_back(htemp_align);
  }

  for (int i_trig = 0; i_trig < num_triggerwords; i_trig++){
    std::stringstream ss_rate, ss_title_rate;
    ss_rate << "rate_trigword"<<i_trig;
    ss_title_rate << "Rate Triggerword "<<i_trig;

    TGraph *graph_trig_rate = new TGraph();
    graph_trig_rate->SetName(ss_rate.str().c_str());
    graph_trig_rate->SetTitle(ss_title_rate.str().c_str());

    if (draw_marker) graph_trig_rate->SetMarkerStyle(20);

    graph_trig_rate->SetMarkerColor(i_trig%4+1);
    graph_trig_rate->SetLineColor(i_trig%4+1);
    graph_trig_rate->SetLineWidth(2);
    graph_trig_rate->SetFillColor(0);
    graph_trig_rate->GetYaxis()->SetTitle("Rate [Hz]");
    graph_trig_rate->GetXaxis()->SetTimeDisplay(1);
    graph_trig_rate->GetXaxis()->SetLabelSize(0.03);
    graph_trig_rate->GetXaxis()->SetLabelOffset(0.03);
    graph_trig_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    
    gr_rate.push_back(graph_trig_rate);

  }

  multi_trig_rate = new TMultiGraph();
  leg_rate = new TLegend(0.7,0.7,0.88,0.88);
  leg_rate->SetLineColor(0);

}

void MonitorTrigger::LoopThroughDecodedEvents(std::map<uint64_t,std::vector<uint32_t>> timetotriggerword){

  Log("MonitorTrigger: LoopThroughDecodedEvents",v_message,verbosity);

  //-------------------------------------------------------
  //-------------LoopThroughDecodedEvents -----------------
  //-------------------------------------------------------

  triggerword_file.clear();
  frequency_file.clear();
  timestamp_file.clear();

  frequency_file.assign(num_triggerwords,0);

  int i_timestamp = 0;
  for (std::map<uint64_t, std::vector<uint32_t>>::iterator it = timetotriggerword.begin(); it != timetotriggerword.end(); it++){

    uint64_t timestamp = it->first;
    uint64_t timestamp_temp = timestamp - utc_to_fermi;
    std::vector<uint32_t> trigword_vec = it->second;
    for (int i_trig=0; i_trig < (int) trigword_vec.size(); i_trig++){
      timestamp_file.push_back(timestamp_temp);
      uint32_t trigword = trigword_vec.at(i_trig)-1;	//Triggerwords in timetotriggerword are index+1, subtract 1 to get index
      triggerword_file.push_back(trigword);
      frequency_file.at(trigword)++;
      i_timestamp++;
    }
  }

}

void MonitorTrigger::WriteToFile(){

  Log("MonitorTrigger: WriteToFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------WriteToFile -------------------------
  //-------------------------------------------------------

  t_file_start = timestamp_file.at(0)/1000000.;
  t_file_end = timestamp_file.at(timestamp_file.size()-1)/1000000.;
  std::string file_start_date = convertTimeStamp_to_Date(t_file_start);
  std::stringstream root_filename;
  root_filename << path_monitoring << "Trigger_" << file_start_date <<".root";

  Log("MonitorTrigger: ROOT filename: "+root_filename.str(),v_message,verbosity);
  
  std::string root_option = "RECREATE";
  if (does_file_exist(root_filename.str())) root_option = "UPDATE";
  TFile *f = new TFile(root_filename.str().c_str(),root_option.c_str());

  ULong64_t t_start, t_end, t_frame;
  std::vector<int> *trigword = new std::vector<int>;
  std::vector<int> *freq_trigword = new std::vector<int>;
  std::vector<double> *rate_trigword = new std::vector<double>;

  TTree *t;
  if (f->GetListOfKeys()->Contains("triggermonitor_tree")) {
    Log("MonitorTrigger: WriteToFile: Tree already exists",v_message,verbosity);
    t = (TTree*) f->Get("triggermonitor_tree");
    t->SetBranchAddress("t_start",&t_start);
    t->SetBranchAddress("t_end",&t_end);
    t->SetBranchAddress("trigword",&trigword);
    t->SetBranchAddress("freq_trigword",&freq_trigword);
    t->SetBranchAddress("rate_trigword",&rate_trigword);
  } else {
    t = new TTree("triggermonitor_tree","Trigger Monitoring tree");
    Log("MonitorTrigger: WriteToFile: Tree is created from scratch",v_message,verbosity);
    t->Branch("t_start",&t_start);
    t->Branch("t_end",&t_end);
    t->Branch("trigword",&trigword);
    t->Branch("freq_trigword",&freq_trigword);
    t->Branch("rate_trigword",&rate_trigword);
  }

  int n_entries = t->GetEntries();
  bool omit_entries = false;
  for (int i_entry = 0; i_entry < n_entries; i_entry++){
    t->GetEntry(i_entry);
    if (t_start == t_file_start) {
      Log("WARNING (MonitorTrigger): WriteToFile: Wanted to write data from file that is already written to DB. Omit entries.",v_warning,verbosity);
      omit_entries = true;
    }
  }

  //If data is already written to DB/File, do not write it again
  if (omit_entries){
    //Don't write file again, but still delete TFile and TTree object!
    f->Close();
    delete trigword;
    delete freq_trigword;
    delete rate_trigword;
    delete f;

    gROOT->cd();

    return;

  }

  trigword->clear();
  freq_trigword->clear();
  rate_trigword->clear();

  t_start = t_file_start;
  t_end = t_file_end;
  t_frame = t_end - t_start;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(t_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_start/MSEC_to_SEC/SEC_to_MIN)%60,int(t_start/MSEC_to_SEC/1000.)%60,t_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  Log("MonitorTrigger: WriteToFile: Writing data to file: "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
    +"..."+std::to_string(endtime_tm.tm_year+1900)+"/"+std::to_string(endtime_tm.tm_mon+1)+"/"+std::to_string(endtime_tm.tm_mday)+"-"+std::to_string(endtime_tm.tm_hour)+":"+std::to_string(endtime_tm.tm_min)+":"+std::to_string(endtime_tm.tm_sec),v_message,verbosity);

  for (int i_trg = 0; i_trg < num_triggerwords; i_trg++){
    trigword->push_back(i_trg);
    freq_trigword->push_back(frequency_file.at(i_trg));
    if (t_frame > 0.) rate_trigword->push_back(frequency_file.at(i_trg)/(t_frame/MSEC_to_SEC));
    else rate_trigword->push_back(0.);
  }

  t->Fill();
  t->Write("",TObject::kOverwrite);
  f->Close();

  delete trigword;
  delete freq_trigword;
  delete rate_trigword;

  delete f;

  gROOT->cd();

}

void MonitorTrigger::ReadFromFile(ULong64_t timestamp_end, double time_frame){

  Log("MonitorTrigger: ReadFromFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------ReadFromFile ------------------------
  //-------------------------------------------------------

  triggerword_plot.clear();
  frequency_plot.clear();
  rate_plot.clear();
  tstart_plot.clear();
  tend_plot.clear();
  labels_timeaxis.clear();

  //take the end time and calculate the start time with the given time_frame
  ULong64_t timestamp_start = timestamp_end - time_frame*MIN_to_HOUR*SEC_to_MIN*MSEC_to_SEC;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(timestamp_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_start/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_start/MSEC_to_SEC/1000.)%60,timestamp_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);

  Log("MonitorTriggerTime: ReadFromFile: Reading in data for time frame "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
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
    root_filename_i << path_monitoring << "Trigger_" << string_date_i <<".root";
    bool tree_exists = true;  

    if (does_file_exist(root_filename_i.str())) {
      TFile *f = new TFile(root_filename_i.str().c_str(),"READ");
      TTree *t;
      if (f->GetListOfKeys()->Contains("triggermonitor_tree")) t = (TTree*) f->Get("triggermonitor_tree");
      else {
        Log("WARNING (MonitorTrigger): File "+root_filename_i.str()+" does not contain triggermonitor_tree. Omit file.",v_warning,verbosity);
        tree_exists = false;
      }
 
      if (tree_exists){

        Log("MonitorTrigger: Tree exists, start reading in data",v_message,verbosity);
  
        ULong64_t t_start, t_end;
        std::vector<int> *trigword = new std::vector<int>;
        std::vector<int> *freq_trigword = new std::vector<int>;
        std::vector<double> *rate_trigword = new std::vector<double>;
  
        int nevents;
        int nentries_tree;
  
        t->SetBranchAddress("t_start",&t_start);
        t->SetBranchAddress("t_end",&t_end);
        t->SetBranchAddress("trigword",&trigword);
        t->SetBranchAddress("freq_trigword",&freq_trigword);
        t->SetBranchAddress("rate_trigword",&rate_trigword);
  
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
            
            triggerword_plot.push_back(*trigword);
            frequency_plot.push_back(*freq_trigword);
            rate_plot.push_back(*rate_trigword);
            tstart_plot.push_back(t_start);
            tend_plot.push_back(t_end);
            boost::posix_time::ptime boost_tend = *Epoch+boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
            struct tm label_timestamp = boost::posix_time::to_tm(boost_tend);
  
            TDatime datime_timestamp(1900+label_timestamp.tm_year,label_timestamp.tm_mon+1,label_timestamp.tm_mday,label_timestamp.tm_hour,label_timestamp.tm_min,label_timestamp.tm_sec);
            labels_timeaxis.push_back(datime_timestamp);
          }
        }
        
        delete trigword;
        delete freq_trigword;
        delete rate_trigword;
  
      }
  
      f->Close();
      delete f;
      gROOT->cd();

    } else {
      Log("MonitorTrigger: ReadFromFile: File "+root_filename_i.str()+" does not exist. Omit file.",v_warning,verbosity);
    }

  }

  //Set the readfromfile time variables to make sure data is not read twice for the same time window
  readfromfile_tend = timestamp_end;
  readfromfile_timeframe = time_frame;

}

void MonitorTrigger::DrawLastFilePlots(){

  Log("MonitorTrigger: DrawLastFilePlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawLastFilePlots -------------------
  //-------------------------------------------------------

  //Draw triggerword frequency & rate histograms
  DrawFrequencyRatePlots(t_file_end,(t_file_end-t_file_start)/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR,"lastFile");

  //Draw time alignment plots
  DrawTimeAlignmentPlots();

}

void MonitorTrigger::UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes){

  Log("MonitorTrigger: UpdateMonitorPlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------UpdateMonitorPlots ------------------
  //-------------------------------------------------------

  //Draw the monitoring plots according to the specifications in the configfiles
  
  for (unsigned int i_time = 0; i_time < timeFrames.size(); i_time++){
    
    ULong64_t zero = 0;
    if (endTimes.at(i_time) == zero) endTimes.at(i_time) = t_file_end;     //set 0 for t_file_end since we did not know what that was at the beginning of initialise

    for (unsigned int i_plot = 0; i_plot < plotTypes.at(i_time).size(); i_plot++){

      if (plotTypes.at(i_time).at(i_plot) == "TriggerRatePlots") this->DrawFrequencyRatePlots(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "TimeEvolution") this->DrawTimeEvolution(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else if (plotTypes.at(i_time).at(i_plot) == "FileHistory") {this->DrawFileHistoryTrig(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time),1);}
      else {
        if (verbosity > 0) std::cout <<"ERROR (MonitorTrigger): UpdateMonitorPlots: Specified plot type -"<<plotTypes.at(i_time).at(i_plot)<<"- does not exist! Omit entry."<<std::endl;
      }
    }

  }

}

void MonitorTrigger::DrawFrequencyRatePlots(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTrigger: DrawFrequencyRatePlots",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //-------------------------------------------------------
  //-------------DrawFrequencyRatePlots -------------------
  //-------------------------------------------------------
  
  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  std::vector<double> overall_rates;
  std::vector<double> overall_freqs;
  ULong64_t current_timeframe;
  ULong64_t overall_timeframe = 0;
  overall_rates.assign(num_triggerwords,0.);
  overall_freqs.assign(num_triggerwords,0);

  for (unsigned int i_file = 0; i_file < tend_plot.size(); i_file++){

    current_timeframe = (tend_plot.at(i_file) - tstart_plot.at(i_file))/MSEC_to_SEC;
    overall_timeframe += current_timeframe;

    for (int i_trg = 0; i_trg < num_triggerwords; i_trg++){
      overall_rates.at(i_trg) += frequency_plot.at(i_file).at(i_trg);
      overall_freqs.at(i_trg) += frequency_plot.at(i_file).at(i_trg);
    }
  }
  
  if (overall_timeframe > 0.){
    for (int i_trg = 0; i_trg < num_triggerwords; i_trg++){
      overall_rates.at(i_trg) /= overall_timeframe;
    }
  }

  gROOT->cd();

  int i_selected=0;
  for (int i_trg = 0; i_trg < num_triggerwords; i_trg++){
    h_triggerword->SetBinContent(i_trg+1,overall_freqs.at(i_trg));
    h_rate_triggerword->SetBinContent(i_trg+1,overall_rates.at(i_trg));
    if (std::find(TriggerMask.begin(),TriggerMask.end(),i_trg)!=TriggerMask.end()){
      h_triggerword_selected->SetBinContent(i_selected+1,overall_freqs.at(i_trg));
      h_rate_triggerword_selected->SetBinContent(i_selected+1,overall_rates.at(i_trg));
      i_selected++;
    }
  }

  std::stringstream ss_title_freq, ss_title_rate;
  ss_title_freq << "Freq Triggerword (last "<<ss_timeframe.str()<<"h) "<<end_time.str()<<std::endl;
  ss_title_rate << "Rate Triggerword (last "<<ss_timeframe.str()<<"h) "<<end_time.str()<<std::endl;
  h_triggerword->SetTitle(ss_title_freq.str().c_str());
  h_rate_triggerword->SetTitle(ss_title_rate.str().c_str());
  h_triggerword_selected->SetTitle(ss_title_freq.str().c_str());
  h_rate_triggerword_selected->SetTitle(ss_title_rate.str().c_str());

  for (int i_sel = 0; i_sel < (int) TriggerMask.size(); i_sel++){
    h_triggerword_selected->GetXaxis()->SetBinLabel(i_sel+1,TriggerWord[TriggerMask.at(i_sel)].c_str());
    h_rate_triggerword_selected->GetXaxis()->SetBinLabel(i_sel+1,TriggerWord[TriggerMask.at(i_sel)].c_str());
  }
  //h_triggerword_selected->LabelsOption("v");
  //h_rate_triggerword_selected->LabelsOption("v");
  h_triggerword_selected->GetXaxis()->SetLabelSize(0.03);
  canvas_triggerword->cd();
  h_triggerword->Draw();
  std::stringstream ss_triggerfreq;
  ss_triggerfreq << outpath << "Trigger_Freq_" << file_ending << "." << img_extension;
  canvas_triggerword->SaveAs(ss_triggerfreq.str().c_str());
  canvas_triggerword->Clear();

  canvas_rate_triggerword->cd();
  h_rate_triggerword->Draw();
  std::stringstream ss_triggerrate;
  ss_triggerrate << outpath << "Trigger_Rate_" << file_ending << "." << img_extension;
  canvas_rate_triggerword->SaveAs(ss_triggerrate.str().c_str());
  canvas_rate_triggerword->Clear();

  canvas_triggerword_selected->cd();
  h_triggerword_selected->Draw();
  std::stringstream ss_triggerfreq_selected;
  ss_triggerfreq_selected << outpath << "Trigger_Freq_selected_" << file_ending << "." << img_extension;
  canvas_triggerword_selected->SaveAs(ss_triggerfreq_selected.str().c_str());
  canvas_triggerword_selected->Clear();

  canvas_rate_triggerword_selected->cd();
  h_rate_triggerword_selected->Draw();
  std::stringstream ss_triggerrate_selected;
  ss_triggerrate_selected << outpath << "Trigger_Rate_selected_" << file_ending << "." << img_extension;
  canvas_rate_triggerword_selected->SaveAs(ss_triggerrate_selected.str().c_str());
  canvas_rate_triggerword_selected->Clear();

}

void MonitorTrigger::DrawTimeAlignmentPlots(){

  Log("MonitorTrigger: DrawTimeAlignmentPlots",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_file_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_file_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_file_end/MSEC_to_SEC/1000.)%60,t_file_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  //-------------------------------------------------------
  //-------------DrawTimeAlignmentPlots -------------------
  //-------------------------------------------------------

  std::stringstream ss_title_hist_timestamp;
  std::stringstream ss_canvas_timestamp;
  std::stringstream ss_title_hist_triggeralign;
  std::stringstream ss_canvas_triggeralign;

  for (int i_trig = 0; i_trig < num_triggerwords; i_trig++){
    if (std::find(TriggerMask.begin(),TriggerMask.end(),i_trig)==TriggerMask.end()) continue;
    h_timestamp.at(i_trig)->SetBins(num_time_bins,(t_file_start+utc_to_t)/MSEC_to_SEC,(t_file_end+utc_to_t)/MSEC_to_SEC);
    h_timestamp.at(i_trig)->GetXaxis()->SetTimeOffset(0.);
    h_timestamp.at(i_trig)->Reset();
    h_timestamp.at(i_trig)->SetMaximum(-1111);
    h_timestamp.at(i_trig)->SetMinimum(-1111);
  }


  uint64_t previous_filestamp = 0;
  int previous_triggerword=-1;
  for (int t=0; t < (int) triggerword_file.size(); t++){
    if (t==0){
      previous_filestamp = timestamp_file.at(t);
      previous_triggerword = triggerword_file.at(t);
    }
    //std::cout <<"Current word: "<<triggerword_file.at(t)<<", previous word: "<<previous_triggerword;
    //std::cout <<", delta t previous filestamp: "<<timestamp_file.at(t)-previous_filestamp<<std::endl;
    int triggerword = triggerword_file.at(t);
    h_timestamp.at(triggerword)->Fill((timestamp_file.at(t)/1000000+utc_to_t)/MSEC_to_SEC);
    for (int i_align=0; i_align < (int) TriggerAlign.size(); i_align++){
      std::vector<int> single_align = TriggerAlign.at(i_align);
      if (triggerword == single_align.at(0) && previous_triggerword == single_align.at(1)) h_triggeralign.at(i_align)->Fill(timestamp_file.at(t)-previous_filestamp);
    }
    previous_filestamp = timestamp_file.at(t);
    previous_triggerword = triggerword_file.at(t);
  }

  for (int i_trig = 0; i_trig < num_triggerwords; i_trig++){
    if (std::find(TriggerMask.begin(),TriggerMask.end(),i_trig)==TriggerMask.end()) continue;
    canvas_timestamp->cd();
    ss_title_hist_timestamp.str("");
    ss_title_hist_timestamp <<TriggerWord.at(i_trig)<<" Timestamps (lastFile) "<< end_time.str();
    h_timestamp.at(i_trig)->SetTitle(ss_title_hist_timestamp.str().c_str());
    if (h_timestamp.at(i_trig)->GetEntries()>0) h_timestamp.at(i_trig)->GetYaxis()->SetRangeUser(0.,1.1*h_timestamp.at(i_trig)->GetMaximum());
    else h_timestamp.at(i_trig)->GetYaxis()->SetRangeUser(0,1.1);
    h_timestamp.at(i_trig)->Draw();
    ss_canvas_timestamp.str("");
    ss_canvas_timestamp << outpath <<"Trigger_Timestamps_"<<TriggerWord.at(i_trig) << "." << img_extension;
    canvas_timestamp->SaveAs(ss_canvas_timestamp.str().c_str());
    canvas_timestamp->Clear();
  }

  for (int i_align = 0; i_align < (int) TriggerAlign.size(); i_align++){
    canvas_triggeralign->cd();
    std::vector<int> single_align = TriggerAlign.at(i_align);
    ss_title_hist_triggeralign.str("");
    ss_title_hist_triggeralign <<"Alignment trgwords "<<single_align.at(0)<<"/"<<single_align.at(1)<<" (lastFile) "<< end_time.str();
    h_triggeralign.at(i_align)->SetTitle(ss_title_hist_triggeralign.str().c_str());
    h_triggeralign.at(i_align)->Draw();
    ss_canvas_triggeralign.str("");
    ss_canvas_triggeralign << outpath << "Trigger_Align_"<<single_align.at(0)<<"_"<<single_align.at(1)<< "." << img_extension;
    canvas_triggeralign->SaveAs(ss_canvas_triggeralign.str().c_str());
    canvas_triggeralign->Clear();
  }

}

void MonitorTrigger::DrawTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorTrigger: DrawTimeEvolution",v_message,verbosity);

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

  //Resetting time evolution graphs
  for (int i_trig=0; i_trig < num_triggerwords; i_trig++){
    gr_rate.at(i_trig)->Set(0);
  }

  for (unsigned int i_file = 0; i_file < tend_plot.size(); i_file++){
    for (int i_trig = 0; i_trig < num_triggerwords; i_trig++){
      gr_rate.at(i_trig)->SetPoint(i_file,labels_timeaxis[i_file].Convert(),rate_plot.at(i_file).at(i_trig));
    }
  }

  // Drawing time evolution plots

  double max_canvas = 0;
  double min_canvas = 999999.;

  int CH_per_CANVAS = 4;	//channels per canvas

  std::stringstream ss_leg_time, ss_trig_rate, ss_title_trig_rate;

  for (int i_trig = 0; i_trig < num_triggerwords; i_trig++){
  
    if (i_trig%CH_per_CANVAS == 0 || i_trig == num_triggerwords-1){
    
      if (i_trig != 0){

        ss_trig_rate.str("");
        int min_triggerword = (i_trig!=num_triggerwords-1)? i_trig-CH_per_CANVAS : i_trig-CH_per_CANVAS+1;
        int max_triggerword = (i_trig==num_triggerwords-1)? i_trig : i_trig-1;
        
        ss_trig_rate << "Rates Triggerwords "<<min_triggerword<<" - "<<max_triggerword<<" (last "<<ss_timeframe.str()<<"h) "<<end_time.str();

        if (i_trig == num_triggerwords-1){
          ss_leg_time.str("");
          ss_leg_time << "trgword "<<i_trig;
          multi_trig_rate->Add(gr_rate.at(i_trig));
          leg_rate->AddEntry(gr_rate.at(i_trig),ss_leg_time.str().c_str(),"l");
        }

        canvas_timeevolution->cd();
        multi_trig_rate->Draw("apl");
        multi_trig_rate->SetTitle(ss_trig_rate.str().c_str());
        multi_trig_rate->GetYaxis()->SetTitle("Rate [Hz]");
        multi_trig_rate->GetXaxis()->SetTimeDisplay(1);
        multi_trig_rate->GetXaxis()->SetLabelSize(0.03);
        multi_trig_rate->GetXaxis()->SetLabelOffset(0.03);
        multi_trig_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
        multi_trig_rate->GetXaxis()->SetTimeOffset(0.);
        leg_rate->Draw();

        ss_title_trig_rate.str("");
        ss_title_trig_rate << outpath << "TrigTimeEvolution_Rate_Trigword"<<i_trig<<"_"<<file_ending<<"."<<img_extension;
        canvas_timeevolution->SaveAs(ss_title_trig_rate.str().c_str());

        for (int i_gr=0; i_gr < CH_per_CANVAS; i_gr++){
          int i_balance = (i_trig == num_triggerwords-1)? 1: 0;
          multi_trig_rate->RecursiveRemove(gr_rate.at(i_trig-CH_per_CANVAS+i_gr+i_balance));
        }
      }

      leg_rate->Clear();
      canvas_timeevolution->Clear();

    }

    if (i_trig != num_triggerwords-1){
      ss_leg_time.str("");
      ss_leg_time << "trgword "<<i_trig;
      multi_trig_rate->Add(gr_rate.at(i_trig));
      if (gr_rate.at(i_trig)->GetMaximum() > max_canvas) max_canvas = gr_rate.at(i_trig)->GetMaximum();
      if (gr_rate.at(i_trig)->GetMinimum() < min_canvas) min_canvas = gr_rate.at(i_trig)->GetMinimum();
      leg_rate->AddEntry(gr_rate.at(i_trig),ss_leg_time.str().c_str(),"l");
    }

  }

}

std::string MonitorTrigger::convertTimeStamp_to_Date(ULong64_t timestamp){
  
  //format of date is YYYY_MM-DD

  boost::posix_time::ptime convertedtime = *Epoch + boost::posix_time::time_duration(int(timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp/MSEC_to_SEC/1000.)%60,timestamp%1000);
  struct tm converted_tm = boost::posix_time::to_tm(convertedtime);
  std::stringstream ss_date;
  ss_date << converted_tm.tm_year+1900 << "_" << converted_tm.tm_mon+1 << "-" <<converted_tm.tm_mday;
  return ss_date.str();

}


bool MonitorTrigger::does_file_exist(std::string filename){

  std::ifstream infile(filename.c_str());
  bool file_good = infile.good();
  infile.close();
  return file_good;

}


void MonitorTrigger::DrawFileHistoryTrig(ULong64_t timestamp_end, double time_frame, std::string file_ending, int _linewidth){

  Log("MonitorTrigger: DrawFileHistoryTrig",v_message,verbosity);

  //-------------------------------------------------------
  //------------------DrawFileHistoryTrig ---------------------
  //-------------------------------------------------------

  //Creates a plot showing the time stamps for all the files within the last time_frame mins
  //The plot is updated with the update_frequency specified in the configuration file (default: 5 mins)

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  timestamp_end += utc_to_t;

  ULong64_t timestamp_start = timestamp_end - time_frame*MSEC_to_SEC*SEC_to_MIN*MIN_to_HOUR;
  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  canvas_logfile_trig->cd();
  log_files_trig->SetBins(num_files_history,timestamp_start/MSEC_to_SEC,timestamp_end/MSEC_to_SEC);
  log_files_trig->GetXaxis()->SetTimeOffset(0.);
  log_files_trig->Draw();

  std::stringstream ss_title_filehistory;
  ss_title_filehistory << "Trigger Files History (last "<<ss_timeframe.str()<<"h)";

  log_files_trig->SetTitle(ss_title_filehistory.str().c_str());

  std::vector<TLine*> file_markers_trig;
  for (unsigned int i_file = 0; i_file < tend_plot.size(); i_file++){
    TLine *line_file = new TLine((tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,0.,(tend_plot.at(i_file)+utc_to_t)/MSEC_to_SEC,1.);
    line_file->SetLineColor(1);
    line_file->SetLineStyle(1);
    line_file->SetLineWidth(_linewidth);
    line_file->Draw("same");
    file_markers_trig.push_back(line_file);
  }

  std::stringstream ss_logfiles;
  ss_logfiles << outpath << "Trigger_FileHistory_" << file_ending << "." << img_extension;
  canvas_logfile_trig->SaveAs(ss_logfiles.str().c_str());
  

  for (unsigned int i_line = 0; i_line < file_markers_trig.size(); i_line++){
    delete file_markers_trig.at(i_line);
  }

  log_files_trig->Reset();
  canvas_logfile_trig->Clear();


}

void MonitorTrigger::PrintFileTimeStamp(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  if (verbosity > 2) std::cout <<"MonitorTrigger: PrintFileTimeStamp"<<std::endl;

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
    time_diff << "#Delta t Last Trigger File: >"<<time_frame<<"h";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  } else {
    ULong64_t timestamp_lastfile = tend_plot.at(tend_plot.size()-1);
    boost::posix_time::ptime filetime = *Epoch + boost::posix_time::time_duration(int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_lastfile/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_lastfile/MSEC_to_SEC/1000.)%60,timestamp_lastfile%1000);
    boost::posix_time::time_duration t_since_file= boost::posix_time::time_duration(endtime - filetime);
    int t_since_file_min = int(t_since_file.total_milliseconds()/MSEC_to_SEC/SEC_to_MIN);
    std::stringstream time_diff;
    time_diff << "#Delta t Last Trigger File: "<<t_since_file_min/60<<"h:"<<t_since_file_min%60<<"min";
    label_timediff = new TLatex(0.04,0.33,time_diff.str().c_str());
    label_timediff->SetNDC(1);
    label_timediff->SetTextSize(0.1);
  }

  canvas_file_timestamp_trig->cd();
  label_lastfile->Draw();
  label_timediff->Draw();
  std::stringstream ss_file_timestamp;
  ss_file_timestamp << outpath << "Trigger_FileTimeStamp_" << file_ending << "." << img_extension;
  canvas_file_timestamp_trig->SaveAs(ss_file_timestamp.str().c_str());

  delete label_lastfile;
  delete label_timediff;

  canvas_file_timestamp_trig->Clear();

}

std::vector<int> MonitorTrigger::LoadTriggerMask(std::string triggermask_file){
  std::vector<int> trigger_mask;
  std::string fileline;
  ifstream myfile(triggermask_file.c_str());
  if (myfile.is_open()){
    while(getline(myfile,fileline)){
      if(fileline.find("#")!=std::string::npos) continue;
      std::vector<std::string> dataline;
      boost::split(dataline,fileline, boost::is_any_of(","), boost::token_compress_on);
      uint32_t triggernum = std::stoul(dataline.at(0));
      if(verbosity>4) std::cout << "Trigger mask will have trigger number " << triggernum << std::endl;
      trigger_mask.push_back(triggernum);
    }
  } else {
    Log("MonitorTrigger Tool: Input trigger mask file not found. "
        " all triggers from CTC will attempt to be paired with PMT/MRD data. ",
        v_warning, verbosity);
  }
  return trigger_mask;
}

std::map<int,std::string> MonitorTrigger::LoadTriggerWord(std::string file_triggerword){
  std::map<int,std::string> trigger_word;
  ifstream myfile(file_triggerword.c_str());
  int triggerword;
  std::string triggerlabel;
  while (!myfile.eof()){
    myfile >> triggerword >> triggerlabel;
    trigger_word.emplace(triggerword,triggerlabel);
    if (myfile.eof()) break;
  }
  return trigger_word;

}

std::vector<std::vector<int>> MonitorTrigger::LoadTriggerAlign(std::string file_triggeralign){
  std::vector<std::vector<int>> trigger_align;
  ifstream myfile(file_triggeralign.c_str());
  int trigword1, trigword2, bin_min, bin_max;
  if (myfile.is_open()){
    while (!myfile.eof()){
      myfile >> trigword1 >> trigword2 >> bin_min >> bin_max;
      std::vector<int> single_trigger_align{trigword1,trigword2,bin_min,bin_max};
      trigger_align.push_back(single_trigger_align);
      if (myfile.eof()) break;
    }
  } else {
    Log("MonitorTrigger Tool: Input trigger alignment file not found. "
	" No alignment plots will be produced.", v_warning,verbosity);
  }
  return trigger_align;
}
