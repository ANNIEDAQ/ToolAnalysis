#include "MonitorDAQ.h"

MonitorDAQ::MonitorDAQ():Tool(){}


bool MonitorDAQ::Initialise(std::string configfile, DataModel &data){

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
  m_variables.Get("PathCompStats",path_compstats);
  m_variables.Get("ImageFormat",img_extension);
  m_variables.Get("ForceUpdate",force_update);
  m_variables.Get("DrawMarker",draw_marker);
  m_variables.Get("UseOnline",online);
  m_variables.Get("SendSlack",send_slack);
  m_variables.Get("Hook",hook);
  m_variables.Get("verbose",verbosity);
  m_variables.Get("TestMode",testmode);

  if (verbosity > 2) std::cout <<"MonitorDAQ: Outpath (temporary): "<<outpath_temp<<std::endl;
  if (outpath_temp == "fromStore") m_data->CStore.Get("OutPath",outpath);
  else outpath = outpath_temp;
  if (verbosity > 2) std::cout <<"MonitorDAQ: Output path for plots is "<<outpath<<std::endl;
  if (update_frequency < 0.1) {
  if (verbosity > 0) std::cout <<"MonitorDAQ: Update Frequency of every "<<update_frequency<<" mins is too high. Setting default value of 5 mins."<<std::endl;
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
  if (online !=0 && online !=1){
    online = 0;
  }
  if (send_slack !=0 && send_slack !=1){
    send_slack = 0;
  }
  if (testmode != 0 && testmode != 1){
    testmode = 0;
  }

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

  //std::cout <<"Hook: "<<hook<<std::endl;
  
  if (online){
    address = "239.192.1.1";
    context=new zmq::context_t(3);
    SD=new ServiceDiscovery(address,port,context,320);
  } 

  if (testmode){
    ifstream testfile("./configfiles/Monitoring/SimMonitoring/TestMonitorDAQ.txt");
    double temp_file, temp_disk;
    bool temp_trig, temp_pmt;
    int temp_vme;
    while (!testfile.eof()){
      testfile >> temp_file >> temp_trig >> temp_pmt >> temp_vme >> temp_disk;
      if (testfile.eof()) break;
      test_filesize.push_back(temp_file);
      test_hastrig.push_back(temp_trig);
      test_haspmt.push_back(temp_pmt);
      test_vme.push_back(temp_vme);
      test_disk.push_back(temp_disk);
    }
    testfile.close();
    testcounter=0;
  }

  return true;
}


bool MonitorDAQ::Execute(){

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

  bool has_file;
  m_data->CStore.Get("HasNewFile",has_file);

  std::string State;
  m_data->CStore.Get("State",State);

  if (has_file){
  if (State == "Wait"){
    //-------------------------------------------------------
    //--------------No tool is executed----------------------
    //-------------------------------------------------------
    
    Log("MonitorDAQ: State is "+State,v_debug,verbosity);
  } else if (State == "DataFile"){

    Log("MonitorDAQ: New data file available.",v_message,verbosity);
    
    //-------------------------------------------------------
    //--------------MonitorDAQ executed------------------
    //-------------------------------------------------------

    file_produced = true;

    //Get parsed trigger information
    this->GetFileInformation();

    //Write the event information to a file
    //TODO: change this to a database later on!
    //Check if data has already been written included in WriteToFile function
    this->WriteToFile();

    //Draw customly defined plots
    this->UpdateMonitorPlots(config_timeframes, config_endtime_long, config_label, config_plottypes);

  } else {
    Log("MonitorDAQ: State not recognized: "+State,v_debug,verbosity);
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
    Log("MonitorDAQ: "+std::to_string(update_frequency)+" mins passed... Updating file history plot.",v_message,verbosity);

    last=current;
    this->GetVMEServices(online);
    this->GetCompStats();
    this->DrawVMEService(current_stamp,24.,"current_24h",1);     //show 24h history of Tank files
    this->PrintInfoBox();

  }

  //Only for debugging memory leaks, otherwise comment out --> test mode
  if (testmode){
    //std::cout <<"List of objects (after Execute step): "<<std::endl;
    //gObjectTable->Print();
    testcounter++;
  }

  return true;
}


bool MonitorDAQ::Finalise(){

  Log("Tool MonitorDAQ: Finalising ....",v_message,verbosity);

  //Deleting things
  
  //Graphs
  delete gr_filesize;
  delete gr_vmeservice;
  delete gr_disk_daq01;
  delete gr_disk_vme01;
  delete gr_disk_vme02;
  delete gr_disk_vme03;
  delete gr_mem_daq01;
  delete gr_mem_vme01;
  delete gr_mem_vme02;
  delete gr_mem_vme03;
  delete gr_cpu_daq01;
  delete gr_cpu_vme01;
  delete gr_cpu_vme02;
  delete gr_cpu_vme03;
  delete multi_disk;
  delete multi_mem;
  delete multi_cpu;

  //Legends
  delete leg_mem;
  delete leg_cpu;
  delete leg_disk;

  //Pie charts
  delete pie_vme;

  //Text blocks
  delete text_summary;
  delete text_vmeservice;
  delete text_currentdate;
  delete text_filesize;
  delete text_filedate;
  delete text_filename;
  delete text_haspmt;
  delete text_hascc;
  delete text_hastrigger;
  delete text_disk_title;
  delete text_disk_daq01;
  delete text_disk_vme01;
  delete text_disk_vme02;
  delete text_disk_vme03;
  delete text_mem_daq01;
  delete text_mem_vme01;
  delete text_mem_vme02;
  delete text_mem_vme03;

  //Canvas
  delete canvas_infobox;
  delete canvas_vmeservice;
  delete canvas_timeevolution_size;
  delete canvas_timeevolution_vme;

  delete canvas_info_diskspace;
  delete canvas_timeevolution_mem;
  delete canvas_timeevolution_cpu;
  delete canvas_timeevolution_disk;

  //Online things
  if (online){
    delete SD;
    //delete context;	//Deleting the context throws an error?
  }

  return true;
}

void MonitorDAQ::ReadInConfiguration(){

  Log("MonitorDAQ: ReadInConfiguration.",v_message,verbosity);

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
        Log("ERROR (MonitorDAQ): ReadInConfiguration: Need at least 4 arguments in one line: TimeFrame - TimeEnd - FileLabel - PlotType1. Please look over the configuration file and adjust it accordingly.",v_error,verbosity);
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
    Log("ERROR (MonitorDAQ): ReadInConfiguration: Could not open file "+plot_configuration+"! Check if path is valid.",v_error,verbosity);
  }
  file.close();

  if (verbosity > 2){
    std::cout <<"---------------------------------------------------------------------"<<std::endl;
    std::cout << "MonitorDAQ: ReadInConfiguration: Read in the following data into configuration variables: "<<std::endl;
    for (unsigned int i_t=0; i_t < config_timeframes.size(); i_t++){
      std::cout <<config_timeframes.at(i_t)<<", "<<config_endtime.at(i_t)<<", "<<config_label.at(i_t)<<", ";
      for (unsigned int i_plot = 0; i_plot < config_plottypes.at(i_t).size(); i_plot++){
        std::cout <<config_plottypes.at(i_t).at(i_plot)<<", ";
      }
      std::cout<<std::endl;
    }
    std::cout <<"-----------------------------------------------------------------------"<<std::endl;
  }

  Log("MonitorDAQ: ReadInConfiguration: Parsing dates: ",v_message,verbosity);
  for (unsigned int i_date = 0; i_date < config_endtime.size(); i_date++){
    if (config_endtime.at(i_date) == "TEND_LASTFILE") {
      Log("MonitorDAQ: TEND_LASTFILE: Starting from end of last read-in file",v_message,verbosity);
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    } else if (config_endtime.at(i_date).size()==15){
        boost::posix_time::ptime spec_endtime(boost::posix_time::from_iso_string(config_endtime.at(i_date)));
        boost::posix_time::time_duration spec_endtime_duration = boost::posix_time::time_duration(spec_endtime - *Epoch);
        ULong64_t spec_endtime_long = spec_endtime_duration.total_milliseconds();
        config_endtime_long.push_back(spec_endtime_long);
    } else {
      Log("MonitorDAQ: Specified end date "+config_endtime.at(i_date)+" does not have the desired format yyyymmddThhmmss. Please change the format in the config file in order to use this tool. Starting from end of last file",v_message,verbosity);
      ULong64_t zero = 0;
      config_endtime_long.push_back(zero);
    }
  }

}

void MonitorDAQ::InitializeHists(){

  Log("MonitorDAQ: Initialize Hists",v_message,verbosity);
  
  //-------------------------------------------------------
  //----------------InitializeHists -----------------------
  //-------------------------------------------------------

  readfromfile_tend = 0;
  readfromfile_timeframe = 0.;

  std::stringstream ss_title_pie, ss_title_evolution_size, ss_title_evolution_vme, ss_title_info;

  ss_title_pie << "VME Services " << title_time.str();
  ss_title_evolution_size << "File sizes " << title_time.str();
  ss_title_evolution_vme << "VME Services " << title_time.str();
  ss_title_info << "DAQ Info " << title_time.str();

  gROOT->cd();

  //VME pie chart
  pie_vme = new TPie("pie_vme",ss_title_pie.str().c_str(),2);	//2 types: running/not running
  pie_vme->GetSlice(0)->SetFillColor(8);
  pie_vme->GetSlice(0)->SetTitle("running");
  pie_vme->GetSlice(1)->SetFillColor(2);
  pie_vme->GetSlice(1)->SetTitle("not running");
  pie_vme->SetCircle(.5,.45,.3);
  pie_vme->GetSlice(1)->SetRadiusOffset(0.05);
  pie_vme->SetLabelFormat("%val");
  pie_vme->SetLabelsOffset(0.02);
  leg_vme = (TLegend*) pie_vme->MakeLegend();
  leg_vme->SetY1(0.78);
  leg_vme->SetY2(0.91);

  //Canvas collection
  canvas_infobox = new TCanvas("canvas_infobox","Infobox",900,600);
  canvas_vmeservice = new TCanvas("canvas_vmeservice","VME services",900,600);
  canvas_timeevolution_size = new TCanvas("canvas_timeevolution_size","Size Time evolution",900,600);
  canvas_timeevolution_vme = new TCanvas("canvas_timeevolution_vme","VME Time evolution",900,600);
  canvas_info_diskspace = new TCanvas("canvas_info_diskspace","DAQ01 Diskspace",900,600);
  canvas_timeevolution_mem = new TCanvas("canvas_timeevolution_mem","Memory Time evolution",900,600);
  canvas_timeevolution_cpu = new TCanvas("canvas_timeevolution_cpu","CPU Time evolution",900,600);
  canvas_timeevolution_disk = new TCanvas("canvas_timeevolution_disk","Diskspace Time evolution",900,600);

  //TGraphs (time evolution plots)
  gr_filesize = new TGraph();
  gr_filesize->SetName("gr_filesize");
  gr_filesize->SetTitle("File sizes");
  if (draw_marker) gr_filesize->SetMarkerStyle(20);
  gr_filesize->SetMarkerColor(1);
  gr_filesize->SetLineColor(1);
  gr_filesize->SetLineWidth(2);
  gr_filesize->SetFillColor(0);
  gr_filesize->GetYaxis()->SetTitle("File size [MB]");
  gr_filesize->GetXaxis()->SetTimeDisplay(1);
  gr_filesize->GetXaxis()->SetLabelSize(0.03);
  gr_filesize->GetXaxis()->SetLabelOffset(0.03);
  gr_filesize->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");

  gr_vmeservice = new TGraph();
  gr_vmeservice->SetName("gr_vmeservice");
  gr_vmeservice->SetTitle("VME services");
  if (draw_marker) gr_vmeservice->SetMarkerStyle(20);
  gr_vmeservice->SetMarkerColor(1);
  gr_vmeservice->SetLineColor(1);
  gr_vmeservice->SetLineWidth(2);
  gr_vmeservice->SetFillColor(0);
  gr_vmeservice->GetYaxis()->SetTitle("# VME services");
  gr_vmeservice->GetXaxis()->SetTimeDisplay(1);
  gr_vmeservice->GetXaxis()->SetLabelSize(0.03);
  gr_vmeservice->GetXaxis()->SetLabelOffset(0.03);
  gr_vmeservice->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");

  gr_mem_daq01 = new TGraph();
  gr_mem_daq01->SetName("gr_mem_daq01");
  gr_mem_daq01->SetTitle("DAQ01 Memory");
  if (draw_marker) gr_mem_daq01->SetMarkerStyle(20);
  gr_mem_daq01->SetMarkerColor(1);
  gr_mem_daq01->SetLineColor(1);
  gr_mem_daq01->SetLineWidth(2);
  gr_mem_daq01->SetFillColor(0);
  gr_mem_daq01->GetYaxis()->SetTitle("Memory [%]");
  gr_mem_daq01->GetXaxis()->SetTimeDisplay(1);
  gr_mem_daq01->GetXaxis()->SetLabelSize(0.03);
  gr_mem_daq01->GetXaxis()->SetLabelOffset(0.03);
  gr_mem_daq01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_mem_vme01 = new TGraph();
  gr_mem_vme01->SetName("gr_mem_vme01");
  gr_mem_vme01->SetTitle("VME01 Memory");
  if (draw_marker) gr_mem_vme01->SetMarkerStyle(20);
  gr_mem_vme01->SetMarkerColor(2);
  gr_mem_vme01->SetLineColor(2);
  gr_mem_vme01->SetLineWidth(2);
  gr_mem_vme01->SetFillColor(0);
  gr_mem_vme01->GetYaxis()->SetTitle("Memory [%]");
  gr_mem_vme01->GetXaxis()->SetTimeDisplay(1);
  gr_mem_vme01->GetXaxis()->SetLabelSize(0.03);
  gr_mem_vme01->GetXaxis()->SetLabelOffset(0.03);
  gr_mem_vme01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_mem_vme02 = new TGraph();
  gr_mem_vme02->SetName("gr_mem_vme02");
  gr_mem_vme02->SetTitle("VME02 Memory");
  if (draw_marker) gr_mem_vme02->SetMarkerStyle(20);
  gr_mem_vme02->SetMarkerColor(4);
  gr_mem_vme02->SetLineColor(4);
  gr_mem_vme02->SetLineWidth(2);
  gr_mem_vme02->SetFillColor(0);
  gr_mem_vme02->GetYaxis()->SetTitle("Memory [%]");
  gr_mem_vme02->GetXaxis()->SetTimeDisplay(1);
  gr_mem_vme02->GetXaxis()->SetLabelSize(0.03);
  gr_mem_vme02->GetXaxis()->SetLabelOffset(0.03);
  gr_mem_vme02->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_mem_vme03 = new TGraph();
  gr_mem_vme03->SetName("gr_mem_vme03");
  gr_mem_vme03->SetTitle("VME03 Memory");
  if (draw_marker) gr_mem_vme03->SetMarkerStyle(20);
  gr_mem_vme03->SetMarkerColor(8);
  gr_mem_vme03->SetLineColor(8);
  gr_mem_vme03->SetLineWidth(2);
  gr_mem_vme03->SetFillColor(0);
  gr_mem_vme03->GetYaxis()->SetTitle("Memory [%]");
  gr_mem_vme03->GetXaxis()->SetTimeDisplay(1);
  gr_mem_vme03->GetXaxis()->SetLabelSize(0.03);
  gr_mem_vme03->GetXaxis()->SetLabelOffset(0.03);
  gr_mem_vme03->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_cpu_daq01 = new TGraph();
  gr_cpu_daq01->SetName("gr_cpu_daq01");
  gr_cpu_daq01->SetTitle("DAQ01 CPU Usage");
  if (draw_marker) gr_cpu_daq01->SetMarkerStyle(20);
  gr_cpu_daq01->SetMarkerColor(1);
  gr_cpu_daq01->SetLineColor(1);
  gr_cpu_daq01->SetLineWidth(2);
  gr_cpu_daq01->SetFillColor(0);
  gr_cpu_daq01->GetYaxis()->SetTitle("CPU [%]");
  gr_cpu_daq01->GetXaxis()->SetTimeDisplay(1);
  gr_cpu_daq01->GetXaxis()->SetLabelSize(0.03);
  gr_cpu_daq01->GetXaxis()->SetLabelOffset(0.03);
  gr_cpu_daq01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_cpu_vme01 = new TGraph();
  gr_cpu_vme01->SetName("gr_cpu_vme01");
  gr_cpu_vme01->SetTitle("VME01 CPU Usage");
  if (draw_marker) gr_cpu_vme01->SetMarkerStyle(20);
  gr_cpu_vme01->SetMarkerColor(2);
  gr_cpu_vme01->SetLineColor(2);
  gr_cpu_vme01->SetLineWidth(2);
  gr_cpu_vme01->SetFillColor(0);
  gr_cpu_vme01->GetYaxis()->SetTitle("CPU [%]");
  gr_cpu_vme01->GetXaxis()->SetTimeDisplay(1);
  gr_cpu_vme01->GetXaxis()->SetLabelSize(0.03);
  gr_cpu_vme01->GetXaxis()->SetLabelOffset(0.03);
  gr_cpu_vme01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_cpu_vme02 = new TGraph();
  gr_cpu_vme02->SetName("gr_cpu_vme02");
  gr_cpu_vme02->SetTitle("VME02 CPU Usage");
  if (draw_marker) gr_cpu_vme02->SetMarkerStyle(20);
  gr_cpu_vme02->SetMarkerColor(4);
  gr_cpu_vme02->SetLineColor(4);
  gr_cpu_vme02->SetLineWidth(2);
  gr_cpu_vme02->SetFillColor(0);
  gr_cpu_vme02->GetYaxis()->SetTitle("CPU [%]");
  gr_cpu_vme02->GetXaxis()->SetTimeDisplay(1);
  gr_cpu_vme02->GetXaxis()->SetLabelSize(0.03);
  gr_cpu_vme02->GetXaxis()->SetLabelOffset(0.03);
  gr_cpu_vme02->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_cpu_vme03 = new TGraph();
  gr_cpu_vme03->SetName("gr_cpu_vme03");
  gr_cpu_vme03->SetTitle("VME03 CPU Usage");
  if (draw_marker) gr_cpu_vme03->SetMarkerStyle(20);
  gr_cpu_vme03->SetMarkerColor(8);
  gr_cpu_vme03->SetLineColor(8);
  gr_cpu_vme03->SetLineWidth(2);
  gr_cpu_vme03->SetFillColor(0);
  gr_cpu_vme03->GetYaxis()->SetTitle("CPU [%]");
  gr_cpu_vme03->GetXaxis()->SetTimeDisplay(1);
  gr_cpu_vme03->GetXaxis()->SetLabelSize(0.03);
  gr_cpu_vme03->GetXaxis()->SetLabelOffset(0.03);
  gr_cpu_vme03->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_disk_daq01 = new TGraph();
  gr_disk_daq01->SetName("gr_disk_daq01");
  gr_disk_daq01->SetTitle("DAQ01 Diskspace");
  if (draw_marker) gr_disk_daq01->SetMarkerStyle(20);
  gr_disk_daq01->SetMarkerColor(1);
  gr_disk_daq01->SetLineColor(1);
  gr_disk_daq01->SetLineWidth(2);
  gr_disk_daq01->SetFillColor(0);
  gr_disk_daq01->GetYaxis()->SetTitle("Disk space [%]");
  gr_disk_daq01->GetXaxis()->SetTimeDisplay(1);
  gr_disk_daq01->GetXaxis()->SetLabelSize(0.03);
  gr_disk_daq01->GetXaxis()->SetLabelOffset(0.03);
  gr_disk_daq01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_disk_vme01 = new TGraph();
  gr_disk_vme01->SetName("gr_disk_vme01");
  gr_disk_vme01->SetTitle("VME01 Diskspace");
  if (draw_marker) gr_disk_vme01->SetMarkerStyle(20);
  gr_disk_vme01->SetMarkerColor(2);
  gr_disk_vme01->SetLineColor(2);
  gr_disk_vme01->SetLineWidth(2);
  gr_disk_vme01->SetFillColor(0);
  gr_disk_vme01->GetYaxis()->SetTitle("Disk space [%]");
  gr_disk_vme01->GetXaxis()->SetTimeDisplay(1);
  gr_disk_vme01->GetXaxis()->SetLabelSize(0.03);
  gr_disk_vme01->GetXaxis()->SetLabelOffset(0.03);
  gr_disk_vme01->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_disk_vme02 = new TGraph();
  gr_disk_vme02->SetName("gr_disk_vme02");
  gr_disk_vme02->SetTitle("VME02 Diskspace");
  if (draw_marker) gr_disk_vme02->SetMarkerStyle(20);
  gr_disk_vme02->SetMarkerColor(4);
  gr_disk_vme02->SetLineColor(4);
  gr_disk_vme02->SetLineWidth(2);
  gr_disk_vme02->SetFillColor(0);
  gr_disk_vme02->GetYaxis()->SetTitle("Disk space [%]");
  gr_disk_vme02->GetXaxis()->SetTimeDisplay(1);
  gr_disk_vme02->GetXaxis()->SetLabelSize(0.03);
  gr_disk_vme02->GetXaxis()->SetLabelOffset(0.03);
  gr_disk_vme02->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  
  gr_disk_vme03 = new TGraph();
  gr_disk_vme03->SetName("gr_disk_vme03");
  gr_disk_vme03->SetTitle("VME03 Diskspace");
  if (draw_marker) gr_disk_vme03->SetMarkerStyle(20);
  gr_disk_vme03->SetMarkerColor(8);
  gr_disk_vme03->SetLineColor(8);
  gr_disk_vme03->SetLineWidth(2);
  gr_disk_vme03->SetFillColor(0);
  gr_disk_vme03->GetYaxis()->SetTitle("Disk space [%]");
  gr_disk_vme03->GetXaxis()->SetTimeDisplay(1);
  gr_disk_vme03->GetXaxis()->SetLabelSize(0.03);
  gr_disk_vme03->GetXaxis()->SetLabelOffset(0.03);
  gr_disk_vme03->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");

  multi_mem = new TMultiGraph();
  multi_cpu = new TMultiGraph();
  multi_disk = new TMultiGraph();
  leg_mem = new TLegend(0.7,0.7,0.88,0.88);
  leg_cpu = new TLegend(0.7,0.7,0.88,0.88);
  leg_disk = new TLegend(0.7,0.7,0.88,0.88);
  leg_mem->SetLineColor(0);
  leg_cpu->SetLineColor(0);
  leg_disk->SetLineColor(0);

  //TText objects for infobox
  text_summary = new TText();
  text_vmeservice = new TText();
  text_filesize = new TText();
  text_filedate = new TText();
  text_haspmt = new TText();
  text_hascc = new TText();
  text_hastrigger = new TText();
  text_currentdate = new TText();
  text_filename = new TText();

  text_disk_title = new TText();
  text_disk_daq01 = new TText();
  text_disk_vme01 = new TText();
  text_disk_vme02 = new TText();
  text_disk_vme03 = new TText();
  text_mem_daq01 = new TText();
  text_mem_vme01 = new TText();
  text_mem_vme02 = new TText();
  text_mem_vme03 = new TText();

  text_summary->SetNDC(1);
  text_vmeservice->SetNDC(1);
  text_filesize->SetNDC(1);
  text_filedate->SetNDC(1);
  text_haspmt->SetNDC(1);
  text_hascc->SetNDC(1);
  text_hastrigger->SetNDC(1);
  text_currentdate->SetNDC(1);
  text_filename->SetNDC(1);

  text_disk_title->SetNDC(1);
  text_disk_daq01->SetNDC(1);
  text_disk_vme01->SetNDC(1);
  text_disk_vme02->SetNDC(1);
  text_disk_vme03->SetNDC(1);
  text_mem_daq01->SetNDC(1);
  text_mem_vme01->SetNDC(1);
  text_mem_vme02->SetNDC(1);
  text_mem_vme03->SetNDC(1);

}

void MonitorDAQ::GetFileInformation(){

  Log("MonitorDAQ: GetFileInformation",v_message,verbosity);

  //-------------------------------------------------------
  //-------------GetFileInformation------------------------
  //-------------------------------------------------------
 
  m_data->CStore.Get("HasTrigData",file_has_trig);
  m_data->CStore.Get("HasCCData",file_has_cc);
  m_data->CStore.Get("HasPMTData",file_has_pmt);
  bool above_hundred = false;
  m_data->CStore.Get("Above100",above_hundred);
  if (above_hundred){
    file_has_trig = true;
    file_has_cc = true;
    file_has_pmt = true;
  }
  m_data->CStore.Get("CurrentFileName",file_name);
  m_data->CStore.Get("CurrentFileSize",file_size_uint);
  file_size = (file_size_uint)/1048576.;	//1MB=1024bytes*1024bytes
  m_data->CStore.Get("CurrentFileTime",file_time);
  file_timestamp = (ULong64_t) file_time*1000;	//cast from time_t to ULong64_t & convert seconds to milliseconds
  file_timestamp -= utc_to_t;	//Correct timestamp to be displayed in Fermilab time


  if (verbosity > 2){
    std::cout <<"////////////////////////////////"<<std::endl;
    std::cout <<"MonitorDAQ: Get File information"<<std::endl;
    std::cout <<"HasTrigData: "<<file_has_trig<<std::endl;
    std::cout <<"HasCCData: "<<file_has_cc<<std::endl;
    std::cout <<"HasPMTData: "<<file_has_pmt<<std::endl;
    std::cout <<"FileName: "<<file_name<<std::endl;
    std::cout <<"FileSize: "<<file_size<<std::endl;
    std::cout <<"FileTime: "<<file_time<<", converted: "<<file_timestamp<<std::endl;
    std::cout <<"//////////////////////////////////"<<std::endl;
  }

  //Get information about the number of VME service processes running
  this->GetVMEServices(online);
  
  //Get statistics about the DAQ & VME computers
  this->GetCompStats();

}

void MonitorDAQ::GetVMEServices(bool is_online){

  if (!is_online){
    num_vme_service=3;
  }
  else {
    num_vme_service=0;

    //What follows is blatantly copied code from Ben's "control" section of the webpage
    std::vector<Store*> RemoteServices;
    
    bool running=true;
    zmq::socket_t Ireceive (*context, ZMQ_DEALER);
    Ireceive.connect("inproc://ServiceDiscovery");

    sleep(7);

    zmq::message_t send(256);
    snprintf ((char *) send.data(), 256 , "%s" ,"All NULL") ;
    Ireceive.send(send);
      
    zmq::message_t receive;
    Ireceive.recv(&receive);
    std::istringstream iss(static_cast<char*>(receive.data()));
    int size;
    iss>>size;
    
    RemoteServices.clear();

    for(int i=0;i<size;i++){

      Store *service = new Store;
      zmq::message_t servicem;
      Ireceive.recv(&servicem);
      std::istringstream ss(static_cast<char*>(servicem.data()));
      service->JsonParser(ss.str());
      RemoteServices.push_back(service);

    }


    for(int i=0;i<(int)RemoteServices.size();i++){

      std::string ip;
      std::string service;
      std::string status;
      std::string colour;

      ip=*((*(RemoteServices.at(i)))["ip"]);
      service=*((*(RemoteServices.at(i)))["msg_value"]);
      status=*((*(RemoteServices.at(i)))["status"]);
      Log("MonitorDAQ: RemoteService "+std::to_string(i)+", service: "+service+", status: "+status,v_debug,verbosity);
      colour="#00FFFF";
      if (status=="Online")colour="#FF00FF";
      else if (status=="Waiting to Initialise ToolChain")colour="#FFFF00";
      else{
        std::stringstream tmpstatus(status);
        tmpstatus>>status;
        if(status=="ToolChain"){
	  tmpstatus>>status;
	  if(status=="running"){
            colour="#00FF00";
            if (service=="VME_service") num_vme_service++;
          }
        }
        status=tmpstatus.str();
      }
    }

    //Some cleanup
    for (int i=0;i<(int)RemoteServices.size();i++){
      delete RemoteServices.at(i);
    }
   
  }
  Log("MonitorDAQ tool: GetVMEServices: Got "+std::to_string(num_vme_service)+" VME services!",v_message,verbosity);
  if (num_vme_service < 3){
    Log("MonitorDAQ tool ERROR: Only "+std::to_string(num_vme_service)+" VME services running! Check DAQ",v_error,verbosity);
  }

}

void MonitorDAQ::GetCompStats(){

  //Get statistics about the DAQ computer (disk space, memory consumption)
  std::stringstream ss_path_daq01;
  ss_path_daq01 << path_compstats << "/compstatsDAQ01";
  if (verbosity > 4) std::cout <<"MonitorDAQ: GetCompStats: Opening file "<<ss_path_daq01.str()<<std::endl;
  Store storeDAQ01;
  storeDAQ01.Initialise(ss_path_daq01.str().c_str());
  long timeDAQ01;
  int diskDAQ01;
  long memtotDAQ01;
  long memfreeDAQ01;
  double idleDAQ01;
  storeDAQ01.Get("Time",timeDAQ01);
  storeDAQ01.Get("Disk",diskDAQ01);
  storeDAQ01.Get("MemTotal",memtotDAQ01);
  storeDAQ01.Get("MemFree",memfreeDAQ01);
  storeDAQ01.Get("CPUidle",idleDAQ01);

  if (verbosity > 4){
    std::cout <<"///////////// Read in DAQ01 file information /////////"<<std::endl;
    std::cout <<"Filepath: "<<ss_path_daq01.str()<<std::endl;
    std::cout <<"Time: "<<timeDAQ01<<std::endl;
    std::cout <<"MemTotal: "<<memtotDAQ01<<std::endl;
    std::cout <<"MemFree: "<<memfreeDAQ01<<std::endl;
    std::cout <<"CPUidle: "<<idleDAQ01<<std::endl;
  }

  timestamp_daq01 = (ULong64_t) timeDAQ01*1000;	//convert to ms
  disk_daq01 = diskDAQ01;
  mem_daq01 = double(memtotDAQ01-memfreeDAQ01)/memtotDAQ01;
  cpu_daq01 = 100.-idleDAQ01;
  timestamp_daq01 -= utc_to_t;   //Correct timestamp to be displayed in Fermilab time
  //Sanity checks
  if (std::isinf(mem_daq01) || std::isnan(mem_daq01)) mem_daq01 = 1.;
  else if (mem_daq01 < 0.) mem_daq01 = 0.;
  else if (mem_daq01 > 1.) mem_daq01 = 1.;
  if (std::isinf(cpu_daq01) || std::isnan(cpu_daq01)) cpu_daq01 = 100.;
  else if (cpu_daq01 < 0.) cpu_daq01 = 0.;
  else if (cpu_daq01 > 100.) cpu_daq01 = 100.;

  std::stringstream ss_path_vme01;
  ss_path_vme01 << path_compstats << "/compstatsVME01";
  if (verbosity > 4) std::cout <<"MonitorDAQ: GetCompStats: Opening file "<<ss_path_vme01.str()<<std::endl;
  Store storeVME01;
  storeVME01.Initialise(ss_path_vme01.str().c_str());
  long timeVME01;
  long memtotVME01;
  long memfreeVME01;
  long memavailVME01;
  double idleVME01;
  int diskVME01;
  storeVME01.Get("Time",timeVME01);
  storeVME01.Get("MemTotal",memtotVME01);
  storeVME01.Get("MemFree",memfreeVME01);
  storeVME01.Get("MemAvailable",memavailVME01);
  storeVME01.Get("CPUidle",idleVME01);
  storeVME01.Get("Disk",diskVME01);
  timestamp_vme01 = (ULong64_t) timeVME01*1000; //convert to ms
  disk_vme01 = diskVME01;
  mem_vme01 = double(memtotVME01-memfreeVME01)/memtotVME01;
  cpu_vme01 = 100.-idleVME01;
  timestamp_vme01 -= utc_to_t;   //Correct timestamp to be displayed in Fermilab time
  //Sanity checks
  if (std::isinf(mem_vme01) || std::isnan(mem_vme01)) mem_vme01 = 1.;
  else if (mem_vme01 < 0.) mem_vme01 = 0.;
  else if (mem_vme01 > 1.) mem_vme01 = 1.;
  if (std::isinf(cpu_vme01) || std::isnan(cpu_vme01)) cpu_vme01 = 100.;
  else if (cpu_vme01 < 0.) cpu_vme01 = 0.;
  else if (cpu_vme01 > 100.) cpu_vme01 = 100.;

  std::stringstream ss_path_vme02;
  ss_path_vme02 << path_compstats << "/compstatsVME02";
  if (verbosity > 4) std::cout <<"MonitorDAQ: GetCompStats: Opening file "<<ss_path_vme02.str()<<std::endl;
  Store storeVME02;
  storeVME02.Initialise(ss_path_vme02.str().c_str());
  long timeVME02;
  long memtotVME02;
  long memfreeVME02;
  long memavailVME02;
  double idleVME02;
  int diskVME02;
  storeVME02.Get("Time",timeVME02);
  storeVME02.Get("MemTotal",memtotVME02);
  storeVME02.Get("MemFree",memfreeVME02);
  storeVME02.Get("MemAvailable",memavailVME02);
  storeVME02.Get("CPUidle",idleVME02);
  storeVME02.Get("Disk",diskVME02);
  timestamp_vme02 = (ULong64_t) timeVME02*1000;  //convert to ms
  disk_vme02 = diskVME02;
  mem_vme02 = double(memtotVME02-memfreeVME02)/memtotVME02;
  cpu_vme02 = 100.-idleVME02;
  timestamp_vme02 -= utc_to_t;   //Correct timestamp to be displayed in Fermilab time
  //Sanity checks
  if (std::isinf(mem_vme02) || std::isnan(mem_vme02)) mem_vme02 = 1.;
  else if (mem_vme02 < 0.) mem_vme02 = 0.;
  else if (mem_vme02 > 1.) mem_vme02 = 1.;
  if (std::isinf(cpu_vme02) || std::isnan(cpu_vme02)) cpu_vme02 = 100.;
  else if (cpu_vme02 < 0.) cpu_vme02 = 0.;
  else if (cpu_vme02 > 100.) cpu_vme02 = 100.;

  std::stringstream ss_path_vme03;
  ss_path_vme03 << path_compstats << "/compstatsVME03";
  if (verbosity > 4) std::cout <<"MonitorDAQ: GetCompStats: Opening file "<<ss_path_vme03.str()<<std::endl;
  Store storeVME03;
  storeVME03.Initialise(ss_path_vme03.str().c_str());
  long timeVME03;
  long memtotVME03;
  long memfreeVME03;
  long memavailVME03;
  double idleVME03;
  int diskVME03;
  storeVME03.Get("Time",timeVME03);
  storeVME03.Get("MemTotal",memtotVME03);
  storeVME03.Get("MemFree",memfreeVME03);
  storeVME03.Get("MemAvailable",memavailVME03);
  storeVME03.Get("CPUidle",idleVME03);
  storeVME03.Get("Disk",diskVME03);
  timestamp_vme03 = (ULong64_t) timeVME03*1000;  //convert to ns
  disk_vme03 = diskVME03;
  mem_vme03 = double(memtotVME03-memfreeVME03)/memtotVME03;
  cpu_vme03 = 100.-idleVME03;
  timestamp_vme03 -= utc_to_t;   //Correct timestamp to be displayed in Fermilab time
  //Sanity checks
  if (std::isinf(mem_vme03) || std::isnan(mem_vme03)) mem_vme03 = 1.;
  else if (mem_vme03 < 0.) mem_vme03 = 0.;
  else if (mem_vme03 > 1.) mem_vme03 = 1.;
  if (std::isinf(cpu_vme03) || std::isnan(cpu_vme03)) cpu_vme03 = 100.;
  else if (cpu_vme03 < 0.) cpu_vme03 = 0.;
  else if (cpu_vme03 > 100.) cpu_vme03 = 100.;

}

void MonitorDAQ::WriteToFile(){

  Log("MonitorDAQ: WriteToFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------WriteToFile -------------------------
  //-------------------------------------------------------

  //Not really possible to distinguish between t_file_start, and t_file_end, set both to the same value
  t_file_start = file_timestamp;
  t_file_end = file_timestamp;
  std::string file_start_date = convertTimeStamp_to_Date(t_file_start);
  std::stringstream root_filename;
  root_filename << path_monitoring << "DAQ_" << file_start_date <<".root";

  Log("MonitorDAQ: ROOT filename: "+root_filename.str(),v_message,verbosity);
  
  std::string root_option = "RECREATE";
  if (does_file_exist(root_filename.str())) root_option = "UPDATE";
  TFile *f = new TFile(root_filename.str().c_str(),root_option.c_str());

  ULong64_t t_start, t_end;
  std::string *file_name_pointer = new std::string;
  bool temp_file_has_trig, temp_file_has_cc, temp_file_has_pmt;
  double temp_file_size;
  ULong64_t temp_file_time;
  int temp_num_vme_service;
  double temp_disk_daq01;
  double temp_disk_vme01;
  double temp_disk_vme02;
  double temp_disk_vme03;
  double temp_mem_daq01;
  double temp_mem_vme01;
  double temp_mem_vme02;
  double temp_mem_vme03;
  double temp_cpu_daq01;
  double temp_cpu_vme01;
  double temp_cpu_vme02;
  double temp_cpu_vme03;
  ULong64_t temp_stamp_daq01;
  ULong64_t temp_stamp_vme01;
  ULong64_t temp_stamp_vme02;
  ULong64_t temp_stamp_vme03;

  TTree *t;
  if (f->GetListOfKeys()->Contains("daqmonitor_tree")) {
    Log("MonitorDAQ: WriteToFile: Tree already exists",v_message,verbosity);
    t = (TTree*) f->Get("daqmonitor_tree");
    t->SetBranchAddress("t_start",&t_start);
    t->SetBranchAddress("t_end",&t_end);
    t->SetBranchAddress("has_trig",&temp_file_has_trig);
    t->SetBranchAddress("has_cc",&temp_file_has_cc);
    t->SetBranchAddress("has_pmt",&temp_file_has_pmt);
    t->SetBranchAddress("file_name",&file_name_pointer);
    t->SetBranchAddress("file_size",&temp_file_size);
    t->SetBranchAddress("file_time",&temp_file_time);
    t->SetBranchAddress("num_vme_services",&temp_num_vme_service);
    t->SetBranchAddress("timestamp_daq01",&temp_stamp_daq01);
    t->SetBranchAddress("timestamp_vme01",&temp_stamp_vme01);
    t->SetBranchAddress("timestamp_vme02",&temp_stamp_vme02);
    t->SetBranchAddress("timestamp_vme03",&temp_stamp_vme03);
    t->SetBranchAddress("disk_daq01",&temp_disk_daq01);
    t->SetBranchAddress("disk_vme01",&temp_disk_vme01);
    t->SetBranchAddress("disk_vme02",&temp_disk_vme02);
    t->SetBranchAddress("disk_vme03",&temp_disk_vme03);
    t->SetBranchAddress("mem_daq01",&temp_mem_daq01);
    t->SetBranchAddress("mem_vme01",&temp_mem_vme01);
    t->SetBranchAddress("mem_vme02",&temp_mem_vme02);
    t->SetBranchAddress("mem_vme03",&temp_mem_vme03);
    t->SetBranchAddress("cpu_daq01",&temp_cpu_daq01);
    t->SetBranchAddress("cpu_vme01",&temp_cpu_vme01);
    t->SetBranchAddress("cpu_vme02",&temp_cpu_vme02);
    t->SetBranchAddress("cpu_vme03",&temp_cpu_vme03);
  } else {
    t = new TTree("daqmonitor_tree","DAQ Monitoring tree");
    Log("MonitorDAQ: WriteToFile: Tree is created from scratch",v_message,verbosity);
    t->Branch("t_start",&t_start);
    t->Branch("t_end",&t_end);
    t->Branch("has_trig",&temp_file_has_trig);
    t->Branch("has_cc",&temp_file_has_cc);
    t->Branch("has_pmt",&temp_file_has_pmt);
    t->Branch("file_name",&file_name_pointer);
    t->Branch("file_size",&temp_file_size);
    t->Branch("file_time",&temp_file_time);
    t->Branch("num_vme_services",&temp_num_vme_service);
    t->Branch("timestamp_daq01",&temp_stamp_daq01);
    t->Branch("timestamp_vme01",&temp_stamp_vme01);
    t->Branch("timestamp_vme02",&temp_stamp_vme02);
    t->Branch("timestamp_vme03",&temp_stamp_vme03);
    t->Branch("disk_daq01",&temp_disk_daq01);
    t->Branch("disk_vme01",&temp_disk_vme01);
    t->Branch("disk_vme02",&temp_disk_vme02);
    t->Branch("disk_vme03",&temp_disk_vme03);
    t->Branch("mem_daq01",&temp_mem_daq01);
    t->Branch("mem_vme01",&temp_mem_vme01);
    t->Branch("mem_vme02",&temp_mem_vme02);
    t->Branch("mem_vme03",&temp_mem_vme03);
    t->Branch("cpu_daq01",&temp_cpu_daq01);
    t->Branch("cpu_vme01",&temp_cpu_vme01);
    t->Branch("cpu_vme02",&temp_cpu_vme02);
    t->Branch("cpu_vme03",&temp_cpu_vme03);
  }

  int n_entries = t->GetEntries();
  bool omit_entries = false;
  for (int i_entry = 0; i_entry < n_entries; i_entry++){
    t->GetEntry(i_entry);
    if ((long) t_start == t_file_start) {
      Log("WARNING (MonitorDAQ): WriteToFile: Wanted to write data from file that is already written to DB. Omit entries.",v_warning,verbosity);
      omit_entries = true;
    }
  }


  //If data is already written to DB/File, do not write it again
  if (omit_entries){
    //Don't write file again, but still delete TFile and TTree object!
    f->Close();
    delete file_name_pointer;
    delete f;

    gROOT->cd();

    return;

  }

  t_start = t_file_start;
  t_end = t_file_end;
  *file_name_pointer = file_name;
  temp_file_has_trig = file_has_trig;
  temp_file_has_cc = file_has_cc;
  temp_file_has_pmt = file_has_pmt;
  temp_file_size = file_size;
  temp_file_time = file_timestamp;
  temp_num_vme_service = num_vme_service;
  temp_stamp_daq01 = (ULong64_t) timestamp_daq01;
  temp_stamp_vme01 = (ULong64_t) timestamp_vme01;
  temp_stamp_vme02 = (ULong64_t) timestamp_vme02;
  temp_stamp_vme03 = (ULong64_t) timestamp_vme03;
  temp_disk_daq01 = disk_daq01;
  temp_disk_vme01 = disk_vme01;
  temp_disk_vme02 = disk_vme02;
  temp_disk_vme03 = disk_vme03;
  temp_mem_daq01 = mem_daq01;
  temp_mem_vme01 = mem_vme01;
  temp_mem_vme02 = mem_vme02;
  temp_mem_vme03 = mem_vme03;
  temp_cpu_daq01 = cpu_daq01;
  temp_cpu_vme01 = cpu_vme01;
  temp_cpu_vme02 = cpu_vme02;
  temp_cpu_vme03 = cpu_vme03;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(t_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_start/MSEC_to_SEC/SEC_to_MIN)%60,int(t_start/MSEC_to_SEC/1000.)%60,t_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  Log("MonitorDAQ: WriteToFile: Writing data to file: "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
    +"..."+std::to_string(endtime_tm.tm_year+1900)+"/"+std::to_string(endtime_tm.tm_mon+1)+"/"+std::to_string(endtime_tm.tm_mday)+"-"+std::to_string(endtime_tm.tm_hour)+":"+std::to_string(endtime_tm.tm_min)+":"+std::to_string(endtime_tm.tm_sec),v_message,verbosity);

  t->Fill();
  t->Write("",TObject::kOverwrite);
  f->Close();

  delete file_name_pointer;
  delete f;

  gROOT->cd();

}

void MonitorDAQ::ReadFromFile(ULong64_t timestamp_end, double time_frame){

  Log("MonitorDAQ: ReadFromFile",v_message,verbosity);

  //-------------------------------------------------------
  //------------------ReadFromFile ------------------------
  //-------------------------------------------------------

  has_trig_plot.clear();
  has_cc_plot.clear();
  has_pmt_plot.clear();
  filename_plot.clear();
  filetimestamp_plot.clear();
  filesize_plot.clear();
  num_vme_plot.clear();
  labels_timeaxis.clear();
  tstart_plot.clear();
  tend_plot.clear();
  disk_daq01_plot.clear();
  disk_vme01_plot.clear();
  disk_vme02_plot.clear();
  disk_vme03_plot.clear();
  mem_daq01_plot.clear();
  mem_vme01_plot.clear();
  mem_vme02_plot.clear();
  mem_vme03_plot.clear();
  cpu_daq01_plot.clear();
  cpu_vme01_plot.clear();
  cpu_vme02_plot.clear();
  cpu_vme03_plot.clear();
  t_daq01_plot.clear();
  t_vme01_plot.clear();
  t_vme02_plot.clear();
  t_vme03_plot.clear();

  //take the end time and calculate the start time with the given time_frame
  ULong64_t timestamp_start = timestamp_end - time_frame*MIN_to_HOUR*SEC_to_MIN*MSEC_to_SEC;

  boost::posix_time::ptime starttime = *Epoch + boost::posix_time::time_duration(int(timestamp_start/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_start/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_start/MSEC_to_SEC/1000.)%60,timestamp_start%1000);
  struct tm starttime_tm = boost::posix_time::to_tm(starttime);
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);

  Log("MonitorDAQTime: ReadFromFile: Reading in data for time frame "+std::to_string(starttime_tm.tm_year+1900)+"/"+std::to_string(starttime_tm.tm_mon+1)+"/"+std::to_string(starttime_tm.tm_mday)+"-"+std::to_string(starttime_tm.tm_hour)+":"+std::to_string(starttime_tm.tm_min)+":"+std::to_string(starttime_tm.tm_sec)
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
    root_filename_i << path_monitoring << "DAQ_" << string_date_i <<".root";
    bool tree_exists = true;  

    if (does_file_exist(root_filename_i.str())) {
      TFile *f = new TFile(root_filename_i.str().c_str(),"READ");
      TTree *t;

	std::string *temp_filename = new std::string;
      if (f->GetListOfKeys()->Contains("daqmonitor_tree")) t = (TTree*) f->Get("daqmonitor_tree");
      else {
        Log("WARNING (MonitorDAQ): File "+root_filename_i.str()+" does not contain daqmonitor_tree. Omit file.",v_warning,verbosity);
        tree_exists = false;
      }
 
      if (tree_exists){

        Log("MonitorDAQ: Tree exists, start reading in data",v_message,verbosity);
  
        ULong64_t t_start, t_end;
	bool has_trig;
	bool has_cc;
	bool has_pmt;
	double temp_filesize;
	ULong64_t temp_filetime;
	int temp_numvme;
  
        int nevents;
        int nentries_tree;

        ULong64_t temp_timestamp_daq01;
        ULong64_t temp_timestamp_vme01;
        ULong64_t temp_timestamp_vme02;
        ULong64_t temp_timestamp_vme03;
        double temp_disk_daq01;
        double temp_disk_vme01;
        double temp_disk_vme02;
        double temp_disk_vme03;
        double temp_mem_daq01;
        double temp_mem_vme01;
        double temp_mem_vme02;
        double temp_mem_vme03;
        double temp_cpu_daq01;
        double temp_cpu_vme01;
        double temp_cpu_vme02;
        double temp_cpu_vme03;
  
        t->SetBranchAddress("t_start",&t_start);
        t->SetBranchAddress("t_end",&t_end);
	t->SetBranchAddress("has_trig",&has_trig);
	t->SetBranchAddress("has_cc",&has_cc);
	t->SetBranchAddress("has_pmt",&has_pmt);
        t->SetBranchAddress("file_name",&temp_filename);
        t->SetBranchAddress("file_size",&temp_filesize);
        t->SetBranchAddress("file_time",&temp_filetime);
	t->SetBranchAddress("num_vme_services",&temp_numvme);
        t->SetBranchAddress("timestamp_daq01",&temp_timestamp_daq01);
        t->SetBranchAddress("timestamp_vme01",&temp_timestamp_vme01);
        t->SetBranchAddress("timestamp_vme02",&temp_timestamp_vme02);
        t->SetBranchAddress("timestamp_vme03",&temp_timestamp_vme03);
        t->SetBranchAddress("disk_daq01",&temp_disk_daq01);
        t->SetBranchAddress("disk_vme01",&temp_disk_vme01);
        t->SetBranchAddress("disk_vme02",&temp_disk_vme02);
        t->SetBranchAddress("disk_vme03",&temp_disk_vme03);
        t->SetBranchAddress("mem_daq01",&temp_mem_daq01);
        t->SetBranchAddress("mem_vme01",&temp_mem_vme01);
        t->SetBranchAddress("mem_vme02",&temp_mem_vme02);
        t->SetBranchAddress("mem_vme03",&temp_mem_vme03);
        t->SetBranchAddress("cpu_daq01",&temp_cpu_daq01);
        t->SetBranchAddress("cpu_vme01",&temp_cpu_vme01);
        t->SetBranchAddress("cpu_vme02",&temp_cpu_vme02);
        t->SetBranchAddress("cpu_vme03",&temp_cpu_vme03);  

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
	    has_trig_plot.push_back(has_trig);
	    has_cc_plot.push_back(has_cc);
	    has_pmt_plot.push_back(has_pmt);
	    filename_plot.push_back(*temp_filename);
	    filesize_plot.push_back(temp_filesize);
	    filetimestamp_plot.push_back(temp_filetime);
	    num_vme_plot.push_back(temp_numvme);
            tstart_plot.push_back(t_start);
            tend_plot.push_back(t_end);
            disk_daq01_plot.push_back(temp_disk_daq01);
            disk_vme01_plot.push_back(temp_disk_vme01);
            disk_vme02_plot.push_back(temp_disk_vme02);
            disk_vme03_plot.push_back(temp_disk_vme03);
            mem_daq01_plot.push_back(temp_mem_daq01);
            mem_vme01_plot.push_back(temp_mem_vme01);
            mem_vme02_plot.push_back(temp_mem_vme02);
            mem_vme03_plot.push_back(temp_mem_vme03);
            cpu_daq01_plot.push_back(temp_cpu_daq01);
            cpu_vme01_plot.push_back(temp_cpu_vme01);
            cpu_vme02_plot.push_back(temp_cpu_vme02);
            cpu_vme03_plot.push_back(temp_cpu_vme03);
            t_daq01_plot.push_back(temp_timestamp_daq01);
            t_vme01_plot.push_back(temp_timestamp_vme01);
            t_vme02_plot.push_back(temp_timestamp_vme02);
            t_vme03_plot.push_back(temp_timestamp_vme03);
            boost::posix_time::ptime boost_tend = *Epoch+boost::posix_time::time_duration(int(t_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(t_end/MSEC_to_SEC/SEC_to_MIN)%60,int(t_end/MSEC_to_SEC/1000.)%60,t_end%1000);
            struct tm label_timestamp = boost::posix_time::to_tm(boost_tend);
  
            TDatime datime_timestamp(1900+label_timestamp.tm_year,label_timestamp.tm_mon+1,label_timestamp.tm_mday,label_timestamp.tm_hour,label_timestamp.tm_min,label_timestamp.tm_sec);
            labels_timeaxis.push_back(datime_timestamp);
          }
        }
      }
  
      f->Close();
      delete temp_filename;
      delete f;
      gROOT->cd();

    } else {
      Log("MonitorDAQ: ReadFromFile: File "+root_filename_i.str()+" does not exist. Omit file.",v_warning,verbosity);
    }

  }

  //Set the readfromfile time variables to make sure data is not read twice for the same time window
  readfromfile_tend = timestamp_end;
  readfromfile_timeframe = time_frame;

}

void MonitorDAQ::UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes){

  Log("MonitorDAQ: UpdateMonitorPlots",v_message,verbosity);

  //-------------------------------------------------------
  //------------------UpdateMonitorPlots ------------------
  //-------------------------------------------------------

  //Draw the monitoring plots according to the specifications in the configfiles
  
  for (unsigned int i_time = 0; i_time < timeFrames.size(); i_time++){
    
    ULong64_t zero = 0;
    if (endTimes.at(i_time) == zero) endTimes.at(i_time) = t_file_end;     //set 0 for t_file_end since we did not know what that was at the beginning of initialise

    for (unsigned int i_plot = 0; i_plot < plotTypes.at(i_time).size(); i_plot++){

      if (plotTypes.at(i_time).at(i_plot) == "DAQTimeEvolution") this->DrawDAQTimeEvolution(endTimes.at(i_time),timeFrames.at(i_time),fileLabels.at(i_time));
      else {
        if (verbosity > 0) std::cout <<"ERROR (MonitorDAQ): UpdateMonitorPlots: Specified plot type -"<<plotTypes.at(i_time).at(i_plot)<<"- does not exist! Omit entry."<<std::endl;
      }
    }

  }

}

void MonitorDAQ::DrawDAQTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending){

  Log("MonitorDAQ: DrawDAQTimeEvolution",v_message,verbosity);

  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  //-------------------------------------------------------
  //-------------DrawDAQTimeEvolution ------------------------
  //-------------------------------------------------------

  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);

  //looping over all files that are in the time interval, each file will be one data point

  std::stringstream ss_timeframe;
  ss_timeframe << round(time_frame*100.)/100.;

  //Resetting time evolution graphs
  gr_filesize->Set(0);
  gr_vmeservice->Set(0);
  gr_mem_daq01->Set(0);
  gr_mem_vme01->Set(0);
  gr_mem_vme02->Set(0);
  gr_mem_vme03->Set(0);
  gr_cpu_daq01->Set(0);
  gr_cpu_vme01->Set(0);
  gr_cpu_vme02->Set(0);
  gr_cpu_vme03->Set(0);
  gr_disk_daq01->Set(0);
  gr_disk_vme01->Set(0);
  gr_disk_vme02->Set(0);
  gr_disk_vme03->Set(0);


  for (int i_file=0; i_file < (int) tend_plot.size(); i_file++){
    gr_filesize->SetPoint(i_file,labels_timeaxis[i_file].Convert(),filesize_plot.at(i_file));
    gr_vmeservice->SetPoint(i_file,labels_timeaxis[i_file].Convert(),num_vme_plot.at(i_file));
    gr_mem_daq01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),mem_daq01_plot.at(i_file)*100.);
    gr_mem_vme01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),mem_vme01_plot.at(i_file)*100.);
    gr_mem_vme02->SetPoint(i_file,labels_timeaxis[i_file].Convert(),mem_vme02_plot.at(i_file)*100.);
    gr_mem_vme03->SetPoint(i_file,labels_timeaxis[i_file].Convert(),mem_vme03_plot.at(i_file)*100.);
    gr_cpu_daq01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),cpu_daq01_plot.at(i_file));
    gr_cpu_vme01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),cpu_vme01_plot.at(i_file));
    gr_cpu_vme02->SetPoint(i_file,labels_timeaxis[i_file].Convert(),cpu_vme02_plot.at(i_file));
    gr_cpu_vme03->SetPoint(i_file,labels_timeaxis[i_file].Convert(),cpu_vme03_plot.at(i_file));
    gr_disk_daq01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),disk_daq01_plot.at(i_file));
    gr_disk_vme01->SetPoint(i_file,labels_timeaxis[i_file].Convert(),disk_vme01_plot.at(i_file));
    gr_disk_vme02->SetPoint(i_file,labels_timeaxis[i_file].Convert(),disk_vme02_plot.at(i_file));
    gr_disk_vme03->SetPoint(i_file,labels_timeaxis[i_file].Convert(),disk_vme03_plot.at(i_file));
  }

  std::stringstream ss_title_filesize, ss_title_vme, ss_title_memory, ss_title_cpu, ss_title_disk;
  ss_title_filesize << "File sizes (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  ss_title_vme << "VME services (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  ss_title_memory << "Memory consumption (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  ss_title_cpu << "CPU load (last "<<ss_timeframe.str()<<"h) "<<end_time.str();
  ss_title_disk << "Disk space (last "<<ss_timeframe.str()<<"h) "<<end_time.str();

  gr_filesize->SetTitle(ss_title_filesize.str().c_str());
  gr_vmeservice->SetTitle(ss_title_vme.str().c_str());
  gr_mem_daq01->SetTitle(ss_title_memory.str().c_str());
  gr_cpu_daq01->SetTitle(ss_title_cpu.str().c_str());
  gr_disk_daq01->SetTitle(ss_title_disk.str().c_str());

  std::stringstream ss_filename_filesize, ss_filename_vme, ss_filename_mem, ss_filename_cpu, ss_filename_disk;
  ss_filename_filesize << outpath << "DAQFileSize_"<<file_ending<<"."<<img_extension;
  ss_filename_vme << outpath << "DAQVMEServices_"<<file_ending<<"."<<img_extension;
  ss_filename_mem << outpath << "DAQMemory_"<<file_ending<<"."<<img_extension;
  ss_filename_cpu << outpath << "DAQCPU_"<<file_ending<<"."<<img_extension;
  ss_filename_disk << outpath << "DAQDisk_"<<file_ending<<"."<<img_extension;

  canvas_timeevolution_size->Clear();
  canvas_timeevolution_size->cd();
  gr_filesize->Draw("apl");
  gr_filesize->GetYaxis()->SetTitle("File size [MB]");
  gr_filesize->GetXaxis()->SetLabelSize(0.03);
  gr_filesize->GetXaxis()->SetLabelOffset(0.03);
  gr_filesize->GetXaxis()->SetTimeDisplay(1);
  gr_filesize->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  gr_filesize->GetXaxis()->SetTimeOffset(0.);
  double max_filesize = TMath::MaxElement(tend_plot.size(),gr_filesize->GetY());
  gr_filesize->GetYaxis()->SetRangeUser(0.001,1.1*max_filesize);
  canvas_timeevolution_size->SaveAs(ss_filename_filesize.str().c_str());

  canvas_timeevolution_vme->Clear();
  canvas_timeevolution_vme->cd();
  gr_vmeservice->Draw("apl");
  gr_vmeservice->GetYaxis()->SetTitle("# VME Services");
  gr_vmeservice->GetXaxis()->SetLabelSize(0.03);
  gr_vmeservice->GetXaxis()->SetLabelOffset(0.03);
  gr_vmeservice->GetXaxis()->SetTimeDisplay(1);
  gr_vmeservice->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  gr_vmeservice->GetXaxis()->SetTimeOffset(0.);
  double max_vme = TMath::MaxElement(tend_plot.size(),gr_vmeservice->GetY());
  gr_vmeservice->GetYaxis()->SetRangeUser(0.001,1.1*max_vme);
  canvas_timeevolution_vme->SaveAs(ss_filename_vme.str().c_str());

  canvas_timeevolution_mem->Clear();
  canvas_timeevolution_mem->cd();
  multi_mem->Add(gr_mem_daq01);
  leg_mem->AddEntry(gr_mem_daq01,"DAQ01","l");
  multi_mem->Add(gr_mem_vme01);
  leg_mem->AddEntry(gr_mem_vme01,"VME01","l");
  multi_mem->Add(gr_mem_vme02);
  leg_mem->AddEntry(gr_mem_vme02,"VME02","l");
  multi_mem->Add(gr_mem_vme03);
  leg_mem->AddEntry(gr_mem_vme03,"VME03","l");
  multi_mem->Draw("apl");
  multi_mem->SetTitle(ss_title_memory.str().c_str());
  multi_mem->GetYaxis()->SetTitle("Memory [%]");
  multi_mem->GetXaxis()->SetTimeDisplay(1);
  multi_mem->GetXaxis()->SetLabelSize(0.03);
  multi_mem->GetXaxis()->SetLabelOffset(0.03);
  multi_mem->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  multi_mem->GetXaxis()->SetTimeOffset(0.);
  leg_mem->Draw();
  canvas_timeevolution_mem->SaveAs(ss_filename_mem.str().c_str());
  multi_mem->RecursiveRemove(gr_mem_daq01);
  multi_mem->RecursiveRemove(gr_mem_vme01);
  multi_mem->RecursiveRemove(gr_mem_vme02);
  multi_mem->RecursiveRemove(gr_mem_vme03);
  leg_mem->Clear();

  canvas_timeevolution_cpu->Clear();
  canvas_timeevolution_cpu->cd();
  multi_cpu->Add(gr_cpu_daq01);
  leg_cpu->AddEntry(gr_cpu_daq01,"DAQ01","l");
  multi_cpu->Add(gr_cpu_vme01);
  leg_cpu->AddEntry(gr_cpu_vme01,"VME01","l");
  multi_cpu->Add(gr_cpu_vme02);
  leg_cpu->AddEntry(gr_cpu_vme02,"VME02","l");
  multi_cpu->Add(gr_cpu_vme03);
  leg_cpu->AddEntry(gr_cpu_vme03,"VME03","l");
  multi_cpu->Draw("apl");
  multi_cpu->SetTitle(ss_title_cpu.str().c_str());
  multi_cpu->GetYaxis()->SetTitle("CPU [%]");
  multi_cpu->GetXaxis()->SetTimeDisplay(1);
  multi_cpu->GetXaxis()->SetLabelSize(0.03);
  multi_cpu->GetXaxis()->SetLabelOffset(0.03);
  multi_cpu->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  multi_cpu->GetXaxis()->SetTimeOffset(0.);
  leg_cpu->Draw();
  canvas_timeevolution_cpu->SaveAs(ss_filename_cpu.str().c_str());
  multi_cpu->RecursiveRemove(gr_cpu_daq01);
  multi_cpu->RecursiveRemove(gr_cpu_vme01);
  multi_cpu->RecursiveRemove(gr_cpu_vme02);
  multi_cpu->RecursiveRemove(gr_cpu_vme03);
  leg_cpu->Clear();
  
  canvas_timeevolution_disk->Clear();
  canvas_timeevolution_disk->cd();
  multi_disk->Add(gr_disk_daq01);
  leg_disk->AddEntry(gr_disk_daq01,"DAQ01","l");
  multi_disk->Add(gr_disk_vme01);
  leg_disk->AddEntry(gr_disk_vme01,"VME01","l");
  multi_disk->Add(gr_disk_vme02);
  leg_disk->AddEntry(gr_disk_vme02,"VME02","l");
  multi_disk->Add(gr_disk_vme03);
  leg_disk->AddEntry(gr_disk_vme03,"VME03","l");
  multi_disk->Draw("apl");
  multi_disk->SetTitle(ss_title_disk.str().c_str());
  multi_disk->GetYaxis()->SetTitle("Disk space [%]");
  multi_disk->GetXaxis()->SetTimeDisplay(1);
  multi_disk->GetXaxis()->SetLabelSize(0.03);
  multi_disk->GetXaxis()->SetLabelOffset(0.03);
  multi_disk->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  multi_disk->GetXaxis()->SetTimeOffset(0.);
  leg_disk->Draw();
  canvas_timeevolution_disk->SaveAs(ss_filename_disk.str().c_str());
  multi_disk->RecursiveRemove(gr_disk_daq01);
  multi_disk->RecursiveRemove(gr_disk_vme01);
  multi_disk->RecursiveRemove(gr_disk_vme02);
  multi_disk->RecursiveRemove(gr_disk_vme03);
  leg_disk->Clear();
}

void MonitorDAQ::DrawVMEService(ULong64_t timestamp_end, double time_frame, std::string file_ending, bool current){

  Log("MonitorDAQ tool: DrawVMEService",v_message,verbosity);
  //-------------------------------------------------------
  //---------------------DrawVMEService--------------------
  //-------------------------------------------------------
  
  if (timestamp_end != readfromfile_tend || time_frame != readfromfile_timeframe) ReadFromFile(timestamp_end, time_frame);
  
  boost::posix_time::ptime endtime = *Epoch + boost::posix_time::time_duration(int(timestamp_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_end/MSEC_to_SEC/1000.)%60,timestamp_end%1000);
  struct tm endtime_tm = boost::posix_time::to_tm(endtime);
  std::stringstream end_time;
  end_time << "Current time: "<<endtime_tm.tm_year+1900<<"/"<<endtime_tm.tm_mon+1<<"/"<<endtime_tm.tm_mday<<"-"<<endtime_tm.tm_hour<<":"<<endtime_tm.tm_min<<":"<<endtime_tm.tm_sec;

  int current_num_vme_service;
  if (!current) current_num_vme_service = num_vme_plot.at(num_vme_plot.size()-1);
  else current_num_vme_service = num_vme_service;

//  int current_num_vme_service = num_vme_service;
  canvas_vmeservice->Clear();
  canvas_vmeservice->cd();

  pie_vme->GetSlice(0)->SetValue(current_num_vme_service);
  pie_vme->GetSlice(1)->SetValue(3-current_num_vme_service);

  std::stringstream ss_title_pie, ss_path_pie;
  ss_title_pie << "# VME services "<< end_time.str();
  ss_path_pie << outpath <<"DAQ_VMEServices_"<<file_ending<<"."<<img_extension;
  pie_vme->SetTitle(ss_title_pie.str().c_str());
  pie_vme->Draw("tsc");
  leg_vme->Draw();
  canvas_vmeservice->SaveAs(ss_path_pie.str().c_str());

}

void MonitorDAQ::PrintInfoBox(){

  Log("MonitorDAQ tool: PrintInfoBox",v_message,verbosity);
  
  double mins_since_file = 0.;
  
  if (file_produced){
  std::stringstream ss_trigger;
  if (file_has_trig) {
    ss_trigger << "TriggerData: True";
    text_hastrigger->SetText(0.06,0.1,ss_trigger.str().c_str());
    text_hastrigger->SetTextColor(1);	//black text if ok
  } else {
    ss_trigger << "TriggerData: False";
    text_hastrigger->SetText(0.06,0.1,ss_trigger.str().c_str());
    text_hastrigger->SetTextColor(2);	//red text if not ok
  }
 
  std::stringstream ss_cc;
  if (file_has_cc) {
    ss_cc << "MRDData: True";
    text_hascc->SetText(0.06,0.2,ss_cc.str().c_str());
    text_hascc->SetTextColor(1);	//black text if ok
  } else {
    ss_cc << "MRDData: False";
    text_hascc->SetText(0.06,0.2,ss_cc.str().c_str());
    text_hascc->SetTextColor(2);	//red text if not ok
  }

  std::stringstream ss_pmt;
  if (file_has_pmt) {
    ss_pmt << "VMEData: True";
    text_haspmt->SetText(0.06,0.3,ss_pmt.str().c_str());
    text_haspmt->SetTextColor(1);	//black text if ok
  } else {
    ss_pmt << "VMEData: False";
    text_haspmt->SetText(0.06,0.3,ss_pmt.str().c_str());
    text_haspmt->SetTextColor(2);	//red text if not ok
  }
  boost::posix_time::ptime currenttime = *Epoch + boost::posix_time::time_duration(int(current_stamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(current_stamp/MSEC_to_SEC/SEC_to_MIN)%60,int(current_stamp/MSEC_to_SEC)%60,current_stamp%1000);
  struct tm currenttime_tm = boost::posix_time::to_tm(currenttime);
  std::stringstream current_time;
  current_time << currenttime_tm.tm_year+1900<<"/"<<currenttime_tm.tm_mon+1<<"/"<<currenttime_tm.tm_mday<<"-"<<currenttime_tm.tm_hour<<":"<<currenttime_tm.tm_min<<":"<<currenttime_tm.tm_sec;
  
  std::stringstream ss_currentdate;
  ss_currentdate << "Current time: "<<current_time.str();
  text_currentdate->SetText(0.06,0.7,ss_currentdate.str().c_str());
  text_currentdate->SetTextColor(1);

  boost::posix_time::ptime lasttime = *Epoch + boost::posix_time::time_duration(int(file_timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(file_timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(file_timestamp/MSEC_to_SEC)%60,file_timestamp%1000);
  struct tm lasttime_tm = boost::posix_time::to_tm(lasttime);
  std::stringstream ss_last_time;
  ss_last_time << lasttime_tm.tm_year+1900<<"/"<<lasttime_tm.tm_mon+1<<"/"<<lasttime_tm.tm_mday<<"-"<<lasttime_tm.tm_hour<<":"<<lasttime_tm.tm_min<<":"<<lasttime_tm.tm_sec;

  boost::posix_time::time_duration time_duration_since_file(currenttime-lasttime);
  long time_since_file=time_duration_since_file.total_milliseconds();
  mins_since_file = (time_since_file/MSEC_to_SEC/SEC_to_MIN);
  mins_since_file = round(mins_since_file*100.)/100.;

  std::stringstream ss_filedate;
  ss_filedate << "Last File time: "<<ss_last_time.str() << " ("<<mins_since_file<<" min ago)";
  text_filedate->SetText(0.06,0.6,ss_filedate.str().c_str());
  if (mins_since_file <= 30.) text_filedate->SetTextColor(1);
  else text_filedate->SetTextColor(2);
  
  std::stringstream ss_filename;
  std::size_t pos = file_name.find("RAWData"); 
  std::string file_name_short = file_name.substr(pos);
  ss_filename << "Last File name: "<<file_name_short;
  text_filename->SetText(0.06,0.5,ss_filename.str().c_str());

  std::stringstream ss_filesize;
  ss_filesize << "Last File size: "<<file_size<<" MB";
  text_filesize->SetText(0.06,0.4,ss_filesize.str().c_str());
  if (file_size<100.) text_filesize->SetTextColor(2);
  else text_filesize->SetTextColor(1);
	  
  if (testmode){
    if (testcounter < (int) test_filesize.size()) file_size = test_filesize.at(testcounter);
    else file_size = test_filesize.at(test_filesize.size()-1);
  }

  if (file_size <= 100.) {
    std::stringstream ss_error_filesize, ss_error_filesize_slack;
    ss_error_filesize << "ERROR (MonitorDAQ tool): Very small filesize < 100 MB for file " << file_name_short << ": Size = " << std::to_string(file_size) << " MB.";
    ss_error_filesize_slack << "payload={\"text\":\"Monitoring: Very small filesize < 100 MB for file " << file_name_short << ": Size = " << std::to_string(file_size) << " MB.\"}";
    bool issue_warning = (!warning_filesize || (warning_filesize_filename!=file_name_short));
    warning_filesize = true;
    warning_filesize_filename = file_name_short;
    if (issue_warning) Log(ss_error_filesize.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_filesize_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_filesize = false;
  }
  
  if (testmode){
    if (testcounter < (int) test_hastrig.size()) file_has_trig = test_hastrig.at(testcounter);
    else file_has_trig = test_hastrig.at(test_hastrig.size()-1);
  }

  if (!file_has_trig) {
    std::stringstream ss_error_trig, ss_error_trig_slack;
    ss_error_trig << "ERROR (MonitorDAQ tool): Did not find Trigger data in last file (" << file_name_short << ")";
    ss_error_trig_slack << "payload={\"text\":\" Monitoring: Did not find Trigger data in last file (" << file_name_short << ")\"}";
    bool issue_warning = (!warning_trigdata || (warning_trigdata_filename!=file_name_short));
    warning_trigdata = true;
    warning_trigdata_filename = file_name_short;
    if (issue_warning) Log(ss_error_trig.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_trig_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_trigdata = false;
  }
  
  if (testmode){
    if (testcounter < (int) test_haspmt.size()) file_has_pmt = test_haspmt.at(testcounter);
    else file_has_pmt = test_haspmt.at(test_haspmt.size()-1);
  }
  
  if (!file_has_pmt) {
    std::stringstream ss_error_pmt, ss_error_pmt_slack;
    ss_error_pmt << "ERROR (MonitorDAQ tool): Did not find VME data in last file (" << file_name_short << ")";
    ss_error_pmt_slack << "payload={\"text\":\"Monitoring: Did not find VME data in last file (" << file_name_short << ")\"}";
    bool issue_warning = (!warning_pmtdata || (warning_pmtdata_filename!=file_name_short));
    warning_pmtdata = true;
    warning_pmtdata_filename = file_name_short;
    if (issue_warning) Log(ss_error_pmt.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_pmt_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_pmtdata = false;
  }	  
  } // end if file_produced

  std::stringstream ss_vmeservice;
  ss_vmeservice << "VME_services: "<<num_vme_service<<" / 3";
  text_vmeservice->SetText(0.06,0.8,ss_vmeservice.str().c_str());
  if (num_vme_service<3) text_vmeservice->SetTextColor(2);
  else text_vmeservice->SetTextColor(1);
  
  if (testmode){
    if (testcounter < (int) test_vme.size()) num_vme_service = test_vme.at(testcounter);
    else num_vme_service = test_vme.at(test_vme.size()-1);
  }

  if (num_vme_service < 3) {
    std::stringstream ss_error_vme, ss_error_vme_slack;
    ss_error_vme << "ERROR (MonitorDAQ tool): Did not find 3 running VME services! Current # of VME_service processes: " << num_vme_service;
    ss_error_vme_slack << "payload={\"text\":\"Monitoring: Did not find 3 running VME services! Current # of VME_service processes: " << num_vme_service << "\"}";
    bool issue_warning = (!warning_vme || (num_vme_service != warning_vme_num));
    warning_vme_num=num_vme_service;
    warning_vme = true;
    if (issue_warning) Log(ss_error_vme.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_vme_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_vme = false;
    warning_vme_num=num_vme_service;
  }


  bool everything_ok=false;
  if (file_produced){
    if (file_has_trig && file_has_pmt && file_size > 100. && num_vme_service==3 && mins_since_file<30.) everything_ok = true;
  }
  else if (num_vme_service==3) everything_ok = true;
  std::stringstream ss_summary;
  if (everything_ok){
    ss_summary << "ANNIE DAQ functioning within normal parameters.";
    text_summary->SetText(0.06,0.9,ss_summary.str().c_str());
    text_summary->SetTextColor(8);
  } else {
    ss_summary << "Something seems to be wrong.";
    text_summary->SetText(0.06,0.9,ss_summary.str().c_str());
    text_summary->SetTextColor(2);
  }

  text_disk_title->SetText(0.06,0.9,"Disk space / Memory");
 
  if (testmode) timestamp_daq01 += (0.3*60*1000*testcounter);

  boost::posix_time::ptime daq01time = *Epoch + boost::posix_time::time_duration(int(timestamp_daq01/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_daq01/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_daq01/MSEC_to_SEC)%60,timestamp_daq01%1000);
  struct tm daq01time_tm = boost::posix_time::to_tm(daq01time);
  std::stringstream ss_daq01_time;
  ss_daq01_time << daq01time_tm.tm_year+1900<<"/"<<daq01time_tm.tm_mon+1<<"/"<<daq01time_tm.tm_mday<<"-"<<daq01time_tm.tm_hour<<":"<<daq01time_tm.tm_min<<":"<<daq01time_tm.tm_sec;

  std::stringstream ss_text_disk_daq01;
  ss_text_disk_daq01 << "DAQ01 Disk space: "<<disk_daq01<<" % ("<<ss_daq01_time.str()<<")";
  text_disk_daq01->SetText(0.06,0.8,ss_text_disk_daq01.str().c_str());
  text_disk_daq01->SetTextColor(1);	//default color
  if (disk_daq01 >= 80.) text_disk_daq01->SetTextColor(kOrange);
  if (disk_daq01 >= 90.) text_disk_daq01->SetTextColor(kRed);

  if (testmode){
    if (testcounter < (int) test_disk.size()) disk_daq01 = test_disk.at(testcounter);
    else disk_daq01 = test_disk.at(test_disk.size()-1);
  }

  if (disk_daq01 >= 80.) {
    std::stringstream ss_error_disk, ss_error_disk_slack;
    ss_error_disk << "ERROR (MonitorDAQ tool): DAQ01 disk space above 80%! (" << disk_daq01 << " %)";
    ss_error_disk_slack << "payload={\"text\":\"Monitoring: DAQ01 disk space above 80%! (" << disk_daq01 << " %)\"}";
    bool issue_warning = (!warning_diskspace_80);
    warning_diskspace_80 = true;
    if (issue_warning) Log(ss_error_disk.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_disk_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_diskspace_80 = false;
  }
  if (disk_daq01 >= 85.) {
    std::stringstream ss_error_disk, ss_error_disk_slack;
    ss_error_disk << "ERROR (MonitorDAQ tool): DAQ01 disk space above 85%! (" << disk_daq01 << " %)";
    ss_error_disk_slack << "payload={\"text\":\"Monitoring: DAQ01 disk space above 85%! (" << disk_daq01 << " %)\"}";
    bool issue_warning = (!warning_diskspace_85);
    warning_diskspace_85 = true;
    if (issue_warning) Log(ss_error_disk.str().c_str(),v_error,verbosity);
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_disk_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_diskspace_85 = false;
  }
  if (disk_daq01 >= 90.) {
    std::stringstream ss_error_disk, ss_error_disk_slack;
    ss_error_disk << "ERROR (MonitorDAQ tool): DAQ01 disk space above 90%! (" << disk_daq01 << " %)";
    ss_error_disk_slack << "payload={\"text\":\"Monitoring: Serious warning! DAQ01 disk space above 90%! (" << disk_daq01 << " %)\"}";
    double t_since_last_warning = 0;
    if (!warning_diskspace_90) timestamp_last_warning_diskspace_90 = daq01time; 
    boost::posix_time::time_duration time_duration_since_warning(daq01time-timestamp_last_warning_diskspace_90);
    long time_since_warning=time_duration_since_warning.total_milliseconds();
    t_since_last_warning = (time_since_warning/MSEC_to_SEC/SEC_to_MIN);
    bool issue_warning = (!warning_diskspace_90 || (t_since_last_warning > 1440));
    warning_diskspace_90 = true;
    if (issue_warning) {
      Log(ss_error_disk.str().c_str(),v_error,verbosity);
      timestamp_last_warning_diskspace_90 = daq01time;  
    }
    if (send_slack && issue_warning){
    try{
      CURL *curl;
      CURLcode res;
      curl_global_init(CURL_GLOBAL_ALL);
      curl=curl_easy_init();
      if (curl){
        curl_easy_setopt(curl,CURLOPT_URL,hook.c_str());
        std::string field = ss_error_disk_slack.str();
        curl_easy_setopt(curl,CURLOPT_POSTFIELDS,field.c_str());
        res=curl_easy_perform(curl);
        if (res != CURLE_OK) Log("MonitorDAQ tool: curl_easy_perform() failed.",v_error,verbosity);
        curl_easy_cleanup(curl);
      }
      curl_global_cleanup();
    }
    catch(...){
      Log("MonitorDAQ tool: Slack send an error",v_warning,verbosity);
    }
    }
  } else {
    warning_diskspace_90 = false;
  }

  std::stringstream ss_text_mem_daq01;
  double mem_round = round(mem_daq01*100.)/100.;
  double cpu_round = round(cpu_daq01*100.)/100.;
  ss_text_mem_daq01 << "DAQ01 Mem: "<<mem_round*100<<" %"<<"  CPU: "<<cpu_round<<" % ("<<ss_daq01_time.str()<<")";
  text_mem_daq01->SetText(0.06,0.4,ss_text_mem_daq01.str().c_str());
  text_mem_daq01->SetTextColor(1);	//default color 
  if (mem_round*100 >= 80.) text_mem_daq01->SetTextColor(kOrange);
  if (mem_round*100 >= 90.) text_mem_daq01->SetTextColor(kRed);
  if (cpu_round >= 80.) text_mem_daq01->SetTextColor(kOrange);
  if (cpu_round >= 90.) text_mem_daq01->SetTextColor(kRed);
 
  boost::posix_time::ptime vme01time = *Epoch + boost::posix_time::time_duration(int(timestamp_vme01/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_vme01/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_vme01/MSEC_to_SEC)%60,timestamp_vme01%1000);
  struct tm vme01time_tm = boost::posix_time::to_tm(vme01time);
  std::stringstream ss_vme01_time;
  ss_vme01_time << vme01time_tm.tm_year+1900<<"/"<<vme01time_tm.tm_mon+1<<"/"<<vme01time_tm.tm_mday<<"-"<<vme01time_tm.tm_hour<<":"<<vme01time_tm.tm_min<<":"<<vme01time_tm.tm_sec;
 
  std::stringstream ss_text_disk_vme01;
  ss_text_disk_vme01 << "VME01 Disk space: "<<disk_vme01<<" % ("<<ss_vme01_time.str()<<")";
  text_disk_vme01->SetText(0.06,0.7,ss_text_disk_vme01.str().c_str());
  text_disk_vme01->SetTextColor(1);     //default color
  if (disk_vme01 >= 80.) text_disk_vme01->SetTextColor(kOrange);
  if (disk_vme01 >= 90.) text_disk_vme01->SetTextColor(kRed);

  std::stringstream ss_text_mem_vme01;  
  mem_round = round(mem_vme01*100.)/100.;
  cpu_round = round(cpu_vme01*100.)/100.;
  ss_text_mem_vme01 << "VME01 Mem: "<<mem_round*100<<" %"<<"  CPU: "<<cpu_round<<" % ("<<ss_vme01_time.str()<<")";
  text_mem_vme01->SetText(0.06,0.3,ss_text_mem_vme01.str().c_str());
  text_mem_vme01->SetTextColor(1);	//default color
  if (mem_round*100 >= 80.) text_mem_vme01->SetTextColor(kOrange);
  if (mem_round*100 >= 90.) text_mem_vme01->SetTextColor(kRed);
  if (cpu_round >= 80.) text_mem_vme01->SetTextColor(kOrange);
  if (cpu_round >= 90.) text_mem_vme01->SetTextColor(kRed);
  
  boost::posix_time::ptime vme02time = *Epoch + boost::posix_time::time_duration(int(timestamp_vme02/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_vme02/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_vme02/MSEC_to_SEC)%60,timestamp_vme02%1000);
  struct tm vme02time_tm = boost::posix_time::to_tm(vme02time);
  std::stringstream ss_vme02_time;
  ss_vme02_time << vme02time_tm.tm_year+1900<<"/"<<vme02time_tm.tm_mon+1<<"/"<<vme02time_tm.tm_mday<<"-"<<vme02time_tm.tm_hour<<":"<<vme02time_tm.tm_min<<":"<<vme02time_tm.tm_sec;

  std::stringstream ss_text_disk_vme02;
  ss_text_disk_vme02 << "VME02 Disk space: "<<disk_vme02<<" % ("<<ss_vme02_time.str()<<")";
  text_disk_vme02->SetText(0.06,0.6,ss_text_disk_vme02.str().c_str());
  text_disk_vme02->SetTextColor(1);     //default color
  if (disk_vme02 >= 80.) text_disk_vme02->SetTextColor(kOrange);
  if (disk_vme02 >= 90.) text_disk_vme02->SetTextColor(kRed);
  
  std::stringstream ss_text_mem_vme02;
  mem_round = round(mem_vme02*100.)/100.;
  cpu_round = round(cpu_vme02*100.)/100.;
  ss_text_mem_vme02 << "VME02 Mem: "<<mem_round*100<<" %"<<"  CPU: "<<cpu_round<<" % ("<<ss_vme02_time.str()<<")";
  text_mem_vme02->SetText(0.06,0.2,ss_text_mem_vme02.str().c_str());
  text_mem_vme02->SetTextColor(1);	//default color
  if (mem_round*100 >= 80.) text_mem_vme02->SetTextColor(kOrange);
  if (mem_round*100 >= 90.) text_mem_vme02->SetTextColor(kRed);
  if (cpu_round >= 80.) text_mem_vme02->SetTextColor(kOrange);
  if (cpu_round >= 90.) text_mem_vme02->SetTextColor(kRed);
  
  boost::posix_time::ptime vme03time = *Epoch + boost::posix_time::time_duration(int(timestamp_vme03/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_vme03/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_vme03/MSEC_to_SEC)%60,timestamp_vme03%1000);
  struct tm vme03time_tm = boost::posix_time::to_tm(vme03time);
  std::stringstream ss_vme03_time;
  ss_vme03_time << vme03time_tm.tm_year+1900<<"/"<<vme03time_tm.tm_mon+1<<"/"<<vme03time_tm.tm_mday<<"-"<<vme03time_tm.tm_hour<<":"<<vme03time_tm.tm_min<<":"<<vme03time_tm.tm_sec;

  std::stringstream ss_text_disk_vme03;
  ss_text_disk_vme03 << "VME03 Disk space: "<<disk_vme03<<" % ("<<ss_vme03_time.str()<<")";
  text_disk_vme03->SetText(0.06,0.5,ss_text_disk_vme03.str().c_str());
  text_disk_vme03->SetTextColor(1);     //default color
  if (disk_vme03 >= 80.) text_disk_vme03->SetTextColor(kOrange);
  if (disk_vme03 >= 90.) text_disk_vme03->SetTextColor(kRed);

  std::stringstream ss_text_mem_vme03;
  mem_round = round(mem_vme03*100.)/100.;
  cpu_round = round(cpu_vme03*100.)/100.;
  ss_text_mem_vme03 << "VME03 Mem: "<<mem_round*100<<" %"<<"  CPU: "<<cpu_round<<" % ("<<ss_vme03_time.str()<<")";
  text_mem_vme03->SetText(0.06,0.1,ss_text_mem_vme03.str().c_str());
  text_mem_vme03->SetTextColor(1);	//default color
  if (mem_round*100 >= 80.) text_mem_vme03->SetTextColor(kOrange);
  if (mem_round*100 >= 90.) text_mem_vme03->SetTextColor(kRed);
  if (cpu_round >= 80.) text_mem_vme03->SetTextColor(kOrange);
  if (cpu_round >= 90.) text_mem_vme03->SetTextColor(kRed);
  
  text_summary->SetTextSize(0.05);
  text_vmeservice->SetTextSize(0.05);
  text_currentdate->SetTextSize(0.05);
  text_filesize->SetTextSize(0.05);
  text_filedate->SetTextSize(0.05);
  text_filename->SetTextSize(0.05);
  text_haspmt->SetTextSize(0.05);
  text_hascc->SetTextSize(0.05);
  text_hastrigger->SetTextSize(0.05);
  text_disk_title->SetTextSize(0.05);
  text_disk_daq01->SetTextSize(0.05);
  text_disk_vme01->SetTextSize(0.05);
  text_disk_vme02->SetTextSize(0.05);
  text_disk_vme03->SetTextSize(0.05);
  text_mem_daq01->SetTextSize(0.05);
  text_mem_vme01->SetTextSize(0.05);
  text_mem_vme02->SetTextSize(0.05);
  text_mem_vme03->SetTextSize(0.05);

  text_summary->SetNDC(1);
  text_vmeservice->SetNDC(1);
  text_filesize->SetNDC(1);
  text_filename->SetNDC(1);
  text_currentdate->SetNDC(1);
  text_filedate->SetNDC(1);
  text_haspmt->SetNDC(1);
  text_hascc->SetNDC(1);
  text_hastrigger->SetNDC(1);
  text_disk_title->SetNDC(1);
  text_disk_daq01->SetNDC(1);
  text_disk_vme01->SetNDC(1);
  text_disk_vme02->SetNDC(1);
  text_disk_vme03->SetNDC(1);
  text_mem_daq01->SetNDC(1);
  text_mem_vme01->SetNDC(1);
  text_mem_vme02->SetNDC(1);
  text_mem_vme03->SetNDC(1);

  canvas_infobox->cd();
  canvas_infobox->Clear();
  text_summary->Draw();
  text_vmeservice->Draw();
  text_currentdate->Draw();
  if (file_produced){ 
    text_filedate->Draw();
    text_filename->Draw();
    text_filesize->Draw();
    text_haspmt->Draw();
    text_hascc->Draw();
    text_hastrigger->Draw();
  }

  std::stringstream ss_path_textinfo;
  ss_path_textinfo << outpath <<"DAQInfo_current."<<img_extension;
  canvas_infobox->SaveAs(ss_path_textinfo.str().c_str());
  canvas_infobox->Clear();

  canvas_info_diskspace->cd();
  canvas_info_diskspace->Clear();
  text_disk_title->Draw();
  text_disk_daq01->Draw();
  text_disk_vme01->Draw();
  text_disk_vme02->Draw();
  text_disk_vme03->Draw();
  text_mem_daq01->Draw();
  text_mem_vme01->Draw();
  text_mem_vme02->Draw();
  text_mem_vme03->Draw();
  std::stringstream ss_path_textinfo_disk;
  ss_path_textinfo_disk << outpath << "DAQInfo_current_Diskspace."<<img_extension;
  canvas_info_diskspace->SaveAs(ss_path_textinfo_disk.str().c_str());
  canvas_info_diskspace->Clear();
}

std::string MonitorDAQ::convertTimeStamp_to_Date(ULong64_t timestamp){
  
  //format of date is YYYY_MM-DD

  boost::posix_time::ptime convertedtime = *Epoch + boost::posix_time::time_duration(int(timestamp/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp/MSEC_to_SEC/1000.)%60,timestamp%1000);
  struct tm converted_tm = boost::posix_time::to_tm(convertedtime);
  std::stringstream ss_date;
  ss_date << converted_tm.tm_year+1900 << "_" << converted_tm.tm_mon+1 << "-" <<converted_tm.tm_mday;
  return ss_date.str();

}


bool MonitorDAQ::does_file_exist(std::string filename){

  std::ifstream infile(filename.c_str());
  bool file_good = infile.good();
  infile.close();
  return file_good;

}

