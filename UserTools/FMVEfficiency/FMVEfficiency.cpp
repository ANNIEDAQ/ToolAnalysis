#include "FMVEfficiency.h"

FMVEfficiency::FMVEfficiency():Tool(){}


bool FMVEfficiency::Initialise(std::string configfile, DataModel &data){

  if(configfile!="") m_variables.Initialise(configfile); // loading config file

  //Default variable values
  outputfile = "fmv_efficiency";
  m_data= &data;

  padPerLayer = 1;

  //Get user configuration  
  m_variables.Get("SinglePEGains",singlePEgains);
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("OutputFile",outputfile);
  m_variables.Get("UseTank",useTank);
  m_variables.Get("IsData",isData);
  m_variables.Get("PadPerLayer",padPerLayer);

  if (useTank != 0 && useTank != 1) useTank = 0;
  if (padPerLayer < 1) padPerLayer = 1;

  //Geometry
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();
  n_veto_pmts = geom->GetNumDetectorsInSet("Veto");

  //PMT gains
  ifstream file_singlepe(singlePEgains.c_str());
  unsigned long temp_chankey;
  double temp_gain;
  while (!file_singlepe.eof()){
    file_singlepe >> temp_chankey >> temp_gain;
    if (file_singlepe.eof()) break;
    pmt_gains.emplace(temp_chankey,temp_gain);
  }
  file_singlepe.close();

  //Define ROOT file & histograms
  std::stringstream root_outputfile;
  root_outputfile << outputfile << ".root";
  file = new TFile(root_outputfile.str().c_str(),"RECREATE");

  //MRD-FMV coincidence
  time_diff_Layer1 = new TH1F("time_diff_Layer1","Time differences FMV-MRD (Layer 1)",2000,-4000,4000);
  time_diff_Layer2 = new TH1F("time_diff_Layer2","Time differences FMV-MRD (Layer 2)",2000,-4000,4000);
  num_paddles_Layer1= new TH1F("num_paddles_Layer1","Number of hit FMV paddles (Layer 1)",14,0,14); 
  num_paddles_Layer2= new TH1F("num_paddles_Layer2","Number of hit FMV paddles (Layer 2)",14,0,14); 
  time_diff_paddle = new TH2F("time_diff_paddle","Time difference paddles",26,0,26,200,-2000,2000); 
  time_diff_paddle_coinc = new TH2F("time_diff_paddle_coinc","Time difference paddles (coincidence)",26,0,26,400,-4000,4000); 
  time_diff_paddle_simplecoinc = new TH2F("time_diff_paddle_simplecoinc","Time difference paddles (simple coincidence)",26,0,26,400,-4000,4000); 

  fmv_layer1 = new TH1F("fmv_layer1","FMV Layer 1",13,0,13);
  fmv_layer2 = new TH1F("fmv_layer2","FMV Layer 2",13,13,26);
  fmv_layer1_layer2 = new TH2F("fmv_layer1_layer2","FMV Layer 1 vs. Layer 2",13,0,13,13,13,26);
 
  fmv_observed_layer1 = new TH1F("fmv_observed_layer1","FMV observed hits (Layer 1)",13,0,13);
  fmv_expected_layer1 = new TH1F("fmv_expected_layer1","FMV expected hits (Layer 1)",13,0,13);
  fmv_observed_layer2 = new TH1F("fmv_observed_layer2","FMV observed hits (Layer 2)",13,13,26);
  fmv_expected_layer2 = new TH1F("fmv_expected_layer2","FMV expected hits (Layer 2)",13,13,26);
  fmv_observed_track_strict_layer1 = new TH1F("fmv_observed_track_strict_layer1","FMV observed hits Layer 1 (Strict MRD track)",13,0,13);
  fmv_expected_track_strict_layer1 = new TH1F("fmv_expected_track_strict_layer1","FMV expected hits Layer 1 (Strict MRD track)",13,0,13);
  fmv_observed_track_strict_layer2 = new TH1F("fmv_observed_track_strict_layer2","FMV observed hits Layer 2 (Strict MRD track)",13,13,26);
  fmv_expected_track_strict_layer2 = new TH1F("fmv_expected_track_strict_layer2","FMV expected hits Layer 2 (Strict MRD track)",13,13,26);
  fmv_observed_track_loose_layer1 = new TH1F("fmv_observed_track_loose_layer1","FMV observed hits Layer 1 (MRD track)",13,0,13);
  fmv_expected_track_loose_layer1 = new TH1F("fmv_expected_track_loose_layer1","FMV expected hits Layer 1 (MRD track)",13,0,13);
  fmv_observed_track_loose_layer2 = new TH1F("fmv_observed_track_loose_layer2","FMV observed hits Layer 2 (MRD track)",13,13,26);
  fmv_expected_track_loose_layer2 = new TH1F("fmv_expected_track_loose_layer2","FMV expected hits Layer 2 (MRD track)",13,13,26);
  track_diff_x_Layer1 = new TH1F("track_diff_x_Layer1","Track difference X (Layer 1)",1000,-5.,5.);
  track_diff_y_Layer1 = new TH1F("track_diff_y_Layer1","Track difference Y (Layer 1)",1000,-5.,5.);
  track_diff_xy_Layer1 = new TH2F("track_diff_xy_Layer1","Track difference X/Y (Layer 1)",200,-5.,5.,200,-5.,5);
  track_diff_x_strict_Layer1 = new TH1F("track_diff_x_strict_Layer1","Track difference X Layer 1 (Strict MRD Track)",1000,-5.,5.);
  track_diff_y_strict_Layer1 = new TH1F("track_diff_y_strict_Layer1","Track difference Y Layer 1 (Strict MRD Track)",1000,-5.,5.);
  track_diff_xy_strict_Layer1 = new TH2F("track_diff_xy_strict_Layer1","Track difference X/Y Layer 1 (Strict MRD Track)",200,-5.,5.,200,-5.,5);
  track_diff_x_loose_Layer1 = new TH1F("track_diff_x_loose_Layer1","Track difference X Layer 1 (MRD Track)",1000,-5.,5.);
  track_diff_y_loose_Layer1 = new TH1F("track_diff_y_loose_Layer1","Track difference Y Layer 1 (MRD Track)",1000,-5.,5.);
  track_diff_xy_loose_Layer1 = new TH2F("track_diff_xy_loose_Layer1","Track difference X/Y Layer 1 (MRD Track)",200,-5.,5.,200,-5.,5);
  track_diff_x_Layer2 = new TH1F("track_diff_x_Layer2","Track difference X (Layer 2)",1000,-5.,5.);
  track_diff_y_Layer2 = new TH1F("track_diff_y_Layer2","Track difference Y (Layer 2)",1000,-5.,5.);
  track_diff_xy_Layer2 = new TH2F("track_diff_xy_Layer2","Track difference X/Y (Layer 2)",200,-5.,5.,200,-5.,5);
  track_diff_x_strict_Layer2 = new TH1F("track_diff_x_strict_Layer2","Track difference X Layer 2 (Strict MRD Track)",1000,-5.,5.);
  track_diff_y_strict_Layer2 = new TH1F("track_diff_y_strict_Layer2","Track difference Y Layer 2 (Strict MRD Track)",1000,-5.,5.);
  track_diff_xy_strict_Layer2 = new TH2F("track_diff_xy_strict_Layer2","Track difference X/Y Layer 2 (Strict MRD Track)",200,-5.,5.,200,-5.,5);
  track_diff_x_loose_Layer2 = new TH1F("track_diff_x_loose_Layer2","Track difference X Layer 2 (MRD Track)",1000,-5.,5.);
  track_diff_y_loose_Layer2 = new TH1F("track_diff_y_loose_Layer2","Track difference Y Layer 2 (MRD Track)",1000,-5.,5.);
  track_diff_xy_loose_Layer2 = new TH2F("track_diff_xy_loose_Layer2","Track difference X/Y Layer 2 (MRD Track)",200,-5.,5.,200,-5.,5);

  //Tank-FMV coincidence
  if (useTank){
    time_diff_tank_Layer1 = new TH1F("time_diff_tank_Layer1","Time differences FMV-Tank (Layer 1)",2000,-4000,4000);
    time_diff_tank_Layer2 = new TH1F("time_diff_tank_Layer2","Time differences FMV-Tank (Layer 2)",2000,-4000,4000);

    fmv_tank_observed_layer1 = new TH1F("fmv_tank_observed_layer1","FMV observed hits Tank Coincidence (Layer 1)",13,0,13);
    fmv_tank_expected_layer1 = new TH1F("fmv_tank_expected_layer1","FMV expected hits Tank Coincidence (Layer 1)",13,0,13);
    fmv_tank_observed_layer2 = new TH1F("fmv_tank_observed_layer2","FMV observed hits Tank Coincidence (Layer 2)",13,13,26);
    fmv_tank_expected_layer2 = new TH1F("fmv_tank_expected_layer2","FMV expected hits Tank Coincidence (Layer 2)",13,13,26);
  }

  //Read in MRD information from geometry
  
  bool is_first_fmv_chankey = true;

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Veto").begin();
                                                    it != Detectors->at("Veto").end();
                                                  ++it){
    //Get geometric properties of FMV paddle
    Detector* amrdpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
    if (is_first_fmv_chankey){
      first_fmv_chankey = chankey;
      first_fmv_detkey = detkey;
      is_first_fmv_chankey=false;
    }
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

    if (amrdpmt->GetDetectorElement()=="Veto") Log("FMVEfficiency tool: Reading in Veto paddle with Chankey: "+std::to_string(chankey)+", x = "+std::to_string(xmin)+"-"+std::to_string(xmax)+", y = "+std::to_string(ymin)+"-"+std::to_string(ymax)+", z = "+std::to_string(zmin)+"-"+std::to_string(zmax),v_debug,verbosity);

    //Define histogram for FMV paddle
    //Count chankeys from 0 to 26
    chankey = chankey - first_fmv_chankey;
    std::stringstream ss_observed_strict, ss_expected_strict, ss_observed_strict_title, ss_expected_strict_title, ss_observed_loose, ss_expected_loose, ss_observed_loose_title, ss_expected_loose_title;
    ss_observed_strict << "hist_observed_strict_chankey" << chankey; 
    ss_expected_strict << "hist_expected_strict_chankey" << chankey; 
    ss_observed_strict_title << "Observed hits (Strict MRD Track) chankey " << chankey; 
    ss_expected_strict_title << "Expected hits (Strict MRD Track) chankey " << chankey; 
    ss_observed_loose << "hist_observed_loose_chankey" << chankey; 
    ss_expected_loose << "hist_expected_loose_chankey" << chankey; 
    ss_observed_loose_title << "Observed hits (MRD Track) chankey " << chankey; 
    ss_expected_loose_title << "Expected hits (MRD Track) chankey " << chankey; 

    TH1F *observed_strict = new TH1F(ss_observed_strict.str().c_str(),ss_observed_strict_title.str().c_str(),10,xmin,xmax);
    TH1F *expected_strict = new TH1F(ss_expected_strict.str().c_str(),ss_expected_strict_title.str().c_str(),10,xmin,xmax);
    TH1F *observed_loose = new TH1F(ss_observed_loose.str().c_str(),ss_observed_loose_title.str().c_str(),10,xmin,xmax);
    TH1F *expected_loose = new TH1F(ss_expected_loose.str().c_str(),ss_expected_loose_title.str().c_str(),10,xmin,xmax);

    if (int(chankey) < n_veto_pmts/2) { 

      //FMV Layer 1
      
      fmv_firstlayer.push_back(chankey); fmv_firstlayer_ymin.push_back(ymin); fmv_firstlayer_ymax.push_back(ymax); fmv_firstlayer_y.push_back(0.5*(ymin+ymax)); fmv_firstlayer_expected.push_back(0); fmv_firstlayer_observed.push_back(0); fmv_firstlayer_expected_track_strict.push_back(0); fmv_firstlayer_observed_track_strict.push_back(0); fmv_firstlayer_expected_track_loose.push_back(0); fmv_firstlayer_observed_track_loose.push_back(0); fmv_firstlayer_z = 0.5*(zmin+zmax); fmv_xmin = xmin; fmv_xmax = xmax; fmv_x = 0.5*(xmin+xmax);
      fmv_tank_firstlayer_expected.push_back(0); fmv_tank_firstlayer_observed.push_back(0);

      vector_observed_strict_layer1.push_back(observed_strict);
      vector_expected_strict_layer1.push_back(expected_strict);
      vector_observed_loose_layer1.push_back(observed_loose);
      vector_expected_loose_layer1.push_back(expected_loose);

    }
    else {
 
      //FMV Layer 2
      
      fmv_secondlayer.push_back(chankey); fmv_secondlayer_ymin.push_back(ymin); fmv_secondlayer_ymax.push_back(ymax); fmv_secondlayer_y.push_back(0.5*(ymin+ymax)); fmv_secondlayer_expected.push_back(0); fmv_secondlayer_observed.push_back(0); fmv_secondlayer_expected_track_strict.push_back(0); fmv_secondlayer_observed_track_strict.push_back(0); fmv_secondlayer_expected_track_loose.push_back(0); fmv_secondlayer_observed_track_loose.push_back(0); fmv_secondlayer_z = 0.5*(zmin+zmax);
      fmv_tank_secondlayer_expected.push_back(0); fmv_tank_secondlayer_observed.push_back(0);

      vector_observed_strict_layer2.push_back(observed_strict);
      vector_expected_strict_layer2.push_back(expected_strict);
      vector_observed_loose_layer2.push_back(observed_loose);
      vector_expected_loose_layer2.push_back(expected_loose);

    }
  }

  return true;
}


bool FMVEfficiency::Execute(){

  //The FMVEfficiency tool looks for coincidences between MRD/Tank clusters and FMV hits and estimates the efficiency of the different FMV paddles 
  //Make sure to execute ClusterFinder & TimeClustering beforehand in the ToolChain, so that m_all_clusters and MrdTimeClusters are available
  //There is a position dependent efficiency resolution for MRD coincidences since a track can be fitted for those
  //For the tank PMT, there is no possibility to fit a track until the LaserBall calibration, here no position resolution on the paddles can be used

  int get_ok;
  
  //  --------------------------------
  //--------Get tank clusters---------
  //  --------------------------------

  if (useTank){
    if (isData){
      get_ok = m_data->CStore.Get("ClusterMap",m_all_clusters);
      if (not get_ok) { Log("FMVEfficiency Tool: Error retrieving ClusterMap from CStore, did you run ClusterFinder beforehand?",v_error,verbosity); return false; }
    } else {
      get_ok = m_data->CStore.Get("ClusterMapMC",m_all_clusters_MC);
      if (not get_ok) { Log("FMVEfficiency Tool: Error retrieving ClusterMapMC from CStore, did you run ClusterFinder beforehand?",v_error,verbosity); return false; }
    } 
    get_ok = m_data->CStore.Get("ClusterMapDetkey",m_all_clusters_detkey);
    if (not get_ok) { Log("FMVEfficiency Tool: Error retrieving ClusterMapDetkey from CStore, did you run ClusterFinder beforehand?",v_error,verbosity); return false; }
  }

  //  --------------------------------
  //--------Get MRD clusters---------
  //  --------------------------------
  
  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not get_ok) { Log("RunValidation Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
  if (MrdTimeClusters.size()!=0){
    get_ok = m_data->CStore.Get("MrdDigitTimes",MrdDigitTimes);
    if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitTimes map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
    get_ok = m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
    if (not get_ok) { Log("EventDisplay Tool: Error retrieving MrdDigitChankeys, did you run TimeClustering beforehand",v_error,verbosity); return false;}
  }

  if (isData) m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);
  else m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData_MC);

  bool pmt_cluster=false;
  int n_pmt_hits=0;
  double max_charge=0.;
  double max_charge_pe=0.;
  double cluster_charge_pe=0;
  double cluster_time=0.; 
 
  //  --------------------------------
  //-----Loop through tank data-------
  //  --------------------------------
  

  if (useTank){
    if (isData){
      if (m_all_clusters){
        int clustersize = m_all_clusters->size();
        if (clustersize != 0){
          for(std::pair<double,std::vector<Hit>>&& apair : *m_all_clusters){
            double cluster_charge = apair.first;
            std::vector<Hit>&Hits = apair.second;
            std::vector<unsigned long> detkeys = m_all_clusters_detkey->at(cluster_charge);
            double global_time=0.;
            double global_charge=0.;
            int nhits=0;
            for (unsigned i_hit = 0; i_hit < Hits.size(); i_hit++){
              unsigned long detkey = detkeys.at(i_hit);
              double time = Hits.at(i_hit).GetTime();
              double charge = Hits.at(i_hit).GetCharge();
              if (pmt_gains[detkey] > 0) charge/=pmt_gains[detkey];
              global_time+=time;
              global_charge+=charge;
              nhits++;
            }
            if (nhits>0) global_time/=nhits;
            if (global_time < 2000. && global_charge > cluster_charge_pe) {
              pmt_cluster=true;
              cluster_charge_pe = global_charge;
              cluster_time = global_time;
            }
          }
        }
      }
    } else {
      if (m_all_clusters_MC){
        int clustersize = m_all_clusters_MC->size();
        if (clustersize != 0){
          for(std::pair<double,std::vector<MCHit>>&& apair : *m_all_clusters_MC){
            double cluster_charge = apair.first;
            std::vector<MCHit>&MCHits = apair.second;
            std::vector<unsigned long> detkeys = m_all_clusters_detkey->at(cluster_charge);
            double global_time=0.;
            double global_charge=0.;
            int nhits=0;
            for (unsigned i_hit = 0; i_hit < MCHits.size(); i_hit++){
              unsigned long detkey = detkeys.at(i_hit);
              double time = MCHits.at(i_hit).GetTime();
              double charge = MCHits.at(i_hit).GetCharge();
              global_time+=time;
              global_charge+=charge;
              nhits++;
            }
            if (nhits>0) global_time/=nhits;
            if (global_time < 2000. && global_charge > cluster_charge_pe) {
              pmt_cluster=true;
              cluster_charge_pe = global_charge;
              cluster_time = global_time;
            }
          }
        }
      }

    }
  }

  //  --------------------------------
  //---------Get FMV data-------------
  //  --------------------------------
  
  std::vector<unsigned long> hit_fmv_detkeys_first, hit_fmv_detkeys_second;
  std::vector<double> vector_fmv_times_first, vector_fmv_times_second;

  if (isData){
    if(TDCData){
      if (TDCData->size()==0){
        Log("FMVEfficiency tool: TDC data is empty in this event.",v_message,verbosity);
      } else {
        for (auto&& anmrdpmt : (*TDCData)){
          unsigned long chankey = anmrdpmt.first;
          Detector* thedetector = geom->ChannelToDetector(chankey);
          unsigned long detkey = thedetector->GetDetectorID();
          if (thedetector->GetDetectorElement()=="Veto") {
            if (int(chankey) < n_veto_pmts/2) {
              hit_fmv_detkeys_first.push_back(detkey);
              double fmv_time=0;
              int nhits_fmv=0;
              for(auto&& hitsonthismrdpmt : anmrdpmt.second){
                fmv_time+=hitsonthismrdpmt.GetTime();
                nhits_fmv++;
              }
              if (nhits_fmv!=0) fmv_time/=nhits_fmv;
              vector_fmv_times_first.push_back(fmv_time);
            }
            else { 
	      /*if (detkey == 14 || detkey == 18 || detkey == 19) continue;
              if (detkey == 23) detkey = 14;
	      else if (detkey == 24) detkey = 18;
	      else if (detkey == 25) detkey = 19;*/
              hit_fmv_detkeys_second.push_back(detkey);
              double fmv_time = 0;
              int nhits_fmv = 0;
              for(auto&&hitsonthismrdpmt : anmrdpmt.second){
                fmv_time+=hitsonthismrdpmt.GetTime();
                nhits_fmv++;
              }
              if (nhits_fmv!=0) fmv_time/=nhits_fmv;
              vector_fmv_times_second.push_back(fmv_time);
            }
          } 
        }
      }
    } else {
      Log("FMVEfficiency tool: No TDC data available in this event.",v_message,verbosity);
    }
  } else {
    if(TDCData_MC){
      if (TDCData_MC->size()==0){
        Log("FMVEfficiency tool: TDC data (MC) is empty in this event.",v_message,verbosity);
      } else {
        for (auto&& anmrdpmt : (*TDCData_MC)){
          unsigned long chankey = anmrdpmt.first;
          Detector* thedetector = geom->ChannelToDetector(chankey);
          unsigned long detkey = thedetector->GetDetectorID();
	  chankey-=first_fmv_chankey;
          detkey-=first_fmv_detkey;
          if (thedetector->GetDetectorElement()=="Veto") {
            if (int(chankey) < n_veto_pmts/2) {
              hit_fmv_detkeys_first.push_back(detkey);
              double fmv_time=0;
              int nhits_fmv=0;
              for(auto&& hitsonthismrdpmt : anmrdpmt.second){
                fmv_time+=hitsonthismrdpmt.GetTime();
                nhits_fmv++;
              }
              if (nhits_fmv!=0) fmv_time/=nhits_fmv;
              vector_fmv_times_first.push_back(fmv_time);
            }
            else { 
              hit_fmv_detkeys_second.push_back(detkey);
              double fmv_time = 0;
              int nhits_fmv = 0;
              for(auto&&hitsonthismrdpmt : anmrdpmt.second){
                fmv_time+=hitsonthismrdpmt.GetTime();
                nhits_fmv++;
              }
              if (nhits_fmv!=0) fmv_time/=nhits_fmv;
              vector_fmv_times_second.push_back(fmv_time);
            }
          } 
        }
      }
    } else {
      Log("FMVEfficiency tool: No TDC data available in this event.",v_message,verbosity);
    }


  }
 
  // How many FMV paddles were hit?
  
  unsigned int npaddles_Layer1 = hit_fmv_detkeys_first.size();
  unsigned int npaddles_Layer2 = hit_fmv_detkeys_second.size();
  num_paddles_Layer1->Fill(npaddles_Layer1);
  num_paddles_Layer2->Fill(npaddles_Layer2);

  //  --------------------------------
  //----Evaluate tank coincidences------
  //  --------------------------------

  if (useTank){
  if (pmt_cluster){
    if (npaddles_Layer1 == 1){
      for (unsigned int i_fmv = 0; i_fmv < hit_fmv_detkeys_first.size(); i_fmv++){
        time_diff_tank_Layer1->Fill(vector_fmv_times_first.at(i_fmv)-cluster_time);  
        if (isData){
          if ((vector_fmv_times_first.at(i_fmv)-cluster_time)<740 || (vector_fmv_times_first.at(i_fmv)-cluster_time)>840) continue;
        } else {
          if ((vector_fmv_times_first.at(i_fmv)-cluster_time)<-100 || (vector_fmv_times_first.at(i_fmv)-cluster_time)>100) continue; 
        }
        unsigned long detkey_first = hit_fmv_detkeys_first.at(i_fmv);
        fmv_tank_secondlayer_expected.at(detkey_first)++;
        if (std::find(hit_fmv_detkeys_second.begin(),hit_fmv_detkeys_second.end(),detkey_first+n_veto_pmts/2)!=hit_fmv_detkeys_second.end()){
          fmv_tank_secondlayer_observed.at(detkey_first)++;
        }
      }
    }

    if (npaddles_Layer2 == 1){
      for (unsigned int i_fmv = 0; i_fmv < hit_fmv_detkeys_second.size(); i_fmv++){
        time_diff_tank_Layer2->Fill(vector_fmv_times_second.at(i_fmv)-cluster_time); 
        if (isData){
          if ((vector_fmv_times_second.at(i_fmv)-cluster_time)<740 || (vector_fmv_times_second.at(i_fmv)-cluster_time)>840) continue;
        } else {
          if ((vector_fmv_times_second.at(i_fmv)-cluster_time)<-100 || (vector_fmv_times_second.at(i_fmv)-cluster_time)>100) continue;
        }
        unsigned long detkey_second = hit_fmv_detkeys_second.at(i_fmv);
        unsigned long detkey_first = detkey_second - n_veto_pmts/2;
        fmv_tank_firstlayer_expected.at(detkey_first)++;
        if (std::find(hit_fmv_detkeys_first.begin(),hit_fmv_detkeys_first.end(),detkey_first)!=hit_fmv_detkeys_first.end()){
          fmv_tank_firstlayer_observed.at(detkey_first)++;
        }
      }
    }
  }
  }

  //  --------------------------------
  //----Evaluate MRD coincidences-----
  //  --------------------------------
  
  //Only use events with one track in the MRD
  if (MrdTimeClusters.size()!=1) return true; 

  //Calculate mean time for MRD track
  double mrd_time = 0;
  int n_pmts = 0; 
  for(unsigned int thiscluster=0; thiscluster<MrdTimeClusters.size(); thiscluster++){
    std::vector<int> single_mrdcluster = MrdTimeClusters.at(thiscluster);
    int numdigits = single_mrdcluster.size();
    if (numdigits >=50) return true;	//reject noise events
    for(int thisdigit=0;thisdigit<numdigits;thisdigit++){
      int digit_value = single_mrdcluster.at(thisdigit);
      unsigned long chankey = mrddigitchankeysthisevent.at(digit_value);
      Detector *thedetector = geom->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement()=="MRD") {
        double mrdtimes=MrdDigitTimes.at(digit_value);
        mrd_time+=mrdtimes;
        n_pmts++;
      }
    }
  }
  if (n_pmts!=0) mrd_time/=n_pmts;

  //Get fitted MRD track
  int numsubevs, numtracksinev;
  m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents",numsubevs);
  m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);
  m_data->Stores["MRDTracks"]->Get("MRDTracks",theMrdTracks);
  BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(0));

  std::vector<int> PMTsHit;
  std::vector<int> LayersHit;
  Position StartVertex;
  Position StopVertex;
  thisTrackAsBoostStore->Get("StartVertex",StartVertex);
  thisTrackAsBoostStore->Get("StopVertex",StopVertex);
  thisTrackAsBoostStore->Get("PMTsHit",PMTsHit);
  thisTrackAsBoostStore->Get("LayersHit",LayersHit);

  //Find intersection of extrapolated MRD track with FMV
  double x_layer1, y_layer1, x_layer2, y_layer2;
  FindPaddleIntersection(StartVertex,StopVertex,x_layer1,y_layer1,fmv_firstlayer_z);
  FindPaddleIntersection(StartVertex,StopVertex,x_layer2,y_layer2,fmv_secondlayer_z);
  unsigned long hit_chankey_layer1 = 999999;
  unsigned long hit_chankey_layer2 = 999999;
  FindPaddleChankey(x_layer1, y_layer1, 1, hit_chankey_layer1);
  FindPaddleChankey(x_layer2, y_layer2, 2, hit_chankey_layer2);

  if (hit_chankey_layer1 == 999999){
    Log("FMVEfficiency tool: Did not find paddle with the desired intersection properties for FMV layer 1.",v_message,verbosity);
  }
  if (hit_chankey_layer2 == 999999){
    Log("FMVEfficiency tool: Did not find paddle with the desired intersection properties for FMV layer 2.",v_message,verbosity);
  }

  Log("FMVEfficiency tool: Intersection properties Layer 1: x_layer1 = "+std::to_string(x_layer1)+", y_layer1 = "+std::to_string(y_layer1),v_debug,verbosity);
  Log("FMVEfficiency tool: Intersection properties Layer 2: x_layer2 = "+std::to_string(x_layer2)+", y_layer2 = "+std::to_string(y_layer2),v_debug,verbosity);
  

  //if (npaddles_Layer1 == 1){
  if (npaddles_Layer1 >= 1 && npaddles_Layer1 <= padPerLayer){
    for (unsigned int i_fmv = 0; i_fmv < hit_fmv_detkeys_first.size(); i_fmv++){

      time_diff_paddle->Fill(hit_fmv_detkeys_first.at(i_fmv),vector_fmv_times_first.at(i_fmv)-mrd_time);
      //Check whether MRD & FMV Layer 1 fired in coincidence (time cut)
      fmv_layer1->Fill(hit_fmv_detkeys_first.at(i_fmv));
      time_diff_Layer1->Fill(vector_fmv_times_first.at(i_fmv)-mrd_time);
      if (fabs(vector_fmv_times_first.at(i_fmv)-mrd_time)>100.) continue;
      
      //Get properties of coincident MRD/FMV Layer 1 hit
      unsigned long detkey_first = hit_fmv_detkeys_first.at(i_fmv);
      track_diff_x_Layer1->Fill(x_layer1-fmv_x);
      track_diff_y_Layer1->Fill(y_layer1-fmv_firstlayer_y.at(detkey_first));
      track_diff_xy_Layer1->Fill(x_layer1-fmv_x,y_layer1-fmv_firstlayer_y.at(detkey_first));
      if (y_layer1-fmv_firstlayer_y.at(detkey_first) > -0.4 && y_layer1-fmv_firstlayer_y.at(detkey_first)< 0.6 && fabs(x_layer1-fmv_x)<1.6){
        vector_expected_loose_layer2.at(detkey_first)->Fill(x_layer1);
        fmv_secondlayer_expected_track_loose.at(detkey_first)++;
        track_diff_x_loose_Layer1->Fill(x_layer1-fmv_x);
        track_diff_y_loose_Layer1->Fill(y_layer1-fmv_firstlayer_y.at(detkey_first));
        track_diff_xy_loose_Layer1->Fill(x_layer1-fmv_x,y_layer1-fmv_firstlayer_y.at(detkey_first));
      }

      //Check second layer whether it saw a coincident hit
      fmv_secondlayer_expected.at(detkey_first)++;
      if (hit_chankey_layer1 == detkey_first) {
        vector_expected_strict_layer2.at(detkey_first)->Fill(x_layer1);
        fmv_secondlayer_expected_track_strict.at(detkey_first)++;
        track_diff_x_strict_Layer1->Fill(x_layer1-fmv_x);
        track_diff_y_strict_Layer1->Fill(y_layer1-fmv_firstlayer_y.at(detkey_first));
        track_diff_xy_strict_Layer1->Fill(x_layer1-fmv_x,y_layer1-fmv_firstlayer_y.at(detkey_first));
      }
      for (int i_second=0; i_second < (int) hit_fmv_detkeys_second.size(); i_second++){
        fmv_layer1_layer2->Fill(detkey_first,hit_fmv_detkeys_second.at(i_second));
      }
      if (std::find(hit_fmv_detkeys_second.begin(),hit_fmv_detkeys_second.end(),detkey_first+n_veto_pmts/2)!=hit_fmv_detkeys_second.end()) {
        if (y_layer1-fmv_firstlayer_y.at(detkey_first) > -0.4 && y_layer1-fmv_firstlayer_y.at(detkey_first)<0.6 && fabs(x_layer1-fmv_x)<1.6){
          fmv_secondlayer_observed_track_loose.at(detkey_first)++;
          vector_observed_loose_layer2.at(detkey_first)->Fill(x_layer1);
	  std::vector<unsigned long>::iterator it = std::find(hit_fmv_detkeys_second.begin(),hit_fmv_detkeys_second.end(),detkey_first+n_veto_pmts/2);
	  int index = std::distance(hit_fmv_detkeys_second.begin(),it);
          time_diff_paddle_coinc->Fill(detkey_first+n_veto_pmts/2,vector_fmv_times_second.at(index)-vector_fmv_times_first.at(i_fmv));
        }
        fmv_secondlayer_observed.at(detkey_first)++;
	std::vector<unsigned long>::iterator it = std::find(hit_fmv_detkeys_second.begin(),hit_fmv_detkeys_second.end(),detkey_first+n_veto_pmts/2);
        int index = std::distance(hit_fmv_detkeys_second.begin(),it);
        time_diff_paddle_simplecoinc->Fill(detkey_first+n_veto_pmts/2,vector_fmv_times_second.at(index)-vector_fmv_times_first.at(i_fmv));
        if (hit_chankey_layer1 == detkey_first){
          vector_observed_strict_layer2.at(detkey_first)->Fill(x_layer1);
          fmv_secondlayer_observed_track_strict.at(detkey_first)++;
        }
      }
    }
  }

  //if (npaddles_Layer2 == 1){
  if (npaddles_Layer2 >= 1 && npaddles_Layer2 <=padPerLayer){
    for (unsigned int i_fmv = 0; i_fmv < hit_fmv_detkeys_second.size(); i_fmv++){

      time_diff_paddle->Fill(hit_fmv_detkeys_second.at(i_fmv),vector_fmv_times_second.at(i_fmv)-mrd_time-20);

      //Check whether MRD & FMV Layer 2 fired in coincidence (time cut)
      time_diff_Layer2->Fill(vector_fmv_times_second.at(i_fmv)-mrd_time);
      fmv_layer2->Fill(hit_fmv_detkeys_second.at(i_fmv));
      if (fabs(vector_fmv_times_second.at(i_fmv)-mrd_time-20.)>100.) continue;
      
      //Get properties of coincident MRD/FMV Layer 2 hit
      unsigned long detkey_second = hit_fmv_detkeys_second.at(i_fmv);
      unsigned long detkey_first = detkey_second - n_veto_pmts/2;
      track_diff_x_Layer2->Fill(x_layer2-fmv_x);
      track_diff_y_Layer2->Fill(y_layer2-fmv_secondlayer_y.at(detkey_first));
      track_diff_xy_Layer2->Fill(x_layer2-fmv_x,y_layer2-fmv_secondlayer_y.at(detkey_first));
      if (y_layer2-fmv_secondlayer_y.at(detkey_first) > -0.4 && y_layer2-fmv_secondlayer_y.at(detkey_first)< 0.6 && fabs(x_layer2-fmv_x)<1.6){
        vector_expected_loose_layer1.at(detkey_first)->Fill(x_layer2);
        fmv_firstlayer_expected_track_loose.at(detkey_first)++;
        track_diff_x_loose_Layer2->Fill(x_layer2-fmv_x);
        track_diff_y_loose_Layer2->Fill(y_layer2-fmv_secondlayer_y.at(detkey_first));
        track_diff_xy_loose_Layer2->Fill(x_layer2-fmv_x,y_layer2-fmv_secondlayer_y.at(detkey_first));
      }

      //Check first layer whether it saw a coincident hit
      fmv_firstlayer_expected.at(detkey_first)++;
      if (hit_chankey_layer2 == detkey_second) {
        vector_expected_strict_layer1.at(detkey_first)->Fill(x_layer2);
        fmv_firstlayer_expected_track_strict.at(detkey_first)++;
        track_diff_x_strict_Layer2->Fill(x_layer2-fmv_x);
        track_diff_y_strict_Layer2->Fill(y_layer2-fmv_secondlayer_y.at(detkey_first));
        track_diff_xy_strict_Layer2->Fill(x_layer2-fmv_x,y_layer2-fmv_secondlayer_y.at(detkey_first));
      }
      if (std::find(hit_fmv_detkeys_first.begin(),hit_fmv_detkeys_first.end(),detkey_first) != hit_fmv_detkeys_first.end()) {
        if (y_layer2-fmv_secondlayer_y.at(detkey_first) > -0.4 && y_layer2-fmv_secondlayer_y.at(detkey_first)<0.6 && fabs(x_layer2-fmv_x)<1.6){
          fmv_firstlayer_observed_track_loose.at(detkey_first)++;
          vector_observed_loose_layer1.at(detkey_first)->Fill(x_layer2);
	  std::vector<unsigned long>::iterator it = std::find(hit_fmv_detkeys_first.begin(),hit_fmv_detkeys_first.end(),detkey_first);
          int index = std::distance(hit_fmv_detkeys_first.begin(),it);
          time_diff_paddle_coinc->Fill(detkey_first,vector_fmv_times_first.at(index)-vector_fmv_times_second.at(i_fmv));
        }
	std::vector<unsigned long>::iterator it = std::find(hit_fmv_detkeys_first.begin(),hit_fmv_detkeys_first.end(),detkey_first);
        int index = std::distance(hit_fmv_detkeys_first.begin(),it);
        time_diff_paddle_simplecoinc->Fill(detkey_first,vector_fmv_times_first.at(index)-vector_fmv_times_second.at(i_fmv));
       for (int i_first=0; i_first < (int) hit_fmv_detkeys_first.size(); i_first++){
        fmv_layer1_layer2->Fill(hit_fmv_detkeys_first.at(i_first),detkey_second);
       }
        fmv_firstlayer_observed.at(detkey_first)++;
        if (hit_chankey_layer2 == detkey_second){
          vector_observed_strict_layer1.at(detkey_first)->Fill(x_layer2);
          fmv_firstlayer_observed_track_strict.at(detkey_first)++;
        }
      }
    }
  }


  return true;
}


bool FMVEfficiency::Finalise(){

  //  --------------------------------
  //---Fill chankey-wise histograms---
  //  --------------------------------

  file->cd();
  fmv_layer1_layer2->Write();
  for (unsigned int i_fmv = 0; i_fmv < fmv_secondlayer.size(); i_fmv++){
    fmv_observed_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_observed.at(i_fmv));
    fmv_expected_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_expected.at(i_fmv));
    fmv_observed_track_strict_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_observed_track_strict.at(i_fmv));
    fmv_expected_track_strict_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_expected_track_strict.at(i_fmv));
    fmv_observed_track_loose_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_observed_track_loose.at(i_fmv));
    fmv_expected_track_loose_layer1->SetBinContent(i_fmv+1,fmv_firstlayer_expected_track_loose.at(i_fmv));
    fmv_observed_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_observed.at(i_fmv));
    fmv_expected_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_expected.at(i_fmv));
    fmv_observed_track_strict_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_observed_track_strict.at(i_fmv));
    fmv_expected_track_strict_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_expected_track_strict.at(i_fmv));
    fmv_observed_track_loose_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_observed_track_loose.at(i_fmv));
    fmv_expected_track_loose_layer2->SetBinContent(i_fmv+1,fmv_secondlayer_expected_track_loose.at(i_fmv));
    if (useTank){
      fmv_tank_observed_layer1->SetBinContent(i_fmv+1,fmv_tank_firstlayer_observed.at(i_fmv));
      fmv_tank_expected_layer1->SetBinContent(i_fmv+1,fmv_tank_firstlayer_expected.at(i_fmv));
      fmv_tank_observed_layer2->SetBinContent(i_fmv+1,fmv_tank_secondlayer_observed.at(i_fmv));
      fmv_tank_expected_layer2->SetBinContent(i_fmv+1,fmv_tank_secondlayer_expected.at(i_fmv));
    }
  }
  file->Write();
  file->Close();
  delete file;


  return true;
}


bool FMVEfficiency::FindPaddleIntersection(Position startpos, Position endpos, double &x, double &y, double z){

        double DirX = endpos.X()-startpos.X();
        double DirY = endpos.Y()-startpos.Y();
        double DirZ = endpos.Z()-startpos.Z();

        if (fabs(DirZ) < 0.001) Log("MrdPaddleEfficiencyPreparer tool: StartVertex = EndVertex! Track was not fitted well",v_error,verbosity);

        double frac = (z - startpos.Z())/DirZ;

        x = startpos.X()+frac*DirX;
        y = startpos.Y()+frac*DirY;

        return true;

}



bool FMVEfficiency::FindPaddleChankey(double x, double y, int layer, unsigned long &chankey){

        bool found_chankey = false;
        for (unsigned int i_channel = 0; i_channel < 13; i_channel++){

                if (found_chankey) break;

                unsigned long chankey_tmp = (unsigned long) i_channel+first_fmv_chankey;
                if (layer == 2) chankey_tmp += 13;
                Detector *mrdpmt = geom->ChannelToDetector(chankey_tmp);
                unsigned long detkey = mrdpmt->GetDetectorID();
                Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
                double xmin = mrdpaddle->GetXmin();
                double xmax = mrdpaddle->GetXmax();
                double ymin = mrdpaddle->GetYmin();
                double ymax = mrdpaddle->GetYmax();

                //Check if expected hit was within the channel or not

                if (xmin <= x && xmax >= x && ymin <= y && ymax >= y){
                        chankey = chankey_tmp;
                        found_chankey = true;
                }
        }

        return true;

}

