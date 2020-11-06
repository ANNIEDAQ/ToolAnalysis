#include "BeamClusterPlots.h"

BeamClusterPlots::BeamClusterPlots():Tool(){}


bool BeamClusterPlots::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load the default threshold settings for finding pulses
  verbosity = 3;
  outputfile = "BeamClusterPlotsDefault_";
  PromptPEMin = 500;
  PromptWindowMin = 0;  //ns
  PromptWindowMax = 2000; //ns
  DelayedWindowMin = 11000;  //ns
  DelayedWindowMax = 70000; //ns
  NeutronPEMin = 5;
  NeutronPEMax = 100;
  PEPerMeV = 12.;

  //Load any configurables set in the config file
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("OutputFile",outputfile);
  m_variables.Get("PromptPEMin",PromptPEMin);
  m_variables.Get("PromptWindowMin",PromptWindowMin);
  m_variables.Get("PromptWindowMax",PromptWindowMax);
  m_variables.Get("DelayedWindowMin",DelayedWindowMin);
  m_variables.Get("DelayedWindowMax",DelayedWindowMax);
  m_variables.Get("NeutronPEMin",NeutronPEMin);
  m_variables.Get("NeutronPEMax",NeutronPEMax);
  m_variables.Get("PEPerMeVInTank",PEPerMeV);

  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap",ChannelKeyToSPEMap);

  //Initialize needed information for occupancy plots; taken from Michael's code 
  auto get_geometry= m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",geom);
  if(!get_geometry){
  	Log("DigitBuilder Tool: Error retrieving Geometry from ANNIEEvent!",v_error,verbosity); 
  	return false; 
  }

  std::string rootfile_out_prefix="_BCA";
  std::string rootfile_out_root=".root";
  std::string rootfile_out_name=outputfile+rootfile_out_prefix+rootfile_out_root;

  bca_file_out=new TFile(rootfile_out_name.c_str(),"RECREATE"); //create one root file for each run to save the detailed plots and fits for all PMTs 

  bca_file_out->cd();

  //You have to do this second, or ROOT breaks
  this->InitializeHistograms();
  this->SetHistogramLabels();

  return true;
}


bool BeamClusterPlots::Execute(){

  //First, get clusters from the BoostStore
  //Clean AmBe triggers with all cluster info first. 
  bool get_clusters = m_data->CStore.Get("ClusterMap",m_all_clusters);
  if(!get_clusters){
    std::cout << "BeamClusterPlots tool: No clusters found!" << std::endl;
    return true;
  }
  bool got_ccp = m_data->Stores.at("ANNIEEvent")->Get("ClusterChargePoints", ClusterChargePoints);
  bool got_ccb = m_data->Stores.at("ANNIEEvent")->Get("ClusterChargeBalanaces", ClusterChargeBalances);
  bool got_cmpe = m_data->Stores.at("ANNIEEvent")->Get("ClusterMaxPEs", ClusterMaxPEs);

  if(verbosity>3) std::cout << "BeamClusterPlots Tool: looping through clusters to get cluster info now" << std::endl;
  if(verbosity>3) std::cout << "BeamClusterPlots Tool: number of clusters: " << m_all_clusters->size() << std::endl;
  
  double max_prompt_clustertime = -1;
  double max_prompt_clusterPE = -1;
  std::vector<double> delayed_cluster_times;
  std::vector<double> delayed_ncandidate_times;

  for (std::pair<double,std::vector<Hit>>&& cluster_pair : *m_all_clusters) {
    double cluster_charge = 0;
    double cluster_time = cluster_pair.first;
    double cluster_PE = 0;
    double num_hits = 0;
    std::vector<Hit> cluster_hits = cluster_pair.second;
    for (int i = 0; i<(int)cluster_hits.size(); i++){
      int hit_ID = cluster_hits.at(i).GetTubeId();
      std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(hit_ID);
      if(it != ChannelKeyToSPEMap.end()){ //Charge to SPE conversion is available
        cluster_charge+=cluster_hits.at(i).GetCharge();
        cluster_PE+=(cluster_hits.at(i).GetCharge() / ChannelKeyToSPEMap.at(hit_ID));
        num_hits += 1;
      } else {
        if(verbosity>2){
          std::cout << "FOUND A HIT FOR CHANNELKEY " << hit_ID << "BUT NO CONVERSION " <<
              "TO PE AVAILABLE.  SKIPPING PE." << std::endl;
        }
      }
    }
    if(verbosity>3) std::cout << "BeamClusterPlots Tool: cluster_time is : " << cluster_time << std::endl;
    if((cluster_time > PromptWindowMin) && (cluster_time < PromptWindowMax)){
      if(cluster_PE > max_prompt_clusterPE){
        max_prompt_clusterPE = cluster_PE;
        max_prompt_clustertime = cluster_time;
      }
      if(verbosity>4) std::cout << "Found a prompt cluster candidate.  Filling prompt info" << std::endl;
      hist_prompt_Time->Fill(cluster_time);
      hist_prompt_PE->Fill(cluster_PE);
      hist_prompt_TimeVPE->Fill(cluster_time,cluster_PE);
      hist_prompt_PEVNHit->Fill(cluster_PE,num_hits);
      double ChargePointZ = ClusterChargePoints.at(cluster_time).Z();
      hist_prompt_ChargePoint->Fill(cluster_PE,ChargePointZ);
      double ChargeBalance = ClusterChargeBalances.at(cluster_time);
      hist_prompt_ChargeBalance->Fill(cluster_PE,ChargeBalance);
      double max_PE = ClusterMaxPEs.at(cluster_time);
      hist_prompt_TotPEVsMaxPE->Fill(cluster_PE,max_PE);
    }
    if((cluster_time > DelayedWindowMin) && (cluster_time < DelayedWindowMax)){
      delayed_cluster_times.push_back(cluster_time);
      if(cluster_PE > NeutronPEMin && cluster_PE < NeutronPEMax) delayed_ncandidate_times.push_back(cluster_time);
      if(verbosity>4) std::cout << "Found a delayed cluster candidate.  Filling delayed info" << std::endl;
      hist_delayed_Time->Fill(cluster_time);
      hist_delayed_PE->Fill(cluster_PE);
      hist_delayed_TimeVPE->Fill(cluster_time,cluster_PE);
      hist_delayed_PEVNHit->Fill(cluster_PE,num_hits);
      double ChargePointZ = ClusterChargePoints.at(cluster_time).Z();
      hist_delayed_ChargePoint->Fill(cluster_PE,ChargePointZ);
      double ChargeBalance = ClusterChargeBalances.at(cluster_time);
      hist_delayed_ChargeBalance->Fill(cluster_PE,ChargeBalance);
      double max_PE = ClusterMaxPEs.at(cluster_time);
      hist_delayed_TotPEVsMaxPE->Fill(cluster_PE,max_PE);
    }
  }

  //Then, make multiplicity histograms
  if(max_prompt_clusterPE>PromptPEMin){
    hist_prompt_delayed_multiplicity->Fill(delayed_cluster_times.size()); 
    for(int j = 0; j< (int) delayed_cluster_times.size(); j++){
      hist_prompt_delayed_deltat->Fill(delayed_cluster_times.at(j) - max_prompt_clustertime);
    }

    hist_prompt_neutron_multiplicity->Fill(delayed_ncandidate_times.size()); 
    hist_prompt_neutron_multiplicityvstankE->Fill(delayed_ncandidate_times.size(),max_prompt_clusterPE/PEPerMeV); 
    for(int j = 0; j< (int) delayed_ncandidate_times.size(); j++){
      hist_prompt_neutron_deltat->Fill(delayed_ncandidate_times.at(j) - max_prompt_clustertime);
    }
  }

  return true;
}


bool BeamClusterPlots::Finalise(){
  this->WriteHistograms();
  bca_file_out->Close();
  std::cout << "BeamClusterPlots tool exitting" << std::endl;
  return true;
}


void BeamClusterPlots::InitializeHistograms(){
  //All Event Histograms
  hist_prompt_Time = new TH1F("hist_prompt_Time","Prompt cluster times",(PromptWindowMax-PromptWindowMin),PromptWindowMin,PromptWindowMax);
  hist_prompt_PE = new TH1F("hist_prompt_PE","Prompt cluster total PEs",6000,0,6000);
  hist_prompt_TimeVPE = new TH2F("hist_prompt_TimeVPE","Hit time as a function of prompt cluster PE",(PromptWindowMax-PromptWindowMin),PromptWindowMin,PromptWindowMax,6000,0,6000);
  hist_prompt_PEVNHit = new TH2F("hist_prompt_PEVNHit","Number of photoelectrons as a function of number of pulses for prompt clusters",6000,0,6000,500,0,500);
  hist_prompt_ChargePoint = new TH2F("hist_prompt_ChargePoint","Number of photoelectrons as a function of charge point z-component for prompt clusters",6000, 0, 6000,100,0,1);
  hist_prompt_ChargeBalance = new TH2F("hist_prompt_ChargeBalance","Number of photoelectrons as a function of charge balance cut for prompt clusters",6000,0,6000,100,0,1);
  hist_prompt_TotPEVsMaxPE = new TH2F("hist_prompt_TotPEVsMaxPE","Number of photoelectrons as a function of the highest photoelectron count for prompt clusters",6000,0,6000,2000,0,2000);
  hist_delayed_Time = new TH1F("hist_delayed_Time","Delayed cluster times",(DelayedWindowMax-DelayedWindowMin),DelayedWindowMin,DelayedWindowMax);
  hist_delayed_PE = new TH1F("hist_delayed_PE","Delayed cluster total PEs",6000,0,6000);
  hist_delayed_TimeVPE = new TH2F("hist_delayed_TimeVPE","Hit time as a function of delayed cluster PE",(DelayedWindowMax-DelayedWindowMin),DelayedWindowMin,DelayedWindowMax,6000,0,6000);
  hist_delayed_PEVNHit = new TH2F("hist_delayed_PEVNHit","Number of photoelectrons as a function of number of pulses for delayed clusters",6000,0,6000,500,0,500);
  hist_delayed_ChargePoint = new TH2F("hist_delayed_ChargePoint","Number of photoelectrons as a function of charge point z-component for delayed clusters",6000, 0, 6000,100,0,1);
  hist_delayed_ChargeBalance = new TH2F("hist_delayed_ChargeBalance","Number of photoelectrons as a function of charge balance cut for delayed clusters",6000,0,6000,100,0,1);
  hist_delayed_TotPEVsMaxPE = new TH2F("hist_delayed_TotPEVsMaxPE","Number of photoelectrons as a function of the highest photoelectron count for delayed clusters",6000,0,6000,2000,0,2000);
  
  //Then, make two bonus histograms; delta t distribution between largest prompt cluster
  hist_prompt_delayed_multiplicity = new TH1F("hist_prompt_delayed_multiplicity","Multiplicity of delayed clusters",20,0,20); 
  hist_prompt_neutron_multiplicity = new TH1F("hist_prompt_neutron_multiplicity","Multiplicity of delayed clusters that are neutron candidates",20,0,20); 
  hist_prompt_delayed_deltat = new TH1F("hist_prompt_delayed_deltat","Time difference between highest PE prompt cluster and delayed cluster",70000,0,70000);
  hist_prompt_neutron_deltat = new TH1F("hist_prompt_neutron_deltat","Time difference between highest PE prompt cluster and neutron cluster",70000,0,70000);
  hist_prompt_neutron_multiplicityvstankE = new TH2F("hist_prompt_neutron_multiplicityvstankE","Delayed cluster multiplicity vs. visible tank energy estimate",20,0,20,500,0,1000);
}

void BeamClusterPlots::WriteHistograms(){
  bca_file_out->cd();
  //All Event Histograms
  hist_prompt_Time->Write();
  hist_prompt_PE->Write();
  hist_prompt_TimeVPE->Write();
  hist_prompt_PEVNHit->Write();
  hist_prompt_ChargePoint->Write();
  hist_prompt_ChargeBalance->Write();
  hist_prompt_TotPEVsMaxPE->Write();
  hist_delayed_Time->Write();
  hist_delayed_PE->Write();
  hist_delayed_TimeVPE->Write();
  hist_delayed_PEVNHit->Write();
  hist_delayed_ChargePoint->Write();
  hist_delayed_ChargeBalance->Write();
  hist_delayed_TotPEVsMaxPE->Write();
  
  //Then, make two bonus histograms->Write(); delta t distribution between largest prompt cluster
  hist_prompt_delayed_multiplicity->Write();
  hist_prompt_neutron_multiplicity->Write();
  hist_prompt_delayed_deltat->Write();
  hist_prompt_neutron_deltat->Write();
  hist_prompt_neutron_multiplicityvstankE->Write();
}

void BeamClusterPlots::SetHistogramLabels(){
  //All Event Histograms
  hist_prompt_Time->GetXaxis()->SetTitle("Prompt cluster time [ns]");
  hist_prompt_PE->GetXaxis()->SetTitle("Prompt PE");
  hist_prompt_TimeVPE->GetXaxis()->SetTitle("Prompt cluster time [ns]");
  hist_prompt_PEVNHit->GetXaxis()->SetTitle("Prompt cluster PE");
  hist_prompt_ChargePoint->GetXaxis()->SetTitle("Prompt cluster PE");
  hist_prompt_ChargeBalance->GetXaxis()->SetTitle("Prompt cluster PE");
  hist_prompt_TotPEVsMaxPE->GetXaxis()->SetTitle("Prompt cluster PE");
  hist_delayed_Time->GetXaxis()->SetTitle("Delayed cluster time [ns]");
  hist_delayed_PE->GetXaxis()->SetTitle("Delayed PE");
  hist_delayed_TimeVPE->GetXaxis()->SetTitle("Delayed cluster time [ns]");
  hist_delayed_PEVNHit->GetXaxis()->SetTitle("Delayed cluster PE");
  hist_delayed_ChargePoint->GetXaxis()->SetTitle("Delayed cluster PE");
  hist_delayed_ChargeBalance->GetXaxis()->SetTitle("Delayed cluster PE");
  hist_delayed_TotPEVsMaxPE->GetXaxis()->SetTitle("Delayed cluster PE");
  
  //Then, make two bonus histograms->GetXaxis()->SetTitle(); delta t distribution between largest prompt cluster
  hist_prompt_delayed_multiplicity->GetXaxis()->SetTitle("Delayed cluster multiplicity");
  hist_prompt_neutron_multiplicity->GetXaxis()->SetTitle("Delayed cluster neutron candidate multiplicity"); 
  hist_prompt_delayed_deltat->GetXaxis()->SetTitle("Delta cluster time [ns]");
  hist_prompt_neutron_deltat->GetXaxis()->SetTitle("Delta cluster time [ns]");
  hist_prompt_neutron_multiplicityvstankE->GetXaxis()->SetTitle("Delayed cluster neutron candidate multiplicity");


  hist_prompt_TimeVPE->GetYaxis()->SetTitle("Prompt cluster PE");
  hist_prompt_PEVNHit->GetYaxis()->SetTitle("Prompt cluster number of pulses");
  hist_prompt_ChargePoint->GetYaxis()->SetTitle("Charge point z-component");
  hist_prompt_ChargeBalance->GetYaxis()->SetTitle("Charge balance parameter");
  hist_prompt_TotPEVsMaxPE->GetYaxis()->SetTitle("Maximum single hit PE");
  hist_delayed_TimeVPE->GetYaxis()->SetTitle("Delayed cluster PE");
  hist_delayed_PEVNHit->GetYaxis()->SetTitle("Delayed cluster number of pulses");
  hist_delayed_ChargePoint->GetYaxis()->SetTitle("Charge point z-component");
  hist_delayed_ChargeBalance->GetYaxis()->SetTitle("Charge balance parameter");
  hist_delayed_TotPEVsMaxPE->GetYaxis()->SetTitle("Maximum single hit PE");
  
  //Then, make two bonus histograms->GetYaxis()->SetTitle(); delta t distribution between largest prompt cluster
  hist_prompt_neutron_multiplicityvstankE->GetYaxis()->SetTitle("Visible tank energy estimate [MeV]");
}
