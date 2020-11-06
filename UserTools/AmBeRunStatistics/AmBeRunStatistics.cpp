#include "AmBeRunStatistics.h"

AmBeRunStatistics::AmBeRunStatistics():Tool(){}


bool AmBeRunStatistics::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load the default threshold settings for finding pulses
  verbosity = 3;
  outputfile = "AmBeStatisticsDefault_";
  SWindowMin = 0;  //ns
  SWindowMax = 4000; //ns
  S1Threshold = 10; //V
  S2Threshold = 10; //V
  ClusterPEMin = 5; //Neighborhood of 5 PE when gain-matched at 7E6
  ClusterPEMax = 150; //Neighborhood of 100 PE when gain-matched at 7E6
  DeltaTimeThreshold = 50; //ns

  //Load any configurables set in the config file
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("OutputFile",outputfile);
  m_variables.Get("SiPMWindowMin",SWindowMin);
  m_variables.Get("SiPMWindowMax",SWindowMax);
  m_variables.Get("SiPM1TriggerThreshold", S1Threshold); 
  m_variables.Get("SiPM2TriggerThreshold", S2Threshold); 
  m_variables.Get("ClusterPEMin", ClusterPEMin); 
  m_variables.Get("ClusterPEMax", ClusterPEMax); 
  m_variables.Get("DeltaTimeThreshold", DeltaTimeThreshold); 


  m_data->CStore.Get("AuxChannelNumToTypeMap",AuxChannelNumToTypeMap);

  std::string rootfile_out_prefix="_AmBeStatistics";
  std::string rootfile_out_root=".root";
  std::string rootfile_out_name=outputfile+rootfile_out_prefix+rootfile_out_root;

  ambe_file_out=new TFile(rootfile_out_name.c_str(),"RECREATE"); //create one root file for each run to save the detailed plots and fits for all PMTs 

  ambe_file_out->cd();

  this->InitializeHistograms();
  this->SetHistogramLabels();

  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap",ChannelKeyToSPEMap);

  std::cout << "AmBeRunStatistics Tool Initialized" << std::endl;
  return true;
}


bool AmBeRunStatistics::Execute(){
  //ambe_file_out->cd();
  NumberOfEvents+=1;

  // Get a pointer to the ANNIEEvent Store
  auto* annie_event = m_data->Stores["ANNIEEvent"];

  annie_event->Get("CalibratedADCAuxData", calibrated_auxwaveform_map);
  annie_event->Get("RecoADCAuxHits", aux_pulse_map);

  ADCPulse SiPM1_MaxPulse;
  ADCPulse SiPM2_MaxPulse;
  int SiPM1_NumPulses;  //Number of pulses crossing the S1Threshold
  int SiPM2_NumPulses;  //Number of pulses crossing the S2Threshold
  this->GetSiPMPulseInfo(SiPM1_MaxPulse,SiPM2_MaxPulse,SiPM1_NumPulses,SiPM2_NumPulses);
  if(verbosity>3){
    std::cout << "AmBeRunStatistics Tool: Number of S1,S2 pulses above threshold: " << 
        SiPM1_NumPulses << "," << SiPM2_NumPulses << std::endl;
  }

  //If there's at least one pulse, fill the non-triggered maxima histograms 
  if(SiPM1_NumPulses>0){ 
    h_SiPM1_Amplitude->Fill(SiPM1_MaxPulse.amplitude());
    h_SiPM1_Charge->Fill(SiPM1_MaxPulse.charge());
    h_SiPM1_PeakTime->Fill(SiPM1_MaxPulse.peak_time());
  }
  if(SiPM2_NumPulses>0){ 
    h_SiPM2_Amplitude->Fill(SiPM2_MaxPulse.amplitude());
    h_SiPM2_Charge->Fill(SiPM2_MaxPulse.charge());
    h_SiPM2_PeakTime->Fill(SiPM2_MaxPulse.peak_time());
  }

  double deltat=0., charge_ratio=0.;
  if(SiPM1_NumPulses>0  && SiPM1_NumPulses>0){ 
    h_S1S2_Amplitudes->Fill(SiPM1_MaxPulse.amplitude(),SiPM2_MaxPulse.amplitude());
    deltat = SiPM1_MaxPulse.peak_time() - SiPM2_MaxPulse.peak_time();
    h_S1S2_Deltat->Fill(deltat);
    charge_ratio = SiPM1_MaxPulse.charge() / SiPM2_MaxPulse.charge();
    h_S1S2_Chargeratio->Fill(charge_ratio);
  }

  if(verbosity>3){
    std::cout << "AmBeRunStatistics Tool: S1-S2 pulse height times (ns): " << 
         deltat << std::endl;
    std::cout << "AmBeRunStatistics Tool: S1/S2 charge ratio: " << 
         charge_ratio << std::endl;
  }

  bool CleanAmBeTrigger = true;
  double MeanSiPMPulseTime = (SiPM1_MaxPulse.peak_time() + SiPM2_MaxPulse.peak_time())/2.;

  //See if any cluster is within DeltaTimeThreshold ns of the mean SiPM pulse time.
  bool get_clusters = m_data->CStore.Get("ClusterMap",m_all_clusters);
  if(!get_clusters){
    std::cout << "AmBeRunStatistics tool: No clusters found!" << std::endl;
    return true;
  }

  bool ClusterNearTrigger = false;
  for (std::pair<double,std::vector<Hit>>&& cluster_pair : *m_all_clusters) {
    double cluster_time = cluster_pair.first;
    if(abs(cluster_time-MeanSiPMPulseTime) < DeltaTimeThreshold) ClusterNearTrigger = true;
  }

  // Check the criteria for a valid clean AmBe trigger event
  if ((SiPM1_NumPulses != 1 || SiPM2_NumPulses != 1) || ClusterNearTrigger || 
      (abs(deltat) > DeltaTimeThreshold) || 
      (SiPM1_MaxPulse.peak_time() < SWindowMin || SiPM1_MaxPulse.peak_time() > SWindowMax) ||
      (SiPM2_MaxPulse.peak_time() < SWindowMin || SiPM2_MaxPulse.peak_time() > SWindowMax)) CleanAmBeTrigger = false;

  if(!CleanAmBeTrigger) return true;
  NumberOfCleanTriggers+=1;

  if(verbosity>3) std::cout << "AmBeRunStatistics Tool: Found a clean AmBe SiPM trigger" << std::endl;
  h_SiPM1_AmplitudeCleanPromptTrig->Fill(SiPM1_MaxPulse.amplitude());
  h_SiPM2_AmplitudeCleanPromptTrig->Fill(SiPM2_MaxPulse.amplitude());
  h_SiPM1_ChargeCleanPromptTrig->Fill(SiPM1_MaxPulse.charge());
  h_SiPM2_ChargeCleanPromptTrig->Fill(SiPM2_MaxPulse.charge());
  h_SiPM1_PeakTimeCleanPromptTrig->Fill(SiPM1_MaxPulse.peak_time());
  h_SiPM2_PeakTimeCleanPromptTrig->Fill(SiPM2_MaxPulse.peak_time());
  h_S1S2_AmplitudesCleanPromptTrig->Fill(SiPM1_MaxPulse.amplitude(),SiPM2_MaxPulse.amplitude());
  h_S1S2_ChargeratioCleanPromptTrig->Fill(charge_ratio);
  h_S1S2_DeltatCleanPromptTrig->Fill(deltat);

  double cluster_charge=0;
  double cluster_time=0;
  double cluster_PE=0;
  if(verbosity>3) std::cout << "AmBeRunStatistics Tool: looping through clusters to get cluster info now" << std::endl;
  if(verbosity>3) std::cout << "AmBeRunStatistics Tool: number of clusters: " << m_all_clusters->size() << std::endl;
  for (std::pair<double,std::vector<Hit>>&& cluster_pair : *m_all_clusters) {
    cluster_charge = 0;
    cluster_time = cluster_pair.first;
    cluster_PE = 0;
    std::vector<Hit> cluster_hits = cluster_pair.second;
    for (int i = 0; i < (int) cluster_hits.size(); i++){
      int hit_ID = cluster_hits.at(i).GetTubeId();
      std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(hit_ID);
      if(it != ChannelKeyToSPEMap.end()){ //Charge to SPE conversion is available
        double hit_charge = cluster_hits.at(i).GetCharge();
        double hit_PE = hit_charge / ChannelKeyToSPEMap.at(hit_ID);
        cluster_charge+=hit_charge;
        cluster_PE+=hit_PE;
      } else {
        if(verbosity>2){
          std::cout << "FOUND A HIT FOR CHANNELKEY " << hit_ID << "BUT NO CONVERSION " <<
              "TO PE AVAILABLE.  SKIPPING PE." << std::endl;
        }
      }
    }
    if(verbosity>3) std::cout << "AmBeRunStatistics Tool: cluster time,charge: : " << cluster_time << "," << cluster_charge << std::endl;
    h_Cluster_ChargeCleanPromptTrig->Fill(cluster_charge);
    h_Cluster_TimeMeanCleanPromptTrig->Fill(cluster_time);
    h_Cluster_MultiplicityCleanPromptTrig->Fill(m_all_clusters->size());
    //See if this cluster is a valid neutron candidate based on input criteria
    if((cluster_PE > ClusterPEMin) && (cluster_PE < ClusterPEMax) && (cluster_time > SWindowMax)){
      h_Cluster_ChargeNeutronCandidate->Fill(cluster_charge);
      h_Cluster_PENeutronCandidate->Fill(cluster_PE);
      h_Cluster_TimeMeanNeutronCandidate->Fill(cluster_time);
      h_Cluster_MultiplicityNeutronCandidate->Fill(m_all_clusters->size());
    }

  }

  bool GoldenNeutronCandidate = true;
  //Now, see if there's only one cluster with a charge relatively consistent with a neutron
  //FIXME: Need to convert to PE count
  if(m_all_clusters->size() != 1 || ((cluster_PE < ClusterPEMin) || (cluster_PE > ClusterPEMax))
          || (cluster_time < SWindowMax)) {
    //Is not a valid single cluster neutron candidate event
    GoldenNeutronCandidate = false;
  }

  if(!GoldenNeutronCandidate) return true;
  if(verbosity>3) std::cout << "AmBeRunStatistics Tool: Found a golden neutron candidate trigger" << std::endl;

  NumberGoldenNeutronCandidates+=1;
  h_Cluster_ChargeGoldenCandidate->Fill(cluster_charge);
  h_Cluster_TimeMeanGoldenCandidate->Fill(cluster_time);
  h_Cluster_MultiplicityGoldenCandidate->Fill(m_all_clusters->size());
  h_Cluster_PEGoldenCandidate->Fill(cluster_PE);
  h_SiPM1_AmplitudeGoldenCandidate->Fill(SiPM1_MaxPulse.amplitude());
  h_SiPM2_AmplitudeGoldenCandidate->Fill(SiPM2_MaxPulse.amplitude());
  h_SiPM1_ChargeGoldenCandidate->Fill(SiPM1_MaxPulse.charge());
  h_SiPM2_ChargeGoldenCandidate->Fill(SiPM2_MaxPulse.charge());
  h_SiPM1_PeakTimeGoldenCandidate->Fill(SiPM1_MaxPulse.peak_time());
  h_SiPM2_PeakTimeGoldenCandidate->Fill(SiPM2_MaxPulse.peak_time());
  h_S1S2_AmplitudesGoldenCandidate->Fill(SiPM1_MaxPulse.amplitude(),SiPM2_MaxPulse.amplitude());
  h_S1S2_ChargeratioGoldenCandidate->Fill(charge_ratio);
  h_S1S2_DeltatGoldenCandidate->Fill(deltat);

  return true;
}


bool AmBeRunStatistics::Finalise(){

  std::cout << "AmBeStatistics tool: Saving ROOT File..." << std::endl;
  
  this->WriteHistograms();
  ambe_file_out->Close();

  std::cout << "AmBeStatistics tool: Number of events: " << NumberOfEvents << std::endl;
  std::cout << "AmBeStatistics tool: Number of clean prompt triggers: " << NumberOfCleanTriggers << std::endl;
  std::cout << "AmBeStatistics tool: Number of clean neutron candidate events: " << NumberGoldenNeutronCandidates << std::endl;
  return true;
}

void AmBeRunStatistics::GetSiPMPulseInfo(ADCPulse& SiPM1_Pulse, 
        ADCPulse& SiPM2_Pulse, int& SiPM1_NumTriggers, int& SiPM2_NumTriggers){
  int NumS2Triggers = 0;
  int NumS1Triggers = 0;
  double S1MaxAmplitude = 0;
  double S2MaxAmplitude = 0;

  //Calibrate the SIPM waveforms
  for (const auto& temp_pair : aux_pulse_map) {
    const auto& channel_key = temp_pair.first;
    //For now, only calibrate the SiPM waveforms
    if(AuxChannelNumToTypeMap->at(channel_key) != "SiPM1" &&
       AuxChannelNumToTypeMap->at(channel_key) != "SiPM2") continue;
    std::string channel_type = AuxChannelNumToTypeMap->at(channel_key);

    if(verbosity>4) {
      std::cout << "AmBeRunStatistics tool: Found SiPM channel " << 
        AuxChannelNumToTypeMap->at(channel_key) << "with channel key " << channel_key << std::endl;
    }
    std::vector< std::vector<ADCPulse>> sipm_minibuffers = temp_pair.second;
    size_t num_minibuffers = sipm_minibuffers.size();  //Should be size 1 in FrankDAQ mode
    for (size_t mb = 0; mb < num_minibuffers; ++mb) {
      std::vector<ADCPulse> thisbuffer_pulses = sipm_minibuffers.at(mb);
      if(verbosity>4) {
        std::cout << "AmBeRunStatistics tool: Number of pulses: " << 
         thisbuffer_pulses.size()  << std::endl;
      }
      for (size_t i = 0; i < thisbuffer_pulses.size(); i++){
        ADCPulse apulse = thisbuffer_pulses.at(i);
        double amplitude = apulse.amplitude();
        if(verbosity>4) {
          std::cout << "AmBeRunStatistics tool: Pulse amplitude is " << 
            amplitude << std::endl;
        }
        if(channel_type == "SiPM1"){
          if (amplitude>S1MaxAmplitude){
            S1MaxAmplitude = amplitude;
            SiPM1_Pulse = apulse;
          }
          if (amplitude*1000>S1Threshold) NumS1Triggers+=1;
        }
        if(channel_type == "SiPM2"){
          if (amplitude>S2MaxAmplitude){
            S2MaxAmplitude = amplitude;
            SiPM2_Pulse = apulse;
          }
          if (amplitude*1000>S2Threshold) NumS2Triggers+=1;
        }
      }
    }
  }
  SiPM1_NumTriggers = NumS1Triggers; 
  SiPM2_NumTriggers = NumS2Triggers; 
  return;
}

void AmBeRunStatistics::InitializeHistograms(){
  //All Event Histograms
  h_SiPM1_Amplitude = new TH1F("h_SiPM1_Amplitude","Max amplitude of SiPM1 in all events",300,0,0.3);
  h_SiPM2_Amplitude = new TH1F("h_SiPM2_Amplitude","Max amplitude of SiPM2 in all events",300,0,0.3);
  h_SiPM1_Charge = new TH1F("h_SiPM1_Charge","Total charge of SiPM1 in all events",300,0,3);
  h_SiPM2_Charge = new TH1F("h_SiPM2_Charge","Total charge of SiPM2 in all events",300,0,3);
  h_S1S2_Amplitudes = new TH2F("h_S1S2_Amplitudes","Max amplitude of SiPM1 and SiPM2 in all events",50,0,0.5,50,0,0.5);
  h_S1S2_Deltat = new TH1F("h_S1S2_Deltat","Time difference between peaks of SiPM1 and SiPM2 in all events (S1-S2)",200,-100,100);
  h_S1S2_Chargeratio = new TH1F("h_S1S2_Chargeratio","SiPM1/SiPM2 charge ratio in all events",50,0.5,1.5);
  h_SiPM1_PeakTime = new TH1F("h_SiPM1_PeakTime","Peak time for largest SiPM1 pulse in all events",4000,0,4000);
  h_SiPM2_PeakTime = new TH1F("h_SiPM2_PeakTime","Peak time for largest SiPM2 pulse in all events",4000,0,4000);
  //Valid AmBe Trigger Histograms
  h_Cluster_ChargeCleanPromptTrig = new TH1F("h_Cluster_ChargeCleanPromptTrig","Cluster charges for events with a clean prompt trigger",1000,0,1);
  h_Cluster_TimeMeanCleanPromptTrig = new TH1F("h_Cluster_TimeMeanCleanPromptTrig","Mean cluster time for events with a clean prompt trigger",80000,0,80000);
  h_Cluster_MultiplicityCleanPromptTrig = new TH1F("h_Cluster_MultiplicityCleanPromptTrig","Multiplicity of clusters for events with a clean prompt trigger",22,-2,20);
  h_SiPM1_AmplitudeCleanPromptTrig = new TH1F("h_SiPM1_AmplitudeCleanPromptTrig","Max amplitude of SiPM1 in clean prompt triggers", 50, 0, 0.5);
  h_SiPM2_AmplitudeCleanPromptTrig = new TH1F("h_SiPM2_AmplitudeCleanPromptTrig","Max amplitude of SiPM2 in clean prompt triggers", 50, 0, 0.5);
  h_SiPM1_ChargeCleanPromptTrig = new TH1F("h_SiPM1_ChargeCleanPromptTrig","Total charge of SiPM1 pulse in clean prompt triggers",300,0,3);
  h_SiPM2_ChargeCleanPromptTrig = new TH1F("h_SiPM2_ChargeCleanPromptTrig","Total charge of SiPM2 pulse in clean prompt triggers",300,0,3);
  h_SiPM1_PeakTimeCleanPromptTrig = new TH1F("h_SiPM1_PeakTimeCleanPromptTrig","Peak time for largest SiPM1 pulse in clean prompt triggers",4000,0,4000);
  h_SiPM2_PeakTimeCleanPromptTrig = new TH1F("h_SiPM2_PeakTimeCleanPromptTrig","Peak time for largest SiPM2 pulse in clean prompt triggers",4000,0,4000);
  h_S1S2_AmplitudesCleanPromptTrig = new TH2F("h_S1S2_AmplitudesCleanPromptTrig","Max amplitude of SiPM1 and SiPM2 in clean prompt triggers",50,0,0.5,50,0,0.5);
  h_S1S2_ChargeratioCleanPromptTrig = new TH1F("h_S1S2_ChargeratioCleanPromptTrig","SiPM1/SiPM2 charge ratio in clean prompt triggers",50,0.5,1.5);
  h_S1S2_DeltatCleanPromptTrig = new TH1F("h_S1S2_DeltatCleanPromptTrig","Time difference between peaks of SiPM1 and SiPM2 (S1-S2) in clean prompt triggers",200,-100,100);
  h_Cluster_ChargeNeutronCandidate = new TH1F("h_Cluster_ChargeNeutronCandidate","Cluster charges for neutron candidate",1000,0,1);
  h_Cluster_TimeMeanNeutronCandidate = new TH1F("h_Cluster_TimeMeanNeutronCandidate","Mean cluster time for neutron candidates",80000,0,80000);
  h_Cluster_MultiplicityNeutronCandidate = new TH1F("h_Cluster_MultiplicityNeutronCandidate","Multiplicity of neutron candidate clusters",22,-2,20);
  h_Cluster_PENeutronCandidate = new TH1F("h_Cluster_PENeutronCandidate","Total PE for neutron candidates",300,0,300);
  //Valid Golden Neutron Candidate Histograms (one cluster only)
  h_Cluster_ChargeGoldenCandidate = new TH1F("h_Cluster_ChargeGoldenCandidate","Cluster charges for events with a golden neutron candidate",1000,0,1);
  h_Cluster_TimeMeanGoldenCandidate = new TH1F("h_Cluster_TimeMeanGoldenCandidate","Mean cluster time for events with a golden neutron candidate",80000,0,80000);
  h_Cluster_MultiplicityGoldenCandidate = new TH1F("h_Cluster_MultiplicityGoldenCandidate","Multiplicity of clusters for events with a golden neutron candidate",22,-2,20);
  h_Cluster_PEGoldenCandidate = new TH1F("h_Cluster_PEGoldenCandidate","Total PE for golden neutron candidates",(ClusterPEMax-ClusterPEMin),ClusterPEMin,ClusterPEMax);
  h_SiPM1_AmplitudeGoldenCandidate = new TH1F("h_SiPM1_AmplitudeGoldenCandidate","Max amplitude of SiPM1 in golden neutron candidates", 300, 0, 0.3);
  h_SiPM2_AmplitudeGoldenCandidate = new TH1F("h_SiPM2_AmplitudeGoldenCandidate","Max amplitude of SiPM2 in golden neutron candidates", 300, 0, 0.3);
  h_SiPM1_ChargeGoldenCandidate = new TH1F("h_SiPM1_ChargeGoldenCandidate","Total charge of SiPM1 pulse in golden neutron candidates",300,0,3);
  h_SiPM2_ChargeGoldenCandidate = new TH1F("h_SiPM2_ChargeGoldenCandidate","Total charge of SiPM2 pulse in golden neutron candidates",300,0,3);
  h_SiPM1_PeakTimeGoldenCandidate = new TH1F("h_SiPM1_PeakTimeGoldenCandidate","Peak time for largest SiPM1 pulse in golden neutron candidates",4000,0,4000);
  h_SiPM2_PeakTimeGoldenCandidate = new TH1F("h_SiPM2_PeakTimeGoldenCandidate","Peak time for largest SiPM2 pulse in golden neutron candidates",4000,0,4000);
  h_S1S2_AmplitudesGoldenCandidate = new TH2F("h_S1S2_AmplitudesGoldenCandidate","Max amplitude of SiPM1 and SiPM2 in golden neutron candidates",50,0,0.5,50,0,0.5);
  h_S1S2_ChargeratioGoldenCandidate = new TH1F("h_S1S2_ChargeratioGoldenCandidate","SiPM1/SiPM2 charge ratio in golden neutron candidates",50,0.5,1.5);
  h_S1S2_DeltatGoldenCandidate = new TH1F("h_S1S2_DeltatGoldenCandidate","Time difference between peaks of SiPM1 and SiPM2 (S1-S2) in golden neutron candidates",200,-100,100);
}

void AmBeRunStatistics::WriteHistograms(){
  ambe_file_out->cd();
  if(verbosity>2) std::cout << "AmBeRunStatistics tool: Making AllHist dir"  << std::endl;
  std::string AllHist = "AllEvents";
  TDirectory* allhistdir = ambe_file_out->mkdir(AllHist.c_str());
  allhistdir->cd();
  if(verbosity>2) std::cout << "AmBeRunStatistics tool: Writing all event histograms"  << std::endl;
  h_SiPM1_Amplitude->Write();
  h_SiPM2_Amplitude->Write();
  h_SiPM1_Charge->Write();
  h_SiPM2_Charge->Write();
  h_SiPM1_PeakTime->Write();
  h_SiPM2_PeakTime->Write();
  h_S1S2_Amplitudes->Write();
  h_S1S2_Chargeratio->Write();
  h_S1S2_Deltat->Write();

  if(verbosity>2) std::cout << "AmBeRunStatistics tool: Writing valid AmBe trigger histograms"  << std::endl;
  std::string AmBePromptHist = "ValidAmBeTrigger";
  TDirectory* ambehistdir = ambe_file_out->mkdir(AmBePromptHist.c_str());
  ambehistdir->cd();
  h_Cluster_ChargeCleanPromptTrig->Write();
  h_Cluster_TimeMeanCleanPromptTrig->Write();
  h_Cluster_MultiplicityCleanPromptTrig->Write();
  h_SiPM1_AmplitudeCleanPromptTrig->Write();
  h_SiPM2_AmplitudeCleanPromptTrig->Write();
  h_SiPM1_ChargeCleanPromptTrig->Write();
  h_SiPM2_ChargeCleanPromptTrig->Write();
  h_SiPM1_PeakTimeCleanPromptTrig->Write();
  h_SiPM2_PeakTimeCleanPromptTrig->Write();
  h_S1S2_AmplitudesCleanPromptTrig->Write();
  h_S1S2_ChargeratioCleanPromptTrig->Write();
  h_S1S2_DeltatCleanPromptTrig->Write();
  h_Cluster_ChargeNeutronCandidate->Write();
  h_Cluster_PENeutronCandidate->Write();
  h_Cluster_TimeMeanNeutronCandidate->Write();
  h_Cluster_MultiplicityNeutronCandidate->Write();
  
  std::string GoldenNeutronCandidHist = "GoldenNeutronCandidate";
  TDirectory* goldneuthistdir = ambe_file_out->mkdir(GoldenNeutronCandidHist.c_str());
  goldneuthistdir->cd();
  h_Cluster_ChargeGoldenCandidate->Write();
  h_Cluster_TimeMeanGoldenCandidate->Write();
  h_Cluster_MultiplicityGoldenCandidate->Write();
  h_Cluster_PEGoldenCandidate->Write();
  h_SiPM1_AmplitudeGoldenCandidate->Write();
  h_SiPM2_AmplitudeGoldenCandidate->Write();
  h_SiPM1_ChargeGoldenCandidate->Write();
  h_SiPM2_ChargeGoldenCandidate->Write();
  h_SiPM1_PeakTimeGoldenCandidate->Write();
  h_SiPM2_PeakTimeGoldenCandidate->Write();
  h_S1S2_AmplitudesGoldenCandidate->Write();
  h_S1S2_ChargeratioGoldenCandidate->Write();
  h_S1S2_DeltatGoldenCandidate->Write();
  ambe_file_out->cd();
  return;
}

void AmBeRunStatistics::SetHistogramLabels(){
  h_SiPM1_Amplitude->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM2_Amplitude->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM1_Charge->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM2_Charge->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM1_PeakTime->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_SiPM2_PeakTime->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_S1S2_Amplitudes->GetXaxis()->SetTitle("SiPM1 Pulse amplitude [V]");
  h_S1S2_Amplitudes->GetYaxis()->SetTitle("SiPM2 Pulse amplitude [V]");
  h_S1S2_Chargeratio->GetXaxis()->SetTitle("S1/S2 charge ratio");
  h_S1S2_Deltat->GetXaxis()->SetTitle("S1-S2 peak time [ns]");

  h_Cluster_ChargeCleanPromptTrig->GetXaxis()->SetTitle("Cluster charge [nC]");
  h_Cluster_TimeMeanCleanPromptTrig->GetXaxis()->SetTitle("Cluster hit time mean [ns]");
  h_Cluster_MultiplicityCleanPromptTrig->GetXaxis()->SetTitle("Cluster multiplicity");
  h_SiPM1_AmplitudeCleanPromptTrig->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM2_AmplitudeCleanPromptTrig->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM1_ChargeCleanPromptTrig->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM2_ChargeCleanPromptTrig->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM1_PeakTimeCleanPromptTrig->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_SiPM2_PeakTimeCleanPromptTrig->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_S1S2_AmplitudesCleanPromptTrig->GetXaxis()->SetTitle("SiPM1 Pulse amplitude [V]");
  h_S1S2_AmplitudesCleanPromptTrig->GetYaxis()->SetTitle("SiPM2 Pulse amplitude [V]");
  h_S1S2_ChargeratioCleanPromptTrig->GetXaxis()->SetTitle("S1/S2 charge ratio");
  h_S1S2_DeltatCleanPromptTrig->GetXaxis()->SetTitle("S1-S2 peak time [ns]");
  h_Cluster_ChargeNeutronCandidate->GetXaxis()->SetTitle("Cluster charge [nC]");
  h_Cluster_TimeMeanNeutronCandidate->GetXaxis()->SetTitle("Cluster hit time mean [ns]");
  h_Cluster_MultiplicityNeutronCandidate->GetXaxis()->SetTitle("Cluster multiplicity");
  
  h_Cluster_ChargeGoldenCandidate->GetXaxis()->SetTitle("Cluster charge [nC]");
  h_Cluster_TimeMeanGoldenCandidate->GetXaxis()->SetTitle("Cluster hit time mean [ns]");
  h_Cluster_MultiplicityGoldenCandidate->GetXaxis()->SetTitle("Cluster multiplicity");
  h_SiPM1_AmplitudeGoldenCandidate->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM2_AmplitudeGoldenCandidate->GetXaxis()->SetTitle("Pulse amplitude [V]");
  h_SiPM1_ChargeGoldenCandidate->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM2_ChargeGoldenCandidate->GetXaxis()->SetTitle("Pulse charge [nC]");
  h_SiPM1_PeakTimeGoldenCandidate->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_SiPM2_PeakTimeGoldenCandidate->GetXaxis()->SetTitle("Pulse peak time [ns]");
  h_S1S2_AmplitudesGoldenCandidate->GetXaxis()->SetTitle("SiPM1 Pulse amplitude [V]");
  h_S1S2_AmplitudesGoldenCandidate->GetYaxis()->SetTitle("SiPM2 Pulse amplitude [V]");
  h_S1S2_ChargeratioGoldenCandidate->GetXaxis()->SetTitle("S1/S2 charge ratio");
  h_S1S2_DeltatGoldenCandidate->GetXaxis()->SetTitle("S1-S2 peak time [ns]");
  return;
}
