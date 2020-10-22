#include "EventSelector.h"

EventSelector::EventSelector():Tool(){}


bool EventSelector::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  fPMTMRDOffset = false;
  fIsMC = true;
  fPMTMRDOffset = 745;

  //Get the tool configuration variables
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("MRDRecoCut", fMRDRecoCut);
  m_variables.Get("MCFVCut", fMCFVCut);
  m_variables.Get("MCPMTVolCut", fMCPMTVolCut);
  m_variables.Get("MCMRDCut", fMCMRDCut);
  m_variables.Get("MCPiKCut", fMCPiKCut);
  m_variables.Get("MCIsMuonCut", fMCIsMuonCut);
  m_variables.Get("MCIsElectronCut", fMCIsElectronCut);
  m_variables.Get("MCIsSingleRingCut", fMCIsSingleRingCut);
  m_variables.Get("MCIsMultiRingCut", fMCIsMultiRingCut);
  m_variables.Get("MCProjectedMRDHit", fMCProjectedMRDHit);
  m_variables.Get("MCEnergyCut", fMCEnergyCut);
  m_variables.Get("Emin",Emin);
  m_variables.Get("Emax",Emax);
  m_variables.Get("MCDSNBlike",fMCDSNBlike);
  m_variables.Get("NHitCut", fNHitCut);
  m_variables.Get("NHitmin", fNHitmin);
  m_variables.Get("PromptTrigOnly", fPromptTrigOnly);
  m_variables.Get("RecoFVCut", fRecoFVCut);
  m_variables.Get("RecoPMTVolCut", fRecoPMTVolCut);
  m_variables.Get("PMTMRDCoincCut",fPMTMRDCoincCut);
  m_variables.Get("PMTMRDOffset",fPMTMRDOffset);
  m_variables.Get("NoVeto",fNoVetoCut);
  m_variables.Get("Veto",fVetoCut);
  m_variables.Get("SaveStatusToStore", fSaveStatusToStore);
  m_variables.Get("IsMC",fIsMC);

  if (!fIsMC){fMCFVCut = false; fMCPMTVolCut = false; fMCMRDCut = false; fMCPiKCut = false; fMCIsMuonCut = false; fMCIsElectronCut = false; fMCIsSingleRingCut = false; fMCIsMultiRingCut = false; fMCProjectedMRDHit = false; fMCEnergyCut = false; fPromptTrigOnly = false;}

  /// Construct the other objects we'll be needing at event level,
  
  // Make the RecoEvent Store if it doesn't exist
  int recoeventexists = m_data->Stores.count("RecoEvent");
  if(recoeventexists==0) m_data->Stores["RecoEvent"] = new BoostStore(false,2);
  //TODO: Have an event selection mask filled based on what cuts are run
  //m_data->Stores.at("RecoEvent")->Set("EvtSelectionMask", fEvtSelectionMask); 
  
  auto get_geometry= m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",fGeometry);
  if(!get_geometry){
    Log("EventSelector Tool: Error retrieving Geometry from ANNIEEvent!",v_error,verbosity); 
    return false; 
  }

  vec_pmtclusters_charge = new std::vector<double>; 
  vec_pmtclusters_time = new std::vector<double>; 
  vec_mrdclusters_time = new std::vector<double>; 

  return true;
}


bool EventSelector::Execute(){
  // Reset everything
  this->Reset();
  
  // see if "ANNIEEvent" exists
  auto get_annieevent = m_data->Stores.count("ANNIEEvent");
  if(!get_annieevent){
  	Log("EventSelector Tool: No ANNIEEvent store!",v_error,verbosity); 
  	return false;
  };

  // ANNIE Event number
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
  
  if (fIsMC){
  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",fMCEventNum);  
  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggernum); 
  
  std::string logmessage = "EventSelector Tool: Processing MCEntry "+to_string(fMCEventNum)+
  ", MCTrigger "+to_string(fMCTriggernum) + ", Event "+to_string(fEventNumber);
  Log(logmessage,v_message,verbosity);
  }

  // Retrive digits from RecoEvent
  auto get_ok = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(not get_ok){
  	Log("EventSelector  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
  	return false;
  }

  // BEGIN CUTS USING TRUTH INFORMATION //

  if (fIsMC){
    // get truth vertex information 
    auto get_truevtx = m_data->Stores.at("RecoEvent")->Get("TrueVertex", fMuonStartVertex);
    if(!get_truevtx){ 
      Log("EventSelector Tool: Error retrieving TrueVertex from RecoEvent!",v_error,verbosity); 
      return false; 
    }
  
    auto get_truestopvtx = m_data->Stores.at("RecoEvent")->Get("TrueStopVertex", fMuonStopVertex);
    if(!get_truestopvtx){ 
      Log("EventSelector Tool: Error retrieving TrueStopVertex from RecoEvent!",v_error,verbosity); 
      return false; 
    }
   
    //Get MC version of MRD hits
    bool get_mrd = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData_MC);
    if (!get_mrd) {
      Log("EventSelector Tool: Error retrieving TDCData, true from ANNIEEvent!",v_error,verbosity);
      return false;
    }
  } else {
  
    //Get data version of MRD hits
    bool get_mrd = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
    if (!get_mrd) {
      Log("EventSelector Tool: Error retrieving TDCData, true from ANNIEEvent!",v_error,verbosity);
      return false;
    }

  }

  bool IsSingleRing = false, IsMultiRing = false, HasProjectedMRDHit = false, passNoPiK = false, passMCFVCut = false, passMCPMTCut = false, passMCMRDCut = false, IsInsideEnergyWindow = false, IsElectron = false, IsMuon = false, isPromptTrigger=false;
 
  if (fIsMC){

    isPromptTrigger = this->PromptTriggerCheck();
    m_data->Stores.at("RecoEvent")->Set("PromptEvent",isPromptTrigger);
    
    IsSingleRing = this->EventSelectionByMCSingleRing();
    m_data->Stores.at("RecoEvent")->Set("MCSingleRingEvent",IsSingleRing);

    IsMultiRing = this->EventSelectionByMCMultiRing();
    m_data->Stores.at("RecoEvent")->Set("MCMultiRingEvent",IsMultiRing);

    HasProjectedMRDHit = this->EventSelectionByMCProjectedMRDHit();
    //information about projected MRD hit already stored in the RecoEvent store by MCRecoEventLoader

    passNoPiK = this->EventSelectionNoPiK();
    m_data->Stores.at("RecoEvent")->Set("MCNoPiK",passNoPiK);

    passMCFVCut = this->EventSelectionByFV(true);
    m_data->Stores.at("RecoEvent")->Set("MCFV",passMCFVCut);

    passMCPMTCut = this->EventSelectionByPMTVol(true);
    m_data->Stores.at("RecoEvent")->Set("MCPMTVol",passMCPMTCut);
 
    passMCMRDCut= this->EventSelectionByMCTruthMRD();
    m_data->Stores.at("RecoEvent")->Set("MCMRDStop",passMCMRDCut);

    IsInsideEnergyWindow = this->EnergyCutCheck(Emin,Emax);
    m_data->Stores.at("RecoEvent")->Set("MCEnergyCut",IsInsideEnergyWindow);

    IsMuon = this->ParticleCheck(13);
    m_data->Stores.at("RecoEvent")->Set("MCIsMuon",IsMuon);

    IsElectron = this->ParticleCheck(11);
    m_data->Stores.at("RecoEvent")->Set("MCIsElectron",IsElectron);

  }


  bool HasEnoughHits = this->NHitCountCheck(fNHitmin);
  m_data->Stores.at("RecoEvent")->Set("NHitCut",HasEnoughHits);  

  bool passPMTMRDCoincCut = this->EventSelectionByPMTMRDCoinc();
  m_data->Stores.at("RecoEvent")->Set("PMTMRDCoinc",passPMTMRDCoincCut);

  bool passVetoCut = this->EventSelectionByVetoCut();
  m_data->Stores.at("RecoEvent")->Set("NoVeto",passVetoCut);

  // Fill the EventSelection mask for the cuts that are supposed to be applied
  if (fMCPiKCut){
    fEventApplied |= EventSelector::kFlagMCPiK; 
    if(!passNoPiK) fEventFlagged |= EventSelector::kFlagMCPiK;
  }  

  if(fMCFVCut || fMCPMTVolCut){
    if (fMCFVCut){
      fEventApplied |= EventSelector::kFlagMCFV; 
      if(!passMCFVCut) fEventFlagged |= EventSelector::kFlagMCFV;
    }
    if (fMCPMTVolCut){
      fEventApplied |= EventSelector::kFlagMCPMTVol; 
      if(!passMCPMTCut) fEventFlagged |= EventSelector::kFlagMCPMTVol;
    }
  }

  if(fMCMRDCut){
    fEventApplied |= EventSelector::kFlagMCMRD; 
    if(!passMCMRDCut) fEventFlagged |= EventSelector::kFlagMCMRD;
  }

  if(fPromptTrigOnly){
    fEventApplied |= EventSelector::kFlagPromptTrig; 
    if(!isPromptTrigger) fEventFlagged |= EventSelector::kFlagPromptTrig;
  }

  if(fNHitCut){
    fEventApplied |= EventSelector::kFlagNHit;
    if(!HasEnoughHits) fEventFlagged |= EventSelector::kFlagNHit;
  }

  if (fMCEnergyCut){
    fEventApplied |= EventSelector::kFlagMCEnergyCut;
    if (!IsInsideEnergyWindow) fEventFlagged |= EventSelector::kFlagMCEnergyCut;
  }

  if (fMCIsMuonCut){
    fEventApplied |= EventSelector::kFlagMCIsMuon;
    if (!IsMuon) fEventFlagged |= EventSelector::kFlagMCIsMuon;
  }

  if (fMCIsElectronCut){
    fEventApplied |= EventSelector::kFlagMCIsElectron;
    if (!IsElectron) fEventFlagged |= EventSelector::kFlagMCIsElectron;
  }

  if (fMCIsSingleRingCut){
    fEventApplied |= EventSelector::kFlagMCIsSingleRing;
    if (!IsSingleRing) fEventFlagged |= EventSelector::kFlagMCIsSingleRing;
  }

  if (fMCIsMultiRingCut){
    fEventApplied |= EventSelector::kFlagMCIsMultiRing;
    if (!IsMultiRing) fEventFlagged |= EventSelector::kFlagMCIsMultiRing;
  }

  if (fMCProjectedMRDHit){
    fEventApplied |= EventSelector::kFlagMCProjectedMRDHit;
    if (!HasProjectedMRDHit) fEventFlagged |= EventSelector::kFlagMCProjectedMRDHit;
  }

  if (fMCDSNBlike){
    bool IsDSNBlike = this->DSNBCheck();
    fEventApplied |= EventSelector::kFlagMCDSNBlike;
    if (!IsDSNBlike) fEventFlagged |= EventSelector::kFlagMCDSNBlike;
  }

  // BEGIN CUTS USING RECONSTRUCTED INFORMATION //
  if(fRecoFVCut || fRecoPMTVolCut){
	// Retrive Reconstructed vertex from RecoEvent 
	auto get_ok = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex",fRecoVertex);  ///> Get reconstructed vertex 
    if(not get_ok){
  	  Log("EventSelector Tool: Error retrieving Extended vertex from RecoEvent!",v_error,verbosity); 
  	  return false;
    }
    if(fRecoFVCut){
      fEventApplied |= EventSelector::kFlagRecoFV; 
      bool passRecoFVCut= this->EventSelectionByFV(false);
      if(!passRecoFVCut) fEventFlagged |= EventSelector::kFlagRecoFV;
    }
    if(fRecoPMTVolCut){
      fEventApplied |= EventSelector::kFlagRecoPMTVol; 
      bool passRecoPMTVolCut= this->EventSelectionByPMTVol(false);
      if(!passRecoPMTVolCut) fEventFlagged |= EventSelector::kFlagRecoPMTVol;
    }
  }

  //FIXME: This isn't working according to Jingbo
  if(fMRDRecoCut){
    fEventApplied |= EventSelector::kFlagRecoMRD; 
    std::cout << "EventSelector Tool: Currently not implemented. Setting to false" << std::endl;
    Log("EventSelector Tool: MRDReco not implemented.  Setting cut bit to false",v_message,verbosity);
    bool passMRDRecoCut = false;
    //bool passMRDRecoCut = this->EventSelectionByMRDReco(); 
    if(!passMRDRecoCut) fEventFlagged |= EventSelector::kFlagRecoMRD; 
  }

  //Fast check whether the times of MRD and PMT clusters are coincident
  if(fPMTMRDCoincCut){
    fEventApplied |= EventSelector::kFlagPMTMRDCoinc;
    if (!passPMTMRDCoincCut) fEventFlagged |= EventSelector::kFlagPMTMRDCoinc;
  }

  if (fNoVetoCut){
    fEventApplied |= EventSelector::kFlagNoVeto;
    if (!passVetoCut) fEventFlagged |= EventSelector::kFlagNoVeto;
  }

  if (fVetoCut){
    fEventApplied |= EventSelector::kFlagVeto;
    if (passVetoCut) fEventFlagged |= EventSelector::kFlagVeto;
  }
  
  if(fEventFlagged != EventSelector::kFlagNone) fEventCutStatus = false;
  if(fEventCutStatus){  
    Log("EventSelector Tool: Event is clean according to current event selection.",v_message,verbosity);
  }
  if(fSaveStatusToStore) m_data->Stores.at("RecoEvent")->Set("EventCutStatus", fEventCutStatus);
  m_data->Stores.at("RecoEvent")->Set("EventFlagApplied", fEventApplied);
  m_data->Stores.at("RecoEvent")->Set("EventFlagged", fEventFlagged);

  if (verbosity >= v_debug){
    std::cout << "EventSelector tool: fEventApplied: "<< fEventApplied << ", fEventFlagged: " << fEventFlagged << std::endl;
    std::cout << "EventSelector tool: Bit representation: fEventApplied: " << std::bitset<32>(fEventApplied) << ", fEventFlagged: " << std::bitset<32>(fEventFlagged) << std::endl;
  }

  return true;
}


bool EventSelector::Finalise(){
  if(verbosity>0) cout<<"EventSelector exitting"<<endl;
  delete vec_pmtclusters_charge;
  delete vec_pmtclusters_time;
  delete vec_mrdclusters_time;
  return true;
}

bool EventSelector::EventSelectionNoPiK() {
  
  // Get the pion and kaon counts from the store.  If any count is greater
  // than zero, the fEventCutStatus is set to "false" for a failed event.
  //Fill in pion counts for this event
  int pi0count, pipcount, pi0mcount, K0count, Kpcount, Kmcount;
  m_data->Stores.at("RecoEvent")->Get("MCPi0Count", pi0count);
  m_data->Stores.at("RecoEvent")->Get("MCPiPlusCount", pipcount);
  m_data->Stores.at("RecoEvent")->Get("MCPiMinusCount", pi0mcount);
  m_data->Stores.at("RecoEvent")->Get("MCK0Count", K0count);
  m_data->Stores.at("RecoEvent")->Get("MCKPlusCount", Kpcount);
  m_data->Stores.at("RecoEvent")->Get("MCKMinusCount", Kmcount);
  int sum = pi0count + pipcount + pi0mcount + K0count + Kpcount + Kmcount;
  if (sum > 0){
    Log("EventSelector: A primary pion or kaon was found. Total count: " +to_string(sum),v_message,verbosity);
    return false;
  } else{
    return true;
  }
}

bool EventSelector::PromptTriggerCheck() {
  /// First, see if this is a delayed trigger in the event
  auto get_mctrigger = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggernum);
  if(!get_mctrigger){ 
    Log("EventSelector Tool: Error retrieving MCTriggernum from ANNIEEvent!",v_error,verbosity); 
    return false; 
  }	
  /// if so, truth analysis is probably not interested in this trigger. Primary muon will not be in the listed tracks.
  if(fMCTriggernum>0){ 
    Log("EventSelector Tool: This event is not a prompt trigger",v_message,verbosity); 
    return false;
  }
  return true;
}

bool EventSelector::NHitCountCheck(int NHitCut) {
  if(verbosity>3){
    std::cout << "SIZE OF DIGIT LIST FOR EVENT IS: " << fDigitList->size() << std::endl;
  }
  
  int size_pmt_digits = 0;

  for (unsigned int i_digit = 0; i_digit < fDigitList->size(); i_digit++){
    RecoDigit thisdigit = fDigitList->at(i_digit);
    int digittype = thisdigit.GetDigitType();
    if (digittype == 0) size_pmt_digits++;
  }

  if(size_pmt_digits<NHitCut) {
    Log("EventSelector Tool: Event has less than 4 digits",v_message,verbosity);
    return false;
  }
  else {
  return true;
  }
}

// This function isn't working now, because the MRDTracks store must have been changed. 
// We have to contact Marcus and ask how we can retieve the MRD track information. 
bool EventSelector::EventSelectionByMRDReco() {
  /// Get number of subentries
  int NumMrdTracks = 0;
  double mrdtracklength = 0.;
  double MaximumMRDTrackLength = 0;
  int mrdlayershit = 0;
  bool mrdtrackisstoped = false;
  int longesttrackEntryNumber = -9999;
  
  /// check if "MRDTracks" store exists
	int mrdtrackexists = m_data->Stores.count("MRDTracks");
	if(mrdtrackexists == 0) {
	  Log("EventSelector Tool: MRDTracks store doesn't exist!",v_message,verbosity);
	  return false;
	} 	
  
  auto getmrdtracks = m_data->Stores.at("MRDTracks")->Get("NumMrdTracks",NumMrdTracks);
  if(!getmrdtracks) {
    Log("EventSelector Tool: Error retrieving MRDTracks Store!",v_error,verbosity); 
  	return false;	
  }
  logmessage = "EventSelector Tool: Found " + to_string(NumMrdTracks) + " MRD tracks";
  Log(logmessage,v_message,verbosity);
  
  // If mrd track isn't found
  if(NumMrdTracks == 0) {
    Log("EventSelector Tool: Found no MRD Tracks",v_message,verbosity);
    return false;
  }
  // If mrd track is found
  // MRD tracks are saved in the "MRDTracks" Store. Each entry is an MRD track
  
  for(int entrynum=0; entrynum<NumMrdTracks; entrynum++) {
    //m_data->Stores["MRDTracks"]->GetEntry(entrynum);
    m_data->Stores.at("MRDTracks")->Get("TrackLength",mrdtracklength);
    m_data->Stores.at("MRDTracks")->Get("LayersHit",mrdlayershit);
    if(MaximumMRDTrackLength<mrdtracklength) {
      MaximumMRDTrackLength = mrdtracklength;
      longesttrackEntryNumber = entrynum;
    }
  }
  // check if the longest track is stopped inside MRD
  //m_data->Stores["MRDTracks"]->GetEntry(longesttrackEntryNumber);
  m_data->Stores.at("MRDTracks")->Get("IsMrdStopped",mrdtrackisstoped);
  if(!mrdtrackisstoped) return false;
  return true;
}


bool EventSelector::EventSelectionByFV(bool isMC) {
  if(isMC && !fMuonStartVertex) return false;
  if(!isMC && !fRecoVertex) return false;
  RecoVertex* checkedVertex;
  if(isMC){
      Log("EventSelector Tool: Checking FV cut for true muon vertex",v_debug,verbosity); 
      checkedVertex=fMuonStartVertex;
  } else {
      Log("EventSelector Tool: Checking FV cut for reconstructed muon vertex",v_debug,verbosity); 
    checkedVertex=fRecoVertex;
  }
  double checkedVtxX, checkedVtxY, checkedVtxZ;
  Position vtxPos = checkedVertex->GetPosition();
  Direction vtxDir = checkedVertex->GetDirection();
  checkedVtxX = vtxPos.X();
  checkedVtxY = vtxPos.Y();
  checkedVtxZ = vtxPos.Z();
  //std::cout<<"checkedVtxX, Y, Z = "<<checkedVtxX<<", "<<checkedVtxY<<", "<<checkedVtxZ<<std::endl;
  double tankradius = ANNIEGeometry::Instance()->GetCylRadius();
  double fidcutradius = 0.8 * tankradius;
  double fidcuty = 50.;
  double fidcutz = 0.;
  if(checkedVtxZ > fidcutz) return false;
  if( (TMath::Sqrt(TMath::Power(checkedVtxX, 2) + TMath::Power(checkedVtxZ,2)) > fidcutradius) 
  	  || (TMath::Abs(checkedVtxY) > fidcuty) 
  	  || (checkedVtxZ > fidcutz) ){
  Log("EventSelector Tool: This event is not contained inside the FV",v_message,verbosity); 
  return false;
  }	
 
  return true;
}

bool EventSelector::EventSelectionByPMTVol(bool isMC) {
  if(isMC && !fMuonStartVertex) return false;
  if(!isMC && !fRecoVertex) return false;
  RecoVertex* checkedVertex;
  if(isMC){
      Log("EventSelector Tool: Checking PMT volume cut for true muon vertex",v_debug,verbosity); 
      checkedVertex=fMuonStartVertex;
  } else {
      Log("EventSelector Tool: Checking PMT volume cut for reconstructed muon vertex",v_debug,verbosity); 
    checkedVertex=fRecoVertex;
  }
  double checkedVtxX, checkedVtxY, checkedVtxZ;
  Position vtxPos = checkedVertex->GetPosition();
  Direction vtxDir = checkedVertex->GetDirection();
  checkedVtxX = vtxPos.X();
  checkedVtxY = vtxPos.Y();
  checkedVtxZ = vtxPos.Z();
  double fidcutradius = fGeometry->GetPMTEnclosedRadius()*100.;
  double fidcuty = fGeometry->GetPMTEnclosedHalfheight()*100.;
  if( (TMath::Sqrt(TMath::Power(checkedVtxX, 2) + TMath::Power(checkedVtxZ,2)) > fidcutradius) 
  	  || (TMath::Abs(checkedVtxY) > fidcuty)){
  Log("EventSelector Tool: This event is not contained within the PMT volume",v_message,verbosity); 
  return false;
  }	
 
  return true;
}

bool EventSelector::EventSelectionByMCTruthMRD() {
  if(!fMuonStopVertex) return false;
  // mrd cut
  double muonStopX, muonStopY, muonStopZ;
  muonStopX = fMuonStopVertex->GetPosition().X();
  muonStopY = fMuonStopVertex->GetPosition().Y();
  muonStopZ = fMuonStopVertex->GetPosition().Z();
  double mrdStartZ = fGeometry->GetMrdStart()*100-168.1;
  double mrdEndZ = fGeometry->GetMrdEnd()*100-168.1;
  double mrdHeightY = fGeometry->GetMrdHeight()*100;                                                                                     
  double mrdWidthX = fGeometry->GetMrdWidth()*100;
  Log("EventSelector tool: Read in MuonStop (X,Y,Z) = ("+std::to_string(muonStopX)+","+std::to_string(muonStopY)+","+std::to_string(muonStopZ)+")");
  if(muonStopZ<mrdStartZ || muonStopZ>mrdEndZ
  	|| muonStopX<-1.0*mrdWidthX || muonStopX>mrdWidthX
  	|| muonStopY<-1.0*mrdHeightY || muonStopY>mrdHeightY) {
    Log("EventSelector Tool: This MC Event's muon does not stop in the MRD",v_message,verbosity); 
    return false;	
  }
  return true;
}

bool EventSelector::EnergyCutCheck(double Emin, double Emax) {

  double energy;
  m_data->Stores.at("RecoEvent")->Get("TrueMuonEnergy",energy);
  if (energy > Emax || energy < Emin) return false;
  return true;

}

bool EventSelector::ParticleCheck(int pdg_number) {

  int pdg;
  m_data->Stores.at("RecoEvent")->Get("PdgPrimary",pdg);
  if (pdg!=pdg_number) return false;
  return true;

}

bool EventSelector::EventSelectionByMCSingleRing() {

  int nrings;
  m_data->Stores.at("RecoEvent")->Get("NRings",nrings);
  if (nrings!=1) return false;
  return true;

}

bool EventSelector::EventSelectionByMCMultiRing() {

  int nrings;
  m_data->Stores.at("RecoEvent")->Get("NRings",nrings);
  if (nrings <=1 ) return false;
  return true;

}

bool EventSelector::EventSelectionByMCProjectedMRDHit() {

  bool has_projected_mrd_hit;
  m_data->Stores.at("RecoEvent")->Get("ProjectedMRDHit",has_projected_mrd_hit);
  return has_projected_mrd_hit;

}

bool EventSelector::DSNBCheck() {
  
  bool is_dsnb_like = false;
  int neutron_count;
  int gamma_count;
  int positron_count;
  m_data->Stores.at("RecoEvent")->Get("NeutronCount",neutron_count);
  m_data->Stores.at("RecoEvent")->Get("GammaCount",gamma_count);
  m_data->Stores.at("RecoEvent")->Get("PositronCount",positron_count);
  if (neutron_count == 1 && (gamma_count >= 1 || positron_count >=1)) is_dsnb_like = true;
  return is_dsnb_like;

}

bool EventSelector::EventSelectionByPMTMRDCoinc() {

  if (fIsMC){
    bool has_clustered_pmt = m_data->CStore.Get("ClusterMapMC",m_all_clusters_MC);
    if (not has_clustered_pmt) { Log("EventSelector Tool: Error retrieving ClusterMapMC from CStore, did you run ClusterFinder beforehand?",v_error,verbosity); return false; }
  } else {
    bool has_clustered_pmt = m_data->CStore.Get("ClusterMap",m_all_clusters);
    if (not has_clustered_pmt) { Log("EventSelector Tool: Error retrieving ClusterMap from CStore, did you run ClusterFinder beforehand?",v_error,verbosity); return false; }
  }

  bool has_clustered_mrd = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not has_clustered_mrd) { Log("EventSelector Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
  if (MrdTimeClusters.size()!=0){
    has_clustered_mrd = m_data->CStore.Get("MrdDigitTimes",MrdDigitTimes);
    if (not has_clustered_mrd) { Log("EventSelector Tool: Error retrieving MrdDigitTimes map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
    has_clustered_mrd = m_data->CStore.Get("MrdDigitChankeys",MrdDigitChankeys);
    if (not has_clustered_mrd) { Log("EventDisplay Tool: Error retrieving MrdDigitChankeys, did you run TimeClustering beforehand",v_error,verbosity); return false;}
  }
  
  int pmt_cluster_size;
  if (fIsMC) pmt_cluster_size = (int) m_all_clusters_MC->size();
  else pmt_cluster_size = (int) m_all_clusters->size();
  m_data->Stores["RecoEvent"]->Set("NumPMTClusters",pmt_cluster_size);
  vec_pmtclusters_charge->clear();
  vec_pmtclusters_time->clear();
  m_data->Stores["RecoEvent"]->Set("PMTClustersCharge",vec_pmtclusters_charge,false);
  m_data->Stores["RecoEvent"]->Set("PMTClustersTime",vec_pmtclusters_time,false);
  vec_mrdclusters_time->clear();
  m_data->Stores["RecoEvent"]->Set("MRDClustersTime",vec_mrdclusters_time);


  bool prompt_cluster = false;
  double pmt_time = 0;


  if (fIsMC){
    if (m_all_clusters_MC->size()){
      double max_charge = 0;
      double cluster_time;
      for(std::pair<double,std::vector<MCHit>>&& apair : *m_all_clusters_MC){
        std::vector<MCHit>&MCHits = apair.second;
        double time_temp = 0;
        double charge_temp = 0;
        for (unsigned int i_hit = 0; i_hit < MCHits.size(); i_hit++){
          time_temp+=MCHits.at(i_hit).GetTime();
          charge_temp+=MCHits.at(i_hit).GetCharge();
        }
        if (MCHits.size()>0) time_temp/=MCHits.size();
        vec_pmtclusters_charge->push_back(charge_temp);
        vec_pmtclusters_time->push_back(time_temp);
        if (time_temp > 2000.) continue;	//not a prompt event
        if (charge_temp > max_charge){
          max_charge = charge_temp;
          prompt_cluster = true;
          pmt_time = time_temp;
        }
      }
    }
  } else {
    if (m_all_clusters->size()){
      double max_charge = 0;
      double cluster_time;
      for(std::pair<double,std::vector<Hit>>&& apair : *m_all_clusters){
        std::vector<Hit>&Hits = apair.second;
        double time_temp = 0;
        double charge_temp = 0;
        for (unsigned int i_hit = 0; i_hit < Hits.size(); i_hit++){
          time_temp+=Hits.at(i_hit).GetTime();
          charge_temp+=Hits.at(i_hit).GetCharge();
        }
        if (Hits.size()>0) time_temp/=Hits.size();
        vec_pmtclusters_charge->push_back(charge_temp);
        vec_pmtclusters_time->push_back(time_temp);
        if (time_temp > 2000.) continue;	//not a prompt event
        if (charge_temp > max_charge){
          max_charge = charge_temp;
          prompt_cluster = true;
          pmt_time = time_temp;
        }
      }
    }
  }

  m_data->Stores["RecoEvent"]->Set("PMTClustersCharge",vec_pmtclusters_charge,false);
  m_data->Stores["RecoEvent"]->Set("PMTClustersTime",vec_pmtclusters_time,false);

  std::vector<double> mrd_meantimes;
  for(unsigned int thiscluster=0; thiscluster<MrdTimeClusters.size(); thiscluster++){
 
    std::vector<int> hitmrd_times;
    std::vector<int> single_mrdcluster = MrdTimeClusters.at(thiscluster);
    int numdigits = single_mrdcluster.size();
    double mrd_meantime = 0.;
    for(int thisdigit=0;thisdigit<numdigits;thisdigit++){
      int digit_value = single_mrdcluster.at(thisdigit);
      unsigned long chankey = MrdDigitChankeys.at(digit_value);
      Detector *thedetector = fGeometry->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement()=="MRD") {
        double mrdtimes=MrdDigitTimes.at(digit_value);
        hitmrd_times.push_back(mrdtimes);
        mrd_meantime += mrdtimes;
      }
      if (hitmrd_times.size()>0) mrd_meantime /= hitmrd_times.size();
      mrd_meantimes.push_back(mrd_meantime);
    }
  }

  vec_mrdclusters_time->clear();
  for (int i=0; i<(int)mrd_meantimes.size(); i++){
    vec_mrdclusters_time->push_back(mrd_meantimes.at(i));
  }
  m_data->Stores["RecoEvent"]->Set("MRDClustersTime",vec_mrdclusters_time);
  
  if (fIsMC){
    if (MrdTimeClusters.size() == 0 || m_all_clusters_MC->size() == 0) return false;
  } else {
    if (MrdTimeClusters.size() == 0 || m_all_clusters->size() == 0) return false;
  }

  double pmtmrd_coinc_min = fPMTMRDOffset - 50;
  double pmtmrd_coinc_max = fPMTMRDOffset + 50;

  bool coincidence = false;
  for (int i_mrd = 0; i_mrd < int(mrd_meantimes.size()); i_mrd++){
    double time_diff = mrd_meantimes.at(i_mrd) - pmt_time;
    if (verbosity > 0) std::cout <<"MRD time: "<<mrd_meantimes.at(i_mrd)<<", PMT time: "<<pmt_time<<", difference: "<<time_diff<<std::endl;
    Log("EventSelector tool: MRD/Tank coincidene candidate "+std::to_string(i_mrd)+ " has time difference: "+std::to_string(time_diff),v_message,verbosity);
    if (time_diff > pmtmrd_coinc_min && time_diff < pmtmrd_coinc_max){
      coincidence = true;
    }
  }

  return coincidence;

}

bool EventSelector::EventSelectionByVetoCut(){

 bool has_veto = false;
 if (fIsMC) {
    if(TDCData_MC){
    if (TDCData_MC->size()==0){
      Log("EventSelector tool: TDC data is empty in this event.",v_message,verbosity);
    } else {
      for (auto&& anmrdpmt : (*TDCData_MC)){
        unsigned long chankey = anmrdpmt.first;
        Detector* thedetector = fGeometry->ChannelToDetector(chankey);
        unsigned long detkey = thedetector->GetDetectorID();
        if (thedetector->GetDetectorElement()=="Veto") has_veto = true;
      }
    }
  } else {
    Log("EventSelector tool: No TDC data available in this event.",v_message,verbosity);
  }
 }
else {
  if(TDCData){
    if (TDCData->size()==0){
      Log("EventSelector tool: TDC data is empty in this event.",v_message,verbosity);
    } else {
      for (auto&& anmrdpmt : (*TDCData)){
        unsigned long chankey = anmrdpmt.first;
        Detector* thedetector = fGeometry->ChannelToDetector(chankey);
        unsigned long detkey = thedetector->GetDetectorID();
        if (thedetector->GetDetectorElement()=="Veto") has_veto = true;
      }
    }
  } else {
    Log("EventSelector tool: No TDC data available in this event.",v_message,verbosity);
  }
  }

  return (!has_veto);	//Successful selection means no veto hit 

}

void EventSelector::Reset() {
  // Reset 
  fEventApplied = EventSelector::kFlagNone;
  fEventFlagged = EventSelector::kFlagNone;
  fEventCutStatus = true; 
} 
