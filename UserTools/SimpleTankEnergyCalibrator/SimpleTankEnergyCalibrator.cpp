#include "SimpleTankEnergyCalibrator.h"

SimpleTankEnergyCalibrator::SimpleTankEnergyCalibrator():Tool(){}


bool SimpleTankEnergyCalibrator::Initialise(std::string configfile, DataModel &data){
  std::cout << "SimpleTankEnergyCalibrator Tool: Initializing..." << std::endl;
  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  verbosity = 0;

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("TankBeamWindowStart",TankBeamWindowStart);
  m_variables.Get("TankBeamWindowEnd",TankBeamWindowEnd);
  m_variables.Get("TankNHitThreshold",TankNHitThreshold);
  m_variables.Get("MinPenetrationDepth",MinPenetrationDepth);
  m_variables.Get("MaxAngle",MaxAngle);
  m_variables.Get("MaxEntryPointRadius",MaxEntryPointRadius);
  
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

  //TODO: Read in a simple CSV file that has channel charge/PE conversions.
  
  std::cout << "SimpleTankEnergyCalibrator Tool: Loading Charge to PE Map" << std::endl;
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap",ChannelKeyToSPEMap);

  std::cout << "SimpleTankEnergyCalibrator Tool Initialized" << std::endl;
  return true;
}


bool SimpleTankEnergyCalibrator::Execute(){

  //All right, so let's do a few things:
  //  - Access the MRDTracks.  return true if the following conditions hold:
  //      Not one track
  //      Angle and minimum track length not satisfied
  //      Entry point not within the circle of consideration
  //  - If we're good here, then Loop through the hits in the time window.  Convert 
  //    each hit's charge to PE count using a ChargeToPE map.  For now, just have the
  //    total PE and estimated Muon energy in MRD printed.
  
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl;}
  
  std::map<unsigned long, std::vector<std::vector<ADCPulse>>> RecoADCHits;
  //First, get the hits from the Store
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  bool got_recoadc = m_data->Stores["ANNIEEvent"]->Get("RecoADCHits",RecoADCHits);

  bool got_hits = m_data->Stores["ANNIEEvent"]->Get("Hits", Hits);
  if (!got_hits){
    std::cout << "No Hits store in ANNIEEvent! " << std::endl;
    return false;
  }
  bool get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);  // a std::map<ChannelKey,vector<TDCHit>>
  if (!got_hits){
    std::cout << "No TDCData store in ANNIEEvent! " << std::endl;
    return false;
  }
  bool PassesCriteria = true;

  //Get beam window hits 
  std::vector<Hit> BeamHits = this->GetInWindowHits();
  if(verbosity>3) std::cout << "SimpleTankEnergyCalculator Tool: Hit count is " << BeamHits.size() << std::endl;
  if((int)BeamHits.size() < TankNHitThreshold){
    if(verbosity>3) std::cout << "SimpleTankEnergyCalculator Tool: Hit count " << BeamHits.size() << " not over threshold" << std::endl;
    PassesCriteria = false;
  }

  //Check muon track
  m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
  m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType",MRDTriggertype);
  
  // In case of a cosmic event, use the fitted track for an efficiency calculation
  
  //if (MRDTriggertype != "Beam"){
  //  if(verbosity>3) std::cout << "SimpleTankEnergyCalculator Tool: MRD trigger " << MRDTriggertype << "not Beam.  Won't be corresponding PMT data." << std::endl;
  //  PassesCriteria = false;
  //}
	
  // Get MRD track information from MRDTracks BoostStore
  m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents",numsubevs);
  m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);
  
  if(numtracksinev!=1) {
    if(verbosity>3) std::cout << "SimpleTankEnergyCalculator Tool: Either zero or more than two reconstructed tracks." << std::endl;
    PassesCriteria = false;
  }
  
  bool facc_hit = false;
  //Check for FACC activity
  if(TDCData->size()==0){
    Log("SimpleTankEnergyCalibrator tool: No TDC hits available!",v_message,verbosity);
    return true;
  } else {
      Log("SimpleTankEnergyCalibrator tool: Looping over FACC/MRD hits... looking for Veto activity",v_message,verbosity);
    for(auto&& anmrdpmt : (*TDCData)){
      unsigned long chankey = anmrdpmt.first;
      Detector* thedetector = geom->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if(thedetector->GetDetectorElement()=="Veto") facc_hit=true; // this is a veto hit, not an MRD hit.
    }
  }

  if(!facc_hit){
    Log("SimpleTankEnergyCalibrator tool: No Veto hits, so not a through-going muon candidate",v_message,verbosity);
    PassesCriteria = false;
  }


  //Check for valid track criteria
  m_data->Stores["MRDTracks"]->Get("MRDTracks",theMrdTracks);
  // Loop over reconstructed tracks
  
  paddlesInTrackReco.clear();
 
  double tracklength = 0;
  TrackAngle = -9999;
  PenetrationDepth = -9999;
  EnergyLoss = -9999;
  MrdEntryPoint.SetX(-9999.);
  MrdEntryPoint.SetY(-9999.);
  MrdEntryPoint.SetZ(-9999.);

  for(int tracki=0; tracki<numtracksinev; tracki++){
    BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
    
    //get track properties that are needed for the through-going muon selection
    thisTrackAsBoostStore->Get("StartVertex",StartVertex);
    thisTrackAsBoostStore->Get("StopVertex",StopVertex);
    thisTrackAsBoostStore->Get("TrackAngle",TrackAngle);
    thisTrackAsBoostStore->Get("TrackAngleError",TrackAngleError);
    thisTrackAsBoostStore->Get("PenetrationDepth",PenetrationDepth);
    thisTrackAsBoostStore->Get("MrdEntryPoint",MrdEntryPoint);
    thisTrackAsBoostStore->Get("NumLayersHit",LayersHit);
    thisTrackAsBoostStore->Get("EnergyLoss",EnergyLoss);
    thisTrackAsBoostStore->Get("EnergyLossError",EnergyLossError);
    tracklength = sqrt(pow((StopVertex.X()-StartVertex.X()),2)+pow(StopVertex.Y()-StartVertex.Y(),2)+pow(StopVertex.Z()-StartVertex.Z(),2));
  }

  double EntryPointRadius = sqrt(pow(MrdEntryPoint.X(),2) + pow(MrdEntryPoint.Y(),2)) * 100.0; // convert to cm
  PenetrationDepth = PenetrationDepth*100.0;

  if(verbosity>3){
    std::cout << "SimpleTankEnergyCalibrator tool: TrackAngle,EntryPointRadius,PenetrationDepth: " << 
        TrackAngle << "," << EntryPointRadius << "," << PenetrationDepth << std::endl;
  }

  if(TrackAngle>MaxAngle || EntryPointRadius>MaxEntryPointRadius || PenetrationDepth<MinPenetrationDepth){
    if(verbosity>3){
      std::cout << "SimpleTankEnergyCalibrator tool: Reco Track criteria doesn't satisfy requirements." << std::endl;
    }
    PassesCriteria = false;
  }

  //We're in business.  We've got a Minimum-ionizing candidate.  Get the total number of PE and
  //print it out.
  //
  if(PassesCriteria){
    double TotalPE = this->GetTotalPE(BeamHits);
    double TotalQ = this->GetTotalQ(BeamHits);
    std::cout << "SimpleTankEnergyCalibrator tool: THROUGH-GOING MUON CANDIDATE FOUND." << std::endl;
    std::cout << "SimpleTankEnergyCalibrator tool: ENTRYX,ENTRYY,ENTRYZ," << MrdEntryPoint.X() << "," <<
        MrdEntryPoint.Y() << "," << MrdEntryPoint.Z() << std::endl;
    std::cout << "SimpleTankEnergyCalibrator tool: TRACKANGLE," << TrackAngle << std::endl;
    std::cout << "SimpleTankEnergyCalibrator tool: PENETRATION," << PenetrationDepth << std::endl;
    std::cout << "SimpleTankEnergyCalibrator tool: CHARGE,PE," << TotalQ << "," << TotalPE << std::endl;
  }

  return true;
}


bool SimpleTankEnergyCalibrator::Finalise(){

  std::cout << "SimpleTankEnergyCalibrator exitting" << std::endl;
  return true;
}

std::vector<Hit> SimpleTankEnergyCalibrator::GetInWindowHits(){
  std::vector<Hit> BeamHits;

  for(std::pair<unsigned long, std::vector<Hit>>&& apair : *Hits){
    unsigned long chankey = apair.first;
    std::vector<Hit> ThisTubeHits = apair.second;
    for(int i=0; i<(int)ThisTubeHits.size(); i++){
      Hit ahit = ThisTubeHits.at(i);
      double hittime = ahit.GetTime();
      if((hittime > TankBeamWindowStart) && (hittime < TankBeamWindowEnd)) BeamHits.push_back(ahit);
    }
  }
  return BeamHits;
}

double SimpleTankEnergyCalibrator::GetTotalQ(std::vector<Hit> AllHits){
  double TotalCharge = 0;
  for (int i = 0; i<(int)AllHits.size(); i++){
    Hit ahit = AllHits.at(i);
    double ThisHitCharge = ahit.GetCharge();
    TotalCharge+=ThisHitCharge;
  }
  return TotalCharge;
}

double SimpleTankEnergyCalibrator::GetTotalPE(std::vector<Hit> AllHits){
  double TotalPE = 0;
  for (int i = 0; i<(int)AllHits.size(); i++){
    Hit ahit = AllHits.at(i);
    double ThisHitCharge = ahit.GetCharge();
    int ThisHitID = ahit.GetTubeId();
    //Search SPEMap for channel ID; if not there, ignore it.
    //Check if this CardData is next in it's sequence for processing
    std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(ThisHitID);
    if(it != ChannelKeyToSPEMap.end()){ //Charge to SPE conversion is available
      double ThisHitPE = ThisHitCharge / ChannelKeyToSPEMap.at(ThisHitID);
      TotalPE+=ThisHitPE;   
    } else {
      if(verbosity>2){
        std::cout << "FOUND A HIT FOR CHANNELKEY " << ThisHitID << "BUT NO CONVERSION " <<
            "TO PE AVAILABLE.  SKIPPING." << std::endl;
      }
    }
  }
  return TotalPE;
}
