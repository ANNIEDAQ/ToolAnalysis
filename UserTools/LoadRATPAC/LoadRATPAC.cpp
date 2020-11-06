#include "LoadRATPAC.h"

LoadRATPAC::LoadRATPAC():Tool(){}


bool LoadRATPAC::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
	
  if(verbosity) std::cout<<"Initializing Tool LoadRATPAC"<<std::endl;
m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  m_variables.Get("InputFile",filename_ratpac);
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("LappdNumStrips", LappdNumStrips);
  m_variables.Get("LappdStripLength", LappdStripLength);           // [mm]
  m_variables.Get("LappdStripSeparation", LappdStripSeparation);   // [mm]
  m_variables.Get("xtankcenter",xtankcenter);
  m_variables.Get("ytankcenter",ytankcenter);
  m_variables.Get("ztankcenter",ztankcenter);
  
  // Initialize the TChain reading event info from our file
  //chain = new TChain("T");
  //chain->Add(filename_ratpac);
  //Initialize the DSReader
  dsReader = new RAT::DSReader(filename_ratpac.c_str()); // cuz stupid dsReader needs to be initialize with a filename..
  
  //Initialize run info
  std::cout<<"Initializing the run information"<<std::endl;
  
  //tri = new TChain("T");
  runtri = new TChain("runT");

  TFile *ftemp = TFile::Open(filename_ratpac.c_str());
  if (ftemp->Get("runT")) {
    runtri->Add(filename_ratpac.c_str());
    RAT::DS::RunStore::SetReadTree(runtri);
  } else {
      std::cout<<"File does not have a run tree!  Not the right format file"<<std::endl;
     return false; // else, no runT, so don't register runtri with RunStore
  }
  delete ftemp;
  
  ds = dsReader->GetEvent(0);
  run = new RAT::DS::Run;
  run = RAT::DS::RunStore::Get()->GetRun(ds);
  RunNumber = static_cast<uint32_t>(run->GetID());
  SubrunNumber = 0; 
  pmtInfo = new RAT::DS::PMTInfo;
  pmtInfo = run->GetPMTInfo();
  lappdInfo = new RAT::DS::LAPPDInfo;
  lappdInfo = run->GetLAPPDInfo();
  
  //Total number of events available to be read
  NbEntries = dsReader->GetTotal();
	
  // things to be saved to the ANNIEEvent Store
  MCParticles = new std::vector<MCParticle>;
  MCHits = new std::map<unsigned long,std::vector<MCHit>>;
  MCLAPPDHits = new std::map<unsigned long,std::vector<MCLAPPDHit>>;
  EventTime = new TimeClass();
  
  this->LoadANNIEGeometry();

  //FIXME: No triggers in RAT-PAC yet, but should be fine without for now
  //TriggerData = new std::vector<TriggerClass>{beamtrigger}; 
  
  //Counter corresponding to Store number
  EventNumber=0;

  //Counter corresponding to event number in MC
  //Without tracking any detector triggers, is the same as EventNum
  MCEventNum=0;

  //Currently, RAT-PAC has no trigger really; just per event.  Keep this zero
  MCTriggernum=0;

  //FIXME: Fill in the run start time as GetUTC from first event?  Or is in runT?
  time_t runstarttime = run->GetStartTime();
  uint64_t runstart = static_cast<uint64_t>(runstarttime);
  TimeClass RunStartTime(runstart);
	
  // use nominal beam values for now.  Eventually, will be loaded in data
  double beaminten=4.777e+12;
  double beampow=3.2545e+16;
  BeamStatus = new BeamStatusClass(RunStartTime, beaminten, beampow, "stable");
  
  if(verbosity) std::cout<<"Initialization of Tool LoadRATPAC Complete"<<std::endl;
  return true;
}


bool LoadRATPAC::Execute(){
  if(verbosity) std::cout<<"Executing Tool LoadRATPAC"<<std::endl;
  // First, reset all counters kept track of per event
  this->Reset();
  

  if (EventNumber >= NbEntries){
    std::cout << "LoadRATPAC tool: end of file reached.  Returning true" << std::endl;
    return true;
  }
  
  if(verbosity) std::cout<<"LoadRATPAC tool: Loading event in DSReader"<<std::endl;
  //chain->GetEntry(EventNumber);
  //Ds takesin a Ulong64_t... this is cludgy but should work
  ULong64_t evnum = static_cast<ULong64_t>(EventNumber);
  ds = dsReader->GetEvent(EventNumber);

  if(verbosity) std::cout<<"LoadRATPAC tool: Loading event time in TimeClass"<<std::endl;
  //FIXME: WCSim gives in trigger data, but we don't have triggers in RAT-PAC yet
  EventTimeNs = 0; //FIXME:  Should be called with chain->GetUTC() when we do have evnet times?
 EventTime->SetNs(EventTimeNs);

  if(verbosity) std::cout<<"LoadRATPAC tool: Begin the PMT loop"<<std::endl;
  //PMT loop
  for( int iPMT = 0; iPMT < ds->GetMC()->GetMCPMTCount(); iPMT++ ){
    if(verbosity>2) std::cout<<"getting PMT "<<iPMT<<std::endl;
    RAT::DS::MCPMT *aPMT = ds->GetMC()->GetMCPMT(iPMT);
    int tubeid = aPMT->GetID();
	  if(verbosity>2) cout<<"tubeid="<<tubeid<<endl;
    if(pmt_tubeid_to_channelkey.count(tubeid)==0){
      cerr<<"LoadRATPAC ERROR: tank PMT with no associated ChannelKey!"<<endl;
      return false;
    }
    unsigned long key = pmt_tubeid_to_channelkey.at(tubeid);
    std::vector<int> hitparents = {}; //FIXME: Get parent particles of hits
    //Loop over photons that hit the PMT for digits
    for(long iPhot = 0; iPhot < aPMT->GetMCPhotonCount(); iPhot++){
      float hitcharge = aPMT->GetMCPhoton(iPhot)->GetCharge();
      float hittime = aPMT->GetMCPhoton(iPhot)->GetHitTime(); //Time relative to event start (ns)
      MCHit thishit(key, hittime, hitcharge, hitparents);
      if(MCHits->count(key)==0) MCHits->emplace(key, std::vector<MCHit>{thishit});
      else MCHits->at(key).push_back(thishit);
      if(verbosity>2) cout<<"digit added"<<endl;
    }
  }
  std::cout << "Done loading event's PMT hits" << std::endl;

  //LAPPD loop
  if(verbosity) std::cout<<"LoadRATPAC tool: Begin the LAPPD loop"<<std::endl;
  for( int iLAPPD = 0; iLAPPD < ds->GetMC()->GetMCLAPPDCount(); iLAPPD++ ){
    if(verbosity>2) std::cout<<"getting LAPPD "<<iLAPPD<<std::endl;
    RAT::DS::MCLAPPD *aLAPPD = ds->GetMC()->GetMCLAPPD(iLAPPD);
    int lappdid = aLAPPD->GetID();
	  if(verbosity>2) cout<<"LAPPDid="<<lappdid<<endl;
    if(lappd_tubeid_to_detectorkey.count(lappdid)==0){
      cerr<<"LoadRATPAC ERROR: LAPPD with no associated ChannelKey!"<<endl;
      return false;
    }
    unsigned int detkey = lappd_tubeid_to_detectorkey.at(lappdid);
    Detector* thedet = anniegeom->GetDetector(detkey);
    unsigned int key = thedet->GetChannels()->begin()->first; // first strip on this LAPPD
    if(verbosity>2) cout<<"Channelkey for LAPPD="<<key<<endl;
    TVector3 lappddir = lappdInfo->GetDirection(lappdid);
    std::cout<<"LAPPD DIRECTION IN RATPAC COORD:: " << lappddir.X() << "," <<
        lappddir.Y() << "," << lappddir.Z() << std::endl;
    TVector3 sensor_direction(lappddir.X(),lappddir.Y(),lappddir.Z());
    //Loop over photons that hit the LAPPD for digits
    for(long iPhot = 0; iPhot < aLAPPD->GetMCPhotonCount(); iPhot++){
      std::vector<int> lappdhitparents = {}; //FIXME: Get parent particles of lappdhits
      double hitcharge = aLAPPD->GetMCPhoton(iPhot)->GetCharge();
      double hittime = aLAPPD->GetMCPhoton(iPhot)->GetHitTime(); //Time relative to event start (ns)
      //Get LAPPD Position in RATPAC coordinates
      double lappdx = lappdInfo->GetPosition(lappdid).X();
      double lappdy = lappdInfo->GetPosition(lappdid).Y();
      double lappdz = lappdInfo->GetPosition(lappdid).Z();
      // Get Local hit position on LAPPD relative to sensor center position
      double localx = aLAPPD->GetMCPhoton(iPhot)->GetPosition().X();
      double localy = aLAPPD->GetMCPhoton(iPhot)->GetPosition().Y();
      double localz = aLAPPD->GetMCPhoton(iPhot)->GetPosition().Z();
      // Rotate local hit position into RATPAC coordinates frame
      TVector3 hit_position(localx,localy,localz);
      TVector3 sensor_position(lappdx,lappdy,lappdz);
      hit_position.RotateUz(sensor_direction);
      //Load hit position into MCLAPPDHit in z-beam coordinates.
      std::vector<double> globalPosition{(hit_position.X()+sensor_position.X())/1000.,
              (hit_position.Z()+sensor_position.Z())/1000.,-(hit_position.Y()+sensor_position.Y())/1000.};
      //std::vector<double> globalPosition{lappdx+trans.X(),lappdy+trans.Y(),lappdz+trans.Z()};
      std::vector<double> localPosition{localx,localy};
	    logmessage = "  LAPPDPosition = ("+to_string(lappdx) + ", " + to_string(lappdy) + ", " + to_string(lappdz) + ") "+ "\n";
      logmessage = "  LAPPDHitGlobalPosition = ("+to_string(globalPosition[0]) + ", " + to_string(globalPosition[1]) + ", " + to_string(globalPosition[2]) + ") "+ "\n";
      //logmessage = "  RotatedLocalPosition = ("+to_string(trans.X()) + ", " + to_string(trans.Y()) + ", " + to_string(trans.Z()) + ") "+ "\n";
      //logmessage += "  LAPPDHitLocalPosition(x,y,z in plane coord) = ("+to_string(localPosition[0]) + ", " + to_string(localPosition[1]) + "," + to_string(localPosition[2]) + ") "+ "\n";
	    Log(logmessage,v_debug,verbosity);
      MCLAPPDHit thishit(key, hittime, hitcharge,globalPosition, localPosition,lappdhitparents);
      if(MCLAPPDHits->count(key)==0){
          MCLAPPDHits->emplace(key, std::vector<MCLAPPDHit>{thishit});
      } else {
        MCLAPPDHits->at(key).push_back(thishit);
      }
      if(verbosity>2) cout<<"LAPPD hit added"<<endl;
    }
  }
  std::cout << "Done loading event's LAPPD hits" << std::endl;

  if(verbosity>1) cout<<"getting "<<ds->GetMC()->GetMCTrackCount()<<" tracks"<<endl;
  //Loop through all the tracks and get particle information.  We'll load it
  //Load MC Parent Particles
	tracktype startstoptype = tracktype::UNDEFINED;
  logmessage = "  Particle count = " +to_string(ds->GetMC()->GetMCTrackCount()) + " \n";
	Log(logmessage,v_debug,verbosity);
  for (long iTrack=0; iTrack<ds->GetMC()->GetMCTrackCount(); iTrack++){
    RAT::DS::MCTrack  *thistrack = ds->GetMC()->GetMCTrack(iTrack);
    int parentid = thistrack->GetParentID();
    int particleid = thistrack->GetID();
    int pdgcode = thistrack->GetPDGCode();
    std::string particlename = thistrack->GetParticleName();
    logmessage = "  Particle name = " +particlename + " \n";
	  Log(logmessage,v_debug,verbosity);
    if(particlename=="opticalphoton") continue; //Do not save optical photons
    if(particlename=="gamma") continue; //Do not save gammas
    Float_t tracklength = thistrack->GetLength();
    RAT::DS::MCTrackStep* firstStep = thistrack->GetMCTrackStep(0);
    RAT::DS::MCTrackStep* lastStep = thistrack->GetMCTrackStep(thistrack->GetMCTrackStepCount()-1); 
    TVector3 endpoint = lastStep->GetEndpoint();
    double EndKE = lastStep->GetKE();
    double endtime = lastStep->GetGlobalTime();
    TVector3 startpoint = firstStep->GetEndpoint();
    double starttime = firstStep->GetGlobalTime();
    logmessage = "  Particle start position = (" +to_string(startpoint.X()) + ", " + to_string(startpoint.Y()) + ", " + to_string(startpoint.Z()) +" \n";
    logmessage += "  Particle end position = (" +to_string(endpoint.X()) + ", " + to_string(endpoint.Y()) + ", " + to_string(endpoint.Z()) +" \n";
    logmessage += "  Particle start time = (" +to_string(starttime) + " \n";
	  Log(logmessage,v_debug,verbosity);
    double StartKE = firstStep->GetKE();
    TVector3 particledir = firstStep->GetMomentum();
    particledir = particledir.Unit();
    logmessage = "  Particle mom. direction = (" +to_string(particledir.X()) + ", " + to_string(particledir.Y()) + ", " + to_string(particledir.Z()) +" \n";
	  Log(logmessage,v_debug,verbosity);
    Direction partdir;
    partdir.SetX(particledir.X());
    partdir.SetY(particledir.Y());
    partdir.SetZ(particledir.Z());
    //Load information into MCParticle 
    //TODO: RATPAC particles are already in z-beam axis coordinates!  We should make particles & detectors match coordinates in RATPAC..
    MCParticle thisparticle(
            pdgcode, StartKE, EndKE, Position((startpoint.X())/1000., 
            (startpoint.Y())/1000., (startpoint.Z())/1000.), Position((endpoint.X())/1000.,
            (endpoint.Y())/1000., (endpoint.Z())/1000.), starttime,
            endtime, partdir,
            tracklength/1000., startstoptype, particleid, parentid, -9999,0);
    MCParticles->push_back(thisparticle);
  }


	// set event level variables
	if(verbosity>1) std::cout<<"setting the store variables"<<std::endl;
	m_data->Stores.at("ANNIEEvent")->Set("RunNumber",RunNumber);
	m_data->Stores.at("ANNIEEvent")->Set("SubrunNumber",SubrunNumber);
	m_data->Stores.at("ANNIEEvent")->Set("EventNumber",EventNumber);
	if(verbosity>2) cout<<"particles"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("MCParticles",MCParticles,true);
	if(verbosity>2) cout<<"hits"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("MCHits",MCHits,true);
	if(verbosity>2) cout<<"LAPPDhits"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("MCLAPPDHits",MCLAPPDHits,true);
	//if(verbosity>2) cout<<"tdcdata"<<endl;
	//m_data->Stores.at("ANNIEEvent")->Set("TDCData",TDCData,true);
	//if(verbosity>2) cout<<"triggerdata"<<endl;
	//m_data->Stores.at("ANNIEEvent")->Set("TriggerData",TriggerData,true);  // FIXME
	if(verbosity>2) cout<<"eventtime"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("EventTime",EventTime,true);
	m_data->Stores.at("ANNIEEvent")->Set("MCEventNum",MCEventNum);
	m_data->Stores.at("ANNIEEvent")->Set("MCTriggernum",MCTriggernum);
	m_data->Stores.at("ANNIEEvent")->Set("MCFile",filename_ratpac);
	m_data->Stores.at("ANNIEEvent")->Set("MCFlag",true);           // constant
	m_data->Stores.at("ANNIEEvent")->Set("BeamStatus",BeamStatus,true);

  EventNumber+=1;
  MCEventNum+=1;
  return true;
}


bool LoadRATPAC::Finalise(){
  Log("LoadRATPACTool: Exitting",v_message,verbosity);
  //delete chain;
  delete dsReader;
  delete run;
  return true;
}

void LoadRATPAC::Reset(){
  // Initialization time
	MCParticles->clear();
	MCHits->clear();
  MCLAPPDHits->clear();
}

void LoadRATPAC::LoadANNIEGeometry(){	
  // Make the ANNIEEvent Store if it doesn't exist
  // =============================================
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(annieeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);

  std::cout<<"LoadRATPAC tool: Loading the ANNIE geometry now"<<std::endl;  
  // construct the Geometry to go in the header Based on the RATPAC Geometry
  // =================================================================
  double RATPACGeometryVersion = 1;                       // TODO pull this from some suitable variable
  
  
  int numtankpmts = pmtInfo->GetPMTCount();
  int numlappds = lappdInfo->GetLAPPDCount();
  int nummrdpmts = 0;
  int numvetopmts = 0;
  
  //TODO: Have these values saved to the RATDS? Hard-coded for now
  //
  double tank_xcenter = (xtankcenter) / 1000.;  // convert [mm] to [m]
  double tank_ycenter = (ytankcenter) / 1000.;
  double tank_zcenter = (ztankcenter) / 1000.;
  Position tank_center(tank_xcenter, tank_ycenter, tank_zcenter);
  double tank_radius = (1524.0) / 1000.;
  double tank_halfheight = (1981.2) / 1000.;
  //Currently hard-coded; estimated with a tape measure on the ANNIE frame :)
  double pmt_enclosed_radius = 1.0;
  double pmt_enclosed_halfheight = 1.45;
  //geometry variables not yet in RATPAC. grabbed from MRDSpecs.hh
  double mrd_width =  (MRDSpecs::MRD_width) / 100.; // convert [cm] to [m]
  double mrd_height = (MRDSpecs::MRD_height) / 100.;
  double mrd_depth =  (MRDSpecs::MRD_depth) / 100.;
  double mrd_start =  (MRDSpecs::MRD_start) / 100.;
  if(verbosity>1) cout<<"we have "<<numtankpmts<<" tank pmts, and "<<numlappds<<" lappds"<<endl;
  
  // construct the goemetry
  anniegeom = new Geometry(RATPACGeometryVersion,
                           tank_center,
                           tank_radius,
                           tank_halfheight,
                           pmt_enclosed_radius,
                           pmt_enclosed_halfheight,
                           mrd_width,
                           mrd_height,
                           mrd_depth,
                           mrd_start,
                           numtankpmts,
                           nummrdpmts,
                           numvetopmts,
                           numlappds,
                           geostatus::FULLY_OPERATIONAL);
  
  if(verbosity>1){ cout<<"constructed anniegom at "<<anniegeom<<" with tank origin "; tank_center.Print(); }
  m_data->Stores.at("ANNIEEvent")->Header->Set("AnnieGeometry",anniegeom,true);
  
  // Construct the Detectors and Channels
  // Taken from LoadWCSim by Marcus
  // ====================================
  // PMTs
  unsigned int ADC_Crate_Num = 0;
  unsigned int ADC_Card_Num  = 0;
  unsigned int ADC_Chan_Num = 0;
  unsigned int MT_Crate_Num = 0;
  unsigned int MT_Card_Num = 0;
  unsigned int MT_Chan_Num = 0;
  // LAPPDs
  unsigned int ACDC_Crate_Num = 0;
  unsigned int ACDC_Card_Num = 0;
  unsigned int ACDC_Chan_Num = 0;
  unsigned int ACC_Crate_Num = 0;
  unsigned int ACC_Card_Num = 0;
  unsigned int ACC_Chan_Num = 0;
  // TDCs (not in RAT-PAC yet)
  //unsigned int TDC_Crate_Num = 0;
  //unsigned int TDC_Card_Num = 0;
  //unsigned int TDC_Chan_Num = 0;
  // HV
  unsigned int CAEN_HV_Crate_Num = 0;
  unsigned int CAEN_HV_Card_Num = 0;
  unsigned int CAEN_HV_Chan_Num = 0;
  // MRD not in RAT-PAC yet , so no LeCroy HV
  //unsigned int LeCroy_HV_Crate_Num = 0;
  //unsigned int LeCroy_HV_Card_Num = 0;
  //unsigned int LeCroy_HV_Chan_Num = 0;
  unsigned int LAPPD_HV_Crate_Num = 0;
  unsigned int LAPPD_HV_Card_Num = 0;
  unsigned int LAPPD_HV_Chan_Num = 0;
  
  
  // lappds
  for(int lappdi=0; lappdi<numlappds; lappdi++){
    TVector3 lappdpos = lappdInfo->GetPosition(lappdi);
    TVector3 lappddir = lappdInfo->GetDirection(lappdi);
    int modelType = lappdInfo->GetModel(lappdi);
    logmessage = "  LAPPD position,RATPAC coord = (" +to_string(lappdpos.X()) + ", " + to_string(lappdpos.Y()) + ", " + to_string(lappdpos.Z()) +" \n";
    logmessage += "  LAPPD direction, RP Coord = (" +to_string(lappddir.X()) + ", " + to_string(lappddir.Y()) + ", " + to_string(lappddir.Z()) +" \n";
    Log(logmessage,v_debug,verbosity);
    // Construct the detector associated with this tile; LAPPD coordinates
    // different than RATPAC global coordinates
    unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
    lappd_tubeid_to_detectorkey.emplace(lappdi,uniquedetectorkey);
    detectorkey_to_lappdid.emplace(uniquedetectorkey,lappdi);
    //Unit conversion for position applied post-hit rotations
    Detector adet(uniquedetectorkey,
                  "LAPPD",
                  "LAPPD",
                  Position( lappdpos.X()/1000.,
                            lappdpos.Z()/1000.,
                            -lappdpos.Y()/1000.),
                  Direction(lappddir.X(),
                            lappddir.Z(),
                            -lappddir.Y()),
                  lappdInfo->GetModelName(modelType),
                  detectorstatus::ON,
                  0.);
    
    // construct all the channels associated with this LAPPD
    for(int stripi=0; stripi<LappdNumStrips; stripi++){
      unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
      lappd_tubeid_to_channelkey.emplace(stripi, uniquechannelkey);
      channelkey_to_lappdid.emplace(uniquechannelkey,stripi);
      int stripside = ((stripi%2)==0);   // StripSide=0 for LHS (x<0), StripSide=1 for RHS (x>0)
      int stripnum = (int)(stripi/2);    // Strip number: add 2 channels per strip as we go
      double xpos = (stripside) ? -LappdStripLength : LappdStripLength;
      double ypos = (stripnum*LappdStripSeparation) - ((LappdNumStrips*LappdStripSeparation)/2.);
      
      // fill up ADC cards and channels monotonically, they're arbitrary for simulation
      ACDC_Chan_Num++;
      if(ACDC_Chan_Num>=RC::ACDC_CHANNELS_PER_CARD)  { ACDC_Chan_Num=0; ACDC_Card_Num++; ACC_Chan_Num++; }
      if(ACDC_Card_Num>=RC::ACDC_CARDS_PER_CRATE)    { ACDC_Card_Num=0; ACDC_Crate_Num++; }
      if(ACC_Chan_Num>=RC::ACC_CHANNELS_PER_CARD)    { ACC_Chan_Num=0; ACC_Card_Num++;    }
      if(ACC_Card_Num>=RC::ACC_CARDS_PER_CRATE)      { ACC_Card_Num=0; ACC_Crate_Num++;   }
      // same for HV
      LAPPD_HV_Chan_Num++;
      if(LAPPD_HV_Chan_Num>=RC::LAPPD_HV_CHANNELS_PER_CARD)    { LAPPD_HV_Chan_Num=0; LAPPD_HV_Card_Num++;  }
      if(LAPPD_HV_Card_Num>=RC::LAPPD_HV_CARDS_PER_CRATE) { LAPPD_HV_Card_Num=0; LAPPD_HV_Crate_Num++; }
      
      Channel lappdchannel(uniquechannelkey,
                           Position(xpos,ypos,0.),
                           stripside,
                           stripnum,
                           ACDC_Crate_Num,
                           ACDC_Card_Num,
                           ACDC_Chan_Num,
                           ACC_Crate_Num,
                           ACC_Card_Num,
                           ACC_Chan_Num,
                           LAPPD_HV_Crate_Num,
                           LAPPD_HV_Card_Num,
                           LAPPD_HV_Chan_Num,
                           channelstatus::ON);
      
      // Add this channel to the geometry
      if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
      adet.AddChannel(lappdchannel);
    }
    if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
    // Add this detector to the geometry
    anniegeom->AddDetector(adet);
    if(verbosity>4) cout<<"printing geometry"<<endl;
    if(verbosity>4) anniegeom->PrintChannels();
  }
  
  // tank PMTs
  for(int pmti=0; pmti<numtankpmts; pmti++){
    TVector3 pmtpos = pmtInfo->GetPosition(pmti);
    TVector3 pmtdir = pmtInfo->GetDirection(pmti);
    int modelType = pmtInfo->GetModel(pmti);
    logmessage = "  PMT position,RATPAC Coord = (" +to_string(pmtpos.X()) + ", " + to_string(pmtpos.Y()) + ", " + to_string(pmtpos.Z()) +" \n";
    Log(logmessage,v_debug,verbosity);
    
    // Construct the detector associated with this PMT; PMT coordinates 
    // Different than RATPAC global coordinates
    unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
    Detector adet(uniquedetectorkey,
                  "Tank",
                  "Tank",
                  Position( pmtpos.X()/1000.,
                            pmtpos.Z()/1000.,
                            -pmtpos.Y()/1000.),
                  Direction(pmtdir.X(),
                            pmtdir.Z(),
                            -pmtdir.Y()),
                  pmtInfo->GetModelName(modelType),
                  detectorstatus::ON,
                  0.);
    
    // construct the channel associated with this PMT
    unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
    pmt_tubeid_to_channelkey.emplace(pmti, uniquechannelkey);
    channelkey_to_pmtid.emplace(uniquechannelkey,pmti);
    
    // fill up ADC cards and channels monotonically, they're arbitrary for simulation
    ADC_Chan_Num++;
    if(ADC_Chan_Num>=RC::ADC_CHANNELS_PER_CARD)  { ADC_Chan_Num=0; ADC_Card_Num++; MT_Chan_Num++; }
    if(ADC_Card_Num>=RC::ADC_CARDS_PER_CRATE)    { ADC_Card_Num=0; ADC_Crate_Num++; }
    if(MT_Chan_Num>=RC::MT_CHANNELS_PER_CARD)    { MT_Chan_Num=0; MT_Card_Num++; }
    if(MT_Card_Num>=RC::MT_CARDS_PER_CRATE)      { MT_Card_Num=0; MT_Crate_Num++; }
    // same for HV
    CAEN_HV_Chan_Num++;
    if(CAEN_HV_Chan_Num>=RC::CAEN_HV_CHANNELS_PER_CARD)    { CAEN_HV_Chan_Num=0; CAEN_HV_Card_Num++;  }
    if(CAEN_HV_Card_Num>=RC::CAEN_HV_CARDS_PER_CRATE) { CAEN_HV_Card_Num=0; CAEN_HV_Crate_Num++; }
    
    Channel pmtchannel( uniquechannelkey,
                        Position(0,0,0.),
                        0, // stripside
                        0, // stripnum
                        ADC_Crate_Num,
                        ADC_Card_Num,
                        ADC_Chan_Num,
                        MT_Crate_Num,
                        MT_Card_Num,
                        MT_Chan_Num,
                        CAEN_HV_Crate_Num,
                        CAEN_HV_Card_Num,
                        CAEN_HV_Chan_Num,
                        channelstatus::ON);
    
    // Add this channel to the geometry
    if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
    adet.AddChannel(pmtchannel);
    
    // Add this detector to the geometry
    if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
    anniegeom->AddDetector(adet);
    if(verbosity>4) cout<<"printing geometry"<<endl;
    if(verbosity>4) anniegeom->PrintChannels();
  }
	// for other RATPAC tools that may need the RATPAC IDs
	m_data->CStore.Set("lappd_tubeid_to_detectorkey",lappd_tubeid_to_detectorkey);
	m_data->CStore.Set("lappd_tubeid_to_channelkey",lappd_tubeid_to_channelkey);
  m_data->CStore.Set("pmt_tubeid_to_channelkey",pmt_tubeid_to_channelkey);
	m_data->CStore.Set("mrd_tubeid_to_channelkey",mrd_tubeid_to_channelkey);
	m_data->CStore.Set("facc_tubeid_to_channelkey",facc_tubeid_to_channelkey);
	// inverse
	m_data->CStore.Set("detectorkey_to_lappdid",detectorkey_to_lappdid);
	m_data->CStore.Set("channelkey_to_lappdid",channelkey_to_lappdid);
	m_data->CStore.Set("channelkey_to_pmtid",channelkey_to_pmtid);
	m_data->CStore.Set("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	m_data->CStore.Set("channelkey_to_faccpmtid",channelkey_to_faccpmtid);
	
}


TMatrixD LoadRATPAC::Rotateatob(TVector3 a, TVector3 b){
  //Get rotation matrix that rotates the z-axis of one coordinate system
  // to vector d.
  //d must be normalized.
  TVector3 au = a.Unit();
  TVector3 bu = b.Unit();
  TVector3 c = au.Cross(bu);
  TVector3 cu = c.Unit();
  double cosa = au.Dot(bu);
  double s = sqrt(1- (cosa*cosa));
  double C = 1-cosa;
  TMatrixD rotn(3,3);
  rotn[0][0]= (cu.X()*cu.X()*C + cosa);
  rotn[0][1]= (cu.X()*cu.Y()*C)-(cu.Z()*s);
  rotn[0][2]= (cu.X()*cu.Z()*C)+(cu.Y()*s);
  rotn[1][0]= (cu.Y()*cu.X()*C)+(cu.Z()*s);
  rotn[1][1]= (cu.Y()*cu.Y()*C)+cosa;
  rotn[1][2]= (cu.Y()*cu.Z()*C)-(cu.X()*s);
  rotn[2][0]= (cu.Z()*cu.X()*C)-(cu.Y()*s);
  rotn[2][1]= (cu.Z()*cu.Y()*C)+(cu.X()*s);
  rotn[2][2]= (cu.Z()*cu.Z()*C)+cosa;
  return rotn;
} 
