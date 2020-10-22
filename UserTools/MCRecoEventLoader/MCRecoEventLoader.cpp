#include "MCRecoEventLoader.h"

MCRecoEventLoader::MCRecoEventLoader():Tool(){}


bool MCRecoEventLoader::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(verbosity) cout<<"Initializing Tool MCRecoEventLoader"<<endl;
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  
  ///////////////////// Defaults for Config ///////////////
  fGetPiKInfo = 1;
  fGetNRings = 1;
  fParticleID = 13;
  fDoParticleSelection = 1;
  xshift = 0.;
  yshift = 14.46469;
  zshift = -168.1;

  /// Get the Tool configuration variables
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("GetPionKaonInfo", fGetPiKInfo);
  m_variables.Get("GetNRings",fGetNRings);
  m_variables.Get("DoParticleSelection",fDoParticleSelection);
  m_variables.Get("ParticleID", fParticleID);
  m_variables.Get("xshift", xshift);
  m_variables.Get("yshift", yshift);
  m_variables.Get("zshift", zshift);

  /// Construct the other objects we'll be setting at event level,
  fMuonStartVertex = new RecoVertex();
  fMuonStopVertex = new RecoVertex();

  // Make the RecoDigit Store if it doesn't exist
  int recoeventexists = m_data->Stores.count("RecoEvent");
  if(recoeventexists==0) m_data->Stores["RecoEvent"] = new BoostStore(false,2);

  // Get particle masses map from CStore (populated by MCParticleProperties tool)
  m_data->CStore.Get("PdgMassMap",pdgcodetomass);

  std::map<int,double>::iterator it_map;
  for (it_map = pdgcodetomass.begin(); it_map != pdgcodetomass.end(); it_map++){
    int pdgcode = it_map->first;
    int pdgmass = it_map->second;
    double cherenkov_thr = GetCherenkovThresholdE(pdgcode);
    pdgcodetocherenkov.emplace(pdgcode,cherenkov_thr);
  }
  
  //std::cout <<"PdgCherenkovMap size (MCRecoEventLoader): "<<pdgcodetocherenkov.size()<<std::endl;
  // Set particle pdg - Cherenkov threshold map to CStore
  m_data->CStore.Set("PdgCherenkovMap",pdgcodetocherenkov);

  return true;
}
bool MCRecoEventLoader::Execute(){

  /// Reset everything
  this->Reset();

  // see if "ANNIEEvent" exists
  auto get_annieevent = m_data->Stores.count("ANNIEEvent");
  if(!get_annieevent){
    Log("DigitBuilder Tool: No ANNIEEvent store!",v_error,verbosity); 
    return false;
  };

  // Load MC Particles information for this event
  auto get_mcparticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",
          fMCParticles);
  if(!get_mcparticles){
    Log("MCRecoEventLoader:: Tool: Error retrieving MCParticles from ANNIEEvent!",
            v_error,verbosity);
    return false;
  }


  ///Get MC Particle information
  this->FindTrueVertexFromMC();
  if (fGetPiKInfo) this->FindPionKaonCountFromMC();
  
  this->PushIBDInfo();

  this->PushTrueVertex(true);
  this->PushTrueStopVertex(true);
  this->PushTrueMuonEnergy(TrueMuonEnergy);
  //std::cout <<"MCRecoEventLoader: Pushing true muon energy "<<TrueMuonEnergy<<std::endl;
  this->PushTrueWaterTrackLength(WaterTrackLength);
  this->PushTrueMRDTrackLength(MRDTrackLength);
  this->PushProjectedMrdHit(projectedmrdhit);

  return true;
}


bool MCRecoEventLoader::Finalise(){
  delete fMuonStartVertex;
  delete fMuonStopVertex;
  if(verbosity>0) cout<<"MCRecoEventLoader exitting"<<endl;
  return true;
}

void MCRecoEventLoader::Reset() {
  // Reset 
  fMuonStartVertex->Reset();
  fMuonStopVertex->Reset();
  TrueMuonEnergy = -9999.;
  WaterTrackLength = -9999.;
  MRDTrackLength = -9999.;
  projectedmrdhit = false;  
}

void MCRecoEventLoader::FindTrueVertexFromMC() {
  
  // loop over the MCParticles to find the highest enery primary muon
  // MCParticles is a std::vector<MCParticle>
  MCParticle primarymuon;  // primary muon
  bool mufound=false;
  if(fMCParticles){
    Log("MCRecoEventLoader::  Tool: Num MCParticles = "+to_string(fMCParticles->size()),v_message,verbosity);
    for(unsigned int particlei=0; particlei<fMCParticles->size(); particlei++){
      MCParticle aparticle = fMCParticles->at(particlei);
      //if(v_debug<verbosity) aparticle.Print();       // print if we're being *really* verbose
      if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
      if (fDoParticleSelection){
        if(aparticle.GetPdgCode()!=fParticleID) continue;       // not a muon
        primarymuon = aparticle;                       // note the particle
        mufound=true;                                  // note that we found it
        m_data->Stores.at("RecoEvent")->Set("PdgPrimary",fParticleID);  //save the primary particle pdg code to the RecoEvent store
        break;                                         // won't have more than one primary muon
      } else {
	//Accept both electrons and muons as primary particles, if no selection is specified
        if( aparticle.GetPdgCode()!=11 && aparticle.GetPdgCode()!=13) continue;
	primarymuon = aparticle;
	mufound=true;
	m_data->Stores.at("RecoEvent")->Set("PdgPrimary",aparticle.GetPdgCode());
	break;
      }
    }
  } else {
    Log("MCRecoEventLoader::  Tool: No MCParticles in the event!",v_error,verbosity);
  }
  if(not mufound){
    Log("MCRecoEventLoader::  Tool: No muon in this event",v_warning,verbosity);
    return;
  }
  
  // retrieve desired information from the particle
  Position muonstartpos = primarymuon.GetStartVertex();    // only true if the muon is primary
  double muonstarttime = primarymuon.GetStartTime();
  Position muonstoppos = primarymuon.GetStopVertex();    // only true if the muon is primary
  double muonstoptime = primarymuon.GetStopTime();
  Direction muondirection = primarymuon.GetStartDirection();
  
  TrueMuonEnergy = primarymuon.GetStartEnergy();
   //std::cout <<"MCRecoEventLoader: FindTrueVertexFromMC: TrueEnergy: "<<TrueMuonEnergy<<std::endl;

  // MCParticleProperties tool fills in MRD track in m, but
  // Water track in cm...
  MRDTrackLength = primarymuon.GetTrackLengthInMrd()*100.;
  WaterTrackLength = primarymuon.GetTrackLengthInTank();

  //std::cout <<"MCRecoEventLoader: Muon start position: ("<<muonstartpos.X()<<","<<muonstartpos.Y()<<","<<muonstartpos.Z()<<")"<<std::endl;
  // set true vertex
  // change unit
  muonstartpos.UnitToCentimeter(); // convert unit from meter to centimeter
  muonstoppos.UnitToCentimeter(); // convert unit from meter to centimeter
  // change coordinate for muon start vertex
  muonstartpos.SetY(muonstartpos.Y()+yshift);
  muonstartpos.SetZ(muonstartpos.Z()+zshift);
  //std::cout <<"MCRecoEventLoader: NEW Muon start position: ("<<muonstartpos.X()<<","<<muonstartpos.Y()<<","<<muonstartpos.Z()<<")"<<std::endl;
  fMuonStartVertex->SetVertex(muonstartpos, muonstarttime);
  fMuonStartVertex->SetDirection(muondirection);
  //  charge coordinate for muon stop vertex
  muonstoppos.SetY(muonstoppos.Y()+yshift);
  muonstoppos.SetZ(muonstoppos.Z()+zshift);
  fMuonStopVertex->SetVertex(muonstoppos, muonstoptime); 
  
  logmessage = "  trueVtx = (" +to_string(muonstartpos.X()) + ", " + to_string(muonstartpos.Y()) + ", " + to_string(muonstartpos.Z()) +", "+to_string(muonstarttime)+ "\n"
            + "           " +to_string(muondirection.X()) + ", " + to_string(muondirection.Y()) + ", " + to_string(muondirection.Z()) + ") " + "\n";
  
  Log(logmessage,v_debug,verbosity);
	logmessage = "  muonStop = ("+to_string(muonstoppos.X()) + ", " + to_string(muonstoppos.Y()) + ", " + to_string(muonstoppos.Z()) + ") "+ "\n";
	Log(logmessage,v_debug,verbosity);

  //get information whether the extended particle trajectory were to hit the MRD
  projectedmrdhit = primarymuon.GetProjectedHitMrd();

}


void MCRecoEventLoader::FindPionKaonCountFromMC() {

  Log("MCRecoEventLoader: Find PionKaonCountFromMC",v_message,verbosity);
  
  // loop over the MCParticles to find the highest enery primary muon
  // MCParticles is a std::vector<MCParticle>
  bool pionfound=false;
  bool kaonfound=false;
  int pi0count = 0;
  int pipcount = 0;
  int pimcount = 0;
  int K0count = 0;
  int Kpcount = 0;
  int Kmcount = 0;

  //set up number of rings to 0 before counting
  int nprimary = 0;
  int nsecondary = 0;
  int nrings = 0;
  std::vector<unsigned int> index_particles_ring;

  if(fMCParticles){
    Log("MCRecoEventLoader::  Tool: Num MCParticles = "+to_string(fMCParticles->size()),v_message,verbosity);
    for(unsigned int particlei=0; particlei<fMCParticles->size(); particlei++){
      MCParticle aparticle = fMCParticles->at(particlei);
      //if(v_debug<verbosity) aparticle.Print();       // print if we're being *really* verbose
      if(aparticle.GetParentPdg()==0) {                //primary particle
        nprimary++;
        if (TMath::Abs(aparticle.GetPdgCode())==11){
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(11)) {nrings++; index_particles_ring.push_back(particlei);}
        } 
        if (TMath::Abs(aparticle.GetPdgCode())==13){
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(13)) {nrings++; index_particles_ring.push_back(particlei);}
        }
        if(aparticle.GetPdgCode()==111){               // is a primary pi0
          pionfound = true;
          pi0count++;
          nrings+=2; 
          index_particles_ring.push_back(particlei);}
        if(aparticle.GetPdgCode()==211){               // is a primary pi+
          pionfound = true;
          pipcount++;
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(211)) {nrings++; index_particles_ring.push_back(particlei);}
        }
        if(aparticle.GetPdgCode()==-211){               // is a primary pi-
          pionfound = true;
          pimcount++;
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-211)) {nrings++; index_particles_ring.push_back(particlei);}
        }
        if(aparticle.GetPdgCode()==311){               // is a primary K0
          kaonfound = true;
          K0count++;
        }
        if(aparticle.GetPdgCode()==321){               // is a primary K+
          kaonfound = true;
          Kpcount++;
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(321)) {nrings++; index_particles_ring.push_back(particlei);}
        }
        if(aparticle.GetPdgCode()==-321){               // is a primary K-
          kaonfound = true;
          Kmcount++;
          if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-321)) {nrings++; index_particles_ring.push_back(particlei);}
        }
      } else {                                          // not a primary particle
        nsecondary++;
        Log("MCRecoEventLoader: Secondary particle with pdg "+std::to_string(aparticle.GetPdgCode()),v_debug,verbosity);
 	//don't count rings from secondary particles for now (should we?)
      }
    }
  } else {
    Log("MCRecoEventLoader::  Tool: No MCParticles in the event!",v_error,verbosity);
  }
  if(not pionfound){
    Log("MCRecoEventLoader::  Tool: No primary pions in this event",v_warning,verbosity);
  }
  if(not kaonfound){
    Log("MCRecoEventLoader::  Tool: No kaons in this event",v_warning,verbosity);
  }
  if (fGetNRings){
    Log("MCRecoEventLoader: Found "+std::to_string(nrings)+" rings in this event, from "+std::to_string(nprimary)+" primary particles and "+std::to_string(nsecondary)+" secondary particles.");
  }
  //Fill in pion counts for this event
  m_data->Stores.at("RecoEvent")->Set("MCPi0Count", pi0count);
  m_data->Stores.at("RecoEvent")->Set("MCPiPlusCount", pipcount);
  m_data->Stores.at("RecoEvent")->Set("MCPiMinusCount", pimcount);
  m_data->Stores.at("RecoEvent")->Set("MCK0Count", K0count);
  m_data->Stores.at("RecoEvent")->Set("MCKPlusCount", Kpcount);
  m_data->Stores.at("RecoEvent")->Set("MCKMinusCount", Kmcount);
  if (fGetNRings) {
    m_data->Stores.at("RecoEvent")->Set("NRings",nrings);
    m_data->Stores.at("RecoEvent")->Set("IndexParticlesRing",index_particles_ring);
  }

}


void MCRecoEventLoader::PushTrueVertex(bool savetodisk) {
  Log("MCRecoEventLoader Tool: Push true vertex to the RecoEvent store",v_message,verbosity);
  m_data->Stores.at("RecoEvent")->Set("TrueVertex", fMuonStartVertex, savetodisk); 
}


void MCRecoEventLoader::PushTrueStopVertex(bool savetodisk) {
  Log("MCRecoEventLoader Tool: Push true stop vertex to the RecoEvent store",v_message,verbosity);
  m_data->Stores.at("RecoEvent")->Set("TrueStopVertex", fMuonStopVertex, savetodisk); 
}

void MCRecoEventLoader::PushTrueMuonEnergy(double MuE) {
	Log("MCRecoEventLoader Tool: Push true muon energy to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("TrueMuonEnergy", MuE);  ///> Add digits to RecoEvent
}

void MCRecoEventLoader::PushTrueWaterTrackLength(double WaterT) {
	Log("MCRecoEventLoader Tool: Push true track length in tank to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("TrueTrackLengthInWater", WaterT);  ///> Add digits to RecoEvent
}

void MCRecoEventLoader::PushTrueMRDTrackLength(double MRDT) {
	Log("MCRecoEventLoader Tool: Push true track length in MRD to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("TrueTrackLengthInMRD", MRDT);  ///> Add digits to RecoEvent
}

void MCRecoEventLoader::PushProjectedMrdHit(bool projectedmrdhit){
  Log("MCRecoEventLoader Tool: Push projected Mrd Hit",v_message,verbosity);
  m_data->Stores.at("RecoEvent")->Set("ProjectedMRDHit", projectedmrdhit);  ///> Add digits to RecoEvent 
}

double MCRecoEventLoader::GetCherenkovThresholdE(int pdg_code) {
  Log("MCRecoEventLoader Tool: GetCherenkovThresholdE",v_message,verbosity);            ///> Calculate Cherenkov threshold energies depending on particle pdg
  double Ethr = pdgcodetomass[pdg_code]*sqrt(1/(1-1/(n*n)));
  return Ethr;
}

void MCRecoEventLoader::PushIBDInfo(){

  Log("MCRecoEventLoader Tool: PushIBDInfo",v_message,verbosity);

  int n_neutrons = 0;
  int n_gammas = 0;
  int n_positrons = 0;

  if(fMCParticles){
    Log("MCRecoEventLoader::  Tool: Num MCParticles = "+to_string(fMCParticles->size()),v_message,verbosity);
    for(unsigned int particlei=0; particlei<fMCParticles->size(); particlei++){
      MCParticle aparticle = fMCParticles->at(particlei);
      if(aparticle.GetParentPdg()==0) {                //primary particle
        int pdg = aparticle.GetPdgCode();
        double energy = aparticle.GetStartEnergy();
        if (pdg == 2112) n_neutrons++;
        if (pdg == -11 && energy < 100) n_positrons++;
        if (pdg == 22 && energy < 100) n_gammas++;
      }
    }
  }

  m_data->Stores.at("RecoEvent")->Set("NeutronCount",n_neutrons);
  m_data->Stores.at("RecoEvent")->Set("PositronCount",n_positrons);
  m_data->Stores.at("RecoEvent")->Set("GammaCount",n_gammas); 
  
}



