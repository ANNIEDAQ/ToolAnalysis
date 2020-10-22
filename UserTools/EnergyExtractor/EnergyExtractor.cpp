#include "EnergyExtractor.h"

EnergyExtractor::EnergyExtractor():Tool(){}

bool EnergyExtractor::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  m_data= &data; //assigning transient data pointer
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  /////////////////////////////////////////////////////////////////
  
  //Get configuration variables
  m_variables.Get("OutputFile",cnn_outpath);
  m_variables.Get("Verbosity",verbosity);
  m_variables.Get("SaveNeutrino",save_neutrino);
  m_variables.Get("SaveElectron",save_elec);
  m_variables.Get("SaveGamma",save_gamma);
  m_variables.Get("SaveMuon",save_muon);
  m_variables.Get("SavePion",save_pion);
  m_variables.Get("SaveKaon",save_kaon);
  m_variables.Get("SaveNeutron",save_neutron);
  m_variables.Get("SaveVisible",save_visible);  

  //Get geometry and mass map
  m_data->CStore.Get("PdgMassMap",pdgcodetomass);
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

  //Emplace filler contents for the Energystore for all particles
  Energystore.emplace("e-",filler);
  Energystore.emplace("e+",filler);
  Energystore.emplace("mu-",filler);
  Energystore.emplace("mu+",filler);
  Energystore.emplace("gamma",filler);
  Energystore.emplace("pi-",filler);
  Energystore.emplace("pi+",filler);
  Energystore.emplace("pi0",filler);
  Energystore.emplace("K+",filler);
  Energystore.emplace("K-",filler);
  Energystore.emplace("K0",filler);

  //Create respective csv-files in which the energy information will be stored
  std::string str_csv = ".csv";
  std::string str_EnergyNeutrino = "_EnergyNeutrino";
  std::string str_EnergyPion = "_EnergyPion";
  std::string str_EnergyKaon = "_EnergyKaon";
  std::string str_EnergyMuon = "_EnergyMuon";
  std::string str_EnergyElectron = "_EnergyElectron";
  std::string str_EnergyGamma = "_EnergyGamma";
  std::string str_NeutronNumber = "_NeutronNumber";
  std::string str_VisibleEnergy = "_VisibleEnergy";

  std::string csvfile_EnergyNeutrino = cnn_outpath + str_EnergyNeutrino + str_csv;
  std::string csvfile_EnergyPion = cnn_outpath + str_EnergyPion + str_csv;
  std::string csvfile_EnergyKaon = cnn_outpath + str_EnergyKaon + str_csv;
  std::string csvfile_EnergyMuon = cnn_outpath + str_EnergyMuon + str_csv;
  std::string csvfile_EnergyElectron = cnn_outpath + str_EnergyElectron + str_csv;
  std::string csvfile_EnergyGamma = cnn_outpath + str_EnergyGamma + str_csv;
  std::string csvfile_NeutronNumber = cnn_outpath + str_NeutronNumber + str_csv;
  std::string csvfile_VisibleEnergy = cnn_outpath + str_VisibleEnergy + str_csv;

  if (save_neutrino) outfile_EnergyNeutrino.open(csvfile_EnergyNeutrino.c_str());
  if (save_pion) outfile_EnergyPion.open(csvfile_EnergyPion.c_str());
  if (save_kaon) outfile_EnergyKaon.open(csvfile_EnergyKaon.c_str());
  if (save_muon) outfile_EnergyMuon.open(csvfile_EnergyMuon.c_str());
  if (save_gamma) outfile_EnergyGamma.open(csvfile_EnergyGamma.c_str());
  if (save_elec) outfile_EnergyElectron.open(csvfile_EnergyElectron.c_str());
  if (save_neutron) outfile_NeutronNumber.open(csvfile_NeutronNumber.c_str());
  if (save_visible) outfile_VisibleEnergy.open(csvfile_VisibleEnergy.c_str());

  return true;
}


bool EnergyExtractor::Execute()
{

  // Was Event Selection cut passed?
  bool passed_eventselection;
  m_data->Stores["RecoEvent"]->Get("EventCutStatus",passed_eventselection);
  std::cout <<"EnergyExtractor tool: passed eventselection: "<<passed_eventselection<<std::endl;
  if (!passed_eventselection) {
    return true;
  }

  // Check if "ANNIEEvent" exists
  auto get_annieevent = m_data->Stores.count("ANNIEEvent");
  if(!get_annieevent){
    Log("DigitBuilder Tool: No ANNIEEvent store!",v_error,verbosity);
    return false;
  }

  // Get ANNIEEvent-related variables
  m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
  m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits", MCLAPPDHits);
  auto get_mcparticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",fMCParticles);
  if(!get_mcparticles){
    Log("EnergyExtractor:: Tool: Error retrieving MCParticles from ANNIEEvent!",v_error,verbosity);
    return false;
  }

  // Evaluate MCHits object
  int vectsize = MCHits->size();
  if (verbosity > 1) std::cout <<"EnergyExtractor tool: MCHits size: "<<vectsize<<std::endl;
  double total_charge_pmts=0;
  for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits)
  {
    unsigned long chankey = apair.first;
    if (verbosity > 3) std::cout <<"chankey: "<<chankey<<std::endl;
    Detector* thistube = geom->ChannelToDetector(chankey);
    unsigned long detkey = thistube->GetDetectorID();
    if (verbosity > 3) std::cout <<"detkey: "<<detkey<<std::endl;
    if (thistube->GetDetectorElement()=="Tank")
    {
      if (thistube->GetTankLocation()=="OD") continue;
      std::vector<MCHit>& Hits = apair.second;
      for (MCHit &ahit : Hits)
      {
        if (ahit.GetTime()>-10. && ahit.GetTime()<40.)  {
          total_charge_pmts += ahit.GetCharge();
        }
      }
    }
  }

  // Evaluate MCLAPPDHits object
  vectsize = MCLAPPDHits->size();
  if (verbosity > 1) std::cout <<"EnergyExtractor tool: MCLAPPDHits size: "<<vectsize<<std::endl;
  double total_charge_lappds=0;
  for(std::pair<unsigned long, std::vector<MCLAPPDHit>>&& apair : *MCLAPPDHits)
  {
    unsigned long chankey = apair.first;
    if (verbosity > 3) std::cout <<"chankey: "<<chankey<<std::endl;
    Detector* thistube = geom->ChannelToDetector(chankey);
    unsigned long detkey = thistube->GetDetectorID();
    if (verbosity > 3) std::cout <<"detkey: "<<detkey<<std::endl;
    if (thistube->GetDetectorElement()=="LAPPD")
    {
      std::vector<MCLAPPDHit>& Hits = apair.second;
      for (MCLAPPDHit &ahit : Hits)
      {
        // Replace with slightly different time cut
        if (ahit.GetTime()>-10. && ahit.GetTime()<40.)  {
          //total_charge_lappds += ahit.GetCharge();	//put back in once charge is correctly implemented for LAPPDs
          total_charge_lappds += 1.;
        }
      }
    }
  }

  // Find true energy of particles from Monte Carlo information (MCParticles object)
  this->FindTrueEnergyFromMC();

  // Save neutrino energy (if desired)
  Float_t neutrino_Energy=-2;
  bool isok;
  MCParticle neutrino;
  isok = m_data->Stores["ANNIEEvent"]->Get("NeutrinoParticle", neutrino);
  if (isok) neutrino_Energy = neutrino.GetStartEnergy();
  if (save_neutrino) outfile_EnergyNeutrino << neutrino_Energy << endl;

  // Save pion energies (if desired)
  if (save_pion){
    outfile_EnergyPion << pimcount<<"," << pipcount << "," << pi0count << ",";
    for (int i=0 ; i < (int) Energystore["pi-"].size() ; i++)
    {
      outfile_EnergyPion << Energystore["pi-"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["pi+"].size() ; i++)
    {
      outfile_EnergyPion << Energystore["pi+"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["pi0"].size() ; i++)
    {
      outfile_EnergyPion << Energystore["pi0"][i] << ",";
    }
    outfile_EnergyPion << endl;
  }

  // Save kaon energies (if desired)
  if (save_kaon){
    outfile_EnergyKaon << Kmcount<<"," << Kpcount << "," << K0count << "," ;
    for (int i=0 ; i < (int) Energystore["K-"].size() ; i++)
    {
      outfile_EnergyKaon << Energystore["K-"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["K+"].size() ; i++)
    {
      outfile_EnergyKaon << Energystore["K+"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["K0"].size() ; i++)
    {
      outfile_EnergyKaon << Energystore["K0"][i] << ",";
    }
    outfile_EnergyKaon << endl;
  }

  // Save muon energies (if desired)
  if (save_muon){
    outfile_EnergyMuon << mcount <<"," << mpcount << ",";
    for (int i=0 ; i < (int) Energystore["mu-"].size() ; i++)
    {
      outfile_EnergyMuon << Energystore["mu-"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["mu+"].size() ; i++)
    {
      outfile_EnergyMuon << Energystore["mu+"][i] << ",";
    }
    outfile_EnergyMuon <<endl;
  }

  // Save electron energies (if desired)
  if (save_elec) {
    outfile_EnergyElectron << ecount<<","<< epcount << ",";
    for (int i=0 ; i < (int) Energystore["e-"].size() ; i++)
    {
      outfile_EnergyElectron << Energystore["e-"][i] << ",";
    }
    for (int i=0 ; i < (int) Energystore["e+"].size() ; i++)
    {
      outfile_EnergyElectron << Energystore["e+"][i] << ",";
    }
    outfile_EnergyElectron << endl;
  }

  // Save gamma energies (if desired)
  if (save_gamma) {
    outfile_EnergyGamma << gammacount << ",";
    for (int i=0; i < (int) Energystore["gamma"].size(); i++){
      outfile_EnergyGamma << Energystore["gamma"][i] << ",";
    }
    outfile_EnergyGamma << endl;
  }

  // Save neutron number (if desired)
  if (save_neutron) outfile_NeutronNumber << ncount << endl;

  // Save visible energy (if desired)
  if (save_visible) outfile_VisibleEnergy << total_charge_pmts << "," << total_charge_lappds << endl;

  // Clear energy store map
  Energystore["pi+"].clear();
  Energystore["pi-"].clear();
  Energystore["K+"].clear();
  Energystore["K-"].clear();
  Energystore["e+"].clear();
  Energystore["e-"].clear();
  Energystore["mu+"].clear();
  Energystore["mu-"].clear();
  Energystore["gamma"].clear();
  Energystore["K0"].clear();
  Energystore["pi0"].clear();
  pi0count = 0;
  pipcount = 0;
  pimcount = 0;
  K0count = 0;
  Kpcount = 0;
  Kmcount = 0;
  ecount = 0;
  epcount = 0;
  mcount = 0;
  mpcount = 0;
  gammacount = 0;
  ncount = 0;

  return true;
}


bool EnergyExtractor::Finalise(){
  
  if (save_neutrino) outfile_EnergyNeutrino.close();
  if (save_pion) outfile_EnergyPion.close();
  if (save_kaon) outfile_EnergyKaon.close();
  if (save_muon) outfile_EnergyMuon.close();
  if (save_elec) outfile_EnergyElectron.close();
  if (save_gamma) outfile_EnergyGamma.close();
  if (save_neutron) outfile_NeutronNumber.close();
  if (save_visible) outfile_VisibleEnergy.close();

  return true;
}


void EnergyExtractor::FindTrueEnergyFromMC() {

  // Find true energies for different particle types (if they occur in the event)
  // Store true energies in a map that connects particle types and a vector of the energies of particles of said type

  if(fMCParticles)
  {
    Log("EnergyExtractor::  Tool: Num MCParticles = "+to_string(fMCParticles->size()),v_message,verbosity);
    for (int particlei=0; particlei < (int) fMCParticles->size(); particlei++)
    {
      MCParticle aparticle = fMCParticles->at(particlei);
      if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
      TrueEnergy = aparticle.GetStartEnergy();
      Log("EnergyExtractor tool: True Energy:  "+std::to_string(TrueEnergy)+" and PDG: "+std::to_string(aparticle.GetPdgCode())+" and Cherenkov Threshold: "+std::to_string(GetCherenkovThresholdE(aparticle.GetPdgCode())),v_message,verbosity);

      // Pions
      if(aparticle.GetPdgCode()==111){               // is a primary pi0
        Energystore["pi0"].push_back(TrueEnergy);  pi0count++;	//Pi0 decays to two gammas --> no Cherenkov threshold check
      }
      if(aparticle.GetPdgCode()==211){               // is a primary pi+
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(211)) {Energystore["pi+"].push_back(TrueEnergy);pipcount++;}
      }
      if(aparticle.GetPdgCode()==-211){               // is a primary pi-
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-211)) {Energystore["pi-"].push_back(TrueEnergy);pimcount++;}
      }
      // Kaons
      if(aparticle.GetPdgCode()==311){               // is a primary K0
        Energystore["K0"].push_back(TrueEnergy);K0count++; // K0 is neutral --> No Cherenkv threshold check
      }
      if(aparticle.GetPdgCode()==321){               // is a primary K+
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(321)) {Energystore["K+"].push_back(TrueEnergy);Kpcount++;}
      }
      if(aparticle.GetPdgCode()==-321){               // is a primary K-
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-321)) {Energystore["K-"].push_back(TrueEnergy);Kmcount++;}
      }
      // Electrons
      if(aparticle.GetPdgCode()==11){               // is a primary e-
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(11)) {Energystore["e-"].push_back(TrueEnergy);ecount++;}
      }
      if(aparticle.GetPdgCode()==-11){               // is a primary e+
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-11)) {Energystore["e+"].push_back(TrueEnergy);epcount++;}
      }
      // Muons
      if(aparticle.GetPdgCode()==13){               // is a primary mu-
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(13)) {Energystore["mu-"].push_back(TrueEnergy);mcount++;}
      }
      if(aparticle.GetPdgCode()==-13){               // is a primary mu+
        if (aparticle.GetStartEnergy() > GetCherenkovThresholdE(-13)) {Energystore["mu+"].push_back(TrueEnergy);mpcount++;}
      }
      // Neutrons
      if(aparticle.GetPdgCode()==2112){               // is a primary neutron
        ncount++;
      }
      // Gammas
      if(aparticle.GetPdgCode() == 22) {		// is a primary gamma
        Energystore["gamma"].push_back(TrueEnergy);  gammacount++;
      }
    }

  } else {
    Log("EnergyExtractor::  Tool: No MCParticles in the event!",v_error,verbosity);
  }
}

double EnergyExtractor::GetCherenkovThresholdE(int pdg_code) {
  Log("EnergyExtractor Tool: GetCherenkovThresholdE",v_message,verbosity);
  Log("EnergyExtractor Tool: PDG code: "+std::to_string(pdg_code)+", mass: "+std::to_string(pdgcodetomass[pdg_code]),v_message,verbosity);             ///> Calculate Cherenkov threshold energies depending on particle pdg
  double Ethr = pdgcodetomass[pdg_code]*sqrt(1./(1-1./(n*n)));
  return Ethr;
}
