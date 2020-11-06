#include "DigitBuilderROOT.h"

DigitBuilderROOT::DigitBuilderROOT():Tool(){}


bool DigitBuilderROOT::Initialise(std::string configfile, DataModel &data){
  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
	m_variables.Get("InputFile",NtupleFile);
	m_variables.Get("verbosity",verbosity);

  Log("---------DigitBuilderROOT Tool initializing ---------- " ,v_message,verbosity);

  int recoeventexists = m_data->Stores.count("RecoEvent");
  if(recoeventexists==0) m_data->Stores["RecoEvent"] = new BoostStore(false,2);

  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(annieeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);

  if(verbosity>4) std::cout << "DigitBuilderROOT Tool: loading ntuple file " << NtupleFile << std::endl;

  P2Chain = new TChain("phaseII");
  P2Chain->Add(NtupleFile.c_str());
 
  fMuonVertex = new RecoVertex();
  fDigitList = new std::vector<RecoDigit>;
  
  fDigitX = new std::vector<double>;
  fDigitY = new std::vector<double>;
  fDigitZ = new std::vector<double>;
  fDigitT = new std::vector<double>;
  fDigitQ = new std::vector<double>;    
  fDigitType  = new std::vector<int>;
  fDigitDetID = new std::vector<int>;
  Log("DigitBuilderROOT Tool: Setting all branch addresses now",v_debug,verbosity);
  //Hit information (PMT and LAPPD)
  //Always output in Phase II Reco Tree
  P2Chain->SetBranchAddress("mcEntryNumber",&fMCEventNum);
  P2Chain->SetBranchAddress("triggerNumber",&fMCTriggerNum);
  P2Chain->SetBranchAddress("eventNumber",&fEventNumber);
  P2Chain->SetBranchAddress("runNumber",&fRunNumber);
  P2Chain->SetBranchAddress("subrunNumber",&fSubRunNumber);
  P2Chain->SetBranchAddress("nhits",&fNhits);
  P2Chain->SetBranchAddress("digitX",&fDigitX);
  P2Chain->SetBranchAddress("digitY",&fDigitY);
  P2Chain->SetBranchAddress("digitZ",&fDigitZ);
  P2Chain->SetBranchAddress("digitT",&fDigitT);
  P2Chain->SetBranchAddress("digitQ",&fDigitQ);
  P2Chain->SetBranchAddress("digitType", &fDigitType);
  P2Chain->SetBranchAddress("digitDetID", &fDigitDetID);

  //Muon Truth Information
  P2Chain->SetBranchAddress("trueVtxX",&fTrueVtxX);
  P2Chain->SetBranchAddress("trueVtxY",&fTrueVtxY);
  P2Chain->SetBranchAddress("trueVtxZ",&fTrueVtxZ);
  P2Chain->SetBranchAddress("trueVtxTime",&fTrueVtxTime);
  P2Chain->SetBranchAddress("trueDirX",&fTrueDirX);
  P2Chain->SetBranchAddress("trueDirY",&fTrueDirY);
  P2Chain->SetBranchAddress("trueDirZ",&fTrueDirZ);
  P2Chain->SetBranchAddress("trueMuonEnergy",&fTrueMuonEnergy);
  P2Chain->SetBranchAddress("trueTrackLengthInWater",&fTrueTrackLengthInWater);
  P2Chain->SetBranchAddress("trueTrackLengthInMRD",&fTrueTrackLengthInMRD);

  //Pion/Kaon Count Truth Information
  P2Chain->SetBranchAddress("Pi0Count",&fPi0Count);
  P2Chain->SetBranchAddress("PiPlusCount",&fPiPlusCount);
  P2Chain->SetBranchAddress("PiMinusCount",&fPiMinusCount);
  P2Chain->SetBranchAddress("K0Count",&fK0Count);
  P2Chain->SetBranchAddress("KPlusCount",&fKPlusCount);
  P2Chain->SetBranchAddress("KMinusCount",&fKMinusCount);

  //Event status information from previous run
  P2Chain->SetBranchAddress("eventStatusApplied",&fEventStatusApplied);
  P2Chain->SetBranchAddress("eventStatusFlagged",&fEventStatusFlagged);

  TotalEntries = P2Chain->GetEntries();
  EntryNum = 0;
  return true;
}


bool DigitBuilderROOT::Execute(){
  this->Reset();
  Log("DigitBuilderROOT Tool: Getting EntryNum " + std::to_string(EntryNum),v_debug,verbosity);
  P2Chain->GetEvent(EntryNum);

  Log("DigitBuilderROOT Tool: loading true muon information from phaseII ntuple",v_debug,verbosity);
  //Load the True Muon Vertex information
  Position muonstartpos(fTrueVtxX, fTrueVtxY, fTrueVtxZ);
  Direction muondirection(fTrueDirX, fTrueDirY, fTrueDirZ);
  double muonstarttime = fTrueVtxTime;
  fMuonVertex->SetVertex(muonstartpos, muonstarttime);
  fMuonVertex->SetDirection(muondirection);
  this->PushTrueVertex(true);

  this->PushTrueWaterTrackLength(fTrueTrackLengthInWater);
  this->PushTrueMRDTrackLength(fTrueTrackLengthInMRD);
  m_data->Stores.at("RecoEvent")->Set("TrueMuonEnergy", fTrueMuonEnergy);

  Log("DigitBuilderROOT Tool: Looping through and loading all digits",v_debug,verbosity);
  //Loop through all Digits and load their information as RecoDigits
  for(int j=0;j<(int)fDigitType->size();j++){
    int region = -999;
    double calT = fDigitT->at(j);
    double calQ = fDigitQ->at(j);
    Position pos_reco(fDigitX->at(j), fDigitY->at(j),fDigitZ->at(j));
    int digitType = fDigitType->at(j);
    int digitID = fDigitDetID->at(j);
    RecoDigit aDigit(region, pos_reco, calT, calQ,digitType,digitID);
    fDigitList->push_back(aDigit); 
  }
  this->PushRecoDigits(true);

  Log("DigitBuilderROOT Tool: Loading all pion/kaon counts",v_debug,verbosity);
  //Fill in pion counts for this event
  m_data->Stores.at("RecoEvent")->Set("MCPi0Count", fPi0Count);
  m_data->Stores.at("RecoEvent")->Set("MCPiPlusCount", fPiPlusCount);
  m_data->Stores.at("RecoEvent")->Set("MCPiMinusCount", fPiMinusCount);
  m_data->Stores.at("RecoEvent")->Set("MCfK0Count", fK0Count);
  m_data->Stores.at("RecoEvent")->Set("MCfKPlusCount", fKPlusCount);
  m_data->Stores.at("RecoEvent")->Set("MCfKMinusCount", fKMinusCount);

  EntryNum++;
  if(EntryNum==TotalEntries){
    m_data->vars.Set("StopLoop",1);
    Log("DigitBuilderROOT Tool: All entries in ROOT file read.",v_message,verbosity);
  }
  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Set("MCEventNum",fMCEventNum);  
  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Set("MCTriggernum",fMCTriggerNum); 
  
  // ANNIE Event number
  m_data->Stores.at("ANNIEEvent")->Set("EventNumber",fEventNumber);



  // Event Status information
  Log("DigitBuilderROOT Tool: Loading Event Status Information",v_debug,verbosity);
  m_data->Stores.at("RecoEvent")->Set("EventFlagApplied",fEventStatusApplied);
  m_data->Stores.at("RecoEvent")->Set("EventFlagged",fEventStatusFlagged);

  if ((fEventStatusApplied & fEventStatusFlagged)==0) fEventCutStatus = true;
  else fEventCutStatus = false;
 
  m_data->Stores.at("RecoEvent")->Set("EventCutStatus", fEventCutStatus);
 
  return true;
}


bool DigitBuilderROOT::Finalise(){
  delete P2Chain;
  delete fMuonVertex;
  delete fDigitList;
  delete fDigitX;
  delete fDigitY;
  delete fDigitZ;
  delete fDigitT;
  delete fDigitQ;
  delete fDigitType;
  delete fDigitDetID;
  return true;
}

void DigitBuilderROOT::Reset(){
  fMuonVertex->Reset();
  fDigitX->clear();
  fDigitY->clear();
  fDigitZ->clear();
  fDigitT->clear();
  fDigitQ->clear();
  fDigitType->clear();
  fDigitDetID->clear();
  fDigitList->clear();
}

void DigitBuilderROOT::PushTrueVertex(bool savetodisk) {
  Log("DigitBuilderROOT Tool: Push true vertex to the RecoEvent store",v_message,verbosity);
  m_data->Stores.at("RecoEvent")->Set("TrueVertex", fMuonVertex, savetodisk); 
}

void DigitBuilderROOT::PushRecoDigits(bool savetodisk) {
	Log("DigitBuilder Tool: Push reconstructed digits to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("RecoDigit", fDigitList, savetodisk);  ///> Add digits to RecoEvent
}

void DigitBuilderROOT::PushTrueWaterTrackLength(double WaterT) {
	Log("DigitBuilderROOT Tool: Push true track length in tank to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("TrueTrackLengthInWater", WaterT);  ///> Add digits to RecoEvent
}

void DigitBuilderROOT::PushTrueMRDTrackLength(double MRDT) {
	Log("DigitBuilderROOT Tool: Push true track length in MRD to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("TrueTrackLengthInMRD", MRDT);  ///> Add digits to RecoEvent
}

