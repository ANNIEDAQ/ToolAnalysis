#include "DigitBuilderDoE.h"

using namespace ROOT::Math;

static DigitBuilderDoE* fgDigitBuilderDoE = 0;
DigitBuilderDoE* DigitBuilderDoE::Instance()
{
  if( !fgDigitBuilderDoE ){
    fgDigitBuilderDoE = new DigitBuilderDoE();
  }

  return fgDigitBuilderDoE;
}

DigitBuilderDoE::DigitBuilderDoE():Tool(){}
DigitBuilderDoE::~DigitBuilderDoE() {
}



bool DigitBuilderDoE::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

	doetreename = "vertextreefiducialmrd";
  fPhotodetectorConfiguration = "All"; //Load all hits by default
  fHistoricOffset = 6670.;

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
	m_variables.Get("InputFile",doefilename);
  m_variables.Get("TreeName",doetreename);
	m_variables.Get("PhotoDetectorConfiguration", fPhotodetectorConfiguration);
	m_variables.Get("HistoricOffset", fHistoricOffset);
  
  // Make the ANNIEEvent Store if it doesn't exist
	int annieeventexists = m_data->Stores.count("ANNIEEvent");
	if(annieeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);
	
  // Make the RecoDigit Store if it doesn't exist
	int recoeventexists = m_data->Stores.count("RecoEvent");
	if(recoeventexists==0) m_data->Stores["RecoEvent"] = new BoostStore(false,2);
  
  //Construct objects
	DigitVertices = new std::vector<LorentzVector<PxPyPzE4D<double>>>;
  DigitWhichDet = new std::vector<std::string>;
	DigitCharges = new std::vector<double>;
	DigitIdArray = new std::vector<int>;
	MuonStartVertex = new TLorentzVector;
	MuonStopVertex = new TLorentzVector;
	MuonDirection = new TVector3;	

  doechain = new TChain(doetreename.c_str());
  doechain->Add(doefilename.c_str());
  NumEvents = doechain->GetEntries();
  std::cout<<"FILE LOADED INTO CHAIN, HAS " << NumEvents << "ENTRIES" << std::endl;

  /// Construct the other objects we'll be setting at event level,
  fMuonStartVertex = new RecoVertex();
  //fMuonStopVertex = new RecoVertex();
	fDigitList = new std::vector<RecoDigit>;

  EventNum = 0;
  MCEventNum = 0;
  //Now, link chain's branch addresses with those defined here
  this->LinkChain();
  return true;
}


bool DigitBuilderDoE::Execute(){
	this->Reset();
  // see if "ANNIEEvent" exists
	auto get_annieevent = m_data->Stores.count("ANNIEEvent");
	if(!get_annieevent){
		Log("DigitBuilder Tool: No ANNIEEvent store!",v_error,verbosity); 
		return false;
	};
	
	/// see if "RecoEvent" exists
 	auto get_recoevent = m_data->Stores.count("RecoEvent");
 	if(!get_recoevent){
  		Log("EventSelector Tool: No RecoEvent store!",v_error,verbosity); 
  		return false;
	};
  
  doechain->GetEntry(EventNum);

  //First, load digits from event into RecoDigit classes
  //Loop through and convert information to the right form for the
  //RecoEvent checks that will be done in EventSelectorDoE
  //
  int region=-999.;
	double calT;
	double calQ = 0.;
	Position pos_reco;
  int DigitId;
  int digitType;
  std::cout << "TOTAL DIGIT ARRAY SIZE: " << DigitIdArray->size() << std::endl;
  for (int i=0;i<(int)DigitIdArray->size();i++){
    if(DigitWhichDet->at(i) == "PMT8inch" && DigitCharges->at(i)>5.0){
      if((fPhotodetectorConfiguration == "PMT_only") || (fPhotodetectorConfiguration == "All")){
        //Load up a PMT RecoDigit's goods
        digitType = 0;
        pos_reco.SetX(DigitVertices->at(i).X());
        pos_reco.SetY(DigitVertices->at(i).Y() + 14.46469);
        pos_reco.SetZ(DigitVertices->at(i).Z() - 168.1);
        calQ = DigitCharges->at(i);
        calT = DigitVertices->at(i).T() - fHistoricOffset;
        DigitId = DigitIdArray->at(i);
  		  RecoDigit recoDigit(region, pos_reco, calT, calQ, digitType, DigitId);
  		  fDigitList->push_back(recoDigit);
      }
    }
    if(DigitWhichDet->at(i) == "lappd_v0"){
      if((fPhotodetectorConfiguration == "LAPPD_only") || (fPhotodetectorConfiguration == "All")){
        DigitId = DigitIdArray->at(i);
        if (DigitId == 236 || DigitId == 203 || DigitId == 268 || DigitId==231 || DigitId == 240){
          digitType = 1;
          pos_reco.SetX(DigitVertices->at(i).X());
          pos_reco.SetY(DigitVertices->at(i).Y() + 14.46469);
          pos_reco.SetZ(DigitVertices->at(i).Z() - 168.1);
          calQ = 1; //Charges not well modeled yet; don't use
          calT = DigitVertices->at(i).T() - fHistoricOffset;
					calT = frand.Gaus(calT, 0.1); // time is smeared with 100 ps time resolution. Harded-coded for now.
    		  RecoDigit recoDigit(region, pos_reco, calT, calQ, digitType, DigitId);
    		  fDigitList->push_back(recoDigit);
        }
      }
    }
  }
  //Now, we'll load in the necessary muon information
  Direction muondirection;
  Position muonstartpos;
  double muonstarttime;
  muondirection.SetX(MuonDirection->X());
  muondirection.SetY(MuonDirection->Y());
  muondirection.SetZ(MuonDirection->Z());
  muonstartpos.SetX(MuonStartVertex->X()); 
  muonstartpos.SetY(MuonStartVertex->Y() + 14.46469);
  muonstartpos.SetZ(MuonStartVertex->Z() - 168.1);
  muonstarttime = MuonStartVertex->T() - fHistoricOffset;
  fMuonStartVertex->SetVertex(muonstartpos, muonstarttime);
  fMuonStartVertex->SetDirection(muondirection);


	/// Push recodigits and muon truth info to RecoEvent
  this->PushRecoDigits(true); 
  this->PushTrueVertex(true);

  //EventSelectorDoE will have a rudimentary check if track went into MRD
  m_data->Stores.at("RecoEvent")->Set("TrackLengthInMRD",TrackLengthInMrd);
  //by default, all events pass the event cut status
  bool EventCutStatus = true;
  m_data->Stores.at("RecoEvent")->Set("EventCutStatus", EventCutStatus);
 
  m_data->Stores.at("ANNIEEvent")->Set("MCEventNum",MCEventNum);
  uint16_t MCTriggernum = 0;
	m_data->Stores.at("ANNIEEvent")->Set("MCTriggernum",MCTriggernum);
	m_data->Stores.at("ANNIEEvent")->Set("EventNumber",EventNum);
  EventNum++;
  MCEventNum++;
  return true;
}


bool DigitBuilderDoE::Finalise(){
  //WE NEED TO DELETE STUFF AND CLEAN UP... 
	delete DigitVertices;
  delete DigitWhichDet; 
  delete DigitCharges;
	delete DigitIdArray;
	delete MuonStartVertex;
	delete MuonStopVertex;
	delete MuonDirection;	
  delete doechain;
  delete fMuonStartVertex;
  //TODO: Track length in MRD deletion here?
  return true;
}

void DigitBuilderDoE::PushTrueVertex(bool savetodisk) {
  Log("EventSelector Tool: Push true vertex to the RecoEvent store",v_message,verbosity);
  m_data->Stores.at("RecoEvent")->Set("TrueVertex", fMuonStartVertex, savetodisk); 
}

void DigitBuilderDoE::PushRecoDigits(bool savetodisk) {
	Log("DigitBuilder Tool: Push reconstructed digits to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("RecoDigit", fDigitList, savetodisk);  ///> Add digits to RecoEvent
}

void DigitBuilderDoE::LinkChain(){
  //Do all the awesome linking here
 doechain->SetBranchAddress("MuonStartVertex", &MuonStartVertex);
 doechain->SetBranchAddress("MuonStopVertex", &MuonStopVertex);
 doechain->SetBranchAddress("TrackLengthInMrd", &TrackLengthInMrd);
 doechain->SetBranchAddress("MuonDirection", &MuonDirection);
 doechain->SetBranchAddress("DigitVertices", &DigitVertices);
 doechain->SetBranchAddress("DigitCharges", &DigitCharges);
 doechain->SetBranchAddress("DigitWhichDet", &DigitWhichDet);
 doechain->SetBranchAddress("DigitPmtId", &DigitIdArray);
}

void DigitBuilderDoE::Reset() {
	// Reset 
  fDigitList->clear();
  fMuonStartVertex->Reset();
}
