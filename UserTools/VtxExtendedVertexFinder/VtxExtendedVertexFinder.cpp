#include "VtxExtendedVertexFinder.h"

VtxExtendedVertexFinder::VtxExtendedVertexFinder():Tool(){}


bool VtxExtendedVertexFinder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  fTmin = -10.0;
  fTmax = 10.0;
  fUseTrueVertexAsSeed = false;
  fSeedGridFits = false;
  /// Get the Tool configuration variables
  m_variables.Get("UseTrueVertexAsSeed",fUseTrueVertexAsSeed);
  m_variables.Get("FitAllOnSeedGrid",fSeedGridFits);
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("FitTimeWindowMin", fTmin);
  m_variables.Get("FitTimeWindowMax", fTmax);
  
  /// Create extended vertex
  /// Note that the objects created by "new" must be added to the "RecoEvent" store. 
  /// The last tool SaveRecoEvent will delete these pointers and free the memory.
  /// If these pointers are not added to any store, the user has to delete the pointers
  /// and free the memory. 
  /// In this tool, the pointer is 
  fExtendedVertex = new RecoVertex();
  
  return true;
}

bool VtxExtendedVertexFinder::Execute(){
  Log("===========================================================================================",v_debug,verbosity);
  Log("VtxExtendedVertexFinder Tool: Executing",v_debug,verbosity);
  // Reset everything
  this->Reset();
  
  // check if event passes the cut
  bool EventCutstatus = false;
  auto get_evtstatus = m_data->Stores.at("RecoEvent")->Get("EventCutStatus",EventCutstatus);
  if(!get_evtstatus) {
    Log("Error: The PhaseITreeMaker tool could not find the Event selection status", v_error, verbosity);
    return false;	
  }

  if(!EventCutstatus) {
  	Log("Message: This event doesn't pass the event selection. ", v_message, verbosity);
    return true;	
  }
  
  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",fMCEventNum);  
  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggerNum); 
  
  // ANNIE Event number
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
  
  std::cout<<"event number = "<<fEventNumber<<std::endl;
	
  // Retrive digits from RecoEvent
  get_ok = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(not get_ok){
    Log("VtxExtendedVertexFinder  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
    return false;
  }
	
  // Load digits to VertexGeometry
  myvtxgeo = VertexGeometry::Instance();
  myvtxgeo->LoadDigits(fDigitList);
  // Do extended vertex (muon track) reconstruction using MC truth information
  if( fUseTrueVertexAsSeed ){
  Log("VtxExtendedVertexFinder Tool: Run vertex reconstruction using MC truth information",v_message,verbosity);
  // get truth vertex information 
  auto get_truevtx = m_data->Stores.at("RecoEvent")->Get("TrueVertex", fTrueVertex);
  if(!get_truevtx){ 
    Log("VtxExtendedVertexFinder Tool: Error retrieving TrueVertex from RecoEvent!",v_error,verbosity); 
    return false; 
  }
  if(verbosity>3){
    Position muonstartpos = fTrueVertex->GetPosition();
    double muonstarttime = fTrueVertex->GetTime();
    Direction muondirection = fTrueVertex->GetDirection();
    std::cout << "VtxExtendedVertexFinder Tool: Printing muon info going into Minuit" << std::endl;
    logmessage = "  trueVtx = (" +to_string(muonstartpos.X()) + ", " + to_string(muonstartpos.Y()) + ", " + to_string(muonstartpos.Z()) +", "+to_string(muonstarttime)+ "\n"
              + "           " +to_string(muondirection.X()) + ", " + to_string(muondirection.Y()) + ", " + to_string(muondirection.Z()) + ") " + "\n";
    
    Log(logmessage,v_debug,verbosity);
  }
  // return vertex
  fExtendedVertex  = (RecoVertex*)(this->FitExtendedVertex(fTrueVertex));
  // Push fitted vertex to RecoEvent store
  this->PushExtendedVertex(fExtendedVertex, true);
  }
  
  else if(fSeedGridFits){
    Log("VtxExtendedVertexFinder Tool: Run vertex reconstruction at all grid seed positions",v_message,verbosity);
  	// Get vertex seed candidates from the store
  	std::vector<RecoVertex>* vSeedVtxList = 0;
  	auto get_seedlist = m_data->Stores.at("RecoEvent")->Get("vSeedVtxList", vSeedVtxList);
  	if(!get_seedlist){ 
      Log("VtxPointPositionFinder Tool: Error retrieving vertex seeds from RecoEvent!",v_error,verbosity);
      Log("VtxPointPositionFinder Tool: Needs to run VtxSeedGenerator first!",v_error,verbosity);
      return false;
    }
    //Now, run FindGridSeeds.
    fExtendedVertex = (RecoVertex*)(this->FitGridSeeds(vSeedVtxList));
    // Push fitted vertex to RecoEvent store
    this->PushExtendedVertex(fExtendedVertex, true);
  }
  
  else {
    Log("VtxExtendedVertexFinder Tool: Run extended vertex reconstruction using point vertex",v_message,verbosity);
    // get point vertex
    RecoVertex* pointvertex = 0; 
    auto get_pointvertex = m_data->Stores.at("RecoEvent")->Get("PointVertex", pointvertex);
    if(!get_pointvertex){ 
      Log("VtxPointVertexFinder Tool: Error retrieving PointVertex from RecoEvent!",v_error,verbosity); 
      return false; 
    }
    // return vertex
    fExtendedVertex  = (RecoVertex*)(this->FitExtendedVertex(pointvertex));
    // Push fitted vertex to RecoEvent store
    this->PushExtendedVertex(fExtendedVertex, true);    	
  }
  return true;
}


bool VtxExtendedVertexFinder::Finalise(){
  // memory has to be freed in the Finalise() function
  delete fExtendedVertex; fExtendedVertex = 0;
  if(verbosity>0) cout<<"VtxExtendedVertexFinder exitting"<<endl;
  return true;
}

RecoVertex* VtxExtendedVertexFinder::FitExtendedVertex(RecoVertex* myVertex) {
  //fit with Minuit
  MinuitOptimizer* myOptimizer = new MinuitOptimizer();
  myOptimizer->SetPrintLevel(-1);
  myOptimizer->SetMeanTimeCalculatorType(1); //Type 1: most probable time
  myOptimizer->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
  myOptimizer->LoadVertex(myVertex); //Load vertex seed
  myOptimizer->SetFitterTimeRange(fTmin, fTmax); //Set time range to fit over 
  myOptimizer->FitExtendedVertexWithMinuit(); //scan the point position in 4D space
  // Fitted vertex must be copied to a new vertex pointer that is created in this class 
  // Once the optimizer is deleted, the fitted vertex is lost. 
  // copy vertex to fExtendedVertex
  RecoVertex* newVertex = new RecoVertex();
  newVertex->CloneVertex(myOptimizer->GetFittedVertex());
  //newVertex->SetFOM(myOptimizer->GetFittedVertex()->GetFOM(),1,1);
  // print vertex
  // ============
  if(verbosity >0) {
  std::cout << "  set extended vertex: " << std::endl
  	        << "  status = "<<newVertex->GetStatus()<<std::endl
            << "     (vx,vy,vz)=(" << newVertex->GetPosition().X() << "," << newVertex->GetPosition().Y() << "," << newVertex->GetPosition().Z() << ") " << std::endl
            << "     vtime=" << newVertex->GetTime() << " itr=" << newVertex->GetIterations() << " fom=" << newVertex->GetFOM() << std::endl;
  }
  delete myOptimizer; myOptimizer = 0;
  return newVertex;
}

RecoVertex* VtxExtendedVertexFinder::FitGridSeeds(std::vector<RecoVertex>* vSeedVtxList) {
  double vtxFOM = -9999.;
  double bestFOM = -1.0;
  int vtxRecoStatus = -1;
  unsigned int nlast = vSeedVtxList->size();
  
  RecoVertex* fSeedPos = 0;
  RecoVertex* fSimpleVertex = new RecoVertex();
  RecoVertex* bestGridVertex = new RecoVertex(); // FIXME: pointer must be deleted by the invoker
  
  for( unsigned int n=0; n<nlast; n++ ){
    //Find best time with Minuit
    MinuitOptimizer* myOptimizer = new MinuitOptimizer();
    myOptimizer->SetPrintLevel(0);
    myOptimizer->SetMeanTimeCalculatorType(1);
    myOptimizer->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
    myOptimizer->SetFitterTimeRange(fTmin, fTmax); //Set time range to fit over 
    fSeedPos = &(vSeedVtxList->at(n));
  	fSimpleVertex= this->FindSimpleDirection(fSeedPos);
    myOptimizer->LoadVertex(fSimpleVertex); //Load vertex seed
    myOptimizer->FitExtendedVertexWithMinuit(); //scan the point position in 4D space
    vtxFOM = myOptimizer->GetFittedVertex()->GetFOM();
    vtxRecoStatus = myOptimizer->GetFittedVertex()->GetStatus();
 
    if((vtxFOM>bestFOM) && (vtxRecoStatus==0)){
      bestGridVertex->CloneVertex(myOptimizer->GetFittedVertex());
      bestFOM = vtxFOM;
    }
    delete myOptimizer; myOptimizer = 0;
  }
  if (verbosity>4){
    std::cout << "Best fit vertex information: " << std::endl;
    std::cout << "bestFOM: " << bestFOM << std::endl;
    std::cout << "best fit reco status: " << bestGridVertex->GetStatus() << std::endl;
    std::cout << "BestVertex info: " << bestGridVertex->Print() << std::endl;
  }
  return bestGridVertex;
}

RecoVertex* VtxExtendedVertexFinder::FindSimpleDirection(RecoVertex* myVertex) {
	
  /// get vertex position
  double vtxX = myVertex->GetPosition().X();
  double vtxY = myVertex->GetPosition().Y();
  double vtxZ = myVertex->GetPosition().Z();
  double vtxTime = myVertex->GetTime();
  
  // current status
  // ==============
  int status = myVertex->GetStatus();
  
  /// loop over digits
  /// ================
  double Swx = 0.0;
  double Swy = 0.0;
  double Swz = 0.0;
  double Sw = 0.0;
  double digitq = 0.;
  double dx, dy, dz, ds, px, py, pz, q;
  
  RecoDigit digit;
  for( int idigit=0; idigit<(int)fDigitList->size(); idigit++ ){
  	digit = fDigitList->at(idigit);
    if( digit.GetFilterStatus() ){ 
      q = digit.GetCalCharge();
      dx = digit.GetPosition().X() - vtxX;
      dy = digit.GetPosition().Y() - vtxY;
      dz = digit.GetPosition().Z() - vtxZ;
      ds = sqrt(dx*dx+dy*dy+dz*dz);
      px = dx/ds;
      py = dx/ds;
      pz = dz/ds;
      Swx += q*px;
      Swy += q*py;
      Swz += q*pz;
      Sw  += q;
    }
  }

  /// average direction
  /// =================
  double dirX = 0.0;
  double dirY = 0.0;
  double dirZ = 0.0;
    
  int itr = 0;
  bool pass = 0; 
  double fom = 0.0;

  if( Sw>0.0 ){
    double qx = Swx/Sw;
    double qy = Swy/Sw;
    double qz = Swz/Sw;
    double qs = sqrt(qx*qx+qy*qy+qz*qz);

    dirX = qx/qs;
    dirY = qy/qs;
    dirZ = qz/qs;

    fom = 1.0;
    itr = 1;
    pass = 1; 
  }  

  // set vertex and direction
  // ========================
  RecoVertex* newVertex = new RecoVertex(); // Note: pointer must be deleted by the invoker
  
  if( pass ){
    newVertex->SetVertex(vtxX,vtxY,vtxZ,vtxTime);
    newVertex->SetDirection(dirX,dirY,dirZ);
    newVertex->SetFOM(fom,itr,pass);
  }

  // set status
  // ==========
  if( !pass ) status |= RecoVertex::kFailSimpleDirection;
  newVertex->SetStatus(status);

  // return vertex
  // =============
  return newVertex;
}

// Add extended vertex to RecoEvent store
void VtxExtendedVertexFinder::PushExtendedVertex(RecoVertex* vtx, bool savetodisk) {  
  // push vertex to RecoEvent store
  Log("VtxExtendedVertexFinder Tool: Push extended vertex to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("ExtendedVertex", fExtendedVertex, savetodisk);
}

void VtxExtendedVertexFinder::Reset() {
	fExtendedVertex->Reset();
}

