#include "VtxPointDirectionFinder.h"

VtxPointDirectionFinder::VtxPointDirectionFinder():Tool(){}


bool VtxPointDirectionFinder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  
  fUseTrueVertexAsSeed = false;
  
  /// Get the Tool configuration variables
	m_variables.Get("UseTrueVertexAsSeed",fUseTrueVertexAsSeed);
	m_variables.Get("verbosity", verbosity);
	
	/// The pointer has to be deleted after usage
	fSimpleDirection = new RecoVertex();
	fPointDirection = new RecoVertex();

  return true;
}

bool VtxPointDirectionFinder::Execute(){
  Log("===========================================================================================",v_debug,verbosity);
	Log("VtxPointDirectionFinder Tool: Executing",v_debug,verbosity);
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
	
	// Retrive digits from RecoEvent
	auto get_recodigit = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(!get_recodigit){
  	Log("VtxPointDirectionFinder  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
  	return false;
  }
	
	// Load digits to VertexGeometry
  VertexGeometry* myvtxgeo = VertexGeometry::Instance();
  myvtxgeo->LoadDigits(fDigitList);
	// Do extended vertex (muon track) reconstruction using MC truth information
	if( fUseTrueVertexAsSeed ){
    Log("VtxPointDirectionFinder Tool: Run direction reconstruction using MC truth information",v_message,verbosity);
    // get truth vertex information 
    auto get_truevertex = m_data->Stores.at("RecoEvent")->Get("TrueVertex", fTrueVertex);
	  if(!get_truevertex){ 
		  Log("VtxPointDirectionFinder Tool: Error retrieving TrueVertex from RecoEvent!",v_error,verbosity); 
		  return false; 
	  }
    // return vertex
    fPointDirection  = (RecoVertex*)(this->FitPointDirection(fTrueVertex));
    // Push fitted vertex to RecoEvent store
    this->PushPointDirection(fPointDirection, true);
  }
  else {
  	Log("VtxPointDirectionFinder Tool: Run direction reconstruction using point position",v_message,verbosity);
    // get point position
    RecoVertex* vtxpos = 0; 
    auto get_pointposition = m_data->Stores.at("RecoEvent")->Get("PointPosition", vtxpos);
	  if(!get_pointposition){ 
		  Log("VtxPointDirectionFinder Tool: Error retrieving PointPosition from RecoEvent!",v_error,verbosity); 
		  return false; 
	  }
	  // Loop over all the digits and find the averaged direction of all photon hits
  	fSimpleDirection = this->FindSimpleDirection(vtxpos);
  	// Push selected point position to RecoEvent store
  	this->PushSimpleDirection(fSimpleDirection, true);
    // return vertex
    fPointDirection  = (RecoVertex*)(this->FitPointDirection(fSimpleDirection));
    // Push fitted vertex to RecoEvent store
    this->PushPointDirection(fPointDirection, true);    	
  }
  return true;
}

bool VtxPointDirectionFinder::Finalise(){
	/// memory has to be freed in the Finalise() function
	/// If the pointer is not delected here, it can be delected in the last tool in the tool chain,
	/// for example the SaveRecoEvent tool
	delete fSimpleDirection; fSimpleDirection = 0;
	delete fPointDirection; fPointDirection = 0;
	if(verbosity>0) cout<<"VtxPointDirectionFinder exitting"<<endl;
  return true;
}

RecoVertex* VtxPointDirectionFinder::FindSimpleDirection(RecoVertex* myVertex) {
	
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
  for( int idigit=0; idigit < (int) fDigitList->size(); idigit++ ){
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

RecoVertex* VtxPointDirectionFinder::FitPointDirection(RecoVertex* myVertex) {
  //fit with Minuit
  MinuitOptimizer* myOptimizer = new MinuitOptimizer();
  myOptimizer->SetPrintLevel(0);
  myOptimizer->SetMeanTimeCalculatorType(1); //Type 1: most probable time
  VertexGeometry* myvtxgeo = VertexGeometry::Instance();
  myOptimizer->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
  myOptimizer->LoadVertex(myVertex); //Load vertex seed
  myOptimizer->FitPointDirectionWithMinuit(); //scan the point position in 4D space
  // Fitted vertex must be copied to a new vertex pointer that is created in this class 
  // Once the optimizer is deleted, the fitted vertex is lost. 
  // copy vertex to fPointDirection
  RecoVertex* newVertex = new RecoVertex();
  newVertex->CloneVertex(myOptimizer->GetFittedVertex());
  // print vertex
  // ============
  if(verbosity >0) {
  std::cout << "  set point direction: " << std::endl
  	        << "  status = "<<newVertex->GetStatus()<<std::endl
            << "    (vx,vy,vz)=(" << newVertex->GetPosition().X() << "," << newVertex->GetPosition().Y() << "," << newVertex->GetPosition().Z() << ") " << std::endl
            << "    (px,py,pz)=(" << newVertex->GetDirection().X() << "," << newVertex->GetDirection().Y() << "," << newVertex->GetDirection().Z() << ") " << std::endl
            << "      angle=" << newVertex->GetConeAngle() << " vtime=" << newVertex->GetTime() << " itr=" << newVertex->GetIterations() << " fom=" << newVertex->GetFOM() << std::endl;
  }
  delete myOptimizer; myOptimizer = 0;
  return newVertex;
}

// Add extended vertex to RecoEvent store
void VtxPointDirectionFinder::PushSimpleDirection(RecoVertex* vtx, bool savetodisk) {  
  // push vertex to RecoEvent store
  Log("VtxPointDirectionFinder Tool: Push simple direction to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("SimpleDirection", fSimpleDirection, savetodisk);
}

// Add extended vertex to RecoEvent store
void VtxPointDirectionFinder::PushPointDirection(RecoVertex* vtx, bool savetodisk) {  
  // push vertex to RecoEvent store
  Log("VtxPointDirectionFinder Tool: Push point direction to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("PointDirection", fPointDirection, savetodisk);
}

void VtxPointDirectionFinder::Reset() {
	fSimpleDirection->Reset();
	fPointDirection->Reset();
}
