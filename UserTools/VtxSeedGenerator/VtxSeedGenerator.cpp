#include "VtxSeedGenerator.h"

VtxSeedGenerator::VtxSeedGenerator():Tool(){}


bool VtxSeedGenerator::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(verbosity) cout<<"Initializing Tool VtxSeedGenerator"<<endl;
  // Initialise variables

  fNumSeeds = 500;
  fThisDigit = 0;
  fLastEntry = 0; 
  fCounter = 0;
  fSeedType = RecoDigit::All;	 	
  	
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  
  /// Get the Tool configuration variables
	m_variables.Get("SeedType",fSeedType);
	m_variables.Get("NumberOfSeeds", fNumSeeds);
	m_variables.Get("verbosity", verbosity);
	m_variables.Get("UseSeedGrid", UseSeedGrid);
  
  // Make the ANNIEEvent Store if it doesn't exist
	// =============================================
	int recoeventexists = m_data->Stores.count("RecoEvent");
	if(recoeventexists==0){
		cerr<<"VtxSeedGenerator tool need to run after DigitBuilder tool!"<<endl;
		return false; 
	}
		
	// Create an object to store the true MC neutrino vertex
	vSeedVtxList = new std::vector<RecoVertex>;
		
  return true;
}

bool VtxSeedGenerator::Execute(){
  Log("===========================================================================================",v_debug,verbosity);
  
  Log("VtxSeedGenerator Tool: Executing",v_debug,verbosity);
  
  // Reset everything
  this->Reset();
  
  auto get_recoevent = m_data->Stores.count("RecoEvent");
  if(!get_recoevent){
  	Log("VtxSeedGenerator Tool: RecoEvent store doesn't exist!",v_error,verbosity);
  	return false;
  };
  
  
  // check if event passes the cut
  Log("VtxSeedGenerator Tool: Loading EventCutStatus",v_debug,verbosity);
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
  
  Log("VtxSeedGenerator Tool: Loading Digits",v_debug,verbosity);
  // Load digits
  auto get_recodigit = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if (!get_recodigit){ 
  Log("VtxSeedGenerator  Tool: Error retrieving RecoDigits,no digit from the RecoEvent store!",v_error,verbosity); 
  return false;
  }

  // Generate vertex candidates and push to "RecoEvent" store
  if (UseSeedGrid){
    Log("VtxSeedGenerator Tool: Generating seed grid",v_debug,verbosity);
    this->GenerateSeedGrid(fNumSeeds);
  } else {
    Log("VtxSeedGenerator Tool: Generating quadfitter seeds",v_debug,verbosity);
    this->GenerateVertexSeeds(fNumSeeds);
  }
  this->PushVertexSeeds(true);
  Log("VtxSeedGenerator Tool: Execution complete",v_debug,verbosity);
  return true;
}

bool VtxSeedGenerator::Finalise(){
  delete vSeedVtxList; vSeedVtxList = 0;
  Log("VtxSeedGenerator exitting", v_debug,verbosity);
  return true;
}

void VtxSeedGenerator::Reset() {
	vSeedVtxList->clear();
}

bool VtxSeedGenerator::GenerateSeedGrid(int NSeeds) {
  // reset list of seeds
  vSeedVtxList->clear();
  
  // always calculate the simple vertex first
  double VtxX1 = 0.0;
  double VtxY1 = 0.0;
  double VtxZ1 = 0.0;
  double VtxTime1 = 0.0;
  this->CalcSimpleVertex(VtxX1,VtxY1,VtxZ1,VtxTime1);
  
  // add this vertex to seed list
  RecoVertex vtxseed;
  vtxseed.SetVertex(VtxX1,VtxY1,VtxZ1,VtxTime1);
  vSeedVtxList->push_back(vtxseed);
  
  // check limit
  if( NSeeds<=1 ) return false;

  // form list of golden digits used in class methods
  // Here, the digit type (PMT or LAPPD or all) is specified.
  Log("VtxSeedGenerator Tool: Getting clean RecoDigits for event", v_debug,verbosity);
  vSeedDigitList.clear();  
  RecoDigit digit;
  for( fThisDigit=0; fThisDigit<(int)fDigitList->size(); fThisDigit++ ){
  	digit = fDigitList->at(fThisDigit);
    if( digit.GetFilterStatus() ){ 
      if( digit.GetDigitType() == fSeedType){ 
        vSeedDigitList.push_back(fThisDigit);
      }
      else if (fSeedType == 2){
        //Use all digits available
        vSeedDigitList.push_back(fThisDigit);
      }
    }
  }
  
  // Recure at least 5 hits to calculate median times 
  if( vSeedDigitList.size()<5 ) {
    Log("VtxSeedGenerator Tool: The number of digits is less than 5! ",v_message,verbosity);
    return false;
  }
  	
  //Now, we generate our grid of position/time guesses.  Position first.
  //We will use Vogel's method to populate disks with equidistant points
  //inside the ANNIE tank.  The z separation for each disk will be approx. the
  //average separation of points on the disk.
  Log("VtxSeedGenerator Tool: Generating grid of positions", v_debug,verbosity);
  //Here, we need to get the radius and height of the ANNIE tank from the geo
  double pmtradius = ANNIEGeometry::Instance()->GetPMTRadius();	
  double pmtlength = ANNIEGeometry::Instance()->GetCylLength();
  pmtlength = pmtlength * (0.75); //FIXME: Can we add a PMTLength to the geometry? 
  //Assuming roughly equal distance, we can calculate the number of
  //layers needed to have close to even spacing, and # points on each disk
  double approx_points = pow((NSeeds*pmtradius/pmtlength),(2.0/3.0));
  int numlayers = (int) (pmtlength/pow(pow(pmtradius,2)/approx_points,(1.0/2.0)));
  int points_ondisk = (int) approx_points;

  //Now, fill a vector of doubles with the xy plane points 
  double ind, phi, radius, z, x; 
  std::vector<double> xpoints;
  std::vector<double> zpoints;
  double increment = TMath::Pi() * (3.0 - sqrt(5.0));
  for (int i=0; i<points_ondisk; i++) {
    ind = (double) i;
    phi = ind * increment;
    radius = pmtradius * sqrt(ind/approx_points);
    z = radius * cos(phi); //z is the beam axis
    x = radius * sin(phi);
    xpoints.push_back(x);
    zpoints.push_back(z);
  }

  //Now, push a vertex for each point in each disk layer
  double diskind,layers,disk_height,mediantime;
  for (int j=0; j<numlayers; j++){
    for (int k=0; k<points_ondisk; k++) {
      diskind = (double) j;
      layers = (double) numlayers;
      disk_height = pmtlength*((diskind+0.5)/layers) - (pmtlength/2.0);
      Position thisgridpos;
      thisgridpos.SetX(xpoints[k]);
      thisgridpos.SetZ(zpoints[k]);
      thisgridpos.SetY(disk_height);
      mediantime = this->GetMedianSeedTime(thisgridpos);
      RecoVertex thisgridseed;
      thisgridseed.SetVertex(thisgridpos,mediantime);
      vSeedVtxList->push_back(thisgridseed);
    }
  }
  Log("VtxSeedGenerator Tool: Grid of positions and median times calculated", v_debug,verbosity);
  return true;
}

double VtxSeedGenerator::GetMedianSeedTime(Position pos){
  double digitx, digity, digitz, digittime;
  double dx,dy,dz,dr;
  double fC, fN;
  double seedtime;
  std::vector<double> extraptimes;
  for (int entry=0; entry<(int)vSeedDigitList.size(); entry++){
    fThisDigit = vSeedDigitList.at(entry);
    digitx = fDigitList->at(fThisDigit).GetPosition().X();
    digity = fDigitList->at(fThisDigit).GetPosition().Y();
    digitz = fDigitList->at(fThisDigit).GetPosition().Z();
    digittime = fDigitList->at(fThisDigit).GetCalTime();
    //Now, find distance to seed position
    dx = digitx - pos.X();
    dy = digity - pos.Y();
    dz = digitz - pos.Z();
    dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

    //Back calculate to the vertex time using speed of light in H20
    //Very rough estimate; ignores muon path before Cherenkov production
    //TODO: add charge weighting?  Kinda like CalcSimpleVertex?
    fC = Parameters::SpeedOfLight();
    fN = Parameters::Index0();
    seedtime = digittime - (dr/(fC/fN));
    extraptimes.push_back(seedtime);
  }
  //return the median of the extrapolated vertex times
  size_t median_index = extraptimes.size() / 2;
  std::nth_element(extraptimes.begin(), extraptimes.begin()+median_index, extraptimes.end());
  return extraptimes[median_index];
}  

bool VtxSeedGenerator::GenerateVertexSeeds(int NSeeds) {
  double VtxX1 = 0.0;
  double VtxY1 = 0.0;
  double VtxZ1 = 0.0;
  double VtxTime1 = 0.0;
  double VtxX2 = 0.0;
  double VtxY2 = 0.0;
  double VtxZ2 = 0.0;
  double VtxTime2 = 0.0;
	// reset list of seeds
  vSeedVtxList->clear();
  // always calculate the simple vertex first
  this->CalcSimpleVertex(VtxX1,VtxY1,VtxZ1,VtxTime1);
  
  // add this vertex to seed list
  RecoVertex vtxseed;
  vtxseed.SetVertex(VtxX1,VtxY1,VtxZ1,VtxTime1);
  vSeedVtxList->push_back(vtxseed);
  // check limit
  if( NSeeds<=1 ) return false;
  
  // form list of golden digits
  // Digit quality cut can be applied here. 
  // However, the HitCleaner class is a better place to do the cuts. 
  // Here only the digit type (PMT or LAPPD or all) is specified. 
  vSeedDigitList.clear();  
  RecoDigit digit;
  for( fThisDigit=0; fThisDigit<(int)fDigitList->size(); fThisDigit++ ){
  	digit = fDigitList->at(fThisDigit);
    if( digit.GetFilterStatus() ){ 
      if( digit.GetDigitType() == fSeedType){ 
        vSeedDigitList.push_back(fThisDigit);
      }
      else if (fSeedType == 2){
        //Use all digits available
        vSeedDigitList.push_back(fThisDigit);
      }
    }
  }
  // check if the digits are enough to calculate the vertex
  if( vSeedDigitList.size()<=4 ) {
    Log("VtxSeedGenerator Tool: The number of digits is less than 4! ",v_message,verbosity);
    return false;
  }
  	
  // generate a new list of seeds
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double t0 = 0.0;

  double x1 = 0.0;
  double y1 = 0.0;
  double z1 = 0.0;
  double t1 = 0.0; 

  double x2 = 0.0;
  double y2 = 0.0;
  double z2 = 0.0;
  double t2 = 0.0;

  double x3 = 0.0;
  double y3 = 0.0;
  double z3 = 0.0;
  double t3 = 0.0;
  
  unsigned int counter = 0;
  unsigned int NSeedsTarget = NSeeds;
  
  while( vSeedVtxList->size()<NSeedsTarget && counter<100*NSeedsTarget ){
    counter++;

    // choose next four digits
    this->ChooseNextQuadruple(x0,y0,z0,t0,
                              x1,y1,z1,t1,
                              x2,y2,z2,t2,
                              x3,y3,z3,t3);
		
    // find common vertex
    ANNIEGeometry::FindVertex(x0,y0,z0,t0,
                              x1,y1,z1,t1,
                              x2,y2,z2,t2,
                              x3,y3,z3,t3,
                              VtxX1,VtxY1,VtxZ1,VtxTime1,
                              VtxX2,VtxY2,VtxZ2,VtxTime2);

//    std::cout << "   result: (x,y,z,t)=(" << VtxX1 << "," << VtxY1 << "," << VtxZ1 << "," << VtxTime1 << ") " << std::endl
//              << "   result: (x,y,z,t)=(" << VtxX2 << "," << VtxY2 << "," << VtxZ2 << "," << VtxTime2 << ") " << std::endl;
//    std::cout << std::endl;

    if(VtxX1==-99999.9 && VtxX2==-99999.9) continue;
    bool inside_det;

    // add first digit
    if( ANNIEGeometry::Instance()->InsideDetector(VtxX1,VtxY1,VtxZ1) ){
      vtxseed.SetVertex(VtxX1,VtxY1,VtxZ1,VtxTime1);
      vSeedVtxList->push_back(vtxseed);
    }
    
    // add second digit
    if( ANNIEGeometry::Instance()->InsideDetector(VtxX2,VtxY2,VtxZ2) ){
			vtxseed.SetVertex(VtxX2,VtxY2,VtxZ2,VtxTime2);
      vSeedVtxList->push_back(vtxseed);
    }
  }
  this->PushVertexSeeds(true);
  return true;
}

void VtxSeedGenerator::CalcSimpleVertex(double& vtxX, double& vtxY, double& vtxZ, double& vtxTime)
{
  // default vertex
  vtxX = 0.0;
  vtxY = 0.0;
  vtxZ = 0.0;
  vtxTime = 0.0;

  // loop over digits
  double Swx = 0.0;
  double Swy = 0.0;
  double Swz = 0.0;
  double Swt = 0.0;
  double Sw = 0.0;
  double digitq = 0.;
  
  int NDigits = fDigitList->size();
  RecoDigit digit;
  for( int idigit=0; idigit<NDigits; idigit++ ){
  	digit = fDigitList->at(idigit);
  	digitq = digit.GetCalCharge();
    if( digit.GetFilterStatus() ){
      Swx += digitq*digit.GetPosition().X();   
      Swy += digitq*digit.GetPosition().Y();
      Swz += digitq*digit.GetPosition().Z();
      Swt += digitq*digit.GetCalTime();
      Sw  += digitq;
    }
  }
  if( Sw>0.0 ){
    vtxX = Swx/Sw;
    vtxY = Swy/Sw;
    vtxZ = Swz/Sw;
    vtxTime = Swt/Sw;
  }   
  return;
}

void VtxSeedGenerator::ChooseNextQuadruple(double& x0, double& y0, double& z0, double& t0, double& x1, double& y1, double& z1, double& t1, double& x2, double& y2, double& z2, double& t2, double& x3, double& y3, double& z3, double& t3) {
  this->ChooseNextDigit(x0,y0,z0,t0);
  this->ChooseNextDigit(x1,y1,z1,t1);
  this->ChooseNextDigit(x2,y2,z2,t2);
  this->ChooseNextDigit(x3,y3,z3,t3);

  return;
}

void VtxSeedGenerator::ChooseNextDigit(double& xpos, double& ypos, double& zpos, double& time) {
  // default
  xpos=0; ypos=0; zpos=0; time=0;

  // random number generator
  int numEntries = vSeedDigitList.size();

  fCounter++;
  if( fCounter>=(int)fDigitList->size() ) fCounter = 0;
  fThisDigit = vSeedDigitList.at(fLastEntry);

  double r = gRandom->Uniform(); 

  fLastEntry = (int)(r*numEntries);

  // return the new digit
  fThisDigit = vSeedDigitList.at(fLastEntry);
  xpos = fDigitList->at(fThisDigit).GetPosition().X();
  ypos = fDigitList->at(fThisDigit).GetPosition().Y();
  zpos = fDigitList->at(fThisDigit).GetPosition().Z();
  time = fDigitList->at(fThisDigit).GetCalTime();
  
  return;
}


void VtxSeedGenerator::PushVertexSeeds(bool savetodisk) {
  m_data->Stores.at("RecoEvent")->Set("vSeedVtxList", vSeedVtxList, savetodisk); 
}

