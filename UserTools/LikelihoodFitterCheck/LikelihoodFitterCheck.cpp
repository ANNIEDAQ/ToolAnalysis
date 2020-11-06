#include "LikelihoodFitterCheck.h"
#include "TVector3.h"

LikelihoodFitterCheck::LikelihoodFitterCheck():Tool(){}


bool LikelihoodFitterCheck::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();
  std::string output_filename;
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFile", output_filename);
  m_variables.Get("ifPlot2DFOM", ifPlot2DFOM);
  m_variables.Get("ShowEvent", fShowEvent);
  fOutput_tfile = new TFile(output_filename.c_str(), "recreate");
  
  // Histograms
  Likelihood2D = new TH2D("Likelihood2D","Figure of merit 2D", 200, -50, 150, 100, -50, 50);
  gr_parallel = new TGraph();
  gr_parallel->SetTitle("Figure of merit parallel to the track direction");
	gr_transverse = new TGraph();
  gr_transverse->SetTitle("Figure of merit transverse to the track direction");
  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool LikelihoodFitterCheck::Execute(){
	Log("===========================================================================================",v_debug,verbosity);
	
	Log("LikelihoodFitterCheck Tool: Executing",v_debug,verbosity);
	
	// Get a pointer to the ANNIEEvent Store
  auto* annie_event = m_data->Stores["RecoEvent"];
  if (!annie_event) {
    Log("Error: The PhaseITreeMaker tool could not find the ANNIEEvent Store",
      0, verbosity);
    return false;
  }
  
  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",fMCEventNum);  
  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggerNum); 
  
  // ANNIE Event number
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
  
  // Only check this event
  if(fShowEvent>0 && (int)fEventNumber!=fShowEvent) return true; 
  	
  logmessage = "Likelihood check for MC Entry Number " + to_string(fMCEventNum) 
               + " , MC Trigger Number" + to_string(fMCTriggerNum) 
               + " and ANNIIE Event Number " + to_string(fEventNumber);
	Log(logmessage,v_message,verbosity);
  
  // Read True Vertex   
  RecoVertex* truevtx = 0;
  auto get_vtx = m_data->Stores.at("RecoEvent")->Get("TrueVertex",fTrueVertex);  ///> Get digits from "RecoEvent" 
  if(!get_vtx){ 
  	Log("LikelihoodFitterCheck  Tool: Error retrieving TrueVertex! ",v_error,verbosity); 
  	return false;
  }
	
	// Retrive digits from RecoEvent
	auto get_digit = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(!get_digit){
  	Log("LikelihoodFitterCheck  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
  	return false;
  }
	
	double recoVtxX, recoVtxY, recoVtxZ, recoVtxT, recoDirX, recoDirY, recoDirZ;
  double trueVtxX, trueVtxY, trueVtxZ, trueVtxT, trueDirX, trueDirY, trueDirZ;
  double seedX, seedY, seedZ, seedT, seedDirX, seedDirY, seedDirZ;
  double ConeAngle = Parameters::CherenkovAngle();

  // Get true Vertex information
  Position vtxPos = fTrueVertex->GetPosition();
	Direction vtxDir = fTrueVertex->GetDirection();
	trueVtxX = vtxPos.X();
  trueVtxY = vtxPos.Y();
  trueVtxZ = vtxPos.Z();
  trueVtxT = fTrueVertex->GetTime();
  trueDirX = vtxDir.X();
  trueDirY = vtxDir.Y();
  trueDirZ = vtxDir.Z();
  
  if(verbosity>0) cout<<"True vertex  = ("<<trueVtxX<<", "<<trueVtxY<<", "<<trueVtxZ<<", "<<trueVtxT<<", "<<trueDirX<<", "<<trueDirY<<", "<<trueDirZ<<")"<<endl;
  
  FoMCalculator * myFoMCalculator = new FoMCalculator();
  VertexGeometry* myvtxgeo = VertexGeometry::Instance();
  myvtxgeo->LoadDigits(fDigitList);
  myFoMCalculator->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
  //parallel direction
  double dl = 1.0; // step size  = 1 cm along the track
  double dx = dl * trueDirX;
  double dy = dl * trueDirY;
  double dz = dl * trueDirZ;
  int nbins = 200;
  double dlpara[200], dlfom[200];
  for(int j=0;j<200;j++) {
    seedX = trueVtxX - 50*dx + j*dx;
    seedY = trueVtxY - 50*dy + j*dy;
    seedZ = trueVtxZ - 50*dz + j*dz;
    seedT = trueVtxT;
    seedDirX = trueDirX;
    seedDirY = trueDirY;
    seedDirZ = trueDirZ;
    myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, 0.0, seedDirX, seedDirY, seedDirZ);
    int nhits = myvtxgeo->GetNDigits();
    double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
    Double_t fom = -999.999*100;
    double timefom = -999.999*100;
    double conefom = -999.999*100;
    myFoMCalculator->TimePropertiesLnL(meantime,timefom);
    myFoMCalculator->ConePropertiesFoM(ConeAngle,conefom);
    fom = timefom*0.5+conefom*0.5;
    cout<<"timeFOM, coneFOM, fom = "<<timefom<<", "<<conefom<<", "<<fom<<endl;
    fom = timefom;
    dlpara[j] = - 50*dl + j*dl;
    dlfom[j] = fom;
    gr_parallel->SetPoint(j, dlpara[j], dlfom[j]);
  } 
  
  //transverse direction
  double dltrans[200];
  //first find the projection vector of x axis on the plane perpendicular to the track direction
  TVector3 n(trueDirX, trueDirY, trueDirZ);
  TVector3 v = n.Orthogonal();
  dl = 1.0; // step size  = 1 cm along the track
  dx = dl * v.X();
  dy = dl * v.Y();
  dz = dl * v.Z();
  for(int j=0;j<100;j++) {
    seedX = trueVtxX - 50*dx + j*dx;
    seedY = trueVtxY - 50*dy + j*dy;
    seedZ = trueVtxZ - 50*dz + j*dz;
    seedT = trueVtxT;
    seedDirX = trueDirX;
    seedDirY = trueDirY;
    seedDirZ = trueDirZ;
    myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, 0.0, seedDirX, seedDirY, seedDirZ);
    int nhits = myvtxgeo->GetNDigits();
    double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
    Double_t fom = -999.999*100;
    double timefom = -999.999*100;
    double conefom = -999.999*100;
    double coneAngle = 42.0;
    myFoMCalculator->TimePropertiesLnL(meantime,timefom);
    myFoMCalculator->ConePropertiesFoM(ConeAngle,conefom);
    fom = timefom*0.5+conefom*0.5;
    //fom = timefom;
    cout<<"timeFOM, coneFOM, fom = "<<timefom<<", "<<conefom<<", "<<fom<<endl;
    dltrans[j] = - 50*dl + j*dl;
    dlfom[j] = fom;
    gr_transverse->SetPoint(j, dlpara[j], dlfom[j]);
  }
  
    if(ifPlot2DFOM) {
      //2D scan around the true vertex position
      double dl_para = 1.0, dl_trans = 1.0;
      double dx_para = dl_para * trueDirX;
      double dy_para = dl_para * trueDirY;
      double dz_para = dl_para * trueDirZ;
      double dx_trans = dl_trans * v.X();
      double dy_trans = dl_trans * v.Y();
      double dz_trans = dl_trans * v.Z();
      for(int k=0; k<100; k++) {
        for(int m=0; m<200; m++) {
        	seedX = trueVtxX - 50*dx_trans + k*dx_trans - 50*dx_para + m*dx_para;
        	seedY = trueVtxY - 50*dy_trans + k*dy_trans - 50*dy_para + m*dy_para;
        	seedZ = trueVtxZ - 50*dz_trans + k*dz_trans - 50*dz_para + m*dz_para;
        	seedT = trueVtxT;
        	seedDirX = trueDirX;
        	seedDirY = trueDirY;
        	seedDirZ = trueDirZ;
        	myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, seedT, seedDirX, seedDirY, seedDirZ);
        	int nhits = myvtxgeo->GetNDigits();
          double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
          Double_t fom = -999.999*100;
          double timefom = -999.999*100;
          double conefom = -999.999*100;
          double coneAngle = 42.0;
          myFoMCalculator->TimePropertiesLnL(meantime,timefom);
          myFoMCalculator->ConePropertiesFoM(coneAngle,conefom);
          fom = timefom*0.5+conefom*0.5;
          //fom = timefom;
          cout<<"k,m, timeFOM, coneFOM, fom = "<<k<<", "<<m<<", "<<timefom<<", "<<conefom<<", "<<fom<<endl;
          Likelihood2D->SetBinContent(m, k, fom);
        }
      }
    }
  delete myFoMCalculator;
  return true;
}


bool LikelihoodFitterCheck::Finalise(){
  fOutput_tfile->cd();
  gr_parallel->Write();
  gr_transverse->Write();
  fOutput_tfile->Write();
  fOutput_tfile->Close();
  Log("LikelihoodFitterCheck exitting", v_debug,verbosity);
  return true;
}
