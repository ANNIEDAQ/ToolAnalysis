#include "CNNImage.h"

CNNImage::CNNImage():Tool(){}


bool CNNImage::Initialise(std::string configfile, DataModel &data){

  if (verbosity >=2) std::cout <<"Initialising tool: CNNImage..."<<std::endl;

  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();
  m_data= &data; //assigning transient data pointer

  //read in configuration file

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("Mode",mode);
  m_variables.Get("Dimension",dimension);
  m_variables.Get("OutputFile",cnn_outpath);
  m_variables.Get("DetectorConf",detector_config);


  if (mode != "Charge" && mode != "Time") mode = "Charge";
  if (verbosity > 2) std::cout <<"Mode: "<<mode<<std::endl;

  //get geometry		

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  tank_radius = geom->GetTankRadius();
  tank_height = geom->GetTankHalfheight();
  double barrel_compression = 0.82;
  if (detector_config == "ANNIEp2v6") tank_height*=barrel_compression;
  if (tank_radius < 1.) tank_radius = 1.37504;
  Position detector_center = geom->GetTankCentre();
  tank_center_x = detector_center.X();
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();
  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  m_data->CStore.Get("channelkey_to_pmtid",channelkey_to_pmtid);
  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  //read in PMT positions

  max_y = -100.;
  min_y = 100.;

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Tank").begin();
                                                    it != Detectors->at("Tank").end();
                                                  ++it){
    Detector* apmt = it->second;
    unsigned long detkey = it->first;
    pmt_detkeys.push_back(detkey);
    unsigned long chankey = apmt->GetChannels()->begin()->first;
    Position position_PMT = apmt->GetDetectorPosition();
    if (verbosity > 2) std::cout <<"detkey: "<<detkey<<std::endl;
    if (verbosity > 2) std::cout <<"chankey: "<<chankey<<std::endl;
    x_pmt.insert(std::pair<int,double>(detkey,position_PMT.X()-tank_center_x));
    y_pmt.insert(std::pair<int,double>(detkey,position_PMT.Y()-tank_center_y));
    z_pmt.insert(std::pair<int,double>(detkey,position_PMT.Z()-tank_center_z));
    if (verbosity > 2) std::cout <<"Detector ID: "<<detkey<<", position: ("<<position_PMT.X()<<","<<position_PMT.Y()<<","<<position_PMT.Z()<<")"<<std::endl;
    if (verbosity > 2) std::cout <<"rho PMT "<<detkey<<": "<<sqrt(x_pmt.at(detkey)*x_pmt.at(detkey)+z_pmt.at(detkey)*z_pmt.at(detkey))<<std::endl;
    if (verbosity > 2) std::cout <<"y PMT: "<<y_pmt.at(detkey)<<std::endl;
    if (y_pmt[detkey]>max_y) max_y = y_pmt.at(detkey);
    if (y_pmt[detkey]<min_y) min_y = y_pmt.at(detkey);
	  
  }

  //define root and csv files to save histograms (root-files temporarily, for cross-checks)

  std::string str_root = ".root";
  std::string str_csv = ".csv";
  std::string rootfile_name = cnn_outpath + str_root;
  std::string csvfile_name = cnn_outpath + str_csv;

  file = new TFile(rootfile_name.c_str(),"RECREATE");
  outfile.open(csvfile_name.c_str());

  return true;
}


bool CNNImage::Execute(){

  if (verbosity >=2) std::cout <<"Executing tool: CNNImage..."<<std::endl;

  //get ANNIEEvent store information
	
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl;}/*return false;*/

  m_data->Stores["ANNIEEvent"]->Get("MCParticles",mcparticles); //needed to retrieve true vertex and direction
  m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  m_data->Stores["ANNIEEvent"]->Get("SubRunNumber",subrunnumber);
  m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);

  //clear variables & containers
  charge.clear();
  time.clear();
  hitpmt_detkeys.clear();

  for (int i_pmt=0; i_pmt<n_tank_pmts;i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    charge.emplace(detkey,0.);
    time.emplace(detkey,0.);
  }

  //make basic selection cuts to only look at clear event signatures

  bool bool_primary=false;
  bool bool_geometry=false;
  bool bool_nhits=false;

  if (verbosity > 1) std::cout <<"Loop through MCParticles..."<<std::endl;
  for(int particlei=0; particlei<(int)mcparticles->size(); particlei++){
    MCParticle aparticle = mcparticles->at(particlei);
    if (verbosity > 1) std::cout <<"particle "<<particlei<<std::endl;
    if (verbosity > 1) std::cout <<"Parent ID: "<<aparticle.GetParentPdg()<<std::endl;
    if (verbosity > 1) std::cout <<"PDG code: "<<aparticle.GetPdgCode()<<std::endl;
    if (verbosity > 1) std::cout <<"Flag: "<<aparticle.GetFlag()<<std::endl;
    if (aparticle.GetParentPdg() !=0 ) continue;
    if (aparticle.GetFlag() !=0 ) continue;
    if (!(aparticle.GetPdgCode() == 11 || aparticle.GetPdgCode() == 13)) continue;    //primary particle for Cherenkov tracking should be muon or electron
    else {
      truevtx = aparticle.GetStartVertex();
      truevtx_x = truevtx.X()-tank_center_x;
      truevtx_y = truevtx.Y()-tank_center_y;
      truevtx_z = truevtx.Z()-tank_center_z;
      double distInnerStr_Hor = tank_radius - sqrt(pow(truevtx_x,2)+pow(truevtx_z,2));
      double distInnerStr_Vert1 = max_y - truevtx_y;
      double distInnerStr_Vert2 = truevtx_y - min_y;
      double distInnerStr_Vert;
      if (distInnerStr_Vert1 > 0 && distInnerStr_Vert2 > 0) {
	if (distInnerStr_Vert1>distInnerStr_Vert2) distInnerStr_Vert=distInnerStr_Vert2;
	else distInnerStr_Vert=distInnerStr_Vert1;
      } else if (distInnerStr_Vert1 <=0) distInnerStr_Vert=distInnerStr_Vert2;
      else distInnerStr_Vert=distInnerStr_Vert1;
      bool_geometry = (distInnerStr_Vert>0.2 && distInnerStr_Hor>0.2);
      bool_primary = true;
    }
  }

  //---------------------------------------------------------------
  //-------------------Iterate over MCHits ------------------------
  //---------------------------------------------------------------

  int vectsize = MCHits->size();
  if (verbosity > 1) std::cout <<"Tool CNNImage: MCHits size: "<<vectsize<<std::endl; 
  total_hits_pmts=0;
  for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits){
    unsigned long chankey = apair.first;
    if (verbosity > 3) std::cout <<"chankey: "<<chankey<<std::endl;
    Detector* thistube = geom->ChannelToDetector(chankey);
    unsigned long detkey = thistube->GetDetectorID();
    if (verbosity > 3) std::cout <<"detkey: "<<detkey<<std::endl;
    if (thistube->GetDetectorElement()=="Tank"){
      hitpmt_detkeys.push_back(detkey);
      std::vector<MCHit>& Hits = apair.second;
      int hits_pmt = 0;
      int wcsim_id;
      for (MCHit &ahit : Hits){
        charge[detkey] += ahit.GetCharge();
        time[detkey] += ahit.GetTime();
        hits_pmt++;
      }
      time[detkey]/=hits_pmt;         //use mean time of all hits on one PMT
      total_hits_pmts++;
    }
  }

  if (total_hits_pmts>=10) bool_nhits=true;

  //---------------------------------------------------------------
  //------------- Determine max+min values ------------------------
  //---------------------------------------------------------------

  maximum_pmts = 0;
  max_time_pmts = 0;
  min_time_pmts = 999.;
  total_charge_pmts = 0;
  for (int i_pmt=0;i_pmt<(int)hitpmt_detkeys.size();i_pmt++){
    unsigned long detkey = hitpmt_detkeys[i_pmt];
    if (charge[detkey]>maximum_pmts) maximum_pmts = charge[detkey];
    total_charge_pmts+=charge[detkey];
    if (time[detkey]>max_time_pmts) max_time_pmts = time[detkey];
    if (time[detkey]<min_time_pmts) min_time_pmts = time[detkey];
  }

  //---------------------------------------------------------------
  //-------------- Create CNN images ------------------------------
  //---------------------------------------------------------------
	
  //define histogram as an intermediate step to the CNN
  std::stringstream ss_cnn, ss_title_cnn;
  ss_cnn<<"hist_cnn"<<evnum;
  ss_title_cnn<<"EventDisplay (CNN), Event "<<evnum;
  TH2F *hist_cnn = new TH2F(ss_cnn.str().c_str(),ss_title_cnn.str().c_str(),dimension,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimension,0.5+min_y/tank_radius*size_top_drawing, 0.5+max_y/tank_radius*size_top_drawing);

  //fill the events into the histogram

  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    double x,y;
    if ((fabs(y_pmt[detkey]-max_y)<0.01) || fabs(y_pmt[detkey]-min_y)<0.01 || fabs(y_pmt[detkey]+1.30912)<0.01) continue; //top/bottom PMTs
    else {
      double phi;
      if (x_pmt[detkey]>0 && z_pmt[detkey]>0) phi = atan(z_pmt[detkey]/x_pmt[detkey])+TMath::Pi()/2;
      else if (x_pmt[detkey]>0 && z_pmt[detkey]<0) phi = atan(x_pmt[detkey]/-z_pmt[detkey]);
      else if (x_pmt[detkey]<0 && z_pmt[detkey]<0) phi = 3*TMath::Pi()/2+atan(z_pmt[detkey]/x_pmt[detkey]);
      else if (x_pmt[detkey]<0 && z_pmt[detkey]>0) phi = TMath::Pi()+atan(-x_pmt[detkey]/z_pmt[detkey]);
      else phi = 0.;
      if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
      phi-=TMath::Pi();
      if (phi < - TMath::Pi()) phi = -TMath::Pi();
      if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<x_pmt[detkey]<<", y="<<y_pmt[detkey]<<", z="<<z_pmt[detkey]<<std::endl;
      x=0.5+phi*size_top_drawing;
      y=0.5+y_pmt[detkey]/tank_height*tank_height/tank_radius*size_top_drawing;
      int binx = hist_cnn->GetXaxis()->FindBin(x);
      int biny = hist_cnn->GetYaxis()->FindBin(y);
      if (verbosity > 2) std::cout <<"binx: "<<binx<<", biny: "<<biny<<", charge fill: "<<charge[detkey]<<", time fill: "<<time[detkey]<<std::endl;

      if (maximum_pmts < 0.001) maximum_pmts = 1.;
      if (fabs(max_time_pmts) < 0.001) max_time_pmts = 1.;

      if (mode == "Charge"){
        double charge_fill = charge[detkey]/maximum_pmts;
	hist_cnn->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+charge_fill);
      }
      if (mode == "Time"){
	double time_fill = time[detkey]/max_time_pmts;
	hist_cnn->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+time_fill);
      }
    }
  }

  //save information from histogram to csv file
  //(1 line corresponds to 1 event, histogram entries flattened out to a 1D array)

  if (bool_primary && bool_geometry && bool_nhits) {

    hist_cnn->Write();
    for (int i_binY=0; i_binY < hist_cnn->GetNbinsY();i_binY++){
      for (int i_binX=0; i_binX < hist_cnn->GetNbinsX();i_binX++){
        outfile << hist_cnn->GetBinContent(i_binX+1,i_binY+1);
        if (i_binX != hist_cnn->GetNbinsX()-1 || i_binY!=hist_cnn->GetNbinsY()-1) outfile<<",";
      }
    }
    outfile << std::endl;
  }

  return true;
}


bool CNNImage::Finalise(){
  
  if (verbosity >=2 ) std::cout <<"Finalising tool: CNNImage..."<<std::endl;
  file->Close();
  outfile.close();

  return true;
}
