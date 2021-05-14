#include "ApplyMRDEff.h"

ApplyMRDEff::ApplyMRDEff():Tool(){}


bool ApplyMRDEff::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  outputfile = "ApplyMRDEff_debug";
  debug_plots = 1;
  mrd_eff_file = "mrd_efficiencies_10-13-2020.txt";
  verbosity = 1;
  use_file_eff = 0;

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("MRDEffFile",mrd_eff_file);
  m_variables.Get("OutputFile",outputfile);
  m_variables.Get("DebugPlots",debug_plots);
  m_variables.Get("MapChankey_WCSimID",file_chankeymap);
  m_variables.Get("UseFileEff",use_file_eff);
  m_variables.Get("FileEff",file_eff);


  // Define histograms for storing general properties of the fitted Mrd tracks
  if (debug_plots){
    std::stringstream rootfilename;
    rootfilename << outputfile << ".root";
    debug_file = new TFile(rootfilename.str().c_str(),"RECREATE");
    mrd_expected = new TH1F("mrd_expected","MRD - Expected hits",310,26,336);
    mrd_observed = new TH1F("mrd_observed","MRD - Observed hits",310,26,336);
    mrd_eff = new TH1F("mrd_eff","MRD - Efficiency",310,26,336);
    tree = new TTree("tree","Tree");
    vec_random = new std::vector<double>;
    dropped_ch = new std::vector<unsigned long>;
    dropped_ch_mc = new std::vector<unsigned long>;
    dropped_ch_time = new std::vector<double>;
    tree->Branch("evnum",&evnum);
    tree->Branch("vec_random",&vec_random);
    tree->Branch("dropped_chkey",&dropped_ch);
    tree->Branch("dropped_chkey_mc",&dropped_ch_mc);
    tree->Branch("dropped_ch_time",&dropped_ch_time);
  }

  // Read in the geometry & get the properties (extents, orientations, etc) of the MRD paddles
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

  //Read in efficiency file
  ifstream eff_file(mrd_eff_file.c_str());
  unsigned long temp_chkey;
  double temp_eff;
  while (!eff_file.eof()){
    eff_file >> temp_chkey >> temp_eff;
    map_eff.emplace(temp_chkey,temp_eff);
    if (eff_file.eof()) break;
  }
  eff_file.close();

  //Read in Channelkey file mapping (data)
  ifstream file_mapping(file_chankeymap);
  unsigned long temp_chankey;
  int temp_wcsimid;
  while (!file_mapping.eof()){
    file_mapping>>temp_chankey>>temp_wcsimid;
    if (file_mapping.eof()) break;
    channelkey_to_mrdpmtid_data.emplace(temp_chankey,temp_wcsimid);
    mrdpmtid_to_channelkey_data.emplace(temp_wcsimid,temp_chankey);
    Log("ApplyMRDEff tool: Emplaced temp_chankey "+std::to_string(temp_chankey)+" with temp_wcsimid "+std::to_string(temp_wcsimid)+"into channelkey_to_mrdpmtid object!",v_debug,verbosity);
  }
  file_mapping.close();

  //Read in Channelkey file mapping (MC)
  m_data->CStore.Get("mrd_tubeid_to_channelkey",mrdpmtid_to_channelkey);
  m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);

  //Initialize objects
  rnd = new TRandom3();
  rnd->SetSeed(0);
  TDCData_mod = new std::map<unsigned long,std::vector<MCHit>>;
  evnum=0;

  //Read-in pre-computed efficiency drop file
  if (use_file_eff){
    File_Eff = new TFile(file_eff.c_str(),"READ");
    Tree_Eff = (TTree*) File_Eff->Get("tree");
    Tree_Eff->SetName("tree_read");
    vec_random_read = new std::vector<double>;
    dropped_ch_read = new std::vector<unsigned long>;
    dropped_ch_mc_read = new std::vector<unsigned long>;
    dropped_ch_time_read = new std::vector<double>;
    Tree_Eff->SetBranchAddress("evnum",&evnum_read);
    Tree_Eff->SetBranchAddress("vec_random",&vec_random_read);
    Tree_Eff->SetBranchAddress("dropped_chkey",&dropped_ch_read);
    Tree_Eff->SetBranchAddress("dropped_chkey_mc",&dropped_ch_mc_read);
    Tree_Eff->SetBranchAddress("dropped_ch_time",&dropped_ch_time_read);   
  }

  return true;
}


bool ApplyMRDEff::Execute(){

  //Get TDCData object from ANNIEStore
  std::cout <<"ApplyMRDEff tool: Data file: Getting TDCData object"<<std::endl;
  bool get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);  // a std::map<unsigned long,vector<MCHit>>

  if(not get_ok){
    Log("TimeClustering Tool: No TDC data in ANNIEEvent!",v_error,verbosity);
    return true;
  }

  std::cout <<"Clear TDCData_mod"<<std::endl;
  //Re-set modified TDCData container
  if (TDCData_mod->size()>0) TDCData_mod->clear();

  vec_random->clear();
  dropped_ch->clear();
  dropped_ch_mc->clear();
  dropped_ch_time->clear();
  if (use_file_eff){
    vec_random_read->clear();
    dropped_ch_read->clear();
    dropped_ch_mc_read->clear();
    dropped_ch_time_read->clear();
    Tree_Eff->GetEntry(evnum);
  }

  std::cout <<"Loop over TDC data hits"<<std::endl;
  //Apply efficiencies to MRD hits
  std::cout <<"TDCData size: "<<TDCData->size()<<std::endl;
  for(auto&& anmrdpmt : (*TDCData)){
    unsigned long chankey = anmrdpmt.first;
    Detector* thedetector = geom->ChannelToDetector(chankey);
    if(thedetector==nullptr){
      Log("TimeClustering Tool: Null detector in TDCData!",v_error,verbosity);
      continue;
    }  
    if(thedetector->GetDetectorElement()!="MRD") continue; // this is a veto hit, not an MRD hit.
    //std::cout <<"chankey: "<<chankey<<std::endl;
    int wcsimid = channelkey_to_mrdpmtid[chankey];
    //std::cout <<"wcsimid: "<<wcsimid<<std::endl;
    unsigned long chankey_data = mrdpmtid_to_channelkey_data[wcsimid-1];
    //std::cout <<"chankey_data: "<<chankey_data<<std::endl;
    double eff = map_eff[chankey_data];
    for(auto&& hitsonthismrdpmt : anmrdpmt.second){
      /*double time = hitsonthismrdpmt.GetTime();
      float charge = hitsonthismrdpmt.GetCharge();
      std::vector<int> parents = hitsonthismrdpmt.GetParents();
      MCHit nexthit(chankey_data, time, charge, parents);*/
      MCHit nexthit = hitsonthismrdpmt;
      double tmp_time = nexthit.GetTime();
      bool omit_hit = false;
      if (use_file_eff){
        if (std::find(dropped_ch_read->begin(),dropped_ch_read->end(),chankey_data)!=dropped_ch_read->end()){
          int index = std::distance(dropped_ch_read->begin(),std::find(dropped_ch_read->begin(),dropped_ch_read->end(),chankey_data));
          double t = dropped_ch_time_read->at(index);
          if (fabs(t-tmp_time)<0.0001) omit_hit = true;
        }
      }
      double tmp = rnd->Uniform();
      vec_random->push_back(tmp);
      if (debug_plots) mrd_expected->Fill(chankey_data);
      if (chankey_data <=52 || chankey_data >=305){ 		//Only apply the efficiency correction to paddles that are not in first or last layer
	if (debug_plots) mrd_observed->Fill(chankey_data);
        if (TDCData_mod->count(chankey)){
          TDCData_mod->at(chankey).push_back(nexthit);
        } else {
          std::vector<MCHit> temp_hits{nexthit};
          TDCData_mod->emplace(chankey,temp_hits);
        }
      } else if (tmp < eff && !use_file_eff){
//      } else if (chankey_data != 77 && chankey_data != 85){
	if (debug_plots) mrd_observed->Fill(chankey_data);
        if (TDCData_mod->count(chankey)){
          TDCData_mod->at(chankey).push_back(nexthit);
        } else {
          std::vector<MCHit> temp_hits{nexthit};
          TDCData_mod->emplace(chankey,temp_hits);
        }
      } else if (use_file_eff && !omit_hit){
	if (debug_plots) mrd_observed->Fill(chankey_data);
        if (TDCData_mod->count(chankey)){
          TDCData_mod->at(chankey).push_back(nexthit);
        } else {
          std::vector<MCHit> temp_hits{nexthit};
          TDCData_mod->emplace(chankey,temp_hits);
        }
      }else {
        dropped_ch->push_back(chankey_data);
        dropped_ch_mc->push_back(chankey);
        dropped_ch_time->push_back(nexthit.GetTime());
      } 
    }
  }

  std::cout <<"Save new TDCData object to ANNIE Store"<<std::endl;
  //Save the new TDCData object in ANNIEStore
  m_data->Stores.at("ANNIEEvent")->Set("TDCData_mod",TDCData_mod);

  tree->Fill();
  evnum++;

  return true;
}


bool ApplyMRDEff::Finalise(){

  //Write histograms to file (if desired), clean up
  if (debug_plots){
    debug_file->cd();
    mrd_expected->Write();
    mrd_observed->Write();
    mrd_eff = (TH1F*) mrd_observed->Clone();
    mrd_eff->Divide(mrd_expected);
    mrd_eff->SetName("mrd_eff");
    mrd_eff->Write();
    tree->Write();
    debug_file->Close();
    delete debug_file;
  }

  return true;
}
