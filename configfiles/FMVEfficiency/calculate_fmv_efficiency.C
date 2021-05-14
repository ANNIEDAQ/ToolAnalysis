void calculate_fmv_efficiency(std::string inputfile, std::string outputfile, bool has_tank){

  TFile *f_out = new TFile(outputfile.c_str(),"RECREATE");

  //Define output histograms
  TH1D *fmv_eff_layer1 = new TH1D("fmv_eff_layer1","FMV Efficiencies Layer 1",13,0,13);
  TH1D *fmv_eff_strict_layer1 = new TH1D("fmv_eff_strict_layer1","FMV Efficiencies Layer 1 (Strict MRD Track)",13,0,13);
  TH1D *fmv_eff_strict_avg_layer1 = new TH1D("fmv_eff_strict_avg_layer1","FMV Average Efficiencies Layer 1 (Strict MRD Track)",13,0,13);
  TH1D *fmv_eff_loose_layer1 = new TH1D("fmv_eff_loose_layer1","FMV Efficiencies Layer 1 (MRD Track)",13,0,13);
  TH1D *fmv_eff_loose_avg_layer1 = new TH1D("fmv_eff_loose_avg_layer1","FMV Average Efficiencies Layer 1 (MRD Track)",13,0,13);
  TH1D *fmv_eff_layer2 = new TH1D("fmv_eff_layer2","FMV Efficiencies Layer 2",13,13,26);
  TH1D *fmv_eff_strict_layer2 = new TH1D("fmv_eff_strict_layer2","FMV Efficiencies Layer 2 (Strict MRD Track)",13,13,26);
  TH1D *fmv_eff_strict_avg_layer2 = new TH1D("fmv_eff_strict_avg_layer2","FMV Average Efficiencies Layer 2 (Strict MRD Track)",13,13,26);
  TH1D *fmv_eff_loose_layer2 = new TH1D("fmv_eff_loose_layer2","FMV Efficiencies Layer 2 (MRD Track)",13,13,26);
  TH1D *fmv_eff_loose_avg_layer2 = new TH1D("fmv_eff_loose_avg_layer2","FMV Average Efficiencies Layer 2 (MRD Track)",13,13,26);
  TH1D *fmv_eff_tank_layer1 = nullptr;
  TH1D *fmv_eff_tank_layer2 = nullptr;
  if (has_tank){
    fmv_eff_tank_layer1 = new TH1D("fmv_eff_tank_layer1","FMV Efficiencies (Tank coincidence) Layer 1",13,0,13);
    fmv_eff_tank_layer2 = new TH1D("fmv_eff_tank_layer2","FMV Efficiencies (Tank coincidence) Layer 2",13,13,26);
  }

  fmv_eff_layer1->SetStats(0);
  fmv_eff_strict_layer1->SetStats(0);
  fmv_eff_strict_avg_layer1->SetStats(0);
  fmv_eff_loose_layer1->SetStats(0);
  fmv_eff_loose_avg_layer1->SetStats(0);
  fmv_eff_layer2->SetStats(0);
  fmv_eff_strict_layer2->SetStats(0);
  fmv_eff_strict_avg_layer2->SetStats(0);
  fmv_eff_loose_layer2->SetStats(0);
  fmv_eff_loose_avg_layer2->SetStats(0);
  fmv_eff_layer1->GetXaxis()->SetTitle("chankey");
  fmv_eff_strict_layer1->GetXaxis()->SetTitle("chankey");
  fmv_eff_strict_avg_layer1->GetXaxis()->SetTitle("chankey");
  fmv_eff_loose_layer1->GetXaxis()->SetTitle("chankey");
  fmv_eff_loose_avg_layer1->GetXaxis()->SetTitle("chankey");
  fmv_eff_layer2->GetXaxis()->SetTitle("chankey");
  fmv_eff_strict_layer2->GetXaxis()->SetTitle("chankey");
  fmv_eff_strict_avg_layer2->GetXaxis()->SetTitle("chankey");
  fmv_eff_loose_layer2->GetXaxis()->SetTitle("chankey");
  fmv_eff_loose_avg_layer2->GetXaxis()->SetTitle("chankey");
  fmv_eff_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_strict_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_strict_avg_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_loose_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_loose_avg_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_strict_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_strict_avg_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_loose_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_loose_avg_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
  fmv_eff_layer1->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_strict_layer1->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_strict_avg_layer1->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_loose_layer1->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_loose_avg_layer1->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_layer2->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_strict_layer2->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_strict_avg_layer2->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_loose_layer2->GetYaxis()->SetRangeUser(0.,1.);
  fmv_eff_loose_avg_layer2->GetYaxis()->SetRangeUser(0.,1.);

  if (has_tank){
    fmv_eff_tank_layer1->SetStats(0);
    fmv_eff_tank_layer2->SetStats(0);
    fmv_eff_tank_layer1->GetXaxis()->SetTitle("chankey");
    fmv_eff_tank_layer2->GetXaxis()->SetTitle("chankey");
    fmv_eff_tank_layer1->GetYaxis()->SetTitle("efficiency #epsilon");
    fmv_eff_tank_layer2->GetYaxis()->SetTitle("efficiency #epsilon");
    fmv_eff_tank_layer1->GetYaxis()->SetRangeUser(0.,1.);
    fmv_eff_tank_layer2->GetYaxis()->SetRangeUser(0.,1.);
  }

  std::vector<TH1D*> eff_hists_strict, eff_hists_loose;

  TFile *f = new TFile(inputfile.c_str(),"READ");

  TH1D *obs_layer1 = (TH1D*) f->Get("fmv_observed_layer1");
  TH1D *exp_layer1 = (TH1D*) f->Get("fmv_expected_layer1");
  TH1D *obs_layer2 = (TH1D*) f->Get("fmv_observed_layer2");
  TH1D *exp_layer2 = (TH1D*) f->Get("fmv_expected_layer2");
  TH1D *obs_track_strict_layer1 = (TH1D*) f->Get("fmv_observed_track_strict_layer1");
  TH1D *exp_track_strict_layer1 = (TH1D*) f->Get("fmv_expected_track_strict_layer1");
  TH1D *obs_track_strict_layer2 = (TH1D*) f->Get("fmv_observed_track_strict_layer2");
  TH1D *exp_track_strict_layer2 = (TH1D*) f->Get("fmv_expected_track_strict_layer2");
  TH1D *obs_track_loose_layer1 = (TH1D*) f->Get("fmv_observed_track_loose_layer1");
  TH1D *exp_track_loose_layer1 = (TH1D*) f->Get("fmv_expected_track_loose_layer1");
  TH1D *obs_track_loose_layer2 = (TH1D*) f->Get("fmv_observed_track_loose_layer2");
  TH1D *exp_track_loose_layer2 = (TH1D*) f->Get("fmv_expected_track_loose_layer2");
 
  TH1D *obs_tank_layer1 = nullptr;
  TH1D *exp_tank_layer1 = nullptr;
  TH1D *obs_tank_layer2 = nullptr;
  TH1D *exp_tank_layer2 = nullptr;
 
  if (has_tank){
    obs_tank_layer1 = (TH1D*) f->Get("fmv_tank_observed_layer1");
    exp_tank_layer1 = (TH1D*) f->Get("fmv_tank_expected_layer1");
    obs_tank_layer2 = (TH1D*) f->Get("fmv_tank_observed_layer2");
    exp_tank_layer2 = (TH1D*) f->Get("fmv_tank_expected_layer2");
  }

  //Layer 1 efficiencies
  for (int i=0; i< 13; i++){
    double obs_counts = obs_layer1->GetBinContent(i+1);
    double exp_counts = exp_layer1->GetBinContent(i+1);
    if (exp_counts == 0) exp_counts = 1;
    fmv_eff_layer1->SetBinContent(i+1,obs_counts/exp_counts);
    double eff_error = sqrt(obs_counts/exp_counts/exp_counts+obs_counts*obs_counts/exp_counts/exp_counts/exp_counts);
    fmv_eff_layer1->SetBinError(i+1,eff_error);
    double obs_counts_strict = obs_track_strict_layer1->GetBinContent(i+1);
    double exp_counts_strict = exp_track_strict_layer1->GetBinContent(i+1);
    if (exp_counts_strict == 0) exp_counts_strict = 1;
    fmv_eff_strict_layer1->SetBinContent(i+1,obs_counts_strict/exp_counts_strict);
    double eff_error_strict = sqrt(obs_counts_strict/exp_counts_strict/exp_counts_strict+obs_counts_strict*obs_counts_strict/exp_counts_strict/exp_counts_strict/exp_counts_strict);
    fmv_eff_strict_layer1->SetBinError(i+1,eff_error_strict);
    double obs_counts_loose = obs_track_loose_layer1->GetBinContent(i+1);
    double exp_counts_loose = exp_track_loose_layer1->GetBinContent(i+1);
    if (exp_counts_loose == 0) exp_counts_loose = 1;
    fmv_eff_loose_layer1->SetBinContent(i+1,obs_counts_loose/exp_counts_loose);
    double eff_error_loose = sqrt(obs_counts_loose/exp_counts_loose/exp_counts_loose+obs_counts_loose*obs_counts_loose/exp_counts_loose/exp_counts_loose/exp_counts_loose);
    fmv_eff_loose_layer1->SetBinError(i+1,eff_error_loose);
    if (has_tank){
      double obs_counts_tank = obs_tank_layer1->GetBinContent(i+1);
      double exp_counts_tank = exp_tank_layer1->GetBinContent(i+1);
      if (exp_counts_tank == 0) exp_counts_tank = 1;
      fmv_eff_tank_layer1->SetBinContent(i+1,obs_counts_tank/exp_counts_tank);
      double eff_error_tank = sqrt(obs_counts_tank/exp_counts_tank/exp_counts_tank+obs_counts_tank*obs_counts_tank/exp_counts_tank/exp_counts_tank/exp_counts_tank);
      fmv_eff_tank_layer1->SetBinError(i+1,eff_error_tank);
    }
  }
  
  //Layer 2 efficiencies
  for (int i=0; i< 13; i++){
    double obs_counts = obs_layer2->GetBinContent(i+1);
    double exp_counts = exp_layer2->GetBinContent(i+1);
    if (exp_counts == 0) exp_counts = 1;
    fmv_eff_layer2->SetBinContent(i+1,obs_counts/exp_counts);
    double eff_error = sqrt(obs_counts/exp_counts/exp_counts+obs_counts*obs_counts/exp_counts/exp_counts/exp_counts);
    fmv_eff_layer2->SetBinError(i+1,eff_error);
    double obs_counts_strict = obs_track_strict_layer2->GetBinContent(i+1);
    double exp_counts_strict = exp_track_strict_layer2->GetBinContent(i+1);
    if (exp_counts_strict == 0) exp_counts_strict = 1;
    fmv_eff_strict_layer2->SetBinContent(i+1,obs_counts_strict/exp_counts_strict);
    double eff_error_strict = sqrt(obs_counts_strict/exp_counts_strict/exp_counts_strict+obs_counts_strict*obs_counts_strict/exp_counts_strict/exp_counts_strict/exp_counts_strict);
    fmv_eff_strict_layer2->SetBinError(i+1,eff_error_strict);
    double obs_counts_loose = obs_track_loose_layer2->GetBinContent(i+1);
    double exp_counts_loose = exp_track_loose_layer2->GetBinContent(i+1);
    if (exp_counts_loose == 0) exp_counts_loose = 1;
    fmv_eff_loose_layer2->SetBinContent(i+1,obs_counts_loose/exp_counts_loose);
    double eff_error_loose = sqrt(obs_counts_loose/exp_counts_loose/exp_counts_loose+obs_counts_loose*obs_counts_loose/exp_counts_loose/exp_counts_loose/exp_counts_loose);
    fmv_eff_loose_layer2->SetBinError(i+1,eff_error_loose);
    if (has_tank){
      double obs_counts_tank = obs_tank_layer2->GetBinContent(i+1);
      double exp_counts_tank = exp_tank_layer2->GetBinContent(i+1);
      if (exp_counts_tank == 0) exp_counts_tank = 1;
      fmv_eff_tank_layer2->SetBinContent(i+1,obs_counts_tank/exp_counts_tank);
      double eff_error_tank = sqrt(obs_counts_tank/exp_counts_tank/exp_counts_tank+obs_counts_tank*obs_counts_tank/exp_counts_tank/exp_counts_tank/exp_counts_tank);
      fmv_eff_tank_layer2->SetBinError(i+1,eff_error_tank);
    }
  }

  //Channel-per-channel efficiencies
  for (int i=0; i< 26; i++){

    //Get observed/expected histograms
    std::stringstream hist_obs_strict_name, hist_exp_strict_name, hist_obs_loose_name, hist_exp_loose_name, hist_eff_strict_name, hist_eff_strict_title, hist_eff_loose_name, hist_eff_loose_title;
    hist_obs_strict_name <<"hist_observed_strict_chankey"<<i;
    hist_exp_strict_name <<"hist_expected_strict_chankey"<<i;
    hist_obs_loose_name <<"hist_observed_loose_chankey"<<i;
    hist_exp_loose_name <<"hist_expected_loose_chankey"<<i;
    hist_eff_strict_name <<"hist_efficiency_strict_chankey"<<i;
    hist_eff_strict_title << "FMV Efficiency (Strict MRD track) Chankey "<<i;
    hist_eff_loose_name << "hist_efficiency_loose_chankey"<<i;
    hist_eff_loose_title << "FMV Efficiency (MRD track) Chankey "<<i;
    TH1F* hist_obs_strict_temp = (TH1F*) f->Get(hist_obs_strict_name.str().c_str());
    TH1F* hist_exp_strict_temp = (TH1F*) f->Get(hist_exp_strict_name.str().c_str());
    TH1F* hist_obs_loose_temp = (TH1F*) f->Get(hist_obs_loose_name.str().c_str());
    TH1F* hist_exp_loose_temp = (TH1F*) f->Get(hist_exp_loose_name.str().c_str());

    //Create efficiency histograms
    f_out->cd();
    TH1D* hist_eff_strict_temp = new TH1D(hist_eff_strict_name.str().c_str(),hist_eff_strict_title.str().c_str(),hist_obs_strict_temp->GetXaxis()->GetNbins()/2,hist_obs_strict_temp->GetXaxis()->GetXmin(),hist_obs_strict_temp->GetXaxis()->GetXmax());
    hist_eff_strict_temp->GetXaxis()->SetTitle("x [m]");
    hist_eff_strict_temp->GetYaxis()->SetTitle("efficiency #epsilon");
    hist_eff_strict_temp->GetYaxis()->SetRangeUser(0.,1.);
    hist_eff_strict_temp->SetStats(0);
    TH1D* hist_eff_loose_temp = new TH1D(hist_eff_loose_name.str().c_str(),hist_eff_loose_title.str().c_str(),hist_obs_loose_temp->GetXaxis()->GetNbins()/2,hist_obs_loose_temp->GetXaxis()->GetXmin(),hist_obs_loose_temp->GetXaxis()->GetXmax());
    hist_eff_loose_temp->GetXaxis()->SetTitle("x [m]");
    hist_eff_loose_temp->GetYaxis()->SetTitle("efficiency #epsilon");
    hist_eff_loose_temp->GetYaxis()->SetRangeUser(0.,1.);
    hist_eff_loose_temp->SetStats(0);
   
    TH1F *hist_obs_strict_rebin = (TH1F*) hist_obs_strict_temp->Rebin(2);
    TH1F *hist_exp_strict_rebin = (TH1F*) hist_exp_strict_temp->Rebin(2);
    TH1F *hist_obs_loose_rebin = (TH1F*) hist_obs_loose_temp->Rebin(2);
    TH1F *hist_exp_loose_rebin = (TH1F*) hist_exp_loose_temp->Rebin(2);

    double avg_efficiency_strict=0;
    double error_avg_efficiency_strict=0;
    double avg_efficiency_loose=0;
    double error_avg_efficiency_loose=0;

    int non_zero_bins_strict=0;
    int non_zero_bins_loose=0;
    for (int bin=0; bin<hist_obs_strict_rebin->GetXaxis()->GetNbins(); bin++){
      double obs_counts = hist_obs_strict_rebin->GetBinContent(bin+1);
      double exp_counts = hist_exp_strict_rebin->GetBinContent(bin+1);
      if (exp_counts == 0) exp_counts = 1;
      else non_zero_bins_strict++;
      hist_eff_strict_temp->SetBinContent(bin+1,obs_counts/exp_counts);
      double eff_error = sqrt(obs_counts/exp_counts/exp_counts+obs_counts*obs_counts/exp_counts/exp_counts/exp_counts);
      hist_eff_strict_temp->SetBinError(bin+1,eff_error);
      double obs_counts_loose = hist_obs_loose_rebin->GetBinContent(bin+1);
      double exp_counts_loose = hist_exp_loose_rebin->GetBinContent(bin+1);
      if (exp_counts_loose == 0) exp_counts_loose = 1;
      else non_zero_bins_loose++;
      hist_eff_loose_temp->SetBinContent(bin+1,obs_counts_loose/exp_counts_loose);
      double eff_loose_error = sqrt(obs_counts_loose/exp_counts_loose/exp_counts_loose+obs_counts_loose*obs_counts_loose/exp_counts_loose/exp_counts_loose/exp_counts_loose);
      hist_eff_loose_temp->SetBinError(bin+1,eff_loose_error);
      if (exp_counts >0){
        avg_efficiency_strict+=(obs_counts/exp_counts);
        error_avg_efficiency_strict+=(eff_error*eff_error);
      }
      if (exp_counts_loose > 0){
        avg_efficiency_loose+=(obs_counts_loose/exp_counts_loose);
        error_avg_efficiency_loose+=(eff_loose_error*eff_loose_error);
      }
    }
    eff_hists_strict.push_back(hist_eff_strict_temp);
    eff_hists_loose.push_back(hist_eff_loose_temp);

    //avg_efficiency_strict/=hist_obs_strict_rebin->GetXaxis()->GetNbins();
    //avg_efficiency_loose/=hist_obs_strict_rebin->GetXaxis()->GetNbins();
    avg_efficiency_strict/=non_zero_bins_strict;
    avg_efficiency_loose/=non_zero_bins_loose;
    error_avg_efficiency_strict=sqrt(error_avg_efficiency_strict);
    //error_avg_efficiency_strict/=hist_obs_strict_rebin->GetXaxis()->GetNbins();
    error_avg_efficiency_strict/=non_zero_bins_strict;
    error_avg_efficiency_loose=sqrt(error_avg_efficiency_loose);
    //error_avg_efficiency_loose/=hist_obs_strict_rebin->GetXaxis()->GetNbins();
    error_avg_efficiency_loose/=non_zero_bins_loose;
    if (i < 13) {
      fmv_eff_strict_avg_layer1->SetBinContent(i+1,avg_efficiency_strict);
      fmv_eff_strict_avg_layer1->SetBinError(i+1,error_avg_efficiency_strict);
      fmv_eff_loose_avg_layer1->SetBinContent(i+1,avg_efficiency_loose);
      fmv_eff_loose_avg_layer1->SetBinError(i+1,error_avg_efficiency_loose);
    }
    else {
      fmv_eff_strict_avg_layer2->SetBinContent(i-13+1,avg_efficiency_strict);
      fmv_eff_strict_avg_layer2->SetBinError(i-13+1,error_avg_efficiency_strict);
      fmv_eff_loose_avg_layer2->SetBinContent(i-13+1,avg_efficiency_loose);
      fmv_eff_loose_avg_layer2->SetBinError(i-13+1,error_avg_efficiency_loose);
    }
  }

  f_out->cd();
  fmv_eff_layer1->Write();
  fmv_eff_strict_layer1->Write();
  fmv_eff_strict_avg_layer1->Write();
  fmv_eff_loose_layer1->Write();
  fmv_eff_loose_avg_layer1->Write();
  fmv_eff_layer2->Write();
  fmv_eff_strict_layer2->Write();
  fmv_eff_strict_avg_layer2->Write();
  fmv_eff_loose_layer2->Write();
  fmv_eff_loose_avg_layer2->Write();
  if (has_tank){
    fmv_eff_tank_layer1->Write();
    fmv_eff_tank_layer2->Write();
  }

  for (int i=0; i< int(eff_hists_strict.size()); i++){
    eff_hists_strict.at(i)->Write();
    eff_hists_loose.at(i)->Write();
  }

  f_out->Close();
  f->Close();

  delete f_out;
  delete f;

}
