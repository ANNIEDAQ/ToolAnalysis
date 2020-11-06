#include "TankCalibrationDiffuser.h"

TankCalibrationDiffuser::TankCalibrationDiffuser():Tool(){}


bool TankCalibrationDiffuser::Initialise(std::string configfile, DataModel &data){

  if(verbose > 0) cout <<"Initialising Tool TankCalibrationDiffuser ..."<<endl;

  //only for debugging memory leaks, otherwise comment out
  /*std::cout <<"List of Objects (beginning of initialise): "<<std::endl;
  gObjectTable->Print();*/

  //---------------------------------------------------------------------------
  //---------------Useful header-----------------------------------------------
  //---------------------------------------------------------------------------

  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  m_data= &data; //assigning transient data pointer

  outputfile = "LEDCalibration";
  outputdir = "./LEDCalibrationFiles/";   //make sure this directory exists! (Or set a different one in the config file)

  //----------------------------------------------------------------------------
  //---------------Get configuration variables for this tool--------------------
  //----------------------------------------------------------------------------

  m_variables.Get("HitStore",HitStoreName);
  m_variables.Get("OutputFile",outputfile);
  m_variables.Get("OutputDirectory",outputdir);
  m_variables.Get("NBinsTimeTotal",nBinsTimeTotal);
  m_variables.Get("TimeTotalMin",timeTotalMin);
  m_variables.Get("TimeTotalMax",timeTotalMax);
  m_variables.Get("NBinsChargeTotal",nBinsChargeTotal);
  m_variables.Get("ChargeTotalMin",chargeTotalMin);
  m_variables.Get("ChargeTotalMax",chargeTotalMax);
  m_variables.Get("NBinsTime",nBinsTime);
  m_variables.Get("TimeMin",timeMin);
  m_variables.Get("TimeMax",timeMax);
  m_variables.Get("NBinsCharge",nBinsCharge);
  m_variables.Get("ChargeMin",chargeMin);
  m_variables.Get("ChargeMax",chargeMax);
  m_variables.Get("NBinsStartTimeTotal",nBinsStartTimeTotal);
  m_variables.Get("NBinsStartTime",nBinsStartTime);
  m_variables.Get("StartTimeMin",startTimeMin);
  m_variables.Get("StartTimeMax",startTimeMax);
  m_variables.Get("NBinsPeakTimeTotal",nBinsPeakTimeTotal);
  m_variables.Get("NBinsPeakTime",nBinsPeakTime);
  m_variables.Get("PeakTimeMin",peakTimeMin);
  m_variables.Get("PeakTimeMax",peakTimeMax);
  m_variables.Get("NBinsBaselineTotal",nBinsBaselineTotal);
  m_variables.Get("NBinsBaseline",nBinsBaseline);
  m_variables.Get("BaselineMin",baselineMin);
  m_variables.Get("BaselineMax",baselineMax);
  m_variables.Get("NBinsSigmaBaselineTotal",nBinsSigmaBaselineTotal);
  m_variables.Get("NBinsSigmaBaseline",nBinsSigmaBaseline);
  m_variables.Get("SigmaBaselineMin",sigmaBaselineMin);
  m_variables.Get("SigmaBaselineMax",sigmaBaselineMax); 
  m_variables.Get("NBinsRawAmplitudeTotal",nBinsRawAmplitudeTotal);
  m_variables.Get("NBinsRawAmplitude",nBinsRawAmplitude);
  m_variables.Get("RawAmplitudeMin",rawAmplitudeMin);
  m_variables.Get("RawAmplitudeMax",rawAmplitudeMax);
  m_variables.Get("NBinsAmplitudeTotal",nBinsAmplitudeTotal);
  m_variables.Get("NBinsAmplitude",nBinsAmplitude);
  m_variables.Get("AmplitudeMin",amplitudeMin);
  m_variables.Get("AmplitudeMax",amplitudeMax);
  m_variables.Get("NBinsRawAreaTotal",nBinsRawAreaTotal);
  m_variables.Get("NBinsRawArea",nBinsRawArea);
  m_variables.Get("RawAreaMin",rawAreaMin);
  m_variables.Get("RawAreaMax",rawAreaMax);
  m_variables.Get("NBinsTimeFit",nBinsTimeFit);
  m_variables.Get("TimeFitMin",timeFitMin);
  m_variables.Get("TimeFitMax",timeFitMax);
  m_variables.Get("NBinsTimeDev",nBinsTimeDev);
  m_variables.Get("TimeDevMin",timeDevMin);
  m_variables.Get("TimeDevMax",timeDevMax);
  m_variables.Get("NBinsChargeFit",nBinsChargeFit);
  m_variables.Get("ChargeFitMin",chargeFitMin);
  m_variables.Get("ChargeFitMax",chargeFitMax);
  m_variables.Get("DiffuserX",diffuser_x);
  m_variables.Get("DiffuserY",diffuser_y);
  m_variables.Get("DiffuserZ",diffuser_z);
  m_variables.Get("ToleranceCharge",tolerance_charge);
  m_variables.Get("ToleranceTime",tolerance_time);
  m_variables.Get("FitMethod",FitMethod);
  m_variables.Get("Gaus1Constant",gaus1Constant);
  m_variables.Get("Gaus1Mean",gaus1Mean);
  m_variables.Get("Gaus1Sigma",gaus1Sigma);
  m_variables.Get("Gaus2Constant",gaus2Constant);
  m_variables.Get("Gaus2Mean",gaus2Mean);
  m_variables.Get("Gaus2Sigma",gaus2Sigma);
  m_variables.Get("ExpConstant",expConstant);
  m_variables.Get("ExpDecay",expDecay);
  m_variables.Get("TApplication",use_tapplication);
  m_variables.Get("verbose",verbose);

  help_file = new TFile("configfiles/TankCalibrationDiffuser/help_file.root","RECREATE");      //let the histograms be associated with this file

  problematic_channels.push_back(333);
  problematic_channels.push_back(346);
  problematic_channels.push_back(349);
  problematic_channels.push_back(352);
//  problematic_channels.push_back(357);
//  problematic_channels.push_back(360);
//  problematic_channels.push_back(404);

  //----------------------------------------------------------------------------
  //---------------Get basic geometry properties -------------------------------
  //----------------------------------------------------------------------------

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  n_lappds = geom->GetNumDetectorsInSet("LAPPD");
  n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
  n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
  Position detector_center=geom->GetTankCentre();
  tank_center_x = detector_center.X();
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();

  std::map<std::string, std::map<unsigned long,Detector*>>* Detectors = geom->GetDetectors();
  std::cout <<"Detectors size: "<<Detectors->size()<<std::endl;  // num detector sets...


  //----------------------------------------------------------------------------
  //---------------Read in geometry properties of PMTs--------------------------
  //----------------------------------------------------------------------------

  //FIXME: put PMT radii in geometry class, read directly and not via PMT type mapping
  map_type_radius.emplace("R5912HQE",0.1016);
  map_type_radius.emplace("R7081",0.127);
  map_type_radius.emplace("R7081HQE",0.127);
  map_type_radius.emplace("D784KFLB",0.1397);
  map_type_radius.emplace("EMI9954KB",0.0508);


  //----------------------------------------------------------------------------
  //---------------calculate expected hit times for PMTs------------------------
  //----------------------------------------------------------------------------

  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Tank").begin(); it != Detectors->at("Tank").end(); ++it){

    Detector* apmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = apmt->GetChannels()->begin()->first;
    pmt_detkeys.push_back(detkey);    
    std::string dettype = apmt->GetDetectorType();

    if (verbose > 1) std::cout <<"PMT with detkey: "<<detkey<<", has type "<<dettype<<std::endl;
    if (verbose > 1) std::cout <<"PMT Print:" <<apmt->Print()<<std::endl;

    Position position_PMT = apmt->GetDetectorPosition();
    if (verbose > 2) std::cout <<"detkey: "<<detkey<<", chankey: "<<chankey<<std::endl;
    if (verbose > 2) std::cout <<"filling PMT position maps"<<std::endl;

    x_PMT.insert(std::pair<unsigned long,double>(detkey,position_PMT.X()-tank_center_x));
    y_PMT.insert(std::pair<unsigned long,double>(detkey,position_PMT.Y()-tank_center_y));
    z_PMT.insert(std::pair<unsigned long,double>(detkey,position_PMT.Z()-tank_center_z));

    double rho = sqrt(pow(x_PMT.at(detkey),2)+pow(y_PMT.at(detkey),2));
    if (x_PMT.at(detkey) < 0) rho*=-1;
    rho_PMT.insert(std::pair<unsigned long, double>(detkey,rho));
    double phi;
    if (x_PMT.at(detkey)>0 && z_PMT.at(detkey)>0) phi = atan(x_PMT.at(detkey)/z_PMT.at(detkey));
    if (x_PMT.at(detkey)>0 && z_PMT.at(detkey)<0) phi = TMath::Pi()/2+atan(-z_PMT.at(detkey)/x_PMT.at(detkey));
    if (x_PMT.at(detkey)<0 && z_PMT.at(detkey)<0) phi = TMath::Pi()+atan(x_PMT.at(detkey)/z_PMT.at(detkey));
    if (x_PMT.at(detkey)<0 && z_PMT.at(detkey)>0) phi = 3*TMath::Pi()/2+atan(z_PMT.at(detkey)/-x_PMT.at(detkey));
    phi_PMT.insert(std::pair<unsigned long,double>(detkey,phi));

    if (verbose > 2) std::cout <<"detectorkey: "<<detkey<<", position: ("<<position_PMT.X()<<","<<position_PMT.Y()<<","<<position_PMT.Z()<<")"<<std::endl;
    if (verbose > 2) std::cout <<"rho PMT "<<detkey<<": "<<rho<<std::endl;
    if (verbose > 2) std::cout <<"y PMT: "<<y_PMT.at(detkey)<<std::endl;
    if (verbose > 2) std::cout <<"phi PMT: "<<phi_PMT.at(detkey)<<std::endl;

    PMT_ishit.insert(std::pair<unsigned long, int>(detkey,0));

    if (dettype.find("R5912") != std::string::npos) radius_PMT[detkey] = map_type_radius["R5912HQE"];
    else if (dettype.find("R7081") != std::string::npos) radius_PMT[detkey] = map_type_radius["R7081"];
    else if (dettype.find("D784KFLB") != std::string::npos) radius_PMT[detkey] = map_type_radius["D784KFLB"];
    else if (dettype.find("EMI9954KB") != std::string::npos) radius_PMT[detkey] = map_type_radius["EMI9954KB"];

    double expectedT = (sqrt(pow(x_PMT.at(detkey)-diffuser_x,2)+pow(y_PMT.at(detkey)-diffuser_y,2)+pow(z_PMT.at(detkey)-diffuser_z,2))-radius_PMT[detkey])/c_vacuum*n_water*1E9;
    expected_time.insert(std::pair<unsigned long,double>(detkey,expectedT));
  }

  if (verbose > 1) std::cout <<"Number of tank PMTs: "<<n_tank_pmts<<std::endl;

  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){
    unsigned long detkey = pmt_detkeys.at(i_pmt);
    if (verbose > 1) std::cout <<"Detkey: "<<detkey<<", Radius PMT "<<i_pmt<<": "<<radius_PMT[detkey]<<", expectedT: "<<expected_time[detkey]<<std::endl;
  }

  std::vector<unsigned long>::iterator it_minkey = std::min_element(pmt_detkeys.begin(),pmt_detkeys.end());
  std::vector<unsigned long>::iterator it_maxkey = std::max_element(pmt_detkeys.begin(),pmt_detkeys.end());

  int min_detkey = std::distance(pmt_detkeys.begin(),it_minkey);
  int max_detkey = std::distance(pmt_detkeys.begin(),it_maxkey);
  int n_detkey_bins = max_detkey-min_detkey;

  //----------------------------------------------------------------------------
  //---------------Initialize Stability histograms -----------------------------
  //----------------------------------------------------------------------------

  std::vector<TH1F> charge_hist;
  hist_charge = new TH1F("hist_charge","Overall charge distribution (all PMTs)",nBinsChargeTotal,chargeTotalMin,chargeTotalMax);
  hist_time = new TH1F("hist_time","Overall time distribution (all PMTs)",nBinsTimeTotal,timeTotalMin,timeTotalMax);
  hist_tubeid = new TH1F("hist_tubeid","Overall Tube ID distribution",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_tubeid_adc = new TH1F("hist_tubeid_adc","Overall Tube ID distribution (ADCReco Store)",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_starttime = new TH1F("hist_starttime","Overall start time distribution (all PMTs)",nBinsStartTimeTotal,startTimeMin,startTimeMax);
  hist_peaktime = new TH1F("hist_peaktime","Overall peak time distribution (all PMTs)",nBinsPeakTimeTotal,peakTimeMin,peakTimeMax);
  hist_baseline = new TH1F("hist_baseline","Overall baseline distribution (all PMTs)",nBinsBaselineTotal,baselineMin,baselineMax);
  hist_sigmabaseline = new TH1F("hist_sigmabaseline","Overall sigma baseline distribution (all PMTs)",nBinsSigmaBaselineTotal,sigmaBaselineMin,sigmaBaselineMax);
  hist_rawamplitude = new TH1F("hist_rawamplitude","Overall raw amplitude distribution (all PMTs)",nBinsRawAmplitudeTotal,rawAmplitudeMin,rawAmplitudeMax);
  hist_amplitude = new TH1F("hist_amplitude","Overall amplitude distribution (all PMTs)",nBinsAmplitudeTotal,amplitudeMin,amplitudeMax);
  hist_rawarea = new TH1F("hist_rawarea","Overall raw area distribution (all PMTs)",nBinsRawAreaTotal,rawAreaMin,rawAreaMax);
  hist_charge_2D_y_phi_fit = new TH2F("hist_charge_2D_y_phi_fit","Spatial distribution of fit mean charges (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_time_2D_y_phi_dev_fit = new TH2F("hist_time_2D_y_phi_dev_fit","Spatial distribution of fit time deviations (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_time_2D_y_phi_fit = new TH2F("hist_time_2D_y_phi_fit","Spatial distribution of fit mean times (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_detkey_charge_fit = new TH1F("hist_detkey_charge_fit","Fit mean charges vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_time_fit = new TH1F("hist_detkey_time_fit","Fit mean times vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_time_dev_fit = new TH1F("hist_detkey_time_dev_fit","Fit time deviations vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_charge_2D_y_phi_mean = new TH2F("hist_charge_2D_y_phi_mean","Spatial distribution of mean charge (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_time_2D_y_phi_dev_mean = new TH2F("hist_time_2D_y_phi_dev_mean","Spatial distribution of mean time deviations (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_time_2D_y_phi_mean = new TH2F("hist_time_2D_y_phi_mean","Spatial distribution of mean time (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_starttime_2D_y_phi_mean = new TH2F("hist_starttime_2D_y_phi_mean","Spatial distribution of mean start time (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_peaktime_2D_y_phi_mean = new TH2F("hist_peaktime_2D_y_phi_mean","Spatial distribution of mean peak time (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_baseline_2D_y_phi_mean = new TH2F("hist_baseline_2D_y_phi_mean","Spatial distribution of mean baseline (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_sigmabaseline_2D_y_phi_mean = new TH2F("hist_sigmabaseline_2D_y_phi_mean","Spatial distribution of mean sigma baseline (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_rawamplitude_2D_y_phi_mean = new TH2F("hist_rawamplitude_2D_y_phi_mean","Spatial distribution of mean raw amplitude (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_amplitude_2D_y_phi_mean = new TH2F("hist_amplitude_2D_y_phi_mean","Spatial distribution of mean amplitude (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_rawarea_2D_y_phi_mean = new TH2F("hist_rawarea_2D_y_phi_mean","Spatial distribution of mean raw area (all PMTs)",100,0,360,25,-2.5,2.5);
  hist_detkey_charge_mean = new TH1F("hist_detkey_charge_mean","Mean charges vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_time_mean = new TH1F("hist_detkey_time_mean","Mean times vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_time_dev_mean = new TH1F("hist_detkey_time_dev_mean","Mean time deviations vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_starttime_mean = new TH1F("hist_detkey_starttime_mean","Mean start times vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_peaktime_mean = new TH1F("hist_detkey_peaktime_mean","Mean peak times vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_baseline_mean = new TH1F("hist_detkey_baseline_mean","Mean baseline vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_sigmabaseline_mean = new TH1F("hist_detkey_sigmabaseline_mean","Mean sigma baseline vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_rawamplitude_mean = new TH1F("hist_detkey_rawamplitude_mean","Mean raw amplitudes vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_amplitude_mean = new TH1F("hist_detkey_amplitude_mean","Mean amplitudes vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_rawarea_mean = new TH1F("hist_detkey_rawarea_mean","Mean raw areas vs. detkey",n_detkey_bins,pmt_detkeys[min_detkey],pmt_detkeys[max_detkey]);
  hist_detkey_2D_y_phi = new TH2F("hist_detkey_2D_y_phi","Spatial distribution of detkeys",100,0,360,25,-2.5,2.5);
  hist_occupied_2D_y_phi = new TH2F("hist_occupied_2D_y_phi","PMT Pixel Occupation in 2D space",100,0,360,25,-2.5,2.5);


  hist_charge->GetXaxis()->SetTitle("charge");
  hist_time->GetXaxis()->SetTitle("time [ns]");
  hist_tubeid->GetXaxis()->SetTitle("detectorkey");
  hist_tubeid_adc->GetXaxis()->SetTitle("detectorkey");
  hist_starttime->GetXaxis()->SetTitle("start time [ns]");
  hist_peaktime->GetXaxis()->SetTitle("peak time [ns]");
  hist_baseline->GetXaxis()->SetTitle("baseline [ADC]");
  hist_sigmabaseline->GetXaxis()->SetTitle("sigma baseline [ADC]");
  hist_rawamplitude->GetXaxis()->SetTitle("raw amplitude [ADC]");
  hist_amplitude->GetXaxis()->SetTitle("amplitude [V]");
  hist_rawarea->GetXaxis()->SetTitle("raw area [ADC x samples]");


  for (int i_tube=0;i_tube<n_tank_pmts;i_tube++){

    unsigned long detkey = pmt_detkeys[i_tube];
    stringstream ss;
    ss<<detkey;
    string detKey=ss.str();
    string name_general_charge="hist_charge_";
    string description_general_charge="Hit charges for detkey ";
    string name_hist_charge=name_general_charge+detKey;
    string description_hist_charge=description_general_charge+detKey;
    hist_charge_singletube[detkey] = new TH1F(name_hist_charge.c_str(),description_hist_charge.c_str(),nBinsCharge,chargeMin,chargeMax);

    string name_general_time="hist_time_";
    string description_general_time="Hit times for detkey ";
    string name_hist_time=name_general_time+detKey;
    string description_hist_time=description_general_time+detKey;
    hist_time_singletube[detkey] = new TH1F(name_hist_time.c_str(),description_hist_time.c_str(),nBinsTime,timeMin,timeMax);

    string name_general_starttime="hist_starttime_";
    string description_general_starttime="Start times for detkey ";
    string name_hist_starttime=name_general_starttime+detKey;
    string description_hist_starttime=description_general_starttime+detKey;
    hist_starttime_singletube[detkey] = new TH1F(name_hist_starttime.c_str(),description_hist_starttime.c_str(),nBinsStartTime,startTimeMin,startTimeMax);  

    string name_general_peaktime="hist_peaktime_";
    string description_general_peaktime="Peak times for detkey ";
    string name_hist_peaktime=name_general_peaktime+detKey;
    string description_hist_peaktime=description_general_peaktime+detKey;
    hist_peaktime_singletube[detkey] = new TH1F(name_hist_peaktime.c_str(),description_hist_peaktime.c_str(),nBinsPeakTime,peakTimeMin,peakTimeMax);  

    string name_general_baseline="hist_baseline_";
    string description_general_baseline="Baselines for detkey ";
    string name_hist_baseline=name_general_baseline+detKey;
    string description_hist_baseline=description_general_baseline+detKey;
    hist_baseline_singletube[detkey] = new TH1F(name_hist_baseline.c_str(),description_hist_baseline.c_str(),nBinsBaseline,baselineMin,baselineMax);  

    string name_general_sigmabaseline="hist_sigmabaseline_";
    string description_general_sigmabaseline="Sigma baselines for detkey ";
    string name_hist_sigmabaseline=name_general_sigmabaseline+detKey;
    string description_hist_sigmabaseline=description_general_sigmabaseline+detKey;
    hist_sigmabaseline_singletube[detkey] = new TH1F(name_hist_sigmabaseline.c_str(),description_hist_sigmabaseline.c_str(),nBinsSigmaBaseline,sigmaBaselineMin,sigmaBaselineMax); 

    string name_general_rawamplitude="hist_rawamplitude_";
    string description_general_rawamplitude="Raw amplitudes for detkey ";
    string name_hist_rawamplitude=name_general_rawamplitude+detKey;
    string description_hist_rawamplitude=description_general_rawamplitude+detKey;
    hist_rawamplitude_singletube[detkey] = new TH1F(name_hist_rawamplitude.c_str(),description_hist_rawamplitude.c_str(),nBinsRawAmplitude,rawAmplitudeMin,rawAmplitudeMax); 

    string name_general_amplitude="hist_amplitude_";
    string description_general_amplitude="Amplitudes for detkey ";
    string name_hist_amplitude=name_general_amplitude+detKey;
    string description_hist_amplitude=description_general_amplitude+detKey;
    hist_amplitude_singletube[detkey] = new TH1F(name_hist_amplitude.c_str(),description_hist_amplitude.c_str(),nBinsAmplitude,amplitudeMin,amplitudeMax);  

    string name_general_rawarea="hist_rawarea_";
    string description_general_rawarea="Raw areas for detkey ";
    string name_hist_rawarea=name_general_rawarea+detKey;
    string description_hist_rawarea=description_general_rawarea+detKey;
    hist_rawarea_singletube[detkey] = new TH1F(name_hist_rawarea.c_str(),description_hist_rawarea.c_str(),nBinsRawArea,rawAreaMin,rawAreaMax);    

    hist_charge_singletube[detkey]->GetXaxis()->SetTitle("charge");
    hist_time_singletube[detkey]->GetXaxis()->SetTitle("time [ns]");
    hist_starttime_singletube[detkey]->GetXaxis()->SetTitle("start time [ns]");
    hist_peaktime_singletube[detkey]->GetXaxis()->SetTitle("peak time [ns]");
    hist_baseline_singletube[detkey]->GetXaxis()->SetTitle("baseline [ADC]");
    hist_sigmabaseline_singletube[detkey]->GetXaxis()->SetTitle("sigma baseline [ADC]");
    hist_rawamplitude_singletube[detkey]->GetXaxis()->SetTitle("raw amplitude [ADC]");
    hist_amplitude_singletube[detkey]->GetXaxis()->SetTitle("amplitude [V]");
    hist_rawarea_singletube[detkey]->GetXaxis()->SetTitle("raw area [ADC x samples]");

  }


  hist_charge_fit = new TH1F("hist_charge_fit","Fit Mean values of detected charges",nBinsChargeFit,chargeFitMin,chargeFitMax);
  hist_time_fit = new TH1F("hist_time_fit","Fit Mean values of detected hit times",nBinsTimeFit,timeFitMin,timeFitMax);
  hist_time_dev_fit = new TH1F("hist_time_dev_fit","Fit Deviation of detected hit times",nBinsTimeDev,timeDevMin,timeDevMax);
  hist_charge_mean = new TH1F("hist_charge_mean","Mean values of detected charges",nBinsChargeFit,chargeFitMin,chargeFitMax);
  hist_time_mean = new TH1F("hist_time_mean","Mean values of detected hit times",nBinsTimeFit,timeFitMin,timeFitMax);
  hist_time_dev_mean = new TH1F("hist_time_dev","Deviation of detected hit times",nBinsTimeDev,timeDevMin,timeDevMax);
  hist_starttime_mean = new TH1F("hist_starttime_mean","Mean values of detected start times",100,startTimeMin,startTimeMax);
  hist_peaktime_mean = new TH1F("hist_peaktime_mean","Mean values of detected peak times",100,peakTimeMin,peakTimeMax);
  hist_baseline_mean = new TH1F("hist_baseline_mean","Mean values of detected baselines",100,baselineMin,baselineMax);
  hist_sigmabaseline_mean = new TH1F("hist_sigmabaseline_mean","Mean values of detected sigma baselines",100,sigmaBaselineMin,sigmaBaselineMax);
  hist_rawamplitude_mean = new TH1F("hist_rawamplitude_mean","Mean values of detected raw amplitudes",100,rawAmplitudeMin,rawAmplitudeMax);
  hist_amplitude_mean = new TH1F("hist_amplitude_mean","Mean values of detected amplitudes",100,amplitudeMin,amplitudeMax);
  hist_rawarea_mean = new TH1F("hist_rawarea_mean","Mean values of detected raw areas",100,rawAreaMin,rawAreaMax);


  //Setting PMT pixel occupation histogram to 0 for all bins
  for (int i_x = 0; i_x < hist_occupied_2D_y_phi->GetNbinsX(); i_x++){
    for (int i_y = 0; i_y < hist_occupied_2D_y_phi->GetNbinsY(); i_y++){
      hist_occupied_2D_y_phi->SetBinContent(i_x+1,i_y+1,0);
    }
  }



  //---------------------------------------------------------------------------------
  //create root-file that will contain all analysis graphs for this calibration run--
  //---------------------------------------------------------------------------------

  //m_data->CStore.Get("FileRunNumber",runnumber);
  //m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  /*std::cout <<"TankCalibrationDiffuser tool: Got RunNumber "<<runnumber<<" from the CStore Boost Store"<<std::endl;
  std::stringstream ss_run;
  ss_run<<runnumber;
  std::string newRunNumber = ss_run.str();
  std::string file_out_prefix="_Run";
  std::string file_out_root=".root";
  std::string file_out_name=outputdir+outputfile+file_out_prefix+newRunNumber+file_out_root;

  file_out=new TFile(file_out_name.c_str(),"RECREATE"); //create one root file for each run to save the detailed plots and fits for all PMTs 
  file_out->cd();*/

  //----------------------------------------------------------------------------
  //---Make histograms appear during execution by declaring TApplication--------
  //----------------------------------------------------------------------------

  if (use_tapplication){
    if (verbose > 0) std::cout <<"Opening TApplication..."<<std::endl;
    int myargc = 0;
    char *myargv[] {(char*) "options"};
    app_stability = new TApplication("app_stability",&myargc,myargv);
    if (verbose > 0) std::cout <<"TApplication initialized..."<<std::endl;
  }

  return true;

}


bool TankCalibrationDiffuser::Execute(){

  if (verbose > 0) std::cout <<"Executing Tool TankCalibrationDiffuser ..."<<endl;

  // get the ANNIEEvent

  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){ std::cerr<<"Error: No ANNIEEvent store!"<<endl; /*return false;*/};

  std::map<unsigned long, std::vector<std::vector<ADCPulse>>> RecoADCHits; 

  //----------------------------------------------------------------------------
  //---------------get the members of the ANNIEEvent----------------------------
  //----------------------------------------------------------------------------
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  m_data->Stores["ANNIEEvent"]->Get("BeamStatus", BeamStatus);
  bool got_recoadc = m_data->Stores["ANNIEEvent"]->Get("RecoADCHits",RecoADCHits);

  if (HitStoreName == "MCHits"){
    bool got_mchits = m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
    if (!got_mchits){
      std::cout << "No MCHits store in ANNIEEvent!" <<  std::endl;
      return false;
    }
  } else if (HitStoreName == "Hits"){
    bool got_hits = m_data->Stores["ANNIEEvent"]->Get("Hits", Hits);
    if (!got_hits){
      std::cout << "No Hits store in ANNIEEvent! " << std::endl;
      return false;
    }
  } else {
    std::cout << "Selected Hits store invalid.  Must be Hits or MCHits" << std::endl;
    return false;
  }
  // Also load hits from the Hits Store, if available

  //----------------------------------------------------------------------------
  //---------------Read out charge and hit values of PMTs-----------------------
  //----------------------------------------------------------------------------
  
  for (int i_pmt = 0; i_pmt < n_tank_pmts; i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    PMT_ishit[detkey] = 0;
  }

  if(HitStoreName=="MCHits"){
    int vectsize = MCHits->size();
    if (verbose > 0) std::cout <<"MCHits size: "<<vectsize<<std::endl; 
    for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits){
      unsigned long chankey = apair.first;
      Detector* thistube = geom->ChannelToDetector(chankey);
      int detectorkey = thistube->GetDetectorID();
      if (thistube->GetDetectorElement()=="Tank"){
        std::vector<MCHit>& ThisPMTHits = apair.second;
        PMT_ishit[detectorkey] = 1;
        for (MCHit &ahit : ThisPMTHits){
          if (verbose > 2) std::cout <<"Charge "<<ahit.GetCharge()<<", time "<<ahit.GetTime()<<std::endl;
          hist_charge->Fill(ahit.GetCharge());
          hist_time->Fill(ahit.GetTime());
          hist_tubeid->Fill(detectorkey);
          hist_charge_singletube[detectorkey]->Fill(ahit.GetCharge());
          hist_time_singletube[detectorkey]->Fill(ahit.GetTime());
        }
      }
    }
  }
  
  if(HitStoreName=="Hits"){
    int vectsize = Hits->size();
    if (verbose > 0) std::cout <<"Hits size: "<<vectsize<<std::endl; 
    for(std::pair<unsigned long, std::vector<Hit>>&& apair : *Hits){
      unsigned long chankey = apair.first;
      Detector* thistube = geom->ChannelToDetector(chankey);
      int detectorkey = thistube->GetDetectorID();
      if (thistube->GetDetectorElement()=="Tank"){
        std::vector<Hit>& ThisPMTHits = apair.second;
        PMT_ishit[detectorkey] = 1;
        for (Hit &ahit : ThisPMTHits){
          if (verbose > 2) std::cout <<"Charge "<<ahit.GetCharge()<<", time "<<ahit.GetTime()<<std::endl;
          hist_charge->Fill(ahit.GetCharge());
          hist_time->Fill(ahit.GetTime());
          hist_tubeid->Fill(detectorkey);
          hist_charge_singletube[detectorkey]->Fill(ahit.GetCharge());
          hist_time_singletube[detectorkey]->Fill(ahit.GetTime());
        }
      }
    }
  }

  //check whether PMT_ishit is filled correctly
  for (int i_pmt = 0; i_pmt < n_tank_pmts ; i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];  
    if (verbose > 2 && PMT_ishit[detkey]) std::cout <<"PMT (detkey "<<detkey<<") is hit."<<std::endl;
  }

  //----------------------------------------------------------------------------
  //---------------Read out RecoADCHits properties of PMTs----------------------
  //---------------------------------------------------------------------------- 

  if (got_recoadc){

    int recoadcsize = RecoADCHits.size();
    int adc_loop = 0;
    if (verbose > 0) std::cout <<"RecoADCHits size: "<<recoadcsize<<std::endl;
    for (std::pair<unsigned long, std::vector<std::vector<ADCPulse>>> apair : RecoADCHits){
      unsigned long chankey = apair.first;
      Detector *thistube = geom->ChannelToDetector(chankey);
      int detectorkey = thistube->GetDetectorID();
      if (thistube->GetDetectorElement()=="Tank"){
        std::vector<std::vector<ADCPulse>> pulses = apair.second;
        for (int i_minibuffer = 0; i_minibuffer < (int) pulses.size(); i_minibuffer++){
          std::vector<ADCPulse> apulsevector = pulses.at(i_minibuffer);
          for (int i_pulse=0; i_pulse < (int) apulsevector.size(); i_pulse++){
            ADCPulse apulse = apulsevector.at(i_pulse);
            double start_time = apulse.start_time();
            double peak_time = apulse.peak_time();
            double baseline = apulse.baseline();
            double sigma_baseline = apulse.sigma_baseline();
            double raw_amplitude = apulse.raw_amplitude();
            double amplitude = apulse.amplitude();
            double raw_area = apulse.raw_area();
            hist_starttime->Fill(start_time);
            hist_peaktime->Fill(peak_time);
            hist_baseline->Fill(baseline);
            hist_sigmabaseline->Fill(sigma_baseline);
            hist_rawamplitude->Fill(raw_amplitude);
            hist_amplitude->Fill(amplitude);
            hist_rawarea->Fill(raw_area);
            hist_tubeid_adc->Fill(detectorkey);
            hist_starttime_singletube[detectorkey]->Fill(start_time);
            hist_peaktime_singletube[detectorkey]->Fill(peak_time);
            hist_baseline_singletube[detectorkey]->Fill(baseline);
            hist_sigmabaseline_singletube[detectorkey]->Fill(sigma_baseline);
            hist_rawamplitude_singletube[detectorkey]->Fill(raw_amplitude);
            hist_amplitude_singletube[detectorkey]->Fill(amplitude);
            hist_rawarea_singletube[detectorkey]->Fill(raw_area);
          }
        }
      }else {
        if (verbose > 0) std::cout <<"TankCalibrationDiffuser: RecoADCHit does not belong to a tank PMT and is ommited. Detector key/element = "<<detectorkey<<" / "<<thistube->GetDetectorElement()<<std::endl;
      }
      adc_loop++;

    }



  } else {

    std::cout <<"TankCalibrationDiffuser: RecoADCHits Store does not exist and is not read out"<<std::endl;
  }

  return true;
}


bool TankCalibrationDiffuser::Finalise(){
  
  if (verbose > 0) std::cout <<"TankCalibrationDiffuser: Finalise"<<std::endl;
  m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  std::cout <<"TankCalibrationDiffuser tool: Got RunNumber "<<runnumber<<" from the CStore Boost Store"<<std::endl;
  std::stringstream ss_run;
  ss_run<<runnumber;
  std::string newRunNumber = ss_run.str();
  std::string file_out_prefix="_Run";
  std::string file_out_root=".root";
  std::string file_out_name=outputdir+outputfile+file_out_prefix+newRunNumber+file_out_root;

  file_out=new TFile(file_out_name.c_str(),"RECREATE"); //create one root file for each run to save the detailed plots and fits for all PMTs 
  file_out->cd();

  //----------------------------------------------------------------------------
  //---------------Create and write single run root files------------------------
  //----------------------------------------------------------------------------

  //cout <<"pointer histogram charge: "<<hist_charge<<endl;
  //cout <<"pointer singletube hist charge: "<<hist_charge_singletube[23]<<endl;
  //hist_charge->Print();
  hist_charge->Write();
  hist_time->Write();
  hist_tubeid->Write();
  hist_tubeid_adc->Write();
  hist_starttime->Write();
  hist_peaktime->Write();
  hist_baseline->Write();
  hist_sigmabaseline->Write();
  hist_rawamplitude->Write();
  hist_amplitude->Write();
  hist_rawarea->Write();

  // create txt file to store calibration data for the specific run
  std::string filename_pre = "_Run";
  std::string filename_post = ".txt";
  std::string filename_complete = outputdir+outputfile+filename_pre+newRunNumber+filename_post;
  std::cout <<"writing to txt file: "<<filename_complete<<std::endl;
  ofstream result_file(filename_complete.c_str());
  result_file << "detkey"<<"    "<<"charge_mean_fit"<<"    "<<"charge_rms_fit"<<"    "<<"charge_mean"<<"    "<<"charge_rms"<<"    "<<"time_mean_fit"<<"    "<<"time_rms_fit"<<"    "<<"time_mean"<<"    "<<"time_rms"<<"    "<<"expected_time"<<"    "<<"time_dev_fit"<<"    "<<"time_dev"<<"    "<<"nentries"<<std::endl;

  //----------------------------------------------------------------------------
  //---------------Perform fits for PMT distributions---------------------------
  //----------------------------------------------------------------------------

  for (int i_tube=0;i_tube<n_tank_pmts;i_tube++){

    unsigned long detkey = pmt_detkeys[i_tube];
    double mean_charge=hist_charge_singletube[detkey]->GetMean();
    double rms_charge = hist_charge_singletube[detkey]->GetRMS();
    double mean_time = hist_time_singletube[detkey]->GetMean();
    double rms_time = hist_time_singletube[detkey]->GetRMS();
    int nentries = hist_charge_singletube[detkey]->GetEntries();

    charge_mean[detkey] = mean_charge;
    charge_rms[detkey] = rms_charge;
    if (std::find(problematic_channels.begin(),problematic_channels.end(),detkey)==problematic_channels.end()){
    time_mean[detkey] = mean_time;
    time_rms[detkey] = rms_time;
    }else {
    time_mean[detkey] = 0.;
    time_rms[detkey] = 0.;
    }
    nentries_hist[detkey] = nentries;

    Double_t par_gaus2exp[8] = {gaus1Constant,gaus1Mean,gaus1Sigma,gaus2Constant,gaus2Mean,gaus2Sigma,expConstant,expDecay};    //old default: {10,0.3,0.1,10,1.0,0.5,1,-1}
    Double_t par_gaus2[6] = {gaus1Constant,gaus1Mean,gaus1Sigma,gaus2Constant,gaus2Mean,gaus2Sigma};  //old default: {10,0.3,0.1,10,1.0,0.5}
    Double_t par_gaus[3] = {gaus1Constant,gaus1Mean,gaus1Sigma};    //old default: {10,1.0,0.5}

    TF1 *total;
    if (FitMethod == "Gaus2Exp") total = new TF1("total","gaus(0)+gaus(3)+expo(6)",chargeMin,chargeMax);
    else if (FitMethod == "Gaus2") total = new TF1("total","gaus(0)+gaus(3)",chargeMin,chargeMax);
    else if (FitMethod == "Gaus") total = new TF1("total","gaus",chargeMin,chargeMax);
    else {
      std::cout <<"ERROR (TankCalibrationDiffuser): FitFunction is not part of the options, please extend the options Using standard Gaus."<<std::endl;
      total = new TF1("total","gaus",chargeMin,chargeMax);
    }
    total->SetLineColor(2);
    if (FitMethod == "Gaus2Exp") total->SetParameters(par_gaus2exp);
    else if (FitMethod == "Gaus2") total->SetParameters(par_gaus2);
    else total->SetParameters(par_gaus2);

    //total = new TF1("total",total_value,0,10,8);
    //total->SetParameters(10,0.3,0.1,10,1.,0.3,1,1);
    TFitResultPtr FitResult = hist_charge_singletube[detkey]->Fit(total,"QR+");
    Int_t FitResult_int = FitResult;

    //hist_charge_singletube[i_tube]->Fit("gaus[0]+gaus[3]+exp[6]");
    hist_charge_singletube[detkey]->Write();
    if (FitResult_int ==0){
      TF1 *fit_result_charge=hist_charge_singletube[detkey]->GetFunction("total");
      if (FitMethod == "Gaus2Exp") {
        charge_mean_fit[detkey] = fit_result_charge->GetParameter(4);
        charge_rms_fit[detkey] = fit_result_charge->GetParameter(5);
      } 
      else if (FitMethod == "Gaus2") {
        charge_mean_fit[detkey]=fit_result_charge->GetParameter(4);
        charge_rms_fit[detkey]=fit_result_charge->GetParameter(5);
      }
      else {
      charge_mean_fit[detkey]=fit_result_charge->GetParameter(1);
      charge_rms_fit[detkey]=fit_result_charge->GetParameter(2);
      }
    } else {
      charge_mean_fit[detkey] = 0.;
      charge_rms_fit[detkey] = 0.;
    }
    //check if fit returned negative standard deviation
    if (charge_rms_fit[detkey] < 0.) charge_rms_fit[detkey] = -charge_rms_fit[detkey];

    hist_charge_fit->Fill(charge_mean_fit[detkey]);  
    hist_charge_mean->Fill(charge_mean[detkey]);
    if (verbose > 2){
      std::cout <<"Fit mean charge for detkey "<<detkey<<": "<<charge_mean_fit[detkey]<<std::endl;
      std::cout <<"Mean charge for detkey "<<detkey<<": "<<charge_mean[detkey]<<std::endl;
    }

    TFitResultPtr FitResultTime = hist_time_singletube[detkey]->Fit("gaus","Q");      //time fit is always the same for now, assume simple Gaussian
    Int_t FitResultTime_int = FitResultTime;
    hist_time_singletube[detkey]->Write();
    if (FitResultTime_int == 0) { //fit was okay and has a result
      TF1 *fit_result_time=hist_time_singletube[detkey]->GetFunction("gaus");
      time_mean_fit[detkey]=fit_result_time->GetParameter(1);
      time_rms_fit[detkey]=fit_result_time->GetParameter(2);
    }
    else {
      time_mean_fit[detkey] = 0.;
      time_rms_fit[detkey] = 0.;
    }
    hist_time_dev_fit->Fill(time_mean_fit[detkey]-expected_time[detkey]);
    hist_time_fit->Fill(time_mean_fit[detkey]);
    hist_time_dev_mean->Fill(time_mean[detkey]-expected_time[detkey]);
    hist_time_mean->Fill(time_mean[detkey]);

    if (verbose > 2){
      std::cout <<"expected hit time (detkey "<<detkey<<") : "<<expected_time[detkey]<<endl;
      std::cout <<"detected hit time (fit) (detkey "<<detkey<<") : "<<time_mean_fit[detkey]<<endl;
      std::cout << "deviation hit time (fit) - expected time (detkey "<<detkey<<") : "<<time_mean_fit[detkey]-expected_time[detkey]<<endl;
      std::cout <<"detected hit time (mean) (detkey "<<detkey<<" ) : "<<time_mean[detkey]<<endl;
      std::cout << "deviation hit time (mean) - expected time (detkey "<<detkey<<") : "<<time_mean[detkey]-expected_time[detkey]<<endl;
    }

    if (charge_mean_fit[detkey] < 0.){        //unphysical charge information --> set to 0
      charge_mean_fit[detkey] = 0.;
      charge_rms_fit[detkey] = 0.;
    }

    //
    //write ADCRecoPulse histograms to file as well
    //


    hist_starttime_singletube[detkey]->Write();
    hist_peaktime_singletube[detkey]->Write();
    hist_baseline_singletube[detkey]->Write();
    hist_sigmabaseline_singletube[detkey]->Write();
    hist_rawamplitude_singletube[detkey]->Write();
    hist_amplitude_singletube[detkey]->Write();
    hist_rawarea_singletube[detkey]->Write();

    
    if (std::find(problematic_channels.begin(),problematic_channels.end(),detkey)==problematic_channels.end()){
    starttime_mean[detkey] = hist_starttime_singletube[detkey]->GetMean();
    peaktime_mean[detkey] = hist_peaktime_singletube[detkey]->GetMean();
    } else {
    starttime_mean[detkey] = 0.;
    peaktime_mean[detkey] = 0.;
    }
    baseline_mean[detkey] = hist_baseline_singletube[detkey]->GetMean();
    sigmabaseline_mean[detkey] = hist_sigmabaseline_singletube[detkey]->GetMean();
    rawamplitude_mean[detkey] = hist_rawamplitude_singletube[detkey]->GetMean();
    amplitude_mean[detkey] = hist_amplitude_singletube[detkey]->GetMean();
    rawarea_mean[detkey] = hist_rawarea_singletube[detkey]->GetMean();
	
    hist_starttime_mean->Fill(starttime_mean[detkey]);
    hist_peaktime_mean->Fill(peaktime_mean[detkey]);
    hist_baseline_mean->Fill(baseline_mean[detkey]);
    hist_sigmabaseline_mean->Fill(sigmabaseline_mean[detkey]);
    hist_rawamplitude_mean->Fill(rawamplitude_mean[detkey]);
    hist_amplitude_mean->Fill(amplitude_mean[detkey]);
    hist_rawarea_mean->Fill(rawarea_mean[detkey]);
	
    //
    //fill spatial detector plots as well
    //

    if (hist_occupied_2D_y_phi->GetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1)<1.){
      hist_time_2D_y_phi_dev_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,fabs(time_mean_fit[detkey]-expected_time[detkey]));  
      hist_time_2D_y_phi_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,time_mean_fit[detkey]);  
      hist_charge_2D_y_phi_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,charge_mean_fit[detkey]);
    } else {
      hist_time_2D_y_phi_dev_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,fabs(time_mean_fit[detkey]-expected_time[detkey]));
      hist_time_2D_y_phi_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,time_mean_fit[detkey]);
      hist_charge_2D_y_phi_fit->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,charge_mean_fit[detkey]);
    }

    if (hist_occupied_2D_y_phi->GetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1)<1.){
      hist_time_2D_y_phi_dev_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,fabs(time_mean[detkey]-expected_time[detkey]));  
      hist_time_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,time_mean[detkey]);  
      hist_charge_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,charge_mean[detkey]);
      hist_starttime_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,starttime_mean[detkey]);
      hist_peaktime_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,peaktime_mean[detkey]);
      hist_baseline_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,baseline_mean[detkey]);
      hist_sigmabaseline_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,sigmabaseline_mean[detkey]);
      hist_rawamplitude_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,rawamplitude_mean[detkey]);
      hist_amplitude_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,amplitude_mean[detkey]);
      hist_rawarea_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,rawarea_mean[detkey]);
      hist_detkey_2D_y_phi->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,detkey);    
      hist_occupied_2D_y_phi->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+1,1);
    }
    else {
      hist_time_2D_y_phi_dev_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,fabs(time_mean[detkey]-expected_time[detkey]));
      hist_time_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,time_mean[detkey]);
      hist_charge_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,charge_mean[detkey]);
      hist_starttime_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,starttime_mean[detkey]);
      hist_peaktime_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,peaktime_mean[detkey]);
      hist_baseline_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,baseline_mean[detkey]);
      hist_sigmabaseline_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,sigmabaseline_mean[detkey]);
      hist_rawamplitude_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,rawamplitude_mean[detkey]);
      hist_amplitude_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,amplitude_mean[detkey]);
      hist_rawarea_2D_y_phi_mean->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,rawarea_mean[detkey]);
      hist_detkey_2D_y_phi->SetBinContent(int(phi_PMT[detkey]/2/TMath::Pi()*100)+1,int((y_PMT[detkey]+2.5)/5.*25)+2,detkey);
    }
    result_file<<detkey<<"  "<<charge_mean_fit[detkey]<<"  "<<charge_rms_fit[detkey]<<"  "<<charge_mean[detkey]<<"  "<<charge_rms[detkey]<<"  "<<time_mean_fit[detkey]<<"  "<<time_rms_fit[detkey]<<"  "<<time_mean[detkey]<<"  "<<time_rms[detkey]<<"  "<<expected_time[detkey]<<"  "<<time_mean_fit[detkey]-expected_time[detkey]<<"  "<<time_mean[detkey]-expected_time[detkey]<<"  "<<nentries_hist[detkey]<<endl;

    //
    //fill detkey 1D histograms as well
    //

    int bin_nr = hist_detkey_charge_mean->FindBin(detkey);
    hist_detkey_charge_fit->SetBinContent(bin_nr,charge_mean_fit[detkey]);
    hist_detkey_time_fit->SetBinContent(bin_nr,time_mean_fit[detkey]);
    hist_detkey_time_dev_fit->SetBinContent(bin_nr,time_mean_fit[detkey]-expected_time[detkey]);
    hist_detkey_charge_mean->SetBinContent(bin_nr,charge_mean[detkey]);
    hist_detkey_time_mean->SetBinContent(bin_nr,time_mean[detkey]);
    hist_detkey_time_dev_mean->SetBinContent(bin_nr,time_mean[detkey]-expected_time[detkey]);
    hist_detkey_starttime_mean->SetBinContent(bin_nr,starttime_mean[detkey]);
    hist_detkey_peaktime_mean->SetBinContent(bin_nr,peaktime_mean[detkey]);
    hist_detkey_baseline_mean->SetBinContent(bin_nr,baseline_mean[detkey]);
    hist_detkey_sigmabaseline_mean->SetBinContent(bin_nr,sigmabaseline_mean[detkey]);
    hist_detkey_rawamplitude_mean->SetBinContent(bin_nr,rawamplitude_mean[detkey]);
    hist_detkey_amplitude_mean->SetBinContent(bin_nr,amplitude_mean[detkey]);
    hist_detkey_rawarea_mean->SetBinContent(bin_nr,rawarea_mean[detkey]);

    vector_tf1.push_back(total);
  }
  result_file<<1000<<"  "<<hist_charge_fit->GetMean()<<"  "<<hist_charge_fit->GetRMS()<<"  "<<hist_charge_mean->GetMean()<<"  "<<hist_charge_mean->GetRMS()<<"  "<<hist_time_fit->GetMean()<<"  "<<hist_time_fit->GetRMS()<<"  "<<hist_time_mean->GetMean()<<"  "<<hist_time_mean->GetRMS()<<"  "<<0<<"  "<<hist_time_dev_fit->GetMean()<<"  "<<hist_time_dev_mean->GetMean()<<"  "<<hist_charge->GetEntries()<<endl; //1000 is identifier key for average value
  result_file.close();

  hist_time_fit->GetXaxis()->SetTitle("t_{arrival} [ns]");
  hist_time_dev_fit->GetXaxis()->SetTitle("t_{arrival} - t_{expected} [ns]");
  hist_charge_fit->GetXaxis()->SetTitle("charge [p.e.]");
  hist_time_mean->GetXaxis()->SetTitle("t_{arrival} [ns]");
  hist_time_dev_mean->GetXaxis()->SetTitle("t_{arrival} - t_{expected} [ns]");
  hist_charge_mean->GetXaxis()->SetTitle("charge [p.e.]");
  hist_starttime_mean->GetXaxis()->SetTitle("start time [ns]");
  hist_peaktime_mean->GetXaxis()->SetTitle("peak time [ns]");
  hist_baseline_mean->GetXaxis()->SetTitle("baseline [ADC]");
  hist_sigmabaseline_mean->GetXaxis()->SetTitle("sigma baseline [ADC]");
  hist_rawamplitude_mean->GetXaxis()->SetTitle("raw amplitude [ADC]");
  hist_amplitude_mean->GetXaxis()->SetTitle("amplitude [V]");
  hist_rawarea_mean->GetXaxis()->SetTitle("raw area [ADC x pulses]");

  hist_detkey_2D_y_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_detkey_2D_y_phi->GetYaxis()->SetTitle("y [m]");
  hist_detkey_2D_y_phi->SetStats(0);
  hist_detkey_2D_y_phi->SetDrawOption("colz text");
  hist_charge_2D_y_phi_fit->GetXaxis()->SetTitle("#phi [deg]");
  hist_charge_2D_y_phi_fit->GetYaxis()->SetTitle("y [m]");
  hist_charge_2D_y_phi_fit->SetStats(0);
  hist_time_2D_y_phi_dev_fit->GetXaxis()->SetTitle("#phi [deg]");
  hist_time_2D_y_phi_dev_fit->GetYaxis()->SetTitle("y [m]");
  hist_time_2D_y_phi_dev_fit->SetStats(0);
  hist_time_2D_y_phi_fit->GetXaxis()->SetTitle("#phi [deg]");
  hist_time_2D_y_phi_fit->GetYaxis()->SetTitle("y [m]");
  hist_time_2D_y_phi_fit->SetStats(0);
  hist_charge_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_charge_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_charge_2D_y_phi_mean->SetStats(0);
  hist_time_2D_y_phi_dev_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_time_2D_y_phi_dev_mean->GetYaxis()->SetTitle("y [m]");
  hist_time_2D_y_phi_dev_mean->SetStats(0);
  hist_time_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_time_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_time_2D_y_phi_mean->SetStats(0);
  hist_starttime_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_starttime_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_starttime_2D_y_phi_mean->SetStats(0);
  hist_peaktime_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_peaktime_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_peaktime_2D_y_phi_mean->SetStats(0);
  hist_baseline_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_baseline_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_baseline_2D_y_phi_mean->SetStats(0);
  hist_sigmabaseline_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_sigmabaseline_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_sigmabaseline_2D_y_phi_mean->SetStats(0);
  hist_rawamplitude_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_rawamplitude_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_rawamplitude_2D_y_phi_mean->SetStats(0);
  hist_amplitude_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_amplitude_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_amplitude_2D_y_phi_mean->SetStats(0);
  hist_rawarea_2D_y_phi_mean->GetXaxis()->SetTitle("#phi [deg]");
  hist_rawarea_2D_y_phi_mean->GetYaxis()->SetTitle("y [m]");
  hist_rawarea_2D_y_phi_mean->SetStats(0);

  hist_detkey_charge_fit->GetXaxis()->SetTitle("detkey");
  hist_detkey_charge_fit->GetYaxis()->SetTitle("Fit mean charge");
  hist_detkey_charge_fit->SetStats(0);
  hist_detkey_time_fit->GetXaxis()->SetTitle("detkey");
  hist_detkey_time_fit->GetYaxis()->SetTitle("Fit mean time [ns]");
  hist_detkey_time_fit->SetStats(0);
  hist_detkey_time_dev_fit->GetXaxis()->SetTitle("detkey");
  hist_detkey_time_dev_fit->GetYaxis()->SetTitle("Fit time deviation [ns]");
  hist_detkey_time_dev_fit->SetStats(0);
  hist_detkey_charge_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_charge_mean->GetYaxis()->SetTitle("Mean charge");
  hist_detkey_charge_mean->SetStats(0);
  hist_detkey_time_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_time_mean->GetYaxis()->SetTitle("Mean time [ns]");
  hist_detkey_time_mean->SetStats(0);
  hist_detkey_time_dev_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_time_dev_mean->GetYaxis()->SetTitle("Time deviation [ns]");
  hist_detkey_time_dev_mean->SetStats(0);
  hist_detkey_starttime_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_starttime_mean->GetYaxis()->SetTitle("Mean start time [ns]");
  hist_detkey_starttime_mean->SetStats(0);
  hist_detkey_peaktime_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_peaktime_mean->GetYaxis()->SetTitle("Mean peak time [ns]");
  hist_detkey_peaktime_mean->SetStats(0);
  hist_detkey_baseline_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_baseline_mean->GetYaxis()->SetTitle("Mean baseline [ADC]");
  hist_detkey_baseline_mean->SetStats(0);
  hist_detkey_sigmabaseline_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_sigmabaseline_mean->GetYaxis()->SetTitle("Mean sigma baseline [ADC]");
  hist_detkey_sigmabaseline_mean->SetStats(0);
  hist_detkey_rawamplitude_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_rawamplitude_mean->GetYaxis()->SetTitle("Mean raw amplitude [ADC]");
  hist_detkey_rawamplitude_mean->SetStats(0);
  hist_detkey_amplitude_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_amplitude_mean->GetYaxis()->SetTitle("Mean amplitude [V]");
  hist_detkey_amplitude_mean->SetStats(0);
  hist_detkey_rawarea_mean->GetXaxis()->SetTitle("detkey");
  hist_detkey_rawarea_mean->GetYaxis()->SetTitle("Mean raw area [ADC x samples]");
  hist_detkey_rawarea_mean->SetStats(0);

  hist_time_fit->Write();
  hist_time_dev_fit->Write();
  hist_charge_fit->Write();
  hist_time_mean->Write();
  hist_time_dev_mean->Write();
  hist_charge_mean->Write();
  hist_starttime_mean->Write();
  hist_peaktime_mean->Write();
  hist_baseline_mean->Write();
  hist_sigmabaseline_mean->Write();
  hist_rawamplitude_mean->Write();
  hist_amplitude_mean->Write();
  hist_rawarea_mean->Write();

  hist_charge_2D_y_phi_fit->Write();
  hist_time_2D_y_phi_fit->Write();
  hist_time_2D_y_phi_dev_fit->Write();
  hist_charge_2D_y_phi_mean->Write();
  hist_time_2D_y_phi_mean->Write();
  hist_time_2D_y_phi_dev_mean->Write();
  hist_starttime_2D_y_phi_mean->Write();
  hist_peaktime_2D_y_phi_mean->Write();
  hist_baseline_2D_y_phi_mean->Write();
  hist_sigmabaseline_2D_y_phi_mean->Write();
  hist_rawamplitude_2D_y_phi_mean->Write();
  hist_amplitude_2D_y_phi_mean->Write();
  hist_rawarea_2D_y_phi_mean->Write();
  hist_detkey_2D_y_phi->Write();

  hist_detkey_charge_fit->Write();
  hist_detkey_time_fit->Write();
  hist_detkey_time_dev_fit->Write();
  hist_detkey_charge_mean->Write();
  hist_detkey_time_mean->Write();
  hist_detkey_time_dev_mean->Write();
  hist_detkey_starttime_mean->Write();
  hist_detkey_peaktime_mean->Write();
  hist_detkey_baseline_mean->Write();
  hist_detkey_sigmabaseline_mean->Write();
  hist_detkey_rawamplitude_mean->Write();
  hist_detkey_amplitude_mean->Write();
  hist_detkey_rawarea_mean->Write();



  //----------------------------------------------------------------------------
  //---------------Create and write stability plots------------------------
  //----------------------------------------------------------------------------

  gr_stability_charge_fit = new TGraphErrors();
  gr_stability_charge_fit->SetName("gr_stability_charge_fit");
  gr_stability_charge_fit->SetTitle("Stability PMT fit mean charges");
  gr_stability_charge_mean = new TGraphErrors();
  gr_stability_charge_mean->SetName("gr_stability_charge_mean");
  gr_stability_charge_mean->SetTitle("Stability PMT mean charges");
  gr_stability_time_fit = new TGraphErrors();
  gr_stability_time_fit->SetName("gr_stability_time_fit");
  gr_stability_time_fit->SetTitle("Stability PMT fit mean times");
  gr_stability_time_mean = new TGraphErrors();
  gr_stability_time_mean->SetName("gr_stability_time_mean");
  gr_stability_time_mean->SetTitle("Stability PMT mean times");
  gr_stability_time_dev_fit = new TGraphErrors();
  gr_stability_time_dev_fit->SetName("gr_stability_time_fit");
  gr_stability_time_dev_fit->SetTitle("Stability PMT fit mean time deviations");
  gr_stability_time_dev_mean = new TGraphErrors();
  gr_stability_time_dev_mean->SetName("gr_stability_time_mean");
  gr_stability_time_dev_mean->SetTitle("Stability PMT mean time deviations");


  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){

    unsigned long detkey = pmt_detkeys[i_pmt];
    stringstream ss;
    ss<<detkey;
    string detKey=ss.str(); 
    string name_general_stability_fit="gr_stability_charge_fit_";
    string description_general_stability_fit="Stability of fit mean charge for detkey ";
    string name_gr_stability_fit=name_general_stability_fit+detKey;
    string description_gr_stability_fit=description_general_stability_fit+detKey;
    gr_stability_charge_fit_singletube[detkey] = new TGraphErrors();
    gr_stability_charge_fit_singletube[detkey]->SetName(name_gr_stability_fit.c_str());
    gr_stability_charge_fit_singletube[detkey]->SetTitle(description_gr_stability_fit.c_str());

    string name_general_stability_mean="gr_stability_charge_mean_";
    string description_general_stability_mean="Stability of mean charge for detkey ";
    string name_gr_stability_mean=name_general_stability_mean+detKey;
    string description_gr_stability_mean=description_general_stability_mean+detKey;
    gr_stability_charge_mean_singletube[detkey] = new TGraphErrors();
    gr_stability_charge_mean_singletube[detkey]->SetName(name_gr_stability_mean.c_str());
    gr_stability_charge_mean_singletube[detkey]->SetTitle(description_gr_stability_mean.c_str());

  }

  //
  //read in information from 100 last runs to produce stability/time evolution plots (if they do not exist, just skip them...)
  //

  const int n_entries=100;
  double run_numbers[n_entries];
  double entries_charge_mean_fit[n_entries];
  double entries_charge_rms_fit[n_entries];
  double entries_charge_mean[n_entries];
  double entries_charge_rms[n_entries];
  double entries_time_mean_fit[n_entries];
  double entries_time_rms_fit[n_entries];
  double entries_time_mean[n_entries];
  double entries_time_rms[n_entries];
  double entries_time_dev_fit[n_entries];
  double entries_time_dev_mean[n_entries];
  double entries_time_expected[n_entries];
  int run_nentries[n_entries];
  bool file_exists[n_entries];
  std::map<unsigned long,std::vector<double>> single_entries_charge_mean_fit;
  std::map<unsigned long,std::vector<double>> single_entries_charge_rms_fit;
  std::map<unsigned long,std::vector<double>> single_entries_charge_mean;
  std::map<unsigned long,std::vector<double>> single_entries_charge_rms;
  std::map<unsigned long,std::vector<double>> single_entries_time_mean_fit;
  std::map<unsigned long,std::vector<double>> single_entries_time_rms_fit;
  std::map<unsigned long,std::vector<double>> single_entries_time_mean;
  std::map<unsigned long,std::vector<double>> single_entries_time_rms;
  std::map<unsigned long,std::vector<double>> single_entries_time_dev_fit;
  std::map<unsigned long,std::vector<double>> single_entries_time_dev_mean;
  std::map<unsigned long,std::vector<int>> single_run_nentries;
  std::vector<int> vec_run_numbers;

  int runnumber_temp;

  //
  //get stability values for time and charge from txt-files
  //

  for (int i_run=n_entries-1;i_run>=0;i_run--){

    std::stringstream ss_runnumber;
    runnumber_temp = runnumber-i_run;
    ss_runnumber<<runnumber_temp;
    std::string filename_mid = ss_runnumber.str();
    std::string filename_temp = outputdir+outputfile+filename_pre+filename_mid+filename_post;
    if (verbose > 1) std::cout <<"Stability: Filename: "<<filename_temp<<", file exists: "<<does_file_exist(filename_temp.c_str())<<std::endl;
    double dummy_temp;
    if (does_file_exist(filename_temp.c_str())){
      ifstream file_temp(filename_temp.c_str());
      file_exists[n_entries-1-i_run]=true;
      unsigned long detkey;
      double temp_charge_mean_fit;
      double temp_charge_rms_fit;
      double temp_charge_mean;
      double temp_charge_rms;
      double temp_time_mean_fit;
      double temp_time_rms_fit;
      double temp_time_mean;
      double temp_time_rms;
      double temp_time_dev_fit;
      double temp_time_dev;
      double temp_time_expected;
      int temp_run_nentries;
      vec_run_numbers.push_back(runnumber-i_run);

      std::cout <<"Removing first line of calibrations txt-file"<<std::endl;
      std::string line;
      file_temp>>line>>line>>line>>line>>line>>line>>line>>line>>line>>line>>line>>line>>line;
      std::cout <<"First line is: "<<line<<std::endl;
      

      for (int i_pmt=0;i_pmt<n_tank_pmts+1;i_pmt++){


        if (i_pmt<n_tank_pmts) {

        file_temp>>detkey>>temp_charge_mean_fit>>temp_charge_rms_fit>>temp_charge_mean>>temp_charge_rms>>temp_time_mean_fit>>temp_time_rms_fit>>temp_time_mean>>temp_time_rms>>temp_time_expected>>temp_time_dev_fit>>temp_time_dev>>temp_run_nentries;  
        
        std::cout <<"detkey = "<<detkey<<", temp_charge_mean_fit = "<<temp_charge_mean_fit<<", temp_charge_mean = "<<temp_charge_mean<<std::endl;

        if (single_entries_charge_mean_fit.find(detkey) != single_entries_charge_mean_fit.end()){

          single_entries_charge_mean_fit.at(detkey).push_back(temp_charge_mean_fit);
          single_entries_charge_rms_fit.at(detkey).push_back(temp_charge_rms_fit);
          single_entries_charge_mean.at(detkey).push_back(temp_charge_mean);
          single_entries_charge_rms.at(detkey).push_back(temp_charge_rms);
          single_entries_time_mean_fit.at(detkey).push_back(temp_time_mean_fit);
          single_entries_time_rms_fit.at(detkey).push_back(temp_time_rms_fit);
          single_entries_time_mean.at(detkey).push_back(temp_time_mean);
          single_entries_time_rms.at(detkey).push_back(temp_time_rms);
          single_entries_time_dev_fit.at(detkey).push_back(temp_time_dev_fit);
          single_entries_time_dev_mean.at(detkey).push_back(temp_time_dev);
          single_run_nentries.at(detkey).push_back(temp_run_nentries);

        } else {

	  std::vector<double> temp_vector;
          temp_vector.push_back(temp_charge_mean_fit);
          single_entries_charge_mean_fit.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_charge_rms_fit);
          single_entries_charge_rms_fit.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_charge_mean);
          single_entries_charge_mean.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_charge_rms);
          single_entries_charge_rms.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_mean_fit);
          single_entries_time_mean_fit.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_rms_fit);
          single_entries_time_rms_fit.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_mean);
          single_entries_time_mean.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_rms);
          single_entries_time_rms.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_dev_fit);
          single_entries_time_dev_fit.emplace(detkey,temp_vector);
          temp_vector.clear();
          temp_vector.push_back(temp_time_dev);
          single_entries_time_dev_mean.emplace(detkey,temp_vector);
          std::vector<int> temp_vector_int;
          temp_vector_int.push_back(temp_run_nentries);
          single_run_nentries.emplace(detkey,temp_vector_int);

       }

        }
        else file_temp>>dummy_temp>>entries_charge_mean_fit[n_entries-1-i_run]>>entries_charge_rms_fit[n_entries-1-i_run]>>entries_charge_mean[n_entries-1-i_run]>>entries_charge_rms[n_entries-1-i_run]>>entries_time_mean_fit[n_entries-1-i_run]>>entries_time_rms_fit[n_entries-1-i_run]>>entries_time_mean[n_entries-1-i_run]>>entries_time_rms[n_entries-1-i_run]>>entries_time_expected[n_entries-1-i_run]>>entries_time_dev_fit[n_entries-1-i_run]>>entries_time_dev_mean[n_entries-1-i_run]>>run_nentries[n_entries-1-i_run];
      }
      file_temp.close();
    }
    else {
      file_exists[n_entries-1-i_run]=false;
    }
  }

  //
  //fill TGraphErrors with stability values for time and charge
  //

  int i_point_graph=0;
  for (int i_run=n_entries-1;i_run>=0;i_run--){
    if (file_exists[n_entries-1-i_run]){
      if (verbose > 1) std::cout <<"Stability, Run "<<n_entries-1-i_run<<": Fit Mean Charge: "<<entries_charge_mean_fit[n_entries-1-i_run]<<", Fit RMS Charge: "<<entries_charge_rms_fit[n_entries-1-i_run]<<", Mean Charge: "<<entries_charge_mean[n_entries-1-i_run]<<", RMS Charge: "<<entries_charge_rms[n_entries-1-i_run]<<", Fit Mean Time: "<<entries_time_mean_fit[n_entries-1-i_run]<<"Fit RMS time: "<<entries_time_rms_fit[n_entries-1-i_run]<<", Mean Time: "<<entries_time_mean[n_entries-1-i_run]<<", RMS Time: "<<entries_time_rms[n_entries-1-i_run]<<", Fit Time deviation: "<<entries_time_dev_fit[n_entries-1-i_run]<<", Time Deviation: "<<entries_time_dev_mean[n_entries-1-i_run]<<std::endl;
      gr_stability_charge_fit->SetPoint(i_point_graph,runnumber-i_run,entries_charge_mean_fit[n_entries-1-i_run]);
      double rms_mean = (run_nentries[n_entries-1-i_run] > 0)? entries_charge_rms_fit[n_entries-1-i_run]/sqrt(run_nentries[n_entries-1-i_run]) : entries_charge_rms_fit[n_entries-1-i_run];
      gr_stability_charge_fit->SetPointError(i_point_graph,0.,rms_mean);
      gr_stability_charge_mean->SetPoint(i_point_graph,runnumber-i_run,entries_charge_mean[n_entries-1-i_run]);
      rms_mean = (run_nentries[n_entries-1-i_run] > 0)? entries_charge_rms[n_entries-1-i_run]/sqrt(run_nentries[n_entries-1-i_run]) : entries_charge_rms[n_entries-1-i_run];
      gr_stability_charge_mean->SetPointError(i_point_graph,0.,rms_mean);
      gr_stability_time_fit->SetPoint(i_point_graph,runnumber-i_run,entries_time_mean_fit[n_entries-1-i_run]);
      rms_mean = (run_nentries[n_entries-1-i_run] > 0)? entries_time_rms_fit[n_entries-1-i_run]/sqrt(run_nentries[n_entries-1-i_run]) : entries_time_rms_fit[n_entries-1-i_run];
      gr_stability_time_fit->SetPointError(i_point_graph,0.,rms_mean);
      gr_stability_time_mean->SetPoint(i_point_graph,runnumber-i_run,entries_time_mean[n_entries-1-i_run]);
      rms_mean = (run_nentries[n_entries-1-i_run] > 0)? entries_time_rms[n_entries-1-i_run]/sqrt(run_nentries[n_entries-1-i_run]) : entries_time_rms[n_entries-1-i_run];
      gr_stability_time_mean->SetPointError(i_point_graph,0.,rms_mean);
      gr_stability_time_dev_fit->SetPoint(i_point_graph,runnumber-i_run,entries_time_dev_fit[n_entries-1-i_run]);
      gr_stability_time_dev_fit->SetPointError(i_point_graph,0.,0.);   // set the error 0 for now
      gr_stability_time_dev_mean->SetPoint(i_point_graph,runnumber-i_run,entries_time_dev_mean[n_entries-1-i_run]);
      gr_stability_time_dev_mean->SetPointError(i_point_graph,0.,0.);  // set the error 0 for now
      i_point_graph++;
    }
  }

  gr_stability_charge_fit->SetMarkerStyle(22);
  gr_stability_charge_fit->SetMarkerSize(1.0);
  gr_stability_charge_fit->SetLineColor(1);
  gr_stability_charge_fit->SetLineWidth(2);
  gr_stability_charge_fit->GetXaxis()->SetTitle("Run Number");
  gr_stability_charge_fit->GetYaxis()->SetTitle("charge [p.e.]");
  gr_stability_time_fit->SetMarkerStyle(22);
  gr_stability_time_fit->SetMarkerSize(1.0);
  gr_stability_time_fit->SetLineColor(1);
  gr_stability_time_fit->SetLineWidth(2);
  gr_stability_time_fit->GetXaxis()->SetTitle("Run Number");
  gr_stability_time_fit->GetYaxis()->SetTitle("hit time [ns]");
  gr_stability_time_dev_fit->SetMarkerStyle(22);
  gr_stability_time_dev_fit->SetMarkerSize(1.0);
  gr_stability_time_dev_fit->SetLineColor(1);
  gr_stability_time_dev_fit->SetLineWidth(2);
  gr_stability_time_dev_fit->GetXaxis()->SetTitle("Run Number");
  gr_stability_time_dev_fit->GetYaxis()->SetTitle("hit time deviation [ns]");
  gr_stability_charge_fit->Write("gr_stability_charge_fit");
  gr_stability_time_fit->Write("gr_stability_time_fit");
  gr_stability_time_dev_fit->Write("gr_stability_time_dev_fit");
  gr_stability_charge_mean->SetMarkerStyle(22);
  gr_stability_charge_mean->SetMarkerSize(1.0);
  gr_stability_charge_mean->SetLineColor(1);
  gr_stability_charge_mean->SetLineWidth(2);
  gr_stability_charge_mean->GetXaxis()->SetTitle("Run Number");
  gr_stability_charge_mean->GetYaxis()->SetTitle("charge [p.e.]");
  gr_stability_time_mean->SetMarkerStyle(22);
  gr_stability_time_mean->SetMarkerSize(1.0);
  gr_stability_time_mean->SetLineColor(1);
  gr_stability_time_mean->SetLineWidth(2);
  gr_stability_time_mean->GetXaxis()->SetTitle("Run Number");
  gr_stability_time_mean->GetYaxis()->SetTitle("hit time [ns]");
  gr_stability_time_dev_mean->SetMarkerStyle(22);
  gr_stability_time_dev_mean->SetMarkerSize(1.0);
  gr_stability_time_dev_mean->SetLineColor(1);
  gr_stability_time_dev_mean->SetLineWidth(2);
  gr_stability_time_dev_mean->GetXaxis()->SetTitle("Run Number");
  gr_stability_time_dev_mean->GetYaxis()->SetTitle("hit time deviation [ns]");
  gr_stability_charge_mean->Write("gr_stability_charge_mean");
  gr_stability_time_mean->Write("gr_stability_time_mean");
  gr_stability_time_dev_mean->Write("gr_stability_time_dev_mean");


  //
  // Fill tube-wise TGraphErrors with stability values for time and charge
  //    

  for (int i_pmt = 0; i_pmt < n_tank_pmts; i_pmt++){

    unsigned long detkey = pmt_detkeys[i_pmt];
    std::cout <<"looping through i_pmt = "<<i_pmt<<", detkey = "<<detkey<<std::endl;
    for (unsigned int i_run = 0; i_run < single_entries_charge_mean_fit.at(detkey).size(); i_run++){
      std::cout <<"i_run = "<<i_run<<", vec_run_numbers.at(i_run) = "<<vec_run_numbers.at(i_run)<<std::endl;
      gr_stability_charge_fit_singletube[detkey]->SetPoint(i_run,vec_run_numbers.at(i_run),single_entries_charge_mean_fit.at(detkey).at(i_run));
      double rms_fit = single_entries_charge_rms_fit.at(detkey).at(i_run);
      if (single_run_nentries.at(detkey).at(i_run) > 0) rms_fit/=sqrt(single_run_nentries.at(detkey).at(i_run));
      gr_stability_charge_fit_singletube[detkey]->SetPointError(i_run,0.,rms_fit);
      std::cout <<"Setting charge_fit = "<<single_entries_charge_mean_fit.at(detkey).at(i_run)<<", charge fit rms = "<<rms_fit<<std::endl;
      gr_stability_charge_mean_singletube[detkey]->SetPoint(i_run,vec_run_numbers.at(i_run),single_entries_charge_mean.at(detkey).at(i_run));
      double rms = single_entries_charge_rms.at(detkey).at(i_run);
      if (single_run_nentries.at(detkey).at(i_run) > 0) rms /=sqrt(single_run_nentries.at(detkey).at(i_run));
      gr_stability_charge_mean_singletube[detkey]->SetPointError(i_run,0.,rms);
      std::cout <<"Setting charge_mean = "<<single_entries_charge_mean.at(detkey).at(i_run)<<", charge rms = "<<rms<<std::endl;
    }

    gr_stability_charge_fit_singletube[detkey]->SetMarkerStyle(22);
    gr_stability_charge_fit_singletube[detkey]->SetMarkerSize(1.0);
    gr_stability_charge_fit_singletube[detkey]->SetLineColor(1);
    gr_stability_charge_fit_singletube[detkey]->SetLineWidth(2);
    gr_stability_charge_fit_singletube[detkey]->GetXaxis()->SetTitle("Run Number");
    gr_stability_charge_fit_singletube[detkey]->GetYaxis()->SetTitle("fit mean charge");

    gr_stability_charge_mean_singletube[detkey]->SetMarkerStyle(22);
    gr_stability_charge_mean_singletube[detkey]->SetMarkerSize(1.0);
    gr_stability_charge_mean_singletube[detkey]->SetLineColor(1);
    gr_stability_charge_mean_singletube[detkey]->SetLineWidth(2);
    gr_stability_charge_mean_singletube[detkey]->GetXaxis()->SetTitle("Run Number");
    gr_stability_charge_mean_singletube[detkey]->GetYaxis()->SetTitle("mean charge");

    std::stringstream ss_detkey;
    ss_detkey << detkey;
    std::stringstream ss_gr_name_fit, ss_gr_name_mean;
    ss_gr_name_fit << "gr_stability_charge_fit_"<<ss_detkey.str();
    ss_gr_name_mean << "gr_stability_charge_mean_"<<ss_detkey.str();

    std::cout <<"Name gr_name_fit: "<<ss_gr_name_fit.str().c_str()<<std::endl;
    gr_stability_charge_mean_singletube[detkey]->Write(ss_gr_name_mean.str().c_str());
    gr_stability_charge_fit_singletube[detkey]->Write(ss_gr_name_fit.str().c_str());

  }

  //
  //Fill new entries into histogram
  //

  std::map<int,double>::iterator it_time = bad_time.begin();
  std::map<int,double>::iterator it_charge = bad_charge.begin();

  std::cout <<"///////////////////////////////////////////////////////////////////"<<std::endl;
  std::cout <<"-------------------------------------------------------------------"<<std::endl;
  std::cout <<"-------------PMTs with deviating properties:-----------------------"<<std::endl;
  std::cout <<"-------------------------------------------------------------------"<<std::endl;  
  std::cout <<"///////////////////////////////////////////////////////////////////"<<std::endl;

  for (int i_tube=0;i_tube<n_tank_pmts;i_tube++){
    unsigned long detkey = pmt_detkeys[i_tube];
    if (fabs(time_mean_fit[detkey]-expected_time[detkey])>tolerance_time) {
      std::cout <<"Abnormally high time deviation for detkey "<<detkey<<": "<<time_mean_fit[detkey]-expected_time[detkey]<<std::endl;
      bad_time.insert(it_time,std::pair<int,double>(detkey,time_mean_fit[detkey]-expected_time[detkey]));    
    }
    if (charge_mean_fit[detkey]<(1.0-tolerance_charge) || charge_mean_fit[detkey]>(1.0+tolerance_charge)){
      std::cout <<"Abnormal charge for detkey "<<detkey<<": "<<charge_mean_fit[detkey]<<std::endl;
      bad_charge.insert(it_charge,std::pair<int,double>(detkey,charge_mean_fit[detkey])); 
    }
  }

  //----------------------------------------------------------------------------
  //---------------Save summary root plots (charge / time fit values)-----------
  //----------------------------------------------------------------------------

  canvas_overview_QT_fit = new TCanvas("canvas_overview_QT_fit","Charge/Time Overview histograms (Fit values)",900,600);
  canvas_overview_QT_fit->Divide(3,2);
  canvas_overview_QT_fit->cd(1);
  hist_charge_2D_y_phi_fit->Draw("colz");
  std::map<int,double>::iterator it_charge_read;                                  //draw red boxes around the misaligned PMTs (in charge as well as in time)
  if (verbose > 0) std::cout << "PMTs with bad charge: "<<std::endl;
  for (it_charge_read=bad_charge.begin();it_charge_read!=bad_charge.end();it_charge_read++){
    int pmt = it_charge_read->first;
    int bin_charge_x = hist_charge_2D_y_phi_fit->GetXaxis()->FindBin(phi_PMT[pmt]/TMath::Pi()*180.);
    int bin_charge_y = hist_charge_2D_y_phi_fit->GetYaxis()->FindBin(y_PMT[pmt]);
    if (verbose > 0) std::cout << "PMT: "<<pmt<<", bin: "<<bin_charge_x<<", "<<bin_charge_y<<std::endl;
    draw_red_box(hist_charge_2D_y_phi_fit,bin_charge_x, bin_charge_y,verbose);
  }
  canvas_overview_QT_fit->cd(2);
  hist_time_2D_y_phi_fit->Draw("colz");
  std::map<int,double>::iterator it_time_read;
  if (verbose > 0) std::cout << "PMTs with bad timing: "<<std::endl;
  for (it_time_read=bad_time.begin();it_time_read!=bad_time.end();it_time_read++){
    int pmt = it_time_read->first;
    //if (verbose > 0) std::cout <<"PMT: "<<pmt<<", phi PMT: "<<phi_PMT[pmt]<<", y PMT: "<<y_PMT[pmt]<<std::endl;
    int bin_time_x = hist_time_2D_y_phi_fit->GetXaxis()->FindBin(phi_PMT[pmt]/TMath::Pi()*180.);
    int bin_time_y = hist_time_2D_y_phi_fit->GetYaxis()->FindBin(y_PMT[pmt]);
    if (verbose > 0) std::cout << "PMT: "<<pmt<<", bin: "<<bin_time_x<<", "<<bin_time_y<<std::endl;
    draw_red_box(hist_time_2D_y_phi_fit,bin_time_x, bin_time_y, verbose);
  }
  canvas_overview_QT_fit->cd(3);
  hist_time_2D_y_phi_dev_fit->Draw("colz");  
  canvas_overview_QT_fit->cd(4);
  hist_charge_fit->Draw();
  canvas_overview_QT_fit->cd(5);
  hist_time_2D_y_phi_fit->Draw("colz");
  canvas_overview_QT_fit->cd(6);
  hist_time_dev_fit->Draw();
  canvas_overview_QT_fit->Modified();
  canvas_overview_QT_fit->Update();
  canvas_overview_QT_fit->Write();

  //----------------------------------------------------------------------------
  //---------------Save summary root plots v2 (charge / time mean values) ------
  //----------------------------------------------------------------------------

  canvas_overview_QT_mean = new TCanvas("canvas_overview_QT_mean","Charge/Time Overview histograms (Mean)",900,600);
  canvas_overview_QT_mean->Divide(3,2);
  canvas_overview_QT_mean->cd(1);
  hist_charge_2D_y_phi_mean->Draw("colz");
  canvas_overview_QT_mean->cd(2);
  hist_time_2D_y_phi_mean->Draw("colz");
  canvas_overview_QT_mean->cd(3);
  hist_time_2D_y_phi_dev_mean->Draw("colz");
  canvas_overview_QT_mean->cd(4);
  hist_charge_mean->Draw();
  canvas_overview_QT_mean->cd(5);
  hist_time_mean->Draw();
  canvas_overview_QT_mean->cd(6);
  hist_time_dev_mean->Draw();
  canvas_overview_QT_mean->Modified();
  canvas_overview_QT_mean->Update();
  canvas_overview_QT_mean->Write();

  //----------------------------------------------------------------------------
  //---------------Save summary root plots v3 (RecoADC Values 2D)----------------
  //----------------------------------------------------------------------------

  canvas_overview_recoadc_2D = new TCanvas("canvas_overview_recoadc_2D","Overview RecoADC values (2D)",900,1000);
  canvas_overview_recoadc_2D->Divide(4,2);
  canvas_overview_recoadc_2D->cd(1);
  hist_starttime_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(2);
  hist_peaktime_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(3);
  hist_baseline_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(4);
  hist_sigmabaseline_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(5);
  hist_rawamplitude_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(6);
  hist_amplitude_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(7);
  hist_rawarea_2D_y_phi_mean->Draw("colz");
  canvas_overview_recoadc_2D->cd(8);
  hist_detkey_2D_y_phi->Draw("colz text");
  canvas_overview_recoadc_2D->Modified();
  canvas_overview_recoadc_2D->Update();
  canvas_overview_recoadc_2D->Write();
  
  //----------------------------------------------------------------------------
  //---------------Save summary root plots v4 (RecoADC Values 1D)---------------
  //----------------------------------------------------------------------------

  canvas_overview_recoadc_1D = new TCanvas("canvas_overview_recoadc_1D","Overview RecoADC values (1D)",900,1000);
  canvas_overview_recoadc_1D->Divide(4,2);
  canvas_overview_recoadc_1D->cd(1);
  hist_detkey_starttime_mean->Draw();
  canvas_overview_recoadc_1D->cd(2);
  hist_detkey_peaktime_mean->Draw();
  canvas_overview_recoadc_1D->cd(3);
  hist_detkey_baseline_mean->Draw();
  canvas_overview_recoadc_1D->cd(4);
  hist_detkey_sigmabaseline_mean->Draw();
  canvas_overview_recoadc_1D->cd(5);
  hist_detkey_rawamplitude_mean->Draw();
  canvas_overview_recoadc_1D->cd(6);
  hist_detkey_amplitude_mean->Draw();
  canvas_overview_recoadc_1D->cd(7);
  hist_detkey_rawarea_mean->Draw();
  canvas_overview_recoadc_1D->Modified();
  canvas_overview_recoadc_1D->Update();
  canvas_overview_recoadc_1D->Write();

  //----------------------------------------------------------------------------
  //---------------Save summary root plots v5 (Stability Q/T mean)--------------
  //----------------------------------------------------------------------------

  canvas_stability_mean = new TCanvas("canvas_stability_mean","Stability canvas (Mean)",900,1000);
  canvas_stability_mean->Divide(2,1);
  canvas_stability_mean->cd(1);
  gr_stability_charge_mean->Draw("AEL");
  canvas_stability_mean->cd(2);
  gr_stability_time_mean->Draw("AEL");
  canvas_stability_mean->Modified();
  canvas_stability_mean->Update();
  canvas_stability_mean->Write();

  //----------------------------------------------------------------------------
  //---------------Save summary root plots v5 (Stability Q/T Fit)--------------
  //----------------------------------------------------------------------------

  canvas_stability_fit = new TCanvas("canvas_stability_fit","Stability canvas (Fit)",900,1000);
  canvas_stability_fit->Divide(2,1);
  canvas_stability_fit->cd(1);
  gr_stability_charge_fit->Draw("AEL");
  canvas_stability_fit->cd(2);
  gr_stability_time_fit->Draw("AEL");
  canvas_stability_fit->Modified();
  canvas_stability_fit->Update();
  canvas_stability_fit->Write();

  file_out->Close();

  //----------------------------------------------------------------------------
  //---------------Prompt user inputs (to avoid abortion of canvas)-------------
  //----------------------------------------------------------------------------

  int valid_calibration;
  if (verbose > 0) {
    std::cout << "//////////////////////////////////////////////////////////////"<<std::endl;
    std::cout << "--------------------------------------------------------------"<<std::endl;
    std::cout << "Calibration run finished. Summary: "<<endl;
    std::cout << "Number of misaligned PMTs (time): "<<bad_time.size()<<endl;
    std::cout << "Number of misaligned PMTs (charge): "<<bad_charge.size()<<endl;
    std::cout << "--------------------------------------------------------------"<<std::endl;
    std::cout << "//////////////////////////////////////////////////////////////"<<std::endl;

    if (use_tapplication){
    std::cout <<"Press enter to close the TApplication window..."<<std::endl;
    std::cin >>valid_calibration;
    app_stability->Terminate();
    }
  }

  help_file->Close();

  //----------------------------------------------------------------------------
  //------------------------Deleting all objects -------------------------------
  //----------------------------------------------------------------------------

  //Histograms already deleted by closing the files ealier

  if (verbose > 0) std::cout <<"Deleting TBoxes..."<<std::endl;
  for (int i_box = 0; i_box < (int) vector_tbox.size();i_box++){
    delete vector_tbox.at(i_box);
  }
  
  if (verbose > 0) std::cout <<"Deleting TF1s..."<<std::endl;
  for (int i_tf = 0; i_tf < (int) vector_tf1.size(); i_tf++){
    delete vector_tf1.at(i_tf);
  } 

  if (verbose > 0) std::cout <<"Deleting canvas, application, files, geom, mchits"<<std::endl;
  delete canvas_overview_QT_fit;
  delete canvas_overview_QT_mean;
  delete canvas_overview_recoadc_1D;
  delete canvas_overview_recoadc_2D;
  delete canvas_stability_mean;
  delete canvas_stability_fit;
  if (use_tapplication) delete app_stability;
  delete help_file;
  delete file_out;
  delete geom;
  delete MCHits;

  //only for debugging memory leaks, otherwise comment out
  /*std::cout <<"List of Objects (end of finalise): "<<std::endl;
  gObjectTable->Print();*/


  return true;

}

//----------------------------------------------------------------------------
//---------------------Helper functions --------------------------------------
//----------------------------------------------------------------------------

bool TankCalibrationDiffuser::does_file_exist(const char *fileName){
  std::ifstream infile(fileName);
  return infile.good();
}

void TankCalibrationDiffuser::draw_red_box(TH2 *hist, int bin_x, int bin_y, int verbosity){

  Int_t bin_global = hist->GetBin(bin_x, bin_y,0);

  double size_x=2.0;
  double size_y=0.05;
  TBox *b = new TBox(hist->GetXaxis()->GetBinLowEdge(bin_x)-size_x,hist->GetYaxis()->GetBinLowEdge(bin_y)-size_y,hist->GetXaxis()->GetBinWidth(bin_x)+hist->GetXaxis()->GetBinLowEdge(bin_x)+size_x, hist->GetYaxis()->GetBinWidth(bin_y)+hist->GetYaxis()->GetBinLowEdge(bin_y)+size_y);
  
  if (verbosity > 2) {
  	std::cout <<"Draw Red Box: lower x: "<<hist->GetXaxis()->GetBinLowEdge(bin_x)-0.1<<std::endl;
  	std::cout <<"Draw Red Box: lower y: "<<hist->GetYaxis()->GetBinLowEdge(bin_y)-0.1<<std::endl;
  	std::cout <<"Draw Red Box: upper x: "<<hist->GetXaxis()->GetBinWidth(bin_x)+hist->GetXaxis()->GetBinLowEdge(bin_x)+0.1<<std::endl;
  	std::cout <<"Draw Red Box: upper y: "<<hist->GetXaxis()->GetBinWidth(bin_y)+hist->GetXaxis()->GetBinLowEdge(bin_y)+0.1<<std::endl;
  }

  b->SetFillStyle(0);
  b->SetLineColor(2);
  b->SetLineWidth(2);
  b->Draw();

  vector_tbox.push_back(b);
}
