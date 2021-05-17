#include <string>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"

// for checking for memory leaks only. Also in $HOME/.rootrc set `Root.ObjectStat` to 1
//#include "TObjectTable.h"

#include <thread>
#include <chrono>

// ToolAnalysis includes
#include "PhaseIIADCCalibrator.h"

PhaseIIADCCalibrator::PhaseIIADCCalibrator() : Tool() {}

bool PhaseIIADCCalibrator::Initialise(std::string config_filename, DataModel& data)
{
  // Load configuration file variables
  if ( !config_filename.empty() ) m_variables.Initialise(config_filename);

  // Assign a transient data pointer
  m_data = &data;
  
  std::cout<<"PhaseIIADCCalibrator Tool: Initializing"<<std::endl;
  
  // get config variables: general
  m_variables.Get("verbosity", verbosity);
 
  // get config variables: Event Building mode?
  m_variables.Get("EventBuilding",eventbuilding_mode);
  if (eventbuilding_mode != false && eventbuilding_mode != true){
    Log("PhaseIIADCCalibrator: Event Building mode not recognized. Default to false",v_warning,verbosity);
    eventbuilding_mode = false;
  }

  //Set defaults in case config file has no entries
  adc_window_db = "none";
  make_led_waveforms = false;
  BEType = "ze3ra";
 
  // algorithm selection
  get_ok = m_variables.Get("BaselineEstimationType", BEType);
  if(BEType != "ze3ra" && BEType != "rootfit" && BEType != "simple" && BEType != "ze3ra_multi"){
    Log("PhaseIIADCCalibrator Tool: Baseline estimation type not recognized!  Default to ze3ra", v_warning, verbosity);
     BEType = "ze3ra";
  }
  Log("PhaseIIADCCalibrator Tool: Configured to use "+BEType+" baseline subtraction method", v_message, verbosity);
  
  //Set defaults in case config file has no entries
  p_critical = 0.01;
  num_sub_waveforms = 6;
  num_baseline_samples = (BEType == "rootfit") ? 980 : 5;
  
  // get baseline variables
  m_variables.Get("NumBaselineSamples", num_baseline_samples);

  // get ze3ra variables 
  m_variables.Get("PCritical", p_critical);
  m_variables.Get("NumSubWaveforms", num_sub_waveforms);

  // Get the Auxiliary channel types; identifies which channels are SiPM channels
  m_data->CStore.Get("AuxChannelNumToTypeMap",AuxChannelNumToTypeMap);
  
  // get ze3ra_multi variables, set defaults
  baseline_rep_samples = 1000;
  baseline_unc_tolerance = 5;
  m_variables.Get("SamplesPerBaselineEstimate", baseline_rep_samples);
  m_variables.Get("BaselineUncertaintyTolerance", baseline_unc_tolerance);

  // get LED waveform-making variables
  m_variables.Get("MakeCalLEDWaveforms",make_led_waveforms);
  m_variables.Get("WindowIntegrationDB", adc_window_db); 
  
  // get ROOT fitting variables
  if(BEType == "rootfit"){
    m_variables.Get("drawBaselineRootFit",draw_baseline_fit);
    if(not get_ok) draw_baseline_fit=false;
    m_variables.Get("BaselineFitStartSample",baseline_start_sample);
    if(not get_ok) baseline_start_sample=0;
    get_ok = m_variables.Get("BaselineFitOrder",baseline_fit_order);
    if(not get_ok) baseline_fit_order = 1; // default to linear fit
    get_ok = m_variables.Get("RedoFitWithoutOutliers",redo_fit_without_outliers);
    if(not get_ok) redo_fit_without_outliers = false; // whether to redo the fit after removing outliers
    get_ok = m_variables.Get("RefitThresholdAdcCounts", refit_threshold); // but only if wfrm range exceeds this
    if(not get_ok) refit_threshold = 5.; // something suitable
    refit_threshold*=ADC_TO_VOLT;
  }
  
  // get or make the ROOT TApplication to show the fit, debug only
  if(draw_baseline_fit){
    int myargc=0;
    intptr_t tapp_ptr=0;
    get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
    if(not get_ok){
      Log("PhaseIIADCCalibrator Tool: making global TApplication",v_debug,verbosity);
      rootTApp = new TApplication("rootTApp",&myargc,0);
      tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
      m_data->CStore.Set("RootTApplication",tapp_ptr);
    } else {
      Log("PhaseIIADCCalibrator Tool: Retrieving global TApplication",v_debug,verbosity);
      rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
    }
    // register ourself as a user
    int tapplicationusers;
    get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
    if(not get_ok) tapplicationusers=1;
    else tapplicationusers++;
    m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
  }


  if(adc_window_db != "none") channel_window_map = this->load_window_map(adc_window_db);

  m_data->CStore.Set("NumBaselineSamples",num_baseline_samples);


  m_data->CStore.Set("NewCalibratedData",false);
  m_data->CStore.Get("TankPMTCrateSpaceToChannelNumMap",TankPMTCrateSpaceToChannelNumMap);
  m_data->CStore.Get("AuxCrateSpaceToChannelNumMap",AuxCrateSpaceToChannelNumMap);

  if (eventbuilding_mode){

    FinishedRawWaveforms = new std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>>;
    FinishedRawWaveformsAux = new std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>>;
    FinishedCalibratedWaveforms = new std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>>;
    FinishedCalibratedWaveformsAux = new std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>>;
    FinishedCalibratedLEDADCData = new std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>>;
    FinishedRawLEDADCData = new std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>>;
    FinishedRawAcqSize = new std::map<uint64_t, std::map<unsigned long, std::vector<int>>>;
    m_data->CStore.Set("FinishedRawWaveforms",FinishedRawWaveforms);
    m_data->CStore.Set("FinishedRawWaveformsAux",FinishedRawWaveformsAux);
    m_data->CStore.Set("FinishedRawAcqSize",FinishedRawAcqSize);
    m_data->CStore.Set("FinishedCalibratedWaveforms",FinishedCalibratedWaveforms);
    m_data->CStore.Set("FinishedCalibratedWaveformsAux",FinishedCalibratedWaveformsAux);
    m_data->CStore.Set("FinishedCalibratedLEDADCData",FinishedCalibratedLEDADCData);
    m_data->CStore.Set("FinishedRawLEDADCData",FinishedRawLEDADCData);
    m_data->CStore.Set("NewCalibratedData",true);

    bool get_executes;
    get_executes = m_variables.Get("ExecutesPerBuild",ExecutesPerBuild);
    if (!get_executes) ExecutesPerBuild = 10;
    ExecuteCount = 0;
  }
  return true;
}

bool PhaseIIADCCalibrator::Execute() {
  
  Log("PhaseIIADCCalibrator Tool: Executing", v_message, verbosity);


  //ANNIEEvent mode
  if (eventbuilding_mode == false){

  // Get a pointer to the ANNIEEvent Store
  auto* annie_event = m_data->Stores["ANNIEEvent"];

  if (!annie_event) {
    Log("Error: The PhaseIIADCCalibrator tool could not find the ANNIEEvent Store", 0,
      verbosity);
    return false;
  }

  // Load the map containing the ADC raw waveform data
  std::map<unsigned long, std::vector<Waveform<unsigned short> > >
    raw_waveform_map;
  // Load the map containing the ADC raw waveform data
  std::map<unsigned long, std::vector<Waveform<unsigned short> > >
    raw_auxwaveform_map;

  bool got_raw_data = annie_event->Get("RawADCData", raw_waveform_map);
  bool got_rawaux_data = annie_event->Get("RawADCAuxData", raw_auxwaveform_map);

  // Check for problems
  if ( !got_raw_data ) {
    Log("Error: The PhaseIIADCCalibrator tool could not find the RawADCData entry", 0,
      verbosity);
    return false;
  }
  else if ( raw_waveform_map.empty() ) {
    Log("Error: The PhaseIIADCCalibrator tool found an empty RawADCData entry", 0,
      verbosity);
    return false;
  }



  // Build the calibrated waveforms
  std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > >
    calibrated_waveform_map;
  // Build the calibrated waveforms
  std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > >
    calibrated_auxwaveform_map;

  // Load the map containing the ADC raw waveform data
  std::map<unsigned long, std::vector<Waveform<unsigned short> > >
    raw_led_waveform_map;

  std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > >
    calibrated_led_waveform_map;

  //Calibrate raw detector waveforms
  for (const auto& temp_pair : raw_waveform_map) {
    const auto& channel_key = temp_pair.first;
    //Default running: raw_waveforms only has one entry.  If we go to a
    //hefty-mode style of running though, this could have multiple minibuffers
    const auto& raw_waveforms = temp_pair.second;
    Log("Making calibrated waveforms for ADC channel " +
      std::to_string(channel_key), 3, verbosity);

    if(BEType == "ze3ra"){
      calibrated_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra(raw_waveforms);
    } else if(BEType == "ze3ra_multi"){
      calibrated_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(raw_waveforms);
    } else if(BEType == "rootfit"){
      calibrated_waveform_map[channel_key] = make_calibrated_waveforms_rootfit(raw_waveforms);
    } else if (BEType == "simple"){
      calibrated_waveform_map[channel_key] = make_calibrated_waveforms_simple(raw_waveforms);
    }

    if(make_led_waveforms){
      Log("Also making LED window waveforms for ADC channel " +
        std::to_string(channel_key), 3, verbosity);
      std::vector<Waveform<unsigned short>> LEDWaveforms;
      this->make_raw_led_waveforms(channel_key,raw_waveforms,LEDWaveforms);
      raw_led_waveform_map.emplace(channel_key,LEDWaveforms);
      if(BEType == "ze3ra"){
        calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra(LEDWaveforms);
      } else if(BEType == "ze3ra_multi"){
        calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(LEDWaveforms);
      } else if(BEType == "rootfit"){
        calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_rootfit(LEDWaveforms);
      } else if(BEType == "simple"){
        calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_simple(LEDWaveforms);
      }
    }
  }
  
  //Calibrate the SIPM waveforms
  for (const auto& temp_pair : raw_auxwaveform_map) {
    const auto& channel_key = temp_pair.first;
    Log("Channel key for Aux channel is " +
      std::to_string(channel_key), 3, verbosity);
    //For now, only calibrate the SiPM waveforms
    Log("Type for Aux channel is " +
      AuxChannelNumToTypeMap->at(channel_key), 3, verbosity);
    if(AuxChannelNumToTypeMap->at(channel_key) != "SiPM1" && 
       AuxChannelNumToTypeMap->at(channel_key) != "SiPM2") continue; 
    //Default running: raw_waveforms only has one entry.  If we go to a
    //hefty-mode style of running though, this could have multiple minibuffers
    const auto& raw_auxwaveforms = temp_pair.second;

    Log("Making calibrated waveforms for Auxiliary channel " +
      std::to_string(channel_key), 3, verbosity);

    if(BEType == "ze3ra"){
      calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_ze3ra(raw_auxwaveforms);
    } else if(BEType == "ze3ra_multi"){
      calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(raw_auxwaveforms);
    } else if(BEType == "rootfit"){
      calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_rootfit(raw_auxwaveforms);
    } else if (BEType == "simple"){
      calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_simple(raw_auxwaveforms);
    }
  }

  Log("PhaseIIADCCalibrator Tool: Setting CalibratedADCData",v_debug,verbosity);
  annie_event->Set("CalibratedADCData", calibrated_waveform_map);
  annie_event->Set("CalibratedADCAuxData", calibrated_auxwaveform_map);
  if(make_led_waveforms){
    std::cout <<"Setting LEDADCData"<<std::endl;
    annie_event->Set("CalibratedLEDADCData", calibrated_led_waveform_map);
    annie_event->Set("RawLEDADCData", raw_led_waveform_map);
  }
  //std::cout <<"Set CalibratedADCData"<<std::endl;
  } else {
    //=========================
    //Event Building mode (new)
    //=========================
  
    ExecuteCount += 1;
    bool toolchain_stopping=false;
    m_data->vars.Get("StopLoop",toolchain_stopping);
    if (ExecuteCount < ExecutesPerBuild && !(toolchain_stopping)) return true;

    bool new_tank_events = false;
    //m_data->CStore.Get("NewTankEvents",new_tank_events);    
    m_data->CStore.Get("NewTankPMTDataAvailable",new_tank_events);
    if (!new_tank_events) return true;	//Don't do anything if no new tank data is available

    //m_data->CStore.Get("FinishedTankEvents",FinishedTankEvents);
    m_data->CStore.Get("InProgressTankEvents",InProgressTankEvents);

    bool new_data = false;
    std::vector<uint64_t> RawTimestampsToDelete;

    //Loop over FinishedTankEvents, fill FinishedCalibratedWaveforms
    //for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : *FinishedTankEvents){
    for(std::pair<uint64_t,std::map<std::vector<int>, std::vector<uint16_t>>> apair : *InProgressTankEvents){
      uint64_t PMTCounterTime = apair.first;
      //std::cout <<"PMTCounterTime: "<<PMTCounterTime<<", waveform size: "<<apair.second.size()<<std::endl;

      //Check if the timestamp was already processed
      //if (FinishedRawWaveforms->count(PMTCounterTime) != 0) continue;
      new_data = true;
      RawTimestampsToDelete.push_back(PMTCounterTime);
      std::map<std::vector<int>, std::vector<uint16_t>> aWaveMap = apair.second;
      std::map<unsigned long, std::vector<Waveform<uint16_t>> > RawADCData;
      std::map<unsigned long, std::vector<Waveform<uint16_t>> > RawADCAuxData;
//      if (FinishedRawWaveforms->count(PMTCounterTime)>0) RawADCData = FinishedRawWaveforms->at(PMTCounterTime);
//      if (FinishedRawWaveformsAux->count(PMTCounterTime)>0) RawADCAuxData = FinishedRawWaveformsAux->at(PMTCounterTime);
      for(std::pair<std::vector<int>, std::vector<uint16_t>> apair : aWaveMap){
        int CardID = apair.first.at(0);
        int ChannelID = apair.first.at(1);
        int CrateNum=-1;
        int SlotNum=-1;
        this->CardIDToElectronicsSpace(CardID, CrateNum, SlotNum);
        std::vector<uint16_t> TheWaveform = apair.second;
        Waveform<uint16_t> TheWave(PMTCounterTime, TheWaveform);
        //Placing waveform in a vector in case we want a hefty-mode minibuffer storage eventually
        std::vector<Waveform<uint16_t>> WaveVec{TheWave};
  
        std::vector<int> CrateSpace{CrateNum,SlotNum,ChannelID};
        unsigned long ChannelKey;
        if(TankPMTCrateSpaceToChannelNumMap.count(CrateSpace)>0){
          ChannelKey = TankPMTCrateSpaceToChannelNumMap.at(CrateSpace);
          RawADCData.emplace(ChannelKey,WaveVec);
        }
        else if (AuxCrateSpaceToChannelNumMap.count(CrateSpace)>0){
          ChannelKey = AuxCrateSpaceToChannelNumMap.at(CrateSpace);
          RawADCAuxData.emplace(ChannelKey,WaveVec);
        } else{
          Log("PhaseIIADCCalibrator:: Cannot find channel key for crate space entry: ",v_error, verbosity);
          Log("PhaseIIADCCalibrator::CrateNum "+to_string(CrateNum),v_error, verbosity);
          Log("PhaseIIADCCalibrator::SlotNum "+to_string(SlotNum),v_error, verbosity);
          Log("PhaseIIADCCalibrator::ChannelID "+to_string(ChannelID),v_error, verbosity);
          Log("PhaseIIADCCalibrator::Passing over the wave; PMT DATA LOST",v_error, verbosity);
          continue;
        }
      }

      //Build raw waveform acquisition size
      std::map<unsigned long, std::vector<int>> waveform_acq_size;
//      if (FinishedRawAcqSize->count(PMTCounterTime)>0) waveform_acq_size = FinishedRawAcqSize->at(PMTCounterTime);
      for (auto& temp_pair : RawADCData) {
        const auto& achannel_key = temp_pair.first;
        auto& araw_waveforms = temp_pair.second;
        for (unsigned int i=0; i< araw_waveforms.size(); i++){
          auto samples = araw_waveforms.at(i).GetSamples();
          int size_sample = samples->size();
          if (waveform_acq_size.count(achannel_key) > 0) waveform_acq_size.at(achannel_key).push_back(size_sample);
          else {
            std::vector<int> temp_size{size_sample};
            waveform_acq_size.emplace(achannel_key,temp_size);
          }
        }
      }

      // Build the calibrated waveforms
      std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > > calibrated_waveform_map;
      // Build the calibrated waveforms
      std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > > calibrated_auxwaveform_map;


      // Load the map containing the ADC raw waveform data
      std::map<unsigned long, std::vector<Waveform<unsigned short> > > raw_led_waveform_map;
      std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > > calibrated_led_waveform_map;
/*
      if (FinishedCalibratedWaveforms->count(PMTCounterTime)>0) calibrated_waveform_map = FinishedCalibratedWaveforms->at(PMTCounterTime);
      if (FinishedCalibratedWaveformsAux->count(PMTCounterTime)>0) calibrated_auxwaveform_map = FinishedCalibratedWaveformsAux->at(PMTCounterTime);

      if (make_led_waveforms){
        if (FinishedCalibratedLEDADCData->count(PMTCounterTime)>0) calibrated_led_waveform_map = FinishedCalibratedLEDADCData->at(PMTCounterTime);
        if (FinishedRawLEDADCData->count(PMTCounterTime)>0) raw_led_waveform_map = FinishedRawLEDADCData->at(PMTCounterTime);
      }*/

      //Calibrate raw detector waveforms
      for (const auto& temp_pair : RawADCData) {
        const auto& channel_key = temp_pair.first;
        //Default running: raw_waveforms only has one entry.  If we go to a
        //hefty-mode style of running though, this could have multiple minibuffers
        const auto& raw_waveforms = temp_pair.second;
        //std::cout <<"Calibrate raw detector waveforms for channel "<<channel_key<<std::endl;
        Log("Making calibrated waveforms for ADC channel " +  std::to_string(channel_key), 3, verbosity);

        if(BEType == "ze3ra"){
          calibrated_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra(raw_waveforms);
        } else if(BEType == "ze3ra_multi"){
          calibrated_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(raw_waveforms);
        } else if(BEType == "rootfit"){
          calibrated_waveform_map[channel_key] = make_calibrated_waveforms_rootfit(raw_waveforms);
        } else if (BEType == "simple"){
          calibrated_waveform_map[channel_key] = make_calibrated_waveforms_simple(raw_waveforms);
        }
        //std::cout <<"temp_calibrated_waveform_map.size(): "<<temp_calibrated_waveform_map[channel_key].size()<<", raw size: "<<raw_waveforms.size()<<std::endl;

/*
	if (raw_waveforms.at(0).Samples().size() != temp_calibrated_waveform_map[channel_key].at(0).Samples().size()){
		std::cout <<"PhaseIIADCCalibrator: raw_waveforms.size(): "<<raw_waveforms.at(0).Samples().size()<<", temp_calibrated_waveforms.size(): "<<temp_calibrated_waveform_map[channel_key].at(0).Samples().size()<<std::endl;
	}

        for (std::pair<unsigned long, std::vector<CalibratedADCWaveform<double>>> apair : temp_calibrated_waveform_map){
          unsigned long chkey = apair.first;
          std::vector<CalibratedADCWaveform<double>> calib_wave = apair.second;
          for (int i_wave=0; i_wave  < calib_wave.size(); i_wave++){
            if (calibrated_waveform_map.count(chkey)>0) calibrated_waveform_map[chkey].push_back(calib_wave.at(i_wave));
            else calibrated_waveform_map.emplace(chkey,std::vector<CalibratedADCWaveform<double>>{calib_wave.at(i_wave)});
          }
        }*/

        if(make_led_waveforms){
          Log("Also making LED window waveforms for ADC channel " + std::to_string(channel_key), 3, verbosity);
          std::vector<Waveform<unsigned short>> LEDWaveforms;
          this->make_raw_led_waveforms(channel_key,raw_waveforms,LEDWaveforms);
          raw_led_waveform_map.emplace(channel_key,LEDWaveforms);
          if(BEType == "ze3ra"){
            calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra(LEDWaveforms);
          } else if(BEType == "ze3ra_multi"){
            calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(LEDWaveforms);
          } else if(BEType == "rootfit"){
            calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_rootfit(LEDWaveforms);
          } else if(BEType == "simple"){
            calibrated_led_waveform_map[channel_key] = make_calibrated_waveforms_simple(LEDWaveforms);
          }
/*
        for (std::pair<unsigned long, std::vector<CalibratedADCWaveform<double>>> apair : temp_calibrated_led_waveform_map){
          unsigned long chkey = apair.first;
          std::vector<CalibratedADCWaveform<double>> calib_wave = apair.second;
          for (int i_wave=0; i_wave  < calib_wave.size(); i_wave++){
            if (calibrated_led_waveform_map.count(chkey)>0) calibrated_led_waveform_map[chkey].push_back(calib_wave.at(i_wave));
            else calibrated_led_waveform_map.emplace(chkey,std::vector<CalibratedADCWaveform<double>>{calib_wave.at(i_wave)});
          }
        }*/
        }
      }
  
      //Calibrate the SIPM waveforms
      for (const auto& temp_pair : RawADCAuxData) {
        const auto& channel_key = temp_pair.first;
        Log("Channel key for Aux channel is " + std::to_string(channel_key), 3, verbosity);
        //For now, only calibrate the SiPM waveforms
        Log("Type for Aux channel is " + AuxChannelNumToTypeMap->at(channel_key), 3, verbosity);
        if(AuxChannelNumToTypeMap->at(channel_key) != "SiPM1" && AuxChannelNumToTypeMap->at(channel_key) != "SiPM2") continue; 
        //Default running: raw_waveforms only has one entry.  If we go to a
        //hefty-mode style of running though, this could have multiple minibuffers
        const auto& raw_auxwaveforms = temp_pair.second;

        Log("Making calibrated waveforms for Auxiliary channel " + std::to_string(channel_key), 3, verbosity);

        if(BEType == "ze3ra"){
          calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_ze3ra(raw_auxwaveforms);
        } else if(BEType == "ze3ra_multi"){
          calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_ze3ra_multi(raw_auxwaveforms);
        } else if(BEType == "rootfit"){
          calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_rootfit(raw_auxwaveforms);
        } else if (BEType == "simple"){
          calibrated_auxwaveform_map[channel_key] = make_calibrated_waveforms_simple(raw_auxwaveforms);
        }
        /*for (std::pair<unsigned long, std::vector<CalibratedADCWaveform<double>>> apair : temp_calibrated_auxwaveform_map){
          unsigned long chkey = apair.first;
          std::vector<CalibratedADCWaveform<double>> calib_wave = apair.second;
          for (int i_wave=0; i_wave  < calib_wave.size(); i_wave++){
            if (calibrated_auxwaveform_map.count(chkey)>0) calibrated_auxwaveform_map[chkey].push_back(calib_wave.at(i_wave));
            else calibrated_auxwaveform_map.emplace(chkey,std::vector<CalibratedADCWaveform<double>>{calib_wave.at(i_wave)});
          }
        }*/

      }

      Log("PhaseIIADCCalibrator Tool: Setting CalibratedADCData",v_debug,verbosity);
      //std::cout <<"RawADCData.size(): "<<RawADCData.size()<<", Calibrated waveforms size: "<<calibrated_waveform_map.size()<<std::endl;

      FinishedRawAcqSize->emplace(PMTCounterTime,waveform_acq_size);
      FinishedRawWaveforms->emplace(PMTCounterTime,RawADCData);
      FinishedRawWaveformsAux->emplace(PMTCounterTime,RawADCAuxData);
      FinishedCalibratedWaveforms->emplace(PMTCounterTime,calibrated_waveform_map);
      FinishedCalibratedWaveformsAux->emplace(PMTCounterTime,calibrated_auxwaveform_map);
      if(make_led_waveforms){
        std::cout <<"Setting LEDADCData"<<std::endl;
        FinishedCalibratedLEDADCData->emplace(PMTCounterTime,calibrated_led_waveform_map);
        FinishedRawLEDADCData->emplace(PMTCounterTime,raw_led_waveform_map);
      }
      //std::cout <<"Set CalibratedADCData"<<std::endl;
    }

    for (int i_del=0; i_del < (int) RawTimestampsToDelete.size(); i_del++){
      //FinishedTankEvents->erase(RawTimestampsToDelete.at(i_del));
      InProgressTankEvents->erase(RawTimestampsToDelete.at(i_del));
      //std::map<std::vector<int>, std::vector<uint16_t>> aRawWaveform = InProgressTankEvents->at(RawTimestampsToDelete.at(i_del));
      //aRawWaveform.clear();
    }

    if (new_data){
      m_data->CStore.Set("FinishedRawWaveforms",FinishedRawWaveforms);
      m_data->CStore.Set("FinishedRawWaveformsAux",FinishedRawWaveformsAux);
      m_data->CStore.Set("FinishedRawAcqSize",FinishedRawAcqSize);
      m_data->CStore.Set("FinishedCalibratedWaveforms",FinishedCalibratedWaveforms);
      m_data->CStore.Set("FinishedCalibratedWaveformsAux",FinishedCalibratedWaveformsAux);
      if (make_led_waveforms){
        m_data->CStore.Set("FinishedCalibratedLEDADCData",FinishedCalibratedLEDADCData);
        m_data->CStore.Set("FinishedRawLEDADCData",FinishedRawLEDADCData);
      }
      m_data->CStore.Set("NewCalibratedData",true);
      }

      ExecuteCount = 0;
    }


  return true;
}


bool PhaseIIADCCalibrator::Finalise() {
  
  if(BEType == "rootfit"){
    Log("PhaseIIADCCalibrator Tool: Cleaning up ROOT fitting objects",v_message,verbosity);
    //std::cout<<"dumping gObjecTable:"<<std::endl;
    //gObjectTable->Print();
    
    if(calibrated_waveform_tgraph) delete calibrated_waveform_tgraph; calibrated_waveform_tgraph=nullptr;
    if(calibrated_waveform_fit) delete calibrated_waveform_fit; calibrated_waveform_fit=nullptr;
    if(gROOT->FindObject("raw_datapoint_hist")!=nullptr) delete raw_datapoint_hist; raw_datapoint_hist=nullptr;
    
    if(draw_baseline_fit || (drawcount>0)){
      if(gROOT->FindObject("baselineFitCanvas")!=nullptr) delete baselineFitCanvas; baselineFitCanvas=nullptr;
      
      int tapplicationusers=0;
      get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
      if(not get_ok || tapplicationusers==1){
        if(rootTApp){
          Log("PhaseIIADCCalibrator Tool: Deleting global TApplication",v_message,verbosity);
          delete rootTApp;
          rootTApp=nullptr;
        }
      } else if(tapplicationusers>1){
        m_data->CStore.Set("RootTApplicationUsers",tapplicationusers-1);
      }
    }
  }
  
  return true;
}

void PhaseIIADCCalibrator::ze3ra_baseline(
  const  Waveform<unsigned short> raw_data,
  double& baseline, double& sigma_baseline, size_t num_baseline_samples,size_t starting_sample)
{

  // Signal ADC means, variances, and F-distribution probability values
  // ("P") for the first num_baseline_samples from each minibuffer
  // (in Hefty mode) or from each sub-minibuffer (in non-Hefty mode)
  std::vector<double> means;
  std::vector<double> variances;
  std::vector<double> Ps;

  // Using the Phase I non-hefty algorithm. Split the early part of the waveform 
  // into sub-minibuffers and compute the mean and variance of each one.
  const auto& data = raw_data.Samples();
  for (size_t sub_mb = 0u; sub_mb < num_sub_waveforms; ++sub_mb) {
    std::vector<unsigned short> sub_mb_data(
      data.cbegin()+ starting_sample + sub_mb * num_baseline_samples,
      data.cbegin() + starting_sample + (1u + sub_mb) * num_baseline_samples);

    double mean, var;
    ComputeMeanAndVariance(sub_mb_data, mean, var, num_baseline_samples);

    means.push_back(mean);
    variances.push_back(var);
  }

  // Compute probabilities for the F-distribution test for each waveform chunk
  for (size_t j = 0; j < variances.size() - 1; ++j) {
    double sigma2_j = variances.at(j);
    double sigma2_jp1 = variances.at(j + 1);
    double F;
    if (sigma2_j > sigma2_jp1) F = sigma2_j / sigma2_jp1;
    else F = sigma2_jp1 / sigma2_j;

    double nu = (num_baseline_samples - 1) / 2.;
    double P = annie_math::Regularized_Beta_Function(1. / (1. + F), nu, nu);

    // Two-tailed hypothesis test (we need to exclude unusually small values
    // as well as unusually large ones). The tails have equal sizes, so we
    // may use symmetry and simply multiply our earlier result by 2.
    P *= 2.;

    // I've never seen this problem (the numerical values for the regularized
    // beta function that I've checked all fall within [0,1]), but the book
    // Numerical Recipes includes this check in a similar block of code,
    // so I'll add it just in case.
    if (P > 1.) P = 2. - P;

    Ps.push_back(P);
  }

  // Compute the mean and standard deviation of the baseline signal
  // for this RawChannel using the mean and standard deviation from
  // each minibuffer whose F-distribution probability falls below
  // the critical value.
  baseline = 0.;
  sigma_baseline = 0.;
  double variance_baseline = 0.;
  size_t num_passing = 0;
  for (size_t k = 0; k < Ps.size(); ++k) {
    if (Ps.at(k) > p_critical) {
      ++num_passing;
      baseline += means.at(k);
      variance_baseline += variances.at(k);
    }
  }

  if (num_passing > 1) {
    baseline /= num_passing;

    variance_baseline *= static_cast<double>(num_baseline_samples - 1)
      / (num_passing*num_baseline_samples - 1);
    // Now that we've combined the sample variances correctly, take the
    // square root to get the standard deviation
    sigma_baseline = std::sqrt( variance_baseline );
  }
  else if (num_passing == 1) {
    // We only have one set of sample statistics, so all we need to
    // do is take the square root of the variance to get the standard
    // deviation.
    sigma_baseline = std::sqrt( variance_baseline );
  }
  else {
    // If none of the minibuffers passed the F-distribution test,
    // choose the one closest to passing (i.e., the one with the largest
    // P-value) and adopt its baseline statistics. For a sufficiently large
    // number of minibuffers (e.g., 40), such a situation should be very rare.
    // TODO: consider changing this approach
    auto max_iter = std::max_element(Ps.cbegin(), Ps.cend());
    int max_index = std::distance(Ps.cbegin(), max_iter);

    baseline = means.at(max_index);
    sigma_baseline = std::sqrt( variances.at(max_index) );
  }

  std::string mb_temp_string = "minibuffer";

  if (verbosity >= 4) {
    for ( size_t x = 0; x < Ps.size(); ++x ) {
      Log("  " + mb_temp_string + " " + std::to_string(x) + ", mean = "
        + std::to_string(means.at(x)) + ", var = "
        + std::to_string(variances.at(x)) + ", p-value = "
        + std::to_string(Ps.at(x)), 4, verbosity);
    }
  }

  Log(std::to_string(num_passing) + " " + mb_temp_string + " pairs passed the"
    " F-test", 3, verbosity);
  Log("Baseline estimate: " + std::to_string(baseline) + " ± "
    + std::to_string(sigma_baseline) + " ADC counts", 3, verbosity);

}

// version based on the ze3bra algorithm; assumes a DC offset is sufficient
std::vector< CalibratedADCWaveform<double> >
PhaseIIADCCalibrator::make_calibrated_waveforms_simple(
  const std::vector< Waveform<unsigned short> >& raw_waveforms)
{

  // Determine the baseline for the set of raw waveforms (assumed to all
  // come from the same readout for the same channel)
  std::vector< CalibratedADCWaveform<double> > calibrated_waveforms;
  for (const auto& raw_waveform : raw_waveforms) {
    double baseline, sigma_baseline;
    std::vector<double> cal_data;
    const std::vector<unsigned short>& raw_data = raw_waveform.Samples();
    ComputeMeanAndVariance(raw_data, baseline, sigma_baseline, num_baseline_samples);
    for (const auto& sample : raw_data) {
      cal_data.push_back((static_cast<double>(sample) - baseline)
        * ADC_TO_VOLT);
    }
    calibrated_waveforms.emplace_back(raw_waveform.GetStartTime(),
      cal_data, baseline, sigma_baseline);
  }
  return calibrated_waveforms;
}

// version based on the ze3bra algorithm; assumes a DC offset is sufficient
std::vector< CalibratedADCWaveform<double> >
PhaseIIADCCalibrator::make_calibrated_waveforms_ze3ra(
  const std::vector< Waveform<unsigned short> >& raw_waveforms)
{

  // Determine the baseline for the set of raw waveforms (assumed to all
  // come from the same readout for the same channel)
  std::vector< CalibratedADCWaveform<double> > calibrated_waveforms;
  for (const auto& raw_waveform : raw_waveforms) {
    double baseline, sigma_baseline;
    ze3ra_baseline(raw_waveform, baseline, sigma_baseline,
      num_baseline_samples, 0);
    std::vector<double> cal_data;
    const std::vector<unsigned short>& raw_data = raw_waveform.Samples();
    for (const auto& sample : raw_data) {
      cal_data.push_back((static_cast<double>(sample) - baseline)
        * ADC_TO_VOLT);
    }
  
    calibrated_waveforms.emplace_back(raw_waveform.GetStartTime(),
      cal_data, baseline, sigma_baseline);
  }
  return calibrated_waveforms;
}

// version based on the ze3bra algorithm; assumes a DC offset is sufficient
std::vector< CalibratedADCWaveform<double> >
PhaseIIADCCalibrator::make_calibrated_waveforms_ze3ra_multi(
  const std::vector< Waveform<unsigned short> >& raw_waveforms)
{

  // Determine the baseline for the set of raw waveforms (assumed to all
  // come from the same readout for the same channel)
  std::vector< CalibratedADCWaveform<double> > calibrated_waveforms;
  for (const auto& raw_waveform : raw_waveforms) {
    std::vector<uint16_t> baselines;
    std::vector<size_t> RepresentationRegion;
    double first_baseline, first_sigma_baseline;
    double baseline, sigma_baseline;
    const size_t nsamples = raw_waveform.Samples().size();
    for(size_t starting_sample = 0; starting_sample < nsamples; starting_sample += baseline_rep_samples){
      double baseline, sigma_baseline;
      ze3ra_baseline(raw_waveform, baseline, sigma_baseline,
        num_baseline_samples,starting_sample);
      if(sigma_baseline<baseline_unc_tolerance){
        RepresentationRegion.push_back(starting_sample + baseline_rep_samples);
        baselines.push_back(baseline);
      } else {
        if(verbosity>4) std::cout << "BASELINE UNCERTAINTY BEYOND SET THRESHOLD.  IGNORING SAMPLE" << std::endl;
      }
    }

    // If NO baselines within tolerance found, just go with the first
    if(baselines.size() == 0){
      if(verbosity>4) std::cout << "NO BASLINE FOUND WITHIN TOLERANCE.  USING FIRST AS BEST ESTIMATE" << std::endl;
      RepresentationRegion.push_back(baseline_rep_samples);
      baselines.push_back(first_baseline);
    }
    std::vector<double> cal_data;
    const std::vector<unsigned short>& raw_data = raw_waveform.Samples();
    for (const auto& asample: raw_data){
      for(int j = 0; j<RepresentationRegion.size(); j++){
        if(asample < RepresentationRegion.at(j)){
          cal_data.push_back((static_cast<double>(asample) - baselines.at(j))
            * ADC_TO_VOLT);
          break;
        } else if (asample >= RepresentationRegion.back()){
          cal_data.push_back((static_cast<double>(asample) - baselines.back())
            * ADC_TO_VOLT);
          break;
        }
      }
    }
    double bl_estimates_mean, bl_estimates_var;
    ComputeMeanAndVariance(baselines, bl_estimates_mean, bl_estimates_var);
    calibrated_waveforms.emplace_back(raw_waveform.GetStartTime(),
      cal_data, bl_estimates_mean, bl_estimates_var);
  }
  return calibrated_waveforms;
}



// version based on a polynomial fit done via ROOT
std::vector< CalibratedADCWaveform<double> >
PhaseIIADCCalibrator::make_calibrated_waveforms_rootfit(
  const std::vector< Waveform<unsigned short> >& raw_waveforms){
  Log("PhaseIIADCCalibrator Tool: Doing ROOT based baseline subtraction", v_debug, verbosity);
  std::vector< CalibratedADCWaveform<double> > calibrated_waveforms;
  
  // Apparently the input waveforms given to us represent all the minibuffers for one channel,
  // obtained in one particular ADC readout. But I don't expect that matters much to us.
  
  // create a TGraph if we don't have one
  //if(gROOT->FindObject("calibrated_waveform_tgraph")==nullptr){  << never finds it.
  if(calibrated_waveform_tgraph==nullptr){
    // we need to know how many points the tgraph will hold
    num_waveform_points = raw_waveforms.front().Samples().size();
    Log("PhaseIIADCCalibrator Tool: Making TGraph with " + to_string(num_waveform_points)
        +" data points", v_debug, verbosity);
    //Log("PhaseIIADCCalibrator Tool: Making new Graph",v_error,verbosity);
    calibrated_waveform_tgraph = new TGraph(num_waveform_points);
    // honestly i have no idea how to stop ROOT interfering with object deletion after drawing
    // see https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html
    // and https://root.cern.ch/root/roottalk/roottalk01/2060.html
    if(not calibrated_waveform_tgraph){
      Log("PhaseIIADCCalibrator Tool: ERROR making TGraph!", v_error, verbosity);
      return calibrated_waveforms;
    }
    calibrated_waveform_tgraph->SetName("calibrated_waveform_tgraph");
    gROOT->Add(calibrated_waveform_tgraph); // make gROOT know about it - doesn't work??
    calibrated_waveform_tgraph->ResetBit(kCanDelete);
    //calibrated_waveform_tgraph->GetHistogram()->SetName("calibrated_waveform_tgraph_histo"); // doesn't help
    //gROOT->Add(calibrated_waveform_tgraph->GetHistogram());
  }
  
  // create a TF1 if we don't have one
  if(calibrated_waveform_fit==nullptr){
    // take waveform range to fit either from minibuffer length or user if less and >0
    if((num_baseline_samples<=0) || (num_baseline_samples>(num_waveform_points-baseline_start_sample)))
        num_baseline_samples = num_waveform_points;
    if(baseline_start_sample<0) baseline_start_sample = 0;
    Log("PhaseIIADCCalibrator Tool: Making fit function of type pol"+to_string(baseline_fit_order)
         + " to fit waveform samples " + to_string(baseline_start_sample) + " to " 
         + to_string(baseline_start_sample+num_waveform_points), v_debug, verbosity);
    //Log("PhaseIIADCCalibrator Tool: Making new Function",v_error,verbosity);
    calibrated_waveform_fit = new TF1("calibrated_waveform_fit",TString::Format("pol%d",baseline_fit_order),
      baseline_start_sample,num_baseline_samples);
    if(not calibrated_waveform_fit){
      Log("PhaseIIADCCalibrator Tool: ERROR making root TF1!", v_error, verbosity);
      return calibrated_waveforms;
    }
  }
  
  // Loop over raw waveforms
  Log("PhaseIIADCCalibrator Tool: Looping over "+to_string(raw_waveforms.size()) + " raw waveforms",
      v_debug, verbosity);
  for (const auto& raw_waveform : raw_waveforms){
    
    // retrieve the raw samples
    const std::vector<unsigned short>& raw_data = raw_waveform.Samples();
    
    // update our TGraph's datapoints with the new datapoints
    Log("PhaseIIADCCalibrator Tool: Setting TGraph datapoints", v_debug, verbosity);
    for(int samplei=0; samplei<raw_data.size(); ++samplei){
      calibrated_waveform_tgraph->SetPoint(samplei,samplei,raw_data.at(samplei));
    }
    
    // fit the graph
    Log("PhaseIIADCCalibrator Tool: Fitting the baseline", v_debug, verbosity);
    TFitResultPtr fit_result = calibrated_waveform_tgraph->Fit(calibrated_waveform_fit,"RCFSQ"); // F
    // R to use range of TF1
    // F option uses minuit fitter for polN... better?
    // C skips calculation of chi2 to save time
    // S prevents conversion of the TFitResultPtr to a plain integer, which can only be used to check status
    // in this case we need to check the status of the fit manually: just by converting the pointer
    // Q prevents Minuit spewing all it's internals
    
    bool fit_succeeded = ((static_cast<Int_t>(fit_result))==0);  // successful fit is 0
    std::vector<double> fitpars(baseline_fit_order+1);
    double baseline=0;
    if(not fit_succeeded){
      Log("PhaseIIADCCalibrator Tool: polynomial fit of baseline failed!",v_warning,verbosity);
    } else {
      Log("PhaseIIADCCalibrator Tool: polynomial fit of baseline succeeded, noting parameters",v_debug,verbosity);
      // make a note of the current fit parameters, in case we re-do the fit later
      for(uint orderi=0; orderi<(baseline_fit_order+1); ++orderi){
        fitpars.at(orderi) = fit_result->Value(orderi);
      }
      logmessage="PhaseIIADCCalibrator Tool: Baseline fit success: fit function was: ";
      for(uint orderi=0; orderi<(baseline_fit_order+1); ++orderi){
        logmessage+= to_string(fitpars.at(orderi))+"*x^"+to_string(orderi);
        if(orderi<baseline_fit_order) logmessage+=" + ";
      }
      Log(logmessage, v_debug, verbosity);
      baseline = fitpars.at(0); // DC offset. FIXME this doesn't fully capture the correction applied
    }
    
    // in either case, draw the data and fit result
    if(draw_baseline_fit){
      Log("PhaseIIADCCalibrator Tool: Drawing initial baseline fit",v_message,verbosity);
      if(gROOT->FindObject("baselineFitCanvas")==nullptr){
        Log("PhaseIIADCCalibrator Tool: Constructing canvas for drawing fit", v_debug, verbosity);
        //Log("PhaseIIADCCalibrator Tool: Making new Canvas",v_error,verbosity);
        baselineFitCanvas = new TCanvas("baselineFitCanvas");
      } else {
        baselineFitCanvas->cd();
      }
      calibrated_waveform_tgraph->Draw("ALP");
      calibrated_waveform_tgraph->GetHistogram()->SetDirectory(0); // canvas should not take ownership
      calibrated_waveform_tgraph->ResetBit(kCanDelete);
      baselineFitCanvas->Modified();
      baselineFitCanvas->Update();
      gSystem->ProcessEvents();
      Log("PhaseIIADCCalibrator Tool: Sleeping while waiting for user to close canvas",v_debug,verbosity);
      while(gROOT->FindObject("baselineFitCanvas")!=nullptr){
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        gSystem->ProcessEvents();
      }
      drawcount++; if(drawcount>10) draw_baseline_fit=false;
    }
    
    // convert samples from ADC count to volts, and subtract the baseline, if fit succeeded
    std::vector<double> cal_data(raw_data.size());
    // keep track of the max and min
    double cal_data_min=std::numeric_limits<double>::max();
    double cal_data_max=std::numeric_limits<double>::min();
    // loop over samples
    Log("PhaseIIADCCalibrator Tool: Subtracting baseline and converting to ADC counts", v_debug, verbosity);
    for(uint samplei=0; samplei<raw_data.size(); ++samplei){
      const unsigned short& sample = raw_data.at(samplei);
      if(gROOT->FindObject("calibrated_waveform_fit")==nullptr){
        Log("PhaseIIADCCalibrator Tool: ERROR! Could not find TF1!", v_error, verbosity);
        return calibrated_waveforms;
      }
      //std::cout<<"Evaluating baseline fit"<<std::endl;
      double baseline_val = (fit_succeeded) ? calibrated_waveform_fit->Eval(samplei) : 0;
      //std::cout<<"Estimated baseline was "<<baseline_val<<std::endl;
      double cal_val = (static_cast<double>(sample) - baseline_val)*ADC_TO_VOLT;
      //std::cout<<"Baseline subtracted val was "<<cal_val<<std::endl;
      cal_data.at(samplei)=cal_val;
      
      if(cal_val<cal_data_min) cal_data_min = cal_val;
      if(cal_val>cal_data_max) cal_data_max = cal_val;
    }
    
    if(draw_baseline_fit){
      Log("PhaseIIADCCalibrator Tool: Drawing baseline subtracted fit",v_message,verbosity);
      
      // update our TGraph's datapoints with the new datapoints
      Log("PhaseIIADCCalibrator Tool: Setting TGraph datapoints", v_debug, verbosity);
      for(int samplei=0; samplei<cal_data.size(); ++samplei){
        calibrated_waveform_tgraph->SetPoint(samplei,samplei,cal_data.at(samplei));
      }
      
      // Redo the fit to guide the eye. Not used.
      TFitResultPtr fit_result = calibrated_waveform_tgraph->Fit(calibrated_waveform_fit,"RCFSQ"); // F
      
      if(gROOT->FindObject("baselineFitCanvas")==nullptr){
        Log("PhaseIIADCCalibrator Tool: Constructing canvas for drawing fit", v_debug, verbosity);
        //Log("PhaseIIADCCalibrator Tool: Making new Canvas",v_error,verbosity);
        baselineFitCanvas = new TCanvas("baselineFitCanvas");
      } else {
        baselineFitCanvas->cd();
      }
      
      calibrated_waveform_tgraph->Draw("ALP");
      calibrated_waveform_tgraph->GetHistogram()->SetDirectory(0); // canvas should not take ownership
      calibrated_waveform_tgraph->ResetBit(kCanDelete);
      baselineFitCanvas->Modified();
      baselineFitCanvas->Update();
      gSystem->ProcessEvents();
      Log("PhaseIIADCCalibrator Tool: Sleeping while waiting for user to close canvas",v_debug,verbosity);
      while(gROOT->FindObject("baselineFitCanvas")!=nullptr){
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        gSystem->ProcessEvents();
      }
    }
    
    // do we need to calculate sigma_baseline? Does anyone use it? TODO
    double sigma_baseline = 0;
    
    // if we're worried about our baseline fit being skewed by the presence of large pulses,
    // remove outliers and do the baseline fit again. We can skip this step if the range is sufficiently
    // small that we know there aren't any pulses to skew the data
    double cal_data_range = cal_data_max - cal_data_min;
    logmessage = "PhaseIIADCCalibrator Tool: calibrated wfrm range was " +to_string(cal_data_range);
    if(redo_fit_without_outliers){
      if(cal_data_range>refit_threshold){
        logmessage+= " which will invoke outlier removal and refit if initial fit succeeded";
      } else {
        logmessage+= " which will not invoke outlier removal and refit";
      }
    }
    Log(logmessage, v_debug, verbosity);
    if((redo_fit_without_outliers) && (cal_data_range>refit_threshold) && (fit_succeeded)){
      Log("PhaseIIADCCalibrator Tool: Removing outliers", v_debug, verbosity);
      
      // find and remove outliers
      std::vector<double> non_outlier_points(cal_data);
      // We'll use interquartile range as a definition of data excluding outliers.
      // First we need to find the interquartile range. Do it the lazy way: with ROOT.
      Log("PhaseIIADCCalibrator Tool: Getting histogram for quantiles", v_debug, verbosity);
      if(gROOT->FindObject("raw_datapoint_hist")==nullptr){
        Log("PhaseIIADCCalibrator Tool: Making histogram for quantile measurement", v_debug, verbosity);
        //Log("PhaseIIADCCalibrator Tool: Making new Histogram",v_error,verbosity);
        raw_datapoint_hist = new TH1D("raw_datapoint_hist","Raw Data Histogram",200,cal_data_min,cal_data_max);
        if(not raw_datapoint_hist){
          Log("PhaseIIADCCalibrator Tool: ERROR! Failed to making raw datapoint histogram!", v_error, verbosity);
          return calibrated_waveforms;
        }
      } else {
        Log("PhaseIIADCCalibrator Tool: Resetting raw data histogram", v_debug, verbosity);
        raw_datapoint_hist->Reset(); raw_datapoint_hist->SetBins(200, cal_data_min, cal_data_max);
      }
      for(double& cal_sample : cal_data){ raw_datapoint_hist->Fill(cal_sample); }
      
      // get the quantile thresholds. We'll take the interquartile range,
      // although we could take arbitrary thresholds
      std::vector<double> threshold_probabilities{0.00,0.95}; // graphs are inverted; only clip top (pulses)
      std::vector<double> threshold_values(threshold_probabilities.size());
      // get the quantile thresholds. Note GetQuantiles asks for it's input arrays backwards...
      Log("PhaseIIADCCalibrator Tool: Getting Quantiles", v_debug, verbosity);
      raw_datapoint_hist->GetQuantiles(threshold_probabilities.size(),
          threshold_values.data(),threshold_probabilities.data());
      
      Log("PhaseIIADCCalibrator Tool: Quantiles were: " + to_string(threshold_values.at(0))
            +" and "+to_string(threshold_values.at(1)),v_debug,verbosity);
      
      // since we know it: XXX note this is before second baseline subtraction!... not really accurate
      Log("PhaseIIADCCalibrator Tool: Getting baseline sigma", v_debug, verbosity);
      sigma_baseline = raw_datapoint_hist->GetStdDev();
      
      // draw the histogram for check
      if(draw_baseline_fit){
        Log("PhaseIIADCCalibrator Tool: Drawing histogrammed data for quantile determination",v_message,verbosity);
        
        if(gROOT->FindObject("baselineFitCanvas")==nullptr){
          Log("PhaseIIADCCalibrator Tool: Constructing canvas for drawing fit", v_debug, verbosity);
          //Log("PhaseIIADCCalibrator Tool: Making new Canvas",v_error,verbosity);
          baselineFitCanvas = new TCanvas("baselineFitCanvas");
        } else {
          baselineFitCanvas->cd();
        }
        
        raw_datapoint_hist->Draw();
        raw_datapoint_hist->SetDirectory(0); // canvas should not take ownership
        baselineFitCanvas->Modified();
        baselineFitCanvas->Update();
        gSystem->ProcessEvents();
        Log("PhaseIIADCCalibrator Tool: Sleeping while waiting for user to close canvas",v_debug,verbosity);
        while(gROOT->FindObject("baselineFitCanvas")!=nullptr){
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          gSystem->ProcessEvents();
        }
      }
      
      // now we can remove any outliers
      Log("PhaseIIADCCalibrator Tool: Erasing outliers", v_debug, verbosity);
      auto newend = std::remove_if(non_outlier_points.begin(), non_outlier_points.end(),
         [&threshold_values](double& dataval){
           return ( (dataval<threshold_values.front()) || (dataval>threshold_values.back()) );
         }
      );
      non_outlier_points.erase(newend, non_outlier_points.end());
      
      // update the contents of the TGraph with the outliers removed
      Log("PhaseIIADCCalibrator Tool: Updating TGraph", v_debug, verbosity);
      for(uint samplei=0; samplei<non_outlier_points.size(); ++samplei){
        calibrated_waveform_tgraph->SetPoint(samplei,samplei,non_outlier_points.at(samplei));
      }
      for(uint samplei=non_outlier_points.size(); samplei<raw_data.size(); ++samplei){
        calibrated_waveform_tgraph->RemovePoint(samplei);
      }
      calibrated_waveform_tgraph->Set(non_outlier_points.size()); // resize, doesn't seem to work, hence remove
      
      // note that this time we have fewer datapoints than before,
      // so we will need to restrict the range of our fit
      Log("PhaseIIADCCalibrator Tool: Setting TF1 range for outlier removed data", v_debug, verbosity);
      calibrated_waveform_fit->SetRange(0, non_outlier_points.size());
      calibrated_waveform_fit->SetMinimum(0);
      calibrated_waveform_fit->SetMaximum(non_outlier_points.size());
      
      // now redo the fit as we did before but with the remaining data
      Log("PhaseIIADCCalibrator Tool: Redoing fit", v_debug, verbosity);
      fit_result = calibrated_waveform_tgraph->Fit(calibrated_waveform_fit,"RCFSQ");
      fit_succeeded = ((static_cast<Int_t>(fit_result))==0);  // successful fit is 0
      
      // Draw the result of the re-fit
      if(draw_baseline_fit){
        Log("PhaseIIADCCalibrator Tool: Drawing outlier-subtracted baseline fit",v_message,verbosity);
        if(gROOT->FindObject("baselineFitCanvas")==nullptr){
          Log("PhaseIIADCCalibrator Tool: Constructing canvas for drawing fit", v_debug, verbosity);
          //Log("PhaseIIADCCalibrator Tool: Making new Canvas",v_error,verbosity);
          baselineFitCanvas = new TCanvas("baselineFitCanvas");
        } else {
          baselineFitCanvas->cd();
        }
        calibrated_waveform_tgraph->Draw("ALP");
        calibrated_waveform_tgraph->GetHistogram()->SetDirectory(0); // canvas should not take ownership
        calibrated_waveform_tgraph->ResetBit(kCanDelete);
        baselineFitCanvas->Modified();
        baselineFitCanvas->Update();
        gSystem->ProcessEvents();
        Log("PhaseIIADCCalibrator Tool: Sleeping while waiting for user to close canvas",v_debug,verbosity);
        while(gROOT->FindObject("baselineFitCanvas")!=nullptr){
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          gSystem->ProcessEvents();
        }
      }
      
      if(not fit_succeeded){
        Log("PhaseIIADCCalibrator Tool: polynomial re-fit of baseline failed!",v_warning,verbosity);
      } else {
        logmessage="PhaseIIADCCalibrator Tool: Baseline re-fit success: fit function was: ";
        for(uint orderi=0; orderi<(baseline_fit_order+1); ++orderi){
          logmessage+= to_string(fit_result->Value(orderi))+"*x^"+to_string(orderi);
          if(orderi<baseline_fit_order) logmessage+=" + ";
        }
        Log(logmessage, v_debug, verbosity);
        
        Log("PhaseIIADCCalibrator Tool: Combining with previous fit for final fit parameters", v_debug, verbosity);
        // get the new fit parameters and add them to the previous ones.
        for(uint orderi=0; orderi<(baseline_fit_order+1); ++orderi){
          double new_parameter_val = fit_result->Value(orderi) + fitpars.at(orderi);
          calibrated_waveform_fit->SetParameter(orderi, new_parameter_val);
        }
        logmessage="PhaseIIADCCalibrator Tool: Final fit function was: ";
        for(uint orderi=0; orderi<(baseline_fit_order+1); ++orderi){
          logmessage+= to_string(fit_result->Value(orderi))+"*x^"+to_string(orderi);
          if(orderi<baseline_fit_order) logmessage+=" + ";
        }
        Log(logmessage,v_debug,verbosity);
        
        baseline = calibrated_waveform_fit->GetParameter(0); // DC offset. FIXME doesn't fully capture correction
        
        // update our calibrated values
        Log("PhaseIIADCCalibrator Tool: Updating calibrated data based on new fit", v_debug, verbosity);
        for(uint samplei=0; samplei<raw_data.size(); ++samplei){
          const unsigned short& sample = raw_data.at(samplei);
          double baseline_val = calibrated_waveform_fit->Eval(samplei);
          double cal_val = (static_cast<double>(sample) - baseline_val)*ADC_TO_VOLT;
          cal_data.at(samplei)=cal_val;
        }
      }
      
      // revert our fit range for the next one before we forget
      Log("PhaseIIADCCalibrator Tool: Resetting TF1 range to default", v_debug, verbosity);
      calibrated_waveform_fit->SetRange(baseline_start_sample,num_baseline_samples);
      calibrated_waveform_fit->SetMinimum(baseline_start_sample);
      calibrated_waveform_fit->SetMaximum(num_baseline_samples);
      
      // Draw the final baseline subtracted data, and a fit, which by defn should be a straight line through 0
      if(draw_baseline_fit){
        // update the contents of the TGraph with the final datapoints
        Log("PhaseIIADCCalibrator Tool: Updating TGraph", v_debug, verbosity);
        calibrated_waveform_tgraph->Set(cal_data.size()); // resize
        for(uint samplei=0; samplei<cal_data.size(); ++samplei){
          calibrated_waveform_tgraph->SetPoint(samplei,samplei,cal_data.at(samplei));
        }
        
        // redo the fit; this time just for check
        Log("PhaseIIADCCalibrator Tool: Redoing fit once more just to check", v_debug, verbosity);
        fit_result = calibrated_waveform_tgraph->Fit(calibrated_waveform_fit,"RCFSQ");
        fit_succeeded = ((static_cast<Int_t>(fit_result))==0);  // successful fit is 0
        
        Log("PhaseIIADCCalibrator Tool: Drawing fit to final baseline-subtracted data",v_message,verbosity);
        if(gROOT->FindObject("baselineFitCanvas")==nullptr){
          Log("PhaseIIADCCalibrator Tool: Constructing canvas for drawing fit", v_debug, verbosity);
          //Log("PhaseIIADCCalibrator Tool: Making new Canvas",v_error,verbosity);
          baselineFitCanvas = new TCanvas("baselineFitCanvas");
        } else {
          baselineFitCanvas->cd();
        }
        calibrated_waveform_tgraph->Draw("ALP");
        calibrated_waveform_tgraph->GetHistogram()->SetDirectory(0); // canvas should not take ownership
        calibrated_waveform_tgraph->ResetBit(kCanDelete);
        baselineFitCanvas->Modified();
        baselineFitCanvas->Update();
        gSystem->ProcessEvents();
        Log("PhaseIIADCCalibrator Tool: Sleeping while waiting for user to close canvas",v_debug,verbosity);
        while(gROOT->FindObject("baselineFitCanvas")!=nullptr){
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          gSystem->ProcessEvents();
        }
      }
    }
    
    // construct the calibrated waveform
    Log("PhaseIIADCCalibrator Tool: Constructing calibrated waveform", v_debug, verbosity);
    calibrated_waveforms.emplace_back(raw_waveform.GetStartTime(), cal_data, baseline, sigma_baseline);
  }
  
  return calibrated_waveforms;
}

void PhaseIIADCCalibrator::make_raw_led_waveforms(unsigned long channel_key,
  const std::vector< Waveform<unsigned short> > raw_waveforms,
  std::vector< Waveform<unsigned short> >& raw_led_waveforms)
{
  //Get the windows for this channel key
  std::vector<std::vector<int>> thispmt_adc_windows;
  thispmt_adc_windows = this->get_db_windows(channel_key);
  for(int j=0; j<raw_waveforms.size(); j++){
    for(int i = 0; i<thispmt_adc_windows.size();i++){
      // Make a waveform out of this subwindow, plus the number of samples
      // used for background estimation
      std::vector<int> awindow = thispmt_adc_windows.at(i);
      int windowmin = awindow.at(0) - ((int)num_baseline_samples);
      if(windowmin<0){
        std::cout << "PhaseIIADCCalibrator Tool WARNING: when making an " <<
            "LED window, there was not enough room prior to the window to " <<
            "include samples for a background estimate.  Don't put your window "<<
           " too close to zero ADC counts." << std::endl;
        windowmin = 0;
      }
      int windowmax = awindow.at(1);

      std::vector<uint16_t> led_waveform;
      //FIXME: want start time of window, not of the full raw waveform
      uint64_t start_time = raw_waveforms.at(j).GetStartTime();
      for (size_t s = windowmin; s < windowmax; ++s) {
        led_waveform.push_back(raw_waveforms.at(j).GetSample(s));
      }
      Waveform<uint16_t> led_rawwave{start_time,led_waveform};
      raw_led_waveforms.push_back(led_rawwave);
    }
  }
  return;
}

std::vector<std::vector<int>> PhaseIIADCCalibrator::get_db_windows(unsigned long channelkey){
    std::vector<std::vector<int>> this_pmt_windows;
  //Look in the map and check if channelkey exists.
  if (channel_window_map.find(channelkey) == channel_window_map.end() ) {
     if (verbosity>3){
       std::cout << "PhaseIIADCHitFinder Warning: no integration windows found" <<
       "for channel_key" << channelkey <<". Not finding pulses." << std::endl;
       }
  } else {
    // gottem
    this_pmt_windows = channel_window_map.at(channelkey);
  }
  return this_pmt_windows;
}

std::map<unsigned long, std::vector<std::vector<int>>> PhaseIIADCCalibrator::load_window_map(std::string window_db)
{
  std::map<unsigned long, std::vector<std::vector<int>>> chanwindowmap;
  std::string fileline;
  ifstream myfile(window_db.c_str());
  if (myfile.is_open()){
    while(getline(myfile,fileline)){
      if(fileline.find("#")!=std::string::npos) continue;
      std::cout << fileline << std::endl; //has our stuff;
      std::vector<std::string> dataline;
      boost::split(dataline,fileline, boost::is_any_of(","), boost::token_compress_on);
      unsigned long chanvalue = std::stoul(dataline.at(0));
      int windowminvalue = std::stoi(dataline.at(1));
      int windowmaxvalue = std::stoi(dataline.at(2));
      std::vector<int> window{windowminvalue,windowmaxvalue};
      if(chanwindowmap.count(chanvalue)==0){ // Place window range into integral windows
        std::vector<std::vector<int>> window_vector{window};
        chanwindowmap.emplace(chanvalue,window_vector);
      }
      else { // Add this window to the windows to integrate in
        chanwindowmap.at(chanvalue).push_back(window);
      }
    }
  } else {
    Log("PhaseIIADCHitFinder Tool: Input integration window DB file not found. "
        " no integration will occur. ",
        1, verbosity);
  }
  return chanwindowmap;
}

void PhaseIIADCCalibrator::CardIDToElectronicsSpace(int CardID,
        int &CrateNum, int &SlotNum)
{
  //CardID = CrateNum * 1000 + SlotNum.  This logic works if we have less than
  // 10 crates and less than 100 Slots (which we do).
  SlotNum = CardID % 100;
  CrateNum = CardID / 1000;
  return;
}
