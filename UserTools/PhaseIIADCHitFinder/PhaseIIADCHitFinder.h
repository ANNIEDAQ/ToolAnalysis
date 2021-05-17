/* vim:set noexpandtab tabstop=4 wrap */

// This tool finds PMT pulses using both raw and calibrated ADC waveforms
// Modified for Phase II by Teal Pershing <tjpershing@ucdavis.edu>
#ifndef PhaseIIADCHitFinder_H
#define PhaseIIADCHitFinder_H

#pragma once

// standard library includes
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>

// ToolAnalysis includes
#include "ADCPulse.h"
#include "CalibratedADCWaveform.h"
#include "Hit.h"

#include "ANNIEconstants.h"
#include "Geometry.h"
#include "Tool.h"
#include "Waveform.h"
#include "Constants.h"
#include "Channel.h"
#include <boost/algorithm/string.hpp>

class PhaseIIADCHitFinder : public Tool {

  public:

    PhaseIIADCHitFinder();
    bool Initialise(const std::string configfile, DataModel& data) override;
    bool Execute() override;
    bool Finalise() override;
    
    // verbosity levels: if 'verbosity' < this level, the message type will be logged.
    int verbosity;
    int v_error=0;
    int v_warning=1;
    int v_message=2;
    int v_debug=3;
  protected:

    Geometry *geom = nullptr;

    //Configurables for HitFinding tool; see README for details
    std::string pulse_finding_approach;
    unsigned short default_adc_threshold;
    std::string threshold_type;
    std::string adc_threshold_db;
    std::string adc_window_db;
    std::string pulse_window_type;
    bool use_led_waveforms;
    int pulse_window_start_shift;
    int pulse_window_end_shift;
    std::map<unsigned long, unsigned short> channel_threshold_map;
    std::map<unsigned long, std::vector<std::vector<int>>> channel_window_map;
    bool eventbuilding_mode; 
   
   
    std::map<int,std::string>* AuxChannelNumToTypeMap;

     // Load the map containing the ADC calibrated waveform data
    std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > >
      calibrated_waveform_map;
    std::map<unsigned long, std::vector<CalibratedADCWaveform<double> > >
      calibrated_aux_waveform_map;
    std::map<unsigned long, std::vector<Waveform<unsigned short> > >
      raw_waveform_map;
    std::map<unsigned long, std::vector<Waveform<unsigned short> > >
      raw_aux_waveform_map;
    
    // Build the map of pulses and Hit Map
    std::map<unsigned long, std::vector< std::vector<ADCPulse>> > pulse_map;
    std::map<unsigned long,std::vector<Hit>>* hit_map;

    std::map<unsigned long, std::vector< std::vector<ADCPulse>> > aux_pulse_map;
    std::map<unsigned long,std::vector<Hit>>* aux_hit_map;

    std::vector<unsigned long> chkey_map;

    // Load a PMT's threshold from the channel_threshold_map. If none, returns default ADC threshold
    unsigned short get_db_threshold(unsigned long channelkey);

    // Load a PMT's integration windows from the channel_window_map. If none, returns an empty vector.
    std::vector<std::vector<int>> get_db_windows(unsigned long channelkey);

    // load a channel threshold map from the source file given
    std::map<unsigned long, unsigned short> load_channel_threshold_map(std::string threshold_db);

    // load an integration window map (CSV file) from the source file given
    std::map<unsigned long, std::vector<std::vector<int>>> load_integration_window_map(std::string window_db);

    void ClearMaps();
    bool build_pulse_and_hit_map(unsigned long ckey,
       std::vector<Waveform<unsigned short> > rawmap, 
      std::vector<CalibratedADCWaveform<double> > calmap,
      std::map<unsigned long, std::vector< std::vector<ADCPulse>> > & pmap,
      std::map<unsigned long,std::vector<Hit>>& hmap);
    // Create a vector of ADCPulse objects using the raw and calibrated signals
    // from a given minibuffer. Note that the vectors of raw and calibrated
    // samples are assumed to be the same size. This function will throw an
    // exception if this assumption is violated.
    std::vector<ADCPulse> find_pulses_bythreshold(
      const Waveform<unsigned short>& raw_minibuffer_data,
      const CalibratedADCWaveform<double>& calibrated_minibuffer_data,
      unsigned short adc_threshold, const unsigned long& channel_key) const;

    std::vector<ADCPulse> find_pulses_bywindow(
      const Waveform<unsigned short>& raw_minibuffer_data,
      const CalibratedADCWaveform<double>& calibrated_minibuffer_data,
      std::vector<std::vector<int>> adc_windows, const unsigned long& channel_key,
      bool MaxHeightPulseOnly) const;

    //Takes the ADC pulse vectors (one per minibuffer) and converts them to a vector of hits
    std::vector<Hit> convert_adcpulses_to_hits(unsigned long channel_key,std::vector<std::vector<ADCPulse>> pulses);

    //EventBuilding mode
    std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>> *FinishedRawWaveforms;      //Key: {MTCTime}, value: map of raw waveforms
    std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>> *FinishedRawWaveformsAux;  //Key: {MTCTime}, value: map of raw waveforms (aux channels)
    std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>> *FinishedCalibratedWaveforms;  //Key: {MTCTime}, value: map of calibrated waveforms
    std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>> *FinishedCalibratedWaveformsAux;  //Key: {MTCTime}, value: map of calibrated waveforms (aux channels)
    std::map<uint64_t, std::map<unsigned long,std::vector<Hit>>*> *InProgressHits;	//Key: {MTCTime}, value: map of  Hit distributions
    std::map<uint64_t, std::map<unsigned long,std::vector<Hit>>*> *InProgressHitsAux;	//Key: {MTCTime}, value: map of  Hit distributions
    std::map<uint64_t, std::map<unsigned long,std::vector<std::vector<ADCPulse>>>> *InProgressRecoADCHits; //Key: {MTCTime}, value: map of found pulses
    std::map<uint64_t, std::map<unsigned long,std::vector<std::vector<ADCPulse>>>> *InProgressRecoADCHitsAux; //Key: {MTCTime}, value: map of found pulses
    std::map<uint64_t, std::vector<unsigned long>> *InProgressChkey;	//Key: {MTCTime}, value: vector of in progress chankeys
    //std::map<uint64_t, std::map<unsigned long,std::vector<Hit>>*> *FinishedHits;	//Key: {MTCTime}, value: map of  Hit distributions
    //std::map<uint64_t, std::map<unsigned long,std::vector<Hit>>*> *FinishedHitsAux;	//Key: {MTCTime}, value: map of  Hit distributions
    //std::map<uint64_t, std::map<unsigned long,std::vector<std::vector<ADCPulse>>>> *FinishedRecoADCHits; //Key: {MTCTime}, value: map of found pulses
    //std::map<uint64_t, std::map<unsigned long,std::vector<std::vector<ADCPulse>>>> *FinishedRecoADCHitsAux; //Key: {MTCTime}, value: map of found pulses

};

#endif
