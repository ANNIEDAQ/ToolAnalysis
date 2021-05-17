// This tool creates a CalibratedWaveform object for each RawWaveform object
// that it finds stored under the "RawADCData" key in the ANNIEEvent store.
// It saves the CalibratedWaveform objects to the ANNIEEvent store using the
// key "CalibratedADCData".
//
// Phase I version by Steven Gardiner <sjgardiner@ucdavis.edu>
// Modified for Phase II by Teal Pershing <tjpershing@ucdavis.edu>
#pragma once

// ToolAnalysis includes
#include "CalibratedADCWaveform.h"
#include "Tool.h"
#include "Waveform.h"
#include "annie_math.h"
#include "ANNIEalgorithms.h"
#include "ANNIEconstants.h"
#include <boost/algorithm/string.hpp>

#include <sstream>

class TApplication;
class TCanvas;
class TGraph;
class TF1;
class TH1D;

class PhaseIIADCCalibrator : public Tool {

  public:

    PhaseIIADCCalibrator();
    bool Initialise(const std::string configfile,DataModel& data) override;
    bool Execute() override;
    bool Finalise() override;
    void CardIDToElectronicsSpace(int CardID, int &CrateNum, int &SlotNum);

  protected:

    /// @brief Compute the baseline for a particular RawChannel
    /// object using a technique taken from the ZE3RA code.
    /// @details See section 2.2 of https://arxiv.org/pdf/1106.0808.pdf for a
    /// description of the algorithm.
    void ze3ra_baseline(const Waveform<unsigned short> raw_data,
      double& baseline, double& sigma_baseline, size_t num_baseline_samples, size_t starting_sample);

    std::vector< CalibratedADCWaveform<double> > make_calibrated_waveforms_ze3ra(
      const std::vector< Waveform<unsigned short> >& raw_waveforms);
    
    std::vector< CalibratedADCWaveform<double> > make_calibrated_waveforms_ze3ra_multi(
      const std::vector< Waveform<unsigned short> >& raw_waveforms);
    
    /// @brief Fit a polynomial to the baseline of each waveform.
    std::vector< CalibratedADCWaveform<double> > make_calibrated_waveforms_rootfit(
      const std::vector<Waveform<short unsigned int> >& raw_waveforms);
    
    /// @brief Calculate mean and standard deviation using num_baseline_samples at beginning of waveform
    std::vector< CalibratedADCWaveform<double> > make_calibrated_waveforms_simple(
      const std::vector<Waveform<short unsigned int> >& raw_waveforms);
    
    bool use_ze3ra_algorithm;
    bool use_root_algorithm;
 
    void make_raw_led_waveforms(unsigned long channel_key,
      const std::vector< Waveform<unsigned short> > raw_waveforms,
      std::vector< Waveform<unsigned short>>& LEDWaveforms);
    // Load a PMT's integration windows from the channel_window_map. If none, returns an empty vector.
    std::vector<std::vector<int>> get_db_windows(unsigned long channelkey);
   
    // load the LED pulse window map (CSV file) from the source file given
    std::map<unsigned long, std::vector<std::vector<int>>> load_window_map(std::string window_db);

    std::map<unsigned long, std::vector<std::vector<int>>> channel_window_map;

    std::string BEType;


    std::map<int,std::string>* AuxChannelNumToTypeMap;
    int verbosity;
    // All F-distribution probabilities above this value will pass the
    // variance consistency test in ze3ra_baseline(). That is, p_critical
    // is the maximum p-value for which we will reject the null hypothesis
    // of equal variances.
    double p_critical;
   
    //ze3ra and ze3ra_multi configurables 
    size_t num_baseline_samples;
    size_t num_sub_waveforms;

    //ze3ra_multi configurables
    size_t baseline_rep_samples;
    size_t baseline_unc_tolerance;

    bool make_led_waveforms;
    std::string adc_window_db; 
    
    size_t num_waveform_points;
    size_t baseline_start_sample;
    
    int baseline_fit_order;
    bool redo_fit_without_outliers;
    double refit_threshold; // V range of the initial baseline subtracted waveform must be > this to trigger refit
 
    //Variables specifically intended for Event Building
    bool eventbuilding_mode = false;
    //std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > > *FinishedTankEvents;  //Key: {MTCTime}, value: map of fully-built waveforms from WaveBank
    std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > > *InProgressTankEvents;  //Key: {MTCTime}, value: map of fully-built waveforms from WaveBank
    std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>> *FinishedRawWaveforms;	//Key: {MTCTime}, value: map of raw waveforms
    std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>> *FinishedRawWaveformsAux;  //Key: {MTCTime}, value: map of raw waveforms (aux channels)
    std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>> *FinishedCalibratedWaveforms;	//Key: {MTCTime}, value: map of calibrated waveforms
    std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>> *FinishedCalibratedWaveformsAux;  //Key: {MTCTime}, value: map of calibrated waveforms (aux channels)
    std::map<uint64_t, std::map<unsigned long,std::vector<CalibratedADCWaveform<double>>>> *FinishedCalibratedLEDADCData;  //Key: {MTCTime}, value: map of calibrated waveforms (LED channels)
    std::map<uint64_t, std::map<unsigned long,std::vector<Waveform<unsigned short>>>> *FinishedRawLEDADCData;	//Key: {MTCTime}, value: map of raw LED waveforms
    std::map<uint64_t, std::map<unsigned long,std::vector<int>>> *FinishedRawAcqSize;   // Key: {MTCTime}, value: map 

    std::map<std::vector<int>,int> TankPMTCrateSpaceToChannelNumMap;
    std::map<std::vector<int>,int> AuxCrateSpaceToChannelNumMap;

    int ExecuteCount;
    int ExecutesPerBuild;

    // ROOT stuff for drawing the fit of the baseline
    bool draw_baseline_fit=false;
    TApplication* rootTApp=nullptr;
    TCanvas* baselineFitCanvas=nullptr;
    TGraph* calibrated_waveform_tgraph=nullptr;
    TF1* calibrated_waveform_fit=nullptr;
    TH1D* raw_datapoint_hist=nullptr;
    int drawcount=0;
    
    // verbosity levels: if 'verbosity' < this level, the message type will be logged.
    int v_error=0;
    int v_warning=1;
    int v_message=2;
    int v_debug=3;
    std::string logmessage;
    int get_ok;
};
