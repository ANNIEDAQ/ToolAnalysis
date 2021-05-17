#ifndef MonitorTrigger_H
#define MonitorTrigger_H

#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "Tool.h"

#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TText.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TObjectTable.h"


/**
 * \class MonitorTrigger
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M.Nieslony $
* $Date: 2020/07/25 19:21:00 $
* Contact: mnieslon@uni-mainz.de
*/
class MonitorTrigger: public Tool {


 public:

  MonitorTrigger(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  //Configuration and initialization functions
  void ReadInConfiguration();
  void InitializeHists();
  void LoopThroughDecodedEvents(std::map<uint64_t,std::vector<uint32_t>> timetotriggerword);
  void WriteToFile();
  void ReadFromFile(ULong64_t timestamp_end, double time_frame);
  void DrawLastFilePlots();
  void UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes);

  void DrawFileHistoryTrig(ULong64_t timestamp_end, double time_frame, std::string file_ending, int _linewidth);
  void PrintFileTimeStamp(ULong64_t timestamp_end, double time_frame, std::string file_ending);

  void DrawFrequencyRatePlots(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawTimeAlignmentPlots();
  std::vector<int> LoadTriggerMask(std::string triggermask_file);
  std::map<int,std::string> LoadTriggerWord(std::string file_triggerword);
  std::vector<std::vector<int>> LoadTriggerAlign(std::string file_triggeralign);

  //helper functions
  std::string convertTimeStamp_to_Date(ULong64_t timestamp);
  bool does_file_exist(std::string filename);

 private:
 
  //Configuration variables
  std::string outpath_temp;
  std::string outpath;
  std::string StartTime;
  double update_frequency;
  std::string plot_configuration;
  std::string path_monitoring;
  std::string img_extension;
  std::string triggermaskfile="none";
  std::string triggerwordfile="none";
  std::string triggeralignfile="none";
  bool force_update;
  bool draw_marker;
  int verbosity;

  //Configuration option for plots
  std::vector<double> config_timeframes;
  std::vector<std::string> config_endtime;
  std::vector<std::string> config_label;
  std::vector<std::vector<std::string>> config_plottypes;
  std::vector<ULong64_t> config_endtime_long;

  //Decoded data variable
  std::map<uint64_t,std::vector<uint32_t>> TimeToTriggerWordMap;
  std::vector<int> TriggerMask;
  std::map<int,std::string> TriggerWord;
  std::vector<std::vector<int>> TriggerAlign;

  //Variables for time calculations
  boost::posix_time::ptime *Epoch;
  boost::posix_time::ptime current;
  boost::posix_time::time_duration period_update;
  boost::posix_time::time_duration duration;
  boost::posix_time::ptime last;
  boost::posix_time::ptime utc;
  boost::posix_time::time_duration current_stamp_duration;
  boost::posix_time::time_duration current_utc_duration;
  time_t t;
  std::stringstream title_time;
  long current_stamp, current_utc;
  ULong64_t readfromfile_tend;
  double readfromfile_timeframe;

  //Variables to convert times
  double MSEC_to_SEC = 1000.;
  double SEC_to_MIN = 60.;
  double MIN_to_HOUR = 60.;
  double HOUR_to_DAY = 24.;
  ULong64_t utc_to_fermi = 2.7e12;  //6h clock delay in ADC clocks (timestamps in UTC time compared to Fermilab time)
  ULong64_t utc_to_t=21600000;  //6h clock delay in millisecons

  //General variables
  int num_triggerwords = 64;
  int num_triggerwords_selected = 5;

  //Live storing variables (for data in current file)
  std::vector<int> frequency_file;
  std::vector<uint32_t> triggerword_file;
  std::vector<uint64_t> timestamp_file;
  long t_file_start, t_file_end;

  //Storing variables for reading in data in given time slot from monitoring file / database
  std::vector<std::vector<int>> triggerword_plot;
  std::vector<std::vector<int>> frequency_plot;
  std::vector<std::vector<double>> rate_plot;
  std::vector<ULong64_t> tstart_plot;
  std::vector<ULong64_t> tend_plot;

  //Histograms
  TH1F *h_triggerword = nullptr;
  TH1F *h_triggerword_selected = nullptr;
  TH1F *h_rate_triggerword = nullptr;
  TH1F *h_rate_triggerword_selected = nullptr;
  std::vector<TH1F*> h_timestamp;
  std::vector<TH1F*> h_triggeralign;
  std::vector<TDatime> labels_timeaxis;
  int num_time_bins=200;
  int num_align_bins=200;

  //TGraphs & related objects
  std::vector<TGraph*> gr_rate;
  TMultiGraph *multi_trig_rate = nullptr;
  TLegend *leg_rate = nullptr;

  //Canvas
  TCanvas *canvas_triggerword = nullptr;
  TCanvas *canvas_triggerword_selected = nullptr;
  TCanvas *canvas_rate_triggerword = nullptr;
  TCanvas *canvas_rate_triggerword_selected = nullptr;
  TCanvas *canvas_timestamp = nullptr;
  TCanvas *canvas_triggeralign = nullptr;
  TCanvas *canvas_timeevolution = nullptr;

  //Histograms showing history (log) of files
  TH1F *log_files_trig = nullptr;
  TCanvas *canvas_logfile_trig = nullptr;
  TCanvas *canvas_file_timestamp_trig = nullptr;
  int num_files_history;

  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;



};


#endif
