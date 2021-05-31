#ifndef MonitorMRDTime_H
#define MonitorMRDTime_H

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "Tool.h"
#include "MRDOut.h"
#include "Geometry.h"
#include "Detector.h"
#include "Paddle.h"

#include "TObjectTable.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TThread.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TDatime.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TPad.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TH2Poly.h"
#include "TPie.h"
#include "TPieSlice.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

class MonitorMRDTime: public Tool {

 public:

  MonitorMRDTime();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  //main functions for organizing the monitoring
  void ReadInConfiguration();
  void InitializeVectors();
  void ReadInData();
  void WriteToFile();
  void ReadFromFile(ULong64_t timestamp_end, double time_frame);
  
  void DrawLastFilePlots();
  void UpdateMonitorPlots(std::vector<double> timeFrames, std::vector<ULong64_t> endTimes, std::vector<std::string> fileLabels, std::vector<std::vector<std::string>> plotTypes);
  void DrawScatterPlots();
  void DrawScatterPlotsTrigger();
  void DrawTDCHistogram();
  void DrawHitMap(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawRatePlotElectronics(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawRatePlotPhysical(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawTimeEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawTriggerEvolution(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawFileHistory(ULong64_t timestamp_end, double time_frame, std::string file_ending, int _linewidth);
  void PrintFileTimeStamp(ULong64_t timestamp_end, double time_frame, std::string file_ending);
  void DrawPieChart(ULong64_t timestamp_end, double time_frame, std::string file_ending);

  //helper functions
  std::string convertTimeStamp_to_Date(ULong64_t timestamp);
  bool does_file_exist(std::string filename);

 private:

  //define MRD stores that contain the data
  BoostStore* MRDdata;
  MRDOut MRDout;
  bool bool_mrddata;

  //MRD store includes the following variables
  unsigned int OutN, Trigger;
  std::vector<unsigned int> Value, Slot, Channel, Crate;
  std::vector<std::string> Type;
  ULong64_t TimeStamp;
  int nentries_MRDData;

  //define configuration file variables
  BoostStore *CCData;
  std::string outpath_temp;
  std::string outpath;
  std::string active_slots;
  std::string inactive_channels;
  std::string loopback_channels;
  std::string StartTime;
  std::string path_monitoring;
  std::string img_extension;
  double update_frequency;
  bool force_update;
  bool draw_marker;
  bool draw_single;
  std::string plot_configuration;
  int verbosity;

  //define variables that contain the configuration option for the plots
  std::vector<double> config_timeframes;
  std::vector<std::string> config_endtime;
  std::vector<std::string> config_label;
  std::vector<std::vector<std::string>> config_plottypes;
  std::vector<ULong64_t> config_endtime_long;

  //define variables for keeping track of time
  boost::posix_time::ptime *Epoch;
  boost::posix_time::ptime current;
  boost::posix_time::ptime utc;
  boost::posix_time::ptime last;
  boost::posix_time::time_duration period_update;
  boost::posix_time::time_duration duration;
  boost::posix_time::time_duration current_stamp_duration;
  boost::posix_time::time_duration current_utc_duration;
  long current_stamp, current_utc;
  time_t t;
  std::stringstream title_time; 
  ULong64_t readfromfile_tend;
  double readfromfile_timeframe;

  //variables to convert times
  double MSEC_to_SEC = 1000.;
  double SEC_to_MIN = 60.;
  double MIN_to_HOUR = 60.;
  double HOUR_to_DAY = 24.;

  //define variables for active slot/channel management
  //MRD electronics is organized in 2 CAMAC crates, which 24 slots each.
  //TDCs in the slots have 32 channels each
  //The crate numbers in the data start from 7, indicating the crate position on rack #7
  static const int num_crates = 2;
  static const int num_slots = 24;
  static const int num_channels = 32;
  static const int min_crate = 7;

  //Only some of the slots have TDCs in them, they are specified in a configuration file
  //Of those TDC slots, not all channels are active. The inactive ones are also read in from a configuration file
  //There are 2 loopback channels duplicating the trigger signals and indicating whether there has been a cosmic or beam trigger. 
  //The positions of those loopback channels are also read in from a configuration file
  int num_active_slots;
  int num_active_slots_cr1, num_active_slots_cr2;
  int active_channel[2][24]={{0}};
  std::vector<unsigned int> active_slots_cr1;
  std::vector<unsigned int> active_slots_cr2;
  std::vector<int> nr_slot;
  std::map<std::vector<unsigned int>,int> CrateSlot_to_ActiveSlot;
  std::map<std::vector<unsigned int>,int> CrateSlotChannel_to_TotalChannel;
  std::map<int,unsigned int> ActiveSlot_to_Crate;
  std::map<int,unsigned int> ActiveSlot_to_Slot;
  std::map<int,unsigned int> TotalChannel_to_Crate;
  std::map<int,unsigned int> TotalChannel_to_Slot;
  std::map<int,unsigned int> TotalChannel_to_Channel;
  std::vector<unsigned int> inactive_ch_crate1, inactive_slot_crate1;
  std::vector<unsigned int> inactive_ch_crate2, inactive_slot_crate2;
  int inactive_crate1, inactive_crate2;
  std::vector<std::string> loopback_name;
  std::vector<unsigned int> loopback_crate, loopback_slot, loopback_channel;
  std::vector<int> mapping_vector_ch;
  int beam_ch=0;
  int cosmic_ch=0;

  //define live storing variables (for data of current file)
  std::vector<std::vector<int>> tdc_file;
  std::vector<std::vector<ULong64_t>> timestamp_file;
  std::vector<std::vector<int>> tdc_file_times;
  std::vector<std::vector<int>> tdc_file_ch;
  ULong64_t t_file_start, t_file_end;
  ULong64_t utc_to_t=21600000;	//6h clock delay (MRD timestamped in UTC time)
  long n_doublehits;
  long n_zerohits;
  long n_noloopback;
  long n_normalhits;

  //define storing variables for reading in data in a given time slot from monitoring file / database
  std::vector<std::vector<double>> tdc_plot;
  std::vector<std::vector<double>> rms_plot;
  std::vector<std::vector<double>> rate_plot;
  std::vector<std::vector<int>> channelcount_plot;
  std::vector<ULong64_t> tstart_plot;
  std::vector<ULong64_t> tend_plot;
  std::vector<double> cosmicrate_plot;
  std::vector<double> beamrate_plot;
  std::vector<double> noloopbackrate_plot;
  std::vector<double> normalhitrate_plot;
  std::vector<double> doublehitrate_plot;
  std::vector<double> zerohitsrate_plot;
  std::vector<int> nevents_plot;

  //labels, lines, etc.
  TPaveText *label_cr1 = nullptr;
  TPaveText *label_cr2 = nullptr;
  std::vector<TDatime> labels_timeaxis;
  TLine *separate_crates = nullptr;
  TF1 *f1 = nullptr;

  //define graphs for time evolution plots
  std::vector<TGraph*> gr_tdc, gr_rms, gr_rate;
  std::vector<TGraph*> gr_trigger;
  TGraph *gr_noloopback = nullptr;
  TGraph *gr_zerohits = nullptr;
  TGraph *gr_doublehits = nullptr;
  TMultiGraph *multi_ch_tdc = nullptr;
  TMultiGraph *multi_ch_rms = nullptr;
  TMultiGraph *multi_ch_rate = nullptr;
  TMultiGraph *multi_trigger = nullptr;
  TMultiGraph *multi_eventtypes = nullptr;
  TLegend *leg_tdc = nullptr;
  TLegend *leg_rms = nullptr;
  TLegend *leg_rate = nullptr;
  TLegend *leg_trigger = nullptr;
  TLegend *leg_noloopback = nullptr;
  TLegend *leg_eventtypes = nullptr;

  //define scatter plot histograms
  std::vector<TH2F*> hist_scatter;
  TLegend *leg_scatter = nullptr;
  TLegend *leg_scatter_trigger = nullptr;
  int n_bins_scatter;

  //define hitmap histograms
  TH1F *hist_hitmap_cr1 = nullptr;
  TH1F *hist_hitmap_cr2 = nullptr;
  std::vector<TH1F*> hist_hitmap_slot;
  std::map<unsigned int,std::vector<TBox*>> vector_box_inactive_hitmap;

  //define TDC histogram
  TH1F *hist_tdc = nullptr;
  TH1F *hist_tdc_cluster = nullptr;
  TH1F *hist_tdc_cluster_20 = nullptr;
  TH1F *hist_tdc_beam = nullptr;
  TH1F *hist_tdc_cluster_beam = nullptr;
  TH1F *hist_tdc_cluster_beam_20 = nullptr;
  TH1F *hist_tdc_cosmic = nullptr;
  TH1F *hist_tdc_cluster_cosmic = nullptr;
  TH1F *hist_tdc_cluster_cosmic_20 = nullptr;
  std::vector<std::vector<int>> overall_mrd_coinc_times;
  std::vector<std::vector<int>> overall_mrd_coinc_times_beam;
  std::vector<std::vector<int>> overall_mrd_coinc_times_cosmic;

  //define rate histograms
  TH2F *rate_crate1=nullptr;
  TH2F *rate_crate2=nullptr;
  TLatex *label_rate_cr1 = nullptr;
  TLatex *label_rate_cr2 = nullptr; 
  TLatex *label_rate_facc = nullptr;
  std::vector<TBox*> vector_box_inactive;

  //geometry conversion table
  Geometry *geom = nullptr;
  double tank_center_x, tank_center_y, tank_center_z;
  std::map<std::vector<int>,int>* CrateSpaceToChannelNumMap = nullptr;
  TH2Poly *rate_top = nullptr;
  TH2Poly *rate_side = nullptr;
  TH2Poly *rate_facc = nullptr;
  double enlargeBoxes = 0.01;
  double shiftSecRow = 0.04;

  //define histogram showing the history (log) of files 
  TH1F *log_files_mrd=nullptr;
  int num_files_history;

  //define TPie objects showing pie charts of the event type composition
  TPie *pie_triggertype = nullptr;
  TPie *pie_weirdevents = nullptr;
  TLegend *leg_triggertype = nullptr;
  TLegend *leg_weirdevents = nullptr;

  //define pointers to canvases
  TCanvas *canvas_hitmap = nullptr;
  TCanvas *canvas_hitmap_slot = nullptr;
  TCanvas *canvas_tdc = nullptr;
  TCanvas *canvas_rate_electronics = nullptr;
  TCanvas *canvas_rate_physical = nullptr;
  TCanvas *canvas_rate_physical_facc = nullptr;
  TCanvas *canvas_logfile_mrd = nullptr;
  TCanvas *canvas_ch_tdc = nullptr;
  TCanvas *canvas_ch_rms = nullptr;
  TCanvas *canvas_ch_rate = nullptr;
  TCanvas *canvas_ch_single = nullptr;
  TCanvas *canvas_trigger = nullptr;
  TCanvas *canvas_trigger_time = nullptr;
  TCanvas *canvas_scatter = nullptr;
  TCanvas *canvas_scatter_single = nullptr;
  TCanvas *canvas_pie = nullptr;
  TCanvas *canvas_file_timestamp = nullptr;

  //define colors for multitude of TGraphs for each channel
  int color_scheme[16] = {1,2,3,4,5,6,7,8,9,13,15,205,94,220,221,225}; 

};


#endif
