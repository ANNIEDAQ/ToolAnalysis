#ifndef FMVEfficiency_H
#define FMVEfficiency_H

#include <string>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"

#include "Tool.h"


/**
 * \class FMVEfficiency
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M.Nieslony $
* $Date: 2020/03/31 11:47:00 $
* Contact: mnieslon@uni-mainz.de
*/
class FMVEfficiency: public Tool {


 public:

  FMVEfficiency(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool FindPaddleIntersection(Position startpos, Position endpos, double &x, double &y, double z);
  bool FindPaddleChankey(double x, double y, int layer, unsigned long &chankey);

 private:

  //configuration variables
  int verbosity;
  std::string singlePEgains;
  std::string outputfile;
  bool useTank;
  bool isData = true;
  int padPerLayer;

  //data objects
  std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
  std::map<unsigned long,std::vector<MCHit>>* TDCData_MC=nullptr;
  std::map<double,std::vector<Hit>>* m_all_clusters;  //from ClusterFinder tool
  std::map<double,std::vector<MCHit>>* m_all_clusters_MC;  //from ClusterFinder tool
  std::map<double,std::vector<unsigned long>>* m_all_clusters_detkey;  //from ClusterFinder tool
  std::vector<std::vector<int>> MrdTimeClusters;  //from TimeClustering tool
  std::vector<double> MrdDigitTimes;  //from TimeClustering tool
  std::vector<unsigned long> mrddigitchankeysthisevent;  //from TimeClustering tool
  std::map<unsigned long,double> pmt_gains; //from txt file
  std::vector<BoostStore>* theMrdTracks;                        // the reconstructed tracks

  //geometry
  Geometry *geom = nullptr;
  int n_veto_pmts;
  double fmv_firstlayer_z, fmv_secondlayer_z, fmv_xmin, fmv_xmax, fmv_x;
  unsigned long first_fmv_chankey=0;
  unsigned long first_fmv_detkey=0;

  //storing containers
  std::vector<double> fmv_firstlayer_ymin, fmv_firstlayer_ymax, fmv_firstlayer_y, fmv_secondlayer_ymin, fmv_secondlayer_ymax, fmv_secondlayer_y;
  std::vector<unsigned long> fmv_firstlayer, fmv_secondlayer;
  std::vector<int> fmv_firstlayer_expected, fmv_firstlayer_observed, fmv_firstlayer_expected_track_strict, fmv_firstlayer_observed_track_strict, fmv_firstlayer_expected_track_loose, fmv_firstlayer_observed_track_loose;
  std::vector<int> fmv_secondlayer_expected, fmv_secondlayer_observed, fmv_secondlayer_expected_track_strict, fmv_secondlayer_observed_track_strict, fmv_secondlayer_expected_track_loose, fmv_secondlayer_observed_track_loose;
  std::vector<int> fmv_tank_firstlayer_expected, fmv_tank_firstlayer_observed, fmv_tank_secondlayer_expected, fmv_tank_secondlayer_observed;

  //Root files, histograms
  TFile *file = nullptr;

  TH1F *time_diff_Layer1 = nullptr;
  TH1F *time_diff_Layer2 = nullptr;
  TH1F *num_paddles_Layer1 = nullptr;
  TH1F *num_paddles_Layer2 = nullptr;
  TH2F *time_diff_paddle = nullptr;
  TH2F *time_diff_paddle_coinc = nullptr;
  TH2F *time_diff_paddle_simplecoinc = nullptr;

  TH1F *fmv_layer1 = nullptr;
  TH1F *fmv_layer2 = nullptr;
  TH2F *fmv_layer1_layer2 = nullptr;

  std::vector<TH1F*> vector_observed_strict_layer1, vector_expected_strict_layer1, vector_observed_loose_layer1, vector_expected_loose_layer1;
  std::vector<TH1F*> vector_observed_strict_layer2, vector_expected_strict_layer2, vector_observed_loose_layer2, vector_expected_loose_layer2;
  TH1F *fmv_observed_layer1 = nullptr;
  TH1F *fmv_expected_layer1 = nullptr;
  TH1F *fmv_observed_layer2 = nullptr;
  TH1F *fmv_expected_layer2 = nullptr;
  TH1F *fmv_observed_track_strict_layer1 = nullptr;
  TH1F *fmv_expected_track_strict_layer1 = nullptr;
  TH1F *fmv_observed_track_strict_layer2 = nullptr;
  TH1F *fmv_expected_track_strict_layer2 = nullptr;
  TH1F *fmv_observed_track_loose_layer1 = nullptr;
  TH1F *fmv_expected_track_loose_layer1 = nullptr;
  TH1F *fmv_observed_track_loose_layer2 = nullptr;
  TH1F *fmv_expected_track_loose_layer2 = nullptr;

  TH1F *track_diff_x_Layer1 = nullptr;
  TH1F *track_diff_y_Layer1 = nullptr;
  TH1F *track_diff_x_strict_Layer1 = nullptr;
  TH1F *track_diff_y_strict_Layer1 = nullptr;
  TH1F *track_diff_x_loose_Layer1 = nullptr;
  TH1F *track_diff_y_loose_Layer1 = nullptr;
  TH1F *track_diff_x_Layer2 = nullptr;
  TH1F *track_diff_y_Layer2 = nullptr;
  TH1F *track_diff_x_strict_Layer2 = nullptr;
  TH1F *track_diff_y_strict_Layer2 = nullptr;
  TH1F *track_diff_x_loose_Layer2 = nullptr;
  TH1F *track_diff_y_loose_Layer2 = nullptr;
  TH2F *track_diff_xy_Layer1 = nullptr;
  TH2F *track_diff_xy_strict_Layer1 = nullptr;
  TH2F *track_diff_xy_loose_Layer1 = nullptr;
  TH2F *track_diff_xy_Layer2 = nullptr;
  TH2F *track_diff_xy_strict_Layer2 = nullptr;
  TH2F *track_diff_xy_loose_Layer2 = nullptr;

  //Tank coincidence
  TH1F *time_diff_tank_Layer1 = nullptr;
  TH1F *time_diff_tank_Layer2 = nullptr;

  TH1F *fmv_tank_observed_layer1 = nullptr;
  TH1F *fmv_tank_expected_layer1 = nullptr;
  TH1F *fmv_tank_observed_layer2 = nullptr;
  TH1F *fmv_tank_expected_layer2 = nullptr;


  //Verbosity settings
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;

};


#endif
