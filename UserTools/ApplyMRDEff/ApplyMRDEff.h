#ifndef ApplyMRDEff_H
#define ApplyMRDEff_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Geometry.h"
#include "Detector.h"
#include "Hit.h"

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TTree.h"

/**
 * \class ApplyMRDEff
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M.Nieslony $
* $Date: 2020/10/30 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/
class ApplyMRDEff: public Tool {


 public:

  ApplyMRDEff(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:

  //Config variables
  bool verbosity;
  std::string mrd_eff_file;
  std::string outputfile;
  bool debug_plots;
  std::string file_chankeymap;
  bool use_file_eff;
  std::string file_eff;

  //ANNIEEvent Store variables
  std::map<unsigned long,std::vector<MCHit>> *TDCData = nullptr;
  std::map<unsigned long,std::vector<MCHit>> *TDCData_mod = nullptr;
  Geometry *geom = nullptr;

  //CStore variables
  std::map<unsigned long, int> channelkey_to_mrdpmtid_data;
  std::map<int, unsigned long> mrdpmtid_to_channelkey_data; 
  std::map<unsigned long, int> channelkey_to_mrdpmtid;
  std::map<int, unsigned long> mrdpmtid_to_channelkey; 

  //ROOT variables
  TFile *debug_file = nullptr;
  TH1F *mrd_expected = nullptr;
  TH1F *mrd_observed = nullptr;
  TH1F *mrd_eff = nullptr;
  TRandom3 *rnd = nullptr;
  TTree *tree = nullptr;
  std::vector<double> *vec_random = nullptr;
  std::vector<unsigned long> *dropped_ch = nullptr;
  std::vector<unsigned long> *dropped_ch_mc = nullptr;
  std::vector<double> *dropped_ch_time = nullptr;
  int evnum;
  TFile *File_Eff = nullptr;
  TTree *Tree_Eff = nullptr;
  std::vector<double> *vec_random_read = nullptr;
  std::vector<unsigned long> *dropped_ch_read = nullptr;
  std::vector<unsigned long> *dropped_ch_mc_read = nullptr;
  std::vector<double> *dropped_ch_time_read = nullptr;
  int evnum_read;

  //Efficiency values
  std::map<unsigned long,double> map_eff;

  //Verbosity variables
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;

};


#endif
