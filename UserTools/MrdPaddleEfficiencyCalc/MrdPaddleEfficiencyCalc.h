#ifndef MrdPaddleEfficiencyCalc_H
#define MrdPaddleEfficiencyCalc_H

#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TROOT.h"
#include "TColor.h"

#include "Tool.h"


/**
 * \class MrdPaddleEfficiencyCalc
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M. Nieslony $
* $Date: 2020/02/25 10:44:00 $
* Contact: mnieslon@uni-mainz.de
*/
class MrdPaddleEfficiencyCalc: public Tool {


 public:

  MrdPaddleEfficiencyCalc(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:

  //Configuration variables
  std::string str_input;
  std::string str_output;
  std::string str_inactive;
  int verbosity;
  bool isData;
  bool layer_plots;
  std::string mc_chankey_path;
  std::map<int,int> chankeymap_MC_data;
  std::map<int,int> chankeymap_data_MC;

  //TFile variables
  TFile *inputfile = nullptr;
  TFile *outputfile = nullptr;

  //Histograms + canvas
  TH1D *eff_chankey = nullptr;
  TH2Poly *eff_top = nullptr;
  TH2Poly *eff_side = nullptr;
  TH2Poly *eff_top_side = nullptr;
  TH2D *eff_crate1 = nullptr;
  TH2D *eff_crate2 = nullptr;
  TCanvas *canvas_elec = nullptr;

  //Channel information (both geometric + electronic)
  int num_slots = 24;
  int num_channels = 32;
  std::map<int,int> map_ch_orientation, map_ch_half, map_ch_side;
  std::map<int,double> map_ch_xmin, map_ch_xmax, map_ch_ymin, map_ch_ymax, map_ch_zmin, map_ch_zmax, map_ch_Crate, map_ch_Slot, map_ch_Channel;
  std::vector<int> active_slots_cr1, active_slots_cr2;
  std::vector<int> inactive_ch_crate1, inactive_slot_crate1, inactive_ch_crate2, inactive_slot_crate2;
  int active_channel[2][24]={{0}};
  std::map<int,std::vector<int>> MRDChannelNumToCrateSpaceMap;
  Geometry *geom = nullptr;

  //Layer-by-layer information
  int channels_per_layer[11] = {26,30,26,34,26,26,26,30,26,30,26};
  int channels_start[11] = {26,52,82,108,142,168,194,220,250,276,306};
  int bins_start[11] = {1,27,57,83,117,143,169,195,225,251,281};
  double extents[11] = {1.318,1.146,1.318,1.299,1.318,1.318,1.318,1.521,1.318,1.521,1.318};
  double zmin[11] = {1.6798,1.8009,1.929,2.0501,2.1782,2.2993,2.4204,2.5415,2.6626,2.7837,2.9048};
  double zmax[11] = {1.6858,1.8139,1.935,2.0631,2.1842,2.3053,2.4264,2.5475,2.6686,2.7897,2.9108};

  //Visualization variables
  double enlargeBoxes = 0.01;
  double shiftSecRow = 0.03;
  std::vector<TBox*> vector_box_inactive;

  //Color palettes (kBird color palette)
  const int n_colors=255;
  double alpha=1.;    //make colors opaque, not transparent
  Int_t Bird_palette[255];
  Int_t Bird_Idx;
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

  //Verbosity levels
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;
  int vv_debug = 4;

};


#endif
