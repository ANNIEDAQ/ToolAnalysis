#ifndef LoadGeometry_H
#define LoadGeometry_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Geometry.h"
#include <boost/algorithm/string.hpp>

class LoadGeometry: public Tool {


 public:

  LoadGeometry();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  bool FileExists(std::string name);
  std::string GetLegendLine(std::string name);
  void InitializeGeometry();
  void LoadFACCMRDDetectors();
  void LoadLAPPDs();
  bool  ParseMRDDataEntry(std::vector<std::string> SpecLine,
                               std::vector<std::string> MRDLegendEntries);

  bool  ParseLAPPDDataEntry(std::vector<std::string> SpecLine,
                              std::vector<std::string> LAPPDLegendEntries);

  void LoadTankPMTDetectors();
  void LoadAuxiliaryChannels();
  bool ParseTankPMTDataEntry(std::vector<std::string> SpecLine,
                               std::vector<std::string> TankPMTLegendEntries);
  bool ParseAuxChannelDataEntry(std::vector<std::string> SpecLine,
                               std::vector<std::string> AuxChannelLegendEntries);
  void LoadTankPMTGains();

  Geometry* AnnieGeometry;


 private:

  int detector_num_store = 0;
  int counter = 0;
  Detector* adet;
  int LAPPD_channel_count;
  std::string fFACCMRDGeoFile;
  std::string fTankPMTGeoFile;
  std::string fTankPMTGainFile;
  std::string fAuxChannelFile;
  std::string fLAPPDGeoFile;
  std::string fDetectorGeoFile;

  //Labels used in Geometry files to mark the legend and data entries
  std::string LegendLineLabel = "LEGEND_LINE";
  std::string DataStartLineLabel = "DATA_START";
  std::string DataEndLineLabel = "DATA_END";

  //Map of channel number to electronics map entry
  std::map<std::vector<int>,int>* MRDCrateSpaceToChannelNumMap;
  std::map<int,std::vector<int>>* MRDChannelNumToCrateSpaceMap;
  std::map<std::vector<int>,int>* TankPMTCrateSpaceToChannelNumMap;
  std::map<std::vector<int>,int>* AuxCrateSpaceToChannelNumMap;
  std::map<int,std::vector<int>>* ChannelNumToTankPMTCrateSpaceMap;
  std::map<int,std::vector<int>>* AuxChannelNumToCrateSpaceMap;
  std::map<int,double>* ChannelNumToTankPMTSPEChargeMap;
  std::map<int,std::string>* AuxChannelNumToTypeMap;
  std::map<std::vector<unsigned int>,int>* LAPPDCrateSpaceToChannelNumMap;

  //Vector of strings indicating variables of interest and their data types in
  //The MRD file.  Used in the LoadFACCMRDDetectors() method
  std::vector<std::string> MRDIntegerValues{"detector_num","channel_num","detector_system","orientation","layer","side","num",
                                   "rack","TDC_slot","TDC_channel","discrim_slot","discrim_ch",
                                   "patch_panel_row","patch_panel_col","amp_slot","amp_channel",
                                   "hv_crate","hv_slot","hv_channel","nominal_HV","polarity"};
  std::vector<std::string> MRDDoubleValues{"x_center","y_center","z_center","x_width","y_width","z_width"};
  std::vector<std::string> MRDStringValues{"PMT_type","cable_label","paddle_label","notes"};


  //Vector of strings indicating variables of interest and their data types in
  //The LAPPD file.  Used in the LoadLAPPDs() method
  std::vector<std::string> LAPPDIntegerValues{"detector_num","channel_strip_side","channel_strip_num"};
  std::vector<std::string> LAPPDDoubleValues{"detector_position_x","detector_position_y","detector_position_z",
                                             "detector_direction_x","detector_direction_y","detector_direction_z",
                                             "channel_position_x","channel_position_y","channel_position_z"};
  std::vector<std::string> LAPPDStringValues{"detector_type","detector_status","channel_status"};
  std::vector<std::string> LAPPDUnIntValues{"channel_signal_crate","channel_signal_card","channel_signal_channel",
                                            "channel_level2_crate","channel_level2_card","channel_level2_channel",
                                            "channel_hv_crate","channel_hv_card","channel_hv_channel","channel_num"};

  //Vector of strings indicating variables of interest and their data types in
  //The TankPMT file.  Used in the LoadTankPMTDetectors() method
  std::vector<std::string> TankPMTIntegerValues{"detector_num","channel_num","panel_number","sb_num","sb_channel",
                                   "signal_crate","signal_slot","signal_channel",
                                   "mt_crate","mt_slot","mt_channel",
                                   "hv_crate","hv_slot","hv_channel","nominal_HV"};
  std::vector<std::string> TankPMTDoubleValues{"x_pos","y_pos","z_pos","x_dir","y_dir","z_dir"};
  std::vector<std::string> TankPMTStringValues{"detector_tank_location","PMT_type","cable_label","detector_status","notes"};

  //Vector of strings indicating variables of interest and their data types in
  //The AuxChannel file.  Used in the LoadAuxiliaryChannels() method
  std::vector<std::string> AuxChannelIntegerValues{"channel_num","signal_crate","signal_slot","signal_channel"};
  std::vector<std::string> AuxChannelStringValues{"channel_type","notes"};

  //verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;


};


#endif
