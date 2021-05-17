#pragma once

// standard library includes
#include <string>
#include <vector>

// ToolAnalysis includes
#include "Tool.h"

class LoadANNIEEvent: public Tool {

  public:

    LoadANNIEEvent();
    bool Initialise(std::string configfile, DataModel& data);
    bool Execute();
    bool Finalise();

  private:
  
    int v_error = 0;
    int v_warning = 1;
    int v_message = 2;
    int v_debug = 3;
    int vv_debug = 4;

  protected:

    /// @brief Integer code that determines the level of logging to show in
    /// the output
    int verbosity_;

    /// @brief FileFormat of input file (combined ANNIEEvent+OrphanStore or separate stores)
    std::string FileFormat;

    /// @brief Boolean whether OrphanStore should be accessed
    bool load_orphan_store;

    /// @brief Vector of filenames for each of the input files
    std::vector<std::string> input_filenames_;

    /// @brief Vector of filenames for each of the orphan input files
    std::vector<std::string> input_filenames_orphan_;
    
    /// @brief The index of the current entry in the ANNIEEvent store
    size_t current_entry_;

    /// @brief Event offset if one wants to ignore the first offset_evnum events
    int offset_evnum;

    /// @brief The index of the current file in this list of input files
    size_t current_file_;

    /// @brief The total number of ANNIEEvent entries in the current file
    size_t total_entries_in_file_;

    /// @brief The total number of OrphanStore entries in the current file
    size_t total_orphans_in_file_;

    /// @brief Flag indicating whether we need to load a new file
    bool need_new_file_;

    bool load_orphanstore;

    std::stringstream logmessage;
};
