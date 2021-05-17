// Tool to download beam information from the Intensity Frontier
// beam database and save it in a BoostStore file for later retrieval
// by the BeamChecker tool
//
// Steven Gardiner <sjgardiner@ucdavis.edu>
#pragma once

// standard library includes
#include <iostream>
#include <string>

// Boost includes
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

// ToolAnalysis includes
#include "Tool.h"

class BeamFetcher: public Tool {

  public:

    BeamFetcher();
    bool Initialise(std::string configfile, DataModel& data);
    bool Execute();
    bool Finalise();

  protected:

    /// @brief The verbosity to use when printing logging messages
    /// @details A larger value corresponds to more verbose output
    int verbosity_;

    /// @brief Transient BoostStore used to save downloaded information
    /// from the Intensity Frontier beam database to disk
    BoostStore beam_db_store_;

    /// @brief Helper function that downloads and processes the beam
    /// information
    bool fetch_beam_data(uint64_t start_ms_since_epoch,
      uint64_t end_ms_since_epoch, uint64_t chunk_step_ms);

    /// @brief Name of the output file in which the beam database information
    /// will be saved
    std::string db_filename_;

    std::string start_date;
    std::string end_date;
    std::string start_timestamp;
    std::string end_timestamp;
    boost::posix_time::ptime current;
    std::string timestamp_mode;
    uint64_t start_ms_since_epoch;
    uint64_t end_ms_since_epoch;
    void ConvertDateToMSec(std::string start_str,std::string end_str,uint64_t &start_ms,uint64_t &end_ms);
    uint64_t TimeZoneShift;
    bool DaylightSavings;
    std::map<int,std::map<std::string,std::string>> RunInfoDB;
    int RunNumber;


};
