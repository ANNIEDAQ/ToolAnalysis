#include "BeamDecoder.h"

//Definitions local to this source file
namespace {

  constexpr double BOGUS_DOUBLE = -9999.;

  std::string make_time_string(uint64_t ms_since_epoch) {
    time_t s_since_epoch = ms_since_epoch / THOUSAND;
    std::string time_string = asctime(gmtime(&s_since_epoch));
    time_string.pop_back(); // remove the trailing \n from the time string
    return time_string;
  }

}

BeamDecoder::BeamDecoder():Tool(),
beam_db_store_(false, BOOST_STORE_MULTIEVENT_FORMAT)
{}



bool BeamDecoder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Retrieve device names needed for beam quality checks from the
  // configuration file
  // 
  // The "first toroid" should be the one farther upstream from the target
  // (E:TOR860 in the current configuration)
  m_variables.Get("FirstToroid", first_toroid);
  m_variables.Get("SecondToroid", second_toroid);
  m_variables.Get("HornCurrentDevice", horn_current_device);

  BeamStatusMap = new std::map<uint64_t,BeamStatus>;

  first_entry = true;

  m_data->CStore.Set("NewBeamDataAvailable",false);

  return true;
}


bool BeamDecoder::Execute(){

  Log("BeamDecoder: First_entry: "+std::to_string(first_entry),4,verbosity_);
  if (first_entry){
    this->initialise_beam_db();
    first_entry = false;
  }

  Log("BeamDecoder: Check if new CTC data available",4,verbosity_);
  //Check if there is new trigger data available
  bool get_ctc = m_data->CStore.Get("NewCTCDataAvailable",NewCTCDataAvailable);
  Log("BeamDecoder: NewCTCDataAvailable: "+std::to_string(NewCTCDataAvailable),4,verbosity_);
  if (!get_ctc){
    Log("BeamDecoder tool: Did not find NewCTCDataAvailable entry in the CStore",0,verbosity_);
    return false;
  }
  
  if (NewCTCDataAvailable){
    m_data->CStore.Get("TimeToTriggerWordMap",TimeToTriggerWordMap);
    for (std::pair<uint64_t,std::vector<uint32_t>> apair : *TimeToTriggerWordMap){
      uint64_t current_timestamp = apair.first;
      std::vector<uint32_t> current_triggerwords = apair.second;
      //Check which minibuffer label we have
      MinibufferLabel current_minibuffer = MinibufferLabel::Unknown;
      if (std::find(current_triggerwords.begin(),current_triggerwords.end(),5)!=current_triggerwords.end()) current_minibuffer = MinibufferLabel::Beam;
      else if (std::find(current_triggerwords.begin(),current_triggerwords.end(),31)!=current_triggerwords.end()) current_minibuffer = MinibufferLabel::LED;
      else if (std::find(current_triggerwords.begin(),current_triggerwords.end(),35)!=current_triggerwords.end()) current_minibuffer = MinibufferLabel::AmBe;
      else if (std::find(current_triggerwords.begin(),current_triggerwords.end(),36)!=current_triggerwords.end()) current_minibuffer = MinibufferLabel::MRDCR; 
      bool is_beam = false;
      //Only add timestamps that are not yet contained in the BeamStatusMap
      if (BeamStatusMap->count(current_timestamp) == 0){
        auto temp_beam_status = get_beam_status(current_timestamp, current_minibuffer); 
        BeamStatusMap->emplace(current_timestamp,temp_beam_status);
        m_data->CStore.Set("NewBeamDataAvailable",true);
      }
    }
  }

  Log("BeamDecoder tool: Set BeamStatus map in CStore",4,verbosity_);
  m_data->CStore.Set("BeamStatusMap",BeamStatusMap);
  Log("BeamDecoder tool: BeamStatusmap->size(): "+std::to_string(BeamStatusMap->size()),4,verbosity_);

  return true;
}


bool BeamDecoder::Finalise(){

  return true;
}

bool BeamDecoder::initialise_beam_db(){

  m_variables.Get("verbosity", verbosity_);

  //Get information about run number
  Log("BeamDecoder tool: Get Postgress information",2,verbosity_);
  Store Postgress;
  int run_nr;
  m_data->CStore.Get("RunInfoPostgress",Postgress);
  Postgress.Get("RunNumber",run_nr);
  std::stringstream ss_db_filename;
  ss_db_filename << run_nr << "_beamdb";
  std::string db_filename = ss_db_filename.str();

  Log("BeamDecoder tool: Check that beam database file exists",2,verbosity_);
  // Check that the beam database file exists using a dummy std::ifstream
 std::ifstream dummy_in_file(db_filename);
  if ( !dummy_in_file.good() ) {
    Log("BeamDecoder tool: Error: Could not open the beam database file \"" + db_filename
     + '\"', 0, verbosity_);
    return false;
  }
  dummy_in_file.close();

  Log("BeamDecoder tool: Initialise BoostStore",2,verbosity_);
  beam_db_store_.Initialise( db_filename );
  Log("BeamDecoder tool: Printing BoostStore",2,verbosity_);
  beam_db_store_.Print(false);

  Log("BeamDecoder tool: Done initialising .. Get BeamDBIndex",4,verbosity_);
  bool got_index = beam_db_store_.Header->Get("BeamDBIndex", beam_db_index_);

  Log("BeamDecoder tool: Got the beam db index",4,verbosity_);
  if ( !got_index ) {
    Log("BeamDecoder tool: Error: Could not find the BeamDBIndex entry in the beam database"
      " BoostStore header stored in the file \"" + db_filename + '\"', 0,
      verbosity_);
    return false;
  }

  bool got_start = beam_db_store_.Header->Get("StartMillisecondsSinceEpoch",
    start_ms_since_epoch_);

  if ( !got_start ) {
    Log("BeamDecoder tool: Error: Could not find the StartMillisecondsSinceEpoch entry in the"
      " beam database BoostStore header stored in the file \"" + db_filename
      + '\"', 0, verbosity_);
    return false;
  }

  bool got_end = beam_db_store_.Header->Get("EndMillisecondsSinceEpoch",
    end_ms_since_epoch_);

  if ( !got_end ) {
    Log("BeamDecoder tool: Error: Could not find the EndMillisecondsSinceEpoch entry in the"
      " beam database BoostStore header stored in the file \"" + db_filename
      + '\"', 0, verbosity_);
    return false;
  }

  if ( !beam_db_store_.Has("BeamDB") ) {
    Log("BeamDecoder tool: Error: Missing BeamDB key in the beam database BoostStore stored"
      " in the file \"" + db_filename
      + '\"', 0, verbosity_);
    return false;
  }

  Log("BeamDecoder tool: Loaded beam database entries for times between "
    + make_time_string(start_ms_since_epoch_) + " and "
    + make_time_string(end_ms_since_epoch_), 1, verbosity_);

  return true;

}


BeamStatus BeamDecoder::get_beam_status(uint64_t ns_since_epoch,
  MinibufferLabel mb_label){

  //Note: Minibuffer nomenclature is from phaseI, but not super important
  //for the BeamStatus labeling --> keep the old classes (Beam/NonBeam Minibuffer)

  // Return immediately if the current minibuffer does not correspond to a beam
  // trigger (we don't want to double-count POT using self-trigger or
  // cosmic trigger minibuffers, for example)
  
  if ( mb_label != MinibufferLabel::Beam ) {
    return BeamStatus(TimeClass( ns_since_epoch ), 0.,
      BeamCondition::NonBeamMinibuffer);
  }

  static std::map<std::string, std::map<uint64_t, BeamDataPoint> >
    beam_data;
  static int current_beam_db_entry = 0;
  static bool got_first_beam_entry = false;

  if ( !got_first_beam_entry ) {
    beam_db_store_.GetEntry(current_beam_db_entry);
    beam_db_store_.Get("BeamDB", beam_data);
    got_first_beam_entry = true;
  }

  // The beam database uses timestamps with ms precision
  uint64_t ms_since_epoch = ns_since_epoch / MILLION;

  Log("Finding beam status information for "
    + make_time_string(ms_since_epoch), 2, verbosity_);

  // Find the beam database entry that contains POT information for the
  // moment of interest
  
  auto iter = std::find_if(beam_db_index_.cbegin(),
    beam_db_index_.cend(),
    [ms_since_epoch](const std::pair<int, std::pair<uint64_t,
      uint64_t> >& pair) -> bool {
      uint64_t start_ms = pair.second.first;
      uint64_t end_ms = pair.second.second;
      return ms_since_epoch >= start_ms && ms_since_epoch <= end_ms;
    }
  );

  // If a suitable entry could not be found, then complain and return
  // a BeamStatus object that indicates that the data were missing
 
  bool found_pot_entry = ( iter != beam_db_index_.cend() );

  if ( !found_pot_entry ) {
    Log("BeamDecoder tool: WARNING: unable to find a suitable entry for "
      + make_time_string(ms_since_epoch) + " (" + std::to_string(ms_since_epoch)
      + " ms since the Unix epoch) in the beam database file", 0, verbosity_);
    return BeamStatus( TimeClass(ns_since_epoch), 0., BeamCondition::Missing );
  }

  // If we need to load a new entry from the beam database, do so.
  // Avoid loading a new entry if you don't have to (the maps stored in
  // each entry are fairly large)   

  int new_entry_number = iter->first;
  if ( new_entry_number != current_beam_db_entry ) {

    beam_db_store_.GetEntry(new_entry_number);
    beam_db_store_.Get("BeamDB", beam_data);

    current_beam_db_entry = new_entry_number;
  }

  // Temporary storage for this function's return value
  BeamStatus beam_status;

  // Retrieve the beam monitoring measurements for the current event, and
  // calculate a POT value
  try {

      // Retrieve a measurement for each of the monitoring devices included
      // in the bundle from the Intensity Frontier beam database
      for (const auto& device_pair : beam_data) {
        const std::string& device_name = device_pair.first;

      // Obtain a map of measurements for the current device.
      // Keys are timestamps (in ms since the Unix epoch),
      // while values are BeamDataPoint objects.
      const std::map<uint64_t, BeamDataPoint>& measurement_map
        = device_pair.second;

      // Find the measurement entry with the closest time to that requested by
      // the user, and add it to the BeamStatus object that will be returned
      std::map<uint64_t, BeamDataPoint>::const_iterator
      low = measurement_map.lower_bound(ms_since_epoch);

      if ( low == measurement_map.cend() ) {
        Log("BeamDecoder tool: WARNING: IF beam database did not have any information"
          " for device " + device_name + " at " + std::to_string(ms_since_epoch)
          + " ms after the Unix epoch (" + make_time_string(ms_since_epoch)
          + ')', 1, verbosity_);
      }

      else if ( low == measurement_map.cbegin() ) {
        beam_status.add_measurement(device_name, low->first, low->second );
      }

      // We're between two time values, so we need to figure out which is
      // closest to the value requested by the user
      else {

        std::map<uint64_t, BeamDataPoint>::const_iterator
          prev = low;
        --prev;

        if ((ms_since_epoch - prev->first) < (low->first - ms_since_epoch)) {
          beam_status.add_measurement(device_name, prev->first, prev->second);
        }
        else {
          beam_status.add_measurement(device_name, low->first, low->second);
        }
      }
    }

    if ( verbosity_ >= 4 ) {
      Log("BeamDecoder tool: Beam monitoring measurements", 4, verbosity_);
      for (const auto& pair : beam_status.data()) {
        Log(pair.first + " at " + std::to_string(pair.second.first)
          + " ms since epoch: " + std::to_string(pair.second.second.value)
          + " " + pair.second.second.unit, 4, verbosity_);
      }
    }

    // Now that we've retrieved all of the needed measurements, calculate a POT
    // value and a timestamp for the BeamStatus object as a whole.
    //
    // We're catching exceptions, so if one of the required devices doesn't
    // exist, things will be handled properly.

    const auto& first_toroid_pair = beam_status.data().at(first_toroid);
    double pot_first_toroid = first_toroid_pair.second.value;
    int64_t ms_since_epoch_first_toroid = first_toroid_pair.first;

    const auto& second_toroid_pair = beam_status.data().at(second_toroid);
    double pot_second_toroid = second_toroid_pair.second.value;
    int64_t ms_since_epoch_second_toroid = second_toroid_pair.first;

    // Get the horn current measurement
    const auto& horn_pair = beam_status.data().at(horn_current_device);
    double peak_horn_current = horn_pair.second.value; // kA
    int64_t ms_since_epoch_horn = horn_pair.first;

    // TODO: perhaps add something more sophisticated here than just asking for
    // a POT value and a time from a single device.
    beam_status.set_pot( pot_second_toroid );
    beam_status.set_time( TimeClass(ms_since_epoch_second_toroid * MILLION) );

    // Wait to set the beam condition until we've applied quality cuts

    // Get beam quality cut parameters
    double pot_min = BOGUS_DOUBLE;
    double pot_max = BOGUS_DOUBLE;
    m_variables.Get("CutPOTMax", pot_max);
    m_variables.Get("CutPOTMin", pot_min);

    double horn_current_min = BOGUS_DOUBLE;
    double horn_current_max = BOGUS_DOUBLE;
    m_variables.Get("CutPeakHornCurrentMax", horn_current_max);
    m_variables.Get("CutPeakHornCurrentMin", horn_current_min);

    double toroid_tol = BOGUS_DOUBLE;
    m_variables.Get("CutToroidAgreement", toroid_tol);

    int64_t t_tol = BOGUS_INT;
    m_variables.Get("CutTimestampAgreement", t_tol);

    // TODO: add beam targeting efficiency cut
    BeamCondition bc = BeamCondition::Ok;

    // POT cut
    beam_status.add_cut("POT in range", (beam_status.pot() >= pot_min)
      && (beam_status.pot() <= pot_max));
    
    // Peak horn current cut
    beam_status.add_cut("peak horn current in range",
      (peak_horn_current >= horn_current_min)
      && (peak_horn_current <= horn_current_max));

    // Toroid disagreement cut
    double tor_diff_frac = 2.*std::abs( pot_second_toroid - pot_first_toroid )
      / ( pot_second_toroid + pot_first_toroid );
    beam_status.add_cut("toroids agree", tor_diff_frac <= toroid_tol);

    // Timestamp cut. Make sure measurements within the timing tolerance
    // were available for all required devices
    int64_t ms_since_epoch = static_cast<int64_t>( ns_since_epoch / MILLION );
    beam_status.add_cut("timestamps agree",
      ( std::abs(ms_since_epoch_horn - ms_since_epoch) <= t_tol )
      &&  ( std::abs(ms_since_epoch_first_toroid - ms_since_epoch) <= t_tol )
      &&  ( std::abs(ms_since_epoch_second_toroid - ms_since_epoch) <= t_tol ));

    Log("BeamDecoder tool: ANNIE DAQ ms since epoch = " + std::to_string(ms_since_epoch), 4,
      verbosity_);
    Log(horn_current_device + " ms since epoch = "
      + std::to_string(ms_since_epoch_horn), 4, verbosity_);
    Log(first_toroid + " ms since epoch = "
      + std::to_string(ms_since_epoch_first_toroid), 4, verbosity_);
    Log(second_toroid + " ms since epoch = "
      + std::to_string(ms_since_epoch_second_toroid), 4, verbosity_);

    // Flag the beam spill as "bad" if it failed any of the cuts above
    if ( !beam_status.passed_all_cuts() )
    {
      bc = BeamCondition::Bad;
      Log("BeamDecoder tool: WARNING: bad beam spill", 1, verbosity_);
    }

    if ( verbosity_ > 3 ) {
      for ( const auto& pair : beam_status.cuts() ) {
        std::string temp_string;
        if (pair.second) temp_string = "PASSED ";
        else temp_string = "FAILED ";
        Log(temp_string + pair.first + " beam quality cut", 4, verbosity_);
      }
    }

    beam_status.set_condition( bc );
  }

  catch (const std::exception& e) {
    Log("BeamDecoder tool: WARNING: problem encountered while querying IF beam"
      " database: " + std::string( e.what() ), 0, verbosity_);

    // There was a problem, so create a BeamStatus object that assumes that
    // the data were missing for some reason.
    return BeamStatus( TimeClass(ns_since_epoch), 0., BeamCondition::Missing );
  }

  return beam_status;
}
