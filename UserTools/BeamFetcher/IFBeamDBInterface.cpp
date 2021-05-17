// standard library includes
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sstream>

// ToolAnalysis includes
#include "ANNIEconstants.h"
#include "IFBeamDBInterface.h"

IFBeamDBInterface::IFBeamDBInterface()
{
  fCurl = curl_easy_init();
  if (!fCurl) throw std::runtime_error("IFBeamDBInterface failed to"
    " initialize libcurl");
}

IFBeamDBInterface::~IFBeamDBInterface()
{
  curl_easy_cleanup(fCurl);
}

const IFBeamDBInterface& IFBeamDBInterface::Instance() {

  // Create the IFBeamDBInterface using a static variable. This ensures that
  // the singleton instance is only created once.
  static std::unique_ptr<IFBeamDBInterface> the_instance(
    new IFBeamDBInterface());

  // Return a reference to the singleton instance
  return *the_instance;
}

std::map<std::string, std::map<uint64_t, BeamDataPoint> >
  IFBeamDBInterface::QueryBeamDB(uint64_t t0, uint64_t t1)
  const
{
  // Temporary storage for the response from the IF beam database
  std::string response;

  int code = QueryBeamDB(t0, t1, response);

  if (code != CURLE_OK) throw std::runtime_error("Error accessing"
    " IF beam database. Please check your internet connection.");

  auto beam_data = ParseDBResponse(response);

  return beam_data;
}

int IFBeamDBInterface::QueryBeamDB(uint64_t t0,
  uint64_t t1, std::string& response_string) const
{
  if (!fCurl) {
    throw std::runtime_error("IFBeamDBInterface::QueryBeamDB() called"
      " without inititalizing libcurl");
    return -1;
  }

  constexpr char BNB_URL_START[] = "http://ifb-data.fnal.gov:8089/ifbeam/"
    "data/data?e=e%2C1d&b=BNBBPMTOR&f=csv&tz=&action=Show+device&t0=";

  std::stringstream url_stream;

  url_stream << BNB_URL_START;
  url_stream << std::fixed << std::setprecision(3) << (t0 - 1)/1000.;
  url_stream << "&t1=" << (t1 + 1)/1000.;

  curl_easy_setopt(fCurl, CURLOPT_URL, url_stream.str().c_str());

  response_string.clear();

  curl_easy_setopt(fCurl, CURLOPT_WRITEFUNCTION,
    static_cast<size_t(*)(char*, size_t, size_t, std::string*)>(
      [](char* ptr, size_t size,
        size_t num_members, std::string* data) -> size_t
      {
        data->append(ptr, size * num_members);
        return size * num_members;
      }
    )
  );
  curl_easy_setopt(fCurl, CURLOPT_WRITEDATA, &response_string);

  int code = curl_easy_perform(fCurl);

  // Check the HTTP response code from the IF beam database server. If
  // it's not 200, then throw an exception (something went wrong).
  // For more information about the possible HTTP status codes, see
  // this Wikipedia article: http://tinyurl.com/8yqvhwf
  long http_response_code;
  constexpr long HTTP_OK = 200;
  curl_easy_getinfo(fCurl, CURLINFO_RESPONSE_CODE, &http_response_code);
  if (http_response_code != HTTP_OK) {
    throw std::runtime_error("HTTP error (code "
      + std::to_string(http_response_code) + ") encountered while querying"
      " the IF beam database");
  }

  return code;
}

// Here's some documentation for some of the parameters stored in the beam
// database. It's taken from the MicroBooNE operations wiki:
// http://tinyurl.com/z3c4mxs
//
// The status page shows the present reading of beamline instrumentation. All
// of this data is being stored to IF beam Database. The "IF Beam DB
// dashboard":http://dbweb4.fnal.gov:8080/ifbeam/app/BNBDash/index provides
// another view of beam data. Some of it is redundant to the status page, but
// it verifies that the data is being stored in the database. At present the
// page shows following devices:
//   * TOR860, TOR875 - two toroids in BNB measuring beam intensity. TOR860 is
//     at the beginning of the beamline, and TOR875 is at the end.
//   * THCURR - horn current
//   * HWTOUT - horn water temperature coming out of the horn.
//   * BTJT2 - target temperature
//   * HP875, VP875 - beam horizontal and vertical positions at the 875
//     location, about 4.5 m upstream of the target center.
//   * HPTG1, VPTG1 - beam horizontal and vertical positions immediately
//     (about 2.5 m) upstream of the target center.
//   * HPTG2, VPTG2 - beam horizontal and vertical positions more immediately
//     (about 1.5 m) upstream of the target center.
//   * Because there are no optics between H/VP875 and H/VPTG2, the movements
//     on these monitors should scale with the difference in distances.
//   * BTJT2 - target air cooling temperature. Four RTD measure the return
//     temperature of the cooling air. This is the one closest to the target.
//   * BTH2T2 - target air cooling temperature. This is the temperature of the
//     air going into the horn.

std::map<std::string, std::map<uint64_t, BeamDataPoint> >
  IFBeamDBInterface::ParseDBResponse(const std::string& response) const
{
  // Create an empty map to store the parsed data.
  std::map<std::string, std::map<uint64_t,
    BeamDataPoint> > beam_data;

  // Use a stringstream to parse the data
  std::istringstream response_stream(response);

  // Temporary storage for the comma-separated fields in the database response
  // string
  uint64_t time_stamp;
  std::string data_type;
  std::string unit;
  double value;

  // Skip the first line (which gives textual column headers)
  std::getline(response_stream, unit, '\n');

  // Parse each line of the response and load the map with the parsed data
  while (response_stream >> time_stamp) {
    response_stream.ignore(1);
    std::getline(response_stream, data_type, ',');
    std::getline(response_stream, unit, ',');
    response_stream >> value;

    // If there were any input problems, give up
    if (!response_stream) break;

    // Otherwise, load the new data into the map
    beam_data[data_type][time_stamp] = BeamDataPoint(value, unit);
  }

  PostprocessParsedResponse(beam_data);

  return beam_data;
}

void IFBeamDBInterface::PostprocessParsedResponse(std::map<std::string,
  std::map<uint64_t, BeamDataPoint> >& parsed_response) const
{
  // Perform postprocessing for the given device names
  // The only two devices currently selected for postprocessing are toroids
  // used to measure the protons-on-target (POT)
  constexpr std::array<const char*, 2> devices = { "E:TOR860", "E:TOR875" };

  std::stringstream ss;
  double value;

  // Check whether each device name is found in the parsed response map
  for (const auto& device : devices) {
    auto iter = parsed_response.find(device);

    // If it is, append the unit (which is known for these devices to be of the
    // form "E12") to the end of the numerical value, then convert the entire
    // string back to a number. Replace the given unit with "POT" (for "protons
    // on target") and update the map with the new value and unit.
    if (iter != parsed_response.end()) {
      for (auto& pair : iter->second) {

        // Clear and reset the string stream
        ss.str("");
        ss.clear();

        // Do the conversion to a numerical value using the exponent from the
        // "E12" unit
        ss << pair.second.value << pair.second.unit;
        ss >> value;

        // Update the map with the new value and unit
        pair.second.value = value;
        pair.second.unit = "POT";
      }
    } else {
      std::cout <<"IFBeamDBInterface: Did not find device name : "<<device<<std::endl;
    }
  }
}

BeamStatus IFBeamDBInterface::GetBeamStatus(uint64_t time) const
{
  constexpr uint64_t TIME_OFFSET = 5000ull; // ms

  // The exact time given by the user might not appear in the database, so
  // look for entries on the interval [time - TIME_OFFSET, time + TIME_OFFSET]
  uint64_t t0 = time - TIME_OFFSET;
  uint64_t t1 = time + TIME_OFFSET;

  // Ask the beam database for information about the time window of interest.
  // If there are any problems, then return a bogus BeamStatus object.
  try {

    auto beam_data = QueryBeamDB(t0, t1);

    // TODO: remove hard-coded device name here
    // Get protons-on-target (POT) information from the parsed data
    const auto& pot_map = beam_data.at("E:TOR875");

    // Find the POT entry with the closest time to that requested by the
    // user, and use it to create the BeamStatus object that will be returned
    const auto low = pot_map.lower_bound(time);

    if (low == pot_map.cend()) {
      std::cerr << "WARNING: IF beam database did not have any information"
        " for " << time << " ms after the Unix epoch\n";

      return BeamStatus(TimeClass(time * MILLION), 0., BeamCondition::Missing);
    }

    else if (low == pot_map.cbegin()) {
      // TODO: add checks for BeamCondition::Bad here?
      return BeamStatus(TimeClass( low->first * MILLION ), low->second.value,
        BeamCondition::Ok);
    }

    // We're between two time values, so we need to figure out which is closest
    // to the value requested by the user
    else {

      auto prev = low;
      --prev;

      if ((time - prev->first) < (low->first - time)) {
        // TODO: add checks for BeamCondition::Bad here?
        return BeamStatus(TimeClass( prev->first * MILLION ),
          prev->second.value, BeamCondition::Ok);
      }
      else {
        // TODO: add checks for BeamCondition::Bad here?
        return BeamStatus(TimeClass( low->first * MILLION ), low->second.value,
          BeamCondition::Ok);
      }
    }

  }

  catch (const std::exception& e) {
    std::cerr << "WARNING: problem encountered while querying IF beam"
      " database:\n  " << e.what() << '\n';

    // Return a default-constructed BeamStatus object since there was a
    // problem. The fOk member will be set to false by default, which will
    // flag the object as problematic.
    return BeamStatus(TimeClass( time ), 0., BeamCondition::Missing);
  }

}
