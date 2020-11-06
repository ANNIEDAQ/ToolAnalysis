#include "LAPPDParseACC.h"
#include "Geometry.h"


LAPPDParseACC::LAPPDParseACC():Tool(){}

using namespace std;



bool LAPPDParseACC::Initialise(std::string configfile, DataModel &data){


  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  //datafilename
  string path;
  m_variables.Get("lappd_data_filepath", path); // does not have a "/" at the end
  string name;
  m_variables.Get("lappd_data_filename", name); // does not have a suffix ".acdc" 
  string filebase = path + "/" + name;

  //open the ACDC filestream
  string line; //dummy line variable
  string datapath = filebase + ".acdc";
  dfs.open(datapath.c_str());

  if(!dfs.is_open())
  {
    cout << "Could not load acdc data file." << endl;
    return false;
  }
  //process the header already
  getline(dfs, line);

  data_stream_pos = dfs.tellg();
  


  //open metadata filestream
  string metapath = filebase + ".meta";
  mfs.open(metapath.c_str());
  if(!mfs.is_open())
  {
    cout << "Could not load acdc data file." << endl;
    return false;
  }
  //process the header already
  getline(mfs, meta_header);


  //event number handling. create event 0
  //at initialization stage (here)
  _event_no = -1;

  return true;
}


bool LAPPDParseACC::Execute(){

  
  //this tool loads one event into the store. 
  //get the event number from the last loaded event
  //and make sure to compare while parsing the
  //raw ACDC files. 
  _event_no++;

  //find the geometry based on a string
  //given by the config file
  string geoname;
  m_variables.Get("geometry_name", geoname);
  string storename; 
  m_variables.Get("store_name", storename);
  bool geomfound;
  Geometry* geom;
  geomfound = m_data->Stores.at(storename)->Header->Get(geoname,geom);
  if(!geomfound)
  {
    cout << "Geometry " << geoname << " in store " << storename << " was not found while attempting to prase ACDC data" << endl;
    return false;
  }


  //isolate the LAPPDs and their channel lists
  //so that when we parse an event in the data file, 
  //we look for that board/channel and re-key the data
  //stream to match the geometry channel keys. 
  map<unsigned long, Channel>* all_lappd_channels = new map<unsigned long, Channel>; 
  GetAllLAPPDChannels(geom, all_lappd_channels);
  
  m_data->Stores.at(storename)->Set("AllLAPPDChannels", all_lappd_channels);

  //create raw lappd data structure
  map<unsigned long, Waveform<double>>* LAPPDWaveforms = new map<unsigned long, Waveform<double>>;


  int bo = 0, ch = 0; //event board and channel parsing vars
  int ev = _event_no;
  string line; //temp parsing var
  int line_counter = 0;
  //key of a certain channel from geometry 
  //associated with the presently parsed channel
  unsigned long geom_key = 0; 
  bool keyfound = false;

  //if you want better precision, take in calibration
  //data, either the ACDC makeLUT table or an output
  //from Whitmer's calibration code
  double adccounts_to_mv = 1.2*1000.0/4096.0; //mv

  //start of this loop is just after
  //the header of the datafile
  Waveform<double> tempwav;


  //in ensuring the data is parsed
  //event by event, we store the start
  //position of the event from the filestream.
  //put the dfs curser at the first line of this event
  dfs.seekg(data_stream_pos); 

  while(getline(dfs, line))
  {
    istringstream iss(line); //the current line in the file
    int temp_bit; //the current integer/bit in the line
    int char_count = -1; //counts the present integer/bit number in the line

    //top of this loop is the start of the line
    //in the data file. 
    while(iss >> temp_bit)
    {
      char_count++;

      //if we are on the event number
      if(char_count == 0)
      {
        ev = temp_bit;
        continue;
      }
      if(char_count == 1)
      {
        //this is the board number
        //not implemented yet
        bo = temp_bit;
        continue;
      }
      if(char_count == 2)
      {
        //this is the channel number
        ch = temp_bit;
        continue;
      }
      
      //if none of the above if statements fire
      //then we are in a sample (0 - 255) of ADC counts
      tempwav.PushSample(temp_bit*adccounts_to_mv);
    }

    //if this is no longer the event we want, 
    //then break out of the loop. if it is still
    //on the event, save this ifstream position
    if(ev == _event_no) data_stream_pos = dfs.tellg();
    else break;

    //find the unique channel key for this
    //board/channel combination given a 
    //geometry class. uses geom_key as a reference
    keyfound = FindChannelKeyFromACDCFile(bo, ch, geom_key, all_lappd_channels);
    //if the channel/board combo is not in the 
    //list of lappd channels in the geometry class,
    //move to the next channel. 
    if(!keyfound) continue;

    //otherwise, add this waveform to the raw lappd waveform
    //structure. 
    LAPPDWaveforms->insert(pair<unsigned int, Waveform<double>>(geom_key, tempwav));

    //reinitialize
    tempwav.ClearSamples();
    line_counter++;
  }

  m_data->Stores.at(storename)->Set("LAPPDWaveforms", LAPPDWaveforms);
  
  /* 

  //Please keep the following comment, it currently does not work and 
  //is not meant to work, but will help in writing the metadata parsing function

  //Metadata 
  //structure metadata[event][board]["key string"] = unsigned int 

  map<int, map<int, map<string, unsigned int>>> metadata;
  stringstream headers(meta_header);
  string header; //temp variable for the particular header column parsing the stringstream headers
  vector<string> header_vec; //a vector containing header column info
  bool first_loop = true; //flag to fill the header_vec only the first time around
  while(getline(mfs, line))
  {
    istringstream iss(line); //the current line in the file
    unsigned int temp_bit;
    int char_count = -1;
    while(iss >> temp_bit)
    {
      headers >> header;
      if(first_loop) header_vec.push_back(header);
      char_count++;

      //if we are on the event number
      if(char_count == 0)
      {
        event = temp_bit;
        continue;
      }
      if(char_count == 1)
      {
        //this is the board number
        //not implemented yet
        board = temp_bit;
        continue;
      }
      //if none of the above if statements fire
      //then we are in a metadata key
      metadata[event][board][header_vec.at(char_count)] = temp_bit;
    }
    first_loop = false;
  }

  m_data->Stores["ANNIEEvent"]->Header->Set("metaData", metadata);


  


  isLoaded = true;
  m_data->Stores["ANNIEEvent"]->Set("isLoaded", isLoaded); //have loaded the entire data file
  */

  return true;
}

//function takes in a board and ch number from
//an ACDC data file and loads a unique channel
//key number found in a map of channels taken 
//from a detector geometry in the Execute function.
//returns 0 if the channel is not found, 1 if the channel is found.
bool LAPPDParseACC::FindChannelKeyFromACDCFile(int acdc_board, int acdc_ch, unsigned long& key, map<unsigned long, Channel>* lappd_channels_from_geom)
{

  //the channel number and board number
  //from the detector class, loaded in the
  //geometry loading tool. 
  unsigned int det_board, det_ch;

  //loop through to find the matching
  //board and ch number. This can seriously
  //be optimized; write a search function to
  //make faster. 
  map<unsigned long, Channel>::iterator itCh;
  for(itCh = lappd_channels_from_geom->begin(); itCh != lappd_channels_from_geom->end(); ++itCh)
  {
    Channel thech = itCh->second;
    unsigned long k = itCh->first;
    det_board = thech.GetSignalCard();
    det_ch = thech.GetSignalChannel();
    //this is the golden find, return if found. 
    if((int)det_board == acdc_board and (int)det_ch == acdc_ch)
    {
      //change the key pointer to match the found channel
      key = k;
      return true;
    }
  }

  //if it never returns, then it doesnt exist in the geometry
  //leave the key reference unaltered.
  return false;
}


//isolate the LAPPDs and their channel lists
//so that when we parse an event in the data file, 
//we look for that board/channel and re-key the data
//stream to match the geometry channel keys. 
void LAPPDParseACC::GetAllLAPPDChannels(Geometry* geom, map<unsigned long, Channel>* all_lappd_channels)
{
  map<string, map<unsigned long,Detector*> >* AllDetectors = geom->GetDetectors();
  map<string, map<unsigned long,Detector*> >::iterator itGeom;
  for(itGeom = AllDetectors->begin(); itGeom != AllDetectors->end(); ++itGeom)
  {
    if(itGeom->first == "LAPPD")
    {
      map<unsigned long,Detector*> LAPPDDetectors = itGeom->second;
      map<unsigned long, Detector*>::iterator itDet;
      for(itDet = LAPPDDetectors.begin(); itDet != LAPPDDetectors.end(); ++itDet)
      {
        //here are the channel objects for this particular LAPPD
        map<unsigned long, Channel>* lappdchannels = itDet->second->GetChannels();
        //now loop through and insert into the all_lappd_channels 
        //map to make a cumulative map of all channels
        map<unsigned long, Channel>::iterator itCh;
        for(itCh = lappdchannels->begin(); itCh != lappdchannels->end(); ++itCh)
        {
          all_lappd_channels->insert(pair<unsigned long, Channel>(itCh->first, itCh->second));
        }
      }
    }
  }
}

bool LAPPDParseACC::Finalise(){
  dfs.close();
  mfs.close();

  return true;
}
