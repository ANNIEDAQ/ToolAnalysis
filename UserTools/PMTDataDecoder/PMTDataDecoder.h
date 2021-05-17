#ifndef PMTDataDecoder_H
#define PMTDataDecoder_H

#include <string>
#include <iostream>
#include <bitset>
#include <deque>

#include "Tool.h"
#include "CardData.h"
#include "TriggerData.h"
#include "BoostStore.h"
#include "Store.h"

#include <boost/algorithm/string.hpp>

/**
 * \class PMTDataDecoder
 *
 This tool is used to decode the binary raw data and construct the ANNIEEvent
 boost stores.
 
 *
* $Author: Teal Pershing $
* $Date: 2019/05/28 10:44:00 $
* Contact: tjpershing@ucdavis.edu 
*/


struct DecodedFrame{
  bool has_recordheader;
  uint32_t frameheader;
  std::vector<uint16_t> samples;
  std::vector<int> recordheader_starts; //Holds indices where a record header starts in samples
  ~DecodedFrame(){
  }
};



class PMTDataDecoder: public Tool {


 public:

  PMTDataDecoder(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  std::vector<DecodedFrame> DecodeFrames(std::vector<uint32_t> bank);

  void ParseFrame(int CardID, DecodedFrame DF);
  void ParseSyncFrame(int CardID, DecodedFrame DF);
  void ParseRecordHeader(int CardID, int ChannelID, std::vector<uint16_t> RH);
  void StoreFinishedWaveform(int CardID, int ChannelID);
  void AddSamplesToWaveBank(int CardID, int ChannelID, std::vector<uint16_t> WaveSlice);
  bool CheckIfCardNextInSequence(CardData aCardData);
  void BuildReadyEvents();


 private:

  std::string InputFile;
  std::string Mode;

  bool NewWavesBuilt;
  int ADCCountsToBuild;  //If a finished wave doesn't have this many ADC counts at least, don't add it for building
  int EntriesPerExecute;
  int PMTDEntryNum = 0; 
  int FileNum = 0;
  int CurrentRunNum;
  int CurrentSubrunNum;
  std::string CurrentFile = "NONE";
  int RECORD_HEADER_SAMPLENUMS = 8;

  //A Record header's first 48-bit word has least significant bits of 0xFFF000.  
  //So, check for this structure within the 12-bit samples being parsed.
  int RECORD_HEADER_LABELPART1 = 0x000;
  int RECORD_HEADER_LABELPART2 = 0xFFF;
  int SYNCFRAME_HEADERID = 10;
  //A record header is made of two 48-bit words, each word in little endian order.  The
  //beginning of the first word has the 0x000 of the Record Header.  Given each 12-bit
  //chunk is stored in a 16-bit word, you want to grab the 7 samples right of 0x000 to
  //get the entire header.
  unsigned int SAMPLES_RIGHTOF_000 = 7;

  BoostStore* PMTData;
  std::vector<CardData>* Cdata = nullptr;
  std::vector<CardData> Cdata_old;

  //Counter used to track the number of entries processed in a PMT file
  int NumPMTDataProcessed = 0;

  //Maps for keeping track of what SequenceID is next for each Card
  std::map<int, int> SequenceMap;  //Key is CardID, Value is next SequenceID 

  //Vector to keep track of fifo errors (type I, type II, for monitoring tools)
  std::vector<int> fifo1;
  std::vector<int> fifo2;
  int FIFOstate = 0;


  //Maps used in decoding frames; specifically, holds record header and record waveform info
  std::map<std::vector<int>, uint64_t> TriggerTimeBank;  //Key: {cardID, channelID}. Value: trigger time associated with wave in WaveBank 
  std::map<std::vector<int>, std::vector<uint16_t>> WaveBank;  //Key: {cardID, channelID}. Value: Waveform being built for this Card and ADC Channel. 
  std::map<int,std::vector<uint64_t>> SyncCounters; //Key: cardID.  Value: vector of sync counters filled in the order they arrive.
 
  //Extra maps used for FIFO overflow info and TimestampsFromTheFuture
  std::map<uint64_t, std::map<std::vector<int>, int> >* FIFOPMTWaves = nullptr;
  std::map<uint64_t,std::map<std::vector<int>,uint64_t>>* TimestampsFromTheFuture = nullptr;
  uint64_t LastGoodTimestamp = 0;

  //Maps that store completed waveforms from cards
  std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > >* FinishedPMTWaves;  //Key: {MTCTime}, value: "WaveMap" with key (CardID,ChannelID), value FinishedWaveform
  std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > > FinishedPMTWaves_old;  //Key: {MTCTime}, value: "WaveMap" with key (CardID,ChannelID), value FinishedWaveform
  std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > > CStorePMTWaves;  //Key: {MTCTime}, value: "WaveMap" with key (CardID,ChannelID), value FinishedWaveform
  bool NewWaveBuilt;
  std::map<uint64_t, std::map<std::vector<int>, std::vector<uint16_t> > > CStoreTankEvents;  //Key: {MTCTime}, value: "WaveMap" with key (CardID,ChannelID), value FinishedWaveform

  // Notes whether DAQ is in lock step running
  // Number of PMTs that must be found in a WaveSet to build the event
  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  int vv_debug=4;
  std::string logmessage;

};


#endif
