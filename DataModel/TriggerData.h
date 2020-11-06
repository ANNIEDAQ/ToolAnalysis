#ifndef TRIGGERDATA_H
#define TRIGGERDATA_H

#include <stdint.h>
#include <iostream>
#include <vector>
#include <zmq.hpp>
#include <SerialisableObject.h>

//#include <boost/iostreams/stream.hpp>
//#include <boost/iostreams/device/back_inserter.hpp>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/serialization/vector.hpp>


class TriggerData : public SerialisableObject{

  friend class boost::serialization::access;

 public:
  
  //  friend class boost::serialization::access;
  int FirmwareVersion;
  int SequenceID;
  int EventSize;
  std::vector<uint16_t> EventIDs;
  std::vector<uint64_t> EventTimes;  // units of nanoseconds.
  //int TriggerSize; 
  int TimeStampSize;
  //std::vector<uint32_t> TriggerMasks;
  //std::vector<uint32_t> TriggerCounters;
  std::vector<uint32_t> TimeStampData; 
  int FIFOOverflow;
  int DriverOverflow;
  ~TriggerData();
  bool Print(){return true;};
  
  void  Send(zmq::socket_t *socket, int flag=0);
  bool Receive(zmq::socket_t *socket);        
  
 private:

  template <class Archive> void serialize(Archive& ar, const unsigned int version){

      ar & FirmwareVersion;
      ar & SequenceID;
      ar & EventSize;
      //ar & TriggerSize;
      ar & TimeStampSize;
      ar & FIFOOverflow;
      ar & DriverOverflow;
      ar & EventIDs;
      ar & EventTimes;
      //ar & TriggerMasks;
      //ar & TriggerCounters;
      ar & TimeStampData; 
    }
    
};

//template <class T>
//std::string serialize_save(const T& obj) {
//  std::string outstr;
//  boost::iostreams::back_insert_device<std::string> inserter(outstr);
//  boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
//  boost::archive::text_oarchive oa(s);
//  oa << obj;
//  s.flush();
//  return outstr;
//}
//
//template <class T>
//void deserialize_load(const std::string data, size_t size) {
//  T obj;
//  boost::iostreams::basic_array_source<char> device(data.data(), size);
//  boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
//  boost::archive::text_iarchive ia(s);
//  ia >> obj;
//}

#endif
