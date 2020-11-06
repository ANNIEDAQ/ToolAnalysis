#ifndef CARDDATA_H
#define CARDDATA_H

//#include <stdint.h>
#include <iostream>
#include <vector>
//#include <array>
#include <stdlib.h>
#include "zmq.hpp"
#include <SerialisableObject.h>

class CardData : public SerialisableObject{

  friend class boost::serialization::access;

 public:

  int CardID;
  int SequenceID;
  int FirmwareVersion;
  std::vector<uint32_t> Data;
  int FIFOstate; 
 
  ~CardData();

  bool Print(){return true;};

  void  Send(zmq::socket_t *socket, int flag=0);
  bool Receive(zmq::socket_t *socket);

 private:

  template <class Archive> void serialize(Archive& ar, const unsigned int version){

    ar & CardID;
    ar & SequenceID;
    ar & FirmwareVersion;
    ar & Data;
    ar & FIFOstate; 

  }
  
};


#endif
