#ifndef RawLoadToRoot_H
#define RawLoadToRoot_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "ChannelKey.h"
#include "Waveform.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

class RawLoadToRoot: public Tool {


 public:

  RawLoadToRoot();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

   TFile* treeoutput;
   TTree* rawtotree;
   TString treeoutfile;
   TH1D* theWaveformHist;

   int lbound;
   int ubound;
   int onoffswitch;
   int RelEventNumber;
   int EventNumber;
   int RunNumber;
   int SubRunNumber;
   std::string key;
   std::string value;

   std::map<ChannelKey, std::vector<Waveform<unsigned short> > >
    RawADCData;

    //ChannelKey ck(subdetector::ADC, pmt_id);
    //std::vector<Waveform<unsigned short> > raw_waveforms;



};


#endif
