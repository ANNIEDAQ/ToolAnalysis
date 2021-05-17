/* vim:set noexpandtab tabstop=4 wrap */
#ifndef PrintANNIEEvent_H
#define PrintANNIEEvent_H

#include <string>
#include <iostream>
#include <sstream>

#include "Tool.h"
#include "Particle.h"
#include "Waveform.h"
#include "Hit.h"
#include "LAPPDHit.h"
#include "TriggerClass.h"
#include "TimeClass.h"
#include "BeamStatus.h"
#include "BeamStatusClass.h"
#include "ADCPulse.h"

class PrintANNIEEvent: public Tool {
	
	public:
	
	PrintANNIEEvent();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	int verbose=1;
	int get_ok;
	std::string inputfile;
	unsigned long NumEvents;
	bool is_mc;	
	bool has_raw;
	bool first_event;

	// contents of ANNIEEvent
	uint32_t RunNumber;
	uint32_t SubrunNumber;
	uint32_t EventNumber;

	//MC variables
	BeamStatusClass* BeamStatusMC = nullptr;
	std::vector<MCParticle>* MCParticles=nullptr;
	std::vector<Particle>* RecoParticles=nullptr;
        std::map<unsigned long,std::vector<MCHit>>* MCHits=nullptr;
	std::map<unsigned long,std::vector<MCLAPPDHit>>* MCLAPPDHits=nullptr;
	std::map<unsigned long,std::vector<MCHit>>* MCTDCData=nullptr;
	bool MCFlag;
	TimeClass* EventTime=nullptr;
	uint64_t MCEventNum;
	uint16_t MCTriggernum;
	std::string MCFile;
	std::vector<TriggerClass>* TriggerData = nullptr;
	
	//Data variables
	std::string MRDTriggertype;
	std::map<std::string,int> MRDLoopbackTDC;
	uint64_t EventTimeMRD;
	uint64_t EventTimeTank;
	uint64_t CTCTimestamp;
	uint32_t TriggerWord;
	std::map<std::string, bool> DataStreams;
	BeamStatus BeamStatusData;
	int PartNumber;
	int TriggerExtended;
	TriggerClass TriggerDataData;
	std::map<unsigned long,std::vector<Hit>>* Hits=nullptr;
	std::map<unsigned long,std::vector<Hit>>* AuxHits=nullptr;
	std::map<unsigned long,std::vector<LAPPDHit>>* LAPPDHits=nullptr;
	std::map<unsigned long,std::vector<Hit>>* TDCData=nullptr;
	std::map<unsigned long,std::vector<Waveform<uint16_t>>> RawADCData;
	std::map<unsigned long,std::vector<Waveform<uint16_t>>> RawAuxADCData;
	std::map<unsigned long,std::vector<Waveform<uint16_t>>> RawLAPPDData;
	std::map<unsigned long,std::vector<Waveform<double>>> CalibratedADCData;
	std::map<unsigned long,std::vector<Waveform<double>>> CalibratedLAPPDData;
	std::map<unsigned long,std::vector<std::vector<ADCPulse>>> RecoADCData;
	std::map<unsigned long,std::vector<std::vector<ADCPulse>>> RecoAuxADCData;
	std::map<unsigned long,std::vector<int>> RawAcqSize;

	std::stringstream logmessage;

	int n_prompt;
	int n_ext;
	int n_ext_cc;
	int n_ext_nc;
};

#endif
