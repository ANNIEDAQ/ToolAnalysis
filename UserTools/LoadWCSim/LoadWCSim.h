/* vim:set noexpandtab tabstop=4 wrap */
#ifndef LoadWCSim_H
#define LoadWCSim_H

#include <string>
#include <iostream>
//#include <boost/algorithm/string.hpp>  // for boost::algorithm::to_lower(string)

#include "Tool.h"
#include "TFile.h"
#include "TTree.h"
#include "wcsimT.h"
#include "Particle.h"
#include "Hit.h"
#include "Waveform.h"
#include "TriggerClass.h"
#include "Geometry.h"
#include "MRDspecs.hh"
#include "Detector.h"
#include "BeamStatus.h"

#include "TObjectTable.h"

namespace{
	//PMTs
	constexpr int ADC_CHANNELS_PER_CARD=4;
	constexpr int ADC_CARDS_PER_CRATE=20;
	constexpr int MT_CHANNELS_PER_CARD=4;
	constexpr int MT_CARDS_PER_CRATE=20;
	//LAPPDs
	constexpr int ACDC_CHANNELS_PER_CARD=30;
	constexpr int ACDC_CARDS_PER_CRATE=20;
	constexpr int ACC_CHANNELS_PER_CARD=8;
	constexpr int ACC_CARDS_PER_CRATE=20;
	//TDCs
	constexpr int TDC_CHANNELS_PER_CARD=32;
	constexpr int TDC_CARDS_PER_CRATE=6;
	//HV
	constexpr int CAEN_HV_CHANNELS_PER_CARD=16;
	constexpr int CAEN_HV_CARDS_PER_CRATE=10;
	constexpr int LECROY_HV_CHANNELS_PER_CARD=16;
	constexpr int LECROY_HV_CARDS_PER_CRATE=16;
	constexpr int LAPPD_HV_CHANNELS_PER_CARD=4; // XXX ??? XXX
	constexpr int LAPPD_HV_CARDS_PER_CRATE=10;  // XXX ??? XXX
}

class LoadWCSim: public Tool {
	
	public:
	
	LoadWCSim();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	std::vector<int> LoadPMTMask(std::string path_to_pmtmask);
	
	private:
	
	// variables from config file
	/////////////////////////////
	int verbosity=1;
	int HistoricTriggeroffset;
	int use_smeared_digit_time;   // digit_time = (T): first photon smeared time, (F): first photon true time
	int LappdNumStrips;           // number of Channels per LAPPD
	double LappdStripLength;      // [mm] for calculating relative x position for dual-ended readout
	double LappdStripSeparation;  // [mm] for calculating relative y position of each stripline
	
	// WCSim variables
	//////////////////
	//TFile* file;
	//TTree* wcsimtree;
	wcsimT* WCSimEntry; // from makeclass
	WCSimRootTrigger* atrigt, *atrigm, *atrigv;
	WCSimRootTrigger* firsttrigt, *firsttrigm, *firsttrigv;
	WCSimRootGeom* wcsimrootgeom;
	WCSimRootOptions* wcsimrootopts;
	int WCSimVersion;   // WCSim version
	
	// Misc Others
	//////////////
	long NumEvents;
	long MaxEntries;   // maximum MC entries to process before stopping the loop. -1 for indefinite
	int numtankpmts;
	int numlappds;
	int nummrdpmts;
	int numvetopmts;
	
	//Masked PMT IDs
	std::vector<int> masked_ids;
        std::string PMTMask;

	// For constructing ToolChain Geometry
	//////////////////////////////////////
	Geometry* ConstructToolChainGeometry();
	std::map<int,unsigned long> lappd_tubeid_to_detectorkey;
	std::map<int,unsigned long> pmt_tubeid_to_channelkey;
	std::map<int,unsigned long> mrd_tubeid_to_channelkey;
	std::map<int,unsigned long> facc_tubeid_to_channelkey;
	// inverse
	std::map<unsigned long,int> detectorkey_to_lappdid;
	std::map<unsigned long,int> channelkey_to_pmtid;
	std::map<unsigned long,int> channelkey_to_mrdpmtid;
	std::map<unsigned long,int> channelkey_to_faccpmtid;
	// alternatively? Better? Save the parentage in each MCHit. Each MCHit will contain
	// the index of it's parent MCParticle in the MCParticles vector
	std::map<int,int>* trackid_to_mcparticleindex=nullptr;
	std::vector<int> GetHitParentIds(WCSimRootCherenkovDigiHit* digihit, WCSimRootTrigger* firstTrig);
	std::map<int,int> timeArrayOffsetMap;
	void BuildTimeArrayOffsetMap(WCSimRootTrigger* firstTrig);
	int triggers_event;	

	// FIXME temporary!! remove me when we have a better way to get FACC paddle origins
	std::vector<double> facc_paddle_yorigins{-198.699875000, -167.999875000, -137.299875000, -106.599875000, -75.899875000, -45.199875000, -14.499875000, 16.200125000, 46.900125000, 77.600125000, 108.300125000, 139.000125000, 169.700125000, -198.064875000, -167.364875000, -136.664875000, -105.964875000, -75.264875000, -44.564875000, -13.864875000, 16.835125000, 47.535125000, 78.235125000, 108.935125000, 139.635125000, 170.335125000}; // taken from geofile.txt
	
	////////////////
	// things that will be filled into the store from this WCSim file.
	// note: filling everything in the right format is a complex process;
	// just do a simple filling here: this will be properly handled by the conversion
	// from WCSim to Raw and the proper RawReader Tools
	// bool MCFlag=true; 
	std::string MCFile;
	uint64_t MCEventNum;
	uint16_t MCTriggernum;
	uint32_t RunNumber;
	uint32_t SubrunNumber;
	uint32_t EventNumber;    // will need to be tracked separately, since we flatten triggers
	TimeClass* EventTime;
	uint64_t RunStartUser;   // sets the run start time, since WCSim doesn't have one
	TimeClass RunStartTime;  // as set from user
	uint64_t EventTimeNs;
	std::vector<MCParticle>* MCParticles;
	std::map<unsigned long,std::vector<MCHit>>* TDCData;
	std::map<unsigned long,std::vector<MCHit>>* MCHits;
	std::vector<TriggerClass>* TriggerData;
	BeamStatusClass* BeamStatus;
	
	int primarymuonindex;
	
	// additional info
	void MakeParticleToPmtMap(WCSimRootTrigger* thisTrig, WCSimRootTrigger* firstTrig, std::map<int,std::map<unsigned long,double>>* ParticleId_to_DigitIds, std::map<int,double>* ChargeFromParticleId, std::map<int,unsigned long> tubeid_to_channelkey);
	std::map<int,std::map<unsigned long,double>>* ParticleId_to_TankTubeIds = nullptr;
	std::map<int,std::map<unsigned long,double>>* ParticleId_to_MrdTubeIds = nullptr;
	std::map<int,std::map<unsigned long,double>>* ParticleId_to_VetoTubeIds = nullptr;
	std::map<int,double>* ParticleId_to_TankCharge = nullptr;
	std::map<int,double>* ParticleId_to_MrdCharge = nullptr;
	std::map<int,double>* ParticleId_to_VetoCharge = nullptr;
	std::map<unsigned long,int> Mrd_Chankey_Layer;
	bool mrd_firstlayer, mrd_lastlayer;
	std::string Triggertype;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
	
};


#endif
