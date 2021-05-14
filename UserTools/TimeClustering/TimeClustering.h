#ifndef TimeClustering_H
#define TimeClustering_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TFile.h"
#include "TObjectTable.h"

class TH1D;
class TH2D;
class TApplication;
class TCanvas;

/**
* \class TimeClustering
*
* This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class TimeClustering: public Tool {
	
	public:
	TimeClustering(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	//Configuration variables
	int minimumdigits=4;                        // a cluster must have at least 4 hits
	double maxsubeventduration=30;              // if all hits within this time, just one subevent
	double minimum_subevent_timeseparation=30;  // minimum empty time to delimit subevents
	bool MakeMrdDigitTimePlot=false;
	bool MakeSingleEventPlots;
	bool LaunchTApplication;
	bool isData;
	std::string output_rootfile;
	std::string file_chankeymap;
	int evnum;
	bool ModifiedTDCData;	
	std::string TimeShiftChannels;

	//Information about which channels are time-shifted
	std::vector<std::pair<unsigned long,unsigned long>> shifted_channels;

	//ANNIEEvent data
	std::map<unsigned long,vector<Hit>>* TDCData;
	std::map<unsigned long,vector<MCHit>>* TDCData_MC;
	Geometry *geom = nullptr;
	
	// From the CStore, for converting WCSim TubeId t channelkey
	std::map<unsigned long,int> channelkey_to_mrdpmtid;
	std::map<unsigned long,int> channelkey_to_faccpmtid;
	std::map<int,unsigned long> mrdpmtid_to_channelkey;	

	// Cluster properties
	std::vector<double> mrddigittimesthisevent;
	std::vector<int> mrddigitpmtsthisevent;
	std::vector<unsigned long> mrddigitchankeysthisevent;
	std::vector<double> mrddigitchargesthisevent;
	std::vector<std::vector<int>> MrdTimeClusters;
	std::vector<std::vector<double>> MrdTimeClusters_Times;
	std::vector<std::vector<double>> MrdTimeClusters_Charges;
	
	// Histograms storing information about the MRD cluster times
	TH1D* mrddigitts_cosmic_cluster = nullptr;
	TH1D* mrddigitts_beam_cluster = nullptr;
	TH1D* mrddigitts_noloopback_cluster = nullptr;
	TH1D* mrddigitts_cluster = nullptr;
	TH1D* mrddigitts_cosmic = nullptr;
	TH1D* mrddigitts_beam = nullptr;
	TH1D* mrddigitts_noloopback = nullptr;
	TH1D* mrddigitts = nullptr;
	TH1D *mrddigitts_cluster_single = nullptr;
	TH1D *mrddigitts_single = nullptr;
	TH1D *mrddigitts_horizontal = nullptr;
	TH1D *mrddigitts_vertical = nullptr;
	TH1D *hist_chankey = nullptr;
	TH1D *hist_chankey_cluster = nullptr;
	TH2D *hist_chankey_time = nullptr;	
	TH2D *hist_chankey_time_cluster = nullptr;	
	TH1D *hist_chankey_multi = nullptr;
	TFile* mrddigitts_file = nullptr;

	//TApplication-related variables
	TCanvas* timeClusterCanvas=nullptr;
	TApplication* rootTApp=nullptr;
	double canvwidth;
	double canvheight;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
};


#endif
