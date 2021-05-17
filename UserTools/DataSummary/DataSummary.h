/* vim:set noexpandtab tabstop=4 wrap */
#ifndef DataSummary_H
#define DataSummary_H

#include <string>
#include <iostream>

#include "Tool.h"

class TFile;
class TTree;
#include "TH2F.h"


/**
* \class DataSummary
*
* A tool to make summary plots of recent run information
*
* $Author: M. O'Flaherty $
* $Date: 2020/07/10 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class DataSummary: public Tool {
	
	public:
	
	DataSummary();   ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();  ///< Execute function used to perform Tool purpose.
	bool Finalise(); ///< Finalise function used to clean up resources.
	
	private:
	
	std::string DataPath;         // input file path
	std::string InputFilePattern; // input file pattern (regex, with some caveats)
	std::string InputFilePatternOrphan; // input file pattern for orphan files (regex, with some caveats)
	std::string FileList;
	int StartRun;                 // limit the range of matched runs to analyse
	int StartSubRun;
	int StartPart;
	int EndRun;
	int EndSubRun;
	int EndPart;
	std::string FileFormat;	
	// list of files to process, in order. key is RunSubrunPart, zero-padded to ensure order
	std::map<std::string,std::string> filelist;
	std::map<std::string,std::string>::iterator nextfile; // the next file to be loaded
	std::map<std::string,std::string> filelist_orphan;
	std::map<std::string,std::string>::iterator nextfile_orphan; // the next file to be loaded

	// BoostStore to read data into
	BoostStore* ProcessedFileStore=nullptr;
	BoostStore* ANNIEEvent=nullptr;
	BoostStore* OrphanStore=nullptr;
	uint64_t globalentry = 0;  // entry number across all processed BoostStores
	uint64_t localentry = 0;   // entry number in the current BoostStore
	uint64_t localentries = 0; // number of entries in the current BoostStore
	uint64_t globalorphan = 0; // same as above but for orphan store
	uint64_t localorphan = 0;
	uint64_t localorphans = 0;
	
	// output file
	std::string OutputFileDir;
	std::string OutputFileName;
	TFile* outfile=nullptr;
	TTree* outtree=nullptr;
	TTree* outtree2=nullptr;
	
	// output branch variables
	// run level
	uint32_t RunNumber;               // RunNumber
	uint32_t SubrunNumber;            // SubrunNumber
	uint32_t PartNumber;              // <need to get this from the filename>
	int RunType;                      // RunType
	std::string RunTypeString;        // convert RunType from int
	uint64_t RunStartTime;            // RunStartTime
	int run,subrun,part;              // extracted from filename
	
	// used for making plots on time axis; first and last timestamp
	double t0;
	double tn;
	
	// event level
	uint32_t EventNumber;             // EventNumber - or use global? local?
	uint64_t CTCtimestamp;            // CTCTimestamp
	uint64_t PMTtimestamp;            // EventTimeTank
	TimeClass mrd_timeclass;          // EventTime in ANNIEEvent
	uint64_t MRDtimestamp;            // <need to convert from mrd_timeclass>
	std::map<std::string,int> MRDLoopbackTDC; // map of loopback type string to TDCVal
	uint64_t BeamLoopbackTimestamp;   // <need to calculate from MRDLoopbackTDC TDCVal>
	uint64_t CosmicLoopbackTimestamp; // <need to calculate from MRDLoopbackTDC TDCVal>
	uint32_t TriggerWord;             // TriggerWord
	std::string TriggerTypeString;    // <need to convert with jonathan's file>
	uint8_t SystemsPresent;           // <int of what readouts were present>
	uint8_t LoopbacksPresent;         // <int of what loopbacks were present>
	TimeClass EventTime;              // <should we make some official TimeClass object? e.g. from CTC?>
	std::map<std::string,bool> datastreams; // DataStreams (which subdetectors saw a coincident timestamp?)
	bool data_ctc;			// Does the event include CTC info?
	bool data_tank;			// Does the event include Tank data?
	bool data_mrd;			// Does the event include MRD data?
	bool window_is_extended = false;   // Was the window an extended window?
	int size_of_window = 2000;	//The size of the acquisition window (in us)
	std::vector<int> window_sizes;	//ADC acquisition sizes for all recorded channels
	std::vector<int> window_chkeys;	//Chankeys for ADC acquisition window sizes
	int CTCWordExtended;  //Trigword information about extended windows

	ULong64_t CTCtimestamp_tree;
	ULong64_t PMTtimestamp_tree;
	ULong64_t MRDtimestamp_tree;	
	ULong64_t BeamLoopbackTimestamp_tree;
	ULong64_t CosmicLoopbackTimestamp_tree;
	double CTCtimestamp_double;
	double PMTtimestamp_double;
	double MRDtimestamp_double;
	double CTCtimestamp_sec;
	double PMTtimestamp_sec;
	double MRDtimestamp_sec;
	ULong64_t orphantimestamp_tree;
	double orphantimestamp_double;
	double orphantimestamp_sec;
	bool trigword_ext;
	bool trigword_ext_cc;
	bool trigword_ext_nc;

	// orphan info
	std::string orphantype;           // EventType
	uint64_t orphantimestamp;         // Timestamp
	std::string orphancause;          // Reason
	int orphannumwaves;		  // Number of waveforms in orphan event
	std::vector<unsigned long> orphanchankeys;	//Chankeys that had a waveform in the orphaned event
	std::vector<int> orphanchankeys_int;	//Chankeys that had a waveform in the orphaned event
	double orphanmintdiff;		  // Time difference to closest CTC timestamp	
	std::vector<std::vector<int>> orphanchannels; 	//Electronics channels that were hit & orphaned
	std::vector<int> orphanchannels_combined;	//Combined electronics channel
	int orphantrigword;		// Trigword of orphaned CTC timestamp

	// functions
	int ScanForFiles(std::string inputdir, std::string filepattern, std::string filepattern_orphan);
	int ReadInFileList(std::string filelist_user);
	bool LoadNextANNIEEventEntry();
	bool LoadNextOrphanStoreEntry();
	bool LoadNextFile();
	bool CreateOutputFile();
	bool CreatePlots();
	bool AddRatePlots(int nbins);
	bool AddEventTypePlots();
	bool AddTDiffPlots();
	
	//TODO: reprocess current data with CTC build type. Standardize filenames including BuildType and .bs extension
	
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
	
};


#endif
