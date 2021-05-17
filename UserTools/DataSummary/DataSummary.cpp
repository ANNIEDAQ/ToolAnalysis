/* vim:set noexpandtab tabstop=4 wrap */
#include "DataSummary.h"
#include "RunTypeLabel.h"
#include "TrigTypeLabel.h"

#include <sys/types.h> // for stat() test to see if file or folder
#include <sys/stat.h>
#include <unistd.h>
#include <memory>
#include <regex>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"

DataSummary::DataSummary():Tool(){}

bool DataSummary::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); // loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// set some defaults
	verbosity=0;
	DataPath=".";
	InputFilePattern="ProcessedRawData";
	InputFilePatternOrphan="OrphanStore";
	StartRun=-1;
	StartSubRun=-1;
	StartPart=-1;
	EndRun=-1;
	EndSubRun=-1;
	EndPart=-1;
	OutputFileDir=".";
	OutputFileName="DataSummary.root";
	FileFormat="CombinedStore";	//Other option: SeparateStores
	FileList="None";	//Option to use filelist instead of regex matching
	
	// read the user's preferences
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("DataPath",DataPath);
	m_variables.Get("InputFilePattern",InputFilePattern);
	m_variables.Get("InputFilePatternOrphan",InputFilePatternOrphan);	//only needed if FileFormat=SeparateStores
	m_variables.Get("FileList",FileList);
	m_variables.Get("StartRun",StartRun);
	m_variables.Get("StartSubRun",StartSubRun);
	m_variables.Get("StartPart",StartPart);
	m_variables.Get("EndRun",EndRun);
	m_variables.Get("EndSubRun",EndSubRun);
	m_variables.Get("EndPart",EndPart);
	m_variables.Get("OutputFileDir",OutputFileDir);
	m_variables.Get("OutputFileName",OutputFileName);
	m_variables.Get("FileFormat",FileFormat);	

	//Scan for wrong file format configuration
	if (FileFormat != "SeparateStores" && FileFormat != "CombinedStore") {
		Log("DataSummary tool: FileFormat option "+FileFormat+" not supported. Use SeparateStores",v_error,verbosity);
		FileFormat = "SeparateStores";
	}

	// scan for matching input files
	int numfilesfound = 0;
	if (FileList == "None") numfilesfound=ScanForFiles(DataPath,InputFilePattern,InputFilePatternOrphan);
	else numfilesfound=ReadInFileList(FileList);
	if(numfilesfound==0){
		Log("No matching files found!",v_error,verbosity);
		return false;
	}

	std::cout << "List of matched filenames: "<<std::endl;
	for (std::map<std::string,std::string>::iterator it=filelist.begin(); it!=filelist.end(); it++){
		std::cout <<it->first<<": "<<it->second<<std::endl;
	}

	if (FileFormat == "SeparateStores"){
		for (std::map<std::string,std::string>::iterator it=filelist_orphan.begin(); it!=filelist_orphan.end(); it++){
			std::cout <<it->first<<": "<<it->second<<std::endl;
		}
	}

	
	// make the output file
	CreateOutputFile();
	
	return true;
}


bool DataSummary::Execute(){
	
	// load next entry. Return of false indicates end of files: StopLoop will be set.
	bool got_annieevent = LoadNextANNIEEventEntry();
	bool got_orphan = LoadNextOrphanStoreEntry();
	if(!got_annieevent && !got_orphan){
		std::cout <<"no annie or orphan event"<<std::endl;
		// i got nothin' ¯\_(ツ)_/¯
		return true;
	}
	
	if(got_annieevent){
		Log("DataSummary Tool: Getting information from ANNIEEvent",v_debug,verbosity);
		
		// variables we can directly retrieve
		ANNIEEvent->Get("EventNumber",EventNumber);
		ANNIEEvent->Get("EventTimeTank",PMTtimestamp);
		ANNIEEvent->Get("CTCTimestamp",CTCtimestamp);
		ANNIEEvent->Get("EventTimeMRD",mrd_timeclass);       // convertme to MRDtimestamp
		ANNIEEvent->Get("TriggerWord",TriggerWord);       // convert to TriggerTypeString
		ANNIEEvent->Get("MRDLoopbackTDC",MRDLoopbackTDC); // convert to LoopbackTimestamp values
		ANNIEEvent->Get("DataStreams",datastreams);  //DataStreams can be used to check which of the subdetectors is included in the data
	        ANNIEEvent->Get("TriggerExtended",CTCWordExtended);

        	// TODO optional sanity checks: consistency of RunNumber and other constants
	
		std::map<unsigned long, std::vector<Waveform<unsigned short>>> raw_waveform_map;
		bool has_raw = ANNIEEvent->Get("RawADCData",raw_waveform_map);
		
		window_is_extended = false;
		size_of_window = 2000;
		window_sizes.clear();
		window_chkeys.clear();
	
		if (has_raw){
			for (auto& temp_pair : raw_waveform_map) {
				const auto& achannel_key = temp_pair.first;
				auto& araw_waveforms = temp_pair.second;
				for (unsigned int i=0; i< (int) araw_waveforms.size(); i++){
					auto samples = araw_waveforms.at(i).GetSamples();
        				int size_sample = 2*samples->size();
					window_sizes.push_back(size_sample);
					window_chkeys.push_back(achannel_key);
					if (size_sample > size_of_window) {
						size_of_window = size_sample;
						window_is_extended = true;
					}
				}
			}
		}

		std::map<unsigned long, std::vector<int>> raw_acqsize_map;
		bool has_raw_acqsize = ANNIEEvent->Get("RawAcqSize",raw_acqsize_map);
		
		if (has_raw_acqsize && !has_raw){
			for (auto& temp_pair : raw_acqsize_map) {
				const auto& achannel_key = temp_pair.first;
				auto& araw_acqsize = temp_pair.second;
				for (unsigned int i=0; i< (int) araw_acqsize.size(); i++){
					int size_sample = 2*araw_acqsize.at(i);
					window_sizes.push_back(size_sample);
					window_chkeys.push_back(achannel_key);
					if (size_sample > size_of_window){
						size_of_window = size_sample;
						window_is_extended = true;
					}						
				}
			}
		}

		if (!has_raw && !has_raw_acqsize) {
   	 		Log("DataSummary tool: Did not find RawADCData or RawAcqSize in ANNIEEvent! Abort",v_error,verbosity);
    			/*return false;*/
  		}

		// calculated variables
		MRDtimestamp = mrd_timeclass.GetNs();

		PMTtimestamp_tree = (ULong64_t) PMTtimestamp;
		CTCtimestamp_tree = (ULong64_t) CTCtimestamp;
		MRDtimestamp_tree = (ULong64_t) MRDtimestamp;
		PMTtimestamp_double = (double) PMTtimestamp_tree;
		CTCtimestamp_double = (double) CTCtimestamp_tree;
		MRDtimestamp_double = (double) MRDtimestamp_tree;
		PMTtimestamp_sec = (PMTtimestamp_double)/(1.E9);
		MRDtimestamp_sec = (MRDtimestamp_double)/(1.E9);
		CTCtimestamp_sec = (CTCtimestamp_double)/(1.E9);
		trigword_ext = false;
		trigword_ext_cc = false;
		trigword_ext_nc = false;
		if (CTCWordExtended > 0) trigword_ext = true;
		if (CTCWordExtended == 1) trigword_ext_cc = true;
		if (CTCWordExtended == 2) trigword_ext_nc = true;

		std::cout <<"PMTtimestamp_tree: "<<PMTtimestamp_tree<<", MRDTimestamp_tree: "<<MRDtimestamp_tree<<std::endl;
		// extract out the TDC vals
		int beamloopbackTDCticks = MRDLoopbackTDC.at("BeamLoopbackTDC");
		int cosmicloopbackTDCticks = MRDLoopbackTDC.at("CosmicLoopbackTDC");
		// convert to ns
		BeamLoopbackTimestamp = 4000. - 4.*(double)beamloopbackTDCticks;
		CosmicLoopbackTimestamp = 4000. - 4.*(double)cosmicloopbackTDCticks;
	
		BeamLoopbackTimestamp_tree = (ULong64_t) BeamLoopbackTimestamp;
		CosmicLoopbackTimestamp_tree = (ULong64_t) CosmicLoopbackTimestamp;
		// FIXME to convert loopback TDC ticks to UTC time, we need to know the time difference
		// between the loopback signal entering the TDC card and the loopback CTC event associated with it
		// for now, i dunno, just neglect this
		BeamLoopbackTimestamp += CTCtimestamp;
		CosmicLoopbackTimestamp += CTCtimestamp;
		// convert int word to enum class TriggerType, then convert enum class to string
		TriggerTypeString = trigtype_to_string(TrigTypeEnum(TriggerWord));
		SystemsPresent=0;
		if(CTCtimestamp!=0) SystemsPresent |= 1;
		if(PMTtimestamp!=0) SystemsPresent |= 2;
		if(MRDtimestamp!=0) SystemsPresent |= 4;
		LoopbacksPresent=0;
		if(BeamLoopbackTimestamp!=0) LoopbacksPresent |= 1;
		if(CosmicLoopbackTimestamp!=0) LoopbacksPresent |= 2;
		
		data_ctc = false;
		data_tank = false;
		data_mrd = false;
		if (datastreams["CTC"]==1) data_ctc = true;
		if (datastreams["Tank"]==1) data_tank = true;
		if (datastreams["MRD"]==1) data_mrd = true;

		Log("DataSummary Tool: Filling Event tree",v_debug,verbosity);
		outtree->Fill();
	}
	
	if(got_orphan){
		Log("DataSummary Tool: Getting information from OrphanStore",v_debug,verbosity);
		// right now this is just a straight translation from BoostStore to ROOT =/
		OrphanStore->Get("EventType",orphantype);
		OrphanStore->Get("Timestamp",orphantimestamp);
		OrphanStore->Get("TriggerWord",orphantrigword);
		OrphanStore->Get("Reason",orphancause);
		OrphanStore->Get("NumWaves",orphannumwaves);
		OrphanStore->Get("WaveformChannels",orphanchannels);
		OrphanStore->Get("WaveformChankeys",orphanchankeys);
		OrphanStore->Get("MinTDiff",orphanmintdiff);	
	
		orphantimestamp_tree = (ULong64_t) orphantimestamp;
		orphantimestamp_double = (double) orphantimestamp_tree;
		orphantimestamp_sec = orphantimestamp_double/(1.E9);
		orphanchankeys_int.clear();
		orphanchannels_combined.clear();
		for (int i_vec=0; i_vec < (int) orphanchankeys.size(); i_vec++){
			orphanchankeys_int.push_back(int(orphanchankeys.at(i_vec)));
		}
		for (int i_vec=0; i_vec < (int) orphanchannels.size(); i_vec++){
			std::vector<int> orphanchannels_single = orphanchannels.at(i_vec);
			orphanchannels_combined.push_back(1000*orphanchannels_single.at(0)+10*orphanchannels_single.at(1)+orphanchannels_single.at(2));
		}

		Log("DataSummary Tool: Filling Orphan tree",v_debug,verbosity);
		outtree2->Fill();
	}
	
	return true;
}


bool DataSummary::Finalise(){
	
	// ensure the ttree is fully written out
	outfile->Write("*",TObject::kOverwrite);
	
	// add the final plots
	CreatePlots();
	
	// ensure all the histograms are fully written out
	outfile->Write("*",TObject::kOverwrite);
	
	// reset ttree addresses
	outtree->ResetBranchAddresses();
	outtree2->ResetBranchAddresses();
	
	// close the file
	outfile->Close();
	
	return true;
}

int DataSummary::ScanForFiles(std::string inputdir, std::string filepattern, std::string filepattern_orphan){
	// Scan the input directory for all files matching the specified pattern and run range
	
	// Step 1:
	// first, find matching files
	// TODO we should give booststores an extension: for now insist R*S*p* comes at the end of the filename
	std::string slash = "/";
	std::string lscommand = std::string("find ") + inputdir +std::string(" -regextype egrep -iregex '.*?")
	   + InputFilePattern + std::string("R([0-9]+)S([0-9]+)[pP]([0-9]+)$'");
	std::cout <<"lscommand: "<<lscommand<<std::endl;
	std::string fileliststring = GetStdoutFromCommand(lscommand);

	std::stringstream ssl;
	ssl << fileliststring;
	std::vector<std::string> flist;
	std::string nextfilestring;
	while(getline(ssl,nextfilestring)){
		flist.push_back(nextfilestring);
	}
	
	std::vector<std::string> flist_orphan;
	if (FileFormat == "SeparateStores"){
		std::string lscommand_orphan = std::string("find ")+inputdir+std::string(" -regextype egrep -iregex '.*?")
	   	+ InputFilePatternOrphan + std::string("R([0-9]+)S([0-9]+)[pP]([0-9]+)$'");
		std::cout <<"lscommand_orphan: "<<lscommand_orphan<<std::endl;
		std::string fileliststring_orphan = GetStdoutFromCommand(lscommand_orphan);
		std::cout <<"fileliststring_orphan: "<<fileliststring_orphan<<std::endl;
	
        	std::stringstream ssl_orphan;
		ssl_orphan << fileliststring_orphan;
		std::string nextfilestring_orphan;
		while(getline(ssl_orphan,nextfilestring_orphan)){
			flist_orphan.push_back(nextfilestring_orphan);
		}
		if (flist.size() != flist_orphan.size()) {
			Log("DataSummary tool error: Did not find the same number of orphan files and ANNIEEvent files. # of ANNIEEvent files: "+std::to_string(flist_orphan.size())+", # of Orphan files: "+std::to_string(flist.size()),v_error,verbosity);
			return 0;
		}
	}


	// now parse the list of files and extract those that the range of run/subrun numbers
	std::cout <<"Filling filelist"<<std::endl;
	filelist.clear();
	int i=0;
	for(auto afile : flist){
		/*//try{
			// extract run, subrun and part using a regex match
			std::cout <<"afile: "<<afile<<std::endl;
			std::smatch submatches;
			std::cout <<"regex theexpression"<<std::endl;
			std::regex theexpression(".*?R([0-9]+)S([0-9]+)[pP]([0-9]+).*",std::regex::extended);
			std::cout <<"regex_match"<<std::endl;
			std::regex_match ((std::string)afile, submatches, theexpression);
			std::cout <<"submatches.size(): "<<submatches.size()<<std::endl;
			std::string runstring = (std::string)submatches[1];
			std::string subrunstring = (std::string)submatches[2];
			std::string partstring = (std::string)submatches[3];
			run = stoi(runstring);
			subrun = stoi(subrunstring);
			part = stoi(partstring);
		//} catch(const std::out_of_range& oor){continue;}
		
		// check within range
		if(    (run>=StartRun || StartRun<0)
			&& (subrun>=StartSubRun || StartSubRun<0)
			&& (part>=StartPart || StartPart<0)
			&& (run<EndRun || EndRun<0)
			&& (subrun<EndSubRun || EndSubRun<0)
			&& (part<EndPart || EndPart<0)
		){
			// ensure files are in the correct order by using a map with suitable key
			// assume no more than 1000 parts per subrun, 1000 subruns per run, 1000000 runs
			char buffer [13];
			snprintf(buffer, 13, "%06d%03d%03d", run, subrun, part);
			filelist.emplace(buffer,afile);
		}*/

		//Get run number
		run = -1;
		subrun = -1;
		part = -1;
		std::stringstream InputFilePatternR;
		InputFilePatternR << InputFilePattern << "R";
		size_t rawdata_pos = afile.find(InputFilePatternR.str().c_str());
		if (rawdata_pos != std::string::npos){
			std::string filenamerun = afile.substr(rawdata_pos+InputFilePattern.length()+1);
			size_t pos_sub = filenamerun.find("S");
			std::string run_str = filenamerun.substr(0,pos_sub);
  			run = std::stoi(run_str);
			std::string filenamesubrun = filenamerun.substr(pos_sub+1);
			size_t pos_part = filenamesubrun.find("p");
			std::string subrun_str = filenamesubrun.substr(0,pos_part);
  			subrun = std::stoi(subrun_str);
			std::string filenamepart = filenamesubrun.substr(pos_part+1);
			part = std::stoi(filenamepart);
		}
		std::cout <<"Found run number: "<<run<<", subrun nr: "<<subrun<<", part nr: "<<part<<std::endl;
		// check within range
		if(    (run>=StartRun || StartRun<0)
			&& (subrun>=StartSubRun || StartSubRun<0)
			&& (part>=StartPart || StartPart<0)
			&& (run<=EndRun || EndRun<0)
			&& (subrun<=EndSubRun || EndSubRun<0)
			&& (part<=EndPart || EndPart<0)
		){
			// ensure files are in the correct order by using a map with suitable key
			// assume no more than 1000 parts per subrun, 1000 subruns per run, 1000000 runs
			char buffer [13];
			snprintf(buffer, 13, "%06d%03d%03d", run, subrun, part);
			filelist.emplace(buffer,afile);
		}
/*
		run=2282;
		subrun=0;
		char buffer [13];
		snprintf(buffer,13,"%06d%03d%03d", 2282, 0, i);
		filelist.emplace(buffer,afile);
		i++;*/
	}
	
	std::cout <<"Filling filelist_orphan"<<std::endl;
	if (FileFormat == "SeparateStores"){
		i=0;
		filelist_orphan.clear();
		for(auto afile : flist_orphan){
		/*	//try{
				// extract run, subrun and part using a regex match
				std::smatch submatches;
				std::regex theexpression (".*?R([0-9]+)S([0-9]+)[pP]([0-9]+).*");
				std::regex_match ((std::string)afile, submatches, theexpression);
				std::string runstring = (std::string)submatches[1];
				std::string subrunstring = (std::string)submatches[2];
				std::string partstring = (std::string)submatches[3];
				run = stoi(runstring);
				subrun = stoi(subrunstring);
				part = stoi(partstring);
			//} catch(const std::out_of_range& oor){continue;}
			
			// check within range
			if(    (run>=StartRun || StartRun<0)
				&& (subrun>=StartSubRun || StartSubRun<0)
				&& (part>=StartPart || StartPart<0)
				&& (run<=EndRun || EndRun<0)
				&& (subrun<=EndSubRun || EndSubRun<0)
				&& (part<=EndPart || EndPart<0)
			){
				// ensure files are in the correct order by using a map with suitable key
				// assume no more than 1000 parts per subrun, 1000 subruns per run, 1000000 runs
				char buffer [13];
				snprintf(buffer, 13, "%06d%03d%03d", run, subrun, part);
				filelist_orphan.emplace(buffer,afile);
			}*/
			
			//Get run number
                	run = -1;
                	subrun = -1;
                	part = -1;
                	std::stringstream InputFilePatternOrphanR;
                	InputFilePatternOrphanR << InputFilePatternOrphan << "R";
                	size_t rawdata_pos = afile.find(InputFilePatternOrphanR.str().c_str());
                	if (rawdata_pos != std::string::npos){
                        	std::string filenamerun = afile.substr(rawdata_pos+InputFilePatternOrphan.length()+1);
                        	size_t pos_sub = filenamerun.find("S");
                        	std::string run_str = filenamerun.substr(0,pos_sub);
                        	run = std::stoi(run_str);
                        	std::string filenamesubrun = filenamerun.substr(pos_sub+1);
                        	size_t pos_part = filenamesubrun.find("p");
                        	std::string subrun_str = filenamesubrun.substr(0,pos_part);
                        	subrun = std::stoi(subrun_str);
                        	std::string filenamepart = filenamesubrun.substr(pos_part+1);
                        	part = std::stoi(filenamepart);
                	}
                	std::cout <<"Found run number: "<<run<<", subrun nr: "<<subrun<<", part nr: "<<part<<std::endl;
               		// check within range
                	if(    (run>=StartRun || StartRun<0)
                        	&& (subrun>=StartSubRun || StartSubRun<0)
                        	&& (part>=StartPart || StartPart<0)
                        	&& (run<=EndRun || EndRun<0)
                       		&& (subrun<=EndSubRun || EndSubRun<0)
                        	&& (part<=EndPart || EndPart<0)
                	){
                        	// ensure files are in the correct order by using a map with suitable key
                        	// assume no more than 1000 parts per subrun, 1000 subruns per run, 1000000 runs
                        	char buffer [13];
                        	snprintf(buffer, 13, "%06d%03d%03d", run, subrun, part);
                        	filelist_orphan.emplace(buffer,afile);
                	}
			/*
			char buffer[13];
			snprintf(buffer,13,"%06d%03d%03d",2282,0,i);
			filelist_orphan.emplace(buffer,afile);
			i++;*/
		}
	}

	nextfile=filelist.begin();
	if (FileFormat == "SeparateStores") nextfile_orphan=filelist_orphan.begin();
	return filelist.size();
}

int DataSummary::ReadInFileList(std::string filelist_user){

	ifstream stream_filelist(filelist_user.c_str());
	int i=0;
	std::string file1, file2;
	while (!stream_filelist.eof()){
		if (FileFormat != "SeparateStores"){
			stream_filelist >> file1;
			if (stream_filelist.eof()) break;
			char buffer[13];
			snprintf(buffer,13,"%06d%03d%03d",0,0,i);
			filelist.emplace(buffer,file1);
			i++;
		} else {
			stream_filelist >> file1 >> file2;
			if (stream_filelist.eof()) break;
			char buffer[13];
			snprintf(buffer,13,"%06d%03d%03d",0,0,i);
			filelist.emplace(buffer,file1);
			filelist_orphan.emplace(buffer,file2);
			i++;	
		}
	}

	nextfile=filelist.begin();
	if (FileFormat == "SeparateStores") nextfile_orphan=filelist_orphan.begin();

	return i;

}

bool DataSummary::LoadNextANNIEEventEntry(){
	// load next ANNIEEvent entry
	localentry++;
	globalentry++;
	std::cout <<"Load localentry: "<<localentry<<" / "<<localentries<<" for ANNIEEvent"<<std::endl;
	if((localentry>=localentries)&&(localorphan>=localorphans)){
		// end of this file, try to load next file
		bool load_ok = LoadNextFile();
		if(not load_ok){
			// failed to load any further files. Terminate.
			m_data->vars.Set("StopLoop",1);
			return false;
		}
		return ANNIEEvent->GetEntry(localentry);
	} else if(localentry<localentries){
		return ANNIEEvent->GetEntry(localentry);
	} else {
          // else no more ANNIEEvents, but we didn't load a new file
	  // because there are still Orphans to process
	  return false;
	}
}

bool DataSummary::LoadNextOrphanStoreEntry(){
	// load next OrphanStore entry
	localorphan++;
	globalorphan++;
	std::cout <<"Load local orphan: "<<localorphan<<" / "<<localorphans<<" for OrphanStore"<<std::endl;
	// if we've run out of entries in both stores, load the next file
	if((localorphan>=localorphans)&&(localentry>=localentries)){
		bool load_ok = LoadNextFile();
		std::cout <<"load_ok: "<<load_ok<<std::endl;
		if(not load_ok){
			// failed to load any further files. Terminate.
			m_variables.Set("StopLoop",true);
			return false;
		}
		return OrphanStore->GetEntry(localorphan);
	} else if(localorphan<localorphans){
		// not a new file, still entries to process
		return OrphanStore->GetEntry(localorphan);
	} else {
	 // else no more orphans, but we didn't load a new file
	  // because there are still ANNIEEvents to process
	return false;
	}
}

bool DataSummary::LoadNextFile(){
	// get next filename
	if(nextfile==filelist.end()){
		// no more files
		return false;
	}
	std::string nextfilename = nextfile->second;
	
	std::string nextfilename_orphan;
	if (FileFormat=="SeparateStores"){
		if (nextfile_orphan==filelist_orphan.end()){
			return false;
		}
		nextfilename_orphan = nextfile_orphan->second;
	}

	// Delete the old file BoostStores if they exist
	if(ProcessedFileStore){
		ProcessedFileStore->Close();
		ProcessedFileStore->Delete();
		delete ProcessedFileStore;
	}
	if(ANNIEEvent){
		ANNIEEvent->Close();
		ANNIEEvent->Delete();
		delete ANNIEEvent;
	}
	if(OrphanStore){
		OrphanStore->Close();
		OrphanStore->Delete();
		delete OrphanStore;
	}
	
	// create a store for the file contents
        if (FileFormat == "CombinedStore") {
		ProcessedFileStore = new BoostStore(false,BOOST_STORE_BINARY_FORMAT);
		// Load the contents from the new input file into it
		std::cout <<"Reading in current file "<<nextfilename<<std::endl;
		int load_ok = ProcessedFileStore->Initialise(nextfilename);
		if(not load_ok) return false;
	}

	// create an ANNIEEvent BoostStore and an OrphanStore BoostStore to load from it
	ANNIEEvent = new BoostStore(false, BOOST_STORE_MULTIEVENT_FORMAT);
	OrphanStore = new BoostStore(false, BOOST_STORE_MULTIEVENT_FORMAT);
	
	if (FileFormat == "CombinedStore"){

	// retrieve the multi-event stores
	ProcessedFileStore->Get("ANNIEEvent",*ANNIEEvent);
	ANNIEEvent->Header->Get("TotalEntries", localentries);
	
	// same for Orphan Store
	ProcessedFileStore->Get("OrphanStore",*OrphanStore);
	OrphanStore->Header->Get("TotalEntries", localorphans);
	}
	else {

		Log("DataSummary tool: Reading in current ANNIEEvent file "+nextfilename,v_message,verbosity);
		ANNIEEvent->Initialise(nextfilename);
		ANNIEEvent->Header->Get("TotalEntries",localentries);

		Log("DataSummary tool: Reading in current Orphan file "+nextfilename_orphan,v_message,verbosity);
		OrphanStore->Initialise(nextfilename_orphan);
		OrphanStore->Header->Get("TotalEntries",localorphans);

	}

	nextfile++;
	localentry = 0;
	localorphan = 0;
	
	if (FileFormat == "SeparateStores") nextfile_orphan++;

	// load run constants and sanity checks
	get_ok = ANNIEEvent->GetEntry(0);
	if(not get_ok){
		Log("DataSummary Tool: Error getting ANNIEEvent entry 0 from new file "+nextfilename,v_error,verbosity);
		return false;
	}
	ANNIEEvent->Get("RunNumber",RunNumber);
	std::cout <<"RunNumber: "<<RunNumber<<", run: "<<run<<std::endl;
	if(RunNumber!=run)
		Log("DataSummary Tool: filename / entry mismatch for RunNumber!",v_error,verbosity);
	ANNIEEvent->Get("SubrunNumber",SubrunNumber);
	if(SubrunNumber!=subrun)
		Log("DataSummary Tool: filename / entry mismatch for SubrunNumber!",v_error,verbosity);
	// TODO... handle these errors? We should have an error log file.
	PartNumber=part; // not stored so assume it's the same
	
	// we can read these just once at the start of the file
	ANNIEEvent->Get("RunStartTime",RunStartTime);
	ANNIEEvent->Get("RunType",RunType);
	RunTypeString = runtype_to_string(RunTypeEnum(RunType));
	
	// TODO sanity checks to ensure these don't change between part files
	
	get_ok = OrphanStore->GetEntry(0);
	if(not get_ok){
		// I guess there aren't necessarily any orphans...
		Log("DataSummary Tool: Warning! No orphans in OrphanStore for file "+nextfilename,v_warning,verbosity);
		//Log("DataSummary Tool: Error getting OrphanStore entry 0 from new file "+nextfilename,v_error,verbosity);
		//return false;
	}
	
	return true;
}

bool DataSummary::CreateOutputFile(){
	outfile = new TFile((OutputFileDir+"/"+OutputFileName).c_str(),"RECREATE");
	outtree = new TTree("EventStats","EventStats Tree");
	
	// set the output branches
	outtree->Branch("RunNumber",&RunNumber);
	outtree->Branch("SubrunNumber",&SubrunNumber);
	outtree->Branch("PartNumber",&PartNumber);
	outtree->Branch("RunType",&RunType);
	outtree->Branch("RunTypeString",&RunTypeString);
//	outtree->Branch("EventNumber",&EventNumber); // processed event num? file event num? run event num?
	outtree->Branch("TriggerWord",&TriggerWord);
	outtree->Branch("TriggerTypeString",&TriggerTypeString);
	outtree->Branch("CTCTimestamp",&CTCtimestamp_tree);
	outtree->Branch("PMTTimestamp",&PMTtimestamp_tree);
	outtree->Branch("MRDTimestamp",&MRDtimestamp_tree);
	outtree->Branch("CTCTimestamp_double",&CTCtimestamp_double);
	outtree->Branch("PMTTimestamp_double",&PMTtimestamp_double);
	outtree->Branch("MRDTimestamp_double",&MRDtimestamp_double);
	outtree->Branch("CTCTimestamp_sec",&CTCtimestamp_sec);
	outtree->Branch("PMTTimestamp_sec",&PMTtimestamp_sec);
	outtree->Branch("MRDTimestamp_sec",&MRDtimestamp_sec);
	outtree->Branch("BeamLoopbackTimestamp",&BeamLoopbackTimestamp_tree);
	outtree->Branch("CosmicLoopbackTimestamp",&CosmicLoopbackTimestamp_tree);
	outtree->Branch("SystemsPresent",&SystemsPresent);
	outtree->Branch("LoopbacksPresent",&LoopbacksPresent);
	outtree->Branch("DataCTC",&data_ctc);
	outtree->Branch("DataTank",&data_tank);
	outtree->Branch("DataMRD",&data_mrd);
	outtree->Branch("ExtendedWindow",&window_is_extended);
	outtree->Branch("WindowSize",&size_of_window);
	outtree->Branch("AllWindowSizes",&window_sizes);
	outtree->Branch("AllWindowChankeys",&window_chkeys);
	outtree->Branch("TrigExtended",&trigword_ext);
	outtree->Branch("TrigExtendedCC",&trigword_ext_cc);
	outtree->Branch("TrigExtendedNC",&trigword_ext_nc);
	outtree->SetAlias("CtcToTankTDiff","CTCTimestamp_double-PMTTimestamp_double");
	outtree->SetAlias("CtcToMrdTDiff","CTCTimestamp_double-MRDTimestamp_double");
	outtree->SetAlias("TankToMrdTDiff","PMTTimestamp_double-MRDTimestamp_double");
	outtree->SetAlias("CtcToBeamLoopbackTDiff","CTCTimestamp_double-BeamLoopbackTimestamp");
	outtree->SetAlias("CtcToCosmicLoopbackTDiff","CTCTimestamp_double-CosmicLoopbackTimestamp");
	
	outtree2 = new TTree("OrpahStats","OrpahStats Tree");
	outtree2->Branch("OrphanedEventType",&orphantype);
	outtree2->Branch("OrphanTimestamp",&orphantimestamp_tree);
	outtree2->Branch("OrphanTimestamp_double",&orphantimestamp_double);
	outtree2->Branch("OrphanCause",&orphancause);
	outtree2->Branch("OrphanNumWaves",&orphannumwaves);
	outtree2->Branch("OrphanChankeys",&orphanchankeys_int);
	outtree2->Branch("OrphanChannels",&orphanchannels_combined);
	outtree2->Branch("OrphanMinTDiff",&orphanmintdiff);
	outtree2->Branch("OrphanTrigWord",&orphantrigword);


	return true;
}

bool DataSummary::CreatePlots(){

	// make plots from the ROOT tree
	// =============================
	// first get the timespan we're going to be plotting
	
	double t0_matched, t0_orphaned;
	double tn_matched, tn_orphaned;

	outtree->GetEntry(0);
	t0_matched = CTCtimestamp_sec;
	outtree->GetEntry(outtree->GetEntries()-1);
	tn_matched = CTCtimestamp_sec;

	outtree2->GetEntry(0);
	t0_orphaned = orphantimestamp_double/1.E9;
	outtree2->GetEntry(outtree2->GetEntries()-1);
	tn_orphaned = orphantimestamp_double/1.E9;
	t0 = (t0_matched < t0_orphaned)? (t0_matched-60) : (t0_orphaned - 60);
	tn = (tn_matched > tn_orphaned)? (tn_matched+60) : (tn_orphaned + 60);
	std::cout <<"CreatePlots: t0: "<<t0<<", tn: "<<tn<<std::endl;
	// we need to set a global "time offset", from which all other time axes will be relative to
	// this should be in seconds
	std::cout <<"Set time offset"<<std::endl;
	//gStyle->SetTimeOffset(t0/1E9);
	gStyle->SetTimeOffset(0.);	

	// Rate plots: should be able to make these just by time-x-axis histograms
	// TODO: Refactor?
	std::cout <<"AddRatePlots"<<std::endl;
	AddRatePlots(500);

	// EventType plots: Rates, fractions & other properties of certain orphan event types
	AddEventTypePlots();
	
	// TODO: Refactor this
	// Time evolution of the timestamp discrepancies
	std::cout <<"Time diff plots"<<std::endl;
	AddTDiffPlots();
	
	return true;

}

bool DataSummary::AddRatePlots(int nbins){

	// -------------------------------------------------------------
	// just histograms of timestamps of a given type with a time axis
	// -------------------------------------------------------------	

	outfile->cd(); // ensure plots get put in the file
	
	// define our histogram min/max relative to global time offset, in seconds
	double t1 = t0;
	double t2 = tn;

//	anything event rate
//	all-systems event rate
//	event with ctc rate
//	event with pmt rate
//	event with mrd rate
//	event with cosmic mrd loopback rate
//	event with beam mrd loopback rate
//	events without any loopback rate
//	orphaned events rate
//	fraction of events orphaned (ratio of bin counts from num orphans rate to (num orphans + all events) rate)
	
	// make the histograms
	TH1D* ctc_rate = new TH1D("ctc_rate","CTC Rate",nbins,t1,t2);
	outtree->Draw("CTCTimestamp_double/(1.E9)>>ctc_rate","CTCTimestamp>0","goff");
	TH1D* tank_rate = new TH1D("tank_rate","Tank Rate",nbins,t1,t2);
	outtree->Draw("PMTTimestamp_double/(1.E9)>>tank_rate","PMTTimestamp>0","goff");
	TH1D* mrd_rate = new TH1D("mrd_rate","MRD Rate",nbins,t1,t2);
	outtree->Draw("MRDTimestamp_double/(1.E9)>>mrd_rate","MRDTimestamp>0","goff");
	TH1D* allsystems_rate = new TH1D("allsystems_rate","All Systems Triggered Rate",nbins,t1,t2);
	outtree->Draw("CTCTimestamp_double/(1.E9)>>allsystems_rate","MRDTimestamp>0&&CTCTimestamp>0&&PMTTimestamp>0","goff");
	TH1D* noloopback_rate = new TH1D("noloopback_rate","No Loopback Triggered Rate",nbins,t1,t2);
	outtree->Draw("CTCTimestamp_double/(1.E9)>>noloopback_rate","LoopbacksPresent==0","goff");
	TH1D* orphan_rate = new TH1D("orphan_rate","Orphan Rate",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate","","goff");
	TH1D* orphan_rate_incomplete = new TH1D("orphan_rate_incomplete","Orphan Rate (incomplete VME)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_incomplete","OrphanCause==\"incomplete_tank_event\"","goff");	
	TH1D* orphan_rate_tank_noctc = new TH1D("orphan_rate_tank_noctc","Orphan Rate (Tank, No CTC)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_tank_noctc","OrphanCause==\"tank_no_ctc\"","goff");	
	TH1D* orphan_rate_pmt = new TH1D("orphan_rate_pmt","Orphan Rate (PMT)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_pmt","OrphanedEventType==\"Tank\"","goff");	
	TH1D* orphan_rate_mrd = new TH1D("orphan_rate_mrd","Orphan Rate (MRD)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_mrd","OrphanedEventType==\"MRD\"","goff");	
	TH1D* orphan_rate_ctc = new TH1D("orphan_rate_ctc","Orphan Rate (CTC)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_ctc","OrphanedEventType==\"CTC\"","goff");	
	TH1D* orphan_rate_ctc_trig5 = new TH1D("orphan_rate_ctc_trig5","Orphan Rate (CTC, trigword 5)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_ctc_trig5","OrphanedEventType==\"CTC\" && OrphanTrigWord == 5","goff");	
	TH1D* orphan_rate_ctc_trig31 = new TH1D("orphan_rate_ctc_trig31","Orphan Rate (CTC, trigword 31)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_ctc_trig31","OrphanedEventType==\"CTC\" && OrphanTrigWord == 31","goff");	
	TH1D* orphan_rate_ctc_trig33 = new TH1D("orphan_rate_ctc_trig33","Orphan Rate (CTC, trigword 33)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_ctc_trig33","OrphanedEventType==\"CTC\" && OrphanTrigWord == 33","goff");	
	TH1D* orphan_rate_ctc_trig36 = new TH1D("orphan_rate_ctc_trig36","Orphan Rate (CTC, trigword 36)",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp_double/(1.E9)>>orphan_rate_ctc_trig36","OrphanedEventType==\"CTC\" && OrphanTrigWord == 36","goff");	

	// then we need to set the axes to time
	ctc_rate->GetXaxis()->SetTimeDisplay(1);
	tank_rate->GetXaxis()->SetTimeDisplay(1);
	mrd_rate->GetXaxis()->SetTimeDisplay(1);
	allsystems_rate->GetXaxis()->SetTimeDisplay(1);
	noloopback_rate->GetXaxis()->SetTimeDisplay(1);
	orphan_rate->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_incomplete->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_tank_noctc->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_pmt->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_mrd->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_ctc->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_ctc_trig5->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_ctc_trig31->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_ctc_trig33->GetXaxis()->SetTimeDisplay(1);
	orphan_rate_ctc_trig36->GetXaxis()->SetTimeDisplay(1);

	// set the timestamp format
	ctc_rate->GetXaxis()->SetLabelSize(0.03);
	ctc_rate->GetXaxis()->SetLabelOffset(0.03);
	tank_rate->GetXaxis()->SetLabelSize(0.03);
	tank_rate->GetXaxis()->SetLabelOffset(0.03);
	mrd_rate->GetXaxis()->SetLabelSize(0.03);
	mrd_rate->GetXaxis()->SetLabelOffset(0.03);
	allsystems_rate->GetXaxis()->SetLabelSize(0.03);
	allsystems_rate->GetXaxis()->SetLabelOffset(0.03);
	noloopback_rate->GetXaxis()->SetLabelSize(0.03);
	noloopback_rate->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate->GetXaxis()->SetLabelSize(0.03);
	orphan_rate->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_incomplete->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_incomplete->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_tank_noctc->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_tank_noctc->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_pmt->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_pmt->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_mrd->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_mrd->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_ctc->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_ctc->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_ctc_trig5->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_ctc_trig5->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_ctc_trig31->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_ctc_trig31->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_ctc_trig33->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_ctc_trig33->GetXaxis()->SetLabelOffset(0.03);
	orphan_rate_ctc_trig36->GetXaxis()->SetLabelSize(0.03);
	orphan_rate_ctc_trig36->GetXaxis()->SetLabelOffset(0.03);
	ctc_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	tank_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	mrd_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	allsystems_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	noloopback_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_incomplete->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_tank_noctc->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_pmt->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_mrd->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_ctc->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_ctc_trig5->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_ctc_trig31->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_ctc_trig33->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	orphan_rate_ctc_trig36->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	
	return true;

}

bool DataSummary::AddEventTypePlots(){


	outfile->cd(); // ensure plots get put in the file

	TH1D* num_waveforms_orphan = new TH1D("num_waveforms_orphan","Number of waveforms (orphaned events)",200,0,200);	
	TH1D* orphan_types = new TH1D("orphan_types","Orphan Event Types",4,0,4);
	TH1D* orphan_types_fractions = new TH1D("orphan_types_fractions","Orphan Event Types (Fractions)",4,0,4);
	TH1D* orphan_types_rates = new TH1D("orphan_types_rates","Orphan Event Types (Rates)",4,0,4);
	TH1D* orphan_reasons = new TH1D("orphan_reasons","Orphan Reasons",4,0,4);
	TH1D* orphan_reasons_fractions = new TH1D("orphan_reasons_fractions","Orphan Reasons (Fractions)",4,0,4);
	TH1D* orphan_reasons_rates = new TH1D("orphan_reasons_rates","Orphan Reasons (Rates)",4,0,4);
	TH1D* waveform_chankeys_orphan = new TH1D("waveform_chankeys_orphan","Waveform channelkeys (orphaned events)",500,0,500);
	TH1D* orphan_mintdiff_tank = new TH1D("orphan_mintdiff_tank","Minimum time difference CTC (orphaned tank events)",500,-1000,1000);
	TH1D* orphan_mintdiff_mrd = new TH1D("orphan_mintdiff_mrd","Minimum time difference CTC (orphaned MRD events)",500,-100000000,100000000);
	TH1D* waveform_channels_orphan = new TH1D("waveform_channels_orphan","Waveform channels (orphaned events)",5000,0,5000);
	TH1D* ctc_orphans_triggerword = new TH1D("ctc_orphans_triggerword","CTC Orphans - Triggerword",64,0,64);
	TH1D* datastreams_present = new TH1D("datastreams_present","Present datastreams",3,0,3);
	TH1D* datastreams_present_fractions = new TH1D("datastreams_present_fractions","Present datastreams (fractions)",3,0,3);
	TH1D* extended_types = new TH1D("extended_types","Extended Triggerwords",4,0,4);
	TH1D* extended_types_fractions = new TH1D("extended_types_fractions","Extended Triggerwords (Fractions)",4,0,4);
	TH1D* extended_types_rates = new TH1D("extended_types_rates","Extended Triggerwords (Rates)",4,0,4);
	TH1D* extended_interplay = new TH1D("extended_interplay","Extended/Prompt Triggerwords",6,0,6);
	TH1D* extended_interplay_fractions = new TH1D("extended_interplay_fractions","Extended/Prompt Triggerwords (Fractions)",6,0,6);
	TH1D* extended_interplay_rates = new TH1D("extended_interplay_rates","Extended/Prompt Triggerwords (Rates)",6,0,6);


	int numwaves_temp;
	std::string *orphantype = new std::string;
	std::string *orphancause = new std::string;
	std::vector<int> *orphanchkeys = new std::vector<int>;
	std::vector<int> *orphanch = new std::vector<int>;
	double mintdiff;
	int trigword;
	outtree2->SetBranchAddress("OrphanNumWaves",&numwaves_temp);
	outtree2->SetBranchAddress("OrphanedEventType",&orphantype);
	outtree2->SetBranchAddress("OrphanCause",&orphancause);
	outtree2->SetBranchAddress("OrphanChankeys",&orphanchkeys);
	outtree2->SetBranchAddress("OrphanChannels",&orphanch);
	outtree2->SetBranchAddress("OrphanMinTDiff",&mintdiff);
	outtree2->SetBranchAddress("OrphanTrigWord",&trigword);
	int entries_tree1 = outtree->GetEntries();
	int entries_tree2 = outtree2->GetEntries();
	for (int i_entry=0; i_entry < entries_tree2; i_entry++){
		outtree2->GetEntry(i_entry);
		orphan_types->Fill(0);
		if (*orphantype == "Tank"){
			orphan_types->Fill(1);
			orphan_mintdiff_tank->Fill(mintdiff);
			std::cout <<"orphanmintdiff: "<<mintdiff<<std::endl;
			if (*orphancause == "incomplete_tank_event"){
				num_waveforms_orphan->Fill(numwaves_temp);
				orphan_reasons->Fill(1);
			} else if (*orphancause == "tank_no_ctc"){
				orphan_reasons->Fill(0);
			}
			for (int i_vec=0; i_vec < (int) orphanchkeys->size(); i_vec++){
				waveform_chankeys_orphan->Fill(orphanchkeys->at(i_vec));
			}
			for (int i_vec=0; i_vec < (int) orphanch->size(); i_vec++){
				waveform_channels_orphan->Fill(orphanch->at(i_vec));
			}
		} else if (*orphantype == "MRD"){
			orphan_types->Fill(2);
			orphan_mintdiff_mrd->Fill(mintdiff);
			std::cout <<"mintdiff (MRD): "<<mintdiff<<std::endl;
			if (*orphancause == "mrd_beam_no_ctc"){
				orphan_reasons->Fill(2);
			}
		} else if (*orphantype == "CTC"){
			orphan_types->Fill(3);
			ctc_orphans_triggerword->Fill(trigword);
			if (*orphancause == "ctc_no_mrd_or_tank"){
				orphan_reasons->Fill(3);
			}
		}
	}

	delete orphantype;
	delete orphancause;

	bool has_ctc, has_tank, has_mrd;
	bool has_extended, has_extended_cc, has_extended_nc, has_extended_vme;
	outtree->SetBranchAddress("DataCTC",&has_ctc);
	outtree->SetBranchAddress("DataTank",&has_tank);
	outtree->SetBranchAddress("DataMRD",&has_mrd);
	outtree->SetBranchAddress("TrigExtended",&has_extended);
	outtree->SetBranchAddress("TrigExtendedCC",&has_extended_cc);
	outtree->SetBranchAddress("TrigExtendedNC",&has_extended_nc);
	outtree->SetBranchAddress("ExtendedWindow",&has_extended_vme);

	for (int i_entry=0; i_entry < entries_tree1; i_entry++){
		outtree->GetEntry(i_entry);
		if (has_ctc && has_tank & has_mrd) datastreams_present->Fill(0);
		else if (has_ctc && has_tank) datastreams_present->Fill(1);
		else if (has_ctc && has_mrd) datastreams_present->Fill(2);
		if (has_extended) {
			extended_types->Fill(0);
			if (has_extended_vme) extended_interplay->Fill(0);
			else extended_interplay->Fill(1);
		}
		if (has_extended_cc && has_extended_nc) extended_types->Fill(3);
		else if (has_extended_cc) {
			extended_types->Fill(1);
			if (has_extended_vme) extended_interplay->Fill(2);
			else extended_interplay->Fill(3);
		}
		else if (has_extended_nc) {
			extended_types->Fill(2);
			if (has_extended_vme) extended_interplay->Fill(4);
			else extended_interplay->Fill(5);
		}
	}

	orphan_types->GetYaxis()->SetTitle("#");
	orphan_types_fractions->GetYaxis()->SetTitle("Fraction");
	const char *label_eventtypes[4] = {"All","Tank","MRD","CTC"};
	orphan_reasons->GetYaxis()->SetTitle("#");
	orphan_reasons_fractions->GetYaxis()->SetTitle("Fraction");
	const char *label_reasons[4] = {"Tank - No CTC","Tank - Incomplete Waveform","MRD - No CTC","CTC - No MRD/Tank"};
	for (int i_bin=0; i_bin<orphan_types->GetXaxis()->GetNbins(); i_bin++){
		orphan_types->GetXaxis()->SetBinLabel(i_bin+1,label_eventtypes[i_bin]);
		orphan_types_fractions->SetBinContent(i_bin+1,orphan_types->GetBinContent(i_bin+1)/(double(entries_tree1)+double(entries_tree2)));
		orphan_types_fractions->GetXaxis()->SetBinLabel(i_bin+1,label_eventtypes[i_bin]);
		orphan_reasons->GetXaxis()->SetBinLabel(i_bin+1,label_reasons[i_bin]);
		orphan_reasons_fractions->SetBinContent(i_bin+1,orphan_reasons->GetBinContent(i_bin+1)/(double(entries_tree1)+double(entries_tree2)));
		orphan_reasons_fractions->GetXaxis()->SetBinLabel(i_bin+1,label_reasons[i_bin]);
	}
	datastreams_present->GetYaxis()->SetTitle("#");
	datastreams_present_fractions->GetYaxis()->SetTitle("Fraction");
	const char *label_datastreams[3]={"CTC+Tank+MRD","CTC+Tank","CTC+MRD"};
	for (int i_bin=0; i_bin < datastreams_present->GetXaxis()->GetNbins(); i_bin++){
		datastreams_present->GetXaxis()->SetBinLabel(i_bin+1,label_datastreams[i_bin]);
		datastreams_present_fractions->SetBinContent(i_bin+1,datastreams_present->GetBinContent(i_bin+1)/(double(entries_tree1)));
		datastreams_present_fractions->GetXaxis()->SetBinLabel(i_bin+1,label_datastreams[i_bin]);
	}
	extended_types->GetYaxis()->SetTitle("#");
	extended_types_fractions->GetYaxis()->SetTitle("Fraction");
	extended_types_rates->GetYaxis()->SetTitle("Rate");
	const char *label_extended[4] = {"Extended","Extended - CC","Extended - NC","Extended - CC+NC"};
	for (int i_bin=0; i_bin < extended_types->GetXaxis()->GetNbins(); i_bin++){
		extended_types->GetXaxis()->SetBinLabel(i_bin+1,label_extended[i_bin]);
		extended_types_fractions->GetXaxis()->SetBinLabel(i_bin+1,label_extended[i_bin]);
		extended_types_fractions->SetBinContent(i_bin+1,extended_types->GetBinContent(i_bin+1)/double(extended_types->GetEntries()));
	}
	extended_interplay->GetYaxis()->SetTitle("#");
	extended_interplay_fractions->GetYaxis()->SetTitle("Fraction");
	extended_interplay_rates->GetYaxis()->SetTitle("Rate");
	const char *label_interplay[6] = {"Extended + VME","Extended No VME","ExtendedCC +VME","ExtendedCC No VME","ExtendedNC +VME","ExtendedNC No VME"};
	for (int i_bin=0; i_bin < extended_interplay->GetXaxis()->GetNbins(); i_bin++){
		extended_interplay->GetXaxis()->SetBinLabel(i_bin+1,label_interplay[i_bin]);
		extended_interplay_fractions->GetXaxis()->SetBinLabel(i_bin+1,label_interplay[i_bin]);
		extended_interplay_fractions->SetBinContent(i_bin+1,extended_interplay->GetBinContent(i_bin+1)/double(extended_interplay->GetEntries()));
	}

	num_waveforms_orphan->GetYaxis()->SetTitle("#");
	num_waveforms_orphan->GetXaxis()->SetTitle("# waveforms");
	waveform_chankeys_orphan->GetXaxis()->SetTitle("Chankey");
	waveform_chankeys_orphan->GetYaxis()->SetTitle("#");
	waveform_channels_orphan->GetXaxis()->SetTitle("Electronic channel");
	waveform_channels_orphan->GetYaxis()->SetTitle("#");
	orphan_mintdiff_tank->GetXaxis()->SetTitle("Minimum #Delta t [ns]");
	orphan_mintdiff_tank->GetYaxis()->SetTitle("#");
	orphan_mintdiff_mrd->GetXaxis()->SetTitle("Minimum #Delta t [ns]");
	orphan_mintdiff_mrd->GetYaxis()->SetTitle("#");
	ctc_orphans_triggerword->GetXaxis()->SetTitle("Triggerword");
	ctc_orphans_triggerword->GetYaxis()->SetTitle("#");

	orphan_types->Write();
	orphan_types_fractions->Write();
	orphan_reasons->Write();
	orphan_reasons_fractions->Write();
	num_waveforms_orphan->Write();
	waveform_chankeys_orphan->Write();
	waveform_channels_orphan->Write();
	orphan_mintdiff_tank->Write();
	orphan_mintdiff_mrd->Write();
	ctc_orphans_triggerword->Write();
	datastreams_present->Write();
	datastreams_present_fractions->Write();
	extended_types->Write();
	extended_types_fractions->Write();
	extended_interplay->Write();
	extended_interplay_fractions->Write();

	return true;
}

// TODO refactor to break up plot types and remove triplets of calls
bool DataSummary::AddTDiffPlots(){
	// make the file the active directory so our histograms get put there
	outfile->cd();
	// TODO we'll do a few plots here:
	// 1. we can make a normal histogram showing the distribution of timestamp differences
	TH1D* all_tank_ctc_tdiffs = new TH1D("all_tank_ctc_tdiffs","All (CTC-Tank) TDiffs",100,-100,100);
	outtree->Draw("CtcToTankTDiff>>all_tank_ctc_tdiffs","","goff");
	// btw we can grab that plotted data for later
	double* datapointer = outtree->GetV1();
	std::vector<double> tank_ctc_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	TH1D* all_tank_ctc_tdiffs_wholerange = new TH1D("all_tank_ctc_tdiffs_wholerange","All (CTC-Tank) TDiffs",500,-1000,1000);
	outtree->Draw("CtcToTankTDiff>>all_tank_ctc_tdiffs_wholerange","","goff");
	TH1D* all_tank_mrd_tdiffs = new TH1D("all_tank_mrd_tdiffs","All (Tank-MRD) TDiffs",100,-1E6,1E6);	//MRD timestamps much cruder, can be up to 1ms away!
	datapointer = outtree->GetV1();
	std::vector<double> mrd_tank_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	outtree->Draw("TankToMrdTDiff>>all_tank_mrd_tdiffs","","goff");
	TH1D* all_mrd_ctc_tdiffs = new TH1D("all_mrd_ctc_tdiffs","All (CTC-MRD) TDiffs",100,-1E6,1E6);
	outtree->Draw("CtcToMrdTDiff>>all_mrd_ctc_tdiffs","","goff");
	datapointer = outtree->GetV1();
	std::vector<double> mrd_ctc_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	TH1D* all_mrd_ctc_tdiffs_wholerange = new TH1D("all_mrd_ctc_tdiffs_wholerange","All (CTC-MRD) TDiffs",500,-100E6,100E6);
	outtree->Draw("CtcToMrdTDiff>>all_mrd_ctc_tdiffs_wholerange","","goff");
	
	// 2. we could do the same for the above, broken down by run
	
	// 3. we could also make 2D scatter versions, where we can see how the differences evolve over time
	// define our histogram min/max - these are relative to the global time offset, in seconds
	double t1 = t0;
	double t2 = tn;
	// TGraph* tank_ctc_diffs = new TGraph("tank_ctc_diffs","All (CTC-Tank) TDiffs",100,t1,t2);
	// apparently there really isn't a way to draw from a TTree straight into a named TGraph
	TH2F *tank_ctc_diffs = new TH2F("tank_ctc_diffs","(Tank - CTC) TDiffs",500,t1,t2,100,-100,100);
	tank_ctc_diffs->GetYaxis()->SetTitle("CTC - Tank [ns]");
	tank_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	tank_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	tank_ctc_diffs->GetXaxis()->SetLabelOffset(0.03);
	tank_ctc_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	outtree->Draw("CtcToTankTDiff:CTCTimestamp_sec>>tank_ctc_diffs","CTCTimestamp>0&&PMTTimestamp>0","goff");
	datapointer = outtree->GetV2();
	std::vector<double> ctc_tvals(datapointer, datapointer + outtree->GetSelectedRows());

	TH2F *mrd_ctc_diffs = new TH2F("mrd_ctc_diffs","(MRD - CTC) TDiffs",500,t1,t2,200,-1E6,1E6);
	mrd_ctc_diffs->GetYaxis()->SetTitle("CTC - MRD [ns]");
	mrd_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	mrd_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	mrd_ctc_diffs->GetXaxis()->SetLabelOffset(0.03);
	mrd_ctc_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	outtree->Draw("CtcToMrdTDiff:CTCTimestamp_sec>>mrd_ctc_diffs","CTCTimestamp>0&&MRDTimestamp>0 && fabs(CtcToMrdTDiff)<10E6","goffcolz");

	TH2F *tank_mrd_diffs = new TH2F("tank_mrd_diffs","(Tank - MRD) TDiffs",500,t1,t2,200,-1E6,1E6);
	tank_mrd_diffs->GetYaxis()->SetTitle("MRD - Tank [ns]");
	tank_mrd_diffs->GetXaxis()->SetTimeDisplay(1);
	tank_mrd_diffs->GetXaxis()->SetLabelSize(0.03);
	tank_mrd_diffs->GetXaxis()->SetLabelOffset(0.03);
	tank_mrd_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	outtree->Draw("TankToMrdTDiff:MRDTimestamp_sec>>tank_mrd_diffs","PMTTimestamp>0&&MRDTimestamp>0 && fabs(TankToMrdTDiff)<10E6","goffcolz");

	//Replaced the TGraphs by TH2 histograms (see above)
	/*TGraph* tank_ctc_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	tank_ctc_diffs->SetName("tank_ctc_diffs");
	tank_ctc_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;CTC - Tank [s]");
	
	outtree->Draw("CTCTimestamp:CtcToMrdTDiff","CTCTimestamp>0&&MRDTimestamp>0","goff");
	TGraph* mrd_ctc_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	mrd_ctc_diffs->SetName("mrd_ctc_diffs");
	mrd_ctc_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;CTC - MRD [s]");
	
	outtree->Draw("PMTTimestamp:TankToMrdTDiff","PMTTimestamp>0&&MRDTimestamp>0","goff");
	TGraph* tank_mrd_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	tank_mrd_diffs->SetName("tank_mrd_diffs");
	tank_mrd_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;Tank - MRD [s]");
	
	// then we need to set the axes to time
	tank_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	mrd_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	tank_mrd_diffs->GetXaxis()->SetTimeDisplay(1);
	// set the timestamp format
	tank_ctc_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	mrd_ctc_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	tank_mrd_diffs->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M:%S}");
	// make labels small enough that we don't overcrowd the axis
	tank_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	mrd_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	tank_mrd_diffs->GetXaxis()->SetLabelSize(0.03);*/
	
	// 4. we may also want to track just the mean and variance of timestamp differences over time,
	// in which case those need to be calculated over some sliding window. That window may be defined
	// either by the last N events, or the last N minutes, and may have some overlap between windows
	// or may use consecutive batches with no overlap.
	// An overlap fraction can define how the averaging window advances between measurements
	// i.e. if it's one, then each averaging window is independent but flush: [1,2,3] -> [4,5,6],...
	// the number of datapoints will be reduced by a factor of the window length
	// if it's 0.5, then the latter half of each window will be the former half of the next window:
	// [1,2,3,4] -> [3,4,5,6] -> [5,6,7,8]...
	// the data reduction factor will be overlap*window
	// this allows you to adjust the averaging (and noise) and data reduction independently.
	// for example
	// number of points to use
	int window_size=10;
	double overlap_fraction=0.5;
	
	// will be calculated, over last N events
	double mean_ctc_to_tank;          // mean time difference between tank and ctc times of matched events
	double mean_ctc_to_mrd;           // variance of time differences between tank and ctc times of matched events
	double mean_tank_to_mrd;          // same for tank to mrd
	double var_tank_to_mrd;           // same for tank to mrd
	double var_ctc_to_tank;           // same for mrd to ctc
	double var_ctc_to_mrd;            // same for mrd to ctc
	
	// calculate mean and variances over this window
	std::vector<double> tank_ctc_means;
	tank_ctc_means.reserve(tank_ctc_diff_vals.size()/(overlap_fraction*window_size));
	std::vector<double> tank_ctc_vars;
	tank_ctc_vars.reserve(tank_ctc_diff_vals.size()/(overlap_fraction*window_size));
	std::vector<double> tank_ctc_ts;
	tank_ctc_ts.reserve(tank_ctc_diff_vals.size()/(overlap_fraction*window_size));
	int start_sample=0;
	int step_size = window_size*overlap_fraction;
	std::cout <<"tank_ctc_diff_vals.size(): "<<tank_ctc_diff_vals.size()<<", ctc_tvals.size(): "<<ctc_tvals.size()<<std::endl;
	for(int i=0; i<(int)tank_ctc_diff_vals.size(); ++i){
		//std::cout <<"i: "<<i<<std::endl;
		ComputeMeanAndVariance(tank_ctc_diff_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size, start_sample);
		tank_ctc_means.push_back(mean_ctc_to_tank);
		tank_ctc_vars.push_back(var_ctc_to_tank);
	//	std::cout <<"mean_ctc_to_tank: "<<mean_ctc_to_tank<<", var_ctc_to_tank: "<<var_ctc_to_tank<<std::endl;
		start_sample += step_size;
		//std::cout <<"start_sample: "<<start_sample<<std::endl;
		if (start_sample < (int) ctc_tvals.size()) tank_ctc_ts.push_back(static_cast<double>(ctc_tvals.at(start_sample)));
	}
	std::vector<double> tank_ctc_binwidths(tank_ctc_diff_vals.size(),step_size);
	std::cout <<"tcerr"<<std::endl;
	TGraphErrors* tcerr = new TGraphErrors(tank_ctc_means.size(), tank_ctc_ts.data(),tank_ctc_means.data(),tank_ctc_binwidths.data(),tank_ctc_vars.data());
	std::string title="(CTC-Tank) Mean and Variances with window size "+std::to_string(window_size)+" and overlap fraction "+std::to_string(int(overlap_fraction));
	//tcerr->SetName(title.c_str());
	tcerr->SetTitle(title.c_str());
	
	outfile->cd();
	tcerr->Draw("apl");
	tcerr->GetYaxis()->SetTitle("#Delta t [ns]");
	tcerr->GetYaxis()->SetTitleSize(0.035);
	tcerr->GetYaxis()->SetTitleOffset(1.3);
	tcerr->GetXaxis()->SetTimeDisplay(1);
	tcerr->GetXaxis()->SetLabelSize(0.03);
	tcerr->GetXaxis()->SetLabelOffset(0.03);
	tcerr->GetXaxis()->SetTimeDisplay(1);
	tcerr->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	//tcerr->GetXaxis()->SetTimeOffset(t1);
	tcerr->GetXaxis()->SetTimeOffset(0);
	tcerr->Write("tcerr");

//	ComputeMeanAndVariance(ctc_to_tank_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size);
//	ComputeMeanAndVariance(ctc_to_tank_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size);
	
	// 5. you could also make a normalized histogram at each step to make a colour band plot
	
	return true;
}
