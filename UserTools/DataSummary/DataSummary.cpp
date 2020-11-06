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
	InputFilePattern="";
	StartRun=-1;
	StartSubRun=-1;
	StartPart=-1;
	EndRun=-1;
	EndSubRun=-1;
	EndPart=-1;
	OutputFileDir=".";
	OutputFileName="DataSummary.root";
	
	// read the user's preferences
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("DataPath",DataPath);
	m_variables.Get("InputFilePattern",InputFilePattern);
	m_variables.Get("StartRun",StartRun);
	m_variables.Get("StartSubRun",StartSubRun);
	m_variables.Get("StartPart",StartPart);
	m_variables.Get("EndRun",EndRun);
	m_variables.Get("EndSubRun",EndSubRun);
	m_variables.Get("EndPart",EndPart);
	m_variables.Get("OutputFileDir",OutputFileDir);
	m_variables.Get("OutputFileName",OutputFileName);
	
	// scan for matching input files
	int numfilesfound = ScanForFiles(DataPath,InputFilePattern);
	if(numfilesfound==0){
		Log("No matching files found!",v_error,verbosity);
		return false;
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
		// i got nothin' ¯\_(ツ)_/¯
		return true;
	}
	
	if(got_annieevent){
		Log("DataSummary Tool: Getting information from ANNIEEvent",v_debug,verbosity);
		
		// variables we can directly retrieve
		ANNIEEvent->Get("EventNumber",EventNumber);
		ANNIEEvent->Get("EventTimeTank",PMTtimestamp);
		ANNIEEvent->Get("CTCTimestamp",CTCtimestamp);
		ANNIEEvent->Get("EventTime",mrd_timeclass);       // convertme to MRDtimestamp
		ANNIEEvent->Get("TriggerWord",TriggerWord);       // convert to TriggerTypeString
		ANNIEEvent->Get("MRDLoopbackTDC",MRDLoopbackTDC); // convert to LoopbackTimestamp values
		// TODO optional sanity checks: consistency of RunNumber and other constants
		
		// calculated variables
		MRDtimestamp = mrd_timeclass.GetNs();
		// extract out the TDC vals
		int beamloopbackTDCticks = MRDLoopbackTDC.at("BeamLoopbackTDC");
		int cosmicloopbackTDCticks = MRDLoopbackTDC.at("CosmicLoopbackTDC");
		// convert to ns
		BeamLoopbackTimestamp = 4000. - 4.*(double)beamloopbackTDCticks;
		CosmicLoopbackTimestamp = 4000. - 4.*(double)cosmicloopbackTDCticks;
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
		
		Log("DataSummary Tool: Filling Event tree",v_debug,verbosity);
		outtree->Fill();
	}
	
	if(got_orphan){
		Log("DataSummary Tool: Getting information from OrphanStore",v_debug,verbosity);
		// right now this is just a straight translation from BoostStore to ROOT =/
		OrphanStore->Get("EventType",orphantype);
		OrphanStore->Get("Timestamp",orphantimestamp);
		OrphanStore->Get("Reason",orphancause);
		
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

int DataSummary::ScanForFiles(std::string inputdir, std::string filepattern){
	// Scan the input directory for all files matching the specified pattern and run range
	
	// Step 1:
	// first, find matching files
	// TODO we should give booststores an extension: for now insist R*S*p* comes at the end of the filename
	std::string lscommand = std::string("find -regextype egrep -iregex '.*?")
	   + InputFilePattern + std::string("R([0-9]+)S([0-9]+)[pP]([0-9]+)$'");
	std::string fileliststring = GetStdoutFromCommand(lscommand);
	std::stringstream ssl;
	ssl << fileliststring;
	std::vector<std::string> flist;
	std::string nextfilestring;
	while(getline(ssl,nextfilestring)){
		flist.push_back(nextfilestring);
	}
	
	// now parse the list of files and extract those that the range of run/subrun numbers
	filelist.clear();
	for(auto afile : flist){
		//try{
			// extract run, subrun and part using a regex match
			std::smatch submatches;
			std::regex theexpression (".*?DataR([0-9]+)S([0-9]+)[pP]([0-9]+).*");
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
			&& (run<EndRun || EndRun<0)
			&& (subrun<EndSubRun || EndSubRun<0)
			&& (part<EndPart || EndPart<0)
		){
			// ensure files are in the correct order by using a map with suitable key
			// assume no more than 1000 parts per subrun, 1000 subruns per run, 1000000 runs
			char buffer [13];
			snprintf(buffer, 13, "%06d%03d%03d", run, subrun, part);
			filelist.emplace(buffer,afile);
		}
	}
	
	nextfile=filelist.begin();
	return filelist.size();
}

bool DataSummary::LoadNextANNIEEventEntry(){
	// load next ANNIEEvent entry
	localentry++;
	globalentry++;
	if((localentry>=localentries)&&(localorphan>=localorphans)){
		// end of this file, try to load next file
		bool load_ok = LoadNextFile();
		if(not load_ok){
			// failed to load any further files. Terminate.
			m_variables.Set("StopLoop",true);
			return false;
		}
		return ANNIEEvent->GetEntry(localentry);
	} else if(localentry<localentries){
		return ANNIEEvent->GetEntry(localentry);
	} // else no more ANNIEEvents, but we didn't load a new file
	  // because there are still Orphans to process
	return false;
}

bool DataSummary::LoadNextOrphanStoreEntry(){
	// load next OrphanStore entry
	localorphan++;
	globalorphan++;
	// if we've run out of entries in both stores, load the next file
	if((localorphan>=localorphans)&&(localentry>=localentries)){
		bool load_ok = LoadNextFile();
		if(not load_ok){
			// failed to load any further files. Terminate.
			m_variables.Set("StopLoop",true);
			return false;
		}
		return ANNIEEvent->GetEntry(localentry);
	} else if(localorphan<localorphans){
		// not a new file, still entries to process
		return OrphanStore->GetEntry(localorphan);
	} // else no more orphans, but we didn't load a new file
	  // because there are still ANNIEEvents to process
	return false;
}

bool DataSummary::LoadNextFile(){
	// get next filename
	if(nextfile==filelist.end()){
		// no more files
		return false;
	}
	std::string nextfilename = nextfile->second;
	
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
	ProcessedFileStore = new BoostStore(false,BOOST_STORE_BINARY_FORMAT);
	// Load the contents from the new input file into it
	std::cout <<"Reading in current file "<<nextfilename<<std::endl;
	int load_ok = ProcessedFileStore->Initialise(nextfilename);
	if(not load_ok) return false;
	
	// create an ANNIEEvent BoostStore and an OrphanStore BoostStore to load from it
	ANNIEEvent = new BoostStore(false, BOOST_STORE_MULTIEVENT_FORMAT);
	OrphanStore = new BoostStore(false, BOOST_STORE_MULTIEVENT_FORMAT);
	
	// retrieve the multi-event stores
	ProcessedFileStore->Get("ANNIEEvent",*ANNIEEvent);
	ANNIEEvent->Header->Get("TotalEntries", localentries);
	
	// same for Orphan Store
	ProcessedFileStore->Get("OrphanStore",*OrphanStore);
	OrphanStore->Header->Get("TotalEntries", localorphans);
	
	nextfile++;
	localentry = 0;
	localorphan = 0;
	
	// load run constants and sanity checks
	get_ok = ANNIEEvent->GetEntry(0);
	if(not get_ok){
		Log("DataSummary Tool: Error getting ANNIEEvent entry 0 from new file "+nextfilename,v_error,verbosity);
		return false;
	}
	ANNIEEvent->Get("RunNumber",RunNumber);
	if((int)RunNumber!=run)
		Log("DataSummary Tool: filename / entry mismatch for RunNumber!",v_error,verbosity);
	ANNIEEvent->Get("SubrunNumber",SubrunNumber);
	if((int)SubrunNumber!=subrun)
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
	outtree->Branch("CTCtimestamp",&CTCtimestamp);
	outtree->Branch("PMTtimestamp",&PMTtimestamp);
	outtree->Branch("MRDtimestamp",&MRDtimestamp);
	outtree->Branch("BeamLoopbackTimestamp",&BeamLoopbackTimestamp);
	outtree->Branch("CosmicLoopbackTimestamp",&CosmicLoopbackTimestamp);
	outtree->Branch("SystemsPresent",&SystemsPresent);
	outtree->Branch("LoopbacksPresent",&LoopbacksPresent);
	outtree->SetAlias("CtcToTankTDiff","CTCtimestamp-PMTtimestamp");
	outtree->SetAlias("CtcToMrdTDiff","CTCtimestamp-MRDtimestamp");
	outtree->SetAlias("TankToMrdTDiff","PMTtimestamp-MRDtimestamp");
	outtree->SetAlias("CtcToBeamLoopbackTDiff","CTCtimestamp-BeamLoopbackTimestamp");
	outtree->SetAlias("CtcToCosmicLoopbackTDiff","CTCtimestamp-CosmicLoopbackTimestamp");
	
	outtree2 = new TTree("OrpahStats","OrpahStats Tree");
	outtree2->Branch("OrphanedEventType",&orphantype);
	outtree2->Branch("OrphanTimestamp",&orphantimestamp);
	outtree2->Branch("OrphanCause",&orphancause);
	return true;
}

bool DataSummary::CreatePlots(){
	// make plots from the ROOT tree
	// =============================
	// first get the timespan we're going to be plotting
	outtree->GetEntry(0);
	t0 = CTCtimestamp;
	outtree->GetEntry(outtree->GetEntriesFast()-1);
	tn = CTCtimestamp;
	// we need to set a global "time offset", from which all other time axes will be relative to
	// this should be in seconds
	gStyle->SetTimeOffset(t0/1E9);
	
	// rate plots: should be able to make these just by time-x-axis histograms
	// TODO refactor?
	AddRatePlots(500);
	
	// TODO refactor this
	// time evolution of the timestamp discrepancies
	AddTDiffPlots();
	//AddTDiffPlots(30,30,0.5,"CtcToMrdTDiff");
	//AddTDiffPlots(30,30,0.5,"TankToMrdTDiff");
	return true;	
}

// just histograms of timestamps of a given type with a time axis
bool DataSummary::AddRatePlots(int nbins){
	outfile->cd(); // ensure plots get put in the file
	
	// define our histogram min/max relative to global time offset, in seconds
	int t1 = 0;
	int t2 = (tn-t0)/1E9;
	
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
	outtree->Draw("CTCtimestamp>>ctc_rate","CTCtimestamp>0");
	TH1D* tank_rate = new TH1D("tank_rate","Tank Rate",nbins,t1,t2);
	outtree->Draw("PMTtimestamp>>tank_rate","PMTtimestamp>0");
	TH1D* mrd_rate = new TH1D("mrd_rate","MRD Rate",nbins,t1,t2);
	outtree->Draw("MRDtimestamp>>mrd_rate","MRDtimestamp>0");
	TH1D* allsystems_rate = new TH1D("allsystems_rate","All Systems Triggered Rate",nbins,t1,t2);
	outtree->Draw("CTCtimestamp>>allsystems_rate","MRDtimestamp>0&&CTCtimestamp>0&&PMTtimestamp>0");
	TH1D* noloopback_rate = new TH1D("noloopback_rate","All Systems Triggered Rate",nbins,t1,t2);
	outtree->Draw("CTCtimestamp>>noloopback_rate","LoopbacksPresent==0");
	TH1D* orphan_rate = new TH1D("orphan_rate","Orphan Rate",nbins,t1,t2);
	outtree2->Draw("OrphanTimestamp>>orphan_rate");
	
	// then we need to set the axes to time
	ctc_rate->GetXaxis()->SetTimeDisplay(1);
	tank_rate->GetXaxis()->SetTimeDisplay(1);
	mrd_rate->GetXaxis()->SetTimeDisplay(1);
	allsystems_rate->GetXaxis()->SetTimeDisplay(1);
	noloopback_rate->GetXaxis()->SetTimeDisplay(1);
	orphan_rate->GetXaxis()->SetTimeDisplay(1);
	// set the timestamp format
	ctc_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	tank_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	mrd_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	allsystems_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	noloopback_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	orphan_rate->GetXaxis()->SetTimeFormat("%d %H:%M");
	// make labels small enough that we don't overcrowd the axis
	ctc_rate->GetXaxis()->SetLabelSize(0.03);
	tank_rate->GetXaxis()->SetLabelSize(0.03);
	mrd_rate->GetXaxis()->SetLabelSize(0.03);
	allsystems_rate->GetXaxis()->SetLabelSize(0.03);
	noloopback_rate->GetXaxis()->SetLabelSize(0.03);
	orphan_rate->GetXaxis()->SetLabelSize(0.03);
	return true;	
}

// TODO refactor to break up plot types and remove triplets of calls
bool DataSummary::AddTDiffPlots(){
	// make the file the active directory so our histograms get put there
	outfile->cd();
	// TODO we'll do a few plots here:
	// 1. we can make a normal histogram showing the distribution of timestamp differences
	TH1D* all_tank_ctc_tdiffs = new TH1D("all_tank_ctc_tdiffs","All (CTC-Tank) TDiffs",100,-100,100);
	outtree->Draw("CtcToTankTDiff>>all_tank_ctc_tdiffs");
	// btw we can grab that plotted data for later
	double* datapointer = outtree->GetV1();
	std::vector<double> tank_ctc_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	TH1D* all_tank_mrd_tdiffs = new TH1D("all_tank_mrd_tdiffs","All (Tank-MRD) TDiffs",100,-100,100);
	datapointer = outtree->GetV1();
	std::vector<double> mrd_ctc_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	outtree->Draw("TankToMrdTDiff>>all_tank_mrd_tdiffs");
	TH1D* all_mrd_ctc_tdiffs = new TH1D("all_ctc_mrd_tdiffs","All (CTC-MRD) TDiffs",100,-100,100);
	outtree->Draw("CtcToMrdTDiff>>all_tank_mrd_tdiffs");
	datapointer = outtree->GetV1();
	std::vector<double> mrd_tank_diff_vals(datapointer, datapointer + outtree->GetSelectedRows());
	
	// 2. we could do the same for the above, broken down by run
	
	// 3. we could also make 2D scatter versions, where we can see how the differences evolve over time
	// define our histogram min/max - these are relative to the global time offset, in seconds
	int t1 = 0;
	int t2 = (tn-t0)/1E9;
	// TGraph* tank_ctc_diffs = new TGraph("tank_ctc_diffs","All (CTC-Tank) TDiffs",100,t1,t2);
	// apparently there really isn't a way to draw from a TTree straight into a named TGraph
	outtree->Draw("CTCtimestamp:CtcToTankTDiff","CTCtimestamp>0&&PMTtimestamp>0");
	datapointer = outtree->GetV2();
	std::vector<double> ctc_tvals(datapointer, datapointer + outtree->GetSelectedRows());
	TGraph* tank_ctc_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	tank_ctc_diffs->SetName("tank_ctc_diffs");
	tank_ctc_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;CTC - Tank [s]");
	
	outtree->Draw("CTCtimestamp:CtcToMrdTDiff","CTCtimestamp>0&&MRDtimestamp>0");
	TGraph* mrd_ctc_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	mrd_ctc_diffs->SetName("mrd_ctc_diffs");
	mrd_ctc_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;CTC - MRD [s]");
	
	outtree->Draw("PMTtimestamp:TankToMrdTDiff","PMTtimestamp>0&&MRDtimestamp>0");
	TGraph* tank_mrd_diffs = (TGraph*)gPad->GetPrimitive("Graph");
	tank_mrd_diffs->SetName("tank_mrd_diffs");
	tank_mrd_diffs->SetTitle("All (CTC-Tank) TDiffs;Time;Tank - MRD [s]");
	
	// then we need to set the axes to time
	tank_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	mrd_ctc_diffs->GetXaxis()->SetTimeDisplay(1);
	tank_mrd_diffs->GetXaxis()->SetTimeDisplay(1);
	// set the timestamp format
	tank_ctc_diffs->GetXaxis()->SetTimeFormat("%d %H:%M");
	mrd_ctc_diffs->GetXaxis()->SetTimeFormat("%d %H:%M");
	tank_mrd_diffs->GetXaxis()->SetTimeFormat("%d %H:%M");
	// make labels small enough that we don't overcrowd the axis
	tank_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	mrd_ctc_diffs->GetXaxis()->SetLabelSize(0.03);
	tank_mrd_diffs->GetXaxis()->SetLabelSize(0.03);
	
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
	for(int i=0; i<(int)tank_ctc_diff_vals.size(); ++i){
		ComputeMeanAndVariance(tank_ctc_diff_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size, start_sample);
		tank_ctc_means.push_back(mean_ctc_to_tank);
		tank_ctc_vars.push_back(var_ctc_to_tank);
		start_sample += step_size;
		tank_ctc_ts.push_back(static_cast<double>(ctc_tvals.at(start_sample)));
	}
	std::vector<double> tank_ctc_binwidths(tank_ctc_diff_vals.size(),step_size);
	TGraphErrors* tcerr = new TGraphErrors(tank_ctc_means.size(), tank_ctc_ts.data(),tank_ctc_means.data(),tank_ctc_binwidths.data(),tank_ctc_vars.data());
	std::string title="(CTC-Tank) Mean and Variances with window size "+std::to_string(window_size)+" and overlap fraction "+std::to_string(int(overlap_fraction));
	//tcerr->SetName(title.c_str());
	tcerr->SetTitle(title.c_str());
	
//	ComputeMeanAndVariance(ctc_to_tank_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size);
//	ComputeMeanAndVariance(ctc_to_tank_vals, mean_ctc_to_tank, var_ctc_to_tank, window_size);
	
	// 5. you could also make a normalized histogram at each step to make a colour band plot
	return true;
	
}
