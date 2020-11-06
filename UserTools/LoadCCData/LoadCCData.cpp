/* vim:set noexpandtab tabstop=4 wrap */

#include "LoadCCData.h"
#include "PMTData.h"
#include "MRDTree.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "THStack.h"

#include <limits>
#include <time.h>
//#include <unordered_map>
//#include <algorithm>

LoadCCData::LoadCCData():Tool(){}

// The MRDTimeStamp is a timestamp of when the TDC readout for that trigger
// ended, which is pretty inaccurate. We need to match it with a more accurate version
// from the ADC data.
// This Tool performs both loading of TDC data, and timestamp alignment.
// See Execute for more information.

namespace {
	constexpr uint64_t MRD_NS_PER_SAMPLE=4;
	constexpr uint64_t MS_TO_NS = 1000000;         // needed for TDC to convert timestamps to ns
	constexpr uint64_t FIVE_HRS_IN_MS = 18000000;  // MRD timezone is off, so timestamps are off by 5 hours
//	constexpr uint64_t NS_TO_SECONDS = 1000000000; // now defined in TimeClass
}

bool LoadCCData::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Usefull header ///////////////////////
	if(configfile!="")  m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	map_version = "v2";
	m_data= &data; //assigning transient data pointer
	// XXX Cannot run Log before we've retrieved m_data!! XXX
	m_variables.Get("verbose",verbosity);
	m_variables.Get("ChannelMapVersion",map_version);	//Variable is introduced to distinguish between the mapping that was present in this tool (v1) and new mapping that seems to work better for the veto efficiency studies (v2)
	Log("LoadCCData Tool: Initializing Tool LoadCCData",v_message,verbosity);

	if (map_version != "v1" && map_version != "v2"){
		Log("LoadCCData Tool: ChannelMapVersion "+map_version+" not recognized. Use v2 as default.",v_error,verbosity);
		map_version = "v2";
	}
	
	// get the files being processed for retrieving the CCData
	std::string events_file;
	get_ok = m_data->CStore.Get("InputFile",events_file);
	if(not get_ok){
		Log("LoadCCData Tool: No InputFile in CStore, expected to be loaded by RawLoader Tool",
			v_error,verbosity);
		return false;
	}
	cout<<"events_file="<<events_file<<endl;
	
	// files used for minibuffer timestamps need to be the same as those loaded by upstream tools
	// to correctly map them.
	std::string timing_file;
	get_ok = m_data->CStore.Get("UsingHeftyMode",useHeftyTimes);
	if(not get_ok){
		Log("LoadCCData Tool: No UsingHeftyMode in CStore, expected to be loaded by RawLoader Tool",
			v_error,verbosity);
		return false;
	}
	if(useHeftyTimes) m_data->CStore.Get("HeftyTimingFile",timing_file);
	else timing_file = events_file;
	cout<<"timing_file="<<timing_file<<endl;
	
	// maximum time difference for a match
	get_ok = m_variables.Get("MaxTimeDiff",maxtimediff);
	if(not get_ok) maxtimediff=std::numeric_limits<uint64_t>::max(); // no limit
	
	// check if making debug files
	m_variables.Get("DrawDebugGraphs",DEBUG_DRAW_TDC_HITS);
	
	// Check the ANNIEEvent Store exists
	// =================================
	int annieeventexists = m_data->Stores.count("ANNIEEvent");
	if(annieeventexists==0){
		Log("LoadCCData Tool: No ANNIEEvent, expected to be loaded by RawLoader Tool",
			v_error,verbosity);
		return false;
	}
	
	// ANNIEEvent TDC hit storage
	// ==========================
	TDCData = new std::map<unsigned long, std::vector<std::vector<Hit>>>;
	
	// Prepare the CC Data files for Reading
	// =====================================
	// 1. Load files from the file list
	std::vector<std::string> filevector;
	// 'filelist' was originally read from LoadCCDataConfig, and is the path to a file in which each line
	// is a file to process. Since we need to load the same files as upstream, it's better not to have
	// to specify files twice, but upstream doesn't provide support for processing multiple files....yet.
//	std::string line;
//	ifstream myfile (filelist.c_str());
//	if (myfile.is_open()){
//		while ( getline (myfile,line) ){
//			filevector.push_back(line);
//		}
//		myfile.close();
//	} else {
//		Log("LoadCCData Tool: Error: Bad CCData filename!",v_error,verbosity);
//		return false;
//	}
	filevector.push_back(events_file);
	
	// 2 : Construct the TChain
	// tree name is "CCData" in RAW files (e.g. DataR900S4p0.root),
	// but "MRDData" in post-processed files (e.g. V6DataR889S0p0T5Sep_25_1650.root).
	// In either case the tree formats are the same, so the Tool should work on both.
	// Use the first file in the list to determine the file type:
	TFile* tempfile = TFile::Open(filevector.front().c_str());
	bool hasccdata = tempfile->GetListOfKeys()->Contains("CCData");
	bool hasmrddata = tempfile->GetListOfKeys()->Contains("MRDData");
	tempfile->Close();
	
	// construct the appropriate TChain
	     if(hasccdata)  { MRDChain = new TChain("CCData");  }
	else if(hasmrddata) { MRDChain = new TChain("MRDData"); }
	else {
		Log("LoadCCData Tool: ERROR: First TDC file contains neither CCData nor MRDData Tree!",v_error,verbosity);
		return false;
	}
	// add the files to the TChain
	for(std::string afile : filevector){
		Log("LoadCCData Tool: Loading file "+afile,v_message,verbosity);
		int filesadded = MRDChain->Add(afile.c_str());
		Log(to_string(filesadded)+" files added to Chain",v_debug,verbosity);
	}
	
	// 3. Construct a Reader class from the TChain
	Log("LoadCCData Tool: MRDChain has "+to_string(MRDChain->GetEntries())+" entries",v_debug,verbosity);
	TTree* testMRD=MRDChain->CloneTree();
	if(testMRD){
		MRDData=new MRDTree(MRDChain->CloneTree());
		MRDData->fChain->SetName("MRDData");
		//m_data->AddTTree("MRDData",MRDData->fChain);
	} else {
		Log("LoadCCData Tool: MRDChain CloneTree failure!",v_error,verbosity);
		return false;
	}
	NumCCDataEntries = MRDChain->GetEntries();
	TDCChainEntry=0;
	
	// Prepare the ADC Data files for reading
	// ======================================
	// 1. Make the TChain
	// we may either use the raw PMTData files, or the heftydb timestamp files
	if(not useHeftyTimes) { ADCTimestampChain = new TChain("PMTData"); }
	else                  { ADCTimestampChain = new TChain("heftydb"); }
	
	// 2. Load files from the file list into the TChain
	Log("LoadCCData Tool: Loading files for timestamp alignment",v_message,verbosity);
	// as per the CCData file list, a list of hefty_timing files was given via an
	// 'alignmentfiles' file, passed in via LoadCCDataConfig
//	myfile.open(alignmentfiles.c_str());
//	if (myfile.is_open()){
//		while ( getline (myfile,line) ){
//			Log("LoadCCData Tool: Loading file "+line+" for PMTData",v_debug,verbosity);
//			int filesadded = ADCTimestampChain->Add(line.c_str());
//			Log(to_string(filesadded)+" files added to alignment Chain",v_debug,verbosity);
//		}
//		myfile.close();
//	} else {
//		Log("LoadCCData Tool: Error: Bad alignment filename "+line+"!",v_error,verbosity);
//		return false;
//	}
	Log("LoadCCData Tool: Loading file "+timing_file+" for PMTData",v_debug,verbosity);
	int filesadded = ADCTimestampChain->Add(timing_file.c_str());
	Log(to_string(filesadded)+" files added to alignment Chain",v_debug,verbosity);
	
	// 3. Construct a Reader class from the TChain
	logmessage = "LoadCCData Tool: Aligning CCData timestamps with ";
	logmessage += ((useHeftyTimes) ? "HeftyData timestamps" : "PMTData timestamps");
	Log(logmessage,v_debug,verbosity);
	if(not useHeftyTimes){ thePMTData   = new                PMTData(ADCTimestampChain); }
	else                 { theHeftyData = new annie::HeftyTreeReader(ADCTimestampChain); }
	NumADCEntries = ADCTimestampChain->GetEntries();
	ADCChainEntry=0;
	
	// Pre-load the first TDC entry
	///////////////////////////////
	Log("LoadCCData Tool: Loading entry "+to_string(TDCChainEntry),v_debug,verbosity);
	Long64_t mentry = MRDData->LoadTree(TDCChainEntry);
	if(mentry < 0){
		Log("LoadCCData Tool: Ran off end of TChain! Stopping ToolChain",v_error,verbosity);
		if(TDCChainEntry==NumCCDataEntries) m_data->vars.Set("StopLoop",1);
		return false;
	} else {
		Log("LoadCCData Tool: Getting TChain localentry " + to_string(mentry),v_debug,verbosity);
		MRDData->GetEntry(mentry);
		TDCChainEntry++;
		if(TDCChainEntry==NumCCDataEntries){
			Log("LoadCCData Tool: Last entry in CCData TChain, stopping chain",v_message,verbosity);
			m_data->vars.Set("StopLoop",1);
		}
	}
	
	// Pre-load the first advance ADC entry
	///////////////////////////////////////
	nextminibufTs.resize(2); // needs at least 2 entries
	LoadPMTDataEntries();
	
	// For making debug Histograms/ROOT output file
	///////////////////////////////////////////////
	if(DEBUG_DRAW_TDC_HITS){
		// get or create the ROOT application to show histograms
		Log("LoadCCData Tool: getting/making TApplication",v_debug,verbosity);
		int myargc=0;
		//char *myargv[] = {(const char*)"mrddist"};
		// get or make the TApplication
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			Log("LoadCCData Tool: Making global TApplication",v_error,verbosity);
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			Log("LoadCCData Tool: Retrieving global TApplication",v_error,verbosity);
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		
		hTDCHitTimes = new TH1D("hTDCHitTimes","TDC Hit Times in Readout",100,0,4200); // phase 1: 4us readout
		hTDCValues = new TH1D("hTDCValues","TDC Tick Values in Readout",100,0,1100);
		hTDCTimeDiffs = new TH1D("hTDCTimeDiffs",
			"Time Diff between TDC Readout Time and Closest Minibuffer Timestamp [ms]",500,-2,5);
		hTDCLastTimeDiffs = new TH1D("hTDCLastTimeDiffs",
			"Time Diff between TDC Readout Time and Last Minibuffer Timestamp [ms]",200,-100,100);
		hTDCNextTimeDiffs = new TH1D("hTDCNextTimeDiffs",
			"Time Diff between TDC Readout Time and Next Minibuffer Timestamp [ms]",200,-100,100);
		
		hVetoL1Times = new TH1D("hVetoL1Times","VetoL1 Hit Times in Readout",100,0,4200);
		hVetoL2Times = new TH1D("hVetoL2Times","VetoL2 Hit Times in Readout",100,0,4200);
		hMrdL1Times = new TH1D("hMrdL1Times","MRDL1 Hit Times in Readout",100,0,4200);
		hMrdL2Times = new TH1D("hMrdL2Times","MRDL2 Hit Times in Readout",100,0,4200);
		hVetoL1Times->SetLineColor(kRed);
		hVetoL2Times->SetLineColor(kViolet);
		hMrdL1Times->SetLineColor(kBlue);
		hMrdL2Times->SetLineColor(kGreen+1);
		
		tdcDebugRootFileOut = new TFile("tdcDebugRootFileOut.root","RECREATE");
		tdcDebugRootFileOut->cd();
		tdcDebugTreeOut = new TTree("tdcDebugTree","Debug Info of TDC Alignment");
		tdcDebugTreeOut->Branch("CamacSlot",&camacslot);
		tdcDebugTreeOut->Branch("CamacChannel",&camacchannel);
		tdcDebugTreeOut->Branch("MRDPMTXNum",&mrdpmtxnum);
		tdcDebugTreeOut->Branch("MRDPMTYNum",&mrdpmtynum);
		tdcDebugTreeOut->Branch("MRDPMTZNum",&mrdpmtznum);
		tdcDebugTreeOut->Branch("MRDTimeInReadout",&mrdtimeinreadout);
		tdcDebugTreeOut->Branch("MRDTicksInReadout",&mrdticksinreadout);
		tdcDebugTreeOut->Branch("MRDReadoutIndex",&mrdreadoutindex);
		tdcDebugTreeOut->Branch("MRDReadoutTime",&mrdreadouttime);
		tdcDebugTreeOut->Branch("ADCReadoutIndex",&adcreadoutindex);
		tdcDebugTreeOut->Branch("ADCMinibufferInReadout",&adcminibufferinreadout);
		tdcDebugTreeOut->Branch("ADCMinibufferIndex",&adcminibufferindex); // flattened minibuffer numbers
		tdcDebugTreeOut->Branch("ADCMinibufferTime",&adcminibuffertime);
		tdcDebugTreeOut->Branch("TDCADCTimeDiff",&tdcadctimediff);
		tdcDebugTreeOut->Branch("TDCHitNum",&tdchitnum);
		gROOT->cd();
		
		if(useHeftyTimes) debugtimesdump.open("debugheftytimes.csv"); // dump minibuf timestamps to file
		else              debugtimesdump.open("debugrawtimes.csv");
		ExecuteIteration=0;
	}
	
	return true;
}

bool LoadCCData::Execute(){
	
	Log("LoadCCData Tool: Executing",v_debug,verbosity);
	
	/*
	On each Execute() call, we should load the ANNIEEvent with any TDC hits matching each of the
	40 ADC minibuffers currently being processed.
	We assume both streams of timestamps increase monotonically, and that we start processing from
	the start of the chain. We then check if a TDC timestamp better matches this ADC timestamp or the next.
	If it better matches the next, it won't better match any previous ones, and we put it off for the next loop.
	If it better matches this timestamp than the next, it won't better match any further future timestamps.
	So we add it to the set of TDC hits for this event.
	*/
	
	// The ADC timestamps of the current event are available from the RawLoader toolchain.
	// For hefty mode, it's in the 'HeftyInfo' saved in the ANNIEEvent,
	// for non-hefty mode they're in the 'MinibufferTimestamps' saved in the ANNIEEvent.
	// Use the presence/absence of HeftyInfo to check if we're reading Hefty data
	
	// Load the ADC info
	std::vector<unsigned long long> currentminibufts;
	HeftyInfo eventheftyinfo;
	get_ok = m_data->Stores["ANNIEEvent"]->Get("HeftyInfo",eventheftyinfo);
	if(not get_ok){
		// Non-Hefty mode
		Log("LoadCCData Tool: ANNIEEvent had no HeftyInfo - assuming non-Hefty data",v_warning,verbosity);
		//m_data->Stores["ANNIEEvent"]->Print(false);
		std::vector<TimeClass> mb_timestamps; // 1-element vector
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MinibufferTimestamps", mb_timestamps);
		if(not get_ok){
			logmessage = "LoadCCData Tool: ANNIEEvent had no HeftyInfo, nor MinibufferTimestamps! ";
			logmessage += "Cannot load trigger timestamps needed for alignment!";
			Log(logmessage,v_error,verbosity);
			//m_data->Stores["ANNIEEvent"]->Print(false);
			return false;
		}
		currentminibufts.emplace_back(mb_timestamps.front().GetNs());
	} else {
		// Hefty Mode
		// CCData tree is 'flattened' relative to Hefty data -
		// each readout corresponds to a single minibuffer, but we process multiple minibuffers
		// simultaneously. We'll need an internal loop to read as many TDCData entries
		// as there are minibuffers in each ADC Readout.
		Log("LoadCCData Tool: Getting HeftyInfo times",v_debug,verbosity);
		currentminibufts = eventheftyinfo.all_times();
		int numminibuffers = eventheftyinfo.num_minibuffers();
		if((int)currentminibufts.size()!=numminibuffers){
			Log("LoadCCData Tool: HeftyInfo mismatch between num minibuffers "
				+ to_string(numminibuffers) + "and size of returned times vector "
				+ to_string(currentminibufts.size()),v_error,verbosity);
		}
	}
	
//	logmessage = "LoadCCData Tool: first minibuffer time = " + to_string(currentminibufts.front());
//	if(currentminibufts.size()>1){
//		logmessage += ", last minibuffer [" + to_string(currentminibufts.size()) + "] at time " +
//						to_string(currentminibufts.back());
//	}
//	Log(logmessage,v_debug,verbosity);
	
	// Also load the next ADCData entry, so we have access to the first timestamp from the upcoming event
	LoadPMTDataEntries();
	
	// Load all matching TDC hits into the TDCData
	TDCData->clear();
	PerformMatching(currentminibufts);
	
	// Update the TDCData in the ANNIEEvent
	m_data->Stores.at("ANNIEEvent")->Set("TDCData",TDCData,true);
	
	ExecuteIteration++;
	return true;
}

bool LoadCCData::Finalise(){
	Log("LoadCCData Tool: Finalizing Tool",v_message,verbosity);
	
	if(DEBUG_DRAW_TDC_HITS){
		tdcDebugRootFileOut->Write("",TObject::kOverwrite);
		
		Double_t canvwidth = 700;
		Double_t canvheight = 600;
		tdcRootCanvas = new TCanvas("tdcRootCanvas","tdcRootCanvas",canvwidth,canvheight);
		tdcRootCanvas->SetWindowSize(canvwidth,canvheight);
		tdcRootCanvas->cd();
		
		hTDCHitTimes->Draw();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("HitTimeInReadout.png");
		
		// breakdown by location
		THStack htimebreakdown("htimebreakdown","Hit Times in Readout Breakdown");
		htimebreakdown.Add(hVetoL1Times);
		htimebreakdown.Add(hVetoL2Times);
		htimebreakdown.Add(hMrdL1Times);
		htimebreakdown.Add(hMrdL2Times);
		htimebreakdown.Draw();
		tdcRootCanvas->BuildLegend();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("HitTimesInReadoutBreakdown.png");
		
		hTDCValues->Draw();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("HitTicksInReadout.png");
		
		hTDCTimeDiffs->Draw();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("TimeDiffToClosestMinibuf.png");
		
		hTDCLastTimeDiffs->Draw();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("TimeDiffToLastMinibuf.png");
		
		hTDCNextTimeDiffs->Draw();
		tdcRootCanvas->Update();
		tdcRootCanvas->SaveAs("TimeDiffToNextMinibuf.png");
		
		// cleanup
		tdcDebugTreeOut->ResetBranchAddresses();
		tdcDebugRootFileOut->Close();
		if(tdcDebugRootFileOut) delete tdcDebugRootFileOut; tdcDebugRootFileOut=nullptr;
		if(hTDCLastTimeDiffs) delete hTDCLastTimeDiffs; hTDCLastTimeDiffs=nullptr;
		if(hTDCHitTimes) delete hTDCHitTimes; hTDCHitTimes=nullptr;
		if(hTDCNextTimeDiffs) delete hTDCNextTimeDiffs; hTDCNextTimeDiffs=nullptr;
		if(hTDCValues) delete hTDCValues; hTDCValues=nullptr;
		if(hTDCTimeDiffs) delete hTDCTimeDiffs; hTDCTimeDiffs=nullptr;
		if(tdcRootCanvas) delete tdcRootCanvas; tdcRootCanvas=nullptr;
		
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				Log("LoadCCData Tool: deleting gloabl TApplication",v_debug,verbosity);
				delete rootTApp;
			}
		} else if (tapplicationusers>1){
			tapplicationusers--;
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		}
		
		debugtimesdump.close();
	}
	
	if(MRDData) delete MRDData;
	if(thePMTData) delete thePMTData;
	if(theHeftyData) delete theHeftyData;
	
	return true;
}

bool LoadCCData::PerformMatching(std::vector<unsigned long long> currentminibufts){
	// scan the TDC data and return all hits matching the current minibuffers
	
	// loop over minibuffers in this Execute() iteration
	for(int minibufi=0; minibufi<(int)currentminibufts.size(); minibufi++){
		// get this minibuffer's timestamp
		uint64_t theminibufferTS = currentminibufts.at(minibufi);
		if(theminibufferTS==0){
			Log("LoadCCData Tool: Event Minibuffer Time Is 0!",v_warning,verbosity);
		}
		
		// get the next minibuffer's timestamp. In the event that this is the last
		// minibuffer of the current readout, we need to use the timestamp of the
		// first minibuffer from the next readout
		uint64_t thenextminibufferTS;
		if((minibufi+1)<(int)currentminibufts.size()){
			thenextminibufferTS = currentminibufts.at(minibufi+1);
			if(thenextminibufferTS==0){
				//Log("LoadCCData Tool: Event Next Minibuffer Time Is 0!",v_warning,verbosity); // don't report twice
			}
		} else {
			thenextminibufferTS = nextreadoutfirstminibufstart;
		}
		
		// check for non-monotonic increase of timestamps. If ADC timestamps go backwards
		// it screws up our check, so we need to discard the sample
		if(thenextminibufferTS<theminibufferTS){
			// get the next next minibuffer, to see where it lies
			uint64_t nextnextminibufferTS;
			if((minibufi+2)<(int)currentminibufts.size()){
				nextnextminibufferTS = currentminibufts.at(minibufi+2);
			} else if(thenextminibufferTS!=nextreadoutfirstminibufstart){
				nextnextminibufferTS=nextreadoutfirstminibufstart;
			} else {
				nextnextminibufferTS=nextminibufTs.at(1);
			}
			
			if(thenextminibufferTS!=0){  // zero timestamps are reported separately
				logmessage="LoadCCData Tool: Non-monotonic increase of ADC timestamps!";
				if((minibufi+1)<(int)currentminibufts.size()){
					logmessage+= " Minibuffer "+to_string(minibufi+1);
				} else {
					logmessage+= " First minibuffer from next readout";
				}
				logmessage+=" is earlier than minibuffer "+to_string(minibufi)+"/"
					  + to_string(currentminibufts.size())+"!";
				Log(logmessage,v_warning,verbosity);
				
				if(not (verbosity<v_debug)){
					// convert to readable format to print some human readable checks
					uint64_t rounded_seconds = static_cast<uint64_t>(floor(theminibufferTS/NS_TO_SECONDS));
					time_t rawtime = rounded_seconds;
					uint64_t mb_remaining_ns = theminibufferTS - rounded_seconds*NS_TO_SECONDS;
					struct tm * ptm;               // time structure from which components can be retrieved
					ptm = gmtime(&rawtime);        // fill timestructure with UTC time
					struct tm tm_mb = tm(*ptm);    // the struct filled by gmtime is internal, so gets overridden
										           // by subsequent calls. We need to make a copy.
					
					rounded_seconds = static_cast<uint64_t>(floor(thenextminibufferTS/NS_TO_SECONDS));
					rawtime = rounded_seconds;
					uint64_t nmb_remaining_ns = thenextminibufferTS - rounded_seconds*NS_TO_SECONDS;
					ptm = gmtime(&rawtime);
					struct tm tm_nmb = tm(*ptm);
					
					rounded_seconds = static_cast<uint64_t>(floor(nextnextminibufferTS/NS_TO_SECONDS));
					rawtime = rounded_seconds;
					uint64_t nnmb_remaining_ns = nextnextminibufferTS - rounded_seconds*NS_TO_SECONDS;
					ptm = gmtime(&rawtime);
					struct tm tm_nnmb = tm(*ptm);
					
					int bufsize=300;
					char logbuffer[bufsize];
					int discardedcharcount = snprintf(logbuffer, bufsize,
						"this minibuffer at %d/%d/%d %d:%d:%d.%09lu (%lu), "
						"next minibuffer at %d/%d/%d %d:%d:%d.%09lu (%lu), "
						"next next minibuffer at %d/%d/%d %d:%d:%d.%09lu (%lu)",
						tm_mb.tm_year+1900,tm_mb.tm_mon,tm_mb.tm_mday,tm_mb.tm_hour,tm_mb.tm_min,
							tm_mb.tm_sec,mb_remaining_ns,theminibufferTS,
						tm_nmb.tm_year+1900,tm_nmb.tm_mon,tm_nmb.tm_mday,tm_nmb.tm_hour,tm_nmb.tm_min,
							tm_nmb.tm_sec,nmb_remaining_ns,thenextminibufferTS,
						tm_nnmb.tm_year+1900,tm_nnmb.tm_mon,tm_nnmb.tm_mday,tm_nnmb.tm_hour,tm_nnmb.tm_min,
							tm_nnmb.tm_sec,nnmb_remaining_ns,nextnextminibufferTS);
				
					Log(logbuffer,v_debug,verbosity);
				} // verbosity >= debug
			} // don't report 0 timestamps as non-monotonic
			
			// now how to handle this...
			if(nextnextminibufferTS<theminibufferTS){
				// maybe *this* timestamp is bad - the next two are both before it...
				// skip trying to match TDC hits to this timestamp...
				Log("LoadCCData Tool: Skipping attempts to match hits to this minibuffer",v_warning,verbosity);
				continue;
			} else {
				// try to skip the next timestamp
				thenextminibufferTS = nextnextminibufferTS;
				Log("LoadCCData Tool: Using nextnextminibufferTS as thenextminibufferTS",v_debug,verbosity);
			}
		}
		
		// scan forward from last matched TDC sample,
		// looking for all TDC samples that match this minibuffer better than the next
		// (we assume both time series increase monotonically)
		bool endoftdctchain=false;
		do{
			// the last TDC entry we loaded didn't match the previous minibuffer
			// so we have a TDC entry pre-loaded
			
			// 1. Get the TDC readout start time
			ULong64_t MRDTimeStamp = MRDData->TimeStamp; // in unix ms - need to convert ms to ns
			TimeClass MRDEventTime = TimeClass(static_cast<uint64_t>((MRDTimeStamp+FIVE_HRS_IN_MS)*MS_TO_NS));
			
			//logmessage = "LoadCCData Tool: TDC readout [" + to_string(MRDData->Trigger)
			//			 + "] at time " + to_string(MRDEventTime.GetNs());
			//Log(logmessage,v_debug,verbosity);
			
			// Compare this TDC readout time with the current and next minibuffer time.
			// if it better fits this minibuffer time, we'll call it a match
			// need to put the result into a signed type before we can take abs!
			int64_t time_to_last_mb = thelastminibufferTS - MRDEventTime.GetNs();
			int64_t time_to_this_mb = theminibufferTS - MRDEventTime.GetNs();
			int64_t time_to_next_mb = thenextminibufferTS - MRDEventTime.GetNs();
			
			if(abs(time_to_this_mb) < abs(time_to_next_mb)){
				if(abs(time_to_this_mb)<(int)maxtimediff){
					// it's a match!
					
					if(not (verbosity<v_debug)){
						// convert to readable format
						uint64_t rounded_seconds = static_cast<uint64_t>(floor(MRDEventTime.GetNs()/NS_TO_SECONDS));
						time_t rawtime = rounded_seconds;
						uint64_t mrd_remaining_ns = MRDEventTime.GetNs() - rounded_seconds*NS_TO_SECONDS;
						struct tm * ptm;               // time structure from which components can be retrieved
						ptm = gmtime(&rawtime);        // fill timestructure with UTC time
						struct tm tm_mrd = tm(*ptm);   // the struct filled by gmtime is internal and overwritten
													   // by subsequent calls. We need to make a copy.
						
						rounded_seconds = static_cast<uint64_t>(floor(theminibufferTS/NS_TO_SECONDS));
						rawtime = rounded_seconds;
						uint64_t mb_remaining_ns = theminibufferTS - rounded_seconds*NS_TO_SECONDS;
						ptm = gmtime(&rawtime);
						struct tm tm_mb = tm(*ptm);
						
						rounded_seconds = static_cast<uint64_t>(floor(thenextminibufferTS/NS_TO_SECONDS));
						rawtime = rounded_seconds;
						uint64_t nmb_remaining_ns = thenextminibufferTS - rounded_seconds*NS_TO_SECONDS;
						ptm = gmtime(&rawtime);
						
						int bufsize=300;
						char logbuffer[bufsize];
						int discardedcharcount = snprintf(logbuffer, bufsize,
							"matched TDC readout at %d/%d/%d %d:%d:%d.%09lu (%lu) "
							"to minibuffer at %d/%d/%d %d:%d:%d.%09lu (%lu), "
							"c.f. next minibuffer at %d/%d/%d %d:%d:%d.%09lu (%lu)",
							tm_mrd.tm_year+1900,tm_mrd.tm_mon,tm_mrd.tm_mday,tm_mrd.tm_hour,tm_mrd.tm_min,
								tm_mrd.tm_sec,mrd_remaining_ns,MRDEventTime.GetNs(),
							tm_mb.tm_year+1900,tm_mb.tm_mon,tm_mb.tm_mday,tm_mb.tm_hour,tm_mb.tm_min,
								tm_mb.tm_sec,mb_remaining_ns,theminibufferTS,
							ptm->tm_year+1900,ptm->tm_mon,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,
								ptm->tm_sec,nmb_remaining_ns,thenextminibufferTS);
						Log(logbuffer,v_debug,verbosity);
					}
					
					// Add all the hits in this readout to the TDC hits matching the minibuffer
					UInt_t numhitsthisevent = MRDData->OutNumber;
					Log("LoadCCData Tool: Looping over "+to_string(numhitsthisevent)+" TDC hits",v_debug,verbosity);
					
					int usedhiti=0; // since trigger hits are skipped
					for(int hiti=0; hiti<(int)numhitsthisevent; hiti++){
						
						// Get the MRD PMT ID plugged into this TDC Card + Channel
						uint32_t tubeid;
						if (map_version == "v1") tubeid = TubeIdFromSlotChannel(MRDData->Slot->at(hiti),MRDData->Channel->at(hiti),1);
						else if (map_version == "v2") tubeid = TubeIdFromSlotChannel(MRDData->Slot->at(hiti),MRDData->Channel->at(hiti),2);
						if(tubeid==std::numeric_limits<uint32_t>::max()){
							// no matching PMT. This was not an MRD / FACC paddle hit.
							// mostly hits on the last channel of the TDC card, to trigger readout
							
							if(MRDData->Channel->at(hiti)!=31){  // last channel is used to trigger TDC readout
								Log("LoadCCData Tool: Ignoring hit on Slot " + to_string(MRDData->Slot->at(hiti))
									+ ", Channel " + to_string(MRDData->Channel->at(hiti))
									+ ", no matching MRD / FACC PMT", v_warning,verbosity);
							}
							continue;
						}
						// FIXME: when we implement loading raw file geometries
						// we'll need to have some way to map from TDC slot+channel to channelkey
						// Then we can probably replace the use of tubeid as this tool interprets it
						// - i.e. as position XXYYZZ - with Channel::GetPosition 
						unsigned long key = tubeid;
						
						//convert time since TDC Timestamp from ticks to ns
						auto hit_time_ticks = MRDData->Value->at(hiti);
						double hit_time_ns = static_cast<double>(hit_time_ticks) * MRD_NS_PER_SAMPLE;
						
						// construct the hit in the ANNIEEvent TDCData
						// Hit nexthit(tubeid, hit_time_ns, -1.); // charge = -1, not recorded
						logmessage  = "LoadCCData Tool: Hit on tube " + to_string(tubeid);
						logmessage += " at " + to_string(hit_time_ns) + " ns relative to TDC readout time";
						Log(logmessage,v_debug,verbosity);
						
						if(TDCData->count(key)==0){
							TDCData->emplace(key, std::vector<std::vector<Hit>>(currentminibufts.size()));
						}
						TDCData->at(key).at(minibufi).emplace_back(tubeid, hit_time_ns, -1.);
						
						if(DEBUG_DRAW_TDC_HITS){
							// fill debug ROOT tree
							camacslot=MRDData->Slot->at(hiti);
							camacchannel=MRDData->Channel->at(hiti);
							int corrected_tubeid = tubeid-1000000;
							mrdpmtxnum=floor(corrected_tubeid/10000);
							mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
							mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
							//printf("LoadCCData: TDC hit on PMT %02d:%02d:%02d in minibuffer %d\n",
							//       mrdpmtxnum,mrdpmtynum,mrdpmtznum,minibufi);
							mrdtimeinreadout=hit_time_ns;
							mrdticksinreadout=hit_time_ticks;
							mrdreadoutindex=TDCChainEntry;
							mrdreadouttime=MRDEventTime.GetNs();
							adcreadoutindex=ADCChainEntry;
							adcminibufferinreadout=minibufi;
							adcminibufferindex=(ADCChainEntry*currentminibufts.size())+minibufi;
							adcminibuffertime=theminibufferTS;
							tdcadctimediff=time_to_this_mb;
							tdchitnum=usedhiti;
							tdcDebugTreeOut->Fill();
							
							// fill debug histograms
							if(usedhiti==0){
								hTDCLastTimeDiffs->Fill(static_cast<double>(time_to_last_mb)/1000000.);
								hTDCTimeDiffs->Fill(static_cast<double>(time_to_this_mb)/1000000.);
								hTDCNextTimeDiffs->Fill(static_cast<double>(time_to_next_mb)/1000000.);
							}
							hTDCHitTimes->Fill(hit_time_ns);
							hTDCValues->Fill(hit_time_ticks);
							if(mrdpmtznum==0){
								if(mrdpmtxnum==0){
									hVetoL1Times->Fill(hit_time_ns);
								} else {
									hVetoL2Times->Fill(hit_time_ns);
								}
							} else if(mrdpmtznum==2){
								hMrdL1Times->Fill(hit_time_ns);
							} else if(mrdpmtznum==3){
								hMrdL2Times->Fill(hit_time_ns);
							}
						}
						usedhiti++;
					} // end loop over hits this TDC readout
				} else {
					Log("LoadCCData Tool: Skipping TDC readout at "+to_string(MRDEventTime.GetNs())
						+" as time to closest minibuffer time is greater than maximum difference",
						v_message,verbosity);
					if(DEBUG_DRAW_TDC_HITS){
						// add them to the debug file, so we can see what we skipped
						UInt_t numhitsthisevent = MRDData->OutNumber;
						int usedhiti=0; // since trigger hits are skipped
						for(int hiti=0; hiti<(int)numhitsthisevent; hiti++){
							uint32_t tubeid=std::numeric_limits<uint32_t>::max();
							if (map_version == "v1") tubeid = TubeIdFromSlotChannel(MRDData->Slot->at(hiti),MRDData->Channel->at(hiti),1);
                                                	else if (map_version == "v2") tubeid = TubeIdFromSlotChannel(MRDData->Slot->at(hiti),MRDData->Channel->at(hiti),2);
							if(tubeid==std::numeric_limits<uint32_t>::max()){ continue; }
							auto hit_time_ticks = MRDData->Value->at(hiti);
							double hit_time_ns = static_cast<double>(hit_time_ticks) * MRD_NS_PER_SAMPLE;
							// fill debug ROOT tree
							camacslot=MRDData->Slot->at(hiti);
							camacchannel=MRDData->Channel->at(hiti);
							int corrected_tubeid = tubeid-1000000;
							mrdpmtxnum=floor(corrected_tubeid/10000);
							mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
							mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
							mrdtimeinreadout=hit_time_ns;
							mrdticksinreadout=hit_time_ticks;
							mrdreadoutindex=TDCChainEntry;
							mrdreadouttime=MRDEventTime.GetNs();
							adcreadoutindex=ADCChainEntry;
							adcminibufferinreadout=minibufi;
							adcminibufferindex=(ADCChainEntry*currentminibufts.size())+minibufi;
							adcminibuffertime=theminibufferTS;
							tdcadctimediff=time_to_this_mb;
							tdchitnum=hiti;
							tdcDebugTreeOut->Fill();
							usedhiti++;
						}
					}
				}
				
				// Load the next TDC readout
				////////////////////////////
				Log("LoadCCData Tool: Loading entry "+to_string(TDCChainEntry),v_debug,verbosity);
				Long64_t mentry = MRDData->LoadTree(TDCChainEntry);
				if(mentry < 0){
					Log("LoadCCData Tool: Ran off end of TChain! Stopping ToolChain",v_warning,verbosity);
					if(TDCChainEntry==NumCCDataEntries) m_data->vars.Set("StopLoop",1);
					return false;
				} else {
					Log("LoadCCData Tool: Getting TChain localentry " + to_string(mentry),v_debug,verbosity);
					MRDData->GetEntry(mentry);
					TDCChainEntry++;
					if(TDCChainEntry==NumCCDataEntries){
						Log("LoadCCData Tool: Last entry in CCData TChain, stopping chain",v_message,verbosity);
						m_data->vars.Set("StopLoop",1);
						endoftdctchain=true;
					}
				} // next entry loaded, loop back round to compare it to this minibuffer timestamp
				
			} else {
				// This TDC readout better matches the next minibuffer. End of matching TDC readouts
				logmessage="next TDC readout better matches next minibuffer time, done matching\n"
						   "mbTS="+to_string(theminibufferTS)+", nextmbTS="+to_string(thenextminibufferTS)
						  +", mrdTS="+to_string(MRDEventTime.GetNs())
						  +"\ntime_to_this_mb="+to_string(time_to_this_mb)
						  +"\ntime_to_next_mb="+to_string(time_to_next_mb);
				Log(logmessage,v_debug,verbosity);
				break;
			}
			
		} while (not endoftdctchain); // end of do loop looking for matching TDC readouts
		
		thelastminibufferTS = theminibufferTS;
	} // end loop over minibuffers in this Execute()
	
	return true;
}

bool LoadCCData::LoadPMTDataEntries(){
	
	if(ADCChainEntry!=NumADCEntries){
		
		// Load more ADC Data
		//Log("LoadCCData Tool: Getting alignment entry "+to_string(ADCChainEntry)
		//	+" in execution "+to_string(ExecuteIteration),v_debug,verbosity);
		if(not useHeftyTimes){
			// if using Raw PMTData tree, we'll use the annie:RawCard class
			// to convert the timestamps to a useable format
			
			Long64_t pentry = thePMTData->LoadTree(ADCChainEntry);
			thePMTData->GetEntry(pentry);
			// n.b. PMTData::TriggerCounts is an array of size PMTData::TriggerNumber
			// but annie::RawCard takes TriggerCounts as a vector
			std::vector<unsigned long long> TriggerCountsAsVector(thePMTData->TriggerCounts,
				thePMTData->TriggerCounts + thePMTData->TriggerNumber);
			//minibufTS = nextminibufTs;
			nextminibufTs = ConvertTimeStamps(	thePMTData->LastSync,
												thePMTData->StartTimeSec,
												thePMTData->StartTimeNSec,
												thePMTData->StartCount,
												TriggerCountsAsVector );
			nextreadoutfirstminibufstart = nextminibufTs.front();
			
		} else {
			// if using the heftydb tree, we can just get the timestamps
			// in a usable format directly
			std::unique_ptr<HeftyInfo> entryinfo = theHeftyData->next();
			if(entryinfo){
				//minibufTS = nextminibufTs;
				//nextminibufTs.clear();
				//for(auto&& anentry: entryinfo->all_times()) nextminibufTs.push_back(anentry);
				nextminibufTs.at(1)=entryinfo->time(1); // for nextnextminibufferTS
				nextreadoutfirstminibufstart = entryinfo->time(0);
			} else {
				Log("LoadCCData Tool: theHeftyData failed to return HeftyInfo on next() call. Last entry?",
					v_warning,verbosity);
			}
			
		}
		
		if(DEBUG_DRAW_TDC_HITS){
			// dump timestamps to file for debug checking
			debugtimesdump << ADCChainEntry <<endl;
			for(uint64_t ats : nextminibufTs){ debugtimesdump << ats << endl; }
			debugtimesdump << endl;
		}
		
		// note that we've read this set of timestamps
		ADCChainEntry++;
		
	} else {
		// no more entries in the TChain.
		nextreadoutfirstminibufstart = std::numeric_limits<uint64_t>::max();
		nextminibufTs.at(1)=std::numeric_limits<uint64_t>::max();
	}
	
	return true;
}

std::vector<uint64_t> LoadCCData::ConvertTimeStamps(unsigned long long LastSync, int StartTimeSec, int StartTimeNSec, unsigned long long StartCount, std::vector<unsigned long long> TriggerCounts){
	// unused members
	const std::vector<unsigned short> Data{};
	const std::vector<unsigned int> Rates{};
	int Channels = 0;                         // must have Channels == Data.size() / BufferSize
	int BufferSize = TriggerCounts.size();    // must have TriggerCounts.size() == (BufferSize / MiniBufferSize)
	int MiniBufferSize = 1;
	int CardID = 0;
	
	// construct a RawCard so we can access the function to obtain the minibuffer TimeStamps
	annie::RawCard tempcard = annie::RawCard(CardID, LastSync, StartTimeSec,
					StartTimeNSec, StartCount, Channels, BufferSize, MiniBufferSize,
					Data, TriggerCounts, Rates);
	std::vector<uint64_t> alltimestamps;
	// loop over all the minibuffers in the readout and retrieve their timestamps
	for(int minibufi=0; minibufi<(int)TriggerCounts.size(); minibufi++){
		unsigned long long thetimestamp = tempcard.trigger_time(minibufi);
		alltimestamps.push_back(thetimestamp);
	}
	return alltimestamps;
}

uint32_t LoadCCData::TubeIdFromSlotChannel(unsigned int slot, unsigned int channel, int version){
	// map TDC Slot + Channel into a MRD PMT Tube ID
	uint32_t tubeid=-1;
	uint16_t slotchan = static_cast<uint16_t>(slot*100)+static_cast<uint16_t>(channel);
	if ((int)slotchan !=1831 && (int)slotchan != 1731) //std::cout <<"LoadCCData: slotchan = "<<slotchan<<std::endl;
	if (version == 1){
		if(slotchantopmtidv1.count(slotchan)){
			std::string stringPMTID = "1"+slotchantopmtidv1.at(slotchan);
			int iPMTID = stoi(stringPMTID);
			tubeid = static_cast<uint32_t>(iPMTID);
		} else if(channel==31){
			// this is the Trigger card output to force the TDCs to always fire
			tubeid = std::numeric_limits<uint32_t>::max();
		} else {
			//Log("LoadCCData Tool: Unknown TDC Slot + Channel combination: "+to_string(slotchan)
			//	+", no matching MRD PMT ID!",v_message,verbosity);  // reported in use
			tubeid = std::numeric_limits<uint32_t>::max();
		}
	} else if (version == 2){
		if(slotchantopmtidv2.count(slotchan)){
			std::string stringPMTID = "1"+slotchantopmtidv2.at(slotchan);
			int iPMTID = stoi(stringPMTID);
			tubeid = static_cast<uint32_t>(iPMTID);
		} else if(channel==31){
			// this is the Trigger card output to force the TDCs to always fire
			tubeid = std::numeric_limits<uint32_t>::max();
		} else {
			//Log("LoadCCData Tool: Unknown TDC Slot + Channel combination: "+to_string(slotchan)
			//	+", no matching MRD PMT ID!",v_message,verbosity);  // reported in use
			tubeid = std::numeric_limits<uint32_t>::max();
		}




	}
	logmessage="LoadCCData Tool: Calculated TubeId from Slot "+to_string(slot)
				+", channel "+to_string(channel)+" = "+to_string(tubeid);
	//Log(logmessage,v_debug,verbosity);
	return tubeid;
}

// PMT ID is constructed from X, Y, Z position:
// e.g. x=0, y=13, z=5 would become '00'+'13'+'05' = 001305.
// since c++ insists on handling integer literals with leading zeros as octal,
// map them as string, then convert the string to int
// phase 1 map:
// Introduce v2 version of the map that seems to fit better with what we see in the veto channels. Old v1 version is kept below.
std::map<uint16_t,std::string> LoadCCData::slotchantopmtidv2{
	std::pair<uint16_t,std::string>{1701u,"000002"},
	std::pair<uint16_t,std::string>{1702u,"000102"},
	std::pair<uint16_t,std::string>{1703u,"000202"},
	std::pair<uint16_t,std::string>{1704u,"000302"},
	std::pair<uint16_t,std::string>{1705u,"000402"},
	std::pair<uint16_t,std::string>{1706u,"000502"},
	std::pair<uint16_t,std::string>{1707u,"000602"},
	std::pair<uint16_t,std::string>{1708u,"000702"},
	std::pair<uint16_t,std::string>{1709u,"000802"},
	std::pair<uint16_t,std::string>{1710u,"000902"},
	std::pair<uint16_t,std::string>{1711u,"001002"},
	std::pair<uint16_t,std::string>{1712u,"001102"},
	std::pair<uint16_t,std::string>{1713u,"001202"},
	std::pair<uint16_t,std::string>{1714u,"010002"},
	std::pair<uint16_t,std::string>{1715u,"010102"},
	std::pair<uint16_t,std::string>{1716u,"010202"},
	std::pair<uint16_t,std::string>{1717u,"010302"},
	std::pair<uint16_t,std::string>{1718u,"010402"},
	std::pair<uint16_t,std::string>{1719u,"010502"},
	std::pair<uint16_t,std::string>{1720u,"010602"},
	std::pair<uint16_t,std::string>{1721u,"010702"},
	std::pair<uint16_t,std::string>{1722u,"010802"},
	std::pair<uint16_t,std::string>{1723u,"010902"},
	std::pair<uint16_t,std::string>{1724u,"011002"},
	std::pair<uint16_t,std::string>{1725u,"011102"},
	std::pair<uint16_t,std::string>{1726u,"011202"},
	std::pair<uint16_t,std::string>{1801u,"000003"},
	std::pair<uint16_t,std::string>{1802u,"010003"},
	std::pair<uint16_t,std::string>{1803u,"020003"},
	std::pair<uint16_t,std::string>{1804u,"030003"},
	std::pair<uint16_t,std::string>{1805u,"040003"},
	std::pair<uint16_t,std::string>{1806u,"050003"},
	std::pair<uint16_t,std::string>{1807u,"060003"},
	std::pair<uint16_t,std::string>{1808u,"070003"},
	std::pair<uint16_t,std::string>{1809u,"080003"},
	std::pair<uint16_t,std::string>{1810u,"090003"},
	std::pair<uint16_t,std::string>{1811u,"100003"},
	std::pair<uint16_t,std::string>{1812u,"110003"},
	std::pair<uint16_t,std::string>{1813u,"120003"},
	std::pair<uint16_t,std::string>{1814u,"130003"},
	std::pair<uint16_t,std::string>{1815u,"000103"},
	std::pair<uint16_t,std::string>{1816u,"010103"},
	std::pair<uint16_t,std::string>{1817u,"020103"},
	std::pair<uint16_t,std::string>{1818u,"030103"},
	std::pair<uint16_t,std::string>{1819u,"040103"},
	std::pair<uint16_t,std::string>{1820u,"050103"},
	std::pair<uint16_t,std::string>{1821u,"060103"},
	std::pair<uint16_t,std::string>{1822u,"070103"},
	std::pair<uint16_t,std::string>{1823u,"080103"},
	std::pair<uint16_t,std::string>{1824u,"090103"},
	std::pair<uint16_t,std::string>{1825u,"100103"},
	std::pair<uint16_t,std::string>{1826u,"110103"},
	std::pair<uint16_t,std::string>{1827u,"120103"},
	std::pair<uint16_t,std::string>{1828u,"130103"},
	std::pair<uint16_t,std::string>{1829u,"140103"},
	std::pair<uint16_t,std::string>{1400u,"000000"},
	std::pair<uint16_t,std::string>{1401u,"000100"},
	std::pair<uint16_t,std::string>{1402u,"000200"},
	std::pair<uint16_t,std::string>{1403u,"000300"},
	std::pair<uint16_t,std::string>{1404u,"000400"},
	std::pair<uint16_t,std::string>{1405u,"000500"},
	std::pair<uint16_t,std::string>{1406u,"000600"},
	std::pair<uint16_t,std::string>{1407u,"000700"},
	std::pair<uint16_t,std::string>{1408u,"000800"},
	std::pair<uint16_t,std::string>{1409u,"000900"},
	std::pair<uint16_t,std::string>{1410u,"001000"},
	std::pair<uint16_t,std::string>{1411u,"001100"},
	std::pair<uint16_t,std::string>{1412u,"001200"},
	std::pair<uint16_t,std::string>{1416u,"010000"},
	std::pair<uint16_t,std::string>{1417u,"010100"},
	std::pair<uint16_t,std::string>{1418u,"010200"},
	std::pair<uint16_t,std::string>{1419u,"010300"},
	std::pair<uint16_t,std::string>{1420u,"010400"},
	std::pair<uint16_t,std::string>{1421u,"010500"},
	std::pair<uint16_t,std::string>{1422u,"010600"},
	std::pair<uint16_t,std::string>{1423u,"010700"},
	std::pair<uint16_t,std::string>{1424u,"010800"},
	std::pair<uint16_t,std::string>{1425u,"010900"},
	std::pair<uint16_t,std::string>{1426u,"011000"},
	std::pair<uint16_t,std::string>{1427u,"011100"},
	std::pair<uint16_t,std::string>{1428u,"011200"}
	};

//Original v1 map: We don't see any veto coincidences between the two layers with this!
std::map<uint16_t,std::string> LoadCCData::slotchantopmtidv1{
        std::pair<uint16_t,std::string>{1701u,"000002"},
        std::pair<uint16_t,std::string>{1702u,"000102"},
        std::pair<uint16_t,std::string>{1703u,"000202"},
        std::pair<uint16_t,std::string>{1704u,"000302"},
        std::pair<uint16_t,std::string>{1705u,"000402"},
        std::pair<uint16_t,std::string>{1706u,"000502"},
        std::pair<uint16_t,std::string>{1707u,"000602"},
        std::pair<uint16_t,std::string>{1708u,"000702"},
        std::pair<uint16_t,std::string>{1709u,"000802"},
        std::pair<uint16_t,std::string>{1710u,"000902"},
        std::pair<uint16_t,std::string>{1711u,"001002"},
        std::pair<uint16_t,std::string>{1712u,"001102"},
        std::pair<uint16_t,std::string>{1713u,"001202"},
        std::pair<uint16_t,std::string>{1714u,"010002"},
        std::pair<uint16_t,std::string>{1715u,"010102"},
        std::pair<uint16_t,std::string>{1716u,"010202"},
        std::pair<uint16_t,std::string>{1717u,"010302"},
        std::pair<uint16_t,std::string>{1718u,"010402"},
        std::pair<uint16_t,std::string>{1719u,"010502"},
        std::pair<uint16_t,std::string>{1720u,"010602"},
        std::pair<uint16_t,std::string>{1721u,"010702"},
        std::pair<uint16_t,std::string>{1722u,"010802"},
        std::pair<uint16_t,std::string>{1723u,"010902"},
        std::pair<uint16_t,std::string>{1724u,"011002"},
        std::pair<uint16_t,std::string>{1725u,"011102"},
        std::pair<uint16_t,std::string>{1726u,"011202"},
        std::pair<uint16_t,std::string>{1801u,"000003"},
        std::pair<uint16_t,std::string>{1802u,"010003"},
        std::pair<uint16_t,std::string>{1803u,"020003"},
        std::pair<uint16_t,std::string>{1804u,"030003"},
        std::pair<uint16_t,std::string>{1805u,"040003"},
        std::pair<uint16_t,std::string>{1806u,"050003"},
        std::pair<uint16_t,std::string>{1807u,"060003"},
        std::pair<uint16_t,std::string>{1808u,"070003"},
        std::pair<uint16_t,std::string>{1809u,"080003"},
        std::pair<uint16_t,std::string>{1810u,"090003"},
        std::pair<uint16_t,std::string>{1811u,"100003"},
        std::pair<uint16_t,std::string>{1812u,"110003"},
        std::pair<uint16_t,std::string>{1813u,"120003"},
        std::pair<uint16_t,std::string>{1814u,"130003"},
        std::pair<uint16_t,std::string>{1815u,"000103"},
        std::pair<uint16_t,std::string>{1816u,"010103"},
        std::pair<uint16_t,std::string>{1817u,"020103"},
        std::pair<uint16_t,std::string>{1818u,"030103"},
        std::pair<uint16_t,std::string>{1819u,"040103"},
        std::pair<uint16_t,std::string>{1820u,"050103"},
        std::pair<uint16_t,std::string>{1821u,"060103"},
        std::pair<uint16_t,std::string>{1822u,"070103"},
        std::pair<uint16_t,std::string>{1823u,"080103"},
        std::pair<uint16_t,std::string>{1824u,"090103"},
        std::pair<uint16_t,std::string>{1825u,"100103"},
        std::pair<uint16_t,std::string>{1826u,"110103"},
        std::pair<uint16_t,std::string>{1827u,"120103"},
        std::pair<uint16_t,std::string>{1828u,"130103"},
        std::pair<uint16_t,std::string>{1829u,"140103"},
	std::pair<uint16_t,std::string>{1401u,"000000"},
	std::pair<uint16_t,std::string>{1402u,"000100"},
	std::pair<uint16_t,std::string>{1403u,"000200"},
	std::pair<uint16_t,std::string>{1404u,"000300"},
	std::pair<uint16_t,std::string>{1405u,"000400"},
	std::pair<uint16_t,std::string>{1406u,"000500"},
	std::pair<uint16_t,std::string>{1407u,"000600"},
	std::pair<uint16_t,std::string>{1408u,"000700"},
	std::pair<uint16_t,std::string>{1409u,"000800"},
	std::pair<uint16_t,std::string>{1410u,"000900"},
	std::pair<uint16_t,std::string>{1411u,"001000"},
	std::pair<uint16_t,std::string>{1412u,"001100"},
	std::pair<uint16_t,std::string>{1413u,"001200"},
	std::pair<uint16_t,std::string>{1414u,"010000"},
	std::pair<uint16_t,std::string>{1415u,"010100"},
	std::pair<uint16_t,std::string>{1416u,"010200"},
	std::pair<uint16_t,std::string>{1417u,"010300"},
	std::pair<uint16_t,std::string>{1418u,"010400"},
	std::pair<uint16_t,std::string>{1419u,"010500"},
	std::pair<uint16_t,std::string>{1420u,"010600"},
	std::pair<uint16_t,std::string>{1421u,"010700"},
	std::pair<uint16_t,std::string>{1422u,"010800"},
	std::pair<uint16_t,std::string>{1423u,"010900"},
	std::pair<uint16_t,std::string>{1424u,"011000"},
	std::pair<uint16_t,std::string>{1425u,"011100"},
	std::pair<uint16_t,std::string>{1426u,"011200"}
	};

// some additional channels frequently have hits,
// even though apparently nothing was connected:
//	18-30
//	14-00
//	14-27
//	14-28
//	17-00
//	17-27
//	17-28
// TODO maybe add these channels, or suppress their messages...
// otherwise it just results in a lot of warnings
// Update 10th of July, 2020
// --> Added 14-00, 14-27, 14-28 in channelmap version 2
// Not sure about 17-00, 17-27, 17-28 (MRD channels might have a similar issue as veto channels!)


/////////////
// MRDTesting TDC setup, part 1:
//	std::map<uint16_t,std::string> slotchantopmtid{
//		std::pair<uint16_t,std::string>{1423u,"000002"},
//		std::pair<uint16_t,std::string>{1416u,"000102"},
//		std::pair<uint16_t,std::string>{1417u,"000202"},
//		std::pair<uint16_t,std::string>{1418u,"000302"},
//		std::pair<uint16_t,std::string>{1419u,"000402"},
//		std::pair<uint16_t,std::string>{1420u,"000502"},
//		std::pair<uint16_t,std::string>{1421u,"000602"},
//		std::pair<uint16_t,std::string>{1422u,"000702"},
//		std::pair<uint16_t,std::string>{1424u,"000802"},
//		std::pair<uint16_t,std::string>{1425u,"000902"},
//		std::pair<uint16_t,std::string>{1426u,"001002"},
//		std::pair<uint16_t,std::string>{1427u,"001102"},
//		std::pair<uint16_t,std::string>{1328u,"001202"},
//		std::pair<uint16_t,std::string>{1823u,"000004"},
//		std::pair<uint16_t,std::string>{1816u,"000104"},
//		std::pair<uint16_t,std::string>{1817u,"000204"},
//		std::pair<uint16_t,std::string>{1818u,"000304"},
//		std::pair<uint16_t,std::string>{1819u,"000404"},
//		std::pair<uint16_t,std::string>{1820u,"000504"},
//		std::pair<uint16_t,std::string>{1821u,"000604"},
//		std::pair<uint16_t,std::string>{1822u,"000704"},
//		std::pair<uint16_t,std::string>{1824u,"000804"},
//		std::pair<uint16_t,std::string>{1825u,"000904"},
//		std::pair<uint16_t,std::string>{1826u,"001004"},
//		std::pair<uint16_t,std::string>{1827u,"001104"},
//		std::pair<uint16_t,std::string>{1828u,"001204"},
//		std::pair<uint16_t,std::string>{1800u,"000006"},
//		std::pair<uint16_t,std::string>{1801u,"000106"},
//		std::pair<uint16_t,std::string>{1802u,"000206"},
//		std::pair<uint16_t,std::string>{1803u,"000306"},
//		std::pair<uint16_t,std::string>{1804u,"000406"},
//		std::pair<uint16_t,std::string>{1805u,"000506"},
//		std::pair<uint16_t,std::string>{1806u,"000606"},
//		std::pair<uint16_t,std::string>{1807u,"000706"},
//		std::pair<uint16_t,std::string>{1808u,"000806"},
//		std::pair<uint16_t,std::string>{1809u,"000906"},
//		std::pair<uint16_t,std::string>{1810u,"001006"},
//		std::pair<uint16_t,std::string>{1811u,"001106"},
//		std::pair<uint16_t,std::string>{1812u,"001206"},
//		std::pair<uint16_t,std::string>{1407u,"000008"},
//		std::pair<uint16_t,std::string>{1400u,"000108"},
//		std::pair<uint16_t,std::string>{1401u,"000208"},
//		std::pair<uint16_t,std::string>{1402u,"000308"},
//		std::pair<uint16_t,std::string>{1403u,"000408"},
//		std::pair<uint16_t,std::string>{1404u,"000508"},
//		std::pair<uint16_t,std::string>{1405u,"000608"},
//		std::pair<uint16_t,std::string>{1406u,"000708"},
//		std::pair<uint16_t,std::string>{1408u,"000808"},
//		std::pair<uint16_t,std::string>{1409u,"000908"},
//		std::pair<uint16_t,std::string>{1410u,"001008"},
//		std::pair<uint16_t,std::string>{1411u,"001108"},
//		std::pair<uint16_t,std::string>{1412u,"001208"},
//		std::pair<uint16_t,std::string>{1700u,"000010"},
//		std::pair<uint16_t,std::string>{1701u,"000110"},
//		std::pair<uint16_t,std::string>{1702u,"000210"},
//		std::pair<uint16_t,std::string>{1703u,"000310"},
//		std::pair<uint16_t,std::string>{1704u,"000410"},
//		std::pair<uint16_t,std::string>{1705u,"000510"},
//		std::pair<uint16_t,std::string>{1706u,"000610"},
//		std::pair<uint16_t,std::string>{1707u,"000710"},
//		std::pair<uint16_t,std::string>{1708u,"000810"},
//		std::pair<uint16_t,std::string>{1709u,"000910"},
//		std::pair<uint16_t,std::string>{1710u,"001010"},
//		std::pair<uint16_t,std::string>{1711u,"001110"},
//		std::pair<uint16_t,std::string>{1712u,"001210"},
//		std::pair<uint16_t,std::string>{1716u,"000012"},
//		std::pair<uint16_t,std::string>{1717u,"000112"},
//		std::pair<uint16_t,std::string>{1718u,"000212"},
//		std::pair<uint16_t,std::string>{1719u,"000312"},
//		std::pair<uint16_t,std::string>{1720u,"000412"},
//		std::pair<uint16_t,std::string>{1721u,"000512"},
//		std::pair<uint16_t,std::string>{1722u,"000612"},
//		std::pair<uint16_t,std::string>{1723u,"000712"},
//		std::pair<uint16_t,std::string>{1724u,"000812"},
//		std::pair<uint16_t,std::string>{1725u,"000912"},
//		std::pair<uint16_t,std::string>{1726u,"001012"},
//		std::pair<uint16_t,std::string>{1727u,"001112"},
//		std::pair<uint16_t,std::string>{1728u,"001212"}
//	};


//		MRDData has the following members:
//		UInt_t                     Trigger;         // TDC readout counter
//		UInt_t                     OutNumber;       // number of hits in this event - size of vectors
//		ULong64_t                  TimeStamp;       // UTC MILLISECONDS since unix epoch
//		std::vector<std::string>*  Type;            // card type - all cards are "TDC"
//		std::vector<unsigned int>* Value;           // time in TDC ticks from TimeStamp
//		std::vector<unsigned int>* Slot;            // slot of hit
//		std::vector<unsigned int>* Channel;         // channel of hit
//		
//		The MRD process is:
//			1. Trigger card sends common start to MRD cards
//			2. A timer is started on all channels.
//			3. When a channel receives a pulse, the timer stops. XXX: only the first pulse is recorded!XXX
//			4. After all channels either record hits or time out (currently 4.2us) everything is read out.
//			   A timestamp is created at time of readout. Note this is CLOSE TO, BUT NOT EQUAL TO the trigger
//			   time. (Probably ~ trigger time + timeout...)
//			5. Channels with a hit will have an entry created with 'Value' = clock ticks between the common
//			   start and when the hit arrived. Channels that timed out have no entry.
//		
//		Timestamp is a UTC [MILLISECONDS] timestamp of when the readout ended.
//		To correctly map Value to an actual time one would need to match the MRD timestamp
//		to the trigger card timestamp (which will be more accurate)
//		then add Value * MRD_NS_PER_SAMPLE to the Trigger time.
//		
//		Only the first pulse will be recorded .... IS THIS OK? XXX
//		What about pre-trigger? - common start issued by trigger is from Beam not on NDigits,
//		so always pre-beam...?
//		TDC records with resolution 4ns = 1 sample per ns? or round times to nearest/round down to 4ns?
//		TDCRes https://github.com/ANNIEDAQ/ANNIEDAQ/blob/master/configfiles/TDCreg
