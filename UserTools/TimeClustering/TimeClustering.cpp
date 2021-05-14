/* vim:set noexpandtab tabstop=4 wrap */
#include "TimeClustering.h"

// for sleeping
#include <thread>  // std::this_thread::sleep_for
#include <chrono>  // std::chrono::seconds

// for ROOT debug plot
#include "TROOT.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"

TimeClustering::TimeClustering():Tool(){}


bool TimeClustering::Initialise(std::string configfile, DataModel &data){
	
	if(verbosity) cout<<"Initializing tool TimeClustering"<<endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	LaunchTApplication = false;
	MakeSingleEventPlots = false;      //very verbose, mostly for debugging
	ModifiedTDCData = false;	
	//Set default values for time-shifted TDC channels
	TimeShiftChannels = "";
	shifted_channels = {std::make_pair(0,51),std::make_pair(82,107),std::make_pair(142,167),std::make_pair(194,219),std::make_pair(250,275),std::make_pair(306,332)};

	m_variables.Get("verbosity",verbosity);
	m_variables.Get("IsData",isData);
	m_variables.Get("MinDigitsForTrack",minimumdigits);
	m_variables.Get("MaxMrdSubEventDuration",maxsubeventduration);
	m_variables.Get("MinSubeventTimeSep",minimum_subevent_timeseparation);
	m_variables.Get("MakeMrdDigitTimePlot",MakeMrdDigitTimePlot);  /// XXX XXX remame
	m_variables.Get("MakeSingleEventPlots",MakeSingleEventPlots);
	m_variables.Get("LaunchTApplication",LaunchTApplication);
	m_variables.Get("OutputROOTFile",output_rootfile);
	m_variables.Get("MapChankey_WCSimID",file_chankeymap);
	m_variables.Get("ModifiedTDCData",ModifiedTDCData);
        m_variables.Get("TimeShiftChannels",TimeShiftChannels);

	if (!MakeMrdDigitTimePlot) LaunchTApplication = false;  //no use launching TApplication when histograms are not produced
	
	if(LaunchTApplication){
		canvwidth = 900;
		canvheight = 600;
		
		// Since we want to make event-wise debug plots, we need a TApplication to show them
		// ======================
		// There may only be one TApplication, so if another tool has already made one
		// register ourself as a user. Otherwise, make one and put a pointer in the CStore for other Tools
		// create the ROOT application to show histograms
		
		int myargc=0;
		//char *myargv[] = {(const char*)"Ahh shark!"};
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			Log("TimeClustering Tool: Making global TApplication",v_error,verbosity);
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			Log("TimeClustering Tool: Retrieving global TApplication",v_error,verbosity);
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}
	
	// setup objects to save clustered hit times to a *.root file
	if (MakeMrdDigitTimePlot){
		std::stringstream ss_rootfilename;
		ss_rootfilename << output_rootfile << ".root";
		Log("TimeClustering tool: Create ROOT-file "+ss_rootfilename.str(),v_message,verbosity);
		//save the cluster time information for different trigger types in different histograms
		mrddigitts_file = new TFile(ss_rootfilename.str().c_str(),"RECREATE");
		mrddigitts_cosmic_cluster = new TH1D("mrddigitts_cosmic_cluster","MRD Cosmic Times Clustered",1000,0,4000);
		mrddigitts_beam_cluster = new TH1D("mrddigitts_beam_cluster","MRD Beam Times Clustered",1000,0,4000);
		mrddigitts_noloopback_cluster = new TH1D("mrddigitts_noloopback_cluster","MRD No Loopback Times Clustered",1000,0,4000);
		mrddigitts_cluster = new TH1D("mrddigitts_cluster","MRD Times Clustered",1000,0,4000);
		mrddigitts_cosmic = new TH1D("mrddigitts_cosmic","MRD Cosmic Times",1000,0,4000);
		mrddigitts_beam = new TH1D("mrddigitts_beam","MRD Beam Times",1000,0,4000);
		mrddigitts_noloopback = new TH1D("mrddigitts_noloopback","MRD No Loopback Times",1000,0,4000);
		mrddigitts = new TH1D("mrddigitts","MRD Times",1000,0,4000);
		mrddigitts_vertical = new TH1D("mrddigitts_vertical","MRD Times (Vertical Layers)",1000,0,4000);
		mrddigitts_horizontal = new TH1D("mrddigitts_horizontal","MRD Times (Horizontal Layers)",1000,0,4000); 
		hist_chankey = new TH1D("hist_chankey","Chankey frequency",340,0,340);
		hist_chankey_cluster = new TH1D("hist_chankey_cluster","Chankey frequency (Cluster)",340,0,340);
		hist_chankey_time = new TH2D("hist_chankey_time","Chankey vs. time",340,0,340,1000,0,4000);
		hist_chankey_time_cluster = new TH2D("hist_chankey_time_cluster","Chankey vs. time (cluster)",340,0,340,1000,0,4000);
		hist_chankey_multi = new TH1D("hist_chankey_multi","Chankey multi-hit frequency",340,0,340);

		if (MakeSingleEventPlots){
			mrddigitts_single = new TH1D("mrddigitts_single","MRD Single Times",1000,0,4000);
			mrddigitts_cluster_single = new TH1D("mrddigitts_cluster_single","MRD Single Times",1000,0,4000);
		}
		gROOT->cd();
	}
	
	// Setup mapping from Channelkeys to WCSim IDs (important for track fitting with old MRD classes in FindMrdTracks)
	if (isData){
		ifstream file_mapping(file_chankeymap);
		unsigned long temp_chankey;
		int temp_wcsimid;
		while (!file_mapping.eof()){
			file_mapping>>temp_chankey>>temp_wcsimid;
			if (file_mapping.eof()) break;
			channelkey_to_mrdpmtid.emplace(temp_chankey,temp_wcsimid);
			mrdpmtid_to_channelkey.emplace(temp_wcsimid,temp_chankey);
			Log("TimeClustering tool: Emplaced temp_chankey "+std::to_string(temp_chankey)+" with temp_wcsimid "+std::to_string(temp_wcsimid)+"into channelkey_to_mrdpmtid object!",v_debug,verbosity);
		}
		file_mapping.close();
		m_data->CStore.Set("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
		m_data->CStore.Set("mrd_tubeid_to_channelkey",mrdpmtid_to_channelkey);
		}
	else {
		get_ok = m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);  //for MC, simply get the sample obtained from the LoadWCSim tool
		if(not get_ok){
			Log("TimeClustering Tool: Error! No channelkey_to_mrdpmtid in CStore!",v_error,verbosity);
			return false;
		}
		get_ok = m_data->CStore.Get("channelkey_to_faccpmtid",channelkey_to_faccpmtid);
		if(not get_ok){
			Log("TimeClustering Tool: Error! No channelkey_to_faccpmtid in CStore!",v_error,verbosity);
			return false;
		}
	}

	//Read in time shifted channels from configuration file
	if (TimeShiftChannels != ""){
		shifted_channels.clear();		
		ifstream file_timeshift(TimeShiftChannels);
		if (file_timeshift.good()){
			unsigned long lower_chkey, upper_chkey;
			while (!file_timeshift.eof()){
				file_timeshift >> lower_chkey >> upper_chkey;
				if (file_timeshift.eof()) break;
				std::cout <<"lower_chkey: "<<lower_chkey<<", upper_chkey: "<<upper_chkey<<std::endl;
				if (lower_chkey > upper_chkey){
					Log("TimeClustering tool: Error when reading in TiemShiftChannels file! Upper chankey "+std::to_string(upper_chkey)+" is lower than lower chankey! "+std::to_string(lower_chkey)+" Abort",v_error,verbosity);
					return false;
				}
				shifted_channels.push_back(std::make_pair(lower_chkey,upper_chkey));
			}
		} else {
			Log("TimeClustering Tool: Error! Timeshift file "+TimeShiftChannels+" does not exist!",v_error,verbosity);
			return false;
		}
	}
	
	// Get Detectors map to divide in horizontal and vertical layers
	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
	
	return true;
}


bool TimeClustering::Execute(){
	
	Log("Tool TimeClustering: finding time clusters in next event",v_message,verbosity);
	
	MrdTimeClusters.clear(); // do not carry over any subevents
	m_data->CStore.Set("MrdTimeClusters",MrdTimeClusters);
	m_data->CStore.Set("NumMrdTimeClusters",0);
	int mrdeventcounter=0;
	
	//Get Trigger type to decide which events are cosmic/beam/etc
	std::string MRDTriggertype;
	m_data->Stores.at("ANNIEEvent")->Get("MRDTriggerType",MRDTriggertype);
	Log("TimeClustering tool: MRDTriggertype is "+MRDTriggertype+" (from ANNIEEvent store)",v_debug,verbosity);
	m_data->Stores.at("ANNIEEvent")->Get("EventNumber",evnum);
	
	// extract the digits from the annieevent and put them into separate vectors used by the track finder
	mrddigitpmtsthisevent.clear();
	mrddigitchankeysthisevent.clear();
	mrddigittimesthisevent.clear();
	mrddigitchargesthisevent.clear();
	
        std::vector<unsigned long> multi_vector;

	if (MakeMrdDigitTimePlot && MakeSingleEventPlots){
		mrddigitts_single->Reset();
		mrddigitts_cluster_single->Reset();
		std::stringstream ss_single, ss_cluster_single;
		ss_single << "mrddigitts_single_" << evnum;
		ss_cluster_single << "mrddigitts_cluster_single_" << evnum;
		mrddigitts_single->SetName(ss_single.str().c_str());
		mrddigitts_cluster_single->SetName(ss_cluster_single.str().c_str());
	}
	
	if (isData){
		std::cout <<"TimeClustering tool: Data file: Getting TDCData object"<<std::endl;
		get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);  // a std::map<unsigned long,vector<Hit>>

		if(not get_ok){
			Log("TimeClustering Tool: No TDC data in ANNIEEvent!",v_error,verbosity);
			return true;
		}
	} else {
		std::cout <<"TimeClustering tool: MC file: Getting TDCData object"<<std::endl;
		if (ModifiedTDCData) get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData_mod",TDCData_MC);  // a std::map<unsigned long,vector<MCHit>>, with artificially applied efficiencies
		else get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData_MC);  // a std::map<unsigned long,vector<MCHit>>, without artifically applied efficiencies
		if(not get_ok){
			Log("TimeClustering Tool: No TDC data in ANNIEEvent!",v_error,verbosity);
			return true;
		}
	}
	
	// check if we have any hits to process
	if (isData){
		if(TDCData->size()==0 || !TDCData){
			Log("TimeClustering tool: No TDC hits to find clusters in / TDCData object does not exist.",v_warning,verbosity);
			return true; // XXX XXX XXX XXX XXX XXX XXX
		}
	} else {
		if(TDCData_MC->size()==0 || !TDCData_MC){
			Log("TimeClustering tool: No TDC hits to find clusters in / TDCData object does not exist.",v_warning,verbosity);
			return true; // XXX XXX XXX XXX XXX XXX XXX
		}
	}
	
	if (isData) Log("TimeClustering Tool: Retrieving digit info from "+to_string(TDCData->size())+" hit pmts",v_debug,verbosity);
	else Log("TimeClustering Tool: Retrieving digit info from "+to_string(TDCData_MC->size())+" hit pmts",v_debug,verbosity);
	
	// just dump all the hit times in this event into a vector. Allows us to sort hit times and search for clusters.
	if (isData){
		for(auto&& anmrdpmt : (*TDCData)){
			unsigned long chankey = anmrdpmt.first;
//			Detector* thedetector = geom->ChannelToDetector(chankey);
//			if(thedetector==nullptr){
//				Log("TimeClustering Tool: Null detector in TDCData!",v_error,verbosity);
//				continue;
//			}
//			if(thedetector->GetDetectorElement()!="MRD") continue; // this is a veto hit, not an MRD hit. XXX keep this cut?
			for(auto&& hitsonthismrdpmt : anmrdpmt.second){
				if (channelkey_to_mrdpmtid.find(chankey) != channelkey_to_mrdpmtid.end()){
					mrddigitpmtsthisevent.push_back(channelkey_to_mrdpmtid[chankey]);
					mrddigitchankeysthisevent.push_back(chankey);
					double time = hitsonthismrdpmt.GetTime();
					//Times of channelkeys in TDC crate 7 (vertical channels) are systematically late by ~20ns with respect to channels in crate 8 --> shift the itmes for those channelkeys manually by 20ns
					//Affected channelkeys are stored in shifted_channels and can be configured in the config file
					for (int i_shift=0; i_shift < (int) shifted_channels.size(); i_shift++){
						std::pair<unsigned long,unsigned long> temp_pair = shifted_channels.at(i_shift);
						if (chankey >= temp_pair.first && chankey <= temp_pair.second) time -= 20.;
					}
					mrddigittimesthisevent.push_back(time);
					mrddigitchargesthisevent.push_back(hitsonthismrdpmt.GetCharge());
					if(MakeMrdDigitTimePlot){  // XXX XXX XXX rename
						// fill the histogram if we're checking
						if (MakeSingleEventPlots) mrddigitts_single->Fill(hitsonthismrdpmt.GetTime());
						mrddigitts->Fill(hitsonthismrdpmt.GetTime());
						if (MRDTriggertype == "Cosmic") mrddigitts_cosmic->Fill(hitsonthismrdpmt.GetTime());
						
						else if (MRDTriggertype == "Beam") mrddigitts_beam->Fill(hitsonthismrdpmt.GetTime());
						else if (MRDTriggertype == "No Loopback") mrddigitts_noloopback->Fill(hitsonthismrdpmt.GetTime());    //this triggertype should not occur if everything is running smoothly, but it can serve as a good cross-check in any case
						Detector* thistube = geom->ChannelToDetector(chankey);
						unsigned long detkey = thistube->GetDetectorID();
						hist_chankey->Fill(detkey);
						hist_chankey_time->Fill(detkey,hitsonthismrdpmt.GetTime());
						if (std::count(mrddigitchankeysthisevent.begin(),mrddigitchankeysthisevent.end(),chankey)>1) hist_chankey_multi->Fill(chankey);
						Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
						int orientation = mrdpaddle->GetOrientation(); // 0 is horizontal, 1 is vertical
						if (orientation == 0) mrddigitts_horizontal->Fill(hitsonthismrdpmt.GetTime());
						else mrddigitts_vertical->Fill(hitsonthismrdpmt.GetTime());
					}
				} else {
					Log("TimeClustering tool: Did not find channelkey "+std::to_string(chankey)+" in chankey_to_mrdpmtid map.",v_warning,verbosity);
				}
			}
		}
	} else {
		for(auto&& anmrdpmt : (*TDCData_MC)){
			unsigned long chankey = anmrdpmt.first;
//			// sanity checks
//			Detector* thedetector = geom->ChannelToDetector(chankey);
//			if(thedetector==nullptr){
//				Log("TimeClustering Tool: Null detector in TDCData_MC!",v_error,verbosity);
//				continue;
//			}
//			if((thedetector->GetDetectorElement()!="MRD") &&
//			   (thedetector->GetDetectorElement()!="Veto") ){
//				Log("TDC hit on channel that is neither MRD nor Veto! Type: " + thedetector->GetDetectorElement(),v_debug,verbosity);
//			}
//			if( (channelkey_to_mrdpmtid.count(chankey) ==0) &&
//				(channelkey_to_faccpmtid.count(chankey)==0)){
//				Log("TimeClustering Tool: MRD PMT with ID not in channelkey_to_mrdpmtid map!",v_error,verbosity);
//				if(true||verbosity>2){
//					std::cerr<<"We have: "<<channelkey_to_mrdpmtid.size()
//							 <<" known MRD mappings and they are: {"<<std::endl;
//					for(auto&& apair : channelkey_to_mrdpmtid){
//						std::cout<<apair.first<<" : "<<apair.second<<std::endl;
//					}
//					std::cout<<"}"<<std::endl;
//					std::cerr<<"We have: "<<channelkey_to_faccpmtid.size()
//							 <<" known FACC mappings and they are: {"<<std::endl;
//					for(auto&& apair : channelkey_to_faccpmtid){
//						std::cout<<apair.first<<" : "<<apair.second<<std::endl;
//					}
//					std::cout<<"}"<<std::endl;
//				}
//				continue;
//			}
//			int wcsimtubeid = -1;
//			wcsimtubeid = (channelkey_to_mrdpmtid.count(chankey)) ? channelkey_to_mrdpmtid.at(chankey) :
//																	channelkey_to_faccpmtid.at(chankey);
//			if(wcsimtubeid==0){
//				Log("TimeClustering Tool: channel with wcsimpmtid 0! IDs should number from 1!",v_error, verbosity);
//				continue;
//			}
			// checking channelkey_to_mrdpmtid (as opposed to channelkey_to_faccpmtid)
			// will filter out MRD PMTs only
			int pmtidwcsim=-1;
			if (channelkey_to_mrdpmtid.count(chankey)){
				pmtidwcsim = channelkey_to_mrdpmtid.at(chankey)-1;
			} /*else if(channelkey_to_faccpmtid.count(chankey)){ /////omit facc hits for MRD clusters
				pmtidwcsim = channelkey_to_faccpmtid.at(chankey)-1;
			}*/
			for(auto&& hitsonthismrdpmt : anmrdpmt.second){
				if(pmtidwcsim>=0){
					mrddigitpmtsthisevent.push_back(pmtidwcsim);
					mrddigitchankeysthisevent.push_back(chankey);
					mrddigittimesthisevent.push_back(hitsonthismrdpmt.GetTime());
					mrddigitchargesthisevent.push_back(hitsonthismrdpmt.GetCharge());
					mrddigitchankeysthisevent.push_back(chankey);
					if(MakeMrdDigitTimePlot){  // XXX XXX XXX rename
						// fill the histogram if we're checking
						if (MakeSingleEventPlots) mrddigitts_single->Fill(hitsonthismrdpmt.GetTime());
						mrddigitts->Fill(hitsonthismrdpmt.GetTime());
						if (MRDTriggertype == "Cosmic") mrddigitts_cosmic->Fill(hitsonthismrdpmt.GetTime());
						
						else if (MRDTriggertype == "Beam") mrddigitts_beam->Fill(hitsonthismrdpmt.GetTime());
						else if (MRDTriggertype == "No Loopback") mrddigitts_noloopback->Fill(hitsonthismrdpmt.GetTime());   //this triggertype should not occur if everything is running smoothly, but it can serve as a good cross-check in any case
						Detector* thistube = geom->ChannelToDetector(chankey);
						unsigned long detkey = thistube->GetDetectorID();
						Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
						int orientation = mrdpaddle->GetOrientation(); // 0 is horizontal, 1 is vertical
						if (orientation == 0) mrddigitts_horizontal->Fill(hitsonthismrdpmt.GetTime());
						else mrddigitts_vertical->Fill(hitsonthismrdpmt.GetTime());
					}
				} else {
					Log("TimeClustering tool: Did not find channelkey "+std::to_string(chankey)+" in chankey_to_mrdpmtid or channelkey_to_faccpmtid maps.",v_warning,verbosity);
				}
			}
		}
	}
	int numdigits = mrddigittimesthisevent.size();
	Log("TimeClustering Tool: Searching "+to_string(numdigits)+" (MRD) hits this event for clusters",v_debug,verbosity);
	
	///////////////////////////
	// now do the time clustering
	
	// First the easy check; do we have enough hits to make a cluster
	if(numdigits<minimumdigits){
		// =====================================================================================
		// NOT ENOUGH DIGITS IN THIS EVENT
		// ======================================
		Log("TimeClustering Tool: Insufficient digits in this event to make any clusters; returning",v_debug,verbosity);
		return true;
		// ======================================================================================
	}
	
	// MEASURE TOTAL EVENT DURATION
	// ============================
	double eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	double eventstarttime = *std::min_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
	double eventduration = (eventendtime - eventstarttime);
	Log("TimeClustering Tool: mrd event start: "+to_string(eventstarttime)
			+", end : "+to_string(eventendtime)+", duration : "+to_string(eventduration),v_debug,verbosity);
	
	// If all hits within maxsubeventduration[ns], just one subevent.
	if((eventduration<maxsubeventduration)&&(numdigits>=minimumdigits)){
	// JUST ONE SUBEVENT
	// =================
		Log("TimeClustering Tool: All hits this event within one subevent.",v_debug,verbosity);
		std::vector<int> digitidsinasubevent(numdigits);    // a vector of indices of the digits in this subevent
		std::iota(digitidsinasubevent.begin(),digitidsinasubevent.end(),0);  // fill with 1-N, as all digits are are in this subevent
		MrdTimeClusters.push_back(digitidsinasubevent);
		if (MakeMrdDigitTimePlot){
			for (unsigned int i_time=0; i_time<mrddigittimesthisevent.size(); i_time++){
				mrddigitts_file->cd();
				hist_chankey_cluster->Fill(mrddigitchankeysthisevent.at(i_time));
				hist_chankey_time_cluster->Fill(mrddigitchankeysthisevent.at(i_time),mrddigittimesthisevent.at(i_time));
				if (MakeSingleEventPlots) mrddigitts_cluster_single->Fill(mrddigittimesthisevent.at(i_time));
				mrddigitts_cluster->Fill(mrddigittimesthisevent.at(i_time));
				if (MRDTriggertype == "Cosmic") mrddigitts_cosmic_cluster->Fill(mrddigittimesthisevent.at(i_time));
				else if (MRDTriggertype == "Beam") mrddigitts_beam_cluster->Fill(mrddigittimesthisevent.at(i_time));
				else if (MRDTriggertype == "No Loopback") mrddigitts_noloopback_cluster->Fill(mrddigittimesthisevent.at(i_time));
			}
		}
		mrdeventcounter++;
		
	} else {
	// MORE THAN ONE SUBEVENT
	// ======================
		// COUNT SUBEVENTS
		// ---------------
		// this event has multiple subevents. We need to split hits into which subevent they belong to.
		// first scan over the times and look for gaps where no digits lie, using these to delimit 'subevents'
		std::vector<float> subeventhittimesv;   // a vector of the starting times of a given subevent
		std::vector<float> subeventendtimesv;
		std::vector<double> sorteddigittimes(mrddigittimesthisevent);
		std::sort(sorteddigittimes.begin(), sorteddigittimes.end());
		subeventhittimesv.push_back(sorteddigittimes.at(0));
		for(unsigned int i=0;i<sorteddigittimes.size()-1;i++){
			float timetonextdigit = sorteddigittimes.at(i+1)-sorteddigittimes.at(i);
			if(timetonextdigit>minimum_subevent_timeseparation){
				subeventendtimesv.push_back(sorteddigittimes.at(i));
				subeventhittimesv.push_back(sorteddigittimes.at(i+1));
				Log("TimeClustering Tool: Setting subevent time threshold at "+to_string(subeventhittimesv.back()),v_debug,verbosity);
			}
		}
		Log("TimeClustering Tool: Found "+to_string(subeventhittimesv.size())+" subevents this event",v_debug,verbosity);
		
		//write subeventhittimesv to CStore for subsequent tools (e.g. FindMrdTracks)
		m_data->CStore.Set("ClusterStartTimes",subeventhittimesv);
		
		// DEBUG CHECK
		// -----------
		// debug check of the timing splitting: draw a histogram of the times
		if(LaunchTApplication){
			// remake the canvas: close it when done viewing (deleted by ROOT)
			if(gROOT->FindObject("timeClusterCanvas")) delete timeClusterCanvas;
			timeClusterCanvas = new TCanvas("timeClusterCanvas","timeClusterCanvas",canvwidth,canvheight);
			timeClusterCanvas->SetWindowSize(canvwidth,canvheight);
			timeClusterCanvas->cd();
			if (MakeSingleEventPlots) mrddigitts_single->Draw();  //the other histograms are simply written to file, otherwise swamped with unnecessary information during execution
			timeClusterCanvas->Modified();
			timeClusterCanvas->Update();
			gSystem->ProcessEvents();
			//rootTApp->Run();
			//std::this_thread::sleep_for (std::chrono::seconds(5));
			while(gROOT->FindObject("timeClusterCanvas")!=nullptr){
				gSystem->ProcessEvents();
				std::this_thread::sleep_for (std::chrono::milliseconds(500));
			}
		}
		
		// SORT HITS INTO SUBEVENTS
		// ------------------------
		// now we need to sort the digits into the subevents they belong to:
		// temporary storage vectors to hold the digits in the present subevent
		std::vector<int> digitidsinasubevent;
		std::vector<int> tubeidsinasubevent;
		std::vector<double> digittimesinasubevent;
		std::vector<double> digitchargesinasubevent;
		// a vector to record the subevent number for each hit, to know if we've allocated it yet.
		std::vector<int> subeventnumthisevent(numdigits,-1);
		// LOOP OVER SUBEVENTS
		// -------------------
		for(unsigned int thissubevent=0; thissubevent<subeventhittimesv.size(); thissubevent++){
			float endtime = (thissubevent<(subeventhittimesv.size()-1)) ?
								//subeventhittimesv.at(thissubevent+1) : (eventendtime+1.);
								subeventendtimesv.at(thissubevent) : (eventendtime+1.);
			Log("TimeClustering Tool: Endtime for subevent "+to_string(thissubevent)+" is "+to_string(endtime), v_debug,verbosity);
			// don't need to worry about search start time as we start from earliest hit and exclude already allocated hits
			// times are not ordered, so scan through all digits for each subevent
			// FIND HITS IN THIS SUBEVENT
			// --------------------------
			for(int thisdigit=0;thisdigit<numdigits;thisdigit++){
				if(subeventnumthisevent.at(thisdigit)<0 && mrddigittimesthisevent.at(thisdigit)<= endtime ){
					// thisdigit is in thissubevent
					if(verbosity>5){
						cout<<"adding digit id "<<thisdigit<<" and digit at "<<mrddigittimesthisevent.at(thisdigit)<<" to subevent "<<thissubevent<<endl;
					}
					digitidsinasubevent.push_back(thisdigit);
					subeventnumthisevent.at(thisdigit)=thissubevent;
					tubeidsinasubevent.push_back(mrddigitpmtsthisevent.at(thisdigit));
					digittimesinasubevent.push_back(mrddigittimesthisevent.at(thisdigit));
					digitchargesinasubevent.push_back(mrddigitchargesthisevent.at(thisdigit));
					if (MakeMrdDigitTimePlot){
						mrddigitts_file->cd();
						if (MakeSingleEventPlots) mrddigitts_single->Fill(mrddigittimesthisevent.at(thisdigit));
					}
				}
			}
			
			// CONSTRUCT THE SUBEVENT
			// -----------------------
			if(int(digitidsinasubevent.size())>=minimumdigits){  // must have enough for a subevent
				Log("TimeClustering Tool: Constructing subevent "+to_string(mrdeventcounter)
					+" with "+to_string(digitidsinasubevent.size())+" digits",v_debug,verbosity);
				//MrdTimeClusters.Add(mrdeventcounter, digitidsinasubevent, tubeidsinasubevent, digittimesinasubevent);
				MrdTimeClusters.push_back(digitidsinasubevent);
				mrdeventcounter++;
				if (MakeMrdDigitTimePlot){
					for (unsigned int i_time = 0; i_time < digittimesinasubevent.size(); i_time++){
						mrddigitts_file->cd();
						mrddigitts_cluster->Fill(digittimesinasubevent.at(i_time));
						hist_chankey_cluster->Fill(mrdpmtid_to_channelkey[tubeidsinasubevent.at(i_time)]);
						hist_chankey_time_cluster->Fill(mrdpmtid_to_channelkey[tubeidsinasubevent.at(i_time)],digittimesinasubevent.at(i_time));
						if (MakeSingleEventPlots) mrddigitts_cluster_single->Fill(digittimesinasubevent.at(i_time));
						if (MRDTriggertype == "Cosmic") mrddigitts_cosmic_cluster->Fill(digittimesinasubevent.at(i_time));
						else if (MRDTriggertype == "Beam") mrddigitts_beam_cluster->Fill(digittimesinasubevent.at(i_time));
						else if (MRDTriggertype == "No Loopback") mrddigitts_noloopback_cluster->Fill(digittimesinasubevent.at(i_time));
					}
				}
			}
			
			// clear the vectors and loop to the next subevent
			digitidsinasubevent.clear();
			tubeidsinasubevent.clear();
			digittimesinasubevent.clear();
			digitchargesinasubevent.clear();
		}
		
		// quick scan to check for any unallocated hits
		for(unsigned int k=0;k<subeventnumthisevent.size();k++){
			if(subeventnumthisevent.at(k)==-1){
				Log("TimeClustering Tool: Found unbinned hit "+to_string(k)+" at "+to_string(mrddigittimesthisevent.at(k)),v_error,verbosity);
			}
		}
		
	}  // end multiple subevents case
	
	// verbose check about what is contained in MrdTimeClusters object
	for (unsigned int i_cluster=0; i_cluster < MrdTimeClusters.size(); i_cluster++){
		if (verbosity >= v_debug) std::cout << "TimeClustering tool: Cluster " << i_cluster+1 << ", MrdTimeClusters.at(i_cluster).size() = " << MrdTimeClusters.at(i_cluster).size() << std::endl;
			for (unsigned int i_entry = 0; i_entry < MrdTimeClusters.at(i_cluster).size(); i_entry++){
				if (verbosity >= v_debug) std::cout <<MrdTimeClusters.at(i_cluster).at(i_entry)<<", ";
			}
			if (verbosity >= v_debug) std::cout << std::endl;
	}
	
	for (unsigned int i=0; i< mrddigittimesthisevent.size(); i++){
		if (verbosity > v_debug) std::cout <<"mrddigitpmts, entry "<<i<<", time: "<<mrddigittimesthisevent.at(i)<<", pmt: "<<mrddigitpmtsthisevent.at(i)<<", charge: "<<mrddigitchargesthisevent.at(i)<<std::endl;
}


	if (MakeMrdDigitTimePlot){
	
		mrddigitts_file->cd();
		if (MakeSingleEventPlots){
			mrddigitts_cluster_single->Write();
			mrddigitts_single->Write();
		}
		gROOT->cd();
	}
	
	// pass the found clusters to the ANNIEEvent
	m_data->CStore.Set("MrdTimeClusters",MrdTimeClusters);
	m_data->CStore.Set("NumMrdTimeClusters",mrdeventcounter);
	m_data->CStore.Set("MrdDigitTimes",mrddigittimesthisevent);
	m_data->CStore.Set("MrdDigitPmts",mrddigitpmtsthisevent);
	m_data->CStore.Set("MrdDigitChankeys",mrddigitchankeysthisevent);
	m_data->CStore.Set("MrdDigitCharges",mrddigitchargesthisevent);
	
	//only for debugging
	//std::cout <<"TimeClustering tool: List of objects (End of Execute): "<<std::endl;
	//gObjectTable->Print();
	
	return true;
}


bool TimeClustering::Finalise(){
	
	// write time cluster histograms to file
	
	if (MakeMrdDigitTimePlot){
		mrddigitts_file->cd();
		mrddigitts_cosmic_cluster->Write();
		mrddigitts_beam_cluster->Write();
		mrddigitts_noloopback_cluster->Write();
		mrddigitts_cluster->Write();
		mrddigitts_cosmic->Write();
		mrddigitts_beam->Write();
		mrddigitts_noloopback->Write();
		mrddigitts->Write();
		mrddigitts_horizontal->Write();
		mrddigitts_vertical->Write();
		hist_chankey->Write();
		hist_chankey_cluster->Write();
		hist_chankey_time->Write();
		hist_chankey_time_cluster->Write();
		hist_chankey_multi->Write();
		mrddigitts_file->Close();
		delete mrddigitts_file;     //histograms get deleted by deleting associated TFile
		mrddigitts_file=0;
		gROOT->cd();
	}
	
	if (LaunchTApplication){
		// see if we're the last user of the TApplication and release it if so,
		// otherwise de-register us as a user since we're done
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				std::cout<<"TimeClustering Tool: Deleting global TApplication"<<std::endl;
				delete rootTApp;
				rootTApp=nullptr;
			}
		} else if(tapplicationusers>1){
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers-1);
		}
	}
	
	return true;
}
