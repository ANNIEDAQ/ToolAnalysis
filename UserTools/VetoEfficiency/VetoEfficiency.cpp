/* vim:set noexpandtab tabstop=4 wrap */
#include "VetoEfficiency.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

#include "ADCPulse.h"
#include "ANNIEconstants.h"
#include "BeamStatus.h"
#include "BoostStore.h"
#include "HeftyInfo.h"
#include "MinibufferLabel.h"
#include "PhaseITreeMaker.h"
#include "TimeClass.h"

// Excluded ADC channel IDs
//  6 = NCV PMT #1 (card 4, channel 1)
// 19 = neutron calibration source trigger input (card 8, channel 2)
// 37 = cosmic trigger input (card 14, channel 0)
// 49 = NCV PMT #2 (card 18, channel 0)
// 61 = summed signals from front veto (card 21, channel 0)
// 62 = summed signals from MRD 2 (card 21, channel 1)
// 63 = RWM (card 21, channel 2)
// 64 = summed signals from MRD 3 (card 21, channel 3)
constexpr std::array<uint32_t, 56> water_tank_pmt_IDs = {
  1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23,
  24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43,
  44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60
};

VetoEfficiency::VetoEfficiency():Tool(){}

bool VetoEfficiency::Initialise(std::string configfile, DataModel &data){
	
	std::cout<<"VetoEfficiency Tool: Initialising"<<std::endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); // loading config file
	//m_variables.Print();
	
	useTApplication=false;
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("drawHistos",drawHistos);
	m_variables.Get("plotDirectory",plotDirectory);
	m_variables.Get("min_unique_water_pmts",min_unique_water_pmts_);
	m_variables.Get("min_tank_charge",min_tank_charge_);
	m_variables.Get("min_mrdl1_pmts",min_mrdl1_pmts_);
	m_variables.Get("min_mrdl2_pmts",min_mrdl2_pmts_);
	m_variables.Get("coincidence_tolerance",coincidence_tolerance_);
	m_variables.Get("pre_trigger_ns",pre_trigger_ns_);
	m_variables.Get("useTApplication",useTApplication);
	get_ok = m_variables.Get("outputfilename",outputfilename);
	// use output name if specified, otherwise try to derive from input filename
	if((not get_ok) || (outputfilename=="auto")){
		get_ok = m_data->CStore.Get("InputFile",outputfilename);
		if(not get_ok){
			// don't have an input filename??
			logmessage="VetoEfficiency Tool: Error! Could not find 'InputFile' in CStore!";
			logmessage+="Will use default filename 'veto_efficiency.root'";
			Log(logmessage,v_error,verbosity);
			outputfilename = "veto_efficiency.root"; // default
		} else {
			// trim the extention
			std::string fname = outputfilename.substr(0,outputfilename.find_last_of("."));
			// trim the path, if it contains one
			uint last_slash = fname.find_last_of("/");
			if(last_slash!=std::string::npos){ fname = fname.substr(last_slash+1,std::string::npos); }
			// build our output name
			outputfilename="veto_efficiency_"+fname+".root";
			debugfilename="veto_debug_"+fname+".root";
		}
		Log("VetoEfficiency Tool: Writing output to "+outputfilename,v_warning,verbosity);
	}
	
	makeOutputFile(plotDirectory+"/"+outputfilename);
	
	if(drawHistos){
		// create the ROOT application to show histograms
		if (useTApplication){
			Log("VetoEfficiency Tool: getting/making TApplication",v_debug,verbosity);
			int myargc=0;
			//char *myargv[] = {(const char*)"mrddist"};
			// get or make the TApplication
			intptr_t tapp_ptr=0;
			get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
			if(not get_ok){
				Log("VetoEfficiency Tool: Making global TApplication",v_error,verbosity);
				rootTApp = new TApplication("rootTApp",&myargc,0);
				tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
				m_data->CStore.Set("RootTApplication",tapp_ptr);
			} else {
				Log("VetoEfficiency Tool: Retrieving global TApplication",v_error,verbosity);
				rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
			}
			int tapplicationusers;
			get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
			if(not get_ok) tapplicationusers=1;
			else tapplicationusers++;
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		
		}
		f_veto = new TFile(debugfilename.c_str(),"RECREATE");

		//The first histograms are built from all events (not restricted to events with veto hits)
		//To evaluate general time differences between the subsystems
		h_all_adc_times = new TH1F("h_all_adc_times","All ADC times",500,0,10000);
		h_all_mrd_times = new TH1F("h_all_mrd_times","All MRD times",500,0,10000);
		h_all_veto_times = new TH1F("h_all_veto_times","All Veto times",100,0,4000);
		h_all_veto_times_layer2 = new TH1F("h_all_veto_times_layer2","All Veto times Layer 2",100,0,4000);
		h_veto_delta_times = new TH1F("h_veto_delta_times","Veto Delta times",500,-2000,2000); //Time difference between the two veto layers
		h_veto_delta_times_coinc = new TH1F("h_veto_delta_times_coinc","Veto Delta times (coincidence)",500,-2000,2000); //Time difference between veto layers for veto-coincidences
		h_veto_2D = new TH2F("h_veto_2D","Veto Chankeys map",13,0,13,13,0,13); //Correlation between veto L1/L2 in-layer channel number
		for (int i=0;i<26;i++){
			std::stringstream temp;
			temp << "h_veto_delta_times_" << i ;
			TH1F *htemp = new TH1F(temp.str().c_str(),temp.str().c_str(),500,-2000,2000);
			vector_all_adc_times.push_back(htemp);
		}
		gROOT->cd();

		// these histograms are filled for events where we had at least one upstream veto hit,
		// but otherwise just plot ALL ADC and TDC times in the minibuffer, with times
		// corrected to be relative to the beam (not the minibuffer)
		h_adc_times = new TH1F("h_adc_times","ADC times",100,0,2000);
		h_tdc_times = new TH1F("h_tdc_times","TDC times",100,0,2000);
		h_tdc_times->SetLineColor(kRed);
		h_adc_delta_times = new TH1F("h_adc_delta_times","ADC Delta times",500,-2000,2000);
		h_adc_delta_times_charge = new TH2F("h_adc_delta_times_charge","ADC Delta times vs charge",200,0,5,200,-2000,2000);
		h_tdc_delta_times = new TH1F("h_tdc_delta_times","TDC Delta times",500,-2000,2000);
		h_tdc_delta_times_L1 = new TH1F("h_tdc_delta_times_L1","TDC Delta times MRD L1",500,-2000,2000);
		h_tdc_delta_times_L2 = new TH1F("h_tdc_delta_times_L2","TDC Delta times MRD L2",500,-2000,2000);
		h_tdc_delta_times->SetLineColor(kRed);
		h_adc_charge = new TH1F("h_adc_charge","ADC charges",500,0,5);
		h_adc_charge_coinc = new TH1F("h_adc_charge_coinc","ADC charges coincidence",500,0,5);
		h_tdc_chankeys_l1 = new TH1F("h_tdc_chankeys_l1","Veto L1 Chankeys",20,1000000,1002000);	//Histogram of occurence of L1 channelkeys
		h_tdc_chankeys_l2 = new TH1F("h_tdc_chankeys_l2","Veto L2 Chankeys",20,1010000,1012000);	//Histogram of occurence of L2 channelkeys
		h_tdc_chankeys_l1_l2 = new TH2F("h_tdc_chankeys_l1_l2","Veto L1/L2 Chankeys",20,1000000,1002000,20,1010000,1012000);	//Correlation between unconverted L1/L2 channelkeys

		// this histogram is also only filled for events where we have an upstream veto layer,
		// and plots ALL ADC hit times relative to the veto event time
		// restricting the timespan to the coincidence required region means this plots the
		// times of alll ADC hits considered as part of the coincidence event
		h_tankhit_time_to_coincidence = new TH1F("h_tankhit_time_to_coincidence", "Tank Hit Times in Coincidence Window", 100,-pre_trigger_ns_, coincidence_tolerance_); // -2000 --> 2000
		// this plots the times of coincidence events - i.e. upstream veto hits minus pre_trigger_ns
		h_coincidence_event_times = new TH1F("h_coincidence_event_times","Start times of Coincidence Events",100,0,4100);
	}
	
	// Get the geometry
	get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
	if(not get_ok){
		Log("VetoEfficiency Tool: Failed to obtain AnnieGeometry from ANNIEEvent header!",v_error,verbosity);
		return false;
	}
	
	// generate the channelkeys for the various MRD/Veto layers, so we can check them for hits
	LoadTDCKeys();
	if(verbosity>v_debug){
		std::cout<<"\nList of veto L1 PMTs: {";
		for(auto&& akey : vetol1keys) std::cout<<akey<<", ";
		std::cout<<"\b\b}\nList of veto L2 PMTs: {";
		for(auto&& akey : vetol2keys) std::cout<<akey<<", ";
		std::cout<<"\b\b}\nList of MRD L1 PMTs: {";
		for(auto&& akey : mrdl1keys) std::cout<<akey<<", ";
		std::cout<<"\b\b}\nList of MRD L2 PMTs: {";
		for(auto&& akey : mrdl2keys) std::cout<<akey<<", ";
		std::cout<<"\b\b}"<<std::endl;
	}
	
	return true;
}


bool VetoEfficiency::Execute(){
	
	Log("VetoEfficiency Tool: Executing",v_debug,verbosity);
	
	// A lot of code here taken from the PhaseITreeMaker tool, since it handles
	// minibuffer validation, POT counting, timestamp fixing etc etc.,
	// but does so in the context of NCV coincidences_, so can't be accessed easily
	// Ultimately what we want out of it is the selection of good data,
	// a measurement of how much data we're analysing,
	// and most critically a measure the tank charge/num unique PMTs hit.
	// This last variable is what we're going to need for the coincidence event
	// for the veto efficiency measurement.
	// ====================================================================
	
	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	// Begin Minibuffer Selection
	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	
	// Get a pointer to the ANNIEEvent Store
	auto* annie_event = m_data->Stores["ANNIEEvent"];
	
	if (!annie_event) {
		Log("Error: The VetoEfficiency tool could not find the ANNIEEvent Store",
			0, verbosity);
		return false;
	}
	
	Log("VetoEfficiency Tool: Getting event information",v_debug,verbosity);
	
	// Load the labels describing the type of data stored in each minibuffer
	std::vector<MinibufferLabel> mb_labels;
	m_data->Stores["ANNIEEvent"]->Get("MinibufferLabels", mb_labels);
	
	get_object_from_store("MinibufferLabels", mb_labels, *annie_event);
	check_that_not_empty("MinibufferLabels", mb_labels);
	
	// Decide whether we're using Hefty vs. non-Hefty data by checking whether
	// the HeftyInfo object is present
	hefty_mode_ = annie_event->Has("HeftyInfo");
	
	// One of these objects will be used to get minibuffer timestamps
	// depending on whether we're using Hefty mode or non-Hefty mode.
	HeftyInfo hefty_info; // used for Hefty mode only
	std::vector<TimeClass> mb_timestamps; // used for non-Hefty mode only
	size_t num_minibuffers = 0u;
	
	if (hefty_mode_) {
		get_object_from_store("HeftyInfo", hefty_info, *annie_event);
		num_minibuffers = hefty_info.num_minibuffers();
		
		if ( num_minibuffers == 0u ) {
			Log("Error: The VetoEfficiency tool found an empty HeftyInfo entry", 0,
				verbosity);
			return false;
		}
		
		// Exclude beam spills (or source triggers) near the end of a full
		// multi-minibuffer readout that included extra self-triggers in the Hefty
		// window that could not be recorded. This is indicated in the heftydb
		// TTree via More[39] == 1 and in the HeftyInfo object by more() == true.
		if ( hefty_info.more() ) {
			// Find the first beam or source minibuffer counting backward from the
			// end of the full readout
			size_t mb = num_minibuffers - 1u;
			for (; mb > 0u; --mb) {
				MinibufferLabel label = mb_labels.at(mb);
				if ( label == MinibufferLabel::Beam
					|| label == MinibufferLabel::Source ) break;
			}
			// Exclude the minibuffers from the incomplete beam or source trigger's
			// Hefty window by setting a new value of num_minibuffers. This will
			// prematurely end the loop over minibuffers below.
			num_minibuffers = mb;
		}
	}
	else {
		// non-Hefty data
		get_object_from_store("MinibufferTimestamps", mb_timestamps, *annie_event);
		check_that_not_empty("MinibufferTimestamps", mb_timestamps);
		num_minibuffers = mb_timestamps.size();
		// Trigger masks are not saved in the tree for non-Hefty mode
		hefty_trigger_mask_ = 0;
	}
	
	// Load the beam status objects for each minibuffer
	std::vector<BeamStatus> beam_statuses;
	get_object_from_store("BeamStatuses", beam_statuses, *annie_event);
	check_that_not_empty("BeamStatuses", beam_statuses);
	
	// Load run, subrun, and event numbers
	get_object_from_store("RunNumber", run_number_, *annie_event);
	get_object_from_store("SubRunNumber", subrun_number_, *annie_event);
	get_object_from_store("EventNumber", event_number_, *annie_event);
	
	// Load the reconstructed ADC hits
	std::map<unsigned long, std::vector< std::vector<ADCPulse> > > adc_hits;
	get_object_from_store("RecoADCHits", adc_hits, *annie_event);
	check_that_not_empty("RecoADCHits", adc_hits);
	
	for (size_t mb = 0; mb < num_minibuffers; ++mb) {
		for (const auto& channel_pair : adc_hits) {
			const unsigned long& ck = channel_pair.first;
			Detector* det = anniegeom->ChannelToDetector(ck);
			if(det==nullptr){
				Log("VetoEfficiency Tool: Warning! ADC hit on channel "+to_string(ck)
				+" with no corresponding Detector obect in Geometry!",v_warning,verbosity);
				continue;
			}
			if (det->GetDetectorElement() != "Tank") continue;
		
			// Skip non-water-tank PMTs
			if ( std::find(std::begin(water_tank_pmt_IDs), std::end(water_tank_pmt_IDs),
			 	uint32_t(ck)) == std::end(water_tank_pmt_IDs) ) continue;
			
			const auto& pulse_vec = channel_pair.second.at(mb);
			//if (pulse_vec.size()>0) std::cout <<"Chankey "<<ck<<", pulse vector size: "<<pulse_vec.size()<<std::endl;
			for (const auto& pulse: pulse_vec){
				size_t pulse_time = pulse.start_time();
				h_all_adc_times->Fill(pulse_time);
			}
		}
	}
	


	// Load the TDC hits
	std::map<unsigned long,std::vector<std::vector<Hit>>>* TDCData=nullptr;
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
	if(not get_ok){
		Log("VetoEfficiency Tool: Error! Failed to get TDCData from ANNIEEvent!",v_error, verbosity);
		return false;
	}
	
	// Flag that vetos minibuffers because the last beam minibuffer failed
	// the quality cuts.
	bool beam_veto_active = false;
	
	Log("VetoEfficiency Tool: Looping over "+to_string(num_minibuffers)+ " minibuffers",v_debug,verbosity);
	
	// loop over minibuffers, for each channel...
	// seems like this channel->minibuffer loop nest should probably be inverted...?
	for (size_t mb = 0; mb < num_minibuffers; ++mb) {
	
		veto_times.clear();
		veto_chankeys.clear();
		veto_times_layer2.clear();
		veto_chankeys_layer2.clear();
		
		//Log("VetoEfficiency Tool: Checking event quality for minibuffer " + to_string(mb),v_debug,verbosity);
		
		minibuffer_number_ = mb;
		
		// Determine the correct label for the events in this minibuffer
		MinibufferLabel event_mb_label = mb_labels.at(mb);
		event_label_ = static_cast<uint8_t>( event_mb_label );
		
		// If this is Hefty mode data, save the trigger mask for the
		// current minibuffer. This is distinct from the MinibufferLabel
		// assigned to the event_label_ variable above. The trigger
		// mask branch isn't used for non-Hefty data.
		if ( hefty_mode_ ) hefty_trigger_mask_ = hefty_info.label(mb);
		else hefty_trigger_mask_ = 0;
		
		// BEAM QUALITY CUTS
		// Skip beam minibuffers with bad or missing beam status information
		// TODO: consider printing a warning message here
		const auto& beam_status = beam_statuses.at(mb);
		
		const auto& beam_condition = beam_status.condition();
		if (beam_condition == BeamCondition::Missing
			|| beam_condition == BeamCondition::Bad)
		{
			// Skip all beam and Hefty window minibuffers until a good-quality beam
			// spill is found again
			beam_veto_active = true;
		}
		if ( beam_veto_active && beam_condition == BeamCondition::Ok ) {
			// We've found a new beam minibuffer that passed the quality check,
			// so disable the beam quality veto
			beam_veto_active = false;
		}
		if (beam_veto_active && (event_mb_label == MinibufferLabel::Hefty
				|| event_mb_label == MinibufferLabel::Beam))
		{
			// Bad beam minibuffers and Hefty window minibuffers belonging to the
			// bad beam spill need to be skipped. Since other minibuffers (e.g.,
			// cosmic trigger minibuffers) are not part of the beam "macroevent,"
			// they may still be processed normally.
			continue;
		}
		
		// NOTE: for minibuffer labels other than "Beam", "Source", and "Hefty"
		// no t_since_beam offset is available (it wouldn't make much sense anyway)
		// so all times are relative to the start of the (arbitrary) minibuffer.
		// We can flag these events as being 'out_of_spill' since they're not beam related
		bool in_spill = false;
		if ( !hefty_mode_ || (event_mb_label == MinibufferLabel::Hefty
			|| event_mb_label == MinibufferLabel::Beam
			|| event_mb_label == MinibufferLabel::Source) )
		{
			in_spill = true;
		}
		else in_spill = false;
		
		// only process beam events, not cosmic, source etc.
		if(not in_spill) continue;
		
		// breakdown of beam condition cuts. Presumably if any of these are false,
		// BeamCondition will be bad and we'll skip the event anyway?
		if ( beam_status.is_beam() ) {
			const auto& cuts = beam_status.cuts();
			bool pot_ok = cuts.at("POT in range");
			bool horn_current_ok = cuts.at("peak horn current in range");
			bool timestamps_ok = cuts.at("toroids agree");
			bool toroids_agree = cuts.at("timestamps agree");
			if( (not pot_ok) ||
				(not horn_current_ok) ||
				(not timestamps_ok) ||
				(not toroids_agree)
			  ){
				Log("VetoEfficiency Tool: Warning! beam_status cuts are not all passed,"
					" but event was to be processed anyway?",v_warning,verbosity);
				continue;
			}
		}
		
		is_cosmic = false;
		// Increment num spills, beam POT count, etc. based on the characteristics of the
		// current minibuffer
		if (beam_condition == BeamCondition::Ok) {
			total_POT_ += beam_status.pot();
			++spill_number_;
		}
		
		// if not a beam event, increment counts of num corresponding trigger types
		// (redundant as we're not processing these)
		if ( event_mb_label == MinibufferLabel::Source ) {
			++num_source_triggers_;
		}
		else if ( event_mb_label == MinibufferLabel::Cosmic ) {
			is_cosmic = true;
			++num_cosmic_triggers_;
		}
		else if ( event_mb_label == MinibufferLabel::Soft ) {
			++num_soft_triggers_;
		}
		
		Log("VetoEfficiency Tool: Minibuffer "+to_string(mb)+" passed event selection",v_debug,verbosity);
		//std::cout<<"hefty t_since_beam for mb "<<mb<<" = "<<hefty_info.t_since_beam(mb)<<std::endl;
		
		/////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////
		// End Minibuffer Selection
		/////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////
		
		// we now need to look for 2-fold coincidences_ between the upstream veto layer
		// and some downstream component - either tank or MRD.
		// since we must have an upstream veto hit, let's find those first,
		// then scan the downstream components for a match on each
		Log("VetoEfficiency Tool: Extracting upstream veto layer hits",v_debug,verbosity);
		
		if(verbosity>v_debug){
			int nvetol1hits=0, nvetol2hits=0, nmrdl1hits=0, nmrdl2hits=0, notherhits=0;
			for(auto&& atdcchannel : *TDCData){
				unsigned long thekey = atdcchannel.first;
				int nhits = atdcchannel.second.at(mb).size();
				ntotaltdchits_+=nhits;
				if(nhits){
					if(std::find(vetol1keys.begin(), vetol1keys.end(), thekey)!=vetol1keys.end()){
						nvetol1hits++; nvetol1hits_++;
					} else if(std::find(vetol2keys.begin(), vetol2keys.end(), thekey)!=vetol2keys.end()){
						nvetol2hits++; nvetol2hits_++;
					} else if(std::find(mrdl1keys.begin(), mrdl1keys.end(), thekey)!=mrdl1keys.end()){
						nmrdl1hits++; nmrdl1hits_++;
					} else if(std::find(mrdl2keys.begin(), mrdl2keys.end(), thekey)!=mrdl2keys.end()){
						nmrdl2hits++; nmrdl2hits_++;
					} else {
						notherhits++; notherhits_++;
					}
					int corrected_tubeid = thekey-1000000;
					int mrdpmtxnum=floor(corrected_tubeid/10000);
					int mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
					int mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
					//printf("VetoEfficiency: TDC hit on PMT %02d:%02d:%02d in minibuffer %d\n",
					//        mrdpmtxnum,mrdpmtynum,mrdpmtznum,mb);
				}
			}
			std::cout<<"hit counts so far: "
					 <<"("<<nvetol1hits<<"/"<<nvetol1hits_<<"), "
					 <<"("<<nvetol2hits<<"/"<<nvetol2hits_<<"), "
					 <<"("<<nmrdl1hits<<"/"<<nmrdl1hits_<<"), "
					 <<"("<<nmrdl2hits<<"/"<<nmrdl2hits_<<"), "
					 <<"("<<notherhits<<"/"<<notherhits_<<"), "
					 <<" == "<<ntotaltdchits_<<" "<<std::endl;
		}
		

		// --------------------------------------------------------------------
		// --------- Fill a lot of debug histograms ---------------------------
		// --------- before looking at veto layer -----------------------------
		// --------------------------------------------------------------------

		if(drawHistos){
			Log("VetoEfficiency Tool: Filling debug histograms",v_debug,verbosity);
			// First get TDC data: loop over TDC channels with a hit in this readout
			// note that not all TDC channels may be present; if a TDC channel
			// didn't have any hits in any of the minibuffers being processed
			// it will not be present
			for(auto&& a_tdc_channel : *TDCData){
				unsigned long tdc_key = a_tdc_channel.first;
				// channelkeys for TDC channels are representative of their location in the MRD system
				h_tdc_chankeys_l1->Fill(tdc_key);
				h_tdc_chankeys_l2->Fill(tdc_key);
				int corrected_tubeid = tdc_key-1000000;
				int mrdpmtxnum=floor(corrected_tubeid/10000);
				int mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
				int mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
				
				// get any hits on this channel in the present minibuffer
				const std::vector<std::vector<Hit>>& tdc_channel_data = a_tdc_channel.second;
				std::vector<Hit> hits_on_this_tdc_channel_in_this_mb = tdc_channel_data.at(mb);
				
				// process the hits
				for(Hit& next_tdc_hit : hits_on_this_tdc_channel_in_this_mb){
					double pulse_start_time_ns = next_tdc_hit.GetTime();
					
					// TDC hit times (like ADC ones) are relative to the start of the minibuffer
					// if using hefty mode, to get times relative to beam,
					// we should add the t_since_beam for the corresponding minibuffer
					if ( hefty_mode_ && event_mb_label == MinibufferLabel::Hefty) {
						pulse_start_time_ns += hefty_info.t_since_beam(mb);
					}
					
					h_all_mrd_times->Fill(pulse_start_time_ns);
				} // end loop over tdc hits on this channel
			} // end loop over tdc channels
			
			// Next the ADC data
			for (const auto& pair : adc_hits) {
				const auto& channel_key = pair.first;
				const auto& minibuffer_pulses = pair.second;
				// Only show ADC channels that belong to water tank PMTs
				Detector* det = anniegeom->ChannelToDetector(channel_key);
				if(det==nullptr){
					Log("VetoEfficiency Tool: Warning! ADC hit on channel "+to_string(channel_key)
						+" with no corresponding Detector obect in Geometry!",v_warning,verbosity);
					continue;
				}
				if (det->GetDetectorElement() != "Tank") {
					continue;
				}
				// Skip non-water-tank PMTs
				if ( std::find(std::begin(water_tank_pmt_IDs), std::end(water_tank_pmt_IDs),
					 uint32_t(channel_key)) == std::end(water_tank_pmt_IDs) )  continue;
				// loop over pulses on this channel, in this minibuffer
				for (const ADCPulse& pulse : minibuffer_pulses.at(mb) ) {
					double pulse_amplitude_ = pulse.amplitude();                  // Volts
					double pulse_charge_ = pulse.charge();                        // nC
					unsigned short pulse_raw_amplitude_ = pulse.raw_amplitude();  // ADC counts
					int64_t pulse_start_time_ns = pulse.start_time();
					if ( hefty_mode_ && event_mb_label == MinibufferLabel::Hefty) {
						pulse_start_time_ns = pulse.start_time() + hefty_info.t_since_beam(mb);
					}
					h_all_adc_times->Fill(pulse_start_time_ns);
				} // end loop over adc hits on this channel
			} // end loop over adc channels
		} // end draw debug histos
		
		// make a vector of the times of upstream veto layer hits
		// and pmts
		std::vector<std::pair<double, unsigned long>> veto_l1_hits;
		// for each upstream veto channel
		for(unsigned long avetol1channel : vetol1keys){
			// did it have any hits this readout
			if(TDCData->count(avetol1channel)){
				// loop over any hits in this minibuffer
				for(auto& ahit : TDCData->at(avetol1channel).at(mb)){
					veto_l1_hits.push_back(std::pair<double,unsigned long>{ahit.GetTime(),avetol1channel});
					veto_times.push_back(ahit.GetTime());
					veto_chankeys.push_back(avetol1channel);
					Log("VetoEfficiency Tool: L1 veto hit on PMT "+to_string(avetol1channel)
						+" at "+to_string(ahit.GetTime()),v_debug+1,verbosity);
				}
			}
		}
		std::vector<std::pair<double, unsigned long>> veto_l2_hits;
		// for each upstream veto channel
		for(unsigned long avetol2channel : vetol2keys){
			// did it have any hits this readout
			if(TDCData->count(avetol2channel)){
				// loop over any hits in this minibuffer
				for(auto& ahit : TDCData->at(avetol2channel).at(mb)){
					veto_l2_hits.push_back(std::pair<double,unsigned long>{ahit.GetTime(),avetol2channel});
					Log("VetoEfficiency Tool: L2 veto hit on PMT "+to_string(avetol2channel)
						+" at "+to_string(ahit.GetTime()),v_debug+1,verbosity);
				}
			}
		}
		for (auto&& this_hit : veto_l2_hits){
			h_all_veto_times_layer2->Fill(this_hit.first);
			veto_times_layer2.push_back(this_hit.first);
			veto_chankeys_layer2.push_back(this_hit.second);
		}
		
		found_coincidence=false;

		//Check very roughly for coincidences of veto L1/L2 channels (mostly for debugging purposes)
		for (int il1=0; il1<(int)veto_times.size();il1++){
			int in_layer_index = std::distance(vetol1keys.begin(),
                        std::find(vetol1keys.begin(), vetol1keys.end(), veto_chankeys.at(il1)));
				for (int il2=0; il2<(int)veto_times_layer2.size();il2++){
					int in_layer_index_layer2 = std::distance(vetol2keys.begin(),
                                                        std::find(vetol2keys.begin(), vetol2keys.end(), veto_chankeys_layer2.at(il2)));
					if (in_layer_index != in_layer_index_layer2) continue;
					if (verbosity >= v_debug){
						std::cout <<"Minibuffer "<<mb<<std::endl;
						std::cout <<"Layer 1 hit, channel "<<veto_chankeys.at(il1)<<", time "<<veto_times.at(il1)<<std::endl;
						std::cout <<"Layer 2 hit, channel "<<veto_chankeys_layer2.at(il2)<<", time "<<veto_times_layer2.at(il2)<<std::endl;
					}
					h_veto_delta_times_coinc->Fill(veto_times_layer2.at(il2)-veto_times.at(il1));	
					found_coincidence = true;
			}

		}

		// --------------------------------------------------------------------
		// --------- End of debug histogram filling ---------------------------
		// --------- Start real coincidence search here------------------------
		// --------------------------------------------------------------------

		bool layer1_hit = true;
		bool layer2_hit = true;
		if(veto_l1_hits.size()==0){
			Log("VetoEfficiency Tool: Found " + to_string(veto_l1_hits.size())
				+" hits on upstream veto layer",v_debug,verbosity);
			layer1_hit = false;
			//continue;
		}
		if(veto_l2_hits.size()==0){
			Log("VetoEfficiency Tool: Found " + to_string(veto_l2_hits.size())
				+" hits on downstream veto layer",v_debug,verbosity);
			layer2_hit = false;
			//continue;
		}

		//Evalute efficiencies for coincidences with the first layer of the veto
		if (layer1_hit){
			Log("VetoEfficiency Tool: Found " + to_string(veto_l1_hits.size())
				+" hits on upstream veto layer, constructing coincidence objects",v_debug,verbosity);
			
			// since we may have multiple vetol1 hits within a close proximity,
			// we should reduce these down before we start looking for matches
			// first sort the hits by time
			std::sort(veto_l1_hits.begin(), veto_l1_hits.end());
			

			// Create coincidences object for coincidence conditions with the first layer of the veto
			// we'll store info about each coincidence within a struct
			coincidences_.clear();
			for(auto&& this_hit : veto_l1_hits){
				// search any existing coincidences_ to see if this hit is close in time
				bool allocated=false;
				h_all_veto_times->Fill(this_hit.first);
				for(CoincidenceInfo& acoincidence : coincidences_){
					int in_layer_index = std::distance(vetol1keys.begin(),
								std::find(vetol1keys.begin(), vetol1keys.end(), this_hit.second));
					if((this_hit.first-acoincidence.event_time_ns)<coincidence_tolerance_){
						if(acoincidence.vetol1hits.count(this_hit.second)){
							acoincidence.vetol1hits.at(this_hit.second).push_back(this_hit.first);
						} else {
							acoincidence.vetol1hits.emplace(this_hit.second,std::vector<double>{this_hit.first});
						}
						allocated=true;
						break;
					}
				}
				// if it wasn't within the time window of any existing coincidences_, make a new one
				if(allocated==false){
					CoincidenceInfo newcoincidence;
					newcoincidence.event_time_ns = this_hit.first - pre_trigger_ns_;
					if(drawHistos) h_coincidence_event_times->Fill(newcoincidence.event_time_ns);
					newcoincidence.vetol1hits.emplace(this_hit.second, std::vector<double>{this_hit.first});
					coincidences_.push_back(newcoincidence);
				}
			}
			
			Log("VetoEfficiency Tool: Reduced to " + to_string(coincidences_.size())
				+ " CoincidenceInfo objects, searching for tank events",v_debug,verbosity);
			
			// Create coincidences object for coincidence conditions with the second layer of the veto
			for (auto&& this_hit : veto_l2_hits){
				for (int i_veto=0; i_veto<(int) veto_times.size(); i_veto++){
					h_veto_delta_times->Fill(this_hit.first-veto_times.at(i_veto));
					int in_layer_index = std::distance(vetol1keys.begin(),std::find(vetol1keys.begin(), vetol1keys.end(), veto_chankeys.at(i_veto)));	 
					int in_layer_index_layer2 = std::distance(vetol2keys.begin(),std::find(vetol2keys.begin(), vetol2keys.end(),this_hit.second));
					in_layer_index_layer2+=13;	 
					
					vector_all_adc_times.at(in_layer_index)->Fill(this_hit.first-veto_times.at(i_veto));
					if (in_layer_index_layer2<26) vector_all_adc_times.at(in_layer_index_layer2)->Fill(this_hit.first-veto_times.at(i_veto));
					h_veto_2D->Fill(in_layer_index,in_layer_index_layer2-13);
					h_tdc_chankeys_l1_l2->Fill(veto_chankeys.at(i_veto),this_hit.second);
				}
			}
			// this provides the set of reference times: now look for coincident events
			// in the rest of the detector. First calculate the tank charge and 
			// num PMTs hit within the tolerance window
			for(CoincidenceInfo& acoincidence : coincidences_){
				Log("VetoEfficiency Tool: Searching for tank charge and hit PMTs from "
					+to_string(acoincidence.event_time_ns)+" to "
					+to_string(acoincidence.event_time_ns+coincidence_tolerance_),v_debug+1,verbosity);
				// calculate tank charge and # PMTs hit within a window around this hit
				int num_unique_water_pmts = BOGUS_INT;
				int num_pmts_above_thr;
				double tank_charge = compute_tank_charge(mb,
					adc_hits, acoincidence.event_time_ns-1000, acoincidence.event_time_ns + coincidence_tolerance_ - 1000,
					num_unique_water_pmts, num_pmts_above_thr);
				// XXX XXX XXX XXX based on VetoEfficiency coincidences, ADC events are 1us
				// earlier than TDC events???? Where does this come from?? XXX XXX XXX XXX
				acoincidence.num_unique_water_pmts = num_unique_water_pmts;
				//acoincidence.num_unique_water_pmts = num_pmts_above_thr;
				acoincidence.tank_charge = tank_charge;
				h_adc_charge->Fill(tank_charge);
				if (found_coincidence) h_adc_charge_coinc->Fill(tank_charge);
				if (found_coincidence && verbosity >= v_message) std::cout <<"ADC result: charge = "<<tank_charge<<", unique PMTs: "<<num_unique_water_pmts<<", cosmic = "<<is_cosmic<<std::endl;
			}
			
			Log("VetoEfficiency Tool: Searching for MRD L1 activity",v_debug,verbosity);
			
			// we might also want to use MRD activity as an alternate selection
			// for through-going events, so find any MRD hits within the windows
			// we separate into L1 and L2, but probably either is fine
			for(const unsigned long& amrd1key : mrdl1keys){
				// did it have any hits this readout
				if(TDCData->count(amrd1key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(amrd1key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						if (found_coincidence && verbosity >= v_message) std::cout <<"MRD L1 hit time "<<ahittime<<std::endl;
						for(CoincidenceInfo& acoincidence : coincidences_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.mrdl1hits.count(amrd1key)){
									acoincidence.mrdl1hits.at(amrd1key).push_back(ahittime);
								} else {
									acoincidence.mrdl1hits.emplace(amrd1key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			Log("VetoEfficiency Tool: Searching for MRD L2 activity",v_debug,verbosity);
			// repeat for MRD L2 PMTs
			for(const unsigned long& amrd2key : mrdl2keys){
				// did it have any hits this readout
				if(TDCData->count(amrd2key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(amrd2key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						if (found_coincidence && verbosity >= v_message) std::cout <<"MRD L2 hit time "<<ahittime<<std::endl;
						for(CoincidenceInfo& acoincidence : coincidences_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.mrdl2hits.count(amrd2key)){
									acoincidence.mrdl2hits.at(amrd2key).push_back(ahittime);
								} else {
									acoincidence.mrdl2hits.emplace(amrd2key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			
			// finally the important bit: seeing if we have any veto L2 activity.
			// whether we had any or not is going to determine our efficiency
			Log("VetoEfficiency Tool: Searching for Veto L2 activity",v_debug,verbosity);
			for(const unsigned long& avetol2key : vetol2keys){
				// did it have any hits this readout
				if(TDCData->count(avetol2key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(avetol2key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						for(CoincidenceInfo& acoincidence : coincidences_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.vetol2hits.count(avetol2key)){
									acoincidence.vetol2hits.at(avetol2key).push_back(ahittime);
								} else {
									acoincidence.vetol2hits.emplace(avetol2key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			
			// write the coincidence events to a ROOT file for easier analysis
			Log("VetoEfficiency Tool: adding coincidences to output tree",v_debug,verbosity);
			int n_coincidences_passed_cuts=0;
			for(CoincidenceInfo& acoincidence : coincidences_){
				// if it passed either the tank or MRD coincidence requirements
				// for a through-going event (upstream veto requirement is always
				// satisfied by construction), then record it to file
				// Require only 1 MRD layer here, can do more stringent cuts
				// when evaulating the output file
				if( ((acoincidence.tank_charge > min_tank_charge_) &&
					 ((uint) acoincidence.num_unique_water_pmts > min_unique_water_pmts_)) ||
					((acoincidence.mrdl1hits.size() > min_mrdl1_pmts_) ||
					 (acoincidence.mrdl2hits.size() > min_mrdl2_pmts_))
					){
					n_coincidences_passed_cuts++;
					
					event_time_ns_ = acoincidence.event_time_ns;
					num_veto_l1_hits_ = acoincidence.vetol1hits.size();
					num_veto_l2_hits_ = acoincidence.vetol2hits.size();
					num_mrdl1_hits_ = acoincidence.mrdl1hits.size();
					num_mrdl2_hits_ = acoincidence.mrdl2hits.size();
					tank_charge_ = acoincidence.tank_charge;
					num_unique_water_pmts_ = acoincidence.num_unique_water_pmts;
					coincidence_info_ = acoincidence;
					coincidence_layer_ = 1;		//First Veto layer used to define coincidences		

					// event cuts have been applied: we have a through-going event.
					// the efficiency measurement:
					// Build a set of the veto PMTs hit this event
					veto_l1_ids_.clear();
					for(auto&& al1hit : acoincidence.vetol1hits){
						// add this pmt if we haven't already
						if(std::find(veto_l1_ids_.begin(),veto_l1_ids_.end(),al1hit.first)==veto_l1_ids_.end()){
							veto_l1_ids_.push_back(al1hit.first);
						}
					}
					veto_l2_ids_.clear();
					for(auto&& al2hit : acoincidence.vetol2hits){
						if(std::find(veto_l2_ids_.begin(),veto_l2_ids_.end(),al2hit.first)==veto_l2_ids_.end()){
							veto_l2_ids_.push_back(al2hit.first);
						}
					}
					
					// record them in the aggregate info
					for(auto&& al1pmt : veto_l1_ids_){
						// convert to index in the layer, i.e. 0-12
						int in_layer_index = std::distance(vetol1keys.begin(),
								std::find(vetol1keys.begin(), vetol1keys.end(), al1pmt));
						// safety check
						if(in_layer_index>=(int)l1_hits_L1_.size()){
							Log("VetoEfficiency Tool: Error! Hit on veto L1 PMT index "
								+to_string(in_layer_index)+" is out of range!",v_error,verbosity);
							return false;
						}
						// note that this paddle saw an event
						l1_hits_L1_.at(in_layer_index)++;
						// check if we had a corresponding hit on the L2 PMT
						// build the coresponding L2 key so we can search for it
						unsigned long l2_key = 1010000 + (100*in_layer_index);
						// see if that key is in our list of L2 keys hit this event
						if (verbosity >= v_message) std::cout <<"In-Layer-Chankey "<<in_layer_index<<" saw a hit! Check if l2_key = "<<l2_key<<" also saw something."<<std::endl;
						for (int i_layer2=0; i_layer2<(int)veto_l2_ids_.size(); i_layer2++){
							if (verbosity >= v_message) std::cout <<"Layer 2 hit with chankey "<<veto_l2_ids_.at(i_layer2)<<std::endl;
						}
						if(std::find(veto_l2_ids_.begin(),veto_l2_ids_.end(),l2_key)!=veto_l2_ids_.end()){
							l2_hits_L1_.at(in_layer_index)++;
							if (verbosity >= v_message) std::cout <<"Layer 2 saw a hit! :)"<<std::endl;
						}
						
					}
					// ultimately the ratio of l2_hits_.at(a_pmt_num) to l1_hits_.at(a_pmt_num)
					// will give us the efficiency of a_pmt_num in L2.
					roottreeout->Fill();
				}
			}
			Log("VetoEfficiency Tool: Found "+to_string(n_coincidences_passed_cuts)
				+" through-going events",v_debug,verbosity);
			if(verbosity>v_debug){
				std::cout<<"VetoEfficiency Tool: Current efficiencies are: {";
				for(int i=0; i<13; ++i){
					double theefficiency = (l1_hits_L1_.at(i)==0) ? 0 : (double(l2_hits_L1_.at(i))/double(l1_hits_L1_.at(i)))*100.;
					logmessage = to_string(l2_hits_L1_.at(i)) + "/"
							   + to_string(l1_hits_L1_.at(i)) + "="
							   + to_string(theefficiency);
					if(i<12) logmessage += ", ";
					std::cout<<logmessage;
				}
				std::cout<<"}"<<std::endl;
			}
			
			/// debug aside: fill some hisograms of all the TDC and ADC hits to check they're on the same footing
			if(drawHistos){
				Log("VetoEfficiency Tool: Filling debug histograms",v_debug,verbosity);
				// First get TDC data: loop over TDC channels with a hit in this readout
				// note that not all TDC channels may be present; if a TDC channel
				// didn't have any hits in any of the minibuffers being processed
				// it will not be present
				for(auto&& a_tdc_channel : *TDCData){
					unsigned long tdc_key = a_tdc_channel.first;
					// channelkeys for TDC channels are representative of their location in the MRD system
					int corrected_tubeid = tdc_key-1000000;
					int mrdpmtxnum=floor(corrected_tubeid/10000);
					int mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
					int mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
					
					// get any hits on this channel in the present minibuffer
					const std::vector<std::vector<Hit>>& tdc_channel_data = a_tdc_channel.second;
					std::vector<Hit> hits_on_this_tdc_channel_in_this_mb = tdc_channel_data.at(mb);
					
					// process the hits
					for(Hit& next_tdc_hit : hits_on_this_tdc_channel_in_this_mb){
						double pulse_start_time_ns = next_tdc_hit.GetTime();
						
						// TDC hit times (like ADC ones) are relative to the start of the minibuffer
						// if using hefty mode, to get times relative to beam,
						// we should add the t_since_beam for the corresponding minibuffer
						if ( hefty_mode_ && event_mb_label == MinibufferLabel::Hefty) {
							pulse_start_time_ns += hefty_info.t_since_beam(mb);
						}
						
						h_tdc_times->Fill(pulse_start_time_ns);
						for(int i_veto=0; i_veto < (int) veto_times.size(); i_veto++){
							if (std::find(vetol1keys.begin(),vetol1keys.end(),tdc_key)!=vetol1keys.end()) continue;
							if (std::find(vetol2keys.begin(),vetol2keys.end(),tdc_key)!=vetol2keys.end()) continue;
							h_tdc_delta_times->Fill(pulse_start_time_ns-veto_times.at(i_veto));
							if (std::find(mrdl1keys.begin(),mrdl1keys.end(),tdc_key)!=mrdl1keys.end()) h_tdc_delta_times_L1->Fill(pulse_start_time_ns-veto_times.at(i_veto));
							if (std::find(mrdl2keys.begin(),mrdl2keys.end(),tdc_key)!=mrdl2keys.end()) h_tdc_delta_times_L2->Fill(pulse_start_time_ns-veto_times.at(i_veto));
						}
					} // end loop over tdc hits on this channel
				} // end loop over tdc channels
				
				// Next the ADC data
				for (const auto& pair : adc_hits) {
					
					const auto& channel_key = pair.first;
					const auto& minibuffer_pulses = pair.second;
					
					// Only show ADC channels that belong to water tank PMTs
					Detector* det = anniegeom->ChannelToDetector(channel_key);
					if(det==nullptr){
						Log("VetoEfficiency Tool: Warning! ADC hit on channel "+to_string(channel_key)
							+" with no corresponding Detector obect in Geometry!",v_warning,verbosity);
						continue;
					}
					if (det->GetDetectorElement() != "Tank") {
						 continue;
					}
					
					// Skip non-water-tank PMTs
					if ( std::find(std::begin(water_tank_pmt_IDs), std::end(water_tank_pmt_IDs),
						 uint32_t(channel_key)) == std::end(water_tank_pmt_IDs) ) {continue;}
					
					// loop over pulses on this channel, in this minibuffer
					for (const ADCPulse& pulse : minibuffer_pulses.at(mb) ) {
						double pulse_amplitude_ = pulse.amplitude();                  // Volts
						double pulse_charge_ = pulse.charge();                        // nC
						unsigned short pulse_raw_amplitude_ = pulse.raw_amplitude();  // ADC counts
						
						int64_t pulse_start_time_ns = pulse.start_time();
						
						// ADCPulse::start_time gives the time relative to the start of the current minibuffer
						// For non-hefty mode only one minibuffer was stored per beam spill,
						// with 10us of pre-beam. 
						// So times relative to the spill are (ADCPulse::start_time - 10us).
						// 
						// For Hefty mode, (labeled by the RawLoader tool with MinibufferLabel::Hefty)
						// multiple minibuffers were recorded for each beam spill,
						// so to get the pulse time relative to the beam spill
						// we need to add the TSinceBeam value for that minibuffer.
						// 
						if ( hefty_mode_ && event_mb_label == MinibufferLabel::Hefty) {
							// (This should be zero for non-hefty events anyway,
							//  since it defaults to that value in the heftydb TTree).
							
							pulse_start_time_ns = pulse.start_time() + hefty_info.t_since_beam(mb);
						}
						
						h_adc_times->Fill(pulse_start_time_ns);
						for(int i_veto=0; i_veto < (int) veto_times.size(); i_veto++){
						  h_adc_delta_times->Fill(pulse_start_time_ns-veto_times.at(i_veto));
						  h_adc_delta_times_charge->Fill(pulse.charge(),pulse_start_time_ns-veto_times.at(i_veto));
						}		
						
						//output_pulse_tree_->Fill();
						
	//					Log("Found pulse on channel " + to_string(channel_key)
	//						+ " in run " + to_string(run_number_) + " subrun "
	//						+ to_string(subrun_number_) + " event "
	//						+ to_string(event_number_) + " in minibuffer "
	//						+ to_string(mb) + " at "
	//						+ to_string(pulse_start_time_ns) + " ns", 3, verbosity);
					} // end loop over adc hits on this channel
				} // end loop over adc channels
			} // end draw debug histos
		} // end veto layer 1 coincidence condition


		// Evaluate efficiencies for layer 2 coincidence conditions
		if (layer2_hit){
			Log("VetoEfficiency Tool: Found " + to_string(veto_l2_hits.size())
				+" hits on downstream veto layer, constructing coincidence objects",v_debug,verbosity);
			
			// since we may have multiple vetol2 hits within a close proximity,
			// we should reduce these down before we start looking for matches
			// first sort the hits by time
			std::sort(veto_l2_hits.begin(), veto_l2_hits.end());
			


			// we'll store info about each coincidence within a struct
			coincidencesl2_.clear();
			for(auto&& this_hit : veto_l2_hits){
				// search any existing coincidencesl2_ to see if this hit is close in time
				bool allocated=false;
				for(CoincidenceInfo& acoincidence : coincidencesl2_){
					int in_layer_index = std::distance(vetol2keys.begin(),
								std::find(vetol2keys.begin(), vetol2keys.end(), this_hit.second));
					if((this_hit.first-acoincidence.event_time_ns)<coincidence_tolerance_){
						if(acoincidence.vetol2hits.count(this_hit.second)){
							acoincidence.vetol2hits.at(this_hit.second).push_back(this_hit.first);
						} else {
							acoincidence.vetol2hits.emplace(this_hit.second,std::vector<double>{this_hit.first});
						}
						allocated=true;
						break;
					}
				}
				// if it wasn't within the time window of any existing coincidencesl2_, make a new one
				if(allocated==false){
					CoincidenceInfo newcoincidence;
					newcoincidence.event_time_ns = this_hit.first - pre_trigger_ns_;
					if(drawHistos) h_coincidence_event_times->Fill(newcoincidence.event_time_ns);
					newcoincidence.vetol2hits.emplace(this_hit.second, std::vector<double>{this_hit.first});
					coincidencesl2_.push_back(newcoincidence);
				}
			}
			
			Log("VetoEfficiency Tool: Reduced to " + to_string(coincidencesl2_.size())
				+ " L2 CoincidenceInfo objects, searching for tank events",v_debug,verbosity);
			
			//TODO
			for (auto&& this_hit : veto_l2_hits){
				for (int i_veto=0; i_veto<(int)veto_times.size(); i_veto++){
					h_veto_delta_times->Fill(this_hit.first-veto_times.at(i_veto));
					int in_layer_index = std::distance(vetol1keys.begin(),std::find(vetol1keys.begin(), vetol1keys.end(), veto_chankeys.at(i_veto)));	 
					int in_layer_index_layer2 = std::distance(vetol2keys.begin(),std::find(vetol2keys.begin(), vetol2keys.end(),this_hit.second));
					in_layer_index_layer2+=13;	 
					
					vector_all_adc_times.at(in_layer_index)->Fill(this_hit.first-veto_times.at(i_veto));
					if (in_layer_index_layer2<26) vector_all_adc_times.at(in_layer_index_layer2)->Fill(this_hit.first-veto_times.at(i_veto));
					h_veto_2D->Fill(in_layer_index,in_layer_index_layer2-13);
					h_tdc_chankeys_l1_l2->Fill(veto_chankeys.at(i_veto),this_hit.second);
				}
			}

			// this provides the set of reference times: now look for coincident events
			// in the rest of the detector. First calculate the tank charge and 
			// num PMTs hit within the tolerance window
			for(CoincidenceInfo& acoincidence : coincidencesl2_){
				Log("VetoEfficiency Tool: Searching for tank charge and hit PMTs from "
					+to_string(acoincidence.event_time_ns)+" to "
					+to_string(acoincidence.event_time_ns+coincidence_tolerance_),v_debug+1,verbosity);
				// calculate tank charge and # PMTs hit within a window around this hit
				int num_unique_water_pmts = BOGUS_INT;
				int num_pmts_above_thr;
				double tank_charge = compute_tank_charge(mb,
					adc_hits, acoincidence.event_time_ns-1000, acoincidence.event_time_ns + coincidence_tolerance_ - 1000,
					num_unique_water_pmts, num_pmts_above_thr);
				// XXX XXX XXX XXX based on VetoEfficiency coincidences, ADC events are 1us
				// earlier than TDC events???? Where does this come from?? XXX XXX XXX XXX
				acoincidence.num_unique_water_pmts = num_unique_water_pmts;
				//acoincidence.num_unique_water_pmts = num_pmts_above_thr;
				acoincidence.tank_charge = tank_charge;
			}
			
			Log("VetoEfficiency Tool: Searching for MRD L1 activity",v_debug,verbosity);
			
			// we might also want to use MRD activity as an alternate selection
			// for through-going events, so find any MRD hits within the windows
			// we separate into L1 and L2, but probably either is fine
			for(const unsigned long& amrd1key : mrdl1keys){
				// did it have any hits this readout
				if(TDCData->count(amrd1key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(amrd1key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						if (found_coincidence && verbosity >= v_message) std::cout <<"MRD L1 hit time "<<ahittime<<std::endl;
						for(CoincidenceInfo& acoincidence : coincidencesl2_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.mrdl1hits.count(amrd1key)){
									acoincidence.mrdl1hits.at(amrd1key).push_back(ahittime);
								} else {
									acoincidence.mrdl1hits.emplace(amrd1key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			Log("VetoEfficiency Tool: Searching for MRD L2 activity",v_debug,verbosity);
			// repeat for MRD L2 PMTs
			for(const unsigned long& amrd2key : mrdl2keys){
				// did it have any hits this readout
				if(TDCData->count(amrd2key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(amrd2key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						if (found_coincidence && verbosity >= v_message) std::cout <<"MRD L2 hit time "<<ahittime<<std::endl;
						for(CoincidenceInfo& acoincidence : coincidencesl2_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.mrdl2hits.count(amrd2key)){
									acoincidence.mrdl2hits.at(amrd2key).push_back(ahittime);
								} else {
									acoincidence.mrdl2hits.emplace(amrd2key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			
			// finally the important bit: seeing if we have any veto L1 activity.
			// whether we had any or not is going to determine our efficiency
			Log("VetoEfficiency Tool: Searching for Veto L1 activity",v_debug,verbosity);
			for(const unsigned long& avetol1key : vetol1keys){
				// did it have any hits this readout
				if(TDCData->count(avetol1key)){
					// loop over any hits in this minibuffer
					for(auto& ahit : TDCData->at(avetol1key).at(mb)){
						double ahittime = ahit.GetTime();
						// scan through our coincidence events and see if this
						// hit lies in any of their windows
						for(CoincidenceInfo& acoincidence : coincidencesl2_){
							if( (ahittime>acoincidence.event_time_ns-100) &&
								(ahittime<(acoincidence.event_time_ns-100+coincidence_tolerance_)) ){
								if(acoincidence.vetol1hits.count(avetol1key)){
									acoincidence.vetol1hits.at(avetol1key).push_back(ahittime);
								} else {
									acoincidence.vetol1hits.emplace(avetol1key,std::vector<double>{ahittime});
								}
							}
						}
					}
				}
			}
			
			// write the coincidence events to a ROOT file for easier analysis
			Log("VetoEfficiency Tool: adding coincidences to output tree",v_debug,verbosity);
			int n_coincidences_passed_cuts=0;
			for(CoincidenceInfo& acoincidence : coincidencesl2_){
				// if it passed either the tank or MRD coincidence requirements
				// for a through-going event (upstream veto requirement is always
				// satisfied by construction), then record it to file
				if( ((acoincidence.tank_charge > min_tank_charge_) &&
					 ((uint) acoincidence.num_unique_water_pmts > min_unique_water_pmts_)) ||
					((acoincidence.mrdl1hits.size() > min_mrdl1_pmts_) ||
					 (acoincidence.mrdl2hits.size() > min_mrdl2_pmts_))
					){
					n_coincidences_passed_cuts++;
					
					event_time_ns_ = acoincidence.event_time_ns;
					num_veto_l1_hits_ = acoincidence.vetol1hits.size();
					num_veto_l2_hits_ = acoincidence.vetol2hits.size();
					num_mrdl1_hits_ = acoincidence.mrdl1hits.size();
					num_mrdl2_hits_ = acoincidence.mrdl2hits.size();
					tank_charge_ = acoincidence.tank_charge;
					num_unique_water_pmts_ = acoincidence.num_unique_water_pmts;
					coincidence_info_ = acoincidence;
					coincidence_layer_ = 2;         //Second Veto layer used to define coincidences          
			
					// event cuts have been applied: we have a through-going event.
					// the efficiency measurement:
					// Build a set of the veto PMTs hit this event
					veto_l1_ids_.clear();
					for(auto&& al1hit : acoincidence.vetol1hits){
						// add this pmt if we haven't already
						if(std::find(veto_l1_ids_.begin(),veto_l1_ids_.end(),al1hit.first)==veto_l1_ids_.end()){
							veto_l1_ids_.push_back(al1hit.first);
						}
					}
					veto_l2_ids_.clear();
					for(auto&& al2hit : acoincidence.vetol2hits){
						if(std::find(veto_l2_ids_.begin(),veto_l2_ids_.end(),al2hit.first)==veto_l2_ids_.end()){
							veto_l2_ids_.push_back(al2hit.first);
						}
					}
					
					// record them in the aggregate info
					for(auto&& al2pmt : veto_l2_ids_){
						// convert to index in the layer, i.e. 0-12
						int in_layer_index = std::distance(vetol2keys.begin(),
								std::find(vetol2keys.begin(), vetol2keys.end(), al2pmt));
						// safety check
						if(in_layer_index>=(int)l2_hits_L2_.size()){
							Log("VetoEfficiency Tool: Error! Hit on veto L2 PMT index "
								+to_string(in_layer_index)+" is out of range!",v_error,verbosity);
							return false;
						}
						// note that this paddle saw an event
						l2_hits_L2_.at(in_layer_index)++;
						// check if we had a corresponding hit on the L2 PMT
						// build the coresponding L1 key so we can search for it
						unsigned long l1_key = 1000000 + (100*in_layer_index);
						// see if that key is in our list of L1 keys hit this event
						if (verbosity >= v_message) std::cout <<"In-Layer-Chankey "<<in_layer_index<<" saw a hit! Check if l1_key = "<<l1_key<<" also saw something."<<std::endl;
						for (int i_layer1=0; i_layer1<(int)veto_l1_ids_.size(); i_layer1++){
							if (verbosity >= v_message) std::cout <<"Layer 1 hit with chankey "<<veto_l1_ids_.at(i_layer1)<<std::endl;
						}
						if(std::find(veto_l1_ids_.begin(),veto_l1_ids_.end(),l1_key)!=veto_l1_ids_.end()){
							l1_hits_L2_.at(in_layer_index)++;
							if (verbosity >= v_message) std::cout <<"Layer 1 saw a hit! :)"<<std::endl;
						}
						
					}
					// ultimately the ratio of l2_hits_.at(a_pmt_num) to l1_hits_.at(a_pmt_num)
					// will give us the efficiency of a_pmt_num in L2.
					roottreeout->Fill();
				}
			}
			Log("VetoEfficiency Tool: Found "+to_string(n_coincidences_passed_cuts)
				+" through-going events",v_debug,verbosity);
			if(verbosity>v_debug){
				std::cout<<"VetoEfficiency Tool: Current efficiencies (Layer1) are: {";
				for(int i=0; i<13; ++i){
					double theefficiency = (l2_hits_L2_.at(i)==0) ? 0 : (double(l1_hits_L2_.at(i))/double(l2_hits_L2_.at(i)))*100.;
					logmessage = to_string(l1_hits_L2_.at(i)) + "/"
							   + to_string(l2_hits_L2_.at(i)) + "="
							   + to_string(theefficiency);
					if(i<12) logmessage += ", ";
					std::cout<<logmessage;
				}
				std::cout<<"}"<<std::endl;
			}
		} //end loop over L2 coincidences
		
	} // end loop over minibuffers in this readout
	Log("VetoEfficiency Tool: Execute done",v_debug,verbosity);
	
	return true;
}


bool VetoEfficiency::Finalise(){
	
	
	// calculate final efficiencies
	std::cout<<"VetoEfficiency Tool: Final efficiencies (Layer 1) are: {";
	for(int i=0; i<13; ++i){
		double theefficiency = (l2_hits_L2_.at(i)==0) ? 0 : (double(l1_hits_L2_.at(i))/double(l2_hits_L2_.at(i)))*100.;
		l1_efficiencies_.at(i)=theefficiency;
		logmessage = to_string(l1_hits_L2_.at(i)) + "/"
				   + to_string(l2_hits_L2_.at(i)) + "="
				   + to_string(theefficiency);
		if(i<12) logmessage += ", ";
		std::cout<<logmessage;
	}
	std::cout<<"}"<<std::endl;
	std::cout<<"VetoEfficiency Tool: Final efficiencies (Layer 2) are: {";
	for(int i=0; i<13; ++i){
		double theefficiency = (l1_hits_L1_.at(i)==0) ? 0 : (double(l2_hits_L1_.at(i))/double(l1_hits_L1_.at(i)))*100.;
		l2_efficiencies_.at(i)=theefficiency;
		logmessage = to_string(l2_hits_L1_.at(i)) + "/"
				   + to_string(l1_hits_L1_.at(i)) + "="
				   + to_string(theefficiency);
		if(i<12) logmessage += ", ";
		std::cout<<logmessage;
	}
	std::cout<<"}"<<std::endl;
	
	// write the aggregate info to the tree
	summarytreeout->Fill();
	
	// write the trees to file
	rootfileout->cd();
	roottreeout->Write();
	summarytreeout->Write();
	
	// cleanup: reset branches and close the file
	roottreeout->ResetBranchAddresses();
	summarytreeout->ResetBranchAddresses();
	rootfileout->Close();
	
	// i can never remember if we need to delete these or not
//	delete roottreeout; roottreeout=nullptr;
//	delete summarytreeout; summarytreeout=nullptr;
	
	// delete the file
	delete rootfileout; rootfileout=nullptr;
	
	if (drawHistos){
		f_veto->cd();
		h_tdc_chankeys_l1->Write();
		h_tdc_chankeys_l2->Write();
		h_tdc_chankeys_l1_l2->Write();
		h_all_veto_times->Write();
		h_all_veto_times_layer2->Write();
		h_veto_2D->Write();
		h_all_mrd_times->Write();
		h_all_adc_times->Write();
		h_tdc_times->Write();
		h_adc_times->Write();
		h_tdc_delta_times->Write();
		h_tdc_delta_times_L1->Write();
		h_tdc_delta_times_L2->Write();
		h_adc_delta_times->Write();
		h_adc_delta_times_charge->Write();
		h_veto_delta_times->Write();
		h_veto_delta_times_coinc->Write();
		h_adc_charge->Write();
		h_adc_charge_coinc->Write();
		for (int i=0;i<26;i++){
			vector_all_adc_times.at(i)->Write();
		}
		f_veto->Close();
		delete f_veto;
	}

	if(drawHistos){
		double canvwidth = 700;
		double canvheight = 600;
		plotCanv = new TCanvas("plotCanv","",canvwidth,canvheight);
		
		plotCanv->cd();
		std::string imgname;
		
		h_tdc_times->Draw();
		imgname=h_tdc_times->GetTitle();
		std::replace(imgname.begin(), imgname.end(), ' ', '_');
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
		
		h_adc_times->Draw();
		imgname=h_adc_times->GetTitle();
		std::replace(imgname.begin(), imgname.end(), ' ', '_');
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
		
		if(h_tdc_times->GetBinContent(h_tdc_times->GetMaximumBin())>h_adc_times->GetBinContent(h_adc_times->GetMaximumBin())){
			h_tdc_times->Draw();
			h_adc_times->Draw("same");
		} else {
			h_adc_times->Draw();
			h_tdc_times->Draw("same");
		}
		//plotCanv->BuildLegend();
		imgname="TDC_ADC_times";
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
		
//		plotCanv->Clear();
//		plotCanv->Modified();
//		plotCanv->Update();
//		gSystem->ProcessEvents();
                h_tdc_times->Scale(1./h_tdc_times->Integral());
                h_adc_times->Scale(1./h_adc_times->Integral());
                if(h_tdc_times->GetBinContent(h_tdc_times->GetMaximumBin())>h_adc_times->GetBinContent(h_adc_times->GetMaximumBin())){
			h_tdc_times->Draw();
			h_adc_times->Draw("same");
		} else {
			h_adc_times->Draw();
			h_tdc_times->Draw("same");
		}
//		plotCanv->Modified();
//		plotCanv->Update();
//		gSystem->ProcessEvents();
		imgname="TDC_ADC_times_scaled";
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
		
		h_coincidence_event_times->Draw();
		imgname=h_coincidence_event_times->GetTitle();
		std::replace(imgname.begin(), imgname.end(), ' ', '_');
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
		
		h_tankhit_time_to_coincidence->Draw();
		imgname=h_tankhit_time_to_coincidence->GetTitle();
		std::replace(imgname.begin(), imgname.end(), ' ', '_');
		plotCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	
		
		// cleanup
		// ~~~~~~~
		Log("VetoEfficiency Tool: deleting reco histograms",v_debug,verbosity);
		if(h_adc_times) delete h_adc_times; h_adc_times=nullptr;
		if(h_tdc_times) delete h_tdc_times; h_tdc_times=nullptr;
		if(h_tankhit_time_to_coincidence) delete h_tankhit_time_to_coincidence; h_tankhit_time_to_coincidence=nullptr;
		if(h_coincidence_event_times) delete h_coincidence_event_times; h_coincidence_event_times=nullptr;
		
		Log("VetoEfficiency Tool: deleting canvas",v_debug,verbosity);
		if(gROOT->FindObject("plotCanv")){ 
			delete plotCanv; plotCanv=nullptr;
		}
		
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				Log("VetoEfficiency Tool: deleting gloabl TApplication",v_debug,verbosity);
				delete rootTApp;
			}
		} else if (tapplicationusers>1){
			tapplicationusers--;
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		}
	}
	
	return true;
}


// function adapted from PhaseITreeMaker tool...
// Returns the integrated tank charge in a given time window.
// Also loads the integer num_unique_water_pmts with the number
// of hit water tank PMTs.
double VetoEfficiency::compute_tank_charge(size_t minibuffer_number,
	const std::map< unsigned long, std::vector< std::vector<ADCPulse> > >& adc_hits,
	uint64_t start_time, uint64_t end_time, int& num_unique_water_pmts, int& num_pmts_above){
	double tank_charge = 0.;
	num_unique_water_pmts = 0;
	num_pmts_above = 0;	

	for (const auto& channel_pair : adc_hits) {
		const unsigned long& ck = channel_pair.first;
		
		// Only consider ADC channels that belong to water tank PMTs
		Detector* det = anniegeom->ChannelToDetector(ck);
		if(det==nullptr){
			Log("VetoEfficiency Tool: Warning! ADC hit on channel "+to_string(ck)
				+" with no corresponding Detector obect in Geometry!",v_warning,verbosity);
			continue;
		}
		if (det->GetDetectorElement() != "Tank") continue;
		
		// Skip non-water-tank PMTs
		if ( std::find(std::begin(water_tank_pmt_IDs), std::end(water_tank_pmt_IDs),
			 uint32_t(ck)) == std::end(water_tank_pmt_IDs) ) continue;
		
		// Get the vector of pulses for the current channel and minibuffer of
		// interest
		const auto& pulse_vec = channel_pair.second.at(minibuffer_number);
		
		bool found_pulse_in_time_window = false;
		
		for ( const auto& pulse : pulse_vec ) {
			size_t pulse_time = pulse.start_time();
			
			if (found_coincidence) //std::cout <<"ADC pulse time: "<<pulse_time<<", pulse charge: "<<pulse.charge()<<std::endl;
			if(drawHistos){
				int64_t tdiff =   (static_cast<int64_t>(pulse_time)
								- (static_cast<int64_t>(start_time)+static_cast<int64_t>(pre_trigger_ns_)));
				//bool in_window = ( pulse_time >= start_time && pulse_time <= end_time );
				// note the time of this pulse within the coincidence event window
				// take off pretrigger so we're plotting time relative to the first L1 veto hit
				h_tankhit_time_to_coincidence->Fill(tdiff);
			}
			if ( pulse_time >= start_time && pulse_time <= end_time) {
				// Increment the unique hit PMT counter if the current pulse is within
				// the given time window (and we hadn't already)
				if (!found_pulse_in_time_window) {
					++num_unique_water_pmts;
					found_pulse_in_time_window = true;
					if (pulse.charge()>0.05) ++num_pmts_above;
				}
				tank_charge += pulse.charge();
			}
		}
	}
	
	Log("Tank charge = " + to_string(tank_charge) + " nC", 4, verbosity);
	Log("Unique PMTs = " + to_string(num_unique_water_pmts), 4, verbosity);
	//std::cout <<"Num PMTs above threshold (0.01nC): "<<num_pmts_above<<std::endl;
	return tank_charge;
}

void VetoEfficiency::LoadTDCKeys(){
	// channelkeys for TDC channels are representative of their location in the MRD system
	// so to find hits on the upstream veto layer we need to generate the corresponding
	// channelkeys
	for(std::string& tdc_string : phase1mrdpmts){
		// strings are of the form "1000102", or "1XXYYZZ"
		// the preceding 1 is required to make channelkeys unique from ADC ones
		int tubeid = stoi(tdc_string);
		int corrected_tubeid = tubeid-1000000;
		int mrdpmtxnum=floor(corrected_tubeid/10000);
		int mrdpmtznum=corrected_tubeid-(100*floor(corrected_tubeid/100));
		int mrdpmtynum=(corrected_tubeid-(10000*mrdpmtxnum)-mrdpmtznum)/100;
		unsigned long tdc_key = static_cast<unsigned long>(tubeid);
		if(mrdpmtznum==0){
			// for no good reason the two veto layers
			// are considered different X layers,
			// rather than two different Z layers
			if(mrdpmtxnum==0){
				vetol1keys.push_back(tdc_key);
			} else {
				vetol2keys.push_back(tdc_key);
			}
		// making this more silly is the fact that
		// the first MRD layer is then layer 2,
		// and there simply isn't a layer 1.
		} else if(mrdpmtznum==2){
			mrdl1keys.push_back(tdc_key);
			allmrdkeys.push_back(tdc_key);
		} else if(mrdpmtznum==3){
			mrdl2keys.push_back(tdc_key);
			allmrdkeys.push_back(tdc_key);
		} else {
			Log("Unknown MRD z number: "+tdc_key,v_error,verbosity);
		}
	}
	std::cout <<"Loaded the following veto-l1 keys:"<<std::endl;
	for (int i=0; i<(int)vetol1keys.size(); i++){
		std::cout <<vetol1keys.at(i)<<std::endl;
	}
	for (int i=0; i<(int)vetol2keys.size(); i++){
		std::cout <<vetol2keys.at(i)<<std::endl;
	}
}

void VetoEfficiency::makeOutputFile(std::string outputfilename){
	Log("VetoEfficiency Tool: Making output file "+outputfilename,v_debug,verbosity);
	rootfileout = new TFile(outputfilename.c_str(),"RECREATE");
	rootfileout->cd();
	
	roottreeout = new TTree("veto_efficiency","PhaseI Veto Coincidences");
	// event information
	roottreeout->Branch("hefty_mode",&hefty_mode_);
	roottreeout->Branch("hefty_trigger_mask", &hefty_trigger_mask_);
	roottreeout->Branch("event_label", &event_label_);
	roottreeout->Branch("run_number", &run_number_);
	roottreeout->Branch("subrun_number",&subrun_number_);
	roottreeout->Branch("event_number",&event_number_);
	roottreeout->Branch("minibuffer_number",&minibuffer_number_);
	roottreeout->Branch("spill_number_",&spill_number_);
	// coincidence information
	roottreeout->Branch("veto_l1_ids",&veto_l1_ids_);
	roottreeout->Branch("veto_l2_ids",&veto_l2_ids_);
	roottreeout->Branch("event_time_ns",&event_time_ns_);
	roottreeout->Branch("num_veto_l1_hits",&num_veto_l1_hits_);
	roottreeout->Branch("num_veto_l2_hits",&num_veto_l2_hits_);
	roottreeout->Branch("num_mrdl1_hits",&num_mrdl1_hits_);
	roottreeout->Branch("num_mrdl2_hits",&num_mrdl2_hits_);
	roottreeout->Branch("tank_charge",&tank_charge_);
	roottreeout->Branch("num_unique_water_pmts",&num_unique_water_pmts_);
	roottreeout->Branch("coincidence_layer",&coincidence_layer_);	//Indicate whether the coincidence condition was done with Veto Layer 1 or Veto Layer 2
//	roottreeout->Branch("coincidence_info",&coincidence_info_pt_);  // further info?
	
	summarytreeout = new TTree("summary_tree","Aggregate Info");
	summarytreeout->Branch("total_POT",&total_POT_);
	summarytreeout->Branch("num_beam_spills",&spill_number_);
	summarytreeout->Branch("num_hits_on_l1_PMTs",&l1_hits_L1_);
	summarytreeout->Branch("coincident_hits_on_l2_PMTs",&l2_hits_L1_);
	summarytreeout->Branch("l2efficiencies",&l2_efficiencies_);
	summarytreeout->Branch("num_hits_on_l2_PMTs",&l2_hits_L2_);
	summarytreeout->Branch("coincident_hits_on_l1_PMTs",&l1_hits_L2_);
	summarytreeout->Branch("l1efficiencies",&l1_efficiencies_);
	
//	summarytreeout->Branch("num_source_triggers",&num_source_triggers_);
//	summarytreeout->Branch("num_cosmic_triggers",&num_cosmic_triggers_);
//	summarytreeout->Branch("num_soft_triggers",&num_soft_triggers_);
	
	gROOT->cd();
}

