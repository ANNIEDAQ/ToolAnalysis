/* vim:set noexpandtab tabstop=4 wrap */
#include "SimulatedWaveformDemo.h"

SimulatedWaveformDemo::SimulatedWaveformDemo():Tool(){}

bool SimulatedWaveformDemo::Initialise(std::string configfile, DataModel &data){
	
	std::cout<<"SimulatedWaveformDemo Tool: Initializing"<<std::endl;
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("WaveformSource",WaveformSource);  // ANNIEEvent::RawADCData=0, pmtDataStore::AllData=1
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// Get the Simulated data BoostStores
	Log("SimulatedWaveformDemo Tool: Getting BoostStores",v_debug,verbosity);
	m_data->Stores.at("SimulatedFileStore")->Get("pmtDataStore",pmtDataStore);
	m_data->Stores.at("SimulatedFileStore")->Get("heftydbStore",heftydbStore);
	m_data->Stores.at("SimulatedFileStore")->Get("CCDataStore",CCDataStore);
	
	// Get the static members from the pmtDataStore
	Log("SimulatedWaveformDemo Tool: Getting pmtDataStore header vals",v_debug,verbosity);
	pmtDataStore->Get("RunStartTimeSec",RunStartTimeSec);
	pmtDataStore->Get("MiniBuffersPerFullBuffer",MiniBuffersPerFullBuffer);
	pmtDataStore->Get("ChannelsPerAdcCard",ChannelsPerAdcCard);
	pmtDataStore->Get("BufferSize",BufferSize);
	SamplesPerMinibuffer = BufferSize/MiniBuffersPerFullBuffer;
	
	Log("SimulatedWaveformDemo Tool: Getting data pointers",v_debug,verbosity);
	// Waveforms are accessible from the pmtDataStore via the 'pmtDataVector' member.
	// This is a vector of pointers to vectors: std::vector< std::vector<uint16_t>* >
	// The vectors pointed to are re-used between events, so these pointers do not change;
	// only the data in the vectors they point at changes.
	// We therefore only need to retrieve the pointers once.
	std::vector<intptr_t> pmtDataVecCopy;
	bool get_ok = pmtDataStore->Get("AllData",pmtDataVecCopy);
	if(not get_ok){
		Log("SimulatedWaveformDemo Tool: Error getting AllData from pmtDataStore!",v_error,verbosity);
		return false;
	}
	for(int i=0; i<(int)pmtDataVecCopy.size(); i++){
		intptr_t ddatapi = pmtDataVecCopy.at(i);
		const std::vector<uint16_t>* ddatap = reinterpret_cast<const std::vector<uint16_t>*>(ddatapi);
		pmtDataVector.push_back(ddatap);
	}
	
	// Similarly, the heftydbStore only holds pointers to arrays on the heap,
	// so we need only retreve them once; the arrays are re-used, so the pointers do not change.
	heftydbStore->Get("minibuffer_starts_unixns",minibuffer_starts_unixns);
	heftydbStore->Get("minibuffer_label",minibuffer_label);
	heftydbStore->Get("minibuffer_tsincebeam",minibuffer_tsincebeam);
	heftydbStore->Get("dropped_delayed_triggers",dropped_delayed_triggers);
	
	// and same for CCDataStore vectors
	CCDataStore->Get("CardType",CardType);
	CCDataStore->Get("Slot",Slot);
	CCDataStore->Get("Channel",Channel);
	CCDataStore->Get("TDCTickValue",TDCTickValue);
	
	// make a numberline to act as the x-axis of the plotting TGraph
	numberline.resize(SamplesPerMinibuffer);
	std::iota(numberline.begin(),numberline.end(),0);
	upcastdata.resize(SamplesPerMinibuffer);
	
	Log("SimulatedWaveformDemo Tool: Acquiring ROOT TApplication",v_debug,verbosity);
	// make the ROOT TApplication so we can display waveforms
	// ======================================================
	// There may only be one TApplication, so if another tool has already made one
	// register ourself as a user. Otherwise, make one and put a pointer in the CStore for other Tools
	int myargc=0;
	//char *myargv[] = {(const char*)"simulatedwfrmdemo"};
	// get or make the TApplication
	intptr_t tapp_ptr=0;
	get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
	if(not get_ok){
		Log("SimulatedWaveformDemo Tool: Making global TApplication",v_error,verbosity);
		rootTApp = new TApplication("rootTApp",&myargc,0);
		tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
		m_data->CStore.Set("RootTApplication",tapp_ptr);
	} else {
		Log("SimulatedWaveformDemo Tool: Retrieving global TApplication",v_error,verbosity);
		rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
	}
	int tapplicationusers;
	get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
	if(not get_ok) tapplicationusers=1;
	else tapplicationusers++;
	m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	canvwidth = 700;
	canvheight = 600;
	
	return true;
}


bool SimulatedWaveformDemo::Execute(){
	
	Log("SimulatedWaveformDemo Tool: Executing",v_debug,verbosity);
	// As described in Initialise, we only need to `Get` basic members that will update between loops
	// =================================================================================================
	
	// The CCDataStore contains one entry per Camac readout.
	// There is one Camac readout for every trigger (minibuffer)
	// -------------------------------------------------------------------------------------------------
	CCDataStore->Get("CCDataReadout",CCDataReadout);
	CCDataStore->Get("CCHitsThisReadout",CCHitsThisReadout);
	CCDataStore->Get("ReadoutTimeUnixMs",ReadoutTimeUnixMs);
	
	// =================================================================================================
	
	// The pmtDataStore contains one entry per ADC readout.
	// There is one ADC readout for each *full buffer* readout, and each full buffer may
	// contain multiple triggers (minibuffers).
	// For the simplest case, configure the PulseSimulation tool with MinibuffersPerFullbuffer = 1.
	// ---------------------------------------------------------------------------------------------------
	pmtDataStore->Get("SequenceID",pmtDataSequenceID);
	
	// =================================================================================================
	
	// The heftydbStore also contains one entry per ADC readout.
	// Each variable in the Store is a c-style array of size MinibuffersPerFullbuffer.
	// -------------------------------------------------------------------------------------------------
	heftydbStore->Get("SequenceID",heftydbSequenceID);
	
	// =================================================================================================
	
	// < NEW >
	// Phase 2 oriented ToolAnalysis approach to analysing waveforms is via ANNIEEvent::RawADCData
	// This is a map<channelkey, vector<Waveform>>, where the vector contains one Waveform per minibuffer.
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("RawADCData",RawADCData);
	if(not get_ok){
		Log("SimulatedWaveformDemo Tool: Failed to get RawADCData from ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	// WORKING WITH THE WAVEFORMS
	// --------------------------
	
	// =================================================================================================
	// METHOD 1
	// =================================================================================================
	if(WaveformSource==0){
		// Method 1: Accessing waveforms via the RawADCData
		Log("SimulatedWaveformDemo Tool: Looping over "+to_string(RawADCData.size())
			 +" Tank PMT channels",v_debug,verbosity);
		for(std::pair<const unsigned long,std::vector<Waveform<uint16_t>>>& achannel : RawADCData){
			const unsigned long channelkey = achannel.first;
			
			// Each Waveform represents one minibuffer on this channel
			Log("SimulatedWaveformDemo Tool: Looping over "+to_string(achannel.second.size())
				 +" minibuffers",v_debug,verbosity);
			for(Waveform<uint16_t>& wfrm : achannel.second){
				std::vector<uint16_t>* samples = wfrm.GetSamples();
				// Note that because these are built directly from the simulated minibuffers,
				// samples will have the Phase 1 interleaving unless it is disabled in the
				// PulseSimulation tool by passing config variable DoPhaseOneRiffle=0
				
				// As a demonstration, let's just plot the waveform on a TGraph
				// To avoid drawing lots of boring plots with no pulses in, only
				// draw the minibuffer if it has a larger pulse than any we've seen before!
				const uint16_t thismax = *std::max_element(samples->begin(), samples->end());
				Log("SimulatedWaveformDemo Tool: checking max "+to_string(thismax)+" against "
					+to_string(maxwfrmamp),v_debug,verbosity);
				if(thismax<maxwfrmamp){ continue; }
				maxwfrmamp = thismax;
				
				// for plotting on a TGraph we need to up-cast the data from uint16_t to int32_t
				Log("SimulatedWaveformDemo Tool: Making TGraph",v_debug,verbosity);
				for(int samplei=0; samplei<SamplesPerMinibuffer; samplei++){
					upcastdata.at(samplei) = samples->at(samplei);
				}
				if(mb_graph){ delete mb_graph; }
				mb_graph = new TGraph(SamplesPerMinibuffer, numberline.data(), upcastdata.data());
				mb_graph->SetName("mb_graph");
				if(gROOT->FindObject("mb_canv")==nullptr) mb_canv = new TCanvas("mb_canv");
				mb_canv->cd();
				mb_canv->Clear();
				mb_graph->Draw("alp");
				mb_canv->Modified();
				mb_canv->Update();
				gSystem->ProcessEvents();
				do{
					gSystem->ProcessEvents();
					std::this_thread::sleep_for(std::chrono::milliseconds(100));
				} while (gROOT->FindObject("mb_canv")!=nullptr); // wait until user closes canvas
				Log("SimulatedWaveformDemo Tool: graph closed, looping",v_debug,verbosity);
			} // end loop over minibuffers
		} // end loop over channelkeys
		
	
	// =================================================================================================
	// METHOD 2
	// =================================================================================================
	} else {
		// Method 2: Accessing waveforms via AllData.
		// This was the initial conversion, based on putting a simplified phase 1 file format into BoostStores.
		// The pmtDataStore AllData vector contains one waveform entry per ADC card.
		Log("SimulatedWaveformDemo Tool: Looping over "+to_string(pmtDataVector.size())
			 +" ADC cards",v_debug,verbosity);
		for(int cardi=0; cardi<(int)pmtDataVector.size(); cardi++){
			Log("SimulatedWaveformDemo Tool: Looping over "+to_string(ChannelsPerAdcCard)
				 +" ADC channels",v_debug,verbosity);
			const std::vector<uint16_t>* card_data = pmtDataVector.at(cardi);
			if(card_data->size()==0){
				Log("SimulatedWaveformDemo Tool: Card data vector is empty!",v_error,verbosity);
				return false;
			}
			// Each waveform entry is a concatenation of waveforms for each channels on that card:
			// [Chan 1][Chan 2][Chan 3][Chan 4]
			for(int chani=0; chani<ChannelsPerAdcCard; chani++){
				int chan_start_index = BufferSize*chani;
				// Furthermore, within each cards subset, all minibuffers for that full readout are concatenated:
				// [Chan 1] == [{chan 1, minibuffer 1}{chan 1, minibuffer 2}...{chan 1, minibuffer N}]
				Log("SimulatedWaveformDemo Tool: Looping over "+to_string(MiniBuffersPerFullBuffer)
					 +" minibuffers",v_debug,verbosity);
				for(int minibufi=0; minibufi<MiniBuffersPerFullBuffer; minibufi++){
					int minibuffer_start_index = chan_start_index + (SamplesPerMinibuffer*minibufi);
					
					// Note that in phase 1 the samples within each waveform were not in order;
					// samples from the first and second half of the waveform were interleaved.
					// To disable this "shuffling" of samples, set DoPhaseOneRiffle = false
					// in the PulseSimulation tool.
					
					// As a demonstration, let's just plot the minibuffer waveform
					// To avoid drawing lots of boring plots with no pulses in, only
					// draw the minibuffer if it has a larger pulse than any we've seen before!
					std::vector<uint16_t>::const_iterator startit = card_data->begin();
					std::advance(startit,minibuffer_start_index);
					std::vector<uint16_t>::const_iterator stopit = startit;
					std::advance(stopit,SamplesPerMinibuffer);
					const uint16_t thismax = *std::max_element(startit, stopit);
					Log("SimulatedWaveformDemo Tool: checking max "+to_string(thismax)+" against "
						+to_string(maxwfrmamp),v_debug,verbosity);
					if(thismax<maxwfrmamp){ continue; }
					maxwfrmamp = thismax;
					
					// for plotting on a TGraph we need to up-cast the data from uint16_t to int32_t
					Log("SimulatedWaveformDemo Tool: Making TGraph",v_debug,verbosity);
					for(int samplei=0; samplei<SamplesPerMinibuffer; samplei++){
						upcastdata.at(samplei) = card_data->at(minibuffer_start_index+samplei);
					}
					if(mb_graph){ delete mb_graph; }
					mb_graph = new TGraph(SamplesPerMinibuffer, numberline.data(), upcastdata.data());
					if(gROOT->FindObject("mb_canv")==nullptr) mb_canv = new TCanvas("mb_canv");
					mb_canv->cd();
					mb_canv->Clear();
					mb_graph->Draw("alp");
					mb_canv->Modified();
					mb_canv->Update();
					gSystem->ProcessEvents();
					do{
						gSystem->ProcessEvents();
						std::this_thread::sleep_for(std::chrono::milliseconds(100));
					} while (gROOT->FindObject("mb_canv")!=nullptr); // wait until user closes canvas
					Log("SimulatedWaveformDemo Tool: graph closed, looping",v_debug,verbosity);
				} // end loop over minibuffers
			} // end loop over channels on this ADC card
		} // end loop over ADC cards
	}
	
	Log("SimulatedWaveformDemo Tool: finished Execute",v_debug,verbosity);
	return true;
}


bool SimulatedWaveformDemo::Finalise(){
	
	Log("SimulatedWaveformDemo Tool: Finalising",v_debug,verbosity);
	
	if(mb_graph){ delete mb_graph; mb_graph=nullptr; }
	if(gROOT->FindObject("mb_canv")!=nullptr){ delete mb_canv; mb_canv=nullptr; }
	
	// see if we're the last user of the TApplication:
	// release it if so, otherwise de-register us as a user since we're done
	int tapplicationusers=0;
	get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
	if(not get_ok || tapplicationusers==1){
		if(rootTApp){
			std::cout<<"SimulatedWaveformDemo Tool: Deleting global TApplication"<<std::endl;
			delete rootTApp;
			rootTApp=nullptr;
		}
	} else if(tapplicationusers>1){
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers-1);
	}
	
	Log("SimulatedWaveformDemo Tool: Finalise done",v_debug,verbosity);
	return true;
}
