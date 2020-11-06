#include <PulseSimulation.h>
/* vim:set noexpandtab tabstop=4 wrap */

// CREATE+OPEN RAW FORMAT OUTPUT FILE
// ==================================
void PulseSimulation::LoadOutputFiles(){
	
	Log("PulseSimulation Tool: Constructing output files",v_debug,verbosity);
	// Use the name of the simulation file to generate the name of the fake raw file
	std::string wcsimfilename;
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCFile", wcsimfilename); // wcsim_0.AAA.B.root
	
	int rawfilerun, rawfilesubrun, rawfilepart;
	// filename is of the form "wcsim_0.AAA.B.root"
	int numargs = 0;
	try{
		std::reverse(wcsimfilename.begin(),wcsimfilename.end());
		// Get the address
		char * p = std::strtok(const_cast<char*>(wcsimfilename.c_str()),".-_");
		while(numargs <3){
			p = std::strtok(NULL,".-_");
			std::string capturedstring(p);
			std::reverse(capturedstring.begin(),capturedstring.end());
			if (numargs == 0) rawfilepart = std::stoi(capturedstring.c_str());
			else if (numargs == 1) rawfilesubrun = std::stoi(capturedstring.c_str());
			else if (numargs == 2) rawfilerun = std::stoi(capturedstring.c_str());
			numargs++;
		}
		std::reverse(wcsimfilename.begin(),wcsimfilename.end());
	} catch (int e){
		cerr<<"unrecognised input file pattern: "<<wcsimfilename
			<<", will set rawfilerun=0, rawfilesubrun=0, rawfilepart=0"<<endl;
		//return;
		rawfilerun=0;
		rawfilesubrun=0;
		rawfilepart = 0;
	}	
	std::string rawfilename="RAWDataR"+to_string(rawfilerun)+
		"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+".root";
	
	// Get the necessary config information on buffer sizes etc
	m_data->CStore.Get("WCSimPreTriggerWindow",pre_trigger_window_ns);
	m_data->CStore.Get("WCSimPostTriggerWindow",post_trigger_window_ns);
	minibuffer_datapoints_per_channel = (abs(pre_trigger_window_ns)+post_trigger_window_ns) / ADC_NS_PER_SAMPLE;
	buffer_size = minibuffer_datapoints_per_channel * minibuffers_per_fullbuffer;
	full_buffer_size = buffer_size * channels_per_adc_card;
	emulated_event_size = minibuffer_datapoints_per_channel / 4.; // the 4 comes from ADC firmware stuff
//	cout<<"pretriggerwindow="<<pre_trigger_window_ns
//		<<", posttriggerwindow="<<post_trigger_window_ns<<endl;
//	cout<<"minibuffer_datapoints_per_channel="<<minibuffer_datapoints_per_channel<<endl;
//	cout<<"buffersize="<<buffer_size<<endl;
//	cout<<"fullbuffersize="<<full_buffer_size<<endl;
//	cout<<"minibuffers_per_fullbuffer="<<minibuffers_per_fullbuffer<<endl;
	
	// Make the static arrays for holding data.
	//*----------------------------------------------------------------------------*
	fileout_TriggerCounts    = new ULong64_t[minibuffers_per_fullbuffer];     // must be >= TriggerNumber
	fileout_Rates            = new UInt_t[channels_per_adc_card];             // must be >= Channels
	fileout_Data             = new UShort_t[full_buffer_size];                // must be >= FullBufferSize
	//*............................................................................*
	// These only store filler data, but the RawLoader tool needs them to exist and have the right sizes
	fileout_EventIDs          = new UShort_t[MAXEVENTSIZE];    // must be >= Eventsize
	fileout_EventTimes        = new ULong64_t[MAXEVENTSIZE];   // must be >= Eventsize
	fileout_TriggerMasks      = new UInt_t[MAXTRIGGERSIZE];    // must be >= TriggerSize
	fileout_TriggerCounters   = new UInt_t[MAXTRIGGERSIZE];    // must be >= TriggerSize
	//*----------------------------------------------------------------------------*
	// Create the file, TTrees, branches etc
	if(GenerateFakeRootFiles){
		Log("PulseSimulation Tool: Creating fake output file " + rawfilename,v_debug,verbosity);
		rawfileout = new TFile(rawfilename.c_str(),"RECREATE");
		//*----------------------------------------------------------------------------*
		tPMTData                 = new TTree("PMTData","");                       // one entry per trigger readout
		TBranch *bLastSync       = tPMTData->Branch("LastSync", &fileout_LastSync);
		TBranch *bSequenceID     = tPMTData->Branch("SequenceID", &fileout_SequenceID);
		TBranch *bStartTimeSec   = tPMTData->Branch("StartTimeSec", &fileout_StartTimeSec);
		TBranch *bStartTimeNSec  = tPMTData->Branch("StartTimeNSec", &fileout_StartTimeNSec);
		TBranch *bStartCount     = tPMTData->Branch("StartCount", &fileout_StartCount);
		TBranch *bTriggerNumber  = tPMTData->Branch("TriggerNumber", &fileout_TriggerNumber);
		TBranch *bTriggerCounts  = tPMTData->Branch("TriggerCounts", &fileout_TriggerCounts,
			"TriggerCounts[TriggerNumber]/l");
		TBranch *bCardID         = tPMTData->Branch("CardID", &fileout_CardID);
		TBranch *bChannels       = tPMTData->Branch("Channels", &fileout_Channels);
		TBranch *bRates          = tPMTData->Branch("Rates", &fileout_Rates, 
			"Rates[Channels]/i");
		TBranch *bBufferSize     = tPMTData->Branch("BufferSize", &fileout_BufferSize);
		TBranch *bEventsize      = tPMTData->Branch("Eventsize", &fileout_Eventsize);
		TBranch *bFullBufferSize = tPMTData->Branch("FullBufferSize", &fileout_FullBufferSize);
		TBranch *bData           = tPMTData->Branch("Data", &fileout_Data, "Data[FullBufferSize]/s");
		//*----------------------------------------------------------------------------*
		tRunInformation = new TTree("RunInformation","");
		TBranch *bInfoTitle = tRunInformation->Branch("InfoTitle", &fileout_InfoTitle);
		TBranch *bInfoMessage = tRunInformation->Branch("InfoMessage", &fileout_InfoMessage);
		//*............................................................................*
		tTrigData                 = new TTree("TrigData","");
		TBranch *bFirmwareVersion = tTrigData->Branch("FirmwareVersion", &fileout_FirmwareVersion);
		TBranch *bSequenceID2     = tTrigData->Branch("SequenceID", &fileout_SequenceID);
		TBranch *bEventsize2      = tTrigData->Branch("EventSize", &fileout_Eventsize);
		TBranch *bTriggerSize     = tTrigData->Branch("TriggerSize", &fileout_TriggerSize);
		TBranch *bFIFOOverflow    = tTrigData->Branch("FIFOOverflow", &fileout_FIFOOverflow);
		TBranch *bDriverOverfow   = tTrigData->Branch("DriverOverfow", &fileout_DriverOverfow);
		TBranch *bEventIDs        = tTrigData->Branch("EventIDs", &fileout_EventIDs,
			"EventIDs[EventSize]/s");
		TBranch *bEventTimes      = tTrigData->Branch("EventTimes", &fileout_EventTimes,
			"EventTimes[EventSize]/l");
		TBranch *bTriggerMasks    = tTrigData->Branch("TriggerMasks", &fileout_TriggerMasks,
			"TriggerMasks[TriggerSize]/I");
		TBranch *bTriggerCounters = tTrigData->Branch("TriggerCounters", &fileout_TriggerCounters, "TriggerCounters[TriggerSize]/I");
		//*----------------------------------------------------------------------------*
		tCCData = new TTree("CCData","");
		TBranch *bTrigger   = tCCData->Branch("Trigger", &fileout_Trigger);         // readout number
		TBranch *bOutNumber = tCCData->Branch("OutNumber", &fileout_OutNumber);     // num hits this event/readout
		TBranch *bTimeStamp = tCCData->Branch("TimeStamp", &fileout_TimeStamp);     // readout time, unix ms
		// branches need to store the address of a pointer to the vector: a pointer to a pointer
		// we should reserve some space as inducing a storage reallocation would invalidate our pointers
		fileout_Type.reserve(500);   // not likely to get > 500 MRD hits in a single event
		fileout_Value.reserve(500);
		fileout_Slot.reserve(500);
		fileout_Channel.reserve(500);
		pfileout_Type       = &fileout_Type;
		pfileout_Value      = &fileout_Value;
		pfileout_Slot       = &fileout_Slot;
		pfileout_Channel    = &fileout_Channel;
		// set branch addresses
		TBranch *bType      = tCCData->Branch("Type", &pfileout_Type);              // card type string, e.g. "ADC"
		TBranch *bValue     = tCCData->Branch("Value", &pfileout_Value);            // see below
		TBranch *bSlot      = tCCData->Branch("Slot", &pfileout_Slot);              // card position in crate
		TBranch *bChannel   = tCCData->Branch("Channel", &pfileout_Channel);        // channel in card
	}
	
	//*----------------------------------------------------------------------------*
	// Since we can't properly synthesize all the timing variables
	// we'll directly create the simplified output Hefty Timing file
	std::string timingfilename="DataR"+to_string(rawfilerun)+
		"S"+to_string(rawfilesubrun)+"p"+to_string(rawfilepart)+"_timing.root";
	//*----------------------------------------------------------------------------*
	timefileout_Time          = new ULong_t[minibuffers_per_fullbuffer];
	timefileout_Label         = new Int_t[minibuffers_per_fullbuffer];
	timefileout_TSinceBeam    = new Long_t[minibuffers_per_fullbuffer];
	timefileout_More          = new Int_t[minibuffers_per_fullbuffer];
	//*----------------------------------------------------------------------------*
	if(GenerateFakeRootFiles){
		Log("PulseSimulation Tool: Creating RAW timing file "+timingfilename,v_warning,verbosity);
		timingfileout = new TFile(timingfilename.c_str(),"RECREATE");
		timingfileout->cd();
		theftydb                  = new TTree("heftydb","");
		theftydb->Branch("SequenceID", &fileout_SequenceID);
		theftydb->Branch("Time", &timefileout_Time, TString::Format("Time[%d]/l",minibuffers_per_fullbuffer));
		theftydb->Branch("Label", &timefileout_Label, TString::Format("Label[%d]/I",minibuffers_per_fullbuffer));
		theftydb->Branch("TSinceBeam", &timefileout_TSinceBeam,
			TString::Format("TSinceBeam[%d]/L",minibuffers_per_fullbuffer));
		theftydb->Branch("More", &timefileout_More, TString::Format("More[%d]/I",minibuffers_per_fullbuffer));
	}
	
	gROOT->cd();
	
	FillInitialFileInfo();
	
}

void PulseSimulation::FillInitialFileInfo(){
	// We can generate some of the output values immediately, as they're just dummy values
	//===============================
	Log("PulseSimulation Tool: Filling run start date and setting constants",v_debug,verbosity);
	FillEmulatedRunInformation();
	
	// StartTimeSec represents the unix seconds of the start of the run
	// Convert "run date", read from the options file, to unixns
	// first convert config file string to parts
	int hh, mm, ss, MM, DD, YYYY;
	sscanf(runStartDateTime.c_str(), "%d/%d/%d_%d:%d:%d", &DD, &MM, &YYYY, &hh, &mm, &ss);
	// combine parts into a time structure
	struct std::tm runstart;
	runstart.tm_year = YYYY;
	runstart.tm_mon = MM;
	runstart.tm_mday = DD;
	runstart.tm_hour = hh;
	runstart.tm_min = mm;
	runstart.tm_sec = ss;
	// finally convert time structure to unix seconds
	time_t runstarttime = mktime(&runstart);
	// use runstarttime = time(NULL);  to get current time.
	RunStartTimeSec = static_cast<Int_t>(runstarttime);
	
	// give proper sizes to the vectors of unused timing variables
	StartCountVals.assign(num_adc_cards,0);
	StartTimeNSecVals.assign(num_adc_cards,0);
	
	// construct the empty data structures (vectors, arrays) needed to hold an ADC readout
	// ===================================================================================
	Log("PulseSimulation Tool: Constructing CardData entries",v_debug,verbosity);
	// fill all the non-minibuffer info and reset the minibuffers
	if((int)emulated_pmtdata_readout.size()<num_adc_cards){
//		cout<<"constructing vector of "<<num_adc_cards<<" CardData objects"
//			<<" and temp databuffer vector with "<<full_buffer_size<<" element data arrays"<<endl;
		emulated_pmtdata_readout= std::vector<MCCardData>(num_adc_cards);
		temporary_databuffers = std::vector< std::vector<uint16_t> >(num_adc_cards);
	}
	// In order for downstream tools to be able to retrieve the waveform vector pointers just
	// once in Initialize, we need to set up those vectors now. This is done in Card::Reset().
	for(auto& acard : emulated_pmtdata_readout){
		acard.Channels = channels_per_adc_card;                   // static
		acard.TriggerNumber = minibuffers_per_fullbuffer;         // static
		acard.FullBufferSize = full_buffer_size;                  // trigger window sizes and sampling rate
		acard.Eventsize = emulated_event_size;                    // PulseSimulation::emulated_event_size etc.
		acard.BufferSize = buffer_size;                           // minibuf_samples_per_ch * minibufs_per_fullbuf
		acard.Reset();                                            // allocate waveform storage
	}
}

void PulseSimulation::FillEmulatedRunInformation(){
	Log("PulseSimulation Tool: Filling the template Run Information tree",v_debug,verbosity);
	// RunInformation tree in RAW file is essentially a map of strings.
	// It stores (at present) 11 entries, each with a 'key' in the InfoTitle branch, 
	// and a corresponding xml object string in the InfoMessage branch - 
	// the xml string may represent multiple variables.
	// Since this is largely irrelevant, we'll just use template values and modify if necessary.
	
	std::vector<std::string> InfoTitles{"ToolChainVariables","InputVariables","PostgresVariables",
	"SlackBotVariables","HVComs","TriggerVariables","NetworkReceiveDataVariables","LoggerVariables",
	"RootDataRecorderVariables","MonitoringVariables","MRDVariables"};
	std::vector<std::string>* TemplateInfoMessages = GetTemplateRunInfo();
	
	if(GenerateFakeRootFiles){
		for(int entryi=0; entryi<(int)InfoTitles.size(); entryi++){
			fileout_InfoTitle=InfoTitles.at(entryi);
			fileout_InfoMessage=TemplateInfoMessages->at(entryi);
			tRunInformation->Fill();
		}
	} else {
		Log("PulseSimulation Tool: Skpping filling of Run Information tree",v_debug,verbosity);
	}
	if(PutOutputsIntoStore){
		// nothing worth doing here
	}
}
