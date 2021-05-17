/* vim:set noexpandtab tabstop=4 wrap */
#include "PrintANNIEEvent.h"

PrintANNIEEvent::PrintANNIEEvent():Tool(){}

bool PrintANNIEEvent::Initialise(std::string configfile, DataModel &data){
	
	if(verbose) cout<<"Initializing tool PrintANNIEEvent"<<endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	is_mc = false;
	has_raw = false;
	first_event = true;

	m_variables.Get("verbose",verbose);
	m_variables.Get("HasRaw",has_raw);
	//verbose=10;
	
	n_prompt = 0;
	n_ext = 0;
	n_ext_cc = 0;
	n_ext_nc = 0;

	return true;
}

bool PrintANNIEEvent::Execute(){
	
	if(verbose) cout<<"executing PrintANNIEEvent"<<endl;
	int annieeventexists = m_data->Stores.count("ANNIEEvent");
	if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl; return false;};
	//m_data->Log->Log(logmessage.str(),1,verbose);
	//logmessage.str(""); ??

        //Print general content overview of ANNIEEvent Store 
        if (first_event){
		Log("PrintANNIEEvent tool: Contents of ANNIEEvent BoostStore: ",0,verbose);
		m_data->Stores["ANNIEEvent"]->Print(false);
		first_event = false;
	}
	
	// Get the contents to be printed
	if(verbose>1) cout<<"getting contents to be printed from ANNIEEvent at "<<m_data->Stores["ANNIEEvent"]<<endl;
	get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber",RunNumber);
	if(get_ok==0) RunNumber=-1;
	get_ok = m_data->Stores["ANNIEEvent"]->Get("SubrunNumber",SubrunNumber);
	if(get_ok==0) SubrunNumber=-1;
	get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
	if(get_ok==0) EventNumber=-1;
	cout<<"RunNumber : "<<RunNumber<<endl;
	cout<<"SubrunNumber : "<<SubrunNumber<<endl;
	cout<<"EventNumber : "<<EventNumber<<endl;
	
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCFlag",MCFlag);
	if (get_ok) {
		if (MCFlag) is_mc = true;
	}

	if (is_mc){
        	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles",MCParticles);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("RecoParticles",RecoParticles);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCHits",MCHits);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits",MCLAPPDHits);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData",MCTDCData);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCFlag",MCFlag);
		if(get_ok==0) MCFlag=0;
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCTriggernum",MCTriggernum);
		if(get_ok==0) MCTriggernum=-1;
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCFile",MCFile);
		if(get_ok==0) MCFile="";
		get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum",MCEventNum);
		if(get_ok==0) MCEventNum=-1;
		get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerData",TriggerData);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("BeamStatus",BeamStatusMC);
		get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);

		if(MCParticles){
			cout<<"Num MCParticles : "<<MCParticles->size()<<endl;
			if(verbose>1){
				cout<<"MCParticles: {"<<endl;
				for(auto&& aparticle : *MCParticles) aparticle.Print();
				cout<<"}"<<endl;
			}
		} else {
			cout<<"No MCParticles"<<endl;
		}
		if(RecoParticles){
			cout<<"Num RecoParticles : "<<RecoParticles->size()<<endl;
			if(verbose>1){
				cout<<"RecoParticles: {"<<endl;
				for(auto&& aparticle : *RecoParticles) aparticle.Print();
				cout<<"}"<<endl;
			}
		} else {
			cout<<"No RecoParticles"<<endl;
		}
		if(MCHits){
			cout<<"Num MCHits : "<<MCHits->size()<<endl;
			if(verbose>2){
				cout<<"MCHits : {"<<endl;
				for(auto&& achannel : *MCHits){
					unsigned long chankey = achannel.first;
					auto& hits = achannel.second;
					cout<<"ChannelKey : "<<chankey<<endl;
					cout<<"Hits : "<<endl;
					for(auto&& ahit : hits) ahit.Print();
					cout<<endl;
				}
				cout<<"}"<<endl;
			}
		} else {
		cout<<"No MCHits"<<endl;
		}
		if(MCLAPPDHits){
			cout<<"Num MCLAPPDHits : "<<MCLAPPDHits->size()<<endl;
			if(verbose>2){
				cout<<"MCLAPPDHits : {"<<endl;
				for(auto&& achannel : *MCLAPPDHits){
					unsigned long chankey = achannel.first;
					auto& hits = achannel.second;
					cout<<"ChannelKey : "<<chankey<<endl;
					cout<<"Hits : "<<endl;
					for(auto&& ahit : hits) ahit.Print();
					cout<<endl;
				}
				cout<<"}"<<endl;
			}
		} else {
			cout<<"No MCLAPPDHits"<<endl;
		}
		if(MCTDCData){
			cout<<"Num TDCData Hits : "<<MCTDCData->size()<<endl;
			if(verbose>2){
				cout<<"MCTDCData : {"<<endl;
				for(auto&& achannel : *MCTDCData){
					unsigned long chankey = achannel.first;
					auto& hits = achannel.second;
					cout<<"ChannelKey : "<<chankey<<endl;
					cout<<"Hits : "<<endl;
					for(auto&& ahit : hits) ahit.Print();
					cout<<endl;
				}
				cout<<"}"<<endl;
			}
		} else {
			cout<<"No TDCData"<<endl;
		}
		if(TriggerData){
			cout<<"Num TriggerData Entries : "<<TriggerData->size()<<endl;
			if(verbose>1){
				cout<<"TriggerData : {"<<endl;
				for(auto&& atrigger : *TriggerData) atrigger.Print();
				cout<<"}"<<endl;
			}
		} else {
			cout<<"No TriggerData"<<endl;
		}
		cout<<"MCFlag : "<<MCFlag<<endl;
		if(EventTime) { cout<<"EventTime : "; EventTime->Print(); }
		else { cout<<"No EventTime"<<endl; }
		cout<<"MCEventNum : "<<MCEventNum<<endl;
		cout<<"MCFile : "<<MCFile<<endl;
		if(BeamStatusMC) { cout<<"BeamStatusMC : "; BeamStatusMC->Print(); }
		else { cout<<"No BeamStatus (MC)"<<endl; }

		//End of MC loop

	} else {
		
		if (has_raw){
			bool get_rawadc, get_rawlappd, get_calibadc, get_caliblappd;
			get_rawadc = m_data->Stores["ANNIEEvent"]->Get("RawADCData",RawADCData);
			get_rawlappd = m_data->Stores["ANNIEEvent"]->Get("RawLAPPDData",RawLAPPDData);
			get_calibadc = m_data->Stores["ANNIEEvent"]->Get("CalibratedADCData",CalibratedADCData);
			get_caliblappd = m_data->Stores["ANNIEEvent"]->Get("CalibratedLAPPDData",CalibratedLAPPDData);
			if(get_rawadc){
				cout<<"Num RawADCData Waveforms : "<<RawADCData.size()<<endl;
				if(verbose>5){
					cout<<"RawADCData : {"<<endl;
					for(auto& achannel : RawADCData){
						unsigned long chankey = achannel.first;
						auto& waveforms = achannel.second;
						cout<<"ChannelKey : "<<chankey<<endl;
						cout<<"Has "<<waveforms.size()<<" waveforms"<<endl;
						if(verbose>2){
							cout<<"Waveforms : "<<endl;
							for(auto& awaveform : waveforms) awaveform.Print();
							cout<<endl;
						}
						cout<<"}"<<endl;
					}
				}
			} else cout<<"No RawADCData in ANNIEEvent"<<endl;
	
			if(get_rawlappd){
				cout<<"Num RawLAPPDData Waveforms : "<<RawLAPPDData.size()<<endl;
				if(verbose>5){
					cout<<"RawLAPPDData : {"<<endl;
					for(auto& achannel : RawLAPPDData){
						unsigned long chankey = achannel.first;
						auto& waveforms = achannel.second;
						cout<<"ChannelKey : "<<chankey<<endl;
						cout<<"Has "<<waveforms.size()<<" waveforms"<<endl;
						if(verbose>2){
							cout<<"Waveforms : "<<endl;
							for(auto& awaveform : waveforms) awaveform.Print();
							cout<<endl;
						}
						cout<<"}"<<endl;
					}
				}
			} else cout<<"No RawLAPPDData"<<endl;

			if(get_calibadc){
				cout<<"Num CalibratedADCData Waveforms : "<<CalibratedADCData.size()<<endl;
				if(verbose>5){
					cout<<"CalibratedADCData : {"<<endl;
					for(auto& achannel : CalibratedADCData){
						unsigned long chankey = achannel.first;
						auto& waveforms = achannel.second;
						cout<<"ChannelKey : "<<chankey<<endl;
						cout<<"Has "<<waveforms.size()<<" waveforms"<<endl;
						if(verbose>2){
							cout<<"Waveforms : "<<endl;
							for(auto& awaveform : waveforms) awaveform.Print();
							cout<<endl;
						}
						cout<<"}"<<endl;
					}
				}
			} else cout<<"No CalibratedADCData in ANNIEEvent"<<endl;

			if(get_caliblappd){
				cout<<"Num CalibratedLAPPDData Waveforms : "<<CalibratedLAPPDData.size()<<endl;
				if(verbose>5){
					cout<<"CalibratedLAPPDData : {"<<endl;
					for(auto& achannel : CalibratedLAPPDData){
						unsigned long chankey = achannel.first;
						auto& waveforms = achannel.second;
						cout<<"ChannelKey : "<<chankey<<endl;
						cout<<"Has "<<waveforms.size()<<" waveforms"<<endl;
						if(verbose>2){
							cout<<"Waveforms : "<<endl;
							for(auto& awaveform : waveforms) awaveform.Print();
							cout<<endl;
						}
						cout<<"}"<<endl;
					}
				}
			} else cout <<"Did not find CalibratedLAPPDData in ANNIEEvent"<<endl;

		} else {
			bool get_hits, get_auxhits, get_recoadc, get_recoauxadc, get_rawacqsize;
			get_hits = m_data->Stores["ANNIEEvent"]->Get("Hits",Hits);
			get_auxhits = m_data->Stores["ANNIEEvent"]->Get("AuxHits",AuxHits);
			get_recoadc = m_data->Stores["ANNIEEvent"]->Get("RecoADCData",RecoADCData);
			get_recoauxadc = m_data->Stores["ANNIEEvent"]->Get("RecoAuxADCData",RecoAuxADCData);
			get_rawacqsize = m_data->Stores["ANNIEEvent"]->Get("RawAcqSize",RawAcqSize);
			if (get_hits){

                        	cout<<"Num Hits : "<<Hits->size()<<endl;
                        	if(verbose>5){
                                	cout<<"Hits : {"<<endl;
                                	for(auto&& achannel : *Hits){ 
                                        	unsigned long chankey = achannel.first;
                                        	auto& hits = achannel.second;
                                        	cout<<"ChannelKey : "<<chankey<<endl;
                                        	cout<<"Hits : "<<endl;
                                        	for(auto&& ahit : hits) ahit.Print();
                                        	cout<<endl;
                                	}
                                	cout<<"}"<<endl;
                        	}
                	} else cout<<"No Hits in ANNIEEvent"<<endl;

			if (get_auxhits){
				cout <<"Num Aux Hits: "<<AuxHits->size()<<std::endl;
				if (verbose > 5){
					cout <<"AuxHits : {"<<endl;
                                	for(auto&& achannel : *AuxHits){ 
                                        	unsigned long chankey = achannel.first;
                                        	auto& hits = achannel.second;
                                        	cout<<"ChannelKey : "<<chankey<<endl;
                                        	cout<<"Hits : "<<endl;
                                        	for(auto&& ahit : hits) ahit.Print();
                                        	cout<<endl;
                                	}
                                	cout<<"}"<<endl;		
				}
			} else cout <<"No AuxHits in ANNIEEvent"<<endl;

			if (get_recoadc){
				cout <<"Num RecoADCData: "<<RecoADCData.size()<<std::endl;
                                if (verbose > 5){
                                        cout <<"RecoADCPulses : {"<<endl;
                                        for(auto& achannel : RecoADCData){
                                                unsigned long chankey = achannel.first;
                                                auto& pulses = achannel.second;
                                                cout<<"ChannelKey : "<<chankey<<endl;
                                                cout<<"Pulses size : "<<pulses.size()<<endl;
						for (int i_pulse=0; i_pulse < pulses.size(); i_pulse++){
							//pulses.at(i_pulse).Print();
						}
                                        }
                                        cout<<"}"<<endl;
                                }
			} else cout <<"No RecoADCHits in ANNIEEvent"<<std::endl;

			if (get_recoauxadc){
			        cout <<"Num RecoAuxADCData: "<<RecoAuxADCData.size()<<std::endl;
                                if (verbose > 5){
                                        cout <<"RecoAuxADCPulses : {"<<endl;
                                        for(auto& achannel : RecoAuxADCData){
                                                unsigned long chankey = achannel.first;
                                                auto& pulses = achannel.second;
                                                cout<<"ChannelKey : "<<chankey<<endl;
                                                cout<<"Pulses size : "<<pulses.size()<<endl;
                                        }       
                                        cout<<"}"<<endl;
                                }
			} else cout <<"No RecoAuxADCData in ANNIEEvent"<<std::endl;

			if (get_rawacqsize){
                                cout <<"Num RawAcqSize: "<<RawAcqSize.size()<<std::endl;
                                if (verbose > 3){
                                        cout <<"RawAcqSize entries : {"<<endl;
                                        for(auto& achannel : RawAcqSize){
                                                unsigned long chankey = achannel.first;
                                                auto& acq = achannel.second;
                                                cout<<"ChannelKey : "<<chankey<<endl;
                                                cout<<"Number of acquisition windows : "<<acq.size()<<endl;
						for (int i_win=0; i_win < (int) acq.size(); i_win++){
							cout <<"Acquisition Window "<<i_win+1<<", Window Size: "<<acq.at(i_win)<<" ADC Counts"<<endl;
						}
                                        }       
                                        cout<<"}"<<endl;
                                }
			} else cout <<"No RawAcqSize in ANNIEEvent"<<std::endl;
		}

		

		bool get_tdc, get_mrdtrig, get_mrdloopback, get_tmrd, get_ttank, get_tctc, get_trigword, get_datastreams;
		bool get_beamstatus, get_partnr, get_trigextended, get_eventtime, get_trigdata;
		get_tdc = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);
		get_mrdtrig = m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType",MRDTriggertype);
		get_mrdloopback = m_data->Stores["ANNIEEvent"]->Get("MRDLoopbackTDC",MRDLoopbackTDC);
		get_tmrd = m_data->Stores["ANNIEEvent"]->Get("EventTimeMRD",EventTimeMRD);
		get_ttank = m_data->Stores["ANNIEEvent"]->Get("EventTimeTank",EventTimeTank);
		get_tctc = m_data->Stores["ANNIEEvent"]->Get("CTCTimestamp",CTCTimestamp);
		get_trigword = m_data->Stores["ANNIEEvent"]->Get("TriggerWord",TriggerWord);
		get_trigdata = m_data->Stores["ANNIEEvent"]->Get("TriggerData",TriggerDataData);
		get_datastreams = m_data->Stores["ANNIEEvent"]->Get("DataStreams",DataStreams);
		get_beamstatus = m_data->Stores["ANNIEEvent"]->Get("BeamStatus",BeamStatusData);
		get_partnr = m_data->Stores["ANNIEEvent"]->Get("PartNumber",PartNumber);
		get_trigextended = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",TriggerExtended);

           

		//General event building statistics
		if (get_partnr) cout <<"Part Number of file: "<<PartNumber<<std::endl;
		else cout <<"Did not find File Part Number in ANNIEEvent"<<endl;
		if (get_tmrd) cout <<"Event Time MRD: "<<EventTimeMRD<<std::endl;
		else cout <<"Did not find MRD Event Time in ANNIEEvent"<<std::endl;
		if (get_ttank) cout <<"Event Time Tank: "<<EventTimeTank<<std::endl;
		else cout <<"Did not find Tank Event Time in ANNIEEvent"<<std::endl;
		if (get_tctc) cout <<"Event Time CTC: "<<CTCTimestamp<<std::endl;
		else cout <<"Did not find CTC Event Time in ANNIEEvent"<<std::endl;
		if (get_trigword) cout <<"Triggerword: "<<TriggerWord<<std::endl;
		else cout <<"Did not find Triggerword in ANNIEEvent"<<std::endl;
		if (get_trigextended) {
			if (TriggerExtended == 0) {cout <<"No extended readout"<<endl; n_prompt++;}
			else if (TriggerExtended == 1) {cout <<"Extended readout (CC)"<<endl; n_ext++; n_ext_cc++;}
			else if (TriggerExtended == 2) {cout <<"Extended readout (Non-CC)"<<endl; n_ext++; n_ext_nc++;}
			else cout <<"Unrecognized TriggerExtended value "<<TriggerExtended<<endl;
		} else {
			cout <<"Did not find TriggerExtended in ANNIEEvent"<<std::endl;
		}
		if (get_trigdata){
			cout <<"TriggerData contents: "<<std::endl;
			TriggerDataData.Print();
			cout << endl;
		} else cout <<"Did not find TriggerData object in ANNIEEvent"<<std::endl;

		//DataStreams Information
		if (get_datastreams) cout <<"DataStreams: Tank = "<<DataStreams["Tank"]<<", MRD = "<<DataStreams["MRD"]<<", CTC = "<<DataStreams["CTC"]<<", LAPPD = "<<DataStreams["LAPPD"]<<std::endl;
		else cout <<"Did not find DataStreams in ANNIEEvent"<<std::endl;
	
		//BeamStatus Information
		if (get_beamstatus) {
			cout <<"BeamStatus information: "<<BeamStatusData.Print();
			cout <<endl;
		} else cout <<"Did not find BeamStatus in ANNIEEvent"<<endl;

		//MRD information
		if (get_mrdtrig) cout <<"MRD Loopback Triggertype: "<<MRDTriggertype<<std::endl;
		else cout <<"Did not find MRD Loopback Triggertype in ANNIEEvent"<<endl;
		if (get_mrdloopback) {

		} else cout <<"Did not find MRDLoopbackTDC in ANNIEEvent"<<std::endl;
		if (get_tdc){
                        cout<<"Num TDCData Hits : "<<TDCData->size()<<endl;
                        if(verbose>5){
                                cout<<"TDCData : {"<<endl;
                                for(auto&& achannel : *TDCData){
                                        unsigned long chankey = achannel.first;
                                        auto& hits = achannel.second;
                                        cout<<"ChannelKey : "<<chankey<<endl;
                                        cout<<"Hits : "<<endl;
                                        for(auto&& ahit : hits) ahit.Print();
                                        cout<<endl;
                                }
                                cout<<"}"<<endl;
                        }
		} else cout <<"Did not find TDCData in ANNIEEvent"<<std::endl;


	}	//End of data loop
		

	
	return true;
}

bool PrintANNIEEvent::Finalise(){
	
       std::cout << "Number of prompt: "<<n_prompt<<std::endl;
       std::cout <<"Number of extended: "<<n_ext<<std::endl;
       std::cout <<"Number of extended (CC): "<<n_ext_cc<<std::endl;
       std::cout <<"Number of extended (NC): "<<n_ext_nc<<std::endl;

	return true;
}
