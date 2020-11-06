/* vim:set noexpandtab tabstop=4 wrap */
#include "FindMrdTracks.h"
#include <numeric>      // std::iota
// for drawing
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
// for sleeping
#include <thread>          // std::this_thread::sleep_for
#include <chrono>          // std::chrono::seconds

FindMrdTracks::FindMrdTracks():Tool(){}

bool FindMrdTracks::Initialise(std::string configfile, DataModel &data){
	
	std::cout << "Initializing tool FindMrdTracks"<<std::endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	outputdir = "./MRDPlots/";
	triggertype = "Cosmic";
	outputfile = "mrdtrackfile";
	
	// get configuration variables for this tool
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("IsData",isData);
	m_variables.Get("OutputDirectory",outputdir);
	m_variables.Get("OutputFile",outputfile);
	m_variables.Get("DrawTruthTracks",DrawTruthTracks);
	m_variables.Get("WriteTracksToFile",writefile);
	m_variables.Get("SelectTriggerType",triggertype_selection);
	m_variables.Get("TriggerType",triggertype);
	
	if (triggertype == "NoLoopback") triggertype = "No Loopback";
	std::cout <<"User Trigger type: "<<triggertype<<std::endl;
	if (isData) {
		Log("FindMrdTracks tool: Selected DrawTruthTracks in configfile for a data file! Setting DrawTruthTracks to false",v_error,verbosity);	
		DrawTruthTracks = false;		//if looking at data, no possiblity to draw truth tracks
	}
	
	// create a store for holding MRD tracks to pass to downstream Tools
	// will be a single entry BoostStore containing a vector of single entry BoostStores
	m_data->Stores["MRDTracks"] = new BoostStore(true,0);
	// the vector of MRD tracks found in the current event
	theMrdTracks = new std::vector<BoostStore>(5);
	// in case we don't find any tracks in the first event, we should populate the entry in the MRDTracks
	// boost store so that other downstream tools that expect it will indeed find it.
	m_data->Stores["MRDTracks"]->Set("MRDTracks",theMrdTracks,false);
	m_data->Stores["MRDTracks"]->Set("NumMrdSubEvents",0);
	m_data->Stores["MRDTracks"]->Set("NumMrdTracks",0);
	
	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geo);
	//numvetopmts = geo->GetNumVetoPMTs();
	
	// create clonesarray for storing the MRD Track details as they're found
	if(SubEventArray==nullptr) SubEventArray = new TClonesArray("cMRDSubEvent");  // string is class name
	// put the pointer in the CStore, so it can be retrieved by MrdTrackPlotter tool Init
	intptr_t subevptr = reinterpret_cast<intptr_t>(SubEventArray);
	m_data->CStore.Set("MrdSubEventTClonesArray",subevptr);
	m_data->CStore.Get("mrd_tubeid_to_channelkey",mrd_tubeid_to_channelkey);
	// pass this to TrackCombiner
	m_data->CStore.Set("DrawMrdTruthTracks",DrawTruthTracks);
	
	return true;
}


bool FindMrdTracks::Execute(){
	
	if(verbosity) cout<<"Tool FindMrdTracks finding tracks in next event."<<endl;
	
	// ensure the previous tracks are cleared so we don't have any carry over
	SubEventArray->Clear("C");
	//m_data->Stores["MRDTracks"]->Delete(); // delete the previous vector of tracks - just clear and refill
	m_data->Stores["MRDTracks"]->Set("NumMrdSubEvents",0);
	m_data->Stores["MRDTracks"]->Set("NumMrdTracks",0);
	
	int lastrunnum=runnum;
	int lastsubrunnum=subrunnum;
	
	// Get the ANNIE event and extract information
	Log("FindMrdTracks tool getting event info from store",v_debug,verbosity);
	m_data->Stores["ANNIEEvent"]->Get("MCFile",MCFile);
	m_data->Stores["ANNIEEvent"]->Get("RunNumber",RunNumber);
	m_data->Stores["ANNIEEvent"]->Get("SubrunNumber",SubrunNumber);
	m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
	m_data->Stores["ANNIEEvent"]->Get("TriggerNumber",MCTriggernum);
	m_data->Stores["ANNIEEvent"]->Get("MCEventNum",MCEventNum);
	if (!isData) m_data->Stores["ANNIEEvent"]->Get("MCParticles",MCParticles);
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MRDTriggerType",MRDTriggertype);
	if (not get_ok){
		Log("FindMrdTracks: Did not find MRDTriggerType in ANNIEEvent. Please check the settings in MRDDataDecoder+BuildANNIEEvent/LoadWCSim?",v_error,verbosity);
	}
	Log("FindMrdTracks tool: MRDTriggertype is "+MRDTriggertype+" (from ANNIEEvent store)",v_debug,verbosity);
	
	// Get time cluster objects from CStore
	// MRDTime Clusters object should be created by previous tool in ToolChain (TimeClustering)
	get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
	if (not get_ok) {
		Log("FindMrdTracks: Did not find MrdTimeClusters in CStore. Did you run the tool TimeClustering beforehand?",v_error,verbosity);
		return true;  //it can happen that an event does not have a TDCData object if MRD is not hit
	}
	int NumMrdTimeClusters;
	get_ok = m_data->CStore.Get("NumMrdTimeClusters",NumMrdTimeClusters);
	if (not get_ok){
		Log("FindMrdTracks: Did not find NumMrdTimeClusters in CStore. Did you run the tool TimeClustering beforehand?",v_error,verbosity);
		return false;
	}
	Log("FindMrdTracks tool: Got MrdTimeClusters of size "+std::to_string(MrdTimeClusters.size())+", with "+std::to_string(NumMrdTimeClusters)+" MrdTimeClusters.",v_message,verbosity);
	
	// FIXME align types until we update MRDTrackClass/MRDSubEventClass
	currentfilestring=MCFile;
	runnum=(int)RunNumber;
	subrunnum=(int)SubrunNumber;
	eventnum=(int)EventNumber;
	triggernum=(int)MCTriggernum;
	
	// TODO: add MCEventNum to MRDTrackClass/MRDSubEventClass
	//
	// Setup the vector to store the information which subevent had a track (for downstream plotting tools)
	std::vector<int> track_subevs;
	
	// make the tree if we want to save the Tracks to file
	if( writefile && (mrdtrackfile==nullptr || (lastrunnum!=runnum) || (lastsubrunnum!=subrunnum)) ){
		Log("FindMrdTracks tool: Creating new file to write MRD tracks to",v_message,verbosity);
		StartNewFile();  // open a new file for each run / subrun
		
	}
	
	
	//verbose check whether the obtained MRD data seems sensible, MRDTimeClusters contains the indices of events belonging to that subevent
	if(verbosity >= v_debug){
		for (unsigned int i_cluster = 0; i_cluster < MrdTimeClusters.size(); i_cluster++){
			std::cout << "FindMrdTracks tool: Cluster " << i_cluster+1 << ", MrdTimeClusters.at(i_cluster).size() = " << MrdTimeClusters.at(i_cluster).size() << std::endl;
			for (unsigned int i_entry = 0; i_entry < MrdTimeClusters.at(i_cluster).size(); i_entry++){
				std::cout <<MrdTimeClusters.at(i_cluster).at(i_entry)<<", ";
			}
			std::cout << std::endl;
		}
	}
	
	// get the time and pmts digits from the CStore and put them into separate vectors used by the track finder, clear those containers
	mrddigittimesthisevent.clear();
	mrddigitpmtsthisevent.clear();
	mrddigitchargesthisevent.clear();
	
	///////////////////////////
	// now do the track finding
	
	Log("FindMrdTracks tool: Searching for MRD tracks in event "+std::to_string(eventnum),v_message,verbosity);
	
/*
if your class contains pointers, use TrackArray.Clear("C"). You MUST then provide a Clear() method in your class that properly performs clearing and memory freeing. (or "implements the reset procedure for pointer objects")
 see https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
*/
	int mrdeventcounter=0;
	
	//select only specific trigger types if prompted to do so
	//this is useful if one wants to look at e.g. only beam events for the analysis or only at cosmic events for the MRD paddle calibration
	
	if ((triggertype_selection && MRDTriggertype == triggertype) || !triggertype_selection ){
		
		// no clusters found
		if (MrdTimeClusters.size() == 0){
			nummrdsubeventsthisevent=0;
			nummrdtracksthisevent=0;
			if(writefile){
				nummrdsubeventsthiseventb->Fill();
				nummrdtracksthiseventb->Fill();
				subeventsinthiseventb->Fill();
				mrdtrackfile->cd();
				mrdtree->SetEntries(nummrdtracksthiseventb->GetEntries());
				mrdtree->Write("",TObject::kOverwrite);
				gROOT->cd();
			}
			Log("FindMrdTracks tool: No MRD digits in this event; returning",v_message,verbosity);
			
			return true;
			// XXX Note for other tools: we've written files and updated BoostStore SubEvent/Track counts to 0.
			// XXX Other tools should check these before processing the MrdTracks BoostStore or MrdSubEventTClonesArray
			// XXX since those are not necessarily cleared!
			// skip remainder
			// ======================================================================================
		}
		
		m_data->CStore.Get("MrdDigitTimes",mrddigittimesthisevent);
		m_data->CStore.Get("MrdDigitPmts",mrddigitpmtsthisevent);
		m_data->CStore.Get("MrdDigitCharges",mrddigitchargesthisevent);
		double eventendtime = *std::max_element(mrddigittimesthisevent.begin(),mrddigittimesthisevent.end());
		
		// LOOP THROUGH MRD SUBEVENTS
		// ===========================
		
		std::vector<int> digitidsinasubevent;
		std::vector<int> tubeidsinasubevent;
		std::vector<double> digitqsinasubevent;
		std::vector<double> digittimesinasubevent;
		std::vector<int> digitnumtruephots;
		std::vector<int> particleidsinasubevent;
		std::vector<double> photontimesinasubevent;
		
		Log("FindMrdTracks tool: "+std::to_string(MrdTimeClusters.size())+" subevents this event!",v_message,verbosity);
		
		// Now we need to sort the digits into the subevents they belong to:
		// Loop over subevents
		
		int mrdtrackcounter=0;   // not all the subevents will have a track
		
		for(unsigned int thiscluster=0; thiscluster<MrdTimeClusters.size(); thiscluster++){
			
			std::vector<int> single_mrdcluster = MrdTimeClusters.at(thiscluster);
			int numdigits = single_mrdcluster.size();
			
			for(int thisdigit=0;thisdigit<numdigits;thisdigit++){
				int digit_value = single_mrdcluster.at(thisdigit); // digit value is index of digit in TDC hits
				// TimeClustering tool builds clusters from both MRD and Veto hits,
				// but we cannot have veto PMTs in the MrdTrackLib reco algorithms
				// (they would be out of bounds in the expected maps...)
				// so convert back to channelkey, get the Detector, and check whether it's MRD or Veto
				int wcsimid = mrddigitpmtsthisevent.at(digit_value);
				if (!isData) wcsimid++;		//mrd_tubeid_to_channelkey map is 1-based in MC
				if(mrd_tubeid_to_channelkey.count(wcsimid)==0){
					Log("FindMrdTracks tool: Error! WCSimID "+to_string(wcsimid)
						+" was not in the mrd_tubeid_to_channelkey map!",v_error,verbosity);
					continue;
				} else {
					unsigned long chankey = mrd_tubeid_to_channelkey.at(wcsimid);
					Detector* thedetector = geo->ChannelToDetector(chankey);
					if(thedetector==nullptr){
						Log("FindMrdTracks Tool: Null detector in TDCData!",v_error,verbosity);
						continue;
					}
					if(thedetector->GetDetectorElement()!="MRD") continue; // this is a veto hit, not an MRD hit
				}
				digitidsinasubevent.push_back(digit_value);
				tubeidsinasubevent.push_back(mrddigitpmtsthisevent.at(digit_value));
				digittimesinasubevent.push_back(mrddigittimesthisevent.at(digit_value));
				digitqsinasubevent.push_back(mrddigitchargesthisevent.at(digit_value));
				
			}
			
			digitnumtruephots.assign(digittimesinasubevent.size(),0);      // FIXME replacement of above
			photontimesinasubevent.assign(digittimesinasubevent.size(),0); // FIXME replacement of above
			particleidsinasubevent.assign(digittimesinasubevent.size(),0); // FIXME replacement of above
			
			// Just initialize MCTruth variables, as before
			std::vector<std::pair<TVector3,TVector3>> truetrackvertices;
			std::vector<Int_t> truetrackpdgs;
			std::vector<std::pair<Position,Position>> truetrackvertices_position;			

			// In case we want to overlay the truth track in MC, also load the MCParticle information
			
			if (DrawTruthTracks){
				int numtracks = MCParticles->size();
				// a vector to record the subevent number for each track, to know if we've allocated it yet.
				std::vector<int> subtrackthisevent(numtracks,-1);
				std::vector<float> cluster_starttimes;
				m_data->CStore.Get("ClusterStartTimes",cluster_starttimes);
				float endtime = (thiscluster<(MrdTimeClusters.size()-1)) ?  cluster_starttimes.at(thiscluster+1) : (eventendtime+1.);
				for(int truetracki=0; truetracki<numtracks; truetracki++){
					MCParticle nextrack = MCParticles->at(truetracki);
					if((subtrackthisevent.at(truetracki)<0)&&(nextrack.GetStartTime()<endtime)
						&&(nextrack.GetPdgCode()!=11)){
						Position startvp = nextrack.GetStartVertex();
						TVector3 startv(startvp.X()*100.,startvp.Y()*100.,startvp.Z()*100.);
						Position stopvp = nextrack.GetStopVertex();
						TVector3 stopv(stopvp.X()*100.,stopvp.Y()*100.,stopvp.Z()*100.);
						truetrackvertices.push_back({startv,stopv});
						truetrackvertices_position.push_back({startvp,stopvp});
						if (verbosity >= v_message) std::cout <<"FindMrdTracks tool: True start vtx: ("<<startvp.X()<<","<<startvp.Y()<<","<<startvp.Z()<<"), true end vtx: ("<<stopvp.X()<<","<<stopvp.Y()<<","<<stopvp.Z()<<")"<<std::endl;
						truetrackpdgs.push_back(nextrack.GetPdgCode());
						subtrackthisevent.at(truetracki)=thiscluster;
					}
				}

				m_data->CStore.Set("TrueTrackVertices",truetrackvertices_position);
				m_data->CStore.Set("TrueTrackPDGs",truetrackpdgs);
			}
			
			Log("FindMrdTracks tool: Constructing subevent "+std::to_string(mrdeventcounter)+" with "+std::to_string(digitidsinasubevent.size())+" digits",v_message,verbosity);
			cMRDSubEvent* currentsubevent = new((*SubEventArray)[mrdeventcounter]) cMRDSubEvent(mrdeventcounter, currentfilestring, runnum, eventnum, triggernum, digitidsinasubevent, tubeidsinasubevent, digitqsinasubevent, digittimesinasubevent, digitnumtruephots, photontimesinasubevent, particleidsinasubevent, truetrackvertices, truetrackpdgs);
			if (currentsubevent->GetTracks()->size() > 0) track_subevs.push_back(mrdeventcounter); //annotate the current subevent number to have a track
			mrdeventcounter++;
			mrdtrackcounter+=currentsubevent->GetTracks()->size();
			Log("FindMrdTracksData: Subevent "+std::to_string(thiscluster)+" found "+std::to_string(currentsubevent->GetTracks()->size())+" tracks",v_message,verbosity);
			
			digitidsinasubevent.clear();
			tubeidsinasubevent.clear();
			digitqsinasubevent.clear();
			digittimesinasubevent.clear();
			particleidsinasubevent.clear();
			photontimesinasubevent.clear();
			digitnumtruephots.clear();
		
		}
		
		nummrdsubeventsthisevent=mrdeventcounter;
		nummrdtracksthisevent=mrdtrackcounter;
		
		if(writefile){
			nummrdsubeventsthiseventb->Fill();
			nummrdtracksthiseventb->Fill();
			subeventsinthiseventb->Fill();
		}
		
	} else {
		//did not pass the triggertype selection cut
		
		nummrdsubeventsthisevent=0;
		nummrdtracksthisevent=0;
		if(writefile){
			nummrdsubeventsthiseventb->Fill();
			nummrdtracksthiseventb->Fill();
			subeventsinthiseventb->Fill();
			mrdtrackfile->cd();
			mrdtree->SetEntries(nummrdtracksthiseventb->GetEntries());
			mrdtree->Write("",TObject::kOverwrite);
			gROOT->cd();
		}
		
		Log("FindMrdTracks tool: Trigger type selection cuts were not passed; returning",v_message,verbosity);
		return true;
		
		// XXX Note for other tools: we've written files and updated BoostStore SubEvent/Track counts to 0.
		// XXX Other tools should check these before processing the MrdTracks BoostStore or MrdSubEventTClonesArray
		// XXX since those are not necessarily cleared!
		// skip remainder
		// =====================================================================================
	
	}
	
	//if(eventnum==735){ assert(false); }
	//if(nummrdtracksthisevent) std::this_thread::sleep_for (std::chrono::seconds(5));
	
	// WRITE+CLOSE OUTPUT FILES
	// ========================
	
	Log("FindMrdTracks tool: Found "+std::to_string(nummrdtracksthisevent)+" MRD tracks in this event, in "+std::to_string(nummrdsubeventsthisevent)+" sub-events. End of finding MRD Tracks in this event",v_message,verbosity);
	
	if(writefile){
		Log("FindMrdTracks tool: Writing update to mrdtree in mrdtrackfile",v_message,verbosity);
		mrdtrackfile->cd();
		mrdtree->SetEntries(nummrdtracksthiseventb->GetEntries());
		mrdtree->Write("",TObject::kOverwrite);
	}
	if(cMRDSubEvent::imgcanvas) cMRDSubEvent::imgcanvas->Update();
	gROOT->cd();
	
	// save MRD tracks to the store for downstream tools
	Log("FindMrdTracks tool: Writing MRD tracks to MRDTracks BoostStore",v_debug,verbosity);
	m_data->Stores["MRDTracks"]->Set("NumMrdSubEvents",nummrdsubeventsthisevent);
	m_data->Stores["MRDTracks"]->Set("NumMrdTracks",nummrdtracksthisevent);
	if(nummrdtracksthisevent>(int)theMrdTracks->size()){
		Log("FindMrdTracks tool: Growing theMrdTracks vector to size "+std::to_string(nummrdtracksthisevent),v_debug,verbosity);
		theMrdTracks->resize(nummrdtracksthisevent);
	}
	
	int itrack_global=0;
	for(int subevi=0; subevi<nummrdsubeventsthisevent; subevi++){
		if (verbosity > v_debug) std::cout <<"FindMrdTracks tool: Looping through subevent "<<subevi<<std::endl;
		cMRDSubEvent* asubev = (cMRDSubEvent*)SubEventArray->At(subevi);
		// let's not save the SubEvent information. It's not much use.
		// N.B. BOOSTSTORE "MRDSubEvents" IS NOT CREATED
//		m_data->Stores["MRDSubEvents"]->Set("SubEventID",asubev->GetSubEventID());
//		m_data->Stores["MRDSubEvents"]->Set("NumDigits",asubev->GetNumDigits());
//		m_data->Stores["MRDSubEvents"]->Set("NumLayersHit",asubev->GetNumLayersHit());
//		m_data->Stores["MRDSubEvents"]->Set("NumPmtsHit",asubev->GetNumPMTsHit());
//		m_data->Stores["MRDSubEvents"]->Set("LayersHit",asubev->GetLayersHit());
//		m_data->Stores["MRDSubEvents"]->Set("LayerEdeps",asubev->GetEdeps());
//		m_data->Stores["MRDSubEvents"]->Set("PMTsHit",asubev->GetPMTsHit());
//		//m_data->Stores["MRDSubEvents"]->Set("TrueTracks",asubev->GetTrueTracks()); // WCSimRootTracks
//		m_data->Stores["MRDSubEvents"]->Set("DigitIds",asubev->GetDigitIds());       // WCSimRootCkvDigi indices
//		m_data->Stores["MRDSubEvents"]->Set("DigitQs",asubev->GetDigitQs());
//		m_data->Stores["MRDSubEvents"]->Set("DigitTs",asubev->GetDigitTs());
//		m_data->Stores["MRDSubEvents"]->Set("NumPhotonsInDigits",asubev->GetDigiNumPhots());
//		// since each digit may have mutiple photons the following are concatenated.
//		// scan through and keep a running total of NumPhots to identify the subarray start/ends
//		m_data->Stores["MRDSubEvents"]->Set("PhotTsInDigits",asubev->GetDigiPhotTs());
//		m_data->Stores["MRDSubEvents"]->Set("PhotParentsInDigits",asubev->GetDigiPhotParents());
		
		// can't nest multi-entry BoostStores
		//m_data->Stores["MRDSubEvents"]->Set("RecoTracks",asubev->GetTracks());       // cMRDTracks
		std::vector<cMRDTrack>* thetracks = asubev->GetTracks();
		for(unsigned int tracki=0; tracki<thetracks->size(); tracki++){
			cMRDTrack* atrack = &thetracks->at(tracki);
			
			// get the BoostStore to hold this track
			Log("FindMrdTracks: Getting Booststore at subevi = "+std::to_string(subevi)+" and tracki = "+std::to_string(tracki)+", global tracki = "+std::to_string(itrack_global),v_debug,verbosity);
			BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(itrack_global));
			itrack_global++;
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// If no tracks are found to match the reconstructed event vertex, the TrackCombiner Tool
			// may perform a cruder reconstruction and append an Mrd Track to this vector.
			// Since we don't remove and re-create the BoostStores each Execute, any members that are
			// not updated by both tools could result in spurious track properties being propagated!
			// !!!!!!!!!!!!! It is therefore crucial to keep these two Tools in sync! !!!!!!!!!!!!!!!
			// XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
			// IF ADDING MEMBERS TO THE MRDTRACK BOOSTSTORE HERE, ADD THEM TO THE TRACKCOMBINER TOO!
			// XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			// fill with this track's data
			thisTrackAsBoostStore->Set("MrdTrackID",atrack->GetTrackID());
			if(subevi!=atrack->GetMrdSubEventID()){
				Log("FindMrdTracks Tool Error: cMRDTrack::GetMrdSubEventID in wrong subevi!",v_error,verbosity);
			}
			thisTrackAsBoostStore->Set("MrdSubEventID",atrack->GetMrdSubEventID());
			thisTrackAsBoostStore->Set("InterceptsTank",atrack->GetInterceptsTank());
			thisTrackAsBoostStore->Set("StartTime",atrack->GetStartTime());
			// convert start posn from TVector3 to ANNIEEVENT Position class
			Position startpos(atrack->GetStartVertex().X() / 100.,
							  atrack->GetStartVertex().Y() / 100.,
							  atrack->GetStartVertex().Z() / 100.);
			// same for endpos
			Position endpos(  atrack->GetStopVertex().X() / 100.,
							  atrack->GetStopVertex().Y() / 100.,
							  atrack->GetStopVertex().Z() / 100.);
			thisTrackAsBoostStore->Set("StartVertex",startpos);
			thisTrackAsBoostStore->Set("StopVertex",endpos);
			thisTrackAsBoostStore->Set("TrackAngle",atrack->GetTrackAngle());
			thisTrackAsBoostStore->Set("TrackAngleError",atrack->GetTrackAngleError()); // TODO!!!
			thisTrackAsBoostStore->Set("LayersHit",atrack->GetLayersHit());
			thisTrackAsBoostStore->Set("TrackLength",atrack->GetTrackLength() / 100.);
			thisTrackAsBoostStore->Set("IsMrdPenetrating",atrack->GetIsPenetrating());
			thisTrackAsBoostStore->Set("EnergyLoss",atrack->GetEnergyLoss());
			thisTrackAsBoostStore->Set("EnergyLossError",atrack->GetEnergyLossError());
			thisTrackAsBoostStore->Set("IsMrdStopped",atrack->GetIsStopped());
			thisTrackAsBoostStore->Set("IsMrdSideExit",atrack->GetIsSideExit());
			thisTrackAsBoostStore->Set("PenetrationDepth",atrack->GetPenetrationDepth() / 100.);
			thisTrackAsBoostStore->Set("HtrackFitChi2",atrack->GetHtrackFitChi2());
			thisTrackAsBoostStore->Set("HtrackFitCov",atrack->GetHtrackFitCov());
			thisTrackAsBoostStore->Set("VtrackFitChi2",atrack->GetVtrackFitChi2());
			thisTrackAsBoostStore->Set("VtrackFitCov",atrack->GetVtrackFitCov());
			thisTrackAsBoostStore->Set("PMTsHit",atrack->GetPMTsHit());
			
			thisTrackAsBoostStore->Set("HtrackOrigin",atrack->GetHtrackOrigin());
			thisTrackAsBoostStore->Set("HtrackOriginError",atrack->GetHtrackOriginError());
			thisTrackAsBoostStore->Set("HtrackGradient",atrack->GetHtrackGradient());
			thisTrackAsBoostStore->Set("HtrackGradientError",atrack->GetHtrackGradientError());
			thisTrackAsBoostStore->Set("VtrackOrigin",atrack->GetVtrackOrigin());
			thisTrackAsBoostStore->Set("VtrackOriginError",atrack->GetVtrackOriginError());
			thisTrackAsBoostStore->Set("VtrackGradient",atrack->GetVtrackGradient());
			thisTrackAsBoostStore->Set("VtrackGradientError",atrack->GetVtrackGradientError());
			
			// convert back projections to Position class
			Position TankExitPoint( atrack->GetTankExitPoint().X() / 100.,
									atrack->GetTankExitPoint().Y() / 100.,
									atrack->GetTankExitPoint().Z() / 100.);
			Position MrdEntryPoint( atrack->GetMrdEntryPoint().X() / 100.,
									atrack->GetMrdEntryPoint().Y() / 100.,
									atrack->GetMrdEntryPoint().Z() / 100.);
			thisTrackAsBoostStore->Set("TankExitPoint",TankExitPoint);
			thisTrackAsBoostStore->Set("MrdEntryPoint",MrdEntryPoint);
			thisTrackAsBoostStore->Set("TrackIndex",tracki);
			
			// this stuff either isn't important or isn't yet implemented, don't store:
//			thisTrackAsBoostStore->Set("NumPMTsHit",atrack->GetNumPMTsHit());
//			thisTrackAsBoostStore->Set("KEStart",atrack->GetKEStart());
//			thisTrackAsBoostStore->Set("KEEnd",atrack->GetKEEnd());
//			thisTrackAsBoostStore->Set("ParticlePID",atrack->GetParticlePID());
//			thisTrackAsBoostStore->Set("NumLayersHit",atrack->GetNumLayersHit());
//			thisTrackAsBoostStore->Set("LayerEdeps",atrack->GetEdeps());
//			thisTrackAsBoostStore->Set("NumDigits",atrack->GetNumDigits());
//			thisTrackAsBoostStore->Set("MrdEntryBoundsX",atrack->GetMrdEntryBoundsX());
//			thisTrackAsBoostStore->Set("MrdEntryBoundsY",atrack->GetMrdEntryBoundsY());
//			//thisTrackAsBoostStore->Set("TrueTrackID",atrack->GetTrueTrackID());
//			//thisTrackAsBoostStore->Set("TrueTrack",atrack->GetTrueTrack());
//			thisTrackAsBoostStore->Set("DigitIds",atrack->GetDigitIds());
//			thisTrackAsBoostStore->Set("DigitQs",atrack->GetDigitQs());
//			thisTrackAsBoostStore->Set("DigitTs",atrack->GetDigitTs());
//			thisTrackAsBoostStore->Set("DigiNumPhots",atrack->GetDigiNumPhots());
//			thisTrackAsBoostStore->Set("DigiPhotTs",atrack->GetDigiPhotTs());
//			thisTrackAsBoostStore->Set("DigiPhotParents",atrack->GetDigiPhotParents());
			
			// differentiate tracks found with this algorithm from those found by the TrackCombiner tool
			thisTrackAsBoostStore->Set("LongTrack",1);
		}
	}
	std::cout <<"Setting MRDTracks"<<std::endl;
	m_data->Stores["MRDTracks"]->Set("MRDTracks",theMrdTracks,false); // the address may be changed by resize()
	//m_data->Stores["MRDTracks"]->Save();   // write to file?
	
	// the TClonesArray doesn't get deleted (for efficiency reasons) so shove a pointer
	// into the CStore in case something wants to use it (e.g. MrdPaddlePlot tool)
	//m_data->CStore.Set("MrdSubEventTClonesArray",SubEventArray,false);
	// XXX shouldn't be necessary as SubEventArray is re-used, so pointer has not changed
	std::cout <<"Setting subevptr in CStore"<<std::endl;
	intptr_t subevptr = reinterpret_cast<intptr_t>(SubEventArray);
	m_data->CStore.Set("MrdSubEventTClonesArray",subevptr);
	
	// Set the information which subevent had a track
	m_data->CStore.Set("TracksSubEvs",track_subevs);
	
	//The following is just for debugging purposes
	//std::cout <<"FindMrdTracks tool: List of objects (End of Execute): "<<std::endl;
	//gObjectTable->Print();
	
	return true;
}

bool FindMrdTracks::Finalise(){
	
	// close output file
	if(mrdtrackfile){
		mrdtrackfile->Close();
		delete mrdtrackfile;
		mrdtrackfile=nullptr;
	}
	
	// clean up any BoostStore items
	m_data->Stores["MRDTracks"]->Delete();
	
	// clean up the space from the no longer needed TClonesArray
	// SubEventArray->Clear("C");
	Log("FindMrdTracks Tool: Calling SubEventArray->Delete()",v_message,verbosity);
	if(SubEventArray){ SubEventArray->Delete(); delete SubEventArray; SubEventArray=0;}
	
	Log("FindMrdTracks tool: Exiting",v_message,verbosity);
	return true;
}

void FindMrdTracks::StartNewFile(){
	TString filenameout = TString::Format("%s/%s.%d.%d.root",outputdir.c_str(),outputfile.c_str(),runnum,subrunnum);
	if (verbosity) std::cout<<"FindMrdTracks tool: Creating mrd output file "<<filenameout.Data()<<std::endl;
	if(mrdtrackfile) mrdtrackfile->Close();
	mrdtrackfile = new TFile(filenameout.Data(),"RECREATE","MRD Tracks file");
	mrdtrackfile->cd();
	mrdtree = new TTree("mrdtree","Tree for reconstruction data");
	mrdeventnumb = mrdtree->Branch("EventID",&eventnum);
	mrdtriggernumb = mrdtree->Branch("TriggerID",&triggernum);
	nummrdsubeventsthiseventb = mrdtree->Branch("nummrdsubeventsthisevent",&nummrdsubeventsthisevent);
	subeventsinthiseventb = mrdtree->Branch("subeventsinthisevent",&SubEventArray, nummrdsubeventsthisevent);
	nummrdtracksthiseventb = mrdtree->Branch("nummrdtracksthisevent",&nummrdtracksthisevent);
	gROOT->cd();
}
