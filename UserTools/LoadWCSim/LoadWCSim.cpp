/* vim:set noexpandtab tabstop=4 wrap */

#include "LoadWCSim.h"

LoadWCSim::LoadWCSim():Tool(){}

bool LoadWCSim::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	
	if(verbosity) cout<<"Initializing Tool LoadWCSim"<<endl;
	
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	
	// Get the Tool configuration variables
	// ====================================
	get_ok = m_variables.Get("verbose", verbosity);
	if(not get_ok) verbosity=1;
	get_ok = m_variables.Get("MaxEntries",MaxEntries);
	if(not get_ok) MaxEntries=-1;
	get_ok = m_variables.Get("InputFile", MCFile);
	if(not get_ok){ Log("LoadWCSim Tool: ERROR: No InputFile defined!",v_error,verbosity); return false; }
	get_ok = m_variables.Get("HistoricTriggeroffset", HistoricTriggeroffset); // [ns]
	if(not get_ok){
		Log("LoadWCSim Tool: ERROR: HistoricTriggeroffset not set in config file!",v_error,verbosity);
		return false;
	}
	get_ok = m_variables.Get("WCSimVersion", WCSimVersion);
		if(not get_ok){
			Log("LoadWCSim Tool: ERROR: WCSimVersion not set in config file!",v_error,verbosity);
			return false;
	}
	get_ok = m_variables.Get("UseDigitSmearedTime", use_smeared_digit_time);
	if(not get_ok){
		Log("LoadWCSim Tool: Assuming to use smeared digit times",v_warning,verbosity);
		use_smeared_digit_time=1;
	}
	get_ok = m_variables.Get("LappdNumStrips", LappdNumStrips);
	if(not get_ok){
		Log("LoadWCSim Tool: Assuming to use 56 LAPPD striplines",v_warning,verbosity);
		LappdNumStrips = 56;
	}
	get_ok = m_variables.Get("LappdStripLength", LappdStripLength);           // [mm]
	if(not get_ok){
		Log("LoadWCSim Tool: Assuming to LAPPD stripline length of 200mm",v_warning,verbosity);
		LappdStripLength = 200;
	}
	get_ok = m_variables.Get("LappdStripSeparation", LappdStripSeparation);   // [mm]
	if(not get_ok){
		Log("LoadWCSim Tool: Assuming LAPPD stripline separation of 7.14mm",v_warning,verbosity);
		LappdStripSeparation = 7.14;
	}
	get_ok = m_variables.Get("RunStartDate", RunStartUser);
	if(not get_ok){
		Log("LoadWCSim Tool: Assuming RunStartDate of 0ns, i.e. unix epoch",v_warning,verbosity);
		RunStartUser = 0;
	}
	get_ok = m_variables.Get("TriggerType",Triggertype);
	if (not get_ok){
		Log("LoadWCSim Tool: No Triggertype specified. Assuming TriggerType = Beam",v_warning,verbosity);
		Triggertype = "Beam";	//other options: Cosmic / No Loopback
	}
	MCEventNum=0;
	get_ok = m_variables.Get("FileStartOffset",MCEventNum);
	
	// put version in the CStore for downstream tools
	m_data->CStore.Set("WCSimVersion", WCSimVersion);
	
	// Short Stores README
	//////////////////////
	// n.b. m_data->vars is a Store (of ben's Store type) that is not saved to disk?
	//      m_data->CStore is a single entry binary BoostStore that is not saved to disk.
	//      m_data->Stores["StoreName"] is a map of binary BoostStores that are saved to disk.
	// If using Stores->BoostStore->Set("MyVariable") it will always be saved to disk
	// Using Stores->BoostStore.Set("MyVariable",&myvar) if myvar is a pointer (to an object on
	// the heap) puts myvar in the Store and it's deletion will be handled by the Store.
	// (provided your class has a suitable destructor.)
	// Is 'BoostStore::Save' needed for single-entry stores?
	// ----------------
	// create a new BoostStore with key "ANNIEEvent" in the Stores std::map
	// BoostStore constructor args: typechecking (bool), m_format (0=binary, 1=ASCII, 2=multievent)
	// A BoostStore has a header where useful constants may be saved. The header is a BoostStore itself,
	// and can be accessed via: Store.Header->Get() and Store.Header->Set().
	// The method 'Store::Save()' writes everything 'Set' since the last 'Save' to the current entry.
	// Store::Clear clears the map of the current entry, to start building a new one.
	// Use 'Store::GetEntry(int entrynum)' to load an entry to then be able to 'Get' it's contents.
	// 'Store->Header->Get("TotalEntries",NumEvents)' will load the num entries into NumEvents
	// ------------------
	// When adding a BoostStore (or class object in general) to a BoostStore (such as ANNIEEvent)
	// you call BoostStore::Set("key",ObjectPointer) - BUT be aware that serialization happens when the
	// 'Set' method is called - although you pass it a pointer, any subsequent changes to the object
	// will NOT get saved! You must call 'Set' AFTER making ALL changes to your object!
	/////////////////////////////////////////////////////////////////
	
	// Make class private members; e.g. the WCSimT and WCSimRootGeom
	// =============================================================
//	file= new TFile(MCFile.c_str(),"READ");
//	wcsimtree= (TTree*) file->Get("wcsimT");
//	WCSimEntry= new wcsimT(wcsimtree);
	WCSimEntry= new wcsimT(MCFile.c_str(),verbosity);
	
	gROOT->cd();
	wcsimrootgeom = WCSimEntry->wcsimrootgeom;
	wcsimrootopts = WCSimEntry->wcsimrootopts;
	int pretriggerwindow=wcsimrootopts->GetNDigitsPreTriggerWindow();
	int posttriggerwindow=wcsimrootopts->GetNDigitsPostTriggerWindow();
	
	// put useful constants into the CStore
	// ====================================
	//m_data->CStore.Set("WCSimEntry",WCSimEntry,false); // pass on the WCSim entry - not possible
	// pass on root options. TODO should these be saved somewhere?
	m_data->CStore.Set("WCSimPreTriggerWindow",pretriggerwindow);
	m_data->CStore.Set("WCSimPostTriggerWindow",posttriggerwindow);
	// store the WCSimRootGeom for LAPPD reader calculation of global coords
	intptr_t geomptr = reinterpret_cast<intptr_t>(wcsimrootgeom);
	if(verbosity>1) cout<<"wcsimrootgeom at "<<wcsimrootgeom<<", geomptr="<<geomptr<<endl;
	m_data->CStore.Set("WCSimRootGeom",geomptr);
//	int wcsimgeomexists = m_data->Stores.count("WCSimRootGeomStore");
//	if(wcsimgeomexists==0){
//		m_data->Stores["WCSimRootGeomStore"] = new BoostStore(false,0);
//		//m_data->Stores.at("WCSimRootGeomStore")->Header->Set("WCSimRootGeom",wcsimrootgeom);
//		m_data->Stores.at("WCSimRootGeomStore")->Set("WCSimRootGeom",&wcsimrootgeom);
//	}
//	
//	// Make a WCSimStore to store additional WCSim info passed between tools
//	// =====================================================================
//	int wcsimstoreexists = m_data->Stores.count("WCSimStore");
//	cout<<"wcsimstoreexists="<<wcsimstoreexists<<endl;
//	if(wcsimstoreexists==0){
//		m_data->Stores["WCSimStore"] = new BoostStore(false,0);
//	}
//	m_data->Stores.at("WCSimStore")->Set("WCSimRootGeom",geomptr);
//	m_data->Stores.at("WCSimStore")->Set("WCSimPreTriggerWindow",pretriggerwindow);
//	m_data->Stores.at("WCSimStore")->Set("WCSimPostTriggerWindow",posttriggerwindow);
	
	// Make the ANNIEEvent Store if it doesn't exist
	// =============================================
	int annieeventexists = m_data->Stores.count("ANNIEEvent");
	if(annieeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);
	
	// Convert WCSimRootGeom into ToolChain Geometry class
	Geometry* anniegeom = ConstructToolChainGeometry();
	
	// Set run-level information in the ANNIEEvent
	// ===========================================
	/*
	At time of writing ANNIEEvent has the following that need to be Set before each Save:
		RunNumber
		SubrunNumber
		EventNumber
		MCParticles
		RecoParticles
		MCHits
		TDCData
		RawADCData
		CalibratedADCData
		RawLAPPDData
		CalibratedLAPPDData
		TriggerData
		MCFlag
		EventTime
		MCEventNum
		MCFile
		BeamStatus
	*/
	
	EventNumber=0;
	MCTriggernum=0;
	// pull the first entry to get the MCFile
	int nbytesread = WCSimEntry->GetEntry(MCEventNum);  // <0 if out of file
	if(nbytesread<=0){
		logmessage = "LoadWCSim Tool had no entry "+to_string(MCEventNum);
		if(nbytesread==-4){
			logmessage+=": Overran end of TChain! Have you specified more iterations than are available in ToolChainConfig?";
		} else if(nbytesread==0){
			logmessage+=": No TChain loaded! Is your filepath correct?";
		}
		Log(logmessage,v_error,verbosity);
		cerr<<"############################"<<endl;
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	MCFile = WCSimEntry->GetCurrentFile()->GetName();
	m_data->Stores.at("ANNIEEvent")->Set("MCFile",MCFile);
	
	// use nominal beam values TODO
	double beaminten=4.777e+12;
	double beampow=3.2545e+16;
	RunStartTime.SetNs(RunStartUser);
	BeamStatus = new BeamStatusClass(RunStartTime, beaminten, beampow, "stable");
	
	// Construct the other objects we'll be setting at event level,
	// pass managed pointers to the ANNIEEvent Store
	MCParticles = new std::vector<MCParticle>;
	MCHits = new std::map<unsigned long,std::vector<MCHit>>;
	TDCData = new std::map<unsigned long,std::vector<MCHit>>;
	EventTime = new TimeClass();
	TriggerClass beamtrigger("beam",5,0,true,0);
	TriggerData = new std::vector<TriggerClass>{beamtrigger}; // FIXME ? one trigger and resetting time is ok?
	
	// we'll put these in the CStore: so don't delete them in Finalise! It'll get handled by the Store
	ParticleId_to_TankTubeIds = new std::map<int,std::map<unsigned long,double>>;
	ParticleId_to_MrdTubeIds = new std::map<int,std::map<unsigned long,double>>;
	ParticleId_to_VetoTubeIds = new std::map<int,std::map<unsigned long,double>>;
	ParticleId_to_TankCharge = new std::map<int,double>;
	ParticleId_to_MrdCharge = new std::map<int,double>;
	ParticleId_to_VetoCharge = new std::map<int,double>;
	trackid_to_mcparticleindex = new std::map<int,int>;
	
	//anniegeom->GetChannel(0); // trigger InitChannelMap
	
	m_data->CStore.Set("UserEvent",false);   //enables the ability for other tools to select a specific event number
	triggers_event = 0;
	
	return true;
}


bool LoadWCSim::Execute(){
	
	// probably not necessary, clears the map for this entry. We're going to re-Set the event entry anyway...
	//m_data->Stores.at("ANNIEEvent")->Clear();
	
	//check if another tool has specified a specific evnumber to load (currently e.g. the EventDisplay has the ability to do that)
	bool user_event;
	m_data->CStore.Get("UserEvent",user_event);
	if (user_event){
		m_data->CStore.Set("UserEvent",false);
		MCTriggernum = 0;   //look at first trigger for user-specified event numbers
		int user_evnum; 
		uint16_t currentTriggernum;
		bool check_further_triggers=false;
		m_data->CStore.Get("LoadEvNr",user_evnum);
		m_data->CStore.Get("CheckFurtherTriggers",check_further_triggers);
		m_data->CStore.Get("CurrentTriggernum",currentTriggernum);
		MCEventNum = user_evnum;
		currentTriggernum++;
		std::cout <<"check_further_triggers = "<<check_further_triggers<<", currentTriggernum = "<<currentTriggernum<<std::endl;
		std::cout <<"Number of Events: "<<triggers_event<<std::endl;
		if (check_further_triggers){
			if (currentTriggernum!=triggers_event){
				//there is a further trigger in the previous event, load the entry
				MCEventNum = user_evnum - 1;
				MCTriggernum=currentTriggernum;
			}
		}
		// Pre-load entry so we can stop the loop if it this was the last one in the chain
		if(MCEventNum>=MaxEntries && MaxEntries>0){
			std::cout<<"LoadWCSim Tool: Reached max entries specified in config file, terminating ToolChain"<<endl;
			m_data->vars.Set("StopLoop",1);
		} else {
			int nbytesread = WCSimEntry->GetEntry(MCEventNum);  // <0 if out of file
			if (verbosity > v_debug) std::cout <<"LoadWCSim tool: Trying to get next event, MCEventNum: "<<MCEventNum<<", nbytesread: "<<nbytesread<<std::endl;
			if(nbytesread<=0){
				Log("LoadWCSim Tool: Reached last entry of WCSim input file, terminating ToolChain",v_warning,verbosity);
				m_data->vars.Set("StopLoop",1);
				return true;
			}
		}
	}
	if(verbosity) cout<<"Executing tool LoadWCSim with MC entry "<<MCEventNum<<", trigger "<<MCTriggernum<<endl;
	int loopstopped=0;
	get_ok = m_data->vars.Get("StopLoop",loopstopped);
	if(get_ok && loopstopped){
		// setting StopLoop doesn't terminate the ToolChain if the number of iterations
		// is specified manually in the ToolChainConfig.
		// This is almost certainly going to result in a segfault somewhere,
		// (e.g. if this tool set it in the last loop iteration because it ran out of entries)
		// but let's do what we can
		Log("WARNING: STOPLOOP HAS BEEN SET. RETURNING",v_error,verbosity);
		return 0;
	}
	MCFile = WCSimEntry->GetCurrentFile()->GetName();
	
	MCHits->clear();
	TDCData->clear();
	mrd_firstlayer=false;
	mrd_lastlayer=false;	

	triggers_event = WCSimEntry->wcsimrootevent->GetNumberOfEvents();
	
	//for(int MCTriggernum=0; MCTriggernum<WCSimEntry->wcsimrootevent->GetNumberOfEvents(); MCTriggernum++){
		if(verbosity>1) cout<<"getting triggers"<<endl;
		// cherenkovhit(times) are all in first trig
		firsttrigt=WCSimEntry->wcsimrootevent->GetTrigger(0);
		firsttrigm=WCSimEntry->wcsimrootevent_mrd->GetTrigger(0);
		firsttrigv=WCSimEntry->wcsimrootevent_facc->GetTrigger(0);
		atrigt = WCSimEntry->wcsimrootevent->GetTrigger(MCTriggernum);
		if(MCTriggernum<(WCSimEntry->wcsimrootevent_mrd->GetNumberOfEvents())){
			atrigm = WCSimEntry->wcsimrootevent_mrd->GetTrigger(MCTriggernum);
		} else { atrigm=nullptr; }
		if(MCTriggernum<(WCSimEntry->wcsimrootevent_facc->GetNumberOfEvents())){
			atrigv = WCSimEntry->wcsimrootevent_facc->GetTrigger(MCTriggernum);
		} else { atrigv=nullptr; }
		if(verbosity>2) cout<<"wcsimrootevent="<<WCSimEntry->wcsimrootevent<<endl;
		if(verbosity>2) cout<<"wcsimrootevent_mrd="<<WCSimEntry->wcsimrootevent_mrd<<endl;
		if(verbosity>2) cout<<"wcsimrootevent_facc="<<WCSimEntry->wcsimrootevent_facc<<endl;
		if(verbosity>2) cout<<"atrigt="<<atrigt<<", atrigm="<<atrigm<<", atrigv="<<atrigv<<endl;
		
		if(verbosity>1) cout<<"getting event date"<<endl;
		RunNumber = atrigt->GetHeader()->GetRun();
		SubrunNumber = 0;
		EventTimeNs = atrigt->GetHeader()->GetDate();
		EventTime->SetNs(EventTimeNs);
		if(verbosity>2) cout<<"EventTime is "<<EventTimeNs<<"ns"<<endl;
		
		// Load ALL MC particles (for all delayed MC triggers) only on MCTrigger 0
		if(MCTriggernum==0){
			MCParticles->clear();
			trackid_to_mcparticleindex->clear();
			primarymuonindex=-1;
			
			std::string geniefilename = firsttrigt->GetHeader()->GetGenieFileName().Data();
			int genieentry = firsttrigt->GetHeader()->GetGenieEntryNum();
			/*if(verbosity>3)*/ cout<<"Genie file is "<<geniefilename<<", genie event num was "<<genieentry<<endl;
			m_data->CStore.Set("GenieFile",geniefilename);
			m_data->CStore.Set("GenieEntry",genieentry);
			
			for(int trigi=0; trigi<WCSimEntry->wcsimrootevent->GetNumberOfEvents(); trigi++){
				
				WCSimRootTrigger* atrigtt = WCSimEntry->wcsimrootevent->GetTrigger(trigi);
				if(verbosity>1) cout<<"getting "<<atrigtt->GetNtrack()<<" tracks from trigger "<<trigi<<endl;
				for(int tracki=0; tracki<atrigtt->GetNtrack(); tracki++){
					if(verbosity>2) cout<<"getting track "<<tracki<<endl;
					WCSimRootTrack* nextrack = (WCSimRootTrack*)atrigtt->GetTracks()->At(tracki);
					/* a WCSimRootTrack has methods:
					Int_t     GetIpnu()             pdg
					Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
					Float_t   GetM()                mass
					Float_t   GetP()                momentum magnitude
					Float_t   GetE()                energy (inc rest mass^2)
					Float_t   GetEndE()             energy on stopping of particle tracking
					Float_t   GetEndP()             momentum on stopping of particle tracking
					Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
					Int_t     GetStopvol()          stopping volume: but these may not be set.
					Float_t   GetDir(Int_t i=0)     momentum unit vector
					Float_t   GetPdir(Int_t i=0)    momentum vector
					Float_t   GetPdirEnd(Int_t i=0) direction vector on stop tracking
					Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
					Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
					Int_t     GetParenttype()       parent pdg, 0 for primary.
					Float_t   GetTime()             trj->GetGlobalTime(); starting time of particle
					Float_t   GetStopTime()
					Int_t     GetId()               wcsim trackid
					*/
					
					tracktype startstoptype = tracktype::UNDEFINED;
					//MC particle times are relative to the trigger time
					if(nextrack->GetFlag()!=0) {
						if (nextrack->GetFlag()==-1){
							MCParticle neutrino(
							nextrack->GetIpnu(), nextrack->GetE(), nextrack->GetEndE(),
							Position(nextrack->GetStart(0) / 100.,
									 nextrack->GetStart(1) / 100.,
									 nextrack->GetStart(2) / 100.),
							Position(nextrack->GetStop(0) / 100.,
									 nextrack->GetStop(1) / 100.,
									 nextrack->GetStop(2) / 100.),
							//MC particle times now stored relative to the trigger time
							(static_cast<double>(nextrack->GetTime()-EventTimeNs)),
							(static_cast<double>(nextrack->GetStopTime()-EventTimeNs)),
							Direction(nextrack->GetDir(0), nextrack->GetDir(1), nextrack->GetDir(2)),
							(sqrt(pow(nextrack->GetStop(0)-nextrack->GetStart(0),2.)+
								 pow(nextrack->GetStop(1)-nextrack->GetStart(1),2.)+
								 pow(nextrack->GetStop(2)-nextrack->GetStart(2),2.))) / 100.,
							startstoptype,
							nextrack->GetId(),
							nextrack->GetParenttype(),
							nextrack->GetFlag(),
							trigi);

							//Set the neutrino as its own particle
							m_data->Stores["ANNIEEvent"]->Set("NeutrinoParticle",neutrino);
							
							continue; // flag 0 only is normal particles: excludes neutrino
						}
					}
					MCParticle thisparticle(
						nextrack->GetIpnu(), nextrack->GetE(), nextrack->GetEndE(),
						Position(nextrack->GetStart(0) / 100.,
								 nextrack->GetStart(1) / 100.,
								 nextrack->GetStart(2) / 100.),
						Position(nextrack->GetStop(0) / 100.,
								 nextrack->GetStop(1) / 100.,
								 nextrack->GetStop(2) / 100.),
						//MC particle times now stored relative to the trigger time
						(static_cast<double>(nextrack->GetTime()-EventTimeNs)),
						(static_cast<double>(nextrack->GetStopTime()-EventTimeNs)),
						Direction(nextrack->GetDir(0), nextrack->GetDir(1), nextrack->GetDir(2)),
						(sqrt(pow(nextrack->GetStop(0)-nextrack->GetStart(0),2.)+
							 pow(nextrack->GetStop(1)-nextrack->GetStart(1),2.)+
							 pow(nextrack->GetStop(2)-nextrack->GetStart(2),2.))) / 100.,
						startstoptype,
						nextrack->GetId(),
						nextrack->GetParenttype(),
						nextrack->GetFlag(),
						trigi);
					// not currently in constructor call, but we now have it in latest WCSim files
					// XXX this will fall over with older WCSim files, whose WCSimLib doesn't have this method!
					thisparticle.SetTankExitPoint(Position(nextrack->GetTankExitPoint(0)/ 100.,
														   nextrack->GetTankExitPoint(1)/ 100.,
														   nextrack->GetTankExitPoint(2)/ 100.));
					if( (nextrack->GetIpnu()==13) &&
						(nextrack->GetParenttype()==0) &&
						(nextrack->GetFlag()==0) &&
						(primarymuonindex<0) ){
							// call this the primary muon. If we have more than one, use the first
							primarymuonindex = MCParticles->size();
					}
					if((abs(nextrack->GetIpnu())==13)||
					   (abs(nextrack->GetIpnu())==211)||
					   (nextrack->GetIpnu()==111)){
							if (verbosity > 0) {
								std::cout<<"Found "<<nextrack->GetIpnu()<<" with parent pdg "
								<<nextrack->GetParenttype()<<", flag "<<nextrack->GetFlag()
								<<" track id "<<nextrack->GetId()
								<< ", start vertex (" + to_string(nextrack->GetStart(0)/100.)
								<< ", " + to_string(nextrack->GetStart(1)/100.)
								<< ", " + to_string(nextrack->GetStart(2)/100.)
								<< "), and end vertex (" + to_string(nextrack->GetStop(0)/100.)
								<< ", " + to_string(nextrack->GetStop(1)/100.)
								<< ", " + to_string(nextrack->GetStop(2)/100.)
								<< ")"
								<<std::endl;
							}
					}
/*
					// Print primary muons or the first track
					if(((nextrack->GetIpnu()==13) && (nextrack->GetParenttype()==0))||
					   (MCParticles->size()==0)){
						if (verbosity) std::cout<<"Found "<<nextrack->GetIpnu()<<" with parent pdg "
								<<nextrack->GetParenttype()<<", flag "<<nextrack->GetFlag()
								<<" track id "<<nextrack->GetId()
								<<" at position "<<MCParticles->size()<<std::endl;
					}
					// print pions
					if((abs(nextrack->GetIpnu())==211)||
					   (nextrack->GetIpnu()==111)){
						if (verbosity) std::cout<<"Found "<<nextrack->GetIpnu()<<" with parent pdg "
								<<nextrack->GetParenttype()<<", flag "<<nextrack->GetFlag()
								<<" track id "<<nextrack->GetId()
								<< ", start vertex (" + to_string(nextrack->GetStart(0)/100.)
								<< ", " + to_string(nextrack->GetStart(1)/100.)
								<< ", " + to_string(nextrack->GetStart(2)/100.)
								<< "), and end vertex (" + to_string(nextrack->GetStop(0)/100.)
								<< ", " + to_string(nextrack->GetStop(1)/100.)
								<< ", " + to_string(nextrack->GetStop(2)/100.)
								<< ")"
								<<std::endl;
					}
					// print muons or gammas from pion decays
					if(((abs(nextrack->GetIpnu())==13)||(nextrack->GetIpnu()==22)) &&
					   ((abs(nextrack->GetParenttype())==211)||(nextrack->GetParenttype()==111))){
						if (verbosity) std::cout<<"Found "<<nextrack->GetIpnu()<<" with parent pdg "
								<<nextrack->GetParenttype()<<", flag "<<nextrack->GetFlag()
								<<" track id "<<nextrack->GetId()
								<<" at position "<<MCParticles->size()<<std::endl;
					}
*/
					if(nextrack->GetIpnu()==13){
						logmessage = "Muon found with flag: "+to_string(nextrack->GetFlag())
							+ ", parent type " + to_string(nextrack->GetParenttype())
							+ ", Id " + to_string(nextrack->GetId())
							+ ", start vertex (" + to_string(nextrack->GetStart(0)/100.)
							+ ", " + to_string(nextrack->GetStart(1)/100.)
							+ ", " + to_string(nextrack->GetStart(2)/100.)
							+ "), and end vertex (" + to_string(nextrack->GetStop(0)/100.)
							+ ", " + to_string(nextrack->GetStop(1)/100.)
							+ ", " + to_string(nextrack->GetStop(2)/100.)
							+ ")";
						Log(logmessage,v_debug,verbosity);
					}
					
					trackid_to_mcparticleindex->emplace(nextrack->GetId(),MCParticles->size());
					MCParticles->push_back(thisparticle);
				}
				if(verbosity>2) cout<<"MCParticles has "<<MCParticles->size()<<" entries"<<endl;
				
			}  // loop over loading particles from all MC triggers on first MC trigger
		} else { 
			// if MCTrigger>0, since particle times are relative to the trigger time,
			// we need to update all the particle times
			double timediff = EventTimeNs - firsttrigt->GetHeader()->GetDate();
			for(MCParticle& aparticle : *MCParticles){
				aparticle.SetStartTime(aparticle.GetStartTime()-timediff);
				aparticle.SetStopTime (aparticle.GetStopTime() -timediff);
			}
		} // end updating particle times
		
		int numtankdigits = atrigt ? atrigt->GetCherenkovDigiHits()->GetEntries() : 0;
		if(verbosity>1) cout<<"looping over "<<numtankdigits<<" tank digits"<<endl;
		for(int digiti=0; digiti<numtankdigits; digiti++){
			if(verbosity>2) cout<<"getting digit "<<digiti<<endl;
			WCSimRootCherenkovDigiHit* digihit =
				(WCSimRootCherenkovDigiHit*)atrigt->GetCherenkovDigiHits()->At(digiti);
			//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ(), GetPhotonIds()
			if(verbosity>2) cout<<"next digihit at "<<digihit<<endl;
			int tubeid = digihit->GetTubeId();  // geometry TubeID->channelkey map is made INCLUDING offset of 1
			if(verbosity>2) cout<<"tubeid="<<tubeid<<endl;
			if(pmt_tubeid_to_channelkey.count(tubeid)==0){
				cerr<<"LoadWCSim ERROR: tank PMT with no associated ChannelKey!"<<endl;
				return false;
			}
			unsigned long key = pmt_tubeid_to_channelkey.at(tubeid);
			if(verbosity>2) cout<<"ChannelKey="<<key<<endl;
			double digittime;
			if(use_smeared_digit_time){
				digittime = static_cast<double>(digihit->GetT()-HistoricTriggeroffset); // relative to trigger
			} else {
				// instead take the true time of the first photon
				std::vector<int> photonids = digihit->GetPhotonIds();   // indices of the digit's photons
				double earliestphotontruetime=999999999999;
				for(int& aphotonindex : photonids){
					WCSimRootCherenkovHitTime* thehittimeobject =
						 (WCSimRootCherenkovHitTime*)firsttrigt->GetCherenkovHitTimes()->At(aphotonindex);
					if(thehittimeobject==nullptr){
						cerr<<"LoadWCSim Tool: ERROR! Retrieval of photon from digit returned nullptr!"<<endl;
						continue;
					}
					double aphotontime = static_cast<double>(thehittimeobject->GetTruetime());
					if(aphotontime<earliestphotontruetime){ earliestphotontruetime = aphotontime; }
				}
				digittime = earliestphotontruetime;
			}
			if(verbosity>2){ cout<<"digittime is "<<digittime<<" [ns] from Trigger"<<endl; }
			float digiq = digihit->GetQ();
			if(verbosity>2) cout<<"digit Q is "<<digiq<<endl;
			// Get hit parent information
			std::vector<int> parents = GetHitParentIds(digihit, firsttrigt);
			
			MCHit nexthit(key, digittime, digiq, parents);
			if(MCHits->count(key)==0) MCHits->emplace(key, std::vector<MCHit>{nexthit});
			else MCHits->at(key).push_back(nexthit);
			if(verbosity>2) cout<<"digit added"<<endl;
		}
		if(verbosity>2) cout<<"done with tank digits"<<endl;
		
		//MRD Hits
		int nummrddigits = atrigm ? atrigm->GetCherenkovDigiHits()->GetEntries() : 0;
		if(verbosity>1) cout<<"adding "<<nummrddigits<<" mrd digits"<<endl;
		for(int digiti=0; digiti<nummrddigits; digiti++){
			if(verbosity>2) cout<<"getting digit "<<digiti<<endl;
			WCSimRootCherenkovDigiHit* digihit =
				(WCSimRootCherenkovDigiHit*)atrigm->GetCherenkovDigiHits()->At(digiti);
			if(verbosity>2) cout<<"next digihit at "<<digihit<<endl;
			int tubeid = digihit->GetTubeId();
			if(verbosity>2) cout<<"tubeid="<<tubeid<<endl;
			if(mrd_tubeid_to_channelkey.count(tubeid)==0){
				cerr<<"LoadWCSim ERROR: MRD PMT with no associated ChannelKey!"<<endl;
				return false;
			}
			unsigned long key = mrd_tubeid_to_channelkey.at(tubeid);
			double digittime;
			if(use_smeared_digit_time){
				digittime = static_cast<double>(digihit->GetT()-HistoricTriggeroffset); // relative to trigger
			} else {
				// instead take the true time of the first photon
				std::vector<int> photonids = digihit->GetPhotonIds();   // indices of the digit's photons
				double earliestphotontruetime=999999999999;
				for(int& aphotonindex : photonids){
					WCSimRootCherenkovHitTime* thehittimeobject =
						 (WCSimRootCherenkovHitTime*)firsttrigm->GetCherenkovHitTimes()->At(aphotonindex);
					if(thehittimeobject==nullptr){
						cerr<<"LoadWCSim Tool: ERROR! Retrieval of photon from digit returned nullptr!"<<endl;
						continue;
					}
					double aphotontime = static_cast<double>(thehittimeobject->GetTruetime());
					if(aphotontime<earliestphotontruetime){ earliestphotontruetime = aphotontime; }
				}
				digittime = earliestphotontruetime;
			}
			if(verbosity>2){ cout<<"digittime is "<<digittime<<" [ns] from Trigger"<<endl; }
			float digiq = digihit->GetQ();
			if(verbosity>2) cout<<"digit Q is "<<digiq<<endl;
			// Get hit parent information
			std::vector<int> parents = GetHitParentIds(digihit, firsttrigm);
			
			MCHit nexthit(key, digittime, digiq, parents);
			if(TDCData->count(key)==0) {TDCData->emplace(key, std::vector<MCHit>{nexthit}); if (Mrd_Chankey_Layer.at(key)==0) mrd_firstlayer=true; if (Mrd_Chankey_Layer.at(key)==10) mrd_lastlayer=true;}
			else TDCData->at(key).push_back(nexthit);
			if(verbosity>2) cout<<"digit added"<<endl;
		}
		if(verbosity>2) cout<<"done with mrd digits"<<endl;
		
		// Veto Hits
		int numvetodigits = atrigv ? atrigv->GetCherenkovDigiHits()->GetEntries() : 0;
		if(verbosity>1) cout<<"adding "<<numvetodigits<<" veto digits"<<endl;
		for(int digiti=0; digiti<numvetodigits; digiti++){
			if(verbosity>2) cout<<"getting digit "<<digiti<<endl;
			WCSimRootCherenkovDigiHit* digihit =
				(WCSimRootCherenkovDigiHit*)atrigv->GetCherenkovDigiHits()->At(digiti);
			if(verbosity>2) cout<<"next digihit at "<<digihit<<endl;
			int tubeid = digihit->GetTubeId();
			if(verbosity>2) cout<<"tubeid="<<tubeid<<endl;
			if(facc_tubeid_to_channelkey.count(tubeid)==0){
				cerr<<"LoadWCSim ERROR: FACC PMT with no associated ChannelKey!"<<endl;
				return false;
			}
			unsigned int key = facc_tubeid_to_channelkey.at(tubeid);
			double digittime;
			if(use_smeared_digit_time){
				digittime = static_cast<double>(digihit->GetT()-HistoricTriggeroffset); // relative to trigger
			} else {
				// instead take the true time of the first photon
				std::vector<int> photonids = digihit->GetPhotonIds();   // indices of the digit's photons
				double earliestphotontruetime=999999999999;
				for(int& aphotonindex : photonids){
					WCSimRootCherenkovHitTime* thehittimeobject =
						 (WCSimRootCherenkovHitTime*)firsttrigv->GetCherenkovHitTimes()->At(aphotonindex);
					if(thehittimeobject==nullptr){
						cerr<<"LoadWCSim Tool: ERROR! Retrieval of photon from digit returned nullptr!"<<endl;
						continue;
					}
					double aphotontime = static_cast<double>(thehittimeobject->GetTruetime());
					if(aphotontime<earliestphotontruetime){ earliestphotontruetime = aphotontime; }
				}
				digittime = earliestphotontruetime;
			}
			if(verbosity>2){ cout<<"digittime is "<<digittime<<" [ns] from Trigger"<<endl; }
			float digiq = digihit->GetQ();
			if(verbosity>2) cout<<"digit Q is "<<digiq<<endl;
			// Get hit parent information
			std::vector<int> parents = GetHitParentIds(digihit, firsttrigv);
			
			MCHit nexthit(key, digittime, digiq, parents);
			if(TDCData->count(key)==0) TDCData->emplace(key, std::vector<MCHit>{nexthit});
			else TDCData->at(key).push_back(nexthit);
			if(verbosity>2) cout<<"digit added"<<endl;
		}
		if(verbosity>2) cout<<"done with veto digits"<<endl;
		
		if(verbosity>2) cout<<"setting triggerdata time to "<<EventTimeNs<<"ns"<<endl;
		TriggerData->front().SetTime(EventTimeNs);
		
		// update the information about tracks and which tank/mrd/veto PMTs they hit
		// this needs updating with each MC trigger, as digits are grouped into MC trigger
		// so these maps will then only contain the digits respective particles create
		// in the active trigger
		// 
		// ParticleId_to_TankTubeIds is a std::map<ParticleId,std::map<ChannelKey,TotalCharge>>
		// where TotalCharge is the total charge from that particle on that tube
		// (in the event that the particle generated several hits on the tube)
		// ParticleId_to_TankCharge is a std::map<ParticleId,TotalCharge> 
		// where TotalCharge is summed over all digits, on all pmts, which contained
		// light from that particle
		ParticleId_to_TankTubeIds->clear();
		ParticleId_to_MrdTubeIds->clear();
		ParticleId_to_VetoTubeIds->clear();
		ParticleId_to_TankCharge->clear();
		ParticleId_to_MrdCharge->clear();
		ParticleId_to_VetoCharge->clear();
		MakeParticleToPmtMap(atrigt, firsttrigt, ParticleId_to_TankTubeIds, ParticleId_to_TankCharge, pmt_tubeid_to_channelkey);
		MakeParticleToPmtMap(atrigm, firsttrigm, ParticleId_to_MrdTubeIds, ParticleId_to_MrdCharge, mrd_tubeid_to_channelkey);
		MakeParticleToPmtMap(atrigv, firsttrigv, ParticleId_to_VetoTubeIds, ParticleId_to_VetoCharge, facc_tubeid_to_channelkey);
		
	//}
	
	//int mrdentries;
	//m_data->Stores.at("TDCData")->Get("TotalEntries",mrdentries); // ??
//	m_data->Stores("WCSimEntries")->Set("wcsimrootevent",WCSimEntry->wcsimrootevent);
//	m_data->Stores("WCSimEntries")->Set("wcsimrootevent_mrd",WCSimEntry->wcsimrootevent_mrd);
//	m_data->Stores("WCSimEntries")->Set("wcsimrootevent_facc",*(WCSimEntry->wcsimrootevent_facc));
	
	
	// set event level variables
	if(verbosity>1) cout<<"setting the store variables"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("RunNumber",RunNumber);
	m_data->Stores.at("ANNIEEvent")->Set("SubrunNumber",SubrunNumber);
	m_data->Stores.at("ANNIEEvent")->Set("EventNumber",EventNumber);
	if(verbosity>2) cout<<"particles"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("MCParticles",MCParticles,true);
	if(verbosity>2) cout<<"hits"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("MCHits",MCHits,true);
	if(verbosity>2) cout<<"tdcdata"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("TDCData",TDCData,true);
	// TODO?
	// right now we have three Time variables:
	// 1. "RunStartTime" which stores just a user passed start time.
	// 2. "TriggerData", a one-element vector of TriggerClass objects, each of which has a time member.
	//    The one used entry has its time member set to the time between MC event start (always 0)
	//     and the MC trigger time.
	// 3. "EventTime" which stores this same time difference between MC Trigger and MC event start
	// This means we simulate many events all happening ~ the unix epoch (time 0).
	// If we want to simulate a 'run' of events, we need to throw some cumulative running time and
	// add this to RunStartTime, and store the Event and Trigger times separately.
	// This is done by PulseSimulation tool (timefileout_Time), but maybe should be here...
	if(verbosity>2) cout<<"triggerdata"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("TriggerData",TriggerData,true);
	if(verbosity>2) cout<<"eventtime"<<endl;
	m_data->Stores.at("ANNIEEvent")->Set("EventTime",EventTime,true);
	m_data->Stores.at("ANNIEEvent")->Set("MCEventNum",MCEventNum);
	m_data->Stores.at("ANNIEEvent")->Set("MCTriggernum",MCTriggernum);
	m_data->Stores.at("ANNIEEvent")->Set("MCFile",MCFile);
	m_data->Stores.at("ANNIEEvent")->Set("MCFlag",true);                   // constant
	m_data->Stores.at("ANNIEEvent")->Set("BeamStatus",BeamStatus,true);
	m_data->CStore.Set("NumTriggersThisMCEvt",WCSimEntry->wcsimrootevent->GetNumberOfEvents());
	// auxilliary information about MC Truth particles
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_TankTubeIds", ParticleId_to_TankTubeIds, false);
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_TankCharge", ParticleId_to_TankCharge, false);
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_MrdTubeIds", ParticleId_to_MrdTubeIds, false);
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_MrdCharge", ParticleId_to_MrdCharge, false);
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_VetoTubeIds", ParticleId_to_VetoTubeIds, false);
	m_data->Stores.at("ANNIEEvent")->Set("ParticleId_to_VetoCharge", ParticleId_to_VetoCharge, false);
	m_data->Stores.at("ANNIEEvent")->Set("TrackId_to_MCParticleIndex",trackid_to_mcparticleindex,false);
	//Change MRD Triggertype if first and last layer saw a hit (Hardware Cosmic Trigger)
	if (mrd_lastlayer && mrd_firstlayer) Triggertype = "Cosmic";
	m_data->Stores.at("ANNIEEvent")->Set("MRDTriggerType",Triggertype);
	m_data->Stores.at("ANNIEEvent")->Set("PrimaryMuonIndex",primarymuonindex);
	
	//Things that need to be set by later tools:
	//RawADCData
	//CalibratedADCData
	//RawLAPPDData
	//CalibratedLAPPDData
	//RecoParticles
	if(verbosity>1) cout<<"done loading event"<<endl;
	
	EventNumber++;
	MCTriggernum++;
	if(verbosity>2) cout<<"checking if we're done on trigs in this event"<<endl;
	bool newentry=false;
	if(MCTriggernum==WCSimEntry->wcsimrootevent->GetNumberOfEvents()){
		MCTriggernum=0;
		MCEventNum++;
		newentry=true;
		if(verbosity>2) cout<<"this is the last trigger in the event: next loop will process a new event"<<endl;
	} else {
		if(verbosity>2) cout<<"there are further triggers in this event: next loop will process the trigger "<<MCTriggernum<<"/"<<WCSimEntry->wcsimrootevent->GetNumberOfEvents()<<endl;
	}
	
	// Pre-load next entry so we can stop the loop if it this was the last one in the chain
	if(newentry){  // if next loop is processing the next trigger in the same entry, no need to re-load it
		if(MCEventNum>=MaxEntries && MaxEntries>0){
			cout<<"LoadWCSim Tool: Reached max entries specified in config file, terminating ToolChain"<<endl;
			m_data->vars.Set("StopLoop",1);
		} else {
			int nbytesread = WCSimEntry->GetEntry(MCEventNum);  // <0 if out of file
			std::cout <<"Trying to get next event, MCEventNum: "<<MCEventNum<<", nbytesread: "<<nbytesread<<std::endl;
			if(nbytesread<=0){
				Log("LoadWCSim Tool: Reached last entry of WCSim input file, terminating ToolChain",v_warning,verbosity);
				m_data->vars.Set("StopLoop",1);
			}
		}
	}
	
	//gObjectTable->Print();
	return true;
}

bool LoadWCSim::Finalise(){
	WCSimEntry->GetCurrentFile()->Close();
	delete WCSimEntry;
	
	// any pointers put in Stores to objects we do not want the Store to clean up
	// must be nullified before in finalise to prevent double free
	// can't just put 0 or nullptr directly as type must be recognisable as a pointer
//	std::map<int,std::map<unsigned long,double>>* ParticleId_to_TankTubeIds_nullptr = nullptr;
//	std::map<int,std::map<unsigned long,double>>* ParticleId_to_MrdTubeIds_nullptr = nullptr;
//	std::map<int,std::map<unsigned long,double>>* ParticleId_to_VetoTubeIds_nullptr = nullptr;
//	std::map<int,double>* ParticleId_to_TankCharge_nullptr = nullptr;
//	std::map<int,double>* ParticleId_to_MrdCharge_nullptr = nullptr;
//	std::map<int,double>* ParticleId_to_VetoCharge_nullptr = nullptr;
	
//	m_data->CStore.Set("ParticleId_to_TankTubeIds", ParticleId_to_TankTubeIds_nullptr, false);
//	m_data->CStore.Set("ParticleId_to_TankCharge", ParticleId_to_TankCharge_nullptr, false);
//	m_data->CStore.Set("ParticleId_to_MrdTubeIds", ParticleId_to_MrdTubeIds_nullptr, false);
//	m_data->CStore.Set("ParticleId_to_MrdCharge", ParticleId_to_MrdCharge_nullptr, false);
//	m_data->CStore.Set("ParticleId_to_VetoTubeIds", ParticleId_to_VetoTubeIds_nullptr, false);
//	m_data->CStore.Set("ParticleId_to_VetoCharge", ParticleId_to_VetoCharge_nullptr, false);
	
	return true;
}

Geometry* LoadWCSim::ConstructToolChainGeometry(){
	// Pull details from the WCSimRootGeom
	// ===================================
	double WCSimGeometryVer = 1;                      // TODO: pull this from some suitable variable
	// some variables are available from wcsimrootgeom.  TODO: put everything into wcsimrootgeom
	numtankpmts = wcsimrootgeom->GetWCNumPMT();
	numlappds = wcsimrootgeom->GetWCNumLAPPD();
	nummrdpmts = wcsimrootgeom->GetWCNumMRDPMT();
	numvetopmts = wcsimrootgeom->GetWCNumFACCPMT();
	double tank_xcentre = (wcsimrootgeom->GetWCOffset(0)) / 100.;  // convert [cm] to [m]
	double tank_ycentre = (wcsimrootgeom->GetWCOffset(1)) / 100.;
	double tank_zcentre = (wcsimrootgeom->GetWCOffset(2)) / 100.;
	Position tank_centre(tank_xcentre, tank_ycentre, tank_zcentre);
	double tank_radius = (wcsimrootgeom->GetWCCylRadius()) / 100.;
	double tank_halfheight = (wcsimrootgeom->GetWCCylLength()) / 100.;
	//Currently hard-coded; estimated with a tape measure on the ANNIE frame :)
	double pmt_enclosed_radius = 1.0;
	double pmt_enclosed_halfheight = 1.45;
	// geometry variables not yet in wcsimrootgeom are in MRDSpecs.hh
	double mrd_width  =  (MRDSpecs::MRD_width)  / 100.;
	double mrd_height =  (MRDSpecs::MRD_height) / 100.;
	double mrd_depth  =  (MRDSpecs::MRD_depth)  / 100.;
	double mrd_start  =  (MRDSpecs::MRD_start)  / 100.;
	if(verbosity>1) cout<<"we have "<<numtankpmts<<" tank pmts, "<<nummrdpmts
					  <<" mrd pmts and "<<numlappds<<" lappds"<<endl;
	
	// construct the ToolChain Goemetry
	// ================================
	Geometry* anniegeom = new Geometry(WCSimGeometryVer,
									   tank_centre,
									   tank_radius,
									   tank_halfheight,
									   pmt_enclosed_radius,
									   pmt_enclosed_halfheight,
									   mrd_width,
									   mrd_height,
									   mrd_depth,
									   mrd_start,
									   numtankpmts,
									   nummrdpmts,
									   numvetopmts,
									   numlappds,
									   geostatus::FULLY_OPERATIONAL);
	if(verbosity>1){
		cout<<"constructed anniegeom at "<<anniegeom<<" with tank origin "; tank_centre.Print();
	}
	m_data->Stores.at("ANNIEEvent")->Header->Set("AnnieGeometry",anniegeom,true);
	
	// Construct the Detectors and Channels
	// ====================================
	// PMTs
	unsigned int ADC_Crate_Num = 0;
	unsigned int ADC_Card_Num  = 0;
	unsigned int ADC_Chan_Num = 0;
	unsigned int MT_Crate_Num = 0;
	unsigned int MT_Card_Num = 0;
	unsigned int MT_Chan_Num = 0;
	// LAPPDs
	unsigned int ACDC_Crate_Num = 0;
	unsigned int ACDC_Card_Num = 0;
	unsigned int ACDC_Chan_Num = 0;
	unsigned int ACC_Crate_Num = 0;
	unsigned int ACC_Card_Num = 0;
	unsigned int ACC_Chan_Num = 0;
	// TDCs
	unsigned int TDC_Crate_Num = 0;
	unsigned int TDC_Card_Num = 0;
	unsigned int TDC_Chan_Num = 0;
	// HV
	unsigned int CAEN_HV_Crate_Num = 0;
	unsigned int CAEN_HV_Card_Num = 0;
	unsigned int CAEN_HV_Chan_Num = 0;
	unsigned int LeCroy_HV_Crate_Num = 0;
	unsigned int LeCroy_HV_Card_Num = 0;
	unsigned int LeCroy_HV_Chan_Num = 0;
	unsigned int LAPPD_HV_Crate_Num = 0;
	unsigned int LAPPD_HV_Card_Num = 0;
	unsigned int LAPPD_HV_Chan_Num = 0;
	
	// lappds
	for(int lappdi=0; lappdi<numlappds; lappdi++){
		WCSimRootPMT anlappd = wcsimrootgeom->GetLAPPD(lappdi);
		
		// Construct the detector associated with this tile
		unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
		lappd_tubeid_to_detectorkey.emplace(anlappd.GetTubeNo(),uniquedetectorkey);
		detectorkey_to_lappdid.emplace(uniquedetectorkey,anlappd.GetTubeNo());
		std::string CylLocString;
		switch (anlappd.GetCylLoc()){
			case 0:  CylLocString = "TopCap";    break;
			case 2:  CylLocString = "BottomCap"; break;
			case 1:  CylLocString = "Barrel";    break;
			case 4:  CylLocString = "MRD";       break;  // TODO set this as H or V paddle? And layer?
			case 5:  CylLocString = "Veto";      break;  // TODO set layer?
			default: CylLocString = "NA";        break;  // unknown
		}
		Detector adet(uniquedetectorkey,
					  "LAPPD",
					  CylLocString,
					  Position( anlappd.GetPosition(0)/100.,
					            anlappd.GetPosition(1)/100.,
					            anlappd.GetPosition(2)/100.),
					  Direction(anlappd.GetOrientation(0),
					            anlappd.GetOrientation(1),
					            anlappd.GetOrientation(2)),
					  anlappd.GetName(),
					  detectorstatus::ON,
					  0.);
		
		// construct all the channels associated with this LAPPD
		for(int stripi=0; stripi<LappdNumStrips; stripi++){
			unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
			
			int stripside = ((stripi%2)==0);   // StripSide=0 for LHS (x<0), StripSide=1 for RHS (x>0)
			int stripnum = (int)(stripi/2);    // Strip number: add 2 channels per strip as we go
			double xpos = (stripside) ? -LappdStripLength : LappdStripLength;
			double ypos = (stripnum*LappdStripSeparation) - ((LappdNumStrips*LappdStripSeparation)/2.);
			
			// fill up ADC cards and channels monotonically, they're arbitrary for simulation
			ACDC_Chan_Num++;
			if(ACDC_Chan_Num>=ACDC_CHANNELS_PER_CARD)  { ACDC_Chan_Num=0; ACDC_Card_Num++; ACC_Chan_Num++; }
			if(ACDC_Card_Num>=ACDC_CARDS_PER_CRATE)    { ACDC_Card_Num=0; ACDC_Crate_Num++; }
			if(ACC_Chan_Num>=ACC_CHANNELS_PER_CARD)    { ACC_Chan_Num=0; ACC_Card_Num++;    }
			if(ACC_Card_Num>=ACC_CARDS_PER_CRATE)      { ACC_Card_Num=0; ACC_Crate_Num++;   }
			// same for HV
			LAPPD_HV_Chan_Num++;
			if(LAPPD_HV_Chan_Num>=LAPPD_HV_CHANNELS_PER_CARD)    { LAPPD_HV_Chan_Num=0; LAPPD_HV_Card_Num++;  }
			if(LAPPD_HV_Card_Num>=LAPPD_HV_CARDS_PER_CRATE) { LAPPD_HV_Card_Num=0; LAPPD_HV_Crate_Num++; }
			
			Channel lappdchannel(uniquechannelkey,
								 Position(xpos,ypos,0.),
								 stripside,
								 stripnum,
								 ACDC_Crate_Num,
								 ACDC_Card_Num,
								 ACDC_Chan_Num,
								 ACC_Crate_Num,
								 ACC_Card_Num,
								 ACC_Chan_Num,
								 LAPPD_HV_Crate_Num,
								 LAPPD_HV_Card_Num,
								 LAPPD_HV_Chan_Num,
								 channelstatus::ON);
			
			// Add this channel to the geometry
			if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
			adet.AddChannel(lappdchannel);
		}
		if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
		// Add this detector to the geometry
		anniegeom->AddDetector(adet);
		if(verbosity>4) cout<<"printing geometry"<<endl;
		if(verbosity>4) anniegeom->PrintChannels();
	}
	
	// tank PMTs
	for(int pmti=0; pmti<numtankpmts; pmti++){
		WCSimRootPMT apmt = wcsimrootgeom->GetPMT(pmti);
		
		// Construct the detector associated with this PMT
		unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
		std::string CylLocString;
		int cylloc = apmt.GetCylLoc();
		if(apmt.GetName().find("EMI9954KB")!= std::string::npos){ cylloc=6; }     // FIXME set OD pmts in WCSim
		switch (cylloc){
			case 0:  CylLocString = "TopCap";    break;
			case 2:  CylLocString = "BottomCap"; break;
			case 1:  CylLocString = "Barrel";    break;
			case 4:  CylLocString = "MRD";       break;  // TODO set this as H or V paddle? And layer?
			case 5:  CylLocString = "Veto";      break;  // TODO set layer?
			case 6:  CylLocString = "OD";        break;
			default: CylLocString = "NA";        break;  // unknown
		}
		Detector adet(uniquedetectorkey,
					  "Tank",
					  CylLocString,
					  Position( apmt.GetPosition(0)/100.,
					            apmt.GetPosition(1)/100.,
					            apmt.GetPosition(2)/100.),
					  Direction(apmt.GetOrientation(0),
					            apmt.GetOrientation(1),
					            apmt.GetOrientation(2)),
					  apmt.GetName(),
					  detectorstatus::ON,
					  0.);
		
		// construct the channel associated with this PMT
		unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
		pmt_tubeid_to_channelkey.emplace(apmt.GetTubeNo(), uniquechannelkey);
		channelkey_to_pmtid.emplace(uniquechannelkey,apmt.GetTubeNo());
		
		// fill up ADC cards and channels monotonically, they're arbitrary for simulation
		ADC_Chan_Num++;
		if(ADC_Chan_Num>=ADC_CHANNELS_PER_CARD)  { ADC_Chan_Num=0; ADC_Card_Num++; MT_Chan_Num++; }
		if(ADC_Card_Num>=ADC_CARDS_PER_CRATE)    { ADC_Card_Num=0; ADC_Crate_Num++; }
		if(MT_Chan_Num>=MT_CHANNELS_PER_CARD)    { MT_Chan_Num=0; MT_Card_Num++; }
		if(MT_Card_Num>=MT_CARDS_PER_CRATE)      { MT_Card_Num=0; MT_Crate_Num++; }
		// same for HV
		CAEN_HV_Chan_Num++;
		if(CAEN_HV_Chan_Num>=CAEN_HV_CHANNELS_PER_CARD)    { CAEN_HV_Chan_Num=0; CAEN_HV_Card_Num++;  }
		if(CAEN_HV_Card_Num>=CAEN_HV_CARDS_PER_CRATE) { CAEN_HV_Card_Num=0; CAEN_HV_Crate_Num++; }
		
		Channel pmtchannel( uniquechannelkey,
							Position(0,0,0.),
							0, // stripside
							0, // stripnum
							ADC_Crate_Num,
							ADC_Card_Num,
							ADC_Chan_Num,
							MT_Crate_Num,
							MT_Card_Num,
							MT_Chan_Num,
							CAEN_HV_Crate_Num,
							CAEN_HV_Card_Num,
							CAEN_HV_Chan_Num,
							channelstatus::ON);
		
		// Add this channel to the geometry
		if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
		adet.AddChannel(pmtchannel);
		
		// Add this detector to the geometry
		if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
		anniegeom->AddDetector(adet);
		if(verbosity>4) cout<<"printing geometry"<<endl;
		if(verbosity>4) anniegeom->PrintChannels();
	}
	
	// mrd PMTs
	for(int mrdpmti=0; mrdpmti<nummrdpmts; mrdpmti++){
		WCSimRootPMT apmt = wcsimrootgeom->GetMRDPMT(mrdpmti);
		
		// Construct the detector associated with this PMT
		unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
		std::string CylLocString;
		switch (apmt.GetCylLoc()){
			case 0:  CylLocString = "TopCap";    break;
			case 2:  CylLocString = "BottomCap"; break;
			case 1:  CylLocString = "Barrel";    break;
			case 4:  CylLocString = "MRD";       break;  // TODO set this as H or V paddle? And layer?
			case 5:  CylLocString = "Veto";      break;  // TODO set layer?
			default: CylLocString = "NA";        break;  // unknown
		}
		Detector adet(uniquedetectorkey,
					  "MRD",
					  CylLocString,
					  Position( apmt.GetPosition(0)/100.,
					            apmt.GetPosition(1)/100.,
					            apmt.GetPosition(2)/100.),
					  Direction(apmt.GetOrientation(0),
					            apmt.GetOrientation(1),
					            apmt.GetOrientation(2)),
					  apmt.GetName(),
					  detectorstatus::ON,
					  0.);
		
		// construct the channel associated with this PMT
		unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
		mrd_tubeid_to_channelkey.emplace(apmt.GetTubeNo(), uniquechannelkey);
		channelkey_to_mrdpmtid.emplace(uniquechannelkey, apmt.GetTubeNo());
		
		// fill up TDC cards and channels monotonically, they're arbitrary for simulation
		TDC_Chan_Num++;
		if(TDC_Chan_Num>=TDC_CHANNELS_PER_CARD)  { TDC_Chan_Num=0; TDC_Card_Num++;  }
		if(TDC_Card_Num>=TDC_CARDS_PER_CRATE)    { TDC_Card_Num=0; TDC_Crate_Num++; }
		// same for HV
		LeCroy_HV_Chan_Num++;
		if(LeCroy_HV_Chan_Num>=LECROY_HV_CHANNELS_PER_CARD)    { LeCroy_HV_Chan_Num=0; LeCroy_HV_Card_Num++;  }
		if(LeCroy_HV_Card_Num>=LECROY_HV_CARDS_PER_CRATE) { LeCroy_HV_Card_Num=0; LeCroy_HV_Crate_Num++; }
		
		Channel pmtchannel( uniquechannelkey,
							Position(0,0,0.),
							0, // stripside
							0, // stripnum
							TDC_Crate_Num,
							TDC_Card_Num,
							TDC_Chan_Num,
							-1,                 // TDC has no level 2 signal handling
							-1,
							-1,
							LeCroy_HV_Crate_Num,
							LeCroy_HV_Card_Num,
							LeCroy_HV_Chan_Num,
							channelstatus::ON);
		
		// Add this channel to the geometry
		if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
		adet.AddChannel(pmtchannel);
		
		// Add this detector to the geometry
		if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
		anniegeom->AddDetector(adet);
		
		// Create a paddle
		// FIXME remove dependency on MRDSpecs
		if(mrdpmti!=(apmt.GetTubeNo()-1)){
			cerr<<"LoadWCSim: mrdpmti!=pmt.GetTubeNo-1!"<<std::endl;
			cerr<<"mrdpmti="<<mrdpmti<<", pmt.GetTubeNo="<<apmt.GetTubeNo()<<endl;
			assert(false);
		}
		// calculate MRD_x_y_z ... MRDSpecs doesn't provide a nice way to do this
		int layernum=0;
		while ((mrdpmti+1) > MRDSpecs::layeroffsets.at(layernum+1)){ layernum++; }
		Mrd_Chankey_Layer.emplace(uniquechannelkey,layernum);
		int in_layer_pmtnum = mrdpmti - MRDSpecs::layeroffsets.at(layernum);
		// paddles in each layer alternate on sides; i.e. paddles 0 and 1 are on opposite sides
		int side = in_layer_pmtnum%2;
		in_layer_pmtnum = std::floor(in_layer_pmtnum/2);
		int orientation = MRDSpecs::paddle_orientations.at(mrdpmti);
		int MRD_x = (orientation) ? in_layer_pmtnum  : side;
		int MRD_y = (orientation) ? side : in_layer_pmtnum;
		int MRD_z = layernum+2;  // first MRD z layer num is 2 (veto are 0,1)
		
		Paddle apaddle( uniquedetectorkey,
						MRD_x,
						MRD_y,
						MRD_z,
						orientation,
						Position( MRDSpecs::paddle_originx.at(mrdpmti)/1000.,
								  MRDSpecs::paddle_originy.at(mrdpmti)/1000.,
								  MRDSpecs::paddle_originz.at(mrdpmti)/1000.),
		std::pair<double,double>{MRDSpecs::paddle_extentsx.at(mrdpmti).first/1000.,
								 MRDSpecs::paddle_extentsx.at(mrdpmti).second/1000.},
		std::pair<double,double>{MRDSpecs::paddle_extentsy.at(mrdpmti).first/1000.,
								 MRDSpecs::paddle_extentsy.at(mrdpmti).second/1000.},
		std::pair<double,double>{MRDSpecs::paddle_extentsz.at(mrdpmti).first/1000.,
								 MRDSpecs::paddle_extentsz.at(mrdpmti).second/1000.});
		if(verbosity>4) cout<<"Setting paddle for detector "<<uniquedetectorkey<<endl;
		anniegeom->SetDetectorPaddle(uniquedetectorkey,apaddle);
		
		if(verbosity>4) cout<<"printing geometry"<<endl;
		if(verbosity>4) anniegeom->PrintChannels();
	}
	
	// veto PMTs
	for(int faccpmti=0; faccpmti<numvetopmts; faccpmti++){
		WCSimRootPMT apmt = wcsimrootgeom->GetMRDPMT(faccpmti);
		
		// Construct the detector associated with this PMT
		unsigned long uniquedetectorkey = anniegeom->ConsumeNextFreeDetectorKey();
		std::string CylLocString;
		switch (apmt.GetCylLoc()){
			case 0:  CylLocString = "TopCap";    break;
			case 2:  CylLocString = "BottomCap"; break;
			case 1:  CylLocString = "Barrel";    break;
			case 4:  CylLocString = "MRD";       break;  // TODO set this as H or V paddle? And layer?
			case 5:  CylLocString = "Veto";      break;  // TODO set layer?
			default: CylLocString = "NA";        break;  // unknown
		}
		Detector adet(uniquedetectorkey,
					  "Veto",
					  CylLocString,
					  Position( apmt.GetPosition(0)/100.,
					            apmt.GetPosition(1)/100.,
					            apmt.GetPosition(2)/100.),
					  Direction(apmt.GetOrientation(0),
					            apmt.GetOrientation(1),
					            apmt.GetOrientation(2)),
					  apmt.GetName(),
					  detectorstatus::ON,
					  0.);
		
		// construct the channel associated with this PMT
		unsigned long uniquechannelkey = anniegeom->ConsumeNextFreeChannelKey();
		facc_tubeid_to_channelkey.emplace(apmt.GetTubeNo(), uniquechannelkey);
		channelkey_to_faccpmtid.emplace(uniquechannelkey, apmt.GetTubeNo());
		
		// fill up TDC cards and channels monotonically, they're arbitrary for simulation
		TDC_Chan_Num++;
		if(TDC_Chan_Num>=TDC_CHANNELS_PER_CARD)  { TDC_Chan_Num=0; TDC_Card_Num++;  }
		if(TDC_Card_Num>=TDC_CARDS_PER_CRATE)    { TDC_Card_Num=0; TDC_Crate_Num++; }
		// same for HV
		LeCroy_HV_Chan_Num++;
		if(LeCroy_HV_Chan_Num>=LECROY_HV_CHANNELS_PER_CARD)    { LeCroy_HV_Chan_Num=0; LeCroy_HV_Card_Num++;  }
		if(LeCroy_HV_Card_Num>=LECROY_HV_CARDS_PER_CRATE) { LeCroy_HV_Card_Num=0; LeCroy_HV_Crate_Num++; }
		
		Channel pmtchannel( uniquechannelkey,
							Position(0,0,0.),
							0, // stripside
							0, // stripnum
							TDC_Crate_Num,
							TDC_Card_Num,
							TDC_Chan_Num,
							-1,                 // TDC has no level 2 signal handling
							-1,
							-1,
							LeCroy_HV_Crate_Num,
							LeCroy_HV_Card_Num,
							LeCroy_HV_Chan_Num,
							channelstatus::ON);
		
		// Add this channel to the geometry
		if(verbosity>4) cout<<"Adding channel "<<uniquechannelkey<<" to detector "<<uniquedetectorkey<<endl;
		adet.AddChannel(pmtchannel);
		
		// Add this detector to the geometry
		if(verbosity>4) cout<<"Adding detector "<<uniquedetectorkey<<" to geometry"<<endl;
		anniegeom->AddDetector(adet);
		
		// Create a paddle
		// FIXME even MRDSpecs can't save us here; instead hard-code the values for now,
		// this will all be replaced eventually anyway
		int MRD_z = (faccpmti>12); // 13 paddles per layer
		int MRD_x = MRD_z;         // i believe PMTs are on LHS for layer 0, RHS for layer 1
		int MRD_y = faccpmti - (13*MRD_z);
		double paddle_zorigin = (MRD_z) ? 0.0728 : 0.0508;  // numbers from geofile.txt
		double paddle_yorigin = facc_paddle_yorigins.at(faccpmti)/100.;
		
		Paddle apaddle( uniquedetectorkey,
						MRD_x,
						MRD_y,
						MRD_z,
						0,  // orientation 0=horizontal, 1=vertical
						Position(0,paddle_yorigin,paddle_zorigin),
		std::pair<double,double>{-1.6,1.6},  // numbers from WCSim source files / measurements
		std::pair<double,double>{paddle_yorigin-0.1525,paddle_yorigin+0.1525},
		std::pair<double,double>{paddle_zorigin-0.01,paddle_zorigin+0.01});
		
		if(verbosity>4) cout<<"Setting paddle for detector "<<uniquedetectorkey<<endl;
		anniegeom->SetDetectorPaddle(uniquedetectorkey,apaddle);
		
		if(verbosity>4) cout<<"printing geometry"<<endl;
		if(verbosity>4) anniegeom->PrintChannels();
	}
	
	// for other WCSim tools that may need the WCSim Tube IDs
	m_data->CStore.Set("lappd_tubeid_to_detectorkey",lappd_tubeid_to_detectorkey);
	m_data->CStore.Set("pmt_tubeid_to_channelkey",pmt_tubeid_to_channelkey);
	m_data->CStore.Set("mrd_tubeid_to_channelkey",mrd_tubeid_to_channelkey);
	m_data->CStore.Set("facc_tubeid_to_channelkey",facc_tubeid_to_channelkey);
	// inverse
	m_data->CStore.Set("detectorkey_to_lappdid",detectorkey_to_lappdid);
	m_data->CStore.Set("channelkey_to_pmtid",channelkey_to_pmtid);
	m_data->CStore.Set("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	m_data->CStore.Set("channelkey_to_faccpmtid",channelkey_to_faccpmtid);
	
	return anniegeom;
}

void LoadWCSim::MakeParticleToPmtMap(WCSimRootTrigger* thistrig, WCSimRootTrigger* firstTrig, std::map<int,std::map<unsigned long,double>>* ParticleId_to_TubeIds, std::map<int,double>* ParticleId_to_Charge, std::map<int,unsigned long> tubeid_to_channelkey){
	if(thistrig==nullptr) return;
	ParticleId_to_TubeIds->clear();
	ParticleId_to_Charge->clear();
	// scan through the parents IDs of the photons contributing to each digit
	// make note of which parent contributes to which digit, and which digits are associated with each parent
	Log("Making Particle to PMT Map",v_message,verbosity);
	// technically the charge will be a lower limit as this sums the charge from all digits
	// that a given particle contributed to, but not all this digit's charge may have been
	// from this particle.
	
	// To match digits to their parent particles we need the corresponding CherenkovHitTimes
	// both CherenkovHits and CherenkovHitTimes are stored in the first trigger
	//--------------------------------------------------------------------------------------
	
	// Loop over all digits 
	int numdigits = thistrig->GetCherenkovDigiHits()->GetEntries();
	for(int digiti=0; digiti<numdigits; digiti++){
		WCSimRootCherenkovDigiHit* digihit =
			(WCSimRootCherenkovDigiHit*)thistrig->GetCherenkovDigiHits()->At(digiti);
		int tubeid = digihit->GetTubeId();
		// loop over the photons in this digit
		std::vector<int> truephotonindices = digihit->GetPhotonIds();
		for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
			int thephotonsid = truephotonindices.at(truephoton);
			// get the index of the photon CherenkovHit object in the TClonesArray
			if(WCSimVersion<2){
				if(timeArrayOffsetMap.size()==0) BuildTimeArrayOffsetMap(firstTrig);
				thephotonsid+=timeArrayOffsetMap.at(tubeid);
			}
			// Get the CherenkovHitTime object that records the photon's Parent ID
			WCSimRootCherenkovHitTime *thehittimeobject = 
				(WCSimRootCherenkovHitTime*)(firstTrig->GetCherenkovHitTimes()->At(thephotonsid));
			if(thehittimeobject==nullptr) cerr<<"HITTIME IS NULL"<<endl;
			// get the parent ID from the CherenkovHitTime
			Int_t thephotonsparenttrackid = (thehittimeobject) ? thehittimeobject->GetParentID() : -1;
			// We'll want a map of particle ID to channel keys, so convert WCSim TubeId to channelkey
			int channelkey = tubeid_to_channelkey.at(tubeid);
			
			// Finally we can record this pairing of Tube to Particle
			if(ParticleId_to_TubeIds->count(thephotonsparenttrackid)==0){
				// we've not recorded any hits for this particle: make an empty map for it
				ParticleId_to_TubeIds->emplace(thephotonsparenttrackid,std::map<unsigned long,double>{});
			}
			if(ParticleId_to_TubeIds->at(thephotonsparenttrackid).count(channelkey)==0){
				// in the map for this particle record that this tube was hit
				ParticleId_to_TubeIds->at(thephotonsparenttrackid).emplace(channelkey,digihit->GetQ());
//				
//				double particletime = -1;
//				int particletrigger = -1;
//				if(trackid_to_mcparticleindex->count(thephotonsparenttrackid)){
//					int particleindex = trackid_to_mcparticleindex->at(thephotonsparenttrackid);
//					MCParticle theparticle = MCParticles->at(particleindex);
//					particletime = theparticle.GetStartTime();
//					particletrigger = theparticle.GetMCTriggerNum();
//				}
//				double hittime = digihit->GetT();
//				double triggertime = thistrig->GetHeader()->GetDate();
//				std::cout<<"Particle "<<thephotonsparenttrackid<<" created at "<<particletime
//						 <<" in trigger "<<particletrigger
//						 <<" produced hit on MRD PMT "<<channelkey<<" at time "<<hittime
//						 <<" relative to trigger time at "<<triggertime
//						 <<std::endl;
//				
			} else {
				// add another hit on this tube from this particle
				ParticleId_to_TubeIds->at(thephotonsparenttrackid).at(channelkey)+=digihit->GetQ();
			}
			if(ParticleId_to_Charge->count(thephotonsparenttrackid)==0){
				// first time seeing this particle
				ParticleId_to_Charge->emplace(thephotonsparenttrackid,digihit->GetQ());
			} else {
				ParticleId_to_Charge->at(thephotonsparenttrackid)+=digihit->GetQ();
			}
		}
	}  // end loop over digits
}

std::vector<int> LoadWCSim::GetHitParentIds(WCSimRootCherenkovDigiHit* digihit, WCSimRootTrigger* firstTrig){
	/* Get the ID of the MCParticle(s) that produced this digit */
	std::vector<int> parentids; // a hit could technically have more than one contrbuting particle
	
	// loop over the photons in this digit
	std::vector<int> truephotonindices = digihit->GetPhotonIds();
	for(int truephoton=0; truephoton<truephotonindices.size(); truephoton++){
		int thephotonsid = truephotonindices.at(truephoton);
		// get the indices of the digit's photon CherenkovHitTime objects
		if(WCSimVersion<2){
			if(timeArrayOffsetMap.size()==0) BuildTimeArrayOffsetMap(firstTrig);
			thephotonsid+=timeArrayOffsetMap.at(digihit->GetTubeId());
		}
		// get the CherenkovHitTime objects themselves, which contain the photon parent IDs
		WCSimRootCherenkovHitTime *thehittimeobject = 
			(WCSimRootCherenkovHitTime*)(firstTrig->GetCherenkovHitTimes()->At(thephotonsid));
		if(thehittimeobject==nullptr) cerr<<"HITTIME IS NULL"<<endl;
		else {
			int theparenttrackid = thehittimeobject->GetParentID();
			// check if this parent track was saved. Not all particles are saved.
			if(trackid_to_mcparticleindex->count(theparenttrackid)){
				parentids.push_back(trackid_to_mcparticleindex->at(theparenttrackid));
			} // else this photon may have come from e.g. an electron or gamma that wasn't recorded
		}
	}
	return parentids;
}

void LoadWCSim::BuildTimeArrayOffsetMap(WCSimRootTrigger* firstTrig){
	if(WCSimVersion<2){
		// The CherenkovHitTimes is a flattened array (over PMTs) of arrays (over photons)
		// For WCSimVersion<2, the PhotonIds available from a digit are the indices 
		// *within the subarray for that PMT*
		// we therefore need we need to offset these indices by the start of the pmt's subarray.
		// This offset may be found by scanning the CherenkovHits array (over PMTs),
		// finding the correct TubeID, and retrieving the 'GetTotalPe(0)' member for this entry.
		int ncherenkovhits=firstTrig->GetCherenkovHits()->GetEntries(); //atrigt->GetNcherenkovhits();
		//int nhittimes = firstTrig->GetCherenkovHitTimes()->GetEntries();
		
		for(int ihit = 0; ihit < ncherenkovhits; ihit++){
			// each WCSimRootCherenkovHit represents a hit PMT
			WCSimRootCherenkovHit* hitobject = 
				(WCSimRootCherenkovHit*)firstTrig->GetCherenkovHits()->At(ihit);
			if(hitobject==nullptr) cerr<<"HITOBJECT IS NULL!"<<endl;
			int tubeNumber = hitobject->GetTubeID();
			int timeArrayOffset = hitobject->GetTotalPe(0);
			timeArrayOffsetMap.emplace(tubeNumber,timeArrayOffset);
		}
	} else {
		Log("LoadWCSim Tool: BuildTimeArrayOffsetMap called with WCSimVersion<2: This is not needed!?",v_error,verbosity);
	}
}
