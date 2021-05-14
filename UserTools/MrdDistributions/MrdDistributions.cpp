/* vim:set noexpandtab tabstop=4 wrap */
#include "MrdDistributions.h"

#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TApplication.h"

#include <sys/types.h> // for stat() test to see if file or folder
#include <sys/stat.h>
//#include <unistd.h>

MrdDistributions::MrdDistributions():Tool(){}

// TODO combine the trees, record all info; muon energy. For all muons - were they primary/secondary,
// how many tank PMTs did they hit, how many tank PMT hits were there, was the muon tank reconstructed,
// how many MRD pmts did they hit, how many MRD hits were there, was the muon long track reconstructed,
// or short track reconstructed? were there reconstructed tracks but not matched to this particle?...etc

// misc function
ROOT::Math::XYZVector PositionToXYZVector(Position posin);
TVector3 PositionToTVector3(Position posin);
double GetClosestApproach(TVector3 start, TVector3 stop, TVector3 point, TVector3* closestapp=nullptr);

bool MrdDistributions::Initialise(std::string configfile, DataModel &data){
	
	cout<<"Initializing tool MrdDistributions"<<endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// get configuration variables for this tool
	//m_variables.Print();
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("printTracks",printTracks);  // from the BoostStore
	m_variables.Get("drawHistos",drawHistos);
	m_variables.Get("plotDirectory",plotDirectory);
	outfilename="MrdDistributions.root";
	m_variables.Get("RootFileOutput",outfilename);
	
	bool isdir=false, plotDirectoryExists=false;
	struct stat s;
	if(stat(plotDirectory.c_str(),&s)==0){
		plotDirectoryExists=true;
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory
			isdir=true;  //it's a directory
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file
			isdir=false; //it's a file
		} else {
			//assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		plotDirectoryExists=false;
		//assert(false&&"stat failed on input path! Is it valid?"); // error
		// errors could also be because this is a file pattern: e.g. wcsim_0.4*.root
		isdir=false;
	}
	
	if(drawHistos && (!plotDirectoryExists || !isdir)){
		Log("MrdDistributions Tool: output directory "+plotDirectory+" does not exist or is not a writable directory; please check and re-run.",v_error,verbosity);
		return false;
	}
	
	// make the ROOT file
	MakeRootFile();
	
	// scaler counters
	// ~~~~~~~~~~~~~~~
	numents=0;    // number of events (~triggers) analysed
	totnumtracks=0;
	numstopped=0;
	numpenetrated=0;
	numsideexit=0;
	numtankintercepts=0;
	numtankmisses=0;
	
	// histograms
	// ~~~~~~~~~~
	if(drawHistos){
		gStyle->SetOptStat(0); // FIXME gStyle mod is too global really....
		// create the ROOT application to show histograms
		Log("MrdDistributions Tool: constructing TApplication",v_debug,verbosity);
		int myargc=0;
		//char *myargv[] = {(const char*)"mrddist"};
		// get or make the TApplication
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			Log("MrdDistributions Tool: Making global TApplication",v_error,verbosity);
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			Log("MrdDistributions Tool: Retrieving global TApplication",v_error,verbosity);
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}
	
	// TODO fix ranges
	gROOT->cd();
	hnumsubevs = new TH1F("hnumsubevs","Number of MRD SubEvents",20,0,10);
	hnumtracks = new TH1F("hnumtracks","Number of MRD Tracks",10,0,10);
	hrun = new TH1F("hrun","Run Number Histogram",20,0,9);
	hevent = new TH1F("hevent","Event Number Histogram",100,0,2000);
	hmrdsubev = new TH1F("hsubevent","MRD SubEvent Number Histogram",20,0,19);
	htrigger = new TH1F("htrigg","Trigger Number Histogram",10,0,9);
	hhangle = new TH1F("hhangle","Track Angle in Top View",100,-TMath::Pi(),TMath::Pi());
	hhangleerr = new TH1F("hhangleerr","Error in Track Angle in Top View",100,0,TMath::Pi());
	hvangle = new TH1F("hvangle","Track Angle in Side View",100,-TMath::Pi(),TMath::Pi());
	hvangleerr = new TH1F("hvangleerr","Error in Track Angle in Side View",100,0,TMath::Pi());
	htotangle = new TH1F("htotangle","Track Angle from Beam Axis",100,0,TMath::Pi());
	htotangleerr = new TH1F("htotangleerr","Error in Track Angle from Beam Axis",100,0,TMath::Pi());
	henergyloss = new TH1F("henergyloss","Track Energy Loss in MRD",100,0,2200);
	henergylosserr = new TH1F("henergylosserr","Error in Track Energy Loss in MRD",100,0,2200);
	htracklength = new TH1F("htracklength","Total Track Length in MRD",100,0,220);
	htrackpen = new TH1F("htrackpen","Track Penetration in MRD",100,0,200);
	htrackpenvseloss = new TH2F("htrackpenvseloss","Track Penetration vs E Loss",100,0,220,100,0,2200);
	htracklenvseloss = new TH2F("htracklenvseloss","Track Length vs E Loss",100,0,220,100,0,2200);
	// htrackstart: this is the start of the reconstructed track. We know the particle got
	// into the MRD, so we also have a back-projected entry point which is likely to be similar,
	// and we have no true equivalent. If done right, it probably isn't much use, but for debug,
	// it might highlight if we have tracks which missed the front few layers and still got reconstructed.
	htrackstart = new TH3D("htrackstart","Reco MRD Track Start Vertices",100,-170,170,100,300,480,100,-230,220);
	htrackstop = new TH3D("htrackstop","Reco MRD Track Stop Vertices",100,-170,170,100,300,480,100,-230,220);
	hpep = new TH3D("hpep","Back Projected Tank Exit",100,-500,500,100,0,480,100,-330,320);
	hmpep = new TH3D("hmpep","Back Projected MRD Entry",100,-500,500,100,0,480,100,-330,320);
	
	// truth versions for comparisons
	hnumsubevstrue = new TH1F("hnumsubevstrue","Number of MRD SubEvents",20,0,10);
	hnumtrackstrue = new TH1F("hnumtrackstrue","Number of MRD Tracks",10,0,10);
	hhangletrue = new TH1F("hangle","Track Angle in Top View",100,-TMath::Pi(),TMath::Pi());
	hvangletrue = new TH1F("vangle","Track Angle in Side View",100,-TMath::Pi(),TMath::Pi());
	htotangletrue = new TH1F("htotangletrue","Track Angle from Beam Axis",100,0,TMath::Pi());
	henergylosstrue = new TH1F("henergylosstrue","Track Energy Loss in MRD",100,0,2200);
	htracklengthtrue = new TH1F("htracklengthtrue","Total Track Length in MRD",100,0,220);
	htrackpentrue = new TH1F("htrackpentrue","Track Penetration in MRD",100,0,200);
	htrackpenvselosstrue = new TH2F("htrackpenvselosstrue","Track Penetration vs E Loss",100,0,220,100,0,2200);
	htracklenvselosstrue = new TH2F("htracklenvselosstrue","Track Length vs E Loss",100,0,220,100,0,2200);
	htrackstoptrue = new TH3D("htrackstoptrue","MRD Track Stop Vertices",100,-170,170,100,300,480,100,-230,220);
	hpeptrue = new TH3D("hpeptrue","Back Projected Tank Exit",100,-500,500,100,0,480,100,-330,320);
	hmpeptrue = new TH3D("hmpeptrue","Back Projected MRD Entry",100,-500,500,100,0,480,100,-330,320);
	
	hnumsubevs->SetLineColor(kBlue);
	hnumsubevstrue->SetLineColor(kRed);
	hnumtracks->SetLineColor(kBlue);
	hnumtrackstrue->SetLineColor(kRed);
	hhangle->SetLineColor(kBlue);
	hhangletrue->SetLineColor(kRed);
	hvangle->SetLineColor(kBlue);
	hvangletrue->SetLineColor(kRed);
	htotangle->SetLineColor(kBlue);
	htotangletrue->SetLineColor(kRed);
	henergyloss->SetLineColor(kBlue);
	henergylosstrue->SetLineColor(kRed);
	htracklength->SetLineColor(kBlue);
	htracklengthtrue->SetLineColor(kRed);
	htrackpen->SetLineColor(kBlue);
	htrackpentrue->SetLineColor(kRed);
	
	htrackpenvseloss->SetMarkerColor(kBlue);
	htrackpenvselosstrue->SetMarkerColor(kRed);
	htracklenvseloss->SetMarkerColor(kBlue);
	htracklenvselosstrue->SetMarkerColor(kRed);
	htrackstart->SetMarkerColor(kBlue);
	htrackstop->SetMarkerColor(kBlue);
	htrackstoptrue->SetMarkerColor(kRed);
	hpep->SetMarkerColor(kBlue);
	hpeptrue->SetMarkerColor(kRed);
	hmpep->SetMarkerColor(kBlue);
	hmpeptrue->SetMarkerColor(kRed);
	
	htrackpenvseloss->SetMarkerStyle(20);
	htrackpenvselosstrue->SetMarkerStyle(20);
	htracklenvseloss->SetMarkerStyle(20);
	htracklenvselosstrue->SetMarkerStyle(20);
	htrackstart->SetMarkerStyle(20);
	htrackstop->SetMarkerStyle(20);
	htrackstoptrue->SetMarkerStyle(20);
	hpep->SetMarkerStyle(20);
	hpeptrue->SetMarkerStyle(20);
	hmpep->SetMarkerStyle(20);
	hmpeptrue->SetMarkerStyle(20);
	
	m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	
	return true;
}


bool MrdDistributions::Execute(){
	
	// Get the ANNIE event and extract general event information
	// Maybe this information could be added to plots for reference
	m_data->Stores.at("ANNIEEvent")->Get("MCFile",MCFile);
	m_data->Stores.at("ANNIEEvent")->Get("RunNumber",RunNumber);
	m_data->Stores.at("ANNIEEvent")->Get("SubrunNumber",SubrunNumber);
	m_data->Stores.at("ANNIEEvent")->Get("EventNumber",EventNumber);
	m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
	m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
	MCEventNumRoot = static_cast<ULong64_t>(MCEventNum);
	TDCData = nullptr;
	m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
	ParticleId_to_MrdTubeIds = nullptr;
	m_data->Stores.at("ANNIEEvent")->Get("ParticleId_to_MrdTubeIds",ParticleId_to_MrdTubeIds);
	
	// retrieving true tracks from the BoostStore
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles",MCParticles);
	if(not get_ok){
		Log("MrdDistributions Tool: No MCParticles in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	
	// retrieving tank reco vertex from the BoostStore
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	recoVtxOrigin = ROOT::Math::XYZVector(-999,-999,-999);
	theExtendedVertex=nullptr;
	get_ok = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex", theExtendedVertex);
	if(get_ok){
		if(theExtendedVertex!=nullptr){
			if(theExtendedVertex->GetFOM()>=0){
				recoVtxOrigin = PositionToXYZVector(theExtendedVertex->GetPosition());
			}
		}
	}	

	// retrieving reconstructed tracks from the BoostStore
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// first the counts of subevents and tracks - do not use tracks beyond these
	get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents",numsubevs);
	if(not get_ok){
		Log("MrdDistributions Tool: No NumMrdSubEvents in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);
	if(not get_ok){
		Log("MrdDistributions Tool: No NumMrdTracks in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	logmessage = "MrdDistributions Tool: Event "+to_string(EventNumber)+" had "
	              +to_string(numtracksinev)+" MRD tracks in "+to_string(numsubevs)+" subevents";
	Log(logmessage,v_debug,verbosity);
	
	// then the actual tracks
	get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", theMrdTracks);
	if(not get_ok){
		Log("MrdDistributions Tool: No MRDTracks member of the MRDTracks BoostStore!",v_error,verbosity);
		Log("MRDTracks store contents:",v_error,verbosity);
		m_data->Stores["MRDTracks"]->Print(false);
		return false;
	}
	// sanity check
	if(theMrdTracks->size()<numtracksinev){
		cerr<<"Too few entries in MRDTracks vector relative to NumMrdTracks!"<<endl;
		// more is fine as we don't shrink for efficiency
	}
	
//	// also get the TClonesArray store, because the objects have some useful functions
//	intptr_t subevptr = nullptr;
//	get_ok = m_data->CStore.Get("MrdSubEventTClonesArray",subevptr);
//	if(not get_ok){
//		Log("MrdDistributions Tool: Failed to get MrdSubEventTClonesArray from CStore! Can't access track objects!",v_warning,verbosity);
//	} else {
//		thesubevptr = reinterpret_cast<TClonesArray*>(subevptr);
//	}
	
	// count the number of mrd hits and hit mrd pmts in the event
	NumHitMrdPMTsInEvent=0;
	NumMrdHitsInEvent=0;
	if(TDCData!=nullptr){
		for(auto&& apmt : *TDCData){
			if(channelkey_to_mrdpmtid.count(apmt.first)){
				NumHitMrdPMTsInEvent++;
				NumMrdHitsInEvent += apmt.second.size();
			}
		}
	}
	if(NumMrdHitsInEvent<0){
		std::cerr<<"ERROR: NEGATIVE NUMBER OF HITS IN EVENT!!!"<<std::endl;
	}
	
	// Try to retrieve matching done by MrdEfficiency Tool
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// TODO These are of course only relevant if analysing MC
	Reco_to_True_Id_Map.clear(); // maps MRD Reco Track ID to MCParticle ID
	True_to_Reco_Id_Map.clear();
	get_ok = m_data->CStore.Get("Reco_to_True_Id_Map",Reco_to_True_Id_Map);
	get_ok &= m_data->CStore.Get("True_to_Reco_Id_Map",True_to_Reco_Id_Map);
	if(not get_ok){
		logmessage  = "MrdDistributions Tool: No Reco_to_True_Id_Map Store in CStore!\n";
		logmessage += "Without running MrdEfficiency tool first, there will be no matching ";
		logmessage += "between MC and true particles, so truth data will be empty!";
		Log(logmessage,v_warning,verbosity);
	}
	// get the map that converts MCParticle ID to index in MCParticles vector,
	// so that we can retrieve the corresponding MCParticle
	// pointers in stores are confusing, best to consider them read-only:
	// XXX do NOT clear the map before hand XXX
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex",trackid_to_mcparticleindex);
	if(not get_ok){
		logmessage  = "MrdDistributions Tool: No TrackId_to_MCParticleIndex in ANNIEEvent!\n";
		logmessage += "Will not be able to retrieve truth information for successfully reconstructed particles!";
		Log(logmessage,v_warning,verbosity);
	}
	
	/////////////////////////////////////
	// RECONSTRUCTED TRACK DISTRIBUTIONS
	/////////////////////////////////////
	
	// Fill counter histograms
	// ~~~~~~~~~~~~~~~~~~~~~~~
	hnumsubevs->Fill(numsubevs);
	hnumsubevstrue->Fill(0);             // TODO calculate the number of true subevents?
	hnumtracks->Fill(numtracksinev);     // num tracks in a given MRD subevent - (a short time window)
	totnumtracks+=numtracksinev;         // total number of tracks analysed in this ToolChain run
	
	// Loop over reconstructed tracks and record their properties
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	ClearBranchVectors();
	Log("MrdDistributions Tool: Looping over theMrdTracks",v_debug,verbosity);
	for(int tracki=0; tracki<numtracksinev; tracki++){
		BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
		if(verbosity>3) cout<<"track "<<tracki<<" at "<<thisTrackAsBoostStore<<endl;
		// Get the track details from the BoostStore
		thisTrackAsBoostStore->Get("MrdTrackID",MrdTrackID);                    // int
		thisTrackAsBoostStore->Get("MrdSubEventID",MrdSubEventID);              // int
		thisTrackAsBoostStore->Get("InterceptsTank",InterceptsTank);            // bool
		thisTrackAsBoostStore->Get("StartTime",StartTime);                      // [ns]
		thisTrackAsBoostStore->Get("StartVertex",StartVertex);                  // [m]
		thisTrackAsBoostStore->Get("StopVertex",StopVertex);                    // [m]
		thisTrackAsBoostStore->Get("TrackAngle",TrackAngle);                    // [rad]
		thisTrackAsBoostStore->Get("TrackAngleError",TrackAngleError);          // dx/dz  TODO
		thisTrackAsBoostStore->Get("HtrackOrigin",HtrackOrigin);                // x posn at global z=0, [cm]
		thisTrackAsBoostStore->Get("HtrackOriginError",HtrackOriginError);      // [cm]
		thisTrackAsBoostStore->Get("HtrackGradient",HtrackGradient);            // dx/dz
		thisTrackAsBoostStore->Get("HtrackGradientError",HtrackGradientError);  // [no units]
		thisTrackAsBoostStore->Get("VtrackOrigin",VtrackOrigin);                // [cm]
		thisTrackAsBoostStore->Get("VtrackOriginError",VtrackOriginError);      // [cm]
		thisTrackAsBoostStore->Get("VtrackGradient",VtrackGradient);            // dy/dz
		thisTrackAsBoostStore->Get("VtrackGradientError",VtrackGradientError);  // [no units]
		thisTrackAsBoostStore->Get("LayersHit",LayersHit);                      // vector<int>
		thisTrackAsBoostStore->Get("TrackLength",TrackLength);                  // [m]
		thisTrackAsBoostStore->Get("EnergyLoss",EnergyLoss);                    // [MeV]
		thisTrackAsBoostStore->Get("EnergyLossError",EnergyLossError);          // [MeV]
		thisTrackAsBoostStore->Get("IsMrdPenetrating",IsMrdPenetrating);        // bool
		thisTrackAsBoostStore->Get("IsMrdStopped",IsMrdStopped);                // bool
		thisTrackAsBoostStore->Get("IsMrdSideExit",IsMrdSideExit);              // bool
		thisTrackAsBoostStore->Get("PenetrationDepth",PenetrationDepth);        // [m]
		thisTrackAsBoostStore->Get("HtrackFitChi2",HtrackFitChi2);              // 
		thisTrackAsBoostStore->Get("HtrackFitCov",HtrackFitCov);                // 
		thisTrackAsBoostStore->Get("VtrackFitChi2",VtrackFitChi2);              // 
		thisTrackAsBoostStore->Get("VtrackFitCov",VtrackFitCov);                // 
		thisTrackAsBoostStore->Get("PMTsHit",PMTsHit);                          // 
		thisTrackAsBoostStore->Get("TankExitPoint",TankExitPoint);              // [m]
		thisTrackAsBoostStore->Get("MrdEntryPoint",MrdEntryPoint);              // [m]
		thisTrackAsBoostStore->Get("LongTrack",LongTrack);                      // int
		
		// some additional histograms are available in the MrdPaddlePlot Tool,
		// since the cMRDTrack classes used for reconstruction contain mrdcluster objects
		// which store information on position of hit paddles within the layer and hit times etc
		// that aren't pushed to the MRDTracks Store
		
		// Fill histograms of reconstructed track property distributions
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if(IsMrdPenetrating) numpenetrated++;
		else if(IsMrdStopped) numstopped++;
		else numsideexit++;
		if(InterceptsTank) numtankintercepts++;
		else numtankmisses++;
		
		hrun->Fill(RunNumber);
		hevent->Fill(EventNumber);
		hmrdsubev->Fill(MrdSubEventID);
		htrigger->Fill(MCTriggernum);
		
		hhangle->Fill(atan(HtrackGradient));
		hhangleerr->Fill(atan(HtrackGradientError));
		hvangle->Fill(atan(VtrackGradient));
		hvangleerr->Fill(atan(VtrackGradientError));
		htotangle->Fill(TrackAngle);
		htotangleerr->Fill(TrackAngleError);
		
		henergyloss->Fill(EnergyLoss);
		//std::cout<<"Reconstructed energy loss is "<<EnergyLoss<<std::endl;
		henergylosserr->Fill(EnergyLossError);
		htracklength->Fill(TrackLength*100.);
		htrackpen->Fill(PenetrationDepth*100.);
		htrackpenvseloss->Fill(PenetrationDepth*100.,EnergyLoss);
		htracklenvseloss->Fill(TrackLength*100.,EnergyLoss);
		
		htrackstart->Fill(StartVertex.X()*100.,StartVertex.Z()*100.,StartVertex.Y()*100.);
		htrackstop->Fill(StopVertex.X()*100.,StopVertex.Z()*100.,StopVertex.Y()*100.);
		hpep->Fill(TankExitPoint.X()*100., TankExitPoint.Z()*100.,TankExitPoint.Y()*100.);
		hmpep->Fill(MrdEntryPoint.X()*100., MrdEntryPoint.Z()*100., MrdEntryPoint.Y()*100.);
		//std::cout<<"Reco MrdEntryPoint("<<MrdEntryPoint.X()*100.<<", "<<MrdEntryPoint.Y()*100.
		//		 <<", "<<MrdEntryPoint.Z()*100.<<")"<<std::endl;
		
		// Print the reconstructed track info
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if(printTracks){
			cout<<"MrdTrack "<<MrdTrackID<<" in subevent "<<MrdSubEventID<<" was reconstructed from "
				<<PMTsHit.size()<<" MRD PMTs hit in "<<LayersHit.size()
				<<" layers for a total penetration depth of "<<PenetrationDepth<<" [m]."<<endl;
			cout<<"The track went from ";
			StartVertex.Print(false);
			cout<<" [m] to ";
			StopVertex.Print(false);
			cout<<" [m] starting at time "<<StartTime<<" [ns], travelling "<<TrackLength
				<<" [m] through the MRD at an angle of "<<TrackAngle
				<<" +/- "<<abs(TrackAngleError)<<" [rads] to the beam axis, before ";
			if(IsMrdStopped) cout<<"stopping in the MRD";
			if(IsMrdPenetrating) cout<<"penetrating the MRD";
			if(IsMrdSideExit) cout<<"leaving through the side of the MRD";
			cout<<endl<<"Back-projection suggests this track ";
			if(InterceptsTank) cout<<"intercepted"; else cout<<"did not intercept";
			cout<<" the tank."<<endl;
		}
		
		// append reco values into the output file vector
		fileout_MrdSubEventID.push_back(MrdSubEventID);
		fileout_HtrackAngle.push_back(atan(HtrackGradient));  // TODO check me, plots looks suspect
		fileout_HtrackAngleError.push_back(atan(HtrackGradientError));
		fileout_VtrackAngle.push_back(atan(VtrackGradient));
		fileout_VtrackAngleError.push_back(atan(VtrackGradientError));
		fileout_TrackAngle.push_back(TrackAngle);
		fileout_TrackAngleError.push_back(TrackAngleError);
		fileout_EnergyLoss.push_back(EnergyLoss);
		fileout_EnergyLossError.push_back(EnergyLossError);
		fileout_TrackLength.push_back(TrackLength*100.);
		fileout_PenetrationDepth.push_back(PenetrationDepth*100.);
		fileout_NumLayersHit.push_back(LayersHit.size());
		fileout_NumPmtsHit.push_back(PMTsHit.size());
		fileout_StartVertex.push_back(PositionToXYZVector(100.*StartVertex));
		fileout_StopVertex.push_back(PositionToXYZVector(100.*StopVertex));
		fileout_TankExitPoint.push_back(PositionToXYZVector(100.*TankExitPoint));
		fileout_MrdEntryPoint.push_back(PositionToXYZVector(100.*MrdEntryPoint));
		fileout_LongTrackReco.push_back(LongTrack);
		
		// See if we had a matching truth track (if MC) and record true details if so
		int mc_particles_index=-1;
		// check if MrdEfficiency tool found a matching true particle
		if(Reco_to_True_Id_Map.count(MrdTrackID)){
			// get its ID
			int true_track_id = Reco_to_True_Id_Map.at(MrdTrackID);
			// confirm LoadWCSim mapped this to an MCParticles index
			if(trackid_to_mcparticleindex->count(true_track_id)>0){
				// retrieve the index
				mc_particles_index = trackid_to_mcparticleindex->at(true_track_id);
			} else {
				// we could scan the MCParticles to find it, but something, somewhere went wrong...
				Log("MrdDistributions Tool: TrackId_to_MCParticleIndex map does not contain Particle ID "
					+ to_string(true_track_id) + " suggested by Reco_to_True_Id_Map",v_error,verbosity);
			}
		} else {
			Log("No MCParticle ID for this track; presumably it was a fake track from noise?",v_debug,verbosity);
		}
		// if we had a true particle, record it's details
		if(mc_particles_index>=0){
			MCParticle* nextparticle = &MCParticles->at(mc_particles_index);
			// Get the track details
			MCTruthParticleID = nextparticle->GetParticleID();
			InterceptsTank = nextparticle->GetEntersTank();
			StartTime = nextparticle->GetStartTime();
			TrackAngleX = nextparticle->GetTrackAngleX();
			TrackAngleY = nextparticle->GetTrackAngleY();
			TrackAngle = nextparticle->GetTrackAngleFromBeam();
			TrackLength = nextparticle->GetTrackLengthInMrd();        // the same as num layers with PMT hits!
			NumLayersHit = nextparticle->GetNumMrdLayersPenetrated(); // XXX ! num layers penetrated is not
			// ParticleId_to_MrdTubeIds is a std::map<ParticleId,std::map<ChannelKey,TotalCharge>>
			if(ParticleId_to_MrdTubeIds->count(MCTruthParticleID)>0){
				NumPmtsHit = ParticleId_to_MrdTubeIds->at(MCTruthParticleID).size();
			} else {
				NumPmtsHit = 0;
			}
			if((NumMrdHitsInEvent==0) && (NumPmtsHit>0)){
				std::cerr<<" ***** ERROR **** ParticleId_to_MrdTubeIds reports "<<NumPmtsHit
						 <<" tubes hit for MCParticle "<<MCTruthParticleID<<", but TDCData"
						 <<" contains no hits!!"<<std::endl;
				assert(false);
			}
			EnergyLoss = nextparticle->GetMrdEnergyLoss();
			StartVertex = nextparticle->GetMrdEntryPoint();
			StopVertex = nextparticle->GetMrdExitPoint();
			TankExitPoint = nextparticle->GetTankExitPoint();
			MrdEntryPoint = nextparticle->GetMrdEntryPoint();
			IsMrdStopped = ((nextparticle->GetEntersMrd())&&(!nextparticle->GetExitsMrd()));
			IsMrdPenetrating = nextparticle->GetPenetratesMrd();
			IsMrdSideExit = ((nextparticle->GetExitsMrd())&&(!IsMrdPenetrating));
			PenetrationDepth = nextparticle->GetMrdPenetration();
			TrueTrackOrigin = nextparticle->GetStartVertex();
			//RecoTrackOrigin = Position(0,0,0); // FIXME get tank reco vertex here
			
			// generate some back-projected point far enough back in the z-plane
			// - e.g 1m before the true track starts
			double zplane = (TrueTrackOrigin.Z() - 1.0)*100.;
			double yval = HtrackOrigin + HtrackGradient*zplane;
			double xval = VtrackOrigin + VtrackGradient*zplane;
			TVector3 pointa(xval, yval, zplane);
			// assume we can use the track endpoint and the corresponding line segment
			// will encompass the point of closest approach
			TVector3 pointb = PositionToTVector3(StopVertex)*100.;
			// find the point and distance of closest approach
			TVector3 trueorigincm = PositionToTVector3(TrueTrackOrigin)*100.;
			TVector3 closestapp(0,0,0);
			ClosestAppDist = GetClosestApproach(pointa, pointb, trueorigincm, &closestapp);
			ClosestAppPoint = Position(closestapp.X(),closestapp.Y(),closestapp.Z());
			
		} else {
			// no matching track, set truth variables to default
			MCTruthParticleID = -1;
			InterceptsTank = 0;
			StartTime = 0;
			TrackAngleX = 0;
			TrackAngleY = 0;
			TrackAngle = 0;
			NumLayersHit = 0;
			NumPmtsHit = 0;
			TrackLength = 0;
			EnergyLoss = 0;
			PenetrationDepth = 0;
			StartVertex = Position(0,0,0);
			StopVertex = Position(0,0,0);
			TankExitPoint = Position(0,0,0);
			MrdEntryPoint = Position(0,0,0);
			IsMrdStopped = 0;
			IsMrdPenetrating = 0;
			IsMrdSideExit = 0;
			TrueTrackOrigin = Position(0,0,0);
			//RecoTrackOrigin = Position(0,0,0);
			ClosestAppDist = 0;
			ClosestAppPoint = Position(0,0,0);
		}
		
		// append true values into the output file vector
		fileout_MCTruthParticleID.push_back(MCTruthParticleID);
		fileout_TrueTrackAngleX.push_back(TrackAngleX);
		fileout_TrueTrackAngleY.push_back(TrackAngleY);
		fileout_TrueTrackAngle.push_back(TrackAngle);
		fileout_TrueEnergyLoss.push_back(EnergyLoss);
		fileout_TrueTrackLength.push_back(TrackLength*100.);
		fileout_TruePenetrationDepth.push_back(PenetrationDepth*100.);
		fileout_TrueStopVertex.push_back(PositionToXYZVector(100.*StopVertex));
		fileout_TrueTankExitPoint.push_back(PositionToXYZVector(100.*TankExitPoint));
		fileout_TrueMrdEntryPoint.push_back(PositionToXYZVector(100.*MrdEntryPoint));
		fileout_TrueInterceptsTank.push_back(InterceptsTank);
		fileout_TrueStartTime.push_back(StartTime);
		fileout_TrueNumLayersHit.push_back(NumLayersHit);
		fileout_TrueNumPmtsHit.push_back(NumPmtsHit);
		fileout_TrueIsMrdStopped.push_back(IsMrdStopped);
		fileout_TrueIsMrdPenetrating.push_back(IsMrdPenetrating);
		fileout_TrueIsMrdSideExit.push_back(IsMrdSideExit);
		fileout_TrueOriginVertex.push_back(PositionToXYZVector(100.*TrueTrackOrigin));
		//fileout_RecoOriginVertex.push_back(PositionToXYZVector(100.*RecoTrackOrigin));
		fileout_TrueClosestApproachPoint.push_back(PositionToXYZVector(ClosestAppPoint));
		fileout_TrueClosestApproachDist.push_back(ClosestAppDist);
		
	}
	
	recotree->Fill();
	outfile->cd();
	recotree->Write("RecoTree",TObject::kOverwrite);
	gROOT->cd();
	ClearBranchVectors();
	
	/////////////////////////////////////
	// TRUE TRACK DISTRIBUTIONS
	/////////////////////////////////////
	
	// Fill counter histograms
	// ~~~~~~~~~~~~~~~~~~~~~~~
	nummrdtracksthisevent=0;
	
	// Loop over true tracks and record their properties
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Log("MrdDistributions Tool: Looping over MCParticles",v_debug,verbosity);
	for(int tracki=0; tracki<MCParticles->size(); tracki++){
		MCParticle* nextparticle = &MCParticles->at(tracki);
		if(fabs(nextparticle->GetPdgCode())!=13)                continue; // only interested in muons 	// EDIT: Added -13 (mu+) as possible option
		if(not (nextparticle->GetEntersMrd()))            continue; // only record tracks that entered the MRD
		if(nextparticle->GetMCTriggerNum()!=MCTriggernum) continue; // only record in their proper trigger
		nummrdtracksthisevent++;
		
		if(verbosity>3) cout<<"track "<<tracki<<endl;
		// Get the track details
		MCTruthParticleID = nextparticle->GetParticleID();
		if(True_to_Reco_Id_Map.count(MCTruthParticleID)){
			MrdTrackID = True_to_Reco_Id_Map.at(MCTruthParticleID);
		} else {
			MrdTrackID = -1;
		}
		InterceptsTank = nextparticle->GetEntersTank();                   // 
		StartTime = nextparticle->GetStartTime();                         // [ns]
		StartVertex = nextparticle->GetMrdEntryPoint();                   // [m] MRD entry vertex!
		StopVertex = nextparticle->GetMrdExitPoint();                     // [m] MRD stop/exit vertex!
		//cout<<"Truth Track goes = ("<<StartVertex.X()<<", "<<StartVertex.Y()<<", "
		//	<<StartVertex.Z()<<") -> ("<<StopVertex.X()<<", "<<StopVertex.Y()<<", "
		//	<<StopVertex.Z()<<")"<<endl;
		TrackAngleX = nextparticle->GetTrackAngleX();                     // [rad]
		TrackAngleY = nextparticle->GetTrackAngleY();                     // [rad]
		//std::cout<<"TrackAngleX="<<TrackAngleX<<", TrackAngleY="<<TrackAngleY<<std::endl;
		TrackAngle = nextparticle->GetTrackAngleFromBeam();               // [rad]
		NumLayersHit = nextparticle->GetNumMrdLayersPenetrated();         // 
		if(ParticleId_to_MrdTubeIds->count(MCTruthParticleID)>0){
			NumPmtsHit = ParticleId_to_MrdTubeIds->at(MCTruthParticleID).size();
		} else {
			NumPmtsHit = 0;
		}
		if((NumMrdHitsInEvent==0) && (NumPmtsHit>0)){
			std::cerr<<" ***** ERROR **** ParticleId_to_MrdTubeIds reports "<<NumPmtsHit
					 <<" tubes hit for MCParticle "<<MCTruthParticleID
					 <<" from trigger "<<nextparticle->GetMCTriggerNum()<<", but TDCData"
					 <<" for trigger num "<<MCTriggernum<<" contains no hits!!"<<std::endl;
			assert(false);
		}
		TrackLength = nextparticle->GetTrackLengthInMrd();                // [m]
		EnergyLoss = nextparticle->GetMrdEnergyLoss();                    // [MeV]
		IsMrdStopped = ((nextparticle->GetEntersMrd())&&(!nextparticle->GetExitsMrd()));    // bool
		IsMrdPenetrating = nextparticle->GetPenetratesMrd();              // bool: fully penetrating
		IsMrdSideExit = ((nextparticle->GetExitsMrd())&&(!IsMrdPenetrating));               // bool
		PenetrationDepth = nextparticle->GetMrdPenetration();             // [m]
		
		TankExitPoint = nextparticle->GetTankExitPoint();                 // [m]
		//cout<<"TankExitPoint: ("<<TankExitPoint.X()<<", "<<TankExitPoint.Y()<<", "<<TankExitPoint.Z()<<")"<<endl;
		MrdEntryPoint = nextparticle->GetMrdEntryPoint();                 // [m]
		
		// Fill histograms of reconstructed track property distributions
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if(IsMrdPenetrating) numpenetratedtrue++;
		else if(IsMrdStopped) numstoppedtrue++;
		else numsideexittrue++;
		if(InterceptsTank) numtankinterceptstrue++;
		else numtankmissestrue++;
		
		hvangletrue->Fill(TrackAngleX);
		hhangletrue->Fill(TrackAngleY);
		htotangletrue->Fill(TrackAngle);
		//std::cout<<"\"True\" Energy loss: "<<EnergyLoss<<std::endl;
		henergylosstrue->Fill(EnergyLoss); // Pull from WCSim
		//cout<<"true track length in MRD: "<<TrackLength*100.<<endl;
		htracklengthtrue->Fill(TrackLength*100.);
		htrackpentrue->Fill(PenetrationDepth*100.);
		htrackpenvselosstrue->Fill(PenetrationDepth*100.,EnergyLoss);
		htracklenvselosstrue->Fill(TrackLength*100.,EnergyLoss);
		
		htrackstoptrue->Fill(StopVertex.X()*100.,StopVertex.Z()*100.,StopVertex.Y()*100.);
		hpeptrue->Fill(TankExitPoint.X()*100., TankExitPoint.Z()*100.,TankExitPoint.Y()*100.);
		hmpeptrue->Fill(MrdEntryPoint.X()*100., MrdEntryPoint.Z()*100., MrdEntryPoint.Y()*100.);
		//std::cout<<"True MrdEntryPoint ("<<MrdEntryPoint.X()*100.<<", "<<MrdEntryPoint.Y()*100.
		//		 <<", "<<MrdEntryPoint.Z()*100.<<")"<<std::endl;
		
		// Print the true track info
		// ~~~~~~~~~~~~~~~~~~~~~~~~~
		if(printTracks){
			cout<<"True MRD-intercepting particle with reco ID "<<MrdTrackID<<" hit "
				<<PMTsHit.size()<<" MRD PMTs hit in "<<NumLayersHit
				<<" layers for a total penetration depth of "<<PenetrationDepth<<" [m]."<<endl;
			cout<<"The track went from ";
			StartVertex.Print(false);
			cout<<" [m] to ";
			StopVertex.Print(false);
			cout<<" [m] starting at time "<<StartTime<<" [ns], travelling "<<TrackLength
				<<" [m] through the MRD at an angle of "<<TrackAngle
				<<" +/- "<<abs(TrackAngleError)<<" [rads] to the beam axis, before ";
			if(IsMrdStopped) cout<<"stopping in the MRD";
			if(IsMrdPenetrating) cout<<"penetrating the MRD";
			if(IsMrdSideExit) cout<<"leaving through the side of the MRD";
			cout<<endl<<"Back-projection suggests this track ";
			if(InterceptsTank) cout<<"intercepted"; else cout<<"did not intercept";
			cout<<" the tank."<<endl;
		}
		
		// append values into the output file vector
		fileout_MCTruthParticleID.push_back(MCTruthParticleID);
		fileout_RecoParticleID.push_back(MrdTrackID);
		fileout_TrueTrackAngleX.push_back(TrackAngleX);
		fileout_TrueTrackAngleY.push_back(TrackAngleY);
		fileout_TrueTrackAngle.push_back(TrackAngle);
		fileout_TrueEnergyLoss.push_back(EnergyLoss);
		fileout_TrueTrackLength.push_back(TrackLength*100.);
		fileout_TruePenetrationDepth.push_back(PenetrationDepth*100.);
		fileout_TrueNumLayersHit.push_back(NumLayersHit);
		fileout_TrueNumPmtsHit.push_back(NumPmtsHit);
		fileout_TrueStopVertex.push_back(PositionToXYZVector(100.*StopVertex));
		fileout_TrueTankExitPoint.push_back(PositionToXYZVector(100.*TankExitPoint));
		fileout_TrueMrdEntryPoint.push_back(PositionToXYZVector(100.*MrdEntryPoint));
		fileout_TrueInterceptsTank.push_back(InterceptsTank);
		fileout_TrueStartTime.push_back(StartTime);
		fileout_TrueIsMrdStopped.push_back(IsMrdStopped);
		fileout_TrueIsMrdPenetrating.push_back(IsMrdPenetrating);
		fileout_TrueIsMrdSideExit.push_back(IsMrdSideExit);
		
	}
	//std::cout<<"MrdDistributions Tool: Filling event "<<EventNumber<<", MCEventNum "<<MCEventNum
	//		 <<", MCTriggerNum "<<MCTriggernum<<std::endl;
	
	truthtree->Fill();
	outfile->cd();
	truthtree->Write("TruthTree",TObject::kOverwrite);
	gROOT->cd();
	ClearBranchVectors();
	
	cout<<"num true particles intercepting MRD this event:" <<nummrdtracksthisevent<<endl;
	hnumtrackstrue->Fill(nummrdtracksthisevent);
	totnumtrackstrue+=nummrdtracksthisevent;
	
	numents++;
	
	return true;
}


bool MrdDistributions::Finalise(){
	
	gROOT->cd();
	cout<<"Analysed "<<numents<<" events, found "<<totnumtracks<<" MRD tracks, of which "
		<<numstopped<<" stopped in the MRD, "<<numpenetrated<<" fully penetrated and the remaining "
		<<numsideexit<<" exited the side."<<endl
		<<"Back projection suggests "<<numtankintercepts<<" tracks would have originated in the tank, while "
		<<numtankmisses<<" would not intercept the tank through back projection."<<endl;
	
	canvwidth = 700;
	canvheight = 600;
	mrdDistCanv = new TCanvas("mrdDistCanv","",canvwidth,canvheight);
	
	mrdDistCanv->cd();
	std::string imgname;

	if(outfile) outfile->cd();
	hnumsubevs->Draw();
	imgname=hnumsubevs->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hnumsubevs->Write();
	hnumtracks->Draw();
	hnumtrackstrue->Draw("same");
	imgname=hnumtracks->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hnumtracks->Write();
	hrun->Draw();
	imgname=hrun->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hrun->Write();
	hevent->Draw();
	imgname=hevent->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hevent->Write();
	hmrdsubev->Draw();
	imgname=hmrdsubev->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hmrdsubev->Write();
	htrigger->Draw();
	imgname=htrigger->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htrigger->Write();
	hhangle->Draw();
	hhangletrue->Draw("same");
	imgname=hhangle->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hhangle->Write();
	if(outfile) hhangletrue->Write();
//	hhangleerr->Draw();
//	imgname=hhangleerr->GetTitle();
//	std::replace(imgname.begin(), imgname.end(), ' ', '_');
//	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
//	if(outfile) hhangleerr->Write();
	hvangle->Draw();
	hvangletrue->Draw("same");
	imgname=hvangle->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hvangle->Write();
	if(outfile) hvangletrue->Write();
//	hvangleerr->Draw();
//	imgname=hvangleerr->GetTitle();
//	std::replace(imgname.begin(), imgname.end(), ' ', '_');
//	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	htotangle->Draw();
	htotangletrue->Draw("same");
	imgname=htotangle->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htotangle->Write();
	if(outfile) htotangletrue->Write();
//	htotangleerr->Draw();
//	imgname=htotangleerr->GetTitle();
//	std::replace(imgname.begin(), imgname.end(), ' ', '_');
//	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
//	if(outfile) htotangleerr->Write();
	henergyloss->Draw();
	henergylosstrue->Draw("same");
	imgname=henergyloss->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) henergyloss->Write();
	if(outfile) henergylosstrue->Write();
//	henergylosserr->Draw();
//	imgname=henergylosserr->GetTitle();
//	std::replace(imgname.begin(), imgname.end(), ' ', '_');
//	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
//	if(outfile) henergylosserr->Write();
	htracklength->Draw();
	imgname=htracklength->GetTitle();
	htracklengthtrue->Draw("same");
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htracklength->Write();
	if(outfile) htracklengthtrue->Write();
	htrackpen->Draw();
	htrackpentrue->Draw("same");
	imgname=htrackpen->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htrackpen->Write();
	if(outfile) htrackpentrue->Write();
	htrackpenvseloss->Draw();
	htrackpenvselosstrue->Draw("same");
	imgname=htrackpenvseloss->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htrackpenvseloss->Write();
	if(outfile) htrackpenvselosstrue->Write();
	htracklenvseloss->Draw();
	htracklenvselosstrue->Draw("same");
	imgname=htracklenvseloss->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htracklenvseloss->Write();
	if(outfile) htracklenvselosstrue->Write();
	htrackstart->Draw(); // has no true equivalent: this is where the reco track starts.
	// we plot where the particle actually enters the MRD separately, and only really
	// have truth information for the latter.
	imgname=htrackstart->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htrackstart->Write();
	htrackstop->Draw();
	htrackstoptrue->Draw("same");
	imgname=htrackstop->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) htrackstop->Write();
	if(outfile) htrackstoptrue->Write();
	hpep->Draw();
	hpeptrue->Draw("same");
	imgname=hpep->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hpep->Write();
	if(outfile) hpeptrue->Write();
	hmpep->Draw();
	hmpeptrue->Draw("same");
	imgname=hmpep->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdDistCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	if(outfile) hmpep->Write();
	if(outfile) hmpeptrue->Write();
	
//	if(drawHistos){ gPad->WaitPrimitive(); }
	
	// cleanup
	// ~~~~~~~
	std::vector<TH1*> histos {hnumsubevs, hnumtracks, hrun, hevent, hmrdsubev, htrigger, hhangle, hhangleerr, hvangle, hvangleerr, htotangle, htotangleerr, henergyloss, henergylosserr, htracklength, htrackpen, htrackpenvseloss, htracklenvseloss, htrackstart, htrackstop, hpep, hmpep};
	Log("MrdDistributions Tool: deleting reco histograms",v_debug,verbosity);
	for(TH1* ahisto : histos){ if(ahisto) delete ahisto; ahisto=nullptr; }
	
	histos = std::vector<TH1*>{hnumsubevstrue, hnumtrackstrue, hhangletrue, hvangletrue, htotangletrue, henergylosstrue, htracklengthtrue, htrackpentrue, htrackpenvselosstrue, htracklenvselosstrue, htrackstoptrue, hpeptrue, hmpeptrue};
	Log("MrdDistributions Tool: deleting truth histograms",v_debug,verbosity);
	for(TH1* ahisto : histos){ if(ahisto) delete ahisto; ahisto=nullptr; }
	
	Log("MrdDistributions Tool: deleting canvas",v_debug,verbosity);
	if(gROOT->FindObject("mrdDistCanv")){ 
		delete mrdDistCanv; mrdDistCanv=nullptr;
	}
	
	if(drawHistos){
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				Log("MrdDistributions Tool: deleting gloabl TApplication",v_debug,verbosity);
				delete rootTApp;
			}
		} else if (tapplicationusers>1){
			tapplicationusers--;
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		}
	}
	
	return true;
}

TFile* MrdDistributions::MakeRootFile(){
	// set the pointers to the vectors first
	// reco track info
	pfileout_MrdSubEventID = &fileout_MrdSubEventID;
	pfileout_HtrackAngle = &fileout_HtrackAngle;
	pfileout_HtrackAngleError = &fileout_HtrackAngleError;
	pfileout_VtrackAngle = &fileout_VtrackAngle;
	pfileout_VtrackAngleError = &fileout_VtrackAngleError;
	pfileout_TrackAngle = &fileout_TrackAngle;
	pfileout_TrackAngleError = &fileout_TrackAngleError;
	pfileout_EnergyLoss = &fileout_EnergyLoss;
	pfileout_EnergyLossError = &fileout_EnergyLossError;
	pfileout_TrackLength = &fileout_TrackLength;
	pfileout_PenetrationDepth = &fileout_PenetrationDepth;
	pfileout_NumLayersHit = &fileout_NumLayersHit;
	pfileout_NumPmtsHit = &fileout_NumPmtsHit;
	pfileout_StartVertex = &fileout_StartVertex;
	pfileout_StopVertex = &fileout_StopVertex;
	pfileout_TankExitPoint = &fileout_TankExitPoint;
	pfileout_MrdEntryPoint = &fileout_MrdEntryPoint;
	pfileout_LongTrackReco = &fileout_LongTrackReco;
	
	// equivalent true track info
	pfileout_MCTruthParticleID = &fileout_MCTruthParticleID;
	pfileout_TrueTrackAngleX = &fileout_TrueTrackAngleX;
	pfileout_TrueTrackAngleY = &fileout_TrueTrackAngleY;
	pfileout_TrueTrackAngle = &fileout_TrueTrackAngle;
	pfileout_TrueEnergyLoss = &fileout_TrueEnergyLoss;
	pfileout_TrueTrackLength = &fileout_TrueTrackLength;
	pfileout_TruePenetrationDepth = &fileout_TruePenetrationDepth;
	pfileout_TrueNumLayersHit = &fileout_TrueNumLayersHit;
	pfileout_TrueNumPmtsHit = &fileout_TrueNumPmtsHit;
	// (n.b. no true start vertex in mrd)
	pfileout_TrueStopVertex = &fileout_TrueStopVertex;
	pfileout_TrueTankExitPoint = &fileout_TrueTankExitPoint;
	pfileout_TrueMrdEntryPoint = &fileout_TrueMrdEntryPoint;
	
	/// additional truth-only track info
	pfileout_TrueStartTime = &fileout_TrueStartTime;
	pfileout_TrueInterceptsTank = &fileout_TrueInterceptsTank;
	pfileout_TrueIsMrdStopped = &fileout_TrueIsMrdStopped;
	pfileout_TrueIsMrdPenetrating = &fileout_TrueIsMrdPenetrating;
	pfileout_TrueIsMrdSideExit = &fileout_TrueIsMrdSideExit;
	pfileout_TrueOriginVertex = &fileout_TrueOriginVertex;
	//pfileout_RecoOriginVertex = &fileout_RecoOriginVertex; // reconstructed tank vertex... TODO add
	pfileout_TrueClosestApproachPoint = &fileout_TrueClosestApproachPoint;
	pfileout_TrueClosestApproachDist = &fileout_TrueClosestApproachDist;
	pfileout_RecoParticleID = &fileout_RecoParticleID;
	
	outfile = new TFile((plotDirectory+"/"+outfilename).c_str(),"RECREATE","MRD Track Distributions");
	outfile->cd();
	recotree = new TTree("RecoTree","Summary of Reconstructed Tracks, with their Truth Values if Matched");
	// event-wise information
	recotree->Branch("MCFile",&MCFile);
	recotree->Branch("RunNumber",&RunNumber);
	recotree->Branch("SubrunNumber",&SubrunNumber);
	recotree->Branch("EventNumber",&EventNumber);
	recotree->Branch("MCEventNum",&MCEventNumRoot);
	recotree->Branch("MCTriggernum",&MCTriggernum);
	recotree->Branch("NumSubEvs",&numsubevs);
	recotree->Branch("NumTracksInSubEv",&numtracksinev);
//	recotree->Branch("NumTrueTracksInSubEv",&nummrdtracksthisevent);	//Makes more sense in truth tree
	recotree->Branch("NumHitMrdPMTsInEvent",&NumHitMrdPMTsInEvent);
	recotree->Branch("NumMrdHitsInEvent",&NumMrdHitsInEvent);
//	recotree->Branch("NumHitMrdPMTsInEvent",&NumHitMrdPMTsInSubEvent);
//	recotree->Branch("NumHitMrdHitsInEvent",&NumMrdHitsInSubEvent);
	recotree->Branch("RecoVtxOrigin",&recoVtxOrigin);
	
	// track-wise information vectors
	recotree->Branch("MrdSubEventID", &pfileout_MrdSubEventID);
	recotree->Branch("TrackAngleX", &pfileout_HtrackAngle);
	recotree->Branch("TrackAngleXError", &pfileout_HtrackAngleError);
	recotree->Branch("TrackAngleY", &pfileout_VtrackAngle);
	recotree->Branch("TrackAngleYError", &pfileout_VtrackAngleError);
	recotree->Branch("TrackAngle", &pfileout_TrackAngle);
	recotree->Branch("TrackAngleError", &pfileout_TrackAngleError);
	recotree->Branch("EnergyLoss", &pfileout_EnergyLoss);
	recotree->Branch("EnergyLossError", &pfileout_EnergyLossError);
	recotree->Branch("TrackLength", &pfileout_TrackLength);
	recotree->Branch("PenetrationDepth", &pfileout_PenetrationDepth);
	recotree->Branch("NumLayersHit",&pfileout_NumLayersHit);
	recotree->Branch("NumPmtsHit",&pfileout_NumPmtsHit);
	recotree->Branch("StartVertex", &pfileout_StartVertex);
	recotree->Branch("StopVertex", &pfileout_StopVertex);
	recotree->Branch("TankExitPoint", &pfileout_TankExitPoint);
	recotree->Branch("MrdEntryPoint", &pfileout_MrdEntryPoint);
	recotree->Branch("LongTrackReco", &pfileout_LongTrackReco);
	
	// information about matched true track
	recotree->Branch("MCTruthParticleID", &pfileout_MCTruthParticleID);
	recotree->Branch("TrueTrackAngleX", &pfileout_TrueTrackAngleX);
	recotree->Branch("TrueTrackAngleY", &pfileout_TrueTrackAngleY);
	recotree->Branch("TrueTrackAngle", &pfileout_TrueTrackAngle);
	recotree->Branch("TrueEnergyLoss", &pfileout_TrueEnergyLoss);
	recotree->Branch("TrueTrackLength", &pfileout_TrueTrackLength);
	recotree->Branch("TruePenetrationDepth", &pfileout_TruePenetrationDepth);
	recotree->Branch("TrueNumLayersHit", &pfileout_TrueNumLayersHit);
	recotree->Branch("TrueNumPmtsHit",&pfileout_TrueNumPmtsHit);
	// no true (MRD reconstruction) start vertex by definition
	recotree->Branch("TrueStopVertex", &pfileout_TrueStopVertex);
	recotree->Branch("TrueTankExitPoint", &pfileout_TrueTankExitPoint);
	recotree->Branch("TrueMrdEntryPoint", &pfileout_TrueMrdEntryPoint);
	
	// truth only info
	recotree->Branch("TrueInterceptsTank", &pfileout_TrueInterceptsTank);
	recotree->Branch("TrueStartTime", &pfileout_TrueStartTime);
	recotree->Branch("TrueIsMrdStopped", &pfileout_TrueIsMrdStopped);
	recotree->Branch("TrueIsMrdPenetrating", &pfileout_TrueIsMrdPenetrating);
	recotree->Branch("TrueIsMrdSideExit", &pfileout_TrueIsMrdSideExit);
	recotree->Branch("TrueTrueOriginVertex", &pfileout_TrueOriginVertex);
	//recotree->Branch("RecoOriginVertex", &pfileout_RecoOriginVertex);
	recotree->Branch("TrueClosestApproachPoint", &pfileout_TrueClosestApproachPoint);
	recotree->Branch("TrueClosestApproachDist", &pfileout_TrueClosestApproachDist);
	
	// we'll have a second tree that will record all the true tracks
	// this ensures we have an un-cut version of the distributions of tracks that went into the MRD,
	// without any selection on whether there was a reco track associated with it
	truthtree = new TTree("TruthTree","Summary of all True Tracks in the Mrd");
	// event-wise information
	truthtree->Branch("MCFile",&MCFile);
	truthtree->Branch("RunNumber",&RunNumber);
	truthtree->Branch("SubrunNumber",&SubrunNumber);
	truthtree->Branch("EventNumber",&EventNumber);
	truthtree->Branch("MCEventNum",&MCEventNumRoot);
	truthtree->Branch("MCTriggernum",&MCTriggernum);
	truthtree->Branch("NumSubEvs",&numsubevs);
	truthtree->Branch("NumHitMrdPMTsInEvent",&NumHitMrdPMTsInEvent);
	truthtree->Branch("NumMrdHitsInEvent",&NumMrdHitsInEvent);
//	truthtree->Branch("NumTracksInSubEv",&numtracksinev);
	truthtree->Branch("NumTrueTracksInSubEv",&nummrdtracksthisevent);	//Makes more sense here compared to recotree
	truthtree->Branch("RecoVtxOrigin",&recoVtxOrigin);
	// track-wise information
	truthtree->Branch("MCTruthParticleID", &pfileout_MCTruthParticleID);
	truthtree->Branch("RecoParticleID",&pfileout_RecoParticleID);
	truthtree->Branch("TrueTrackAngleX", &pfileout_TrueTrackAngleX);
	truthtree->Branch("TrueTrackAngleY", &pfileout_TrueTrackAngleY);
	truthtree->Branch("TrueTrackAngle", &pfileout_TrueTrackAngle);
	truthtree->Branch("TrueEnergyLoss", &pfileout_TrueEnergyLoss);
	truthtree->Branch("TrueTrackLength", &pfileout_TrueTrackLength);
	truthtree->Branch("TruePenetrationDepth", &pfileout_TruePenetrationDepth);
	truthtree->Branch("TrueNumLayersHit", &pfileout_TrueNumLayersHit);
	truthtree->Branch("TrueNumPmtsHit", &pfileout_TrueNumPmtsHit);
	truthtree->Branch("TrueStopVertex", &pfileout_TrueStopVertex);
	truthtree->Branch("TrueTankExitPoint", &pfileout_TrueTankExitPoint);
	truthtree->Branch("TrueMrdEntryPoint", &pfileout_TrueMrdEntryPoint);
	truthtree->Branch("TrueInterceptsTank", &pfileout_TrueInterceptsTank);
	truthtree->Branch("TrueStartTime", &pfileout_TrueStartTime);
	truthtree->Branch("TrueIsMrdStopped", &pfileout_TrueIsMrdStopped);
	truthtree->Branch("TrueIsMrdPenetrating", &pfileout_TrueIsMrdPenetrating);
	truthtree->Branch("TrueIsMrdSideExit", &pfileout_TrueIsMrdSideExit);
	
	gROOT->cd();
	return outfile;
}

void MrdDistributions::ClearBranchVectors(){
	fileout_MrdSubEventID.clear();
	fileout_HtrackAngle.clear();
	fileout_HtrackAngleError.clear();
	fileout_VtrackAngle.clear();
	fileout_VtrackAngleError.clear();
	fileout_TrackAngle.clear();
	fileout_TrackAngleError.clear();
	fileout_EnergyLoss.clear();
	fileout_EnergyLossError.clear();
	fileout_TrackLength.clear();
	fileout_PenetrationDepth.clear();
	fileout_NumLayersHit.clear();
	fileout_NumPmtsHit.clear();
	fileout_StartVertex.clear();
	fileout_TankExitPoint.clear();
	fileout_MrdEntryPoint.clear();
	fileout_LongTrackReco.clear();
	
	fileout_MCTruthParticleID.clear();
	fileout_TrueTrackAngleX.clear();
	fileout_TrueTrackAngleY.clear();
	fileout_TrueTrackAngle.clear();
	fileout_TrueEnergyLoss.clear();
	fileout_TrueTrackLength.clear();
	fileout_TruePenetrationDepth.clear();
	fileout_TrueStopVertex.clear();
	fileout_TrueTankExitPoint.clear();
	fileout_TrueMrdEntryPoint.clear();
	fileout_TrueStartTime.clear();
	fileout_TrueInterceptsTank.clear();
	fileout_TrueNumLayersHit.clear();
	fileout_TrueNumPmtsHit.clear();
	fileout_TrueIsMrdStopped.clear();
	fileout_TrueIsMrdPenetrating.clear();
	fileout_TrueIsMrdSideExit.clear();
	fileout_TrueOriginVertex.clear();
	//fileout_RecoOriginVertex.clear();
	fileout_TrueClosestApproachPoint.clear();
	fileout_TrueClosestApproachDist.clear();
	fileout_RecoParticleID.clear();
}

ROOT::Math::XYZVector PositionToXYZVector(Position posin){
	return ROOT::Math::XYZVector(posin.X(),posin.Y(),posin.Z());
}

TVector3 PositionToTVector3(Position posin){
	return TVector3(posin.X(),posin.Y(),posin.Z());
}

double GetClosestApproach(TVector3 start, TVector3 stop, TVector3 point, TVector3* closestapp){
	/*
	READ CODE FIRST
	(reinterpret based on using dot product to determine if moving away from point,
	and if directional distance to point is greater than distance to end)
	
	The closest distance to the point may either be the segment perpendicular to the point,
	if the line segment passes through it, or the distance from a segment end to the point,
	if it does not.
	We can determine which we want using the dot product of the line segment and the vector
	from the line end to the point:
	If the dot product is <0 the angle between the vectors is obtuse and the line
	segment does not pass through the perpendicular. Otherwise, it does.
	If we call the vectors from the segment endpoints to the point P, v0 and v1,
	and the vector of the line segment s, then the two values we need are (v0.s) and (v1.s).
	Writing v1 = (v0+s), we have (v1.s) = (v0.s) + (s.s), so we can determine which
	distance to use after calculating just (v0.s) and (s.s).
	*/
	
	TVector3 segment = stop - start;
	TVector3 start_to_point = point - start;
	
	double c1 = start_to_point.Dot(segment);  // angle between (directional) segment and from start to P
	if(c1<=0) return start_to_point.Mag();    // track moves away from point, closest approach is start
	
	double c2 = segment.Mag2();
	if(c2<=c1) return (point-stop).Mag();     // distance from start to point > distance from start to
	                                          // end of segment - closest approach is end
	
	double b = c1/c2;                         // neither above condition is met: the point must lie
	TVector3 Pb = start + b*segment;          // between endpoints, along the segment.
	if(closestapp!=nullptr) (*closestapp)=Pb;
	return (point-Pb).Mag();                  // Project from start to perpendicular (parametric base)
}
