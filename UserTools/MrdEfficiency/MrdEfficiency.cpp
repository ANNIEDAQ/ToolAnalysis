/* vim:set noexpandtab tabstop=4 wrap */
#include "MrdEfficiency.h"
#include "TH1.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TFile.h"
#include <thread>
#include <chrono>
#include <sys/types.h> // for stat() test to see if file or folder
#include <sys/stat.h>
//#include <unistd.h>

MrdEfficiency::MrdEfficiency():Tool(){}

bool MrdEfficiency::Initialise(std::string configfile, DataModel &data){
	
	if(verbosity) cout<<"Initializing tool MrdEfficiency"<<endl;
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile);  // loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////
	
	// get configuration variables for this tool
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("plotDirectory",plotDirectory);
	m_variables.Get("drawHistos",drawHistos);
	std::string outputfile="";
	get_ok = m_variables.Get("outputfile",outputfile);
	if(not get_ok){ outputfile="MrdEfficiency.root"; }
	
	// check the output directory exists and is suitable
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
	
	if((!plotDirectoryExists) || (!isdir)){
		Log("MrdEfficiency Tool: output directory "+plotDirectory+" does not exist or is not a writable directory; please check and re-run.",v_error,verbosity);
		return false;
	}
	
	// histograms
	// ~~~~~~~~~~
	if(drawHistos){ // if showing plots as well as saving, we need a TApplication
		// create the ROOT application to show histograms
		int myargc=0;
		//char *myargv[] = {(const char*)"mrdeff"};
		// get or make the TApplication
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			Log("MrdEfficiency Tool: Making global TApplication",v_error,verbosity);
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			Log("MrdEfficiency Tool: Retrieving global TApplication",v_error,verbosity);
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}

	std::string outputfilepath = plotDirectory+"/"+outputfile;
	fileout = new TFile(outputfilepath.c_str(),"RECREATE");
	fileout->cd();
	canvwidth = 700;
	canvheight = 600;
	
	hnumcorrectlymatched = new TH1F("hnumtruematched","Number of True Tracks Matched",20,0,10);
	hnumtruenotmatched = new TH1F("hnumtruenotmatched","Number of True Tracks Not Matched",20,0,10);
	hnumreconotmatched = new TH1F("hnumreconotmatched","Number of Reconstructed Tracks Not Matched",20,0,10);
	
	// distributions of properties for primary muons that were successfully reconstructed
	hhangle_recod = new TH1F("hhangle_recod","Track Angle in Top View, Reconstructed",20,-TMath::Pi()/2.,TMath::Pi()/2.);
	hvangle_recod = new TH1F("hvangle_recod","Track Angle in Side View, Reconstructed",20,-TMath::Pi()/2.,TMath::Pi()/2.);
	htotangle_recod = new TH1F("htotangle_recod","Track Angle from Beam Axis, Reconstructed",20,0,TMath::Pi()/2.);
	henergyloss_recod = new TH1F("henergyloss_recod","Track Energy Loss in MRD, Reconstructed",20,0,2200);
	htracklength_recod = new TH1F("htracklength_recod","Total Track Length in MRD, Reconstructed",20,0,220);
	htrackpen_recod = new TH1F("htrackpen_recod","Track Penetration in MRD, Reconstructed",20,0,200);
	hnummrdpmts_recod = new TH1F("hnummrdpmts_recod","Number of True MRD PMTs Hit by Particle, Reconstucted",20,0,20);
	hq2_recod = new TH1F("hq2_recod","True Q2 of Events, Reconstructed",20,0,2000);
	htrackstart_recod = new TH3D("htrackstart_recod","MRD Track Start Vertices, Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	htrackstop_recod = new TH3D("htrackstop_recod","MRD Track Stop Vertices, Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	hpep_recod = new TH3D("hpep_recod","Back Projected Tank Exit, Reconstructed", 100,-500,500,100,0,480,100,-330,320);
	hmpep_recod = new TH3D("hmpep_recod","Back Projected MRD Entry, Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	
	// distributions of properties for primary muons that were not reconstructed
	hhangle_nrecod = new TH1F("hhangle_nrecod","Track Angle in Top View, Not Reconstructed",20,-TMath::Pi()/2.,TMath::Pi()/2.);
	hvangle_nrecod = new TH1F("hvangle_nrecod","Track Angle in Side View, Not Reconstructed",20,-TMath::Pi()/2.,TMath::Pi()/2.);
	htotangle_nrecod = new TH1F("htotangle_nrecod","Track Angle from Beam Axis, Not Reconstructed", 20,0,TMath::Pi()/2.);
	henergyloss_nrecod = new TH1F("henergyloss_nrecod","Track Energy Loss in MRD, Not Reconstructed", 20,0,2200);
	htracklength_nrecod = new TH1F("htracklength_nrecod","Total Track Length in MRD, Not Reconstructed", 20,0,220);
	htrackpen_nrecod = new TH1F("htrackpen_nrecod","Track Penetration in MRD, Not Reconstructed",20,0,200);
	hnummrdpmts_nrecod = new TH1F("hnummrdpmts_nrecod","Number of True MRD PMTs Hit by Particle, Not Reconstucted",20,0,20);
	hq2_nrecod = new TH1F("hq2_nrecod","True Q2 of Events, Not Reconstructed",20,0,2000);
	htrackstart_nrecod = new TH3D("htrackstart_nrecod","MRD Track Start Vertices, Not Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	htrackstop_nrecod = new TH3D("htrackstop_nrecod","MRD Track Stop Vertices, Not Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	hpep_nrecod = new TH3D("hpep_nrecod","Back Projected Tank Exit, Not Reconstructed", 100,-500,500,100,0,480,100,-330,320);
	hmpep_nrecod = new TH3D("hmpep_nrecod","Back Projected MRD Entry, Not Reconstructed", 100,-170,170,100,300,480,100,-230,220);
	
	hhangle_recod->SetLineColor(kRed);
	hvangle_recod->SetLineColor(kRed);
	htotangle_recod->SetLineColor(kRed);
	henergyloss_recod->SetLineColor(kRed);
	htracklength_recod->SetLineColor(kRed);
	htrackpen_recod->SetLineColor(kRed);
	hnummrdpmts_recod->SetLineColor(kRed);
	hq2_recod->SetLineColor(kRed);
	htrackstart_recod->SetMarkerColor(kRed);
	htrackstart_recod->SetMarkerStyle(20);
	htrackstop_recod->SetMarkerColor(kRed);
	htrackstop_recod->SetMarkerStyle(20);
	hpep_recod->SetMarkerColor(kRed);
	hpep_recod->SetMarkerStyle(20);
	hmpep_recod->SetMarkerColor(kRed);
	hmpep_recod->SetMarkerStyle(20);
	
	hhangle_nrecod->SetLineColor(kBlue);
	hvangle_nrecod->SetLineColor(kBlue);
	htotangle_nrecod->SetLineColor(kBlue);
	henergyloss_nrecod->SetLineColor(kBlue);
	htracklength_nrecod->SetLineColor(kBlue);
	htrackpen_nrecod->SetLineColor(kBlue);
	hnummrdpmts_nrecod->SetLineColor(kBlue);
	hq2_nrecod->SetLineColor(kBlue);
	htrackstart_nrecod->SetMarkerColor(kBlue);
	htrackstart_nrecod->SetMarkerStyle(20);
	htrackstop_nrecod->SetMarkerColor(kBlue);
	htrackstop_nrecod->SetMarkerStyle(20);
	hpep_nrecod->SetMarkerColor(kBlue);
	hpep_nrecod->SetMarkerStyle(20);
	hmpep_nrecod->SetMarkerColor(kBlue);
	hmpep_nrecod->SetMarkerStyle(20);
	
	gROOT->cd();
	
	// to match true and reconstructed tracks we match the ids of paddles they contain
	// But ParticleId_to_MrdTubeIds matches to paddle ChannelKeys, whereas reconstructed
	// tracks contain a list of WCSim TubeIds. We need to map from one to the other
	get_ok = m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	if(not get_ok){
		Log("MrdEfficiency Tool: Could not find channelkey_to_mrdpmtid map in the CStore!",v_error,verbosity);
		return false;
	}
	
	// XXX Remove me, make less global
	gStyle->SetOptStat(0);
	
	// for later tools, in case we don't find any, put an empty pair of maps in now
	m_data->CStore.Set("Reco_to_True_Id_Map",Reco_to_True_Id_Map);
	m_data->CStore.Set("True_to_Reco_Id_Map",True_to_Reco_Id_Map);
	
	return true;
}


bool MrdEfficiency::Execute(){
	
	if(verbosity) cout<<"Executing tool MrdEfficiency"<<endl;
	// To measure efficiency we need to try to match true tracks with reconstructed tracks
	// The eaiest way to do this is by comparing the hit paddles that make up the track.
	
	// before anything else let's clear the mapping from the last execution
	Reco_to_True_Id_Map.clear();
	True_to_Reco_Id_Map.clear();
	m_data->CStore.Set("Reco_to_True_Id_Map",Reco_to_True_Id_Map);
	m_data->CStore.Set("True_to_Reco_Id_Map",True_to_Reco_Id_Map);
	
	// retrieve the collections of PMTs hit by true particles
	get_ok = m_data->Stores["ANNIEEvent"]->Get("ParticleId_to_MrdTubeIds", ParticleId_to_MrdTubeIds);
	if(not get_ok){
		Log("MrdEfficiency Tool: Could not find ParticleId_to_MrdTubeIds in CStore!",v_error,verbosity);
		return false;
	}
	m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
	
	// find the primary muon to check if it was reconstructed
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// loop over the MCParticles to find the highest enery primary muon
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
	if(not get_ok){
		Log("MrdEfficiency Tool: Error retrieving MCParticles,true from ANNIEEvent!",v_error,verbosity);
		return false;
	}
	// MCParticles is a std::vector<MCParticle>
	MCParticle primarymuon;  // primary muon
	bool mufound=false;
	int primarymuonid = 0;
	if(MCParticles){
		Log("MrdEfficiency Tool: Num MCParticles = "+to_string(MCParticles->size()),v_message,verbosity);
		for(int particlei=0; particlei<(int)MCParticles->size(); particlei++){
			MCParticle aparticle = MCParticles->at(particlei);
			if(aparticle.GetPdgCode()==13){
				logmessage = "True muon found with parent type " + to_string(aparticle.GetParentPdg())
					+ ", Id " + to_string(aparticle.GetParticleID())
					+ ", start vertex (" + to_string(aparticle.GetStartVertex().X())
					+ ", " + to_string(aparticle.GetStartVertex().Y())
					+ ", " + to_string(aparticle.GetStartVertex().Z())
					+ "), and end vertex (" + to_string(aparticle.GetStopVertex().X())
					+ ", " + to_string(aparticle.GetStopVertex().Y())
					+ ", " + to_string(aparticle.GetStopVertex().Z())
					+ ")";
				Log(logmessage,v_debug,verbosity);
			}
			if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
			if(aparticle.GetPdgCode()!=13) continue;       // not a muon
			primarymuon = aparticle;                       // note the particle
			primarymuonid = primarymuon.GetParticleID();   // note the ID
			mufound=true;                                  // note that we found it
			break;                                         // XXX assume we don't have more than one primary muon
		}
	} else {
		Log("MrdEfficiency Tool: No MCParticles in the event!",v_warning,verbosity);
	}
	if(not mufound){
		Log("MrdEfficiency Tool: No muon in this event",v_debug,verbosity);
		return true;
	} else {
		//Log("MrdEfficiency Tool: found muon",v_debug,verbosity);
	}
	
	// Get the true particles and which PMTs they hit
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::map<int,std::vector<int>> paddlesInTrackTrue;
	std::vector<int> trueIds;
	int npaddleshitbyprimarymuon=0;
	for(std::pair<const int,std::map<unsigned long,double>>& aparticle : *ParticleId_to_MrdTubeIds){
		//if(aparticle.first==theprimarymuonid){
		//	std::cout<<"found MrdTubeIds hit by primarymuon: {";
		//}
		std::map<unsigned long,double>* pmtshit = &aparticle.second;
		std::vector<int> tempvector;
		for(std::pair<const unsigned long,double>& apmt : *pmtshit){
			unsigned long channelkey = apmt.first;
			// map to WCSim PMT ID
			int pmtsid = channelkey_to_mrdpmtid.at(channelkey);
			tempvector.push_back(pmtsid - 1); // -1 to align with MrdTrackLib
			//if(aparticle.first==theprimarymuonid){
			//	std::cout<<(pmtsid-1)<<", ";
			//}
		}
		if(tempvector.size()){
			paddlesInTrackTrue.emplace(aparticle.first, tempvector);
			trueIds.push_back(aparticle.first);
		}
		if(aparticle.first==primarymuonid) npaddleshitbyprimarymuon = tempvector.size();
	}
	//std::cout<<"}"<<std::endl;
	
	// Get the reconstructed particles from the BoostStore
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents",numsubevs);
	if(not get_ok){
		Log("MrdEfficiency Tool: No NumMrdSubEvents in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);
	if(not get_ok){
		Log("MrdEfficiency Tool: No NumMrdTracks in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", theMrdTracks);
	if(not get_ok){
		Log("MrdEfficiency Tool: No MRDTracks member of the MRDTracks BoostStore!",v_error,verbosity);
		Log("MRDTracks store contents:",v_error,verbosity);
		m_data->Stores["MRDTracks"]->Print(false);
		return false;
	}
	if((int)theMrdTracks->size()<numtracksinev){
		cerr<<"Too few entries in MRDTracks vector relative to NumMrdTracks!"<<endl;
		// more is fine as we don't shrink for efficiency
	}
	
	std::cout<<"MRDEfficiency tool reports "<<numtracksinev<<" true tracks found"<<std::endl;
	
	// Loop over reconstructed tracks and which PMTs they hit
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::map<int,std::vector<int>> paddlesInTrackReco;
	std::vector<int> recoIds;
	for(int tracki=0; tracki<numtracksinev; tracki++){
		BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
		// Get the track details from the BoostStore
		PMTsHit.clear();
		thisTrackAsBoostStore->Get("PMTsHit",PMTsHit);
		thisTrackAsBoostStore->Get("MrdTrackID",MrdTrackID);
		if(PMTsHit.size()) paddlesInTrackReco.emplace(MrdTrackID,PMTsHit);
		recoIds.push_back(MrdTrackID);
	}
	
	logmessage = "MrdEfficiency Tool: Event had "+to_string(paddlesInTrackReco.size()) 
				 + " reconstructed tracks and " + to_string(paddlesInTrackTrue.size())
				 + " true particles that hit MRD PMTs";
	Log(logmessage,v_debug,verbosity);
	
	// Produce a matrix of figures of merit for how well each reco track matches each true particle
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<std::vector<double> > matchmerits (paddlesInTrackReco.size(), std::vector<double>(paddlesInTrackTrue.size(),0));
	for(int recotracki=0; recotracki<(int)paddlesInTrackReco.size(); recotracki++){
		// get the next reco track
		std::map<int,std::vector<int>>::iterator arecotrack = paddlesInTrackReco.begin();
		std::advance(arecotrack,recotracki);
		int RecoTrackId = arecotrack->first;
		std::vector<int> reco_tracks_paddles = arecotrack->second;
		
		logmessage = "MrdEfficiency Tool: Scanning matches against recotrack "+to_string(recotracki)
						+" which hit "+to_string(reco_tracks_paddles.size())+" tubes";
		Log(logmessage,v_debug,verbosity);
		
		// scan through the true tracks and calculate FOM against this reco track
		for(int truetracki=0; truetracki<(int)paddlesInTrackTrue.size(); truetracki++){
			// get the next true track
			std::map<int,std::vector<int>>::iterator atruetrack = paddlesInTrackTrue.begin();
			std::advance(atruetrack,truetracki);
			int TrueTrackId = atruetrack->first;
			std::vector<int> true_tracks_paddles = atruetrack->second;
			
			logmessage = "MrdEfficiency Tool: Making figure-of-merit with true track "+to_string(truetracki)
							+", which hit "+to_string(true_tracks_paddles.size())+" tubes";
			Log(logmessage,v_debug,verbosity);
			// measure how well these two sets of paddles match and assign a merit
			int num_true_paddles_matched = 0;
			int num_true_paddles_unmatched = 0;
			int num_extra_reco_paddles = 0;
			
			// scan over the reconstructed paddles
			for(int arecopaddlei=0; arecopaddlei<(int)reco_tracks_paddles.size(); arecopaddlei++){
				int therecopaddleid = reco_tracks_paddles.at(arecopaddlei);
				// see if this paddle is also in the true paddles
				bool in_true = std::count(true_tracks_paddles.begin(),true_tracks_paddles.end(),therecopaddleid);
				if(in_true) num_true_paddles_matched++;
				else num_extra_reco_paddles++;
			}
			// scan over the true paddles
			for(int atruepaddlei=0; atruepaddlei<(int)true_tracks_paddles.size(); atruepaddlei++){
				int thetruepaddleid = true_tracks_paddles.at(atruepaddlei);
				bool in_reco = std::count(reco_tracks_paddles.begin(),reco_tracks_paddles.end(),thetruepaddleid);
				if(not in_reco) num_true_paddles_unmatched++;
			}
			
			// combine the three values and the number of paddles in each track to form an FOM
			double thisfom = 0;
			thisfom +=   num_true_paddles_matched;
			thisfom += - num_true_paddles_unmatched;
			thisfom += - num_extra_reco_paddles;
			thisfom += - abs(static_cast<double>(reco_tracks_paddles.size()) -
							 static_cast<double>(true_tracks_paddles.size()));
			
			// add the merit to the matrix
			logmessage = "MrdEfficiency Tool: Final figure-of-merit for recotrack "+to_string(recotracki)
				+" with true track "+to_string(truetracki) + ": " + to_string(thisfom);
			Log(logmessage,v_debug,verbosity);
			matchmerits.at(recotracki).at(truetracki)=thisfom;
		}
	}
	
	Log("MrdEfficiency Tool: Matching reco vs true tracks",v_debug,verbosity);
	// having assigned figures of merit for how well each true track matches each reconstructed track,
	// scan over this matrix and pull out the best matches, stopping once we've run out or hit a FOM
	// lower limit
	double fomthreshold=2.0; // a pair of tracks must have at least this figure-of-merit to be matched
	
	// start of track pairing
	//~~~~~~~~~~~~~~~~~~~~~~~
	std::map<int,int> Reco_to_True_Index_Map;
	while(true){
		double currentmax=0.;
		int maxrow=-1, maxcolumn=-1;
		for(int rowi=0; rowi<(int)matchmerits.size(); rowi++){
			std::vector<double> therow = matchmerits.at(rowi);
			std::vector<double>::iterator thisrowsmaxit = std::max_element(therow.begin(), therow.end());
			if((thisrowsmaxit!=therow.end())&&((*thisrowsmaxit)>currentmax)){
				currentmax=*thisrowsmaxit; 
				maxrow = rowi;
				maxcolumn = std::distance(therow.begin(),thisrowsmaxit);
			}
		}
		
		logmessage="MrdEfficiency Tool: Next max FOM in the matrix is "+to_string(currentmax);
		Log(logmessage,v_debug,verbosity);
		// find the maximum of the matchmerits matrix, allocate that pair, mask off that row/column, and repeat
		if(currentmax>fomthreshold){
			bool makethepair=true;
			std::map<int,std::vector<int>>::iterator aparticle = paddlesInTrackReco.begin();
			std::advance(aparticle,maxrow);
			std::vector<int>* recotrackpmts = &(aparticle->second);
			aparticle = paddlesInTrackTrue.begin();
			std::advance(aparticle,maxcolumn);
			std::vector<int>* truetrackpmts = &(aparticle->second);
			
			// print the matching
			logmessage = "MrdEfficiency Tool: Candidate track "+to_string(Reco_to_True_Index_Map.size())
				+" reco has "+to_string(recotrackpmts->size())+" PMTs; {";
			for(auto&& apmt : (*recotrackpmts)){ logmessage+=to_string(apmt)+", "; }
			logmessage += 
				"}, true track has "+to_string(truetrackpmts->size()) +" PMTs; {";
			for(auto&& apmt : (*truetrackpmts)){ logmessage+=to_string(apmt)+", "; }
			logmessage += "}\n";
			Log(logmessage,v_debug,verbosity);
			
			// add this matching to the list of pairs
			if(makethepair){
				logmessage = "MrdEfficiency Tool: Matching reco track "
					+ to_string(maxrow) +" to true track "+to_string(maxcolumn)
					+ ", match FOM "+to_string(currentmax)+"\n";
				Log(logmessage,v_debug,verbosity);
				Reco_to_True_Index_Map.emplace(std::make_pair(maxrow, maxcolumn));  // TODO needs to get Ids
			}
			// remove the matched tracks from the matchmerits matrix
			matchmerits.at(maxrow).assign(matchmerits.at(maxrow).size(),-1);
			for(auto&& avec : matchmerits) avec.at(maxcolumn) = -1;
			
		} else {
			break;  // the maximum FOM between any remaining tracks is insufficient. We're done here.
		}
	}
	// end of track pairing
	//~~~~~~~~~~~~~~~~~~~~~
	
	// the map generated matches *indices* of tracks in the vectors of reco or true tracks
	// we need to convert these indices to IDs
	for(const std::pair<int,int>& amapping : Reco_to_True_Index_Map){
		Reco_to_True_Id_Map.emplace(recoIds.at(amapping.first),trueIds.at(amapping.second));
		True_to_Reco_Id_Map.emplace(trueIds.at(amapping.second),recoIds.at(amapping.first));
	}
	
	// put the matching map into the store
	// maps MC particle ids, from MCParticle::GetParticleID(),
	// to MRD track ids, from thisTrackAsBoostStore->Get("MrdTrackID",MrdTrackID);
	m_data->CStore.Set("Reco_to_True_Id_Map",Reco_to_True_Id_Map);
	m_data->CStore.Set("True_to_Reco_Id_Map",True_to_Reco_Id_Map);
	
	// now measure the efficiency by comparing how many tracks were correctly matched,
	// how many true tracks were not reconstructed, and how many reconstructed tracks had
	// no true particle (perhaps several true particles created hits that looked like a track)
	int num_correctly_matched_tracks=0;
	int num_true_tracks_not_reconstructed=0;
	int num_reco_tracks_without_match=0;
	for(std::pair<const int,std::vector<int>>& apaddleset : paddlesInTrackReco){
		int reco_track_id = apaddleset.first;
		bool is_matched = Reco_to_True_Id_Map.count(reco_track_id);
		if(is_matched) num_correctly_matched_tracks++;
		else num_reco_tracks_without_match++;
	}
	for(std::pair<const int,std::vector<int>>& apaddleset : paddlesInTrackTrue){
		int true_track_id = apaddleset.first;
		bool is_matched = True_to_Reco_Id_Map.count(true_track_id);
		if(not is_matched) num_true_tracks_not_reconstructed++;
	}
	
	logmessage = "MrdEfficiency Tool: This event had "+to_string(num_correctly_matched_tracks)
		+ " matched tracks, " + to_string(num_true_tracks_not_reconstructed)
		+ " true tracks not reconstructed and " + to_string(num_reco_tracks_without_match)
		+ " reconstructed tracks that had no recorded true particle";
	Log(logmessage,v_message,verbosity);
	
	// update the histograms
	hnumcorrectlymatched->Fill(num_correctly_matched_tracks);
	hnumtruenotmatched->Fill(num_true_tracks_not_reconstructed);
	hnumreconotmatched->Fill(num_reco_tracks_without_match);
	
	// scan the vector of MRD tubes hit by the primary muon, if any.
	// if none, this will not count toward the efficiency
	bool primarymuonhitmrd=false;
	if(ParticleId_to_MrdTubeIds->count(primarymuon.GetParticleID())==0){
		primarymuonhitmrd = false;
	} else if(ParticleId_to_MrdTubeIds->at(primarymuon.GetParticleID()).size()==0){
		primarymuonhitmrd = false;
	} else {
		primarymuonhitmrd = true;
	}
	if(not primarymuonhitmrd){
		//Log("MrdEfficiency Tool: primary muon did not hit the MRD",v_debug,verbosity);
		num_primary_muons_that_missed_MRD++;
	} else {
		//Log("MrdEfficiency Tool: muon hit the MRD",v_debug,verbosity);
		num_primary_muons_that_hit_MRD++;
		// check if we reconstructed it
		if(True_to_Reco_Id_Map.count(primarymuon.GetParticleID())){
			num_primary_muons_reconstructed++;
			
			// update the histos
			Log("MrdEfficiency Tool: filling the recod histos",v_debug,verbosity);
			hhangle_recod->Fill(primarymuon.GetTrackAngleX());
			hvangle_recod->Fill(primarymuon.GetTrackAngleY());
			htotangle_recod->Fill(primarymuon.GetTrackAngleFromBeam());
			henergyloss_recod->Fill(primarymuon.GetMrdEnergyLoss());
			htracklength_recod->Fill(primarymuon.GetTrackLengthInMrd()*100.);
			htrackpen_recod->Fill(primarymuon.GetMrdPenetration()*100.);
			hnummrdpmts_recod->Fill(npaddleshitbyprimarymuon);
			hq2_recod->Fill(0/*primarymuon.GetQ2()*/);
			// truth tank exit point
			hpep_recod->Fill(primarymuon.GetTankExitPoint().X()*100.,primarymuon.GetTankExitPoint().Z()*100.,primarymuon.GetTankExitPoint().Y()*100.);
			// truth mrd entry point
			hmpep_recod->Fill(primarymuon.GetMrdEntryPoint().X()*100.,primarymuon.GetMrdEntryPoint().Z()*100.,primarymuon.GetMrdEntryPoint().Y()*100.);
			//cout<<"back projected mrd entry recod: "; primarymuon.GetMrdEntryPoint().Print();
			// truth track endpoint (if in MRD) or MRD exit point
			htrackstop_recod->Fill(primarymuon.GetMrdExitPoint().X()*100., primarymuon.GetMrdExitPoint().Z()*100., primarymuon.GetMrdExitPoint().Y()*100.);
			//cout<<"trackstop recod: "; primarymuon.GetMrdExitPoint().Print();
			
		} else {
			num_primary_muons_not_reconstructed++;
			
			// update the histos
			Log("MrdEfficiency Tool: filling the nrecod histos",v_debug,verbosity);
			hhangle_nrecod->Fill(primarymuon.GetTrackAngleX());
			hvangle_nrecod->Fill(primarymuon.GetTrackAngleY());
			htotangle_nrecod->Fill(primarymuon.GetTrackAngleFromBeam());
			henergyloss_nrecod->Fill(primarymuon.GetMrdEnergyLoss());
			htracklength_nrecod->Fill(primarymuon.GetTrackLengthInMrd()*100.);
			htrackpen_nrecod->Fill(primarymuon.GetMrdPenetration()*100.);
			hnummrdpmts_nrecod->Fill(npaddleshitbyprimarymuon);
			//cout<<"q2"<<endl;
			hq2_nrecod->Fill(0/*primarymuon.GetQ2()*/);
			// truth tank exit point
			//cout<<"hpep"<<endl;
			hpep_nrecod->Fill(primarymuon.GetTankExitPoint().X()*100.,primarymuon.GetTankExitPoint().Z()*100.,primarymuon.GetTankExitPoint().Y()*100.);
			// truth mrd entry point
			//cout<<"hmpep"<<endl;
			hmpep_nrecod->Fill(primarymuon.GetMrdEntryPoint().X()*100.,primarymuon.GetMrdEntryPoint().Z()*100.,primarymuon.GetMrdEntryPoint().Y()*100.);
			//cout<<"back projected mrd entry not recod: "; primarymuon.GetMrdEntryPoint().Print();
			// truth track endpoint (if in MRD) or MRD exit point
			htrackstop_nrecod->Fill(primarymuon.GetMrdExitPoint().X()*100.,primarymuon.GetMrdExitPoint().Z()*100., primarymuon.GetMrdExitPoint().Y()*100.);
			//cout<<"trackstop not recod: "; primarymuon.GetMrdExitPoint().Print();
		}
	}
	
	return true;
}

bool MrdEfficiency::Finalise(){
	
	logmessage = "MrdEfficiency Tool: Counted "+to_string(num_primary_muons_that_hit_MRD)
		+ " primary muons that entered the MRD, of which "+to_string(num_primary_muons_reconstructed)
		+ " (" + to_string(100.*double(num_primary_muons_reconstructed)/double(num_primary_muons_that_hit_MRD))
		+ "%) were successfully reconstructed while "+to_string(num_primary_muons_not_reconstructed)
		+ " (" + to_string(100.*double(num_primary_muons_not_reconstructed)/double(num_primary_muons_that_hit_MRD))
		+ "%) were not reconstructed. "+to_string(num_primary_muons_that_missed_MRD)
		+" muons missed the MRD";
	Log(logmessage,v_message,verbosity);
	
	///////////////////////
	// generate the plots of Efficiency (y-axis) vs various metrics (x-axis)
	// the binned (in metric) efficiency is calculated by comparing the corresponding bin contents
	// of the 'Recod' and 'Not Recod' distributions
	
	Log("Generating efficiency plots",v_debug,verbosity);
	
	// Efficiency vs Track Length
	std::vector<double> eff_vs_mrdlength_yvec;
	std::vector<double> eff_vs_mrdlength_xvec;
	std::vector<double> eff_vs_mrdlength_yerrvec;
	std::vector<double> eff_vs_mrdlength_xerrvec;
	std::string title="Efficiency vs Track Length in MRD";
	for(int bini=1; bini<(htracklength_recod->GetNbinsX()+1); bini++){
		int recod_events = htracklength_recod->GetBinContent(bini);
		int nrecod_events = htracklength_nrecod->GetBinContent(bini);
		
		double the_efficiency = (recod_events) ? (double(recod_events)/double(recod_events+nrecod_events)) : 0.;
		eff_vs_mrdlength_yvec.push_back(the_efficiency);
		eff_vs_mrdlength_xvec.push_back(htracklength_recod->GetBinCenter(bini));
		eff_vs_mrdlength_xerrvec.push_back(htracklength_recod->GetBinWidth(bini)/2.);
		eff_vs_mrdlength_yerrvec.push_back(Efficiency_Error(recod_events,nrecod_events));
	}
	heff_vs_mrdlength = new TGraphErrors(eff_vs_mrdlength_yvec.size(), eff_vs_mrdlength_xvec.data(), eff_vs_mrdlength_yvec.data(), eff_vs_mrdlength_xerrvec.data(), eff_vs_mrdlength_yerrvec.data());
	heff_vs_mrdlength->SetName("heff_vs_mrdlength");
	heff_vs_mrdlength->SetTitle(title.c_str());
	
	//Efficiency vs Penetration
	std::vector<double> eff_vs_penetration_yvec;
	std::vector<double> eff_vs_penetration_xvec;
	std::vector<double> eff_vs_penetration_yerrvec;
	std::vector<double> eff_vs_penetration_xerrvec;
	title="Efficiency vs Track Penetration";
	for(int bini=1; bini<(htrackpen_recod->GetNbinsX()+1); bini++){
		int recod_events = htrackpen_recod->GetBinContent(bini);
		int nrecod_events = htrackpen_nrecod->GetBinContent(bini);
		
		double the_efficiency = (recod_events) ? (double(recod_events)/double(recod_events+nrecod_events)) : 0.;
		eff_vs_penetration_yvec.push_back(the_efficiency);
		eff_vs_penetration_xvec.push_back(htrackpen_recod->GetBinCenter(bini));
		eff_vs_penetration_xerrvec.push_back(htrackpen_recod->GetBinWidth(bini)/2.);
		eff_vs_penetration_yerrvec.push_back(Efficiency_Error(recod_events,nrecod_events));
	}
	heff_vs_penetration = new TGraphErrors(eff_vs_penetration_yvec.size(), eff_vs_penetration_xvec.data(), eff_vs_penetration_yvec.data(), eff_vs_penetration_xerrvec.data(), eff_vs_penetration_yerrvec.data());
	heff_vs_penetration->SetName("heff_vs_penetration");
	heff_vs_penetration->SetTitle(title.c_str());
	
	//Efficiency vs Q2
	std::vector<double> eff_vs_q2_yvec;
	std::vector<double> eff_vs_q2_xvec;
	std::vector<double> eff_vs_q2_yerrvec;
	std::vector<double> eff_vs_q2_xerrvec;
	title="Efficiency vs Event Q^2";
	for(int bini=1; bini<(hq2_recod->GetNbinsX()+1); bini++){
		int recod_events = hq2_recod->GetBinContent(bini);
		int nrecod_events = hq2_nrecod->GetBinContent(bini);
		
		double the_efficiency = (recod_events) ? (double(recod_events)/double(recod_events+nrecod_events)) : 0.;
		if(bini==1){
			std::cout<<"eff vs q2 says: nrecod events="<<recod_events<<", nrecod_events="<<nrecod_events
			<<" so total events="<<(recod_events+nrecod_events)<<" and eff="<<(double(recod_events)/double(recod_events+nrecod_events))<<std::endl;
		}
		eff_vs_q2_yvec.push_back(the_efficiency);
		eff_vs_q2_xvec.push_back(hq2_recod->GetBinCenter(bini));
		eff_vs_q2_xerrvec.push_back(hq2_recod->GetBinWidth(bini)/2.);
		eff_vs_q2_yerrvec.push_back(Efficiency_Error(recod_events,nrecod_events));
	}
	heff_vs_q2 = new TGraphErrors(eff_vs_q2_yvec.size(), eff_vs_q2_xvec.data(), eff_vs_q2_yvec.data(), eff_vs_q2_xerrvec.data(), eff_vs_q2_yerrvec.data());
	heff_vs_q2->SetName("heff_vs_q2");
	heff_vs_q2->SetTitle(title.c_str());
	
	//Efficiency vs Angle From Beam Axis
	std::vector<double> eff_vs_angle_yvec;
	std::vector<double> eff_vs_angle_xvec;
	std::vector<double> eff_vs_angle_yerrvec;
	std::vector<double> eff_vs_angle_xerrvec;
	title="Efficiency vs Track Angle From Beam";
	for(int bini=1; bini<(htotangle_recod->GetNbinsX()+1); bini++){
		int recod_events = htotangle_recod->GetBinContent(bini);
		int nrecod_events = htotangle_nrecod->GetBinContent(bini);
		
		double the_efficiency = (recod_events) ? (double(recod_events)/double(recod_events+nrecod_events)) : 0.;
		eff_vs_angle_yvec.push_back(the_efficiency);
		eff_vs_angle_xvec.push_back(htotangle_recod->GetBinCenter(bini));
		eff_vs_angle_xerrvec.push_back(htotangle_recod->GetBinWidth(bini)/2.);
		eff_vs_angle_yerrvec.push_back(Efficiency_Error(recod_events,nrecod_events));
	}
	heff_vs_angle = new TGraphErrors(eff_vs_angle_yvec.size(), eff_vs_angle_xvec.data(), eff_vs_angle_yvec.data(), eff_vs_angle_xerrvec.data(), eff_vs_angle_yerrvec.data());
	heff_vs_angle->SetName("heff_vs_angle");
	heff_vs_angle->SetTitle(title.c_str());
	
	//Efficiency vs Num MRD PMTs Hit By True Particle
	std::vector<double> eff_vs_npmts_yvec;
	std::vector<double> eff_vs_npmts_xvec;
	std::vector<double> eff_vs_npmts_yerrvec;
	std::vector<double> eff_vs_npmts_xerrvec;
	title="Efficiency vs Num MRD PMTs Hit";
	for(int bini=1; bini<(hnummrdpmts_recod->GetNbinsX()+1); bini++){
		int recod_events = hnummrdpmts_recod->GetBinContent(bini);
		int nrecod_events = hnummrdpmts_nrecod->GetBinContent(bini);
		
		double the_efficiency = (recod_events) ? (double(recod_events)/double(recod_events+nrecod_events)) : 0.;
		eff_vs_npmts_yvec.push_back(the_efficiency);
		eff_vs_npmts_xvec.push_back(hnummrdpmts_recod->GetBinCenter(bini));
		eff_vs_npmts_xerrvec.push_back(hnummrdpmts_recod->GetBinWidth(bini)/2.);
		eff_vs_npmts_yerrvec.push_back(Efficiency_Error(recod_events,nrecod_events));
	}
	heff_vs_npmts = new TGraphErrors(eff_vs_npmts_yvec.size(), eff_vs_npmts_xvec.data(), eff_vs_npmts_yvec.data(), eff_vs_npmts_xerrvec.data(), eff_vs_npmts_yerrvec.data());
	heff_vs_npmts->SetName("heff_vs_npmts");
	heff_vs_npmts->SetTitle(title.c_str());
	
	///////////////////////
	
	// now draw all the histograms and save images
	mrdEffCanv = new TCanvas("mrdEffCanv","",canvwidth,canvheight);
	mrdEffCanv->cd();
	std::string imgname;
	TLegend* leg = new TLegend(0.5,0.85,0.98,0.93,NULL,"brNDC");
	TLegendEntry* leg_entry=nullptr;
	leg->SetBorderSize(1);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	
	hhangle_nrecod->Draw();
	hhangle_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("hhangle_recod",hhangle_recod->GetTitle(),"lpf");
	leg->AddEntry("hhangle_nrecod",hhangle_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=hhangle_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	hvangle_nrecod->Draw();
	hvangle_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("hvangle_recod",hvangle_recod->GetTitle(),"lpf");
	leg->AddEntry("hvangle_nrecod",hvangle_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=hvangle_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	htotangle_nrecod->Draw();
	htotangle_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("htotangle_recod",htotangle_recod->GetTitle(),"lpf");
	leg->AddEntry("htotangle_nrecod",htotangle_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=htotangle_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	henergyloss_nrecod->Draw();
	henergyloss_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("henergyloss_recod",henergyloss_recod->GetTitle(),"lpf");
	leg->AddEntry("henergyloss_nrecod",henergyloss_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=henergyloss_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	hnummrdpmts_nrecod->Draw();
	hnummrdpmts_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("hnummrdpmts_recod",hnummrdpmts_recod->GetTitle(),"lpf");
	leg->AddEntry("hnummrdpmts_nrecod",hnummrdpmts_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=htracklength_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	hq2_nrecod->Draw();
	hq2_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("hq2_recod",hq2_recod->GetTitle(),"lpf");
	leg->AddEntry("hq2_nrecod",hq2_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=htracklength_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	htracklength_nrecod->Draw();
	htracklength_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("htracklength_recod",htracklength_recod->GetTitle(),"lpf");
	leg->AddEntry("htracklength_nrecod",htracklength_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=htracklength_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	htrackpen_nrecod->Draw();
	htrackpen_recod->Draw("same");
	leg->Clear();
	leg->AddEntry("htrackpen_recod",htrackpen_recod->GetTitle(),"lpf");
	leg->AddEntry("htrackpen_nrecod",htrackpen_nrecod->GetTitle(),"lpf");
	leg->Draw();
	imgname=htrackpen_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	hpep_recod->Draw();
	hpep_nrecod->Draw("same");
	leg->Clear();
	leg->AddEntry("hpep_recod",hpep_recod->GetTitle(),"p");
	leg->AddEntry("hpep_nrecod",hpep_nrecod->GetTitle(),"p");
	leg->Draw();
	imgname=hpep_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	hmpep_recod->Draw();
	hmpep_nrecod->Draw("same");
	leg->Clear();
	leg->AddEntry("hmpep_recod",hmpep_recod->GetTitle(),"p");
	leg->AddEntry("hmpep_nrecod",hmpep_nrecod->GetTitle(),"p");
	leg->Draw();
	imgname=hmpep_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
//	while(gROOT->FindObject("mrdEffCanv")!=nullptr){
//		std::this_thread::sleep_for (std::chrono::milliseconds(500));
//		gSystem->ProcessEvents();
//	}
//	mrdEffCanv = new TCanvas("mrdEffCanv","",canvwidth,canvheight);
//	mrdEffCanv->cd();
	htrackstop_recod->Draw();
	htrackstop_nrecod->Draw("same");
	leg->Clear();
	leg->AddEntry("htrackstop_recod",htrackstop_recod->GetTitle(),"p");
	leg->AddEntry("htrackstop_nrecod",htrackstop_nrecod->GetTitle(),"p");
	leg->Draw();
	imgname=htrackstop_recod->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	// TGraphs of efficiency vs metric
	mrdEffCanv->Clear();
	heff_vs_mrdlength->Draw("ap");
	imgname=heff_vs_mrdlength->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	mrdEffCanv->Clear();
	heff_vs_penetration->Draw("ap");
	imgname=heff_vs_penetration->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	mrdEffCanv->Clear();
	heff_vs_q2->Draw("ap");
	imgname=heff_vs_q2->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	mrdEffCanv->Clear();
	heff_vs_angle->Draw("ap");
	imgname=heff_vs_angle->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	mrdEffCanv->Clear();
	heff_vs_npmts->Draw("ap");
	imgname=heff_vs_npmts->GetTitle();
	std::replace(imgname.begin(), imgname.end(), ' ', '_');
	mrdEffCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),imgname.c_str()));
	
	// cleanup
	if(fileout) fileout->cd();
	std::vector<TObject*> histos {hhangle_recod, hvangle_recod, htotangle_recod, henergyloss_recod, htracklength_recod, htrackpen_recod, hnummrdpmts_recod, hq2_recod, hpep_recod, htrackstop_recod, hhangle_nrecod, hvangle_nrecod, htotangle_nrecod, henergyloss_nrecod, htracklength_nrecod, htrackpen_nrecod, hnummrdpmts_nrecod, hq2_nrecod, hpep_nrecod, htrackstop_nrecod, heff_vs_mrdlength, heff_vs_penetration, heff_vs_q2, heff_vs_angle, heff_vs_npmts};
	std::cout<<"deleting histograms"<<std::endl;
	for(TObject* ahisto : histos){ if(ahisto && fileout) ahisto->Write(); delete ahisto; ahisto=nullptr; }
	if(fileout!=nullptr){ fileout->Close(); delete fileout; fileout=nullptr; }
	
	Log("MrdEfficiency Tool: deleting canvas",v_debug,verbosity);
	if(gROOT->FindObject("mrdEffCanv")){ delete mrdEffCanv; mrdEffCanv=nullptr; }
	
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

double MrdEfficiency::Efficiency_Error(double n_recod_events,double n_nrecod_events){
	// variance of efficiency, from Ullrich and Xu: arXiv:physics/0701199v1
	// via Calculating Efficiency Uncertainties - Louise Heelan;
	// https://indico.cern.ch/event/66256/contribution/1/attachments/1017176/1447814/EfficiencyErrors.pdf
	double total_events       = n_recod_events+n_nrecod_events;
	double first_numerator    = ((n_recod_events+1.)*(n_recod_events+2.));
	double first_denominator  = ( (total_events+2.) * (total_events+3.) );
	double first_term         = first_numerator/first_denominator;
	double second_numerator   = pow((n_recod_events+1.),2.);
	double second_denominator = pow((total_events+2.),2.);
	double second_term        = second_numerator/second_denominator;
	return  sqrt(first_term - second_term);
}
