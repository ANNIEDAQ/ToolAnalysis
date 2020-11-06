/* vim:set noexpandtab tabstop=4 wrap */

#include "PlotLAPPDTimesFromStore.h"

// for drawing
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TPolyMarker3D.h"
#include "TSystem.h"

PlotLAPPDTimesFromStore::PlotLAPPDTimesFromStore():Tool(){}

bool PlotLAPPDTimesFromStore::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	
	if(verbosity) cout<<"Initializing Tool PlotLAPPDTimesFromStore"<<endl;
	
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// Get the Tool configuration variables
	// ====================================
	m_variables.Get("verbose",verbosity);
	
	// Get trigger window parameters from CStore
	// =========================================
	int pretriggerwindow;
	int posttriggerwindow;
	m_data->CStore.Get("WCSimPreTriggerWindow",pretriggerwindow);
	m_data->CStore.Get("WCSimPostTriggerWindow",posttriggerwindow);
	logmessage = "WCSimPreTriggerWindow="+to_string(pretriggerwindow)
				+", WCSimPostTriggerWindow="+to_string(posttriggerwindow);
	Log(logmessage,v_message,verbosity);
	
	// create the ROOT application to show histograms
	int myargc=0;
	//char *myargv[] = {(const char*)"somestring"};
	lappdRootDrawApp = new TApplication("lappdRootStoreTimesDrawApp",&myargc,0);
	digitime = new TH1D("lappdstoretimes","lappd digitimes",100,-10,50);
	digitimewithmuon = new TH1D("lappdstoretimes_withmuon","lappd digitimes for events w/ a muon",100,-10,50);
	mutime = new TH1D("muonstoretimes","muon start times from store",100,-10,50);
//	digitime = new TH1D("lappdstoretimes","lappd digitimes",100,6660,6720);
//	digitimewithmuon = new TH1D("lappdstoretimes_withmuon","lappd digitimes for events w/ a muon",100,6660,6720);
//	mutime = new TH1D("muonstoretimes","muon start times from store",100,6660,6720);
	
	// Get the anniegeom
	get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not get the AnnieGeometry from ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	return true;
}

bool PlotLAPPDTimesFromStore::Execute(){
	
	Log("Executing tool PlotLAPPDTimesFromStore",v_message,verbosity);
	
	bool allok=true;
	MCLAPPDHits=nullptr; MCParticles=nullptr;
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not find MCEventNum in ANNIEEvent!",v_error,verbosity);
		allok=false;
	}
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not find MCTriggernum in ANNIEEvent!",v_error,verbosity);
		allok=false;
	}
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventTime",EventTime);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not find EventTime in ANNIEEvent!",v_error,verbosity);
		allok=false;
	}
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not find MCLAPPDHits in ANNIEEvent!",v_error,verbosity);
		allok=false;
	}
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
	if(not get_ok){
		Log("PlotLAPPDTimesFromStore: Could not find MCParticles in ANNIEEvent!",v_error,verbosity);
		allok=false;
	}
	if(not allok){ return false; }
	
	double wcsimtriggertime = static_cast<double>(EventTime->GetNs()); // returns a uint64_t
	logmessage = "LAPPDTool: Trigger time for event "+to_string(MCEventNum)
				+", trigger "+to_string(MCTriggernum)+" was "+to_string(wcsimtriggertime)+"ns";
	Log(logmessage,v_message,verbosity);
	
	// check times relative to the time of the primary muon
	MCParticle primarymuon;  // primary muon
	bool mufound=false;
	if(MCParticles){
		// MCParticles is a std::vector<MCParticle>*
		for(int particlei=0; particlei<(int)MCParticles->size(); particlei++){
			MCParticle aparticle = MCParticles->at(particlei);
			if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
			if(aparticle.GetPdgCode()!=13) continue;       // not a muon
			primarymuon = aparticle;                       // note the particle
			mufound=true;                                  // note that we found it
			Log("primary muon at "+to_string(aparticle.GetStartTime()),v_message,verbosity);
			break;                                         // won't have more than one primary muon
		}
	} else {
		Log("WCSimDemo Tool: No MCParticles in the event!",v_message,verbosity);
	}
	
	if(MCLAPPDHits){
		Log("MCLAPPDHits has "+to_string(MCLAPPDHits->size())+" entries",v_message,verbosity);
		// MCLAPPDHits is a std::map<unsigned long,std::vector<LAPPDHit>>*
		//std::cout<<"LAPPD hits at: ";
		for(auto&& achannel : *MCLAPPDHits){
			unsigned long chankey = achannel.first;
			Detector* thedet = anniegeom->ChannelToDetector(chankey);
			int detid = thedet->GetDetectorID();
			auto& hits = achannel.second;
			logmessage = "tile "+to_string(detid)
						+" has "+to_string(hits.size())+" hits";
			Log(logmessage,v_message,verbosity);
			for(auto&& ahit : hits){
				digitime->Fill(ahit.GetTime());
				//std::cout<<ahit.GetTime()<<", ";
				if(mufound){
					digitimewithmuon->Fill(ahit.GetTime());
					mutime->Fill(primarymuon.GetStartTime());  // or GetStopTime
				}
			}
		}
		//std::cout<<std::endl;
	} else {
		Log("No MCLAPPDHits",v_warning,verbosity);
	}
	
	return true;
}


bool PlotLAPPDTimesFromStore::Finalise(){
	
	Log("plotlappdtimesfromstore finalising",v_message,verbosity);
	Double_t canvwidth = 700;
	Double_t canvheight = 600;
	lappdRootCanvas = new TCanvas("lappdTimesFromStoreCanvas","lappdTimesFromStoreCanvas",canvwidth,canvheight);
	Log("canvas created",v_message,verbosity);
	lappdRootCanvas->SetWindowSize(canvwidth,canvheight);
	lappdRootCanvas->cd();
	Log("canvas selected",v_message,verbosity);
	
	digitime->Draw();
	Log("histo drawn",v_message,verbosity);
	//lappdRootCanvas->SetLogy(1);  // spits errors about negative axis
	//gPad->SetLogy(1);
	//Log("logy set",v_message,verbosity);
	lappdRootCanvas->Update();
	lappdRootCanvas->SaveAs("lappddigitimes_instore.png");
	Log("saved, deleting hist",v_message,verbosity);
	
	digitimewithmuon->Draw();
	Log("histo 2 drawn",v_message,verbosity);
	lappdRootCanvas->Update();
	lappdRootCanvas->SaveAs("lappddigitimeswithmuon_instore.png");
	Log("saved, deleting hist",v_message,verbosity);
	
	mutime->Draw();
	Log("histo 3 drawn",v_message,verbosity);
	lappdRootCanvas->Update();
	lappdRootCanvas->SaveAs("lappdmuontimes_instore.png");
	Log("saved, deleting hist",v_message,verbosity);
	
	//lappdRootDrawApp->Run();
	//std::this_thread::sleep_for (std::chrono::seconds(5));
	//lappdRootDrawApp->Terminate(0);
	
	lappdRootCanvas->Clear();
	
	if(digitime){ delete digitime; digitime=nullptr; }
	if(digitimewithmuon){ delete digitimewithmuon; digitimewithmuon=nullptr; }
	if(mutime){ delete mutime; mutime=nullptr; }
	
	Log("deleting canvas",v_message,verbosity);
	delete lappdRootCanvas;
	Log("deleting app",v_message,verbosity);
	delete lappdRootDrawApp;
	Log("exiting",v_message,verbosity);
	return true;
}
