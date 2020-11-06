/* vim:set noexpandtab tabstop=4 wrap */
#include "HitResiduals.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1.h"
#include "TLegend.h"

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

HitResiduals::HitResiduals():Tool(){}

const double SPEED_OF_LIGHT=29.9792458; // cm/ns - match time and position units of PMT
const double REF_INDEX_WATER=1.42535;   // maximum Rindex gives maximum travel time

bool HitResiduals::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Usefull header ///////////////////////
	if(configfile!="")	m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	Double_t canvwidth = 700;
	Double_t canvheight = 600;
	// create the ROOT application to show histograms
	int myargc=0;
	//char *myargv[] = {(const char*)"somestring"};
	hitResidualApp = new TApplication("HitResidualApp",&myargc,0);
	htimeresidpmt = new TH1D("htimeresidpmt","PMT Digit Time residual",300,-10,25);
	htimeresidlappd = new TH1D("htimeresidlappd","LAPPD Digit Time residual",300,-10,25);
	hitResidCanv = new TCanvas("hitResidCanv","hitResidCanv",canvwidth,canvheight);
	
	get_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",anniegeom);
	if(get_ok==0){ cerr<<"No AnnieGeometry in ANNIEEVENT!"<<endl; return false; }
	
	return true;
}


bool HitResiduals::Execute(){
	
	int annieeventexists = m_data->Stores.count("ANNIEEvent");
	if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl; return false;};
	get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);
	if(get_ok==0){ cerr<<"No EventTime in ANNIEEVENT!"<<endl; return false; }
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles",MCParticles);
	if(get_ok==0){ cerr<<"No MCParticles in ANNIEEVENT!"<<endl; return false; }
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCHits",MCHits);
	if(get_ok==0){ cerr<<"No MCHits in ANNIEEVENT!"<<endl; return false; }
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits",MCLAPPDHits);
	if(get_ok==0){ cerr<<"No MCLAPPDHits in ANNIEEVENT!"<<endl; return false; }
	
	// FIRST GET THE TRIGGER TIME. IF WE HAVE A PROMPT TRIGGER, THIS SHOULD BE 0
	double triggertime = static_cast<double>(EventTime->GetNs());
	//if(triggertime!=0) cout<<"Non-zero trigger time: "<<EventTime->GetNs()<<endl;
	
	// NEXT GET THE PRIMARY MUON TIME
	bool mufound=false;
	MCParticle primarymuon;
	double muontime;
	Position muonvertex;
	if(MCParticles){
		for(int particlei=0; particlei<(int)MCParticles->size(); particlei++){
			MCParticle aparticle = MCParticles->at(particlei);
			//if(v_debug<verbosity) aparticle.Print();     // print if we're being *really* verbose
			if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
			if(aparticle.GetPdgCode()!=13) continue;       // not a muon
			primarymuon = aparticle;                       // note the particle
			mufound=true;                                  // note that we found it
			break;                                         // won't have more than one primary muon
		}
		if(mufound){
			muontime = primarymuon.GetStartTime();
			muonvertex = primarymuon.GetStartVertex();
		} else {
			cerr<<"No Primary Muon"<<endl;
			return true;
		}
	} else {
		cerr<<"No MCParticles"<<endl;
		return false;
	}
	if(muontime!=0) cout<<"Non-zero muon time: "<<muontime<<endl;
	
	// NOW LOOP OVER HITS AND CALCULATE RESIDUALS
	if(MCHits){
		cout<<"Num MCHits : "<<MCHits->size()<<endl;
		// iterate over the map of sensors with a measurement
		for(std::pair<unsigned long,std::vector<Hit>>&& apair : *MCHits){
			unsigned long chankey = apair.first;
			Detector* thistube = anniegeom->ChannelToDetector(chankey);
			Position tubeposition = thistube->GetDetectorPosition();
			double distance = sqrt(pow((tubeposition.X()-muonvertex.X()),2.)+
								   pow((tubeposition.Y()-muonvertex.Y()),2.)+
								   pow((tubeposition.Z()-muonvertex.Z()),2.));
			if(thistube->GetDetectorElement()=="Tank"){  // ADC, LAPPD, TDC
				std::vector<Hit>& hits = apair.second;
				for(Hit& ahit : hits){
					//if(v_message<verbosity) ahit.Print(); // << VERY verbose
					double hittime = ahit.GetTime();
					double expectedarrivaltime = muontime + distance*REF_INDEX_WATER/SPEED_OF_LIGHT;
					double timeresidual = hittime-expectedarrivaltime;
					htimeresidpmt->Fill(timeresidual);
					cout<<"Hittime "<<hittime<<", expected arrival time "<<expectedarrivaltime
						<<", residual "<<timeresidual<<endl;
				}
			}
		} // end loop over MCHits
	} else {
		cout<<"No MCHits"<<endl;
	}
	if(MCLAPPDHits){
		cout<<"Num MCLAPPDHits : "<<MCLAPPDHits->size()<<endl;
		// iterate over the map of sensors with a measurement
		for(std::pair<unsigned long,std::vector<LAPPDHit>>&& apair : *MCLAPPDHits){
			unsigned long chankey = apair.first;
			Detector* thistube = anniegeom->ChannelToDetector(chankey);
			Position tubeposition = thistube->GetDetectorPosition();
			double distance = sqrt(pow((tubeposition.X()-muonvertex.X()),2.)+
								   pow((tubeposition.Y()-muonvertex.Y()),2.)+
								   pow((tubeposition.Z()-muonvertex.Z()),2.));
			if(thistube->GetDetectorElement()=="LAPPD"){  // ADC, LAPPD, TDC
				std::vector<LAPPDHit>& hits = apair.second;
				for(LAPPDHit& ahit : hits){
					//if(v_message<verbosity) ahit.Print(); // << VERY verbose
					double hittime = ahit.GetTime();
					double expectedarrivaltime = muontime + distance*REF_INDEX_WATER/SPEED_OF_LIGHT;
					double timeresidual = hittime-expectedarrivaltime;
					htimeresidlappd->Fill(timeresidual);
					cout<<"LAPPDHittime "<<hittime<<", expected arrival time "<<expectedarrivaltime
						<<", residual "<<timeresidual<<endl;
				}
			}
		} // end loop over MCLAPPDHits
	} else {
		cout<<"No MCLAPPDHits"<<endl;
	}
	
	return true;
}


bool HitResiduals::Finalise(){
	
	gSystem->ProcessEvents();
	hitResidCanv->cd();
	htimeresidpmt->SetLineColor(kBlue);
	htimeresidlappd->SetLineColor(kRed);
	htimeresidpmt->Draw();
	htimeresidlappd->Draw("same");
	hitResidCanv->BuildLegend();
	hitResidCanv->SaveAs("HitTimeResiduals.png");
	while(gROOT->FindObject("hitResidCanv")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	
	// cleanup track drawing TApplication
	if(htimeresidpmt) delete htimeresidpmt;
	if(htimeresidlappd) delete htimeresidlappd;
	if(gROOT->FindObject("hitResidCanv")!=nullptr)delete hitResidCanv;
	delete hitResidualApp;
	
	return true;
}
