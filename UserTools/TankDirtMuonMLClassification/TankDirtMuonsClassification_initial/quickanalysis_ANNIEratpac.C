///////////  Reads a RAT-PAC output file and does a quick analysis ////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TVector3.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

// #include "rootstart.h"

// Header file for the classes stored in the TTree if any.
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Calib.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DSReader.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

using namespace std;
using namespace TMath;

void quick_analyzer(const char* filename_ratpac) {
	
	// ------------------------------------------------------------------------------------------- //
	// Load RAT libraries (for dsReader)
	gSystem->Load("$(RATPAC_PATH)/lib/libRATEvent.so");
	
	// Initialization
	RAT::DSReader *dsReader;
	RAT::DS::Root *ds;
	RAT::TrackNav *nav;
	RAT::TrackCursor *cursor;
	RAT::TrackNode *node;
	TChain* tri;
	TChain* runtri;
	
	RAT::DS::Run* run;
	RAT::DS::PMTInfo* pmtInfo;
	
	std::clock_t start;
	double duration;
	
	ULong64_t NbEntries;
	
	TVector3 muon_momentum;
	TVector3 unit_vect(0.,0.,1.);
	
	TVector3 init_pos;
	TVector3 fin_pos;
	Double_t init_time, fin_time;
	TString nucl_cap_pdg_code;
	
	TH1::SetDefaultSumw2();
// 	SetMyStyle();
	
	// ------------------------------------------------------------------------------------------- //
	// Starts the timer
	start = clock();
	
	gRandom = new TRandom3();
	
	// Output file
	TFile f_output("Analyzer_output.root","RECREATE");
	ofstream csvfile;
        csvfile.open ("tank_muons_data_1000evts.csv");
	//csvfile.open ("dirt_muons_data_1000evts.csv");

        int maxhits0=10000; 
      //--- write to file: ---//
	csvfile<<"particlePDGCode"<<",";
        csvfile<<"trueKE"<<",";
        csvfile<<"vtxX"<<",";
        csvfile<<"vtxY"<<",";
        csvfile<<"vtxZ"<<",";
        csvfile<<"numberofPEs"<<",";
        for (int i=0; i<133;++i){
             stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string pmt_id0= "pmt_id_";
             pmt_id0.append(temp_str);
             const char * pmtid = pmt_id0.c_str();
             csvfile<<pmtid<<",";
         }
         for (int i=0; i<133;++i){
	     stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string pmt_name0= "pmt_charge_";
             pmt_name0.append(temp_str);
             const char * pmtname = pmt_name0.c_str();
             csvfile<<pmtname<<",";
	 }
         for (int i=0; i<maxhits0;++i){
             stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string X_name= "pmtx_";
             X_name.append(temp_str);
             const char * xname = X_name.c_str();
             csvfile<<xname<<",";
         }
         for (int i=0; i<maxhits0;++i){
             stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string Y_name= "pmty_";
             Y_name.append(temp_str);
             const char * yname = Y_name.c_str();
             csvfile<<yname<<",";
         }
         for (int i=0; i<maxhits0;++i){
             stringstream strs;
             strs << i;
             string temp_str = strs.str();
             string Z_name= "pmtz_";
             Z_name.append(temp_str);
             const char * zname = Z_name.c_str();
             csvfile<<zname<<",";
         }
         for (int ii=0; ii<maxhits0;++ii){
             stringstream strs4;
             strs4 << ii;
             string temp_str4 = strs4.str();
             string T_name= "T_";
             T_name.append(temp_str4);
             const char * tname = T_name.c_str();
             csvfile<<tname<<",";
         }
         csvfile<<'\n';


	// Load the files
	dsReader = new RAT::DSReader(filename_ratpac);
	NbEntries = dsReader->GetTotal();
	
	// Load the trees 
	// load the ratpac trees and DS
	tri = new TChain("T");
	runtri = new TChain("runT");
	
	if (TString(filename_ratpac).MaybeWildcard()) {
		// Assume there is a runT in all files
		runtri->Add(filename_ratpac);
		RAT::DS::RunStore::SetReadTree(runtri);
	} else {
		// In single file case, we can check
		TFile *ftemp = TFile::Open(filename_ratpac);
		if (ftemp->Get("runT")) {
			runtri->Add(filename_ratpac);
			RAT::DS::RunStore::SetReadTree(runtri);
		} // else, no runT, so don't register runtri with RunStore
		
		delete ftemp;
	}
	
	RAT::DS::Root *branchDS = new RAT::DS::Root();
	tri->SetBranchAddress("ds", &branchDS);
	RAT::DS::RunStore::GetRun(branchDS);
	
	TH1::SetDefaultSumw2(kTRUE);
	
	// Some global variables
	Double_t prompt_window_time_low = 0, prompt_window_time_high = 200;
	Double_t delayed_window_time_low = 1000, delayed_window_time_high = 1000000;
	vector<Double_t> v_prompt_hits_times, v_delayed_hits_times;
	
	TVector3 Interaction_Vertex; TVector3 photon_pos; TVector3 Particle_Dir;
	//std::vector<double> *digitX=0; std::vector<double> *digitY=0;  std::vector<double> *digitZ=0;
        //std::vector<double> *digitt=0;

	// ------------------------------------------------------------------------------------------- //
	// Analysis loop over all the events
	for (ULong64_t entry=0; entry<NbEntries; ++entry) {
        //for (ULong64_t entry=0; entry<2; ++entry) {
		
		//cout << "Entry: " << entry << endl;
		ds = dsReader->GetEvent(entry);    
		run = RAT::DS::RunStore::Get()->GetRun(ds);
		
 		pmtInfo = run->GetPMTInfo();
		
              	// Some initializations
		double digitt[10000]={0.}; double digitX[10000]={0.}; double digitY[10000]={0.}; double digitZ[10000]={0.};
                double pmt_charge[132] = {0.}; double pmt_id[132] = {0.};
		TVector3 sensor_position; sensor_position.SetXYZ(0.,0.,0.);
		TVector3 hit_position; hit_position.SetXYZ(0.,0.,0.);
		TVector3 hit_position_global; hit_position_global.SetXYZ(0.,0.,0.);
		TVector3 sensor_direction; sensor_direction.SetXYZ(0.,0.,0.);		
		int total_hits=0; int allhits=0;	

		Interaction_Vertex = ds->GetMC()->GetMCParticle(0)->GetPosition();
		//Particle_Dir = ds->GetMC()->GetMCParticle(0)->GetDirection();
		//cout<<"Particle_Dir: "<<Particle_Dir[0]<<","<<Particle_Dir[1]<<","<<Particle_Dir[2]<<endl;
		
		// ---------- PMT loop ---------------- //
   		//cout<<" ds->GetMC()->GetMCPMTCount(): "<< ds->GetMC()->GetMCPMTCount()<<endl;
		for(long iPMT = 0; iPMT < ds->GetMC()->GetMCPMTCount(); iPMT++ ){
 	           //cout<<"iPMT= "<<iPMT<<" PMT ID "<<ds->GetMC()->GetMCPMT(iPMT)->GetID()<<endl;
                   pmt_id[iPMT] = ds->GetMC()->GetMCPMT(iPMT)->GetID();
	           pmt_charge[iPMT] = ds->GetMC()->GetMCPMT(iPMT)->GetCharge();
 		   //cout<<"pmt_charge: "<<pmt_charge[iPMT]<<endl;
		   sensor_direction.SetXYZ(pmtInfo->GetDirection(ds->GetMC()->GetMCPMT(iPMT)->GetID()).X(), 
			    		   pmtInfo->GetDirection(ds->GetMC()->GetMCPMT(iPMT)->GetID()).Y(), 
					   pmtInfo->GetDirection(ds->GetMC()->GetMCPMT(iPMT)->GetID()).Z());
                   sensor_position.SetXYZ(pmtInfo->GetPosition(ds->GetMC()->GetMCPMT(iPMT)->GetID()).X(), 
					  pmtInfo->GetPosition(ds->GetMC()->GetMCPMT(iPMT)->GetID()).Y(), 
					  pmtInfo->GetPosition(ds->GetMC()->GetMCPMT(iPMT)->GetID()).Z());
			//cout<<"sensor_position: "<<sensor_position[0]<<","<<sensor_position[1]<<","<<sensor_position[2]<<endl;
                        
			//cout<<"number of hits: "<<ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount()<<endl;
                        total_hits += ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount();
			for(long iPhot = 0; iPhot < ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount(); iPhot++){
			    allhits++;
     			    if(allhits>10000){cout<<"-----> ALARM: increase the array length!! "<<allhits<<endl;}
			    digitt[allhits] = ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime();
			    //---- Hit position in local PMT coordinates
			    photon_pos = ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetPosition();
 		            //the best position resolution you have is the PMT position so I'm addign this for hit position:
			    digitX[allhits] = sensor_position.X(); //photon_pos[0];
			    digitY[allhits] = sensor_position.Y(); //photon_pos[1];
			    digitZ[allhits] = sensor_position.Z(); //photon_pos[2];
			    //cout<<allhits<<" digit: "<<digitX[allhits]<<","<<digitY[allhits]<<","<<digitZ[allhits]<<endl;

 			    hit_position.SetXYZ(ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetPosition().X(),
			  			ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetPosition().Y(),
				         	ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetPosition().Z());

			    hit_position.RotateUz(sensor_direction);	
			    hit_position_global.SetXYZ(hit_position.X() + sensor_position.X(),
						    hit_position.Y() + sensor_position.Y(),
		  				    hit_position.Z() + sensor_position.Z());
   			    //cout<<"hit_position: "<<hit_position_global[0]<<","<<hit_position_global[1]<<","<<hit_position_global[2]<<endl;
			}
		}
		//cout<<"Total number of hits per evt is: "<<total_hits<<endl;
		// -------------------------------------- //
							
		// kind of progress bar...
		if (NbEntries > 10) {
		   if ( entry%(NbEntries/10) == 0 ) { 
		      cout<<"Evt "<<entry<<" out of "<<NbEntries <<" events ===>"<< Floor(Double_t(entry)/Double_t(NbEntries)*100.)<<" %\n";
		   }
		}
		  		  
		//----------- write to csv file ---------------------
		csvfile<<entry<<",";
		csvfile<<ds->GetMC()->GetMCParticle(0)->GetPDGCode()<<",";
		csvfile<<ds->GetMC()->GetMCParticle(0)->GetKE()<<",";
    		csvfile<<(Interaction_Vertex[0])<<",";
    		csvfile<<(Interaction_Vertex[1])<<",";
    		csvfile<<(Interaction_Vertex[2])<<",";
		csvfile<<(ds->GetMC()->GetNumPE())<<",";
                for(int i=0; i<133;++i){       
                   csvfile<<pmt_id[i]<<",";
                }   
		for(int i=0; i<133;++i){              
                   csvfile<<pmt_charge[i]<<",";
                }	       
	       for(int i=0; i<maxhits0;++i){
        	   csvfile<<digitX[i]<<",";
        	}
		for(int i=0; i<maxhits0;++i){
        	   csvfile<<digitY[i]<<",";
        	}
		for(int i=0; i<maxhits0;++i){
        	   csvfile<<digitZ[i]<<",";
        	}
	       for(int i=0; i<maxhits0;++i){
        	   csvfile<<digitt[i]<<",";
        	}
               	csvfile<<'\n';

 		//--------------------------------

		
	}
	// ------------------------------------------------------------------------------------------- //
	
	delete run;
	delete tri, runtri, branchDS;
	delete dsReader;
	
	f_output.cd();
	f_output.Close();
	
		// Ends the timer
	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Execution time: " << duration << " seconds\n";
	
}
