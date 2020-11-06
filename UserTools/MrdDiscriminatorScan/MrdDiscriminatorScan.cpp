/* vim:set noexpandtab tabstop=4 wrap */
#include "MrdDiscriminatorScan.h"

#include <boost/algorithm/string.hpp>   // for trim
#include <sys/types.h>     // for stat() test to see if file or folder
#include <sys/stat.h>
#include <thread>          // std::this_thread::sleep_for
#include <chrono>          // std::chrono::seconds
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TString.h"
#include "TApplication.h"


MrdDiscriminatorScan::MrdDiscriminatorScan():Tool(){}

bool MrdDiscriminatorScan::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// Get the Tool configuration variables
	// ====================================
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("plotDirectory",plotDirectory);
	m_variables.Get("drawHistos",drawHistos);
	m_variables.Get("filelist",filelist);
	m_variables.Get("filedir",filedir);
	m_variables.Get("Graph_X_variable",x_variable);
	
	// check the output directory exists and is suitable
	Log("MrdDiscriminatorScan Tool: Checking output directory "+plotDirectory,v_message,verbosity);
	bool isdir=false, plotDirectoryExists=false;
	struct stat s;
	if(stat(plotDirectory.c_str(),&s)==0){
		plotDirectoryExists=true;
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory
			isdir=true;                   //it's a directory
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file
			isdir=false;                  //it's a file
		} else {
			isdir=false;                  // neither file nor folder??
		}
	} else {
		plotDirectoryExists=false;
		//assert(false&&"stat failed on input path! Is it valid?"); // error
		// errors could also be because this is a file pattern: e.g. wcsim_0.4*.root
		isdir=false;
	}
	
	if(!plotDirectoryExists || !isdir){
		Log("MrdDiscriminatorScan Tool: output directory "+plotDirectory+" does not exist or is not a writable directory; please check and re-run.",v_error,verbosity);
		return false;
	} else {
		Log("MrdDiscriminatorScan Tool: output directory OK",v_debug,verbosity);
	}
	
	// If we wish to show the histograms during running, we need a TApplication
	// There may only be one TApplication, so see if another tool has already made one
	// and register ourself as a user if so. Otherwise, make one and put a pointer in the
	// CStore for other Tools
	if(drawHistos){
		Log("MrdDiscriminatorScan Tool: Checking TApplication",v_debug,verbosity);
		// create the ROOT application to show histograms
		int myargc=0;
		//char *myargv[] = {(const char*)"mrdeff"};
		// get or make the TApplication
		intptr_t tapp_ptr=0;
		get_ok = m_data->CStore.Get("RootTApplication",tapp_ptr);
		if(not get_ok){
			Log("MrdDiscriminatorScan Tool: Making global TApplication",v_error,verbosity);
			rootTApp = new TApplication("rootTApp",&myargc,0);
			tapp_ptr = reinterpret_cast<intptr_t>(rootTApp);
			m_data->CStore.Set("RootTApplication",tapp_ptr);
		} else {
			Log("MrdDiscriminatorScan Tool: Retrieving global TApplication",v_error,verbosity);
			rootTApp = reinterpret_cast<TApplication*>(tapp_ptr);
		}
		int tapplicationusers;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok) tapplicationusers=1;
		else tapplicationusers++;
		m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
	}
	
	// try to get the list of files
	Log("MrdDiscriminatorScan Tool: Opening list of files to process: "+filelist,v_message,verbosity);
	fin.open(filelist.c_str());
	if(not fin.is_open()){
		Log("MrdDiscriminatorScan Tool: could not open file list "+filelist,v_error,verbosity);
		return false;
	}
	
	// scan through the file list for the next file name
	filei=0;
	Log("MrdDiscriminatorScan Tool: Getting first file name",v_debug,verbosity);
	get_ok = GetNextFile();
	if(not get_ok){
		Log("MrdDiscriminatorScan Tool: No files found in file list!",v_error,verbosity);
		return false;
	}
	
	Log("MrdDiscriminatorScan Tool: End of Initialise",v_debug,verbosity);
	return true;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

bool MrdDiscriminatorScan::Execute(){
	
	Log("MrdDiscriminatorScan Tool: Executing",v_debug,verbosity);
	
	// Open the next file
	// ==================
	BoostStore* indata=new BoostStore(false,0);
	Log("MrdDiscriminatorScan Tool: Initializing BoostStore from file "+filepath,v_message,verbosity);
	indata->Initialise(filepath);
	
	BoostStore* MRDData= new BoostStore(false,2);
	Log("MrdDiscriminatorScan Tool: Getting CCData Multi-Entry BoostStore from loaded file",v_debug,verbosity);
	get_ok = indata->Get("CCData",*MRDData);
	if(not get_ok){
		Log("MrdDiscriminatorScan Tool: No CCData in file BoostStore!",v_error,verbosity);
		return false;
	}
	
	// Get the pulse rates for each channel
	// ====================================
	// we'll need to count how many hits we have on each channel in this file
	// make a set of maps from crate:card:channel to number of hits
	std::map<int,std::map<int,std::map<int,int>>> hit_counts_on_channels;
	TimeClass first_timestamp, last_timestamp;
	Log("MrdDiscriminatorScan Tool: Scanning hits in this file",v_debug,verbosity);
	bool success = CountChannelHits(MRDData, hit_counts_on_channels, first_timestamp, last_timestamp);
	if(not success) return false;
	// each run has 10k readouts, each readout records a 4us window (single word mode with resolution set at 4ns)
	// so in total we record 10e3 x 4e-6 = 40e-3ns of actual time
	//double total_run_seconds = 40E-3;
	unsigned long num_readouts_this_run = 0;
        get_ok = MRDData->Header->Get("TotalEntries",num_readouts_this_run);
	double total_run_seconds = 4E-6 * num_readouts_this_run;
	// get relative time of the start of this run....
	if(very_first_timestamp==0){
		very_first_timestamp = first_timestamp.GetNs();
		Log("MrdDiscriminatorScan Tool: Found start time of first run to be "+to_string(very_first_timestamp),v_debug,verbosity);
	}
	unsigned long relative_runstart_ns = 
first_timestamp.GetNs() - very_first_timestamp;
	Log("MrdDiscriminatorScan Tool: Total duration of this run was: "
		+std::to_string(total_run_seconds)+" seconds",v_debug,verbosity);
	cout<<"MrdDiscriminatorScan Tool: Run ran from " << first_timestamp.AsString() <<" to "
	    <<last_timestamp.AsString()<<endl;
	cout<<"This run started at "<<relative_runstart_ns<<" nanoseconds since the first"<<endl;
	// discriminator threshold for this run
	current_threshold = filei+15;  // in... mV, maybe?
	
	double current_x_value = (x_variable=="threshold") ? current_threshold : (relative_runstart_ns/(86400E9));
	
	// Loop over the TGraphs and set the next datapoint
	Log("MrdDiscriminatorScan Tool: Converting hit counts to rates",v_debug,verbosity);
	for( auto&& acrate : hit_counts_on_channels){
		int cratei=acrate.first;
		//cout<<"crate "<<cratei<<endl;
		for(auto&& acard : acrate.second){
			int sloti=acard.first;
			//cout<<"	slot "<<sloti<<endl;
			for(auto&& achannel : acard.second){
				// Add the rate for this channel to the TGraph of rate vs threshold for this card
				int channeli=achannel.first;
				if(cratei==7&&sloti==11&&channeli==15) continue; // RWM
				int num_hits_on_this_channel = achannel.second;
				//cout<<"		channel "<<channeli<<" had "<<num_hits_on_this_channel<<" hits"<<endl;
				double hitrate = double(num_hits_on_this_channel) / total_run_seconds;
				
				// check if this graph key exists and make it if not
				if(rategraphs.count(cratei)==0){
					rategraphs.emplace(cratei,std::map<int,std::map<int,TGraphErrors*>>{});
				}
				if(rategraphs.at(cratei).count(sloti)==0){
					rategraphs.at(cratei).emplace(sloti,std::map<int,TGraphErrors*>{});
				}
				if(rategraphs.at(cratei).at(sloti).count(channeli)==0){
					TGraphErrors* thisgraph = new TGraphErrors();
					//thisgraph->SetName(std::to_string(cratei)+"_"+std::to_string(sloti)+"_"+std::to_string(channeli));
					thisgraph->SetName(std::to_string(channeli).c_str());
					rategraphs.at(cratei).at(sloti).emplace(channeli,thisgraph);
				}
				
				// set the next datapoint on the tgraph for this channel
				logmessage = "MrdDiscriminatorScan Tool: Setting rate for "
					+to_string(cratei)+"_"+to_string(sloti)+"_"+to_string(channeli)
					+" based on "+to_string(num_hits_on_this_channel)+ " hits to "+to_string(hitrate)+"Hz";
				//Log(logmessage,v_debug,verbosity);
				TGraphErrors* thisgraph = rategraphs.at(cratei).at(sloti).at(channeli);
				if(thisgraph==nullptr){
					Log("MrdDiscriminatorScan Tool: Null TGraph pointer setting datapoint!?",v_error,verbosity);
				} else {
					thisgraph->SetPoint(thisgraph->GetN(), current_x_value, hitrate);  // XXX add errors
				}
			}
			Log("MrdDiscriminatorScan Tool: Looping to next slot",v_debug,verbosity);
		}
		Log("MrdDiscriminatorScan Tool: Looping to next crate",v_debug,verbosity);
	}
	Log("MrdDiscriminatorScan Tool: Done setting rates for this threshold",v_debug,verbosity);
	
	// Load the next file
	// ==================
	// Do this at the end of the loop so that we can break the loop now if there are no more files
	Log("MrdDiscriminatorScan Tool: Getting next file name",v_debug,verbosity);
	get_ok = GetNextFile();
	if(not get_ok){
		Log("MrdDiscriminatorScan Tool: reached end of file list.",v_error,verbosity);
		fin.close();
		m_data->vars.Set("StopLoop",1);
	} else {
		Log("MrdDiscriminatorScan Tool: Next file will be: "+filepath,v_message,verbosity);
		filei++;
	}
	
	// Cleanup??
	// ========
	Log("MrdDiscriminatorScan Tool: deleting indata and MRDData...",v_debug,verbosity);
	delete indata;
	delete MRDData;
	
	Log("MrdDiscriminatorScan Tool: Execute complete",v_debug,verbosity);
	return true;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

bool MrdDiscriminatorScan::Finalise(){
	
	Log("MrdDiscriminatorScan Tool: Processed "+std::to_string(filei+1)+" files",v_message,verbosity);
	
	// make a canvas
	int canvwidth = 700;
	int canvheight = 600;
	mrdScanCanv = new TCanvas("mrdScanCanv","",canvwidth,canvheight);
	mrdScanCanv->cd();
	
	// loop over crates
	Log("MrdDiscriminatorScan Tool: Drawing TGraphs",v_debug,verbosity);
	for(auto&& acrate : rategraphs){
		int cratei = acrate.first;
		//cout<<"crate : "<<cratei<<endl;
		// loop over card in this crate
		for(auto&& aslot : acrate.second){
			int sloti = aslot.first;
			//cout<<"	slot : "<<sloti<<endl;
			// clear the canvas (just in case), we'll plot one tdc per graph
			mrdScanCanv->Clear();
			TMultiGraph* allgraphsforthiscard = new TMultiGraph();
			for(auto&& achannel : aslot.second){
				int channeli=achannel.first;
				//cout<<"		channel : "<<channeli<<endl;
				TGraphErrors* thisgraph = rategraphs.at(cratei).at(sloti).at(channeli);
				if(thisgraph==nullptr){
					logmessage = "MrdDiscriminatorScan Tool: No TGraph to draw for crate " + to_string(cratei)
					+ ", slot " + to_string(sloti) + ", channel " + to_string(channeli);
					Log(logmessage,v_error,verbosity);
				} else {
					thisgraph->SetMarkerStyle(20);  // filled circles
					thisgraph->SetMarkerSize(0.7);  // big enough to see the fill colour
					thisgraph->SetLineStyle(0);     // no lines, makes it too busy
					thisgraph->SetFillStyle(0);
					thisgraph->SetFillColor(0); // even with no fill, this makes a border around the legend entry
					auto thecolour = gStyle->GetColorPalette((float)gStyle->GetNumberOfColors()/32*channeli);
					thisgraph->SetLineColor(thecolour);
					thisgraph->SetMarkerColor(thecolour);
					allgraphsforthiscard->Add(thisgraph,"P"); // ,"PL");
					// for sufficiently new ROOT, using TMultiGraph::Add(graph,"PL") and the option
					// TMultiGraph::Draw("same PLC PMC") will automatically pick unique colors for multiple TH1s
					// (or for THStack just "PFC nostack" when drawing the stack)
					// but we don't have sufficiently new ROOT
					std::string graphname = TString::Format("%d_%d_%d",cratei,sloti,channeli).Data();
					std::string xaxistitle = (x_variable=="threshold") ? "Threshold [mV]" : "Run Start [days]";
					std::string alltitles = graphname+";"+xaxistitle+";Pulse Rate [Hz]";
					thisgraph->GetHistogram()->SetTitle(alltitles.c_str());
					thisgraph->Draw("AP");
					//mrdScanCanv->BuildLegend();
					mrdScanCanv->Modified();
					mrdScanCanv->SaveAs(TString::Format("%s/%s.png",plotDirectory.c_str(),graphname.c_str()));

				}
				//std::cout<<"Channel "<<channeli<<" points:"<<std::endl;
				//thisgraph->Print();
			}
			Log("MrdDiscriminatorScan Tool: Drawing TMultiGraph for crate "+to_string(cratei)
				+", slot "+to_string(sloti),v_debug,verbosity);
			allgraphsforthiscard->Draw("AP");
			TLegend* leg = mrdScanCanv->BuildLegend();  // XXX position and resize
			leg-> SetNColumns(4);
			leg->SetX1NDC(0.01);
			leg->SetX2NDC(0.4);
			leg->SetY1NDC(0.25);
			leg->SetY2NDC(0.5);
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);
			// set options for the entries
			TList *legentries = leg->GetListOfPrimitives();
			TIter next(legentries);
			TObject *obj;
			TLegendEntry *alegentry;
			while ((obj = next())) {
				alegentry = (TLegendEntry*)obj;
				alegentry->SetFillStyle(0);
				alegentry->SetLineStyle(0);
			}
			// apply all settings
			gPad->Update();
			gPad->Modified();
			
			Log("MrdDiscriminatorScan Tool: Saving TMultiGraph",v_debug,verbosity);
			mrdScanCanv->SaveAs(TString::Format("%s/Crate_%i_%i.C",plotDirectory.c_str(),cratei, sloti));
			mrdScanCanv->SaveAs(TString::Format("%s/Crate_%i_%i.png",plotDirectory.c_str(),cratei, sloti));
//			// wait and allow the user to inspect
//			while(gROOT->FindObject("mrdScanCanv")!=nullptr){
//				std::this_thread::sleep_for(std::chrono::milliseconds(500));
//				gSystem->ProcessEvents();
//			}
//			mrdScanCanv = new TCanvas("mrdScanCanv","",canvwidth,canvheight);
			Log("MrdDiscriminatorScan Tool: graph cleanup",v_debug,verbosity);
			delete leg; leg=nullptr;
			delete allgraphsforthiscard; allgraphsforthiscard=nullptr;
			// the TMultiGraph owns it's contents: do not delete them individually
			Log("MrdDiscriminatorScan Tool: Looping to next slot",v_debug,verbosity);
		}
		Log("MrdDiscriminatorScan Tool: Looping to next crate",v_debug,verbosity);
	}
	
	
	// Cleanup
	// =======
	Log("MrdDiscriminatorScan Tool: Canvas cleanup",v_debug,verbosity);
	delete mrdScanCanv; mrdScanCanv=nullptr;
	
	Log("MrdDiscriminatorScan Tool: TApplication cleanup",v_debug,verbosity);
	// Cleanup TApplication if we have used it and we are the last Tool
	if(drawHistos){
		int tapplicationusers=0;
		get_ok = m_data->CStore.Get("RootTApplicationUsers",tapplicationusers);
		if(not get_ok || tapplicationusers==1){
			if(rootTApp){
				Log("MrdDiscriminatorScan Tool: deleting gloabl TApplication",v_debug,verbosity);
				delete rootTApp;
			}
		} else if (tapplicationusers>1){
			tapplicationusers--;
			m_data->CStore.Set("RootTApplicationUsers",tapplicationusers);
		}
	}
	
	Log("MrdDiscriminatorScan Tool: Finalise done",v_debug,verbosity);
	return true;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

bool MrdDiscriminatorScan::GetNextFile(){
	std::string filename="";
	bool got_file=false;
	while(getline(fin, filename)){
		boost::trim(filename);             // trim any surrounding whitespace
		if(filename.empty()){
			continue;
		} else if (filename[0] == '#'){
			continue;
		} else {
			Log("MrdDiscriminatorScan Tool: Next file will be "+filename,v_message,verbosity);
			got_file=true;
			break;
		}
	}
	if(not got_file) return false;
	filepath = filedir + "/" + filename;
	return true;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

bool MrdDiscriminatorScan::CountChannelHits(BoostStore* MRDData, std::map<int,std::map<int,std::map<int,int>>> &hit_counts_on_channels, TimeClass& first_timestamp, TimeClass& last_timestamp){
	
	// we should have 10k events per threshold
	unsigned long total_number_entries=0;
	get_ok = MRDData->Header->Get("TotalEntries",total_number_entries);
	if(not get_ok){
		Log("MrdDiscriminatorScan Tool: No TotalEntries member in the header of MRDData!",v_error,verbosity);
		return false;
	}
	logmessage = "MrdDiscriminatorScan Tool: file " + to_string(filei)
				  +" has " + std::to_string(total_number_entries) + " readouts";
	Log(logmessage,v_debug,verbosity);
	
	// loop over the readouts
	for (int readouti = 0; readouti <  (int) total_number_entries; readouti++){
		
		MRDData->GetEntry(readouti);       // MRDData is a multi-entry BoostStore: load the next event
		MRDOut mrdReadout;                 // class defined in the DataModel to contain an MRD readout
		get_ok = MRDData->Get("Data",mrdReadout);   // get the readout
		if(not get_ok){
			Log("MrdDiscriminatorScan Tool: Could not get Data member from MRDData entry "+to_string(readouti),v_error,verbosity);
			return false;
		}
		
		// each mrdReadout has a trigger number, a readout number (usually the same), a trigger timestamp,
		// and then vectors of crate, slot, channel, tdc value and a type string, one entry for each hit
		ULong64_t timestamp = mrdReadout.TimeStamp;
		TimeClass timestampclass(timestamp*1000*1000);                            // [ms] to [ns]
		if(readouti==0) first_timestamp = timestampclass;                         // [ms] to [s]
		if(readouti==(int)(total_number_entries-1)) last_timestamp = timestampclass;   // [ms] to [s]
		
		// we don't actually care about the times of the hits, we just want to count how many there were
		// so just for interest
		if(readouti == 0){
			Log("First entry of next file had timestamp: " + timestampclass.AsString(),v_debug,verbosity);
			//Log("There were " + to_string(mrdReadout.Channel.size())+" hits in this readout",v_debug,verbosity);
		}
		
		// loop over all hits in this readout
		for (int hiti = 0; hiti < (int)mrdReadout.Slot.size(); hiti++){
			
			// get the channel info
			int crate_num = mrdReadout.Crate.at(hiti);
			int slot_num = mrdReadout.Slot.at(hiti);
			int channel_num = mrdReadout.Channel.at(hiti);
			// n.b. crate 8 slot 9 channel ? is the trigger? XXX
			
			// increment the map, creating the necessary keys first if necessary
			if(hit_counts_on_channels.count(crate_num)==0){
				hit_counts_on_channels.emplace(crate_num,std::map<int,std::map<int,int>>{});
			}
			if(hit_counts_on_channels.at(crate_num).count(slot_num)==0){
				hit_counts_on_channels.at(crate_num).emplace(slot_num,std::map<int,int>{});
			}
			if(hit_counts_on_channels.at(crate_num).at(slot_num).count(channel_num)==0){
				hit_counts_on_channels.at(crate_num).at(slot_num).emplace(channel_num,0);
			}
			hit_counts_on_channels.at(crate_num).at(slot_num).at(channel_num)++;
			
//			// XXX all the following is unneeded, merely for demonstration XXX
//			unsigned int hit_time_ticks = mrdReadout.Value.at(hiti);
//			// combine the crate and slot number to a unique card index
//			int tdc_index = slot_num+(crate_num-min_crate)*100;
//			// check if this TDC contains any currently channels - 
//			// see if this tdc index is in our list of active TDC cards supplied by config file
//			std::vector<int>::iterator it = std::find(nr_slot.begin(), nr_slot.end(), tdc_index);
//			if (it == nr_slot.end()){
//				std::cout <<"Read-out Crate/Slot/Channel number not active according to configuration file."
//									<<" Check the configfile to process the data..."<<std::endl
//									<<"Crate: "<<crate_num<<", Slot: "<<slot_num<<std::endl;
//				continue;
//			}
//			// get the reduced tdc number (the n'th active TDC)
//			int tdc_num = std::distance(nr_slot.begin(),it);
//			// get a reduced channel number (the n'th active channel)
//			int abs_channel = tdc_num*num_channels+channel_num;
			
		}  // end of loop over hits this readout
		
		//clear mrdReadout vectors afterwards.... do we need to do this? XXX
		mrdReadout.Value.clear();
		mrdReadout.Slot.clear();
		mrdReadout.Channel.clear();
		mrdReadout.Crate.clear();
		mrdReadout.Type.clear();
	}  // loop to next readout
	
	return true;
}
