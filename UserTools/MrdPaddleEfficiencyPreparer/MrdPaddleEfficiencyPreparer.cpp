#include "MrdPaddleEfficiencyPreparer.h"


MrdPaddleEfficiencyPreparer::MrdPaddleEfficiencyPreparer():Tool(){}


bool MrdPaddleEfficiencyPreparer::Initialise(std::string configfile, DataModel &data){

	// Get configuration variables

	if(configfile!="")  m_variables.Initialise(configfile); //loading config file
	m_data= &data; //assigning transient data pointer

	isData = true;
	outputfile = "MRDTestFile_observed";
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("OutputFile",outputfile);
	m_variables.Get("IsData",isData);
	m_variables.Get("UseTrueTrack",usetruetrack);

	if (isData) {
		Log("MrdPaddleEfficiencyPreparer tool: Cannot use true track in data mode! Switch to non-true-track mode...",v_error,verbosity);
		usetruetrack = false;
	}

	// Define histograms for storing general properties of the fitted Mrd tracks

	std::stringstream rootfilename;
	rootfilename << outputfile << ".root";
	hist_file = new TFile(rootfilename.str().c_str(),"RECREATE");

	// Read in the geometry & get the properties (extents, orientations, etc) of the MRD paddles

	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

	std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  	for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("MRD").begin();
                                                    it != Detectors->at("MRD").end();
                                                  ++it){
	    Detector* amrdpmt = it->second;
	    unsigned long detkey = it->first;
	    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
	    Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);

	    double xmin = mrdpaddle->GetXmin();
	    double xmax = mrdpaddle->GetXmax();
	    double xmean = 0.5*(xmin+xmax);
	    double ymin = mrdpaddle->GetYmin();
	    double ymax = mrdpaddle->GetYmax();
	    double ymean = 0.5*(ymin+ymax);
	    double zmin = mrdpaddle->GetZmin();
	    double zmax = mrdpaddle->GetZmax();
	    double zmean = 0.5*(zmin+zmax);
	    int orientation = mrdpaddle->GetOrientation();    //0 is horizontal, 1 is vertical
	    int half = mrdpaddle->GetHalf();                  //0 or 1
	    int side = mrdpaddle->GetSide();
	    int layer = mrdpaddle->GetLayer();

	    map_chkey_half.emplace(chankey,half);

	    //MRD layer enumeration starts at 2, set to 0
	    layer-=2;

	    double layermin, layermax;
	    if (orientation==0) {layermin = ymin; layermax = ymax;}
	    else {layermin = xmin; layermax = xmax;}

	    if (zLayers.count(layer)==0) {
	    	zLayers.emplace(layer,zmean);
	    	orientationLayers.emplace(layer,orientation);
		std::vector<unsigned long> vec_chankey;
		vec_chankey.push_back(chankey);
		channelsLayers.emplace(layer,vec_chankey);
	    } else {
		channelsLayers.at(layer).push_back(chankey);
	    }
	    if (expected_MRDHits.count(layer)==0) {
	    	std::stringstream ss_expectedhist, ss_expectedhist_layer;
	    	std::stringstream ss_layer, ss_layer_layer;
	    	ss_expectedhist << "expectedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer << "Expected Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_expected = new TH1D(ss_expectedhist.str().c_str(),ss_layer.str().c_str(),50,layermin, layermax);
		if (orientation == 0) hist_expected->GetXaxis()->SetTitle("y [m]");
		else hist_expected->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_expected;
		map_chkey_hist_expected.emplace(chankey,hist_expected);
		expected_MRDHits.emplace(layer,map_chkey_hist_expected);
	    	ss_expectedhist_layer << "expectedhits_wholelayer"<<layer<<"_chkey"<<chankey;
	    	ss_layer_layer << "Expected Hits Whole Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_expected_layer = new TH1D(ss_expectedhist_layer.str().c_str(),ss_layer_layer.str().c_str(),200,-extents[layer],extents[layer]);
		if (orientation == 0) hist_expected_layer->GetXaxis()->SetTitle("y [m]");
		else hist_expected_layer->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_expected_layer;
		map_chkey_hist_expected_layer.emplace(chankey,hist_expected_layer);
		expected_MRDHits_layer.emplace(layer,map_chkey_hist_expected_layer);
	    } else {
		std::stringstream ss_expectedhist, ss_expectedhist_layer;
		std::stringstream ss_layer, ss_layer_layer;
		ss_expectedhist <<"expectedhits_layer"<<layer<<"_chkey"<<chankey;
		ss_layer << "Expected Hits Layer "<<layer<<" Chkey "<<chankey;
		TH1D *hist_expected = new TH1D(ss_expectedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_expected->GetXaxis()->SetTitle("y [m]");
		else hist_expected->GetXaxis()->SetTitle("x [m]");
		expected_MRDHits.at(layer).emplace(chankey,hist_expected);	
		ss_expectedhist_layer <<"expectedhits_wholelayer"<<layer<<"_chkey"<<chankey;
		ss_layer_layer << "Expected Hits Whole Layer "<<layer<<" Chkey "<<chankey;
		TH1D *hist_expected_layer = new TH1D(ss_expectedhist_layer.str().c_str(),ss_layer_layer.str().c_str(),200,-extents[layer],extents[layer]);
		if (orientation == 0) hist_expected_layer->GetXaxis()->SetTitle("y [m]");
		else hist_expected_layer->GetXaxis()->SetTitle("x [m]");
		expected_MRDHits_layer.at(layer).emplace(chankey,hist_expected_layer);	
	    }
	    if (observed_MRDHits.count(layer)==0) {
	    	std::stringstream ss_observedhist, ss_observedhist_layer;
	    	std::stringstream ss_layer, ss_layer_layer;
	    	ss_observedhist << "observedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer <<"Observed Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed = new TH1D(ss_observedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_observed->GetXaxis()->SetTitle("y [m]");
		else hist_observed->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_observed;
	    	map_chkey_hist_observed.emplace(chankey,hist_observed);
	    	observed_MRDHits.emplace(layer,map_chkey_hist_observed);
	    	ss_observedhist_layer << "observedhits_wholelayer"<<layer<<"_chkey"<<chankey;
	    	ss_layer_layer <<"Observed Hits Whole Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed_layer = new TH1D(ss_observedhist_layer.str().c_str(),ss_layer_layer.str().c_str(),200,-extents[layer],extents[layer]);
		if (orientation == 0) hist_observed_layer->GetXaxis()->SetTitle("y [m]");
		else hist_observed_layer->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_observed_layer;
	    	map_chkey_hist_observed_layer.emplace(chankey,hist_observed_layer);
	    	observed_MRDHits_layer.emplace(layer,map_chkey_hist_observed_layer);
	    } else {
	    	std::stringstream ss_observedhist, ss_observedhist_layer;
	    	std::stringstream ss_layer, ss_layer_layer;
	    	ss_observedhist << "observedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer <<"Observed Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed = new TH1D(ss_observedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_observed->GetXaxis()->SetTitle("y [m]");
		else hist_observed->GetXaxis()->SetTitle("x [m]");
	    	observed_MRDHits.at(layer).emplace(chankey,hist_observed);
	    	ss_observedhist_layer << "observedhits_wholelayer"<<layer<<"_chkey"<<chankey;
	    	ss_layer_layer <<"Observed Hits Whole Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed_layer = new TH1D(ss_observedhist_layer.str().c_str(),ss_layer_layer.str().c_str(),200,-extents[layer],extents[layer]);
		if (orientation == 0) hist_observed_layer->GetXaxis()->SetTitle("y [m]");
		else hist_observed_layer->GetXaxis()->SetTitle("x [m]");
	    	observed_MRDHits_layer.at(layer).emplace(chankey,hist_observed_layer);
	    }

	}

	hist_chankey = new TH1F("hist_chankey","MRD Chankeys - Efficiency",340,0,340);
	//hist_timediff = new TH1F("hist_timediff","Timediff first & last layer",2000,-2000,2000);
	tree_trackfit = new TTree("tree_trackfit","tree_trackfit");
	tree_trackfit->Branch("EvNum",&EventNumber);
	tree_trackfit->Branch("MRDTriggerType",&MRDTriggertype);
	tree_trackfit->Branch("NumMrdTracks",&numtracksinev);
	tree_trackfit->Branch("MrdTrackID",&MrdTrackID);
	tree_trackfit->Branch("NPMTsHit",&NPMTsHit);
	tree_trackfit->Branch("NLayersHit",&NLayersHit);
	tree_trackfit->Branch("StartVertexX",&StartVertexX);
	tree_trackfit->Branch("StartVertexY",&StartVertexY);
	tree_trackfit->Branch("StartVertexZ",&StartVertexZ);
	tree_trackfit->Branch("StopVertexX",&StopVertexX);
	tree_trackfit->Branch("StopVertexY",&StopVertexY);
	tree_trackfit->Branch("StopVertexZ",&StopVertexZ);

	m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	m_data->CStore.Get("mrd_tubeid_to_channelkey",mrdpmtid_to_channelkey);

	numtracksinev = 0;

	return true;

}

bool MrdPaddleEfficiencyPreparer::Execute(){

	// Get all relevant data from ANNIEEvent BoostStore

	m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
	m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType",MRDTriggertype);

	// In case of a cosmic event, use the fitted track for an efficiency calculation

	if (MRDTriggertype == "Cosmic"){
		
		// Get MRD track information from MRDTracks BoostStore

		m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents",numsubevs);
		m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);

		if(verbosity>2) std::cout<<"MrdPaddleEfficiencyPreparer tool: Event "<<EventNumber<<" had "<<numtracksinev<<" tracks in "<<numsubevs<<" subevents"<<endl;
			
		// Only look at events with a single track (cut may be relaxed in the future)		

//		if (numtracksinev >= 1) {
		if (numtracksinev == 1) {

			// Get reconstructed tracks

			m_data->Stores["MRDTracks"]->Get("MRDTracks",theMrdTracks);

			// Loop over reconstructed tracks
			
			for(int tracki=0; tracki<numtracksinev; tracki++){
				
				BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
				PMTsHit.clear();
				LayersHit.clear();

				//get track properties that are needed for the efficiency analysis
				
				thisTrackAsBoostStore->Get("StartVertex",StartVertex);
				thisTrackAsBoostStore->Get("StopVertex",StopVertex);
				thisTrackAsBoostStore->Get("PMTsHit",PMTsHit);
				thisTrackAsBoostStore->Get("LayersHit",LayersHit);
				thisTrackAsBoostStore->Get("MrdTrackID",MrdTrackID);

				StartVertexX = StartVertex.X();
				StartVertexY = StartVertex.Y();
				StartVertexZ = StartVertex.Z();
				StopVertexX = StopVertex.X();
				StopVertexY = StopVertex.Y();
				StopVertexZ = StopVertex.Z();
				NPMTsHit = (int) PMTsHit.size();
				NLayersHit = (int) LayersHit.size();
				tree_trackfit->Fill();

				if (usetruetrack){
					std::vector<std::pair<Position,Position>> truetrackvertices;
					std::vector<Int_t> truetrackpdgs;
					m_data->CStore.Get("TrueTrackVertices",truetrackvertices);
					m_data->CStore.Get("TrueTrackPDGs",truetrackpdgs);
					for (int i=0; i< truetrackvertices.size(); i++){
						if (i==0) continue; 	//first vertex is not saved correctly
						Position startvertex = truetrackvertices.at(i).first;
						Position stopvertex = truetrackvertices.at(i).second;
						StartVertex.SetX(startvertex.X());
						StartVertex.SetY(startvertex.Y());
						StartVertex.SetZ(startvertex.Z());
						StopVertex.SetX(stopvertex.X());
						StopVertex.SetY(stopvertex.Y());
						StopVertex.SetZ(stopvertex.Z());
					}
				}

				numpmtshit = PMTsHit.size();
				numlayershit = LayersHit.size();
				tracklength = sqrt(pow((StopVertex.X()-StartVertex.X()),2)+pow(StopVertex.Y()-StartVertex.Y(),2)+pow(StopVertex.Z()-StartVertex.Z(),2));

				// Minor selection cuts to only select well fit tracks

				if (PMTsHit.size() < 50 ){

					for (int i_pmt=0; i_pmt < (int) PMTsHit.size(); i_pmt++){
						int mrdid = PMTsHit.at(i_pmt);
						unsigned long chankey = mrdpmtid_to_channelkey[mrdid];
						hist_chankey->Fill(chankey);
						/*if (chankey < 52){
							for (int j_pmt = i_pmt; j_pmt < (int) PMTsHit.size(); j_pmt++){
								unsigned long j_chankey = mrdpmtid_to_channelkey[PMTsHit.at(j_pmt)];
								if (j_chankey > 305){
									
								}
							}
						}*/
					}

					for (int i_layer = 0; i_layer < zLayers.size(); i_layer++){
						
						// Exclude first and layer from efficiency determination						
						if (i_layer == 0 || i_layer == int(zLayers.size()) -1) continue;

							double x_layer, y_layer;
							FindPaddleIntersection(StartVertex, StopVertex, x_layer, y_layer, zLayers.at(i_layer));
							if (verbosity > 2) std::cout <<"MrdPaddleEfficiencyPreparer tool: FindPaddleIntersection found x_layer = "<<x_layer<<"& y_layer = "<<y_layer<<" for track with start position ("<<StartVertex.X()<<","<<StartVertex.Y()<<","<<StartVertex.Z()<<"), stop position ("<<StopVertex.X()<<","<<StopVertex.Y()<<","<<StopVertex.Z()<<") and z intersection point "<<zLayers.at(i_layer)<<std::endl;

							unsigned long hit_chankey=99999;
							FindPaddleChankey(x_layer, y_layer, i_layer,hit_chankey);
							if (hit_chankey == 99999) {
								std::cout <<"FindMrdPaddleEfficiencyPreparer: Did not find paddle with the desired intersection properties for this MRD layer, abort. "<<std::endl;
								continue;
							}
	
							if (orientationLayers.at(i_layer) == 0) {
								expected_MRDHits.at(i_layer).at(hit_chankey)->Fill(y_layer);
								expected_MRDHits_layer.at(i_layer).at(hit_chankey)->Fill(y_layer);
								int half_expected_ch = map_chkey_half[hit_chankey];
								int mrdid;
								if (isData) mrdid = channelkey_to_mrdpmtid[hit_chankey];
								else mrdid = channelkey_to_mrdpmtid[hit_chankey]-1;
								if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end()) {
									observed_MRDHits.at(i_layer).at(hit_chankey)->Fill(y_layer);
									observed_MRDHits_layer.at(i_layer).at(hit_chankey)->Fill(y_layer);
									
								}
								for (unsigned long i = channels_start[i_layer]; i < channels_start[i_layer+1]; i++){
									int half_channel = map_chkey_half[i];
									if (i == hit_chankey) continue;
									if (half_channel == half_expected_ch){
										if (isData) mrdid = channelkey_to_mrdpmtid[i];
                                                                		else mrdid = channelkey_to_mrdpmtid[i]-1;
										expected_MRDHits_layer.at(i_layer).at(i)->Fill(y_layer);
										if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)==PMTsHit.end()) {
                                                                		if (fabs(i-hit_chankey)<=1 || fabs(hit_chankey-i) <=1){ //Only look at neighboring channels
										if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end())  observed_MRDHits_layer.at(i_layer).at(i)->Fill(y_layer);		}								
										}
									}
								}
							}
							else {
								expected_MRDHits.at(i_layer).at(hit_chankey)->Fill(x_layer);
								expected_MRDHits_layer.at(i_layer).at(hit_chankey)->Fill(x_layer);
								int half_expected_ch = map_chkey_half[hit_chankey];
								int mrdid;
								if (isData) mrdid = channelkey_to_mrdpmtid[hit_chankey];
								else mrdid = channelkey_to_mrdpmtid[hit_chankey]-1;
								if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end()){
									observed_MRDHits.at(i_layer).at(hit_chankey)->Fill(x_layer);
									observed_MRDHits_layer.at(i_layer).at(hit_chankey)->Fill(x_layer);
								}
								for (unsigned long i = channels_start[i_layer]; i < channels_start[i_layer+1]; i++){
									int half_channel = map_chkey_half[i];
									if (i==hit_chankey) continue;
									if (half_channel == half_expected_ch){
										if (isData) mrdid = channelkey_to_mrdpmtid[i];
                                                                		else mrdid = channelkey_to_mrdpmtid[i]-1;
										expected_MRDHits_layer.at(i_layer).at(i)->Fill(x_layer);
										if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)==PMTsHit.end()){
                                                                		if (fabs(i-hit_chankey)<=1 || fabs(hit_chankey-i) <=1){
										if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end())  observed_MRDHits_layer.at(i_layer).at(i)->Fill(x_layer);										
										}
									}
									}
								}
							}
					}
				}
			}
		} else {
			tree_trackfit->Fill();
		}
	} 

	return true;

}

bool MrdPaddleEfficiencyPreparer::Finalise(){

	hist_file->cd();

	hist_chankey->Write();
        tree_trackfit->Write();
	for (unsigned int i_layer = 0; i_layer < zLayers.size(); i_layer++){
		for (unsigned int i_ch = 0; i_ch < channelsLayers.at(i_layer).size(); i_ch++){
			unsigned long temp_chkey = channelsLayers.at(i_layer).at(i_ch);
			observed_MRDHits.at(i_layer).at(temp_chkey)->Write();
			expected_MRDHits.at(i_layer).at(temp_chkey)->Write();
			observed_MRDHits_layer.at(i_layer).at(temp_chkey)->Write();
			expected_MRDHits_layer.at(i_layer).at(temp_chkey)->Write();
		}
	}

	hist_file->Close();
	delete hist_file;

	return true;


}

bool MrdPaddleEfficiencyPreparer::FindPaddleIntersection(Position startpos, Position endpos, double &x, double &y, double z){

	double DirX = endpos.X()-startpos.X();
	double DirY = endpos.Y()-startpos.Y();
	double DirZ = endpos.Z()-startpos.Z();

	if (fabs(DirZ) < 0.001) Log("MrdPaddleEfficiencyPreparer tool: StartVertex = EndVertex! Track was not fitted well",v_error,verbosity);

	double frac = (z - startpos.Z())/DirZ;

	x = startpos.X()+frac*DirX;
	y = startpos.Y()+frac*DirY;

	return true;

}

bool MrdPaddleEfficiencyPreparer::FindPaddleChankey(double x, double y, int layer, unsigned long &chankey){

	bool found_chankey = false;
	for (unsigned int i_channel = 0; i_channel < channelsLayers.at(layer).size(); i_channel++){

		if (found_chankey) break;

		unsigned long chankey_tmp = channelsLayers.at(layer).at(i_channel);
		Detector *mrdpmt = geom->ChannelToDetector(chankey_tmp);
		unsigned long detkey = mrdpmt->GetDetectorID();
		Paddle *mrdpaddle = (Paddle*) geom->GetDetectorPaddle(detkey);
		double xmin = mrdpaddle->GetXmin();
		double xmax = mrdpaddle->GetXmax();
		double ymin = mrdpaddle->GetYmin();
		double ymax = mrdpaddle->GetYmax();

		//Check if expected hit was within the channel or not
		if (xmin <= x && xmax >= x && ymin <= y && ymax >= y){
			chankey = chankey_tmp;
			found_chankey = true;
		}
	}	

	return true;

} 
