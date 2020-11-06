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
	rootfilename << outputfile << "_mrdefficiency.root";
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
	    	std::stringstream ss_expectedhist;
	    	std::stringstream ss_layer;
	    	ss_expectedhist << "expectedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer << "Expected Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_expected = new TH1D(ss_expectedhist.str().c_str(),ss_layer.str().c_str(),50,layermin, layermax);
		if (orientation == 0) hist_expected->GetXaxis()->SetTitle("y [m]");
		else hist_expected->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_expected;
		map_chkey_hist_expected.emplace(chankey,hist_expected);
		expected_MRDHits.emplace(layer,map_chkey_hist_expected);
	    } else {
		std::stringstream ss_expectedhist;
		std::stringstream ss_layer;
		ss_expectedhist <<"expectedhits_layer"<<layer<<"_chkey"<<chankey;
		ss_layer << "Expected Hits Layer "<<layer<<" Chkey "<<chankey;
		TH1D *hist_expected = new TH1D(ss_expectedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_expected->GetXaxis()->SetTitle("y [m]");
		else hist_expected->GetXaxis()->SetTitle("x [m]");
		expected_MRDHits.at(layer).emplace(chankey,hist_expected);	
	    }
	    if (observed_MRDHits.count(layer)==0) {
	    	std::stringstream ss_observedhist;
	    	std::stringstream ss_layer;
	    	ss_observedhist << "observedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer <<"Observed Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed = new TH1D(ss_observedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_observed->GetXaxis()->SetTitle("y [m]");
		else hist_observed->GetXaxis()->SetTitle("x [m]");
	    	std::map<unsigned long, TH1D*> map_chkey_hist_observed;
	    	map_chkey_hist_observed.emplace(chankey,hist_observed);
	    	observed_MRDHits.emplace(layer,map_chkey_hist_observed);
	    } else {
	    	std::stringstream ss_observedhist;
	    	std::stringstream ss_layer;
	    	ss_observedhist << "observedhits_layer"<<layer<<"_chkey"<<chankey;
	    	ss_layer <<"Observed Hits Layer "<<layer<<" Chkey "<<chankey;
	    	TH1D *hist_observed = new TH1D(ss_observedhist.str().c_str(),ss_layer.str().c_str(),50,layermin,layermax);
		if (orientation == 0) hist_observed->GetXaxis()->SetTitle("y [m]");
		else hist_observed->GetXaxis()->SetTitle("x [m]");
	    	observed_MRDHits.at(layer).emplace(chankey,hist_observed);
	    }

	}

	m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);

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

				if (usetruetrack){
					std::vector<std::pair<Position,Position>> truetrackvertices;
					std::vector<Int_t> truetrackpdgs;
					m_data->CStore.Get("TrueTrackVertices",truetrackvertices);
					m_data->CStore.Get("TrueTrackPDGs",truetrackpdgs);
					for (int i=0; i< (int) truetrackvertices.size(); i++){
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

					for (int i_layer = 0; i_layer < (int) zLayers.size(); i_layer++){
						
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
								int mrdid;
								if (isData) mrdid = channelkey_to_mrdpmtid[hit_chankey];
								else mrdid = channelkey_to_mrdpmtid[hit_chankey]-1;
								if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end())  observed_MRDHits.at(i_layer).at(hit_chankey)->Fill(y_layer);
							}
							else {
								expected_MRDHits.at(i_layer).at(hit_chankey)->Fill(x_layer);
								int mrdid;
								if (isData) mrdid = channelkey_to_mrdpmtid[hit_chankey];
								else mrdid = channelkey_to_mrdpmtid[hit_chankey]-1;
								if (std::find(PMTsHit.begin(),PMTsHit.end(),mrdid)!=PMTsHit.end()) observed_MRDHits.at(i_layer).at(hit_chankey)->Fill(x_layer);
							}
					}
				}
			}
		}
	}

	return true;

}

bool MrdPaddleEfficiencyPreparer::Finalise(){

	hist_file->cd();

	for (unsigned int i_layer = 0; i_layer < zLayers.size(); i_layer++){
		for (unsigned int i_ch = 0; i_ch < channelsLayers.at(i_layer).size(); i_ch++){
			unsigned long temp_chkey = channelsLayers.at(i_layer).at(i_ch);
			observed_MRDHits.at(i_layer).at(temp_chkey)->Write();
			expected_MRDHits.at(i_layer).at(temp_chkey)->Write();
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
