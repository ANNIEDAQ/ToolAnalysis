/* vim:set noexpandtab tabstop=4 wrap */
#include "WCSimDemo.h"

WCSimDemo::WCSimDemo():Tool(){}


bool WCSimDemo::Initialise(std::string configfile, DataModel &data){

	/////////////////// Usefull header ///////////////////////
	//loading config file
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	//assign transient data pointer
	m_data= &data; 
	
	// Get the Tool configuration variables
	m_variables.Get("verbosity",verbosity);
	
	// Get the geometry
	int get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
	if(not get_ok){
		Log("No AnnieGeometry in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	Log("WCSimDemo Tool: Initializing",v_message,verbosity);
	return true;
}


bool WCSimDemo::Execute(){
	Log("WCSimDemo Tool: Executing",v_debug,verbosity);
	get_ok = m_data->Stores.count("ANNIEEvent");
	if(!get_ok){
		Log("WCSimDemo Tool: No ANNIEEvent store!",v_error,verbosity);
		return false;
	};
	
	// First, see if this is a delayed trigger in the event
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCTriggernum from ANNIEEvent!",v_error,verbosity); return false; }
	// if so, truth analysis is probably not interested in this trigger. Primary muon will not be in the listed tracks.
	if(MCTriggernum>0){ Log("WCSimDemo Tool: Skipping delayed trigger",v_debug,verbosity); return true;}
	
	// Retrieve the info from the ANNIEEvent
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCFile",MCFile);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCFile from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCEventNum from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber",EventNumber);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving EventNumber from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCParticles,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits",MCHits);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving MCLAPPDHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving TDCData,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventTime",EventTime);
	if(not get_ok){ Log("WCSimDemo Tool: Error retrieving EventTime,true from ANNIEEvent!",v_error,verbosity); return false; }
	
	// Print the information. Jingbo, do what you do here. :) 
	logmessage = "WCSimDemo Tool: Processing truth tracks and digits for "+MCFile
				+", MCEvent "+to_string(MCEventNum)+", MCTrigger "+to_string(MCTriggernum);
	Log(logmessage,v_debug,verbosity);
	
	// loop over the MCParticles to find the highest enery primary muon
	// MCParticles is a std::vector<MCParticle>
	MCParticle primarymuon;  // primary muon
	bool mufound=false;
	if(MCParticles){
		Log("WCSimDemo Tool: Num MCParticles = "+to_string(MCParticles->size()),v_message,verbosity);
		for(int particlei=0; particlei<(int)MCParticles->size(); particlei++){
			MCParticle aparticle = MCParticles->at(particlei);
			//if(v_debug<verbosity) aparticle.Print();       // print if we're being *really* verbose
			if(aparticle.GetParentPdg()!=0) continue;      // not a primary particle
			if(aparticle.GetPdgCode()!=13) continue;       // not a muon
			primarymuon = aparticle;                       // note the particle
			mufound=true;                                  // note that we found it
			break;                                         // won't have more than one primary muon
		}
	} else {
		Log("WCSimDemo Tool: No MCParticles in the event!",v_error,verbosity);
	}
	if(not mufound){
		Log("WCSimDemo Tool: No muon in this event",v_warning,verbosity);
		return true;
	}
	
	// retrieve desired information from the particle
	const Position neutrinovtx = primarymuon.GetStartVertex();    // only true if the muon is primary
	const Direction muondirection = primarymuon.GetStartDirection();
	double muonenergy = primarymuon.GetStartEnergy();
	logmessage = "WCSimDemo Tool: Interaction Vertex is at ("+to_string(neutrinovtx.X())
		+", "+to_string(neutrinovtx.Y())+", "+to_string(neutrinovtx.Z())+")\n"
		+"Primary muon has energy "+to_string(muonenergy)+"GeV and direction ("
		+to_string(muondirection.X())+", "+to_string(muondirection.Y())+", "+to_string(muondirection.Z())+")";
	Log(logmessage,v_debug,verbosity);
	// see Particle.h for other information in the MCParticle class
	
	// now move to digit retrieval
	// MCHits is a std::map<unsigned long,std::vector<Hit>>
	if(MCHits){
		Log("WCSimDemo Tool: Num PMT Digits = "+to_string(MCHits->size()),v_message,verbosity);
		// iterate over the map of sensors with a measurement
		for(std::pair<unsigned long,std::vector<Hit>>&& apair : *MCHits){
			unsigned long chankey = apair.first;
			Detector* thedet = anniegeom->ChannelToDetector(chankey);
			
			if(thedet->GetDetectorElement()=="Tank"){
				std::vector<Hit>& hits = apair.second;
				for(Hit& ahit : hits){
					//if(v_message<verbosity) ahit.Print(); // << VERY verbose
					// a Hit has tubeid, charge and time
				}
			}
		} // end loop over MCHits
	} else {
		cout<<"No MCHits"<<endl;
	}
	
	// repeat for LAPPD hits
	// MCLAPPDHits is a std::map<unsigned long,std::vector<LAPPDHit>>
	if(MCLAPPDHits){
		Log("WCSimDemo Tool: Num LAPPD Digits = "+to_string(MCLAPPDHits->size()),v_message,verbosity);
		// iterate over the map of sensors with a measurement
		for(std::pair<unsigned long,std::vector<LAPPDHit>>&& apair : *MCLAPPDHits){
			unsigned long chankey = apair.first;
			Detector* thedet = anniegeom->ChannelToDetector(chankey);
			
			if(thedet->GetDetectorElement()=="LAPPD"){ // redundant
				std::vector<LAPPDHit>& hits = apair.second;
				for(LAPPDHit& ahit : hits){
					//if(v_message<verbosity) ahit.Print(); // << VERY verbose
					// an LAPPDHit has adds (global x-y-z) position, (in-tile x-y) local position
					// and time psecs
				}
			}
		} // end loop over MCLAPPDHits
	} else {
		cout<<"No MCLAPPDHits"<<endl;
	}
	
	return true;
}


bool WCSimDemo::Finalise(){
	
	return true;
}
