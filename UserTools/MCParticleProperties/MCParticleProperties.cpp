/* vim:set noexpandtab tabstop=4 wrap */
#include "MCParticleProperties.h"
#include "MRDSubEventClass.hh"      // a class for defining subevents
#include "MRDTrackClass.hh"         // a class for defining MRD tracks  - needed for cMRDTrack EnergyLoss TF1

#include "Position.h"
#include "TF1.h"

MCParticleProperties::MCParticleProperties():Tool(){}

bool MCParticleProperties::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Usefull header ///////////////////////
	if(configfile!="")  m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	m_variables.Get("verbosity",verbosity);
	
	// Get the anniegeom
	get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
	if(not get_ok){
		Log("MCParticleProperties Tool: Could not get the AnnieGeometry from ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	fidcutradius = anniegeom->GetFiducialCutRadius()*100.; // convert to cm
	fidcutz = anniegeom->GetFiducialCutZ()*100.;
	fidcuty = anniegeom->GetFiducialCutY()*100.;
	tank_radius = anniegeom->GetTankRadius()*100.;
	tank_start = (anniegeom->GetTankCentre().Z()-tank_radius)*100.;
	tank_yoffset = (anniegeom->GetTankCentre().Y())*100.;
	tank_halfheight = (anniegeom->GetTankHalfheight())*100.;
	
	// XXX override as anniegeom versions are wrong
	tank_radius = MRDSpecs::tank_radius;
	tank_start = MRDSpecs::tank_start;
	tank_yoffset = MRDSpecs::tank_yoffset;
	tank_halfheight = MRDSpecs::tank_halfheight;

	// generate pdg code to name / mass maps & save them to the CStore for further use by downstream tools
	this->GeneratePdgMap();
	this->GeneratePdgMassMap();
	m_data->CStore.Set("PdgNameMap",pdgcodetoname);
	m_data->CStore.Set("PdgMassMap",pdgcodetomass);
	
	return true;
}


bool MCParticleProperties::Execute(){
	
	Log("Tool MCParticleProperties Executing",v_message,verbosity);
	
	// first check if this is a delayed MCTrigger; if it is, we've already run
	// on this set of particles and we don't need to do it again
	uint16_t MCTriggernum;
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCTriggernum",MCTriggernum);
	if(not get_ok){
		Log("MCParticleProperties Tool: No MCTriggernum in ANNIEEvent!",v_error,verbosity);
		return false;
	} else if(MCTriggernum>0){
		return true; // nothing to do
	}
	
	// retrieve the tracks from the BoostStore
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles",MCParticles);
	if(not get_ok){
		Log("MCParticleProperties Tool: No MCParticles in ANNIEEvent!",v_error,verbosity);
		return false;
	}
	
	// Loop over reconstructed tracks and calculate additional properties
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Log("MCParticleProperties Tool: Looping over MCParticles",v_debug,verbosity);
	for(int tracki=0; tracki<(int)MCParticles->size(); tracki++){
		Log("MCParticleProperties Tool: Processing MCParticle "+to_string(tracki),v_debug,verbosity);
		// Get the track details
		MCParticle* nextparticle = &MCParticles->at(tracki);
		Log("MCParticleProperties Tool: Printing existing particle stats",v_debug,verbosity);
		if(verbosity>v_debug){
			cout<<"<<<<<<<<<<<<<<<<"<<endl;
			nextparticle->Print();
			cout<<"<<<<<<<<<<<<<<<<"<<endl;
		}
		Position startvertex = nextparticle->GetStartVertex();
		startvertex.UnitToCentimeter();
		Position stopvertex = nextparticle->GetStopVertex();
		stopvertex.UnitToCentimeter();
		
		Position differencevector = (stopvertex-startvertex);
		double atracklengthtotal = differencevector.Mag();
		// XXX DEBUG THIS, why isn't this returning the correct answer?
		//double atrackangle = differencevector.Angle(Position(0,0,1));
		
		// check if it starts in the fiducial volume
		bool isinfiducialvol=false;
		if( (sqrt(pow(startvertex.X(), 2) 
				+ pow(startvertex.Z()-tank_start-tank_radius,2)) < fidcutradius) && 
			(abs(startvertex.Y()-tank_yoffset) < fidcuty) && 
			((startvertex.Z()-tank_start-tank_radius) < fidcutz) ){
				isinfiducialvol=true;
		}
		nextparticle->SetStartsInFiducialVolume(isinfiducialvol);
		
		//====================================================================================================
		// Estimate MRD penetration TODO add mrd entry/exit points to WCSim Tracks?
		//====================================================================================================
		
		// The mrd is as wide as the tank. We can ensure a track enters the mrd by projecting
		// the track forward from the start vertex at the angle between start and stop vertices,
		// and requiring:
		// 1) at z=MRD_start, x is between +/-(MRD_width/2);
		// 2) z endpoint is > MRD_start
		// For at least 2 layers of penetration, as above but with requirement on x @ z=MRD_layer2
		// For range-out, as above with requirement on x @ z=MRD_start+MRD_depth && 
		//   z endpoint is > MRD_start+MRD_depth
		
		float oppx = stopvertex.X() - startvertex.X();
		float adj = stopvertex.Z() - startvertex.Z();
		double avgtrackgradx = (adj!=0) ? (oppx/adj) : 1000000;
		float avgtrackanglex = atan(avgtrackgradx);
		float oppy = stopvertex.Y() - startvertex.Y();
		double avgtrackgrady = (adj!=0) ? (oppy/adj) : 1000000;
		float avgtrackangley = atan(avgtrackgrady);
		// XXX debug override until we fix version above
		double atrackangle = sqrt(pow(avgtrackanglex,2.)+pow(avgtrackangley,2.));
		
		float xatmrd = startvertex.X() + (MRDSpecs::MRD_start-startvertex.Z())*avgtrackgradx;
		float yatmrd = startvertex.Y() + (MRDSpecs::MRD_start-startvertex.Z())*avgtrackgrady;
		
		if(verbosity>3){
			cout<<"primary start vertex: "; startvertex.Print();
			cout<<"primary stop vertex: "; stopvertex.Print();
			cout<<"oppx = "<<oppx<<endl;
			cout<<"adj = "<<adj<<endl;
			cout<<"anglex = "<<avgtrackanglex<<endl;
			cout<<"tan(anglex) = "<<avgtrackgradx<<endl;
			cout<<"projected x at z="<<MRDSpecs::MRD_start<<" is "<<xatmrd<<endl;
			cout<<"xatmrd="<<xatmrd<<", MRD_width="<<MRDSpecs::MRD_width<<endl;
			cout<<"yatmrd="<<yatmrd<<", MRD_height="<<MRDSpecs::MRD_height<<endl;
		}
		
		// variables to be saved to EventDistributions
		Position MRDentrypoint(0,0,0), MRDexitpoint(0,0,0), MuTrackInMRD(0,0,0);
		double atracklengthinmrd, mrdpenetrationcm;
		int mrdpenetrationlayers;
		bool atrackprojectedhitmrd;
		bool atrackentersmrd, atrackstopsinmrd, atrackpenetratesmrd;
		// check for intercept and record entry point
		bool checkboxlinerror=false;

		//======================================================================================================
		// Calculate whether the extended (projected) particle trajectory would hit the MRD (first layer)
		//======================================================================================================

		atrackprojectedhitmrd  =  CheckProjectedMRDHit(startvertex, stopvertex, MRDSpecs::MRD_width, MRDSpecs::MRD_height,MRDSpecs::MRD_start);
		
		//new version based on external function calls
		///////////////////////////////////////////////
		// bool CheckLineBox( Position L1, Position L2, Position B1, Position B2, 
		//					  Position &Hit, Position &Hit2, bool &error)
		// returns true if line (L1, L2) intersects with the box (B1, B2)
		// returns intersection with smaller Z in Hit
		// if 2 interceptions are found, returns interception with larger Z,
		// if 1 interception is found, returns L2 (stopping point).
		// error returns true if >2 intercepts are found, or other error.
		atrackentersmrd  =  CheckLineBox(startvertex, stopvertex, 
										Position(-MRDSpecs::MRD_width,-MRDSpecs::MRD_height,MRDSpecs::MRD_start), 
										Position(MRDSpecs::MRD_width,MRDSpecs::MRD_height,MRDSpecs::MRD_end),
										MRDentrypoint, MRDexitpoint, checkboxlinerror);
		// sanity check: DISABLE TO ALLOW TRACKS STARTING IN THE MRD
		//assert(MRDentrypoint!=startvertex&&"track starts in MRD!?");
		// check if MRD stops in the MRD
		atrackstopsinmrd = ( abs(stopvertex.X())<MRDSpecs::MRD_width&&
							 abs(stopvertex.Y())<MRDSpecs::MRD_height&&
							 stopvertex.Z()>MRDSpecs::MRD_start&&
							 stopvertex.Z()<MRDSpecs::MRD_end );
		if(atrackentersmrd){
			atrackpenetratesmrd = ( (MRDentrypoint.Z()==MRDSpecs::MRD_start)&&
									(MRDexitpoint.Z()==MRDSpecs::MRD_end) );
			MuTrackInMRD = (MRDexitpoint-MRDentrypoint);
			atracklengthinmrd = MuTrackInMRD.Mag();
			mrdpenetrationcm = MuTrackInMRD.Z();
			mrdpenetrationlayers=0;
			for(auto layerzval : MRDSpecs::mrdscintlayers){
				if(MRDexitpoint.Z()<layerzval) break;
				mrdpenetrationlayers++;
			}
		} else {
			atrackpenetratesmrd=false;
			mrdpenetrationcm=0.;
			mrdpenetrationlayers=0;
			atracklengthinmrd=0.;
		}
		
		double muxtracklength=MuTrackInMRD.X();
		double muytracklength=MuTrackInMRD.Y();
		double muztracklength=MuTrackInMRD.Z();
		if(verbosity>3){
			cout<<"particle travels "<<atracklengthinmrd<<"cm before exiting MRD bounds"<<endl
				<<"particle travels "<<muxtracklength<<"cm before leaving X bounds"<<endl
				<<"particle travles "<<muytracklength<<"cm before leaving Y bounds"<<endl;
		}
		double maxtracklengthinMRD = sqrt(
				pow(MRDSpecs::MRD_width*2.,2.) +
				pow(MRDSpecs::MRD_height*2.,2.) +
				pow(MRDSpecs::MRD_depth,2.) );
		if((atracklengthinmrd>maxtracklengthinMRD)||
			(muxtracklength>((MRDSpecs::MRD_width*2)*1.01)) || (muytracklength>((MRDSpecs::MRD_height*2)*1.01))||
			(mrdpenetrationcm>(MRDSpecs::MRD_depth*1.01))   || checkboxlinerror ){ // float inaccuracies...
			cerr<<"MRD track length is bad!"<<endl
				<<"Max track length is "<<maxtracklengthinMRD<<"cm"
				<<", this track is "<<atracklengthinmrd<<"cm, "
				<<"distances are ("<<muxtracklength<<", "
				<<muytracklength<<", "<<mrdpenetrationcm<<")"<<endl
				<<"compare to maximum extents ("<<(2*MRDSpecs::MRD_width)
				<<", "<<(2*MRDSpecs::MRD_height)<<", "<<MRDSpecs::MRD_depth<<")"<<endl
				<<"MRD Track goes = ("<<MRDentrypoint.X()<<", "<<MRDentrypoint.Y()<<", "
				<<MRDentrypoint.Z()<<") -> ("<<MRDexitpoint.X()<<", "<<MRDexitpoint.Y()<<", "
				<<MRDexitpoint.Z()<<") and ";
				if(atrackstopsinmrd) cerr<<"stops in the MRD"<<endl;
				else if(atrackpenetratesmrd) cerr<<"penetrates the MRD"<<endl;
				else cerr<<"exits the side of the MRD"<<endl;
			cerr<<"MRD width is "<<MRDSpecs::MRD_width<<", MRD height is "<<MRDSpecs::MRD_height
				<<", MRD start is "<<MRDSpecs::MRD_start<<", MRD end is "<<(MRDSpecs::MRD_end)<<endl
				<<"total path is from ("<<startvertex.X()<<", "<<startvertex.Y()
				<<", "<<startvertex.Z()<<") -> ("<<stopvertex.X()<<", "
				<<stopvertex.Y()<<", "<<stopvertex.Z()<<")"<<endl
				<<"avgtrackangley="<<avgtrackangley<<", avgtrackanglex="<<avgtrackanglex<<endl
				<<"CheckLineBox error is "<<checkboxlinerror<<endl;
			
			//return false;
		}
		
		// transfer to ANNIEEvent
		// ~~~~~~~~~~~~~~~~~~~~~~
		nextparticle->SetTrackAngleX(avgtrackanglex);
		nextparticle->SetTrackAngleY(avgtrackangley);
		nextparticle->SetTrackAngleFromBeam(atrackangle);
		nextparticle->SetProjectedHitMrd(atrackprojectedhitmrd);
		nextparticle->SetEntersMrd(atrackentersmrd);
		MRDentrypoint.UnitToMeter();
		nextparticle->SetMrdEntryPoint(MRDentrypoint);
		MRDexitpoint.UnitToMeter();
		nextparticle->SetExitsMrd(!atrackstopsinmrd);
		nextparticle->SetMrdExitPoint(MRDexitpoint);
		nextparticle->SetPenetratesMrd(atrackpenetratesmrd);
		nextparticle->SetMrdPenetration(mrdpenetrationcm/100.);
		nextparticle->SetNumMrdLayersPenetrated(mrdpenetrationlayers);
		nextparticle->SetTrackLengthInMrd(atracklengthinmrd/100.);
		double dEdx = cMRDTrack::MRDenergyvspenetration.Eval(atrackangle);
		double energylossinmrd = (atracklengthinmrd) ? (atracklengthinmrd*dEdx) : 0;
		nextparticle->SetMrdEnergyLoss(energylossinmrd);
		if(verbosity>3){
			cout<<"atrackangle: "<<atrackangle<<endl
				<<"atrackprojectedhitmrd: "<<atrackprojectedhitmrd<<endl
				<<"atrackentersmrd: "<<atrackentersmrd<<endl
				<<"atrackstopsinmrd: "<<atrackstopsinmrd<<endl
				<<"mrdpenetrationcm: "<<mrdpenetrationcm<<endl
				<<"mrdpenetrationlayers: "<<mrdpenetrationlayers<<endl
				<<"atracklengthinmrd: "<<atracklengthinmrd<<endl
				<<"energylossinmrd: "<<(energylossinmrd)<<endl
				<<"MRD track goes = ("<<MRDentrypoint.X()<<", "<<MRDentrypoint.Y()<<", "
				<<MRDentrypoint.Z()<<") -> ("<<MRDexitpoint.X()<<", "<<MRDexitpoint.Y()<<", "
				<<MRDexitpoint.Z()<<")"<<endl;
		}
		
		//====================================================================================================
		// calculate the track length in water
		//====================================================================================================
		// to calculate track length _in water_ find distance from start vertex to the point
		// where it intercepts the tank. if this length > total track length, use total track length
		// otherwise use this length. 
		
		// first check if the start and endpoints are in the tank - if so, there is no tank exit point
		// and track length in tank is total length. 
		Log("MCParticleProperties Tool: Estimating track length in tank",v_debug,verbosity);
		double atracklengthintank=0;
		bool interceptstank=true;
		bool hastrueexitvtx=false;
		Position tankentryvtx, tankexitvtx, truetankexitvtx;
		
		bool trackstartsintank = ( sqrt(pow(startvertex.X(),2.) +
										pow(startvertex.Z()-tank_start-tank_radius,2.))
										 < tank_radius ) &&
								 ( abs(startvertex.Y()-tank_yoffset)
										 < tank_halfheight );
		bool trackstopsintank = ( sqrt(pow(stopvertex.X(),2.) +
										pow(stopvertex.Z()-tank_start-tank_radius,2.))
										 < tank_radius ) &&
								 ( abs(stopvertex.Y()-tank_yoffset)
										 < tank_halfheight );
		
//		// if TankExitPoint is not close to either the tank caps or tank walls,
//		// we probably don't have one. So project it instead.
//		double tankexitradius = sqrt(pow(nextparticle->GetTankExitPoint().X(),2.) +
//		                             pow(nextparticle->GetTankExitPoint().Z(),2.));
//		double tankexitheight = abs(nextparticle->GetTankExitPoint().Y());
//		if((tankexitradius<(tank_radius*0.8))&&(tankexitheight<(tank_height*0.8))){
//			
		if(nextparticle->GetTankExitPoint().Mag()>0.2){ // cover rounding errors
			// Later versions of WCSim record the tank exit point explicitly
			truetankexitvtx = nextparticle->GetTankExitPoint();
			truetankexitvtx.UnitToCentimeter();
			hastrueexitvtx=true;
		}
		
		if(trackstartsintank&&(trackstopsintank||hastrueexitvtx)){
			// If we have know true start and endpoints in tank, we can calculate track length directly
			logmessage="MCParticleProperties Tool: Track endpoints in tank are both known";
			Log(logmessage,v_debug,verbosity);
			tankentryvtx = startvertex;
			tankexitvtx = stopvertex;
			
		} else {
			logmessage="MCParticleProperties Tool: Estimating tank track length";
			Log(logmessage,v_debug,verbosity);
			// Three potential cases here:
			// 1. Start is outside the tank (entry unknown), exit point is known (either inside tank or recorded)
			//     - we need to estimate entry point
			// 2. Start is inside tank, but endpoint is outside the tank and exit was not recorded by WCSim
			//     - we need to estimate the exit point
			// 3. Both start and endpoints are outside the tank and unknown
			//     - we need to check the track intercepts the tank, and estimate both entry & exit points
			// If we have a recorded truth exit point we can use it as the endpoint for CheckTankIntercepts
			// to get a more accurate estimate
			
			Position theendvertex = (hastrueexitvtx) ? truetankexitvtx : stopvertex;
			interceptstank = CheckTankIntercepts(startvertex, theendvertex, trackstartsintank, trackstopsintank, tankexitvtx, tankentryvtx);
			if(verbosity>3){
				cout<<"checktankintercepts returned "<<interceptstank<<endl;
				cout<<"and set tankexitvtx to: ("<<tankexitvtx.X()
					<<", "<<tankexitvtx.Y()<<", "<<tankexitvtx.Z()<<")"<<endl;
				cout<<"and set tankentryvtx to ("<<tankentryvtx.X()<<", "<<tankentryvtx.Y()
					<<", "<<tankentryvtx.Z()<<")"<<endl;
			}
			if(trackstartsintank) tankentryvtx=startvertex;
			if(trackstopsintank) tankexitvtx=stopvertex;
		}
		if(hastrueexitvtx) tankexitvtx = truetankexitvtx;
		
		// we're now able to determine particle track length in the tank:
		atracklengthintank = (interceptstank) ? (tankexitvtx-tankentryvtx).Mag() : 0;
		
		// Sanity checks
		double maxtanktracklength = 
		sqrt(pow(tank_radius*2.,2.)+pow(tank_halfheight*2.,2.));
		if(verbosity>3){
			cout<<"TRACK SUMMARY"<<endl;
			cout<<"track start point : ("<<startvertex.X()<<", "<<startvertex.Y()
				<<", "<<startvertex.Z()<<")"<<endl
				<<"track stop point : ("<<stopvertex.X()<<", "<<stopvertex.Y()
				<<", "<<stopvertex.Z()<<")"<<endl;
			cout<<"track starts "<<((trackstartsintank) ? "inside" : "outside")<<" tank, and ends "
				<<((trackstopsintank) ? "inside" : "outside")<<" the tank"<<endl;
			cout<<"tank entry point : ("<<tankentryvtx.X()<<", "<<tankentryvtx.Y()
				<<", "<<tankentryvtx.Z()<<")"<<endl;
			cout<<"tank exit point: ("<<tankexitvtx.X()<<", "<<tankexitvtx.Y()<<", "<<tankexitvtx.Z()<<") "
				<<endl;
			cout<<"tank "<<((hastrueexitvtx) ? "had" : "didn't have")<<" a truth recorded exit point"<<endl;
			cout<<"tank track length: ("<<(tankexitvtx.X()-tankentryvtx.X())
				<<", "<<(tankexitvtx.Y()-tankentryvtx.Y())<<", "
				<<(tankexitvtx.Z()-tankentryvtx.Z())<<") = "<<atracklengthintank<<"cm total"<<endl;
			cout<<"c.f. max possible tank track length is "<<maxtanktracklength<<endl;
		}
		if(atracklengthintank > maxtanktracklength){
			cerr<<"MCParticleProperties Tool: Track length is impossibly long!"<<endl;
			//return false;
		}
		if(atracklengthintank > differencevector.Mag()){
			cerr<<"MCParticleProperties Tool: Track length in tank is greater than total track length"<<endl;
			//return false;
		}
		if(TMath::IsNaN(atracklengthintank)){
			cerr<<"MCParticleProperties Tool: NaN RESULT FROM MU TRACK LENGTH IN TANK?!"<<endl;
			//return false;
		}
		
		// We can also use a slightly modified version of the CheckTankIntercepts to find
		// the projected tank exit point for particles that started in the tank, regardless
		// of whether they left it. This is useful for e.g. showing pion decay gamma vertices,
		// even if the light is not produced by the gammas themselves
		Log("MCParticleProperties Tool: Calculating projected tank exit",v_debug,verbosity);
		Position projectedexitvertex;
		if(hastrueexitvtx||(interceptstank&&(!trackstopsintank))){
			Log("MCParticleProperties Tool: Using true/estimated tank exit",v_debug,verbosity);
		} else {
			// exit projection only makes sense for tracks starting in tank
			if(trackstartsintank && trackstopsintank){
			bool projectok = ProjectTankIntercepts(startvertex, stopvertex, projectedexitvertex);
				if(not projectok){
					cerr<<"MCParticleProperties Tool: ProjectTankIntercepts returned false?!"<<std::endl;
				} else {
					Log("MCParticleProperties Tool: Projected track exit OK",v_debug,verbosity);
				}
			} else {
				// TODO we could relax the requirement that it started in the tank;
				// provided it either started in the tank *or* had a tankentryvtx
				// we can project an exit. Leave that for now though. FIXME
				Log("MCParticleProperties Tool: Track did not start in tank; no projected exit",v_debug,verbosity);
			}
		}
		
		// transfer to ANNIEEvent
		// ~~~~~~~~~~~~~~~~~~~~~~
		nextparticle->SetEntersTank(interceptstank);
		tankentryvtx.UnitToMeter();
		nextparticle->SetTankEntryPoint(tankentryvtx);
		nextparticle->SetExitsTank(interceptstank&&!trackstopsintank);
		// TankExitPoint is either a projected or true exit point
		if(hastrueexitvtx){
			// do nothing: the tankexitvtx is already populated by LoadWCSim
		}
		else if(interceptstank&&(!trackstopsintank)){
			tankexitvtx.UnitToMeter();
			nextparticle->SetTankExitPoint(tankexitvtx);
		} else {
			projectedexitvertex.UnitToMeter();
			nextparticle->SetTankExitPoint(projectedexitvertex);
		}
		nextparticle->SetTrackLengthInTank(atracklengthintank/100.);
		if(verbosity>3){
			cout<<"tankentryvtx: ("<<tankentryvtx.X()<<", "<<tankentryvtx.Y()<<", "<<tankentryvtx.Z()<<")"
				<<", tankexitvtx: ("<<tankexitvtx.X()<<", "<<tankexitvtx.Y()<<", "<<tankexitvtx.Z()<<")"
				<<", atracklengthintank: "<<atracklengthintank<<endl;
		}
		
	} // end of loop over MCParticles
	
	return true;
}


bool MCParticleProperties::Finalise(){
	
	return true;
}


//============================================================================

// a test to see if a projected point in a plane is within a box in that plane
int MCParticleProperties::InBox( Position Hit, Position B1, Position B2, const int Axis) {
	if ( Axis==1 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	if ( Axis==2 && Hit.Z() > B1.Z() && Hit.Z() < B2.Z() && Hit.X() > B1.X() && Hit.X() < B2.X()) return 1;
	if ( Axis==3 && Hit.X() > B1.X() && Hit.X() < B2.X() && Hit.Y() > B1.Y() && Hit.Y() < B2.Y()) return 1;
	return 0;
}

// projects the hitpoint by adding a scaled vector to the start point
int MCParticleProperties::GetIntersection( float fDst1, float fDst2, Position P1, Position P2, Position &Hit) {
	if ( (fDst1 * fDst2) >= 0.0f) return 0;
	if ( fDst1 == fDst2) return 0;
	Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
	return 1;
}

// returns true if line (L1, L2) intersects with the box (B1, B2)
// returns intersection point in Hit
bool MCParticleProperties::CheckLineBox( Position L1, Position L2, Position B1, Position B2, Position &Hit, Position &Hit2, bool &error, int verbose){
	error=false;
	
	if(verbose){
		cout<<"CheckLineBox called with start:"; L1.Print(); cout<<" and stop "; L2.Print();
		cout<<", c.f. box corners are "; B1.Print(); cout<<" and "; B2.Print(); cout<<endl;
	}
	
	// check if it misses the box entirely by being on one side of a plane over entire track
	if (L2.X() <= B1.X() && L1.X() <= B1.X()) return false;
	if (L2.X() >= B2.X() && L1.X() >= B2.X()) return false;
	if (L2.Y() <= B1.Y() && L1.Y() <= B1.Y()) return false;
	if (L2.Y() >= B2.Y() && L1.Y() >= B2.Y()) return false;
	if (L2.Z() <= B1.Z() && L1.Z() <= B1.Z()) return false;
	if (L2.Z() >= B2.Z() && L1.Z() >= B2.Z()) return false;
	
	if(verbose){
		cout<<"passes initial tests as not missing box"<<endl;
	}
	
	// check if it's inside the box to begin with (classed as an interception at start vtx)
	if (L1.X() > B1.X() && L1.X() < B2.X() &&
		L1.Y() > B1.Y() && L1.Y() < B2.Y() &&
		L1.Z() > B1.Z() && L1.Z() < B2.Z()){
		if(verbose){
			cout<<"starts in the box"<<endl;
		}
		Hit = L1;
	}
	// check if it's inside the box at endpoint (classed as an interception at end vtx)
	if (L2.X() > B1.X() && L2.X() < B2.X() &&
		L2.Y() > B1.Y() && L2.Y() < B2.Y() &&
		L2.Z() > B1.Z() && L2.Z() < B2.Z()){
		Hit2 = L2;
		if(verbose){
			cout<<"ends in the box"<<endl;
		}
	}
	if(Hit==L1&&Hit2==L2) return true;
	
	// check for an interception in X, Y then Z.
	//if ( (GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
	//  || (GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
	//  || (GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
	//  || (GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
	//	return true;
	
	// Above seems to assume there will only be one interception!!
	// e.g. if X has an interception, there are no checks for Z interception - if it enters
	// the front face and exits the side, only the side exit will be returned. 
	// Instead, note all interception points and return the first (and second if it exists)
	std::vector<Position> interceptions;
	bool anyinterception=false;
	bool thisinterception;
	
	if(verbose){
		cout<<"finding interceptions"<<endl;
	}
	Position temphit;
	thisinterception=
	GetIntersection( L1.X()-B1.X(), L2.X()-B1.X(), L1, L2, temphit) && InBox(temphit, B1, B2, 1);
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B1.Y(), L2.Y()-B1.Y(), L1, L2, temphit) && InBox( temphit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B1.Z(), L2.Z()-B1.Z(), L1, L2, temphit) && InBox( temphit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.X()-B2.X(), L2.X()-B2.X(), L1, L2, temphit) && InBox( temphit, B1, B2, 1 );
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Y()-B2.Y(), L2.Y()-B2.Y(), L1, L2, temphit) && InBox( temphit, B1, B2, 2 );
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	thisinterception=
	GetIntersection( L1.Z()-B2.Z(), L2.Z()-B2.Z(), L1, L2, temphit) && InBox( temphit, B1, B2, 3 );
	if(thisinterception){ interceptions.push_back(temphit); anyinterception=true; }
	
	if(verbose){
		cout<<"found "<<interceptions.size()<<" interceptions";
		if(interceptions.size()>0) { cout<<" at "; interceptions.at(0).Print(); }
		if(interceptions.size()>1) { cout<<" and "; interceptions.at(1).Print(); }
		cout<<endl;
	}
	
	if(interceptions.size()>2){
		cerr<<"CheckLineBox found more than two intercepts?! They are at:"<<endl;
		for(auto&& avec : interceptions)
			cerr<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")"<<endl;
		error=true;
		return false;
	} else if(interceptions.size()==2){
		auto vec1 = interceptions.at(0);
		auto vec2 = interceptions.at(1);
		if(vec1.Z()<vec2.Z()){
			Hit=vec1;
			Hit2=vec2;
		} else {
			Hit=vec2;
			Hit2=vec1;
		}
		return true;
	} else if(interceptions.size()==1){
		if(Hit==L1){ // track starts in mrd - found intercept is exit point
			Hit2=interceptions.at(0);
			return true;
		} else {     // track starts outside mrd - found intercept is entry point
			Hit=interceptions.at(0);
			return true;
		}
	} else {
		return false;
	}
}

bool MCParticleProperties::CheckProjectedMRDHit(Position startvertex, Position stopvertex, double mrdwidth, double mrdheight, double mrdstart){


	bool projected_mrdhit = true;
	Position diffvertex = stopvertex - startvertex;

	//calculate projected x and y positions of the extended particle track from the track to the MRD
	double ProjMRDHit = (mrdstart - startvertex.Z())/(diffvertex.Z());
	double ProjXMRD = startvertex.X() + ProjMRDHit*diffvertex.X();
	double ProjYMRD = startvertex.Y() + ProjMRDHit*diffvertex.Y();

	//check if the projected x and y positions at the z-position of the first MRD layer are contained within the first MRD layers
	if (fabs(ProjXMRD) > mrdwidth) projected_mrdhit = false;
	if (fabs(ProjYMRD) > mrdheight) projected_mrdhit = false;
        if (ProjMRDHit < 0.) projected_mrdhit = false;			//particle going in wrong direction

	return projected_mrdhit;

}

//=========================================================================================
// A test to see whether a projected line within a plane intersects a circle within the plane

bool MCParticleProperties::CheckTankIntercepts( Position startvertex, Position stopvertex, bool trackstartsintank, bool trackstopsintank, Position &Hit, Position &Hit2, int verbose){
	if(verbose){
		cout<<"CheckTankIntercepts with startvertex=("<<startvertex.X()<<", "<<startvertex.Y()
			<<", "<<startvertex.Z()<<"), stopvertex=("<<stopvertex.X()<<", "
			<<stopvertex.Y()<<", "<<stopvertex.Z()<<")"<<endl;
	}
	
	// simple checks to save time
	if( 
		( (startvertex.X() > tank_radius) && (stopvertex.X() > tank_radius) ) ||
		( (startvertex.X() < -tank_radius) && (stopvertex.X() < -tank_radius) ) ||
		( ((startvertex.Y()-tank_yoffset) > tank_halfheight) && ((stopvertex.Y()-tank_yoffset) > tank_halfheight) ) ||
		( ((startvertex.Y()-tank_yoffset) < -tank_halfheight) && ((stopvertex.Y()-tank_yoffset) < -tank_halfheight) ) ||
		( (startvertex.Z() < tank_start) && (stopvertex.Z() < tank_start) ) ||
		( (startvertex.Z() > (tank_start+2.*tank_radius)) && (stopvertex.Z() > (tank_start+2.*tank_radius)) )
	){
		logmessage = "MCParticleProperties Tool::CheckTankIntercepts breaking early as start and endpoints"
					 " are both on the same side outside of the tank bounds";
		Log(logmessage,v_debug,verbosity);
		
		if(verbose){
		if( (startvertex.X() > tank_radius) && (stopvertex.X() > tank_radius) )
			cout<<"track entirely on right of tank"<<endl;
		if( (startvertex.X() < -tank_radius) && (stopvertex.X() < -tank_radius) )
			cout<<"track entirely on left of tank"<<endl;
		if( ((startvertex.Y()-tank_yoffset) > tank_halfheight) && ((stopvertex.Y()-tank_yoffset) > tank_halfheight) )
			cout<<"track entirely above tank"<<endl;
		if( ((startvertex.Y()-tank_yoffset) < -tank_halfheight) && ((stopvertex.Y()-tank_yoffset) < -tank_halfheight) )
			cout<<"track entirely below tank"<<endl;
		if( (startvertex.Z() < tank_start) && (stopvertex.Z() < tank_start) )
			cout<<"track entirely before tank"<<endl;
		if( (startvertex.Z() > (tank_start+2.*tank_radius)) && (stopvertex.Z() > (tank_start+2.*tank_radius)) )
			cout<<"track entirely after tank"<<endl;
		}
		
		return false;
	}
	
	// calculate average gradients
	double oppx = stopvertex.X() - startvertex.X();
	double adj = stopvertex.Z() - startvertex.Z();
	double avgtrackgradx = (adj!=0) ? (oppx/adj) : 1000000;
	double oppy = stopvertex.Y() - startvertex.Y();
	double avgtrackgrady = (adj!=0) ? (oppy/adj) : 1000000;
	
	// make it unit so we don't get issues with very short tracks
	Position differencevector = (stopvertex-startvertex).Unit();
	
	// first check for the track being in the z plane, as this will produce infinite gradients
	if(abs(differencevector.Z())<0.1){
		// we have a simpler case: the tank is simply a box of height tankheight
		// and width = length of a chord at the given z
		if(verbose){
			cout<<"track is in z plane, using alt method"<<endl;
		}
		double base = abs(startvertex.Z()-tank_start-tank_radius);
		double chordlen = sqrt(pow(tank_radius,2.)-pow(base,2.)); // half the chord length, strictly. 
		if(verbose){
			cout<<"half width of tank in plane z="<<startvertex.Z()<<" is "<<chordlen<<endl;
		}
		
		// simple checks to save time
		if(base>tank_radius) return false; // track misses tank in z plane
		if( ((startvertex.Y()-tank_yoffset)<-tank_halfheight) &&
		    ((stopvertex.Y() -tank_yoffset)<-tank_halfheight) ){
		    return false; // track is entirely below tank
		} else 
		if( ((startvertex.Y()-tank_yoffset)>tank_halfheight) &&
		    ((stopvertex.Y() -tank_yoffset)>tank_halfheight) ){
		    return false; // track is entirely above tank
		} else 
		if( (startvertex.X()<-chordlen) &&
		    (stopvertex.Z( )<-chordlen) ){
		    return false; // track is entirely to left of tank
		} else 
		if( (startvertex.X()>chordlen) &&
		    (stopvertex.X() >chordlen) ){
		    return false; // track is entirely to right of tank
		}
		
		if(verbose){
			cout<<"track intercepted tank somewhere"<<endl;
		}
		// second possibility: track is parallel to x-axis: then d*/dz=inf AND dx/dy=inf
		if(abs(differencevector.Y())<0.1){
			if(verbose){
				cout<<"track ran parallel to x axis"<<endl;
			}
			// even simpler: entry and exit points are just the tank walls in that z plane
			bool isrightgoing = (stopvertex.X()>startvertex.X());
			if(!trackstartsintank){
				(isrightgoing) ? Hit2.SetX(-chordlen) : Hit2.SetX(chordlen);
				Hit2.SetY(startvertex.Y()); Hit2.SetZ(startvertex.Z());
			}
			if(!trackstopsintank){
				(isrightgoing) ? Hit.SetX(chordlen) : Hit.SetX(-chordlen);
				Hit.SetY(startvertex.Y()); Hit.SetZ(startvertex.Z());
			}
			
		} else {
			
			// track intercepted the tank at least once. find where
			double trackdxdy = (stopvertex.X()-startvertex.X())/(stopvertex.Y()-startvertex.Y());
			if(verbose){
				cout<<"trackdxdy="<<trackdxdy<<endl;
			}
			// there are 4 possible intersection points (2x, 2y), of which only 1 or 2 may be relevant
			double xentry=-1., xexit=-1., yentry=-1., yexit=-1.;
			bool entryfound=false, exitfound=false;
			
			// entry point
			if(!trackstartsintank){ // if startvol is not in tank we must have had an entry point
				if(verbose){
					cout<<"track started outside tank, finding entry point"<<endl;
				}
				// find y entry - track must start outside tank y bounds to have a cap entry
				if(abs(startvertex.Y()-tank_yoffset)>tank_halfheight){
					bool isupgoing = (stopvertex.Y()>startvertex.Y());
					(isupgoing) ? yentry=-tank_halfheight : yentry=tank_halfheight;
					double distancetoyentry =  yentry - startvertex.Y(); // signed
					double xatyentry = startvertex.X()+(trackdxdy*distancetoyentry);
					if(abs(xatyentry)<chordlen){
						Hit2.SetX(xatyentry);
						Hit2.SetY(yentry+tank_yoffset);
						entryfound=true;
					}
					// else track entered via a wall
					if(verbose){
						cout<<"track started outside caps ";
						if(entryfound){
							cout<<"and entered at ("<<Hit2.X()<<", "<<(Hit2.Y()-tank_yoffset)<<")"<<endl;
						} else {
							cout<<"but entered via the wall, as xatyentry="<<xatyentry<<endl;
						}
					}
				} // else track did not start outside tank y bounds: no cap entry
				
				// find wall entry - track must start outside tank x bounds to have a wall entry
				if((!entryfound) && abs(startvertex.X())>chordlen){
					bool isrightgoing = (stopvertex.X()>startvertex.X());
					(isrightgoing) ? xentry=-chordlen : xentry=chordlen;
					double distancetoxentry =  xentry - (startvertex.X()); // signed
					double yatxentry = startvertex.Y()+(distancetoxentry/trackdxdy);
					if(abs(yatxentry)<tank_halfheight){ Hit2.SetX(xentry); Hit2.SetY(yatxentry); entryfound=true; }
					// else track entered via a cap
					if(verbose){
						cout<<"track started outside walls ";
						if(entryfound){
							cout<<"and entered at ("<<Hit2.X()<<", "<<(Hit2.Y()-tank_yoffset)<<")"<<endl;
						} else {
							cout<<"but yatxentry="<<yatxentry<<endl;
						}
					}
				} // else track did not start outside tank x bounds: no wall entry
				
				if(!entryfound) cerr<<"could not find track entry point!?"<<endl;
				else Hit2.SetZ(startvertex.Z());
				if(verbose){
					if(entryfound) cout<<"setting entry Z to "<<(Hit2.Z()-tank_start-tank_radius)<<endl;
				}
			}
			
			// repeat for exit point
			if(!trackstopsintank){
				if(verbose){
					cout<<"track ended outside tank, finding exit point"<<endl;
				}
				// find y exit - track must stop outside tank y bounds to have a cap exit
				if(abs(stopvertex.Y()-tank_yoffset)>tank_halfheight){
					bool isupgoing = (stopvertex.Y()>startvertex.Y());
					(isupgoing) ? yexit=tank_halfheight : yexit=-tank_halfheight;
					double distancetoyexit =  yexit - startvertex.Y(); // signed
					double xatyexit = startvertex.X()+(trackdxdy*distancetoyexit);
					if(abs(xatyexit)<chordlen){ Hit.SetX(xatyexit); Hit.SetY(yexit+tank_yoffset); exitfound=true; }
					// else track entered via a wall
					if(verbose){
						cout<<"track stopped outside caps ";
						if(exitfound) cout<<"and exited at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
						else cout<<"but exited via the wall, as xatyexit="<<xatyexit<<endl;
					}
				} // else track did not start outside tank y bounds: no cap exit
				
				// find wall exit - track must start outside tank x bounds to have a wall exit
				if((!exitfound) && abs(stopvertex.X())>chordlen){
					bool isrightgoing = (stopvertex.X()>startvertex.X());
					(isrightgoing) ? xexit=chordlen : xexit=-chordlen;
					double distancetoxexit =  xexit - (startvertex.X()); // signed
					double yatxexit = startvertex.Y()+(distancetoxexit/trackdxdy);
					if(abs(yatxexit)<tank_halfheight){ Hit.SetX(xexit); Hit.SetY(yatxexit); exitfound=true; }
					// else track entered via a cap
					if(verbose){
						cout<<"track stopped outside walls ";
						if(entryfound) cout<<"and exited at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
						else cout<<"but yatxexit="<<yatxexit<<endl;
					}
				} // else track did not start outside tank y bounds: no wall exit
				
				if(!exitfound) cerr<<"could not find track exit point!?"<<endl;
				else Hit.SetZ(startvertex.Z());
				if(verbose){
					if(exitfound) cout<<"setting exit Z to "<<(Hit.Z()-tank_start-tank_radius)<<endl;
				}
			}
			if( ((!trackstopsintank)&&(!exitfound)) || ((!trackstartsintank)&&(!entryfound)) ) return false;
		
		}
		
	} else {  // track was not in z plane: use old method, based on dx/dz, dy/dz
		
		// parameterize the track by it's z position, and find the 
		// z value where the track would leave the radius of the tank
		double xattankcentre = startvertex.X() - (startvertex.Z()-tank_start-tank_radius)*avgtrackgradx;
		double firstterm = -avgtrackgradx*xattankcentre;
		double thirdterm = 1+pow(avgtrackgradx,2.);
		double secondterm = (pow(tank_radius,2.)*thirdterm) - pow(xattankcentre,2.);
		if(secondterm<=0){
			if(verbose){ cout<<"Tank miss!"<<endl; }
			return false; // Doens't hit the tank.
		}
		double solution1 = (firstterm + sqrt(secondterm))/thirdterm;
		double solution2 = (firstterm - sqrt(secondterm))/thirdterm;
		double tankstartpointz, tankendpointz;
		if(stopvertex.Z() > startvertex.Z()){
			tankendpointz = solution1;    //forward going track
			tankstartpointz = solution2;  // more upstream intercept is entry point
		} else {
			tankendpointz = solution2;    // backward going track
			tankstartpointz = solution1;  // more upstream intercept is exit point
		}
		// correct for tank z offset
		tankendpointz += tank_start+tank_radius;
		tankstartpointz += tank_start+tank_radius;
		if(verbose){
			// for sanity check:
			cout<<"z'1 = "<<solution1<<", z1 = "<<(solution1+startvertex.Z())
				<<", x1 = "<<(xattankcentre+(avgtrackgradx*solution1))
				<<", r1 = "<<sqrt(pow(solution1,2.)+pow(xattankcentre+(avgtrackgradx*solution1),2.))<<endl;
			cout<<"z'2 = "<<solution2<<", z2 = "<<(solution2+startvertex.Z())
				<<", x2 = "<<(xattankcentre+(avgtrackgradx*solution2))
				<<", r2 = "<<sqrt(pow(solution2,2.)+pow(xattankcentre+(avgtrackgradx*solution2),2.))<<endl;
		}
		
		// calculate x by projecting track along parameterization
		double tankendpointx = startvertex.X() + (tankendpointz-startvertex.Z())*avgtrackgradx;
		double tankstartpointx = startvertex.X() + (tankstartpointz-startvertex.Z())*avgtrackgradx;
		
		// and same for y
		double tankendpointy = startvertex.Y() + (tankendpointz-startvertex.Z())*avgtrackgrady;
		double tankstartpointy = startvertex.Y() + (tankstartpointz-startvertex.Z())*avgtrackgrady;
		
		if(verbose){
			cout<<"tank start: "<<tank_start<<endl;
			cout<<"tank end: "<<(tank_start+2*tank_radius)<<endl;
			cout<<"tank radius: "<<tank_radius<<endl;
			Position trackdir = (stopvertex - startvertex).Unit();
			cout<<"start dir =("<<trackdir.X()<<", "<<trackdir.Y()
				<<", "<<trackdir.Z()<<")"<<endl;
			cout<<"avgtrackgradx="<<avgtrackgradx<<endl;
			cout<<"avgtrackgrady="<<avgtrackgrady<<endl;
			cout<<"xattankcentre="<<xattankcentre<<endl;
			cout<<"firstterm="<<firstterm<<endl;
			cout<<"thirdterm="<<thirdterm<<endl;
			cout<<"secondterm="<<secondterm<<endl;
			cout<<"tank intercept z solution1="<<solution1<<endl;
			cout<<"tank intercept z solution2="<<solution2<<endl<<endl;
			
			cout<<"values before cap exit check"<<endl;
			cout<<"tankstartpointz="<<tankstartpointz<<endl;
			cout<<"tankstartpointx="<<tankstartpointx<<endl;
			cout<<"tankstartpointy="<<tankstartpointy<<endl;
			
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
			
			cout<<"TMath::Abs(startvertex.Y()-tank_yoffset)="
				<<TMath::Abs(startvertex.Y()-tank_yoffset)
				<<"(tank_halfheight)="<<(tank_halfheight)<<endl;
		}
		
		// now check if the particle would have exited through one of the caps before reaching this point
		// if projected y value at the radial exit is outside the tank, and the track started inside the tank,
		// then the tank must have left prior to this point via one of the caps
		if( TMath::Abs(tankendpointy-tank_yoffset)>(tank_halfheight) && 
			TMath::Abs(startvertex.Y()-tank_yoffset)<(tank_halfheight)){
			// ^ second condition due to slight (mm) differences in tank y height in WCSim.
			// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
			if(stopvertex.Y()>startvertex.Y()){
				tankendpointy = tank_halfheight+tank_yoffset;  // by definition of leaving through cap
			} else {
				tankendpointy = -tank_halfheight+tank_yoffset;
			}
			tankendpointz = 
			startvertex.Z() + (tankendpointy-startvertex.Y())/avgtrackgrady;
			tankendpointx = 
			startvertex.X() + (tankendpointz-startvertex.Z())*avgtrackgradx;
		} // else this trajectory exited the tank by a side point; existing value is valid
		
		// repeat for the tankstartpoint
		if( TMath::Abs(tankstartpointy-tank_yoffset)>(tank_halfheight) && 
			TMath::Abs(startvertex.Y()-tank_yoffset)<(tank_halfheight)){
			// ^ second condition due to slight (mm) differences in tank y height in WCSim.
			// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
			if(stopvertex.Y()>startvertex.Y()){
				tankstartpointy = tank_halfheight+tank_yoffset;	// by definition of leaving through cap
			} else {
				tankstartpointy = -tank_halfheight+tank_yoffset;
			}
			tankstartpointz = 
			startvertex.Z() + (tankstartpointy-startvertex.Y())/avgtrackgrady;
			tankstartpointx = 
			startvertex.X() + (tankstartpointz-startvertex.Z())*avgtrackgradx;
		} // else this trajectory entered the tank by a side point; existing value is valid
		
		if(verbose){
			cout<<"values after cap exit check"<<endl;
			cout<<"tankstartpointz="<<tankstartpointz<<endl;
			cout<<"tankstartpointx="<<tankstartpointx<<endl;
			cout<<"tankstartpointy="<<tankstartpointy<<endl;
			
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
		}
		
		// return only the appropriate vertices
		if(!trackstopsintank) Hit=Position(tankendpointx,tankendpointy,tankendpointz);
		if(!trackstartsintank) Hit2=Position(tankstartpointx,tankstartpointy,tankstartpointz);
	}
	
	return true;
}

//=========================================================================================
//=========================================================================================
// Find the particle would have left the tank, had it done so

// XXX remove me, hacked these to static so we can make the following function static and usable by other tools
double MCParticleProperties::tank_radius=MRDSpecs::tank_radius;
double MCParticleProperties::tank_start=MRDSpecs::tank_start;
double MCParticleProperties::tank_yoffset=MRDSpecs::tank_yoffset;
double MCParticleProperties::tank_halfheight=MRDSpecs::tank_halfheight;

bool MCParticleProperties::ProjectTankIntercepts(Position startvertex, Position stopvertex, Position &Hit, int verbose){
	// this is based on CheckTankIntercepts which estimates the tank entry and exit points,
	// but only if the track did indeed enter or exit the tank
	// ProjectTankIntercepts cuts out all the checks, and
	// 1. assumes the particle started in the tank
	// 2. projects an exit point regardless of whether it is beyond the track end
	
	if(verbose){
		cout<<"ProjectTankIntercepts with startvertex=("<<startvertex.X()<<", "<<startvertex.Y()
			<<", "<<startvertex.Z()<<"), stopvertex=("<<stopvertex.X()<<", "
			<<stopvertex.Y()<<", "<<stopvertex.Z()<<")"<<endl;
	}
	
	// calculate average gradients
	double oppx = stopvertex.X() - startvertex.X();
	double adj = stopvertex.Z() - startvertex.Z();
	double avgtrackgradx = (adj!=0) ? (oppx/adj) : 1000000;
	double oppy = stopvertex.Y() - startvertex.Y();
	double avgtrackgrady = (adj!=0) ? (oppy/adj) : 1000000;
	
	// make it unit so we don't get rounding issues with very short tracks
	Position differencevector = (stopvertex-startvertex).Unit();
	
	// first check for the track being in the z plane, as this will produce infinite gradients
	if(abs(differencevector.Z())<0.1){
		// we have a simpler case: the tank is simply a box of height tankheight
		// and width = length of a chord at the given z
		if(verbose){
			cout<<"track is in z plane, using alt method"<<endl;
		}
		double base = abs(startvertex.Z()-tank_start-tank_radius);
		double chordlen = sqrt(pow(tank_radius,2.)-pow(base,2.)); // half the chord length, strictly. 
		if(verbose){
			cout<<"half width of tank in plane z="<<startvertex.Z()<<" is "<<chordlen<<endl;
		}
		// second possibility: track is parallel to x-axis: then d*/dz=inf AND dx/dy=inf
		if(abs(differencevector.Y())<0.1){
			if(verbose){
				cout<<"track ran parallel to x axis"<<endl;
			}
			// even simpler: entry and exit points are just the tank walls in that z plane
			bool isrightgoing = (stopvertex.X()>startvertex.X());
			(isrightgoing) ? Hit.SetX(chordlen) : Hit.SetX(-chordlen); // exit point is tank extent in that plane
			Hit.SetY(startvertex.Y()); Hit.SetZ(startvertex.Z());      // parallel to y and z axes: no change
			
		} else {
			// track is in z-plane, but not along the X axis: line through a box scenario
			double trackdxdy = (stopvertex.X()-startvertex.X())/(stopvertex.Y()-startvertex.Y());
			if(verbose){
				cout<<"trackdxdy="<<trackdxdy<<endl;
			}
			// there are 4 possible intersection points (2x, 2y), of which only 1 or 2 may be relevant
			double xentry=-1., xexit=-1., yentry=-1., yexit=-1.;
			bool entryfound=true, exitfound=false;
			
			// find exit point
			if(verbose){
				cout<<"finding projected exit point"<<endl;
			}
			// find y exit - calculate points for both barrel and cap exit
			bool isupgoing = (stopvertex.Y()>startvertex.Y());
			(isupgoing) ? yexit=tank_halfheight : yexit=-tank_halfheight;
			double distancetoyexit =  yexit - startvertex.Y(); // signed
			double xatyexit = startvertex.X()+(trackdxdy*distancetoyexit);
			if(abs(xatyexit)<chordlen){ Hit.SetX(xatyexit); Hit.SetY(yexit+tank_yoffset); exitfound=true; }
			// else track entered via a wall
			if(verbose){
				if(exitfound) cout<<"track exited via cap at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
				else cout<<"track exited via the wall, as xatyexit="<<xatyexit<<endl;
			}
			
			// find wall exit
			if(!exitfound){
				bool isrightgoing = (stopvertex.X()>startvertex.X());
				(isrightgoing) ? xexit=chordlen : xexit=-chordlen;
				double distancetoxexit =  xexit - (startvertex.X()); // signed
				double yatxexit = startvertex.Y()+(distancetoxexit/trackdxdy);
				if(abs(yatxexit)<tank_halfheight){ Hit.SetX(xexit); Hit.SetY(yatxexit); exitfound=true; }
				// else track entered via a cap
				if(verbose){
					if(entryfound) cout<<"track exited through the barrel at ("<<Hit.X()<<", "<<(Hit.Y()-tank_yoffset)<<")"<<endl;
					else cerr<<"but yatxexit="<<yatxexit
							 <<">tank_halfheight suggests we should have already found a cap exit!"<<endl;
				}
			} // else we already found a cap exit
			
			if(!exitfound) cerr<<"MCParticleProperties::ProjectTankExit could not find track exit point!?"<<endl;
			else Hit.SetZ(startvertex.Z());
			if(verbose){
				if(exitfound) cout<<"setting exit Z to "<<(Hit.Z()-tank_start-tank_radius)<<endl;
			}
		} // end else track did not run along x-axis
		
	} else {  // track was not in z plane: use old method, based on dx/dz, dy/dz
		
		// parameterize the track by it's z position, and find the 
		// z value where the track would leave the radius of the tank
		double xattankcentre = startvertex.X() - (startvertex.Z()-tank_start-tank_radius)*avgtrackgradx;
		double firstterm = -avgtrackgradx*xattankcentre;
		double thirdterm = 1+pow(avgtrackgradx,2.);
		double secondterm = (pow(tank_radius,2.)*thirdterm) - pow(xattankcentre,2.);
		if(secondterm<=0){
			if(verbose){ cout<<"MCParticleProperties::ProjectTankExit says Tank miss?!"<<endl; }
			return false; // Doens't hit the tank.
		}
		double solution1 = (firstterm + sqrt(secondterm))/thirdterm;
		double solution2 = (firstterm - sqrt(secondterm))/thirdterm;
		double tankendpointz;
		if(stopvertex.Z() > startvertex.Z()){
			tankendpointz = solution1;    //forward going track
		} else {
			tankendpointz = solution2;    // backward going track
		}
		// correct for tank z offset
		tankendpointz += tank_start+tank_radius;
		if(verbose){
			// for sanity check:
			cout<<"z'1 = "<<solution1<<", z1 = "<<(solution1+startvertex.Z())
				<<", x1 = "<<(xattankcentre+(avgtrackgradx*solution1))
				<<", r1 = "<<sqrt(pow(solution1,2.)+pow(xattankcentre+(avgtrackgradx*solution1),2.))<<endl;
			cout<<"z'2 = "<<solution2<<", z2 = "<<(solution2+startvertex.Z())
				<<", x2 = "<<(xattankcentre+(avgtrackgradx*solution2))
				<<", r2 = "<<sqrt(pow(solution2,2.)+pow(xattankcentre+(avgtrackgradx*solution2),2.))<<endl;
		}
		
		// calculate x by projecting track along parameterization
		double tankendpointx = startvertex.X() + (tankendpointz-startvertex.Z())*avgtrackgradx;
		
		// and same for y
		double tankendpointy = startvertex.Y() + (tankendpointz-startvertex.Z())*avgtrackgrady;
		
		if(verbose){
			cout<<"tank start: "<<tank_start<<endl;
			cout<<"tank end: "<<(tank_start+2*tank_radius)<<endl;
			cout<<"tank radius: "<<tank_radius<<endl;
			Position trackdir = (stopvertex - startvertex).Unit();
			cout<<"start dir =("<<trackdir.X()<<", "<<trackdir.Y()
				<<", "<<trackdir.Z()<<")"<<endl;
			cout<<"avgtrackgradx="<<avgtrackgradx<<endl;
			cout<<"avgtrackgrady="<<avgtrackgrady<<endl;
			cout<<"xattankcentre="<<xattankcentre<<endl;
			cout<<"firstterm="<<firstterm<<endl;
			cout<<"thirdterm="<<thirdterm<<endl;
			cout<<"secondterm="<<secondterm<<endl;
			cout<<"tank intercept z solution1="<<solution1<<endl;
			cout<<"tank intercept z solution2="<<solution2<<endl<<endl;
			
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
			
			cout<<"TMath::Abs(startvertex.Y()-tank_yoffset)="
				<<TMath::Abs(startvertex.Y()-tank_yoffset)
				<<"(tank_halfheight)="<<(tank_halfheight)<<endl;
		}
		
		// now check if the particle would have exited through one of the caps before reaching this point
		// if projected y value at the radial exit is outside the tank, and the track started inside the tank,
		// then the tank must have left prior to this point via one of the caps
		if( TMath::Abs(tankendpointy-tank_yoffset)>(tank_halfheight) && 
			TMath::Abs(startvertex.Y()-tank_yoffset)<(tank_halfheight)){
			// ^ second condition due to slight (mm) differences in tank y height in WCSim.
			// this trajectory exits through the cap. Need to recalculate x, z exiting points...!
			if(stopvertex.Y()>startvertex.Y()){
				tankendpointy = tank_halfheight+tank_yoffset; // by definition of leaving through cap
			} else {
				tankendpointy = -tank_halfheight+tank_yoffset;
			}
			tankendpointz = 
			startvertex.Z() + (tankendpointy-startvertex.Y())/avgtrackgrady;
			tankendpointx = 
			startvertex.X() + (tankendpointz-startvertex.Z())*avgtrackgradx;
		} // else this trajectory exited the tank by a side point; existing value is valid
		
		if(verbose){
			cout<<"values after cap exit check"<<endl;
			cout<<"tankendpointz="<<tankendpointz<<endl;
			cout<<"tankendpointx="<<tankendpointx<<endl;
			cout<<"tankendpointy="<<tankendpointy<<endl;
		}
		Hit.SetX(tankendpointx);
		Hit.SetY(tankendpointy);
		Hit.SetZ(tankendpointz);
	}
	
	return true;
}

//=========================================================================================



std::string MCParticleProperties::PdgToString(int code){
	if(pdgcodetoname.size()==0) GeneratePdgMap();
	if(pdgcodetoname.count(code)!=0){
		return pdgcodetoname.at(code);
	} else {
		cerr<<"unknown pdg code "<<code<<endl;
		return std::to_string(code);
	}
}

std::map<int,std::string>* MCParticleProperties::GeneratePdgMap(){
	if(pdgcodetoname.size()!=0) return &pdgcodetoname;
	pdgcodetoname.emplace(2212,"Proton");
	pdgcodetoname.emplace(-2212,"Anti Proton");
	pdgcodetoname.emplace(11,"Electron");
	pdgcodetoname.emplace(-11,"Positron");
	pdgcodetoname.emplace(12,"Electron Neutrino");
	pdgcodetoname.emplace(-12,"Anti Electron Neutrino");
	pdgcodetoname.emplace(22,"Gamma");
	pdgcodetoname.emplace(2112,"Neutron");
	pdgcodetoname.emplace(-2112,"Anti Neutron");
	pdgcodetoname.emplace(-13,"Muon+");
	pdgcodetoname.emplace(13,"Muon-");
	pdgcodetoname.emplace(130,"Kaonlong");
	pdgcodetoname.emplace(211,"Pion+");
	pdgcodetoname.emplace(-211,"Pion-");
	pdgcodetoname.emplace(321,"Kaon+");
	pdgcodetoname.emplace(-321,"Kaon-");
	pdgcodetoname.emplace(3122,"Lambda");
	pdgcodetoname.emplace(-3122,"Antilambda");
	pdgcodetoname.emplace(310,"Kaonshort");
	pdgcodetoname.emplace(3112,"Sigma-");
	pdgcodetoname.emplace(3222,"Sigma+");
	pdgcodetoname.emplace(3212,"Sigma0");
	pdgcodetoname.emplace(111,"Pion0");
	pdgcodetoname.emplace(311,"Kaon0");
	pdgcodetoname.emplace(-311,"Antikaon0");
	pdgcodetoname.emplace(14,"Muon Neutrino");
	pdgcodetoname.emplace(-14,"Anti Muon Neutrino");
	pdgcodetoname.emplace(-3222,"Anti Sigma-");
	pdgcodetoname.emplace(-3212,"Anti Sigma0");
	pdgcodetoname.emplace(-3112,"Anti Sigma+");
	pdgcodetoname.emplace(3322,"Xsi0");
	pdgcodetoname.emplace(-3322,"Anti Xsi0");
	pdgcodetoname.emplace(3312,"Xsi-");
	pdgcodetoname.emplace(-3312,"Xsi+");
	pdgcodetoname.emplace(3334,"Omega-");
	pdgcodetoname.emplace(-3334,"Omega+");
	pdgcodetoname.emplace(-15,"Tau+");
	pdgcodetoname.emplace(15,"Tau-");
	pdgcodetoname.emplace(100,"OpticalPhoton");
	pdgcodetoname.emplace(3328,"Alpha");
	pdgcodetoname.emplace(3329,"Deuteron");
	pdgcodetoname.emplace(3330,"Triton");
	pdgcodetoname.emplace(3351,"Li7");
	pdgcodetoname.emplace(3331,"C10");
	pdgcodetoname.emplace(3345,"B11");
	pdgcodetoname.emplace(3332,"C12");
	pdgcodetoname.emplace(3350,"C13");
	pdgcodetoname.emplace(3349,"N13");
	pdgcodetoname.emplace(3340,"N14");
	pdgcodetoname.emplace(3333,"N15");
	pdgcodetoname.emplace(3334,"N16");
	pdgcodetoname.emplace(3335,"O16");
	pdgcodetoname.emplace(3346,"Al27");
	pdgcodetoname.emplace(3341,"Fe54");
	pdgcodetoname.emplace(3348,"Mn54");
	pdgcodetoname.emplace(3342,"Mn55");
	pdgcodetoname.emplace(3352,"Mn56");
	pdgcodetoname.emplace(3343,"Fe56");
	pdgcodetoname.emplace(3344,"Fe57");
	pdgcodetoname.emplace(3347,"Fe58");
	pdgcodetoname.emplace(3353,"Eu154");
	pdgcodetoname.emplace(3336,"Gd158");
	pdgcodetoname.emplace(3337,"Gd156");
	pdgcodetoname.emplace(3338,"Gd157");
	pdgcodetoname.emplace(3339,"Gd155");
	return &pdgcodetoname;
}



double MCParticleProperties::PdgToMass(int code){
	if(pdgcodetomass.size()==0) GeneratePdgMassMap();
	if(pdgcodetomass.count(code)!=0){
		return pdgcodetomass.at(code);
	} else {
		cerr<<"unknown pdg code "<<code<<endl;
		return double(code);
	}
}

std::map<int,double>* MCParticleProperties::GeneratePdgMassMap(){
	//all masses in MeV
	if(pdgcodetomass.size()!=0) return &pdgcodetomass;
	pdgcodetomass.emplace(2212,938.272);
	pdgcodetomass.emplace(-2212,938.272);
	pdgcodetomass.emplace(11,0.511);
	pdgcodetomass.emplace(-11,0.511);
	pdgcodetomass.emplace(12,0.000);
	pdgcodetomass.emplace(-12,0.000);
	pdgcodetomass.emplace(22,0.000);
	pdgcodetomass.emplace(2112,939.565);
	pdgcodetomass.emplace(-2112,939.485);
	pdgcodetomass.emplace(-13,105.658);
	pdgcodetomass.emplace(13,105.658);
	pdgcodetomass.emplace(211,139.571);
	pdgcodetomass.emplace(-211,139.571);
	pdgcodetomass.emplace(321,493.666);
	pdgcodetomass.emplace(-321,493.666);
	pdgcodetomass.emplace(3122,1115.683);
	pdgcodetomass.emplace(-3122,1115.683);
	pdgcodetomass.emplace(3112,1197.449);
	pdgcodetomass.emplace(3222,1189.37);
	pdgcodetomass.emplace(3212,1192.642);
	pdgcodetomass.emplace(111,134.977);
	pdgcodetomass.emplace(311,497.611);
	pdgcodetomass.emplace(-311,497.611);
	pdgcodetomass.emplace(14,0.000);
	pdgcodetomass.emplace(-14,0.000);
	pdgcodetomass.emplace(-3222,1189.37);
	pdgcodetomass.emplace(-3212,1192.642);
	pdgcodetomass.emplace(-3112,1197.449);
	pdgcodetomass.emplace(3322,1314.86);
	pdgcodetomass.emplace(-3322,1314.86);
	pdgcodetomass.emplace(3312,1321.71);
	pdgcodetomass.emplace(-3312,1321.71);
	pdgcodetomass.emplace(3334,1672.45);
	pdgcodetomass.emplace(-3334,1672.45);
	pdgcodetomass.emplace(-15,1776.86);
	pdgcodetomass.emplace(15,1776.86);
	pdgcodetomass.emplace(100,0.000);
	//masses of nuclei not specified since they won't be able to generate Cherenkov light
	return &pdgcodetomass;
}

