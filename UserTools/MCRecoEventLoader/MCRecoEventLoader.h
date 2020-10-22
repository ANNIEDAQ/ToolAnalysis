#ifndef MCRecoEventLoader_H
#define MCRecoEventLoader_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "ANNIEGeometry.h"
#include "Detector.h"
#include "TMath.h"

class MCRecoEventLoader: public Tool {


 public:

  MCRecoEventLoader();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  // Config file inputs
  int verbosity=1;
  bool fGetPiKInfo;
  bool fGetNRings;
  bool fDoParticleSelection;
  int fParticleID;
  double xshift;
  double yshift;
  double zshift;

  RecoVertex* fMuonStartVertex = nullptr; 	 ///< true muon start vertex
  RecoVertex* fMuonStopVertex = nullptr; 	 ///< true muon stop vertex
  std::vector<MCParticle>* fMCParticles=nullptr;  ///< truth tracks
  double TrueMuonEnergy = -9999.;
  double WaterTrackLength = -9999.;
  double MRDTrackLength = -9999.;
  bool projectedmrdhit = false;
  std::map<int,double> pdgcodetomass;     ///< map that links pdg code to particle mass in MeV, provided by MCParticleProperties
  std::map<int,std::string> pdgcodetoname;
  std::map<int,double> pdgcodetocherenkov;

  double n = 1.333;   //refractive index of water, should probably rather be defined somewhere else


  /// \brief Find true neutrino vertex
  ///
  /// Loop over all MC particles and find the particle with highest energy. 
  /// This particle is the primary muon. The muon start position, time and 
  /// the muon direction are used to initise the true neutrino vertex 
  void FindTrueVertexFromMC();
  
  /// \brief Find PionKaon Count 
  ///
  /// Loop over all MC particles and find any particles with PDG codes
  /// Consistent with Pions or Kaons of any charges. Racks up a count
  /// of the number of each type of particle
  /// In addition: Loop over MC particles and count the number of rings that should be produced
  /// by those particles. The particle needs to be above Cherenkov threshold to
  /// produce a ring. Neutrally charged particles like the Pi0 produce 2 rings
 	
  void FindPionKaonCountFromMC();

  /// \brief GetCherenkovThresholdE
  /// Get Cherenkov threshold energy for a given particle PDG number code

  double GetCherenkovThresholdE(int pdg_code);
  
  /// \brief Save true neutrino vertex
  ///
  /// Push true muon vertex to "RecoVertex"
  /// \param[in] bool savetodisk: save object to disk if savetodisk=true
  void PushTrueVertex(bool savetodisk);
  
  /// \brief Save true neutrino vertex
  ///
  /// Push true muon stop vertex to "RecoVertex"
  /// \param[in] bool savetodisk: save object to disk if savetodisk=true
  void PushTrueStopVertex(bool savetodisk);

  /// \brief Push muon track lengths to RecoEvent Store
  void PushTrueMuonEnergy(double MuE);
  void PushTrueWaterTrackLength(double WaterT);
  void PushTrueMRDTrackLength(double MRDT);

  /// \brief Push projected particle MRD hit boolean to RecoEvent store
  void PushProjectedMrdHit(bool projectedmrdhit);

  /// \brief Push IBD/IBD-like true information to RecoEvent store
  void PushIBDInfo();

  /// \brief Reset initialized classes
  ///
  /// Clear True Vertices 
  void Reset();

  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;

 private:





};


#endif
