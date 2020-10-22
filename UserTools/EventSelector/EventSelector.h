//This tool reads information from both the Monte Carlo and
//Event reconstruction and flags each event with the
//chosen flags in the EventSelector config file.

#ifndef EventSelector_H
#define EventSelector_H

#include <string>
#include <iostream>
#include <bitset>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "Tool.h"
// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "ANNIEGeometry.h"
#include "TMath.h"

class EventSelector: public Tool {


 public:

  EventSelector();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  typedef enum EventFlags {
   kFlagNone         = 0x00, //0
   kFlagMCFV         = 0x01, //1
   kFlagMCPMTVol     = 0x02, //2
   kFlagMCMRD        = 0x04, //4
   kFlagMCPiK        = 0x08, //8
   kFlagRecoMRD      = 0x10, //16
   kFlagPromptTrig   = 0x20, //32
   kFlagNHit         = 0x40, //64
   kFlagRecoFV       = 0x80, //128
   kFlagRecoPMTVol   = 0x100, //256
   kFlagMCIsMuon       = 0x200, //512
   kFlagMCIsElectron   = 0x400, //1024
   kFlagMCIsSingleRing = 0x800, //2048
   kFlagMCIsMultiRing  = 0x1000, //4096
   kFlagMCProjectedMRDHit = 0x2000, //8192
   kFlagMCEnergyCut   = 0x4000, //16384
   kFlagPMTMRDCoinc   = 0x8000, //32768
   kFlagNoVeto        = 0x10000, //65536
   kFlagVeto        = 0x20000, //131072
   kFlagMCDSNBlike = 0x40000 //262144
  } EventFlags_t;

 private:
 	
  /// Clear reconstruction info.
  void Reset();

  /// \brief Event selection by MRD reconstructed information
  ///
  /// Loop over all the MRC tracks. Find the track with the longest track
  /// length. If the longest track is stoped inside the MRD, the event is 
  /// selected
  bool EventSelectionByMRDReco();

  /// \brief Event selection by trigger number
  ///
  /// The selection is based on the trigger number for the event
  /// in the store "ANNIEEvent".  Events are selected if they have
  /// an MCTriggernum == 0 (i.e. they are a prompt trigger) 
  bool PromptTriggerCheck();

  /// \brief Event selection by number of digits hit
  ///
  /// Require a minimum amount of digits to be hit in the event.
  /// If the criteria is not met, the event is flagged.
  bool NHitCountCheck(int NHitCount);

  /// \brief Event selection by fidicual volume
  ///
  /// The selection is based on the muon interaction point. 
  /// If isMC is true, checks the Event using muon truth info. 
  /// If False, the reconstructed vertex is used. 
  bool EventSelectionByFV(bool isMC);

  /// \brief Event selection by PMT Volume 
  ///
  /// The selection is based on the muon interaction point. 
  /// If isMC is true, checks the Event using muon truth info. 
  /// If False, the reconstructed vertex is used. 
 bool EventSelectionByPMTVol(bool isMC);

  /// \brief Event selection by Muon MRD stop position
  /////
  /// The selection is based on the true vertex stop position from MC. 
  /// If the true muon vertex stops inside the MRD, the event 
  /// is selected. 
  bool EventSelectionByMCTruthMRD();

  /// \brief Event selection by Pion Kaon count
  /////
  /// This event selection criteria requires that no pions or 
  /// kaons are parent particles in the event.  This will help
  /// Select the CC0Pi events when testing reconstruction.
  bool EventSelectionNoPiK();

  /// \brief Event selection by requiring the primary to be a certain particle
  ////
  /// This event selection criteria selects events with a 
  /// certain primary particle (to be selected by pdg number)
  bool ParticleCheck(int pdg_number);

  /// \brief Event selection by Energy cut
  /////
  /// This event selection criteria requires that the energy 
  /// of the primary electron/muon is between Emin and Emax [MeV]
  /// 
  bool EnergyCutCheck(double Emin, double Emax);

  /// \brief Event selection by Single Rings
  ////
  /// This event selection criteria requires events to have only
  /// one ring. Number of rings are counted by looping through primary
  /// particles and selecting charged particles above Cherenkov threshold
  bool EventSelectionByMCSingleRing();

  /// \brief Event selection by Multi Rings
  ////
  /// This event selection criteria requires events to have multiple
  /// rings. Number of rings are counted by looping through primary
  /// particles and selecting charged particles above Cherenkov threshold
  bool EventSelectionByMCMultiRing();

  /// \brief Event selection by Projected MRD Hits
  ////
  /// This event selection criteria requires events to have a projected
  /// hit in the MRD. This event selection includes events that do not actually
  /// enter the MRD but whose projected trajectory from the tank intersects with
  /// the MRD.
  bool EventSelectionByMCProjectedMRDHit();

  /// \brief Event selection by DSNB-like
  ////
  ///This event selection criteria requires a positron+neutron coincidence or a 
  ///gamma+neutron coincidence that might be misinterpreted as an IBD signature
  bool DSNBCheck();


  /// \brief Event selection by PMT/MRD time coincidence
  ////
  /// This event selection criteria requires clustered events in tank & MRD
  /// to have coincidicent time activity and therefore correspond to a single
  /// event.
  bool EventSelectionByPMTMRDCoinc();

  /// \brief Event selection by rejecting veto hits
  ////
  /// This event selection criteria requires that no veto paddles
  /// of the Front Muon Veto fired during the event
  bool EventSelectionByVetoCut();
  
  /// \brief MC entry number
  uint64_t fMCEventNum;
  
  /// \brief trigger number
  uint16_t fMCTriggernum;
  
  /// \brief ANNIE event number
  uint32_t fEventNumber;

  // \brief Event Status bitwords
  int fEventApplied; //Integer indicates what event cleaning flags were checked for the event
  int fEventFlagged; //Integer indicates what evt. cleaning flags the event was flagged with
  
  Geometry *fGeometry = nullptr;    ///< ANNIE Geometry
  RecoVertex* fMuonStartVertex = nullptr; 	 ///< true muon start vertex
  RecoVertex* fMuonStopVertex = nullptr; 	 ///< true muon stop vertex
  std::vector<RecoDigit>* fDigitList;		///< Reconstructed Hits including both LAPPD hits and PMT hits
  RecoVertex* fRecoVertex = nullptr; 	 ///< Reconstructed Vertex 
  std::map<double,std::vector<Hit>>* m_all_clusters;   ///< clustered PMT hits
  std::map<double,std::vector<MCHit>>* m_all_clusters_MC;   ///< clustered PMT hits
  std::vector<std::vector<int>> MrdTimeClusters;      ///< clustered MRD hits
  std::vector<double> MrdDigitTimes;          ///< clustered MRD times
  std::vector<unsigned long> MrdDigitChankeys;          ///< clustered MRD chankeys
  std::map<unsigned long,std::vector<MCHit>>* TDCData_MC;	///< MRD hits (MC)
  std::map<unsigned long,std::vector<Hit>>* TDCData;	///< MRD hits (data)
  std::vector<double> *vec_pmtclusters_charge = nullptr;
  std::vector<double> *vec_pmtclusters_time = nullptr;
  std::vector<double> *vec_mrdclusters_time = nullptr;

  //verbosity initialization
  int verbosity=1;
  
  std::string fInputfile;
  bool fMRDRecoCut = false;
  bool fMCFVCut = false;
  bool fMCPMTVolCut = false;
  bool fMCMRDCut = false;
  bool fMCPiKCut = false;
  bool fMCIsMuonCut = false;
  bool fMCIsElectronCut = false;
  bool fMCIsSingleRingCut = false;
  bool fMCIsMultiRingCut = false;
  bool fMCProjectedMRDHit = false;
  bool fMCEnergyCut = false;
  bool fNHitCut = true;
  int  fNHitmin = 4;
  double Emin = 0.;
  double Emax = 10000.;
  bool fMCDSNBlike = false;
  bool fRecoPMTVolCut = false;
  bool fRecoFVCut = false;
  bool fPMTMRDCoincCut = false;
  double fPMTMRDOffset = 745;
  bool fPromptTrigOnly = true;
  bool fNoVetoCut = false;
  bool fVetoCut = false;
  bool fEventCutStatus;
  bool fIsMC; 

  
  bool fSaveStatusToStore = true;
  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  
};


#endif
