#ifndef EnergyExtractor_H
#define EnergyExtractor_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class EnergyExtractor
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class EnergyExtractor: public Tool {


 public:

  EnergyExtractor(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  void FindTrueEnergyFromMC(); ///< Find true energies of particles from MCParticles object
  double GetCherenkovThresholdE(int pdg_code); ///< Get Cherenkov threshold for particles

 private:

  //User configuration variables
  int verbosity;
  bool save_neutrino = 0;
  bool save_elec = 0;
  bool save_gamma = 0;
  bool save_muon = 0;
  bool save_pion = 0;
  bool save_kaon = 0;
  bool save_neutron;
  bool save_visible = 0;

  //Variables stored in the ANNIEEvent or CStore BoostStores
  double TrueEnergy;
  ofstream outfile_EnergyNeutrino, outfile_EnergyPion, outfile_EnergyKaon, outfile_EnergyMuon, outfile_EnergyElectron, outfile_EnergyGamma, outfile_NeutronNumber ,outfile_VisibleEnergy;
  std::map<int,double> pdgcodetomass;
  string cnn_outpath;
  std::vector<MCParticle>* fMCParticles=nullptr;
  std::map<unsigned long, std::vector<MCHit>>* MCHits = nullptr;
  std::map<unsigned long, std::vector<MCLAPPDHit>>* MCLAPPDHits = nullptr;
  Geometry *geom = nullptr;
  
  //Refractive index
  double n=1.33;

  //Energy map related variables
  std::map<string, vector<double>> Energystore;
  std::vector<double> filler;
  int pi0count = 0;
  int pipcount = 0;
  int pimcount = 0;
  int K0count = 0;
  int Kpcount = 0;
  int Kmcount = 0;
  int ecount = 0;
  int epcount = 0;
  int mcount = 0;
  int mpcount = 0;
  int ncount = 0;
  int gammacount = 0;

  //Verbosity variables
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  int vv_debug=4;


};


#endif
