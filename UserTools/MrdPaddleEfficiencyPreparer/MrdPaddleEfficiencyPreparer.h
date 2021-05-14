#ifndef MrdPaddleEfficiencyPreparer_H
#define MrdPaddleEfficiencyPreparer_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "Tool.h"
#include "Geometry.h"
#include "Paddle.h"
#include "Detector.h"
#include "Hit.h"
#include "MRDSubEventClass.hh"      // a class for defining subevents
#include "MRDTrackClass.hh"         // a class for defining MRD tracks

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

/**
 * \class MrdPaddleEfficiencyPreparer
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M.Nieslony $
* $Date: 2020/01/07 18:16:00 $
* Contact: mnieslon@uni-mainz.de
*/

class MrdPaddleEfficiencyPreparer: public Tool {

 public:

  MrdPaddleEfficiencyPreparer(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  bool FindPaddleIntersection(Position startpos, Position endpos, double &x, double &y, double z);
  bool FindPaddleChankey(double x, double y, int layer, unsigned long &chankey);


 private:

 	int verbosity;
 	std::vector<BoostStore>* theMrdTracks;                        // the reconstructed tracks
 	std::string MRDTriggertype;
	std::string outputfile; 
	bool usetruetrack;
	Geometry *geom = nullptr;
	bool isData;

	std::map<unsigned long,int> channelkey_to_mrdpmtid;
	std::map<int,unsigned long> mrdpmtid_to_channelkey;

 	int numsubevs;
 	int numtracksinev;
 	int EventNumber;
 	Position StartVertex;
 	Position StopVertex;
 	std::vector<int> PMTsHit;
	int numpmtshit;
 	int MrdTrackID;
 	double HtrackFitChi2;
 	double VtrackFitChi2;
 	std::vector<int> LayersHit;
 	int numlayershit;
	double tracklength;
	int NPMTsHit = 0;
	int NLayersHit = 0;
	double StartVertexX = 0.;
	double StartVertexY = 0.;
	double StartVertexZ = 0.;
	double StopVertexX = 0.;
	double StopVertexY = 0.;
	double StopVertexZ = 0.;

 	TFile *hist_file = nullptr;
	TH1F *hist_chankey = nullptr;
	TTree *tree_trackfit = nullptr;

 	std::map<int,std::map<unsigned long,TH1D*>> observed_MRDHits;
 	std::map<int,std::map<unsigned long,TH1D*>> expected_MRDHits;
 	std::map<int,double> zLayers;
 	std::map<int,int> orientationLayers;
	std::map<int,std::vector<unsigned long>> channelsLayers;
	std::map<unsigned long, int> map_chkey_half;

        double extents[11] = {1.318,1.146,1.318,1.299,1.318,1.318,1.318,1.521,1.318,1.521,1.318};
 	int channels_start[11] = {26,52,82,108,142,168,194,220,250,276,306};
	std::map<int,std::map<unsigned long,TH1D*>> observed_MRDHits_layer;
 	std::map<int,std::map<unsigned long,TH1D*>> expected_MRDHits_layer;

 	int v_error = 0;
 	int v_warning = 1;
 	int v_message = 2;
 	int v_debug = 3;


 };

 #endif
