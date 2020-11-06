/*
 * LAPPDDisplay.cpp
 *
 *  Created on: April 25, 2019
 *      Author: stenderm
 */

#include "LAPPDDisplay.h"
#include <thread>

/**
 * Constructor LAPPDDisplay: Initialises the class with a TApplication, the canvases and the output file.
 * @param filePath           Path and name of the output root file.
 */
LAPPDDisplay::LAPPDDisplay(std::string filePath, int confignumber):_LAPPD_sim_app(nullptr),_LAPPD_MC_all_canvas(nullptr),_LAPPD_MC_canvas(nullptr),
_LAPPD_MC_time_canvas(nullptr),_LAPPD_all_waveforms_canvas(nullptr),_LAPPD_waveform_canvas(nullptr),_all_hits(nullptr),_output_file(nullptr),
_output_file_name(filePath),_config_number(confignumber)
{
  //TApplication
	int myargc = 0;
	char *myargv[] =
	{ (const char*) "somestring" };
	_LAPPD_sim_app = new TApplication("LAPPDSimApp", &myargc, myargv);
  //Canvases only need to be created, if one wants to display the plots while the program is running
	if (_config_number == 2)
	{
    //size of all canvases is 1280x720
		Double_t canvwidth = 1280;
		Double_t canvheight = 720;

    //Canvas for the histogram of the MC hits of all LAPPDs in one plot
		_LAPPD_MC_all_canvas = new TCanvas("AllLAPPDMCCanvas", "All LAPPDs at once", canvwidth, canvheight);

    //Canvas for the histogram: transverse versus parallel coordinate with arrival time as colour code for each LAPPD
		_LAPPD_MC_canvas = new TCanvas("LAPPDMCCanvas", "2D temperature for one LAPPD", canvwidth, canvheight);
		_LAPPD_MC_canvas->SetRightMargin(0.15);

    //Canvas for the histograms: transverse coordinate versus time and parallel coordinate versus time with number of events as colour code
		_LAPPD_MC_time_canvas = new TCanvas("LAPPDMCTimeCanvas", "LAPPD time versus parallel and transvere coordinate", canvwidth, canvheight);
		_LAPPD_MC_time_canvas->Divide(2, 1, 0.01, 0.01, 0);
    _LAPPD_MC_time_canvas->GetPad(1)->SetRightMargin(0.15);
    _LAPPD_MC_time_canvas->GetPad(2)->SetRightMargin(0.15);

    //Canvas for all waveforms in one plot: Strip number versus time with voltage as colour code
		_LAPPD_all_waveforms_canvas = new TCanvas("LAPPDAllWaveformCanvas", "All waveforms", canvwidth, canvheight);
		_LAPPD_all_waveforms_canvas->Divide(2, 1, 0.015, 0.01, 0);
		_LAPPD_all_waveforms_canvas->GetPad(1)->SetRightMargin(0.15);
		_LAPPD_all_waveforms_canvas->GetPad(2)->SetRightMargin(0.15);

    //Canvas for single waveforms: Voltage versus time
    _LAPPD_waveform_canvas = new TCanvas("LAPPDWaveformCanvas", "One waveform", canvwidth, canvheight);
    _LAPPD_waveform_canvas->Divide(2, 1, 0.015, 0.01, 0);
	}

  //initial opening of the output file
  const char *filename = filePath.c_str();
	_output_file = new TFile(filename, "Recreate");
}

/**
 * Destructor LAPPDDisplay: Closes the output file, deletes the canvas and the TApplication.
 */

LAPPDDisplay::~LAPPDDisplay()
{
	if (_config_number == 2)
	{
		if (gROOT->FindObject("LAPPDSimCanvas") != nullptr)
		{
			delete _LAPPD_MC_time_canvas;
		}
		if (gROOT->FindObject("AllLAPPDMCCanvas") != nullptr)
		{
			delete _LAPPD_MC_all_canvas;
		}
		if (gROOT->FindObject("LAPPDMCCanvas") != nullptr)
		{
			delete _LAPPD_MC_canvas;
		}
		if (gROOT->FindObject("LAPPDAllWaveformCanvas") != nullptr)
		{
			delete _LAPPD_all_waveforms_canvas;
		}
    if (gROOT->FindObject("LAPPDWaveformCanvas") != nullptr)
    {
      delete _LAPPD_waveform_canvas;
    }
	}
	delete _LAPPD_sim_app;
	_output_file->Close();
}

/**
 * [LAPPDDisplay::InitialiseHistoAllLAPPDs description]
 * @param eventNumber [description]
 */
void LAPPDDisplay::InitialiseHistoAllLAPPDs(int eventNumber){
	//Create the name for the histogram for all LAPPDs
	std::string eventnumber = boost::lexical_cast < std::string > (eventNumber);
	string allHitsName = "event" + eventnumber + "AllLAPPDs";
	const char *allHitsNamec = allHitsName.c_str();
	//Initialisation of the histogram for all LAPPDs
	_all_hits = new TH2D(allHitsNamec, allHitsNamec, 200, 0, 180, 200, -1.1, 1.1);
}

/**
 * Method OpenNewFile: This method is needed to avoid the creation of too big .root files.
 *                     The splitting can be adjusted in the LAPPDSim tool.
 * @param filenumber: The number of the current file
 */
void LAPPDDisplay::OpenNewFile(int filenumber){

  //Close the original file
  _output_file->Close();
  std::string number = boost::lexical_cast < std::string > (filenumber);

  //adjust the name for sensible numbering
  _output_file_name.erase(_output_file_name.size()-7);
  if(filenumber < 10)
  {
    _output_file_name = _output_file_name + "0" + number + ".root";
  }
  else
  {
    _output_file_name = _output_file_name + number + ".root";
  }
  const char *filename = _output_file_name.c_str();
  _output_file = new TFile(filename, "Recreate");
}



/**
 * Method MCTruthDrawing: Draws the MCtruth hits as 2D histogram with transverse coordinate as x-axis, parallel coordinate as y-axis and hit time as colour code.
 *                        This is drawn in one canvas, whereas the following two histograms share a canvas.
 *                        Draws also the MCtruth hits as 2D histogram with time as x-axis and transverse coordinate as y-axis and number of events as colour code.
 *                        Draws also the MCtruth hits as 2D histogram with time as x-axis and parallel coordinate as y-axis and number of events as colour code.
 *                        This is done for every LAPPD individually and also for all LAPPDs in an event.
 * @param eventNumber The number of the event, needed for the naming of the histograms.
 * @param lappdmchits Map of all LAPPD hits
 */

void LAPPDDisplay::MCTruthDrawing(int eventNumber, unsigned long actualTubeNo, vector <MCLAPPDHit> mchits)
{
		std::string eventnumber = boost::lexical_cast < std::string > (eventNumber);

		//Create the name for the histogram for individual LAPPDs
		std::string lappdnumber = boost::lexical_cast < std::string > (actualTubeNo);
		string individualName = "event" + eventnumber + "lappd" + lappdnumber;
		const char *individualNamec = individualName.c_str();

		//Initialisation of the histogram for individual LAPPDs
		TH2D* LAPPDMCHits = new TH2D(individualNamec, individualNamec, 60, -0.12, 0.12, 60, -0.12, 0.12);

    //Find the minimal time to create the histograms in a decent range
		double mintime = 1000000;
		for (int k = 0; k < (int) mchits.size(); k++)
		{
			LAPPDHit ahit = mchits.at(k);
			double atime = ahit.GetTime();
			if (mintime > atime)
			{
				mintime = atime;
			}
		}

    //Name of the histograms
		std::string nameTrans = "event" + eventnumber + "lappd" + lappdnumber + "trans";
		std::string namePara = "event" + eventnumber + "lappd" + lappdnumber + "para";
		const char * charNameTrans = nameTrans.c_str();
		const char * charNamePara = namePara.c_str();

    //Initialisation of the histograms
		TH2D* LAPPDTrans = new TH2D(charNameTrans, charNameTrans, 256, mintime-1, mintime + 25.6, 60, -0.12, 0.12);
		TH2D* LAPPDPara = new TH2D(charNamePara, charNamePara, 256, mintime-1, mintime + 25.6, 60, -0.12, 0.12);

    //loop over all MC hits
    for (int i = 0; i < (int) mchits.size(); i++)
		{
      //Getting the hit information
			LAPPDHit ahit = mchits.at(i);
			double atime = ahit.GetTime(); //*1000.;
			vector<double> localpos = ahit.GetLocalPosition();
			double trans = localpos.at(1);
			double para = localpos.at(0);

      //Filling of the histograms
			LAPPDMCHits->Fill(trans, para, atime);
			LAPPDTrans->Fill(atime, trans, 1);
			LAPPDPara->Fill(atime, para, 1);

      //Calculate the hits in the 2D representation of the tank wall and fill them in the histogram for all LAPPDs
      Position globalHit(ahit.GetPosition().at(0), ahit.GetPosition().at(1), ahit.GetPosition().at(2));
			double phi = TMath::ATan2(globalHit.Z(), globalHit.X()) * (180 / TMath::Pi());
			_all_hits->Fill(phi, globalHit.Y(), atime);
		}
    //Cosmetics for the histograms
		LAPPDMCHits->GetXaxis()->SetTitle("Transverse coordinate [m]");
		LAPPDMCHits->GetYaxis()->SetTitle("Parallel coordinate [m]");
		LAPPDMCHits->GetYaxis()->SetTitleOffset(1.4);
		LAPPDMCHits->GetZaxis()->SetTitle("Arrival time [ns]");
		LAPPDMCHits->GetZaxis()->SetTitleOffset(1.4);
		LAPPDMCHits->Write();

		LAPPDTrans->GetYaxis()->SetTitle("Transverse coordinate [m]");
		LAPPDTrans->GetYaxis()->SetTitleOffset(1.4);
		LAPPDTrans->GetXaxis()->SetTitle("Arrival time [ns]");
		LAPPDTrans->GetZaxis()->SetTitle("Events");
    LAPPDTrans->GetZaxis()->SetTitleOffset(1.4);
		LAPPDTrans->Write();

		LAPPDPara->GetYaxis()->SetTitle("Parallel coordinate [m]");
		LAPPDPara->GetYaxis()->SetTitleOffset(1.4);
		LAPPDPara->GetXaxis()->SetTitle("Arrival time [ns]");
		LAPPDPara->GetZaxis()->SetTitle("Events");
    LAPPDPara->GetZaxis()->SetTitleOffset(1.4);
		LAPPDPara->Write();

    //Adjust the canvases and draw the histograms in the right canvases
		if(_config_number == 2)
		{
		_LAPPD_MC_canvas->cd();
		LAPPDMCHits->SetStats(0);
		LAPPDMCHits->Draw("COLZ");
		_LAPPD_MC_canvas->Modified();
		_LAPPD_MC_canvas->Update();

		_LAPPD_MC_time_canvas->cd(1);
		LAPPDTrans->SetMarkerStyle(21);
		LAPPDTrans->SetMarkerSize(5);
		LAPPDTrans->SetStats(0);
		LAPPDTrans->Draw("COLZ");

		_LAPPD_MC_time_canvas->cd(2);
		LAPPDPara->SetMarkerStyle(21);
		LAPPDPara->SetMarkerSize(5);
		LAPPDPara->SetStats(0);
		LAPPDPara->Draw("COLZ");
		_LAPPD_MC_time_canvas->Modified();
		_LAPPD_MC_time_canvas->Update();
		}

		LAPPDMCHits->Clear();
		LAPPDTrans->Clear();
		LAPPDPara->Clear();
}

/**
 * Method FinaliseHistoAllLAPPDs: Cosmetics and drawing of the histogram for all LAPPDs
 */
void LAPPDDisplay::FinaliseHistoAllLAPPDs(){
  //Cosmetics
	_all_hits->GetXaxis()->SetTitle("Radius [m]");
	_all_hits->GetYaxis()->SetTitle("Height [m]");
	_all_hits->GetZaxis()->SetTitle("Arrival time [ns]");
	_all_hits->GetZaxis()->SetTitleOffset(1.4);
	_all_hits->Write();
  //Canvas adjustments
	if(_config_number == 2)
	{
	_LAPPD_MC_all_canvas->cd();
	_all_hits->SetStats(0);
	_all_hits->Draw("COLZ");
	_LAPPD_MC_all_canvas->Modified();
	_LAPPD_MC_all_canvas->Update();
	}
	_all_hits->Clear();
}

/**
 * Method RecoDrawing:   This method draws the waveforms. One histogram for the left and the right side of each strip:
 *                       Strip number as y-axis and time as x-axis and voltage as colour code.
 * @param eventCounter   Number of the event used for the names of the histograms
 * @param tubeNumber     The detector ID of the LAPPD also used for the names of the histograms.
 * @param waveformVector The vector, from which the waveforms can be retrieved
 */
void LAPPDDisplay::RecoDrawing(int eventCounter, unsigned long tubeNumber, std::vector<Waveform<double>> waveformVector)
{
  //Creation of the histogram names
  std::string eventnumber = boost::lexical_cast < std::string > (eventCounter);
	std::string lappdnumber = boost::lexical_cast < std::string > (tubeNumber);
	std::string nameleft = "event" + eventnumber + "lappd" + lappdnumber + "left";
	const char *d = nameleft.c_str();
  std::string nameright = "event" + eventnumber + "lappd" + lappdnumber + "right";
	const char *e = nameright.c_str();
  double nbinsx = 25.6;

  //Initialisation of the histograms
  TH2D* leftAllWaveforms = new TH2D(d, d, 256, 0.0, nbinsx, 30, 0.0, 30.0);
	TH2D* rightAllWaveforms = new TH2D(e, e, 256, 0.0, nbinsx, 30, 0.0, 30.0);

  //loop over all strips
	for (int i = 0; i < 30; i++)
	{
    //Get the samples from the waveforms. For the right side 59-i is used,
    //because in the LAPPDSim tool the stripes range from -30 to 30, with
    //-30 is the left side of the strip and 30 the right side. The waveforms are
    //saved in a way, that 0 equals -30 and 59 equals 30. So for the right side
    //one needs to start with the last entry.
		std::vector<double>* samplesleft = waveformVector[i].GetSamples();
		std::vector<double>* samplesright = waveformVector[59 - i].GetSamples();

    //Creation of histogram names for every strip
    std::string stripNumber = boost::lexical_cast < std::string > (i);
    std::string nameWaveformLeft = "event" + eventnumber + "lappd" + lappdnumber + "strip" + stripNumber + "left";
    std::string nameWaveformRight = "event" + eventnumber + "lappd" + lappdnumber + "strip" + stripNumber + "right";
    const char *waveformLeftChar = nameWaveformLeft.c_str();
    const char *waveformRightChar = nameWaveformRight.c_str();

    //Initialisation of the histograms
    TH1D* waveformLeft = new TH1D(waveformLeftChar, waveformLeftChar, 256, 0, 25.6);
    TH1D* waveformRight = new TH1D(waveformRightChar, waveformRightChar, 256, 0, 25.6);

    //loop over the samples
		for (int j = 0; j < (int) samplesleft->size(); j++)
		{
      //The samples range from 0 to 256, which equals 0 to 25.6 ns.
      //Therefore one needs the sample number divided with 10 to get ns.
      //0.00001 is added to the time to avoid binning effects.
			double time = ((double)j/10) + 0.00001;

      //Filling of the histograms.
      waveformLeft->Fill(time, -samplesleft->at(j));
      waveformRight->Fill(time, -samplesright->at(j));

			leftAllWaveforms->Fill(time, i, -samplesleft->at(j));
			rightAllWaveforms->Fill(time, i, -samplesright->at(j));
		}
    //Cosmetics
    waveformLeft->GetXaxis()->SetTitle("Time [ns]");
    waveformLeft->GetYaxis()->SetTitle("Voltage [mV]");
    waveformLeft->GetYaxis()->SetTitleOffset(1.4);
    waveformLeft->Write();

    waveformRight->GetXaxis()->SetTitle("Time [ns]");
    waveformRight->GetYaxis()->SetTitle("Voltage [mV]");
    waveformRight->GetYaxis()->SetTitleOffset(1.4);
    waveformRight->Write();

    //Canvas adjustments
    if(_config_number == 2){
      _LAPPD_waveform_canvas->cd(1);
      waveformLeft->SetStats(0);
      waveformLeft->Draw("HIST");

      _LAPPD_waveform_canvas->cd(2);
      waveformRight->SetStats(0);
      waveformRight->Draw("HIST");

      _LAPPD_waveform_canvas->Modified();
    	_LAPPD_waveform_canvas->Update();
      }
	}

  //Cosmetics
	leftAllWaveforms->GetXaxis()->SetTitle("Time [ns]");
	leftAllWaveforms->GetYaxis()->SetTitle("Strip number");
	leftAllWaveforms->GetZaxis()->SetTitle("Voltage [mV]");
	//leftAllWaveforms->GetZaxis()->SetTitleOffset(1.4);
	leftAllWaveforms->Write();

	rightAllWaveforms->GetXaxis()->SetTitle("Time [ns]");
	rightAllWaveforms->GetYaxis()->SetTitle("Strip number");
	rightAllWaveforms->GetZaxis()->SetTitle("Voltage [mV]");
	//rightAllWaveforms->GetZaxis()->SetTitleOffset(1.4);
	rightAllWaveforms->Write();

  //Canvas adjustments
	if(_config_number == 2)
	{
	_LAPPD_all_waveforms_canvas->cd(1);
	leftAllWaveforms->SetStats(0);
	leftAllWaveforms->Draw("COLZ");

	_LAPPD_all_waveforms_canvas->cd(2);
	rightAllWaveforms->SetStats(0);
	rightAllWaveforms->Draw("COLZ");

	_LAPPD_all_waveforms_canvas->Modified();
	_LAPPD_all_waveforms_canvas->Update();
	}

	leftAllWaveforms->Clear();
	rightAllWaveforms->Clear();

}
