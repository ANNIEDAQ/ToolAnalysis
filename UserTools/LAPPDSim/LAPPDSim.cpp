#include "LAPPDSim.h"
#include <unistd.h>

LAPPDSim::LAPPDSim():Tool(),myTR(nullptr),_tf(nullptr),_event_counter(0),_file_number(0),_display_config(0),_is_artificial(false),_display(nullptr),_geom(nullptr),LAPPDWaveforms(nullptr)
{
}

bool LAPPDSim::Initialise(std::string configfile, DataModel &data)
{

	/////////////////// Usefull header ///////////////////////
	if (configfile != "")
	m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	m_data = &data; //assigning transient data pointer

	//Get the config parameters and print them

	//File path to pulsecharacteristics.root
  std::string pulsecharacteristicsFile;
	m_variables.Get("PathToPulsecharacteristics", pulsecharacteristicsFile);
	std::cout << "Path to pulsecharacteristics.root: " << pulsecharacteristicsFile << std::endl;
	const char * pulsecharacteristicsFileChar = pulsecharacteristicsFile.c_str();

	//Config number for the displaying
	m_variables.Get("EventDisplay", _display_config);
	std::cout << "DisplayNumber " << _display_config << std::endl;

	//Whether artifical events or MC events are used
	m_variables.Get("ArtificialEvent", _is_artificial);
	if(_is_artificial)
	{
		std::cout << "Artifical events will be used." << std::endl;
	}
	else
	{
		std::cout << "MC events will be used." << std::endl;
	}

	//Path to the file to write the output to
	//Since one file is not enough to save all events of a regular WCSim file, there will be several .root files in the specified path.
	std::string outputFile;
	m_variables.Get("OutputFile", outputFile);
	std::cout << "OutputFile " << outputFile << std::endl;
	outputFile.erase(outputFile.size()-5);
	outputFile = outputFile + "00.root";

	//Get the Geometry information
	bool testgeom = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", _geom);
	if (not testgeom)
	{
		std::cerr << "LAPPDSim Tool: Could not find Geometry in the ANNIEEvent!" << std::endl;
		return false;
	}

	// This quantity should be set to false if we are working with real data later
	//bool isSim = true;
	//m_data->Stores["ANNIEEvent"]->Header->Set("isSim",isSim);

	// initialize the ROOT random number generator
	myTR = new TRandom3();
	_tf = new TFile(pulsecharacteristicsFileChar, "READ");

	if (_display_config > 0)
	{
		_display = new LAPPDDisplay(outputFile, _display_config);
	}
	return true;
}

bool LAPPDSim::Execute()
{
	std::cout << "Executing LAPPDSim; event counter " << _event_counter << std::endl;

	//The files become too large, if one tries to save all WCSim events into one file.
	//Every 100 events get a new file.
	if(_event_counter == (20 * (_file_number + 1)))
	{
		_display->OpenNewFile(_file_number);
		_file_number++;
	}
	//Initialise the histogram for displaying all LAPPDs at once
	if (_display_config > 0)
	{
		_display->InitialiseHistoAllLAPPDs(_event_counter);
	}

	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	//Artifical Events: This is only used for test purposes and is just a quick and dirty solution to try some things.
	//This is the reason, why there are no comments

	if(_is_artificial)
	{
		vector<MCLAPPDHit> artificialHits;

	std::map<std::string, std::map<unsigned long,Detector*> >* AllDetectors = _geom->GetDetectors();
	std::map<std::string, std::map<unsigned long,Detector*> >::iterator itGeom;
	for(itGeom = AllDetectors->begin(); itGeom != AllDetectors->end(); ++itGeom){
		if(itGeom->first == "LAPPD"){
			std::map<unsigned long,Detector*> LAPPDDetectors = itGeom->second;
  		std::map<unsigned long, Detector*>::iterator itDet;
					for(itDet = LAPPDDetectors.begin(); itDet != LAPPDDetectors.end(); ++itDet){

					LAPPDresponse response;
					response.Initialise(_tf);
					artificialHits.clear();
					int detectorID = itDet->second->GetDetectorID();
					Position LAPPDPosition = itDet->second->GetDetectorPosition();
					Position LAPPDDirection(itDet->second->GetDetectorDirection().X(),itDet->second->GetDetectorDirection().Y(),itDet->second->GetDetectorDirection().Z());
					Position normalHeight(0,1,0);
					Position side = normalHeight.Cross(LAPPDDirection);
					std::vector<int> parents{0,0};
					double charge = 1.0;
					for(int i = 0; i < 2; i++){
						double timeNs = 1.0 + 10 * i;
						double paraMeter = 0.0;
						double transMeter = 0.0;
						// if(detectorID == 0){
						// 	timeNs = 1.0;
						// 	paraMeter = 0.1;
						// 	transMeter = 0.0;
						// }

						std::vector<double> localPosition{paraMeter, transMeter};
						std::vector<double> globalPosition{LAPPDPosition.X()+paraMeter*side.X(), LAPPDPosition.Y()+transMeter, LAPPDPosition.Z()+paraMeter*side.Z()};
						MCLAPPDHit firstHit(detectorID, timeNs, charge, globalPosition, localPosition, parents);

						artificialHits.push_back(firstHit);
						double trans = transMeter * 1000;
						double para = paraMeter * 1000;
						double time = timeNs * 1000;
						response.AddSinglePhotonTrace(trans, para, time);
					}
					if (_display_config > 0)
					{
						_display->MCTruthDrawing(_event_counter, detectorID, artificialHits);
					}
					vector<Waveform<double>> Vwavs;

					for (int i = -30; i < 31; i++)
					{
						if (i == 0)
						{
							continue;
						}
						Waveform<double> awav = response.GetTrace(i, 0.0, 100, 256, 1.0);
						Vwavs.push_back(awav);
					}

					if (_display_config > 0)
					{
						_display->RecoDrawing(_event_counter, detectorID, Vwavs);

						if (_display_config == 2)
						{
							do
							{
								std::cout << "Press a key to continue..." << std::endl;
							} while (cin.get() != '\n');

							std::cout << "Continuing" << std::endl;
						}
					}
					Vwavs.clear();
					if(detectorID > 1){
						break;
					}
				}
			}
		}
	}
//---------------------------------------------------------------------------------------------------------------------
//MC events: Here is the implementation for the MC events
	else
	{
		//storage for the waveforms
		LAPPDWaveforms = new std::map<unsigned long, Waveform<double> >;
		LAPPDWaveforms->clear();
		// get the MC Hits
		std::map<unsigned long, std::vector<MCLAPPDHit> >* lappdmchits;
		bool testval = m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits", lappdmchits);
		if (not testval)
		{
			std::cerr << "LAPPDSim Tool: Could not find MCLAPPDHits in the ANNIEEvent!" << std::endl;
			return false;
		}

		// loop over the number of lappds
		std::map<unsigned long, std::vector<MCLAPPDHit> >::iterator itr;
		for (itr = lappdmchits->begin(); itr != lappdmchits->end(); ++itr)
		{
			//Get the Channelkey
			unsigned long tubeno = itr->first;

			//Retrieve the detector object with the Channelkey
			Detector* thelappd = _geom->ChannelToDetector(tubeno);

			//Use the detector object to get the detector ID
			unsigned long actualTubeNo = thelappd->GetDetectorID();

			//Get the Hits on the LAPPD
			vector<MCLAPPDHit> mchits = itr->second;

			//If display is active, draw histograms showing the MC hits on the LAPPDs
			if (_display_config > 0)
			{
				_display->MCTruthDrawing(_event_counter, actualTubeNo, mchits);
			}

			//Create an object of the LAPPDresponse class, which is used for the electronics simulation
			LAPPDresponse response;
			response.Initialise(_tf);

			//loop over the hits on each lappd
			for (int j = 0; j < (int) mchits.size(); j++)
			{
				LAPPDHit ahit = mchits.at(j);
				//Time is in [ns], we need [ps] for the LAPPDrespnse class' methods.
				double atime = ahit.GetTime()*1000.;
				//local position is in [m], we need [mm] for the LAPPDresponse class' methods.
				vector<double> localpos = ahit.GetLocalPosition();
				double trans = localpos.at(1) * 1000;
				double para = localpos.at(0) * 1000;
				//Add the traces to retrieve them later
				response.AddSinglePhotonTrace(trans, para, atime);
			}

			vector<Waveform<double>> Vwavs;
			Vwavs.clear();
			//loop over the channels on each LAPPD
			//Positive numbers describe one side, negative numbers the other side in a way, that left number = -right numbers
			for (int i = -30; i < 31; i++)
			{
				if (i == 0)
				{
					continue;
				}
				//Retrive the traces, which were stored with the AddSinglePhotonTrace method
				Waveform<double> awav = response.GetTrace(i, 0.0, 100, 256, 1.0);
				Vwavs.push_back(awav);
			}

			//Get the channels of each LAPPD
			std::map<unsigned long, Channel>* lappdchannel = thelappd->GetChannels();
			int numberOfLAPPDChannels = lappdchannel->size();
			std::map<unsigned long, Channel>::iterator chitr;
			//Loop over all channels for the assignment of the waveforms to the channels for storing the waveforms
			for (chitr = lappdchannel->begin(); chitr != lappdchannel->end(); ++chitr)
			{
				Channel achannel = chitr->second;
				//achannel->Print();
				//This assignment uses the following numbering scheme:
				//Channelkey 0-29 is the one side, Channelkey 30-59 is the other side in a way that 0 is the left side of the strip, where 30 denotes the right side.
				if (achannel.GetStripSide() == 0)
				{
					LAPPDWaveforms->insert(pair<unsigned long, Waveform<double>>(achannel.GetChannelID(), Vwavs[achannel.GetStripNum()]));
				}
				else
				{
					LAPPDWaveforms->insert(pair<unsigned long, Waveform<double>>(achannel.GetChannelID(), Vwavs[numberOfLAPPDChannels - achannel.GetStripNum() - 1]));
				}

			}
			//Waveforms are drawn
			if (_display_config > 0)
			{
				_display->RecoDrawing(_event_counter, actualTubeNo, Vwavs);
				//This command is used to display a set of histograms and then wait for input to display the next set.
				if (_display_config == 2)
				{
					do
					{
						std::cout << "Press a key to continue..." << std::endl;
					} while (cin.get() != '\n');

					std::cout << "Continuing" << std::endl;
				}
			}

		}				//end loop over LAPPDs
	} //end else of if(_is_artificial)

	if (_display_config > 0)
	{
		_display->FinaliseHistoAllLAPPDs();
	}

	//The waveforms are only saved if MC events are used.
	//The artifical events are not meant to be saved, because they cannot be used in any other tool,
	//since there won't be any hit information in the MCHits or MCLAPPDHits
	if(!_is_artificial)
	{
		std::cout << "Saving waveforms to store" << std::endl;

		m_data->Stores.at("ANNIEEvent")->Set("LAPPDWaveforms", LAPPDWaveforms, true);
	}
	_event_counter++;

	return true;
}

bool LAPPDSim::Finalise()
{
	_tf->Close();
	_display->~LAPPDDisplay();
	return true;
}
