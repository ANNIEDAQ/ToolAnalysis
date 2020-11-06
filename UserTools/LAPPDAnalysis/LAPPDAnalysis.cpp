#include "LAPPDAnalysis.h"

LAPPDAnalysis::LAPPDAnalysis():Tool(){}


bool LAPPDAnalysis::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool LAPPDAnalysis::Execute(){



  //BaselineSubtract
    Waveform<double> bwav;
    // get raw lappd data
    map<int,vector<Waveform<double>>> rawlappddata;
    m_data->Stores["ANNIEEvent"]->Get("RawLAPPDData",rawlappddata);
    

    // the filtered Waveform
    std::map<int,vector<Waveform<double>>> blsublappddata;

    map <int, vector<Waveform<double>>> :: iterator itr;
    for (itr = rawlappddata.begin(); itr != rawlappddata.end(); ++itr){
      int channelno = itr->first;
      vector<Waveform<double>> Vwavs = itr->second;
      vector<Waveform<double>> Vfwavs;

      //loop over all Waveforms
      for(int i=0; i < (int) Vwavs.size(); i++){

          Waveform<double> bwav = Vwavs.at(i);
          Waveform<double> blswav = SubtractSine(bwav);
          Vfwavs.push_back(blswav);
        }

        blsublappddata.insert(pair <int,vector<Waveform<double>>> (channelno,Vfwavs));
      }

    m_data->Stores["ANNIEEvent"]->Set("BLsubtractedLAPPDData",blsublappddata);








  return true;
}




























Waveform<double> LAPPDAnalysis::SubtractSine(Waveform<double> iwav) {

  Waveform<double> subWav;

  int nbins = iwav.GetSamples()->size();
  double starttime=0.;
  double endtime = starttime + ((double)nbins)*100.;
  TH1D* hwav_raw = new TH1D("hwav_raw","hwav_raw",nbins,starttime,endtime);

  for(int i=0; i<nbins; i++){
    hwav_raw->SetBinContent(i+1,iwav.GetSample(i));
    hwav_raw->SetBinError(i+1,0.1);

  }

  TF1* sinit = new TF1("sinit","([0]*sin([2]*x+[1]))",0,DimSize*Deltat);
  sinit->SetParameter(0,0.4);
  sinit->SetParameter(1,0.0);
  sinit->SetParameter(2,0.0);
  sinit->SetParameter(2,0.00055);
  sinit->SetParLimits(2,0.0003,0.0008);
  sinit->SetParLimits(0,0.,1.0);

  hwav_raw->Fit("sinit","QNO","",LowBLfitrange,HiBLfitrange);
  //cout<<"Parameters: "<< sinit->GetParameter(3)<<" "<<LowBLfitrange<<" "<<HiBLfitrange<<endl;

  for(int j=0; j<nbins; j++){

    subWav.PushSample((hwav_raw->GetBinContent(j+1))-(sinit->Eval(hwav_raw->GetBinCenter(j+1))));
  }

  delete hwav_raw;
  delete sinit;
  return subWav;
}















bool LAPPDAnalysis::Finalise(){

  return true;
}
