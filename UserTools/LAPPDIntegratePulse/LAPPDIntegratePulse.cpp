#include "LAPPDIntegratePulse.h"

LAPPDIntegratePulse::LAPPDIntegratePulse():Tool(){}

// THIS CODE INTEGRATES AN LAPPD WAVEFORM OVER A FIXED RANGE DEFINED IN THE CONFIG FILE

bool LAPPDIntegratePulse::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Make note that the pulse has been integrated
  bool isIntegrated = true;
  m_data->Stores["ANNIEEvent"]->Header->Set("isIntegrated",isIntegrated);

  // Get from the configuration file the parameters of the waveform
  m_variables.Get("Nsamples", DimSize);
  m_variables.Get("SampleSize",Deltat);

  // Get from the config file, the range over which to integrate
  m_variables.Get("IntegLow",lowR);
  m_variables.Get("IntegHi",hiR);

  return true;
}


bool LAPPDIntegratePulse::Execute(){

  // get raw lappd data
  std::map<int,vector<Waveform<double>>> rawlappddata;
  bool testval =  m_data->Stores["ANNIEEvent"]->Get("RawLAPPDData",rawlappddata);

  std::map<int,vector<double>> thecharge;

  map <int, vector<Waveform<double>>> :: iterator itr;
  for (itr = rawlappddata.begin(); itr != rawlappddata.end(); ++itr){
    int channelno = itr->first;
    vector<Waveform<double>> Vwavs = itr->second;

    //loop over all Waveforms
    vector<double> acharge;
    for(int i=0; i<(int)Vwavs.size(); i++){

        Waveform<double> bwav = Vwavs.at(i);

        // integrate the pulse from the low range to high range in units of mV*psec
        double Qmvpsec = CalcIntegral(bwav,lowR,hiR);
        // convert to coulomb
        double Qcoulomb = Qmvpsec/(1000.*50.*1e12);
        //convert to number of electrons
        double Qelectrons = Qcoulomb/(1.60217733e-19);
        acharge.push_back(Qelectrons);
        //cout<<channelno<<" "<<Qelectrons<<endl;
      }

      // store the charge vector by channel
      thecharge.insert(pair <int,vector<double>> (channelno,acharge));
    }

    // add the charge information to the Boost Store
    m_data->Stores["ANNIEEvent"]->Set("theCharges",thecharge);


  return true;
}


bool LAPPDIntegratePulse::Finalise(){

  return true;
}

double LAPPDIntegratePulse::CalcIntegral(Waveform<double> hwav, double lowR, double hiR){

  double sT=0.; // currently hard coded

  int lowb = (lowR - sT)/Deltat;
  int hib = (hiR - sT)/Deltat;

  double tQ=0.;

  if( (lowb>=0) && (hib<(int)hwav.GetSamples()->size()) ){
    for(int i=lowb; i<hib; i++){
      tQ+=((-hwav.GetSample(i))*Deltat);
    }
  } else std::cout<<"OUT OF RANGE!!!!";

  return tQ;
}
