#include "LAPPDcfd.h"

LAPPDcfd::LAPPDcfd():Tool(){}


bool LAPPDcfd::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  TString CIWL;
  m_variables.Get("CFDInputWavLabel",CIWL);
  CFDInputWavLabel = CIWL;
  // Get the CFD threshold from the config file
  m_variables.Get("Fraction_CFD", Fraction_CFD);
  //std::cout<<"Fraction_CFD="<<Fraction_CFD<<std::endl;

  isSim=false;
  // Check in the Boost Store whether this is a simulated event or not
  m_data->Stores["ANNIEEvent"]->Header->Get("isSim",isSim);

  return true;
}


bool LAPPDcfd::Execute(){

  Waveform<double> bwav;

  // get raw lappd data from the Boost Store
  std::map<int,vector<Waveform<double>>> rawlappddata;
  bool testval =  m_data->Stores["ANNIEEvent"]->Get(CFDInputWavLabel,rawlappddata);

  // get first-level pulse reco from the Boost Store
  std::map<int,vector<LAPPDPulse>> SimpleRecoLAPPDPulses;
  m_data->Stores["ANNIEEvent"]->Get("SimpleRecoLAPPDPulses",SimpleRecoLAPPDPulses);

  // get the sim-level information

  std::map<int,vector<LAPPDHit>> lappdmchits;
  if(isSim){
    m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHit",lappdmchits);
  }

  // Place to store the reconstructed pulses
  std::map<int,vector<LAPPDPulse>> CFDRecoLAPPDPulses;

  // Loop over all channels
  map <int, vector<Waveform<double>>> :: iterator itr;
  for (itr = rawlappddata.begin(); itr != rawlappddata.end(); ++itr){

    // Get the channel number and a vector of Waveforms
    int channelno = itr->first;
    vector<Waveform<double>> Vwavs = itr->second;

    // get the vector of pulses correseponding to the channel
    map<int, vector<LAPPDPulse>>::iterator p;
    p = SimpleRecoLAPPDPulses.find(channelno);
    vector<LAPPDPulse> Vpulses = p->second;

//    std::cout<<"************************************************"<<std::endl;
//    std::cout<<"IN LAPPDCFD:: channel: "<<channelno<<std::endl;

    if(isSim){
      // If the data is simulated data, we loop over the true hits
      // and print them to screen
//      std::cout<<"simulated pulse: ";
      map<int, vector<LAPPDHit>>::iterator p;
      p = lappdmchits.find(0);
      vector<LAPPDHit> Vhits = p->second;
/*
      std::cout<<Vhits.size()<<" simulated pulses. ";
      for(int np=0; np<Vhits.size(); np++){
        std::cout<<"p"<<np<<"=(t="<<(Vhits.at(np)).GetTpsec()<<") ";
      }
      std::cout<<" "<<std::endl;
*/
    }

    std::vector<LAPPDPulse> thepulses;

    //loop over all Waveforms
    for(int i=0; i<(int)Vwavs.size(); i++){

        Waveform<double> bwav = Vwavs.at(i);

        //std::cout<<"reconstructed pulses: "<<std::endl;;

        // loop over all candidate pulses on each waveform, as determined by the LAPPDFindPeak Tool
        for(int j=0; j<(int)Vpulses.size(); j++){

          // for each pulse on the Waveform find the time using CFD1 algorithm
          double cfdtime = CFD_Discriminator1(bwav.GetSamples(),Vpulses.at(j));
          //std::cout<<"for pulse #"<<j<<" (Q="<<(Vpulses.at(j)).GetCharge()<<",Amp="<<(Vpulses.at(j)).GetPeak()<<",LowRange="<<(Vpulses.at(j)).GetLowRange()<<",HiRange="<<(Vpulses.at(j)).GetHiRange()<<") "<<"  cfd_time="<<cfdtime<<std::endl;

          // Store the reconstructed time in a new LAPPDPulse
          LAPPDPulse apulse(0,channelno,(cfdtime/1000.),(Vpulses.at(j)).GetCharge(),(Vpulses.at(j)).GetPeak(),(Vpulses.at(j)).GetLowRange(),(Vpulses.at(j)).GetHiRange());
          thepulses.push_back(apulse);
        }
      }
        //std::cout<<" "<<std::endl;
        // Put the newly reconsructed LAPPDPulses into a map, by channel
        CFDRecoLAPPDPulses.insert(pair <int,vector<LAPPDPulse>> (channelno,thepulses));
    }

    // Add the CFD reconstructed information to the Boost Store
    m_data->Stores["ANNIEEvent"]->Set("CFDRecoLAPPDPulses",CFDRecoLAPPDPulses);



  return true;
}


bool LAPPDcfd::Finalise(){

  return true;
}


double	LAPPDcfd::CFD_Discriminator1(std::vector<double>* trace, LAPPDPulse pulse) {

  double amp = pulse.GetPeak();
  double time = 0;
  // this should be a global variable that gets set from a config...
  double th = Fraction_CFD * amp;
  double eps = 1e-2;
  double FitWindow_min = (pulse.GetLowRange()-5)*100.;
  double FitWindow_max = (pulse.GetHiRange()+5)*100.;

  double PointsPerSpline = 5;

  // this is kludgy...I'm converting the trace to a histo
  // (many parameters are hard coded)
  int nbins = trace->size();
  double starttime=0.;
  double endtime = starttime + ((double)nbins)*100.;
  TH1D* hwav = new TH1D("hwav","hwav",nbins,starttime,endtime);

  for(int i=0; i<nbins; i++){
    hwav->SetBinContent(i+1,-trace->at(i));
  }

//	int bin = pulse.GetTpsec();
//		while(fsplinewav->V(hwav->GetBinCenter(bin))<=th) {bin--;}
	double xlow = FitWindow_min;
//	double xhigh = hwav->GetBinCenter(bin);
  double xhigh = FitWindow_max;
	double xmid = (xlow+xhigh)/2;

  if(hwav->Interpolate(xmid)-th==0) time = xmid;
  //std::cout<<"^^^ "<<xhigh<<" "<<xlow<<" "<<th<<" "<<hwav->Interpolate(xmid)<<std::endl;

  // gradually moving the high and low range of the search towards the point where the
  // pulse crosses the threshold. Stop when the low and high range are closer than eps
  while ((xhigh-xlow) >= eps ){

    xmid = (xlow + xhigh)/2.;
    if(hwav->Interpolate(xmid)-th==0) time = xmid;
    if ((hwav->Interpolate(xmid)-th) > 0)
    xhigh = xmid;
    else
    xlow = xmid;
  }

  time = xlow;
  delete hwav;
	return time;
}


double LAPPDcfd::CFD_Discriminator2(std::vector<double>* trace, LAPPDPulse pulse){

  double time=0;

//  sprintf(CFDSplineName, "CFDSpline_ch%d_%d", chno, evtno);
//  sprintf(CFDSplineTitle, "CFDSpline | CFDSpline_ch%d_%d", chno, evtno);


//  if(amp <= AmpThreshold || npeaks==0) { time = 0; } //TDC threshold
//	else {
//    if(IfDynamicWindow == 1) Calculate_fitrange();
/*

  int Delay_CFDBin = (int)(Delay_CFD/Deltat);
  for(int i=0;i<DimSize;i++) {hcfd->SetBinContent(i+1, hwav->GetBinContent(i+1));} //copy
  TH1D *hinv = new TH1D("hinv","hinv",DimSize,0,DimSize*Deltat);
  for(int i=0;i<DimSize;i++) {
  		if(i+Delay_CFDBin<DimSize) hinv->SetBinContent(i+1, hwav->GetBinContent(i+1+Delay_CFDBin)); //copy and delay
  		else hinv->SetBinContent(i+1,0);
  }
  hinv->Scale(-1. * Fraction_CFD); //inverse and attenuate
  hcfd->Add(hinv);
  fsplinecfd = new TSplineFit(CFDSplineName,CFDSplineTitle,50,PointsPerSpline, hcfd,FitWindow_min, FitWindow_max);
  //		fsplinecfd->UpdateFile(true); //SplineFit database:  Skip this line for fast analysis
  fsplinecfd->ReduceMemory();

  //solve equation
  double eps = 1e-4;
  double xhigh = hcfd->GetMinimumBin()*Deltat;
  //		double xlow = hcfd->GetMaximumBin()*Deltat;
  double xlow = FitWindow_min;
  //		cout<<xlow<<"\t"<<xhigh<<endl;
  double xmid = (xlow+xhigh)/2;

  if(fsplinecfd->V(xmid)==0) time = xmid;

  while ((xhigh-xlow) >= eps) {
  	xmid = (xlow + xhigh) / 2;
  	if (fsplinecfd->V(xmid) == 0)
  			time =  xmid;
  		if (fsplinecfd->V(xlow)*fsplinecfd->V(xmid) < 0)
  			xhigh = xmid;
  		else
  			xlow = xmid;
  	}
  	time = xlow;
  	delete hinv;
  	hinv = NULL;
  }
  */
  return time;
}
