#include "LAPPDresponse.h"
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <cmath>


//ClassImp(LAPPDresponse)

LAPPDresponse::LAPPDresponse()
{
/************************************************************************************************************************
  Had to move this part to an own function because otherwise the constructor opens the TFile for every WCSim event.
  This means one gets segmentation faults because of too many open instances of the file.
  @author: Malte Stender
*************************************************************************************************************************

  TFile* tf = new TFile("/nashome/m/mstender/ToolAnalysis-MrdEfficiency2/UserTools/LAPPDSim/pulsecharacteristics.root","READ");

  the shape of a typical pulse
 _templatepulse = (TH1D*) tf->Get("templatepulse");
  variations in the peak signal on the central strip
 _PHD = (TH1D*) tf->Get("PHD");

  charge spreading of a pulse in the transverse direction (in mm)
  as a function of nearness to strip center. The charge tends to
  spread more in the transverse direction if the centroid of the
  signal is between two striplines
 _pulsewidth = (TH1D*) tf->Get("pulsewidth");

  structure to store the pulses, count them, and organize them by channel
  _pulseCluster = new LAPPDpulseCluster()  This is no longer needed, kept for reference for now

  random numbers for generating noise
  mrand = new TRandom3();
  */
}
LAPPDresponse::~LAPPDresponse()
{

}

void LAPPDresponse::Initialise(TFile* tf){
  // the shape of a typical pulse
  _templatepulse = (TH1D*) tf->Get("templatepulse");
  // variations in the peak signal on the central strip
  _PHD = (TH1D*) tf->Get("PHD");

  // charge spreading of a pulse in the transverse direction (in mm)
  // as a function of nearness to strip center. The charge tends to
  // spread more in the transverse direction if the centroid of the
  // signal is between two striplines
  _pulsewidth = (TH1D*) tf->Get("pulsewidth");

  // structure to store the pulses, count them, and organize them by channel
  //_pulseCluster = new LAPPDpulseCluster()  This is no longer needed, kept for reference for now

  // random numbers for generating noise
  mrand = new TRandom3();
}

void LAPPDresponse::AddSinglePhotonTrace(double trans, double para, double time)
{
  // Draw a random value for the peak signal peak
  double peak = (_PHD->GetRandom())/10.;

  // find nearest strip
  int neareststripnum = this->FindStripNumber(trans);

  // calculate distance from nearest strip center
  double offcenter = fabs(trans - (this->StripCoordinate(neareststripnum))) ; //trans - striptrans;

  //std::cout<<"THE PEAK "<<peak<<std::endl;

  // width of the charge sharing
  double thesigma =  _pulsewidth->Interpolate(offcenter);

  //std::cout << "/* message */" << '\n';std::cout<<"nearest stripnum: "<<neareststripnum<<" off center: "<<offcenter<<" thesigma "<<thesigma<<std::endl;

  TF1* theChargeSpread = new TF1("theChargeSpread","gaus",-100,100);
  theChargeSpread->SetParameter(0,peak);
  theChargeSpread->SetParameter(1,0.0);
  theChargeSpread->SetParameter(2,thesigma);

  // calculate distances and times in the parallel direction
  double leftdistance = fabs(-114.554 - para); // annode is 229.108 mm in parallel direction
  double rightdistance = fabs(114.554 - para);

  if(leftdistance+rightdistance!=229.108) std::cout<<"WHAT!? "<<(leftdistance+rightdistance)<<std::endl;;

  double lefttime = leftdistance/(0.53*(0.299792458)); // 53% speed of light (picoseconds per mm) on transmission lines
  double righttime = rightdistance/(0.53*(0.299792458)); // 53% speed of light (picoseconds per mm) on transmission lines

  //std::cout<<leftdistance<<" "<<rightdistance<<" "<<lefttime<<" "<<righttime<<std::endl;

  //loop over five-strip cluster about the central strip
  for(int i=0; i<5; i++){

    int wstrip = (neareststripnum-2)+i;
    double wtrans = this->StripCoordinate(wstrip);
    double wspeak = theChargeSpread->Eval(trans-wtrans);

    //signal has to be larger than 0.5 mV
    if( (wspeak>0.5) && (wstrip>0) && (wstrip<31) ) {
      int tubeid = 0;
      double charge =0;
      double low = 0;
      double hi = 0;

      //std::cout<<"which strip "<<wstrip<<" peakvalue"<<wspeak<<std::endl;


      LAPPDPulse pulse(tubeid, wstrip, (time + righttime)/1000., charge, wspeak, low, hi);  //SD
      if(LAPPDPulseCluster.count(wstrip)==1){
        std::vector<LAPPDPulse> tempVector;
        tempVector = LAPPDPulseCluster.at(wstrip);
        tempVector.push_back (pulse);
        LAPPDPulseCluster.at(wstrip)=tempVector;
        //add pulse to already existing vector at key  SD
      }
      else{
        std::vector<LAPPDPulse> PulseVector ;
        PulseVector.push_back (pulse);
        LAPPDPulseCluster.insert (std::pair<int,vector<LAPPDPulse>>(wstrip,PulseVector) );//SD
        //create vector at key and add pulse into that vector SD
      }

      pulse.SetChannelID(-1.0*wstrip); //SD
      pulse.SetTime((time +lefttime)/1000.); //SD
      if(LAPPDPulseCluster.count(-wstrip)==1){
        std::vector<LAPPDPulse> tempVector;
        tempVector = LAPPDPulseCluster.at(-wstrip);
        tempVector.push_back (pulse);
        LAPPDPulseCluster.at(-wstrip)=tempVector;
        //add pulse to already existing vector at key  SD
      }
      else{
        std::vector<LAPPDPulse> PulseVector ;
        PulseVector.push_back (pulse);
        LAPPDPulseCluster.insert (std::pair<int,vector<LAPPDPulse>>(-wstrip,PulseVector) );   //SD
        //create vector at key and add pulse into that vector SD
      }
    }
  }

  //std::cout<<"Done Adding Pulse"<<std::endl;
  delete theChargeSpread;
}


Waveform<double> LAPPDresponse::GetTrace(int CHnumber, double starttime, double samplesize, int numsamples, double thenoise)
{

  // parameters for the histogram of the scope trace
  double lowend = (starttime-(samplesize/2.));
  double upend = lowend + samplesize*((double)numsamples);
  TString tracename;
  tracename += "trace_";
  tracename += CHnumber;
  //tracename += "_";
  //if(parity==1) tracename+="right";
//  else tracename+="left";
  //create said histogram

  TH1D* trace = new TH1D(tracename,tracename,numsamples,lowend,upend);
  Waveform<double> wav_trace;
  //if there are no pulses on the strip, just generate white noise
  if(LAPPDPulseCluster.count(CHnumber)==0) {   //SD
    for(int j=0; j<numsamples; j++){

      double mnoise = thenoise*(mrand->Rndm()-0.5);

      trace->SetBinContent(j+1, mnoise);

    }
  }
  else{

    //if there are pulses on the strip, loop over the N pulses on that strip
	std::vector<LAPPDPulse> tempoVector = LAPPDPulseCluster.at(CHnumber);   //SD
    for(int k=0; k<(int)tempoVector.size(); k++){           //SD
      //  for(int k=0; k<4; k++){


      //looks up the index number for pulse "k" on strip "CHnumber"
      //int wPulse = _pulseCluster->GetPulseNum(CHnumber,k);
      //gets the pulse with that index number
      //LAPPDpulse* mpulse = _pulseCluster->GetPulse(wPulse);
      //peak value of the signal on that strip
      //double peakv = mpulse->Getpeakvalue();
      double peakv=tempoVector.at(k).GetPeak();
      //arrival time of the pulse
      //double ptime = mpulse->Getpulsetime();
      double thetime=0;
      double ptime= thetime;
      //transit time of the pulse along the strip
      double stime=tempoVector.at(k).GetTime()*1000.;
      //looking at the signal to the left (parity=-1) or right (parity=1)?

      //if(parity<0) stime = mpulse->Getlefttime();
      //else stime = mpulse->Getrighttime();
      //	 if(parity<0) stime=0.5;    //SD out
      //     else stime = 0.3;    //SD out

      //sum the pulse arrival time with transit time on the strip to determine
      //when the signal will arrive
      double tottime = ptime + stime;

      //loop over number of samples

      for(int j=0; j<numsamples; j++){

        //get the global time when each sample is acquired
        double bcent = trace->GetBinCenter(j+1);
        double mbincontent=0.0;
        double mnoise=0.0;

        //only add the noise on ONCE
        if(k==0) mnoise = thenoise*(mrand->Rndm()-0.5);
        mbincontent+=mnoise;

        //if the sample time actually falls in the window for when the pulse
        //should arrive, evaluate the pulse value at that sample point
        if( (bcent > tottime) && (bcent< tottime+3000) ) mbincontent+=(peakv*(_templatepulse->Interpolate(bcent-tottime)));

        //add this on to the contributions to the trace from previous pulses
        double obincontent = trace->GetBinContent(j+1);
        trace->SetBinContent(j+1,obincontent+mbincontent);

      }
    }
  }

  for (int i = 1; i <= numsamples; i++){
      wav_trace.PushSample(trace->GetBinContent(i));
  }

  delete trace;
  return wav_trace;
}


int LAPPDresponse::FindStripNumber(double trans){

  double newtrans = trans + 101.6;

  // the first and last strips have a different width
  int stripnum=-1;
  if(newtrans<5.765) stripnum = 1;
  if(newtrans>197.435) stripnum = 30;

  double stripdouble;
  if(stripnum==-1){
    // divide the 28 remaining strips into the remaining area
    double stripdouble = 28.0*((newtrans-5.765)/(203.2 - 11.53));
    stripnum = 2 + floor(stripdouble);
  }

  return stripnum;
}


double LAPPDresponse::StripCoordinate(int stripnumber){

  double coor = -55555.;
  // the first and last strips have a different width
  if(stripnumber==1) coor = (2.31-101.6);
  if(stripnumber==30) coor = (101.6-2.31);

  if( stripnumber>1 && stripnumber<30 ){
    // remaining 28 strips have the same spacing
    coor= (5.765-101.6+3.455) + (stripnumber-2)*6.91;
  }

  return coor;

}
