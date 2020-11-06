#include "LAPPDlasertestHitFinder.h"


const double c = 29.98; //cm/ns
const int side = 1;



LAPPDlasertestHitFinder::LAPPDlasertestHitFinder():Tool(){}


bool LAPPDlasertestHitFinder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("stripID1",stripID1);
  m_variables.Get("stripID2",stripID2);
  m_variables.Get("stripID3",stripID3);
  m_variables.Get("MaxTimeWindow", MaxTimeWindow);
  m_variables.Get("MinTimeWindow", MinTimeWindow);
  m_variables.Get("TwoSided", TwoSided);
  m_variables.Get("CenterChannel", CenterChannel);
  m_variables.Get("PTRange", PTRange);
  return true;
}


bool LAPPDlasertestHitFinder::Execute(){
  stripID.push_back(stripID1);
  stripID.push_back(stripID2);
  stripID.push_back(stripID3);
std::vector<double> AbsPosition;
AbsPosition.push_back(0);
AbsPosition.push_back(0);
AbsPosition.push_back(0);

std::vector<double> LocalPosition;
std::vector<LAPPDPulse> NeighboursPulses;
std::vector<LAPPDPulse> HitPulses;
double ParaPosition;
double PerpPosition;
double HitTime;

std::map<int,vector<LAPPDPulse>> LAPPDPulseCluster;
m_data->Stores["ANNIEEvent"]->Get("CFDRecoLAPPDPulses",LAPPDPulseCluster);
//m_data->Stores["ANNIEEvent"]->Get("SimpleRecoLAPPDPulses",LAPPDPulseCluster);

   LAPPDPulse MaxPulse;
double MaxAmp =-1.0;

//Find the biggest pulse inside the time window
  for (int i = 0; i < 3; i++){
    std::vector<LAPPDPulse> TempVector;
    TempVector = LAPPDPulseCluster.at(i);
     std::cout<<"Temp Vector Size "<<TempVector.size()<<endl;
     //     std::cout<<"Temp Vector Channel "<<TempVector.at(i).GetChannelID()<<endl;
    for (int j = 0; j < (int)TempVector.size(); j++){
      if (TempVector.at(j).GetTime()*1000. >= MinTimeWindow && TempVector.at(j).GetTime()*1000. <= MaxTimeWindow){
        if (TempVector.at(j).GetPeak() > MaxAmp){
          MaxAmp = TempVector.at(j).GetPeak();
          MaxPulse = LAPPDPulseCluster.at(i).at(j);
        }
      }
    }
  }
  std::cout<<"LargestPulse found"<<endl;
  int MPC=MaxPulse.GetChannelID();

//ParallelPosition and time
LAPPDPulse OpposPulse;
  if (TwoSided) {

    // what is the strip ID of the max channel?
    int cid = MaxPulse.GetChannelID();
    int sid = stripID.at(cid);

    bool hasopposite=false;
    // look to see if we have an opposite strip...if so, we grab it
    std::vector<LAPPDPulse> negaVector;
    for(int k=0; k<3; k++){

      if(sid == -stripID.at(k)){
        hasopposite=true;
        negaVector = LAPPDPulseCluster.at(k);
      }
    }

    // is there a pulse in the time window, on the opposite side of the strip? (I hope so)
    LAPPDPulse OpposPulse;
    for (int j = 0; j < (int)negaVector.size(); j++){
    //Find pulse inside this vector that has the same peak and thetime
      if (fabs(negaVector.at(j).GetTime()*1000. - MaxPulse.GetTime()*1000.) < PTRange && (fabs(negaVector.at(j).GetPeak() - MaxAmp) / MaxAmp) <= 0.1){
        OpposPulse = negaVector.at(j);
        Deltatime = (MaxPulse.GetTime()*1000. - OpposPulse.GetTime()*1000.);
      }
      ParaPosition = Deltatime*c;
      LocalPosition.push_back(ParaPosition);
      HitTime = (MaxPulse.GetTime()*1000. + OpposPulse.GetTime()*1000.)/2;
    }
  }
  else {
    ParaPosition = 0;
    HitTime = MaxPulse.GetTime()*1000.;
    LocalPosition.push_back(ParaPosition);
  }


//Perpendicular position
  double SumAbove=0;
  double SumBelow=0;

  if(!TwoSided){
  std::cout<<"Not twosided"<<endl;
    for (int i =0;i<3; i++){
      std::vector<LAPPDPulse> neighbourVector;
       if(i!=MaxPulse.GetChannelID()){

      neighbourVector = LAPPDPulseCluster.at(i);

      for (int j = 0; j < (int)neighbourVector.size(); j++){
        if (fabs(neighbourVector.at(j).GetTime()*1000. - MaxPulse.GetTime()*1000.) < PTRange){

	   NeighboursPulses.push_back(neighbourVector.at(j));
        }
	std::cout<<"NeighboursPulses size  "<<NeighboursPulses.size()<<endl;
      }
    }
}
    if(NeighboursPulses.size()>0){
  for (int i=0; i<3; i++){
    int j=0;
    if(j<(int)NeighboursPulses.size()){
      if(i!=MaxPulse.GetChannelID())
	{
	SumAbove +=  (stripID[i] * NeighboursPulses.at(j).GetPeak());
	SumBelow +=  (NeighboursPulses.at(j).GetPeak());
	j++;
	}
      else
	{
	SumAbove += (stripID[i] * MaxPulse.GetPeak());
	SumBelow += MaxPulse.GetPeak();
	}
    }
  }
  PerpPosition = (SumAbove / SumBelow);
    }
    else {
      PerpPosition =stripID.at(MaxPulse.GetChannelID());
    }
  }
  else
    {
      for(int i =0; i<3; i++){
	if(i!=MaxPulse.GetChannelID() && i!=OpposPulse.GetChannelID()){
	  for (int j = 0; j < (int)LAPPDPulseCluster.at(i).size(); j++){
	    if (fabs(LAPPDPulseCluster.at(i).at(j).GetTime()*1000. - MaxPulse.GetTime()*1000.) < PTRange){
	       NeighboursPulses.push_back(LAPPDPulseCluster.at(i).at(j));
	     }
	   }
	}
      }
      PerpPosition=stripID.at(MaxPulse.GetChannelID());
	}
 LocalPosition.push_back(PerpPosition);

 std::cout<<stripID.at(0)<<endl;
 std::cout<<stripID.at(1)<<endl;
 std::cout<<stripID.at(2)<<endl;

  //Store data in LAPPDHit
  LAPPDHit kHit(MaxPulse.GetChannelID(), MaxPulse.GetTime(), MaxPulse.GetPeak(), AbsPosition, LocalPosition);
  std::cout<<"Hit created "<<endl;
  int PulseNum;
  m_data->Stores["ANNIEEvent"]->Set("RecoLaserTestHit",kHit);
  std::cout<<"Hit saved  "<<endl;
  m_data->Stores["ANNIEEvent"]->Set("isTwoSided",TwoSided);
  std::cout<<"TwoSided Saved  "<<endl;
  if(TwoSided){
    if(NeighboursPulses.size()>0){
  HitPulses.push_back(NeighboursPulses.at(0));
  HitPulses.push_back(MaxPulse);
  HitPulses.push_back(OpposPulse);
  PulseNum=3;
  std::cout<<"Stored twosided"<<endl;
    }
    else {
      LAPPDPulse FillerPulse;
      FillerPulse.SetChannelID(7);
      HitPulses.push_back(FillerPulse);
      HitPulses.push_back(MaxPulse);
      HitPulses.push_back(OpposPulse);
      PulseNum=2;
        std::cout<<"Stored twosided +filler"<<endl;
}
  }
else
{
  if(NeighboursPulses.size()>1){
  std::cout<<NeighboursPulses.size()<<std::endl;
  HitPulses.push_back(NeighboursPulses.at(0));
  HitPulses.push_back(MaxPulse);
  HitPulses.push_back(NeighboursPulses.at(1));
  PulseNum=3;
  std::cout<<"Stored Onesided"<<endl;
  }
  else if(NeighboursPulses.size()>0){
      LAPPDPulse FillerPulse;
      FillerPulse.SetChannelID(7);
      HitPulses.push_back(NeighboursPulses.at(0));
      HitPulses.push_back(MaxPulse);
      HitPulses.push_back(FillerPulse);
      PulseNum=2;
    }
    else
    {
      LAPPDPulse FillerPulse;
      FillerPulse.SetChannelID(7);
      HitPulses.push_back(FillerPulse);
      HitPulses.push_back(MaxPulse);
      HitPulses.push_back(FillerPulse);
      PulseNum=1;
  std::cout<<"Stored Onesided + filler"<<endl;
}

  m_data->Stores["ANNIEEvent"]->Set("HitPulses",HitPulses);
 }
  std::cout<<"first pulse "<<HitPulses.at(0).GetChannelID()<<endl;
  std::cout<<"second pulse "<<HitPulses.at(1).GetChannelID()<<endl;
  std::cout<<"third pulse "<<HitPulses.at(2).GetChannelID()<<endl;
  m_data->Stores["ANNIEEvent"]->Set("PulseNum",PulseNum);
  return true;
}


bool LAPPDlasertestHitFinder::Finalise(){

  return true;
}
