#include "RawLoadToRoot.h"

RawLoadToRoot::RawLoadToRoot():Tool(){}


bool RawLoadToRoot::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("treeoutfile", treeoutfile);
  m_variables.Get("lowerbound", lbound);
  m_variables.Get("upperbound", ubound);
  m_variables.Get("onoffswitch", onoffswitch);
  m_variables.Get("onoffswitchII", onoffswitchII);
  m_variables.Get("lowerboundII", lboundII);
  m_variables.Get("upperboundII", uboundII);
  m_variables.Get("timewindow",ttot);

  treeoutput = new TFile(treeoutfile, "RECREATE");
  //treeoutput->mkdir("fullwaveforms");
  //treeoutput->mkdir("pulses");
  //treeoutput->mkdir("starttime");
  //treeoutput->mkdir("bsubwaveforms");
  treeoutput->mkdir("pulsewaveforms");
  //treeoutput->mkdir("movingbaseline");

  rawtotree = new TTree("RawLoadTree", "RawLoadTree");

  rawtotree->Branch("EventNumber", &EventNumber);
  rawtotree->Branch("RunNumber", &RunNumber);
  rawtotree->Branch("SubRunNumber", &SubRunNumber);

  AllStartTimesHist = new TH1D("AllStarts","AllStarts",(ttot/500),0,ttot);
  ChanPulseHist = new TH1D("ChanPulses", "ChanPulses",66,0,65);

  RelEventNumber = 0;

  return true;
}


bool RawLoadToRoot::Execute(){

  auto* annieevent = m_data->Stores["ANNIEEvent"];

  annieevent->Get("EventNumber", EventNumber);
  annieevent->Get("RunNumber", RunNumber);
  annieevent->Get("SubRunNumber", SubRunNumber);
  annieevent->Get("RawADCData", RawADCData);
  annieevent->Get("TrigEvents",trigev);
  annieevent->Get("CalibratedADCData", caladcdata);

  //restricts tool to range of trigger events specified in config file
  if((lbound<=RelEventNumber && RelEventNumber<=ubound && onoffswitch==1)||onoffswitch==0){

    map<unsigned long,vector<Waveform<unsigned short>>> :: iterator itr;
    map<unsigned long,std::vector<CalibratedADCWaveform<double>>> :: iterator ijk;

    int chancount=0;
    for(itr=RawADCData.begin(),ijk=caladcdata.begin(); itr!=RawADCData.end(),ijk!=caladcdata.end(); ++itr,++ijk){  //loop through channels

      unsigned long ck = itr->first;
      vector<Waveform<unsigned short>> TheWaveforms = itr->second;
      unsigned long chank = ijk->first;
      vector<CalibratedADCWaveform<double>> calwaves = ijk->second;

      for(int mmm=0; mmm<(int)TheWaveforms.size(); mmm++){  //loop through minibuffers

        CalibratedADCWaveform<double> acalwave = calwaves.at(mmm);
        double basline = acalwave.GetBaseline();

        int pulsize = trigev.at(mmm).at(ck).size();  // TODO check if (ck-1) was still needed!!

        //TString wname;
        //wname+="Waveform_CH";
        //wname+=chancount;
        //wname+="_Minibuffer";
        //wname+=mmm;
        //wname+="_Read";
        //wname+=EventNumber;
        //TString bwname;
        //bwname+="BaseSubWaveform_CH";
        //bwname+=chancount;
        //bwname+="_Minibuffer";
        //bwname+=mmm;
        //bwname+="_Read";
        //bwname+=EventNumber;
        TString puname;
        puname+="PulseWaveform_CH";
        puname+=chancount;
        puname+="_Minibuffer";
        puname+=mmm;
        puname+="_Read";
        puname+=EventNumber;
        //theWaveformHist = new TH1D(wname,wname,(ttot/2),0,ttot);
        //BSubWaveformHist = new TH1D(bwname,bwname,(ttot/2),0,ttot);
        Waveform<unsigned short> aWaveform = TheWaveforms.at(mmm);
        //Excludes TDC cards and Channel 57 which is faulty
        if(chancount<60&&chancount!=57){
        //for(int nnn=0; nnn<aWaveform.Samples().size(); nnn++){
          //theWaveformHist->SetBinContent(nnn+1,(double)aWaveform.GetSample(nnn));
          //BSubWaveformHist->SetBinContent(nnn+1,((double)aWaveform.GetSample(nnn))-basline);
        //}

        double bassline = 0.;
        if(pulsize>0){
          thePulseformHist = new TH1D(puname,puname,(ttot/2),0,ttot);
          for(int vvv=0; vvv<(int)aWaveform.Samples().size(); vvv++){
            //115-206 establish the dynamic baseline
            if(vvv==0){
              bassline=aWaveform.GetSample(vvv);
            }

            if(vvv==1
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1))/2;
            }

            if(vvv==2
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2))/3;
            }

            if(vvv==3
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3))/4;
            }

            if(vvv==4
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4))/5;
            }

            if(vvv==5
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-5)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4)+aWaveform.GetSample(vvv-5))/6;
            }

            if(vvv==6
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-5)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-6)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4)+aWaveform.GetSample(vvv-5)+aWaveform.GetSample(vvv-6))/7;
            }

            if(vvv==7
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-5)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-6)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-7)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4)+aWaveform.GetSample(vvv-5)+aWaveform.GetSample(vvv-6)+aWaveform.GetSample(vvv-7))/8;
            }

            if(vvv==8
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-5)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-6)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-7)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-8)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4)+aWaveform.GetSample(vvv-5)+aWaveform.GetSample(vvv-6)+aWaveform.GetSample(vvv-7)+aWaveform.GetSample(vvv-8))/9;
            }

            if(vvv>=9
            &&fabs(aWaveform.GetSample(vvv)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-1)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-2)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-3)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-4)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-5)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-6)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-7)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-8)-bassline)<5
            &&fabs(aWaveform.GetSample(vvv-9)-bassline)<5){
              bassline=(aWaveform.GetSample(vvv)+aWaveform.GetSample(vvv-1)+aWaveform.GetSample(vvv-2)+aWaveform.GetSample(vvv-3)+aWaveform.GetSample(vvv-4)+aWaveform.GetSample(vvv-5)+aWaveform.GetSample(vvv-6)+aWaveform.GetSample(vvv-7)+aWaveform.GetSample(vvv-8)+aWaveform.GetSample(vvv-9))/10;
            }

            thePulseformHist->SetBinContent(vvv+1,(double)aWaveform.GetSample(vvv)-bassline);
          }
          treeoutput->cd("pulsewaveforms");
          thePulseformHist->Write();
          delete thePulseformHist;
        }
      }
        //treeoutput->cd("fullwaveforms");
        //theWaveformHist->Write();
        //    mycounter++;
        //delete theWaveformHist;
        //treeoutput->cd("bsubwaveforms");
        //BSubWaveformHist->Write();
        //delete BSubWaveformHist;
      }
      chancount++;

    }

    rawtotree->Fill();

  }else{

  }

  RelEventNumber++;

  if((lboundII<=EventNumber && EventNumber<=uboundII && onoffswitchII==1)||onoffswitchII==0){

    //TString aname;
    //aname += "StartTime_Read";
    //aname += EventNumber;

    //TH1D* StartTimeHist = new TH1D(aname,aname,(ttot/20),0,ttot);

    map<int,map<int,std::vector<ADCPulse>>> :: iterator iter;
    for(iter=trigev.begin(); iter!=trigev.end(); ++iter){

      int minibuffernum = iter->first;
      map<int,std::vector<ADCPulse>> chanpulses = iter->second;
      map<int,std::vector<ADCPulse>> :: iterator itera;
      for(itera=chanpulses.begin(); itera!=chanpulses.end(); ++itera){

        int channelnum = itera ->first;
        std::vector<ADCPulse> somepulses = itera->second;

        if(channelnum<59&&channelnum!=57){
        ADCPulse prevpulse;
          for(int ccc=0; ccc<(int)somepulses.size(); ccc++){

            ADCPulse apulse = somepulses.at(ccc);
            if(ccc>0){prevpulse = somepulses.at(ccc-1);}

            //StartTimeHist->Fill(apulse.start_time());
            AllStartTimesHist->Fill(apulse.start_time());
            ChanPulseHist->Fill(channelnum);
          }
        }

      }

    }

    //treeoutput->cd("starttime");
    //StartTimeHist->Write();
    //delete StartTimeHist;

  }

  return true;
}

bool RawLoadToRoot::Finalise(){

  treeoutput->cd();
  rawtotree->Write();
  AllStartTimesHist->Write();
  ChanPulseHist->Write();
  treeoutput->Close();

  return true;
}
