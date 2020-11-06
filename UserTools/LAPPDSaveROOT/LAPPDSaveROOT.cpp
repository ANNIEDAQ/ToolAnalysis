#include "LAPPDSaveROOT.h"

LAPPDSaveROOT::LAPPDSaveROOT():Tool(){}


bool LAPPDSaveROOT::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // setup the output files
  TString OutFile = "lappdout.root";	 //default input file
  m_variables.Get("outfile", OutFile);
  tf = new TFile(OutFile,"RECREATE");
  tf->mkdir("wavs");
  tf->mkdir("filteredwavs");
  tf->mkdir("blswavs");

  isFiltered=false;
  isBLsub=false;

  // have the pulses been filtered? integrated? is this real or simulated data?
  bool ifilt = m_data->Stores["ANNIEEvent"]->Header->Get("isFiltered",isFiltered);
  bool ibl = m_data->Stores["ANNIEEvent"]->Header->Get("isBLsubtracted",isBLsub);
  bool iI = m_data->Stores["ANNIEEvent"]->Header->Get("isIntegrated",isIntegrated);
  bool isim = m_data->Stores["ANNIEEvent"]->Header->Get("isSim",isSim);
  cout<<"PP "<<ifilt<<" "<<ibl<<" "<<iI<<" "<<isim<<endl;
  // from the config file, get the number of channels, the trigger channel, num histos, samplesize
  m_variables.Get("NChannels", NChannel);
  m_variables.Get("TrigChannel", TrigChannel);
  m_variables.Get("NHistos", NHistos);
  m_variables.Get("SampleSize",Deltat);

  // keep count of the loop number (starting from 0)
  miter=0;

  // set up the branches of the ROOT tree
  outtree = new TTree("LAPPDdata","LAPPDdata");
  outtree->Branch("iteration",&miter);
  outtree->Branch("channel",&chno);
  outtree->Branch("cfdtime",&cfdtime);
  outtree->Branch("amplitude",&amp);
  outtree->Branch("Twidth",&twidth);


LAPPDTree = new TTree("ISULAPPDData","ISULAPPDData");
LAPPDTree->Branch("Channel1Charge",&chrgCH1);
LAPPDTree->Branch("Channel2Charge",&chrgCH2);
LAPPDTree->Branch("Channel3Charge",&chrgCH3);
LAPPDTree->Branch("Channel1Time",&tpsecCH1);
LAPPDTree->Branch("Channel2Time",&tpsecCH2);
LAPPDTree->Branch("Channel3Time",&tpsecCH3);
LAPPDTree->Branch("Channel1Amp",&ampCH1);
LAPPDTree->Branch("Channel2Amp",&ampCH2);
LAPPDTree->Branch("Channel3Amp",&ampCH3);
LAPPDTree->Branch("ParallelPosition",&ParaPos);
LAPPDTree->Branch("Transverseposition",&TransPos);
LAPPDTree->Branch("PulseNum",&PulseNum);
  // declare the histograms
  hAmp = new TH1D*[NChannel];
  hTime = new TH1D*[NChannel];
  hQ = new TH1D*[NChannel];

  // initialize the histograms
  for(int i=0; i<NChannel; i++){
    TString AmpName;
    AmpName+="Amplitudes_CH";
    AmpName+=i;
    hAmp[i] = new TH1D(AmpName,AmpName,1000,0.,50.);

    TString QName;
    QName+="Charge_CH";
    QName+=i;
    hQ[i] = new TH1D(QName,QName,8800,-1e7,10e7);

    TString TimeName;
    TimeName+="Time_CH";
    TimeName+=i;
    hTime[i] = new TH1D(TimeName,TimeName,10000,0.,100000.);
  }

  return true;
}


  bool LAPPDSaveROOT::Execute(){
  vector<LAPPDPulse> LAPPDHitPulses;
  LAPPDHit ISUHits;
  m_data->Stores["ANNIEEvent"]->Get("RecoLaserTestHit",ISUHits);
  m_data->Stores["ANNIEEvent"]->Get("HitPulses",LAPPDHitPulses);
  m_data->Stores["ANNIEEvent"]->Get("PulseNum",PulseNum);
  std::cout<<LAPPDHitPulses.size()<<endl;
  for(int i=0; i<(int)LAPPDHitPulses.size(); i++){
    LAPPDPulse tempPulse = LAPPDHitPulses.at(i);
    std::cout<<tempPulse.GetChannelID()<<endl;
    if(tempPulse.GetChannelID()==0)
    {
      ampCH1=LAPPDHitPulses.at(i).GetPeak();
      tpsecCH1=LAPPDHitPulses.at(i).GetTime()*1000.;
      chrgCH1=LAPPDHitPulses.at(i).GetCharge();
       std::cout<<"Amp 1 "<<ampCH1<<endl;
       std::cout<<"Tpsec 1 "<<tpsecCH1<<endl;
       std::cout<<"Charge 1 "<<chrgCH1<<endl;
    }
     if(tempPulse.GetChannelID()==1)
    {
      ampCH2=LAPPDHitPulses.at(i).GetPeak();
      tpsecCH2=LAPPDHitPulses.at(i).GetTime()*1000.;
      chrgCH2=LAPPDHitPulses.at(i).GetCharge();
      std::cout<<"Amp 2 "<<ampCH2<<endl;
 std::cout<<"Tpsec 2 " <<tpsecCH2<<endl;
 std::cout<<"Charge 2 "<<chrgCH2<<endl;
    }
    if(tempPulse.GetChannelID()==2)
    {
      ampCH3=LAPPDHitPulses.at(i).GetPeak();
      tpsecCH3=LAPPDHitPulses.at(i).GetTime()*1000.;
      chrgCH3=LAPPDHitPulses.at(i).GetCharge();
      std::cout<<"Amp 3 "<<ampCH3<<endl;
 std::cout<<"Tpsec 3 "<<tpsecCH3<<endl;
 std::cout<<"Charge 3 "<<chrgCH3<<endl;
    }
    
    
    ParaPos=ISUHits.GetLocalPosition().at(0);
    TransPos=ISUHits.GetLocalPosition().at(1);
    std::cout<<"ParaPos "<<ParaPos<<endl;
    std::cout<<"TransPos "<<TransPos<<endl;

}
 LAPPDTree->Fill(); 
      // get raw lappd data
  std::map<int,vector<Waveform<double>>> rawlappddata;
  m_data->Stores["ANNIEEvent"]->Get("RawLAPPDData",rawlappddata);
  // get filtered data
  std::map<int,vector<Waveform<double>>> filteredlappddata;
  if(isFiltered) m_data->Stores["ANNIEEvent"]->Get("FiltLAPPDData",filteredlappddata);
  std::map<int,vector<Waveform<double>>> BLsubtractedlappddata;
  if(isBLsub)    m_data->Stores["ANNIEEvent"]->Get("BLsubtractedLAPPDData",BLsubtractedlappddata);
  // get charge information
  std::map<int, vector<double>> TheCharges;
  if(isIntegrated) m_data->Stores["ANNIEEvent"]->Get("theCharges",TheCharges);
  // get reconstructed LAPPD pulses
  std::map<int,vector<LAPPDPulse>> CFDRecoLAPPDPulses;
  m_data->Stores["ANNIEEvent"]->Get("CFDRecoLAPPDPulses",CFDRecoLAPPDPulses);

  // loop over all channels
  std::map<int, vector<Waveform<double>>> :: iterator itr;
  for (itr = rawlappddata.begin(); itr != rawlappddata.end(); ++itr){

    int channelno = itr->first;
    vector<Waveform<double>> Vwavs = itr->second;

    vector<Waveform<double>> Vfwavs;
    if(isFiltered){
      std::map<int, vector<Waveform<double>>>::iterator fp;
      fp = filteredlappddata.find(channelno);
      Vfwavs = fp->second;
    }

    vector<Waveform<double>> Vblswavs;
    if(isBLsub){
      std::map<int, vector<Waveform<double>>>::iterator bp;
      bp = BLsubtractedlappddata.find(channelno);
      Vblswavs = bp->second;
    }

    // get the pulses for each channel
    map<int, vector<LAPPDPulse>>::iterator p;
    p = CFDRecoLAPPDPulses.find(channelno);
    vector<LAPPDPulse> Vpulses = p->second;

    // get the vector of charges for each channel
    map<int, vector<double>>::iterator pq;
    pq = TheCharges.find(channelno);
    vector<double> Vqs = pq->second;

    // loop over the charges and store to a hist
    for(int pqi=0; pqi<(int)Vqs.size(); pqi++){
      double theq = Vqs.at(pqi);
      if(channelno>=0&&channelno<NChannel){
        if(channelno!=TrigChannel) hQ[channelno]->Fill(theq);
      } else {cout<<"TOO MANY CHANNELS"<<endl;}
    }

    //loop over all pulses fill histos and fill tree with pulse information
    for(int i=0; i<(int)Vpulses.size(); i++){
      LAPPDPulse thepulse = Vpulses.at(i);
      if(channelno!=TrigChannel){

        if(channelno>=0&&channelno<NChannel){
            hAmp[channelno]->Fill(thepulse.GetPeak());
            hTime[channelno]->Fill(thepulse.GetTime()*1000.);
        }

        chno=channelno;
        cfdtime=thepulse.GetTime()*1000.;
        amp=thepulse.GetPeak();
        twidth=(thepulse.GetHiRange()-thepulse.GetLowRange())*Deltat;
        outtree->Fill();
      }

    }

    //loop over all Waveforms, make histrograms
    std::vector<LAPPDPulse> thepulses;
    for(int i=0; i<(int)Vwavs.size(); i++){

        Waveform<double> bwav = Vwavs.at(i);
        Waveform<double> bfwav;
        if(isFiltered) bfwav = Vfwavs.at(i);
        Waveform<double> blswav;
        if(isBLsub) blswav = Vblswavs.at(i);

        int nbins = bwav.GetSamples()->size();
        double starttime=0.;
        double endtime = starttime + ((double)nbins)*100.;

        TString hname;
        hname+="wav_ch";
        hname+=channelno;
        hname+="_wav";
        hname+=i;
        hname+="_evt";
        hname+=miter;

        TH1D* hwav = new TH1D(hname,hname,nbins,starttime,endtime);
        TH1D* fwav = new TH1D("filtered"+hname,"filtered"+hname,nbins,starttime,endtime);
        TH1D* blwav = new TH1D("blssub"+hname,"blssub"+hname,nbins,starttime,endtime);

        for(int i=0; i<nbins; i++){
          hwav->SetBinContent(i+1,-bwav.GetSamples()->at(i));
          if(isFiltered) fwav->SetBinContent(i+1,-bfwav.GetSamples()->at(i));
          if(isBLsub)  blwav->SetBinContent(i+1,-blswav.GetSamples()->at(i));

        }

        // sotre wavform hists in the "wav" folder
        tf->cd("wavs");
        if(miter<NHistos) hwav->Write();

        if(isFiltered){
          tf->cd("filteredwavs");
          if(miter<NHistos) fwav->Write();
        }

        if(isBLsub){
          tf->cd("blswavs");
          if(miter<NHistos) blwav->Write();
        }

        delete blwav;
        delete hwav;
        delete fwav;
    }
  }

  miter++;
  return true;
}


bool LAPPDSaveROOT::Finalise(){

  // go to the top level of the output file
  tf->cd();

  // write the summary histos to file
  for(int i=0; i<NChannel; i++){
    hAmp[i]->Write();
    hTime[i]->Write();
    hQ[i]->Write();
  }
  // write the output tree to file
  outtree->Write();
  LAPPDTree->Write();
  tf->Close();

  return true;
}
