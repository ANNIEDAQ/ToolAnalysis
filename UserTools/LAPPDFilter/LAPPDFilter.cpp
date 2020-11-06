#include "LAPPDFilter.h"

LAPPDFilter::LAPPDFilter():Tool(){}


bool LAPPDFilter::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  bool isFiltered = true;
  m_data->Stores["ANNIEEvent"]->Header->Set("isFiltered",isFiltered);

  TString FIWL;
  //FilterInputWavLabel;
  m_variables.Get("FilterInputWavLabel",FIWL);
  FilterInputWavLabel = FIWL;
  m_variables.Get("Nsamples", DimSize);
  m_variables.Get("CutoffFrequency", CutoffFrequency);
  m_variables.Get("SampleSize",Deltat);
  return true;
}


bool LAPPDFilter::Execute(){

  Waveform<double> bwav;

  // get raw lappd data
  std::map<int,vector<Waveform<double>>> rawlappddata;

  //m_data->Stores["ANNIEEvent"]->Get("RawLAPPDData",rawlappddata);
  m_data->Stores["ANNIEEvent"]->Get(FilterInputWavLabel,rawlappddata);

  // the filtered Waveform
  std::map<int,vector<Waveform<double>>> filteredlappddata;

  map <int, vector<Waveform<double>>> :: iterator itr;
  for (itr = rawlappddata.begin(); itr != rawlappddata.end(); ++itr){
    int channelno = itr->first;
    vector<Waveform<double>> Vwavs = itr->second;
    vector<Waveform<double>> Vfwavs;

    //loop over all Waveforms
    for(int i=0; i<(int)Vwavs.size(); i++){

        Waveform<double> bwav = Vwavs.at(i);
        Waveform<double> filtwav = Waveform_FFT(bwav);
        Vfwavs.push_back(filtwav);
      }

      filteredlappddata.insert(pair <int,vector<Waveform<double>>> (channelno,Vfwavs));
    }

  m_data->Stores["ANNIEEvent"]->Set("FiltLAPPDData",filteredlappddata);


  return true;
}


bool LAPPDFilter::Finalise(){

  return true;
}



Waveform<double> LAPPDFilter::Waveform_FFT(Waveform<double> iwav) {

  // Hack !!!!! Converting Waveform to Hist
  int nbins = iwav.GetSamples()->size();
  double starttime=0.;
  double endtime = starttime + ((double)nbins)*100.;
  TH1D* hwav = new TH1D("hwav","hwav",nbins,starttime,endtime);

  for(int i=0; i<nbins; i++){
    hwav->SetBinContent(i+1,iwav.GetSample(i));
  }

  //Look at the real part of the output
	TH1 *hr = 0;
	hr = hwav->FFT(hr, "RE");
	hr->SetName("hr");
	hr->SetTitle("Real part of the 1st transform");
	hr->GetXaxis()->Set(DimSize,0,1e10); //Hz

	//Look at the imaginary part of the output
	TH1 *him = 0;
	him = hwav->FFT(him, "IM");
	him->SetName("him");
	him->SetTitle("Imaginary part of the 1st transform");
	him->GetXaxis()->Set(DimSize,0,1e10);

	//Look at the DC component and the Nyquist harmonic:
	double re, im;
	//That's the way to get the current transform object:
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	//Use the following method to get the full output:
	double *re_full = new double[DimSize];
	double *im_full = new double[DimSize];
	fft->GetPointsComplex(re_full,im_full);
	//filter
	TH1D* freqq = new TH1D("frequency","frequency",DimSize,0,DimSize);
	TH1D* freqq2 = new TH1D("frequency2","frequency2",DimSize,0,DimSize);
	for(int i=0;i<DimSize;i++) {
	  //		double f = Waveform_Filter1(5e8, i*1.0e10/DimSize);
	  double f = Waveform_Filter2(CutoffFrequency, 4, i*1.0e10/DimSize);

    re_full[i] = re_full[i]*f;
    im_full[i] = im_full[i]*f;
	  freqq->SetBinContent(i+1,re_full[i]);
	}

	//Now let's make a backward transform:
  Waveform<double> owav;
	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &DimSize, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();
	TH1 *hb = 0;
	//Let's look at the output
	hb = TH1::TransformHisto(fft_back,hb,"Re");
	hb->SetTitle("The backward transform result");
	hb->Scale(1.0/DimSize);
	hb->GetXaxis()->Set(DimSize,0,DimSize*Deltat); //ps
	for(int i=0;i<DimSize;i++) {
    owav.PushSample(hb->GetBinContent(i+1));
	  //	im_fft.push_back(him->GetBinContent(i+1));
	  //	vol_fft[i] = hb->GetBinContent(i+1);
	  //	re_fft[i] = hr->GetBinContent(i+1);
	  //	im_fft[i] = him->GetBinContent(i+1);
	  //  freqq2->SetBinContent(i+1,abs(re_fft[i]));
	}

  delete freqq;
  delete freqq2;
  delete hwav;
	delete hr;
  delete him;
	delete hb;
	delete re_full;
	delete im_full;
	delete fft;
	delete fft_back;
	hr=0;
  him=0;
	hb=0;
	re_full=0;
	im_full=0;
	fft=0;
	fft_back=0;

  return owav;
}

float LAPPDFilter::Waveform_Filter2(float CutoffFrequency, int T, float finput) {
	double f = 1.0/(1+TMath::Power(finput/CutoffFrequency, 2*T));
	return f;
}
