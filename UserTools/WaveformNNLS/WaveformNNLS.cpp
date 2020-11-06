/*
This tool is a waveform processing tool that uses the Lawson and Hanson 
non-negative least squares deconvolution algorithm to find signals 
with a known "template" in a noisy waveform. The libraries found
in this tool directory are taken from: http://suvrit.de/work/soft/nnls.html 

The algorithm numerically solves the matrix equation Ax = b by
minimizing 0.5*|Ax - b|^2 . In the context of a template fit, 
the matrix A represents information about the template, b is 
a vector form of the noisy waveform where the vector index is a time
index, and x is a "solution vector". Each row of A is a time series
(indexed by the column) that is a template waveform. Subsequent rows of
A have the template waveform shifted in time by one sample. Therefore, each 
element of solution vector "x" represents how many template waveforms should
be placed at a particular index in time. The sum of all such components is 
the final fitted solution, A*x . 

Expects raw data, named in the config file, of the form map<int, vector<Waveform<double>>>.
Expects a template that is a TH1D coming from a root file. The template often has a much
higher sampling rate than the raw data. This algorithm can handle that and can use it as
extra information. There is a parameter "samplingFactor" that controls how to up/downsample
the template and raw data so that they match timescales. A list of sample times is also
needed to place the solution components in the correct position relative to the rawdata
time-series. 

-Evan Angelico, 7/11/2019
*/


#include "WaveformNNLS.h"
#include <TCanvas.h>
#include <math.h>
#include <string>



WaveformNNLS::WaveformNNLS():Tool(){}


bool WaveformNNLS::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool WaveformNNLS::Execute(){

	

	//rawdata loaded in format of LAPPDSim
	string RawDataName;
	m_variables.Get("rawdataname", RawDataName);

	map<int,vector<Waveform<double>>> rawData;
	m_data->Stores["ANNIEEvent"]->Get(RawDataName,rawData);

	//This variable controls the sampling depth of the
	//nnls algorithm and allows the template to have
	//more samples than the signal waveform. The template
	//waveform is reduced to (nsamples template)/samplingFactor
	//number of samples. The timestep between template
	//and signal waveform is made equal by expanding the
	//signal waveform to a larger number of total samples. 
	double samplingFactor;
	m_variables.Get("sampling_factor", samplingFactor);

	//you can use any template you want for the 
	//algorithm. It assumes that it is a histogram
	//in a root file named ""
	TString tempfilename = "dummy.root";
	m_variables.Get("tempfilename", tempfilename);
	TString temphistname = "dummy";
	m_variables.Get("temphistname", temphistname);

	//get the template waveform. This assumes
	//that the file is .root file with a TH1D
	//named templatehistname
	TFile* tempfile = new TFile(tempfilename, "READ");
	TH1D* temphist = (TH1D*)tempfile->Get(temphistname);
	Waveform<double> tempwave; //signal of template
	vector<double> temptimes;  //times of template
	int nbins = temphist->GetNbinsX(); 
	double binwidth = temphist->GetBinWidth(0); //should be in ns
	double starttime = temphist->GetBinLowEdge(0); //first time in template histogram, used to set template to start at 0ps
	double endtime = temphist->GetBinLowEdge(nbins-1);

	//this new timestep will be the new timestep for the 
	//raw data. Even if the raw data is sampled at a non-constant
	//rate (for example, well calibrated LAPPD electronics), the raw
	//waveform is re-sampled (not interpolated) at this new timestep
	double newtimestep = binwidth*samplingFactor; //in ps

	// check units and such
	/*
	cout << "nbins = " << nbins << endl;
	cout << "binwidth = " << binwidth << endl;
	cout << "start = " << starttime << endl;
	cout << "endtime = " << endtime << endl;
	cout << "newtimestep = " << newtimestep << endl;
	*/
	

	//fill the template waveform and time vectors 
	//at a downsampled/upsampled time sampling using
	//linear interpolation. Interpolation not necessary, it
	//can just be re-sampled. One method may perform better than another
	for(double t = starttime; t <= endtime; t+=newtimestep)
	{
		tempwave.PushSample(temphist->Interpolate(t)); 

		//zero the time vector, in ps just like signal wave
		temptimes.push_back((t - starttime)); 
	}


	//It is necessary to make sure that the template
	//and the rawdata are on the same time-scale. 
	//What follows is code to set the timescale of the
	//raw data. It later is adjusted/resampled to match
	//the template times. 
	vector<float> sampletimes;
    
	//different data types will have different sample
	//times. For example, LAPPDs will eventually have
	//non-constant sampling rates and times associated with
	//each sample of each channel. 
	//Possible TODO: include a times vector for your desired
	//datatype into the ANNIE store?. For now, do logic on
	//the string that indexes the rawdata type
	if(RawDataName == "RawLAPPDData")
	{
		//later, we need to include real sample time data format
	    //based on non-zero aperture jitter of boards on an lappd/channel
	    //by channel basis. For now, just assume constant sampling rate
	    double dt = (1.0/(256*40*1e6))*1e12; //ps
	    for(int i = 0; i < 256; i++)
	    {
	    	sampletimes.push_back(i*dt);
	    }
	}
	//default to avoid crashing
	else
	{
		Waveform<double> example_waveform = rawData[0].front();
		for(int i = 0; i < (int) example_waveform.GetSamples()->size(); i++)
		{
			sampletimes.push_back(i);
		}
	}

	//create a new signal times list based on the 
	//new timestep that was determined from the template.
	//this is really only used to (1) return the time at which
	//a solved nnls component waveform should appear in the rawdata
	//and (2) to count the number of rows and columns in the nnls matrix
	vector<double> newsignaltimes;
	for(double t = sampletimes.front(); t <= sampletimes.back(); t+=newtimestep)
	{
		newsignaltimes.push_back(t);
	}

	//the number of rows for the matrix
	//and the number of samples in the signal waveform
	//is set to match the timesteps of the NEW template
	//waveform and the signal waveform. Here, we get
	//the signal waveform and expand it to that nrows sampling rate
	map <int, vector<Waveform<double>>> :: iterator itr;
    Waveform<double> signalwave; 
    size_t nrows = newsignaltimes.size();
    int flag; //error flag on nnls solver

    cout << "doing a " << nrows << " x " << nrows << " = " << nrows*nrows << " matrix " << endl;

    //***Assumes all waveforms are same number of samples
    //a number of things in this tool would change if that
    //were not the case

    //initialize the main algebra variables
    //for the nnls algo

    //maximum number of nnls iterations
	int maxiter;
	m_variables.Get("maxiter", maxiter);

    nnlsvector* b = new nnlsvector(nrows);
    nnlsvector* x = new nnlsvector(nrows);
    nnlsmatrix* A = new denseMatrix(nrows, nrows);
    BuildTemplateMatrix(A, tempwave, nrows);

    //The solution to the NNLS algorithm is 
    //stored in a class defined in the DataModel. 
    //Check out that class for more info.
    map<int, NnlsSolution> soln;


    for(itr = rawData.begin(); itr != rawData.end(); ++itr)
    {
    	int ch = itr->first; //channel of this waveform in the loop
    	cout << "On channel " << ch << endl;

    	//in the solution object, save the template
    	//waveform for each channel in the map
    	soln[ch].SetTemplate(tempwave, temptimes);


    	//I've found that LAPPD has only first element
    	//of the vector<Waveform<double>> component of the rawData
    	//populated. Is this true with PMTs? If not, you will need 
    	//another loop here. 
    	signalwave = itr->second.front(); 
    	
		//pass the signal vector pointer that gets
		//modified in the following function
		BuildWaveformVector(b, signalwave, sampletimes, newtimestep);
    	
    	nnls* solver = new nnls(A, b, maxiter);
    	flag = solver->optimize();
    	if(flag < 0)
    	{
    		cout << "NNLS solver terminated with an error flag" << endl;
    	}
    	x = solver->getSolution();

    	SaveNNLSOutput(&soln[ch], A, x, newsignaltimes);
    }


    m_data->Stores["ANNIEEvent"]->Set("nnls_solution", soln);



	return true;
}

//makes the nnls matrix A given a root template file. 
//The algorithm requires a careful accounting of the
//number of samples in both the template and the signal waveform.
//Furthermore, to preserve conclusions about the resulting signals,
//one needs to maintain that the signal and the template matrix
//have the same time-step. The variable nrows is the number of rows in the end matrix, 
//and is determined by the samplingFactor from Execute and the number of samples in the template. 
//Both here and in BuildWaveformVector, the template and signal waveform
//are over/undersampled to match timestamps and number of samples.
void WaveformNNLS::BuildTemplateMatrix(nnlsmatrix* A, Waveform<double> tempwave, size_t nrows)
{

	int temp_size = tempwave.GetSamples()->size();
	

	//the nnls matrix consists of:
	//each row is a version of the template waveform
	//inserted starting at a sample indexed by the column number. 
	//the next column then has another template waveform inserted
	//one sample further in time. The elements below the diagonal
	//are 0, and any elements that are not part of the template
	//are 0. 
	for(int row = 0; row < (int) nrows; row++)
	{
		for(int col = 0; col < (int) nrows; col++)
		{
			//elements below diagonal are zero
			if(col < row)
			{
				A->set(row, col, 0.0);
			}
			//on the upper diagonal, the template
			//waveform is written to a row, starting
			//with the column = row
			else if((col - row) < temp_size)
			{
				A->set(row, col, tempwave.GetSample(col - row));
			}
			//otherwise (outside bounds of template)
			//set to 0. 
			else {A->set(row, col, 0.0);}
		}
	}


} 

//formats the waveform into the vector format expected by nnls algo.
//see comment above BuildTemplateMatrix for explanation of nrows
void WaveformNNLS::BuildWaveformVector(nnlsvector* b, Waveform<double> wave, vector<float> times, double template_timestep)
{

	//make a new signal vector that is
	//NOT interpolated, but has more samples
	// so that the timesteps of the template
	// and signal waveform are equal. 
	vector<double> newwave;
	float t0;
	float t1;
	double current_time;
	for(int j = 0; j < (int) b->length(); j++)
	{
		//current time iterating through 
		//times scaled for nnls matrix
		current_time = times.front() + j*template_timestep;

		//find value of waveform at this time by 
		//finding closest sample times. Assumes the
		//times vector is ordered
		for(int i = 0; i < (int) times.size() - 1; i++)
		{
			t0 = times.at(i);
			t1 = times.at(i+1);
			if(current_time >= t0 and current_time < t1)
			{
				b->set(j, wave.GetSample(i));
				break;
			}
		}
	}

} 


//want to save the nnls output such that you have access to the
//full comprehensive fit (A*x = bsolution) as well as the 
//individual template components contained in the vector x. 
//Each element of x represents a template waveform scaled by
//the magnitude of the element and placed at a time "t" 
//(read off of the index of the element)
void WaveformNNLS::SaveNNLSOutput(NnlsSolution* soln, nnlsmatrix* A, nnlsvector* x, vector<double> signaltimes)
{

	//first, the fully composed (full fit) solution
	nnlsvector* bsolv = new nnlsvector(x->length());
	A->dot(false, x, bsolv); //now bsolv is the fitted vector waveform
	//turn vector into waveform
	Waveform<double> ff; //waveform version of nnls full (summed) solution
	for(int i = 0; i < (int)bsolv->length(); i++)
	{
		ff.PushSample(bsolv->get(i));

		//now, the decomposed pulses, one for every
		//element of x that is non-zero.

		//(1e-5 is for possible floating point error)
		//ignore scales that are close to 0
		if(x->get(i) <= 1e-6)
		{
			continue;
		}
		//store the scale and time for each component
		//of the nnls-solved waveform
		else
		{
			soln->AddComponent(signaltimes.at(i), x->get(i));
		}

	}
	soln->SetFullSoln(ff);
}




bool WaveformNNLS::Finalise(){

  return true;
}
