/*
This class is a data structure to store the already solved 
nnls algorithm outputs. It is meant to reduce the memory overhead
used by the nnls matrices/vectors and just store info on 
the full solution and components of the solution. Originally, 
components of the solution were stored as full waveforms, but
this would crash due to memory overusage. This class stores
components as a tuple of time index and scale factor applied to template. 
-Evan Angelico, 7/11/2019
*/


#ifndef NnlsSolution_H
#define NnlsSolution_H

#include <SerialisableObject.h>
#include "Waveform.h"

using namespace std;

class NnlsSolution : public SerialisableObject
{
	friend class boost::serialization::access;

public:
	NnlsSolution() : component_times(vector<double>{}), component_scales(vector<double>{}), templateWaveform(Waveform<double>{}),templateTimes(vector<double>{}),fullNnlsWaveform(Waveform<double>{}) {serialise=true;}
	~NnlsSolution();
	void AddComponent(double t, double s);
	void SetFullSoln(Waveform<double> fullsoln); //sum of all components, full waveform fit
	void SetTemplate(Waveform<double> tmpwfm, vector<double> temptimes);
	inline int GetNumberOfComponents() {return component_times.size();} //how many components are populated in components vector
	inline double GetComponentTime(int index) {return component_times.at(index);}
	inline double GetComponentScale(int index) {return component_scales.at(index);}

	bool Print() {
		int verbose=0;
		bool filled = false;
		if(component_times.size() != 0) filled = true;
		cout << "Solution filled?: " << filled << endl;
		if(filled) cout << "Number of components: " << component_times.size() << endl;
		if(verbose){
			cout<<"Times and scales:" << endl;
			for(int i = 0; i < (int) component_times.size(); i++){
				cout << component_times.at(i) << ", " << component_scales.at(i) << endl;
			}
		}

		return true;
	}

	

protected:
	//components of the solution are pairs of doubles
	//first element is the time, second element is the scale
	//at which a template waveform appears relative to the
	//raw waveform in the time units of the raw waveform
	vector<double> component_times;
	vector<double> component_scales;
	Waveform<double> templateWaveform; //template used in the nnls algo, at the time binning used in the algorithm
	vector<double> templateTimes; //time series of length of templateWaveform
	Waveform<double> fullNnlsWaveform; //the full solution, sum of all components

	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & component_times;
			ar & component_scales;
			ar & templateWaveform;
			ar & templateTimes;
			ar & fullNnlsWaveform;
		}
	}
};


#endif
