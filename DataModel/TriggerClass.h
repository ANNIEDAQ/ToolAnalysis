/* vim:set noexpandtab tabstop=4 wrap */
#ifndef TRIGGERCLASS_H
#define TRIGGERCLASS_H

#include<SerialisableObject.h>
#include "TimeClass.h"

using namespace std;

class TriggerClass : public SerialisableObject{
	
	friend class boost::serialization::access;
	
	public:
  TriggerClass() : TriggerName(""), TriggerWord(0), ExtendedWindow(0), TriggerOccurred(false), TriggerTime(){serialise=true;}
  TriggerClass(std::string trigname, uint32_t trigword, int extended, bool trigocc, TimeClass trigtime) : TriggerName(trigname), TriggerWord(trigword), ExtendedWindow(extended), TriggerOccurred(trigocc), TriggerTime(trigtime){serialise=true;}
	
	inline std::string GetName(){return TriggerName;}
	inline TimeClass GetTime(){return TriggerTime;}
	inline bool GetOccurred(){return TriggerOccurred;}
	inline uint32_t GetWord(){return TriggerWord;}
	inline int GetExtended(){return ExtendedWindow;}
	
	inline void SetName(std::string namein){TriggerName=namein;}
	inline void SetOccurred(bool trigoc){TriggerOccurred=trigoc;}
	inline void SetTime(TimeClass trigt){TriggerTime=trigt;}
	inline void SetWord(uint32_t wordin){TriggerWord=wordin;}
	inline void SetExtended(int extin){ExtendedWindow=extin;}
	

	bool Print() {
		cout<<"TriggerName : "<<TriggerName<<endl;
		cout<<"TriggerWord : "<<TriggerWord<<endl;
		cout<<"ExtendedWindow : "<<ExtendedWindow<<endl;
		cout<<"TriggerTime : "; TriggerTime.Print();
		cout<<"TriggerOccurred : "<<TriggerOccurred<<endl;
		
		return true;
	}
	
	private:
	std::string TriggerName;
	uint32_t TriggerWord;
	int ExtendedWindow;
	bool TriggerOccurred;
	TimeClass TriggerTime;
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & TriggerName;
			ar & TriggerWord;
			ar & ExtendedWindow;
			ar & TriggerTime;
			ar & TriggerOccurred;
		}
	}
};

#endif
