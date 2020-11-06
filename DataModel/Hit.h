/* vim:set noexpandtab tabstop=4 wrap */
#ifndef HITCLASS_H
#define HITCLASS_H

#include<SerialisableObject.h>

#include <iostream>

using namespace std;

class Hit : public SerialisableObject{
	
	friend class boost::serialization::access;
	
	public:
	Hit() : TubeId(0), Time(0), Charge(0){serialise=true;}
	Hit(int thetubeid, double thetime, double thecharge) : TubeId(thetubeid), Time(thetime), Charge(thecharge){serialise=true;}
	virtual ~Hit(){};
	
	inline int GetTubeId() const {return TubeId;}
	inline double GetTime() const {return Time;}
	inline double GetCharge() const {return Charge;}
	
	inline void SetTubeId(int tubeid){TubeId=tubeid;}
	inline void SetTime(double tc){Time=tc;}
	inline void SetCharge(double chg){Charge=chg;}
	
	bool Print() {
	  std::cout<<"TubeId : "<<TubeId<<endl;
	  std::cout<<"Time : "<<Time<<endl;
	  std::cout<<"Charge : "<<Charge<<endl;
		return true;
	}
	
	protected:
	int TubeId;
	double Time;
	double Charge;
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & TubeId;
			ar & Time;
			ar & Charge;
		}
	}
};

//  Derived classes

class MCHit : public Hit {
	// XXX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XXX
	// XXX ~~~~~~~~~~~~~~~~~~~~~~~~ UPDATING THIS CLASS? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XXX
	// XXX ~~~~~ Everything added in this class must be duplicated in MCLAPPDHit!~~~~ XXX
	// XXX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XXX
	
	friend class boost::serialization::access;
	
	public:
	MCHit() : Hit(), Parents(std::vector<int>{}) {serialise=true;}
	MCHit(int tubeid, double thetime, double thecharge, std::vector<int> theparents) : Hit(tubeid, thetime, thecharge), Parents(theparents) {serialise=true;}
	virtual ~MCHit(){};
	
	const std::vector<int>* GetParents() const { return &Parents; }
	void SetParents(std::vector<int> parentsin){ Parents = parentsin; }
	
	bool Print(){
	  std::cout<<"TubeId : "<<TubeId<<endl;
	  std::cout<<"Time : "<<Time<<endl;
		std::cout<<"Charge : "<<Charge<<endl;
		if(Parents.size()){
			std::cout<<"Parent MCPartice indices: {";
			for(int parenti=0; parenti<(int)Parents.size(); ++parenti){
				std::cout<<Parents.at(parenti);
				if((parenti+1)<(int)Parents.size()) std::cout<<", ";
			}
			std::cout<<"}"<<endl;
		} else {
			std::cout<<"No recorded parents"<<endl;
		}
		return true;
	}
	
	protected:
	std::vector<int> Parents;
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & TubeId;
			ar & Time;
			ar & Charge;
			// do not serialize parents; the indices by themselves are not meaningful
		}
	}
};

/*
class RecoHit : public Hit {
	public:
	RecoHit(double thetime, double thecharge) : Time(thetime), Charge(thecharge){};

	inline double GetCharge(){return Charge;}
	inline void SetCharge(double chg){Charge=chg;}

	protected:
	double Charge;
};
*/

#endif
