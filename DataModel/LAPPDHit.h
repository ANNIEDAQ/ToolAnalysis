/* vim:set noexpandtab tabstop=4 wrap */
#ifndef LAPPDHITCLASS_H
#define LAPPDHITCLASS_H

#include<Hit.h>
#include<SerialisableObject.h>

using std::cout;
using std::endl;

class LAPPDHit : public Hit{
	
	friend class boost::serialization::access;
	
	public:
	LAPPDHit() : Hit(), Position(0), LocalPosition(0) {serialise=true;}
	LAPPDHit(int thetubeid, double thetime, double thecharge, std::vector<double> theposition, std::vector<double> thelocalposition) : Hit(thetubeid,thetime,thecharge), Position(theposition), LocalPosition(thelocalposition) {serialise=true;}
	virtual ~LAPPDHit(){};
	
	inline std::vector<double> GetPosition() const {return Position;}
	inline std::vector<double> GetLocalPosition() const {return LocalPosition;}
	inline void SetPosition(std::vector<double> pos){Position=pos;}
	inline void SetLocalPosition(std::vector<double> locpos){LocalPosition=locpos;}
	
	bool Print() {
		cout<<"TubeId : "<<TubeId<<endl;
		cout<<"Time : "<<Time<<endl;
		cout<<"X Pos : "<<Position.at(0)<<endl;
		cout<<"Y Pos : "<<Position.at(1)<<endl;
		cout<<"Z Pos : "<<Position.at(2)<<endl;
		cout<<"Parallel Pos : "<<LocalPosition.at(0)<<endl;
		cout<<"Transverse Pos : "<<LocalPosition.at(1)<<endl;
		cout<<"Charge : "<<Charge<<endl;
		return true;
	}
	
	protected:
	std::vector<double> Position;
	std::vector<double> LocalPosition;
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & TubeId;
			ar & Time;
			ar & Position;
			ar & LocalPosition;
			ar & Charge;
		}
	}
};

//  Derived classes
class MCLAPPDHit : public LAPPDHit{
	
	friend class boost::serialization::access;
	
	public:
	MCLAPPDHit() : LAPPDHit(), Parents(std::vector<int>{}) {serialise=true;}
	MCLAPPDHit(int thetubeid, double thetime, double thecharge, std::vector<double> theposition, std::vector<double> thelocalposition, std::vector<int> theparents) : LAPPDHit(thetubeid, thetime, thecharge,theposition,thelocalposition), Parents(theparents) {serialise=true;}
	
	const std::vector<int>* GetParents() const { return &Parents; }
	void SetParents(std::vector<int> parentsin){ Parents = parentsin; }
	
	bool Print() {
		cout<<"TubeId : "<<TubeId<<endl;
		cout<<"Time : "<<Time<<endl;
		cout<<"X Pos : "<<Position.at(0)<<endl;
		cout<<"Y Pos : "<<Position.at(1)<<endl;
		cout<<"Z Pos : "<<Position.at(2)<<endl;
		cout<<"Parallel Pos : "<<LocalPosition.at(0)<<endl;
		cout<<"Transverse Pos : "<<LocalPosition.at(1)<<endl;
		cout<<"Charge : "<<Charge<<endl;
		if(Parents.size()){
			cout<<"Parent MCPartice indices: {";
			for(int parenti=0; parenti<(int)Parents.size(); ++parenti){
				cout<<Parents.at(parenti);
				if((parenti+1)<(int)Parents.size()) cout<<", ";
			}
			cout<<"}"<<endl;
		} else {
			cout<<"No recorded parents"<<endl;
		}
		return true;
	}
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & TubeId;
			ar & Time;
			ar & Position;
			ar & LocalPosition;
			ar & Charge;
			// n.b. at time of writing MCHit stores no additional persistent members
			// - it only adds parent MCParticle indices, and these aren't saved... 
		}
	}
	
	protected:
	std::vector<int> Parents;
};

/*
class TDCHit : public Hit {
	public:

	private:
}

class RecoHit : public Hit {
	public:
	RecoHit(double thetime, double thecharge) : Time(thetime), Charge(thecharge){};

	inline double GetCharge(){return Charge;}
	inline void SetCharge(double chg){Charge=chg;}

	private:
	double Charge;
};
*/

#endif
