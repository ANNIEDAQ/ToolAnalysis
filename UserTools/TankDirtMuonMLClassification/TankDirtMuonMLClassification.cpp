#include "TankDirtMuonMLClassification.h"

TankDirtMuonMLClassification::TankDirtMuonMLClassification():Tool(){}


bool TankDirtMuonMLClassification::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool TankDirtMuonMLClassification::Execute(){

  return true;
}


bool TankDirtMuonMLClassification::Finalise(){

  return true;
}
