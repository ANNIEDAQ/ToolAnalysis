#include "LoadRunInfo.h"

LoadRunInfo::LoadRunInfo():Tool(){}


bool LoadRunInfo::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("RunInfoFile",runinfofile);
  m_variables.Get("verbosity",verbosity);

  std::map<int,std::map<std::string,std::string>> RunInfoDB;
  std::vector<int> RunsWithoutEndTime;

  //Loop through txt-file copy of RunInfo database
  //This txt file needs to be updated manually whenever a new run is recorded
  ifstream f(runinfofile);
  std::string line;
  while (getline(f,line)){
    std::istringstream iss(line);
    std::string entry;
    std::vector<std::string> runinfo_vector;
    std::map<std::string,std::string> RunInfoDB_Single;
    if (iss.str().length() == 0) break;	//Don't read out lines without entries
    while (std::getline(iss, entry, '\t')) runinfo_vector.push_back(entry);
    if (runinfo_vector.size() < 8) continue;
    Log("LoadRunInfo tool: Run: "+runinfo_vector.at(1)+", Run Type: "+runinfo_vector.at(6)+", Start Time: "+runinfo_vector.at(3)+", Num Events: "+runinfo_vector.at(7),v_message,verbosity);
    RunInfoDB_Single.emplace("StartTime",runinfo_vector.at(3));
    RunInfoDB_Single.emplace("RunType",runinfo_vector.at(5));
    RunInfoDB_Single.emplace("RunStatus",runinfo_vector.at(6));
    RunInfoDB_Single.emplace("NumEvents",runinfo_vector.at(7));
    if (runinfo_vector.at(4).length() > 0) {
      RunInfoDB_Single.emplace("EndTime",runinfo_vector.at(4));
    } else {
      RunsWithoutEndTime.push_back(std::stoi(runinfo_vector.at(1)));
    }
    RunInfoDB.emplace(std::stoi(runinfo_vector.at(1)),RunInfoDB_Single);
  }


  //For those runs with no end time, use start time of next run as end time
  //Might not be the best solution in all cases (TODO: better solution)
  for (int i_vec=0; i_vec < (int) RunsWithoutEndTime.size(); i_vec++){
    int Run = RunsWithoutEndTime.at(i_vec);
    if (RunInfoDB.count(Run+1)>0){
      RunInfoDB.at(Run).emplace("EndTime",RunInfoDB.at(Run+1)["StartTime"]);
    }
  }


  //Emplace RunInfoDB to CStore for other tools to use
  m_data->CStore.Set("RunInfoDB",RunInfoDB);

  return true;
}


bool LoadRunInfo::Execute(){

  return true;
}


bool LoadRunInfo::Finalise(){

  return true;
}
