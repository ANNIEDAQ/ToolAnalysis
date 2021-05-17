// standard library includes
#include <fstream>

// ToolAnalysis includes
#include "LoadANNIEEvent.h"

LoadANNIEEvent::LoadANNIEEvent():Tool() {}

bool LoadANNIEEvent::Initialise(std::string config_filename, DataModel &data) {

  // Load settings from the configuration file
  if ( !config_filename.empty() ) m_variables.Initialise(config_filename);

  // Assign transient data pointer
  m_data= &data;
  offset_evnum = 0;
  FileFormat = "SeparateStores";	//Other option: "CombinedStore"
  load_orphan_store = false;

  m_variables.Get("verbose", verbosity_);
  m_variables.Get("EventOffset", offset_evnum);
  m_variables.Get("FileFormat",FileFormat);
  m_variables.Get("LoadOrphanStore",load_orphan_store);

  std::string input_list_filename;
  bool got_input_file_list = m_variables.Get("FileForListOfInputs",
    input_list_filename);

  if ( !got_input_file_list ) {
    Log("Error: Missing input list file in the configuration for the"
      " LoadANNIEEvent tool", 0, verbosity_);
    return false;
  }

  std::ifstream list_file(input_list_filename);
  if ( !list_file.good() ) {
    Log("Error: Could not open the input list file for the LoadANNIEEvent tool",
      0, verbosity_);
    return false;
  }

  std::string temp_str;
  while ( list_file >> temp_str ) input_filenames_.push_back( temp_str );

  if (load_orphan_store && FileFormat == "SeparateStores"){
    std::string input_list_filename_orphan;
    bool got_input_file_list_orphan = m_variables.Get("FileForListOfInputsOrphan", input_list_filename_orphan);
    if (!got_input_file_list_orphan){
      Log("Error: Missing orphan input list file in the configuration for the"
      " LoadANNIEEvent tool", 0, verbosity_);
      return false;
    }
    std::ifstream list_file_orphan(input_list_filename_orphan);
    if ( !list_file_orphan.good() ) {
      Log("Error: Could not open the orphaninput list file for the LoadANNIEEvent tool",
        0, verbosity_);
      return false;
    }

    std::string temp_str_orphan;
    while ( list_file_orphan >> temp_str_orphan ) input_filenames_orphan_.push_back( temp_str_orphan );


  }

  current_entry_ = 0u;
  current_file_ = 0u;
  need_new_file_ = true;

  m_data->CStore.Set("UserEvent",false);

  current_entry_ += offset_evnum;
 
  return true;
}


bool LoadANNIEEvent::Execute() {

  int stop_the_loop = -1;
  m_data->vars.Get("StopLoop", stop_the_loop);
  if ( stop_the_loop == 1 ) return false;

  if (need_new_file_) {
    need_new_file_=false;

    if (FileFormat == "CombinedStore"){
      // Delete the old file BoostStore if there is one
      if(m_data->Stores.count("ProcessedFileStore")){
        BoostStore* ProcessedFileStore = m_data->Stores.at("ProcessedFileStore");
        //ProcessedFileStore->Close();   // XXX should we be calling these?
        //ProcessedFileStore->Delete();  // XXX
        delete ProcessedFileStore;
      }
    }

    // also close and cleanup the associated contained stores
    if ( m_data->Stores.count("ANNIEEvent") ) {
      auto* annie_event = m_data->Stores.at("ANNIEEvent");
      if (annie_event){
        //annie_event->Close();
        //annie_event->Delete();
        delete annie_event;
      }
    }

    if (load_orphan_store){
      if(m_data->Stores.count("OrphanStore")){
        auto* annie_orphans = m_data->Stores.at("OrphanStore");
        if (annie_orphans){
          //annie_orphans->Close();
          //annie_orphans->Delete();
        delete annie_orphans;
        }
      }
    }

    if (FileFormat == "CombinedStore"){
      // create a store for the file contents
      BoostStore* ProcessedFileStore = new BoostStore(false,BOOST_STORE_BINARY_FORMAT);
      // Load the contents from the new input file into it
      std::string input_filename = input_filenames_.at(current_file_);
      std::cout <<"Reading in current file "<<current_file_<<std::endl;
      ProcessedFileStore->Initialise(input_filename);
      m_data->Stores["ProcessedFileStore"]=ProcessedFileStore;
    
      // create an ANNIEEvent BoostStore and an OrphanStore BoostStore to load from it
      BoostStore* theANNIEEvent = new BoostStore(false,
        BOOST_STORE_MULTIEVENT_FORMAT);
    
      // retrieve the multi-event stores
      ProcessedFileStore->Get("ANNIEEvent",*theANNIEEvent);
      // set a pointer into the Stores map
      m_data->Stores["ANNIEEvent"]=theANNIEEvent;
      // get the number of entries
      m_data->Stores.at("ANNIEEvent")->Header->Get("TotalEntries",
        total_entries_in_file_);
    
      if (load_orphan_store){
        // same for Orphan Store
        BoostStore* theOrphanStore = new BoostStore(false,
          BOOST_STORE_MULTIEVENT_FORMAT);
        ProcessedFileStore->Get("OrphanStore",*theOrphanStore);
        m_data->Stores["OrphanStore"] = theOrphanStore;
        m_data->Stores.at("OrphanStore")->Header->Get("TotalEntries",
        total_orphans_in_file_);
      }
    } else if (FileFormat == "SeparateStores"){
      BoostStore *theANNIEEvent = new BoostStore(false,
        BOOST_STORE_MULTIEVENT_FORMAT);
      std::string input_filename = input_filenames_.at(current_file_);
      theANNIEEvent->Initialise(input_filename);
      m_data->Stores["ANNIEEvent"] = theANNIEEvent;
      m_data->Stores.at("ANNIEEvent")->Header->Get("TotalEntries",total_entries_in_file_);

      if (load_orphan_store){
        BoostStore *theOrphanStore = new BoostStore(false,
          BOOST_STORE_MULTIEVENT_FORMAT);
        std::string input_filename_orphan = input_filenames_orphan_.at(current_file_);
        theOrphanStore->Initialise(input_filename_orphan);
        m_data->Stores["OrphanStore"] = theOrphanStore;
        m_data->Stores.at("OrphanStore")->Header->Get("TotalEntries",
          total_orphans_in_file_);
      }
    }
  }

   bool user_event=false;
   m_data->CStore.Get("UserEvent",user_event);
   if (user_event){
     m_data->CStore.Set("UserEvent",false);
     int user_evnum;
     m_data->CStore.Get("LoadEvNr",user_evnum);
     if (user_evnum < total_entries_in_file_ && user_evnum >=0) current_entry_ = user_evnum;
   }

  Log("ANNIEEvent store has "+std::to_string(total_entries_in_file_)+" entries",v_debug,verbosity_);
  Log("Loading entry " + std::to_string(current_entry_) + " from the"
    " ANNIEEvent input file \"" + input_filenames_.at(current_file_)
    + '\"', 1, verbosity_);
 
  if (current_entry_ != offset_evnum) m_data->Stores["ANNIEEvent"]->Delete();	//ensures that we can access pointers without problems

  m_data->Stores["ANNIEEvent"]->GetEntry(current_entry_);  
  ++current_entry_;
  
  if ( current_entry_ >= total_entries_in_file_ ) {
    ++current_file_;
    if ( current_file_ >= input_filenames_.size() ) {
      m_data->vars.Set("StopLoop", 1);
    }
    else {
      current_entry_ = 0u;
      need_new_file_ = true;
    }
  }


  return true;
}


bool LoadANNIEEvent::Finalise() {
  return true;
}
