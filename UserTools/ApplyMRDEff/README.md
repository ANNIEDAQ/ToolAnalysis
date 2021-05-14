# ApplyMRDEff

The `ApplyMRDEff` tool artificially removes hits from the `TDCData` object according to the measured paddle efficiencies. 

## Data

`ApplyMRD` efficiency takes the `TDCData` object and decides for each hit whether it should be kept or rejected according to the paddle efficiency value. The modified MRD data object containing the removed paddle hits is then stored in the `ANNIEEvent` store. 

**TDCData** `map<unsigned long, vector<MCHit>>`
* Takes this data from the `ANNIEEvent` store and evaluates each paddle hit

**TDCData_mod** `<map<unsigned long, vector<MCHit>>`
* Writes this modified `TDCData` object to the `ANNIEEvent` store. 

The specifics about which hits were rejected can be stored in a separate ROOT file if requested in the configuration file (`DebugPlots 1`). This ROOT file will have a tree with the following variables:
* `evnum (int)`: Event number
* `vec_random (vector<double>)`: Random variables used for the hits in this event.
* `dropped_chkey (vector<unsigned long>)`: Dropped channelkeys (data)
* `dropped_chkey_mc (vector<unsigned long>)`: Dropped channelkeys (MC)
* `dropped_ch_time (vector<double>)`: Time values of dropped hits.


## Configuration

`ApplyMRDEff` has the following configuration variables:

```
verbosity 1   #How verbose the tool is
MRDEffFile mrd_efficiencies_10-13-2020.txt  #The file specifying the MRD efficiency values
OutputFile ApplyMRDEff_tree   #Output file name for a ROOT tree that contains specifics on the removed paddle hits
DebugPlots 1  #Boolean to decide whether the debug ROOT tree should be filled
MapChankey_WCSimID  #File that contains a map of channelkeys to WCSim IDs
UseFileEff 0  #Should pre-computed values be used for the removal of paddle hits?
FileEff ApplyMRDEff_tree.root   #If pre-computed values should be used, look for them in this file
```
