# PrintANNIEEvent toolchain

***********************
# Description
**********************

The `PrintANNIEEvent` toolchain loops through the contents of an ANNIEEvent BoostStore file and summarizes the information that is contained inside the file. It is supposed to be an overview tool that one can use to get a better feeling of the data that is stored in the file.

It is used in combination with the `LoadANNIEEvent` tool:
* `LoadANNIEEvent`: Loads the ANNIEEvent BoostStore into the shared `m_data->Stores` BoostStore space
* `PrintANNIEEvent`: Accesses the contents of the ANNIEEvent BoostStore and prints the information.

************************
# Usage
************************

The PrintANNIEEvent toolchain can be configured in the following ways:
* The verbosity can be set with the "verbose" keyword
* The data/MC nature of the ANNIEEvent store can be selected with the "IsMC" keyword
* The availability of raw PMT data (waveforms) in the ANNIEEvent BoostStore can be specified with the "HasRaw" keyword.

```
verbose 3
IsMC 0          # Is it a MC (1) or data (0) file?
HasRaw 0        # Does file contain raw data objects?
```
