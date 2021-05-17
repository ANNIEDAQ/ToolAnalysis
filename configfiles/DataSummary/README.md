# DataSummary toolchain

***********************
# Overview
**********************

The `DataSummary` toolchain provides an overview of the event building process in ANNIE. Processed data files are investigated with respect to the number of orphaned events, the reasons for orphaned events and the properties of the merged data streams in general. It saves the information about merged events and orphaned events in two separate ROOT trees that can be accessed for further analysis. In addition, some diagnostic histograms are also produced and saved in the same output file.


************************
# Configuration
************************

The toolchain consists of only one tool, i.e. the `DataSummary` tool. One has to provide the range of runs, subruns and part files that one wants to look at in the configuration file. Furthermore, it needs to be configured whether both the `ANNIEEvent` and the `OrphanStore` data was saved in a combined BoostStore output file or whether there are  two separate BoostStore files for the two data streams. This can be set with the `FileFormat` configuration variable.

Exemplary configuration file:

```
# DataSummary config file

verbosity 10
FileList None       #FileList provides a separate way of feeding a simple list of files as inputs instead of the run/subrun/part file ranges
DataPath .

# input files will be found by a regex search. All subdirectories within DataPath will be searched.
# Filenames must at minimum end with 'R**S**p**', which is used to extract the run, subrun and part numbers
# To restrict the search further a regular expression 'InputFilePattern' may be given.
# filenames must then be of the form [anything][InputFilePattern][RXXXSYYYpZZZ]
InputFilePattern ProcessedRawData_TankAndMRDAndCTC_
InputFilePatternOrphan OrphanStore_TankAndMRDAndCTC_
FileFormat SeparateStores           #Other option: CombinedStore

# specify range of data to analyse. Use -1 to match anything.
StartRun 2282
StartSubRun 0
StartPart 4
EndRun -1
EndSubRun -1
EndPart 4

# ROOT file output
OutputFileDir .
OutputFileName DataSummary_R2282S0p4_trigoverlap.root
```
