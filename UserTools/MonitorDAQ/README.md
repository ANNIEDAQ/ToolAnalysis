# MonitorDAQ

The `MonitorDAQ` tool checks whether the currently produced data files have reasonable sizes and whether all VME_service processes are running in the DAQ.

## Data

`MonitorDAQ` uses information provided by the `MonitorReceive` tool about the most recent file timestamp and file size. Furthermore it checks the currently running processes in the DAQ and evaluates the number of running `VME_service` processes.

## Configuration

The configuration parameters for MonitorDAQ are similar to the other `Monitor*`-Tools.

```
verbose 5               #Verbosity setting
OutputPath fromStore    #Either fromStore or explicit path
StartTime 1970/1/1 #start time for the conversion of timestamps into dates
PlotConfiguration configfiles/Monitoring/DAQTimePlotConfig.txt #file containing instructions what to plot
PathMonitoring /monitoringfiles/	#path at which the monitoring plots are going to be saved
ImageFormat png	#format in which monitoring plots are saved. Options: png, jpg, jpeg
UpdateFrequency 1. #specify frequency for the file history plot, in mins
ForceUpdate 0	#force monitor plots to be produced even if there was no new data file available
DrawMarker 0	#specify whether to use markers for the time evolution graphs or not
UseOnline 1	#Should the VME services be obtained online? Only possible in network
```
