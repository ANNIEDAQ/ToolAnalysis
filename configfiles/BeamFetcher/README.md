# BeamFetcher toolchain

***********************
# Description
**********************

The `BeamFetcher` toolchain gets beam status information from the IF beam database and stores the information in a BoostStore file.

Tools in ToolChain:
* `LoadRunInfo`: Tool to load start and end times of single runs in case one needs to load the beam status information for the duration of a single run
* `BeamFetcher`: Tool that loads the information about POT and horn currents for the specified time frame

************************
#Usage
************************

The `BeamFetcher` tool can be configured in multiple ways: either via a start and end time or via the duration of a specific run. This is configured with the `TimeStampMode` variable:
* Option "MSEC": Start and end times have to be specified in milliseconds since 1970/1/1 via the "StartMillisecondsSinceEpoch" and "EndMillisecondsSinceEpoch" variables. Please note that those timestamps are in UTC time. Output filename will be set according to the `OutputFile` configuration variable.
* Option "LOCALDATE: Start and end times have to be specified as a string that is contained in the "StartDate" and "EndDate" files. One can copy the start and end times for specific runs from the database. Note that those times are in Fermilab time, they are converted to UTC within the tool. Output filename will be set according to the `OutputFile` configuration variable.
* Option "DB": Start and end times are taken from the Run Information database text-file, one only needs to specify the run number with the "RunNumber" variable. The output file name will be adjusted automatically according to the syntax ${RunNumber}_beamdb (the `OutputFileName` configuration variable is ignored in this case)

```
verbose 5
OutputFile ./2369_beamdb
TimestampMode DB        #Options: MSEC / LOCALDATE / DB
DaylightSavings 0
RunNumber 2295
StartDate ./configfiles/BeamFetcher/my_start_date.txt
EndDate ./configfiles/BeamFetcher/my_end_date.txt
#StartMillisecondsSinceEpoch 1491132659000 # 6:30:49 AM 2 April 2017 (FNAL time)
#EndMillisecondsSinceEpoch   1491164001000 # 3:13:21 PM 2 April 2017 (FNAL time)
TimeChunkStepInMilliseconds       7200000 # two hours
```

