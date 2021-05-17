# BeamFetcher

The `BeamFetcher` tool obtains information about the status of the BNB from the IF database in a time interval given by the user. It stores the information in a BoostStore that can be read and accessed by other tools. The `IFBeamDBInterface` class is a helper class that handles the details of the communication with the IF database.

## Data

The following objects will be saved in the BeamStatus BoostStore:
* "BeamDBIndex" (Header) `map<int,pair<uint64_t,uint64_t>`
  * Designates the time interval that is stored in the given BoostStore entry
* "StartMillisecondsSinceEpoch" (Header) `uint64_t`
  * Designates the overall start time of database entries stored in the BoostStore
* "EndMillisecondsSinceEpoch" (Header) `uint64_t`
  * Designates the overall end time of database entries stored in the BoostStore
* "BeamDB" `map<string,map<uint64_t,BeamDataPoint>>`
  * Actual beam status information, keys are the device names (e.g. E:TOR875)

## Configuration

The main configuration variable for the `BeamFetcher` tool is the time frame it is supposed to look at. This time frame can be specified in terms of a start and an end time, either given in milliseconds or in a string format. The preferred format can be chosen with the `TimestampMode` variable (LOCALDATE/MSEC). The format for the string timestamps can be inferred from the timestamps that are stored in the ANNIE PSQL Run database (You can just copy the string for the respective date from the SQL monitoring page).

The `TimeChunkInMilliseconds` variable determines in what time frames the data is saved as an entry to the BeamStatus BoostStore.

```
# BeamFetcher config file
verbose 5
OutputFile ./1604_beamdb
TimestampMode LOCALDATE
DaylightSavings 0
StartDate ./configfiles/BeamFetcher/my_start_date.txt #String form of start date stored in a file
EndDate ./configfiles/BeamFetcher/my_end_date.txt #String form of end date stored in a file
#StartMillisecondsSinceEpoch 1491132659000 # 6:30:49 AM 2 April 2017 (FNAL time) #msec format of start time
#EndMillisecondsSinceEpoch   1491164001000 # 3:13:21 PM 2 April 2017 (FNAL time) #msec format of end time
TimeChunkStepInMilliseconds       7200000 # two hours
```
