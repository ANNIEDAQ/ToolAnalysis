# LoadRunInfo

The `LoadRunInfo` tool loads run-specific information from a txt-file and stores the information in a map. The txt-file can be a simply copy of the SQL RunInfo page and is needed because the database is not accessible everywhere. The map with the information is stored in the CStore, so other tools can access the run information.

## Data

The tool fills the following map with run-related information:

**RunInfoDB** `map<int, map<string,string>>`
* Keys are the run numbers, values are maps of strings that includes the detailed run information.

The following information is saved:
* "StartTime": The start time of the run
* "EndTime": The end time of the run
* "RunType": The run type
* "NumEvents": Number of (VME) events in run. Does not correspond to number of physics-quality events.


## Configuration

LoadRunInfo primarily needs the path of the RunInfo txt file in order to access the information.

```
verbosity 1
/pnfs/annie/persistent/users/mnieslon/data/runinfo/ANNIE_RunInformation_PSQL.txt
```
