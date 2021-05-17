# BeamDecoder

The `BeamDecoder` tool is part of the Event Building chain in ANNIE and forwards information about the Beam Status to the `ANNIEEventBuilder` tool. It uses information that was previously retrieved from the beam database with the `BeamFetcher` tool. Note that the basic functionality of this tool is a blatant copy of Steven's `BeamChecker` tool and was changed in a way to integrate this beam status information in the Event Building process. 

The information is saved in the form of the `BeamStatus` class. This class contains some basic information like the POT for the timestamp in question and some more detailed information about the horn currents. It is possible to already choose some beam quality cuts for the timestamp tolerance, horn currents, and POT values. However, the full information will be stored in the object, so it will always be possible to use slightly different cuts when analyzing the data later.

## Data

The Beam Status information is stored in the `BeamStatusMap` object and put in the CStore. The `ANNIEEventBuilder` tool can access the object in the CStore and write the information to the ANNIEEvent BoostStore.

The `BeamDecoder` tool goes through the decoded trigger timestamps and searches for the beam status at each of these trigger timestamps (in case there was a beam trigger). The properties of the beam are then saved in the BeamStatus object and put into the `BeamStatusMap`.

**BeamStatusMap** `map<uint64_t, BeamStatus>`
* Beam status for the trigger timestamps

The `BeamStatusMap` is stored in the form of a pointer, and the `ANNIEEventBuilder` will delete already built entries from the map to free up memory.

## Configuration

BeamDecoder has the following configuration variables:

```
# BeamDecoder config file
verbosity 2
# Names of devices needed for beam quality cuts
HornCurrentDevice E:THCURR
# The "first" toroid is the one farther upstream from the target
FirstToroid E:TOR860
SecondToroid E:TOR875
# POT window
CutPOTMin 5e11
CutPOTMax 8e12
# Peak horn current window (in kA)
CutPeakHornCurrentMin 172
CutPeakHornCurrentMax 176
# Toroid agreement tolerance (fractional error)
CutToroidAgreement 0.05
# DB vs DAQ timestamp agreement tolerance (ms)
CutTimestampAgreement 100
```
