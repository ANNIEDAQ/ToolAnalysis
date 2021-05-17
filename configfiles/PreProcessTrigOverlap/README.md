# PreProcessTrigOverlap

***********************
# Overview
**********************

The `PreProcessTrigOverlap` toolchain needs to be run before DataDecoding jobs on single part files can be run.

If the raw data in a part file `RawDataRxSypz` is supposed to be processed, the raw data file itself doesn't actually contain all the necessary trigger information for a complete decoding. The previous data file `RawDataRxSyp(z-1)` contains record headers that are necessary for decoding the first couple of entries in the current TriggerData file. The subsequent data file `RawDataRxSyp(z+1)` contains necessary trigger information to pair the last couple of VME and MRD timestamps in the current file.

All this necessary information for a rawdata file is saved in a dedicated BoostStore file by running this toolchain. This BoostStore file will be named `TrigOverlap_RxSypz`.

When the decoding chain is run afterwards for this very same rawdata file, the `TriggerDataDecoder` and `LoadRawData` tool will load the information stored in the `TrigOverlapRxSypz` BoostStore and use the information to match all available VME and MRD data in the file to trigger timestamps.

************************
# Usage
************************

Both the `LoadRawData` tool and the `TriggerDataDecoder` tool have two important configuration variables that are used to access the trigger overlap information. Those variables are called `ReadTrigOverlap` and `StoreTrigOverlap`. The `StoreTrigOverlap` variable should be set to 1 when running this `PreProcessTrigOverlap` toolchain to generate the `TrigOverlap` BoostStore files, while the `ReadTrigOverlap` variable should be set to 1 when running the subsequent `DataDecoder` toolchain that will use the previously created `TrigOverlap` BoostStores.

The two variables would hence be set like this in both `TriggerDataDecoderConfig` and `LoadRawDataConfig` for the `PreProcessTrigOverlap` toolchain:

```
ReadTrigOverlap 0
StoreTrigOverlap 1
```

And like this in the `DataDecoder` toolchain:

```
ReadTrigOverlap 1
StoreTrigOverlap 0
```
