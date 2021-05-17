# DataDecoderTankCTC

***********************
# Overview
**********************

The `DataDecoderTankCTC` toolchain enables the Event Building process that merges the timestamps from the tank VME system with timestamps of selected triggerwords from the CTC.


************************
# Usage
************************

The `DataDecoderTankCTC` toolchain needs the following tools in the `ToolsConfig` file:

* `LoadGeometry`
* `LoadRawData`
* `PMTDataDecoder`
* `TriggerDataDecoder`
* `ANNIEEventBuilder`

Both the `LoadRawData` tool and the `ANNIEEventBuilder` tool need to set the `BuildType` configuration variable to `TankAndCTC`.

