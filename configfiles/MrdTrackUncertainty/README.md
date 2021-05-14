# MrdTrackUncertainy

## Description

The `MrdTrackUncertainty` toolchain can be used to evaluate the MRD track reconstruction performance with or without applied MRD paddle efficiency properties. In case the paddle efficiency values should be applied to the data, the tool `ApplyMRDEff` is needed and the `TimeClustering` tool needs to be configured to read in the modified `TDCData` object.

## Tools in the ToolChain

A typical toolchain to evaluate the MRD reconstruction performance will have the following sequence of tools. The `ApplyMRDEff` tool should be used if one wants to evaluate the effects of inefficient MRD paddles on the reconstruction.

* `LoadWCSim`: Loads the data from the MC file
* `MCParticleProperties`: Extends the MC particle information
* `ApplyMRDEff`: Applies MRD efficiencies to the simulation data
* `TimeClustering`: Clusters MRD hits in time
* `FindMrdTracks`: Finds MRD tracks in clustered MRD data.
* `MrdPaddlePlot`: Plots MRD paddles and reconstructed tracks
* `MrdEfficiency`: Evaluates the efficiency of how many true tracks were successfully reconstructed
* `MrdDistributions`: Plots distributions showing the properties of the reconstructed MRD tracks. 

