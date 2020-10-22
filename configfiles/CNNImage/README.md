# CNNImage toolchain

The `CNNImage` toolchain creates csv-files representing pictures of the ANNIE detector in charge & time. The picture can either be a rolled-out version of the tank where the top and bottom caps are represented by circles. The other option is a rolled-out version of the tank where the top and bottom caps are also stretched out over the whole length of the tank.

************************
# CNNImage_Data
************************

The `CNNImage_Data` toolchain operates on data and hence uses the `LoadANNIEEvent` tool to load in the data. The recommended sequence of tools is

* `LoadGeometry`: Load the geometry of ANNIE.
* `LoadANNIEEvent`: Load the event information from the `ANNIEEvent` BoostStore.
* `PhaseIIADCCalibrator`: Calibrates the baseline of the PMT ADC signal.
* `PhaseIIADCHitFinder`: Finds hits in the PMT signal.
* `ClusterFinder`: Finds (time) clusters in the obtained PMT hits.
* `DigitBuilder`: Builds reco digits from hits in found clusters.
* `TimeClustering`: Finds (time) clusters in MRD hits.
* `EventSelector`: Choose event selection cuts.
* `CNNImage`: Construct charge/time detector images from the clustered hits objects.

************************
# CNNImage_MC
************************

The `CNNImage_MC` toolchain operates on MC and hence uses the `LoadWCSim`/`LoadWCSimLAPPD` tools to load the data into the `ANNIEEvent` BoostStore. The recommended sequence of tools is

* `LoadWCSim`: Load the PMT information from the WCSim simulation file, as well as geometry information.
* `LoadWCSimLAPPD`: Load the LAPPD information from the WCSim simulation file.
* `MCRecoEventLoader`: Monte Carlo properties are loaded into the RecoEvent BoostStore.
* `DigitBuilder`: Build reco digits from prompt MC hits.
* `TimeClustering`: Finds (time) clusters in MRD hits.
* `EnergyExtractor`: Extracts energies of primary particles into separate csv-files. (For performance evaluations of the network as a function of those energies)
* `EventSelector`: Choose event selection cuts.
* `CNNImage`: Construct charge/time detector images from the MC hits objects.
