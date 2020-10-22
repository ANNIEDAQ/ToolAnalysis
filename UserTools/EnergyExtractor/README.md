# EnergyExtractor

The `EnergyExtractor` tool can be used to extract information about the true energy of particles into csv files. This is particularly useful if one wants to analyze the behavior of classification algorithms as a function of the true energy of particles. In this case it is useful to run the `EnergyExtractor` tool e.g. in combination with the `CNNImage` tool in the same toolchain.

## Data

The `EnergyExtractor` tool primarily uses the information from the `MCParticles` object about the true energies of the involved particles. In addition, it also uses the `MCHits` and `MCLAPPDHits` objects to evaluate the visible charge of the event. It also uses the number of neutrons that were produced in the event.

**MCParticles** `vector<MCParticle>`
* Takes this vector of MCParticles and extracts the true energies of the particles

**MCHits/MCLAPPDHits** `map<unsigned long,vector<MCHit>>/map<unsigned long,vector<MCLAPPDHit>>`
* Takes these maps and extracts the visible charge

**EventCutStatus** `bool`
* Checks this variable to make sure specified event selection cuts were passed in the `EventSelector` tool.


## Configuration

It can be specified which particles the energies should be extracted from as well as the general name of the output file.

```
OutputFile beam_0
Verbosity 1
SaveNeutrino 0/1
SaveElectron 0/1
SaveGamma 0/1
SaveMuon 0/1
SavePion 0/1
SaveKaon 0/1
SaveNeutron 0/1
SaveVisible 0/1
```
