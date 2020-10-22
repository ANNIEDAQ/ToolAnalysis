# CNNImage

CNNImage creates ANNIE event display information in a csv-file format that can directly be loaded into Machine Learning classifier frameworks like CNNs for image classification purposes. Currently only PMT information from the side PMTs is loaded, no information from the top/bottom PMTs and LAPPDs is used.

## Data

CNNImage creates multiple `.csv`-files and one `.root`-file. Each csv-file contains the event display information in single rows  per event for a specific event type (PMTs/LAPPDs, charge/time), whereas the root-file provides the same information in a 2D histogram format.

The PMT data can be written PMT-wise (every pixel is reserved for just 1 PMT, no spaces in between) or geometrically (every pixel is defined by its geometric location, if a PMT happens to be within the pixel its value is added to the pixel value). The configuration variables `DimensionX` and `DimensionY` hence only apply when using the `Geometric` save mode and not for the `PMT-wise` save mode. The dimensions for the `PMT-wise` save mode are simply given by the number of PMTs in x- and y-direction. When choosing `PMT-wise`, the LAPPD information is also saved per LAPPD in a 20x20cm grid.

The `DataMode` configuration specifies whether the average hit time per PMT is calculated as a simple mean (`Normal`) or weighted by the charge (`Charge-Weighted`).


## Configuration

CNNImage uses the following configuration variables:

```
verbosity 1     
DataMode Normal			#options: Normal / Charge-Weighted
SaveMode Geometric		#options: Geometric / PMT-wise
DimensionX 16			#choose something suitable (32/64/...)
DimensionY 20			#choose something suitable (32/64/...)
OutputFile beam_test0		#csv file name
DetectorConf ANNIEp2v6		#specify the detector version used in simulation
```
