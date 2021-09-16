**GFA v.3.0**
----------

### Geomorphic Flood Area
----------
<p align="left"><img src="https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area/blob/master/GeomorphicFloodIndex/icona.png"   width="150" height="150" //></p>

**GFA - tool** is an open-source QGIS plug-in to realize a fast and cost-effective delineation of the floodplains in the contexts where the available data is scarce to carry out hydrological/hydraulic analyses. This algorithm uses the linear binary classification technique based on the recently proposed Geomorphic Flood Index (GFI), (Manfreda et al., 2015; Samela et al., 2016; Samela et al., 2017)


## Table of Contents

* [**Team - HydroLAB**](#team)

* [**Project Details**](#project-details)  
    * [Scope](#scope)
    * [Background](#background)
    * [License](#license)
    * [Metadata](#metadata)

* [**GFA Installation Instruction**](#Installation-Instruction)  

* [**GFA Quick Start**](#quick-start)
    * [Input data](#input)
    * [Benchmark Input DataSet](#inputdata)
    * [output data](#input)

* [**Reference**](#reference)

* [**Bug Report**](#bug-report)

## Team - HydroLAB


- [Raffaele Albano](http://www2.unibas.it/raffaelealbano) (Co-Founder of [Wat-TUBE](http://wat-tube.it) spin-off UNIBAS, Research Associate at University of Basilicata (UNIBAS))
- Samela Caterina (Research Associate at University of Basilicata (UNIBAS))
- Aurelia Sole (Professor at University of Basilicata  (UNIBAS))
- [Salvatore Manfreda](http://www2.unibas.it/manfreda/HydroLab) (Professor at University of Basilicata, Head of HydroLAB, Co-Founder of [Wat-TUBE](http://wat-tube.it) spin-off UNIBAS)

<p align="left"><img src="https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area/blob/master/GeomorphicFloodIndex/HydroLAB.PNG" width="300" height="100" /></p>

## Project Details
Learn more about the **GFA** project: scope, background and licensing.

### Scope
Nowadays, the most used approach to modelling fluvial hydraulics and to obtain hydraulic hazard maps, is to make use of hydrological and hydraulic models. 
Generally, more is sophisticated the flood inundation model and higher is the accuracy.
However, we pay this accuracy with a price: they are expensive and time consuming; moreover, they require an extensive input dataset not readily available for all areas. The scarcity of adequate data for flood hazard studies is generally most pronounced in developing counties. This gap has stimulated a lot of research in this area to overcome these limitations and several researchers have recently shown that the delineation of flood-prone areas can be carried out using simplified methods that rely on the analysis of the basin morphology (Manfreda et al., 2011; Mafreda et al., 2014; Mafreda et al., 2015; Samela et al., 2016; Samela et al., 2017). These methods are based on the following assumption: in a drainage basin, a mutual causal relationship exists between flooding and the shape and extension of floodplains.

### Background
Samela et al. (2018) proposed a practical and cost-effective procedure for the preliminarily delineation of the flood-prone areas in data poor environments and for large-scale analyses based on information easily available worldwide.  The potentials of a geomorphic classifier, based on the recently proposed Geomorphic Flood Index, have been implemented in the open-source geographic information system Quantum GIS in the form of a new plugin named Geomorphic Flood Area – tool (GFA – tool) designed to create a user-friendly interface for the detection of the food-prone areas. Moreover, the tool allows also to generate a number of complementary information, such as the GFI, that may be used as ancillary data for remote sensing detection of inundated areas. 


### License
This project is completely licensed [GPL v2+](https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area/blob/master/GeomorphicFloodIndex/LICENSE.txt).

### Metadata
| Metadata Title | Description |
|----------------|-------------|
|  Current Code version              |    v. 1.0         |
|   Permanent link to code / repository used of this code version             |       https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area       |
|       Software Code Language used          |        Python (v. 2.7)       |
|    Compilation requirements, Operating environments & dependencies             |    QGIS 2 or higher         |
|   Computing platform / Operating System             |   Linux, OS X, Microsoft Windows          |
|     License           |      GNU GPL v.2       |

## GFA Installation Instruction
**GFA - tool** can be installed in [QGIS](https://qgis.org) using the Plugin Manager (see [here](http://docs.qgis.org/2.8/en/docs/user_manual/plugins/plugins.html#managing-plugins)). The plugin was tested on the 2.14 ltr release of QGIS. It uses the python library of QGIS core. 

## GFA Quick Start


### Input data
- DEM: Digital Elevation Model grid
- Filled DEM: Depressionless DEM grid
- Flow Direction D8: eight direction flow model grid
- Flow Accumulation: Flow accumulation grid

### Benchmark Input DataSet
A set of sample data to test the plugin is provided 
[here](https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area-doc).

### Output data
- GFI raster: the Geomorphic Flood Index map (Samela et al., 2017)
- GFI normalized raster: the Geomorphic Flood Index map, with values normalized in the range -1:1
- GFI derived flood-prone areas map: the binary raster of the flood prone areas identified by the GFI classifier
- GFI performance metrics: a text file which stores the calibration threshold, the false positive rate, the false negative rate, the area under the ROC curve
- Create intermediate files: locally store all the main intermediate file utilized to produce the above results

## Reference
- Greenlee, D., (1987). Raster and vector processing for scanned linework. Photogrammetric Engineering and Remote Sensing, 53(10), 1383–1387.
- Cazorzi, F. (2002), Software GIS HydroGrid2002 (HyGrid2k2), 1-38, 16.
- Garbrecht, J. and Martz, L. W., (1997). The Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models, Journal of Hydrology, 193: 204-213.
- Manfreda, S., Di Leo, M., Sole, A. (2011). Detection of Flood Prone Areas using Digital Elevation Models, Journal of Hydrologic Engineering, 16(10), 781-790.
- Manfreda, S., Nardi, F., Samela, C., Grimaldi, S., Taramasso, A. C., Roth, G., and Sole, A., (2014). Investigation on the Use of Geomorphic Approaches for the Delineation of Flood Prone Areas, Journal of Hydrology, 517, 863-876. 
- Manfreda, S., Samela, C., Gioia, A., Consoli, G., Iacobellis, V.,  Giuzio, L., Cantisani, A., Sole, A., (2015). Flood-Prone Areas Assessment Using Linear Binary Classifiers based on flood maps obtained from 1D and 2D hydraulic models,  Natural Hazards, 79 (2), 735-754.
- Samela, C., Manfreda, S., De Paola, F., Giugni, M., Sole, A., Fiorentino, M., (2016). DEM-based approaches for the delineation of flood prone areas in an ungauged basin in Africa,  Journal of Hydrologic Engineering, 21(2).
- Samela, C., Troy, T.J., Manfreda, S., (2017). Geomorphic classifiers for flood-prone areas delineation for data-scarce environments, Advances in Water Resources, 102: 13-28 
- **Samela, C., Albano, R., Sole, A., Manfreda, S. (2018). Geomorphic Flood Area (GFA): a QGIS tool for a cost-effective delineation of the flood-prone areas, Computers, Environment and Urban Systems, (doi: 10.1016/j.compenvurbsys.2018.01.013)**

## Bug Report
The best place to file bug reports is at the [Bug Tracker](https://github.com/HydroLAB-UNIBAS/GFA-Geomorphic-Flood-Area/issues), this requires a free [GitHub](https://github.com/) account.

Please ensure that bug reports clearly describe the bug and if possible provide a simple script that can reproduce the problem. If in doubt contact [Raffaele Albano](http://www2.unibas.it/raffaelealbano/?page_id=115).
