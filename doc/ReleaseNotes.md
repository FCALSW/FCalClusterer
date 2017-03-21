# v04-00

## PR #14 from andresailer 
 
* Replace deprecated calls from dd4hep
* Fix warnings for GCC 

## PR #13 from andresailer 
 
* BeamCalReco: Fix bug causing printouts of integration errors from ROOT  

## PR #11 from andresailer 
 
* LumiCalReco: allow the use or arbitrary processor names 

## PR #6 from andresailer 
 
* LumiCalReco: corrected the cell neighbor selection to pick the correct neighbors in phi and r 

## PR #5 from andresailer 
 
* LumiCalReco can now use geometry from DD4hep 

## PR #3 from andresailer 
 
* LumiCalReco: fix bug in LCCluster setZ function 

## PR #2 from pawlikb 
 
* Fix cellID Issue in LumiCal 

BeamCalReco:
* Implemented compatibility with DD4hep Geometry


# v0-03

BeamCalReco (A.Sailer):
* Reduce memory footprint and fix memory leaks for Parameterized background method
* Silence warnigns from root integral and getrandom, enabled for DEBUG0
* Compatible with ROOT 6

LumiCalReco (B.Pawlik):
* Massive fixes, mainly coordinate system, cluster merging. Creating ClustersClass Map and root output is optional now
* Fix LCIO interface to pass hit positions, geo info. Create LucasGear for Marlin reco, and more...


# v0-02 

LumiCalReco:
* Bugfixes/improvements by Bogdan Pawlik
  - Read Parameters from Gear
  - More steering parameters from xml
  - Fix to create the cluster on the proper side

BeamCalReco:
* Bugfixes/improvements from Andrey Sapronov
  - BackgroundMethod Parametrised renamed to Gaussian
  - new BackgroundMethod Parametrised with Gaus/x parametrisation

* More examples in the example folder

* Bugfixes from Andre Sailer


# v0-01

First release of FCalClusterer

* BeamCalReco
* LumiCalReco
