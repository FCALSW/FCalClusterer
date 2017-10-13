# v00-06

* 2017-07-12 Andre Sailer ([PR#19](https://github.com/FCalSW/FCalClusterer/pull/19))
  - LumiCalClusterer: wrap decoder object in unique_ptr to guarantee cleanup, fixes memory leak

* 2017-07-17 Andre Sailer ([PR#20](https://github.com/FCalSW/FCalClusterer/pull/20))
  - LumiCalClusterer: fix memory leak if there are not enough hits to attempt clustering

* 2017-07-18 Andre Sailer ([PR#22](https://github.com/FCalSW/FCalClusterer/pull/22))
  - LumiCalReco: use lambda function in TF1 instead of string
  - BeamCalReco: GeometryDD layer starts at 0 if simulated with ddsim
  - BeamCalReco: fix infinite loop  in background subtraction if nBX = 0 (no background), fixes #21

* 2017-10-02 Andre Sailer ([PR#26](https://github.com/FCalSW/FCalClusterer/pull/26))
  - Adapt to changes in AidaSoft/DD4hep#238
  - CMake add C language to PROJECT for macs and cmake 3.3.2

* 2017-07-26 Andre Sailer ([PR#24](https://github.com/FCalSW/FCalClusterer/pull/24))
  - LumiCalReco: disable Fiducial Volume cuts, cuts on geometry should be done in later step of analysis

* 2017-07-26 Andre Sailer ([PR#23](https://github.com/FCalSW/FCalClusterer/pull/23))
  - BeamCalReco: correct the phi position of cluster reconstructed in the backward direction, which was incorrectly calculated due to the rotation of the backward BeamCal detector

* 2017-10-06 Andre Sailer ([PR#27](https://github.com/FCalSW/FCalClusterer/pull/27))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

# v00-05

* 2017-06-29 Andre Sailer ([PR#18](https://github.com/FCALSW/FCalClusterer/pull/18))
  - LumiCalReco: parameters are now a shared_ptr

* 2017-04-03 Andre Sailer ([PR#15](https://github.com/FCALSW/FCalClusterer/pull/15))
  * Fixed warnings for llvm

* 2017-06-20 Andre Sailer ([PR#17](https://github.com/FCALSW/FCalClusterer/pull/17))
  - Adapt to changes in namespaces and LCDD -->  Detector in DD4hep

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
