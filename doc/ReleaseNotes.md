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
