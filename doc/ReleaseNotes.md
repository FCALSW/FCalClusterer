# v01-01

* 2024-04-03 Andre Sailer ([PR#74](https://github.com/fcalsw/fcalclusterer/pull/74))
  - CI: update underlying stacks

* 2024-04-03 scott snyder ([PR#73](https://github.com/fcalsw/fcalclusterer/pull/73))
  - Always record the output collection, even they are empty.  Avoids crashes from PodioOutput.

# v01-00-06

* 2023-07-04 Andre Sailer ([PR#72](https://github.com/FCALSW/FCalClusterer/pull/72))
  - MarlinLumiCalClusterer: if used with SimCalorimeterHitCollection as input, the CalorimeterHitCollection will always be created correctly and stored in the event.

# v01-00-05

* 2023-03-07 Andre Sailer ([PR#71](https://github.com/fcalsw/fcalclusterer/pull/71))
  - LumiCalReco: always write out collections, even if they are empty. Needed for EDM4hep compatibility, collections always have to be present

# v01-00-04

* 2021-11-10 Placido Fernandez Declara ([PR#69](https://github.com/FCALSW/FCalClusterer/pull/69))
  * BeamCalCluster Processor:  Always write the output collections, even when they are empty. Fixes issue with EDM4hep always expecting collections to be present
  * CMake: add LCIO include directory early to the list to avoid picking up outdated LCIO from the stack, mostly an issue for the CI tests

# v01-00-03

* 2020-09-18 Joseph C Wang ([PR#67](https://github.com/FCALSW/FCalClusterer/pull/67))
  - Fix compiler error for llvm (static constexpr) in BeamCalBackgroundPregen

# v01-00-02

* 2020-09-04 Andre Sailer ([PR#65](https://github.com/FCALSW/FCalClusterer/pull/65))
  - CMake: Fix hardcoding of CMAKE_INSTALL_PREFIX, now other installation locations can be set as expected

* 2020-04-23 Andre Sailer ([PR#64](https://github.com/FCALSW/FCalClusterer/pull/64))
  - Update CI to lcg97, cc7, clang10, gcc9 (from #63)
  - Fix Warnings from clang
  - Fix issues with c++17 (random_shuffle), also make shuffling of background files deterministic

# v01-00-01

* 2018-05-16 Shaojun Lu ([PR#59](https://github.com/fcalsw/fcalclusterer/pull/59))
  - added "FIND_PACKAGE( Doxygen REQUIRED )" for doc/CMakeLists.txt.
    - to build document with available Doxygen with "INSTALL_DOC=ON".

* 2018-05-15 Andre Sailer ([PR#58](https://github.com/fcalsw/fcalclusterer/pull/58))
  -  CI: Add build on CC7
  - CMAKE: fix linker problem on CC7 (newer ld version)

* 2018-05-15 Andre Sailer ([PR#57](https://github.com/fcalsw/fcalclusterer/pull/57))
  - Doc: add doxygen file

* 2018-05-14 Andre Sailer ([PR#56](https://github.com/fcalsw/fcalclusterer/pull/56))
  - CMAKE: use TARGET_[INCLUDE_DIRECTORIES|LINK_LIBRARIES|COMPILE_DEFINITIONS] instead of global include or linking

* 2018-05-08 Andre Sailer ([PR#55](https://github.com/fcalsw/fcalclusterer/pull/55))
  - CMake: let include directories be decided via dependencies

# v01-00

* 2017-11-16 Yorgos Voutsinas ([PR#30](https://github.com/FCALSW/FCalClusterer/pull/30))
  * LumiCalReco: Correcting the units when reading the thickness of the LumiCal layer
  * LumiCalReco: Correcting position reconstruction of the ZCell

* 2017-11-17 Andre Sailer ([PR#29](https://github.com/FCALSW/FCalClusterer/pull/29))
  - LumiCalClusterer: Assign hits to the created cluster, fixes #8 
    - Instead of the simulated hits collection a Digitized Hits collection for the LumiCal should be passed. If a SimCalorimterHit collection is passed dummy digitization is done and the CalorimeterHit collection is written out, controlled by Parameter `LumiCal_Hits` with default value `LumiCalHits`

* 2017-11-20 Andre Sailer ([PR#31](https://github.com/FCALSW/FCalClusterer/pull/31))
  - LumiCalReco: base all calculations on calibrated hits, allows "external" calibration
  - LumiCalReco: correct the cluster position calculation for optimal theta resolutions, ensure consistency between theta calculation and cluster position (x,y,z). Theta calculation for theta in the rootOutput was correct, the theta for SLCIO output is now much better reconstructed.
  - LumiCalReco: add RCellOffset to calculate radial cell coordinate for DD4hep segmentation, automatically derived from segmentation (affected the ROOT output for Theta)
  - LumiCalReco: re-factoring for weight calculations, cleanup, debug outputs

* 2017-11-22 Andre Sailer ([PR#32](https://github.com/FCALSW/FCalClusterer/pull/32))
  - LumiCalClusterer: scale weights by RMin/RCell to reduce polar angle bias
  - LumiCalClusterer: fix azimuthal angle offset when taking segmentation from DD4hep (actually no longer used, all (I hope...) position calculations are now based on the position of the hits)
  - LumiCalClusterer: Refactor ClusterClass, now using MCInfo and LCCluster to hold information, remove duplicated cluster position calculations

* 2017-11-27 Marko Petric ([PR#34](https://github.com/FCALSW/FCalClusterer/pull/34))
  - Add missing `#include <functional>` for `std::placeholders`

* 2018-01-19 Andre Sailer ([PR#35](https://github.com/FCALSW/FCalClusterer/pull/35))
  - Add automatic test for LumiCal reconstruction, some re-factoring to run tests without instantiating geometry
  - Fix bug for cluster at negative Z, where cluster at 0/180 degrees might have been split into 2
  - Optionally add cling-tidy and fix some issues it pointed out
  - Fix a bug in the finding of the max neighbouring pad that was introduced long time ago in the initial import and clean-up
  - Fix bugs in the calculation of the position at the LumiCal front used for MCParticle <-> cluster matching

* 2018-01-24 Andre Sailer ([PR#36](https://github.com/FCALSW/FCalClusterer/pull/36))
  - Add linking of clusters and hits, fixes #9 , one can use CalorimeterHit or SimCalorimeterHit collections
  - Added the possibility to run the BeamCal reconstruction on the LumiCal, which gives more reasonable results in the presence of background and multiple clusters.
  - New parameters for BeamCalClusterReco: 
     * MaxPadDistance: limit the distance between primary tower and other towers for clustering),  
     * LogWeightingConstant: if larger 0 uses logarithmic weighting for calculating shower position, like the LumiCal reco, 
     * DetectorStartingLayerID: The BeamCal driver starts layers at 1 , the LumiCal driver at 0
     * DetectorName: the name of the detector in the compactXML file this processor instance is for
  - New background method: Empty, allows running without (pair) background file
  - Fixes for position calculation of cluster in BeamCalClusterReco

* 2018-01-27 Andre Sailer ([PR#38](https://github.com/FCALSW/FCalClusterer/pull/38))
  - Fix DrawBeamCalFromDD4hep processor, fixes #37 , removed superfluous and incompatible geometry checks

* 2018-01-29 StrahinjaLukic ([PR#39](https://github.com/FCALSW/FCalClusterer/pull/39))
  - Optional feature was added to DrawBeamCalFromDD4hep to draw hit energy density per unit area rather than hit energy. This option is steerable via the new parameter `DrawDensities`
  - The z-axis range is calculated automatically in the polar segmentation plot of DrawBeamCalFromDD4hep.
  - A simple legend is drawn in the polar segmentation plot of DrawBeamCalFromDD4hep to indicate the dynamic range of the z-axis.
  - Minor style changes were made to the polar segmentation plot of DrawBeamCalFromDD4hep.

* 2018-02-01 Frank Gaede ([PR#41](https://github.com/FCALSW/FCalClusterer/pull/41))
  - add missing header cmath to TestUtilities.hh for gcc49 and gcc54

* 2018-02-02 Andre Sailer ([PR#43](https://github.com/FCALSW/FCalClusterer/pull/43))
  - BeamCalReco: fill subclusterEnergy to given ID via the ProcessorParameter: `SubClusterEnergyID`
  - BeamCalReco/LumiCalReco: set charge to 0 so to not confuse the RecoMCTruth linker

* 2018-02-02 Andre Sailer ([PR#42](https://github.com/FCALSW/FCalClusterer/pull/42))
  - Cleanup includes with "Include what you use", remove unused ones, explicitly add all used ones
  - Remove some pragma warnings for implemented features
  - Replace NULL with nullptr
  - Adapt to changes proposed in ilcsoft/LCIO#35, use LCIO namespaces

* 2018-02-06 Marko Petric ([PR#44](https://github.com/FCALSW/FCalClusterer/pull/44))
  - Include `DD4hepUnits.h` from `DD4hep` and not from `DDParsers` (which was renamed to `Evaluator`)

* 2018-03-08 Andre Sailer ([PR#46](https://github.com/FCALSW/FCalClusterer/pull/46))
  - ReadBeamCal: The processor can now also be used to read backgrounds for the LumiCal detector. New processor Parameters `DetectorName` and `DetectorStartingLayerID`

* 2018-03-16 Andre Sailer ([PR#50](https://github.com/FCALSW/FCalClusterer/pull/50))
  -  BeamCalClusterReco: fix problem when using CalorimeterHit collection as input
       * add optional parameter ReadoutName to be set in that case
  - BeamCalClusterReco: fix setting of CalorimeterHit collection flags , fixes #49 , where the cluster was only having nullptrs when the collection was written out
  - BeamCalClusterReco: only add hits to cluster if they don't come from background
  - BeamCalClusterReco: Pregenerated background, remove explicit gear dependency

* 2018-03-27 Andre Sailer ([PR#52](https://github.com/FCALSW/FCalClusterer/pull/52))
  - BeamCalClusterReco: add tree for efficiency post processing, and efficiency objects to get fake rates for different fake cluster energy ranges
  - BeamCalClusterReco: Protect against crash when input collection is empty, fixes #51 
  - BeamCalClusterReco PrintThisEvent option: make compatible with LumiCal detector
  - ReadBeamCal: fix for using LumiCal (only half the LumiCal was properly read, pads with sector<0 were rejected)

* 2018-03-28 Andre Sailer ([PR#54](https://github.com/FCALSW/FCalClusterer/pull/54))
  - LumiCalReco: fix setting of flags for CalorimeterHit collection (ID1, LONG) instead of copying from SimCalorimeterHit input collection, whose bits mean different things

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
