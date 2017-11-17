
/* =========================================================================
   GENERAL NOTES:
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION HERE ......
   ============================================================================ */
#include "LumiCalClusterer.h"

#include <IMPL/CalorimeterHitImpl.h>


namespace EVENT{
  class LCEvent;
}

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>



/* ============================================================================
   Constructor
   ========================================================================= */
LumiCalClustererClass::LumiCalClustererClass(std::string const& lumiNameNow)
    : _superClusterIdToCellId(),
      _superClusterIdToCellEngy(),
      _superClusterIdClusterInfo(),
      _lumiName(lumiNameNow),
      _clusterMinNumHits(15),
      _hitMinEnergy(5 * 1e-6),
      // global variables
      _numEventsPerTree(0),
      _resetRootTrees(0),
      _maxLayerToAnalyse(0),
      _zFirstLayer(0),
      _zLayerThickness(0.0),
      _zLayerPhiOffset(0.0),
      _rMin(0.0),
      _rMax(0.0),
      _rCellLength(0.0),
      _phiCellLength(0.0),
      _beamCrossingAngle(0.),
      _elementsPercentInShowerPeakLayer(0.03),
      _logWeightConst(0.0),
      _nNearNeighbor(6),
      _cellRMax(0),
      _cellPhiMax(0),
      _middleEnergyHitBoundFrac(0.01),
      _methodCM(GlobalMethodsClass::LogMethod),
      _moliereRadius(),
      _thetaContainmentBounds(),
      _minSeparationDistance(),
      _minClusterEngyGeV(),
      _totEngyArm(),
      _numHitsInArm(),
      _mydecoder(),
      _gmc(),
      _useDD4hep(false),
      RotMat() {}

/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void LumiCalClustererClass::init( GlobalMethodsClass const& gmc ){

  _gmc = gmc;
  /* --------------------------------------------------------------------------
     constants specific to this class
  _armsToCluster.clear();
  _armsToCluster.push_back(-1);
  _armsToCluster.push_back(1);
     -------------------------------------------------------------------------- */
  _methodCM = gmc.getMethod(gmc.GlobalParamS.at(GlobalMethodsClass::WeightingMethod));  // GlobalMethodsClass::LogMethod
  _clusterMinNumHits			= gmc.GlobalParamI.at(GlobalMethodsClass::ClusterMinNumHits); // = 15
  _hitMinEnergy				= gmc.GlobalParamD.at(GlobalMethodsClass::MinHitEnergy); // = 5e-6
  _zLayerThickness			= gmc.GlobalParamD.at(GlobalMethodsClass::ZLayerThickness); // = 4.5
  _zLayerPhiOffset			= gmc.GlobalParamD.at(GlobalMethodsClass::ZLayerPhiOffset); // = 3.75 [deg]
  _elementsPercentInShowerPeakLayer	= gmc.GlobalParamD.at(GlobalMethodsClass::ElementsPercentInShowerPeakLayer); // = 0.03  //APS 0.04;
  _nNearNeighbor			= gmc.GlobalParamI.at(GlobalMethodsClass::NumOfNearNeighbor); // = 6; // number of near neighbors to consider
  _beamCrossingAngle                    = gmc.GlobalParamD.at(GlobalMethodsClass::BeamCrossingAngle)/2.;
  RotMat[-1]["cos"]                     = cos( - _beamCrossingAngle );
  RotMat[-1]["sin"]                     = sin( - _beamCrossingAngle );
  RotMat[ 1]["cos"]                     = cos( _beamCrossingAngle );
  RotMat[ 1]["sin"]                     = sin( _beamCrossingAngle );

  // the minimal energy to take into account in the initial clustering pass is
  // defined as _middleEnergyHitBoundFrac of the minimal energy that is taken into
  // account when computing weighted averages in the log' weighting method
  _middleEnergyHitBoundFrac = gmc.GlobalParamD.at(GlobalMethodsClass::MiddleEnergyHitBoundFrac);                     // =.01;


  /* --------------------------------------------------------------------------
     constants set by: GlobalMethodsClass
     -------------------------------------------------------------------------- */
  _logWeightConst = gmc.GlobalParamD.at(GlobalMethodsClass::LogWeightConstant);
  _moliereRadius  = gmc.GlobalParamD.at(GlobalMethodsClass::MoliereRadius);

  // minimal separation distance and energy (of either cluster) to affect a merge
  _minSeparationDistance = gmc.GlobalParamD.at(GlobalMethodsClass::MinSeparationDist);
  _minClusterEngyGeV = gmc.GlobalParamD.at(GlobalMethodsClass::MinClusterEngyGeV);


  _thetaContainmentBounds[0] = gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMin);
  _thetaContainmentBounds[1] = gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMax);

  _maxLayerToAnalyse = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsZ);
  _cellRMax	   = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsR);
  _cellPhiMax	   = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsPhi);

  _zFirstLayer = gmc.GlobalParamD.at(GlobalMethodsClass::ZStart);
  _rMin	     = gmc.GlobalParamD.at(GlobalMethodsClass::RMin);
  _rMax	     = gmc.GlobalParamD.at(GlobalMethodsClass::RMax);

  _rCellLength   = gmc.GlobalParamD.at(GlobalMethodsClass::RCellLength);
  _phiCellLength = gmc.GlobalParamD.at(GlobalMethodsClass::PhiCellLength);

  _useDD4hep = gmc.isUsingDD4hep();

  /* --------------------------------------------------------------------------
     Print out Parameters
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  streamlog_out(MESSAGE) << std::endl << "Global parameters for LumiCalClustererClass:"  << std::endl;
  streamlog_out(MESSAGE) << " _cellRMax: "			    << _cellRMax			 << std::endl
			 << " _cellPhiMax: "			    << _cellPhiMax			 << std::endl
			 << " _zFirstLayer: "		            << _zFirstLayer			 << std::endl
			 << " _zLayerThickness: "		    << _zLayerThickness			 << std::endl
			 << " _zLayerPhiOffset[deg]: "	            << _zLayerPhiOffset*180./M_PI	 << std::endl
			 << " _rMin: "			            << _rMin				 << std::endl
			 << " _rMax: "			            << _rMax				 << std::endl
			 << " _rCellLength [mm]: "		    << _rCellLength			 << std::endl
			 << " _phiCellLength [rad]:"		    << _phiCellLength			 << std::endl
			 << " _methodCM: "			    << _methodCM			 << std::endl
			 << " _logWeightConst: "		    << _logWeightConst			 << std::endl
			 << " _elementsPercentInShowerPeakLayer: "  << _elementsPercentInShowerPeakLayer << std::endl
			 << " _moliereRadius: "		            << _moliereRadius			 << std::endl
			 << " _minSeparationDistance: "	            << _minSeparationDistance		 << std::endl
			 << " _minClusterEngy - GeV: "	            << _minClusterEngyGeV		 << std::endl
			 << " _hitMinEnergy: "		            << _hitMinEnergy		         << std::endl
			 << " _thetaContainmentBounds[0]: "	    << _thetaContainmentBounds[0]	 << std::endl
			 << " _thetaContainmentBounds[1]: "	    << _thetaContainmentBounds[1]	 << std::endl
			 << " _middleEnergyHitBoundFrac: "	    << _middleEnergyHitBoundFrac	 << std::endl
			 << " Clustering Options : " << std::endl
			 << "           _CLUSTER_MIDDLE_RANGE_ENGY_HITS    " <<  _CLUSTER_MIDDLE_RANGE_ENGY_HITS    << std::endl
			 << "           _MOLIERE_RADIUS_CORRECTIONS        " <<  _MOLIERE_RADIUS_CORRECTIONS        << std::endl
			 << "           _CLUSTER_MIXING_ENERGY_CORRECTIONS " <<  _CLUSTER_MIXING_ENERGY_CORRECTIONS << std::endl
	    << std::endl;
#endif



}


/* ============================================================================
   main actions in each event:
   ========================================================================= */
int LumiCalClustererClass::processEvent( EVENT::LCEvent * evt ) {
  int OK = 1;
  int NOK = 0;
  // increment / initialize global variables
  _totEngyArm[-1] = _totEngyArm[1] = 0.;
  _numHitsInArm[-1] = _numHitsInArm[1] = 0;

  MapIntMapIntVCalHit calHits;
  MapIntMapIntLCCluster superClusterCM;

  MapIntMapIntCalHit calHitsCellIdGlobal;

  _superClusterIdToCellId.clear();
  _superClusterIdToCellEngy.clear();
  _superClusterIdClusterInfo.clear();

  /* --------------------------------------------------------------------------
     Loop over all hits in the LCCollection and write the hits into std::vectors
     of IMPL::CalorimeterHitImpl. Hits are split in two std::vectors, one for each arm
     of LumiCal.
     -------------------------------------------------------------------------- */
  if ( !getCalHits(evt , calHits) ) return NOK;


  /* --------------------------------------------------------------------------
     cccccccccccccc
     --------------------------------------------------------------------------
  const int numArmsToCluster = _armsToCluster.size();
    int armNow = _armsToCluster[armToClusterNow];
 */
  for(int armNow = -1; armNow < 2; armNow += 2) {
    streamlog_out(DEBUG6) << std::endl
	      << "ARM = " << armNow << " : " << std::endl << std::endl;
    /* --------------------------------------------------------------------------
       Construct clusters for each arm
       -------------------------------------------------------------------------- */
    streamlog_out(DEBUG6) << "\tRun LumiCalClustererClass::buildClusters()" << std::endl;
    streamlog_out(DEBUG5) << "\tEnergy deposit: "<< _totEngyArm[-1] << "\t" << _totEngyArm[1] <<"\n"
			  << "\tNumber of hits: "<< _numHitsInArm[-1] << "\t" << _numHitsInArm[1] << std::endl;
    buildClusters( calHits[armNow],
		   calHitsCellIdGlobal[armNow],
		   _superClusterIdToCellId[armNow],
		   _superClusterIdToCellEngy[armNow],
		   superClusterCM[armNow],
		   armNow);


    /* --------------------------------------------------------------------------
       Merge superClusters according the minDistance nad minEngy rules
       -------------------------------------------------------------------------- */
    streamlog_out(DEBUG6) << "\tRun LumiCalClustererClass::clusterMerger()" << std::endl;
    clusterMerger(		_superClusterIdToCellEngy[armNow],
				_superClusterIdToCellId[armNow],
				superClusterCM[armNow],
				calHitsCellIdGlobal[armNow] );


    /* --------------------------------------------------------------------------
       Perform fiducial volume cuts
       -------------------------------------------------------------------------- */
    streamlog_out(DEBUG6) << "\tRun LumiCalClustererClass::fiducialVolumeCuts()" << std::endl;
    fiducialVolumeCuts(		_superClusterIdToCellId[armNow],
				_superClusterIdToCellEngy[armNow],
				superClusterCM[armNow] );


    /* --------------------------------------------------------------------------
       Perform energy correction for inter-mixed superClusters
       -------------------------------------------------------------------------- */
#if _CLUSTER_MIXING_ENERGY_CORRECTIONS == 1
    if(superClusterCM[armNow].size() == 2) {
      streamlog_out(DEBUG6) << "Run LumiCalClustererClass::energyCorrections()" << std::endl;
      energyCorrections( _superClusterIdToCellId[armNow],
			 _superClusterIdToCellEngy[armNow],
			 superClusterCM[armNow],
			 calHitsCellIdGlobal[armNow] );
     }
#endif

  }


  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
  streamlog_out(DEBUG6) << "Final clusters:" << std::endl;
  for (int armNow = -1; armNow < 2; armNow += 2) {
    for (MapIntLCCluster::iterator superClusterCMIterator = superClusterCM[armNow].begin();
         superClusterCMIterator != superClusterCM[armNow].end(); ++superClusterCMIterator) {
      const int superClusterId = (int)superClusterCMIterator->first;

      streamlog_out(DEBUG6) << "  Arm:"    << std::setw(4) << armNow
                            << "  Id:"     << std::setw(4) << superClusterId
                            << superClusterCMIterator->second
                            << std::endl;
      _superClusterIdClusterInfo[armNow][superClusterId] = superClusterCMIterator->second;
    }
  }

  return OK;

}
