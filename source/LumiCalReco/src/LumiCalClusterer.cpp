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
LumiCalClustererClass::LumiCalClustererClass(std::string const& lumiNameNow):
  _superClusterIdToCellId(),
  _superClusterIdToCellEngy(),
  _superClusterIdClusterInfo(),
  _lumiName( lumiNameNow ),
  _clusterMinNumHits(15),
  _hitMinEnergy(5*1e-6),
  // global variables
  _numEventsPerTree(0), _resetRootTrees(0),
  _maxLayerToAnalyse(0),
  _zFirstLayer(0), _zLayerThickness(0.0), _rMin(0.0), _rMax(0.0), _rCellLength(0.0), _phiCellLength(0.0),
  _elementsPercentInShowerPeakLayer(0.03),
  _logWeightConst(0.0),
  _nNearNeighbor (6),
  _cellRMax(0), _cellPhiMax (0),
  _middleEnergyHitBoundFrac(0.01),
  _methodCM("LogMethod"),
  _moliereRadius(),
  _thetaContainmentBounds(),
  _minSeparationDistance(), _minClusterEngyGeV(), _minClusterEngySignal(),
  _totEngyArm(),
  //  _armsToCluster(),
  _mydecoder(NULL)
{
}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void LumiCalClustererClass::init( GlobalMethodsClass::ParametersString const& GlobalParamS,
				  GlobalMethodsClass::ParametersInt    const& GlobalParamI,
				  GlobalMethodsClass::ParametersDouble const& GlobalParamD ){


  /* --------------------------------------------------------------------------
     constants specific to this class
  _armsToCluster.clear();
  _armsToCluster.push_back(-1);
  _armsToCluster.push_back(1);
     -------------------------------------------------------------------------- */
  _methodCM				= GlobalParamS.at(GlobalMethodsClass::WeightingMethod);                  // GlobalMethodsClass::LogMethod
  _clusterMinNumHits			= GlobalParamI.at(GlobalMethodsClass::ClusterMinNumHits);                // = 15
  _hitMinEnergy				= GlobalParamD.at(GlobalMethodsClass::MinHitEnergy);                     // = 5e-6
  _zLayerThickness			= GlobalParamD.at(GlobalMethodsClass::ZLayerThickness);                  // = 4.5
  _elementsPercentInShowerPeakLayer	= GlobalParamD.at(GlobalMethodsClass::ElementsPercentInShowerPeakLayer); // = 0.03  //APS 0.04;
  _nNearNeighbor			= GlobalParamI.at(GlobalMethodsClass::NumOfNearNeighbor);                // = 6; // number of near neighbors to consider

  // the minimal energy to take into account in the initial clustering pass is
  // defined as _middleEnergyHitBoundFrac of the minimal energy that is taken into
  // account when computing weighted averages in the log' weighting method
  _middleEnergyHitBoundFrac = GlobalParamD.at(GlobalMethodsClass::MiddleEnergyHitBoundFrac);                     // =.01;


  /* --------------------------------------------------------------------------
     constants set by: GlobalMethodsClass
     -------------------------------------------------------------------------- */
  _logWeightConst = GlobalParamD.at(GlobalMethodsClass::LogWeightConstant);
  _moliereRadius  = GlobalParamD.at(GlobalMethodsClass::MoliereRadius);

  // minimal separation distance and energy (of either cluster) to affect a merge
  _minSeparationDistance = GlobalParamD.at(GlobalMethodsClass::MinSeparationDist);
  _minClusterEngyGeV = GlobalParamD.at(GlobalMethodsClass::MinClusterEngyGeV);
  _minClusterEngySignal = GlobalParamD.at(GlobalMethodsClass::MinClusterEngySignal);


  _thetaContainmentBounds[0] = GlobalParamD.at(GlobalMethodsClass::ThetaMin);
  _thetaContainmentBounds[1] = GlobalParamD.at(GlobalMethodsClass::ThetaMax);

  _maxLayerToAnalyse = GlobalParamI.at(GlobalMethodsClass::NumCellsZ);
  _cellRMax	   = GlobalParamI.at(GlobalMethodsClass::NumCellsR);
  _cellPhiMax	   = GlobalParamI.at(GlobalMethodsClass::NumCellsPhi);

  _zFirstLayer = GlobalParamD.at(GlobalMethodsClass::ZStart);
  _rMin	     = GlobalParamD.at(GlobalMethodsClass::RMin);
  _rMax	     = GlobalParamD.at(GlobalMethodsClass::RMax);

  _rCellLength = (_rMax - _rMin) / _cellRMax;
  _phiCellLength = 2*M_PI / _cellPhiMax;



  /* --------------------------------------------------------------------------
     Print out Parameters
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  streamlog_out(MESSAGE) << std::endl << "Global parameters for LumiCalClustererClass:"  << std::endl;
  streamlog_out(MESSAGE) << " _cellRMax: "			    << _cellRMax			 << std::endl
	    << " _cellPhiMax: "			    << _cellPhiMax			 << std::endl
	    << " _zFirstLayer: "		    << _zFirstLayer			 << std::endl
	    << " _zLayerThickness: "		    << _zLayerThickness			 << std::endl
	    << " _rMin: "			    << _rMin				 << std::endl
	    << " _rMax: "			    << _rMax				 << std::endl
	    << " _rCellLength [mm]: "		    << _rCellLength			 << std::endl
	    << " _phiCellLength [rad]:"		    << _phiCellLength			 << std::endl
	    << " _methodCM: "			    << _methodCM			 << std::endl
	    << " _logWeightConst: "		    << _logWeightConst			 << std::endl
	    << " _elementsPercentInShowerPeakLayer: " << _elementsPercentInShowerPeakLayer	<< std::endl
	    << " _moliereRadius: "		    << _moliereRadius			 << std::endl
	    << " _minSeparationDistance: "	    << _minSeparationDistance		 << std::endl
	    << " _minClusterEngy - GeV: "	    << _minClusterEngyGeV		 << std::endl
	    << " _minClusterEngy - Signal: "	    << _minClusterEngySignal             << std::endl
	    << " _hitMinEnergy: "		    << _hitMinEnergy		         << std::endl
	    << " _thetaContainmentBounds[0]: "	    << _thetaContainmentBounds[0]	 << std::endl
	    << " _thetaContainmentBounds[1]: "	    << _thetaContainmentBounds[1]	 << std::endl
	    << " _middleEnergyHitBoundFrac: "	    << _middleEnergyHitBoundFrac	 << std::endl
	    << std::endl;
#endif



}


/* ============================================================================
   main actions in each event:
   ========================================================================= */
void LumiCalClustererClass::processEvent( EVENT::LCEvent * evt ) {

  // increment / initialize global variables
  _totEngyArm[-1] = _totEngyArm[1] = 0.;

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
  getCalHits(evt , calHits);


  /* --------------------------------------------------------------------------
     cccccccccccccc
     --------------------------------------------------------------------------
  const int numArmsToCluster = _armsToCluster.size();
    int armNow = _armsToCluster[armToClusterNow];
 */
  for(int armNow = -1; armNow < 2; armNow += 2) {
#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << std::endl
	      << "ARM = " << armNow << " : " << std::endl << std::endl;
#endif

    /* --------------------------------------------------------------------------
       Construct clusters for each arm
       -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << "\tRun LumiCalClustererClass::buildClusters()" << std::endl;
#endif

    buildClusters( calHits[armNow],
		   calHitsCellIdGlobal[armNow],
		   _superClusterIdToCellId[armNow],
		   _superClusterIdToCellEngy[armNow],
		   superClusterCM[armNow],
		   armNow);


    /* --------------------------------------------------------------------------
       Merge superClusters according the minDistance nad minEngy rules
       -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << "\tRun LumiCalClustererClass::clusterMerger()" << std::endl;
#endif

    clusterMerger(		_superClusterIdToCellEngy[armNow],
				_superClusterIdToCellId[armNow],
				superClusterCM[armNow],
				calHitsCellIdGlobal[armNow] );


    /* --------------------------------------------------------------------------
       Perform fiducial volume cuts
       -------------------------------------------------------------------------- */
    fiducialVolumeCuts(		_superClusterIdToCellId[armNow],
				_superClusterIdToCellEngy[armNow],
				superClusterCM[armNow] );


    /* --------------------------------------------------------------------------
       Perform energy correction for inter-mixed superClusters
       -------------------------------------------------------------------------- */
#if _CLUSTER_MIXING_ENERGY_CORRECTIONS == 1
    if(superClusterCM[armNow].size() == 2) {
#if _GENERAL_CLUSTERER_DEBUG == 1
      streamlog_out( DEBUG ) << "\tRun LumiCalClustererClass::energyCorrections()" << std::endl;
#endif

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
#if _GENERAL_CLUSTERER_DEBUG == 1
  streamlog_out( DEBUG4 ) << "Final clusters:" << std::endl;
  /*
 for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
     const int armNow = _armsToCluster[armToClusterNow];
  */
 for(int armNow = -1; armNow < 2; armNow += 2) {
     for(MapIntLCCluster::iterator superClusterCMIterator = superClusterCM[armNow].begin();
	superClusterCMIterator != superClusterCM[armNow].end();
	++superClusterCMIterator) {
      const int superClusterId = (int)superClusterCMIterator->first;

      streamlog_out( DEBUG4 ) << "  Arm:"    << std::setw(4)  << armNow
			     << "  Id:"     << std::setw(4)  << superClusterId
			     << superClusterCMIterator->second
			     << std::endl;

      //Store information of clusters
      _superClusterIdClusterInfo[armNow][superClusterId] = superClusterCMIterator->second;
    }
  }

#endif



  // cleanUp
  for (std::map<int, std::map < int, std::vector <IMPL::CalorimeterHitImpl*> > >::iterator it = calHits.begin(); it != calHits.end(); ++it) {
    std::map < int, std::vector <IMPL::CalorimeterHitImpl*> >& mapVecHits = it->second;
    for (std::map < int, std::vector <IMPL::CalorimeterHitImpl*> >::iterator it2 = mapVecHits.begin(); it2 != mapVecHits.end(); ++it2) {
      std::vector <IMPL::CalorimeterHitImpl*>& vecHits = it2->second;
      for (std::vector<IMPL::CalorimeterHitImpl*>::iterator it3 = vecHits.begin(); it3 != vecHits.end(); ++it3) {
	delete (*it3);
      }
    }
  }

  calHits.clear();
  superClusterCM.clear();
  calHitsCellIdGlobal.clear();


}
