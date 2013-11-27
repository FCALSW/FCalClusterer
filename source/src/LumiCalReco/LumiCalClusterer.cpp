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
  _methodCM(GlobalMethodsClass::LogMethod),
  _moliereRadius(),
  _thetaContainmentBouds(),
  _minSeparationDistance(), _minClusterEngyGeV(),
  _totEngyArm(),
  _armsToCluster(),
  _mydecoder(NULL)
{
}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void LumiCalClustererClass::init( GlobalMethodsClass::ParametersInt    const& GlobalParamI,
				  GlobalMethodsClass::ParametersDouble const& GlobalParamD ){


  /* --------------------------------------------------------------------------
     constants specific to this class
     -------------------------------------------------------------------------- */
  _armsToCluster.clear();
  _armsToCluster.push_back(-1);
  _armsToCluster.push_back(1);
#pragma message ("Make Parameters steerable")
  _clusterMinNumHits			= 15;
  _hitMinEnergy				= 5*1e-6;
  _hitMinEnergy				= 1e-16;
  _zLayerThickness			= 4.5;
  _methodCM				= GlobalMethodsClass::LogMethod;	// "Energy";
  _elementsPercentInShowerPeakLayer	= 0.03;//APS 0.04;
  _nNearNeighbor				= 6;		// number of near neighbors to consider

  // the minimal energy to take into account in the initial clustering pass is
  // defined as _middleEnergyHitBoundFrac of the minimal energy that is taken into
  // account when computing weighted averages in the log' weighting method
  _middleEnergyHitBoundFrac = .01;


  /* --------------------------------------------------------------------------
     constants set by: GlobalMethodsClass
     -------------------------------------------------------------------------- */
  _logWeightConst = GlobalParamD.at(GlobalMethodsClass::LogWeightConstant);
  _moliereRadius  = GlobalParamD.at(GlobalMethodsClass::MoliereRadius);

  // minimal separation distance and energy (of either cluster) to affect a merge
  _minSeparationDistance = GlobalParamD.at(GlobalMethodsClass::MinSeparationDist);
  _minClusterEngyGeV = GlobalMethodsClass::SignalGevConversion( GlobalMethodsClass::Signal_to_GeV,
								GlobalParamD.at(GlobalMethodsClass::MinClusterEngy) );


  _thetaContainmentBouds[0] = GlobalParamD.at(GlobalMethodsClass::ThetaMin);
  _thetaContainmentBouds[1] = GlobalParamD.at(GlobalMethodsClass::ThetaMax);

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
  std::cout << std::endl << "Global parameters for LumiCalClustererClass:"  << std::endl;
  std::cout << " _cellRMax: "			    << _cellRMax			 << std::endl
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
	    << " _minClusterEngy - Signal: "	    
	    << GlobalMethodsClass::SignalGevConversion(GlobalMethodsClass::GeV_to_Signal, _minClusterEngyGeV)
	    << std::endl
	    << " _hitMinEnergy: "		    << _hitMinEnergy		         << std::endl
	    << " _thetaContainmentBouds[0]: "	    << _thetaContainmentBouds[0]	 << std::endl
	    << " _thetaContainmentBouds[1]: "	    << _thetaContainmentBouds[1]	 << std::endl
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


  int	numArmsToCluster, numSuperClusters;

  std::map < int , std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > >	calHits;

  std::map < int , std::map < int , LCCluster > >	superClusterCM;
  std::map < int , LCCluster > :: iterator	superClusterCMIterator;

  std::map < int , std::map < int , IMPL::CalorimeterHitImpl* > > calHitsCellIdGlobal;
  std::map < int , IMPL::CalorimeterHitImpl* > :: iterator	calHitsCellIdGlobalIterator;


  _superClusterIdToCellId.clear();
  _superClusterIdToCellEngy.clear();
  _superClusterIdClusterInfo.clear();

  /* --------------------------------------------------------------------------
     Loop over al hits in the LCCollection and write the hits into std::vectors
     of IMPL::CalorimeterHitImpl. Hits are split in two std::vectors, one for each arm
     of LumiCal.
     -------------------------------------------------------------------------- */
  getCalHits(evt , calHits);


  /* --------------------------------------------------------------------------
     cccccccccccccc
     -------------------------------------------------------------------------- */
  numArmsToCluster = _armsToCluster.size();
  for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
    int armNow = _armsToCluster[armToClusterNow];

#if _GENERAL_CLUSTERER_DEBUG == 1
    std::cout << std::endl
	      << "ARM = " << armNow << " : " << std::endl << std::endl;
#endif

    /* --------------------------------------------------------------------------
       Construct clusters for each arm
       -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
    std::cout << "\tRun LumiCalClustererClass::buildClusters()" << std::endl;
#endif

    buildClusters( calHits[armNow],
		   (calHitsCellIdGlobal[armNow]),
		   (_superClusterIdToCellId[armNow]),
		   (_superClusterIdToCellEngy[armNow]),
		   (superClusterCM[armNow]) );


    /* --------------------------------------------------------------------------
       Merge superClusters according the minDistance nad minEngy rules
       -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
    std::cout << "\tRun LumiCalClustererClass::clusterMerger()" << std::endl;
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
    numSuperClusters = superClusterCM[armNow].size();
    if(numSuperClusters == 2) {
#if _GENERAL_CLUSTERER_DEBUG == 1
      std::cout << "\tRun LumiCalClustererClass::energyCorrections()" << std::endl;
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
  std::cout << std::endl << "Final clusters:" << std::endl;

  for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
    int armNow = _armsToCluster[armToClusterNow];

    superClusterCMIterator = superClusterCM[armNow].begin();
    numSuperClusters       = superClusterCM[armNow].size();
    for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
      int superClusterId = (int)(*superClusterCMIterator).first;

      std::cout << "  Arm:"    << std::setw(4)  << armNow
		<< "  Id:"     << std::setw(4)  << superClusterId
		<< superClusterCM[armNow][superClusterId]
		<< std::endl;

      //Store information of clusters
      _superClusterIdClusterInfo[armNow][superClusterId] = superClusterCM[armNow][superClusterId];
    }
  }

  std::cout << std::endl;
#endif



  // cleanUp
  numArmsToCluster = _armsToCluster.size();
  for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
    int armNow = _armsToCluster[armToClusterNow];

    calHitsCellIdGlobalIterator = calHitsCellIdGlobal[armNow].begin();
    int numHitsLayer            = calHitsCellIdGlobal[armNow].size();
    for(int hitNow = 0; hitNow < numHitsLayer; hitNow++, calHitsCellIdGlobalIterator++) {
      int cellId = (int)(*calHitsCellIdGlobalIterator).first;  // Id of cell

      delete calHitsCellIdGlobal[armNow][cellId];
    }
  }

  superClusterCM.clear();
  calHitsCellIdGlobal.clear();


}
