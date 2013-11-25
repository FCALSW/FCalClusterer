/* =========================================================================
   GENERAL NOTES:
============================================================================
  (1). Description:
--------------------------------
   - SOME DESCRIPTION HERE ......
============================================================================ */
// clustering options
#define	_CLUSTER_MIDDLE_RANGE_ENGY_HITS 1
#define	_MOLIERE_RADIUS_CORRECTIONS 1
#define	_CLUSTER_MIXING_ENERGY_CORRECTIONS 1


// verbisity
#define _GENERAL_CLUSTERER_DEBUG 1
#define _CLUSTER_BUILD_DEBUG 0
#define _VIRTUALCLUSTER_BUILD_DEBUG 0
#define _MOL_RAD_CORRECT_DEBUG 0

#include "MarlinLumiCalClusterer.h"

#include "LumiCalClusterer.h"
#include "LumiCalClusterer_auxiliary.h"
#include "LumiCalClusterer_getCalHits.h"
#include "LumiCalClusterer_buildClusters.h"
#include "LumiCalClusterer_buildClusters_auxiliary.h"
#include "LumiCalClusterer_clusterMerger.h"
#include "LumiCalClusterer_fiducialVolumeCuts.h"
#include "LumiCalClusterer_energyCorrections.h"



/* ============================================================================
   Constructor
========================================================================= */
LumiCalClustererClass::LumiCalClustererClass(string lumiNameNow){

	_lumiName = lumiNameNow;

}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
========================================================================= */
void LumiCalClustererClass::init(	map < TString , int >		GlobalParamI,
					map < TString , double >	GlobalParamD ){


	/* --------------------------------------------------------------------------
	   constants specific to this class
	-------------------------------------------------------------------------- */
	_armsToCluster.clear();
	_armsToCluster.push_back(-1);
	_armsToCluster.push_back(1);

	_clusterMinNumHits			= 15;
	_hitMinEnergy				= 5*1e-6;
	_zLayerThickness			= 4.5;
	_methodCM				= "Log";	// "Energy";
	_elementsPercentInShowerPeakLayer	= 0.04;
	_hitMinEnergy				= 1e-16;
	_nNearNeighbor				= 6;		// number of near neighbors to consider

	// the minimal energy to take into account in the initial clustering pass is
	// defined as _middleEnergyHitBoundFrac of the minimal energy that is taken into
	// account when computing weighted averages in the log' weighting method
	_middleEnergyHitBoundFrac = .01;


	/* --------------------------------------------------------------------------
	   constants set by: GlobalMethodsClass
	-------------------------------------------------------------------------- */
	_logWeightConst = GlobalParamD["LogWeightConstant"];
	_moliereRadius = GlobalParamD["MoliereRadius"];

	// minimal separation distance and energy (of either cluster) to affect a merge
	_minSeparationDistance = GlobalParamD["MinSeparationDist"];
	_minClusterEngyGeV = engySignalGeV(GlobalParamD["MinClusterEngy"], "SignalToGeV");


	_thetaContainmentBouds[0] = GlobalParamD["ThetaMin"];
	_thetaContainmentBouds[1] = GlobalParamD["ThetaMax"];

	_maxLayerToAnalyse = GlobalParamI["NumCellsZ"];
	_cellRMax = GlobalParamI["NumCellsR"];;
	_cellPhiMax = GlobalParamI["NumCellsPhi"];

	_zFirstLayer = GlobalParamD["ZStart"];
	_rMin = GlobalParamD["RMin"];
	_rMax = GlobalParamD["RMax"];

	_rCellLength = (_rMax - _rMin) / _cellRMax;
	_phiCellLength = TwoPi() / _cellPhiMax;



	/* --------------------------------------------------------------------------
	   Print out Parameters
	-------------------------------------------------------------------------- */
	#if _GENERAL_CLUSTERER_DEBUG == 1
	cout	<< endl << coutUnderLine << coutBlue << "Global parameters for LumiCalClustererClass:" << coutDefault << endl;
	cout	<< " _cellRMax: "	 		<< coutRed << _cellRMax 			<< coutDefault << endl
		<< " _cellPhiMax: "			<< coutRed << _cellPhiMax 			<< coutDefault << endl
		<< " _zFirstLayer: "			<< coutRed << _zFirstLayer 			<< coutDefault << endl
		<< " _zLayerThickness: "		<< coutRed << _zLayerThickness 			<< coutDefault << endl
		<< " _rMin: "				<< coutRed << _rMin 				<< coutDefault << endl
		<< " _rMax: "				<< coutRed << _rMax 				<< coutDefault << endl
		<< " _rCellLength [mm]: "		<< coutRed << _rCellLength 			<< coutDefault << endl
		<< " _phiCellLength [rad]:"		<< coutRed << _phiCellLength 			<< coutDefault << endl
		<< " _methodCM: "			<< coutRed << _methodCM  			<< coutDefault << endl
		<< " _logWeightConst: "			<< coutRed << _logWeightConst			<< coutDefault << endl
		<< " _elementsPercentInShowerPeakLayer: "
							<< coutRed << _elementsPercentInShowerPeakLayer	<< coutDefault << endl
		<< " _moliereRadius: "			<< coutRed << _moliereRadius 			<< coutDefault << endl
		<< " _minSeparationDistance: "		<< coutRed << _minSeparationDistance 		<< coutDefault << endl
		<< " _minClusterEngy - GeV: "		<< coutRed << _minClusterEngyGeV		<< coutDefault << endl
		<< " _minClusterEngy - Signal: "	<< coutRed << engySignalGeV(_minClusterEngyGeV, "GeVToSignal")
													<< coutDefault << endl
		<< " _hitMinEnergy: "			<< coutRed << _hitMinEnergy			<< coutDefault << endl
		<< " _thetaContainmentBouds[0]: "	<< coutRed << _thetaContainmentBouds[0]		<< coutDefault << endl
		<< " _thetaContainmentBouds[1]: "	<< coutRed << _thetaContainmentBouds[1]		<< coutDefault << endl
		<< " _middleEnergyHitBoundFrac: "	<< coutRed << _middleEnergyHitBoundFrac 	<< coutDefault << endl
		<< endl;
	#endif



}


/* ============================================================================
   main actions in each event:
========================================================================= */
void LumiCalClustererClass::processEvent( LCEvent * evt ) {

	// increment / initialize global variables
	_totEngyArm[-1] = _totEngyArm[1] = 0.;


	int	numArmsToCluster, numSuperClusters;

	map < int , map < int , vector <CalorimeterHitImpl*> > >	calHits;

	map < int , map < int , vector<double> > > 	superClusterCM;
	map < int , vector<double> > :: iterator	superClusterCMIterator;

	map < int , map < int , CalorimeterHitImpl* > > calHitsCellIdGlobal;
	map < int , CalorimeterHitImpl* > :: iterator 	calHitsCellIdGlobalIterator;


	_superClusterIdToCellId.clear();
	_superClusterIdToCellEngy.clear();

	/* --------------------------------------------------------------------------
	   Loop over al hits in the LCCollection and write the hits into vectors
	   of CalorimeterHitImpl. Hits are split in two vectors, one for each arm
	   of LumiCal. 
	-------------------------------------------------------------------------- */
	getCalHits(evt , &(calHits));


	/* --------------------------------------------------------------------------
	   cccccccccccccc
	-------------------------------------------------------------------------- */
	numArmsToCluster = _armsToCluster.size();
	for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
		int armNow = _armsToCluster[armToClusterNow];

		#if _GENERAL_CLUSTERER_DEBUG == 1
		cout	<< endl << coutUnderLine << coutGreen << "ARM = " << armNow << " : " << coutDefault << endl << endl;
		#endif

		/* --------------------------------------------------------------------------
		   Construct clusters for each arm
		-------------------------------------------------------------------------- */
		#if _GENERAL_CLUSTERER_DEBUG == 1
		cout	<< coutBlue << "\tRun LumiCalClustererClass::buildClusters()" << coutDefault << endl;
		#endif

		buildClusters( calHits[armNow],
			       &(calHitsCellIdGlobal[armNow]),
			       &(_superClusterIdToCellId[armNow]),
			       &(_superClusterIdToCellEngy[armNow]),
			       &(superClusterCM[armNow]) );


		/* --------------------------------------------------------------------------
		   Merge superClusters according the minDistance nad minEngy rules
		-------------------------------------------------------------------------- */
		#if _GENERAL_CLUSTERER_DEBUG == 1
		cout	<< coutBlue << "\tRun LumiCalClustererClass::clusterMerger()" << coutDefault << endl;
		#endif

		clusterMerger( 		&(_superClusterIdToCellEngy[armNow]),
					&(_superClusterIdToCellId[armNow]),
					&(superClusterCM[armNow]),
					calHitsCellIdGlobal[armNow] );


		/* --------------------------------------------------------------------------
		   Perform fiducial volume cuts
		-------------------------------------------------------------------------- */
		fiducialVolumeCuts( 	&(_superClusterIdToCellId[armNow]),
					&(_superClusterIdToCellEngy[armNow]),
					&(superClusterCM[armNow]) );


		/* --------------------------------------------------------------------------
		   Perform energy correction for inter-mixed superClusters
		-------------------------------------------------------------------------- */
		#if _CLUSTER_MIXING_ENERGY_CORRECTIONS == 1
		numSuperClusters = superClusterCM[armNow].size();
		if(numSuperClusters == 2) {
			#if _GENERAL_CLUSTERER_DEBUG == 1
			cout	<< coutBlue << "\tRun LumiCalClustererClass::energyCorrections()" << coutDefault << endl;
			#endif

			energyCorrections( &(_superClusterIdToCellId[armNow]),
					   &(_superClusterIdToCellEngy[armNow]),
					   &(superClusterCM[armNow]),
					   calHitsCellIdGlobal[armNow] );
		}
		#endif

	}


	/* --------------------------------------------------------------------------
	   verbosity
	-------------------------------------------------------------------------- */
	#if _GENERAL_CLUSTERER_DEBUG == 1
	cout	<< endl << coutGreen << "Final clusters:" << coutDefault << endl;

	for(int armToClusterNow = 0; armToClusterNow < numArmsToCluster; armToClusterNow++) {
		int armNow = _armsToCluster[armToClusterNow];

		superClusterCMIterator = superClusterCM[armNow].begin();
		numSuperClusters       = superClusterCM[armNow].size();
		for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
			int superClusterId = (int)(*superClusterCMIterator).first;

			cout	<< "\t arm , Id " << armNow << "  ,  " << superClusterId
				<< "  \t energy " << superClusterCM[armNow][superClusterId][0] 
				<< "     \t pos(x,y) =  ( " << superClusterCM[armNow][superClusterId][1]
				<< " , " << superClusterCM[armNow][superClusterId][2] << " )"
				<< "     \t pos(theta,phi) =  ( " << superClusterCM[armNow][superClusterId][6]
				<< " , " << superClusterCM[armNow][superClusterId][7] << " )"
				<< coutDefault << endl;
		}
	}

	cout	<<  endl;
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






