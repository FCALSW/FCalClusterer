#ifndef LumiCalClusterer_h
#define LumiCalClusterer_h 1


// clustering options
#define	_CLUSTER_MIDDLE_RANGE_ENGY_HITS 1
#define	_MOLIERE_RADIUS_CORRECTIONS 1
#define	_CLUSTER_MIXING_ENERGY_CORRECTIONS 1


// verbosity
#define _GENERAL_CLUSTERER_DEBUG 0
#define _CREATE_CLUSTERS_DEBUG 0
#define _CLUSTER_BUILD_DEBUG 0
#define _VIRTUALCLUSTER_BUILD_DEBUG 0
#define _MOL_RAD_CORRECT_DEBUG 0

#include "Global.hh"

#include "GlobalMethodsClass.h"
#include "LCCluster.hh"

#include <UTIL/CellIDDecoder.h>

#include <string>
#include <memory>

enum RETVAL {
  NOK = 0,
  OK = 1,
};

namespace EVENT {
  class CalorimeterHit;
  class LCEvent;
  class LCCollection;
}

class LumiCalClustererClass {

public:

  // Constructor
  LumiCalClustererClass( std::string const& lumiNameNow ) ;

  // initialization routine - Called at the begining of the job.
  void init( GlobalMethodsClass const& gmc );

  /// set the cutOnFiducialVolume flag
  void setCutOnFiducialVolume( bool cutFlag ) { _cutOnFiducialVolume = cutFlag; }

  // main actions in each event -Called for every event - the working horse.
  RETVAL processEvent( EVENT::LCEvent * evt ) ;

  MapIntMapIntVInt       _superClusterIdToCellId;
  MapIntMapIntVDouble    _superClusterIdToCellEngy;
  MapIntMapIntLCCluster  _superClusterIdClusterInfo;

  void setLumiCollectionName(std::string const& lumiNameNow) { _lumiName = lumiNameNow; }
  void setLumiOutCollectionName(std::string const& name) { _lumiOutName = name; }

protected:

  LumiCalClustererClass(LumiCalClustererClass const& rhs);
  LumiCalClustererClass& operator=(LumiCalClustererClass const& rhs);

  // Processor Parameters
  std::string	_lumiName;
  std::string   _lumiOutName = "";
  int		_clusterMinNumHits;
  double	_hitMinEnergy;


  // global variables
  int	_numEventsPerTree, _resetRootTrees;
  int	_maxLayerToAnalyse;
  double	_zFirstLayer, _zLayerThickness, _zLayerPhiOffset, _rMin, _rMax, _rCellLength, _phiCellLength;
  double        _beamCrossingAngle;
  double	_elementsPercentInShowerPeakLayer;
  double	_logWeightConst;
  int	_nNearNeighbor ;
  int	_cellRMax, _cellPhiMax ;
  double	_middleEnergyHitBoundFrac;
  GlobalMethodsClass::WeightingMethod_t _methodCM;
  double	_moliereRadius;
  double	_thetaContainmentBounds[2];
  double	_minSeparationDistance, _minClusterEngyGeV;

  MapIntDouble _totEngyArm;
  MapIntInt    _numHitsInArm;
  //  VInt _armsToCluster;

  std::unique_ptr<UTIL::CellIDDecoder<EVENT::CalorimeterHit>> _mydecoder{};

  GlobalMethodsClass _gmc;

  bool _useDD4hep;
  bool _cutOnFiducialVolume=false;

  // methods:
  int	getCalHits( EVENT::LCEvent * evt,
		    MapIntMapIntVCalHit & calHits );

  EVENT::LCCollection* createCaloHitCollection(EVENT::LCCollection* simCaloHitCollection) const;

  int	buildClusters(	MapIntVCalHit const& calHits,
			MapIntCalHit & calHitsCellIdGlobal,
			MapIntVInt & superClusterIdToCellId,
			MapIntVDouble & superClusterIdToCellEngy,
			MapIntLCCluster & superClusterCM, 
			const int detectorArm);

  int	initialClusterBuild( MapIntCalHit const& calHitsCellId,
			     MapIntInt			  & cellIdToClusterId,
			     MapIntVInt	  & clusterIdToCellId,
			     MapIntLCCluster & clusterCM,
			     VInt const& controlVar );

  int	initialLowEngyClusterBuild( MapIntCalHit const& calHitsSmallEngyCellId,
				    MapIntCalHit & calHitsCellId,
				    MapIntInt			 & cellIdToClusterId,
				    MapIntVInt		 & clusterIdToCellId,
				    MapIntLCCluster	 & clusterCM );


  int	virtualCMClusterBuild( MapIntCalHit const&	  calHitsCellId,
			       MapIntInt	&			 cellIdToClusterId,
			       MapIntVInt & clusterIdToCellId,
			       MapIntLCCluster		& clusterCM,
			       MapIntVirtualCluster const& virtualClusterCM );

  int	virtualCMPeakLayersFix(	MapIntCalHit const&	calHitsCellId,
				MapIntInt				& cellIdToClusterId,
				MapIntVInt		& clusterIdToCellId,
				MapIntLCCluster		& clusterCM,
				MapIntVirtualCluster virtualClusterCM );

  int	buildSuperClusters ( MapIntCalHit & calHitsCellIdGlobal,
			     VMapIntCalHit const&	calHitsCellId,
			     VMapIntVInt const&	clusterIdToCellId,
			     VMapIntLCCluster const&	clusterCM,
			     VMapIntVirtualCluster const& virtualClusterCM,
			     MapIntInt & cellIdToSuperClusterId,
			     MapIntVInt & superClusterIdToCellId,
			     MapIntLCCluster & superClusterCM );

  int	engyInMoliereCorrections ( MapIntCalHit const& calHitsCellIdGlobal,
				   MapIntVCalHit const& calHits,
				   VMapIntCalHit const& calHitsCellIdLayer,
				   VMapIntVInt & clusterIdToCellId,
				   VMapIntLCCluster & clusterCM,
				   VMapIntInt & cellIdToClusterId,
				   MapIntInt & cellIdToSuperClusterId,
				   MapIntVInt & superClusterIdToCellId,
				   MapIntLCCluster & superClusterCM,
				   double middleEnergyHitBound,
				   int detectorArm );



  void	energyCorrections (	MapIntVInt & superClusterIdToCellId,
				MapIntVDouble & superClusterIdToCellEngy,
				MapIntLCCluster & superClusterCM,
				MapIntCalHit const& calHitsCellIdGlobal ) ;

  void clusterMerger(MapIntVDouble& clusterIdToCellEngy, MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                     MapIntCalHit& calHitsCellIdGlobal);

  void	fiducialVolumeCuts (	MapIntVInt & superClusterIdToCellId,
				MapIntVDouble & superClusterIdToCellEngy,
				MapIntLCCluster & superClusterCM ) ;


  void	getThetaPhiZCluster( MapIntCalHit const& calHitsCellId,
			     VInt const& clusterIdToCellId,
			     double totEngy,
			     double * output );

  int	getNeighborId( int	cellId,
		       int	neighborIndex );

  template <class T> double posWeight(T const& calHit, GlobalMethodsClass::WeightingMethod_t method) const;

  template <class T>
  double posWeightTrueCluster(T const& calHit, double cellEngy, GlobalMethodsClass::WeightingMethod_t method) const;

  template <class T> double posWeight(T const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method) const;

  template <class T>
  double posWeight(T const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method,
                   double logWeightConstNow) const;

  double	distance2DPolar( double * pos1,
				 double * pos2 );

  LCCluster calculateEngyPosCM( VInt const& cellIdV,
				MapIntCalHit const& calHitsCellId,
				GlobalMethodsClass::WeightingMethod_t method );

  void	calculateEngyPosCM_EngyV( VInt const& cellIdV,
				  VDouble const& cellEngyV,
				  MapIntCalHit const& calHitsCellId,
				  MapIntLCCluster & clusterCM,
				  int clusterId,
				  GlobalMethodsClass::WeightingMethod_t method );

  template <class T> void updateEngyPosCM(T const& calHit, LCCluster& clusterCM);

  int	checkClusterMergeCM( int clusterId1,
			     int clusterId2,
			     MapIntVInt const& clusterIdToCellId,
			     MapIntCalHit const& calHitsCellId,
			     double				distanceAroundCM,
			     double				percentOfEngyAroungCM,
			     GlobalMethodsClass::WeightingMethod_t method );

  double	getDistanceAroundCMWithEnergyPercent( LCCluster const& clusterCM,
						      VInt const& clusterIdToCellId,
						      MapIntCalHit const& calHitsCellId,
						      double engyPercentage );

  double	getMoliereRadius( MapIntCalHit	const& calHitsCellId,
				  VInt const& clusterIdToCellId,
				  LCCluster const& clusterCM );

  double	getEngyInMoliereFraction( MapIntCalHit const& calHitsCellId,
					  VInt const& clusterIdToCellId,
					  LCCluster const&	clusterCM,
					  double  moliereFraction );

  double	getEngyInMoliereFraction( MapIntCalHit	const& calHitsCellId,
					  VInt const& clusterIdToCellId,
					  LCCluster const& clusterCM,
					  double moliereFraction,
					  MapIntInt & flag );

  void dumpClusters( MapIntLCCluster const& clusterCM );

  std::string printClusters(const int armNow, MapIntMapIntLCCluster const& superClusterCM) const;
  std::string printClusters(MapIntLCCluster const& superClusterCM) const;
};

#endif // LumiCalClusterer_h
