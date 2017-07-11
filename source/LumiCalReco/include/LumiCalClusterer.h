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
#include "VirtualCluster.hh"
#include "ProjectionInfo.hh"

#include <IMPL/SimCalorimeterHitImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <string>
#include <map>
#include <memory>
#include <vector>

namespace EVENT {
  class LCEvent;
}

namespace IMPL{
  class CalorimeterHitImpl;
}


class LumiCalClustererClass {

  typedef std::vector <IMPL::CalorimeterHitImpl*> VecCalHit;
  typedef std::vector < double >                  VDouble;
  typedef std::vector < int >                     VInt;

  typedef std::map < int , IMPL::CalorimeterHitImpl* >  MapIntCalHit;
  typedef std::map < int , LCCluster >                  MapIntLCCluster;

  typedef std::map < int , VecCalHit >                  MapIntVCalHit;
  typedef std::map < int , VDouble >                    MapIntVDouble;
  typedef std::map < int , VInt >                       MapIntVInt;

  typedef std::map < int , MapIntLCCluster >            MapIntMapIntLCCluster;
  typedef std::map < int , MapIntVCalHit >              MapIntMapIntVCalHit;
  typedef std::map < int , MapIntVDouble >              MapIntMapIntVDouble; 
  typedef std::map < int , MapIntVInt >                 MapIntMapIntVInt; 
  typedef std::map < int , ProjectionInfo >             MapIntProjectionInfo;
  typedef std::map < int , VirtualCluster >             MapIntVirtualCluster;
  typedef std::map < int , double >                     MapIntDouble;
  typedef std::map < int , int >                        MapIntInt;

  typedef std::map < int , MapIntCalHit > MapIntMapIntCalHit;


  typedef std::vector < MapIntCalHit >         VMapIntCalHit;
  typedef std::vector < MapIntInt >            VMapIntInt;
  typedef std::vector < MapIntLCCluster >      VMapIntLCCluster;
  typedef std::vector < MapIntVInt >           VMapIntVInt;
  typedef std::vector < MapIntVirtualCluster > VMapIntVirtualCluster;

  typedef std::vector < VDouble >              VVDouble;

public:

  // Constructor
  LumiCalClustererClass( std::string const& lumiNameNow ) ;

  // initialization routine - Called at the begining of the job.
  void init( GlobalMethodsClass const& gmc );


  // main actions in each event -Called for every event - the working horse.
  int processEvent( EVENT::LCEvent * evt ) ;

  MapIntMapIntVInt       _superClusterIdToCellId;
  MapIntMapIntVDouble    _superClusterIdToCellEngy;
  MapIntMapIntLCCluster  _superClusterIdClusterInfo;

  void setLumiCollectionName(std::string const& lumiNameNow) { _lumiName = lumiNameNow; }

protected:

  LumiCalClustererClass(LumiCalClustererClass const& rhs);
  LumiCalClustererClass& operator=(LumiCalClustererClass const& rhs);

  // Processor Parameters
  std::string	_lumiName;
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
  double	_minSeparationDistance, _minClusterEngyGeV, _minClusterEngySignal;

  MapIntDouble _totEngyArm;
  MapIntInt    _numHitsInArm;
  //  VInt _armsToCluster;

  std::unique_ptr<CellIDDecoder<SimCalorimeterHit>> _mydecoder{};

  GlobalMethodsClass _gmc;
  bool _useDD4hep;

  // global to local rotations mtx elements
  std::map < int , std::map< std::string, double> > RotMat;

  // methods:
  int	getCalHits( EVENT::LCEvent * evt,
		    MapIntMapIntVCalHit & calHits );


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


  void	clusterMerger (	      MapIntVDouble & clusterIdToCellEngy,
			      MapIntVInt & clusterIdToCellId,
			      MapIntLCCluster & clusterCM,
			      MapIntCalHit calHitsCellIdGlobal ) ;


  void	fiducialVolumeCuts (	MapIntVInt & superClusterIdToCellId,
				MapIntVDouble & superClusterIdToCellEngy,
				MapIntLCCluster & superClusterCM ) ;


  void	getThetaPhiZCluster( MapIntCalHit const& calHitsCellId,
			     VInt const& clusterIdToCellId,
			     double totEngy,
			     double * output );

  int	getNeighborId( int	cellId,
		       int	neighborIndex );

  double	posWeight( IMPL::CalorimeterHitImpl const* calHit ,
			   GlobalMethodsClass::WeightingMethod_t method );

  double	posWeightTrueCluster( IMPL::CalorimeterHitImpl const* calHit,
				      double			cellEngy,
				      GlobalMethodsClass::WeightingMethod_t method );

  double	posWeight( IMPL::CalorimeterHitImpl const* calHit,
			   double		totEngy,
			   GlobalMethodsClass::WeightingMethod_t method );

  double	posWeight( IMPL::CalorimeterHitImpl const* calHit,
			   double		totEngy,
			   GlobalMethodsClass::WeightingMethod_t method,
			   double		logWeightConstNow );


  double	distance2DPolar( double * pos1,
				 double * pos2 );

  double	thetaPhiCell( int	cellId,
			      GlobalMethodsClass::Coordinate_t	output );

  LCCluster calculateEngyPosCM( VInt const& cellIdV,
				MapIntCalHit const& calHitsCellId,
				GlobalMethodsClass::WeightingMethod_t method );

  void	calculateEngyPosCM_EngyV( VInt const& cellIdV,
				  VDouble const& cellEngyV,
				  MapIntCalHit const& calHitsCellId,
				  MapIntLCCluster & clusterCM,
				  int clusterId,
				  GlobalMethodsClass::WeightingMethod_t method );

  void	updateEngyPosCM( IMPL::CalorimeterHitImpl	* calHit,
			 LCCluster & clusterCM );

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

};

#endif // LumiCalClusterer_h
