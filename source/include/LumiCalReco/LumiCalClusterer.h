#ifndef LumiCalClusterer_h
#define LumiCalClusterer_h 1

// clustering options
#define	_CLUSTER_MIDDLE_RANGE_ENGY_HITS 1
#define	_MOLIERE_RADIUS_CORRECTIONS 1
#define	_CLUSTER_MIXING_ENERGY_CORRECTIONS 1


// verbosity
#define _GENERAL_CLUSTERER_DEBUG 1
#define _CLUSTER_BUILD_DEBUG 1
#define _VIRTUALCLUSTER_BUILD_DEBUG 0
#define _MOL_RAD_CORRECT_DEBUG 0


#include "GlobalMethodsClass.h"
#include "LCCluster.hh"
#include "VirtualCluster.hh"

#include <IMPL/SimCalorimeterHitImpl.h>
#include <UTIL/CellIDDecoder.h>


#include <string>
#include <map>
#include <vector>

namespace EVENT {
  class LCEvent;
}

namespace IMPL{
  class CalorimeterHitImpl;
}

class LumiCalClustererClass {


public:

  // Constructor
  LumiCalClustererClass( std::string const& lumiNameNow ) ;

  // initialization routine - Called at the begining of the job.
  void init( GlobalMethodsClass::ParametersInt	  const& GlobalParamI,
	     GlobalMethodsClass::ParametersDouble const& GlobalParamD );


  // main actions in each event -Called for every event - the working horse.
  void processEvent( EVENT::LCEvent * evt ) ;

  std::map < int , std::map < int , std::vector<int> > >	_superClusterIdToCellId;
  std::map < int , std::map < int , std::vector<double> > >	_superClusterIdToCellEngy;
  std::map < int , std::map < int , LCCluster > > _superClusterIdClusterInfo;

  void setLumiCollectionName(std::string const& lumiNameNow) { _lumiName = lumiNameNow; }

  /* --------------------------------------------------------------------------
     calculate the distance between two poins in 2D (in cartezian coordinates)
     (first index of arrays is X and the second is Y coordinates {or the other way around})
     -------------------------------------------------------------------------- */
  template <class T, class U>
  static double distance2D(const T *vec1, const U *vec2) {
    const double diff0 = vec1[0]-vec2[0];
    const double diff1 = vec1[1]-vec2[1];
    const double distance = sqrt( diff0*diff0 + diff1*diff1 );
    //assert (distance >= 0);
    return distance;
  }



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
  double	_zFirstLayer, _zLayerThickness, _rMin, _rMax, _rCellLength, _phiCellLength;
  double	_elementsPercentInShowerPeakLayer;
  double	_logWeightConst;
  int	_nNearNeighbor ;
  int	_cellRMax, _cellPhiMax ;
  double	_middleEnergyHitBoundFrac;
  GlobalMethodsClass::WeightingMethod_t _methodCM;
  double	_moliereRadius;
  double	_thetaContainmentBouds[2];
  double	_minSeparationDistance, _minClusterEngyGeV;

  std::map < int , double >	_totEngyArm;
  std::vector < int >		_armsToCluster;

  CellIDDecoder<SimCalorimeterHit> * _mydecoder;

  // methods:
  void	getCalHits( EVENT::LCEvent * evt,
		    std::map < int , std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > > & calHits );


  int	buildClusters(	std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > const& calHits,
			std::map < int , IMPL::CalorimeterHitImpl* > & calHitsCellIdGlobal,
			std::map < int , std::vector<int> > & superClusterIdToCellId,
			std::map < int , std::vector<double> > & superClusterIdToCellEngy,
			std::map < int , LCCluster > & superClusterCM );

  int	initialClusterBuild( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
			     std::map < int , int >			  & cellIdToClusterId,
			     std::map < int , std::vector<int> >	  & clusterIdToCellId,
			     std::map < int , LCCluster > & clusterCM,
			     std::vector < int > const& controlVar );

  int	initialLowEngyClusterBuild( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsSmallEngyCellId,
				    std::map < int , IMPL::CalorimeterHitImpl* > & calHitsCellId,
				    std::map < int , int >			 & cellIdToClusterId,
				    std::map < int , std::vector<int> >		 & clusterIdToCellId,
				    std::map < int , LCCluster >	 & clusterCM );


  int	virtualCMClusterBuild( std::map < int , IMPL::CalorimeterHitImpl* > const&	  calHitsCellId,
			       std::map < int , int >	&			 cellIdToClusterId,
			       std::map < int , std::vector<int> > & clusterIdToCellId,
			       std::map < int , LCCluster >		& clusterCM,
			       std::map < int , VirtualCluster > const& virtualClusterCM );

  int	virtualCMPeakLayersFix(	std::map < int , IMPL::CalorimeterHitImpl* > const&	calHitsCellId,
				std::map < int , int >				& cellIdToClusterId,
				std::map < int , std::vector<int> >		& clusterIdToCellId,
				std::map < int , LCCluster >		& clusterCM,
				std::map < int , VirtualCluster > virtualClusterCM );

  int	buildSuperClusters ( std::map <int , IMPL::CalorimeterHitImpl* > & calHitsCellIdGlobal,
			     std::vector < std::map < int , IMPL::CalorimeterHitImpl* > > const&	calHitsCellId,
			     std::vector < std::map < int , std::vector<int> > >	const&	clusterIdToCellId,
			     std::vector < std::map < int , LCCluster > >	const&	clusterCM,
			     std::vector < std::map < int , VirtualCluster > > const& virtualClusterCM,
			     std::map < int , int > & cellIdToSuperClusterId,
			     std::map < int , std::vector<int> > & superClusterIdToCellId,
			     std::map < int , LCCluster > & superClusterCM );

  int	engyInMoliereCorrections ( std::map <int , IMPL::CalorimeterHitImpl* > const& calHitsCellIdGlobal,
				   std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > const& calHits,
				   std::vector < std::map < int , IMPL::CalorimeterHitImpl* > > const& calHitsCellIdLayer,
				   std::vector < std::map < int , std::vector<int> > > & clusterIdToCellId,
				   std::vector < std::map < int , LCCluster > >        & clusterCM,
				   std::vector < std::map < int , int > >	       & cellIdToClusterId,
				   std::map < int , int >			       & cellIdToSuperClusterId,
				   std::map < int , std::vector<int> >			       & superClusterIdToCellId,
				   std::map < int , LCCluster >			& superClusterCM,
				   double					middleEnergyHitBound,
				   int							detectorArm );



  void	energyCorrections (	std::map < int , std::vector<int> >		& superClusterIdToCellId,
				std::map < int , std::vector<double> >			& superClusterIdToCellEngy,
				std::map < int , LCCluster >		& superClusterCM,
				std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellIdGlobal ) ;


  void	clusterMerger (	      std::map < int , std::vector<double> >	& clusterIdToCellEngy,
			      std::map < int , std::vector<int> >	& clusterIdToCellId,
			      std::map < int , LCCluster > & clusterCM,
			      std::map < int , IMPL::CalorimeterHitImpl* >	calHitsCellIdGlobal ) ;


  void	fiducialVolumeCuts (	std::map < int , std::vector<int> >		& superClusterIdToCellId,
				std::map < int , std::vector<double> >			& superClusterIdToCellEngy,
				std::map < int , LCCluster > & superClusterCM ) ;


  void	getThetaPhiZCluster( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
			     std::vector <int> const& clusterIdToCellId,
			     double totEngy,
			     double * output );

  int	getNeighborId( int	cellId,
		       int	neighborIndex );

  double	posWeight( IMPL::CalorimeterHitImpl const* calHit ,
			   GlobalMethodsClass::WeightingMethod_t method );

  double	posWeightTureCluster( IMPL::CalorimeterHitImpl	* calHit,
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

  LCCluster getEngyPosCMValues( std::vector <int> const& cellIdV,
				std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
				GlobalMethodsClass::WeightingMethod_t method );

  void	calculateEngyPosCM( std::vector <int> const& cellIdV,
			    std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
			    //std::map < int , LCCluster > & clusterCM, int clusterId,
			    LCCluster & cluster,
			    GlobalMethodsClass::WeightingMethod_t method );

  void	calculateEngyPosCM_EngyV(      std::vector <int> const& cellIdV,
				       std::vector <double> const& cellEngyV,
				       std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
				       std::map < int , LCCluster > & clusterCM,
				       int clusterId,
				       GlobalMethodsClass::WeightingMethod_t method );

  void	updateEngyPosCM( IMPL::CalorimeterHitImpl	* calHit,
			 LCCluster & clusterCM );

  int	checkClusterMergeCM( int clusterId1,
			     int clusterId2,
			     std::map < int , std::vector<int> > const& clusterIdToCellId,
			     std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
			     double				distanceAroundCM,
			     double				percentOfEngyAroungCM,
			     GlobalMethodsClass::WeightingMethod_t method );

  double	getDistanceAroundCMWithEnergyPercent( LCCluster const& clusterCM,
						      std::vector < int > const& clusterIdToCellId,
						      std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
						      double engyPercentage );

  double	getMoliereRadius( std::map < int , IMPL::CalorimeterHitImpl* >	const& calHitsCellId,
				  std::vector <int> const& clusterIdToCellId,
				  LCCluster const& clusterCM );

  double	getEngyInMoliereFraction( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
					  std::vector <int> const& clusterIdToCellId,
					  LCCluster const&	clusterCM,
					  double				moliereFraction );

  double	getEngyInMoliereFraction( std::map < int , IMPL::CalorimeterHitImpl* >	const& calHitsCellId,
					  std::vector < int > const& clusterIdToCellId,
					  LCCluster const& clusterCM,
					  double moliereFraction,
					  std::map < int , int > & flag );

};

#endif // LumiCalClusterer_h
