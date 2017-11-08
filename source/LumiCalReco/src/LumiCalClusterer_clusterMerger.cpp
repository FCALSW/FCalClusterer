// Local
#include "SuperTrueClusterWeights.hh"
#include "LumiCalClusterer.h"
// Stdlib
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
//Forward Declaration
namespace IMPL{
  class CalorimeterHitImpl;
}

void LumiCalClustererClass::clusterMerger(MapIntVDouble& clusterIdToCellEngy,
                                          MapIntVInt& clusterIdToCellId,
                                          MapIntLCCluster& clusterCM, MapIntCalHit& calHitsCellIdGlobal) {
  int clusterId, clusterId1, clusterId2;
  int numClusters, numClusters2, numElementsInCluster, cellId;
  double	engyNow;

  std::map < int , std::vector<int> > :: iterator	clusterIdToCellIdIterator;

  std::vector < int >	allClusterIds, cellIdV;
  std::vector < double >	cellEngyV;

  std::map < int , double >			idToCellEngy;
  std::map < int , double > :: iterator		idToCellEngyIterator;

  std::vector < SuperTrueClusterWeights * >	clusterPairWeightsV;





  while(1) {
    /* --------------------------------------------------------------------------
       make a std::vector of all of the cluster Ids, and than write the pair-wise
       properties for all pair combinations
       -------------------------------------------------------------------------- */
    allClusterIds.clear();
    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();

    for(int clusterNow=0; clusterNow<numClusters; clusterNow++, clusterIdToCellIdIterator++){

      clusterId   = (int)(*clusterIdToCellIdIterator).first;  // Id of cluster
      allClusterIds.push_back(clusterId);

    }


    numClusters2 =  allClusterIds.size();
    for(int clusterNow1=0; clusterNow1<numClusters2; clusterNow1++) {
      clusterId1 = allClusterIds[clusterNow1];
      //(BP) avoid creation of "reverse order" pairs
      //     for(int clusterNow2=0; clusterNow2<numClusters2; clusterNow2++) {
      for(int clusterNow2=clusterNow1+1; clusterNow2<numClusters2; clusterNow2++) {
	clusterId2 = allClusterIds[clusterNow2];

	if(clusterId1 == clusterId2) continue;

	SuperTrueClusterWeights * clusterPairWeightsNow = new SuperTrueClusterWeights( clusterId1,
										       clusterId2,
										       clusterCM[clusterId1],
										       clusterCM[clusterId2] );
	/* BP. in the following only distance seems to be important
	 * so make it simpler    
	clusterPairWeightsNow->setWeight( "minEngyDistance",
					  _minSeparationDistance,
					  _minClusterEngySignal );

	if(clusterPairWeightsNow->weight < 0) {         // at least one condition is not true
	  clusterPairWeightsNow->setWeight("distance");
	  if( clusterPairWeightsNow->weight <= minSeparationDistance ){
	    clusterPairWeightsV.push_back(clusterPairWeightsNow);
	  }else{
	    delete clusterPairWeightsNow;
	  }
	} else {                                        // both conditions are false
	  delete clusterPairWeightsNow;
	}
	*/
	clusterPairWeightsNow->setWeight( "distance" );
	if( clusterPairWeightsNow->weight < _minSeparationDistance ) clusterPairWeightsV.push_back(clusterPairWeightsNow);
	else                                                        delete clusterPairWeightsNow;
      }
    }

    // if all the pairs made the cut than finish the loop
    if( clusterPairWeightsV.empty() ) break;

    // choose the pair with the shortest weight (distance)
    sort(clusterPairWeightsV.begin(), clusterPairWeightsV.end(), SuperTrueClusterWeights::Compare);


    clusterId1 = clusterPairWeightsV[0]->superClusterId;
    clusterId2 = clusterPairWeightsV[0]->trueClusterId;

    /* --------------------------------------------------------------------------
       go over all hits in the discarded cluster write a new hit-energy std::map
       -------------------------------------------------------------------------- */
    for(int clusterNow = 0; clusterNow < 2; clusterNow++) {
      if(clusterNow == 0) clusterId = clusterId1;
      if(clusterNow == 1) clusterId = clusterId2;

      numElementsInCluster = clusterIdToCellId[clusterId].size();
      for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
	cellId = clusterIdToCellId[clusterId][hitNow];

	engyNow = clusterIdToCellEngy[clusterId][hitNow];
	idToCellEngy[cellId] += engyNow;
      }
    }

    //cleanUp both clusters
    clusterIdToCellId.erase(clusterId1);  clusterCM.erase(clusterId1);  clusterIdToCellEngy.erase(clusterId1);
    clusterIdToCellId.erase(clusterId2);  clusterCM.erase(clusterId2);  clusterIdToCellEngy.erase(clusterId2);

    /* --------------------------------------------------------------------------
       write the new merged cluster with the Id of one of the clusters
       -------------------------------------------------------------------------- */
    clusterId = clusterId1;
    idToCellEngyIterator = idToCellEngy.begin();
    numClusters          = idToCellEngy.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, idToCellEngyIterator++){
      cellId = (int)(*idToCellEngyIterator).first;

      engyNow = idToCellEngy[cellId];

      clusterIdToCellId[clusterId].push_back(cellId);
      clusterIdToCellEngy[clusterId].push_back(engyNow);
    }


    /* --------------------------------------------------------------------------
       compute the total energy and center of mass of the new merged cluster
       -------------------------------------------------------------------------- */
    cellIdV   = clusterIdToCellId[clusterId];
    cellEngyV = clusterIdToCellEngy[clusterId];

    // initialize the energy/position std::vector for new clusters only
    clusterCM[clusterId] = LCCluster();

    // calculate/update the energy/position of the CM
    calculateEngyPosCM_EngyV( cellIdV, cellEngyV, calHitsCellIdGlobal,
			      clusterCM, clusterId, _methodCM);

    // cleanUp
    cellIdV.clear();  cellEngyV.clear();

    /* --------------------------------------------------------------------------
       cleanUp virtual memory allocations
       -------------------------------------------------------------------------- */
    int numClustersToMerge = clusterPairWeightsV.size();
    for(int clusterNow1=0; clusterNow1<numClustersToMerge; clusterNow1++)
      delete clusterPairWeightsV[clusterNow1];

    clusterPairWeightsV.clear();
    idToCellEngy.clear();
  }




  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  if ( not clusterIdToCellId.empty() ) {
    streamlog_out( MESSAGE3 ) << "Clusters:"  << std::endl;
  } else {
    streamlog_out( MESSAGE3 ) << "No Clusters on this side" << std::endl;
  }

  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++) {
    clusterId = (int)(*clusterIdToCellIdIterator).first;

    streamlog_out( MESSAGE3 ) << "Id " << std::setw(4) << clusterId
			      << clusterCM[clusterId] << std::endl;

  }

#endif



  // cleanUp
  allClusterIds.clear();

}
