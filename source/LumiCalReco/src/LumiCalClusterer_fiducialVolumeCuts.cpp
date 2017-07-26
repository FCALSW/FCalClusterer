//Local
#include "Global.hh"
#include "LumiCalClusterer.h"
// Stdlib
#include <map>
#include <vector>
#include <cmath>
#include <iostream>

void LumiCalClustererClass:: fiducialVolumeCuts( std::map < int , std::vector<int> > & superClusterIdToCellId,
						 std::map < int , std::vector<double> > & superClusterIdToCellEngy,
						 std::map < int , LCCluster > & superClusterCM ) {
  if( not _cutOnFiducialVolume ) {
    streamlog_out( DEBUG ) << "Skipping fiducial volume cuts" << std::endl;
    return;
  }

  int	numSuperClusters, superClusterId;
  double	thetaSuperCluster;

  std::map < int , std::vector<int> > ::iterator superClusterIdToCellIdIterator;

  std::vector < int >	clusterIdToErase;

  /* --------------------------------------------------------------------------
     discard true/reconstructed clusters that are outside the fiducial volume
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters               = superClusterIdToCellId.size();
  for(int superClusterNow=0; superClusterNow<numSuperClusters; superClusterNow++, superClusterIdToCellIdIterator++){
    superClusterId   = (int)(*superClusterIdToCellIdIterator).first;  // Id of cluster

    thetaSuperCluster = fabs(superClusterCM[superClusterId].getTheta());
    if(thetaSuperCluster < _thetaContainmentBounds[0] || thetaSuperCluster > _thetaContainmentBounds[1])
      clusterIdToErase.push_back(superClusterId);

  }

  numSuperClusters = clusterIdToErase.size();
  for(int superClusterNow=0; superClusterNow<numSuperClusters; superClusterNow++) {
    superClusterId = clusterIdToErase[superClusterNow];

#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << "\tErase cluster " << superClusterId << std::endl;
#endif

    superClusterIdToCellId.erase(superClusterId);
    superClusterIdToCellEngy.erase(superClusterId);
    superClusterCM.erase(superClusterId);
  }
  clusterIdToErase.clear();

  return ;
}
