// Local
#include "LumiCalClusterer.h"
#include "LCCluster.hh"
#include "SortingFunctions.hh"
#include "Global.hh"
//LCIO
#include <IMPL/CalorimeterHitImpl.h>
// Stdlib
#include <map>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <stdexcept>



/* =========================================================================
   LumiCalClustererClass :: initialClusterBuild
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */

int LumiCalClustererClass::initialClusterBuild(	std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
						std::map < int , int > & cellIdToClusterId,
						std::map < int , std::vector<int> > & clusterIdToCellId,
						std::map < int , LCCluster > & clusterCM,
						std::vector < int > const& controlVar  ) {

  // general variables
  int	numClusters, numElementsInCluster;
  int	cellIdHit, cellIdNeighbor, nNeighborsConectedToMe;

  const int numElementsInLayer = calHitsCellId.size();

  /* --------------------------------------------------------------------------
     layer parameters
     -------------------------------------------------------------------------- */
  // (*). Small clusters are merged with large clusters if they are close. A weight is computed for each
  //	close cluster candidate, and the max-weighted cluster is chosen. This is done twice. First small
  //	clusters try to be merged into large ones, and then large clusters try to merge small ones into themselvs.
  // ----------------------------------------------------------------------------------------------------------------
  //	(1). definitions of 'small' and 'large' sizes
  int	numElementsSmallClusterToMerge = int(.1 * numElementsInLayer); //5;  // max number of elements in small cluster
  int	numElementsLargeClusterToMerge = int(.15 * numElementsInLayer); //10; // min number of elements in a large cluster

  //	(2). decide to merge if the two clusters are close to each other, and if the merged cluster's CM
  //	isnt shifted to an area where there are few cal hits (small percentage of CM energy)
  double	mergeScanDistanceFromLargeClusterCM = _moliereRadius ;  // max distance of small cluster from CM of large one
  double	distanceToCollectEngyAroundCM = .5*_moliereRadius ;  // distance from the CM to collect energy
  double	percentOfEngyAroungCM = .4;  // percenatge of cluster CM energy

  //	(3). the weight for merging two clusters is a Power() function of the energy of the large cluster and the
  //	distance to its CM. two different weight functions for the two pases.
  int	smallClusterEngyCMPower = 1;
  int	smallClusterDistanceCMPower = -3;

  int	largeClusterEngyCMPower = 1;
  int	largeClusterDistanceCMPower = -1;


  // (*). Remaining small clusters are force-merged into large clusters. The choise of who to merge with is done through
  //	weights.
  // -------------------------------------------------------------------------------------------------------------------
  //	(1). definitions of 'small' and 'large' sizes
  int	numElementsSmallClusterToMergeForced = int(.08 * numElementsInLayer); //5;  // max number of elements in small cluster

  //	(2). the weight for merging two clusters is a Power() function of the energy of the
  //		large cluster and the distance to its CM.
  int	smallClusterEngyCMPowerForced = 0;
  int	smallClusterDistanceCMPowerForced = -1;


  // control options for enabeling the different rutines (1 is for active)
  int	mergeOneHitClusters		= controlVar.at(0);
  int	mergeSmallToLargeClusters	= controlVar.at(1);
  int	mergeLargeToSmallClusters	= controlVar.at(2);
  int	forceMergeSmallToLargeClusters	= controlVar.at(3);


  /* --------------------------------------------------------------------------
     maps for................
     -------------------------------------------------------------------------- */
  // map param: (1). cellId of cal hit , (2). Id of cluster the cal hit belongs to
  // std::map < int , int >		cellIdToClusterId;

  // map param: (1). Id of cluster the cal hit belongs to , (2). std::vector of cellIds of cal hit
  // std::map < int , std::vector<int> >		clusterIdToCellId;
  std::map < int , std::vector<int> > :: iterator	clusterIdToCellIdIterator  ;

  // map param: (1). Id of cluster , (2). std::vector of: (A). total energy , (B). CM coordinats in the order: x,y
  //std::map < int , LCCluster > clusterCM;

  // std::map param: (1). cellId of cal hit , (2). ponter to call hit

  // std::vector for holding the Ids of cells/clusters
  //std::vector <int>	cellIdV;
  std::vector <int>	clusterIdV;

  // // temporary std::vectors for sorting purposes
  // std::vector < std::vector <double> >	clusterIdEngyV2;
  // std::vector <double>			clusterIdEngyV;

  // counter for cluster Id
  int	clusterId = 0 ;

  /* --------------------------------------------------------------------------
     maps for storing for each cal hit the highest-energy near neighbor, and for
     each cal hit to keep track of the other cal hits for which it is the highest
     energy nearest neighbor.
     -------------------------------------------------------------------------- */
  // map param: (1). cellId of cal hit , (2). cellId of highest energy near neighbor
  std::map <int , int >	isConnectedToNeighbor ;
  // map param: (1). cellId of cal hit , (2). cellIds of neighbors which are connected to the cal hit
  std::map <int , std::vector <int> >	neighborsConectedToMe ;

  // copy hits in this layer to a cal hit std::vector
  std::vector <IMPL::CalorimeterHitImpl*>	calHitsLayer ;

  std::map <int , IMPL::CalorimeterHitImpl* > :: const_iterator calHitsCellIdIterator= calHitsCellId.begin(),
    calHitsEnd = calHitsCellId.end();
  for(; calHitsCellIdIterator != calHitsEnd; ++calHitsCellIdIterator){
    cellIdHit = (int)(*calHitsCellIdIterator).first;
    try {
      calHitsLayer.push_back( calHitsCellIdIterator->second );
      // initialization
      isConnectedToNeighbor[cellIdHit] = 0;
      cellIdToClusterId[cellIdHit] = 0;
    } catch (std::out_of_range &e) {
      streamlog_out( ERROR ) << "Out Of Range " << cellIdHit
			     << "in calHitsCellID" << std::endl;
    }
  }

  // sort acording to energy in ascending order (lowest energy is first)
  sort( calHitsLayer.begin(), calHitsLayer.end(), HitEnergyCmpAsc );

  /* --------------------------------------------------------------------------
     connect each cal hit to it's highest-energy near neighbor.
     -------------------------------------------------------------------------- */
  for(int j=0; j<(int)calHitsLayer.size(); j++) {

    cellIdHit  = (int)calHitsLayer[j]->getCellID0();
    const double engyCalHit   = calHitsLayer[j]->getEnergy();

    // go on to next cal hit if this hit has already been registered
    if(isConnectedToNeighbor[cellIdHit] != 0) continue;

    double maxEngyNeighbor = 0.;

    // get the Ids of the neighbor
    for(int neighborIndex=0; neighborIndex < _nNearNeighbor ; neighborIndex++) {
      // find cellId of neighbor
      cellIdNeighbor = getNeighborId(cellIdHit, neighborIndex);
      if(cellIdNeighbor == 0) continue;
      try {
	// if the neighbor has a cal hit...
	std::map <int , IMPL::CalorimeterHitImpl* > :: const_iterator neighborIt =
	  calHitsCellId.find(cellIdNeighbor);
	if( neighborIt != calHitsCellId.end() ) {
	  IMPL::CalorimeterHitImpl* neighbor = neighborIt->second;

	  //if(tmpFlag==1) cout << "neighbor " << cellIdNeighbor
	  //	<< " " << neighborIndex <<endl;
	  double engyNeighbor = neighbor->getEnergy();

	  //if(tmpFlag==1) cout << "\tmax/me/neighbor \t"<<maxEngyNeighbor <<"\t"
	  //	<< engyCalHit<< "\t"<<engyNeighbor <<endl;
	  if((maxEngyNeighbor < engyNeighbor) && (engyNeighbor >= engyCalHit)) {
	    // register the neighbor at the cal hit
	    isConnectedToNeighbor[cellIdHit] = cellIdNeighbor;
	    // update highest-energy counter
	    maxEngyNeighbor = engyNeighbor;
	  }
	}//if
      } catch (std::out_of_range & e) {
	streamlog_out( ERROR ) << "Out of range " << cellIdNeighbor
			       << " in calHitsCellId:" <<  std::endl;
      }

    }//For all neighbor

    if(maxEngyNeighbor > 0) {
      // check if the neighbor has already been registered
      cellIdNeighbor = isConnectedToNeighbor[cellIdHit];
      // modify the conditional control variable

      //APS: Does this make sense: When the possible
      //neighbor already has neighbors we don't quit
      //here, actually we will always quit, because of the else statement
      // if(isConnectedToNeighbor[cellIdNeighbor] != 0) neighborFound = true;

      // register the cal hit at the neighbor
      neighborsConectedToMe[cellIdNeighbor].push_back(cellIdHit);
      // in the next iteration work with the cal hit's highest-energy near neighbor
      cellIdHit = isConnectedToNeighbor[cellIdHit];
    }

  }//for all hits in layer


  /* --------------------------------------------------------------------------
     create clusters from each connected bunch of cal hits.
     -------------------------------------------------------------------------- */
  // sort according to energy in descending order (highest energy is first)
  sort( calHitsLayer.begin(), calHitsLayer.end(), HitEnergyCmpDesc );

  for(int j=0; j<(int)calHitsLayer.size(); j++) {

    std::vector <int> neighborFoundId ;

    cellIdHit  = (int)calHitsLayer[j]->getCellID0();

    // if the cal hit has already been registered in a cluster continue to the next one
    if(cellIdToClusterId[cellIdHit] > 0) continue;

    // add the hit to a new cluster
    clusterId++;
    cellIdToClusterId[cellIdHit] = clusterId;
    clusterIdToCellId[clusterId].push_back(cellIdHit);

    bool neighborFound = true ;
    while(neighborFound == true) {
      // get the Ids of the connected neighbor
      nNeighborsConectedToMe = (int)neighborsConectedToMe[cellIdHit].size();
      for(int k=0; k<nNeighborsConectedToMe; k++) {

	// register the neighbor in the cluster
	cellIdNeighbor = neighborsConectedToMe[cellIdHit][k];
	cellIdToClusterId[cellIdNeighbor] = clusterId ;
	clusterIdToCellId[clusterId].push_back(cellIdNeighbor);

	// store the connected cal hit Id - it will be registered in a cluster next...
	neighborFoundId.push_back(cellIdNeighbor);
      }

      // in case of a loop in connections, make sure that if a cell has already been
      // considered, then it wont be reanalyzed again...
      neighborsConectedToMe[cellIdHit].clear();

      // every neighbor found is taken off from the neighborFoundId std::vector
      // once this std::vector is empty, end the while(neighborFound == 1) loop

      if(neighborFoundId.size() == 0) {
	neighborFound = false;
      } else {
	cellIdHit = neighborFoundId.back();
	neighborFoundId.pop_back();
      }
    }
  }
  // cleanUp
  isConnectedToNeighbor.clear(); neighborsConectedToMe.clear();

  /* --------------------------------------------------------------------------
     merge clusters that have one hit only to the nearest cluster which is
     a near neighbor (choose the neighbor with the highest energy).
     -------------------------------------------------------------------------- */
  if(mergeOneHitClusters == 1) {
    // go over all cellIds which have been registered in a cluster (first parameter of cellIdToClusterId).
    std::map < int , int > :: iterator
      cellIdToClusterIdIterator = cellIdToClusterId.begin(),
      cEnd = cellIdToClusterId.begin();
    for(; cellIdToClusterIdIterator != cEnd; ++cellIdToClusterIdIterator){
      cellIdHit = (int)(*cellIdToClusterIdIterator).first;	// cellId of cal hit
      int clusterIdHit = (int)(*cellIdToClusterIdIterator).second;	// cluster Id of cal hit

      if(clusterIdToCellId[clusterIdHit].size() == 1) {
	int	maxEngyNeighborId = 0;
	double	maxEngyNeighbor   = 0.;

	// get the Ids of the neighbor
	for(int neighborIndex=0; neighborIndex < _nNearNeighbor ; neighborIndex++) {
	  // find cellId of neighbor
	  cellIdNeighbor = getNeighborId(cellIdHit, neighborIndex);
	  if(cellIdNeighbor == 0) continue;

	  // if the neighbor has a cal hit...
	  if( calHitsCellId.find(cellIdNeighbor) != calHitsCellId.end() ) {
	    double engyNeighbor = calHitsCellId.at(cellIdNeighbor)->getEnergy();

	    // find the neighbor with the highest energy
	    if(maxEngyNeighbor < engyNeighbor) {
	      // store the cellId of the neighbor
	      maxEngyNeighborId = cellIdNeighbor;
	      // update highest-energy counter
	      maxEngyNeighbor = engyNeighbor;
	    }

	  }
	}

	if(maxEngyNeighbor > 0) {
	  // connect to the neighbor with the highest energy
	  int clusterIdNeighbor = cellIdToClusterId[maxEngyNeighborId];
	  cellIdToClusterId[cellIdHit] = clusterIdNeighbor;
	  clusterIdToCellId[clusterIdNeighbor].push_back(cellIdHit);

	  // cleanUp the discarded cluster
	  clusterIdToCellId.erase(clusterIdHit);
	}
      }
    }
  }


  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each cluster
     -------------------------------------------------------------------------- */
  clusterCM.clear();
  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
    clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
    std::vector<int> const& cellIdV = clusterIdToCellIdIterator->second;	// cal-hit-Ids in cluster

    // initialize the energy/position std::vector for this cluster
    //clusterCM[clusterId].clear();
    // calculate the energy/position of the CM
    calculateEngyPosCM(cellIdV, calHitsCellId, clusterCM[clusterId], _methodCM);

  }
  //printMapVector(__func__,__LINE__,clusterCM);

  /* --------------------------------------------------------------------------
     merge clusters with close CM positions -
     (1). first go over small clusters and merge with large clusters.
     -------------------------------------------------------------------------- */
  if(mergeSmallToLargeClusters == 1){
    // sort the clusterIds according to clusterCM energy in ascending order (lowest energy is first)
    std::vector< std::vector<double> > clusterIdEngyV2;

    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();
    for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
      clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
      std::vector<double> clusterIdEngyV(2, 0.0);
      clusterIdEngyV[0] = clusterCM[clusterId].getE();
      clusterIdEngyV[1] = (double)clusterId;
      clusterIdEngyV2.push_back( clusterIdEngyV );
    }

    sort(clusterIdEngyV2.begin(),clusterIdEngyV2.end(),clusterCMEnergyCmpAsc);

    // copy the Ids that are now in order to a std::vector
    numClusters = clusterIdEngyV2.size();
    for(int j=0; j<numClusters; j++)
      clusterIdV.push_back( (int)clusterIdEngyV2[j][1] );

    // cleanUp
    clusterIdEngyV2.clear();

    // merge close neighboring clusters
    while(clusterIdV.size()>1){

      std::vector <double> weightedDistanceV, engyCloseCluster;
      std::vector <int> closeClusterId, closeClusterIdPositionInVec;

      // start with the highest-energy cluster
      if((int)clusterIdToCellId[clusterIdV[0]].size() > numElementsSmallClusterToMerge) {
	clusterIdV.erase( clusterIdV.begin() );
	continue;
      }

      // compute distance to neighboring clusters
      numClusters = clusterIdV.size();
      for(int j=1; j<numClusters; j++){
	clusterId = clusterIdV[j];

	const double distanceCM = distance2D(clusterCM[clusterIdV[0]].getPosition(),clusterCM[clusterId].getPosition());

	int considerCloseCluster = 0;
	// choose close clusters
	if((int)clusterIdToCellId[clusterId].size() >= numElementsLargeClusterToMerge)
	  if(distanceCM < mergeScanDistanceFromLargeClusterCM)
	    considerCloseCluster = checkClusterMergeCM(  clusterIdV[0],clusterIdV[j],
							 clusterIdToCellId,calHitsCellId,
							 distanceToCollectEngyAroundCM,
							 percentOfEngyAroungCM,
							 _methodCM  );

	// if all the conditions regarding the close cluster were satisfied, add
	// it to the list of possible merging partner, one of which will be chosen next
	if(considerCloseCluster == 1) {
	  // store the weight of the close cluster
	  const double engyCM = clusterCM[clusterId].getE();
	  double closeClusterWeight = pow(engyCM,smallClusterEngyCMPower)
	    * pow(distanceCM,smallClusterDistanceCMPower);
	  weightedDistanceV.push_back(closeClusterWeight);

	  // store the energy, Id and position in the std::vector
	  // clusterIdV of the close cluster
	  engyCloseCluster.push_back(engyCM);
	  closeClusterId.push_back(clusterId);
	  closeClusterIdPositionInVec.push_back(j);
	}
      }

      // find the close cluster with the higest weightedDistance
      int	maxWDClusterId(0);//, positionInVec = 0;
      double	weightedDistance = 0.;

      int numCloseClusters = weightedDistanceV.size();
      if(numCloseClusters>0){
	for(int j=0; j<numCloseClusters; j++){
	  if(weightedDistance < weightedDistanceV[j]) {
	    weightedDistance = weightedDistanceV[j];
	    maxWDClusterId   = closeClusterId[j];
	    //positionInVec    = closeClusterIdPositionInVec[j];
	  }

	}
	// merge the two clusters
	const int clusterId1 = clusterIdV[0];
	const int clusterId2 = maxWDClusterId;

	numElementsInCluster = clusterIdToCellId[clusterId1].size();
	for(int j=0; j<numElementsInCluster; j++){
	  // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
	  cellIdHit = clusterIdToCellId[clusterId1][j];
	  clusterIdToCellId[clusterId2].push_back(cellIdHit);
	  cellIdToClusterId[cellIdHit] = clusterId2;

	  //  update the totalEnergy counter and CM position of the cluster
	  updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
	}

	// cleanUp
	clusterIdToCellId.erase(clusterId1);
	clusterCM.erase(clusterId1);
      }

      // erase the merged cluster Id from the clusterIdV std::vector.
      // if there are no close neighbors than positionInVec still has the initial zero
      // value, and so the initial cluster is erased from clusterIdV.
      clusterIdV.erase( clusterIdV.begin() );
    }
    // cleanUp
    clusterIdV.clear();
  }


  /* --------------------------------------------------------------------------
     merge clusters with close CM positions
     (2). go over all clusters and try to merge.
     -------------------------------------------------------------------------- */
  if(mergeLargeToSmallClusters == 1) {
    // sort the clusterIds according to clusterCM energy in
    // descending order (highest energy is first)

    std::vector< std::vector<double> > clusterIdEngyV2;

    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();
    for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
      clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
      std::vector<double> clusterIdEngyV(2, 0.0);
      clusterIdEngyV[0] = clusterCM[clusterId].getE();
      clusterIdEngyV[1] = (double)clusterId;
      clusterIdEngyV2.push_back( clusterIdEngyV );
    }

    sort(clusterIdEngyV2.begin(),clusterIdEngyV2.end(),clusterCMEnergyCmpDesc);

    // copy the Ids that are now in order to a std::vector
    numClusters = clusterIdEngyV2.size();
    for(int j=0; j<numClusters; j++)
      clusterIdV.push_back( (int)clusterIdEngyV2[j][1] );

    // cleanUp
    clusterIdEngyV2.clear();

    // merge close neighboring clusters
    while(clusterIdV.size()>1){

      double	CM1[3], CM2[3], distanceCM, engyCM;
      int	clusterId1, clusterId2;
      std::vector <double> weightedDistanceV, engyCloseCluster;
      std::vector <int> closeClusterId, closeClusterIdPositionInVec;

      // start with the highest-energy cluster
      clusterId = clusterIdV[0];
      CM1[0] = clusterCM[clusterId].getX();
      CM1[1] = clusterCM[clusterId].getY();

      // compute distance to neighboring clusters
      numClusters = clusterIdV.size();

      for(int j=1; j<numClusters; j++){

	int clusterIdJ = clusterIdV[j];
	CM2[0] = clusterCM[clusterIdJ].getX();
	CM2[1] = clusterCM[clusterIdJ].getY();

	distanceCM = distance2D(CM1,CM2);

	int considerCloseCluster = 0;
	if(distanceCM < mergeScanDistanceFromLargeClusterCM)
	  considerCloseCluster = checkClusterMergeCM(  clusterIdV[0],clusterIdV[j],
						       clusterIdToCellId,calHitsCellId,
						       distanceToCollectEngyAroundCM,
						       percentOfEngyAroungCM,
						       _methodCM  );

	// if all the conditions regarding the close cluster were satisfied, add
	// it to the list of possible merging partner, one of which will be chosen next
	if(considerCloseCluster == 1) {
	  // store the weight of the close cluster
	  engyCM = clusterCM[clusterIdJ].getE();
	  assert ( distanceCM >= 0 && engyCM > 0 );//APS:: assert >= now so 0.0 passes

	  double closeClusterWeight = pow(engyCM,largeClusterEngyCMPower)
	    * pow(distanceCM,largeClusterDistanceCMPower);
	  weightedDistanceV.push_back(closeClusterWeight);

	  // store the energy, Id and position in the std::vector
	  // clusterIdV of the close cluster
	  engyCloseCluster.push_back(engyCM);
	  closeClusterId.push_back(clusterIdJ);
	  closeClusterIdPositionInVec.push_back(j);
	}
      }// for all other clusters

      // find the close cluster with the higest weightedDistance
      int	maxWDClusterId(0), positionInVec = 0;
      double	weightedDistance = 0.;

      int numCloseClusters = weightedDistanceV.size();
      if(numCloseClusters>0){
	for(int j=0; j<numCloseClusters; j++){
	  if(weightedDistance < weightedDistanceV[j]) {
	    weightedDistance = weightedDistanceV[j];
	    maxWDClusterId   = closeClusterId[j];
	    positionInVec    = closeClusterIdPositionInVec[j];
	  }

	}
	// merge the two clusters
	clusterId2 = clusterIdV[0];
	clusterId1 = maxWDClusterId;

	numElementsInCluster = clusterIdToCellId[clusterId1].size();
	for(int j=0; j<numElementsInCluster; j++){
	  // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
	  cellIdHit = clusterIdToCellId[clusterId1][j];
	  clusterIdToCellId[clusterId2].push_back(cellIdHit);
	  cellIdToClusterId[cellIdHit] = clusterId2;

	  //  update the totalEnergy counter and CM position of the cluster
	  updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
	}

	// cleanUp
	clusterIdToCellId.erase(clusterId1);
	clusterCM.erase(clusterId1);
      }

      // erase the merged cluster Id from the clusterIdV std::vector.
      // if there are no close neighbors than positionInVec still has the initial zero
      // value, and so the initial cluster is erased from clusterIdV.
      clusterIdV.erase( clusterIdV.begin() + positionInVec );
    } //while
  }

  /* --------------------------------------------------------------------------
     merge small clusters with large ones
     -------------------------------------------------------------------------- */
  if(forceMergeSmallToLargeClusters == 1) {
    std::vector <int> smallClusterIdV, largeClusterIdV;
    int numSmallClusters, numLargeClusters;

    // sort clusters into two std::vectors according to the number of elements in each cluster
    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();
    for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
      clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
      numElementsInCluster = clusterIdToCellId[clusterId].size();

      // ???????? DECIDE/FIX - improve this number ????????
      if(numElementsInCluster <= numElementsSmallClusterToMergeForced)
	smallClusterIdV.push_back(clusterId);
      else
	largeClusterIdV.push_back(clusterId);

    }

    numSmallClusters = smallClusterIdV.size();
    for(int j=0; j<numSmallClusters; j++) {

      double	CM1[3], CM2[3], distanceCM, engyCM;
      int	clusterId1, clusterId2;
      std::vector <double> weightedDistanceV, engyCloseCluster;
      std::vector <int> closeClusterId;

      // start with the highest-energy cluster
      clusterId = smallClusterIdV[j];
      CM1[0] = clusterCM[clusterId].getX();
      CM1[1] = clusterCM[clusterId].getY();

      // compute distance to neighboring clusters
      numLargeClusters = largeClusterIdV.size();
      for(int k=0; k<numLargeClusters; k++){
	clusterId = largeClusterIdV[k];
	CM2[0] = clusterCM[clusterId].getX();
	CM2[1] = clusterCM[clusterId].getY();

	distanceCM = distance2D(CM1,CM2);

	int considerCloseCluster = 1;
	// if all the conditions regarding the close cluster were satisfied, add
	// it to the list of possible merging partner, one of which will be chosen next
	if(considerCloseCluster == 1) {
	  // store the weight of the close cluster
	  engyCM = clusterCM[clusterId].getE();

	  assert (distanceCM > 0 && engyCM > 0);

	  double closeClusterWeight = pow(engyCM,smallClusterEngyCMPowerForced)
	    * pow(distanceCM,smallClusterDistanceCMPowerForced);
	  weightedDistanceV.push_back(closeClusterWeight);
	  closeClusterId.push_back(clusterId);
	}
      }

      // find the close cluster with the higest weightedDistance
      int	maxWDClusterId(0);
      double	weightedDistance = 0.;

      int numCloseClusters = weightedDistanceV.size();
      if(numCloseClusters>0){
	for(int k=0; k<numCloseClusters; k++){
	  if(weightedDistance < weightedDistanceV[k]) {
	    weightedDistance = weightedDistanceV[k];
	    maxWDClusterId   = closeClusterId[k];
	  }

	}
	// merge the two clusters
	clusterId1 = smallClusterIdV[j];
	clusterId2 = maxWDClusterId;

	numElementsInCluster = clusterIdToCellId[clusterId1].size();
	for(int k=0; k<numElementsInCluster; k++){
	  // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
	  cellIdHit = clusterIdToCellId[clusterId1][k];
	  clusterIdToCellId[clusterId2].push_back(cellIdHit);
	  cellIdToClusterId[cellIdHit] = clusterId2;

	  //  update the totalEnergy counter and CM position of the cluster
	  updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
	}

	// cleanUp
	clusterIdToCellId.erase(clusterId1);
	clusterCM.erase(clusterId1);
      }
    }
  }

  return 1;
}


/* =========================================================================
   LumiCalClustererClass :: initialLowEngyClusterBuild
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */

int LumiCalClustererClass::initialLowEngyClusterBuild( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsSmallEngyCellId,
						       std::map < int , IMPL::CalorimeterHitImpl* >	& calHitsCellId,
						       std::map < int , int >			& cellIdToClusterId,
						       std::map < int , std::vector<int> >		& clusterIdToCellId,
						       std::map < int , LCCluster > & clusterCM ) {

  int	cellIdHit, clusterId, numClusters, numElementsInLayer;

  std::map < int , int > :: iterator	cellIdToClusterIdIterator  ;
  std::map < int , std::vector<int> > :: iterator	clusterIdToCellIdIterator  ;
  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator calHitsCellIdIterator;
  std::vector < int > cellIdV, clusterIdV;

  double	CM1[3], CM2[3], distanceCM;
  std::map < int , double >			weightedDistanceV;
  std::map < int , double > :: iterator	weightedDistanceVIterator;

  /* --------------------------------------------------------------------------
     merge the unclustered cal hits with the existing clusters
     -------------------------------------------------------------------------- */
  calHitsCellIdIterator    = calHitsSmallEngyCellId.begin();
  numElementsInLayer       = calHitsSmallEngyCellId.size();
  for(int hitNow = 0; hitNow < numElementsInLayer; hitNow++, calHitsCellIdIterator++){
    cellIdHit = (int)(*calHitsCellIdIterator).first;

    // add the small energy hits that have now been clustred to the cal hit list at calHitsCellId
    calHitsCellId[cellIdHit] = calHitsSmallEngyCellId.at(cellIdHit);

    weightedDistanceV.clear();

    // position of the cal hit
    CM1[0] = calHitsSmallEngyCellId.at(cellIdHit) -> getPosition()[0];
    CM1[1] = calHitsSmallEngyCellId.at(cellIdHit) -> getPosition()[1];

    // compute weight for the cal hit and each cluster
    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
      clusterId = (int)(*clusterIdToCellIdIterator).first;

      CM2[0]     = clusterCM[clusterId].getX();
      CM2[1]     = clusterCM[clusterId].getY();
      distanceCM = distance2D(CM1,CM2);

      if(distanceCM > 0)
	weightedDistanceV[clusterId] = 1./distanceCM;
      else
	weightedDistanceV[clusterId] = 1e10;
    }

    // decide to which cluster to merge the hit according to a proper weight
    double maxClusterWeight = 0, clusterWeight;
    int   maxWeightClusterId(0);
    weightedDistanceVIterator = weightedDistanceV.begin();
    numClusters               = weightedDistanceV.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
      clusterId     = (int)(*weightedDistanceVIterator).first;
      clusterWeight = (double)(*weightedDistanceVIterator).second;

      if(maxClusterWeight < clusterWeight) {
	maxWeightClusterId = clusterId;
	maxClusterWeight   = clusterWeight;
      }
    }

    // add the hit to the cluster with the max weight
    if(maxClusterWeight > 0){
      clusterId = maxWeightClusterId;
      // add the hit to the chosen cluster
      clusterIdToCellId[clusterId].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = clusterId;

      //  update the totalEnergy counter and CM position of the cluster
      updateEngyPosCM(calHitsCellId[cellIdHit], clusterCM[clusterId]);

    }
  }

  return 1;
}





/* =========================================================================
   LumiCalClustererClass :: virtualCMClusterBuild
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */
int LumiCalClustererClass::virtualCMClusterBuild( std::map < int , IMPL::CalorimeterHitImpl* >	const& calHitsCellId,
						  std::map < int , int >			& cellIdToClusterId,
						  std::map < int , std::vector<int> >		& clusterIdToCellId,
						  std::map < int , LCCluster >	& clusterCM,
						  std::map < int , VirtualCluster >	const& virtualClusterCM ) {

  int	numElementsInLayer, cellIdHit;
  int	numClusters, clusterId, numUnClusteredElements;

  std::map < int , std::vector<int> > :: iterator		clusterIdToCellIdIterator  ;

  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator	calHitsCellIdIterator;

  std::map < int , VirtualCluster > :: const_iterator	virtualClusterCMIterator;

  std::vector < int >	cellIdV,
    unClusteredCellId;

  double	CM1[3], CM2[3], distanceCM;
  std::map < int , double >			weightedDistanceV;
  std::map < int , double > :: iterator	weightedDistanceVIterator;

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
  cout	<< endl << coutBlue << "Hits inside virtualRadius of virtualClusters" << coutDefault << endl << endl;
#endif

  /* --------------------------------------------------------------------------
     form clusters by gathering hits inside the virtual cluster radius
     -------------------------------------------------------------------------- */
  calHitsCellIdIterator = calHitsCellId.begin();
  numElementsInLayer    = (int)calHitsCellId.size();
  for(int hitNow = 0; hitNow < numElementsInLayer; hitNow++, calHitsCellIdIterator++){
    cellIdHit = (int)(*calHitsCellIdIterator).first;

    weightedDistanceV.clear();

    // compute the distance of the cal hit from the virtual cluster CMs, and keep score of
    // the clusters that the hit is in range of
    const IMPL::CalorimeterHitImpl *thisHit = calHitsCellId.at(cellIdHit);
    CM1[0] = thisHit -> getPosition()[0];
    CM1[1] = thisHit -> getPosition()[1];

    virtualClusterCMIterator = virtualClusterCM.begin();
    numClusters              = virtualClusterCM.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, virtualClusterCMIterator++){
      clusterId = (int)(*virtualClusterCMIterator).first;

      CM2[0]     = virtualClusterCM.at(clusterId).getX();
      CM2[1]     = virtualClusterCM.at(clusterId).getY();
      distanceCM = distance2D(CM1,CM2);

      if(distanceCM <= virtualClusterCM.at(clusterId).getZ()) {
	if(distanceCM > 0)
	  weightedDistanceV[clusterId] = 1./distanceCM;
	else
	  weightedDistanceV[clusterId] = 1e10;
      }

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< "\tclusterId: " << clusterId << " \t hit(x,y): " << CM1[0] << " \t " << CM1[1]
		<< " \t virtualCM(x,y):" << CM2[0] << " \t " << CM2[1] << endl
		<< "\t\tvirtualRadius: " << virtualClusterCM[clusterId][0] << " \t distance: "
		<< distanceCM << endl << endl;
#endif
    }

    // decide to which cluster to merge the hit according to a proper weight
    double	maxClusterWeight = 0, clusterWeight;
    int	maxWeightClusterId(0);
    weightedDistanceVIterator = weightedDistanceV.begin();
    numClusters                      = weightedDistanceV.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
      clusterId     = (int)(*weightedDistanceVIterator).first;
      clusterWeight = (double)(*weightedDistanceVIterator).second;

      if(maxClusterWeight < clusterWeight) {
	maxWeightClusterId = clusterId;
	maxClusterWeight   = clusterWeight;
      }
    }

    // create clusters from the hits which made the cut, and keep score of those which didnt
    if(maxClusterWeight > 0){
      clusterId = maxWeightClusterId;
      // add the hit to the chosen cluster
      clusterIdToCellId[clusterId].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = clusterId;

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< coutRed << "\tnum possible -  " << weightedDistanceV.size()
		<< " \t chosen - " << clusterId << coutDefault << endl << endl;
#endif
    } else	{
      unClusteredCellId.push_back( cellIdHit );

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< coutRed << "\tno cluster is within range ..." << coutDefault << endl << endl;
#endif
    }
  }


  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each cluster
     -------------------------------------------------------------------------- */
  clusterCM.clear();
  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
    clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
    cellIdV = (std::vector<int>)(*clusterIdToCellIdIterator).second;	// cal-hit-Ids in cluster
    //clusterCM[clusterId] = LCCluster();
    // initialize the energy/position std::vector for this cluster
    //for(int k=0; k<8; k++) clusterCM[clusterId].push_back(0.);
    // calculate the energy/position of the CM
    calculateEngyPosCM(cellIdV, calHitsCellId, clusterCM[clusterId], _methodCM);

  }
  // cleanUp
  cellIdV.clear();


  /* --------------------------------------------------------------------------
     if no clusters were created, make them into temporary clusters with
     zero energy (and no cal hits), so that hits will have a chance to merge with them
     -------------------------------------------------------------------------- */
  virtualClusterCMIterator = virtualClusterCM.begin();
  numClusters              = virtualClusterCM.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, virtualClusterCMIterator++){
    clusterId = (int)(*virtualClusterCMIterator).first;

    if(clusterIdToCellId[clusterId].size() == 0){
      clusterCM[clusterId] = LCCluster(virtualClusterCM.at(clusterId));
      // clusterCM[clusterId].push_back( 0. );
      // clusterCM[clusterId].push_back( virtualClusterCM[clusterId][1] );
      // clusterCM[clusterId].push_back( virtualClusterCM[clusterId][2] );
      // clusterCM[clusterId].push_back( virtualClusterCM[clusterId][3] );
      //APS: HUGE ERROR HERE???
      // clusterCM[clusterId].push_back( 0. );
      // clusterCM[clusterId].push_back( 0. );
    }
  }


  /* --------------------------------------------------------------------------
     merge unclustered cal hits with the existing clusters
     -------------------------------------------------------------------------- */
  numUnClusteredElements = unClusteredCellId.size();
  for(int hitNow = 0; hitNow < numUnClusteredElements; hitNow++){

    cellIdHit = unClusteredCellId[hitNow];

    weightedDistanceV.clear();

    // position of the cal hit
    CM1[0] = calHitsCellId.at(cellIdHit)->getPosition()[0];
    CM1[1] = calHitsCellId.at(cellIdHit)->getPosition()[1];

    // compute weight for the cal hit and each cluster
    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters = clusterIdToCellId.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
      clusterId = (int)(*clusterIdToCellIdIterator).first;

      CM2[0] = clusterCM[clusterId].getX();
      CM2[1] = clusterCM[clusterId].getY();
      distanceCM = distance2D(CM1,CM2);

      if(distanceCM > 0)
	weightedDistanceV[clusterId] = 1./distanceCM;
      else
	weightedDistanceV[clusterId] = 1e10;
    }

    // decide to which cluster to merge the hit according to a proper weight
    double	maxClusterWeight = 0, clusterWeight;
    int	maxWeightClusterId(0);
    weightedDistanceVIterator = weightedDistanceV.begin();
    numClusters               = weightedDistanceV.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
      clusterId     = (int)(*weightedDistanceVIterator).first;
      clusterWeight = (double)(*weightedDistanceVIterator).second;

      if(maxClusterWeight < clusterWeight) {
	maxWeightClusterId = clusterId;
	maxClusterWeight   = clusterWeight;
      }
    }

    // add the hit to the cluster with the max weight
    if(maxClusterWeight > 0){
      clusterId = maxWeightClusterId;
      // add the hit to the chosen cluster
      clusterIdToCellId[clusterId].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = clusterId;

      //  update the totalEnergy counter and CM position of the cluster
      updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId]);

      cellIdV.clear();
    }
  }

  return 1;
}




/* =========================================================================
   LumiCalClustererClass :: virtualCMPeakLayersFix
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */
int LumiCalClustererClass::virtualCMPeakLayersFix( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
						   std::map < int , int > & cellIdToClusterId,
						   std::map < int , std::vector<int> > & clusterIdToCellId,
						   std::map < int , LCCluster > & clusterCM,
						   std::map < int , VirtualCluster > virtualClusterCM ) {

  // general variables
  int	numElementsInLayer, cellIdHit;
  int	numClusters, clusterId, numElementsInCluster;
  int	virtualClusterId, numVirtualClusters;

  std::map < int , std::vector<int> > :: const_iterator		clusterIdToCellIdIterator;

  std::map <int , IMPL::CalorimeterHitImpl* > :: const_iterator	calHitsCellIdIterator;

  std::map < int , VirtualCluster > :: iterator	virtualClusterCMIterator;

  std::vector < int >	cellIdV;

  std::vector < std::vector <double> >	unClusteredCellId;

  std::map < int , int >		virtualToRealClusterId,
    oldVirtualClusterCM,
    oldClusterCM,
    oldToNewVirtualClusterIds;
  std::map < int , int > :: iterator	oldToNewVirtualClusterIdsIterator;

  double	CM1[3], CM2[3], distanceCM;
  std::map < int , double >			weightedDistanceV;
  std::map < int , double > :: iterator	weightedDistanceVIterator;

  /* --------------------------------------------------------------------------
     make sure that all the virtual cluster Ids are different than the
     existing real cluster Ids.
     -------------------------------------------------------------------------- */
  virtualClusterCMIterator = virtualClusterCM.begin();
  numVirtualClusters       = virtualClusterCM.size();
  for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, virtualClusterCMIterator++){
    virtualClusterId = (int)(*virtualClusterCMIterator).first;
    oldVirtualClusterCM[virtualClusterId] = 1;
  }

  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
    clusterId = (int)(*clusterIdToCellIdIterator).first;
    oldClusterCM[clusterId] = 1;
  }


  virtualClusterCMIterator = virtualClusterCM.begin();
  numVirtualClusters       = virtualClusterCM.size();
  for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, virtualClusterCMIterator++){
    virtualClusterId = (int)(*virtualClusterCMIterator).first;

    int newVirtualClusterId = virtualClusterId;
    int loopFlag = 1;
    while(loopFlag == 1){
      if( (oldClusterCM[newVirtualClusterId] == 1) || (oldVirtualClusterCM[newVirtualClusterId] == 1) )
	newVirtualClusterId++;
      else {
	loopFlag = 0;
	oldVirtualClusterCM[newVirtualClusterId] = 1;
      }
    }
    oldToNewVirtualClusterIds[virtualClusterId] = newVirtualClusterId;
  }

  oldToNewVirtualClusterIdsIterator = oldToNewVirtualClusterIds.begin();
  numVirtualClusters                = oldToNewVirtualClusterIds.size();
  for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, oldToNewVirtualClusterIdsIterator++){
    virtualClusterId = (int)(*oldToNewVirtualClusterIdsIterator).first;

    int newVirtualClusterId = oldToNewVirtualClusterIds[virtualClusterId];

    virtualClusterCM[newVirtualClusterId] = virtualClusterCM[virtualClusterId];
    virtualClusterCM.erase(virtualClusterId);

    // initialization for a later stage
    virtualToRealClusterId[newVirtualClusterId] = 0;
  }

  oldToNewVirtualClusterIds.clear(); oldVirtualClusterCM.clear(); oldClusterCM.clear();


  /* --------------------------------------------------------------------------
     associate the real clusters with virtual clusters and thus find the
     virtual clusters which dont have a real cluster of their own.
     -------------------------------------------------------------------------- */
  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
    clusterId = (int)(*clusterIdToCellIdIterator).first;

    weightedDistanceV.clear();

    CM1[0] = clusterCM[clusterId].getX();
    CM1[1] = clusterCM[clusterId].getY();

    virtualClusterCMIterator = virtualClusterCM.begin();
    numVirtualClusters       = virtualClusterCM.size();
    for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, virtualClusterCMIterator++){
      virtualClusterId = (int)(*virtualClusterCMIterator).first;

      CM2[0] = virtualClusterCM[virtualClusterId].getX();
      CM2[1] = virtualClusterCM[virtualClusterId].getY();

      distanceCM = distance2D(CM1,CM2);

      if(distanceCM > 0)
	weightedDistanceV[virtualClusterId] = 1./distanceCM;
      else
	weightedDistanceV[virtualClusterId] = 1e10;
    }

    // decide which virtualCluster to associate with the real cluster
    double	maxClusterWeight = 0, clusterWeight;
    int	maxWeightClusterId(0);
    weightedDistanceVIterator = weightedDistanceV.begin();
    numVirtualClusters        = weightedDistanceV.size();
    for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, weightedDistanceVIterator++){
      virtualClusterId = (int)(*weightedDistanceVIterator).first;
      clusterWeight    = (double)(*weightedDistanceVIterator).second;

      if(maxClusterWeight < clusterWeight) {
	maxWeightClusterId = virtualClusterId;
	maxClusterWeight   = clusterWeight;
      }
    }

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
    cout	<< "cluster " << clusterId << " is associated with virtualCluster " << maxWeightClusterId << endl;
#endif

    virtualToRealClusterId[maxWeightClusterId] = 1;
  }


  /* --------------------------------------------------------------------------
     go over all of the hits, and change the cluster id of the hits inside
     the empty virtual clusters to that of the virtual cluster
     -------------------------------------------------------------------------- */
  calHitsCellIdIterator = calHitsCellId.begin();
  numElementsInLayer    = (int)calHitsCellId.size();
  for(int j=0; j<numElementsInLayer; j++, calHitsCellIdIterator++){
    cellIdHit = (int)(*calHitsCellIdIterator).first;

    weightedDistanceV.clear();

    // compute the distance of the cal hit from the virtual cluster CMs
    CM1[0] = calHitsCellId.at(cellIdHit) -> getPosition()[0];
    CM1[1] = calHitsCellId.at(cellIdHit) -> getPosition()[1];

    virtualClusterCMIterator = virtualClusterCM.begin();
    numVirtualClusters       = virtualClusterCM.size();
    for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, virtualClusterCMIterator++){
      virtualClusterId = (int)(*virtualClusterCMIterator).first;

      if(virtualToRealClusterId[virtualClusterId] == 1) continue;

      CM2[0] = virtualClusterCM[virtualClusterId].getX();
      CM2[1] = virtualClusterCM[virtualClusterId].getY();
      distanceCM = distance2D(CM1,CM2);

      if(distanceCM <= virtualClusterCM[virtualClusterId].getZ()) {
	if(distanceCM > 0)
	  weightedDistanceV[virtualClusterId] = 1./distanceCM;
	else
	  weightedDistanceV[virtualClusterId] = 1e10;
      }
    }

    // decide to which virtualCluster to merge the hit according to a proper weight
    double	maxClusterWeight = 0, clusterWeight;
    int	maxWeightClusterId = -1;
    weightedDistanceVIterator = weightedDistanceV.begin();
    numVirtualClusters        = weightedDistanceV.size();
    for(int virtualClusterNow = 0; virtualClusterNow < numVirtualClusters; virtualClusterNow++, weightedDistanceVIterator++){
      virtualClusterId = (int)(*weightedDistanceVIterator).first;
      clusterWeight    = (double)(*weightedDistanceVIterator).second;

      if(maxClusterWeight < clusterWeight) {
	maxWeightClusterId = virtualClusterId;
	maxClusterWeight   = clusterWeight;
      }
    }

    if(numVirtualClusters == 0) continue;

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
    cout	<< " - num possible -  " << weightedDistanceV.size() << " \t id of chosen - (" << maxWeightClusterId
		<< ") \t virtualClusterId/hit(x,y): " <<  " \t " << virtualClusterCM[maxWeightClusterId][1] << " \t "
		<< virtualClusterCM[maxWeightClusterId][2]  <<endl;
#endif

    // add the hit to the chosen cluster
    clusterIdToCellId[maxWeightClusterId].push_back(cellIdHit);
    cellIdToClusterId[cellIdHit] = maxWeightClusterId;
  }


  /* --------------------------------------------------------------------------
     initialize and rewrite the clusterIdToCellId[] std::vectors for all clusters
     -------------------------------------------------------------------------- */
  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
    clusterId = (int)(*clusterIdToCellIdIterator).first;
    clusterIdToCellId[clusterId].clear();
  }

  calHitsCellIdIterator = calHitsCellId.begin();
  numElementsInLayer    = (int)calHitsCellId.size();
  for(int j=0; j<numElementsInLayer; j++, calHitsCellIdIterator++){
    cellIdHit = (int)(*calHitsCellIdIterator).first;

    clusterId = cellIdToClusterId[cellIdHit];
    clusterIdToCellId[clusterId].push_back(cellIdHit);
  }


  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each new cluster
     -------------------------------------------------------------------------- */
  clusterCM.clear();
  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int j=0; j<numClusters; j++, clusterIdToCellIdIterator++){
    clusterId  = (int)(*clusterIdToCellIdIterator).first;	// Id of cluster
    cellIdV = (std::vector<int>)(*clusterIdToCellIdIterator).second;	// cal-hit-Ids in cluster

    numElementsInCluster = clusterIdToCellId[clusterId].size();
    if(numElementsInCluster == 0) continue;
    // initialize the energy/position std::vector
    //clusterCM[clusterId] = LCCluster();
    // calculate/update the energy/position of the CM
    calculateEngyPosCM(cellIdV, calHitsCellId, clusterCM[clusterId], _methodCM);
  }
  cellIdV.clear();

  return 1;
}




/* =========================================================================
   LumiCalClustererClass :: buildSuperClusters
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */
int LumiCalClustererClass::buildSuperClusters ( std::map <int , IMPL::CalorimeterHitImpl* > & calHitsCellIdGlobal,
						std::vector < std::map < int , IMPL::CalorimeterHitImpl* > > const& calHitsCellIdLayer,
						std::vector < std::map < int , std::vector<int> > > const& clusterIdToCellId,
						std::vector < std::map < int , LCCluster > > const& clusterCM,
						std::vector < std::map < int , VirtualCluster > > const& virtualClusterCM,
						std::map < int , int > & cellIdToSuperClusterId,
						std::map < int , std::vector<int> > & superClusterIdToCellId,
						std::map < int , LCCluster > & superClusterCM ) {


  // general variables
  int	numClusters, numSuperClusters, numVirtualClusters, clusterId, superClusterId;
  int	numElementsInCluster;

  std::map < int , std::vector<int> > :: const_iterator	superClusterIdToCellIdIterator;


  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator	calHitsCellIdIterator;

  std::map < int , std::vector<int> > :: const_iterator		clusterIdToCellIdIterator;
  std::map < int , VirtualCluster > :: const_iterator	virtualClusterCMIterator;

  std::vector < int > cellIdV;

  double	CM1[3], CM2[3], distanceCM;
  std::map < int , double >			weightedDistanceV;
  std::map < int , double > :: iterator	weightedDistanceVIterator;

  std::vector < int >	reClusterHits;


  /* --------------------------------------------------------------------------
     merge all layer-clusters into global superClusters according to distance
     from the virtual clusters' CM
     -------------------------------------------------------------------------- */
#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
  cout	<< coutBlue << "Assign each cluster to a virtual cluster:" << coutDefault << endl;
#endif

  for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++ ){
    reClusterHits.clear();

    clusterIdToCellIdIterator = clusterIdToCellId[layerNow].begin();
    numClusters               = clusterIdToCellId[layerNow].size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
      clusterId = (int)(*clusterIdToCellIdIterator).first;

      weightedDistanceV.clear();

      CM1[0] = clusterCM.at(layerNow).at(clusterId).getX();
      CM1[1] = clusterCM.at(layerNow).at(clusterId).getY();

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< " - layer " << layerNow << " \t real cluster id,x,y,engy : "
		<< clusterId << " \t " << CM1[0] << " \t " << CM1[1]
		<< " \t " << clusterCM[layerNow][clusterId][0] << endl;
#endif

      virtualClusterCMIterator = virtualClusterCM.at(layerNow).begin();
      numVirtualClusters         = virtualClusterCM.at(layerNow).size();
      for(int superClusterNow = 0; superClusterNow < numVirtualClusters; superClusterNow++, virtualClusterCMIterator++){
	superClusterId = (int)(*virtualClusterCMIterator).first;

	CM2[0] = virtualClusterCM.at(layerNow).at(superClusterId).getX();
	CM2[1] = virtualClusterCM.at(layerNow).at(superClusterId).getY();

	distanceCM = distance2D(CM1,CM2);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
	cout	<<  "\t\t virtual cluster id,x,y : "
		<< superClusterId << " \t " << CM2[0] << " \t " << CM2[1]
		<< " \t distanceCM : " << distanceCM << endl;
#endif

	if(distanceCM > 0)
	  weightedDistanceV[superClusterId] = 1./distanceCM;
	else
	  weightedDistanceV[superClusterId] = 1e10;
      }

      // decide to which superCluster to merge the cluster according to a proper weight
      double	maxClusterWeight = 0, clusterWeight;
      int	maxWeightClusterId;
      weightedDistanceVIterator = weightedDistanceV.begin();
      numVirtualClusters        = weightedDistanceV.size();
      for(int superClusterNow = 0; superClusterNow < numVirtualClusters; superClusterNow++, weightedDistanceVIterator++){
	superClusterId = (int)(*weightedDistanceVIterator).first;
	clusterWeight  = (double)(*weightedDistanceVIterator).second;

	if(maxClusterWeight < clusterWeight) {
	  maxWeightClusterId = superClusterId;
	  maxClusterWeight   = clusterWeight;
	}
      }

      // make sure that all clusters have been added to superclusters
      if(maxClusterWeight <= 0)	return 0;

      /* --------------------------------------------------------------------------
	 in the case where there are more than 1 virtualClusters and the real
	 cluster is not close to any virtualCM, than instead of assigning the whole
	 real cluster to a virtualCluster, the real cluster will be dismanteled
	 and each hit will be assigned to a virtualCluster separatly
	 -------------------------------------------------------------------------- */
      numVirtualClusters = virtualClusterCM[layerNow].size();
      distanceCM         = 1./maxClusterWeight;
      if(distanceCM > _moliereRadius && numVirtualClusters > 1) {
	reClusterHits.push_back(clusterId);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
	cout	<< coutRed << "\t\t no virtual cluster is chosen... " << coutDefault << endl << endl;;
#endif
	continue;
      }

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< coutRed << "\t\t chosen virtual cluster " <<  maxWeightClusterId << coutDefault << endl << endl;;
#endif

      // go over all hits in the cluster and add to the chosen superCluster
      numElementsInCluster = clusterIdToCellId.at(layerNow).at(clusterId).size();
      for(int j=0; j<numElementsInCluster; j++){
	// add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
	const int cellIdHit = clusterIdToCellId.at(layerNow).at(clusterId)[j];

	superClusterIdToCellId[maxWeightClusterId].push_back(cellIdHit);
	cellIdToSuperClusterId[cellIdHit] = maxWeightClusterId;
      }
    }


    /* --------------------------------------------------------------------------
       real clusters far away from virtualCMs are dismanteled and each hit is
       assigned to a virtualCluster separatly according to distance
       -------------------------------------------------------------------------- */
    numClusters = reClusterHits.size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++){
      clusterId = reClusterHits[clusterNow];

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout	<< coutBlue << "\t - dismantel cluster " << clusterId
		<< " and assign each hit to a virtualCluster "<< coutDefault << endl ;
#endif

      numElementsInCluster = clusterIdToCellId.at(layerNow).at(clusterId).size();
      for(int j=0; j<numElementsInCluster; j++){
	// add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
	const int cellIdHit = clusterIdToCellId.at(layerNow).at(clusterId).at(j);

	weightedDistanceV.clear();

	// position of the cal hit
	CM1[0] = calHitsCellIdLayer.at(layerNow).at(cellIdHit)->getPosition()[0];
	CM1[1] = calHitsCellIdLayer.at(layerNow).at(cellIdHit)->getPosition()[1];

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
	cout	<< "\t cellId " << cellIdHit << " \t x,y : " << CM1[0] << " \t " << CM1[1] << endl;
#endif

	virtualClusterCMIterator = virtualClusterCM[layerNow].begin();
	numVirtualClusters       = virtualClusterCM[layerNow].size();
	for(int superClusterNow = 0; superClusterNow < numVirtualClusters; superClusterNow++, virtualClusterCMIterator++){
	  superClusterId = (int)(*virtualClusterCMIterator).first;

	  CM2[0] = virtualClusterCM.at(layerNow).at(superClusterId).getX();
	  CM2[1] = virtualClusterCM.at(layerNow).at(superClusterId).getY();

	  distanceCM = distance2D(CM1,CM2);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
	  cout	<<  "\t\t virtual cluster id,x,y : "
		<< superClusterId << " \t " << CM2[0] << " \t " << CM2[1]
		<< " \t distanceCM : " << distanceCM << endl;
#endif

	  if(distanceCM > 0)
	    weightedDistanceV[superClusterId] = 1./distanceCM;
	  else
	    weightedDistanceV[superClusterId] = 1e10;
	}

	// decide to which superCluster to merge the cluster according to a proper weight
	double	maxClusterWeight = 0, clusterWeight;
	int	maxWeightClusterId;
	weightedDistanceVIterator = weightedDistanceV.begin();
	numVirtualClusters        = weightedDistanceV.size();
	for(int superClusterNow = 0; superClusterNow < numVirtualClusters; superClusterNow++, weightedDistanceVIterator++){
	  superClusterId = (int)(*weightedDistanceVIterator).first;
	  clusterWeight  = (double)(*weightedDistanceVIterator).second;

	  if(maxClusterWeight < clusterWeight) {
	    maxWeightClusterId = superClusterId;
	    maxClusterWeight   = clusterWeight;
	  }
	}

	// make sure that all clusters have been added to superclusters
	assert	(maxClusterWeight > 0);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
	cout	<< coutRed << "\t\t chosen virtual cluster " <<  maxWeightClusterId << coutDefault << endl << endl;;
#endif

	superClusterIdToCellId[maxWeightClusterId].push_back(cellIdHit);
	cellIdToSuperClusterId[cellIdHit] = maxWeightClusterId;
      }
    }


    // write all the layer cal hits from calHitsCellIdLayer to the global calHitsCellId map
    std::map < int , IMPL::CalorimeterHitImpl* > const& layerHits = calHitsCellIdLayer.at(layerNow);
    calHitsCellIdIterator = layerHits.begin();
    std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator	calHitsEnd = layerHits.end();
    for(; calHitsCellIdIterator  != calHitsEnd; ++calHitsCellIdIterator){
      const int cellIdHit = calHitsCellIdIterator->first;
      calHitsCellIdGlobal[cellIdHit] = calHitsCellIdIterator->second;
    }
  }


  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each new superCluster
     -------------------------------------------------------------------------- */
  superClusterCM.clear();
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters               = superClusterIdToCellId.size();
  for(int j=0; j<numSuperClusters; j++, superClusterIdToCellIdIterator++){
    superClusterId  = (int)(*superClusterIdToCellIdIterator).first;  // Id of cluster
    cellIdV = (std::vector<int>)(*superClusterIdToCellIdIterator).second; // cal-hit-Ids in cluster

    // initialize the energy/position std::vector for new clusters only
    //superClusterCM[superClusterId] = LCCluster();

    // calculate/update the energy/position of the CM
    calculateEngyPosCM(cellIdV, calHitsCellIdGlobal, superClusterCM[superClusterId], _methodCM);
  }
  cellIdV.clear();


#if _CLUSTER_BUILD_DEBUG == 1
  streamlog_out( DEBUG ) << " - superClusters:"  << std::endl;

  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters               = superClusterIdToCellId.size();
  for(int superClusterNow=0; superClusterNow<numSuperClusters; superClusterNow++, superClusterIdToCellIdIterator++){
    superClusterId  = (int)(*superClusterIdToCellIdIterator).first;	// Id of cluster

    streamlog_out( DEBUG ) << "\t Id " << superClusterId << "  \t  energy " << superClusterCM[superClusterId].getE()
	      << "     \t pos(x,y) =  ( " << superClusterCM[superClusterId].getX()
	      << " , " << superClusterCM[superClusterId].getY() << " )"
	      << "     \t pos(theta,phi) =  ( " << superClusterCM[superClusterId].getTheta()
	      << " , " << superClusterCM[superClusterId].getPhi() << " )"
	      << std::endl;

  }
#endif

  return 1;

}


/* =========================================================================
   LumiCalClustererClass :: engyInMoliereCorrections
   ============================================================================
   (1). Description:
   --------------------------------
   - SOME DESCRIPTION ......
   ============================================================================ */
int LumiCalClustererClass::engyInMoliereCorrections ( std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellIdGlobal,
						      std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > const&	calHits,
						      std::vector < std::map < int , IMPL::CalorimeterHitImpl* > > const& ,//	calHitsCellIdLayer,
						      std::vector < std::map < int , std::vector<int> > > & clusterIdToCellId,
						      std::vector < std::map < int , LCCluster > > & clusterCM,
						      std::vector < std::map < int , int > > & cellIdToClusterId,
						      std::map < int , int > & cellIdToSuperClusterId,
						      std::map < int , std::vector<int> > & superClusterIdToCellId,
						      std::map < int , LCCluster > & superClusterCM,
						      double middleEnergyHitBound,
						      int // detectorArm
						      )  {

  // ???????? DECIDE/FIX - incorparate the parameter given here better in the code ????????
  int	engyHitBoundMultiply	= 1;
  double	molRadPercentage	= .4;
  double	engyPercentInMolFrac	= .8;
  double	minEngyPercentInMol	= .5;
  double	midEngyPercentInMol	= .8;
  double	baseEngyPercentInMol	= .9;

  // general variables
  int	numClusters, cellIdHit, superClusterId, numElementsInSuperCluster;
  int	numElementsInArm, numHitsProjection, numSuperClusters, rejectFlag;
  double	totEngyArmAboveMin, totEngyInAllMol;
  double	superClusterMolRatio = 0., superClusterMolRatio_Tmp = 0., projectionClusterMolRatio = 0.;

  std::map < int , IMPL::CalorimeterHitImpl* > calHitsCellIdProjection, calHitsCellIdProjectionFull;

  std::map < int , int > cellIdToSuperClusterId_Tmp;
  std::map < int , std::vector <int> > superClusterIdToCellId_Tmp;
  std::map < int , std::vector <int> > :: const_iterator	superClusterIdToCellIdIterator;

  std::map < int , std::vector <int> > :: iterator	clusterIdToCellIdIterator  ;

  std::map < int , LCCluster > superClusterCM_Tmp;
  std::map < int , LCCluster > :: iterator	superClusterCMIterator;

  std::map < int , LCCluster > :: iterator	clusterCMIterator;

  std::map < int , double >	superClusterEngyInMoliere,
    superClusterEngyInMoliere_Tmp,
    projectionClusterEngyInMoliere;

  std::map < int , std::vector<double> >			calHitsProjection,
    calHitsProjectionFull;
  std::map < int , std::vector<double> > :: iterator	calHitsProjectionIterator;

  std::map < int , int >	projectionFlag;

  std::map < int , double >			weightedDistanceV;
  std::map < int , double > :: iterator	weightedDistanceVIterator;

  std::vector < int >	initialClusterControlVar(4),
    superClusterRejected,
    superClusterAccepted,
    cellIdV;

#if _MOL_RAD_CORRECT_DEBUG == 1
  cout	<< coutBlue << "Before Moliere corrections:" << coutDefault << endl;
#endif

  totEngyArmAboveMin = 0.;
  totEngyInAllMol    = 0.;

  /* --------------------------------------------------------------------------
     find the percentage of energy for each cluster within _moliereRadius
     -------------------------------------------------------------------------- */
  superClusterCMIterator = superClusterCM.begin();
  numSuperClusters       = superClusterCM.size();
  for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
    superClusterId = (int)(*superClusterCMIterator).first;

    superClusterEngyInMoliere[superClusterId]
      = getEngyInMoliereFraction(	calHitsCellIdGlobal,
					superClusterIdToCellId[superClusterId],
					superClusterCM[superClusterId],
					1.  );

    totEngyInAllMol       += superClusterEngyInMoliere[superClusterId];
    totEngyArmAboveMin    += superClusterCM[superClusterId].getE();

#if _MOL_RAD_CORRECT_DEBUG == 1
    double engyPercentInMol = superClusterEngyInMoliere[superClusterId]
      / superClusterCM[superClusterId].getE() ;

    cout	<< "\t Id " << superClusterId << "  \t energy " << superClusterCM[superClusterId].getE()
		<< "     \t pos(x,y) =  ( " << superClusterCM[superClusterId].getX()
		<< " , " << superClusterCM[superClusterId].getY() << " )"
		<< "     \t pos(theta,phi) =  ( " << superClusterCM[superClusterId].getTheta()
		<< " , " << superClusterCM[superClusterId].getPhi() << " )" << endl
		<< "\t\t engy in _moliereRadius  \t=   " << superClusterEngyInMoliere[superClusterId]
		<< coutRed << " \t -> totEngy in Moliere = \t " << engyPercentInMol << " %" << coutDefault << endl << endl;
#endif
  }

  superClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

#if _MOL_RAD_CORRECT_DEBUG == 1
  cout	<< "\t (-) Tot engy in arm = " << totEngyArmAboveMin
	<< " , tot engy in all super clusters in MolRad = " << totEngyInAllMol
	<< coutRed << "  -> their ratio = \t " << superClusterMolRatio << coutDefault << endl << endl;
#endif

  /* --------------------------------------------------------------------------
     set a global variable that warns about an unstable clustering, in
     the case that not enough energy is found within _moliereRadius around
     the superClusters:
     (A). reject if below baseEngyPercentInMol of the energy is in moliere
     for all clusters combined
     (B). reject if below midEngyPercentInMol of the energy is in moliere
     for each cluster seperatly
     (*). also store the id's of clusters with very below minEngyPercentInMol
     energy in moliere, for force merging later, in case the profile
     clusters arent accepted.
     -------------------------------------------------------------------------- */
  if(superClusterMolRatio < baseEngyPercentInMol)	rejectFlag = 1;
  else						rejectFlag = 0;

  superClusterCMIterator = superClusterCM.begin();
  numSuperClusters       = superClusterCM.size();
  for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
    superClusterId = (int)(*superClusterCMIterator).first;

    double engyPercentInMol = superClusterEngyInMoliere[superClusterId]
      / superClusterCM[superClusterId].getE() ;
    if(engyPercentInMol < midEngyPercentInMol) rejectFlag = 1;
    if(engyPercentInMol < minEngyPercentInMol)
      superClusterRejected.push_back(superClusterId);
    else
      superClusterAccepted.push_back(superClusterId);
  }


  /* --------------------------------------------------------------------------
     in the case that not enough energy is found within _moliereRadius around
     each one of the superClusters, try to build projectionClusters.
     a projection layer is created out of all the hits of the real layers.
     clusters are built from high energy hits and the energy distribution is
     measured around the projectionCluster CMs.
     -------------------------------------------------------------------------- */
  if(rejectFlag == 1) {
    /* --------------------------------------------------------------------------
       sup up the energy for each Phi/R cell for all Z layers
       -------------------------------------------------------------------------- */
    //	  for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) {
    std::map < int , std::vector <IMPL::CalorimeterHitImpl*> >::const_iterator calHitsIt = calHits.begin(),
      calHitsEnd = calHits.end();
    //	    for ( auto const& layerHits: calHits) {
    for (; calHitsIt != calHitsEnd; ++calHitsIt) {
      std::pair < int , std::vector < IMPL::CalorimeterHitImpl*> > const& layerHits = (*calHitsIt);
      int numElementsInLayer = (int)layerHits.second.size();
      //const int layerNow =  layerHits.first;
      for(int j=0; j<numElementsInLayer; j++){
	IMPL::CalorimeterHitImpl const* thisCalHit = layerHits.second.at(j);
	cellIdHit = (int)thisCalHit->getCellID0();

	double cellEngy = (double)thisCalHit->getEnergy();
	///APS: This encoding needs to be fixed, now using
	// get the new cellId for the projection
	int	cellIdHitZ   = _maxLayerToAnalyse + 1;
	int	cellIdHitPhi = GlobalMethodsClass::CellIdZPR(cellIdHit, GlobalMethodsClass::COP);
	int	cellIdHitR   = GlobalMethodsClass::CellIdZPR(cellIdHit, GlobalMethodsClass::COR);

	int	cellIdProjection = GlobalMethodsClass::CellIdZPR(cellIdHitZ, cellIdHitPhi, cellIdHitR);

	// the original hit's layer number is stored in (the previously unused) CellID1
	cellIdHitZ = (cellIdHit >> 0 ) & 0xff;
	// only high energy projection hits will be considered in the first stage
	if(cellEngy > middleEnergyHitBound * engyHitBoundMultiply ) {
	  std::vector<double> & thisVector = calHitsProjection[cellIdProjection];
	  if(thisVector.empty()) {
	    thisVector.reserve(5);
	    thisVector.push_back( thisCalHit -> getEnergy() );
	    thisVector.push_back( thisCalHit -> getPosition()[0] );
	    thisVector.push_back( thisCalHit -> getPosition()[1] );
	    thisVector.push_back( thisCalHit -> getPosition()[2] );
	    thisVector.push_back( double(cellIdHitZ) );
	  } else {
	    thisVector[0] += thisCalHit -> getEnergy();
	  }
	}
	// store all hits above _hitMinEnergy energy
	if(cellEngy > _hitMinEnergy ) {
	  std::vector<double> & thisVector = calHitsProjectionFull[cellIdProjection];
	  if(thisVector.empty()){
	    thisVector.reserve(5);
	    thisVector.push_back( thisCalHit -> getEnergy() );
	    thisVector.push_back( thisCalHit -> getPosition()[0] );
	    thisVector.push_back( thisCalHit -> getPosition()[1] );
	    thisVector.push_back( thisCalHit -> getPosition()[2] );
	    thisVector.push_back( double(cellIdHitZ) );
	  } else {
	    thisVector[0] += thisCalHit -> getEnergy();
	  }
	}
      }
    }

    /* --------------------------------------------------------------------------
       input the results into new cal hit objects
       -------------------------------------------------------------------------- */
    IMPL::CalorimeterHitImpl * calHitNew;

    calHitsProjectionIterator = calHitsProjection.begin();
    numHitsProjection     = calHitsProjection.size();
    for(int hitNow = 0; hitNow < numHitsProjection; hitNow++, calHitsProjectionIterator++ ){
      int cellIdProjection = (int)(*calHitsProjectionIterator).first;

      double	engyCellProjection = calHitsProjection[cellIdProjection][0];

      float	hitPosV[3];
      hitPosV[0] = calHitsProjection[cellIdProjection][1];
      hitPosV[1] = calHitsProjection[cellIdProjection][2];
      hitPosV[2] = calHitsProjection[cellIdProjection][3];

      int	cellIdHitZ = (int)calHitsProjection[cellIdProjection][4];

      calHitNew = new IMPL::CalorimeterHitImpl();
      calHitNew -> setCellID0(cellIdProjection);
      calHitNew -> setCellID1(cellIdHitZ);
      calHitNew -> setEnergy(engyCellProjection);
      calHitNew -> setPosition(hitPosV);

      calHitsCellIdProjection[cellIdProjection] = calHitNew;
    }

    /* --------------------------------------------------------------------------
       build clusters out of the projection hits
       -------------------------------------------------------------------------- */
    // set the control std::vector for the initialClusterBuild clustering options
    initialClusterControlVar[0] = 1;  // mergeOneHitClusters
    initialClusterControlVar[1] = 1;  // mergeSmallToLargeClusters
    initialClusterControlVar[2] = 1;  // mergeLargeToSmallClusters
    initialClusterControlVar[3] = 1;  // forceMergeSmallToLargeClusters

    initialClusterBuild(	calHitsCellIdProjection,
				cellIdToClusterId[_maxLayerToAnalyse],
				clusterIdToCellId[_maxLayerToAnalyse],
				clusterCM[_maxLayerToAnalyse],
				initialClusterControlVar   );

    /* --------------------------------------------------------------------------
       find the percentage of energy for each cluster within _moliereRadius
       -------------------------------------------------------------------------- */
#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << coutBlue << "Projection clusters (initial clustering with all hits):" << coutDefault << endl;
#endif

    int engyInMoliereFlag = 0;
    clusterCMIterator = clusterCM[_maxLayerToAnalyse].begin();
    numClusters       = clusterCM[_maxLayerToAnalyse].size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterCMIterator++) {
      int clusterId = (int)(*clusterCMIterator).first;

      projectionClusterEngyInMoliere[clusterId]
	= getEngyInMoliereFraction(	calHitsCellIdProjection,
					clusterIdToCellId[_maxLayerToAnalyse][clusterId],
					clusterCM[_maxLayerToAnalyse][clusterId],
					molRadPercentage  );

      double engyPercentInMol = projectionClusterEngyInMoliere[clusterId]
	/ clusterCM[_maxLayerToAnalyse][clusterId].getE();

      if(engyPercentInMol < engyPercentInMolFrac) engyInMoliereFlag = 1;

#if _MOL_RAD_CORRECT_DEBUG == 1
      cout	<< "\tprojection " << clusterId << "  at (x,y) = (" << clusterCM[_maxLayerToAnalyse][clusterId][1]
		<< " , " << clusterCM[_maxLayerToAnalyse][clusterId][2] << ")   \t engy in (" << molRadPercentage
		<< " * _moliereRadius)  =   " << projectionClusterEngyInMoliere[clusterId]
		<< coutRed << " \t-> % totEngy = \t " << engyPercentInMol << coutDefault << endl;
#endif
    }


    /* --------------------------------------------------------------------------
       optimise the projection clusters (if needed).
       stop when there is a percentage of engyPercentInMolFrac total energy within
       a percentage molRadPercentage od moliere radius (each iteration the total
       energy decreases due to a higher cut on minHitEnergy).
       -------------------------------------------------------------------------- */
#if _MOL_RAD_CORRECT_DEBUG == 1
    if(engyInMoliereFlag == 1)
      cout	<< endl << coutBlue << "Optimizing ...  -  ignore hits with energy below    " << coutDefault;
#endif

    while(engyInMoliereFlag == 1){

      engyInMoliereFlag     = 0;
      engyHitBoundMultiply += 30;

#if _MOL_RAD_CORRECT_DEBUG == 1
      cout	<< "  ->  " << middleEnergyHitBound * engyHitBoundMultiply ;
#endif

      // remove hits that are of low energy
      std::vector <int> idsToErase;
      std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator calHitsCellIdIterator
	= calHitsCellIdProjection.begin();
      numHitsProjection = calHitsCellIdProjection.size();
      for(int hitNow = 0; hitNow < numHitsProjection; hitNow++, calHitsCellIdIterator++ ){
	int	cellIdProjection = (int)(*calHitsCellIdIterator).first;

	double engyHit = (double)calHitsCellIdProjection[cellIdProjection]->getEnergy();
	if(engyHit < middleEnergyHitBound * engyHitBoundMultiply)
	  idsToErase.push_back(cellIdProjection);
      }

      int numIdsToErase = idsToErase.size();
      int numHitsRemaining = (int)calHitsCellIdProjection.size() - numIdsToErase;
      if(numHitsRemaining < 5) {
#if _MOL_RAD_CORRECT_DEBUG == 1
	cout << endl << coutRed << "  -- optimization of the projection clusters has failed ... --" << coutDefault << endl;
#endif
	break;
      }

      for(int hitNow = 0; hitNow < numIdsToErase; hitNow++){
	int idsToEraseNow = idsToErase[hitNow];

	// cleanUp dynamicly allocated memory
	delete calHitsCellIdProjection[idsToEraseNow];
	// erase entry from map
	calHitsCellIdProjection.erase(idsToEraseNow);
      }


      /* --------------------------------------------------------------------------
	 build clusters out of the projection hits
	 -------------------------------------------------------------------------- */
      // clean up the clustering results from the previous run
      cellIdToClusterId[_maxLayerToAnalyse].clear();
      clusterIdToCellId[_maxLayerToAnalyse].clear();
      clusterCM[_maxLayerToAnalyse].clear();

      // set the control std::vector for the initialClusterBuild clustering options
      initialClusterControlVar[0] = 1;  // mergeOneHitClusters
      initialClusterControlVar[1] = 1;  // mergeSmallToLargeClusters
      initialClusterControlVar[2] = 1;  // mergeLargeToSmallClusters
      initialClusterControlVar[3] = 1;  // forceMergeSmallToLargeClusters

      initialClusterBuild(	calHitsCellIdProjection,
				cellIdToClusterId[_maxLayerToAnalyse],
				clusterIdToCellId[_maxLayerToAnalyse],
				clusterCM[_maxLayerToAnalyse],
				initialClusterControlVar );


      /* --------------------------------------------------------------------------
	 find the percentage of energy for each cluster within _moliereRadius
	 -------------------------------------------------------------------------- */
      clusterCMIterator = clusterCM[_maxLayerToAnalyse].begin();
      numClusters       = clusterCM[_maxLayerToAnalyse].size();
      for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterCMIterator++) {
	int clusterIdHit = (int)(*clusterCMIterator).first;

	projectionClusterEngyInMoliere[clusterIdHit]
	  =  getEngyInMoliereFraction(	calHitsCellIdProjection,
					clusterIdToCellId[_maxLayerToAnalyse][clusterIdHit],
					clusterCM[_maxLayerToAnalyse][clusterIdHit],
					molRadPercentage);

	double engyPercentInMol = projectionClusterEngyInMoliere[clusterIdHit]
	  / clusterCM[_maxLayerToAnalyse][clusterIdHit].getE();

	if(engyPercentInMol < engyPercentInMolFrac) engyInMoliereFlag = 1;
      }
    }

    /* --------------------------------------------------------------------------
       since the optimization of the projection clusters involved deleting of low
       energy cal hits, these need to be re-registered so as to calculate the
       total energy around the CM within _moliereRadius.
       -------------------------------------------------------------------------- */
    calHitsProjectionIterator = calHitsProjectionFull.begin();
    numHitsProjection         = calHitsProjectionFull.size();
    for(int hitNow = 0; hitNow < numHitsProjection; hitNow++, calHitsProjectionIterator++ ){
      int	cellIdProjection = (int)(*calHitsProjectionIterator).first;

      double	engyCellProjection = calHitsProjectionFull[cellIdProjection][0];
      float	hitPosV[3];
      hitPosV[0] = calHitsProjectionFull[cellIdProjection][1];
      hitPosV[1] = calHitsProjectionFull[cellIdProjection][2];
      hitPosV[2] = calHitsProjectionFull[cellIdProjection][3];

      calHitNew = new IMPL::CalorimeterHitImpl();
      calHitNew -> setCellID0(cellIdProjection);
      calHitNew -> setEnergy(engyCellProjection);
      calHitNew -> setPosition(hitPosV);

      calHitsCellIdProjectionFull[cellIdProjection] = calHitNew;
      // a flag map for avoiding double counting of clustered cells
      projectionFlag[cellIdProjection] = 0;
    }


#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << endl << coutBlue << "Projection clusters (all hits):" << coutDefault << endl;
#endif

    /* --------------------------------------------------------------------------
       find the percentage of energy for each cluster within _moliereRadius
       -------------------------------------------------------------------------- */
    totEngyInAllMol = 0.;

    clusterCMIterator = clusterCM[_maxLayerToAnalyse].begin();
    numClusters       = clusterCM[_maxLayerToAnalyse].size();
    for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterCMIterator++) {
      int clusterIdHit = (int)(*clusterCMIterator).first;

      projectionClusterEngyInMoliere[clusterIdHit]
	= getEngyInMoliereFraction( calHitsCellIdProjectionFull,
				    clusterIdToCellId[_maxLayerToAnalyse][clusterIdHit],
				    clusterCM[_maxLayerToAnalyse][clusterIdHit],
				    1.,
				    projectionFlag    );

      totEngyInAllMol += projectionClusterEngyInMoliere[clusterIdHit];

#if _MOL_RAD_CORRECT_DEBUG == 1
      cout	<< "\tfull-Projection " << clusterIdHit << " \tat (x,y) = ("
		<< clusterCM[_maxLayerToAnalyse][clusterIdHit][1] << " , "
		<< clusterCM[_maxLayerToAnalyse][clusterIdHit][2] << ")   \t engy in _moliereRadius  \t=   "
		<< projectionClusterEngyInMoliere[clusterIdHit] << " (possibly under-counted...)"<<endl;
#endif
    }

    projectionClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << "\t (-) Tot engy in arm = " << totEngyArmAboveMin << " , tot engy in all projection clusters  in MolRad = "
		<< totEngyInAllMol << coutRed << "  -> their ratio = \t " << projectionClusterMolRatio << coutDefault << endl << endl;
#endif
  }


  /* --------------------------------------------------------------------------
     if the projectionClusters (after possible optimization) have enough energy
     within _moliereRadius around their CMs, and the superClusters don't,
     than dump the superClusters and re-cluster around the projectionClusters.
     this time a very rudementary clustering algorithm is used.
     -------------------------------------------------------------------------- */
  if(superClusterMolRatio < projectionClusterMolRatio) {
    // create new clusters around the projection CMs
    std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator
      calHitsCellIdIterator = calHitsCellIdGlobal.begin();
    numElementsInArm      = calHitsCellIdGlobal.size();
    for(int hitNow = 0; hitNow < numElementsInArm; hitNow++, calHitsCellIdIterator++) {
      cellIdHit = (int)(*calHitsCellIdIterator).first;

      weightedDistanceV.clear();
      const IMPL::CalorimeterHitImpl * thisHit = calHitsCellIdGlobal.at(cellIdHit);

      clusterCMIterator = clusterCM[_maxLayerToAnalyse].begin();
      numClusters       = clusterCM[_maxLayerToAnalyse].size();
      for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterCMIterator++) {
	int clusterId = (int)(*clusterCMIterator).first;

	const double distanceCM = distance2D(thisHit->getPosition(),clusterCM[_maxLayerToAnalyse][clusterId].getPosition());

	if(distanceCM > 0)
	  weightedDistanceV[clusterId] = 1./distanceCM;
	else
	  weightedDistanceV[clusterId] = 1e10;
	//				cout << "cellId/hit(x,y): " << cellIdHit << " \t " << CM1[0] << " \t " << CM1[1]
	//					<< " \t clusterId/hit(x,y): " << clusterId << " \t " << CM2[0] << " \t " << CM2[1]
	//					<< " \t distance: " << distanceCM <<endl;
      }

      // decide to which superCluster to merge the cluster according to a proper weight
      double maxClusterWeight = 0, clusterWeight;
      int   maxWeightClusterId(0);
      weightedDistanceVIterator = weightedDistanceV.begin();
      numClusters               = weightedDistanceV.size();
      for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
	int clusterId     = (int)(*weightedDistanceVIterator).first;
	clusterWeight = (double)(*weightedDistanceVIterator).second;

	if(maxClusterWeight < clusterWeight) {
	  maxWeightClusterId = clusterId;
	  maxClusterWeight   = clusterWeight;
	}
      }
      //			cout	<< " - num possible -  " << clusterCM[_maxLayerToAnalyse].size()
      //				<< " \t id of chosen - " << maxWeightClusterId << endl <<endl;

      // make sure that all hits have have been added to clusters
      assert	(maxClusterWeight > 0);

      // add the hit to the chosen superCluster
      superClusterId = maxWeightClusterId;
      superClusterIdToCellId_Tmp[superClusterId].push_back(cellIdHit);
      cellIdToSuperClusterId_Tmp[cellIdHit] = superClusterId;
    }

    /* --------------------------------------------------------------------------
       compute total energy and center of mass of each new superCluster
       -------------------------------------------------------------------------- */
    superClusterCM_Tmp.clear();
    superClusterIdToCellIdIterator  = superClusterIdToCellId_Tmp.begin();
    numSuperClusters               = superClusterIdToCellId_Tmp.size();
    for(int j=0; j<numSuperClusters; j++, superClusterIdToCellIdIterator++){
      superClusterId  = (int)(*superClusterIdToCellIdIterator).first;  // Id of cluster
      cellIdV = (std::vector<int>)(*superClusterIdToCellIdIterator).second; // cal-hit-Ids in cluster

      // initialize the energy/position std::vector for new clusters only
      //superClusterCM_Tmp[superClusterId] = LCCluster();

      // calculate/update the energy/position of the CM
      calculateEngyPosCM(cellIdV, calHitsCellIdGlobal, superClusterCM_Tmp[superClusterId], _methodCM);
    }
    cellIdV.clear();


    /* --------------------------------------------------------------------------
       find the percentage of energy for each new superCluster within _moliereRadius
       -------------------------------------------------------------------------- */
    totEngyArmAboveMin = 0.;
    totEngyInAllMol    = 0.;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << coutBlue << "New superClusters:" << coutDefault << endl;
#endif

    superClusterCMIterator = superClusterCM_Tmp.begin();
    numSuperClusters       = superClusterCM_Tmp.size();
    for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
      superClusterId = (int)(*superClusterCMIterator).first;

      superClusterEngyInMoliere_Tmp[superClusterId]
	= getEngyInMoliereFraction(	calHitsCellIdGlobal,
					superClusterIdToCellId_Tmp[superClusterId],
					superClusterCM_Tmp[superClusterId],
					1.  );

      totEngyInAllMol       += superClusterEngyInMoliere_Tmp[superClusterId];
      totEngyArmAboveMin    += superClusterCM_Tmp[superClusterId].getE();

#if _MOL_RAD_CORRECT_DEBUG == 1
      double engyPercentInMol = superClusterEngyInMoliere_Tmp[superClusterId]
	/ superClusterCM_Tmp[superClusterId][0] ;

      cout << "superCluster " << superClusterId << " \tat (x,y) = (" << superClusterCM_Tmp[superClusterId][1]
	   << " , " << superClusterCM_Tmp[superClusterId][2] << ")   \t engy in _moliereRadius  \t=   "
	   << superClusterEngyInMoliere_Tmp[superClusterId]
	   << coutRed << " \t-> % totEngy = \t " << engyPercentInMol << coutDefault << endl;
#endif
    }

    superClusterMolRatio_Tmp = totEngyInAllMol / totEngyArmAboveMin;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
		<< " , tot engy in all super clusters in MolRad = " << totEngyInAllMol
		<< coutRed << "  -> their ratio = \t " << superClusterMolRatio_Tmp << coutDefault << endl << endl;
#endif

    if(superClusterMolRatio < superClusterMolRatio_Tmp) {
      superClusterIdToCellId		= superClusterIdToCellId_Tmp;
      cellIdToSuperClusterId		= cellIdToSuperClusterId_Tmp;
      superClusterCM			= superClusterCM_Tmp;
      superClusterEngyInMoliere	= superClusterEngyInMoliere_Tmp;

      superClusterMolRatio = superClusterMolRatio_Tmp;

#if _MOL_RAD_CORRECT_DEBUG == 1
      cout	<< coutRed << "\t -- ACCEPTED new superCluster(s) -- " << coutDefault <<  endl << endl;
#endif
      rejectFlag = 0;
    } else {
#if _MOL_RAD_CORRECT_DEBUG == 1
      cout	<< coutRed << "\t -- REJECTED new superCluster(s) -- " << coutDefault <<  endl << endl;
#endif
      rejectFlag = 1;
    }

    // cleanUp
    superClusterIdToCellId_Tmp.clear();
    cellIdToSuperClusterId_Tmp.clear();
    superClusterCM_Tmp.clear();
    superClusterEngyInMoliere_Tmp.clear();

  } else {
#if _MOL_RAD_CORRECT_DEBUG == 1
    if(projectionClusterMolRatio > 0)
      cout	<< coutRed << " \t -- the projrction is not an improvement on the SuperCluster(s) -- "
		<< coutDefault << endl << endl;
#endif
    rejectFlag = 1;
  }


  /* --------------------------------------------------------------------------
     in the case that the projection cluster was not accepted, and some of the
     superClusters have a low energy content in their moliere redius, than
     merge these superClusters with the rest
     -------------------------------------------------------------------------- */
  if(rejectFlag == 1 && (int)superClusterRejected.size() > 0 && (int)superClusterAccepted.size() > 0){

    int numSuperClustersBad  = superClusterRejected.size();
    for(int superClusterNowBad = 0; superClusterNowBad < numSuperClustersBad; superClusterNowBad++){
      int superClusterIdBad = superClusterRejected[superClusterNowBad];

      numElementsInSuperCluster = superClusterIdToCellId[superClusterIdBad].size();
      for(int hitNow = 0; hitNow < numElementsInSuperCluster; hitNow++) {
	cellIdHit = superClusterIdToCellId[superClusterIdBad][hitNow];

	weightedDistanceV.clear();
	const IMPL::CalorimeterHitImpl * thisHit = calHitsCellIdGlobal.at(cellIdHit);

	int numSuperClustersGood = superClusterAccepted.size();
	for(int superClusterNowGood = 0; superClusterNowGood < numSuperClustersGood; superClusterNowGood++) {
	  const int superClusterIdGood = superClusterAccepted[superClusterNowGood];

	  const double distanceCM = distance2D(thisHit->getPosition(),
					       superClusterCM[superClusterIdGood].getPosition());

	  if(distanceCM > 0)
	    weightedDistanceV[superClusterIdGood] = 1./distanceCM;
	  else
	    weightedDistanceV[superClusterIdGood] = 1e10;

	  //					cout << "cellId/hit(x,y): " << cellIdHit << " \t " << CM1[0] << " \t " << CM1[1]
	  //						<< " \t superClusterId/hit(x,y): " << superClusterIdGood << " \t " << CM2[0] << " \t " << CM2[1]
	  //						<< " \t distance: " << distanceCM <<endl;
	}

	// decide to which superCluster to merge the hit according to a proper weight
	double maxSuperClusterWeight = 0, superClusterWeight;
	int   maxWeightSuperClusterId(0);
	weightedDistanceVIterator = weightedDistanceV.begin();
	numClusters               = weightedDistanceV.size();
	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
	  superClusterId     = (int)(*weightedDistanceVIterator).first;
	  superClusterWeight = (double)(*weightedDistanceVIterator).second;

	  if(maxSuperClusterWeight < superClusterWeight) {
	    maxWeightSuperClusterId = superClusterId;
	    maxSuperClusterWeight   = superClusterWeight;
	  }
	}
	//				cout	<< " - num possible -  " << weightedDistanceV.size()
	//					<< " \t id of chosen - " << maxWeightSuperClusterId << endl <<endl;

	// make sure that all hits have have been added to clusters
	assert	(maxSuperClusterWeight > 0);

	// add the hit to the chosen superCluster
	superClusterId = maxWeightSuperClusterId;
	superClusterIdToCellId[superClusterId].push_back(cellIdHit);
	cellIdToSuperClusterId[cellIdHit] = superClusterId;

	//  update the totalEnergy counter and CM position of the superCluster
	updateEngyPosCM(calHitsCellIdGlobal.at(cellIdHit), superClusterCM[superClusterId]);
      }
      // cleanUp
      superClusterIdToCellId.erase(superClusterIdBad);
      superClusterCM.erase(superClusterIdBad);
    }

    // find the percentage of energy for each cluster within _moliereRadius
    totEngyArmAboveMin = 0.;
    totEngyInAllMol    = 0.;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << coutBlue << "Merged superCluster(s):" << coutDefault << endl;
#endif

    superClusterCMIterator = superClusterCM.begin();
    numSuperClusters       = superClusterCM.size();
    for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
      superClusterId = (int)(*superClusterCMIterator).first;

      superClusterEngyInMoliere[superClusterId]
	= getEngyInMoliereFraction(	calHitsCellIdGlobal,
					superClusterIdToCellId[superClusterId],
					superClusterCM[superClusterId],
					1.  );

      totEngyInAllMol       += superClusterEngyInMoliere[superClusterId];
      totEngyArmAboveMin    += superClusterCM[superClusterId].getE();

#if _MOL_RAD_CORRECT_DEBUG == 1
      double engyPercentInMol = superClusterEngyInMoliere[superClusterId]
	/ superClusterCM[superClusterId][0] ;
      cout	<< "superCluster " << superClusterId << " \tat (x,y) = (" << superClusterCM[superClusterId][1]
		<< " , " << superClusterCM[superClusterId][2] << ")   \t engy in _moliereRadius  \t=   "
		<< superClusterEngyInMoliere[superClusterId]
		<< coutRed << " \t-> % totEngy = \t " << engyPercentInMol << coutDefault << endl;
#endif
    }

    superClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout	<< endl << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
		<< " , tot engy in all super clusters in MolRad = " << totEngyInAllMol
		<< coutRed << "  -> their ratio = \t " << superClusterMolRatio << coutDefault << endl << endl;
#endif
  }


  // cleanUp dynamicly allocated memory
  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator
    calHitsCellIdIterator = calHitsCellIdProjection.begin();
  numHitsProjection     = calHitsCellIdProjection.size();
  for(int hitNow = 0; hitNow < numHitsProjection; hitNow++, calHitsCellIdIterator++ ){
    int	cellIdProjection = (int)(*calHitsCellIdIterator).first;
    delete calHitsCellIdProjection[cellIdProjection];
  }
  calHitsCellIdIterator = calHitsCellIdProjectionFull.begin();
  numHitsProjection = calHitsCellIdProjectionFull.size();
  for(int hitNow = 0; hitNow < numHitsProjection; hitNow++, calHitsCellIdIterator++ ){
    int	cellIdProjection = (int)(*calHitsCellIdIterator).first;
    delete calHitsCellIdProjectionFull[cellIdProjection];
  }



  /* --------------------------------------------------------------------------
     re-compute total energy and center of mass of each superCluster (just in case...)
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters               = superClusterIdToCellId.size();
  for(int j=0; j<numSuperClusters; j++, superClusterIdToCellIdIterator++){
    superClusterId  = (int)(*superClusterIdToCellIdIterator).first;  // Id of cluster
    cellIdV = (std::vector<int>)(*superClusterIdToCellIdIterator).second; // cal-hit-Ids in cluster

    // initialize the energy/position std::vector for new clusters only
    //superClusterCM[superClusterId] = LCCluster();
    // superClusterCM[superClusterId].clear();
    // for(int k=0; k<8; k++) superClusterCM[superClusterId].push_back(0.);

    // calculate/update the energy/position of the CM
    calculateEngyPosCM(cellIdV, calHitsCellIdGlobal, superClusterCM[superClusterId], _methodCM);
  }
  cellIdV.clear();



  // write out the results

  return 1;

}
