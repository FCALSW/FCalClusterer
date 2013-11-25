



void LumiCalClustererClass::clusterMerger(	map < int , vector<double> >		* clusterIdToCellEngyP,
						map < int , vector<int> >		* clusterIdToCellIdP,
						map < int , vector<double> >		* clusterCMP,
						map < int , CalorimeterHitImpl* >	calHitsCellIdGlobal ){


  int	detectorArm, clusterId, clusterId1, clusterId2, numClustersToMerge;
  int	numClusters, numClusters2, numElementsInCluster, cellId;
  double	engyNow, zPos;

  map < int , vector<int> >		clusterIdToCellId;
  map < int , vector<int> > :: iterator	clusterIdToCellIdIterator;
  clusterIdToCellId	= * clusterIdToCellIdP;

  map < int , vector<double> >	clusterCM;
  clusterCM =	* clusterCMP;

  map < int , vector<double>  >	clusterIdToCellEngy;
  clusterIdToCellEngy = * clusterIdToCellEngyP;

  vector < int >	allClusterIds,
    cellIdV;

  vector < double >	cellEngyV;

  map < int , double >				idToCellEngy;
  map < int , double > :: iterator	idToCellEngyIterator;

  vector < SuperTrueClusterWeights * >	clusterPairWeightsV;

  SuperTrueClusterWeights	* clusterPairWeightsNow;




  while(1) {
    /* --------------------------------------------------------------------------
       make a vector of all of the cluster Ids, and than write the pair-wise
       properties for all pair combinations
       -------------------------------------------------------------------------- */
    allClusterIds.clear();
    clusterIdToCellIdIterator = clusterIdToCellId.begin();
    numClusters               = clusterIdToCellId.size();
    for(int clusterNow=0; clusterNow<numClusters; clusterNow++, clusterIdToCellIdIterator++){
      clusterId   = (int)(*clusterIdToCellIdIterator).first;  // Id of cluster

      // discount backscattered hits
      zPos = clusterCM[clusterId][3];
      detectorArm = int( zPos/Abs(zPos) );

      allClusterIds.push_back(clusterId);
    }

    numClusters2 =  allClusterIds.size();
    for(int clusterNow1=0; clusterNow1<numClusters2; clusterNow1++) {
      clusterId1 = allClusterIds[clusterNow1];

      for(int clusterNow2=0; clusterNow2<numClusters2; clusterNow2++) {
	clusterId2 = allClusterIds[clusterNow2];

	if(clusterId1 == clusterId2) continue;

	clusterPairWeightsNow = new SuperTrueClusterWeights(	clusterId1,
								clusterId2,
								clusterCM[clusterId1],
								clusterCM[clusterId2] );

	clusterPairWeightsNow->setWeight(	"minEngyDistance",
						_minSeparationDistance,
						engySignalGeV(_minClusterEngyGeV, "GeVToSignal") );

	if(clusterPairWeightsNow->weight < 0) {
	  clusterPairWeightsNow->setWeight("distance");
	  clusterPairWeightsV.push_back(clusterPairWeightsNow);
	}
      }
    }

    // if all the pairs made the cut than finish the loop
    numClustersToMerge = clusterPairWeightsV.size();
    if(numClustersToMerge == 0) break;

    // choose the cluster pair with the shortest weight (distance)
    sort(clusterPairWeightsV.begin(), clusterPairWeightsV.end(), SuperTrueClusterWeightsCompare);

    clusterId1 = clusterPairWeightsV[0]->superClusterId;
    clusterId2 = clusterPairWeightsV[0]->trueClusterId;

    /* --------------------------------------------------------------------------
       go over all hits in the discarded cluster write a new hit-energy map
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

    // initialize the energy/position vector for new clusters only
    clusterCM[clusterId].clear();
    for(int k=0; k<8; k++) clusterCM[clusterId].push_back(0.);

    // calculate/update the energy/position of the CM
    calculateEngyPosCM_EngyV(	cellIdV, cellEngyV, calHitsCellIdGlobal,
				&clusterCM, clusterId, _methodCM);

    // cleanUp
    cellIdV.clear();  cellEngyV.clear();

    /* --------------------------------------------------------------------------
       cleanUp virtual memory allocations
       -------------------------------------------------------------------------- */
    numClustersToMerge = clusterPairWeightsV.size();
    for(int clusterNow1=0; clusterNow1<numClustersToMerge; clusterNow1++)
      delete clusterPairWeightsV[clusterNow1];

    clusterPairWeightsV.clear();
    idToCellEngy.clear();
  }




  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  cout	<< coutRed << "\tClusters:" << coutDefault << endl;

  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters               = clusterIdToCellId.size();
  for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++) {
    clusterId = (int)(*clusterIdToCellIdIterator).first;

    cout	<< "\t Id "  << clusterId
		<< "  \t energy " << clusterCM[clusterId][0]
		<< "     \t pos(x,y) =  ( " << clusterCM[clusterId][1]
		<< " , " << clusterCM[clusterId][2] << " )"
		<< "     \t pos(theta,phi) =  ( " << clusterCM[clusterId][6]
		<< " , " << clusterCM[clusterId][7] << " )"
		<< coutDefault << endl;
  }

  cout	<<  endl;
#endif



  // cleanUp
  allClusterIds.clear();

  * clusterIdToCellIdP	   = clusterIdToCellId;
  * clusterCMP		   = clusterCM;
  * clusterIdToCellEngyP     = clusterIdToCellEngy;
}
