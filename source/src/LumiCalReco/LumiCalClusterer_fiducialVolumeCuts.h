
void LumiCalClustererClass:: fiducialVolumeCuts(	map < int , vector<int> > 		* superClusterIdToCellIdP,
							map < int , vector<double> > 		* superClusterIdToCellEngyP,
							map < int , vector<double> >		* superClusterCMP ) {


	map < int , vector<int> > 	superClusterIdToCellId   = * superClusterIdToCellIdP;
	map < int , vector<double> > 	superClusterIdToCellEngy = * superClusterIdToCellEngyP;
	map < int , vector<double> >	superClusterCM           = * superClusterCMP;

	int	numSuperClusters, superClusterId;
	double	thetaSuperCluster;

	map < int , vector<int> > ::iterator superClusterIdToCellIdIterator;

	vector < int >	clusterIdToErase;

	/* --------------------------------------------------------------------------
	   discard true/reconstructed clusters that are outside the fiducial volume
	-------------------------------------------------------------------------- */
	superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
	numSuperClusters               = superClusterIdToCellId.size();
	for(int superClusterNow=0; superClusterNow<numSuperClusters; superClusterNow++, superClusterIdToCellIdIterator++){
		superClusterId   = (int)(*superClusterIdToCellIdIterator).first;  // Id of cluster

		thetaSuperCluster = Abs(superClusterCM[superClusterId][6]);
		if(thetaSuperCluster < _thetaContainmentBouds[0] || thetaSuperCluster > _thetaContainmentBouds[1])
			clusterIdToErase.push_back(superClusterId);
		
	}

	numSuperClusters = clusterIdToErase.size();
	for(int superClusterNow=0; superClusterNow<numSuperClusters; superClusterNow++) {
		superClusterId = clusterIdToErase[superClusterNow];
		
		#if _GENERAL_CLUSTERER_DEBUG == 1
		cout	<< coutRed << "\tErase cluster " << superClusterId << coutDefault << endl;
		#endif

		superClusterIdToCellId.erase(superClusterId);
		superClusterIdToCellEngy.erase(superClusterId);
		superClusterCM.erase(superClusterId);
	}
	clusterIdToErase.clear();


	* superClusterIdToCellIdP   = superClusterIdToCellId;
	* superClusterIdToCellEngyP = superClusterIdToCellEngy;
	* superClusterCMP           = superClusterCM;

	return ;
}

