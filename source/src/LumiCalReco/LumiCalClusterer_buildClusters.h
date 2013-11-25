	
/* =========================================================================
   LumiCalClustererClass :: buildClusters
============================================================================
  (1). Description:
--------------------------------
   - SOME DESCRIPTION ......
============================================================================ */

int LumiCalClustererClass::buildClusters(	map < int , vector <CalorimeterHitImpl*> > 	& calHits,
				 	 	map < int , CalorimeterHitImpl* > 		* calHitsCellIdGlobalP,
					 	map < int , vector<int> > 			* superClusterIdToCellIdP,
				 	 	map < int , vector<double> > 			* superClusterIdToCellEngyP,
				 	 	map < int , vector<double> > 			* superClusterCMP ) {

	
	int	detectorArm, numClusters;
	int	numElementsInArm , minNumElementsInShowerPeakLayer, numElementsInLayer;
	int	numClustersMax, numClustersMajority, numClustersAllLayers;
	int	maxEngyLayerN, engyInMoliereFlag;
	double	maxEngyLayer;
	double	middleEnergyHitBound;
	double	CM1[2], CM2[2], distanceCM;

	int	numSuperClusters, superClusterId, numElementsInSuperCluster, cellIdHit;
	double	engyHit;

	TString	hisName;
	TF1	* fitFunc;

	vector < map <int , CalorimeterHitImpl* > > 	calHitsCellId(_maxLayerToAnalyse),
	       						calHitsSmallEngyCellId(_maxLayerToAnalyse);

	map < int , CalorimeterHitImpl* > 		calHitsCellIdGlobal,
	    						calHitsCellIdProjection,
							calHitsCellIdProjectionFull;
	map < int , CalorimeterHitImpl* > :: iterator 	calHitsCellIdIterator;
	
	vector < map < int , int > > 	cellIdToClusterId(_maxLayerToAnalyse+1);
	
	vector < map < int , vector<int> > > 	clusterIdToCellId(_maxLayerToAnalyse+1);
	map < int , vector<int> > :: iterator 	clusterIdToCellIdIterator;

	vector < map < int , vector<double> > > 	clusterCM(_maxLayerToAnalyse+1);
	map < int , vector<double> > :: iterator	clusterCMIterator;

	map < int , int > 		cellIdToSuperClusterId;
	map < int , int > :: iterator 	cellIdToSuperClusterIdIterator;
	
	map < int , vector<int> > 	superClusterIdToCellId;
	map < int , vector<double> > 	superClusterIdToCellEngy;
	
	map < int , vector<double> > 			superClusterCM;
	map < int , vector<double> > :: iterator	superClusterCMIterator;

	vector < int > 	initialClusterControlVar(4),
	       		isShowerPeakLayer(_maxLayerToAnalyse);

	map < int , vector < vector<double> > >			engyPosCMLayer;
	map < int , vector < vector<double> > >	:: iterator	engyPosCMLayerIterator;

	map < int , vector <int> >	thisLayer;

	vector < vector <double> >	avrgCM;
	
	map < int , int > 		numClustersCounter;
	map < int , int > :: iterator 	numClustersCounterIterator;

	map <int , double>		weightedDistanceV;
	map <int , double> :: iterator	weightedDistanceVIterator;

	vector < map < int , vector<double> > >	virtualClusterCM(_maxLayerToAnalyse);
	vector < double >			virtualClusterCMV(3);

	map < int , TH1F * > 		xLineFitCM, yLineFitCM;
	map < int , vector<double> > 	fitParamX, fitParamY;
	
	map < int , double > 			layerToPosX, layerToPosY, layerToEngy;
	map < int , double > :: iterator	layerToPosXYIterator;


	/* --------------------------------------------------------------------------
	   determine the detector arm
	-------------------------------------------------------------------------- */
	detectorArm = 0;
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++){
		if( (int)calHits[layerNow].size() > 0) {
			detectorArm  = int( (double)calHits[layerNow][0]->getPosition()[2] );
			detectorArm /= Abs( detectorArm );
			break;
		}
	}


	/* --------------------------------------------------------------------------
	   determine the total energy of the hits in the arm
	-------------------------------------------------------------------------- */
	#if _GENERAL_CLUSTERER_DEBUG == 1
	string detectorArmName;
	if(detectorArm > 0) detectorArmName = "positive detector arm";
	if(detectorArm < 0) detectorArmName = "negative detector arm";
	cout 	<< coutRed << "\tTotal arm energy =  " <<  coutDefault
		<< _totEngyArm[detectorArm] << endl << endl;
	#endif

	
	/* --------------------------------------------------------------------------
	   calculate the minimal energy to take into account in the initial
	   clustering pass
	-------------------------------------------------------------------------- */
	middleEnergyHitBound = Exp(-1*_logWeightConst) * _totEngyArm[detectorArm] * _middleEnergyHitBoundFrac;

	
	/* --------------------------------------------------------------------------
	   determine the min number of hits that makeup a showerPeak layer and
	   flag the layers that make the cut
	-------------------------------------------------------------------------- */
	numElementsInArm = 0;
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++){
		numElementsInArm += (int)calHits[layerNow].size();
	}
	minNumElementsInShowerPeakLayer = int(numElementsInArm * _elementsPercentInShowerPeakLayer);

	#if _CLUSTER_BUILD_DEBUG == 1
	cout	<< coutBlue << "Shower peak layers:" << coutDefault << endl;
	cout 	<< "\t min # of hits, min energy/hit([signal],[GeV]) : " << minNumElementsInShowerPeakLayer
		<< "\t(" <<  middleEnergyHitBound << " , " << engySignalGeV(middleEnergyHitBound, "SignalToGeV")
		<< ")" << endl << "\t layers chosen : ";
	#endif

	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++){
		numElementsInLayer = (int)calHits[layerNow].size();
		if(numElementsInLayer >= minNumElementsInShowerPeakLayer) {
			isShowerPeakLayer[layerNow] = 1;
			#if _CLUSTER_BUILD_DEBUG == 1
			cout << layerNow << " , " ;
			#endif

		}
		else
			isShowerPeakLayer[layerNow] = 0;
	}
	#if _CLUSTER_BUILD_DEBUG == 1
	cout << endl;
	#endif


	/* --------------------------------------------------------------------------
	   fill calHitsCellId with the layer taged cal hits. for showerPeak layers
	   separate cal hits with energy above/below the middleEnergyHitBound.
	   in any case only choose hits with energy above the _hitMinEnergy cut
	-------------------------------------------------------------------------- */
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) {
		numElementsInLayer = (int)calHits[layerNow].size();
		for(int j=0; j<numElementsInLayer; j++){
			int	cellIdHit = (int)calHits[layerNow][j]->getCellID0();
			double	cellEngy = (double)calHits[layerNow][j]->getEnergy();

			if(cellEngy < _hitMinEnergy) continue;

			if(cellEngy > middleEnergyHitBound)
				calHitsCellId[layerNow][cellIdHit] = calHits[layerNow][j];

			#if _CLUSTER_MIDDLE_RANGE_ENGY_HITS == 1
			if(isShowerPeakLayer[layerNow] == 1) {
				if(cellEngy <= middleEnergyHitBound)
					calHitsSmallEngyCellId[layerNow][cellIdHit] = calHits[layerNow][j];
			} else
				calHitsCellId[layerNow][cellIdHit] = calHits[layerNow][j];
			#endif
		}
	}
	

	/* --------------------------------------------------------------------------
	   form initial clusters for the shower-peak layers
	-------------------------------------------------------------------------- */	
	#if _CLUSTER_BUILD_DEBUG == 1
	cout	<< endl << coutBlue
		<< "run initialClusterBuild and initialLowEngyClusterBuild:" << coutDefault << endl;
	#endif

	// set the control vector for the initialClusterBuild clustering options
	initialClusterControlVar[0] = 1;  // mergeOneHitClusters
	initialClusterControlVar[1] = 1;  // mergeSmallToLargeClusters
	initialClusterControlVar[2] = 1;  // mergeLargeToSmallClusters
	initialClusterControlVar[3] = 1;  // forceMergeSmallToLargeClusters

	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) if(isShowerPeakLayer[layerNow] == 1) {
		// run the initial clustering algorithm for the high energy hits
		initialClusterBuild( calHitsCellId[layerNow],
				     &(cellIdToClusterId[layerNow]),
				     &(clusterIdToCellId[layerNow]),
				     &(clusterCM[layerNow]),
				     initialClusterControlVar );

		#if _CLUSTER_MIDDLE_RANGE_ENGY_HITS == 1
		// cluster the low energy hits
		initialLowEngyClusterBuild( calHitsSmallEngyCellId[layerNow],
					    &(calHitsCellId[layerNow]),
					    &(cellIdToClusterId[layerNow]),
					    &(clusterIdToCellId[layerNow]),
					    &(clusterCM[layerNow]) );
		#endif

		#if _CLUSTER_BUILD_DEBUG == 1
		cout << "\t layer " << layerNow << endl;
		clusterCMIterator = clusterCM[layerNow].begin();
		numClusters       = clusterCM[layerNow].size();
		for(int clusterNow1 = 0; clusterNow1 < numClusters; clusterNow1++, clusterCMIterator++){
			int clusterId = (int)(*clusterCMIterator).first;
			cout	<< "\t\t cluster Id, pos(x,y), engy: "
				<< clusterId << "\t (" << clusterCM[layerNow][clusterId][1] << " , "
				<< clusterCM[layerNow][clusterId][2] << ") \t " << clusterCM[layerNow][clusterId][0] << endl;
		}
		cout << endl;
		#endif
	}


	/* --------------------------------------------------------------------------
	   decide how many global clusters there are
	-------------------------------------------------------------------------- */
	// find the number of clusters in the majority of layers
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) if(isShowerPeakLayer[layerNow] == 1) {
		numClusters = clusterCM[layerNow].size();
		numClustersCounter[numClusters]++;
	}

	numClustersMax = 0;
	numClustersCounterIterator = numClustersCounter.begin();
	numClustersAllLayers       = numClustersCounter.size();
	for(int i=0; i<numClustersAllLayers; i++,numClustersCounterIterator++){
		int numClustersNow  = (int)(*numClustersCounterIterator).first;
		int numClustersSize = (int)(*numClustersCounterIterator).second;

		if(numClustersMax < numClustersSize) {
			numClustersMax      = numClustersSize;
			numClustersMajority = numClustersNow;
		}
	}
	numClustersCounter.clear();

	#if _CLUSTER_BUILD_DEBUG == 1
	cout << coutRed << "\t -> Assume thet there are " << numClustersMajority
		<< " global clusters" << coutDefault << endl << endl;
	#endif


	/* --------------------------------------------------------------------------
	   choose only layers which have a numClustersMajority number of clusters,
	   and add their CM to engyPosCMLayer[clusterNow] the decision to choose a
	   certain clusterNow is made according to the distance of the given
	   cluster CM from an averaged CM
	-------------------------------------------------------------------------- */
	// find the layer with the most energy which has numClustersMajority clusters
	maxEngyLayer = 0.;
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) if(isShowerPeakLayer[layerNow] == 1) {
		clusterCMIterator = clusterCM[layerNow].begin();
		numClusters       = clusterCM[layerNow].size();

		double	engyLayerNow = 0.;
		if(numClusters != numClustersMajority) continue;
		for(int clusterNow1 = 0; clusterNow1 < numClustersMajority; clusterNow1++, clusterCMIterator++){
			int clusterId = (int)(*clusterCMIterator).first;

			engyLayerNow += clusterCM[layerNow][clusterId][0];
		}

		if(maxEngyLayer < engyLayerNow) {
			maxEngyLayer  = engyLayerNow;
			maxEngyLayerN = layerNow;
		}
	}

	// for the layer with the most energy which has numClustersMajority clusters,
	// initialize the averageCM vector
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) if(isShowerPeakLayer[layerNow] == 1) {
		clusterCMIterator = clusterCM[layerNow].begin();
		numClusters       = clusterCM[layerNow].size();
		
		if(numClusters != numClustersMajority)  continue;
		if(layerNow != maxEngyLayerN)		continue;

		for(int clusterNow = 0; clusterNow < numClustersMajority; clusterNow++, clusterCMIterator++){
			int clusterId = (int)(*clusterCMIterator).first;
			
			// store the CM energy/position vector for each cluster
			engyPosCMLayer[clusterNow].push_back( clusterCM[layerNow][clusterId] );
			thisLayer[clusterNow].push_back( layerNow );

			// initialize the averaged CM position vector
			avrgCM.push_back( clusterCM[layerNow][clusterId] );
		}
		break;
	}
	
	// for all layers but the layer with the most energy which has numClustersMajority clusters,
	// update the averageCM vector
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++) if(isShowerPeakLayer[layerNow] == 1) {
		clusterCMIterator = clusterCM[layerNow].begin();
		numClusters       = clusterCM[layerNow].size();
		
		if(numClusters != numClustersMajority)	continue;
		if(layerNow == maxEngyLayerN)		continue;

		for(int clusterNow1 = 0; clusterNow1 < numClustersMajority; clusterNow1++, clusterCMIterator++){
			int clusterId = (int)(*clusterCMIterator).first;
			
			CM1[0] = clusterCM[layerNow][clusterId][1];
			CM1[1] = clusterCM[layerNow][clusterId][2];

			// compare the position of the CM to the averaged CM positions
			for(int clusterNow2 = 0; clusterNow2 < numClustersMajority; clusterNow2++){
				CM2[0] = avrgCM[clusterNow2][1];
				CM2[1] = avrgCM[clusterNow2][2];
				
				distanceCM = distance2D(CM1,CM2);
				if(distanceCM > 0)
					weightedDistanceV[clusterNow2] = Power(distanceCM,-1);
				else
					weightedDistanceV[clusterNow2] = 1e10;
			}

			double	maxClusterWeight = -1;
			int	maxWeightClusterId;
			weightedDistanceVIterator = weightedDistanceV.begin();
			numClusters               = weightedDistanceV.size();
			for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, weightedDistanceVIterator++){
				int	clusterIdNow  = (int)(*weightedDistanceVIterator).first;
				double	clusterWeight = (double)(*weightedDistanceVIterator).second;

				if(maxClusterWeight < clusterWeight) {
					maxWeightClusterId = clusterIdNow;
					maxClusterWeight   = clusterWeight;
				}
			}

			// add the CM to the right vector
			engyPosCMLayer[maxWeightClusterId].push_back(clusterCM[layerNow][clusterId]);			
			thisLayer[maxWeightClusterId].push_back( layerNow );
		
			// update the multi-layer CM position
			avrgCM[maxWeightClusterId][0] += clusterCM[layerNow][clusterId][0];
			avrgCM[maxWeightClusterId][1]  = (CM1[0]+CM2[0])/2.;
			avrgCM[maxWeightClusterId][2]  = (CM1[1]+CM2[1])/2.;
		}	
	}


	/* --------------------------------------------------------------------------
	   fit a stright line through each cluster from engyPosCMLayer. results
	   (line parametrizations) are stored in fitParamX(Y)
	-------------------------------------------------------------------------- */
	#if _CLUSTER_BUILD_DEBUG == 1
	cout 	<< coutBlue << "Fit lines through the averaged CM" << coutDefault << endl << endl;
	#endif

	fitFunc = new TF1("fitFunc","[0]+[1]*x",-3000,-2000);

	engyPosCMLayerIterator = engyPosCMLayer.begin();
	numClusters            = engyPosCMLayer.size();
	for(int clusterNow=0; clusterNow<numClusters; clusterNow++, engyPosCMLayerIterator++){
		int	clusterId = (int)(*engyPosCMLayerIterator).first;

		int	numLayers = engyPosCMLayer[clusterId].size();
		if(numLayers < 3) {
			#if _CLUSTER_BUILD_DEBUG == 1
			cout	<< coutRed << "\tdecrease the global cluster number by 1"
				<< coutDefault << endl << endl;
			#endif
			numClustersMajority--;
			continue;
		}

		#if _CLUSTER_BUILD_DEBUG == 1
		cout << coutRed << "clusterId " << clusterId << coutDefault << endl;
		#endif

		hisName = "_xLineFitCM_Cluster"; hisName += clusterNow;
		xLineFitCM[clusterNow] = new TH1F(  hisName,hisName,_maxLayerToAnalyse*10,0,_maxLayerToAnalyse);
		
		hisName = "_yLineFitCM_Cluster"; hisName += clusterNow;
		yLineFitCM[clusterNow] = new TH1F(  hisName,hisName,_maxLayerToAnalyse*10,0,_maxLayerToAnalyse);


		/* --------------------------------------------------------------------------
		   since more than on cluster may have choosen the same averagedCM in a
		   given layer, some layers may have more than one entry in the engyPosCMLayer
		   map. therefore an averaging is performed for each layer
		-------------------------------------------------------------------------- */
		// initialize the layer-averaging position/energy maps
		numLayers = engyPosCMLayer[clusterId].size();
		for(int layerN = 0; layerN < numLayers; layerN++) {
			int	layerNow = thisLayer[clusterId][layerN];
			layerToPosX[layerNow] = layerToPosY[layerNow] = layerToEngy[layerNow] = 0.;
		}

		// fill the maps with energy-weighted positions
		for(int layerN = 0; layerN < numLayers; layerN++) {
			int	layerNow = thisLayer[clusterId][layerN];

			layerToPosX[layerNow] += engyPosCMLayer[clusterNow][layerN][1] * engyPosCMLayer[clusterNow][layerN][0];
			layerToPosY[layerNow] += engyPosCMLayer[clusterNow][layerN][2] * engyPosCMLayer[clusterNow][layerN][0];
			layerToEngy[layerNow] += engyPosCMLayer[clusterNow][layerN][0];
		}

		// fill historgams of x(z) and y(z) of the CM positions
		layerToPosXYIterator = layerToPosX.begin();
		numLayers            = layerToPosX.size();
		for(int layerN = 0; layerN < numLayers; layerN++, layerToPosXYIterator++) {
			int	layerNow = (int)(*layerToPosXYIterator).first;

			// get back to units of position
			layerToPosX[layerNow] /= layerToEngy[layerNow];
			layerToPosY[layerNow] /= layerToEngy[layerNow];

			xLineFitCM[clusterNow] -> Fill(layerNow , layerToPosX[layerNow]);
			yLineFitCM[clusterNow] -> Fill(layerNow , layerToPosY[layerNow]);

			#if _CLUSTER_BUILD_DEBUG == 1
			cout 	<< "\tlayer , avPos(x,y) : " << layerNow << " \t (" << layerToPosX[layerNow]
				<< " , " << layerToPosY[layerNow] << ")" << endl;
			#endif
		}

		// fit a stight line for each histogram, and store the fit results
		xLineFitCM[clusterNow]->Fit("fitFunc","+CQ0"); 
		fitParamX[clusterNow].push_back( fitFunc->GetParameter(0) );
		fitParamX[clusterNow].push_back( fitFunc->GetParameter(1) );

		#if _CLUSTER_BUILD_DEBUG == 1
		cout 	<< "\t -> xFitPar 0,1:  " 
			<< fitFunc->GetParameter(0) << " (+-) " << fitFunc->GetParError(0)
		 	<< " \t,\t " << fitFunc->GetParameter(1) << " (+-) " << fitFunc->GetParError(1) << endl;
		#endif

		yLineFitCM[clusterNow] -> Fit("fitFunc","+CQ0");
		fitParamY[clusterNow].push_back( fitFunc->GetParameter(0) );
		fitParamY[clusterNow].push_back( fitFunc->GetParameter(1) );

		#if _CLUSTER_BUILD_DEBUG == 1
		cout 	<< "\t -> yFitPar 0,1:  " 
			<< fitFunc->GetParameter(0) << " (+-) " << fitFunc->GetParError(0)
		 	<< " \t,\t " << fitFunc->GetParameter(1) << " (+-) " << fitFunc->GetParError(1) << endl << endl;
		#endif

		// cleanUp
		delete xLineFitCM[clusterNow]; delete yLineFitCM[clusterNow];
		layerToPosX.clear();  layerToPosY.clear();  layerToEngy.clear();
	}
	// cleanUp
	xLineFitCM.clear(); yLineFitCM.clear(); delete fitFunc;


	/* --------------------------------------------------------------------------
	   extrapolate the CM positions in all layers (and form clusters in the
	   non shower-peak layers)
	-------------------------------------------------------------------------- */	
	#if _CLUSTER_BUILD_DEBUG == 1
	cout 	<< coutBlue << "Extrapolate virtual cluster CMs" << coutDefault << endl << endl;
	#endif

	// fill virtual cluster CM vectors for all the layers
	for(int layerNow = 0; layerNow < _maxLayerToAnalyse; layerNow++ ){

		numElementsInLayer = calHitsCellId[layerNow].size();
		if(numElementsInLayer == 0) continue;

		for(int clusterNow=0; clusterNow<numClustersMajority; clusterNow++){
			int	maxLayerToRaiseVirtualClusterSize = int(0.75*_maxLayerToAnalyse);
			double	fitPar0, fitPar1, hitLayerRatio;

			// extrapolated x/y positions
			virtualClusterCMV[1] =  fitParamX[clusterNow][0] + fitParamX[clusterNow][1] * layerNow; // x position
			virtualClusterCMV[2] =  fitParamY[clusterNow][0] + fitParamY[clusterNow][1] * layerNow; // y position

			// ???????? DECIDE/FIX - incorparate the parameters given here better in the code ????????
			// ???????? DECIDE/FIX - consider a different middle layer for the else condition ????????
			// extrapolated cluster radius around CM position
			if(avrgCM[clusterNow][0] > 1)	{ fitPar0 = 236.7; fitPar1 = 9.11; hitLayerRatio = 22/2618.; }
			else				{ fitPar0 = 226.5; fitPar1 = 10.3; hitLayerRatio = 22/2570.; }
			
			if(layerNow < maxLayerToRaiseVirtualClusterSize)
				virtualClusterCMV[0] = Exp( (fitPar0 + fitPar1 * layerNow) * hitLayerRatio );
			else
				virtualClusterCMV[0] = Exp( (fitPar0 + fitPar1 * maxLayerToRaiseVirtualClusterSize) * hitLayerRatio );
			
			// The numbers above for fitPar0/1 were derived for a detector with moliereRadius=18.2
			// they must, therefore, they must be corrected for according to the _moliereRadius used now
			virtualClusterCM[layerNow][clusterNow] = virtualClusterCMV;
		}

		// form clusters for the non shower-peak layers in the non shower-peak layers only.
		if(isShowerPeakLayer[layerNow] == 0) { 
			virtualCMClusterBuild( calHitsCellId[layerNow],
					       &(cellIdToClusterId[layerNow]),
					       &(clusterIdToCellId[layerNow]),
					       &(clusterCM[layerNow]),
					       virtualClusterCM[layerNow] );
			#if _CLUSTER_BUILD_DEBUG == 1
			cout	<< coutRed << "\tbuild cluster around a virtual CM in layer " << layerNow 
				<< coutDefault << endl;
			#endif
		}

		// only do something if there are less real clusters than virtual clusters in the layer
		else {
			int numRealClusters    = clusterIdToCellId[layerNow].size();
			int numVirtualClusters = virtualClusterCM[layerNow].size();
			if(numRealClusters < numVirtualClusters) {
				#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
				cout	<< endl << coutRed << "\tin layer " << layerNow << " there are real/virtual clusters: "
					<< numRealClusters << "  ,  " << numVirtualClusters <<  coutDefault << endl;
				#endif

				virtualCMPeakLayersFix(	calHitsCellId[layerNow],
							&(cellIdToClusterId[layerNow]),
							&(clusterIdToCellId[layerNow]),
							&(clusterCM[layerNow]),
							virtualClusterCM[layerNow] );

				#if _CLUSTER_BUILD_DEBUG == 1
				cout	<< coutRed << "\tre-cluster in layer " << layerNow <<  coutDefault << endl;
				#endif
			}
		}
	}
	// cleanUp
	avrgCM.clear(); fitParamX.clear(); fitParamY.clear();


	/* --------------------------------------------------------------------------
	   merge all the clusters from the different layers into superClusters.
	   the number of final clusters will be numClustersMajority (the majority
	   found in the shower-peak layers), and each cluster will be merged into
	   a superCluster according to its distance from it (stored as virtual clusters
	   in virtualClusterCM).
	-------------------------------------------------------------------------- */
	#if _CLUSTER_BUILD_DEBUG == 1
	cout 	<< endl << coutBlue << "Build superClusters" << coutDefault << endl << endl;
	#endif

	int	buildSuperClustersFlag = buildSuperClusters(	&(calHitsCellIdGlobal),
								calHitsCellId,
								clusterIdToCellId,
								clusterCM,
								virtualClusterCM,
								&(cellIdToSuperClusterId),
								&(superClusterIdToCellId),
								&(superClusterCM) );

	if(buildSuperClustersFlag == 0)
		return 0;

	// cleanUp
	virtualClusterCM.clear();




	#if _MOLIERE_RADIUS_CORRECTIONS == 1
	#if _CLUSTER_BUILD_DEBUG == 1 || _MOL_RAD_CORRECT_DEBUG == 1
	cout 	<< endl << coutPurple << "RUN engyInMoliereCorrections() ..." << coutDefault << endl << endl;
	#endif

	engyInMoliereFlag = engyInMoliereCorrections( calHitsCellIdGlobal,
						      calHits,
						      calHitsCellId,
						      &(clusterIdToCellId),
						      &(clusterCM),
		    				      &(cellIdToClusterId),
						      &(cellIdToSuperClusterId),
		    				      &(superClusterIdToCellId),
		    				      &(superClusterCM),
		    				      middleEnergyHitBound,
		    				      detectorArm );

	if(engyInMoliereFlag == 0) 
		return 0;
	
	#endif  // #if _MOLIERE_RADIUS_CORRECTIONS == 1






	/* --------------------------------------------------------------------------
	   --------------------------------------------------------------------------
	   NOTE:
	   --------------------------------------------------------------------------
	     FROM THIS POINT ON ENERGY OF HITS BELONGING TO A SUPERCLUSTER
	     MUST BE ACCESED BY:	superClusterIdToCellEngy
	     AND NOT BY:	   	calHitsCellIdGlobal[cellIdHit]->getEnergy()
	     AS THE ENERGY OF SINGLE HITS MAY BE DIVIDED BY SEVERAL CLUSTERS !
	   --------------------------------------------------------------------------
	-------------------------------------------------------------------------- */

	
	/* --------------------------------------------------------------------------
	   fill the superClusterIdToCellEngy container
	-------------------------------------------------------------------------- */
	superClusterCMIterator = superClusterCM.begin();
	numSuperClusters       = superClusterCM.size();
	for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
		superClusterId = (int)(*superClusterCMIterator).first;

		numElementsInSuperCluster = superClusterIdToCellId[superClusterId].size();
		for(int cellNow = 0; cellNow < numElementsInSuperCluster; cellNow++) {
			cellIdHit = superClusterIdToCellId[superClusterId][cellNow];

			engyHit = calHitsCellIdGlobal[cellIdHit] -> getEnergy();

			superClusterIdToCellEngy[superClusterId].push_back(engyHit);
		}
	}

	/* --------------------------------------------------------------------------
	   verbosity
	-------------------------------------------------------------------------- */
	#if _GENERAL_CLUSTERER_DEBUG == 1
	cout	<< coutRed << "\tClusters:" << coutDefault << endl;

	superClusterCMIterator = superClusterCM.begin();
	numSuperClusters       = superClusterCM.size();
	for(int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
		superClusterId = (int)(*superClusterCMIterator).first;

		cout	<< "\t Id "  << superClusterId
			<< "  \t energy " << superClusterCM[superClusterId][0] 
			<< "     \t pos(x,y) =  ( " << superClusterCM[superClusterId][1]
			<< " , " << superClusterCM[superClusterId][2] << " )"
			<< "     \t pos(theta,phi) =  ( " << superClusterCM[superClusterId][6]
			<< " , " << superClusterCM[superClusterId][7] << " )"
			<< coutDefault << endl;
	}

	cout	<<  endl;
	#endif

	* superClusterIdToCellIdP   = superClusterIdToCellId;
	* superClusterIdToCellEngyP = superClusterIdToCellEngy;
	* superClusterCMP           = superClusterCM;
	* calHitsCellIdGlobalP      = calHitsCellIdGlobal;

	return 1;
}


