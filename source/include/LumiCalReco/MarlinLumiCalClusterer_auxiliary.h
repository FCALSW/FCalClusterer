namespace MarlinLumiCalClusterer {


/* --------------------------------------------------------------------------
-------------------------------------------------------------------------- */
void MarlinLumiCalClusterer::TryMarlinLumiCalClusterer(  LCEvent * evt  ){

	try{
  
		TString hisName, optionName;
		int	numMCParticles, particleId;
		double	engyNow, totEngyNow, thetaNow, phiNow;

		map < int , map < int , ClusterClass * > >	clusterClassMap;
		map < int , ClusterClass * > :: iterator	clusterClassMapIterator;

		map < int , map < int , vector<int> > > 	clusterIdToCellId;
		map < int , map < int , vector<double> > >	cellIdToCellEngy;

		int	numClusters, clusterId;
		double	weightNow;
	
		vector < SortingClass * >	sortingClassV;
		SortingClass			* sortingClass;


		/* --------------------------------------------------------------------------
		   create clusters using: LumiCalClustererClass
		-------------------------------------------------------------------------- */
		// instantiate a clusterClass object for each mcParticle which was created
		// infront of LumiCal and was destroyed after lumical.
		LumiCalClusterer -> processEvent(evt);
		
		clusterIdToCellId = LumiCalClusterer -> _superClusterIdToCellId;
		cellIdToCellEngy  = LumiCalClusterer -> _superClusterIdToCellEngy;

		CreateClusters( clusterIdToCellId , cellIdToCellEngy , &clusterClassMap );

		
		/* --------------------------------------------------------------------------
		   flag the highest-energy cluster in each arm
		-------------------------------------------------------------------------- */
		for(int armNow = -1; armNow < 2; armNow += 2) {
			sortingClassV.clear();
			
			numClusters             = clusterClassMap[armNow].size();
			clusterClassMapIterator = clusterClassMap[armNow].begin();
			for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
				clusterId = (int)(*clusterClassMapIterator).first;

				if(clusterClassMap[armNow][clusterId]->OutsideFlag == 1)	continue;

				weightNow   = clusterClassMap[armNow][clusterId]->Engy;
				weightNow   = Power(weightNow , -1);

				// instantiate the sorting class with the weight and insert to the vector
				sortingClass = new SortingClass(clusterId , weightNow);
				sortingClassV.push_back(sortingClass);	
			}

			// arrange the vector according to inverse of the energy (lowest first)
			sort(sortingClassV.begin(),sortingClassV.end(),cmpRuleDesc);

			// flag the highest-energy cluster
			numClusters = sortingClassV.size();
			for(int clusterNow = 0; clusterNow < numClusters; clusterNow++) {
				clusterId = sortingClassV[clusterNow]->Id;

				if(clusterNow == 0)
					clusterClassMap[armNow][clusterId]->HighestEnergyFlag = 1;
				else
					clusterClassMap[armNow][clusterId]->HighestEnergyFlag = 0;
			}


			// cleanUp
			numClusters = sortingClassV.size();
			for(int clusterNow = 0; clusterNow < numClusters; clusterNow++) {
				delete sortingClassV[clusterNow];
			}	

		}



		/* --------------------------------------------------------------------------
		   histograming
		-------------------------------------------------------------------------- */		
		for(int armNow = -1; armNow < 2; armNow += 2) {	
			totEngyNow = 0.;

			clusterClassMapIterator = clusterClassMap[armNow].begin();
			numMCParticles          = clusterClassMap[armNow].size();
			for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
				particleId = (int)(*clusterClassMapIterator).first;

				// only take into account particles inside LumiCal
				if(clusterClassMap[armNow][particleId]->OutsideFlag == 1)	continue;

				engyNow = clusterClassMap[armNow][particleId] -> Engy;
				engyNow = GlobalMethods->SignalGevConversion("Signal_to_GeV" , engyNow);

				totEngyNow += engyNow;
			}
			// total energy inside LumiCal
			if(totEngyNow > 0) {
				hisName = "totEnergyIn";	OutputManager->HisMap1D[hisName] -> Fill (totEngyNow);
			}
		}



		for(int armNow = -1; armNow < 2; armNow += 2) {	
			totEngyNow = 0.;
			
			clusterClassMapIterator = clusterClassMap[armNow].begin();
			numMCParticles          = clusterClassMap[armNow].size();
			for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
				particleId = (int)(*clusterClassMapIterator).first;

				// only take into account particles outside LumiCal
				if(clusterClassMap[armNow][particleId]->OutsideFlag == 0)	continue;
				// only count in particles from the current arm
				if(clusterClassMap[armNow][particleId]->SignMC != armNow)	continue;

				// energy/theta dependance of particles outside LumiCal the fiducial volume
				if(clusterClassMap[armNow][particleId]->OutsideReason == "Reconstructed outside the fiducial volume") {
					engyNow  = clusterClassMap[armNow][particleId] -> Engy;
					engyNow  = GlobalMethods->SignalGevConversion("Signal_to_GeV" , engyNow);
					thetaNow = clusterClassMap[armNow][particleId] -> Theta;

					hisName = "thetaEnergyOut_DepositedEngy";
					OutputManager->HisMap2D[hisName] -> Fill (engyNow , thetaNow);

				} else {
					engyNow  = clusterClassMap[armNow][particleId] -> EngyMC;
					thetaNow = clusterClassMap[armNow][particleId] -> ThetaMC;
				}

				totEngyNow += engyNow;
			}

			// total energy outside LumiCal
			if(totEngyNow > 0) {
				hisName = "totEnergyOut";	OutputManager->HisMap1D[hisName] -> Fill (totEngyNow);
			}
		}


		for(int armNow = -1; armNow < 2; armNow += 2) {	
			clusterClassMapIterator = clusterClassMap[armNow].begin();
			numMCParticles          = clusterClassMap[armNow].size();
			for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
				particleId = (int)(*clusterClassMapIterator).first;

				// only take into account the highest-energy particle in the arm
				if(clusterClassMap[armNow][particleId]->HighestEnergyFlag == 0)	continue;

				engyNow  = clusterClassMap[armNow][particleId] -> Engy;
				engyNow  = GlobalMethods->SignalGevConversion("Signal_to_GeV" , engyNow);
				thetaNow = clusterClassMap[armNow][particleId] -> Theta;

				hisName = "higestEngyParticle_Engy";	OutputManager->HisMap1D[hisName] -> Fill (engyNow);
				hisName = "higestEngyParticle_Theta";	OutputManager->HisMap1D[hisName] -> Fill (thetaNow);
			}
		}


		/* --------------------------------------------------------------------------
		   write criteria for selection cuts for Bhabha events into a tree
		-------------------------------------------------------------------------- */
		int	clusterInFlag = 0;
		for(int armNow = -1; armNow < 2; armNow += 2) {	
			clusterClassMapIterator = clusterClassMap[armNow].begin();
			numMCParticles          = clusterClassMap[armNow].size();
			for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
				particleId = (int)(*clusterClassMapIterator).first;

				// only take into account the highest-energy particle in the arm
				// (by default this means that this particle's shower is contained)
				if(clusterClassMap[armNow][particleId]->HighestEnergyFlag == 0)	continue;
				clusterInFlag++;

				engyNow  = clusterClassMap[armNow][particleId] -> Engy;
				engyNow  = GlobalMethods->SignalGevConversion("Signal_to_GeV" , engyNow);
				thetaNow = clusterClassMap[armNow][particleId] -> Theta;
				phiNow   = clusterClassMap[armNow][particleId] -> Phi;


				OutputManager -> TreeIntV["nEvt"]	= EvtNumber;
				OutputManager -> TreeIntV["sign"]	= armNow;
				OutputManager -> TreeDoubleV["engy"]	= engyNow;
				OutputManager -> TreeDoubleV["theta"]	= thetaNow;
				OutputManager -> TreeDoubleV["phi"]	= phiNow;

				OutputManager -> TreeMap["bhabhaSelectionTree"] -> Fill();
			}
		}

		// fil in a flag entry (all values = -1) in the tree in the case that no arm has any clusters
		if(clusterInFlag == 0) {
		
			OutputManager -> TreeIntV["nEvt"]	= EvtNumber;
			OutputManager -> TreeDoubleV["engy"]	= -1;
			OutputManager -> TreeDoubleV["theta"]	= -1;
			OutputManager -> TreeDoubleV["phi"]	= -1;

			OutputManager -> TreeIntV["sign"]	= 1;
			OutputManager -> TreeMap["bhabhaSelectionTree"] -> Fill();

			OutputManager -> TreeIntV["sign"]	= -1;
			OutputManager -> TreeMap["bhabhaSelectionTree"] -> Fill();
		}


		
		#if _GLOBAL_COUNTERS_UPDATE_DEBUG == 1
		// write out the counter map
		int numCounters = OutputManager->Counter.size();
	
		if(numCounters > 0)
			cout 	<< endl << coutRed << "Global counters:"  << coutDefault << endl;

		OutputManager->CounterIterator = OutputManager->Counter.begin();
		for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager->CounterIterator++) {
			TString counterName = (TString)(*OutputManager->CounterIterator).first;
				cout << "\t" << OutputManager->Counter[counterName] << "  \t <->  " << counterName << endl;
		}
		#endif
		


		/* --------------------------------------------------------------------------
		   cleanUp
		-------------------------------------------------------------------------- */
		for(int armNow = -1; armNow < 2; armNow += 2) {	
			clusterClassMapIterator = clusterClassMap[armNow].begin();
			numMCParticles          = clusterClassMap[armNow].size();
			for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
				particleId = (int)(*clusterClassMapIterator).first;

				delete	clusterClassMap[armNow][particleId];
			}
		}


		/* --------------------------------------------------------------------------
		   write to the root tree
		-------------------------------------------------------------------------- */
		OutputManager->WriteToRootTree("" , NumEvt);

	} // try

	// if an !E!9exception has been thrown (no *col for this event) than do....
	catch( DataNotAvailableException &e ){
		#ifdef _LC_DEBUG
		cout << "Event " << NumEvt << " has an exception"<< endl;
		#endif
	}

	return;
}




void MarlinLumiCalClusterer::CreateClusters(	map < int , map < int , vector<int> > > 	clusterIdToCellId,
						map < int , map < int , vector<double> > > 	cellIdToCellEngy,
						map < int , map < int , ClusterClass * > >	* clusterClassMapP ) {


	map < int , map < int , ClusterClass * > >		clusterClassMap = * clusterClassMapP;
	map < int , ClusterClass * > :: iterator		clusterClassMapIterator;

	map < int , vector<int> >  ::iterator	clusterIdToCellIdIterator;

	int	numClusters, clusterId, numElementsInCluster, cellId;
	double	engyHit;

	for(int armNow = -1; armNow < 2; armNow += 2) {	
		
		clusterIdToCellIdIterator = clusterIdToCellId[armNow].begin();
		numClusters               = clusterIdToCellId[armNow].size();
		for(int superClusterNow = 0; superClusterNow < numClusters; superClusterNow++, clusterIdToCellIdIterator++) {
			clusterId = (int)(*clusterIdToCellIdIterator).first;

			clusterClassMap[armNow][clusterId] = new ClusterClass(clusterId);
			clusterClassMap[armNow][clusterId] -> SetStatsMC();
			clusterClassMap[armNow][clusterId] -> SignMC = armNow;

			numElementsInCluster = clusterIdToCellId[armNow][clusterId].size();
			for(int cellNow = 0; cellNow < numElementsInCluster; cellNow++) {

				cellId  = clusterIdToCellId[armNow][clusterId][cellNow];
				engyHit = cellIdToCellEngy [armNow][clusterId][cellNow];

				clusterClassMap[armNow][clusterId] -> FillHit(cellId , engyHit);
			}
		}
	}



	#if _CLUSTER_RESET_STATS_DEBUG == 1
	cout	<< endl << coutUnderLine << coutPurple << "Transfering information into ClusterClass objects ......  "
		<< coutDefault << endl;
	#endif

	/* --------------------------------------------------------------------------
	   calculate the energy and position of each cluster
	-------------------------------------------------------------------------- */
	for(int armNow = -1; armNow < 2; armNow += 2) {
		#if _CLUSTER_RESET_STATS_DEBUG == 1
		cout << endl << coutGreen << "Initial Reset Stats for LumiCal arm = " << armNow << "  ...... " << coutDefault << endl;
		#endif

		clusterClassMapIterator = clusterClassMap[armNow].begin();
		numClusters             = clusterClassMap[armNow].size();
		for (int MCParticleNow = 0; MCParticleNow < numClusters; MCParticleNow++, clusterClassMapIterator++) {
			clusterId = (int)(*clusterClassMapIterator).first;

			int	resetStatsFlag = clusterClassMap[armNow][clusterId] -> ResetStats();

			#if _CLUSTER_RESET_STATS_DEBUG == 1
			if(resetStatsFlag == 1) continue;		
			if(clusterClassMap[armNow][clusterId] -> SignMC != armNow) continue;		

			cout	<< coutRed << "\tParticle Out ("
				<< clusterClassMap[armNow][clusterId] -> OutsideReason << "):   " << clusterId << endl
				<< "\t\t pdg, parentId , NumMCDaughters = "
				<< "\t" << clusterClassMap[armNow][clusterId] -> Pdg
				<< "\t" << clusterClassMap[armNow][clusterId] -> ParentId
				<< "\t" << clusterClassMap[armNow][clusterId] -> NumMCDaughters << endl
				<< "\t\t vertex -> endPoint X, Y, Z = "
				<< "\t" << clusterClassMap[armNow][clusterId] -> VtxX
				<< "\t" << clusterClassMap[armNow][clusterId] -> VtxY
				<< "\t" << clusterClassMap[armNow][clusterId] -> VtxZ << "   ->   "
				<< "\t" << clusterClassMap[armNow][clusterId] -> EndPointX
				<< "\t" << clusterClassMap[armNow][clusterId] -> EndPointY
				<< "\t" << clusterClassMap[armNow][clusterId] -> EndPointZ << endl
				<< "\t\t engy, theta, phi (mc):\t "
				<< clusterClassMap[armNow][clusterId] -> EngyMC  << "\t"
				<< clusterClassMap[armNow][clusterId] -> ThetaMC << "\t"
				<< clusterClassMap[armNow][clusterId] -> PhiMC
				<< coutDefault << endl << endl;
			#endif
		}
	}


	* clusterClassMapP = clusterClassMap;

	return;

}




}  // namespace MarlinLumiCalClusterer


