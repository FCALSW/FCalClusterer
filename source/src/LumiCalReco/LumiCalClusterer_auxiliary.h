
/* =========================================================================
   Auxiliary functions
----------------------------------------------------------------------------
  Iftach Sadeh - ???? 2007
============================================================================ */


/* --------------------------------------------------------------------------
sorting of clusterCM[id] (for a cluster with Id 'id') with respect to the cluster CM energy
-------------------------------------------------------------------------- */
//in descending order (highest energy is first)
bool clusterCMEnergyCmpDesc( vector<double> a, vector<double> b ) {
	return a[0] > b[0];
}
//in ascending order (lowest energy is first)
bool clusterCMEnergyCmpAsc( vector<double> a, vector<double> b ) {
	return a[0] < b[0];
}


/* --------------------------------------------------------------------------
sorting of hits with respect to their energies
-------------------------------------------------------------------------- */
//in descending order (highest energy is first)
bool HitEnergyCmpDesc( CalorimeterHitImpl* a, CalorimeterHitImpl* b ) {
	return a->getEnergy() > b->getEnergy();
}

//in ascending order (lowest energy is first)
bool HitEnergyCmpAsc( CalorimeterHitImpl* a, CalorimeterHitImpl* b ) {
	return a->getEnergy() < b->getEnergy();
}


/* --------------------------------------------------------------------------
sorting of hits with respect to their distance from the CM of their cluster
-------------------------------------------------------------------------- */
//in ascending order (shortest distance is first)
bool HitDistanceCMCmpAsc(  vector<double> a, vector<double> b ) {
	return a[3] < b[3];
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
-------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(CalorimeterHitImpl* calHit ,string method) {

	double	posWeightHit = 0.;
	int	detectorArm = int((double)calHit->getPosition()[2] 
					/ Abs((double)calHit->getPosition()[2]));

	if(method == "Energy") posWeightHit = calHit->getEnergy();
	

	// ???????? DECIDE/FIX - improve the log weight constants ????????
	if(method == "Log") {
		posWeightHit = Log(calHit->getEnergy() / _totEngyArm[detectorArm]) + _logWeightConst;
		
		if(posWeightHit < 0) posWeightHit = 0. ;
	}

	return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
-------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeightTureCluster(CalorimeterHitImpl* calHit , double cellEngy, string method) {

	double	posWeightHit = 0.;
	int	detectorArm = int((double)calHit->getPosition()[2] 
					/ Abs((double)calHit->getPosition()[2]));

	if(method == "Energy") posWeightHit = cellEngy;
	

	// ???????? DECIDE/FIX - improve the log weight constants ????????
	if(method == "Log") {
		posWeightHit = Log(cellEngy / _totEngyArm[detectorArm]) + _logWeightConst;
		
		if(posWeightHit < 0) posWeightHit = 0. ;
	}

	return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   - overloaded version with a given energy normalization
-------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(CalorimeterHitImpl* calHit, double totEngy, string method) {

	double	posWeightHit = 0.;
	double	detectorArm = calHit->getPosition()[2]; detectorArm /= Abs(detectorArm);

	if(method == "Energy") posWeightHit = calHit->getEnergy();
	

	// ???????? DECIDE/FIX - improve the log weight constants ????????
	if(method == "Log") {
		posWeightHit = Log(calHit->getEnergy() / totEngy) + _logWeightConst;
		
		if(posWeightHit < 0) posWeightHit = 0. ;
	}

	return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   - overloaded version with a given energy normalization and a logWeightConst
-------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(CalorimeterHitImpl* calHit, double totEngy, string method, double logWeightConstNow) {

	double	posWeightHit = 0.;
	double	detectorArm = calHit->getPosition()[2]; detectorArm /= Abs(detectorArm);

	if(method == "Energy") posWeightHit = calHit->getEnergy();
	

	// ???????? DECIDE/FIX - improve the log weight constants ????????
	if(method == "Log") {
		posWeightHit = Log(calHit->getEnergy() / totEngy) + logWeightConstNow;
		
		if(posWeightHit < 0) posWeightHit = 0. ;
	}

	return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate the distance between two poins in 2D (in cartezian coordinates)
   (first index of arrays is X and the second is Y coordinates {or the other way around})
-------------------------------------------------------------------------- */
double LumiCalClustererClass::distance2D(double *pos1, double *pos2) {

	double	distance = Sqrt( Power(pos1[0]-pos2[0],2) + Power(pos1[1]-pos2[1],2) );

	assert (distance >= 0);

	return distance;
}


/* --------------------------------------------------------------------------
   calculate the distance between two poins in 2D (in polar coordinates)
   (first index of arrays is R and the second is PHI coordinates)
-------------------------------------------------------------------------- */
double LumiCalClustererClass::distance2DPolar(double *pos1, double *pos2) {


        double  distance = Sqrt( pos1[0]*pos1[0] + pos2[0]*pos2[0] - 2 * pos1[0]*pos2[0] * Cos(pos1[1] - pos2[1]) );

	assert (distance >= 0);

	return distance;

}


/* --------------------------------------------------------------------------
   ZPRToId:	return a cellId for a given Z (layer), R (cylinder) and Phi (sector)
   idToZPR:	return Z (layer), R (cylinder) and Phi (sector) for a given cellId
-------------------------------------------------------------------------- */
int LumiCalClustererClass::ZPRToId(int cellZ, int cellPhi, int cellR) {

	int cellId = 0;

	cellId |= ( cellZ   << 0  ) ;
        cellId |= ( cellPhi << 8  ) ;
        cellId |= ( cellR   << 16 ) ;

	return cellId;
}

int LumiCalClustererClass::idToZPR(int cellId, string ZPR) {

	int cellZ, cellPhi, cellR;
	int result = 0;

	// compute Z,Phi,R coordinets acording to the cellId
	cellZ   = (cellId >> 0 ) & 0xff ;
	cellPhi = (cellId >> 8 ) & 0xff ;
	cellR   = (cellId >> 16 ) & 0xff ;

	if(ZPR == "Z") result = cellZ;
	if(ZPR == "R") result = cellR;
	if(ZPR == "P") result = cellPhi;

	return result;
}



/* --------------------------------------------------------------------------
   compute the nearest neigbour cellId
   - in case of the next neighbor being outside the detector, return 0
-------------------------------------------------------------------------- */
int LumiCalClustererClass::getNeighborId(int cellId, int neighborIndex) {

	int cellZ, cellPhi, cellR ;

	// compute Z,Phi,R coordinets acording to the cellId
	cellZ   = idToZPR(cellId,"Z") ;
	cellPhi = idToZPR(cellId,"P") ;
	cellR   = idToZPR(cellId,"R") ;


	// change R cells according to the neighborIndex
	if(neighborIndex  == 0) { 
		cellR += 1;
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 1) { 
		cellR += 2; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 2) {
		cellR -= 1 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 3) {
		cellR -= 2 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 4) {
		cellPhi++ ;
		if(cellPhi == _cellPhiMax) cellPhi = 0;
	}
	if(neighborIndex  == 5) {
		cellPhi-- ;
		if(cellPhi < 0) cellPhi = _cellPhiMax - 1;
	}

	
	if(neighborIndex  == 6) { 
		cellR += 3; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 7) {
		cellR -= 3 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 8) { 
		cellR += 4; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 9) {
		cellR -= 4 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 10) { 
		cellR += 5; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 11) {
		cellR -= 5 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 12) { 
		cellR += 6; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 13) {
		cellR -= 6 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 14) { 
		cellR += 7; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 15) {
		cellR -= 7 ;
		if(cellR   < 0) return 0;
	}
	if(neighborIndex  == 16) { 
		cellR += 8; 
		if(cellR   == _cellRMax) return 0;
	}
	if(neighborIndex  == 17) {
		cellR -= 8 ;
		if(cellR   < 0) return 0;
	}

	// compute the new cellId according to the new Z,Phi,R coordinets
	cellId = ZPRToId(cellZ, cellPhi, cellR) ;

	return cellId ;
}




/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   	(3). calculate the map clusterCM from scratch
-------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM(	vector <int> cellIdV, map < int , CalorimeterHitImpl* > calHitsCellId,
						map < int , vector<double> > * clusterCMp, int clusterId, string method) {

	map < int , vector<double> > clusterCM;
	clusterCM = * clusterCMp;

	int 	numElementsInCluster, cellIdHit;
	double	totEngy, xHit, yHit, thetaHit, phiHit, weightHit, weightSum;
	totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;

	int loopFlag = 1;
	while(loopFlag == 1) {
		numElementsInCluster = cellIdV.size();
		for(int k=0; k<numElementsInCluster; k++) {
			cellIdHit  = cellIdV[k];
				
			weightHit  = posWeight(calHitsCellId[cellIdHit],method);
			weightSum += weightHit;

			xHit      += calHitsCellId[cellIdHit]->getPosition()[0] * weightHit;
			yHit      += calHitsCellId[cellIdHit]->getPosition()[1] * weightHit;

			thetaHit  += thetaPhiCell(cellIdHit, "theta")  * weightHit;
			phiHit    += thetaPhiCell(cellIdHit, "phi")    * weightHit;

			totEngy   += calHitsCellId[cellIdHit]->getEnergy();
		}	
		if(weightSum > 0) {
			xHit     /= weightSum;   yHit   /= weightSum;
			thetaHit /= weightSum;   phiHit /= weightSum;
			loopFlag = 0;
		} else {
			// initialize counters and recalculate with the Energy-weights method
			method = "Energy";
			totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;
		}
	}
	
	clusterCM[clusterId][0] = totEngy;
	clusterCM[clusterId][1] = xHit;
	clusterCM[clusterId][2] = yHit;
	clusterCM[clusterId][3] = calHitsCellId[cellIdHit]->getPosition()[2];
	clusterCM[clusterId][4] = weightSum;
	if(method == "Energy")	clusterCM[clusterId][5] =  1.;
	if(method == "Log")	clusterCM[clusterId][5] = -1.;
	clusterCM[clusterId][6] = thetaHit;
	clusterCM[clusterId][7] = phiHit;

	* clusterCMp = clusterCM;
	clusterCM.clear();
}


/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   	(3). calculate the map clusterCM from scratch
-------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM_EngyV(	vector <int> cellIdV, vector <double> cellEngyV,
							map < int , CalorimeterHitImpl* > calHitsCellId,
							map < int , vector<double> > * clusterCMp, int clusterId, string method) {

	map < int , vector<double> > clusterCM;
	clusterCM = * clusterCMp;

	int 	numElementsInCluster, cellIdHit;
	double	totEngy, xHit, yHit, thetaHit, phiHit, weightHit, weightSum;
	totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;

	int loopFlag = 1;
	while(loopFlag == 1) {
		numElementsInCluster = cellIdV.size();
		for(int k=0; k<numElementsInCluster; k++) {
			cellIdHit  = cellIdV[k];
				
			weightHit  = posWeightTureCluster(calHitsCellId[cellIdHit],cellEngyV[k],method);
			weightSum += weightHit;

			xHit      += calHitsCellId[cellIdHit]->getPosition()[0] * weightHit;
			yHit      += calHitsCellId[cellIdHit]->getPosition()[1] * weightHit;

			thetaHit  += thetaPhiCell(cellIdHit, "theta")  * weightHit;
			phiHit    += thetaPhiCell(cellIdHit, "phi")    * weightHit;

			totEngy   += cellEngyV[k];
		}	
		if(weightSum > 0) {
			xHit     /= weightSum;   yHit   /= weightSum;
			thetaHit /= weightSum;   phiHit /= weightSum;
			loopFlag = 0;
		} else {
			// initialize counters and recalculate with the Energy-weights method
			method = "Energy";
			totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;
		}
	}
	
	clusterCM[clusterId][0] = totEngy;
	clusterCM[clusterId][1] = xHit;
	clusterCM[clusterId][2] = yHit;
	clusterCM[clusterId][3] = calHitsCellId[cellIdHit]->getPosition()[2];
	clusterCM[clusterId][4] = weightSum;
	if(method == "Energy")	clusterCM[clusterId][5] =  1.;
	if(method == "Log")	clusterCM[clusterId][5] = -1.;
	clusterCM[clusterId][6] = thetaHit;
	clusterCM[clusterId][7] = phiHit;

	* clusterCMp = clusterCM;
	clusterCM.clear();
}

/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   	(2). update the map clusterCM with the new cal hit
-------------------------------------------------------------------------- */
void LumiCalClustererClass::updateEngyPosCM(CalorimeterHitImpl* calHit, vector<double> * clusterCMp) {

	vector <double> clusterCM;
	clusterCM = * clusterCMp;
	string method;

	double	engyHit = (double)calHit->getEnergy();
	clusterCM[0] += engyHit;

	if(clusterCM[5] > 0) method = "Energy";
	if(clusterCM[5] < 0) method = "Log";

	double weightHit = posWeight(calHit,method);
	
	if(weightHit > 0){
		double 	xCM     = clusterCM[1];
		double 	yCM     = clusterCM[2];
		double	thetaCM = clusterCM[6];
		double	phiCM   = clusterCM[7];
		double 	weightCM = clusterCM[4];

		xCM *= weightCM; yCM *= weightCM; thetaCM *= weightCM; phiCM *= weightCM;

		int	cellIdHit = (int)  calHit->getCellID0();
		double 	xHit      = (double)calHit->getPosition()[0];
		double 	yHit      = (double)calHit->getPosition()[1];
		double 	thetaHit  = thetaPhiCell(cellIdHit, "theta");
		double 	phiHit    = thetaPhiCell(cellIdHit, "phi");
	
		xHit *= weightHit; yHit *= weightHit; thetaHit *= weightHit;  phiHit *= weightHit;
	
		weightCM += weightHit;
		xCM      += xHit;	xCM 	/= weightCM;
		yCM      += yHit;	yCM	/= weightCM;
		thetaCM  += thetaHit;	thetaCM	/= weightCM;
		phiCM    += phiHit;	phiCM	/= weightCM;

		clusterCM[1] = xCM;
		clusterCM[2] = yCM;
		clusterCM[4] = weightCM;
		clusterCM[6] = thetaCM;
		clusterCM[7] = phiCM;
	}


	* clusterCMp = clusterCM;
	clusterCM.clear();
}

/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   	(3). compute the values without updating clusterCM
-------------------------------------------------------------------------- */
void LumiCalClustererClass::getEngyPosCMValues(	vector <int> cellIdV, map < int , CalorimeterHitImpl* > calHitsCellId,
						double * engyPosCM, string method) {

	int	numElementsInCluster, cellIdHit;
	double	totEngy, xHit, yHit, thetaHit, phiHit, weightHit, weightSum;
	totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;

	int loopFlag = 1;
	while(loopFlag == 1) {
		numElementsInCluster = cellIdV.size();
		for(int k=0; k<numElementsInCluster; k++) {
			cellIdHit  = cellIdV[k];
				
			weightHit  = posWeight(calHitsCellId[cellIdHit],method);
			weightSum += weightHit;

			xHit      += calHitsCellId[cellIdHit]->getPosition()[0] * weightHit;
			yHit      += calHitsCellId[cellIdHit]->getPosition()[1] * weightHit;

			thetaHit  += thetaPhiCell(cellIdHit, "theta")  * weightHit;
			phiHit    += thetaPhiCell(cellIdHit, "phi")    * weightHit;

			totEngy   += calHitsCellId[cellIdHit]->getEnergy();
		}	
		if(weightSum > 0) {
			xHit /= weightSum;   yHit /= weightSum;
			thetaHit /= weightSum;   phiHit /= weightSum;
			loopFlag = 0;
		} else {
			// initialize counters and recalculate with the Energy-weights method
			method = "Energy";
			totEngy = xHit = yHit = thetaHit = phiHit = weightHit = weightSum = 0.;
		}
	}

	engyPosCM[0] = totEngy;
	engyPosCM[1] = xHit;
	engyPosCM[2] = yHit;
	engyPosCM[3] = calHitsCellId[cellIdHit]->getPosition()[2];
	engyPosCM[4] = thetaHit;
	engyPosCM[5] = phiHit;
}


/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
-------------------------------------------------------------------------- */
int LumiCalClustererClass::checkClusterMergeCM(  int clusterId1, int clusterId2, map < int , vector<int> > clusterIdToCellId, 
					    map <int , CalorimeterHitImpl* > calHitsCellId, double distanceAroundCM,
					    double percentOfEngyAroungCM, string method   ){
 	
	// vector for holding the Ids of clusters
	vector <int> cellIdV;
	int	numElementsInCluster, cellIdHit;
	double	engyPosCM[6], CM1[3], CM2[3], distanceCM, engyCM;
	double	totEngyAroundCM = 0.;
			
	// add to cellIdV hits from both clusters, and sum up each cluster's energy
	numElementsInCluster = clusterIdToCellId[clusterId1].size();
	for(int i=0; i<numElementsInCluster; i++) {
		cellIdHit = clusterIdToCellId[clusterId1][i];
		cellIdV.push_back(cellIdHit);
	}
					
	numElementsInCluster = clusterIdToCellId[clusterId2].size();
	for(int i=0; i<numElementsInCluster; i++) {
		cellIdHit = clusterIdToCellId[clusterId2][i];
		cellIdV.push_back(cellIdHit);
	}

	getEngyPosCMValues(cellIdV, calHitsCellId, engyPosCM, method);
	CM1[0] = engyPosCM[1];
	CM1[1] = engyPosCM[2];
	engyCM = engyPosCM[0];

	// check that the CM position has a significant amount of the energy around it
	numElementsInCluster = cellIdV.size();
	for(int i=0; i<numElementsInCluster; i++) {
		cellIdHit  = cellIdV[i];
		CM2[0] = calHitsCellId[cellIdHit]->getPosition()[0];
		CM2[1] = calHitsCellId[cellIdHit]->getPosition()[1];
						
		distanceCM = distance2D(CM1,CM2); 
		if(distanceCM < distanceAroundCM) {
			double engyHit = calHitsCellId[cellIdHit]->getEnergy();
			totEngyAroundCM += engyHit;
		}
	}
	
	if(totEngyAroundCM > (engyCM * percentOfEngyAroungCM))
		return 1;
	else
		return 0;
		
}

/* --------------------------------------------------------------------------
   get the energy around a cluster CM within a disranceToScan raduis
-------------------------------------------------------------------------- */
double LumiCalClustererClass::getEngyInMoliereFraction(	map < int , CalorimeterHitImpl* >  calHitsCellId, 
							vector <int> clusterIdToCellId, 
							vector <double> clusterCM,
							double moliereFraction   ){

	int	cellIdHit, numElementsInCluster;
	double	CM1[3], CM2[3], distanceCM;

	double	disranceToScan = _moliereRadius * moliereFraction;
	double	engyAroundCM = 0., engyHit;

	CM1[0] = clusterCM[1];
	CM1[1] = clusterCM[2];

	cellIdHit = clusterIdToCellId[0];

	numElementsInCluster = clusterIdToCellId.size();
	for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
		cellIdHit = clusterIdToCellId[hitNow];
			
		CM2[0] = calHitsCellId[cellIdHit]->getPosition()[0];
		CM2[1] = calHitsCellId[cellIdHit]->getPosition()[1];

		engyHit = calHitsCellId[cellIdHit]->getEnergy();

		distanceCM = distance2D(CM1,CM2);
		if(distanceCM < disranceToScan)
			engyAroundCM += engyHit;
	}

	return engyAroundCM;

}


// overloaded with different variables and functionality...
double LumiCalClustererClass::getEngyInMoliereFraction(  	map < int , CalorimeterHitImpl* >  calHitsCellId, 
							vector <int> clusterIdToCellId, 
							vector <double> clusterCM,
							double moliereFraction,
							map < int , int >  * flagP  ){

	int	cellIdHit, numElementsInLayer;
	double	engyAroundCM = 0.;
	double	CM1[3], CM2[3], distanceCM;
	map < int , int >  flag = * flagP;
	map < int , CalorimeterHitImpl* > :: iterator calHitsCellIdIterator;

	double	disranceToScan = _moliereRadius * moliereFraction;

	CM1[0] = clusterCM[1];
	CM1[1] = clusterCM[2];

	calHitsCellIdIterator = calHitsCellId.begin();
	numElementsInLayer    = calHitsCellId.size();
	for(int hitNow = 0; hitNow < numElementsInLayer; hitNow++, calHitsCellIdIterator++) {
		cellIdHit = (int)(*calHitsCellIdIterator).first;

		CM2[0] = calHitsCellId[cellIdHit]->getPosition()[0];
		CM2[1] = calHitsCellId[cellIdHit]->getPosition()[1];

		distanceCM = distance2D(CM1,CM2);

		if(distanceCM < disranceToScan && flag[cellIdHit] == 0){
			engyAroundCM += calHitsCellId[cellIdHit]->getEnergy();
			flag[cellIdHit] = 1;
		}
	}

	* flagP = flag;
	return engyAroundCM;
}


/* --------------------------------------------------------------------------
   compute the weighed generated-theta / weighed zPosition of a cluster
-------------------------------------------------------------------------- */
void LumiCalClustererClass::getThetaPhiZCluster(	map < int , CalorimeterHitImpl* >  calHitsCellId,
						vector <int> clusterIdToCellId,
						double totEngy, double * output   ) {
	
	int	numElementsInCluster, cellIdHit;
	double	zCell, zCluster = 0., thetaCell, thetaCluster = 0., phiCell, phiCluster = 0.;
	double	weightHit, weightSum = -1., logWeightConstFactor = 0.;
	string	method = "Log";

	while(weightSum < 0){
		numElementsInCluster = clusterIdToCellId.size();
		for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
			cellIdHit  = clusterIdToCellId[hitNow];
			thetaCell  = thetaPhiCell(cellIdHit, "theta");
			phiCell    = thetaPhiCell(cellIdHit, "phi");
			zCell      = calHitsCellId[cellIdHit]->getPosition()[2];
			weightHit  = posWeight(calHitsCellId[cellIdHit],totEngy,method,(_logWeightConst + logWeightConstFactor));

			if(weightHit > 0){
				if(weightSum < 0) weightSum = 0.;
				weightSum    += weightHit;
				thetaCluster += weightHit * thetaCell;
				phiCluster   += weightHit * phiCell;
				zCluster     += zCell * weightHit;
			}
		}
		if(weightSum < 0) logWeightConstFactor -= .5;
	}
	thetaCluster /= weightSum;  phiCluster /= weightSum;  zCluster /= weightSum;
	
	output[0] = thetaCluster;
	output[1] = phiCluster;
	output[2] = zCluster;

	return ;
}


/* --------------------------------------------------------------------------
   compute the theta/phi of a cell with a given cellId
-------------------------------------------------------------------------- */
double LumiCalClustererClass::thetaPhiCell(int cellId, string output) {

	int	cellIdR, cellIdZ, cellIdPhi;
	double	rCell, zCell, thetaCell, phiCell, outputVal;

	if(output == "theta") {
		cellIdR    = idToZPR(cellId,"R");
		rCell      = _rMin + (cellIdR + .5) * _rCellLength;
		cellIdZ    = idToZPR(cellId,"Z");
		zCell      = Abs(_zFirstLayer) + _zLayerThickness * (cellIdZ - 1);
		thetaCell  = ATan(rCell / zCell);
		outputVal  = thetaCell;
	}

	if(output == "phi") {
		cellIdPhi = idToZPR(cellId,"P");
		phiCell   = TwoPi() * (double(cellIdPhi) + .5) / _cellPhiMax;
		outputVal = phiCell;
	}

	return outputVal;

}



/* --------------------------------------------------------------------------
   get the energy around a cluster CM within a disranceToScan raduis
-------------------------------------------------------------------------- */
double LumiCalClustererClass::getMoliereRadius(	map < int , CalorimeterHitImpl* >  calHitsCellId, 
						vector <int> clusterIdToCellId, 
						vector <double> clusterCM    ){

	vector < vector<double> > clusterHitsEngyPos;
	vector < double > oneHitEngyPos(4);

	int	cellIdHit, numElementsInCluster;
	double	CM1[3], CM2[3], distanceCM;
	double	engyPercentage = .9;

//	double	disranceToScan = _moliereRadius * moliereFraction;
	double	engyAroundCM = 0.;

	CM1[0] = clusterCM[1];
	CM1[1] = clusterCM[2];

	cellIdHit = clusterIdToCellId[0];

	#if _IMPROVE_PROFILE_LAYER_POS == 1
	double thetaZClusterV[3];

	getThetaPhiZCluster(calHitsCellId, clusterIdToCellId, clusterCM[0],thetaZClusterV);
	double	thetaCluster = thetaZClusterV[0];
	double	zCluster     = Abs(thetaZClusterV[2]);

	int	layerMiddle  = int( ( zCluster - Abs(_zFirstLayer) ) / _zLayerThickness ) + 1;
	#endif

	numElementsInCluster = clusterIdToCellId.size();
	for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
		cellIdHit = clusterIdToCellId[hitNow];
			
		CM2[0] = calHitsCellId[cellIdHit]->getPosition()[0];
		CM2[1] = calHitsCellId[cellIdHit]->getPosition()[1];

		#if _IMPROVE_PROFILE_LAYER_POS == 1
		int	layerHit   = idToZPR(cellIdHit,"Z");

		// projection hits have an incoded layer number of (_maxLayerToAnalyse + 1)
		// the original hit's layer number is stored in (the previously unused) CellID1
		if(layerHit == (_maxLayerToAnalyse + 1))
			layerHit = (int)calHitsCellId[cellIdHit]->getCellID1();

		if(layerHit != layerMiddle) {
			double	tngPhiHit = CM2[1] / CM2[0];
			double	deltaZ    = (layerHit - layerMiddle) * _zLayerThickness;
			double	deltaR    = Tan(thetaCluster) * deltaZ;
			int	xSignFix  = int(CM2[0] / Abs(CM2[0])) * int(deltaR / Abs(deltaR));

			double	deltaX  = Sqrt( Abs(deltaR) / (1 + Power(tngPhiHit,2) ) );
				deltaX *= xSignFix;
			double	deltaY  = deltaX * tngPhiHit;

			CM2[0] += deltaX;
			CM2[1] += deltaY;
		}
		#endif

		oneHitEngyPos[0] = calHitsCellId[cellIdHit]->getEnergy();
		oneHitEngyPos[1] = CM2[0];
		oneHitEngyPos[2] = CM2[1];

		oneHitEngyPos[3] = distance2D(CM1,CM2);

		clusterHitsEngyPos.push_back( oneHitEngyPos );
	}

	// sort the vector of the cal hits according to distance from CM in ascending order (shortest distance is first)
	sort (clusterHitsEngyPos.begin(), clusterHitsEngyPos.end(), HitDistanceCMCmpAsc );

	numElementsInCluster = clusterIdToCellId.size();
	for(int i=0; i<numElementsInCluster; i++) {
		if (engyAroundCM < (engyPercentage * clusterCM[0])) {
			engyAroundCM += clusterHitsEngyPos[i][0];
			distanceCM    = clusterHitsEngyPos[i][3];
		} else
			break;
	}

	return distanceCM;
}




/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
-------------------------------------------------------------------------- */
double LumiCalClustererClass::getDistanceAroundCMWithEnergyPercent(	vector <double> clusterCM,
								vector <int> clusterIdToCellId,
								map <int , CalorimeterHitImpl* > calHitsCellId,
								double engyPercentage   ){
	
	vector < vector<double> > clusterHitsEngyPos;
	vector < double > oneHitEngyPos(4);
	double CM1[2], CM2[2], engyAroundCM = 0., distanceCM;
	int numElementsInCluster;
	
	// position of cluster CM
	CM1[0] = clusterCM[1];
	CM1[1] = clusterCM[2];

	engyPercentage *= clusterCM[0];

	// fill a vector with energy, position and distance from CM of every cal hit
	numElementsInCluster = clusterIdToCellId.size();
	for(int i=0; i<numElementsInCluster; i++) {
		int cellIdHit  = clusterIdToCellId[i];

		oneHitEngyPos[0] = calHitsCellId[cellIdHit]->getEnergy();
		oneHitEngyPos[1] = calHitsCellId[cellIdHit]->getPosition()[0];
		oneHitEngyPos[2] = calHitsCellId[cellIdHit]->getPosition()[1];

		CM2[0] = oneHitEngyPos[1];
		CM2[1] = oneHitEngyPos[2];
	
		oneHitEngyPos[3] = distance2D(CM1,CM2); 
		
		clusterHitsEngyPos.push_back( oneHitEngyPos );
	}

	// sort the vector of the cal hits according to distance from CM in ascending order (shortest distance is first)
	sort (clusterHitsEngyPos.begin(), clusterHitsEngyPos.end(), HitDistanceCMCmpAsc );

	numElementsInCluster = clusterIdToCellId.size();
	for(int i=0; i<numElementsInCluster; i++)
		if (engyAroundCM < engyPercentage) {
			engyAroundCM += clusterHitsEngyPos[i][0];
			distanceCM    = clusterHitsEngyPos[i][3];
		}

	
	return distanceCM;

}


// ???????? DECIDE/FIX - make sure these are the correct callibration constants ????????
/* --------------------------------------------------------------------------
   get the parent particle Id of the parent which contributed the most
   energy to the cluster
-------------------------------------------------------------------------- */
double LumiCalClustererClass::engySignalGeV( double engy , TString transformMethod ) {

	if(transformMethod == "SignalToGeV")
		engy = (engy) / .0105;
//		engy = (engy - .0013) / .0105;

	if(transformMethod == "GeVToSignal")
		engy = engy * .0105;
//		engy = engy * .0105 + .0013;

	return engy;
}




/* --------------------------------------------------------------------------
   class and sort rule for computing weights required for assignement of
   reconstructed to true clusters
-------------------------------------------------------------------------- */
class SuperTrueClusterWeights {

public:
		SuperTrueClusterWeights(int superClusterIdNow,
					int trueClusterIdNow,
					vector <double> superClusterCM,
					vector <double> trueClusterCM);

	double	distance2D(double *pos1, double *pos2);
	void	setWeight(TString weightMethod);
	void	setWeight(TString weightMethod, double minSeparationDistance, double minClusterEngyGeV);
		
	int	superClusterId, trueClusterId;
	double	distance, deltaEngy, minEngy, weight;
};

SuperTrueClusterWeights::SuperTrueClusterWeights(	int superClusterIdNow,
							int trueClusterIdNow,
							vector <double> superClusterCM,
							vector <double> trueClusterCM) {
	
	double	pos1[2], pos2[2];
	
	superClusterId = superClusterIdNow;
	trueClusterId  = trueClusterIdNow;

	pos1[0] = superClusterCM[1];	pos1[1] = superClusterCM[2];
	pos2[0] = trueClusterCM[1];	pos2[1] = trueClusterCM[2];

	distance = distance2D(pos1,pos2);

	deltaEngy = Abs(superClusterCM[0] - trueClusterCM[0]);
	
	minEngy = min(superClusterCM[0] , trueClusterCM[0]);
}

void SuperTrueClusterWeights::setWeight(TString weightMethod) {

	if(weightMethod == "distance") 	weight = distance;

	if(weightMethod == "deltaEngy") weight = deltaEngy;
}

void SuperTrueClusterWeights::setWeight(TString weightMethod, double minSeparationDistance, double minClusterEngyGeV) {

	if(weightMethod == "minEngyDistance") {
		if(distance > minSeparationDistance && minEngy > minClusterEngyGeV)
			weight = 1.;
		else
			weight = -1.;
	}
}

double SuperTrueClusterWeights::distance2D(double *pos1, double *pos2) {

	double	distanceNow = Sqrt( Power(pos1[0]-pos2[0],2) + Power(pos1[1]-pos2[1],2) );

	assert (distanceNow >= 0);

	return distanceNow;
}

bool SuperTrueClusterWeightsCompare(SuperTrueClusterWeights * a, SuperTrueClusterWeights * b) {
	return a->weight < b->weight;
}


















