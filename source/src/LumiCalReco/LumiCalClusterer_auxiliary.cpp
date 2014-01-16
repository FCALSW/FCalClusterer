//Local
#include "SortingFunctions.hh"
#include "LumiCalClusterer.h"
//LCIO
#include <IMPL/CalorimeterHitImpl.h>
// stdlib
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <iomanip>



/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight( IMPL::CalorimeterHitImpl const* calHit , GlobalMethodsClass::WeightingMethod_t method) {

  double	posWeightHit = 0.;
  int	detectorArm = ((calHit->getPosition()[2] < 0) ? -1 : 1 );

  if(method == GlobalMethodsClass::EnergyMethod) posWeightHit = calHit->getEnergy();


  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if(method == GlobalMethodsClass::LogMethod) {
    posWeightHit = log(calHit->getEnergy() / _totEngyArm[detectorArm]) + _logWeightConst;
    if(posWeightHit < 0) posWeightHit = 0. ;
  }

  return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeightTureCluster(IMPL::CalorimeterHitImpl const* calHit , double cellEngy, GlobalMethodsClass::WeightingMethod_t method) {

  double	posWeightHit = 0.;
  int	detectorArm = ((calHit->getPosition()[2] < 0) ? -1 : 1 );

  if(method == GlobalMethodsClass::EnergyMethod) posWeightHit = cellEngy;


  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if(method == GlobalMethodsClass::LogMethod) {
    posWeightHit = log(cellEngy / _totEngyArm[detectorArm]) + _logWeightConst;

    if(posWeightHit < 0) posWeightHit = 0. ;
  }

  return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   - overloaded version with a given energy normalization
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(IMPL::CalorimeterHitImpl const* calHit, double totEngy,
					GlobalMethodsClass::WeightingMethod_t method) {

  double	posWeightHit = 0.;

  if(method == GlobalMethodsClass::EnergyMethod) posWeightHit = calHit->getEnergy();


  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if(method == GlobalMethodsClass::LogMethod) {
    posWeightHit = log(calHit->getEnergy() / totEngy) + _logWeightConst;

    if(posWeightHit < 0) posWeightHit = 0. ;
  }

  return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   - overloaded version with a given energy normalization and a logWeightConst
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(IMPL::CalorimeterHitImpl const* calHit, double totEngy,
					GlobalMethodsClass::WeightingMethod_t method, double logWeightConstNow) {

  double	posWeightHit = 0.;

  if(method == GlobalMethodsClass::EnergyMethod ) posWeightHit = calHit->getEnergy();


  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if(method == GlobalMethodsClass::LogMethod ) {
    posWeightHit = log(calHit->getEnergy() / totEngy) + logWeightConstNow;

    if(posWeightHit < 0) posWeightHit = 0. ;
  }

  return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate the distance between two poins in 2D (in polar coordinates)
   (first index of arrays is R and the second is PHI coordinates)
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::distance2DPolar(double *pos1, double *pos2) {


  double  distance = sqrt( pos1[0]*pos1[0] + pos2[0]*pos2[0] - 2 * pos1[0]*pos2[0] * cos(pos1[1] - pos2[1]) );

  assert (distance >= 0);

  return distance;

}

/* --------------------------------------------------------------------------
   compute the nearest neigbour cellId
   - in case of the next neighbor being outside the detector, return 0
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::getNeighborId(int cellId, int neighborIndex) {

  int cellZ, cellPhi, cellR ;

  // compute Z,Phi,R coordinets acording to the cellId
  cellZ   = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COZ) ;
  cellPhi = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COP) ;
  cellR   = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COR) ;


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

  // compute the new cellId according to the new Z,Phi,R coordinates
  cellId = GlobalMethodsClass::CellIdZPR(cellZ, cellPhi, cellR) ;

  return cellId ;
}




/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (3). calculate the map clusterCM from scratch
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM(	std::vector <int> const& cellIdV,
						std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
						LCCluster & cluster,
						GlobalMethodsClass::WeightingMethod_t method) {

  double totEngy, xHit, yHit, zHit, thetaHit, weightHit, weightSum;
  totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;
  int loopFlag = 1;
  while(loopFlag == 1) {
    for (std::vector<int>::const_iterator it = cellIdV.begin(); it != cellIdV.end(); ++it) {
      int cellIdHit  = *it;
      const IMPL::CalorimeterHitImpl* calHit = calHitsCellId.at(cellIdHit);
      weightHit  = posWeight(calHit,method);
      weightSum += weightHit;

      const float* position = calHit->getPosition();
      xHit      += position[0] * weightHit;
      yHit      += position[1] * weightHit;
      zHit      += position[2] * weightHit;
      totEngy   += calHit->getEnergy();
      thetaHit  += thetaPhiCell(cellIdHit, GlobalMethodsClass::COTheta)  * weightHit;

    }
    if(weightSum > 0) {
      xHit     /= weightSum;   yHit   /= weightSum; zHit /= weightSum;
      thetaHit /= weightSum;
      loopFlag = 0;
    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = GlobalMethodsClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;
    }
  }

  cluster = LCCluster(totEngy, xHit, yHit, zHit, weightSum, method, thetaHit, 0.0);


}


/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (3). calculate the map clusterCM from scratch
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM_EngyV(	std::vector <int> const& cellIdV, std::vector <double> const& cellEngyV,
							std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
							std::map < int , LCCluster > & clusterCM, int clusterId,
							GlobalMethodsClass::WeightingMethod_t method) {

  int	numElementsInCluster, cellIdHit;
  double	totEngy, xHit, yHit, zHit, thetaHit, weightHit, weightSum;
  totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;

  int loopFlag = 1;
  while(loopFlag == 1) {
    numElementsInCluster = cellIdV.size();
    for(int k=0; k<numElementsInCluster; k++) {
      cellIdHit  = cellIdV[k];
      const IMPL::CalorimeterHitImpl * thisHit = calHitsCellId.at(cellIdHit);
      weightHit  = posWeightTureCluster(thisHit,cellEngyV[k],method);
      weightSum += weightHit;

      xHit      += thisHit->getPosition()[0] * weightHit;
      yHit      += thisHit->getPosition()[1] * weightHit;
      zHit      += thisHit->getPosition()[2] * weightHit;

      thetaHit  += thetaPhiCell(cellIdHit, GlobalMethodsClass::COTheta)  * weightHit;

      totEngy   += cellEngyV[k];
    }
    if(weightSum > 0) {
      xHit     /= weightSum;   yHit   /= weightSum; zHit /= weightSum;
      thetaHit /= weightSum;
      loopFlag = 0;
    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = GlobalMethodsClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;
    }
  }

  clusterCM[clusterId] = LCCluster(totEngy, xHit, yHit, zHit,
				   weightSum, method, thetaHit, 0.0);

}

/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (2). update the map clusterCM with the new cal hit
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::updateEngyPosCM(IMPL::CalorimeterHitImpl* calHit, LCCluster & clusterCM) {

  GlobalMethodsClass::WeightingMethod_t method;

  double	engyHit = (double)calHit->getEnergy();
  clusterCM.addToEnergy(engyHit);

  method = clusterCM.getMethod();

  double weightHit = posWeight(calHit,method);

  if(weightHit > 0){
    double xCM      = clusterCM.getX();
    double yCM      = clusterCM.getY();
    double thetaCM  = clusterCM.getTheta();
    double phiCM    = clusterCM.getPhi();
    double weightCM = clusterCM.getWeight();

    xCM *= weightCM; yCM *= weightCM; thetaCM *= weightCM; phiCM *= weightCM;

    int	cellIdHit = (int)  calHit->getCellID0();
    double xHit      = (double)calHit->getPosition()[0];
    double yHit      = (double)calHit->getPosition()[1];
    double thetaHit  = thetaPhiCell(cellIdHit, GlobalMethodsClass::COTheta);
    double phiHit    = thetaPhiCell(cellIdHit, GlobalMethodsClass::COPhi);

    xHit *= weightHit; yHit *= weightHit; thetaHit *= weightHit;  phiHit *= weightHit;

    weightCM += weightHit;
    xCM      += xHit;	xCM	/= weightCM;
    yCM      += yHit;	yCM	/= weightCM;
    thetaCM  += thetaHit;	thetaCM	/= weightCM;
    phiCM    += phiHit;	phiCM	/= weightCM;

    clusterCM.setX(xCM);
    clusterCM.setY(yCM);
    clusterCM.setWeight(weightCM);
    clusterCM.setTheta(thetaCM);
    clusterCM.setPhi(phiCM);
  }

}

/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (3). compute the values without updating clusterCM
   -------------------------------------------------------------------------- */
LCCluster LumiCalClustererClass::getEngyPosCMValues(	std::vector <int> const& cellIdV,
							std::map < int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
							GlobalMethodsClass::WeightingMethod_t method) {

  int	numElementsInCluster, cellIdHit;
  double	totEngy, xHit, yHit, zHit, thetaHit, weightHit, weightSum;
  totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;

  int loopFlag = 1;
  while(loopFlag == 1) {
    numElementsInCluster = cellIdV.size();
    for(int k=0; k<numElementsInCluster; k++) {
      cellIdHit  = cellIdV[k];
      IMPL::CalorimeterHitImpl * thisHit = calHitsCellId.at(cellIdHit);
      weightHit  = posWeight(thisHit,method);
      weightSum += weightHit;

      xHit      += thisHit->getPosition()[0] * weightHit;
      yHit      += thisHit->getPosition()[1] * weightHit;
      zHit      += thisHit->getPosition()[2] * weightHit;
      thetaHit  += thetaPhiCell(cellIdHit, GlobalMethodsClass::COTheta)  * weightHit;
      totEngy   += thisHit->getEnergy();
    }
    if(weightSum > 0) {
      xHit /= weightSum; yHit /= weightSum; zHit /= weightSum;
      thetaHit /= weightSum;
      loopFlag = 0;
    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = GlobalMethodsClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightHit = weightSum = 0.;
    }
  }

  return LCCluster(totEngy, xHit, yHit, zHit,
		   weightSum, method, thetaHit, 0.0);

}


/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::checkClusterMergeCM(  int clusterId1, int clusterId2,
						 std::map < int , std::vector<int> > const& clusterIdToCellId,
						 std::map <int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
						 double distanceAroundCM, double percentOfEngyAroungCM,
						 GlobalMethodsClass::WeightingMethod_t method ){

  // std::vector for holding the Ids of clusters
  std::vector <int> cellIdV;
  int	numElementsInCluster, cellIdHit;
  double	CM1[3], CM2[3], distanceCM, engyCM;
  double	totEngyAroundCM = 0.;

  // add to cellIdV hits from both clusters, and sum up each cluster's energy
  numElementsInCluster = clusterIdToCellId.at(clusterId1).size();
  for(int i=0; i<numElementsInCluster; i++) {
    cellIdHit = clusterIdToCellId.at(clusterId1).at(i);
    cellIdV.push_back(cellIdHit);
  }

  numElementsInCluster = clusterIdToCellId.at(clusterId2).size();
  for(int i=0; i<numElementsInCluster; i++) {
    cellIdHit = clusterIdToCellId.at(clusterId2).at(i);
    cellIdV.push_back(cellIdHit);
  }

  LCCluster engyPosCM( getEngyPosCMValues(cellIdV, calHitsCellId, method) );
  CM1[0] = engyPosCM.getX();
  CM1[1] = engyPosCM.getY();
  engyCM = engyPosCM.getE();

  // check that the CM position has a significant amount of the energy around it
  numElementsInCluster = cellIdV.size();
  for(int i=0; i<numElementsInCluster; i++) {
    cellIdHit  = cellIdV[i];
    CM2[0] = calHitsCellId.at(cellIdHit)->getPosition()[0];
    CM2[1] = calHitsCellId.at(cellIdHit)->getPosition()[1];

    distanceCM = distance2D(CM1,CM2);
    if(distanceCM < distanceAroundCM) {
      double engyHit = calHitsCellId.at(cellIdHit)->getEnergy();
      totEngyAroundCM += engyHit;
    }
  }

  if(totEngyAroundCM > (engyCM * percentOfEngyAroungCM))
    return 1;
  else
    return 0;

}

/* --------------------------------------------------------------------------
   get the energy around a cluster CM within a distanceToScan raduis
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::getEngyInMoliereFraction(	std::map < int , IMPL::CalorimeterHitImpl* >  const& calHitsCellId,
							std::vector <int> const&, //clusterIdToCellId,
							LCCluster const& clusterCM,
							double moliereFraction   ){

  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator calHitsCellIdIterator, calEnd;

  const double distanceToScan = _moliereRadius * moliereFraction;
  double	engyAroundCM = 0.0;

  calHitsCellIdIterator = calHitsCellId.begin();
  calEnd    = calHitsCellId.end();

  for(; calHitsCellIdIterator != calEnd; ++calHitsCellIdIterator) {
    const IMPL::CalorimeterHitImpl * calHit = calHitsCellIdIterator->second;
    const double distanceCM = distance2D(clusterCM.getPosition(),calHit->getPosition());
    if(distanceCM < distanceToScan)
      engyAroundCM += calHit->getEnergy();
  }

  return engyAroundCM;

}


// overloaded with different variables and functionality...
double LumiCalClustererClass::getEngyInMoliereFraction( std::map < int , IMPL::CalorimeterHitImpl* >  const& calHitsCellId,
							std::vector <int> const&,//clusterIdToCellId,
							LCCluster const& clusterCM,
							double moliereFraction,
							std::map < int , int >  & flag  ){

  double	engyAroundCM = 0.;
  std::map < int , IMPL::CalorimeterHitImpl* > :: const_iterator calHitsCellIdIterator, calEnd;

  double	distanceToScan = _moliereRadius * moliereFraction;

  calHitsCellIdIterator = calHitsCellId.begin();
  calEnd    = calHitsCellId.end();
  for(; calHitsCellIdIterator != calEnd; ++calHitsCellIdIterator) {
    const int cellIdHit = calHitsCellIdIterator->first;
    const IMPL::CalorimeterHitImpl * calHit = calHitsCellIdIterator->second;

    const double distanceCM = distance2D(clusterCM.getPosition(),calHit->getPosition());

    if(distanceCM < distanceToScan && flag[cellIdHit] == 0){
      engyAroundCM += calHit->getEnergy();
      flag[cellIdHit] = 1;
    }
  }

  return engyAroundCM;
}


/* --------------------------------------------------------------------------
   compute the weighed generated-theta / weighed zPosition of a cluster
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::getThetaPhiZCluster( std::map < int , IMPL::CalorimeterHitImpl* >  const& calHitsCellId,
						 std::vector <int> const& clusterIdToCellId,
						 double totEngy, double * output   ) {

  int	numElementsInCluster, cellIdHit;
  double	zCell, zCluster = 0., thetaCell, thetaCluster = 0., phiCell, phiCluster = 0.;
  double	weightHit, weightSum = -1., logWeightConstFactor = 0.;
  GlobalMethodsClass::WeightingMethod_t method = GlobalMethodsClass::LogMethod;

  while(weightSum < 0){
    numElementsInCluster = clusterIdToCellId.size();
    for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
      cellIdHit  = clusterIdToCellId[hitNow];
      thetaCell  = thetaPhiCell(cellIdHit, GlobalMethodsClass::COTheta);
      phiCell    = thetaPhiCell(cellIdHit, GlobalMethodsClass::COPhi);
      zCell      = calHitsCellId.at(cellIdHit)->getPosition()[2];
      weightHit  = posWeight(calHitsCellId.at(cellIdHit),totEngy,method,(_logWeightConst + logWeightConstFactor));

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
double LumiCalClustererClass::thetaPhiCell(int cellId, GlobalMethodsClass::Coordinate_t output) {

  int	cellIdR, cellIdZ, cellIdPhi;
  double	rCell, zCell, thetaCell, phiCell, outputVal(-1);

  if(output == GlobalMethodsClass::COTheta) {
    cellIdR    = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COR);
    rCell      = _rMin + (cellIdR + .5) * _rCellLength;
    cellIdZ    = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COZ);
    zCell      = fabs(_zFirstLayer) + _zLayerThickness * (cellIdZ - 1);
    thetaCell  = atan(rCell / zCell);
    outputVal  = thetaCell;
  }
  else if(output == GlobalMethodsClass::COPhi) {
    cellIdPhi = GlobalMethodsClass::CellIdZPR(cellId,GlobalMethodsClass::COP);
    phiCell   = 2*M_PI * (double(cellIdPhi) + .5) / _cellPhiMax;
    outputVal = phiCell;
  }

  return outputVal;

}



/* --------------------------------------------------------------------------
   get the energy around a cluster CM within a disranceToScan raduis
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::getMoliereRadius(	std::map < int , IMPL::CalorimeterHitImpl* >  const& calHitsCellId,
						std::vector <int> const& clusterIdToCellId,
						LCCluster const& clusterCM    ){

  std::vector < double > oneHitEngyPos(2);

  const int numElementsInCluster = clusterIdToCellId.size();
  std::vector < std::vector<double> > clusterHitsEngyPos(numElementsInCluster, oneHitEngyPos);

  double   distanceCM(0.0);
  double	engyPercentage = .9;

  //	double	disranceToScan = _moliereRadius * moliereFraction;
  double	engyAroundCM = 0.;


#if _IMPROVE_PROFILE_LAYER_POS == 1
  double thetaZClusterV[3];
  getThetaPhiZCluster(calHitsCellId, clusterIdToCellId, clusterCM[0],thetaZClusterV);
  double	thetaCluster = thetaZClusterV[0];
  double	zCluster     = Abs(thetaZClusterV[2]);

  int	layerMiddle  = int( ( zCluster - Abs(_zFirstLayer) ) / _zLayerThickness ) + 1;
#endif

  for(int hitNow = 0; hitNow < numElementsInCluster; hitNow++) {
    const int cellIdHit = clusterIdToCellId[hitNow];
    const IMPL::CalorimeterHitImpl * thisHit = calHitsCellId.at(cellIdHit);
    double CM2[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

#if _IMPROVE_PROFILE_LAYER_POS == 1
    int	layerHit   = GlobalMethodsClass::CellIdZPR(cellIdHit,GlobalMethodsClass::COZ);

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

    clusterHitsEngyPos[hitNow][0] = thisHit->getEnergy();
    clusterHitsEngyPos[hitNow][1] = distance2D(clusterCM.getPosition(),CM2);
  }

  // sort the std::vector of the cal hits according to distance from CM in ascending order (shortest distance is first)
  std::sort (clusterHitsEngyPos.begin(), clusterHitsEngyPos.end(), HitDistanceCMCmpAsc<1> );

  for(int i=0; i<numElementsInCluster; i++) {
    if (engyAroundCM < (engyPercentage * clusterCM.getE())) {
      engyAroundCM += clusterHitsEngyPos[i][0];
      distanceCM    = clusterHitsEngyPos[i][1];
    } else
      break;
  }

  return distanceCM;
}




/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::getDistanceAroundCMWithEnergyPercent( LCCluster const& clusterCM,
								    std::vector <int> const& clusterIdToCellId,
								    std::map <int , IMPL::CalorimeterHitImpl* > const& calHitsCellId,
								    double engyPercentage   ){
  engyPercentage *= clusterCM.getE();
  double engyAroundCM = 0., distanceCM(0.0);

  const int numElementsInCluster(clusterIdToCellId.size());

  std::vector < double > oneHitEngyPos(2);
  std::vector < std::vector<double> > clusterHitsEngyPos (numElementsInCluster, oneHitEngyPos);

  // fill a std::vector with energy, position and distance from CM of every cal hit
  for(int i=0; i<numElementsInCluster; i++) {
    int cellIdHit  = clusterIdToCellId[i];

    IMPL::CalorimeterHitImpl *thisHit = calHitsCellId.at(cellIdHit);

    clusterHitsEngyPos[i][0] = thisHit->getEnergy();

    clusterHitsEngyPos[i][1] = distance2D(clusterCM.getPosition(),thisHit->getPosition());

  }

  // sort the std::vector of the cal hits according to distance from CM in ascending order (shortest distance is first)
  sort (clusterHitsEngyPos.begin(), clusterHitsEngyPos.end(), HitDistanceCMCmpAsc<1> );

  for(int i=0; i<numElementsInCluster; i++)
    if (engyAroundCM < engyPercentage) {
      engyAroundCM += clusterHitsEngyPos[i][0];
      distanceCM    = clusterHitsEngyPos[i][1];
    }


  return distanceCM;

}
