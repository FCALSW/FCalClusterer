//Local
#include "SortingFunctions.hh"
#include "LumiCalClusterer.h"
#include "Distance2D.hh"
using LCHelper::distance2D;

// stdlib
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <iomanip>


#ifdef FCAL_WITH_DD4HEP
#include <DD4hep/DD4hepUnits.h>
#endif



/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   -------------------------------------------------------------------------- */
template <> double LumiCalClustererClass::posWeight(CalHit const& calHit, GlobalMethodsClass::WeightingMethod_t method) {
  double posWeightHit = -1.;

  if( method == GlobalMethodsClass::EnergyMethod ) {
 
      return ( calHit->getEnergy() ) ;

  }else if (method == GlobalMethodsClass::LogMethod)  {            // ???????? DECIDE/FIX - improve the log weight constants ????????

    int	detectorArm = ((calHit->getPosition()[2] < 0) ? -1 : 1 );
    posWeightHit = log(calHit->getEnergy() / _totEngyArm[detectorArm]) + _logWeightConst;
    if(posWeightHit < 0) posWeightHit = 0. ;
    return posWeightHit;

  }
  return posWeightHit;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   -------------------------------------------------------------------------- */
template <>
double LumiCalClustererClass::posWeightTrueCluster(CalHit const& calHit, double cellEngy,
                                                   GlobalMethodsClass::WeightingMethod_t method) {
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
template <>
double LumiCalClustererClass::posWeight(CalHit const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method) {
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
template <>
double LumiCalClustererClass::posWeight(CalHit const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method,
                                        double logWeightConstNow) {
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
   - get neighbors in Phi +/- 1 (neighborIndex < 2), and neighbors in R +/- n ( n = neighborIndex/2 )
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::getNeighborId(int cellId, int neighborIndex) {

  int cellZ, cellPhi, cellR, arm;
  // compute Z,Phi,R coordinates according to the cellId
  GlobalMethodsClass::CellIdZPR( cellId, cellZ, cellPhi, cellR, arm);

  // change iRho cell index  according to the neighborIndex
  if( neighborIndex >= 2 ){
    if ( neighborIndex%2 ) cellR -= neighborIndex/2;
    else cellR += neighborIndex/2;
  }

  // change iPhi cell index
  if ( neighborIndex == 0 )     { cellPhi += 1; }
  else if( neighborIndex == 1 ) { cellPhi -= 1; }

  if ( cellR == _cellRMax || cellR < 0 ) return 0; 

  if( cellPhi == _cellPhiMax ) cellPhi = 0;
  else if ( cellPhi < 0 ) cellPhi = _cellPhiMax-1 ;
  // compute neighbor cellId according to the new Z,Phi,R coordinates
  cellId = GlobalMethodsClass::CellIdZPR(cellZ, cellPhi, cellR, arm) ;

  return cellId ;
}




/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (3). calculate the map clusterCM from scratch
   -------------------------------------------------------------------------- */
LCCluster LumiCalClustererClass::calculateEngyPosCM( VInt const& cellIdV,
                                                     MapIntCalHit const& calHitsCellId,
                                                     GlobalMethodsClass::WeightingMethod_t method) {

  double totEngy(0.0), xHit(0.0), yHit(0.0), zHit(0.0), thetaHit(0.0), weightSum(0.0);
  int loopFlag = 1;
  VecCalHit caloHits;
  while(loopFlag == 1) {
    caloHits.clear();
    for (VInt::const_iterator it = cellIdV.begin(); it != cellIdV.end(); ++it) {
      auto& calHit = calHitsCellId.at(*it);
      caloHits.push_back(calHit);

      const double weightHit = posWeight(calHit,method);
      weightSum += weightHit;

      const double* position = calHit->getPosition();

      xHit      += position[0] * weightHit;
      yHit      += position[1] * weightHit;
      zHit      += position[2] * weightHit;
      totEngy   += calHit->getEnergy();

    }
    if(weightSum > 0.) {
      xHit     /= weightSum;   yHit   /= weightSum; zHit /= weightSum;
      thetaHit  = atan( sqrt( xHit*xHit + yHit*yHit)/fabs( zHit ));
      loopFlag = 0;

    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = GlobalMethodsClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightSum = 0.;
    }
  }

  return LCCluster(totEngy, xHit, yHit, zHit, weightSum, method, thetaHit, 0.0, caloHits);
}


/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (3). calculate the map clusterCM from scratch
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM_EngyV( VInt const& cellIdV,
                                                      VDouble const& cellEngyV,
                                                      MapIntCalHit const& calHitsCellId,
                                                      MapIntLCCluster & clusterCM, int clusterId,
                                                      GlobalMethodsClass::WeightingMethod_t method) {

  double totEngy(0.0), xHit(0.0), yHit(0.0), zHit(0.0), thetaHit(0.0), weightSum(0.0);
  int loopFlag = 1;
  VecCalHit caloHits;
  while(loopFlag == 1) {
    caloHits.clear();
    for (VInt::const_iterator it = cellIdV.begin(); it != cellIdV.end(); ++it) {
      const int k = it - cellIdV.begin();
      const auto& calHit = calHitsCellId.at(*it);
      caloHits.push_back(calHit);
      const double weightHit = posWeightTrueCluster(calHit,cellEngyV[k],method);
      weightSum += weightHit;
      const double* position = calHit->getPosition();

      xHit      += position[0] * weightHit;
      yHit      += position[1] * weightHit;
      zHit      += position[2] * weightHit;
      totEngy   += cellEngyV[k];
    }
    if(weightSum > 0) {
      xHit     /= weightSum;   yHit   /= weightSum; zHit /= weightSum;
      thetaHit  = atan( sqrt( xHit*xHit + yHit*yHit)/fabs( zHit ));
      //      thetaHit /= weightSum;
      loopFlag = 0;

    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = GlobalMethodsClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightSum = 0.;
    }
  }

  clusterCM[clusterId] = LCCluster(totEngy, xHit, yHit, zHit, weightSum, method, thetaHit, 0.0, caloHits);
}

/* --------------------------------------------------------------------------
   compute center of mass of each cluster
   (2). update the map clusterCM with the new cal hit
   -------------------------------------------------------------------------- */
template <> void LumiCalClustererClass::updateEngyPosCM(CalHit const& calHit, LCCluster& clusterCM) {
  double	engyHit = (double)calHit->getEnergy();
  GlobalMethodsClass::WeightingMethod_t method =  clusterCM.getMethod();

  clusterCM.addToEnergy(engyHit);

  double weightHit = posWeight(calHit,method);

  if(weightHit > 0){
    double weightCM = clusterCM.getWeight();

    double xCM = clusterCM.getX() * weightCM;
    double yCM = clusterCM.getY() * weightCM;
    double zCM = clusterCM.getZ() * weightCM;

    double xHit = double(calHit->getPosition()[0]) * weightHit;
    double yHit = double(calHit->getPosition()[1]) * weightHit;
    double zHit = double(calHit->getPosition()[2]) * weightHit;

    weightCM += weightHit;
    xCM      += xHit;	xCM	/= weightCM;
    yCM      += yHit;	yCM	/= weightCM;
    zCM      += zHit;	zCM	/= weightCM;
    clusterCM.setPosition( xCM, yCM, zCM );
    clusterCM.setWeight(weightCM);

  }

}

/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::checkClusterMergeCM(  int clusterId1, int clusterId2,
						 MapIntVInt const& clusterIdToCellId,
						 MapIntCalHit const& calHitsCellId,
						 double distanceAroundCM, double percentOfEngyAroungCM,
						 GlobalMethodsClass::WeightingMethod_t method ){

  // std::vector for holding the Ids of clusters
  VInt cellIdV;
  double totEngyAroundCM = 0.;

  // add to cellIdV hits from both clusters, and sum up each cluster's energy
  for(VInt::const_iterator cellIt = clusterIdToCellId.at(clusterId1).begin();
      cellIt != clusterIdToCellId.at(clusterId1).end();
      ++cellIt) {
    cellIdV.push_back(*cellIt);
  }

  for(VInt::const_iterator cellIt = clusterIdToCellId.at(clusterId2).begin();
      cellIt != clusterIdToCellId.at(clusterId2).end();
      ++cellIt) {
    cellIdV.push_back(*cellIt);
  }

  LCCluster engyPosCM( calculateEngyPosCM(cellIdV, calHitsCellId, method) );
  double CM1[3] = { engyPosCM.getX(), engyPosCM.getY(), 0.0 };
  const double engyCM = engyPosCM.getE();

  // check that the CM position has a significant amount of the energy around it
  for( VInt::iterator cellIt = cellIdV.begin();
       cellIt != cellIdV.end();
       ++cellIt ) {
    const auto&  hit        = calHitsCellId.at(*cellIt);
    const double distanceCM = distance2D(CM1,hit->getPosition());
    if(distanceCM < distanceAroundCM) {
      double engyHit = hit->getEnergy();
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
double LumiCalClustererClass::getEngyInMoliereFraction(	MapIntCalHit  const& calHitsCellId,
							VInt const&, //clusterIdToCellId,
							LCCluster const& clusterCM,
							double moliereFraction   ){

  const double distanceToScan = _moliereRadius * moliereFraction;
  double engyAroundCM = 0.0;

  for(MapIntCalHit::const_iterator calHitsCellIdIterator = calHitsCellId.begin();
      calHitsCellIdIterator != calHitsCellId.end();
      ++calHitsCellIdIterator) {
    const auto&  calHit     = calHitsCellIdIterator->second;
    const double distanceCM = distance2D(clusterCM.getPosition(),calHit->getPosition());
    if(distanceCM < distanceToScan)
      engyAroundCM += calHit->getEnergy();
  }

  return engyAroundCM;

}


// overloaded with different variables and functionality...
double LumiCalClustererClass::getEngyInMoliereFraction( MapIntCalHit  const& calHitsCellId,
							VInt const&,//clusterIdToCellId,
							LCCluster const& clusterCM,
							double moliereFraction,
							MapIntInt & flag  ){

  const double distanceToScan = _moliereRadius * moliereFraction;
  double engyAroundCM = 0.;

  for(MapIntCalHit::const_iterator calHitsCellIdIterator = calHitsCellId.begin();
      calHitsCellIdIterator != calHitsCellId.end();
      ++calHitsCellIdIterator) {
    const auto&  calHit     = calHitsCellIdIterator->second;
    const double distanceCM = distance2D(clusterCM.getPosition(),calHit->getPosition());

    const int cellIdHit = calHitsCellIdIterator->first;
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
void LumiCalClustererClass::getThetaPhiZCluster( MapIntCalHit  const& calHitsCellId,
                                                 VInt const& clusterIdToCellId,
                                                 double totEngy, double * output   ) {

  double zCluster = 0., xCluster = 0., yCluster = 0.;
  double weightSum = -1., logWeightConstFactor = 0.;

  // (BP) hardwired method ?!
  GlobalMethodsClass::WeightingMethod_t method = GlobalMethodsClass::LogMethod;

  while(weightSum < 0){
    for( VInt::const_iterator cellIt = clusterIdToCellId.begin();
	 cellIt != clusterIdToCellId.end();
	 ++cellIt ) {
      const auto& hit = calHitsCellId.at(*cellIt);
      //(BP) posWeight returns always weight >= 0.
      const double weightHit = posWeight(hit, totEngy, method, (_logWeightConst + logWeightConstFactor) );
      if(weightHit > 0){
	if(weightSum < 0) weightSum = 0.;
	weightSum += weightHit;
	xCluster += double(hit->getPosition()[0]) * weightHit;
	yCluster += double(hit->getPosition()[1]) * weightHit;
	zCluster += double(hit->getPosition()[2]) * weightHit;
      }
    }
    // (BP) even one hit with positive weight make this false, is it what we want ?
    // since weightHit is always >= 0  this will happen if all hits weights are 0 !
    // and reducing logWeightConstFactor does not help in case weightSum < 0 (?!).....
    if(weightSum < 0) logWeightConstFactor -= .5;
  }

  xCluster /= weightSum; yCluster /= weightSum; zCluster /= weightSum;

  output[0] = atan( sqrt( xCluster*xCluster + yCluster*yCluster )/fabs( zCluster ));
  output[1] = atan2( yCluster, xCluster );
  output[2] = zCluster;

  return ;
}


/* --------------------------------------------------------------------------
   compute the theta/phi of a cell with a given cellId
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::thetaPhiCell(int cellId, GlobalMethodsClass::Coordinate_t output) {

  int	cellIdR, cellIdZ, cellIdPhi, arm;
  GlobalMethodsClass::CellIdZPR(cellId, cellIdZ, cellIdPhi, cellIdR, arm);
  double	rCell, zCell, thetaCell, phiCell, outputVal(-1);

  if(output == GlobalMethodsClass::COTheta) {
    rCell      = _rMin + (cellIdR + .5) * _rCellLength;
    zCell      = fabs(_zFirstLayer) + _zLayerThickness * (cellIdZ - 1);
    thetaCell  = atan(rCell / zCell);
    outputVal  = thetaCell;
  }
  else if(output == GlobalMethodsClass::COPhi) {
    //(BP) account for possible layer offset
    //    phiCell   = 2*M_PI * (double(cellIdPhi) + .5) / _cellPhiMax;
    phiCell   = (double(cellIdPhi) + .0) * _phiCellLength + double( (cellIdZ)%2 ) * _zLayerPhiOffset;
    phiCell = ( phiCell > M_PI ) ? phiCell - 2.*M_PI : phiCell;
    outputVal = phiCell;
  }

  return outputVal;

}



/* --------------------------------------------------------------------------
   get the energy around a cluster CM within a distanceToScan raduis
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::getMoliereRadius( MapIntCalHit const& calHitsCellId,
                                                VInt const& clusterIdToCellId,
                                                LCCluster const& clusterCM ) {

  const double engyPercentage = .9;

  VVDouble clusterHitsEngyPos(clusterIdToCellId.size(), VDouble(2, 0.0));

#if _IMPROVE_PROFILE_LAYER_POS == 1
  double thetaZClusterV[3];
  getThetaPhiZCluster(calHitsCellId, clusterIdToCellId, clusterCM[0],thetaZClusterV);
  double	thetaCluster = thetaZClusterV[0];
  double	zCluster     = Abs(thetaZClusterV[2]);

  int	layerMiddle  = int( ( zCluster - Abs(_zFirstLayer) ) / _zLayerThickness ) + 1;
#endif

  for( VInt::const_iterator cellIt = clusterIdToCellId.begin();
       cellIt != clusterIdToCellId.end();
       ++cellIt ){
    const int cellIdHit = *cellIt;
    const int hitNow = cellIt - clusterIdToCellId.begin();
    const auto& thisHit   = calHitsCellId.at(cellIdHit);
    double CM2[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

#if _IMPROVE_PROFILE_LAYER_POS == 1
    int	layerHit   = GlobalMethodsClass::CellIdZPR(cellIdHit,GlobalMethodsClass::COZ);

    // projection hits have an incoded layer number of (_maxLayerToAnalyse + 1)
    // the original hit's layer number is stored in (the previously unused) CellID1
    if(layerHit == (_maxLayerToAnalyse + 1))
      layerHit = thisHit->getCellID1();

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
  double distanceCM(0.0), engyAroundCM(0.0);
  for( VVDouble::iterator energyIt = clusterHitsEngyPos.begin();
       energyIt != clusterHitsEngyPos.end();
       ++energyIt) {

    if (engyAroundCM < (engyPercentage * clusterCM.getE())) {
      engyAroundCM += (*energyIt)[0];
      distanceCM    = (*energyIt)[1];
    } else {
      break;
    }

  }

  return distanceCM;
}




/* --------------------------------------------------------------------------
   make sure that the CM of two merged clusters is where most of the merged
   cluster's energy is deposited
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::getDistanceAroundCMWithEnergyPercent( LCCluster const& clusterCM,
								    VInt const& clusterIdToCellId,
								    MapIntCalHit const& calHitsCellId,
								    double engyPercentage ) {
  VVDouble clusterHitsEngyPos (clusterIdToCellId.size(), VDouble(2, 0.0));

  // fill a std::vector with energy, position and distance from CM of every cal hit
  for( VInt::const_iterator cellIt = clusterIdToCellId.begin();
       cellIt != clusterIdToCellId.end();
       ++cellIt) {

    const int clusterId = cellIt - clusterIdToCellId.begin();
    const auto& thisHit   = calHitsCellId.at(*cellIt);

    clusterHitsEngyPos[clusterId][0] = thisHit->getEnergy();
    clusterHitsEngyPos[clusterId][1] = distance2D(clusterCM.getPosition(),thisHit->getPosition());

  }

  // sort the std::vector of the cal hits according to distance from CM in ascending order (shortest distance is first)
  sort (clusterHitsEngyPos.begin(), clusterHitsEngyPos.end(), HitDistanceCMCmpAsc<1> );

  double engyAroundCM(0.0), distanceCM(0.0);
  for( VVDouble::iterator energyIt = clusterHitsEngyPos.begin();
       energyIt != clusterHitsEngyPos.end();
       ++energyIt ) {
    if (engyAroundCM < engyPercentage * clusterCM.getE()) {
      engyAroundCM += (*energyIt)[0];
      distanceCM    = (*energyIt)[1];
    }
  }

  return distanceCM;

}


void LumiCalClustererClass::dumpClusters( MapIntLCCluster const& clusterCM ) {
  for( MapIntLCCluster::const_iterator clusterCMIterator = clusterCM.begin();
       clusterCMIterator != clusterCM.end();
       ++clusterCMIterator ) {
    int clusterId = clusterCMIterator->first;
    streamlog_out(DEBUG5) << "\t\t cluster Id, pos(x,y,z), engy(weight): "
			  << clusterId << "\t ("
			  << clusterCMIterator->second.getX() << " , "
			  << clusterCMIterator->second.getY() << " , "
			  << clusterCMIterator->second.getZ() << ") \t "
			  << clusterCMIterator->second.getE() <<" ( "
			  << clusterCMIterator->second.getWeight() <<" ) " <<std::endl;
  }

}
