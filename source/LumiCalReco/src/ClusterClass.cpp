
#include "Global.hh"

#include "ClusterClass.h"

#include <EVENT/MCParticle.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <cmath>
#include <iomanip>

using namespace streamlog;

/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
ClusterClass::ClusterClass(int idNow):
  GlobalMethodsClass(),
  Id(idNow), Pdg(0), SignMC(0), ParentId(idNow), NumMCDaughters(0),
  OutsideFlag(0), MergedFlag(0), HighestEnergyFlag(0), ModifiedFlag(0),NumHits(0),
  Engy(0.0), Theta(0.0), Phi(0.0), RZStart(0.0),
  VtxX(0.0), VtxY(0.0), VtxZ(0.0), EndPointX(0.0), EndPointY(0.0), EndPointZ(0.0),
  EngyMC(0.0), ThetaMC(0.0), PhiMC(0.0),
  MergedV(0),
  OutsideReason(""),
  Hit()
{
  // inherited method from GlobalMethods class, to set all global constants
  SetConstants();
  clusterPosition[0] = 0.0;
  clusterPosition[1] = 0.0;
  clusterPosition[2] = 0.0;
}

ClusterClass::~ClusterClass() {

}

void ClusterClass::SetStatsMC(EVENT::MCParticle * mcParticle) {

  Pdg = (int)mcParticle->getPDG();
  Id  = (int)mcParticle -> id();
  NumMCDaughters = (int)mcParticle -> getDaughters().size();

  VtxX = (double)mcParticle->getVertex()[0];
  VtxY = (double)mcParticle->getVertex()[1];
  VtxZ = (double)mcParticle->getVertex()[2];

  EndPointX = (double)mcParticle->getEndpoint()[0];
  EndPointY = (double)mcParticle->getEndpoint()[1];
  EndPointZ = (double)mcParticle->getEndpoint()[2];

  EngyMC        = (double)mcParticle->getEnergy();

  const double MCmomX = (double)mcParticle->getMomentum()[0];
  const double MCmomY = (double)mcParticle->getMomentum()[1];
  const double MCmomZ = (double)mcParticle->getMomentum()[2];

  SignMC = int(MCmomZ / fabs(MCmomZ));
  // unboost momentum
  double betgam = tan( GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2. );
  double gama = sqrt( 1.+ betgam*betgam);
  double locPX = betgam*EngyMC - gama*MCmomX;

  double MCr = sqrt( locPX*locPX  + MCmomY*MCmomY);
  ThetaMC = atan(MCr / fabs(MCmomZ));                  // Theta in local LCal coordinates

  // correction for the polar angle for cases where the particle does not
  // come from the IP (according to a projection on the face of LumiCal)
  if( fabs(VtxZ) > _VERY_SMALL_NUMBER ) {
    //    cout <<"   ---   "<< Id << "\t" << ThetaMC << "\t -> \t" ;
    MCr = tan(ThetaMC) * (GlobalParamD[ZStart] - fabs(VtxZ));
    ThetaMC = atan( MCr / GlobalParamD[ZStart] );
    // cout << ThetaMC << "\t r,z:  " << MCr <<"\t" << Abs(vtxZ) <<  "\t"
    //      <<  mcParticle->getDaughters().size() << endl ;
    // EVENT::MCParticle *mcParticleParent = mcParticle -> getDaughters()[0];
    // cout << "\t\t\t"<< mcParticleParent->getVertex()[2]
    //      <<  "\t"<< mcParticleParent->id() <<  "\t"<< endl ;
  }
  // (BP)get phi angle in the range (0.; +- M_PI) - default in LC  
  PhiMC = atan2( MCmomY, MCmomX );
  /*(BP) replaced with above 
  PhiMC = atan(fabs(MCmomY / MCmomX));
  if(MCmomZ > 0) {
    // if(MCmomX>0 && MCmomY>0)   PhiMC = PhiMC;//no change, unneeded
    if(MCmomX<0 && MCmomY>0)      PhiMC = M_PI   - PhiMC;
    else if(MCmomX<0 && MCmomY<0) PhiMC = M_PI   + PhiMC;
    else if(MCmomX>0 && MCmomY<0) PhiMC = 2*M_PI - PhiMC;
  } else {
    if(MCmomX>0 && MCmomY>0)      PhiMC = M_PI   - PhiMC;
    //if(MCmomX<0 && MCmomY>0)    PhiMC = PhiMC;//no change, unneeded
    else if(MCmomX<0 && MCmomY<0) PhiMC = 2*M_PI - PhiMC;
    else if(MCmomX>0 && MCmomY<0) PhiMC = M_PI   + PhiMC;
  }
  */


  // find the original parent of the particle
  EVENT::MCParticle *mcParticleParent = mcParticle;
  ParentId = Id;
  while(1) {
    if( !(mcParticleParent -> isCreatedInSimulation()) )
      break;

    mcParticleParent = mcParticleParent -> getParents()[0];
    ParentId = (int)mcParticleParent -> id();
  }


  return;
}



void ClusterClass::SetStatsMC() {

  ParentId = -1;
  Pdg = Id = SignMC = NumHits = NumMCDaughters = 0;
  VtxX = VtxY = VtxZ = 0.;
  EndPointX = EndPointY = EndPointZ = 0.;
  ThetaMC = PhiMC = EngyMC = 0.;


  return;
}



void ClusterClass::FillHit(int cellNow, double engyNow) {

  Hit[cellNow] += engyNow;
  Engy         += engyNow;

  return;
}


int ClusterClass::ResetStats() {

  double thetaSum(0.0), weightSum(0.0), engySum(0.0);
  std::map < GlobalMethodsClass::Coordinate_t , double > thetaPhiCellV;

  // initialize modification flag
  ModifiedFlag = 0;

  // re-sum energy from all cells (just in case)
  NumHits = 0;
  for( std::map < int , double > :: iterator hitIterator = Hit.begin();
       hitIterator != Hit.end(); ++hitIterator) {
    engySum += hitIterator->second;
    NumHits++;
  }
  Engy = engySum;

  // if the particle has no deposits in LumiCal
  // (BP) this already done in getCalHits method !?
  if(Engy < _VERY_SMALL_NUMBER) {
    OutsideFlag = 1;
    OutsideReason = "No energy deposits at all";
    return 0;
  }

  if( NumHits < GlobalParamI[ClusterMinNumHits] ){
    OutsideFlag = 1;
    OutsideReason = " Number of hits below minimum";
    return 0;
  }

  double xtemp(0.0), ytemp(0.0), ztemp(0.0);
  for( std::map < int , double > :: iterator hitIterator = Hit.begin();
       hitIterator != Hit.end(); ++hitIterator) {
    int cellId = hitIterator->first;
    ThetaPhiCell(cellId , thetaPhiCellV);

    const double engyNow = hitIterator->second;
    const double weightNow = GlobalParamD[LogWeightConstant] + log(engyNow/engySum);

    if(weightNow < 0) continue;

    weightSum += weightNow;
    thetaSum  += thetaPhiCellV[GlobalMethodsClass::COTheta] * weightNow;
    const double phi = thetaPhiCellV[GlobalMethodsClass::COPhi];
    ztemp += thetaPhiCellV[GlobalMethodsClass::COZ] * weightNow;
    xtemp += cos(phi) * weightNow;
    ytemp += sin(phi) * weightNow;
  }


  Theta    = thetaSum / weightSum;
  double xCOG = xtemp / weightSum;
  double yCOG = ytemp / weightSum;
  double zCOG = ztemp / weightSum;
  RZStart  = atan(fabs(Theta)) * GlobalParamD[ZStart];
  //  RZStart  = atan(fabs(Theta)) * zCOG;
  Phi = atan2(yCOG, xCOG);
  // (BP) Phi = ( Phi >= 0. ) ? Phi : Phi + 2.0 * M_PI;
  
  clusterPosition[0] = RZStart*cos(Phi);
  clusterPosition[1] = RZStart*sin(Phi);
  clusterPosition[2] = zCOG;
  
  if( Theta < GlobalParamD[ThetaMin] || Theta > GlobalParamD[ThetaMax] ){
    OutsideFlag = 1;
    OutsideReason = "Reconstructed outside the fiducial volume";
    return 0;
  }

  if(engySum < GlobalParamD[MinClusterEngySignal]) {
    OutsideFlag = 1;
    OutsideReason = "Cluster energy below minimum";
    return 0;
  }

#if _CLUSTER_RESET_STATS_DEBUG == 1
  streamlog_out( MESSAGE4 ) << std::endl;
    << "Cluster Information:   " << Id << std::endl
    << std::setw(30) << "sign, engyHits, engyMC:  "
    << std::setw(13) << SignMC
    << std::setw(13) << SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , Engy)
    << std::setw(13) << EngyMC << std::endl
    << std::setw(30) << "theta mc,rec:"
    << std::setw(13) << ThetaMC
    << std::setw(13) << Theta << std::endl
    << std::setw(30) << "phi mc,rec:"
    << std::setw(13) << PhiMC
    << std::setw(13) << Phi << std::endl
    << std::setw(30) << "vertex X,Y,Z: "
    << std::setw(13) << VtxX
    << std::setw(13) << VtxY
    << std::setw(13) << VtxZ
    << std::endl;
 #endif

  return 1;
}

