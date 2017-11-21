
#include "Global.hh"

#include "ClusterClass.h"

#include <EVENT/MCParticle.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <cmath>
#include <iomanip>

using namespace streamlog;

/* --------------------------------------------------------------------------
   This class holds information about the clusters
   -------------------------------------------------------------------------- */
ClusterClass::ClusterClass(GlobalMethodsClass& _gmc, int idNow, int armNow, VInt const& cells, VDouble const& cellEnergies)
    : Id(idNow),
      Pdg(0),
      SignMC(armNow),
      ParentId(-1),
      NumMCDaughters(0),
      OutsideFlag(0),
      MatchFlag(0),
      HighestEnergyFlag(0),
      ModifiedFlag(0),
      NumHits(0),
      Engy(0.0),
      Theta(0.0),
      Phi(0.0),
      RZStart(0.0),
      VtxX(0.0),
      VtxY(0.0),
      VtxZ(0.0),
      EndPointX(0.0),
      EndPointY(0.0),
      EndPointZ(0.0),
      EngyMC(0.0),
      ThetaMC(0.0),
      PhiMC(0.0),
      DiffTheta(0.),
      DiffPosXY(0.),
      MergedV(0),
      OutsideReason(""),
      Hit(),
      gmc(_gmc) {
  int cellNow = 0;
  for (VInt::const_iterator cellIt = cells.begin(); cellIt != cells.end(); ++cellIt, ++cellNow) {
    FillHit(*cellIt, cellEnergies.at(cellNow));  //Adds to Engy!
  }

  clusterPosition[0] = 0.0;
  clusterPosition[1] = 0.0;
  clusterPosition[2] = 0.0;
  mcpPosition[0] = 0.0;
  mcpPosition[1] = 0.0;
  mcpPosition[2] = 0.0;

  ResetStats();  // calculate energy, position for the cluster
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

  SignMC = int (MCmomZ / fabs(MCmomZ));
  // unboost momentum
  double betgam = tan( gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2. );
  double gama = sqrt( 1.+ betgam*betgam);
  double localPX =-betgam*EngyMC + gama*MCmomX;
  // position at face of Lcal
  mcpPosition[0] = VtxX + double(SignMC)*gmc.GlobalParamD[GlobalMethodsClass::ZStart]*tan(localPX /MCmomZ);
  mcpPosition[1] = VtxY + double(SignMC)*gmc.GlobalParamD[GlobalMethodsClass::ZStart]*tan(MCmomY/MCmomZ);
  mcpPosition[2] = gmc.GlobalParamD[GlobalMethodsClass::ZStart]*double(SignMC);

  double MCr = sqrt( localPX*localPX  + MCmomY*MCmomY);
  ThetaMC = atan(MCr / fabs(MCmomZ));                  // Theta in local LCal coordinates
  // (BP)get phi angle in the range (0.; +- M_PI) - default in LC  
  PhiMC = atan2( mcpPosition[1], mcpPosition[0]);

  // correction for the polar angle for cases where the particle does not
  // come from the IP (according to a projection on the face of LumiCal)
  if( fabs(VtxZ) > _VERY_SMALL_NUMBER ) {
    MCr = tan(ThetaMC) * (gmc.GlobalParamD[GlobalMethodsClass::ZStart] - fabs(VtxZ));
    ThetaMC = atan( MCr / gmc.GlobalParamD[GlobalMethodsClass::ZStart] );
  }


  // find the original parent of the particle
  EVENT::MCParticle *mcParticleParent = mcParticle;
  ParentId = Id;
  while(1) {
    if( !(mcParticleParent -> isCreatedInSimulation()) or mcParticleParent -> getParents().empty() )
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

  // initialize modification flag
  ModifiedFlag = 0;

  // re-sum energy from all cells (just in case)
  NumHits = Hit.size();
  Engy    = 0.0;
  for (auto const& pairIDEng : Hit) {
    Engy += pairIDEng.second;
  }

  // if the particle has no deposits in LumiCal
  // (BP) this already done in getCalHits method !?
  if(Engy < _VERY_SMALL_NUMBER) {
    OutsideFlag = 1;
    OutsideReason = "No energy deposits at all";
    return 0;
  }

  if( NumHits < gmc.GlobalParamI[GlobalMethodsClass::ClusterMinNumHits] ){
    OutsideFlag = 1;
    OutsideReason = " Number of hits below minimum";
  }

  const double logConstant(gmc.GlobalParamD[GlobalMethodsClass::LogWeightConstant]);
  const double rMin(gmc.GlobalParamD[GlobalMethodsClass::RMin]);

  double       xtemp(0.0), ytemp(0.0), ztemp(0.0), thetaSum(0.0), weightSum(0.0);
  std::map<GlobalMethodsClass::Coordinate_t, double> thetaPhiCellV;

  for (auto const& pairIDEng : Hit) {
    const int cellId = pairIDEng.first;
    gmc.ThetaPhiCell(cellId, thetaPhiCellV);
    const double phi   = thetaPhiCellV[GlobalMethodsClass::COPhi];
    const double rcell = thetaPhiCellV[GlobalMethodsClass::COR];

    const double weightNow =
        GlobalMethodsClass::posWeight(pairIDEng.second, Engy, GlobalMethodsClass::LogMethod, logConstant) * rMin / rcell;

    if (not(weightNow > 0))
      continue;

    weightSum += weightNow;
    thetaSum  += thetaPhiCellV[GlobalMethodsClass::COTheta] * weightNow;
    xtemp += rcell*cos(phi) * weightNow;
    ytemp += rcell*sin(phi) * weightNow;
    ztemp += thetaPhiCellV[GlobalMethodsClass::COZ] * weightNow;
  }


  Theta  = thetaSum / weightSum;
  xtemp /= weightSum;
  ytemp /= weightSum;
  ztemp /= weightSum;
  RZStart  = atan(fabs(Theta)) * gmc.GlobalParamD[GlobalMethodsClass::ZStart];
  Phi = atan2(ytemp, xtemp);

  const double r     = sqrt(xtemp * xtemp + ytemp * ytemp + ztemp * ztemp);
  clusterPosition[0] = r * sin(Theta) * cos(Phi);
  clusterPosition[1] = r * sin(Theta) * sin(Phi);
  clusterPosition[2] = r * cos(Theta);

  if( Theta < gmc.GlobalParamD[GlobalMethodsClass::ThetaMin] || Theta > gmc.GlobalParamD[GlobalMethodsClass::ThetaMax] ){
    OutsideFlag = 1;
    OutsideReason = "Reconstructed outside the fiducial volume";
  }

  if (Engy < gmc.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV]) {
    OutsideFlag = 1;
    OutsideReason = "Cluster energy below minimum";
  }


  return 1;
}

void ClusterClass::PrintInfo(){
  // clang-format off
  streamlog_out( MESSAGE4 ) << std::endl
    << "ClusterClass Information:   " << Id << std::endl
    << std::setw(30) << "sign, engyHits, engyMC, nHits:"
    << std::setw(13) << SignMC
    << std::setw(13) << Engy
    << std::setw(13) << EngyMC
    << std::setw(13) << Hit.size() << std::endl
    << std::setw(30) << "theta mc,rec:  "
    << std::setw(13) << ThetaMC
    << std::setw(13) << Theta << std::endl
    << std::setw(30) << "phi mc,rec:  "
    << std::setw(13) << PhiMC
    << std::setw(13) << Phi << std::endl
    << std::setw(30) << "vertex X,Y,Z:  "
    << std::setw(13) << VtxX
    << std::setw(13) << VtxY
    << std::setw(13) << VtxZ << std::endl
    << std::setw(30) << "CLstart  X,Y,Z:  "
    /*
    << std::setw(13) << clusterPosition[0]
    << std::setw(13) << clusterPosition[1]
    << std::setw(13) << clusterPosition[2] << std::endl
    */
			    << std::setw(13) << RZStart*cos(Phi)
			    << std::setw(13) << RZStart*sin(Phi)
			    << std::setw(13) << gmc.GlobalParamD[GlobalMethodsClass::ZStart]*double(SignMC) << std::endl
    << std::setw(30) << "MCstart  X,Y,Z:  "
    << std::setw(13) << mcpPosition[0]
    << std::setw(13) << mcpPosition[1]
    << std::setw(13) << mcpPosition[2] << std::endl
    << std::endl;
  // clang-format on
}
