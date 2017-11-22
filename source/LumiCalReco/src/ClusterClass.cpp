
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
ClusterClass::ClusterClass(GlobalMethodsClass& _gmc, int idNow, LCCluster const* thisCluster)
    : m_clusterInfo(thisCluster), Id(idNow), gmc(_gmc) {
  ResetStats();
}

void ClusterClass::SetStatsMC(SMCInfo const& mcInfo) {
  m_mcInfo = mcInfo;

  return;
}



void ClusterClass::SetStatsMC() {
  Id       = 0;
  m_mcInfo = std::make_shared<MCInfo>();
  return;
}

int ClusterClass::ResetStats() {
  // if the cluster has no deposits in LumiCal
  if (m_clusterInfo->getEnergy() < _VERY_SMALL_NUMBER) {
    OutsideFlag = 1;
    OutsideReason = "No energy deposits at all";
    return 0;
  }

  if (int(m_clusterInfo->getCaloHits().size()) < gmc.GlobalParamI[GlobalMethodsClass::ClusterMinNumHits]) {
    OutsideFlag = 1;
    OutsideReason = " Number of hits below minimum";
  }

  if (m_clusterInfo->getTheta() < gmc.GlobalParamD[GlobalMethodsClass::ThetaMin] ||
      m_clusterInfo->getTheta() > gmc.GlobalParamD[GlobalMethodsClass::ThetaMax]) {
    OutsideFlag = 1;
    OutsideReason = "Reconstructed outside the fiducial volume";
  }

  if (m_clusterInfo->getEnergy() < gmc.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV]) {
    OutsideFlag = 1;
    OutsideReason = "Cluster energy below minimum";
  }


  return 1;
}

void ClusterClass::PrintInfo(){
  // clang-format off
  streamlog_out(DEBUG6) << std::endl
                        << "ClusterClass Information:   " << Id << std::endl
                        << std::setw(30) << "sign, engyHits, engyMC, nHits:"
                        << std::setw(13) << m_mcInfo->sign
                        << std::setw(13) << m_clusterInfo->getEnergy()
                        << std::setw(13) << m_mcInfo->engy
                        << std::setw(13) << m_clusterInfo->getCaloHits().size() << std::endl
                        << std::setw(30) << "theta mc,rec:  "
                        << std::setw(13) << m_mcInfo->theta
                        << std::setw(13) << m_clusterInfo->getTheta() << std::endl
                        << std::setw(30) << "phi mc,rec:  "
                        << std::setw(13) << m_mcInfo->phi
                        << std::setw(13) << m_clusterInfo->getPhi() << std::endl
                        << std::setw(30) << "vertex X,Y,Z:  "
                        << std::setw(13) << m_mcInfo->m_vtxX
                        << std::setw(13) << m_mcInfo->m_vtxY
                        << std::setw(13) << m_mcInfo->m_vtxZ << std::endl
                        << std::setw(30) << "CLstart  X,Y,Z:  "
                        << std::setw(13) << m_clusterInfo->getPositionAtFront()[0]
                        << std::setw(13) << m_clusterInfo->getPositionAtFront()[1]
                        << std::setw(13) << m_clusterInfo->getPositionAtFront()[2] << std::endl
                        << std::setw(30) << "MCstart  X,Y,Z:  "
                        << std::setw(13) << m_mcInfo->mcPosition[0]
                        << std::setw(13) << m_mcInfo->mcPosition[1]
                        << std::setw(13) << m_mcInfo->mcPosition[2] << std::endl
                        << std::endl;
  // clang-format on
}

std::ostream& operator<<(std::ostream& o, const ClusterClass& rhs) {
  // clang-format off
  o << std::setw(20) << "X, Y, Z:" << std::endl << std::fixed << std::setprecision(3)
    << std::setw(13) << rhs.getPosition()[0]
    << std::setw(13) << rhs.getPosition()[1]
    << std::setw(13) << rhs.getPosition()[2] << std::endl
    << std::setw(20) << "Energy, Theta, Phi: " << std::endl
    << std::setw(13) << rhs.getEnergy()
    << std::setw(13) << rhs.getTheta()
    << std::setw(13) << rhs.getPhi() << std::endl
    << std::endl;
  return o;
  // clang-format on
}
