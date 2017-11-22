#ifndef ClusterClass_h
#define ClusterClass_h 1

#include "GlobalMethodsClass.h"
#include "LCCluster.hh"
#include "MCInfo.h"

#include <map>
#include <string>
#include <vector>

namespace EVENT {
  class MCParticle;
}

/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
class ClusterClass {

public:
  ClusterClass(GlobalMethodsClass& gmc, int idNow, LCCluster const* clusterInfo);
  ~ClusterClass() = default;
  ClusterClass(const ClusterClass&) = default;
  ClusterClass& operator=(const ClusterClass&) = default;

  void	FillHit(int cellNow, double engyNow);
  int	ResetStats();
  void SetStatsMC(SMCInfo const& mcParticle);
  void	SetStatsMC();
  void  PrintInfo();

private:
  LCCluster const* m_clusterInfo;
  SMCInfo          m_mcInfo = std::make_shared<MCInfo>();

public:
  int         Id        = 0;
  double      DiffTheta = 0.0, DiffPosXY = 0.0;
  int         OutsideFlag = 0, HighestEnergyFlag = 0;
  std::string OutsideReason{};

  const double* getPosition() const { return m_clusterInfo->getPosition(); }
  const double* getPositionAtFront() const { return m_clusterInfo->getPositionAtFront(); }

  const double* getMCPosition() const { return m_mcInfo->mcPosition; }

  friend std::ostream& operator<<(std::ostream & o, const ClusterClass& rhs);
  
  GlobalMethodsClass& gmc;

  double getEnergy() const { return m_clusterInfo->getEnergy(); }
  double getTheta() const { return m_clusterInfo->getTheta(); }
  double getPhi() const { return m_clusterInfo->getPhi(); }
  double getRZStart() const { return m_clusterInfo->getRZStart(); }
  int    getNHits() const { return m_clusterInfo->getCaloHits().size(); }

  double getEnergyMC() const { return m_mcInfo->engy; }
  double getThetaMC() const { return m_mcInfo->theta; }
  double getPhiMC() const { return m_mcInfo->phi; }

  int getPDG() const { return m_mcInfo->pdg; }
  int getSign() const { return m_mcInfo->sign; }
  int getParentId() const { return m_mcInfo->m_parentID; }
  int getNumMCDaughters() const { return m_mcInfo->m_numMCDaughters; }
};

#endif // ClusterClass_h
