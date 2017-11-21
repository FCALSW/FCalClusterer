#ifndef ClusterClass_h
#define ClusterClass_h 1

#include "GlobalMethodsClass.h"

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
  ClusterClass(GlobalMethodsClass& gmc, int idNow, int armNow, VInt const& cells, VDouble const& cellEnergies);
  ~ClusterClass() = default;

  void	FillHit(int cellNow, double engyNow);
  int	ResetStats();
  void	SetStatsMC(EVENT::MCParticle * mcParticle);
  void	SetStatsMC();
  void  PrintInfo();

  int   	Id, Pdg, SignMC, ParentId, NumMCDaughters;
  int   	OutsideFlag, MatchFlag, HighestEnergyFlag, ModifiedFlag;
  int           NumHits;
  double	Engy, Theta, Phi, RZStart;
  double	VtxX, VtxY, VtxZ, EndPointX, EndPointY, EndPointZ;
  double	EngyMC, ThetaMC, PhiMC;
  double        DiffTheta, DiffPosXY;
  double        clusterPosition[3];           // cluster position at Zbeg
  double        mcpPosition[3];               // MCparticle position at Zbeg

  std::vector <int> MergedV;

  std::string	OutsideReason;

  std::map < int , double >	Hit;

  const double * getPosition() const { return clusterPosition; }
  
  friend std::ostream& operator<<(std::ostream & o, const ClusterClass& rhs);
  
  GlobalMethodsClass& gmc;

};

#endif // ClusterClass_h
