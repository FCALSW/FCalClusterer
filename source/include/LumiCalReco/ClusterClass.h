#ifndef ClusterClass_h
#define ClusterClass_h 1

#include "GlobalMethodsClass.h"

#include <map>
#include <string>
#include <vector>
#include <iomanip>

namespace EVENT {
  class MCParticle;
}

/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
class ClusterClass : public GlobalMethodsClass {

public:
  ClusterClass(int idNow);
  ~ClusterClass();

  void	FillHit(int cellNow, double engyNow);
  int	ResetStats();
  void	SetStatsMC(EVENT::MCParticle * mcParticle);
  void	SetStatsMC();

  int	Id, Pdg, SignMC, ParentId, NumMCDaughters;
  int	OutsideFlag, MergedFlag, HighestEnergyFlag, ModifiedFlag;
  double	Engy, Theta, Phi, RZStart;
  double	VtxX, VtxY, VtxZ, EndPointX, EndPointY, EndPointZ;
  double	EngyMC, ThetaMC, PhiMC;

  std::vector <int> MergedV;

  std::string	OutsideReason;

  std::map < int , double >	Hit;

};

#endif // ClusterClass_h
