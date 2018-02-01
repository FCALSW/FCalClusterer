#ifndef ProjectionInfo_hh
#define ProjectionInfo_hh 1

#include <memory> // IWYU pragma: keep
#include <set>

// IWYU pragma: no_include <bits/shared_ptr.h>

namespace EVENT{
  class CalorimeterHit;
}

class LumiCalHit;

using CalHit = std::shared_ptr<LumiCalHit>;
using VecLCIOCalHit = std::set<EVENT::CalorimeterHit*>;

class ProjectionInfo {

public:
  ProjectionInfo();
  ProjectionInfo(CalHit const& calHit, int cellIdHitZ);

  void addHit(CalHit const& calHit );

  const double* getPosition() const { return position; }
  double getEnergy() const { return energy; }
  int getCellIdHitZ() const { return cellIdHitZ; }
  VecLCIOCalHit const& getCaloHits() const { return hits; }

private:
  double energy;
  double position[3];
  int cellIdHitZ;
  VecLCIOCalHit hits{};
public:
  bool newObject;
};


#endif // ProjectionInfo_hh
