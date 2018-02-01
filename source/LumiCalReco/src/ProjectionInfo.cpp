#include "ProjectionInfo.hh"
#include "LumiCalHit.hh"

#include <memory>

ProjectionInfo::ProjectionInfo (): energy (0.0), cellIdHitZ(0), newObject(true) {
  position[0] = 0.0;
  position[1] = 0.0;
  position[2] = 0.0;
}

ProjectionInfo::ProjectionInfo(CalHit const& calHit, int cellIdZ)
    : energy(calHit->getEnergy()), cellIdHitZ(cellIdZ), newObject(false) {
  position[0] = calHit->getPosition()[0];
  position[1] = calHit->getPosition()[1];
  position[2] = calHit->getPosition()[2];
  hits.insert(calHit->beginHits(), calHit->endHits());
}

void ProjectionInfo::addHit(CalHit const& calHit) {
  hits.insert(calHit->beginHits(), calHit->endHits());
  energy += calHit->getEnergy();
}
