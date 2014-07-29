#include "ProjectionInfo.hh"

//LCIO
#include <IMPL/CalorimeterHitImpl.h>


ProjectionInfo::ProjectionInfo (): energy (0.0), cellIdHitZ(0), newObject(true) {
  position[0] = 0.0;
  position[1] = 0.0;
  position[2] = 0.0;
}


ProjectionInfo::ProjectionInfo ( IMPL::CalorimeterHitImpl const* calHit, int cellIdZ):
  energy ( calHit->getEnergy() ), cellIdHitZ(cellIdZ), newObject(false){
  position[0] = calHit->getPosition()[0];
  position[1] = calHit->getPosition()[1];
  position[2] = calHit->getPosition()[2];
}
