#include "Distance2D.hh"
#include "SuperTrueClusterWeights.hh"

#include <cmath>
#include <cassert>


SuperTrueClusterWeights::SuperTrueClusterWeights(	int superClusterIdNow,
							int trueClusterIdNow,
							LCCluster superClusterCM,
							LCCluster trueClusterCM):

  superClusterId(superClusterIdNow),
  trueClusterId(trueClusterIdNow),
  distance(LCHelper::distance2D(superClusterCM.getPosition(),trueClusterCM.getPosition())),
  deltaEngy(fabs(superClusterCM.getE() - trueClusterCM.getE())),
  minEngy(std::min(superClusterCM.getE() , trueClusterCM.getE())),
  weight(-1)
{
}

void SuperTrueClusterWeights::setWeight(std::string weightMethod) {

  weight  = ( weightMethod == "distance" ) ? distance : deltaEngy ;

}

void SuperTrueClusterWeights::setWeight(std::string weightMethod, double minSeparationDistance, double minClusterEngyGeV) {

  if(weightMethod == "minEngyDistance") {

    weight = (distance > minSeparationDistance && minEngy > minClusterEngyGeV) ?  1. : -1. ;

  }
}
