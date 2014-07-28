#include "LumiCalClusterer.h"
#include "SuperTrueClusterWeights.hh"

#include <cmath>
#include <cassert>


SuperTrueClusterWeights::SuperTrueClusterWeights(	int superClusterIdNow,
							int trueClusterIdNow,
							LCCluster superClusterCM,
							LCCluster trueClusterCM):

  superClusterId(superClusterIdNow),
  trueClusterId(trueClusterIdNow),
  distance(LumiCalClustererClass::distance2D(superClusterCM.getPosition(),trueClusterCM.getPosition())),
  deltaEngy(fabs(superClusterCM.getE() - trueClusterCM.getE())),
  minEngy(std::min(superClusterCM.getE() , trueClusterCM.getE())),
  weight(-1)
{
}

void SuperTrueClusterWeights::setWeight(std::string weightMethod) {

  if(weightMethod == "distance")	weight = distance;

  if(weightMethod == "deltaEngy") weight = deltaEngy;
}

void SuperTrueClusterWeights::setWeight(std::string weightMethod, double minSeparationDistance, double minClusterEngyGeV) {

  if(weightMethod == "minEngyDistance") {
    if(distance > minSeparationDistance && minEngy > minClusterEngyGeV)
      weight = 1.;
    else
      weight = -1.;
  }
}
