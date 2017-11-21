#include "LCCluster.hh"
#include "VirtualCluster.hh"

#include <iostream>
#include <iomanip>

LCCluster::LCCluster()
    : _position{0.0, 0.0, 0.0},
      _energy(0.0),
      _weight(0.0),
      _method(GlobalMethodsClass::LogMethod),
      _theta(0.0),
      _phi(0.0),
      _caloHits{} {}

LCCluster::LCCluster(const VirtualCluster& vc)
    : _position{vc.getX(), vc.getY(), vc.getZ()},
      _energy(0.0),
      _weight(0.0),
      _method(GlobalMethodsClass::LogMethod),
      _theta(0.0),
      _phi(0.0),
      _caloHits{} {}

LCCluster::LCCluster(double energy, double x, double y, double z, double weight,
                     GlobalMethodsClass::WeightingMethod_t method, double theta, double phi, VecCalHit const& caloHitVector)
    : _position{x, y, z},
      _energy(energy),
      _weight(weight),
      _method(method),
      _theta(theta),
      _phi(phi),
      _caloHits(caloHitVector) {
  CalculatePhi();
}

void LCCluster::clear() {
  _energy = 0.0;
  _position[0] = 0.0;
  _position[1] = 0.0;
  _position[2] = 0.0;
  _weight = 0.0;
  _method      = GlobalMethodsClass::LogMethod;
  _theta = 0.0;
  _phi = 0.0;
  _caloHits.clear();
}

std::ostream& operator<<(std::ostream & o, const LCCluster& rhs) {
  o << "  Energy "              << std::setw(10) << rhs._energy
    << "  Method "              << std::setw(4)  << rhs._method
    << "  Weight "              << std::setw(10) << rhs._weight
    << "  N Calo Hits "         << std::setw(7) <<  rhs._caloHits.size()
    << "  pos(x,y,z) =  ( "
    << std::setw(10) << rhs._position[0] << " , "
    << std::setw(10) << rhs._position[1] << " , "
    << std::setw(10) << rhs._position[2] << " )"
    << "  pos(theta,phi) =  ( " << std::setw(10) << rhs._theta << " , " << std::setw(10) << rhs._phi << " )";
  return o;
}

/** recalculate the position of the cluster based on the theta and phi averages
 *
 * Resolution in Theta (R) is better than in RPhi so averaging theta gives better results
 */
void LCCluster::recalculatePositionFromHits(GlobalMethodsClass const& gmc) {
  const double logConstant(gmc.GlobalParamD.at(GlobalMethodsClass::LogWeightConstant));

  //re-set new cluster energy
  _energy = 0.0;
  for (auto const& calHit : _caloHits) {
    _energy += calHit->getEnergy();
  }

  const double rMin = gmc.GlobalParamD.at(GlobalMethodsClass::RMin);

  double thetaTemp(0.0), weightsTemp(0.0), xTemp(0.0), yTemp(0.0), zTemp(0.0);
  for (auto const& calHit : _caloHits) {
    const auto*  pos    = calHit->getPosition();
    const double rCell  = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    //cell area scales with radius, reduce weight for cells at larger radii
    const double weight = GlobalMethodsClass::posWeight(calHit->getEnergy(), _energy, _method, logConstant) * rMin / rCell;

    if (not(weight > 0))
      continue;

    const double theta = atan(rCell / pos[2]);

    thetaTemp += theta * weight;

    xTemp += pos[0] * weight;
    yTemp += pos[1] * weight;
    zTemp += pos[2] * weight;

    weightsTemp += weight;
  }

  xTemp /= weightsTemp;
  yTemp /= weightsTemp;
  zTemp /= weightsTemp;

  thetaTemp /= weightsTemp;

  _phi   = atan2(yTemp, xTemp);
  _theta = thetaTemp;

  const double r = sqrt(xTemp * xTemp + yTemp * yTemp + zTemp * zTemp);

  _position[0] = r * sin(_theta) * cos(_phi);
  _position[1] = r * sin(_theta) * sin(_phi);
  _position[2] = r * cos(_theta);
}
