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
