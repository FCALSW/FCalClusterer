#include "VirtualCluster.hh"

#include <iostream>
#include <iomanip>


VirtualCluster::VirtualCluster():
  _x(0.0),
  _y(0.0),
  _z(0.0)
{}

VirtualCluster::VirtualCluster(double x, double y, double z):
  _x(x),
  _y(y),
  _z(z)
{}

void VirtualCluster::clear() {

  _x = 0.0;
  _y = 0.0;
  _z = 0.0;

}

std::ostream& operator<<(std::ostream & o, const VirtualCluster& rhs) {
  o << std::setw(13) << rhs._x
    << std::setw(13) << rhs._y
    << std::setw(13) << rhs._z;
  return o;
}
