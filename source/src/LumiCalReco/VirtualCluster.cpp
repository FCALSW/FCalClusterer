#include "VirtualCluster.hh"

#include <iostream>
#include <iomanip>


VirtualCluster::VirtualCluster()
{
  _position[0] = 0.0;
  _position[1] = 0.0;
  _position[2] = 0.0;
}

VirtualCluster::VirtualCluster(double x, double y, double z)
{
  _position[0] = x;
  _position[1] = y;
  _position[2] = z;
}

void VirtualCluster::clear() {
  _position[0] = 0.0;
  _position[1] = 0.0;
  _position[2] = 0.0;
}

std::ostream& operator<<(std::ostream & o, const VirtualCluster& rhs) {
  o << std::setw(13) << rhs._position[0]
    << std::setw(13) << rhs._position[1]
    << std::setw(13) << rhs._position[2];
  return o;
}
