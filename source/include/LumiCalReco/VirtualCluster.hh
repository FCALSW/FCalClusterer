#ifndef VirtualCluster_hh
#define VirtualCluster_hh 1


#include "GlobalMethodsClass.h"


class VirtualCluster
{

public:
  VirtualCluster();
  VirtualCluster(double x, double y, double z);

  ~VirtualCluster(){}


  void clear();

  inline double getX() const { return _x;}
  inline double getY() const { return _y;}
  inline double getZ() const { return _z;}

  inline void setX( double x ) { _x = x; }
  inline void setY( double y ) { _y = y; }
  inline void setZ( double z ) { _z = z; }

  friend std::ostream& operator<<(std::ostream & o, const VirtualCluster& rhs);

private:

  double _x, _y, _z;//Not sure the first parameter is Z or E, or Distance??




};


#endif // VirtualCluster_hh
