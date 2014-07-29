#ifndef VirtualCluster_hh
#define VirtualCluster_hh 1

#include <iosfwd>

class VirtualCluster
{

public:
  VirtualCluster();
  VirtualCluster(double x, double y, double z);

  ~VirtualCluster(){}


  void clear();

  inline double getX() const { return _position[0]; }
  inline double getY() const { return _position[1]; }
  inline double getZ() const { return _position[2]; }

  inline void setX( double x ) { _position[0] = x; }
  inline void setY( double y ) { _position[1] = y; }
  inline void setZ( double z ) { _position[2] = z; }

  inline const double* getPosition() const { return _position; }

  friend std::ostream& operator<<(std::ostream & o, const VirtualCluster& rhs);

private:

  double _position[3];//Not sure the first parameter is Z or E, or Distance??
  
};


#endif // VirtualCluster_hh
