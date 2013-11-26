#ifndef LCCluster_hh
#define LCCluster_hh 1


#include "GlobalMethodsClass.h"
#include <cmath>

class VirtualCluster;

class LCCluster
{

public:
  LCCluster();
  LCCluster( const VirtualCluster& vc );
  LCCluster( double energy, double x, double y, double z, double weight,
	     GlobalMethodsClass::WeightingMethod_t method, double theta, double phi );

  ~LCCluster(){}


  void clear();

  inline double getX() const { return _position[0];}
  inline double getY() const { return _position[1];}
  inline double getZ() const { return _position[2];}
  inline double getE() const { return _energy ;}
  inline double getWeight() const { return _weight; }
  inline GlobalMethodsClass::WeightingMethod_t getMethod() const { return _method; }
  inline double getEnergy() const { return getE(); }
  inline double getTheta() const { return _theta; }
  inline double getPhi() const { return _phi; }

  inline const double * getPosition() const { return _position; }

  inline void setX( double x) { _position[0] = x; CalculatePhi(); }
  inline void setY( double y) { _position[1] = y; CalculatePhi(); }
  inline void setWeight (double w) { _weight = w; }
  inline void setTheta  (double t) { _theta = t; }
  inline void setPhi    (double p) { _phi = p; CalculatePhi(); }

  void addToEnergy( double E) { _energy += E; }
  friend std::ostream& operator<<(std::ostream & o, const LCCluster& rhs);


private:

  void CalculatePhi();
  double _position[3];
  double _energy, _weight;
  GlobalMethodsClass::WeightingMethod_t _method;
  double _theta, _phi;

};


inline void LCCluster::CalculatePhi() {
  _phi = atan2(_position[1],_position[0]);
  return;
}


#endif // LCCluster_hh
