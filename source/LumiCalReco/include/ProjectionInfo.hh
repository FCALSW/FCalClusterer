#ifndef ProjectionInfo_hh
#define ProjectionInfo_hh 1

namespace IMPL {
  class CalorimeterHitImpl;
}

class ProjectionInfo {

public:
  ProjectionInfo ();
  ProjectionInfo ( IMPL::CalorimeterHitImpl const* calHit, int cellIdHitZ);

  const double* getPosition() const { return position; }

  double energy;
  double position[3];
  int cellIdHitZ;
  bool newObject;

};


#endif // ProjectionInfo_hh
