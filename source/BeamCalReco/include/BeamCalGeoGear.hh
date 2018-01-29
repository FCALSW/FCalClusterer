#ifndef BeamCalGeoGear_hh
#define BeamCalGeoGear_hh 1

#include "BeamCalGeo.hh"

namespace gear{
  class CalorimeterParameters;
  class GearMgr;
}

#include <vector>

class BeamCalGeoGear : public BeamCalGeo {     

public:
  BeamCalGeoGear(gear::GearMgr* gearMgr);

  virtual ~BeamCalGeoGear() {}

  // virtual int getPadsPerBeamCal() const;
  // virtual int getPadsPerLayer() const;
  // virtual int getLayer(int padIndex) const;
  // virtual int getRing(int padIndex) const;
  // virtual int getLocalPad(int padIndex) const;
  // virtual void getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const;
 
  //  virtual void countNumberOfPadsInRing();
  //virtual int getPadIndex(int layer, int ring, int pad) const throw(std::out_of_range);

  // virtual double getPadPhi(int ring, int pad) const;
  // virtual double getPadPhi(int globalPandIndex) const;
  // virtual double getThetaFromRing(double averageRing) const;

  virtual double                getBCInnerRadius()   const;
  virtual double                getBCOuterRadius()   const;
  virtual int                   getBCLayers()        const;
  virtual int                   getBCRings()         const;
  virtual std::vector<double> const&  getPhiSegmentation() const;
  virtual std::vector<double> const&  getRadSegmentation() const;
  virtual std::vector<int>    const&  getNSegments()       const;
  virtual double                getCutout()          const;
  virtual double                getBCZDistanceToIP() const;
  virtual double                getDeadAngle()       const;

  //virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const ;
  // virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const;


  //virtual int getPadsInRing( int ring ) const;
  //we have 8 full segments
  virtual int getSymmetryFold() const { return 8; }


private:
  const gear::CalorimeterParameters& m_BCPs;





};



#endif // BeamCalGeoGear_hh
