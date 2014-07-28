#ifndef BeamCalGeoDD_hh
#define BeamCalGeoDD_hh 1

#include "BeamCalGeo.hh"

namespace DD4hep{
  namespace Geometry {
    struct DetElement;
    struct LCDD;
  }
}

class BeamCalGeoDD : public BeamCalGeo {     

public:
  BeamCalGeoDD( DD4hep::Geometry::LCDD const& lcdd );
  BeamCalGeoDD();


  virtual int getPadsPerBeamCal() const;
  virtual int getPadsPerLayer() const;
  virtual int getLayer(int padIndex) const;
  virtual int getRing(int padIndex) const;
  virtual int getLocalPad(int padIndex) const;
  virtual void getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const;
  virtual bool arePadsNeighbours(int padIndex1, int padIndex2, bool mustBeInSameLayer = false) const;

 
  virtual void countNumberOfPadsInRing();
  virtual int getPadIndex(int layer, int ring, int pad) const throw(std::out_of_range);
  virtual double getPadPhi(int ring, int pad) const;
  virtual double getPadPhi(int globalPandIndex) const;

  virtual double getThetaFromRing(double averageRing) const;


  virtual double                getBCInnerRadius()   const;
  virtual double                getBCOuterRadius()   const;
  virtual int                   getBCLayers()        const;
  virtual int                   getBCRings()         const;
  virtual std::vector<double>   getSegmentation()    const;
  virtual std::vector<int>      getNSegments()       const;
  virtual double                getCutout()          const;
  virtual double                getBCZDistanceToIP() const;

private:
  const DD4hep::Geometry::DetElement& m_BeamCal;

};



#endif // BeamCalGeoDD_hh
