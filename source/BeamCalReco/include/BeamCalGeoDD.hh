#ifndef BeamCalGeoDD_hh
#define BeamCalGeoDD_hh 1

#include "BeamCalGeo.hh"

#include <DD4hep/LCDD.h>
#include <DD4hep/Detector.h>

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

  virtual double                getBCInnerRadius()   const;
  virtual double                getBCOuterRadius()   const;
  virtual int                   getBCLayers()        const;
  virtual int                   getBCRings()         const;
  virtual std::vector<double> const&  getPhiSegmentation() const;
  virtual std::vector<double> const&  getRadSegmentation() const;
  virtual std::vector<int> const& getNSegments()       const;
  virtual double                getCutout()          const;
  virtual double                getBCZDistanceToIP() const;

  virtual double                getLayerZDistanceToIP(const int lr) const;
  virtual double                getDeadAngle()       const;

  virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const ;
  virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const;

  //we have 8 full segments
  virtual int getSymmetryFold() const { return 8; }

  virtual int getPadsInRing( int ring ) const;


private:
  DD4hep::Geometry::DetElement m_BeamCal;
  DD4hep::Geometry::Segmentation m_segmentation;

  ///Read From DD4hep/DDRec/Segmentation
  double		m_innerRadius;
  double		m_outerRadius;
  int			m_layers;
  int			m_rings;
  std::vector<double>	m_phiSegmentation;
  std::vector<double>	m_radSegmentation;
  std::vector<int>	m_nPhiSegments;
  std::vector<double>   m_layerDistanceToIP;
  double		m_beamCalZPosition;
  double                m_deadAngle;

  //calculated
  std::vector<int>	m_padsPerRing;
  std::vector<int>	m_padsBeforeRing;
  int			m_padsPerLayer;
  int			m_padsPerBeamCal;
  double                m_crossingAngle;
  double		m_cutOut;
  int                   m_firstFullRing;

  void setPadsInRing();
  void setPadsBeforeRing();
  void setPadsPerLayer();
  void setPadsPerBeamCal();
  void setCrossingAngle();
  void setFirstFullRing();
  void setCutOut();

};



#endif // BeamCalGeoDD_hh
