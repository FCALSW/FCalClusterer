#ifndef BeamCalGeoDD_hh
#define BeamCalGeoDD_hh 1

#include "BeamCalGeo.hh"

#include <DD4hep/Detector.h>

class BeamCalGeoDD : public BeamCalGeo {     

public:
  BeamCalGeoDD( dd4hep::Detector const& theDetector, std::string const& detectorName="BeamCal",
                std::string const& colName="BeamCalCollection");
  BeamCalGeoDD();


  virtual int getPadsPerBeamCal() const;
  virtual int getPadsPerLayer() const;
  virtual int getLayer(int padIndex) const;
  virtual int getRing(int padIndex) const;
  virtual int getLocalPad(int padIndex) const;

  virtual void getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const;
  virtual int getPadIndex(int layer, int ring, int pad) const;

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
  virtual double                getPhiOffset()       const;

  virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const ;
  virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const;

  //we have 8 full segments
  virtual int getSymmetryFold() const { return m_symmetryFold; }

  virtual int getPadsInRing( int ring ) const;


private:
  dd4hep::DetElement m_BeamCal;
  dd4hep::Segmentation m_segmentation;

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
  double	        m_phiOffset;

  //hard coded
  int                   m_symmetryFold = 8;

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

  void readPolarGridRPhi2();
  void readPolarGridRPhi();
  
};



#endif // BeamCalGeoDD_hh
