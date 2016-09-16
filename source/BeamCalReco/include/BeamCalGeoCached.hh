#ifndef BeamCalGeoCached_hh
#define BeamCalGeoCached_hh 1

#include "BeamCalGeo.hh"

namespace gear{
  class CalorimeterParameters;
  class GearMgr;
}

class BeamCalGeoCached : public BeamCalGeo {     

public:
  BeamCalGeoCached(gear::GearMgr* gearMgr);

  virtual ~BeamCalGeoCached() {}

  virtual int getPadsPerBeamCal() const;
  virtual int getPadsPerLayer() const;

  virtual double                getBCInnerRadius()   const;
  virtual double                getBCOuterRadius()   const;
  virtual int                   getBCLayers()        const;
  virtual int                   getBCRings()         const;
  virtual std::vector<double> const&  getPhiSegmentation() const;
  virtual std::vector<double> const&  getRadSegmentation() const;
  virtual std::vector<int>    const&  getNSegments()       const;
  virtual double                getCutout()          const;
  virtual double                getBCZDistanceToIP() const;
  virtual double                getLayerZDistanceToIP(const int lr) const;
  virtual double                getDeadAngle()       const;

  virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const ;
  virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const;

  virtual int getPadsInRing( int ring ) const;
  //we have 8 full segments
  virtual int getSymmetryFold() const { return 8; }


private:
  const gear::CalorimeterParameters& m_BCPs;


  ///Read From GearFile
  double		m_innerRadius;
  double		m_outerRadius;
  int			m_layers;
  int			m_rings;
  std::vector<double>	m_phiSegmentation;
  std::vector<double>	m_radSegmentation;
  std::vector<int>	m_nPhiSegments;
  double		m_cutOut;
  int                   m_firstFullRing;
  double		m_beamCalZPosition;
  double                m_deadAngle;
  double                m_crossingAngle;
  //calculated
  std::vector<int>	m_padsPerRing;
  std::vector<int>	m_padsBeforeRing;
  int			m_padsPerLayer;
  int			m_padsPerBeamCal;

  void setPadsInRing();
  void setPadsBeforeRing();
  void setPadsPerLayer();
  void setPadsPerBeamCal();
  void setFirstFullRing();
};



#endif // BeamCalGeoCached_hh
