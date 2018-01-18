#include "BeamCalGeoCached.hh"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearMgr.h>


#include <vector>
#include <algorithm>
#include <cmath>


BeamCalGeoCached::BeamCalGeoCached(gear::GearMgr* gearMgr): m_BCPs(gearMgr->getBeamCalParameters()),
							    m_innerRadius(m_BCPs.getExtent()[0]),
							    m_outerRadius(m_BCPs.getExtent()[1]),
							    m_layers(m_BCPs.getLayerLayout().getNLayers()),
							    m_rings(m_BCPs.getDoubleVals("phi_segmentation").size()),
							    m_phiSegmentation(m_BCPs.getDoubleVals("phi_segmentation")),
							    m_radSegmentation(0),
			                                    m_nPhiSegments(0),
                                                            m_cutOut(m_BCPs.getDoubleVal("dead_area_outer_r")+0.1),
  m_firstFullRing(0),
  m_beamCalZPosition(m_BCPs.getExtent()[2]+100.0), //add distance from graphite shield
  m_deadAngle(2.0 * ( m_BCPs.getDoubleVal("cylinder_starting_phi") - M_PI )),
  m_phiOffset(180.0 / M_PI * 0.5 * m_deadAngle),
  m_crossingAngle(m_BCPs.getDoubleVal("beam_crossing_angle")),
  m_padsPerRing(m_rings+1),
  m_padsBeforeRing(m_rings+1),
  m_padsPerLayer(-1),
  m_padsPerBeamCal(-1)
{

  try {
    m_radSegmentation = m_BCPs.getDoubleVals("rad_segmentation");
  } catch (gear::UnknownParameterException &e) {
    const double step = (m_outerRadius - m_innerRadius) / double(m_rings);
    for (int i = 0; i < m_rings ;++i) {
      m_radSegmentation.push_back(m_innerRadius+i*step);
    }
  }

  try {
    m_nPhiSegments = m_BCPs.getIntVals("nPhi_segmentation");
  } catch (gear::UnknownParameterException &e) {
    for (int i = 0; i < m_rings;++i) {
      m_nPhiSegments.push_back(int(2*M_PI/m_phiSegmentation[i]+0.5)/9);
    }
  }

  if ( m_radSegmentation.size()-1 != m_phiSegmentation.size()) {
    m_radSegmentation.push_back(m_outerRadius);
  }

  setFirstFullRing();

  setPadsInRing();
  setPadsBeforeRing();
  setPadsPerLayer();
  setPadsPerBeamCal();

}

//Wrappers around Gear Interface:
inline double                BeamCalGeoCached::getBCInnerRadius() const {
  return m_innerRadius;
}

inline double                BeamCalGeoCached::getBCOuterRadius() const {
  return m_outerRadius;
}

inline int                   BeamCalGeoCached::getBCLayers()      const {
  return m_layers;
}

inline int                   BeamCalGeoCached::getBCRings()       const {
  return m_rings;
}

inline std::vector<double> const& BeamCalGeoCached::getPhiSegmentation()  const {
  return m_phiSegmentation;
}

inline std::vector<double> const& BeamCalGeoCached::getRadSegmentation()  const {
  return m_radSegmentation;
}

inline std::vector<int> const& BeamCalGeoCached::getNSegments()     const {
  return m_nPhiSegments;
}

inline double                BeamCalGeoCached::getCutout()        const {
  return m_cutOut;
}

inline double                BeamCalGeoCached::getBCZDistanceToIP() const {
  return m_beamCalZPosition;
}


inline double                BeamCalGeoCached::getLayerZDistanceToIP(int layer) const {
#pragma message "FIXME: make thickness of graphite shield to be read from GEAR"
  const double graphiteShield_dZ = 100.;
  double lr_zdist(graphiteShield_dZ);
  lr_zdist+=m_BCPs.getLayerLayout().getDistance(layer);
  return lr_zdist;
}


int BeamCalGeoCached::getPadsInRing( int ring ) const {
  return m_padsPerRing[ring];
}


void BeamCalGeoCached::setPadsInRing()  {
  for (int ring = 0; ring < m_rings; ++ring) {
    m_padsPerRing[ring] = BeamCalGeo::getPadsInRing(ring);
  }//for all rings
}//setPadsRings

inline int BeamCalGeoCached::getPadsBeforeRing( int ring ) const {
  //  std::cout << "Ring " << ring << " val " << m_padsBeforeRing[ring]  << std::endl;
  return m_padsBeforeRing[ring];
}


void BeamCalGeoCached::setPadsBeforeRing() {
  int nPads = 0;
  for (int ring = 0; ring <= m_rings; ++ring) {//we want the number of pads _before_ the ring
    m_padsBeforeRing[ring] = nPads;
    nPads += getPadsInRing(ring);
  }
  return;
}


inline int BeamCalGeoCached::getPadsPerLayer() const {
  return m_padsPerLayer;
}


inline int BeamCalGeoCached::getPadsPerBeamCal() const {
  return m_padsPerBeamCal;
}


void BeamCalGeoCached::setPadsPerLayer() {
  m_padsPerLayer = getPadsBeforeRing( getBCRings() );
}


void BeamCalGeoCached::setPadsPerBeamCal() {
  m_padsPerBeamCal = getPadsPerLayer() * getBCLayers();
}


double BeamCalGeoCached::getDeadAngle() const {
  return m_deadAngle;
}


 double BeamCalGeoCached::getFullKeyHoleCutoutAngle() const {
  return getDeadAngle();
}

 double BeamCalGeoCached::getPhiOffset() const { return m_phiOffset; }

 double BeamCalGeoCached::getCrossingAngle() const {
  return m_crossingAngle;
}

int BeamCalGeoCached::getFirstFullRing() const {
  return m_firstFullRing;
}

void BeamCalGeoCached::setFirstFullRing() {
  int ring = 0;
  while ( m_cutOut > getRadSegmentation()[ring] ) { ++ring; }
  m_firstFullRing = ring;
}
