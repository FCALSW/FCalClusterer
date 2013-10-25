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
							    m_radSegmentation(m_BCPs.getDoubleVals("rad_segmentation")),
			                                    m_nPhiSegments(m_BCPs.getIntVals("nPhi_segmentation")),
                                                            m_cutOut(m_BCPs.getDoubleVal("dead_area_outer_r")+0.1),
                                                            m_beamCalZPosition(m_BCPs.getExtent()[2]),
  m_deadAngle(2.0 * ( m_BCPs.getDoubleVal("cylinder_starting_phi") - M_PI )),
  m_crossingAngle(m_BCPs.getDoubleVal("beam_crossing_angle")),
  m_padsPerRing(m_rings),
  m_padsBeforeRing(m_rings+1),
  m_padsPerLayer(-1),
  m_padsPerBeamCal(-1)
{
  if ( m_radSegmentation.size()-1 != m_phiSegmentation.size()) {
    m_radSegmentation.push_back(m_outerRadius);
  }
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
  // static const double globalBeamCalDist = 3350; return globalBeamCalDist;
  return m_beamCalZPosition;
}


// void BeamCalGeoCached::countNumberOfPadsInRing(){
//   //  m_PadsBeforeRing.clear();

//   for (std::vector<int>::iterator it = getNSegments().begin(); it != getNSegments().end(); ++it) {
//     m_PadsBeforeRing.push_back(getPadsPerLayer());//first ring has 0 entries! No it doesn't, yes it does have 0 BEFORE it
//     getPadsPerLayer()+= *it;
//   }

//   //  if(getNSegments().size() != m_PadsBeforeRing.size()) throw std::logic_error("this should not have happened");

//   m_nPadsPerBeamCal = getBCLayers() * getPadsPerLayer();
//   //m_PadEnergies.resize(m_nPadsPerBeamCal, 0.0);
// }


int BeamCalGeoCached::getPadsInRing( int ring ) const {
  return m_padsPerRing[ring];
}


void BeamCalGeoCached::setPadsInRing()  {
  for (int ring = 0; ring < m_rings; ++ring) {
    if ( ring < getFirstFullRing() ) {
      m_padsPerRing[ring] = getNSegments()[ring] * getSymmetryFold();
    } else {
      int segmentsInRing = getSymmetryFold()*getNSegments()[ring];
      double deltaPhi = ( 2.0*M_PI - getDeadAngle())/double(segmentsInRing);
      const int additionalSegments = (int)(getDeadAngle()/deltaPhi);
      m_padsPerRing[ring] = segmentsInRing + additionalSegments;
    }
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

 double BeamCalGeoCached::getCrossingAngle() const {
  return m_crossingAngle;
}
