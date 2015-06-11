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
  m_beamCalZPosition(m_BCPs.getExtent()[2]+100.0), //add distance from graphite shield
  m_deadAngle(2.0 * ( m_BCPs.getDoubleVal("cylinder_starting_phi") - M_PI )),
  m_crossingAngle(m_BCPs.getDoubleVal("beam_crossing_angle")),
  m_padsPerRing(m_rings),
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

inline double                BeamCalGeoCached::getLayerZDistanceToIP(int layer) const {
  // static const double globalBeamCalDist = 3350; return globalBeamCalDist;
#pragma message "FIXME: make thickness of graphite shield to be read from GEAR"
  const double graphiteShield_dZ = 100.;
  double lr_zdist(graphiteShield_dZ);
  lr_zdist+=m_BCPs.getLayerLayout().getDistance(layer);
  return lr_zdist;
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

 double BeamCalGeoCached::getCrossingAngle() const {
  return m_crossingAngle;
}

void BeamCalGeoCached::getPadExtentsById(int globalPadIndex, double *extents) const
{
  // pad index within layer
  int padIndex = globalPadIndex%this->getPadsPerLayer();

  std::vector<int>::const_iterator it_cylinder = 
    std::upper_bound(m_padsBeforeRing.begin(), m_padsBeforeRing.end(), padIndex)-1;
  int cylinder = int(it_cylinder - m_padsBeforeRing.begin());
  int sector = padIndex - m_padsBeforeRing.at(cylinder);

  this->getPadExtents(cylinder, sector, extents);
  //std::cout << m_padsBeforeRing.at(cylinder) << "\t" << cylinder << "\t" << sector <<"\t" << extents[4] << "\t" << extents[5]<< std::endl;
}


double BeamCalGeoCached::getPadsDistance(int padIndex1, int padIndex2) const
{
  const double DEGRAD=M_PI/180.;
  static double e1[6], e2[6];
  static int prev_pi1(-1), prev_pi2(-1);
  if ( padIndex1 != prev_pi1) getPadExtentsById(padIndex1, e1), prev_pi1 = padIndex1;
  if ( padIndex2 != prev_pi2) getPadExtentsById(padIndex2, e2), prev_pi2 = padIndex2;
  //std::cout << padIndex1 << "\t" << padIndex2 << "\t" << e2[4] << "\t" << e2[5] << std::endl;

  double dx = e1[4]*cos(e1[5]*DEGRAD) - e2[4]*cos(e2[5]*DEGRAD);
  double dy = e1[4]*sin(e1[5]*DEGRAD) - e2[4]*sin(e2[5]*DEGRAD);
  double d = sqrt(dx*dx+dy*dy);
  return d;
}
