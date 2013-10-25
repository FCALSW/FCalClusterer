#include "BeamCalGeoGear.hh"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearMgr.h>


#include <vector>
#include <algorithm>
#include <cmath>

BeamCalGeoGear::BeamCalGeoGear(gear::GearMgr* gearMgr): m_BCPs(gearMgr->getBeamCalParameters()) {
}

//Wrappers around Gear Interface:
inline double                BeamCalGeoGear::getBCInnerRadius() const { 
  return m_BCPs.getExtent()[0]; 
}

inline double                BeamCalGeoGear::getBCOuterRadius() const { 
  return m_BCPs.getExtent()[1]; 
}

inline int                   BeamCalGeoGear::getBCLayers()      const { 
  return m_BCPs.getLayerLayout().getNLayers(); 
} 

inline int                   BeamCalGeoGear::getBCRings()       const { 
  return m_BCPs.getDoubleVals("phi_segmentation").size(); 
}

inline std::vector<double> const&  BeamCalGeoGear::getPhiSegmentation()  const { 
  return m_BCPs.getDoubleVals("phi_segmentation"); 
}

inline std::vector<double> const&  BeamCalGeoGear::getRadSegmentation()  const { 
  return m_BCPs.getDoubleVals("rad_segmentation"); 
}


inline std::vector<int> const&     BeamCalGeoGear::getNSegments()     const {
  return m_BCPs.getIntVals("nPhi_segmentation"); 
}

inline double                BeamCalGeoGear::getCutout()        const {
  return m_BCPs.getDoubleVal("dead_area_outer_r")+0.1; 
}

inline double                BeamCalGeoGear::getBCZDistanceToIP() const { 
  // static const double globalBeamCalDist = 3350; return globalBeamCalDist; 
  return m_BCPs.getExtent()[2];
}


// void BeamCalGeoGear::countNumberOfPadsInRing(){
//   //  m_PadsBeforeRing.clear();

//   for (std::vector<int>::iterator it = getNSegments().begin(); it != getNSegments().end(); ++it) {
//     m_PadsBeforeRing.push_back(getPadsPerLayer());//first ring has 0 entries! No it doesn't, yes it does have 0 BEFORE it
//     getPadsPerLayer()+= *it;
//   }

//   //  if(getNSegments().size() != m_PadsBeforeRing.size()) throw std::logic_error("this should not have happened");

//   m_nPadsPerBeamCal = getBCLayers() * getPadsPerLayer();
//   //m_PadEnergies.resize(m_nPadsPerBeamCal, 0.0);
// }


 double BeamCalGeoGear::getFullKeyHoleCutoutAngle() const {
  return getDeadAngle();
}

double BeamCalGeoGear::getDeadAngle() const {
  return 2.0 * ( m_BCPs.getDoubleVal("cylinder_starting_phi") - M_PI );
}


 double BeamCalGeoGear::getCrossingAngle() const {
  return m_BCPs.getDoubleVal("beam_crossing_angle");
}
