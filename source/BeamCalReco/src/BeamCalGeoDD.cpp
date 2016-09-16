#include "BeamCalGeoDD.hh"

#include <DD4hep/LCDD.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/DetectorData.h>
#include <DDRec/API/IDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>


//#include <DD4hep/UserExtension/BeamCalInfo.h>
// Use the DDRec CalorimeterExtension, I guess?
// Use DDSegmentation RPphi parameters? This knows about tghe inner radius, too


class BeamCalInfo;

#include <vector>
#include <algorithm>

BeamCalGeoDD::BeamCalGeoDD(DD4hep::Geometry::LCDD const& lcdd): m_BeamCal(lcdd.detector("BeamCal")),
								m_segmentation(lcdd.readout("BeamCalCollection").segmentation()),
								m_innerRadius(0.0),
								m_outerRadius(0.0),
								m_layers(0),
								m_rings(0),
								m_phiSegmentation(0),
								m_radSegmentation(0),
								m_nPhiSegments(0),
								m_layerDistanceToIP(0),
								m_beamCalZPosition(0),
								m_deadAngle(0),
								m_padsPerRing(0),
								m_padsBeforeRing(0),
								m_padsPerLayer(-1),
								m_padsPerBeamCal(-1),
								m_crossingAngle(0),
								m_cutOut(0),
								m_firstFullRing(0)
{

  const DD4hep::DDRec::LayeredCalorimeterData * theExtension = m_BeamCal.extension<DD4hep::DDRec::LayeredCalorimeterData>();
  m_innerRadius = theExtension->extent[0]/dd4hep::mm;
  m_outerRadius = theExtension->extent[1]/dd4hep::mm;
  m_beamCalZPosition = theExtension->extent[2]/dd4hep::mm;
  const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer>& layers= theExtension->layers;

  m_layers = layers.size();
  for (size_t i = 0; i < layers.size(); i++) {
    const DD4hep::DDRec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);
    //no inner thickness means no sensitive detector, like the old graphite layer
    if ( theLayer.inner_thickness == 0.0 ) {
      m_layers--;
      continue;
    }
    //Distance to center of sensitive element
    m_layerDistanceToIP.push_back((theLayer.distance+theLayer.inner_thickness)/dd4hep::mm);
  }

  streamlog_out(MESSAGE) << "Segmentation Type" << m_segmentation.type()  << std::endl;
  streamlog_out(MESSAGE) <<"FieldDef: " << m_segmentation.segmentation()->fieldDescription()  << std::endl;

  if (m_segmentation.type() != "PolarGridRPhi2" ){
    streamlog_out(ERROR) << "Cannot use this segmentation type" 
			 << "(" << m_segmentation.type() << ")"
			 <<  "for BeamCalReco at the moment"  << std::endl;
    throw std::runtime_error( "BeamCalReco: Incompatible segmentation type" );
  }

  typedef DD4hep::DDSegmentation::TypedSegmentationParameter< std::vector<double> > ParVec;
  ParVec* rPar = dynamic_cast<ParVec*>(m_segmentation.segmentation()->parameter("grid_r_values"));
  ParVec* pPar = dynamic_cast<ParVec*>(m_segmentation.segmentation()->parameter("grid_phi_values"));
  std::vector<double> rValues = rPar->typedValue();
  std::vector<double> pValues = pPar->typedValue();

  m_rings = rValues.size()-1;
  m_radSegmentation = rValues;
  m_phiSegmentation = pValues;

  //DeadAngle Calculations
  typedef DD4hep::DDSegmentation::TypedSegmentationParameter< double > ParDou;
  ParDou* oPPar = dynamic_cast<ParDou*>(m_segmentation.segmentation()->parameter("offset_phi"));
  double offsetPhi = oPPar->typedValue();
  m_deadAngle = 2.0 * ( M_PI + offsetPhi/dd4hep::radian );

  for (int i = 0; i < m_rings;++i) {
    const int NUMBER_OF_SENSOR_SEGMENTS = 8;
    m_nPhiSegments.push_back(int((2*M_PI-m_deadAngle)/m_phiSegmentation[i]+0.5)/NUMBER_OF_SENSOR_SEGMENTS);
    m_radSegmentation[i] = m_radSegmentation[i]/dd4hep::mm;

  }
  // the outer radius needs to be fixed as well
  m_radSegmentation[m_rings] = m_radSegmentation[m_rings]/dd4hep::mm;


  m_padsPerRing.resize(m_rings+1);
  m_padsBeforeRing.resize(m_rings+1);

  setCrossingAngle();

  setCutOut();
  setFirstFullRing();

  setPadsInRing();
  setPadsBeforeRing();
  setPadsPerLayer();
  setPadsPerBeamCal();

}

//Wrappers around DD4hep Interface:
inline double                BeamCalGeoDD::getBCInnerRadius() const { 
  return m_innerRadius;
}

inline double                BeamCalGeoDD::getBCOuterRadius() const { 
  return m_outerRadius;
}

inline int                   BeamCalGeoDD::getBCLayers()      const { 
  return m_layers;
} 

inline int                   BeamCalGeoDD::getBCRings()       const { 
  return m_rings;
}

inline std::vector<int> const& BeamCalGeoDD::getNSegments()     const {
  return m_nPhiSegments;
}

inline double                BeamCalGeoDD::getCutout()        const {
  return m_cutOut;
}

inline double                BeamCalGeoDD::getBCZDistanceToIP() const { 
  return m_beamCalZPosition;
}


inline double BeamCalGeoDD::getLayerZDistanceToIP(int layer) const {
  return m_layerDistanceToIP.at(layer);
}

inline std::vector<double> const& BeamCalGeoDD::getPhiSegmentation() const {
  return m_phiSegmentation;
}

inline int BeamCalGeoDD::getPadsBeforeRing(int ring) const {
  return m_padsBeforeRing.at(ring);
}

void BeamCalGeoDD::setPadsBeforeRing() {
  int nPads = 0;
  for (int ring = 0; ring <= m_rings; ++ring) {//we want the number of pads _before_ the ring
    m_padsBeforeRing[ring] = nPads;
    nPads += getPadsInRing(ring);
  }
  return;
}

void BeamCalGeoDD::setPadsPerLayer() {
  m_padsPerLayer = getPadsBeforeRing( getBCRings() );
}

void BeamCalGeoDD::setPadsPerBeamCal() {
  m_padsPerBeamCal = getPadsPerLayer() * getBCLayers();
}

int BeamCalGeoDD::getPadsInRing( int ring ) const {
  return m_padsPerRing[ring];
}

void BeamCalGeoDD::setPadsInRing()  {
  for (int ring = 0; ring < m_rings; ++ring) {
    m_padsPerRing[ring] = BeamCalGeo::getPadsInRing(ring);
  }//for all rings
}//setPadsRings


inline double BeamCalGeoDD::getFullKeyHoleCutoutAngle() const {
  return getDeadAngle();
}

inline double BeamCalGeoDD::getDeadAngle() const {
  return m_deadAngle;
}

inline double BeamCalGeoDD::getCrossingAngle() const {
  return m_crossingAngle;
}

inline std::vector<double> const& BeamCalGeoDD::getRadSegmentation() const {
  return m_radSegmentation;
}

int BeamCalGeoDD::getFirstFullRing() const {
  return m_firstFullRing;
}

int BeamCalGeoDD::getPadsPerBeamCal() const{
  return m_padsPerBeamCal;
}

int BeamCalGeoDD::getPadsPerLayer() const {
  return m_padsPerLayer;
}

int BeamCalGeoDD::getRing(int padIndex) const {
  std::vector<int>::const_iterator element =  std::upper_bound(m_padsBeforeRing.begin()+1, m_padsBeforeRing.end(),
							       padIndex % m_padsPerLayer);
  return element - m_padsBeforeRing.begin() - 1 ;
}

int BeamCalGeoDD::getLocalPad(int padIndex) const {
  int layer, ring, pad;
  getLayerRingPad(padIndex, layer, ring, pad);
  return pad;
}

int BeamCalGeoDD::getLayer(int padIndex) const {
  //layer starts at 1
  return   ( padIndex / m_padsPerLayer ) + 1;
}

void BeamCalGeoDD::setFirstFullRing() {
  int ring = 0;
  while ( m_cutOut > getRadSegmentation()[ring] ) { ++ring; }
  m_firstFullRing = ring;
}

void BeamCalGeoDD::setCutOut() {
  const double cutOutRadius =  tan(m_crossingAngle/1000.0)*getLayerZDistanceToIP(m_layers-1)+3.5;
  m_cutOut = cutOutRadius;
}

/// get the full crossing angle from the position of the BeamCal
void BeamCalGeoDD::setCrossingAngle() {
  DD4hep::Geometry::DetElement::Children children = m_BeamCal.children();
  for ( DD4hep::Geometry::DetElement::Children::const_iterator it = children.begin();
	it != children.end();
	++it ) {
    DD4hep::Geometry::Position loc(0.0, 0.0, 0.0);
    DD4hep::Geometry::Position glob(0.0, 0.0, 0.0);
    it->second.localToWorld( loc, glob );
    m_crossingAngle = 2.0*fabs( atan( glob.x() / glob.z() ) / dd4hep::mrad );
    return;
  }
  throw std::runtime_error( "Cannot obtain crossing angle from this BeamCal, update lcgeo?" );
}
