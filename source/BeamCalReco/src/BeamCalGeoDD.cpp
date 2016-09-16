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

  for (int i = 0; i < m_rings;++i) {
    const int NUMBER_OF_SENSOR_SEGMENTS = 8;
    m_nPhiSegments.push_back(int((2*M_PI-m_deadAngle)/m_phiSegmentation[i]+0.5)/NUMBER_OF_SENSOR_SEGMENTS);
    m_radSegmentation[i] = m_radSegmentation[i]/dd4hep::mm;

  }
  // the outer radius needs to be fixed as well
  m_radSegmentation[m_rings] = m_radSegmentation[m_rings]/dd4hep::mm;

  //DeadAngle Calculations
  typedef DD4hep::DDSegmentation::TypedSegmentationParameter< double > ParDou;
  ParDou* oPPar = dynamic_cast<ParDou*>(m_segmentation.segmentation()->parameter("offset_phi"));
  double offsetPhi = oPPar->typedValue();
  m_deadAngle = 2.0 * ( M_PI + offsetPhi/dd4hep::radian );


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

// [ VERBOSE "BCReco"] PadsBeforeRing 0 0
// [ VERBOSE "BCReco"] PadsBeforeRing 1 32
// [ VERBOSE "BCReco"] PadsBeforeRing 2 72
// [ VERBOSE "BCReco"] PadsBeforeRing 3 120
// [ VERBOSE "BCReco"] PadsBeforeRing 4 168
// [ VERBOSE "BCReco"] PadsBeforeRing 5 224
// [ VERBOSE "BCReco"] PadsBeforeRing 6 288

// [ VERBOSE "BCReco"] PadsBeforeRing 7 360
// [ VERBOSE "BCReco"] PadsBeforeRing 8 432
// [ VERBOSE "BCReco"] PadsBeforeRing 9 512
// [ VERBOSE "BCReco"] PadsBeforeRing 10 600
// [ VERBOSE "BCReco"] PadsBeforeRing 11 696
// [ VERBOSE "BCReco"] PadsBeforeRing 12 792
// [ VERBOSE "BCReco"] PadsBeforeRing 13 896
// [ VERBOSE "BCReco"] PadsBeforeRing 14 1008
// [ VERBOSE "BCReco"] PadsBeforeRing 15 1128

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


// void BeamCalGeoDD::countNumberOfPadsInRing(){
//   //  m_PadsBeforeRing.clear();

//   for (std::vector<int>::iterator it = m_Segments.begin(); it != m_Segments.end(); ++it) {
//     m_PadsBeforeRing.push_back(m_nPadsPerLayer);//first ring has 0 entries! No it doesn't, yes it does have 0 BEFORE it
//     m_nPadsPerLayer+= *it;
//   }

//   //  if(m_Segments.size() != m_PadsBeforeRing.size()) throw std::logic_error("this should not have happened");

//   m_nPadsPerBeamCal = m_nLayers * m_nPadsPerLayer;
//   //m_PadEnergies.resize(m_nPadsPerBeamCal, 0.0);
// }



int BeamCalGeoDD::getPadsPerBeamCal() const{
  return m_padsPerBeamCal;
}

int BeamCalGeoDD::getPadsPerLayer() const {
  return m_padsPerLayer;
}

//layers start at 1, ring and pad start at 0
int BeamCalGeoDD::getPadIndex(int layer, int ring, int pad) const throw(std::out_of_range){
  if( layer < 1 || m_layers < layer) {//starting at 1 ending at nLayers
    throw std::out_of_range("Layer out of range:");
  } else if(ring < 0 || m_rings <= ring) {//starting at 0, last entry is nRings-1
    throw std::out_of_range("Ring out of range:");
  } else if( pad < 0 || m_nPhiSegments[ring] <= pad ) {//starting at 0
    throw std::out_of_range("Pad out of range:");
  }

  return (layer-1) * (m_padsPerLayer) + m_padsBeforeRing[ring] + (pad);
}


void BeamCalGeoDD::getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const{

  //how often does nPadsPerLayer fit into padIndex;
  //layer starts at 1!
  layer = getLayer(padIndex);// ( padIndex / m_nPadsPerLayer ) + 1;
  ring = getRing(padIndex);

  const int ringIndex(padIndex % m_padsPerLayer);
  pad = ringIndex - m_padsBeforeRing[ring];

#ifdef DEBUG
  if (padIndex != this->getPadIndex(layer, ring, pad) ) {
    std::stringstream error;
    error << "PadIndex " 
  	  << std::setw(7) << padIndex
  	  << std::setw(7) << this->getPadIndex(layer, ring, pad);
    throw std::logic_error(error.str());
  }
#endif

  return;

}//getLayerRingPad

int BeamCalGeoDD::getRing(int padIndex) const {
  //int ringIndex = padIndex % m_nPadsPerLayer;
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

//ARG
// double BeamCalGeoDD::getPadPhi(int ring, int pad) const {
//   double phi = 0;
//   if( ring < m_firstFullRing) {
//     phi += m_fullKeyholeCutoutAngle/2.0;
//     double deltaPhi = (360.0 - m_fullKeyholeCutoutAngle/2.0)/double(m_Segments[ring]);
//     phi +=  deltaPhi * double(pad) + deltaPhi/2.0;
//   } else {
//     double deltaPhi = 360.0 / double(m_Segments[ring]);
//     //Pads still start at the top of the cutout! which is negative -180 + half cutoutangle
//     phi +=  m_fullKeyholeCutoutAngle/2.0 + deltaPhi * double(pad) + deltaPhi/2.0;
//   }

//   phi  += 180.0;
//   if( phi > 360.0 ) {
//     phi -= 360.0;
//   }

//   //  std::cout << std::setw(10) << ring << std::setw(10) << pad << std::setw(10) << phi  << std::endl;
//   return phi;
// }


// double BeamCalGeoDD::getPadPhi(int globalPandIndex) const {
//   int layer, ring, pad;
//   getLayerRingPad(globalPandIndex, layer, ring, pad);
//   return getPadPhi(ring, pad);
// }

double BeamCalGeoDD::getThetaFromRing(int layer, double averageRing) const  {
  const double radiusStep = getBCOuterRadius() - getBCInnerRadius();
  const double radius = getBCInnerRadius() + ( averageRing + 0.5 )  * ( ( radiusStep  ) / double(m_rings) );
  return atan( radius / getLayerZDistanceToIP(layer) );
}



/**
 * Returns true if the pads are in adjacent cells in the same ring
 * or if they are in adjacent rings, then it will return true if the overlap in phi is smaller than
 * some amount of degrees...
 * do the layers have to be the same? --> use boolean flag
 * What about adjacent layers and the same index? we are going to write a second function for those, if we have to
 * Pads are Neighbours if they are in the same Ring, but have Pads with a difference of one
 * What about the keyhole cutout? They will not be considered to be neighbours
 */

bool BeamCalGeoDD::arePadsNeighbours(int globalPadIndex1, int globalPadIndex2, bool mustBeInSameLayer) const {
  int layerIndex1, layerIndex2, ringIndex1, ringIndex2, padIndex1, padIndex2;
  getLayerRingPad(globalPadIndex1, layerIndex1, ringIndex1, padIndex1);
  getLayerRingPad(globalPadIndex2, layerIndex2, ringIndex2, padIndex2);

  if( mustBeInSameLayer and ( not ( layerIndex1 == layerIndex2 ) ) ) {
    //    std::cout << "Not Neighbours: Layers differ"  << std::endl;
    return false;
  }
  //If they are in the same ring
  if( ringIndex1 == ringIndex2 ) {

    if( 
       //The index of neighbouring pads differs by 1
       ( abs( padIndex1 - padIndex2 ) <= 1 ) or
       // or one is zero and the other is m_segments[ringIndex1] - 1
       (  ringIndex1 >= getFirstFullRing() && ( ( abs( padIndex1 - padIndex2 ) + 1 )  == m_nPhiSegments[ringIndex1] )  )
	) {
      // std::cout << "Neighbours: Adjacent pads in same ring"  
      // 		<< std::setw(5) << padIndex1
      // 		<< std::setw(5) << padIndex2
      // 		<< std::endl;
      return true;
    } else { //if they are in the same ring, but not neighbours, we can bail out here
      // std::cout << "Not Neighbours: Same ring but not adjacent pads" 
      // 		<< std::setw(5) << padIndex1
      // 		<< std::setw(5) << padIndex2
      // 		<< std::endl;
      return false;
    }
  }//if in the same ring

  //adjacent rings' indices differ by 1
  if( abs( ringIndex1 - ringIndex2 ) == 1 ) {
    //Angles are in Degrees!
    double phi1 = getPadPhi(ringIndex1, padIndex1);
    double phi2 = getPadPhi(ringIndex2, padIndex2);

    if( fabs( phi1 - phi2 ) < 18.0 ) {
      // std::cout << "Neighbours: Adjacent pads in different ring"  
      // 		<< std::setw(10) << phi1
      // 		<< std::setw(10) << phi2
      // 		<< std::endl;
      return true;
    } else {
      // std::cout << "Not Neighbours: Non adjacent pads in different ring"  
      // 		<< std::setw(10) << phi1
      // 		<< std::setw(10) << phi2
      // 		<< std::endl;
      return false;
    }
  }

  return false;

}



void BeamCalGeoDD::getPadExtentsById(int globalPadIndex, double *extents) const
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


double BeamCalGeoDD::getPadsDistance(int padIndex1, int padIndex2) const
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

void BeamCalGeoDD::setFirstFullRing() {
  int ring = 0;
  while ( m_cutOut > getRadSegmentation()[ring] ) { ++ring; }
  m_firstFullRing = ring;
}

void BeamCalGeoDD::setCutOut() {
  const double cutOutRadius =  tan(m_crossingAngle/1000.0)*getLayerZDistanceToIP(m_layers-1)+3.5;
  m_cutOut = cutOutRadius;
}

