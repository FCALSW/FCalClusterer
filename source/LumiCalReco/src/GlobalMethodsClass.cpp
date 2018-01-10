
#include "GlobalMethodsClass.h"
#include "MarlinLumiCalClusterer.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearMgr.h>

#ifdef FCAL_WITH_DD4HEP

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/DetectorData.h>

#endif


#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

using streamlog::MESSAGE;

#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

// utility copied from marlin 
template <class T>
bool convert(std::string input, T &value) {
 std::istringstream stream(input);
 return ( ! (stream >> std::setbase(0) >> value).fail() ) && stream.eof();
}

GlobalMethodsClass :: GlobalMethodsClass() :
  _useDD4hep(false),
  _backwardRotationPhi(0.0),
  GlobalParamI(),
  GlobalParamD(),
  GlobalParamS(),
  _forwardCalo(NULL),
  _backwardCalo(NULL)
{
}

GlobalMethodsClass :: ~GlobalMethodsClass(){
}

/* --------------------------------------------------------------------------
   (1):	return a cellId for a given Z (layer), R (cylinder) and Phi (sector)
   (2):	return Z (layer), R (cylinder) and Phi (sector) for a given cellId
   -------------------------------------------------------------------------- */

#define SHIFT_I_32Fcal 0  // I = 10 bits  ( ring )
#define SHIFT_J_32Fcal 10 // J = 10 bits  ( sector)
#define SHIFT_K_32Fcal 20 // K = 10 bits  ( layer )
#define SHIFT_S_32Fcal 30 // S =  2 bits  ( side/arm ) 

#define MASK_I_32Fcal (unsigned int) 0x000003FF
#define MASK_J_32Fcal (unsigned int) 0x000FFC00
#define MASK_K_32Fcal (unsigned int) 0x3FF00000
#define MASK_S_32Fcal (unsigned int) 0xC0000000

int GlobalMethodsClass::CellIdZPR(int cellZ, int cellPhi, int cellR, int arm) {

  int cellId = 0;
  int side = ( arm < 0 ) ? 0 : arm;
  cellId  = (
          ( (side    << SHIFT_S_32Fcal) & MASK_S_32Fcal) |
          ( (cellR   << SHIFT_I_32Fcal) & MASK_I_32Fcal) |
          ( (cellPhi << SHIFT_J_32Fcal) & MASK_J_32Fcal) |
          ( (cellZ   << SHIFT_K_32Fcal) & MASK_K_32Fcal)
               );
  return cellId;
}

void GlobalMethodsClass::CellIdZPR(int cellID, int& cellZ, int& cellPhi, int& cellR, int& arm) {

  // compute Z,Phi,R indices according to the cellId

        cellR   = ((((unsigned int)cellID)&MASK_I_32Fcal) >> SHIFT_I_32Fcal);
        cellPhi = ((((unsigned int)cellID)&MASK_J_32Fcal) >> SHIFT_J_32Fcal);
        cellZ   = ((((unsigned int)cellID)&MASK_K_32Fcal) >> SHIFT_K_32Fcal);
        arm     = ((((unsigned int)cellID)&MASK_S_32Fcal) >> SHIFT_S_32Fcal);
  return;
}

int GlobalMethodsClass::CellIdZPR(int cellId, GlobalMethodsClass::Coordinate_t ZPR) {

 
  int cellZ, cellPhi, cellR, arm;
  CellIdZPR(cellId, cellZ, cellPhi, cellR, arm);
  arm = ( arm == 0 ) ? -1 : 1;
  if(ZPR == GlobalMethodsClass::COZ) return cellZ;
  else if(ZPR == GlobalMethodsClass::COR) return cellR;
  else if(ZPR == GlobalMethodsClass::COP) return cellPhi;
  else if(ZPR == GlobalMethodsClass::COA) return arm;

  return 0;
}


void GlobalMethodsClass::SetConstants( marlin::Processor* procPTR ) {

  
  std::shared_ptr<marlin::StringParameters> _lcalRecoPars = procPTR->parameters();

  //SetGeometryConstants
  if( SetGeometryDD4HEP() ) {
    _useDD4hep=true;
  } else {
    SetGeometryGear();
  }

  //------------------------------------------------------------------------ 
  // Processor Parameters 
  // Clustering/Reco parameters
  //(BP) layer relative phi offset - must go sometimes to GEAR params
  const std::string parname = "ZLayerPhiOffset";
  double val = 0.;
    if ( _lcalRecoPars->isParameterSet(parname) ){ 
      val = _lcalRecoPars->getFloatVal( parname );
    }else {
      marlin::ProcessorParameter* par = marlin::CMProcessor::instance()->getParam( procPTR->type(), parname );
      if ( convert( par->defaultValue(), val ) ){
	streamlog_out(WARNING)<<"\tParameter <"<< parname <<"> not set default value : "<< val << "\t is used"<<"\n";
      }else{ 
	streamlog_out(WARNING)<<"\tParameter <"<< parname <<"> not set default value : "<< 0.  << "\t is used"<<"\n";
      }
   }
   
  // check units just in case ( convert to rad as needed )
  val = ( val <= GlobalParamD[PhiCellLength] ) ? val : val*M_PI/180.;
  GlobalParamD[ZLayerPhiOffset] = val; 

  GlobalParamD[Signal_to_GeV] = 1. / _lcalRecoPars->getFloatVal("EnergyCalibConst");

  // logarithmic constant for position reconstruction
  GlobalParamD[LogWeightConstant] = _lcalRecoPars->getFloatVal("LogWeigthConstant");
  
  GlobalParamD[MinHitEnergy] = toGev(_lcalRecoPars->getFloatVal("MinHitEnergy"));
  GlobalParamD[MiddleEnergyHitBoundFrac] = _lcalRecoPars->getFloatVal(  "MiddleEnergyHitBoundFrac" );
  GlobalParamD[ElementsPercentInShowerPeakLayer] = _lcalRecoPars->getFloatVal(  "ElementsPercentInShowerPeakLayer" );
  GlobalParamI[ClusterMinNumHits]   = _lcalRecoPars->getIntVal(  "ClusterMinNumHits" );

  // Moliere radius of LumiCal [mm]
  GlobalParamD[MoliereRadius] = _lcalRecoPars->getFloatVal(  "MoliereRadius" );

  // Geometrical fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  // (BP) Note, this in local LumiCal Reference System ( crossing angle not accounted )
  // quite large - conservative, further reco-particles selection can be done later if desired
  GlobalParamD[ThetaMin] = (GlobalParamD[RMin] + GlobalParamD[MoliereRadius])/GlobalParamD[ZEnd];
  GlobalParamD[ThetaMax] = (GlobalParamD[RMax] - GlobalParamD[MoliereRadius])/GlobalParamD[ZStart];
 
  // minimal separation distance between any pair of clusters [mm]
  GlobalParamD[MinSeparationDist] = GlobalParamD[MoliereRadius];

  // minimal energy of a single cluster
  GlobalParamD[MinClusterEngyGeV] = _lcalRecoPars->getFloatVal(  "MinClusterEngy" );  // value in GeV

  // hits positions weighting method 
  GlobalParamS[WeightingMethod] = _lcalRecoPars->getStringVal(  "WeightingMethod" );
  GlobalParamD[ElementsPercentInShowerPeakLayer] = _lcalRecoPars->getFloatVal("ElementsPercentInShowerPeakLayer");
  GlobalParamI[NumOfNearNeighbor] = _lcalRecoPars->getIntVal("NumOfNearNeighbor");
  // IO
  GlobalParamS[LumiInColName] = _lcalRecoPars->getStringVal(  "LumiCal_Collection" );

  initializeAdditionalParameters();
}

double GlobalMethodsClass::toSignal(double value) const { return value * getCalibrationFactor(); }

double GlobalMethodsClass::toGev(double value) const { return value / getCalibrationFactor(); }

void GlobalMethodsClass::ThetaPhiCell(const int cellId,
                                      std::map<GlobalMethodsClass::Coordinate_t, double>& thetaPhiCell) const {
  // compute Z,Phi,R coordinates according to the cellId
  // returned Phi is in the range (-M_PI, M_PI )
 
  int cellIdZ, cellIdPhi, cellIdR, arm;
  CellIdZPR(cellId, cellIdZ, cellIdPhi, cellIdR, arm);

  // theta
  double rCell = GlobalParamD.at(RMin) + (cellIdR + 0.5) * GlobalParamD.at(RCellLength) - GlobalParamD.at(RCellOffset);
  double zCell =
      fabs(GlobalParamD.at(ZStart)) + GlobalParamD.at(ZLayerThickness) * (cellIdZ) + GlobalParamD.at(ZLayerZOffset);
  double thetaCell  = atan(rCell / zCell);

  // phi
  //(BP) use phiCell size and account for possible layers relative offset/stagger
  // double phiCell   = 2 * M_PI * (double(cellIdPhi) + .5) / double(GlobalParamI[NumCellsPhi]) + double( cellIdZ % 2 ) * GlobalParamD[;
  double phiCell = (double(cellIdPhi)) * GlobalParamD.at(PhiCellLength) +
                   double((cellIdZ) % 2) * GlobalParamD.at(ZLayerPhiOffset) + GlobalParamD.at(PhiCellOffset);
  phiCell = ( phiCell > M_PI ) ? phiCell-2.*M_PI : phiCell;
  // fill output container
  thetaPhiCell[GlobalMethodsClass::COTheta] = thetaCell;
  thetaPhiCell[GlobalMethodsClass::COPhi]   = phiCell;
  thetaPhiCell[GlobalMethodsClass::COR]     = rCell;
  thetaPhiCell[GlobalMethodsClass::COZ]     = zCell;
  return;
}


std::string GlobalMethodsClass::GetParameterName ( Parameter_t par ){

  switch (par) {
  case ZStart:                           return "ZStart";
  case ZEnd:                             return "ZEnd";
  case RMin:		                 return "RMin";
  case RMax:		                 return "RMax";
  case NumCellsR:	                 return "NumCellsR";
  case NumCellsPhi:	                 return "NumCellsPhi";
  case NumCellsZ:	                 return "NumCellsZ";
  case RCellLength:	                 return "RCellLength";
  case RCellOffset:	                 return "RCellOffset";
  case PhiCellLength:	                 return "PhiCellLength";
  case PhiCellOffset:	                 return "PhiCellOffset";
  case ZLayerThickness:	                 return "ZLayerThickness";
  case ZLayerPhiOffset:	                 return "ZLayerPhiOffset";
  case ZLayerZOffset:	                 return "ZLayerZOffset";
  case ThetaMin:	                 return "ThetaMin";
  case ThetaMax:	                 return "ThetaMax";
  case LogWeightConstant:                return "LogWeightConstant";
  case MoliereRadius:	                 return "MoliereRadius";
  case MinSeparationDist:                return "MinSeparationDist";
  case ElementsPercentInShowerPeakLayer: return "ElementsPercentInShowerPeakLayer";
  case NumOfNearNeighbor:                return "NumOfNearNeighbor";
  case ClusterMinNumHits:                return "ClusterMinNumHits";
  case MinHitEnergy:                     return "MinHitEnergy";
  case MinClusterEngyGeV:	         return "MinClusterEngyGeV";
  case MiddleEnergyHitBoundFrac:         return "MiddleEnergyHitBoundFrac";
  case WeightingMethod:                  return "WeightingMethod";
  case Signal_to_GeV:                    return "Signal_to_GeV";
  case BeamCrossingAngle:                return "BeamCrossingAngle";
  case LumiInColName:                    return "LumiInColName";
  case BetaGamma:                        return "BetaGamma";
  case Gamma:                            return "Gamma";
  default: return "Unknown Parameter";
  }


}


void GlobalMethodsClass::PrintAllParameters() const {
  streamlog_out(MESSAGE) << "------------------------------------------------------------------" << std::endl;
  streamlog_out(MESSAGE) << "********* LumiCalReco Parameters set in GlobalMethodClass ********" << std::endl;

  for (ParametersInt::const_iterator it = GlobalParamI.begin();it != GlobalParamI.end() ;++it) {
    streamlog_out(MESSAGE) << " - (int)     " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

  for (ParametersDouble::const_iterator it = GlobalParamD.begin();it != GlobalParamD.end() ;++it) {
    streamlog_out(MESSAGE) << " - (double)  " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

  for (ParametersString::const_iterator it = GlobalParamS.begin();it != GlobalParamS.end() ;++it) {
    streamlog_out(MESSAGE) << " - (string)  " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

  streamlog_out(MESSAGE) << "Using DD4hep based geometry and hit enconding? " << std::boolalpha << _useDD4hep  << std::endl;

  streamlog_out(MESSAGE) << "---------------------------------------------------------------" << std::endl;

}

void GlobalMethodsClass::SetGeometryGear() {

  // GEAR access
  // BP. access gear file
  gear::GearMgr* gearMgr =  marlin::Global::GEAR;
  const gear::CalorimeterParameters& _lcalGearPars = gearMgr->getLcalParameters();
  // geometry parameters of LumiCal
  // inner and outer radii [mm]
  GlobalParamD[RMin] = _lcalGearPars.getExtent()[0];
  GlobalParamD[RMax] = _lcalGearPars.getExtent()[1];

  // starting/end position [mm]
  GlobalParamD[ZStart] = _lcalGearPars.getExtent()[2];
  GlobalParamD[ZEnd]   = _lcalGearPars.getExtent()[3];

  // cell division numbers
  GlobalParamD[RCellLength]   = _lcalGearPars.getLayerLayout().getCellSize0(0);
  GlobalParamD[RCellOffset]   = 0.0;
  GlobalParamD[PhiCellLength] = _lcalGearPars.getLayerLayout().getCellSize1(0);
  GlobalParamD[PhiCellOffset] = 0.0;
  GlobalParamI[NumCellsR]   = (int)( (GlobalParamD[RMax]-GlobalParamD[RMin]) / GlobalParamD[RCellLength]);
  GlobalParamI[NumCellsPhi] = (int)(2.0 * M_PI / GlobalParamD[PhiCellLength] + 0.5);
  GlobalParamI[NumCellsZ] = _lcalGearPars.getLayerLayout().getNLayers();

 // beam crossing angle ( convert to rad )
  GlobalParamD[BeamCrossingAngle] = _lcalGearPars.getDoubleVal("beam_crossing_angle") / 1000.; 

  // layer thickness
  GlobalParamD[ZLayerThickness] = _lcalGearPars.getLayerLayout().getThickness(0);
  GlobalParamD[ZLayerZOffset]   = 0.0;
}


bool GlobalMethodsClass::SetGeometryDD4HEP() {
#ifdef FCAL_WITH_DD4HEP

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  if( theDetector.detectors().count("LumiCal") == 0 ) return false;

  dd4hep::DetElement lumical(theDetector.detector("LumiCal"));
  dd4hep::Segmentation readout(theDetector.readout("LumiCalCollection").segmentation());

  streamlog_out(MESSAGE) << "Segmentation Type" << readout.type()  << std::endl;
  streamlog_out(MESSAGE) <<"FieldDef: " << readout.segmentation()->fieldDescription()  << std::endl;

  const dd4hep::rec::LayeredCalorimeterData * theExtension = lumical.extension<dd4hep::rec::LayeredCalorimeterData>();
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers= theExtension->layers;

  GlobalParamD[RMin] = theExtension->extent[0]/dd4hep::mm;
  GlobalParamD[RMax] = theExtension->extent[1]/dd4hep::mm;

  // starting/end position [mm]
  GlobalParamD[ZStart] = theExtension->extent[2]/dd4hep::mm;
  GlobalParamD[ZEnd]   = theExtension->extent[3]/dd4hep::mm;

  // cell division numbers
  typedef dd4hep::DDSegmentation::TypedSegmentationParameter< double > ParDou;
  ParDou* rPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_r"));
  ParDou* rOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_r"));
  ParDou* pPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_phi"));
  ParDou* pOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_phi"));

  if (rPar == NULL or pPar == NULL ) {
    throw std::runtime_error( "Could not obtain parameters from segmentation" );
  }

  GlobalParamD[RCellLength]   = rPar->typedValue()/dd4hep::mm;
  GlobalParamD[RCellOffset]   = 0.5 * GlobalParamD[RCellLength] + GlobalParamD[RMin] - rOff->typedValue() / dd4hep::mm;
  GlobalParamD[PhiCellLength] = pPar->typedValue()/dd4hep::radian;
  GlobalParamD[PhiCellOffset] = pOff->typedValue() / dd4hep::radian;
  GlobalParamI[NumCellsR]     = (int)( (GlobalParamD[RMax]-GlobalParamD[RMin]) / GlobalParamD[RCellLength]);
  GlobalParamI[NumCellsPhi]   = (int)(2.0 * M_PI / GlobalParamD[PhiCellLength] + 0.5);
  GlobalParamI[NumCellsZ]     = layers.size();

  // beam crossing angle ( convert to rad )
  dd4hep::DetElement::Children children = lumical.children();

  if( children.empty() ) {
    throw std::runtime_error( "Cannot obtain crossing angle from this LumiCal, update lcgeo?" );
  }

  for( dd4hep::DetElement::Children::const_iterator it = children.begin(); it != children.end(); ++it ) {
    dd4hep::Position loc(0.0, 0.0, 0.0);
    dd4hep::Position glob(0.0, 0.0, 0.0);
    it->second.nominal().localToWorld( loc, glob );
    GlobalParamD[BeamCrossingAngle] = 2.0*fabs( atan( glob.x() / glob.z() ) / dd4hep::rad );
    if( glob.z() > 0.0 ) {
      _forwardCalo = &it->second.nominal().worldTransformation();
    } else {
      _backwardCalo = &it->second.nominal().worldTransformation();

      //get phi rotation from global to local transformation
      TGeoHMatrix *tempMat = (TGeoHMatrix*) _backwardCalo->Clone();
      double nulltr[] = { 0.0, 0.0, 0.0 };
      // undo backward and crossing angle rotation
      tempMat->SetTranslation( nulltr );
      // root matrices need degrees as argument
      tempMat->RotateY( GlobalParamD[BeamCrossingAngle]/2.0 * 180/M_PI );
      tempMat->RotateY( -180.0 );
      double local[] =  { 0.0, 1.0, 0.0 };
      double global[] = { 0.0, 0.0, 0.0 };

      tempMat->LocalToMaster( local, global );

      _backwardRotationPhi = atan2( local[1], local[0] ) - atan2( global[1], global[0] ) ;
      if(_backwardRotationPhi < M_PI) _backwardRotationPhi += 2*M_PI;

      delete tempMat;

    }
  }

  // layer thickness
  GlobalParamD[ZLayerThickness] = (layers[0].inner_thickness + layers[0].outer_thickness)/dd4hep::mm;
  GlobalParamD[ZLayerZOffset] = (layers[0].inner_thickness)/dd4hep::mm;
  
  //successfully created geometry from DD4hep
  return true;
#endif
  //no dd4hep geometry
  return false;
}

GlobalMethodsClass::WeightingMethod_t GlobalMethodsClass::getMethod(std::string const& methodName) const {
  if (methodName == "LogMethod") {
    return LogMethod;
  }
  if (methodName == "EnergyMethod") {
    return EnergyMethod;
  }
  throw std::runtime_error("Unkown weighting method");
}

double GlobalMethodsClass::posWeight(double cellEngy, double totEngy, GlobalMethodsClass::WeightingMethod_t method,
                                     double logWeightConstNow) {
  if (method == GlobalMethodsClass::EnergyMethod)
    return cellEngy;
  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if (method == GlobalMethodsClass::LogMethod) {
    return std::max(0.0, log(cellEngy / totEngy) + logWeightConstNow);
  }
  return -1;
}

void GlobalMethodsClass::initializeAdditionalParameters() {
  // Lorentz boost params
  const double beta       = tan(GlobalParamD[BeamCrossingAngle] / 2.0);
  GlobalParamD[BetaGamma] = beta;
  GlobalParamD[Gamma]     = sqrt(1. + beta * beta);

  m_armCosAngle[-1] = cos(-GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armCosAngle[1]  = cos(GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armSinAngle[-1] = sin(-GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armSinAngle[1]  = sin(GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
}
