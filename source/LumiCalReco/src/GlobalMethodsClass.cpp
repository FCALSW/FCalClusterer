
#include "GlobalMethodsClass.h"
#include "MarlinLumiCalClusterer.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearMgr.h>


#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

using streamlog::MESSAGE;

#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <sstream>

// optional DD4hep
#ifdef FCAL_WITH_DD4HEP
#include <DD4hep/LCDD.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/DetectorData.h>
#include <DDRec/API/IDDecoder.h>
#include <XML/XMLChildValue.h>
#endif

// utility copied from marlin 
template <class T>
bool convert(std::string input, T &value) {
 std::istringstream stream(input);
 return ( ! (stream >> std::setbase(0) >> value).fail() ) && stream.eof();
}

GlobalMethodsClass::WeightingMethod_t GlobalMethodsClass::LogMethod = "LogMethod";
GlobalMethodsClass::WeightingMethod_t GlobalMethodsClass::EnergyMethod = "EnergyMethod";
double GlobalMethodsClass::EnergyCalibrationFactor = 0.0105;

GlobalMethodsClass :: GlobalMethodsClass() :
  _procName( "MarlinLumiCalClusterer" ),
  GlobalParamI(),
  GlobalParamD(),
  GlobalParamS()
{
}
GlobalMethodsClass :: GlobalMethodsClass(const std::string &name) :
  _procName( name ),
  GlobalParamI(),
  GlobalParamD(),
  GlobalParamS()
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


void GlobalMethodsClass::SetConstants() {


  
  // BP. access processor parameters
  marlin::StringParameters* _lcalRecoPars = NULL;
  marlin::ProcessorMgr* procMgr= marlin::ProcessorMgr::instance();
  marlin::Processor* procPTR = procMgr->getActiveProcessor( _procName );
  //procMgr->dumpRegisteredProcessorsXML();


  if ( procPTR )
    _lcalRecoPars = procPTR->parameters();
  else{                                                // try with My+_procName
    _procName = "My" + _procName;
    procPTR = procMgr->getActiveProcessor( _procName ); 
    _lcalRecoPars = procPTR -> parameters();
  }
  if( !_lcalRecoPars ) {
    streamlog_out( MESSAGE ) << " GlobalMethodsClass::SetConstants : Non of processor names ("<<_procName<< ", My"<<_procName<< ") found !" << std::endl;
    exit( EXIT_FAILURE );
  }

  // BP. GEAR access
  // geometry parameters of LumiCal
  // inner and outer radii [mm]
  // 
#ifdef FCAL_WITH_DD4HEP

  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  DD4hep::Geometry::DetElement mLCal = lcdd.detector("LumiCal");
  const DD4hep::DDRec::LayeredCalorimeterData * theExtension = mLCal.extension<DD4hep::DDRec::LayeredCalorimeterData>();
  const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer>& mLayers = theExtension->layers;
 
  GlobalParamD[RMin]   = theExtension->extent[0];
  GlobalParamD[RMax]   = theExtension->extent[1];
  GlobalParamD[ZStart] = theExtension->extent[2];
  GlobalParamD[ZEnd]   = theExtension->extent[3];

  // all layers LumiCal are identical, take first one
  const DD4hep::DDRec::LayeredCalorimeterStruct::Layer & oneLayer = mLayers.at(0); 
  GlobalParamD[RCellLength]   = oneLayer.cellSize0;
  GlobalParamD[PhiCellLength] = oneLayer.cellSize1;
  GlobalParamI[NumCellsR]     = (int)( (GlobalParamD[RMax]-GlobalParamD[RMin]) / GlobalParamD[RCellLength]);
  GlobalParamI[NumCellsPhi]   = (int)(2.0 * M_PI / GlobalParamD[PhiCellLength] + 0.5);
  GlobalParamI[NumCellsZ]     = mLayers.size();
  GlobalParamD[ZLayerThickness]   = oneLayer.inner_thickness + oneLayer.outer_thickness;
  GlobalParamD[BeamCrossingAngle] = lcdd.constantAsDouble("crossingangle") / 1000.;

#else
  gear::GearMgr* gearMgr =  marlin::Global::GEAR;
  const gear::CalorimeterParameters& _lcalPars = gearMgr->getLcalParameters();

  GlobalParamD[RMin] = _lcalPars.getExtent()[0];
  GlobalParamD[RMax] = _lcalPars.getExtent()[1];

  // starting/end position [mm]
  GlobalParamD[ZStart] = _lcalPars.getExtent()[2];
  GlobalParamD[ZEnd]   = _lcalPars.getExtent()[3];

  // cell division numbers
  GlobalParamD[RCellLength]   = _lcalPars.getLayerLayout().getCellSize0(0);
  GlobalParamD[PhiCellLength] = _lcalPars.getLayerLayout().getCellSize1(0);
  GlobalParamI[NumCellsR]   = (int)( (GlobalParamD[RMax]-GlobalParamD[RMin]) / GlobalParamD[RCellLength]);
  GlobalParamI[NumCellsPhi] = (int)(2.0 * M_PI / GlobalParamD[PhiCellLength] + 0.5);
  GlobalParamI[NumCellsZ] = _lcalPars.getLayerLayout().getNLayers();

 // beam crossing angle ( convert to rad )
  GlobalParamD[BeamCrossingAngle] = _lcalPars.getDoubleVal("beam_crossing_angle") / 1000.; 

  // layer thickness
  GlobalParamD[ZLayerThickness] = _lcalPars.getLayerLayout().getThickness(0);

  
#endif

  //------------------------------------------------------------------------ 
  // Processor Parameters 
  // Clustering/Reco parameters
  //(BP) layer relative phi offset - must go sometimes to GEAR params
  const std::string parname = "ZLayerPhiOffset";
  double val = 0.;
    if ( _lcalRecoPars->isParameterSet(parname) ){ 
      val = _lcalRecoPars->getFloatVal( parname );
    }else {
      marlin::ProcessorParameter* par = marlin::CMProcessor::instance()->getParam( _procName, parname );
      //      val = (double)std::stof( par->defaultValue() );
      if ( convert( par->defaultValue(), val ) ){
	streamlog_out(WARNING)<<"\tParameter <"<< parname <<"> not set default value : "<< val << "\t is used"<<"\n";
      }else{ 
	streamlog_out(WARNING)<<"\tParameter <"<< parname <<"> not set default value : "<< 0.  << "\t is used"<<"\n";
      }
   }
   
  // check units just in case ( convert to rad as needed )
  val = ( val <= GlobalParamD[PhiCellLength] ) ? val : val*M_PI/180.;
  GlobalParamD[ZLayerPhiOffset] = val; 

  EnergyCalibrationFactor = _lcalRecoPars->getFloatVal("EnergyCalibConst");

  // logarithmic constant for position reconstruction
  GlobalParamD[LogWeightConstant] = _lcalRecoPars->getFloatVal("LogWeigthConstant");
  
  GlobalParamD[MinHitEnergy] = _lcalRecoPars->getFloatVal("MinHitEnergy");
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
  GlobalParamD[MinClusterEngySignal] = SignalGevConversion(GeV_to_Signal , GlobalParamD[MinClusterEngyGeV]); 
  // conversion factor "detector Signal to primary particle energy"
  GlobalParamD[Signal_to_GeV] = SignalGevConversion( Signal_to_GeV, 1. );             
  // hits positions weighting method 
  GlobalParamS[WeightingMethod] = _lcalRecoPars->getStringVal(  "WeightingMethod" );
  GlobalParamD[ElementsPercentInShowerPeakLayer] = _lcalRecoPars->getFloatVal("ElementsPercentInShowerPeakLayer");
  GlobalParamI[NumOfNearNeighbor] = _lcalRecoPars->getIntVal("NumOfNearNeighbor");
  // IO
  GlobalParamS[LumiInColName] = _lcalRecoPars->getStringVal(  "LumiInColName" );

}


/* --------------------------------------------------------------------------
   ccccccccc
   -------------------------------------------------------------------------- */
double GlobalMethodsClass::SignalGevConversion( Parameter_t optName , double valNow ){

#pragma message("FIXME: SignalToGeV conversion")
  double	returnVal = -1;
  //  const double conversionFactor(0.0105*1789/1500*(1488/1500.0));
  //  const double conversionFactor(0.0105);
  //const double conversionFactor(0.0105);

  if(optName == GeV_to_Signal)
    returnVal = valNow * EnergyCalibrationFactor;
  //		returnVal = valNow * .0105 + .0013;

  if(optName == Signal_to_GeV)
    returnVal = valNow / EnergyCalibrationFactor;
  //		returnVal = (valNow-.0013) / .0105;

  return returnVal;
}



void GlobalMethodsClass::ThetaPhiCell(int cellId , std::map <GlobalMethodsClass::Coordinate_t , double> &thetaPhiCell) {

  // compute Z,Phi,R coordinates according to the cellId
  // returned Phi is in the range (-M_PI, M_PI )
 
  int cellIdZ, cellIdPhi, cellIdR, arm;
  CellIdZPR(cellId, cellIdZ, cellIdPhi, cellIdR, arm);

  // theta
  double rCell      = GlobalParamD[RMin] + (cellIdR + .5) * GlobalParamD[RCellLength];
  double zCell      = fabs(GlobalParamD[ZStart]) + GlobalParamD[ZLayerThickness] * (cellIdZ - 1);
  double thetaCell  = atan(rCell / zCell);

  // phi
  //(BP) use phiCell size and account for possible layers relative offset/stagger
  // double phiCell   = 2 * M_PI * (double(cellIdPhi) + .5) / double(GlobalParamI[NumCellsPhi]) + double( cellIdZ % 2 ) * GlobalParamD[;
  double phiCell   = (double(cellIdPhi) + .0) * GlobalParamD[PhiCellLength] + double( (cellIdZ) % 2 ) * GlobalParamD[ZLayerPhiOffset];
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
  case PhiCellLength:	                 return "PhiCellLength";
  case ZLayerThickness:	                 return "ZLayerThickness";
  case ZLayerPhiOffset:	                 return "ZLayerPhiOffset";
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
  case MinClusterEngySignal:	         return "MinClusterEngySignal";
  case MiddleEnergyHitBoundFrac:         return "MiddleEnergyHitBoundFrac";
  case WeightingMethod:                  return "WeightingMethod";
  case GeV_to_Signal:	                 return "GeV_to_Signal";
  case Signal_to_GeV:                    return "Signal_to_GeV";
  case BeamCrossingAngle:                return "BeamCrossingAngle";
  case LumiInColName:                    return "LumiInColName";
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
  streamlog_out(MESSAGE) << "---------------------------------------------------------------" << std::endl;

}
