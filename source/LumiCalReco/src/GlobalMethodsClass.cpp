
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
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>

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
int GlobalMethodsClass::CellIdZPR(int cellZ, int cellPhi, int cellR, int arm) {

  int cellId = 0;

  cellId |= ( cellZ   <<  0 ) ;
  cellId |= ( cellPhi << 10 ) ;
  cellId |= ( cellR   << 20 ) ;
  cellId |= ( arm     << 30 ) ;

  return cellId;
}

void GlobalMethodsClass::CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm) {

  // compute Z,Phi,R indices according to the cellId
  cellZ   = (cellId >> 0  ) & (int)(( 1 << 10 ) -1) ; 
  cellPhi = (cellId >> 10 ) & (int)(( 1 << 10 ) -1) ;
  cellR   = (cellId >> 20 ) & (int)(( 1 << 10 ) -1) ;
  arm     = (cellId >> 30 ) & (int)(( 1 <<  2 ) -1) ; arm = ( arm == 0 ) ? -1 : 1;
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


  // BP. access gear file
  gear::GearMgr* gearMgr =  marlin::Global::GEAR;
  const gear::CalorimeterParameters& _lcalGearPars = gearMgr->getLcalParameters();
  
  // BP. access processor parameters
  marlin::StringParameters* _lcalRecoPars = NULL;
  marlin::ProcessorMgr* procMgr= marlin::ProcessorMgr::instance();
  marlin::Processor* procPTR = procMgr->getActiveProcessor( _procName );
  //procMgr->dumpRegisteredProcessorsXML();


  if ( procPTR )
    _lcalRecoPars = procPTR->parameters();
  else{                                                // try with My+_procName
    std::string my_procName = "My" + _procName; 
    _lcalRecoPars = procMgr->getActiveProcessor( my_procName )->parameters();
  }
  if( !_lcalRecoPars ) {
    streamlog_out( MESSAGE ) << " GlobalMethodsClass::SetConstants : Non of processor names ("<<_procName<< ", My"<<_procName<< ") found !" << std::endl;
    exit( EXIT_FAILURE );
  }

  // geometry parameters of LumiCal
  // inner and outer radii [mm]
  GlobalParamD[RMin] = _lcalGearPars.getExtent()[0];
  GlobalParamD[RMax] = _lcalGearPars.getExtent()[1];

  // starting/end position [mm]
  GlobalParamD[ZStart] = _lcalGearPars.getExtent()[2];
  GlobalParamD[ZEnd]   = _lcalGearPars.getExtent()[3];

  // cell division numbers
  GlobalParamD[RCellLength]   = _lcalGearPars.getLayerLayout().getCellSize0(0);
  GlobalParamD[PhiCellLength] = _lcalGearPars.getLayerLayout().getCellSize1(0);
  GlobalParamI[NumCellsR]   = (int)( (GlobalParamD[RMax]-GlobalParamD[RMin]) / GlobalParamD[RCellLength]);
  GlobalParamI[NumCellsPhi] = (int)(2.0 * M_PI / GlobalParamD[PhiCellLength] + 0.5);
  GlobalParamI[NumCellsZ] = _lcalGearPars.getLayerLayout().getNLayers();

 // beam crossing angle
  GlobalParamD[BeamCrossingAngle] = _lcalGearPars.getDoubleVal("beam_crossing_angle") / 1000.; 

  // layer thickness
  GlobalParamD[ZLayerThickness] = _lcalGearPars.getLayerLayout().getThickness(0);
  //(BP) layer relative phi offset - must go sometimes to GEAR params
  const std::string parname = "ZLayerStagger";
  double val = 0.;
    if ( _lcalRecoPars->isParameterSet(parname) ){ 
      val = _lcalRecoPars->getFloatVal( parname );
    }else {
      marlin::ProcessorParameter* par = marlin::CMProcessor::instance()->getParam( _procName, parname );
      val = (double)std::stof( par->defaultValue() );
      streamlog_out(MESSAGE)<<"\tParameter <"<< parname <<"> not set default value : "<< val << "\t is used"<<"\n"; 
    }
   
  // check units just in case
  val = ( val < GlobalParamD[PhiCellLength] ) ? val : val*M_PI/180.;
  GlobalParamD[ZLayerPhiOffset] = val;  

  // Clustering/Reco parameters
  // logarithmic constant for position reconstruction
  EnergyCalibrationFactor = _lcalRecoPars->getFloatVal("EnergyCalibConst");

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
