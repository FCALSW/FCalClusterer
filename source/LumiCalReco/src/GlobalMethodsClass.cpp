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
#include <iostream>


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

int GlobalMethodsClass::CellIdZPR(int cellId, GlobalMethodsClass::Coordinate_t ZPR) {

  int cellZ, cellPhi, cellR, arm;
  CellIdZPR(cellId, cellZ, cellPhi, cellR, arm);

  if(ZPR == GlobalMethodsClass::COZ) return cellZ;
  if(ZPR == GlobalMethodsClass::COR) return cellR;
  if(ZPR == GlobalMethodsClass::COP) return cellPhi;
  if(ZPR == GlobalMethodsClass::COA) return arm;

  return 0;
}

void GlobalMethodsClass::CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm) {

  // compute Z,Phi,R coordinates according to the cellId
  cellZ   = (cellId >> 0  ) & (int)(( 1 << 10 ) -1) ; 
  cellPhi = (cellId >> 10 ) & (int)(( 1 << 10 ) -1) ;
  cellR   = (cellId >> 20 ) & (int)(( 1 << 10 ) -1) ;
  arm     = (cellId >> 30 ) & (int)(( 1 <<  2 ) -1) ;
  return;

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

  // layer thickness
  GlobalParamD[ZLayerThickness] = _lcalGearPars.getLayerLayout().getThickness(0);


  // beam crossing angle
  GlobalParamD[BeamCrossingAngle] = _lcalGearPars.getDoubleVal("beam_crossing_angle"); 


  // Clustering/Reco parameters
  // logarithmic constant for position reconstruction
  GlobalParamD[LogWeightConstant] = _lcalRecoPars->getFloatVal("LogWeigthConstant");
  
  GlobalParamD[MinHitEnergy] = _lcalRecoPars->getFloatVal("MinHitEnergy");
  GlobalParamD[MiddleEnergyHitBoundFrac] = _lcalRecoPars->getFloatVal(  "MiddleEnergyHitBoundFrac" );
  GlobalParamD[ElementsPercentInShowerPeakLayer] = _lcalRecoPars->getFloatVal(  "ElementsPercentInShowerPeakLayer" );
  GlobalParamI[ClusterMinNumHits]   = _lcalRecoPars->getIntVal(  "ClusterMinNumHits" );

  // Moliere radius of LumiCal [mm]
  GlobalParamD[MoliereRadius] = _lcalRecoPars->getFloatVal(  "MoliereRadius" );

   // Geometrical fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  // BP Note this in local LumiCal Reference System ( crossing angle not accounted )
  GlobalParamD[ThetaMin] = (GlobalParamD[RMin] + GlobalParamD[MoliereRadius])/GlobalParamD[ZStart];
  GlobalParamD[ThetaMax] = (GlobalParamD[RMax] - GlobalParamD[MoliereRadius])/GlobalParamD[ZEnd];
 
  // minimal separation distance between any pair of clusters [mm]
  GlobalParamD[MinSeparationDist] = GlobalParamD[MoliereRadius];

  // minimal energy of a single cluster
  GlobalParamD[MinClusterEngy] = _lcalRecoPars->getFloatVal(  "MinClusterEngy" );  // value in GeV
  GlobalParamD[MinClusterEngy] = SignalGevConversion(GeV_to_Signal , GlobalParamD[MinClusterEngy]);  // conversion to "detector Signal"
  // hits positions weighting method 
  GlobalParamI[WeightingMethod] = _lcalRecoPars->getIntVal(  "WeightingMethod" );
  GlobalParamD[ElementsPercentInShowerPeakLayer] = _lcalRecoPars->getFloatVal("ElementsPercentInShowerPeakLayer");
  GlobalParamI[NumOfNearNeighbor] = _lcalRecoPars->getIntVal("NumOfNearNeighbor");
}


/* --------------------------------------------------------------------------
   ccccccccc
   -------------------------------------------------------------------------- */
double GlobalMethodsClass::SignalGevConversion( Parameter_t optName , double valNow ){

#pragma message("FIXME: SignalToGeV conversion")
  double	returnVal = -1;
  const double conversionFactor(0.0105*1789/1500*(1488/1500.0));
  //const double conversionFactor(0.0105);

  if(optName == GeV_to_Signal)
    returnVal = valNow * conversionFactor;
  //		returnVal = valNow * .0105 + .0013;

  if(optName == Signal_to_GeV)
    returnVal = valNow / conversionFactor;
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
  double phiCell   = 2 * M_PI * (double(cellIdPhi) + .5) / double(GlobalParamI[NumCellsPhi]);

  // fill output container
  thetaPhiCell[GlobalMethodsClass::COTheta] = thetaCell;
  thetaPhiCell[GlobalMethodsClass::COPhi]   = phiCell;
  thetaPhiCell[GlobalMethodsClass::COR]     = rCell;

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
  case ThetaMin:	                 return "ThetaMin";
  case ThetaMax:	                 return "ThetaMax";
  case LogWeightConstant:                return "LogWeightConstant";
  case MoliereRadius:	                 return "MoliereRadius";
  case MinSeparationDist:                return "MinSeparationDist";
  case ElementsPercentInShowerPeakLayer: return "ElementsPercentInShowerPeakLayer";
  case NumOfNearNeighbor:                return "NumOfNearNeighbor";
  case ClusterMinNumHits:                return "ClusterMinNumHits";
  case MinHitEnergy:                     return "MinHitEnergy";
  case MinClusterEngy:	                 return "MinClusterEngy";
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
