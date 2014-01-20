#include "GlobalMethodsClass.h"

#include <marlin/VerbosityLevels.h>

#include <map>
#include <string>
#include <cmath>
#include <iostream>


GlobalMethodsClass :: GlobalMethodsClass() :
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
int GlobalMethodsClass::CellIdZPR(int cellZ, int cellPhi, int cellR) {

  int cellId = 0;

  cellId |= ( cellZ   << 0  ) ;
  cellId |= ( cellPhi << 8  ) ;
  cellId |= ( cellR   << 16 ) ;

  return cellId;
}

int GlobalMethodsClass::CellIdZPR(int cellId, GlobalMethodsClass::Coordinate_t ZPR) {

  int cellZ, cellPhi, cellR;
  CellIdZPR(cellId, cellZ, cellPhi, cellR);

  if(ZPR == GlobalMethodsClass::COZ) return cellZ;
  if(ZPR == GlobalMethodsClass::COR) return cellR;
  if(ZPR == GlobalMethodsClass::COP) return cellPhi;

  return 0;
}

void GlobalMethodsClass::CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR) {

  // compute Z,Phi,R coordinates according to the cellId
  cellZ   = (cellId >> 0 ) & 0xff ;
  cellPhi = (cellId >> 8 ) & 0xff ;
  cellR   = (cellId >> 16 ) & 0xff ;
  return;

}

void GlobalMethodsClass::SetConstants() {
#pragma message ("FIXME: Pick geometry up from GearFile")
  // fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  GlobalParamD[ThetaMin] = 40e-3;
  GlobalParamD[ThetaMax] = 200e-3;///APS this is large, because we have the crossing angle not taken into account

  // geometry parameters of LumiCal

  // starting position [mm]
  GlobalParamD[ZStart]   = 2654.2;
  // inner and outer radii [mm]
  GlobalParamD[RMin] = 100;
  GlobalParamD[RMax] = 290;
  // cell division numbers
  GlobalParamI[NumCellsR] = 64;
  GlobalParamI[NumCellsPhi] = 48;
  GlobalParamI[NumCellsZ] = 40;

  GlobalParamD[RCellLength] = (GlobalParamD[RMax] - GlobalParamD[RMin]) / double(GlobalParamI[NumCellsR]);
  // layer thickness
  GlobalParamD[ZLayerThickness] = 3.5+0.2+0.32+0.25;

  // logarithmic constant for position reconstruction
  GlobalParamD[LogWeightConstant] = 6.;

  // Moliere radius of LumiCal [mm]
  GlobalParamD[MoliereRadius] = 14;

  // minimal separation distance between any pair of clusters [mm]
  GlobalParamD[MinSeparationDist] = GlobalParamD[MoliereRadius];

  // minimal energy of a single cluster
  GlobalParamD[MinClusterEngy] = 2;  // value in GeV
  GlobalParamD[MinClusterEngy] = SignalGevConversion(GeV_to_Signal , GlobalParamD[MinClusterEngy]);  // conversion to "detector Signal"


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
  int cellIdZ, cellIdPhi, cellIdR;
  CellIdZPR(cellId, cellIdZ, cellIdPhi, cellIdR);

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
  case ZStart:             return "ZStart";
  case RMin:		   return "RMin";
  case RMax:		   return "RMax";
  case NumCellsR:	   return "NumCellsR";
  case NumCellsPhi:	   return "NumCellsPhi";
  case NumCellsZ:	   return "NumCellsZ";
  case RCellLength:	   return "RCellLength";
  case ZLayerThickness:	   return "ZLayerThickness";
  case ThetaMin:	   return "ThetaMin";
  case ThetaMax:	   return "ThetaMax";
  case LogWeightConstant:  return "LogWeightConstant";
  case MoliereRadius:	   return "MoliereRadius";
  case MinSeparationDist:  return "MinSeparationDist";
  case MinClusterEngy:	   return "MinClusterEngy";
  case GeV_to_Signal:	   return "GeV_to_Signal";
  case Signal_to_GeV:      return "Signal_to_GeV";
  default: return "Unknown Parameter";
  }


}


void GlobalMethodsClass::PrintAllParameters() const {


  for (ParametersInt::const_iterator it = GlobalParamI.begin();it != GlobalParamI.end() ;++it) {
    streamlog_out(MESSAGE) << " - (int)     " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

  for (ParametersDouble::const_iterator it = GlobalParamD.begin();it != GlobalParamD.end() ;++it) {
    streamlog_out(MESSAGE) << " - (double)  " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

  for (ParametersString::const_iterator it = GlobalParamS.begin();it != GlobalParamS.end() ;++it) {
    streamlog_out(MESSAGE) << " - (string)  " << GetParameterName(it->first) << "  =  " << it->second<< std::endl;
  }

}
