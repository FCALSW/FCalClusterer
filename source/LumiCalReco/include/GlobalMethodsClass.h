#ifndef GlobalMethodsClass_H
#define GlobalMethodsClass_H 1

#include "Global.hh"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <gear/GearMgr.h>

#include <marlin/ProcessorMgr.h>
#include <marlin/CMProcessor.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorParameter.h>
#include <marlin/StringParameters.h>


#include <map>
#include <string>


class GlobalMethodsClass {

 private:
  std::string _procName;
 public:

  enum Parameter_t{
    ZStart,
    ZEnd,
    RMin,
    RMax,
    NumCellsR,
    NumCellsPhi,
    NumCellsZ,
    RCellLength,
    PhiCellLength,
    ZLayerThickness,
    ZLayerPhiOffset,
    ThetaMin,
    ThetaMax,
    LogWeightConstant,
    MoliereRadius,
    MinSeparationDist,
    ElementsPercentInShowerPeakLayer,
    NumOfNearNeighbor,
    ClusterMinNumHits,
    MinHitEnergy,
    MinClusterEngyGeV,
    MinClusterEngySignal,
    MiddleEnergyHitBoundFrac,
    WeightingMethod,
    GeV_to_Signal,
    Signal_to_GeV,
    BeamCrossingAngle,
    LumiInColName
};
  //
  
 enum Coordinate_t {
    COTheta,
    COPhi,
    COZ,
    COR,
    COP,
    COA
  }; 
  
  static std::string GetParameterName ( Parameter_t par );

  typedef std::string                           WeightingMethod_t;
  typedef std::map < Parameter_t, int >         ParametersInt;
  typedef std::map < Parameter_t, double >      ParametersDouble;
  typedef std::map < Parameter_t, std::string > ParametersString;

  GlobalMethodsClass();
  GlobalMethodsClass( const std::string &procType );
  virtual ~GlobalMethodsClass();
 
  void SetConstants();
  static WeightingMethod_t LogMethod;
  static WeightingMethod_t EnergyMethod;
  static double EnergyCalibrationFactor;
  ParametersInt    GlobalParamI;
  ParametersDouble GlobalParamD;
  ParametersString GlobalParamS;

  double SignalGevConversion( Parameter_t optName , double valNow );
  void	ThetaPhiCell(int cellId , std::map <GlobalMethodsClass::Coordinate_t , double> & thetaPhiCell);

  static void CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int  CellIdZPR(int cellZ, int cellPhi, int cellR, int arm);
  static int  CellIdZPR(int cellId, Coordinate_t ZPR);

  void PrintAllParameters() const;

};



#endif
