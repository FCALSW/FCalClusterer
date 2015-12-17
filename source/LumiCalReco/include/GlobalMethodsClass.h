#ifndef GlobalMethodsClass_H
#define GlobalMethodsClass_H 1

#include "Global.hh"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <gear/GearMgr.h>

#include <marlin/ProcessorMgr.h>
#include <marlin/Processor.h>
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
    ThetaMin,
    ThetaMax,
    LogWeightConstant,
    MoliereRadius,
    MinSeparationDist,
    ElementsPercentInShowerPeakLayer,
    NumOfNearNeighbor,
    ClusterMinNumHits,
    MinHitEnergy,
    MinClusterEngy,
    MiddleEnergyHitBoundFrac,
    WeightingMethod,
    GeV_to_Signal,
    Signal_to_GeV,
    BeamCrossingAngle
};
  typedef short WeightingMethod_t;
  static const WeightingMethod_t LogMethod = -1;
  static const WeightingMethod_t EnergyMethod = 1;
  /*  
  enum {
    LogMethod=-1,
    EnergyMethod=1
  };
  */
  enum Coordinate_t {
    COTheta,
    COPhi,
    COZ,
    COR,
    COP,
    COA
  }; 
  
  static std::string GetParameterName ( Parameter_t par );

  typedef std::map < Parameter_t, int >         ParametersInt;
  typedef std::map < Parameter_t, double >      ParametersDouble;
  typedef std::map < Parameter_t, std::string > ParametersString;

  GlobalMethodsClass();
  GlobalMethodsClass( const std::string &procType );
  ~GlobalMethodsClass();
 
  void SetConstants();

  ParametersInt    GlobalParamI;
  ParametersDouble GlobalParamD;
  ParametersString GlobalParamS;

  static double SignalGevConversion( Parameter_t optName , double valNow );

  void	ThetaPhiCell(int cellId , std::map <GlobalMethodsClass::Coordinate_t , double> & thetaPhiCell);

  static void CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int  CellIdZPR(int cellZ, int cellPhi, int cellR, int arm);
  static int  CellIdZPR(int cellId, Coordinate_t ZPR);

  void PrintAllParameters() const;

};



#endif
