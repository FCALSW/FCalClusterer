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

class TGeoHMatrix;


class GlobalMethodsClass {
private:
  std::string _procName;
  bool _useDD4hep;

public:
  enum WeightingMethod_t { LogMethod = -1, EnergyMethod = 1 };

  enum Parameter_t{
    ZStart,
    ZEnd,
    RMin,
    RMax,
    NumCellsR,
    NumCellsPhi,
    NumCellsZ,
    RCellLength,
    RCellOffset,
    PhiCellLength,
    ZLayerThickness,
    ZLayerPhiOffset,
    ZLayerZOffset,
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
    MiddleEnergyHitBoundFrac,
    WeightingMethod,
    Signal_to_GeV,
    BeamCrossingAngle,
    LumiInColName,
    BetaGamma,
    Gamma
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

  typedef std::map < Parameter_t, int >         ParametersInt;
  typedef std::map < Parameter_t, double >      ParametersDouble;
  typedef std::map < Parameter_t, std::string > ParametersString;

  GlobalMethodsClass();
  GlobalMethodsClass( const std::string &procType );
  GlobalMethodsClass( const GlobalMethodsClass &rhs );
  GlobalMethodsClass& operator=( const GlobalMethodsClass &rhs );

  virtual ~GlobalMethodsClass();
 
  void SetConstants( marlin::Processor* procPTR );
  WeightingMethod_t getMethod(std::string const& methodName) const;

  WeightingMethod_t method = LogMethod;

  double _backwardRotationPhi;
  ParametersInt    GlobalParamI;
  ParametersDouble GlobalParamD;
  ParametersString GlobalParamS;

  double toSignal(double valNow) const;
  double toGev(double valNow) const;

  void	ThetaPhiCell(int cellId , std::map <GlobalMethodsClass::Coordinate_t , double> & thetaPhiCell);

  static void CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int  CellIdZPR(int cellZ, int cellPhi, int cellR, int arm);
  static int  CellIdZPR(int cellId, Coordinate_t ZPR);

  void PrintAllParameters() const;

  inline bool isUsingDD4hep() const { return _useDD4hep; }

  static double posWeight(double cellEngy, double totEngy, GlobalMethodsClass::WeightingMethod_t method,
                          double logWeightConstNow);

  inline double getCalibrationFactor() const { return GlobalParamD.at(Signal_to_GeV); }

private:
  void SetGeometryGear();
  bool SetGeometryDD4HEP();

  const TGeoHMatrix *_forwardCalo;
  const TGeoHMatrix *_backwardCalo;

};



#endif
