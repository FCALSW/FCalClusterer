#ifndef GlobalMethodsClass_H
#define GlobalMethodsClass_H 1

#include <map>
#include <string>
#include <tuple>

class LCCluster;

class TGeoHMatrix;

namespace marlin {
  class Processor;
}

namespace IMPL {
  class ClusterImpl;
  class ReconstructedParticleImpl;
}

class GlobalMethodsClass {
private:
  bool _useDD4hep;

public:
  enum WeightingMethod_t { LogMethod = -1, EnergyMethod = 1 };

  enum Parameter_t {
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
    PhiCellOffset,
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
  GlobalMethodsClass(const GlobalMethodsClass& rhs) = delete;
  GlobalMethodsClass& operator=(const GlobalMethodsClass& rhs) = default;

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

  void ThetaPhiCell(const int cellId, std::map<GlobalMethodsClass::Coordinate_t, double>& thetaPhiCell) const;

  static void CellIdZPR(int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int  CellIdZPR(int cellZ, int cellPhi, int cellR, int arm);
  static int  CellIdZPR(int cellId, Coordinate_t ZPR);

  void PrintAllParameters() const;
  void initializeAdditionalParameters();

  inline bool isUsingDD4hep() const { return _useDD4hep; }

  static double posWeight(double cellEngy, double totEngy, GlobalMethodsClass::WeightingMethod_t method,
                          double logWeightConstNow);

  inline double getCalibrationFactor() const { return GlobalParamD.at(Signal_to_GeV); }

  template <class T, class U> inline void rotateToGlobal(const T* loc, U* glob) const;
  template <class T, class U> inline void rotateToLumiCal(const T* glob, U* loc) const;

  std::tuple<IMPL::ClusterImpl*, IMPL::ReconstructedParticleImpl*> getLCIOObjects(LCCluster const& thisClusterInfo,
                                                                                  double minClusterEnergy,
                                                                                  bool   cutOnFiducialVolume) const;

private:
  void SetGeometryGear();
  bool SetGeometryDD4HEP();

  const TGeoHMatrix *_forwardCalo;
  const TGeoHMatrix *_backwardCalo;

  // LumiCal rotations angles ( local->Global )
  std::map<int, double> m_armCosAngle{};
  std::map<int, double> m_armSinAngle{};
};

template <class T, class U> void GlobalMethodsClass::rotateToGlobal(const T* loc, U* glob) const {
  const int armNow = (loc[2] < 0) ? -1 : 1;
  glob[0]          = +m_armCosAngle.at(armNow) * loc[0] + m_armSinAngle.at(armNow) * loc[2];
  glob[1]          = loc[1];
  glob[2]          = -m_armSinAngle.at(armNow) * loc[0] + m_armCosAngle.at(armNow) * loc[2];
}

template <class T, class U> void GlobalMethodsClass::rotateToLumiCal(const T* glob, U* loc) const {
  const int armNow = (glob[2] < 0) ? -1 : 1;
  loc[0]           = +m_armCosAngle.at(armNow) * glob[0] - m_armSinAngle.at(armNow) * glob[2];
  loc[1]           = glob[1];
  loc[2]           = +m_armSinAngle.at(armNow) * glob[0] + m_armCosAngle.at(armNow) * glob[2];
}

#endif
