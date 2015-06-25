#ifndef GlobalMethodsClass_H
#define GlobalMethodsClass_H 1

#include <map>
#include <string>

class GlobalMethodsClass {

public:

  enum Parameter_t{
    ZStart,
    RMin,
    RMax,
    NumCellsR,
    NumCellsPhi,
    NumCellsZ,
    RCellLength,
    ZLayerThickness,
    ThetaMin,
    ThetaMax,
    LogWeightConstant,
    MoliereRadius,
    MinSeparationDist,
    MinClusterEngy,
    GeV_to_Signal,
    Signal_to_GeV
};

  enum WeightingMethod_t{
    LogMethod=-1,
    EnergyMethod=1
  };

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
