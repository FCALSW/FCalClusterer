/**
* @file BeamCalBackground.hh
* @brief File contains manager for BeamCal background
* @author Andrey Sapronov
* @version 0.0.1
* @date 2015-02-18
*/

#pragma once

#include <string>
#include <vector>

#include "BCPadEnergies.hh"

class TChain;
class TFile;
class TH1;
class TRandom3;
class TTree;

class BeamCalGeo;

using std::vector;
using std::string;

typedef struct {
  double zero_rate;
  double mean;
  double stdev;
  double sum;
  double minm;
  double maxm;
  double chi2;
  double par0;
  double par1;
  double par2;
} PadEdepRndPar_t;
 
class BeamCalBackground {
 public:
  BeamCalBackground(const string &bg_method_name, const BeamCalGeo* BCG);
  ~BeamCalBackground();

 public: 
  enum BackgroundMethod_t { kPregenerated, kParametrised, kAveraged };

 private:

  BackgroundMethod_t m_bgMethod;

  int m_nBX;

  vector<double> *m_BeamCalDepositsLeft;
  vector<double> *m_BeamCalDepositsRight;
  BCPadEnergies* m_BeamCalAverageLeft;
  BCPadEnergies* m_BeamCalAverageRight;

  BCPadEnergies* m_BeamCalErrorsLeft;
  BCPadEnergies* m_BeamCalErrorsRight;

  vector<BCPadEnergies*> m_listOfBunchCrossingsLeft;
  vector<BCPadEnergies*> m_listOfBunchCrossingsRight;

  vector<double>* m_TowerErrorsLeft;
  vector<double>* m_TowerErrorsRight;

  TRandom3 *m_random3;
  TChain* m_backgroundBX;

  vector<PadEdepRndPar_t> *m_padParLeft;
  vector<PadEdepRndPar_t> *m_padParRight;

  const BeamCalGeo *m_BCG;
  int m_startLayer;
  int m_numberForAverage;

 public:
  int init(vector<string> &bg_files, const int n_bx);
  void setRandom3Seed(const int seed);
  void setStartLayer(const int sl) { m_startLayer = sl; }
  void setNumberForAverage(const int nav) { m_numberForAverage = nav; }

  void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  void getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  void getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

  int getTowerErrorsBG(int padIndex, const BCPadEnergies::BeamCalSide_t bc_side, 
        double &tower_sigma);

  int getPadsCovariance(vector<int> &pad_list, vector<double> &covinv, 
        const BCPadEnergies::BeamCalSide_t &bc_side) const;

 private:
  int initPregenerated(vector<string> &bg_files);
  int initParametrised(vector<string> &bg_files);
  int initAveraged(vector<string> &bg_files);

  void getEventPregeneratedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  void getEventParametrisedBG(BCPadEnergies &pe,
        const BCPadEnergies::BeamCalSide_t bc_side);
  void getEventAveragedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

  void setTowerErrors(const std::vector<BCPadEnergies*> singles, 
         const BCPadEnergies::BeamCalSide_t bc_side);

  BCPadEnergies* getBeamCalErrors(const BCPadEnergies *averages, 
                   const std::vector<BCPadEnergies*> singles );

  void readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side);

 public:
  BeamCalBackground(const BeamCalBackground&);
  BeamCalBackground& operator=(const BeamCalBackground&);

};

