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
  enum BackgroundMethod_t { kPregenerated, kParametrised };

 private:

  BackgroundMethod_t m_bgMethod;

  int m_nBX;

  vector<double> *m_BeamCalDepositsLeft;
  vector<double> *m_BeamCalDepositsRight;
  BCPadEnergies* m_BeamCalAverageLeft;
  BCPadEnergies* m_BeamCalAverageRight;

  BCPadEnergies* m_BeamCalErrorsLeft;
  BCPadEnergies* m_BeamCalErrorsRight;

  TRandom3 *m_random3;
  TChain* m_backgroundBX;

  vector<PadEdepRndPar_t> *m_padParLeft;
  vector<PadEdepRndPar_t> *m_padParRight;

  const BeamCalGeo *m_BCG;

 public:
  int init(vector<string> &bg_files, const int n_bx);
  void setRandom3Seed(const int seed);
  int getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  int getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  int getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

 private:
  int initPregenerated(vector<string> &bg_files);
  int initParametrised(vector<string> &bg_files);

  int getEventPregeneratedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  int getEventParametrisedBG(BCPadEnergies &pe, 
        const BCPadEnergies::BeamCalSide_t bc_side);

  BCPadEnergies* getBeamCalErrors(const BCPadEnergies *averages, 
                   const std::vector<BCPadEnergies*> singles, int numberForAverage );

  int readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side);

 public:
  BeamCalBackground(const BeamCalBackground&);
  BeamCalBackground& operator=(const BeamCalBackground&);

};

