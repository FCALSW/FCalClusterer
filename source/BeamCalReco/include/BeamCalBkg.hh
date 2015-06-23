/**
* @file BeamCalBkg.hh
* @brief File contains interface for BeamCal background
* @author Andrey Sapronov
* @version 0.0.1
* @date 2015-06-11
*/

#pragma once

#include <string>
#include <vector>

#include "BCPadEnergies.hh"

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
 
class BeamCalBkg {
 public:
  BeamCalBkg(const string &bg_method_name, const BeamCalGeo* BCG);
  virtual ~BeamCalBkg();

 public: 
  enum BackgroundMethod_t { kPregenerated, kParametrised, kAveraged };

 protected:

  BackgroundMethod_t m_bgMethod;

  int m_nBX;

  vector<double> *m_BeamCalDepositsLeft;
  vector<double> *m_BeamCalDepositsRight;
  BCPadEnergies* m_BeamCalAverageLeft;
  BCPadEnergies* m_BeamCalAverageRight;

  BCPadEnergies* m_BeamCalErrorsLeft;
  BCPadEnergies* m_BeamCalErrorsRight;

  vector<double>* m_TowerErrorsLeft;
  vector<double>* m_TowerErrorsRight;

  TRandom3 *m_random3;

  const BeamCalGeo *m_BCG;
  int m_startLayer;

 public:
  virtual void init(const int n_bx);
  virtual void init(vector<string>& bgfiles, const int n_bx) = 0;
  void setRandom3Seed(const int seed);
  void setStartLayer(const int sl) { m_startLayer = sl; }

  virtual void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight) = 0;
  virtual void getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);
  virtual void getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

//  commented out for now
//  virtual int getPadsCovariance(vector<int> &pad_list, vector<double> &covinv, 
//        const BCPadEnergies::BeamCalSide_t &bc_side) const;

  virtual int getTowerErrorsBG(int padIndex, const BCPadEnergies::BeamCalSide_t bc_side, 
        double &tower_sigma);

 protected:
  virtual void setTowerErrors(const BCPadEnergies::BeamCalSide_t bc_side);

  public:
  BeamCalBkg(const BeamCalBkg&);
  BeamCalBkg& operator=(const BeamCalBkg&);

};

