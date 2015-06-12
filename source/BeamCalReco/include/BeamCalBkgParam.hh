/**
* @file BeamCalBkgParam.hh
* @brief File contains manager for BeamCal parametrised background
* @author Andrey Sapronov
* @version 0.0.1
* @date 2015-02-18
*/

#pragma once

#include <string>
#include <vector>

#include "BeamCalBkg.hh"
#include "BCPadEnergies.hh"

class TFile;
class TRandom3;
class TTree;

class BeamCalGeo;

using std::vector;
using std::string;

class BeamCalBkgParam : public BeamCalBkg {
 public:
  BeamCalBkgParam(const string &bg_method_name, const BeamCalGeo* BCG);
  ~BeamCalBkgParam();

 private:
  vector<PadEdepRndPar_t> *m_padParLeft;
  vector<PadEdepRndPar_t> *m_padParRight;

 public:
  void init(vector<string> &bg_files, const int n_bx);

  void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

 private:
  void setTowerErrors(const std::vector<BCPadEnergies*> singles, 
         const BCPadEnergies::BeamCalSide_t bc_side);

  void readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side);

 public:
  BeamCalBkgParam(const BeamCalBkgParam&);
  BeamCalBkgParam& operator=(const BeamCalBkgParam&);

};

