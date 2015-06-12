/**
* @file BeamCalBkgAverage.hh
* @brief File contains manager for BeamCal averaged background
*
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

class BeamCalBkgAverage : public BeamCalBkg {
 public:
  BeamCalBkgAverage(const string &bg_method_name, const BeamCalGeo* BCG);
  ~BeamCalBkgAverage();

 public:
  void init(vector<string> &bg_files, const int n_bx);

  void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

 public:
  BeamCalBkgAverage(const BeamCalBkgAverage&);
  BeamCalBkgAverage& operator=(const BeamCalBkgAverage&);

};

