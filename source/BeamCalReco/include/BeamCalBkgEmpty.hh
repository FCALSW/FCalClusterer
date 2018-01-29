/**
* @file BeamCalBkgEmpty.hh
* @brief File contains manager for BeamCal empty background if no background is used
*
*/

#pragma once

#include "BeamCalBkg.hh"

#include <string>
#include <vector>

class BeamCalGeo;
class BCPadEnergies;

class BeamCalBkgEmpty : public BeamCalBkg {
 public:
  BeamCalBkgEmpty(const string &bg_method_name, const BeamCalGeo* BCG);
  ~BeamCalBkgEmpty();
  BeamCalBkgEmpty(const BeamCalBkgEmpty&) = delete;
  BeamCalBkgEmpty& operator=(const BeamCalBkgEmpty&) = delete;

 public:
  void init(vector<string> &bg_files, const int n_bx);

  void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);


};

