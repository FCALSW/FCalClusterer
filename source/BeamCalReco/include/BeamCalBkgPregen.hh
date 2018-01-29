/**
* @file BeamCalBkgPregen.hh
* @brief File contains manager for BeamCal pregenerated background
* @author Andrey Sapronov
* @version 0.0.1
* @date 2015-02-18
*/

#pragma once

#include <string>
#include <vector>

#include "BeamCalBkg.hh"

class TChain;

class BCPadEnergies;
class BeamCalGeo;

using std::vector;
using std::string;

class BeamCalBkgPregen : public BeamCalBkg {
 public:
  BeamCalBkgPregen(const string &bg_method_name, const BeamCalGeo* BCG);
  ~BeamCalBkgPregen();

 private:
  vector<BCPadEnergies*> m_listOfBunchCrossingsLeft;
  vector<BCPadEnergies*> m_listOfBunchCrossingsRight;

  TChain* m_backgroundBX;

  int m_numberForAverage;

 public:
  void init(vector<string> &bg_files, const int n_bx);
  void setNumberForAverage(const int nav) { m_numberForAverage = nav; }

  void getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight);

//  commented out for now
//  int getPadsCovariance(vector<int> &pad_list, vector<double> &covinv, 
//        const BCPadEnergies::BeamCalSide_t &bc_side) const;

 private:
  BCPadEnergies* getBeamCalErrors(const BCPadEnergies *averages, 
                   const std::vector<BCPadEnergies*> singles );

 public:
  BeamCalBkgPregen(const BeamCalBkgPregen&);
  BeamCalBkgPregen& operator=(const BeamCalBkgPregen&);

};

