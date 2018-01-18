/**
* @file BeamCalBkgEmpty.cpp
* @brief Implementation for BeamCalBkgEmpty methods
*/
#include "BeamCalBkgEmpty.hh"
#include "BCPadEnergies.hh"

BeamCalBkgEmpty::BeamCalBkgEmpty(const string& bg_method_name, const BeamCalGeo* BCG) : BeamCalBkg(bg_method_name, BCG) {}

BeamCalBkgEmpty::~BeamCalBkgEmpty() {}

void BeamCalBkgEmpty::init(vector<string>&, const int) {
  this->BeamCalBkg::init(0);

  m_BeamCalAverageLeft  = new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight = new BCPadEnergies(m_BCG);
  m_BeamCalErrorsLeft   = new BCPadEnergies(m_BCG);
  m_BeamCalErrorsRight  = new BCPadEnergies(m_BCG);

  return;
}

void BeamCalBkgEmpty::getEventBG(BCPadEnergies&, BCPadEnergies&) { return; }
