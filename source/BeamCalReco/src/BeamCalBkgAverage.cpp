/**
* @file BeamCalBkgAverage.cpp
* @brief Implementation for BeamCalBkgAverage methods
* @author Andre Sailer <andre.philippe.sailer@cern.ch>
* @version 0.0.1
* @date 2015-02-18
* 
* Modified by Andrey Sapronov <andrey.sapronov@cern.ch>
* to add method for parametrised background.
*
*/
#include "BeamCalBkgAverage.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"
#include "BeamCalGeo.hh"

#include <TRandom3.h>

#include <cmath>

using std::vector;
using std::string;


BeamCalBkgAverage::BeamCalBkgAverage(const string& bg_method_name, 
                     const BeamCalGeo *BCG) : BeamCalBkg(bg_method_name, BCG)
{}

BeamCalBkgAverage::~BeamCalBkgAverage()
{}

void BeamCalBkgAverage::init(vector<string> &bg_files, const int n_bx)
{
  this->BeamCalBkg::init(n_bx);

  //We don't have any average energy in some cases, and we don't really need it
  //because we just subtract it after adding it in some cases anyway....
  m_BeamCalAverageLeft  = new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight = new BCPadEnergies(m_BCG);

  std::vector<BCPadEnergies> backgroundBeamCals(2, m_BCG);

  BCUtil::ReadBecasFile(bg_files[0], backgroundBeamCals, 
    "tBcDensAverage", "sEdepErr", true);

  m_BeamCalErrorsLeft  = new BCPadEnergies(backgroundBeamCals[0]);
  m_BeamCalErrorsRight = new BCPadEnergies(backgroundBeamCals[1]);

  // calculate st.dev. of tower energies
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kLeft);
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kRight);
}


void BeamCalBkgAverage::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  const int nBCpads = m_BCG->getPadsPerBeamCal();

  // generate at once with sigma = dE*sqrts(m_nBX)
  const double rd_coef = sqrt(m_nBX);
  for (int i = 0; i < nBCpads ;++i) { //Add gaussian randomisation of background to each cell
    peRight.addEnergy(i, m_random3->Gaus(0.0, rd_coef*m_BeamCalErrorsRight->getEnergy(i)));
    peLeft.addEnergy(i, m_random3->Gaus(0.0, rd_coef*m_BeamCalErrorsLeft->getEnergy(i)));
  }//for all pads
} // getEventBG
