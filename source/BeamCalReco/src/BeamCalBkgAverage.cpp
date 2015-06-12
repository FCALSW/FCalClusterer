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
#include "BeamCalGeoCached.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"


// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

// ROOT
#include <TChain.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom3.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

using std::vector;
using std::string;
using std::map;

using marlin::Global;

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
