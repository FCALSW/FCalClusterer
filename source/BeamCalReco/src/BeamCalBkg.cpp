/**
* @file BeamCalBkg.cpp
* @brief Implementation for BeamCalBkg methods
* @author Andre Sailer <andre.philippe.sailer@cern.ch>
* @version 0.0.1
* @date 2015-02-18
* 
* Modified by Andrey Sapronov <andrey.sapronov@cern.ch>
* to split into abstract classes and implementations for specific
* background methods.
*
*/
#include "BeamCalBkg.hh"
#include "BCPCuts.hh"
#include "BCPadEnergies.hh"
#include "BeamCalGeo.hh"

// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

// ROOT
#include <TRandom3.h>

#include <cmath>
#include <iostream>
#include <map>

using std::vector;
using std::string;
using std::map;

BeamCalBkg::BeamCalBkg(const string& bg_method_name, const BeamCalGeo* BCG)
    : m_bgMethod(kPregenerated),
      m_nBX(0),
      m_BeamCalDepositsLeft(nullptr),
      m_BeamCalDepositsRight(nullptr),
      m_BeamCalAverageLeft(nullptr),
      m_BeamCalAverageRight(nullptr),
      m_BeamCalErrorsLeft(nullptr),
      m_BeamCalErrorsRight(nullptr),
      m_TowerErrorsLeft(nullptr),
      m_TowerErrorsRight(nullptr),
      m_random3(nullptr),
      m_BCG(BCG),
      m_bcpCuts(nullptr) {
  streamlog_out(MESSAGE) << "Initialising BeamCal background with \""
			 << bg_method_name << "\" method" << std::endl;
}

BeamCalBkg::~BeamCalBkg()
{
  delete m_random3;

  delete m_BeamCalAverageLeft;
  delete m_BeamCalAverageRight;

  delete m_BeamCalDepositsLeft;
  delete m_BeamCalDepositsRight;
  
  delete m_BeamCalErrorsLeft;
  delete m_BeamCalErrorsRight;

  delete m_TowerErrorsLeft;
  delete m_TowerErrorsRight;
}

void BeamCalBkg::init(const int n_bx)
{
  m_random3 = new TRandom3();
  m_nBX = n_bx;

  m_TowerErrorsLeft  =  new vector<double>;
  m_TowerErrorsRight =  new vector<double>;
}

/**
* @brief Calculates background stdev for energies projected along the tower
*
* For a given side it computes the projected energy error from pad errors 
* in the tower. The pad-to-pad correlatins are neglected, since they are below 0.3
* for layers deeper than 5.
*
* @param bc_side BeamCal side, Left or Right
*/
void BeamCalBkg::setTowerErrors(const BCPadEnergies::BeamCalSide_t bc_side)
{
  BCPadEnergies * BC_errors = (BCPadEnergies::kLeft == bc_side 
    ? m_BeamCalErrorsLeft : m_BeamCalErrorsRight );

  // variance of tower energies
  vector<double>* te_var = (BCPadEnergies::kLeft == bc_side 
    ? m_TowerErrorsLeft : m_TowerErrorsRight );
  te_var->clear();

  const int ppl = m_BCG->getPadsPerLayer();
  const int start_layer = m_bcpCuts->getStartingLayer();
  const int end_layer = m_bcpCuts->getStartingLayer() + m_bcpCuts->getCountingLayers();

  // loop over pads in one layer == towers in BC
  for (int ip = 0; ip < ppl; ip++){
    te_var->push_back(0.);
    // loop over pads in a tower
    for (int jp = ip+start_layer*ppl; jp < ip+(end_layer)*ppl; jp+=ppl){
      te_var->back() += pow( BC_errors->getEnergy(jp), 2);
    }
    te_var->back() = sqrt(te_var->back());
    if (te_var->back() == 0.0) {
      te_var->back() = 0.00001;
    }
  }

} // setTowerErrors

void BeamCalBkg::getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalAverageLeft);
  peRight.setEnergies(*m_BeamCalAverageRight);
}

void BeamCalBkg::getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalErrorsLeft);
  peRight.setEnergies(*m_BeamCalErrorsRight);
}

int BeamCalBkg::getTowerErrorsBG(int padIndex, 
      const BCPadEnergies::BeamCalSide_t bc_side, double &tower_sigma)
{
  tower_sigma = (BCPadEnergies::kLeft == bc_side ? m_TowerErrorsLeft->at(padIndex) 
    : m_TowerErrorsRight->at(padIndex));
  
  return tower_sigma;
}


void BeamCalBkg::setRandom3Seed(int seed)
{ 
  m_random3->SetSeed(seed); 
}
