/**
* @file BeamCalBkgPregen.cpp
* @brief Implementation for BeamCalBkgPregen methods
* @author Andre Sailer <andre.philippe.sailer@cern.ch>
* @version 0.0.1
* @date 2015-02-18
* 
* Modified by Andrey Sapronov <andrey.sapronov@cern.ch>
* to split into separate class for pregenerated background
*
*/
#include "BeamCalBkgPregen.hh"
#include "BCPadEnergies.hh"
#include "BeamCalBkg.hh"
#include "BeamCalGeo.hh"

// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

// ROOT
#include <TChain.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

using std::vector;
using std::string;
using std::map;

/// A UniformRandomBitGenerator based on TRandom3
class Random3_URBG {
private:
  TRandom3 * const r3;
public:
  using result_type = UInt_t;
  Random3_URBG() = delete;
  explicit Random3_URBG(TRandom3 * const random3): r3(random3) {}
  UInt_t operator()() { return r3->Integer(kMaxUInt); }
  constexpr UInt_t min() { return 0; }
  constexpr UInt_t max() { return kMaxUInt-1; }
};

BeamCalBkgPregen::BeamCalBkgPregen(const string& bg_method_name, const BeamCalGeo* BCG)
    : BeamCalBkg(bg_method_name, BCG),
      m_listOfBunchCrossingsLeft(vector<BCPadEnergies*>()),
      m_listOfBunchCrossingsRight(vector<BCPadEnergies*>()),
      m_backgroundBX(nullptr),
      m_numberForAverage(1) {
  streamlog_out(MESSAGE) << "Initialising BeamCal background with \""
			 << bg_method_name << "\" method" << std::endl;
}

BeamCalBkgPregen::~BeamCalBkgPregen()
{
  delete m_backgroundBX;

  vector<BCPadEnergies*>::iterator it_pe = m_listOfBunchCrossingsRight.begin();
  for(; it_pe != m_listOfBunchCrossingsRight.end(); it_pe++)
    delete *it_pe;

  it_pe = m_listOfBunchCrossingsLeft.begin();
  for(; it_pe != m_listOfBunchCrossingsLeft.end(); it_pe++)
    delete *it_pe;
}

void BeamCalBkgPregen::init(vector<string> &bg_files, const int n_bx)
{
  this->BeamCalBkg::init(n_bx);

  m_numberForAverage = 10;

  //Open the Files given as the list into a TChain...
  m_backgroundBX = new TChain("bcTree");

  //mix up the files, because the random numbers are ordered to avoid repeating
  std::shuffle(bg_files.begin(), bg_files.end(), Random3_URBG(m_random3));

  for (std::vector<std::string>::iterator file = bg_files.begin(); file != bg_files.end(); ++file) {
    streamlog_out(DEBUG1) << *file << std::endl;
    m_backgroundBX->Add((*file).c_str());
  }

  //Ready the energy deposit vectors for the tree
  m_BeamCalDepositsLeft  = nullptr;
  m_BeamCalDepositsRight = nullptr;

  m_backgroundBX->SetBranchAddress("vec_left" , &m_BeamCalDepositsLeft);
  m_backgroundBX->SetBranchAddress("vec_right", &m_BeamCalDepositsRight);

  streamlog_out(DEBUG2) << "We have " << m_backgroundBX->GetEntries() << " background BXs" << std::endl;

  m_BeamCalAverageLeft  =  new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight =  new BCPadEnergies(m_BCG);

  // Create the average BCPadEnergies objects used to subtract average backgrounds
  m_listOfBunchCrossingsRight.clear();
  m_listOfBunchCrossingsLeft.clear();

  for (int i = 0; i < 10; ++i) {
    m_listOfBunchCrossingsLeft.push_back( new BCPadEnergies(m_BCG) );
    m_listOfBunchCrossingsRight.push_back( new BCPadEnergies(m_BCG) );
  }

  std::set<int> randomNumbers;
  const unsigned int nBackgroundBX = m_backgroundBX->GetEntries();

  //Check that we have
  if( int(nBackgroundBX) < m_nBX*10 ) {
    streamlog_out(ERROR) << "There are not enough BeamCal " \
     " Background files to calculate a proper average!" << std::endl;
    throw std::runtime_error( "Not enough BeamCal Background bunch crossings available");
  }

  //we use a set so no duplication occurs
  while( int(randomNumbers.size()) < m_nBX*m_numberForAverage ){//do it ten times as often
    randomNumbers.insert( int(m_random3->Uniform(0, nBackgroundBX)) );
  }

  int counter = 0;

  for (std::set<int>::iterator it = randomNumbers.begin(); it != randomNumbers.end();++it) {
    streamlog_out(DEBUG1) << std::setw(5) << *it << std::flush;
    m_backgroundBX->GetEntry(*it);
    m_BeamCalAverageLeft->addEnergies(*m_BeamCalDepositsLeft);
    m_BeamCalAverageRight->addEnergies(*m_BeamCalDepositsRight);
    m_listOfBunchCrossingsLeft.at(counter/m_nBX)-> addEnergies(*m_BeamCalDepositsLeft);
    m_listOfBunchCrossingsRight.at(counter/m_nBX)-> addEnergies(*m_BeamCalDepositsRight);
    ++counter;
  }

  //Now divide by ten, to get the average distributions...
  m_BeamCalAverageLeft ->scaleEnergies(1./double(m_numberForAverage));
  m_BeamCalAverageRight->scaleEnergies(1./double(m_numberForAverage));
  double totalEnergyMean = m_BeamCalAverageLeft->getTotalEnergy();
  double varEn(0.0);
  for (int l = 0; l < m_numberForAverage;++l) {
    varEn += (m_listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean) * (m_listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean);
  }//histograms
  varEn /= double(m_numberForAverage);
  varEn = sqrt(varEn);
  streamlog_out(MESSAGE4) << "Total Energy " << totalEnergyMean << " +- "
			  << varEn << " GeV/" << m_nBX <<"BX" << std::endl;

  // m_h3BeamCalAverageLeft = bc.getBeamCalHistogram("AverageLeft");
  // m_h3BeamCalAverageRight = bc.getBeamCalHistogram("AverageRight");

  //And now calculate the error for every bin....
  m_BeamCalErrorsLeft  = getBeamCalErrors(m_BeamCalAverageLeft,  m_listOfBunchCrossingsLeft);
  m_BeamCalErrorsRight = getBeamCalErrors(m_BeamCalAverageRight, m_listOfBunchCrossingsRight);

  //Add one sigma to the averages -- > just do it once here
  //m_BeamCalAverageLeft ->addEnergies( m_BeamCalErrorsLeft );
  //m_BeamCalAverageRight->addEnergies( m_BeamCalErrorsRight);

  // calculate st.dev. of tower energies
  this->setTowerErrors(BCPadEnergies::kLeft);
  this->setTowerErrors(BCPadEnergies::kRight);

  streamlog_out(DEBUG1) << std::endl;

  /*
  for (int i = 0; i < m_numberForAverage;++i) {
    delete m_listOfBunchCrossingsLeft[i];
    delete m_listOfBunchCrossingsRight[i];
  }
  */
}


/// Calculate the one sigma errors for the given selected background bunch crossings
BCPadEnergies* BeamCalBkgPregen::getBeamCalErrors(const BCPadEnergies *averages, const std::vector<BCPadEnergies*> singles) {

  BCPadEnergies * BCPErrors = new BCPadEnergies(m_BCG);
  for (int i = 0; i < m_BCG->getPadsPerBeamCal()  ;++i) {
    double mean(averages->getEnergy(i));
    double variance(0);
    double nHistos(m_numberForAverage);
    for (int l = 0; l < nHistos;++l) {
      double energy(singles[l]->getEnergy(i));
      variance += ( energy - mean ) * ( energy - mean );
    }//histograms
    variance /= nHistos;
    variance = sqrt(variance);
    BCPErrors->setEnergy(i, variance);

  }//for all pads

  return BCPErrors;

}//getBeamCalErrors


/*
int BeamCalBkgPregen::getPadsCovariance(vector<int> &pad_list, vector<double> &covinv, 
      const BCPadEnergies::BeamCalSide_t &bc_side) const
{
  // local copy of bunchcrossings
  const vector<BCPadEnergies*> bxl = (BCPadEnergies::kLeft == bc_side 
    ? m_listOfBunchCrossingsLeft : m_listOfBunchCrossingsRight );
  
  // average Edep
  BCPadEnergies * pe_aver = (BCPadEnergies::kLeft == bc_side 
    ? m_BeamCalAverageLeft : m_BeamCalAverageRight );

  // number of elements
  const int nm = pad_list.size();

  // new covariance matrix, to be inverted and deleted
  double *ma = new double [nm*nm];

  vector<double> te_aver;
  for (int i = 0; i<nm; i++){
    te_aver.push_back(pe_aver->getTowerEnergy(pad_list[i],m_startLayer));
  }

  for (int i = 0; i<nm; i++){
    for (int j = 0; j<nm; j++){
      ma[i*nm+j] = 0.;
      // loop over averaged entries to calculate covariance
      for (int ibx = 0; ibx < m_numberForAverage; ibx++){
      
        ma[i*nm+j] += (bxl.at(ibx)->getTowerEnergy(pad_list[i], m_startLayer) - te_aver[i])
	         *(bxl.at(ibx)->getTowerEnergy(pad_list[j],m_startLayer) - te_aver[j])
		 /(m_numberForAverage-1);
      }
    }
  }

  double det(0.);
  TMatrixD *mtx= new TMatrixD(nm, nm, ma); 
  TMatrixD mtxi(mtx->InvertFast(&det));
  if ( 0. == det) {
    streamlog_out(WARNING5) << "Unable to calculate inverse of covariance matrix" \
         " for requested pads" << std::endl;;

    delete mtx;

    return -1;
  }

  double *mtxarr = mtxi.GetMatrixArray();
  covinv.assign(mtxarr, mtxarr+nm*nm);

  delete mtx;
  delete ma;

  return covinv.size();
}
*/

void BeamCalBkgPregen::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  ////////////////////////////////////////////////////////
  // Prepare the randomly chosen Background BeamCals... //
  ////////////////////////////////////////////////////////
  std::set<int> randomNumbers;
  unsigned int nBackgroundBX = m_backgroundBX->GetEntries();
  while( int(randomNumbers.size()) < m_nBX ){
    randomNumbers.insert( int(m_random3->Uniform(0, nBackgroundBX)) );
  }

  ////////////////////////
  // Sum them all up... //
  ////////////////////////
  for (std::set<int>::iterator it = randomNumbers.begin(); it != randomNumbers.end();++it) {
    m_backgroundBX->GetEntry(*it);
    peRight.addEnergies(*m_BeamCalDepositsRight);
    peLeft.addEnergies(*m_BeamCalDepositsLeft);
  }

}
