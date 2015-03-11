/**
* @file BeamCalBackground.cpp
* @brief Implementation for BeamCalBackground methods
* @author Andre Sailer <andre.philippe.sailer@cern.ch>
* @version 0.0.1
* @date 2015-02-18
* 
* Modified by Andrey Sapronov <andrey.sapronov@cern.ch>
* to add method for parametrised background.
*
*/
#include "BeamCalBackground.hh"
#include "BeamCalGeoCached.hh"
#include "BCPadEnergies.hh"

// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

// ROOT
#include <TChain.h>
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
using std::endl;

using marlin::Global;

BeamCalBackground::BeamCalBackground(const string& bg_method_name, 
                     const BeamCalGeo *BCG) : 
                                           m_bgMethod(kParametrised),
					   m_nBX(0),
                                           m_BeamCalDepositsLeft(NULL),
                                           m_BeamCalDepositsRight(NULL),
                                           m_BeamCalAverageLeft(NULL),
                                           m_BeamCalAverageRight(NULL),
                                           m_BeamCalErrorsLeft(NULL),
                                           m_BeamCalErrorsRight(NULL),
                                           m_random3(NULL),
					   m_backgroundBX(NULL),
					   m_padParLeft(NULL),
					   m_padParRight(NULL),
                                           m_BCG(BCG)
{
  if( string("Pregenerated") == bg_method_name ) m_bgMethod = kPregenerated;
  else if( string("Parametrised") == bg_method_name ) m_bgMethod = kParametrised;
  else {
    streamlog_out(ERROR) <<"== Error From BeamCalBackground == " \
            "Unknown BeamCal method of background esitmation \""+bg_method_name+"\"!" << endl;;
  }

  streamlog_out(MESSAGE) << "initialising BeamCal background "\
    " with \"" << bg_method_name << "\" method" << endl;
}

BeamCalBackground::~BeamCalBackground()
{
  delete m_random3;
  if (kPregenerated == m_bgMethod ) delete m_backgroundBX;

  delete m_BeamCalAverageLeft;
  delete m_BeamCalAverageRight;

  delete m_BeamCalDepositsLeft;
  delete m_BeamCalDepositsRight;
  
  delete m_BeamCalErrorsLeft;
  delete m_BeamCalErrorsRight;
}

int BeamCalBackground::init(vector<string> &bg_files, const int n_bx)
{
  m_random3 = new TRandom3();
  m_nBX = n_bx;

  switch (m_bgMethod) {
    case kPregenerated     : initPregenerated(bg_files); break;
    case kParametrised     : initParametrised(bg_files); break;
  }

  return 0;
}


int BeamCalBackground::initPregenerated(vector<string> &bg_files)
{

  //Open the Files given as the list into a TChain...
  m_backgroundBX = new TChain("bcTree");

  //mix up the files, because the random numbers are ordered to avoid repeating
  std::random_shuffle(bg_files.begin(), bg_files.end());

  for (std::vector<std::string>::iterator file = bg_files.begin(); file != bg_files.end(); ++file) {
    streamlog_out(DEBUG1) << *file << std::endl;
    m_backgroundBX->Add(TString(*file));
  }

  //Ready the energy deposit vectors for the tree
  m_BeamCalDepositsLeft=NULL;
  m_BeamCalDepositsRight=NULL;

  m_backgroundBX->SetBranchAddress("vec_left" , &m_BeamCalDepositsLeft);
  m_backgroundBX->SetBranchAddress("vec_right", &m_BeamCalDepositsRight);

  streamlog_out(DEBUG2) << "We have " << m_backgroundBX->GetEntries() << " background BXs" << std::endl;

  //Create an Average BeamCal, needed to setup the BCPadEnergies
  m_BCG = new BeamCalGeoCached(marlin::Global::GEAR);

  m_BeamCalAverageLeft  =  new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight =  new BCPadEnergies(m_BCG);

  // Create the average BCPadEnergies objects used to subtract average backgrounds
  std::vector<BCPadEnergies*> listOfBunchCrossingsRight, listOfBunchCrossingsLeft;

  for (int i = 0; i < 10; ++i) {
    listOfBunchCrossingsLeft.push_back( new BCPadEnergies(m_BCG) );
    listOfBunchCrossingsRight.push_back( new BCPadEnergies(m_BCG) );
  }


  std::set<int> randomNumbers;
  const unsigned int nBackgroundBX = m_backgroundBX->GetEntries();

  //Check that we have
  if( int(nBackgroundBX) < m_nBX*10 ) {
    streamlog_out(ERROR) << "There are not enough BeamCal " \
     " Background files to calculate a proper average!" << endl;
  }

  //we use a set so no duplication occurs
  const int numberForAverage = 10;
  while( int(randomNumbers.size()) < m_nBX*numberForAverage ){//do it ten times as often
    randomNumbers.insert( int(m_random3->Uniform(0, nBackgroundBX)) );
  }

  int counter = 0;

  for (std::set<int>::iterator it = randomNumbers.begin(); it != randomNumbers.end();++it) {
    streamlog_out(DEBUG1) << std::setw(5) << *it << std::flush;
    m_backgroundBX->GetEntry(*it);
    m_BeamCalAverageLeft->addEnergies(*m_BeamCalDepositsLeft);
    m_BeamCalAverageRight->addEnergies(*m_BeamCalDepositsRight);
    listOfBunchCrossingsLeft.at(counter/m_nBX)-> addEnergies(*m_BeamCalDepositsLeft);
    listOfBunchCrossingsRight.at(counter/m_nBX)-> addEnergies(*m_BeamCalDepositsRight);
    ++counter;
  }

  //Now divide by ten, to get the average distributions...
  m_BeamCalAverageLeft ->scaleEnergies(1./double(numberForAverage));
  m_BeamCalAverageRight->scaleEnergies(1./double(numberForAverage));
  Double_t totalEnergyMean = m_BeamCalAverageLeft->getTotalEnergy();
  Double_t varEn(0.0);
  for (int l = 0; l < numberForAverage;++l) {
    varEn += (listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean) * (listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean);
  }//histograms
  varEn /= double(numberForAverage);
  varEn = sqrt(varEn);
  streamlog_out(MESSAGE4) << "Total Energy " << totalEnergyMean << " +- "
			  << varEn << " GeV/" << m_nBX <<"BX" << std::endl;

  // m_h3BeamCalAverageLeft = bc.getBeamCalHistogram("AverageLeft");
  // m_h3BeamCalAverageRight = bc.getBeamCalHistogram("AverageRight");

  //And now calculate the error for every bin....
  m_BeamCalErrorsLeft  = getBeamCalErrors(m_BeamCalAverageLeft,  listOfBunchCrossingsLeft, numberForAverage);
  m_BeamCalErrorsRight = getBeamCalErrors(m_BeamCalAverageRight, listOfBunchCrossingsRight, numberForAverage);

  //Add one sigma to the averages -- > just do it once here
  m_BeamCalAverageLeft ->addEnergies( m_BeamCalErrorsLeft );
  m_BeamCalAverageRight->addEnergies( m_BeamCalErrorsRight);


  streamlog_out(DEBUG1) << std::endl;

  for (int i = 0; i < numberForAverage;++i) {
    delete listOfBunchCrossingsLeft[i];
    delete listOfBunchCrossingsRight[i];
  }
  
  return 0;
}

/// Calculate the one sigma errors for the given selected background bunch crossings
BCPadEnergies* BeamCalBackground::getBeamCalErrors(const BCPadEnergies *averages, const std::vector<BCPadEnergies*> singles, int numberForAverage ) {

  BCPadEnergies * BCPErrors = new BCPadEnergies(m_BCG);
  for (int i = 0; i < m_BCG->getPadsPerBeamCal()  ;++i) {
    Double_t mean(averages->getEnergy(i));
    Double_t variance(0);
    Double_t nHistos(numberForAverage);
    for (int l = 0; l < nHistos;++l) {
      Double_t energy( singles[l]->getEnergy(i) );
      variance += ( energy - mean ) * ( energy - mean );
    }//histograms
    variance /= nHistos;
    variance = sqrt(variance);
    BCPErrors->setEnergy(i, variance);

  }//for all pads

  return BCPErrors;

}//getBeamCalErrors

int BeamCalBackground::initParametrised(vector<string> &bg_files)
{
  // too many checks
  if (1 != bg_files.size() ){
    streamlog_out(ERROR) << "Single background file should be specified"\
      "for parametrised BG simulation method." << std::endl;
  }

  TTree *bg_par_tree;
  TString bgfname(bg_files.at(0).c_str());
  TFile *bgfile = TFile::Open(bgfname);
  if ( !bgfile ) {
    streamlog_out(ERROR) << "Background file " << bg_files.at(0) << " not found" << std::endl;;
  }

  bgfile->GetObject("bc_bg_fitpars", bg_par_tree);
  if ( !bg_par_tree ) {
    streamlog_out(ERROR) << "Root tree with background parameters" \
         "bc_bg_fitpars is not found in the file "+bgfname << std::endl;;

    bgfile->Close();
    delete bgfile;
  }


  const int nBCpads = m_BCG->getPadsPerBeamCal();
  m_padParLeft  = new vector<PadEdepRndPar_t>(nBCpads); 
  m_padParRight = new vector<PadEdepRndPar_t>(nBCpads); 

  m_BeamCalAverageLeft  =  new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight =  new BCPadEnergies(m_BCG);
  m_BeamCalErrorsLeft   =  new BCPadEnergies(m_BCG);
  m_BeamCalErrorsRight  =  new BCPadEnergies(m_BCG);

  readBackgroundPars(bg_par_tree, BCPadEnergies::kLeft);
  readBackgroundPars(bg_par_tree, BCPadEnergies::kRight);

  bgfile->Close();

  return 0;
}

int BeamCalBackground::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{

  switch (m_bgMethod) {
    case kPregenerated: getEventPregeneratedBG(peLeft, peRight); break;
    case kParametrised: {
      getEventParametrisedBG(peLeft, BCPadEnergies::kLeft); 
      getEventParametrisedBG(peRight, BCPadEnergies::kRight); 
      break;
    }
  }

  return 0;
}

int BeamCalBackground::getEventPregeneratedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
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

  return 0;
}

int BeamCalBackground::getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalAverageLeft);
  peRight.setEnergies(*m_BeamCalAverageRight);

  return 0;
}

int BeamCalBackground::getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalErrorsLeft);
  peRight.setEnergies(*m_BeamCalErrorsRight);
  
  return 0;
}


int BeamCalBackground::getEventParametrisedBG(BCPadEnergies &pe, 
      const BCPadEnergies::BeamCalSide_t bc_side)
{
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<double> vedep(nBCpads, 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = (BCPadEnergies::kLeft == bc_side ? m_padParLeft->at(ip) : m_padParRight->at(ip));
    //std::cout << ip << "\t" << pep.zero_rate << "\t" << pep.mean << "\t" << pep.stdev<< "\t" << pep.par0 << "\t" <<  pep.par1 << "\t" <<  pep.par2 << "\t" << pep.chi2 << endl;

    // put the average to 0 at once by generating fluctiations with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    vedep.at(ip) = m_random3->Gaus(0., pep.stdev*sqrt(m_nBX));
    /*
    for (int ibx=0; ibx<m_nBX; ibx++){
      // check if zero
      if (m_random3->Uniform(1.) < pep.zero_rate ) continue ; // == {vedep.at(ip)+=0.};
      // do gauss from mean, stdev
      vedep.at(ip) += m_random3->Gaus(pep.mean, pep.stdev);
    }
    */
  }

  pe.setEnergies(vedep);

  streamlog_out(DEBUG) << "BeamCalBackground: total energy generated with parametrised method" \
    " for " << (BCPadEnergies::kLeft == bc_side ? string("Left") : string("Right") ) <<
    " BeamCal = " << pe.getTotalEnergy() << endl;

  return 0;
}


int BeamCalBackground::readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side)
{
  vector<double> *zero_rate(NULL) , *mean(NULL) , *stdev(NULL);
  vector<double> *sum(NULL)       , *minm(NULL) , *maxm(NULL);
  vector<double> *chi2(NULL)      , *par0(NULL) , *par1(NULL)   , *par2(NULL);
  string side_name = ( BCPadEnergies::kLeft == bc_side ? "left_" : "right_" );

  // map the branch names to containers
  typedef map<string, vector<double> *> MapStrVec_t;
  MapStrVec_t br_cont_map;

  br_cont_map[side_name+"zero_rate"] = zero_rate;
  br_cont_map[side_name+"mean"]      = mean;
  br_cont_map[side_name+"stdev"]     = stdev;
  br_cont_map[side_name+"sum"]       = sum;
  br_cont_map[side_name+"minm"]      = minm;
  br_cont_map[side_name+"maxm"]      = maxm;
  br_cont_map[side_name+"chi2"]      = chi2;
  br_cont_map[side_name+"par0"]      = par0;
  br_cont_map[side_name+"par1"]      = par1;
  br_cont_map[side_name+"par2"]      = par2;

  // check the branch presence in the tree
  MapStrVec_t::iterator im = br_cont_map.begin();
  for (; im != br_cont_map.end(); im++){
    TBranch* br = dynamic_cast<TBranch*> (bg_par_tree->GetListOfBranches()->FindObject((im->first).c_str()));
    if (! br) {
      streamlog_out(ERROR7) << "BeamCalBackground: Missing " << im->first << 
      " branch in the tree with background parameters." << std::endl;
    }
  }

  // set the branch addresses
  for (im = br_cont_map.begin(); im != br_cont_map.end(); im++){
    bg_par_tree->SetBranchAddress((im->first).c_str(), &(im->second));
  }

  bg_par_tree->GetEntry(0);

  // check that number of pads in read vectors is equal one in current geometry
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  if ( nBCpads != int(br_cont_map[side_name+"zero_rate"]->size()) ){
    streamlog_out(ERROR7) << "Number of BeaCal pads in the background "\
      "file is not equal with one in current geometry." << std::endl;
  }

  vector<double> *pad_sigma = new vector<double>; 
  vector<double> *pad_mean = new vector<double>; 

  PadEdepRndPar_t pad_par;
  for(int ip=0; ip< nBCpads; ip++){
    pad_par.zero_rate = br_cont_map[side_name+"zero_rate"]->at(ip);
    pad_par.mean      = br_cont_map[side_name+"mean"]->at(ip);
    pad_par.stdev     = br_cont_map[side_name+"stdev"]->at(ip);
    pad_par.sum       = br_cont_map[side_name+"sum"]->at(ip);
    pad_par.minm      = br_cont_map[side_name+"minm"]->at(ip);
    pad_par.maxm      = br_cont_map[side_name+"maxm"]->at(ip);
    pad_par.chi2      = br_cont_map[side_name+"chi2"]->at(ip);
    pad_par.par0      = br_cont_map[side_name+"par0"]->at(ip);
    pad_par.par1      = br_cont_map[side_name+"par1"]->at(ip);
    pad_par.par2      = br_cont_map[side_name+"par2"]->at(ip);

    if (BCPadEnergies::kLeft == bc_side )  m_padParLeft->at(ip) = pad_par;
    else m_padParRight->at(ip) = pad_par;

    pad_sigma->push_back(pad_par.stdev*sqrt(m_nBX));
    pad_mean->push_back(pad_par.mean);
  }


  if (BCPadEnergies::kLeft == bc_side )  {
    m_BeamCalErrorsLeft->setEnergies(*pad_sigma);
    //m_BeamCalAverageLeft->setEnergies(*pad_mean);
    m_BeamCalAverageLeft->resetEnergies();
  } else {
    m_BeamCalErrorsRight->setEnergies(*pad_sigma);
    //m_BeamCalAverageRight->setEnergies(*pad_mean);
    m_BeamCalAverageRight->resetEnergies();
  }

  delete pad_sigma;
  delete pad_mean;

  // drop the branch addresses
  bg_par_tree->ResetBranchAddresses();
}

void BeamCalBackground::setRandom3Seed(int seed)
{ 
  m_random3->SetSeed(seed); 
}
