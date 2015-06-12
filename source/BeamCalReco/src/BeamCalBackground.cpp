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
  if(      string("Pregenerated") == bg_method_name ) {m_bgMethod = kPregenerated;}
  else if( string("Parametrised") == bg_method_name ) {m_bgMethod = kParametrised;}
  else if( string("Averaged")     == bg_method_name ) {m_bgMethod = kAveraged;}
  else {
    streamlog_out(ERROR) <<"== Error From BeamCalBackground == "
			 << "Unknown BeamCal method of background esitmation \""
			 << bg_method_name << "\"" << std::endl;
    throw std::runtime_error("Unknown BeamCal background method");
  }

  streamlog_out(MESSAGE) << "initialising BeamCal background with \""
			 << bg_method_name << "\" method" << std::endl;
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

  delete m_TowerErrorsLeft;
  delete m_TowerErrorsRight;

  vector<BCPadEnergies*>::iterator it_pe = m_listOfBunchCrossingsRight.begin();
  for(; it_pe != m_listOfBunchCrossingsRight.end(); it_pe++)
    delete *it_pe;

  it_pe = m_listOfBunchCrossingsLeft.begin();
  for(; it_pe != m_listOfBunchCrossingsLeft.end(); it_pe++)
    delete *it_pe;
}

int BeamCalBackground::init(vector<string> &bg_files, const int n_bx)
{
  m_random3 = new TRandom3();
  m_nBX = n_bx;
  m_numberForAverage = 10;

  m_TowerErrorsLeft  =  new vector<double>;
  m_TowerErrorsRight =  new vector<double>;

  switch (m_bgMethod) {
    case kPregenerated     : return initPregenerated(bg_files);
    case kParametrised     : return initParametrised(bg_files);
    case kAveraged         : return initAveraged(bg_files);
  }

  return 1;
}


int BeamCalBackground::initAveraged(vector<string> &bg_files) {


  //We don't have any average energy in some cases, and we don't really need it
  //because we just subtract it after adding it in some cases anyway....
  m_BeamCalAverageLeft  = new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight = new BCPadEnergies(m_BCG);

  std::vector<BCPadEnergies> backgroundBeamCals(2, m_BCG);

  BCUtil::ReadBecasFile(bg_files[0], backgroundBeamCals, "tBcDensAverage", "sEdepErr", true);

  m_BeamCalErrorsLeft  = new BCPadEnergies(backgroundBeamCals[0]);
  m_BeamCalErrorsRight = new BCPadEnergies(backgroundBeamCals[1]);


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
  Double_t totalEnergyMean = m_BeamCalAverageLeft->getTotalEnergy();
  Double_t varEn(0.0);
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
  this->setTowerErrors(m_listOfBunchCrossingsLeft, BCPadEnergies::kLeft);
  this->setTowerErrors(m_listOfBunchCrossingsRight, BCPadEnergies::kRight);

  streamlog_out(DEBUG1) << std::endl;

  /*
  for (int i = 0; i < m_numberForAverage;++i) {
    delete m_listOfBunchCrossingsLeft[i];
    delete m_listOfBunchCrossingsRight[i];
  }
  */
  
  return 0;
}

/// Calculate the one sigma errors for the given selected background bunch crossings
BCPadEnergies* BeamCalBackground::getBeamCalErrors(const BCPadEnergies *averages, const std::vector<BCPadEnergies*> singles) {

  BCPadEnergies * BCPErrors = new BCPadEnergies(m_BCG);
  for (int i = 0; i < m_BCG->getPadsPerBeamCal()  ;++i) {
    Double_t mean(averages->getEnergy(i));
    Double_t variance(0);
    Double_t nHistos(m_numberForAverage);
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

void BeamCalBackground::setTowerErrors(const std::vector<BCPadEnergies*> singles, 
       const BCPadEnergies::BeamCalSide_t bc_side)
{
  BCPadEnergies * pe_aver = (BCPadEnergies::kLeft == bc_side 
    ? m_BeamCalAverageLeft : m_BeamCalAverageRight );

  // variance of tower energies
  vector<double>* te_var = (BCPadEnergies::kLeft == bc_side 
    ? m_TowerErrorsLeft : m_TowerErrorsRight );
  te_var->clear();

  // loop over pads in one layer == towers in BC
  for (int ip = 0; ip < m_BCG->getPadsPerLayer(); ip++){
    double te_aver = pe_aver->getTowerEnergy(ip,m_startLayer);
    te_var->push_back(0.);
    // loop over averaged entries
    for (int ibx = 0; ibx < m_numberForAverage; ibx++){
      te_var->back() += pow ( singles.at(ibx)->getTowerEnergy(ip, m_startLayer) - te_aver, 2);
    }
  }

  vector<double>::iterator it_tv = te_var->begin();
  for (; it_tv != te_var->end(); it_tv++){
    *it_tv = sqrt(*it_tv/m_numberForAverage);
  }
} // setTowerErrors

int BeamCalBackground::getPadsCovariance(vector<int> &pad_list, vector<double> &covinv, 
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
      
      /*
      if (i == j ) {
      std::cout << bxl.at(ibx)->getTowerEnergy(pad_list[i], m_startLayer)
      -   te_aver[j]
      << std::endl;
      }
      */
      
        ma[i*nm+j] += (bxl.at(ibx)->getTowerEnergy(pad_list[i], m_startLayer) - te_aver[i])
	         *(bxl.at(ibx)->getTowerEnergy(pad_list[j],m_startLayer) - te_aver[j])
		 /(m_numberForAverage-1);
      }
      //std::cout << ma[i*nm+j] << "\t";
    }
    //std::cout  << std::endl;
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
  /*
  for (int i = 0; i<nm; i++)
    for (int j = 0; j<nm; j++)
      std::cout << pad_list[i]<< "\t" <<pad_list[j]<< "\t" << ma[i*nm+j]<< std::endl;
      */

  double *mtxarr = mtxi.GetMatrixArray();
  covinv.assign(mtxarr, mtxarr+nm*nm);

  delete mtx;
  delete ma;

  return covinv.size();
}

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
} // initParametrised

void BeamCalBackground::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{

  switch (m_bgMethod) {
  case kPregenerated: {
      getEventPregeneratedBG(peLeft, peRight);
      break;
  }
  case kParametrised: {
    getEventParametrisedBG(peLeft, BCPadEnergies::kLeft);
    getEventParametrisedBG(peRight, BCPadEnergies::kRight);
    break;
  }
  case kAveraged: {
    getEventAveragedBG(peLeft, peRight);
    break;
  }
  }
}

void BeamCalBackground::getEventAveragedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight) {

  const int nBCpads = m_BCG->getPadsPerBeamCal();

  for (int j = 0; j < m_nBX;++j) {  //Add one for each BunchCrossing we want???
    for (int i = 0; i < nBCpads ;++i) { //Add gaussian randomisation of background to each cell
      peRight.addEnergy(i, m_random3->Gaus(0.0, m_BeamCalErrorsRight->getEnergy(i)));
      peLeft.addEnergy(i, m_random3->Gaus(0.0, m_BeamCalErrorsLeft->getEnergy(i)));
    }//for all BXs
  }//for all pads
} // getEventBG

void BeamCalBackground::getEventPregeneratedBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
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
} // getEventPregeneratedBG

void BeamCalBackground::getAverageBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalAverageLeft);
  peRight.setEnergies(*m_BeamCalAverageRight);
}

void BeamCalBackground::getErrorsBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  peLeft.setEnergies(*m_BeamCalErrorsLeft);
  peRight.setEnergies(*m_BeamCalErrorsRight);
}


int BeamCalBackground::getTowerErrorsBG(int padIndex, 
      const BCPadEnergies::BeamCalSide_t bc_side, double &tower_sigma)
{
  tower_sigma = (BCPadEnergies::kLeft == bc_side ? m_TowerErrorsLeft->at(padIndex) 
    : m_TowerErrorsRight->at(padIndex));
  
  return tower_sigma;
}

void BeamCalBackground::getEventParametrisedBG(BCPadEnergies &pe,
      const BCPadEnergies::BeamCalSide_t bc_side)
{
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<double> vedep(nBCpads, 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = (BCPadEnergies::kLeft == bc_side ? m_padParLeft->at(ip) : m_padParRight->at(ip));
    //std::cout << ip << "\t" << pep.zero_rate << "\t" << pep.mean << "\t" << pep.stdev<< "\t" << pep.par0 << "\t" <<  pep.par1 << "\t" <<  pep.par2 << "\t" << pep.chi2 << endl;

    // generating fluctiations at once with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
  }

  pe.setEnergies(vedep);

  streamlog_out(DEBUG) << "BeamCalBackground: total energy generated with parametrised method for "
		       << (BCPadEnergies::kLeft == bc_side ? string("Left") : string("Right") )
		       << " BeamCal = " << pe.getTotalEnergy()
		       << std::endl;

  return;
}


void BeamCalBackground::readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side)
{
  vector<double> *mean(NULL) , *stdev(NULL);
  string side_name = ( BCPadEnergies::kLeft == bc_side ? "left_" : "right_" );

  // map the branch names to containers
  typedef map<string, vector<double> *> MapStrVec_t;
  MapStrVec_t br_cont_map;

  br_cont_map[side_name+"mean"]      = mean;
  br_cont_map[side_name+"stdev"]     = stdev;

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
  if ( nBCpads != int(br_cont_map[side_name+"mean"]->size()) ){
    streamlog_out(ERROR7) << "Number of BeaCal pads in the background "\
      "file is not equal with one in current geometry." << std::endl;
  }

  vector<double> *pad_sigma = new vector<double>; 
  vector<double> *pad_mean = new vector<double>; 

  PadEdepRndPar_t pad_par;
  for(int ip=0; ip< nBCpads; ip++){
    pad_par.mean      = br_cont_map[side_name+"mean"]->at(ip);
    pad_par.stdev     = br_cont_map[side_name+"stdev"]->at(ip);

    if (BCPadEnergies::kLeft == bc_side )  m_padParLeft->at(ip) = pad_par;
    else m_padParRight->at(ip) = pad_par;

    pad_sigma->push_back(pad_par.stdev*sqrt(m_nBX));
    pad_mean->push_back(pad_par.mean*m_nBX);
  }


  if (BCPadEnergies::kLeft == bc_side )  {
    m_BeamCalErrorsLeft->setEnergies(*pad_sigma);
    m_BeamCalAverageLeft->setEnergies(*pad_mean);
    //m_BeamCalAverageLeft->resetEnergies();
  } else {
    m_BeamCalErrorsRight->setEnergies(*pad_sigma);
    m_BeamCalAverageRight->setEnergies(*pad_mean);
    //m_BeamCalAverageRight->resetEnergies();
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
