/**
* @file BeamCalBkgParam.cpp
* @brief Methods of derived class for parametrised background
* @author Andrey Sapronov
* @version v0.01
* @date 2015-06-12
*/

#include "BeamCalBkgParam.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"
#include "BeamCalBkg.hh"
#include "BeamCalGeo.hh"

// ----- include for verbosity dependent logging ---------
#include <streamlog/baselevels.h>
#include <streamlog/loglevels.h>
#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>


// ROOT
#include <RVersion.h>
#include <TBranch.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

using std::vector;
using std::string;
using std::map;

// IWYU pragma: no_include <ext/alloc_traits.h>

BeamCalBkgParam::BeamCalBkgParam(const string& bg_method_name, const BeamCalGeo* BCG)
    : BeamCalBkg(bg_method_name, BCG),
      m_padParLeft(nullptr),
      m_padParRight(nullptr),
      m_unuransLeft(vector<TF1*>()),
      m_unuransRight(vector<TF1*>()) {}

BeamCalBkgParam::~BeamCalBkgParam()
{
  vector<TF1*>::iterator iun = m_unuransLeft.begin();
  for (;iun!=m_unuransLeft.end(); iun++){
    delete *iun;
  }

  iun = m_unuransRight.begin();
  for (;iun!=m_unuransRight.end(); iun++){
    delete *iun;
  }

  delete m_padParLeft;
  delete m_padParRight;

}

void BeamCalBkgParam::init(vector<string> &bg_files, const int n_bx)
{
  this->BeamCalBkg::init(n_bx);

  // too many checks
  if (1 != bg_files.size() ){
    streamlog_out(ERROR) << "Single background file should be specified"\
      "for parametrised BG simulation method." << std::endl;
    throw std::runtime_error("Need exactly 1 (one) background fill for this BeamCalReco Background");
  }

  TTree *bg_par_tree;
  TString bgfname(bg_files.at(0).c_str());
  TFile *bgfile = TFile::Open(bgfname);
  if ( !bgfile ) {
    streamlog_out(ERROR) << "Background file " << bg_files.at(0) << " not found" << std::endl;;
    throw std::runtime_error("Could not find background file for BeamCalReco");
  }

  bgfile->GetObject("bc_bg_fitpars", bg_par_tree);
  if ( !bg_par_tree ) {
    streamlog_out(ERROR) << "Root tree with background parameters" \
         "bc_bg_fitpars is not found in the file "+bgfname << std::endl;;

    bgfile->Close();
    delete bgfile;
    throw std::runtime_error("Could not find background tree in file for BeamCalReco");
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

  // set background distributions
  setBkgDistr(BCPadEnergies::kLeft);
  setBkgDistr(BCPadEnergies::kRight);

  // calculate st.dev. of tower energies
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kLeft);
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kRight);

  bgfile->Close();

}


void BeamCalBkgParam::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  BCUtil::IgnoreRootError ire( not streamlog_level(DEBUG0) );

  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<double> vedep(nBCpads, 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParLeft->at(ip);

    		//std::cout << ip<< "\t" <<m_unuransLeft.at(ip) << std::endl;
    if (m_unuransLeft.at(ip)){
      for (int ibx=0; ibx<m_nBX; ibx++){
        // check if zero
        if (m_random3->Uniform(1.) < pep.zero_rate ) continue ; // == {vedep.at(ip)+=0.};
	// add value
        vedep.at(ip) += m_unuransLeft.at(ip)->GetRandom();
      }
    } else  {
      // if unuran is null, than it's just gaus
      // generating fluctuations at once with stdev*sqrt(nBX)
      // otherwise the time to generate each event grows too much
      vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
    }
  }

  peLeft.setEnergies(vedep);
  std::fill(vedep.begin(), vedep.end(), 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParRight->at(ip);

    if (m_unuransRight.at(ip)){
      for (int ibx=0; ibx<m_nBX; ibx++){
        // check if zero
        if (m_random3->Uniform(1.) < pep.zero_rate ) continue ; // == {vedep.at(ip)+=0.};
	// add value
        vedep.at(ip) += m_unuransRight.at(ip)->GetRandom();
      }
    } else  {
      // if unuran is null, than it's just gaus
      // generating fluctiations at once with stdev*sqrt(nBX)
      // otherwise the time to generate each event grows too much
      vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
    }
  }

  peRight.setEnergies(vedep);

  streamlog_out(DEBUG) << "BeamCalBkgParam: total energy generated with parametrised method for "
		       << "Left and Right BeamCal = " << peLeft.getTotalEnergy() << "\t" 
		       << peRight.getTotalEnergy() << std::endl;

}

void BeamCalBkgParam::readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side)
{
  string side_name = ( BCPadEnergies::kLeft == bc_side ? "left_" : "right_" );

  // map the branch names to containers
  typedef map<string, vector<double> *> MapStrVec_t;
  MapStrVec_t br_cont_map;
  //Set to nullptr so root deals with the memory itself
  br_cont_map[side_name + "zero_rate"] = nullptr;
  br_cont_map[side_name + "mean"]      = nullptr;
  br_cont_map[side_name + "stdev"]     = nullptr;
  br_cont_map[side_name + "sum"]       = nullptr;
  br_cont_map[side_name + "minm"]      = nullptr;
  br_cont_map[side_name + "maxm"]      = nullptr;
  br_cont_map[side_name + "chi2"]      = nullptr;
  br_cont_map[side_name + "par0"]      = nullptr;
  br_cont_map[side_name + "par1"]      = nullptr;
  br_cont_map[side_name + "par2"]      = nullptr;

  // check the branch presence in the tree
  bool errorGettingBranches = false;
  for (MapStrVec_t::iterator im = br_cont_map.begin(); im != br_cont_map.end(); im++){
    TBranch* br = dynamic_cast<TBranch*> (bg_par_tree->GetListOfBranches()->FindObject((im->first).c_str()));
    if (! br) {
      streamlog_out(ERROR7) << "BeamCalBkgParam: Missing " << im->first << 
      " branch in the tree with background parameters." << std::endl;
      errorGettingBranches = true;
    }
  }

  if( errorGettingBranches ) {
    throw std::runtime_error("Failed to get branches from tree in background file");
  }

  // set the branch addresses
  for (MapStrVec_t::iterator im = br_cont_map.begin(); im != br_cont_map.end(); im++){
    bg_par_tree->SetBranchAddress((im->first).c_str(), &(im->second));
  }

  bg_par_tree->GetEntry(0);

  // check that number of pads in read vectors is equal one in current geometry
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  if ( nBCpads != int(br_cont_map[side_name+"mean"]->size()) ){
    streamlog_out(ERROR7) << "Number of BeamCal pads in the background "\
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
  for (MapStrVec_t::iterator im = br_cont_map.begin(); im != br_cont_map.end(); im++){
    delete im->second;
  }

  // drop the branch addresses
  bg_par_tree->ResetBranchAddresses();
}

double gausOverX(double *x, double* p) {
  return p[0]*TMath::Gaus(x[0], p[1], p[2], false)/x[0];
}


int BeamCalBkgParam::setBkgDistr(const BCPadEnergies::BeamCalSide_t bc_side)
{
  //this in highest debug level only
  BCUtil::IgnoreRootError ire( not streamlog_level(DEBUG0) );

  streamlog_out( MESSAGE0 ) << "Creating Background Distributions: " << bc_side << std::endl;

  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<TF1*> &vunr = (BCPadEnergies::kLeft == bc_side ? m_unuransLeft : m_unuransRight);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = (BCPadEnergies::kLeft == bc_side ? 
                           m_padParLeft->at(ip) : m_padParRight->at(ip));
    vunr.push_back(nullptr);
    //TF1 *func_pad_edep = new TF1("fedep", "gaus(0)/x", pep.minm, pep.maxm);
    // if chi2 is good enough we want to initialise our unuran randomiser
    if (pep.chi2 <= 200. && pep.par1 >= 2*pep.par2) {
      TF1 *func_pad_edep = new TF1("fedep", gausOverX, pep.minm, pep.maxm, 3);
      func_pad_edep->SetParameters(pep.par0, pep.par1, pep.par2);
		/*
    		std::cout << ip << "\t" 
		          << pep.mean << "\t" << pep.stdev<< "\t" 
			  << pep.par0 << "\t" <<  pep.par1 << "\t" 
			  <<  pep.par2 << "\t" << pep.chi2 << std::endl;
			  */

      double funcparam[3] ={pep.par0, pep.par1, pep.par2};

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
      double integral = func_pad_edep->Integral(pep.minm, pep.maxm, funcparam, 0.0001);
#else
      func_pad_edep->SetParameters(funcparam);
      double integral = func_pad_edep->Integral(pep.minm, pep.maxm, 0.0001);
#endif
      if (integral > 0.001) {
	vunr.back() = func_pad_edep;
      } else {
	streamlog_out( DEBUG1 ) << "Failed to create gaus/x background distribution for this pad: " << ip << std::endl;
	delete func_pad_edep;
      }
    }
  }

  return 0;
}
