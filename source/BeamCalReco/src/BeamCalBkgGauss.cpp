/**
* @file BeamCalBkgGauss.cpp
* @brief Methods of derived class for gaussian background
* @author Andrey Sapronov
* @version v0.01
* @date 2015-06-12
*/

#include "BeamCalBkgGauss.hh"
#include "BCPadEnergies.hh"
#include "BeamCalBkg.hh"
#include "BeamCalGeo.hh"

// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

// ROOT
#include <TBranch.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

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

BeamCalBkgGauss::BeamCalBkgGauss(const string& bg_method_name, const BeamCalGeo* BCG)
    : BeamCalBkg(bg_method_name, BCG), m_padParLeft(nullptr), m_padParRight(nullptr) {}

BeamCalBkgGauss::~BeamCalBkgGauss() {
  delete m_padParLeft;
  delete m_padParRight;

}

void BeamCalBkgGauss::init(vector<string> &bg_files, const int n_bx)
{
  this->BeamCalBkg::init(n_bx);

  // too many checks
  if (1 != bg_files.size() ){
    streamlog_out(ERROR) << "Single background file should be specified"\
      "for gaussian BG simulation method." << std::endl;
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

  // calculate st.dev. of tower energies
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kLeft);
  this->BeamCalBkg::setTowerErrors(BCPadEnergies::kRight);

  bgfile->Close();
  delete bgfile;

}


void BeamCalBkgGauss::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<double> vedep(nBCpads, 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParLeft->at(ip);

    // generating fluctiations at once with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    //vedep.at(ip) = m_random3->Gaus(0., pep.stdev*sqrt(m_nBX));
    vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
  }

  peLeft.setEnergies(vedep);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParRight->at(ip);

    // generating fluctiations at once with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    //vedep.at(ip) = m_random3->Gaus(0., pep.stdev*sqrt(m_nBX));
    vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
  }

  peRight.setEnergies(vedep);

  streamlog_out(DEBUG) << "BeamCalBkgGauss: total energy generated with gaussian method for "
		       << "Left and Right BeamCal = " << peLeft.getTotalEnergy() << "\t" 
		       << peRight.getTotalEnergy() << std::endl;

}

void BeamCalBkgGauss::readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side)
{

  string side_name = ( BCPadEnergies::kLeft == bc_side ? "left_" : "right_" );

  // map the branch names to containers
  typedef map<string, vector<double> *> MapStrVec_t;
  MapStrVec_t br_cont_map;

  br_cont_map[side_name + "mean"]  = nullptr;
  br_cont_map[side_name + "stdev"] = nullptr;

  // check the branch presence in the tree
  MapStrVec_t::iterator im = br_cont_map.begin();
  for (; im != br_cont_map.end(); im++){
    TBranch* br = dynamic_cast<TBranch*> (bg_par_tree->GetListOfBranches()->FindObject((im->first).c_str()));
    if (! br) {
      streamlog_out(ERROR7) << "BeamCalBkgGauss: Missing " << im->first << 
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
    streamlog_out(ERROR7) << "Number of BeamCal pads in the background "\
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

  delete br_cont_map[side_name+"mean"];
  delete br_cont_map[side_name+"stdev"];

}
