/**
* @file BeamCalBkgParam.cpp
* @brief Methods of derived class for parametrised background
* @author Andrey Sapronov
* @version v0.01
* @date 2015-06-12
*/


#include "BeamCalBkg.hh"
#include "BeamCalBkgParam.hh"
#include "BeamCalGeoCached.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"


// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

// ROOT
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

BeamCalBkgParam::BeamCalBkgParam(const string& bg_method_name, 
                     const BeamCalGeo *BCG) : BeamCalBkg(bg_method_name, BCG),
					   m_padParLeft(NULL),
					   m_padParRight(NULL)
{}

BeamCalBkgParam::~BeamCalBkgParam()
{}

void BeamCalBkgParam::init(vector<string> &bg_files, const int n_bx)
{
  this->BeamCalBkg::init(n_bx);

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

}


void BeamCalBkgParam::getEventBG(BCPadEnergies &peLeft, BCPadEnergies &peRight)
{
  const int nBCpads = m_BCG->getPadsPerBeamCal();
  vector<double> vedep(nBCpads, 0.);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParLeft->at(ip);

    // generating fluctiations at once with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    vedep.at(ip) = m_random3->Gaus(0., pep.stdev*sqrt(m_nBX));
    //vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
  }

  peLeft.setEnergies(vedep);

  for (int ip=0; ip< nBCpads; ip++){
    // Parameters of energy deposition in a pad:
    PadEdepRndPar_t pep = m_padParRight->at(ip);

    // generating fluctiations at once with stdev*sqrt(nBX)
    // otherwise the time to generate each event grows too much
    vedep.at(ip) = m_random3->Gaus(0., pep.stdev*sqrt(m_nBX));
    //vedep.at(ip) = m_random3->Gaus(pep.mean*m_nBX, pep.stdev*sqrt(m_nBX));
  }

  peRight.setEnergies(vedep);

  streamlog_out(DEBUG) << "BeamCalBkgParam: total energy generated with parametrised method for "
		       << "Left and Right BeamCal = " << peLeft.getTotalEnergy() << "\t" 
		       << peRight.getTotalEnergy() << std::endl;
}

void BeamCalBkgParam::readBackgroundPars(TTree *bg_par_tree, const BCPadEnergies::BeamCalSide_t bc_side)
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
      streamlog_out(ERROR7) << "BeamCalBkgParam: Missing " << im->first << 
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
