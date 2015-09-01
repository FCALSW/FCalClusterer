/**
* @file BCBackgroundPar.cpp
* @brief The program to parametrise beamstrahlung background
* @author Andrey Sapronov <andrey.sapronov@cern.ch>
* @version 1.0.0
* @date 2015-03-04
*
* The code taks a list of BeamCal background energy depositons
* and calculates mean, stdev and other parameters of the 
* bx to bx distributions. It includes fitting the distribution
* with gaus()/x shape, with storing fit results for each pad.
* The fit parameters can potentially be used for more precise
* background representation at the reconstruction stage, 
* while in the simplest case it is Gauss-shaped with mean and stdev.
*
* Usage example:
* > find /path/to/BX/files/ -name *.root | ./BCBackgroundPar
*
*/


#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <numeric>

#include "TTree.h"
#include "TFile.h"

#include "BackgroundFitter.hh"

using namespace std;

typedef vector<vector<double> > Tvvd;

// reads pad energy vectors from bgfname and pushes them to 
// vbgleft, vbgright
int read_bg_vec(string &bgfname, Tvvd &vbgleft, Tvvd &vbgright);

// slicing functions to through-select corresponding array:
// puts ip'th pad from all BXs to vout
int slice_pad(Tvvd &vvbg, int ip, vector<double>& vout);
// puts il'th layer to vvout
//int slice_layer(Tvvd &vvbg, int il, Tvvd& vvout);
// puts tower under ip'th pad to vvout
//int slice_tower(Tvvd &vvbg_left, Tvvd &vvbg_right, int ip, Tvvd& vvout);

// function to estimate parameters of background distribution
int estimate_pars(vector<double> &vpad, double& zr, 
            double& mean, double &stdev, double &sum, double &minm, double &maxm);

int main(int argc, char **argv) {
  if (argc < 2 ) {
    std::cerr << "No input files provided"  << std::endl;
    std::cerr << "BCPackgroundPar background.root [[background2.root] ...]"  << std::endl;
    return 1;
  }


  Tvvd vbgleft, vbgright;
  //first entry is program name, start at 1
  for (int i = 1; i < argc; ++i) {
    std::cout << "Reading background file " << argv[i] << std::endl;
    std::string bgfname = argv[i];
    read_bg_vec(bgfname, vbgleft, vbgright);
  }

  unsigned int nbcpads = vbgleft.at(0).size();

  cout << "Left: " << nbcpads << endl;
  BackgroundFitter *fmleft = new BackgroundFitter(nbcpads);

  for (unsigned int ip = 0; ip<nbcpads; ip++){
    vector<double> vpad;
    int npads = slice_pad(vbgleft, ip, vpad);
    if (0 == npads ) cout << "No pads in slice " << ip << endl;
    fmleft->Fit(ip, vpad); 
//    if (chi2 > 150 ) cout << ip<< ":\t" << chi2 << endl;
  }

  if (vbgright.at(0).size() != nbcpads){
    cout << "Warning, number of pads is not equal in left and right" << endl;
  }

  nbcpads = vbgright.at(0).size();
  cout << "Right: " << nbcpads << endl;
  BackgroundFitter *fmright = new BackgroundFitter(nbcpads);

  for (unsigned int ip = 0; ip<nbcpads; ip++){
    vector<double> vpad;
    int npads = slice_pad(vbgright, ip, vpad);
    if (0 == npads ) cout << "No pads in slice " << ip << endl;
    fmright->Fit(ip, vpad); 
//    cout << ip<< ":\t" << chi2 << endl;
  }

  TFile *fout = new TFile("BeamCal_bg.root", "RECREATE");
  TTree *tree_fitpars = new TTree("bc_bg_fitpars", "bc_bg_fitpars");
  fmright->WriteFitPars(tree_fitpars, 1);
  fmleft->WriteFitPars(tree_fitpars, -1);
  tree_fitpars->Fill();

  tree_fitpars->Print();
  fout->Write();
  fout->Close();

  delete fmleft; delete fmright;
//  delete cm;
}


int read_bg_vec(string &bgfname, Tvvd &vbgleft, Tvvd &vbgright){
  TTree* tree;
  vector<double> *depLeft=NULL;
  vector<double> *depRight=NULL;

  TFile* file = TFile::Open(bgfname.c_str());
  if ( not file ) {
    std::cerr << "File not found: "<< bgfname  << std::endl;
  }

  file->GetObject("bcTree", tree);
  if ( not tree ) {
    std::cerr << "Tree not found in file " << bgfname  << std::endl;
    file->Close();
    delete file;
    return 0;
  }
  
  tree->SetBranchAddress("vec_left" , &depLeft);
  tree->SetBranchAddress("vec_right", &depRight);

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    vbgleft.push_back(*depLeft);
    vbgright.push_back(*depRight);
  }
  file->Close();

  return 1;
}


int slice_pad(Tvvd &vvbg, int ip, vector<double>& vout){
  vout.clear();
  Tvvd::iterator ivv = vvbg.begin();
  for (;ivv != vvbg.end(); ivv++){
    vout.push_back(ivv->at(ip));
  }
  return vout.size();
}

int estimate_pars(vector<double> &vpad, double& zr, 
            double& mean, double &stdev, double &sum, double &minm, double &maxm){
  zr = 0.; mean = 0.; stdev = 0.; sum = 0.; minm = 0.; maxm = 0.;
  unsigned int nz = count(vpad.begin(), vpad.end(), 0.);
  if ( nz == vpad.size() ){
    zr = 1.;
    return 0;
  }

  zr = double(nz)/vpad.size();
  sum = accumulate(vpad.begin(), vpad.end(), 0.0);
  mean = sum/(vpad.size()-nz);

  double sq_sum = inner_product(vpad.begin(), vpad.end(), vpad.begin(), 0.0);
  stdev = sqrt(sq_sum / (vpad.size()-nz) - mean * mean);

  minm = (double)*min_element(vpad.begin(), vpad.end());
  maxm = (double)*max_element(vpad.begin(), vpad.end());

  return vpad.size()-nz;
}
