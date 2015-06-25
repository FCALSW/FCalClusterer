#include "BackgroundFitter.hh"

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <numeric>

#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;

extern int slice_pad(vector<vector<double> > &vvbg, int ip, vector<double>& vout);
extern int estimate_pars(vector<double> &vpad, double& zr, 
            double& mean, double &stdev, double &sum, double &minm, double &maxm);

BackgroundFitter::BackgroundFitter(int npads) : _mean(NULL), _stdev(NULL) 
{
  _mean = new vector<double>(npads,0.);
  _stdev = new vector<double>(npads,0.);
}

BackgroundFitter::BackgroundFitter(const BackgroundFitter&bf) : 
                  _mean(bf._mean),
		  _stdev(bf._stdev)
{}

BackgroundFitter &
BackgroundFitter::operator=(const BackgroundFitter&bf)
{
  _mean = bf._mean;
  _stdev = bf._stdev;

  return *this;
}

BackgroundFitter::~BackgroundFitter(){
  delete _mean;
  delete _stdev;
}

double 
BackgroundFitter::Fit(int ip, vector<double> &vpad){

  double zr(0.), mean(0.), stdev(0.), sum(0.), minm(0.), maxm(0.);
  estimate_pars(vpad, zr, mean, stdev, sum, minm, maxm);

  if ( 1. == zr ) return 0.;

  // trying to calculate number of bins in the histogram to make
  // a reasonable fit.
  const int n_entries_per_bin = 10;
  int nbins = vpad.size()*(1.-zr)/n_entries_per_bin; 

  // set fit, parameters and limits
  TF1 *bgfit = new TF1("bgfit", "gaus(0)/x", 0., maxm*1.04);

  //cout << vpad.size() << "\t" << zr << "\t" <<  mean << "\t" <<  stdev << "\t" <<  maxm<<  endl;
  bgfit->SetParameter(0, 0.5*double(vpad.size()*(1.-zr)));
  bgfit->SetParLimits(0, 0., double(vpad.size()*(1.-zr)));
  bgfit->SetParameter(1, min(mean*(1.+pow(zr,2)), mean+3*stdev));
  bgfit->SetParLimits(1, mean-3*stdev, mean+3*stdev);
  bgfit->SetParameter(2, 2.*stdev);
  bgfit->SetParLimits(2, 0.1*stdev, 5*stdev);

  // fill a histogram with non-zero values
  TH1F *h = new TH1F("h", "h", nbins, 0., maxm*1.04);
  for (size_t ie=0; ie<vpad.size(); ie++){
    if (vpad.at(ie) != 0. ) {
      h->Fill(vpad.at(ie));
    }
  }

  // get fitted parameters
  TFitResultPtr r = h->Fit("bgfit","QS");

  _mean->at(ip)      = mean;
  _stdev->at(ip)     = stdev;

  delete h;
  delete bgfit;

  return r->Chi2();
}

void
BackgroundFitter::WriteFitPars(TTree* tree, int dir){
  string lr("left");
  if(dir >0. ) lr.assign("right");
    
  tree->Branch((lr+"_mean").c_str(),      "std::vector<double>", &_mean);
  tree->Branch((lr+"_stdev").c_str(),     "std::vector<double>", &_stdev);

  return 0;
}

