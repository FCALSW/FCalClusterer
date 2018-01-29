#include "BackgroundFitter.hh"

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TTree.h>

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

extern int slice_pad(vector<vector<double> > &vvbg, int ip, vector<double>& vout);
extern int estimate_pars(vector<double> &vpad, double& zr, 
            double& mean, double &stdev, double &sum, double &minm, double &maxm);

BackgroundFitter::BackgroundFitter(int npads) 
                 : 
		   _zero_rate(new vector<double>(npads,0.)), 
		   _chi2(new vector<double>(npads,0.)), 
		   _par0(new vector<double>(npads,0.)), 
		   _par1(new vector<double>(npads,0.)), 
		   _par2(new vector<double>(npads,0.)), 
		   _mean(new vector<double>(npads,0.)), 
		   _sum(new vector<double>(npads,0.)), 
		   _minm(new vector<double>(npads,0.)), 
		   _maxm(new vector<double>(npads,0.)), 
		   _stdev(new vector<double>(npads,0.)) 
{}

BackgroundFitter::BackgroundFitter(const BackgroundFitter&bf)
		 : _zero_rate(NULL), 
		   _chi2(NULL), 
		   _par0(NULL), 
		   _par1(NULL), 
		   _par2(NULL), 
		   _mean(NULL), 
		   _sum(NULL), 
		   _minm(NULL), 
		   _maxm(NULL), 
		   _stdev(NULL) 
{
  _zero_rate = new vector<double>(*(bf._zero_rate));
  _chi2 = new vector<double>(*(bf._chi2));
  _par0 = new vector<double>(*(bf._par0));
  _par1 = new vector<double>(*(bf._par1));
  _par2 = new vector<double>(*(bf._par2));
  _mean = new vector<double>(*(bf._mean));
  _sum = new vector<double>(*(bf._sum));
  _minm = new vector<double>(*(bf._minm));
  _maxm = new vector<double>(*(bf._maxm));
  _stdev = new vector<double>(*(bf._stdev));
}

BackgroundFitter &
BackgroundFitter::operator=(const BackgroundFitter&bf)
{
  _zero_rate = new vector<double>(*(bf._zero_rate));
  _chi2 = new vector<double>(*(bf._chi2));
  _par0 = new vector<double>(*(bf._par0));
  _par1 = new vector<double>(*(bf._par1));
  _par2 = new vector<double>(*(bf._par2));
  _mean = new vector<double>(*(bf._mean));
  _sum = new vector<double>(*(bf._sum));
  _maxm = new vector<double>(*(bf._maxm));
  _minm = new vector<double>(*(bf._minm));
  _stdev = new vector<double>(*(bf._stdev));

  return *this;
}

BackgroundFitter::~BackgroundFitter(){
  delete _zero_rate;
  delete _chi2;
  delete _par0;
  delete _par1;
  delete _par2;
  delete _mean;
  delete _sum;
  delete _minm;
  delete _maxm;
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

  _zero_rate->at(ip) = zr;
  _chi2->at(ip)      = r->Chi2();
  _par0->at(ip)      = r->Parameter(0);
  _par1->at(ip)      = r->Parameter(1);
  _par2->at(ip)      = r->Parameter(2);
  _mean->at(ip)      = mean;
  _sum->at(ip)       = sum;
  _minm->at(ip)      = minm;
  _maxm->at(ip)      = maxm;
  _stdev->at(ip)     = stdev;

//  std::cout << _chi2->at(ip) << std::endl;

  delete h;
  delete bgfit;

  return r->Chi2();
}

void
BackgroundFitter::WriteFitPars(TTree* tree, int dir){
  string lr("left");
  if(dir >0. ) lr.assign("right");
    
  tree->Branch((lr+"_zero_rate").c_str(), "std::vector<double>", &_zero_rate);
  tree->Branch((lr+"_chi2").c_str(),      "std::vector<double>", &_chi2);
  tree->Branch((lr+"_par0").c_str(),      "std::vector<double>", &_par0);
  tree->Branch((lr+"_par1").c_str(),      "std::vector<double>", &_par1);
  tree->Branch((lr+"_par2").c_str(),      "std::vector<double>", &_par2);
  tree->Branch((lr+"_mean").c_str(),      "std::vector<double>", &_mean);
  tree->Branch((lr+"_sum").c_str(),       "std::vector<double>", &_sum);
  tree->Branch((lr+"_minm").c_str(),      "std::vector<double>", &_minm);
  tree->Branch((lr+"_maxm").c_str(),      "std::vector<double>", &_maxm);
  tree->Branch((lr+"_stdev").c_str(),     "std::vector<double>", &_stdev);
}

