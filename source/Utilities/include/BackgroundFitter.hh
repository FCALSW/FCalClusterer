#pragma once

#include <vector>

class TTree;

using std::vector;

class BackgroundFitter {
 public:
  BackgroundFitter(int npads);
  ~BackgroundFitter();

 public:
  double Fit(int ip, vector<double> &vpad);
  int WriteFitPars(TTree* tree, int dir);

 private:
  vector<double> *_zero_rate;
  vector<double> *_chi2;
  vector<double> *_par0;
  vector<double> *_par1;
  vector<double> *_par2;
  vector<double> *_mean;
  vector<double> *_sum;
  vector<double> *_minm;
  vector<double> *_maxm;
  vector<double> *_stdev;
};
