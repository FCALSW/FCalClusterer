#pragma once

#include <vector>

class TTree;

using std::vector;

class BackgroundFitter {
 public:
  BackgroundFitter(int npads);
  ~BackgroundFitter();

 private:
  BackgroundFitter(const BackgroundFitter&);
  BackgroundFitter &operator=(const BackgroundFitter&);

 public:
  double Fit(int ip, vector<double> &vpad);
  void WriteFitPars(TTree* tree, int dir);

 private:
  vector<double> *_mean;
  vector<double> *_stdev;
};
