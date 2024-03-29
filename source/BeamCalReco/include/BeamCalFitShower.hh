/**
* @file BeamCalFitShower.hh
* @brief A simple shower fitting class for BeamCal simulation
* @author Andrey Sapronov <sapronov@cern.ch>
* @version 0.0.1
* @date 2015-03-30
*/

#pragma once

#include "BCPadEnergies.hh"

#include <vector>
#include <map>

class BeamCalGeo;
class BeamCalBkg;
class BeamCalPadGeometry;

/**
* @brief Segment parameters for profile of the calorimeter energy deposition
*/
struct EdepProfile_t {
  int id=0;
  double towerChi2=0;
  double totalEdep=0;
  double bkgEdep=0;
  double bkgSigma=0;
  BeamCalPadGeometry *padGeom=nullptr;
  std::map<int, double> padIDs{};
} ;


class BeamCalFitShower {
 public:
  /**
  * @brief Initialises the shower fitter
  *
  * @param vep profile of energy deposit in the BeamCal as a vector of 
  * segment parameters.
  */
  BeamCalFitShower(std::vector<EdepProfile_t*> &vep, const BCPadEnergies::BeamCalSide_t bc_side);
  ~BeamCalFitShower(){};

  /**
   * @brief Copy constructor for shower fitter
   *
   * @param fs reference to BeamCalFitShower object to be copied
   */
  BeamCalFitShower(const BeamCalFitShower& fs);

  /**
   * @brief Assignment operator
   *
   * @param fs reference to BeamCalFitShower object to be copied
   */
  BeamCalFitShower& operator=(const BeamCalFitShower&fs);

  /**
  * @brief Fits the shower with 2d Gauss
  *
  * The method takes maximum energy pad from profile, 
  * fits shower with 2d Gauss around it, estimates probability 
  * that it is actually a shower, calculates its parameters and returns all that
  *
  * Independing on the probability, the fitted shower energy 
  * is subtracted from the profile.
  *
  * @param theta 
  * @param phi
  * @param en_shower 
  *
  * @return probability value that we have an actuall shower
  */
  double fitShower(double &theta, double &phi, double &en_shwr, double &chi2, std::map<int, double>& padIDsInCluster);

  double operator()(const double *par);
  //double showerChi2(double *par);
  
  void setGeometry(const BeamCalGeo *BCG) { m_BCG = BCG; }
  void setBackground(const BeamCalBkg *BCbg) { m_BCbackground = BCbg; }
  void setStartLayer(const int sl) { m_startLayer = sl; }
  void setCountingLayers(const int cl) { m_countingLayer = cl; }
  void setEshwrLimit(double elimit) { m_enTowerLimit = elimit; }
  void setTowerChi2Limit(double tchi2lim) { m_towerChi2Limit = tchi2lim;}

 private:
  void estimateShowerPars(double &rc, double &phic, double &A0, double &sig0);
  int selectSpotPads(std::vector<int> &pad_ids);
  int calcCovar();
  void deleteSpotPads();

 private:
  std::vector<EdepProfile_t*> m_vep;
  /**
  * @brief List of profile segments falling into the shower spot
  */
  std::vector<EdepProfile_t*> m_spotPads;

  /**
  * @brief Inverse covariance matrix for spot pads
  */
  std::vector<double> m_covInv;

  /**
  * @brief integral of gaus-distributed energy in spot pads
  */
  std::vector<double> m_spotEint;

  const BeamCalGeo* m_BCG;
  const BeamCalBkg *m_BCbackground;
  const BCPadEnergies::BeamCalSide_t m_BCside;

  const double m_rhom;
  double m_enTowerLimit;
  double m_towerChi2Limit;
  int m_startLayer;
  int m_countingLayer;

  bool m_flagUncorr;
};


