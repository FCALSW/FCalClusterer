/**
* @file BeamCalFitShower.hh
* @brief A simple shower fitting class for BeamCal simulation
* @author Andrey Sapronov <sapronov@cern.ch>
* @version 0.0.1
* @date 2015-03-30
*/

#pragma once

#include <vector>
#include <utility>
#include "BCPadEnergies.hh"

using std::vector;

class BeamCalGeo;
class BeamCalBackground;
class PadGeometry;

/**
* @brief Segment parameters for profile of the calorimeter energy deposition
*/
typedef struct {
  int id; 
  int ndf;
  double chi2;
  double en;
  double bg;
  double er;
  PadGeometry *pg;
} EdepProfile_t;


class BeamCalFitShower {
 public:
  /**
  * @brief Initialises the shower fitter
  *
  * @param vep profile of energy deposit in the BeamCal as a vector of 
  * segment parameters.
  */
  BeamCalFitShower(vector<EdepProfile_t*> &vep, const BCPadEnergies::BeamCalSide_t bc_side);
  ~BeamCalFitShower(){};

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
  double fitShower(double &theta, double &phi, double &en_shwr, double &chi2);

  double operator()(const double *par);
  //double showerChi2(double *par);
  
  void setGeometry(const BeamCalGeo *BCG) { m_BCG = BCG; }
  void setBackground(const BeamCalBackground *BCbg) { m_BCbackground = BCbg; }
  void setStartLayer(const int sl) { m_startLayer = sl; }
  void setEshwrLimit(double elimit) { m_enTowerLimit = elimit; }

 private:
  void estimateShowerPars(double &rc, double &phic, double &A0, double &sig0);
  int selectSpotPads(vector<int> &pad_ids);
  int calcCovar();
  void deleteSpotPads();

 private:
  vector<EdepProfile_t*> m_vep;
  /**
  * @brief List of profile segments falling into the shower spot
  */
  vector<EdepProfile_t*> m_spotPads;

  /**
  * @brief Inverse covariance matrix for spot pads
  */
  vector<double> m_covInv;

  /**
  * @brief integral of gaus-distributed energy in spot pads
  */
  vector<double> m_spotEint;

  const BeamCalGeo* m_BCG;
  const BeamCalBackground *m_BCbackground;
  const BCPadEnergies::BeamCalSide_t m_BCside;

  const double m_rhom;
  double m_enTowerLimit;
  int m_startLayer;

  bool m_flagUncorr;
};


class PadGeometry {
 public:
  /**
  * @brief The constructor for pad geometry
  *
  * Takes coordinates and size of the pad in polar coordinates relative 
  * to calorimeter center.
  *
  * @param r_pad center R-coordinate
  * @param phi_pad center phi-coordinate
  * @param dr_pad radial extent
  * @param dphi_pad polar angle extent
  */
  PadGeometry(double r_pad, double phi_pad, double dr_pad, double dphi_pad);
  ~PadGeometry();

  /**
  * @brief Length of the arc in the pad
  *
  * Calculates arc of a circle of radius r0 with center at (R,Phi)
  * which lies within the pad. Returns the length if any, returns 0 if
  * doesn't intersect.
  *
  * @param R circle center R-coordinate
  * @param Phi circle center phi-coordinate
  * @param r0 circle radius
  *
  * @return Length of the arc within the pad.
  */
  double getArcWithin(double &r0);

  /**
  * @brief Sets pad coordinates to local
  *
  * Coordinates of pad vertices are calculated relative to a point at R, 
  * Phi. This is needed to simplify integration.
  *
  * @param R
  * @param Phi
  */
  void setLocalCoords(const double &R, const double &Phi);

  double m_R;
  double m_phi;
  double m_dR;
  double m_dphi;

  /**
  * @brief Pad center XY coordinate in showere center reference frame
  */
  double m_pcX;
  double m_pcY;

  /**
  * @brief A flag telling wether this pad is one with largest energy
  *
  * If it's true the shower center is within it.
  */
  bool m_isCentral;

 private:
  typedef struct {
    double x1, x2;
    double y1, y2;
    double a, b;
  } PadSide_t;

  bool arcOpenClose(PadSide_t &ps, double &xi, double &yi, double &r0,
         std::pair<double,bool> &ang);

  vector<PadSide_t*> m_sides;

};

