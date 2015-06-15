/**
* @file BeamCalPadGeometry.hh
* @brief BeamCal pad geometry class
* @author Andrey Sapronov <sapronov@cern.ch>
* @version 0.0.1
* @date 2015-03-30
*/

#pragma once

#include <vector>
#include <utility>

class BeamCalPadGeometry {
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
  BeamCalPadGeometry(double r_pad, double phi_pad, double dr_pad, double dphi_pad);
  ~BeamCalPadGeometry();

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
  double getArcWithin(const double &r0);

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

  bool arcOpenClose(PadSide_t &ps, double &xi, double &yi, const double &r0,
         std::pair<double,bool> &ang);

  std::vector<PadSide_t> m_sides;

};

