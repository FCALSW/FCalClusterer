/**
* @file BeamCalPadGeometry.cpp
* @brief Pad geometry class to calculate arc of a circle which
* lays within a pad.
* @author Andrey sapronov <sapronov@cern.ch>
* @version 0.0.1
* @date 2015-03-31
*/

#include "BeamCalPadGeometry.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

using std::vector;
using std::pair;

//=================================================================//
//                      BeamCalPadGeometry methods                        //
//=================================================================//

BeamCalPadGeometry::BeamCalPadGeometry(double R, double phi, double dR,  double dphi): 
                         m_R(R),
                         m_phi(phi),
		         m_dR(dR),
		         m_dphi(dphi),
			 m_pcX(0.),
			 m_pcY(0.),
		         m_isCentral(false),
			 m_sides(vector<PadSide_t>())
{
  m_sides.clear();
  m_sides.assign(4, PadSide_t());

  //for (int iside = 0; iside <4 ; iside++) m_sides.push_back(new PadSide_t);
}

BeamCalPadGeometry::~BeamCalPadGeometry()
{
  /*
  while (m_sides.size() != 0 ){
    delete m_sides.back();
    m_sides.pop_back();
  }
  */
}

void BeamCalPadGeometry::setLocalCoords(const double &R, const double &Phi)
{
  // side counter
  int iside = 0;

  // loop over pad vertices
  for (int i = -1; i<=1 ; i+=2){
    for (int j= -1; j<=1 ; j+=2){
      // calculate coordinates
      double r = m_R + i*0.5*m_dR;
      double phi = m_phi - i*j*0.5*m_dphi; // lol
      double x = r*cos(phi) - R*cos(Phi);
      double y = r*sin(phi) - R*sin(Phi);
      // assign them to the ends of pad sides
      m_sides.at(iside).x1 = x;
      m_sides.at(iside).y1 = y;
      if (iside == 0){
        m_sides.at(m_sides.size()-1).x2 = x;
        m_sides.at(m_sides.size()-1).y2 = y;
      } else {
        m_sides.at(iside-1).x2 = x;
        m_sides.at(iside-1).y2 = y;
      }

      iside++;
    }
  }
  
  // calculate line parameters y = a*x + b for every side
  vector<PadSide_t>::iterator it_ps = m_sides.begin();
  for (;it_ps != m_sides.end(); it_ps++){
    if ( fabs(it_ps->x1 - it_ps->x2 ) < 1.e-10)  it_ps->x1+= 1.e-8;
    if ( fabs(it_ps->y1 - it_ps->y2 ) < 1.e-10)  it_ps->y1+= 1.e-8;
    double x1 = it_ps->x1, y1 = it_ps->y1;
    double x2 = it_ps->x2, y2 = it_ps->y2;
    // make small offset, if x-corrdinate accidentally coinside
    it_ps->a = (y2-y1)/(x2-x1);
    it_ps->b = (x2*y1-x1*y2)/(x2-x1);
//    std::cout << (*it_ps)->a<< "\t" <<(*it_ps)->b<< std::endl;
  }

  // also shift pad center to the local coords:
  m_pcX = m_R*cos(m_phi) - R*cos(Phi);
  m_pcY = m_R*sin(m_phi) - R*sin(Phi);
}

double BeamCalPadGeometry::getArcWithin(const double &r0)
{
  // intersections with pad sides
  vector< pair<double, bool> > isects;
  //std::cout << "---------------------" << std::endl;
  double r_vertex_min(1000.);

  // loop over pad sides
  vector<PadSide_t>::iterator it_ps = m_sides.begin();
  for (;it_ps != m_sides.end(); it_ps++){
    // line parameters
    PadSide_t ps = *it_ps;
    double a = ps.a, b = ps.b;

    // check if determinant is > 0. otherwise there is no intersection
    double det = r0*r0*(1+a*a)-b*b;
    if ( det <= 0. ) continue;

    // check if both ends are within the circle, and skip if yes
    double r1 = sqrt(pow(ps.x1,2)+pow(ps.y1,2));
    double r2 = sqrt(pow(ps.x2,2)+pow(ps.y2,2));
    if (r1 < r_vertex_min ) r_vertex_min = r1;
    if (r2 < r_vertex_min ) r_vertex_min = r2;
    		//std::cout << r1 << "\t" <<r2<< "\t" <<r0<<"\n" ;
    if ( r1 <= r0 - 0.01 && r2 <= r0 - 0.01 ) continue;

    // get the roots ( = intersection with line )
    double xi1(0.), xi2(0.), yi1(0.), yi2(0.);
    xi1 = (-a*b - sqrt(det))/(1.+a*a);
    yi1 = a*xi1+b;
    xi2 = (-a*b + sqrt(det))/(1.+a*a);
    yi2 = a*xi2+b;

    pair<double, bool> ang;
    if (arcOpenClose(ps, xi1, yi1, r0, ang)) isects.push_back(ang);
    //std::cout << a<< "\t" <<b<< "\t" <<det << "\t" <<xi1<< "\t" <<yi1<< "\t" << ang.first << "\t" <<ang.second<< std::endl;
    if (arcOpenClose(ps, xi2, yi2, r0, ang)) isects.push_back(ang);
    //std::cout << a<< "\t" <<b<< "\t" <<det << "\t" <<xi2<< "\t" <<yi2<< "\t" << ang.first << "\t" <<ang.second<< std::endl;
  }

  if ( 0 == isects.size() )
    if ( m_isCentral && r0 < r_vertex_min ) return 2.*M_PI*r0;
  if ( 0 == isects.size() ) return 0.;

  // every pad must have even number of intersections
  if ( 0 != isects.size() % 2 ) {
    std::cout << "Warning in BeamCalFitShower: pad shower integration algorithm misbehaved." << std::endl;
    
    
    /*
    for (it_ps = m_sides.begin();it_ps != m_sides.end(); it_ps++){
      PadSide_t ps = *it_ps;
      std::cout << ps.x1<< "\t" <<ps.y1<< "\t" <<ps.x2<< "\t" <<ps.y2<< "\t" <<ps.a<< "\t" <<ps.b<< "\t" <<r0 << std::endl;
    }
    
    std::cout << isects.size();
    vector<pair<double,bool> >::iterator it_is = isects.begin();
    for (; it_is!=isects.end(); it_is++){
      std::cout  << "\t" <<it_is->first-2*M_PI<< "/" << it_is->second ;
    }
    std::cout << std::endl;
    std::cout << m_sides[0].a<< "\t" <<m_sides[0].b << std::endl;
    */
    
  }

  // sort intersections and, if necessary,
  // shift so that first element is the arc opening 
  vector<pair<double, bool> > isects_srt;
  std::sort(isects.begin(), isects.end());
  if (isects.begin()->second != true ) {
    isects_srt.assign(isects.begin()+1, isects.end());
    isects_srt.push_back(*isects.begin());
    isects.assign(isects_srt.begin(),isects_srt.end());
  }

  // calculate total arc length within the pad
  double tot_arc(0.);
  double phi_open(0.);
  vector<pair<double,bool> >::iterator it_is = isects.begin();
  for (; it_is!=isects.end(); it_is++){
    if(it_is->second) phi_open = it_is->first;
    double dphi = it_is->first - phi_open;
    // deal with case when arc opens at pasitive phi, and closes at negative
    if ( dphi < 0. ) dphi+=2.*M_PI; 
    if(! it_is->second) tot_arc += dphi;
  }

  return tot_arc*r0;

}

bool BeamCalPadGeometry::arcOpenClose(PadSide_t &ps, double &xi, 
                               double &yi,const double &r0, pair<double,bool> &ang)
{
    const double phi_oft = 0.001;
    double a = ps.a, b = ps.b;
    ang.first=0., ang.second=false;

    // find if intersection is within the line segment
    if ( (xi > std::min (ps.x1, ps.x2) && xi < std::max (ps.x1, ps.x2)) &&
         (yi > std::min (ps.y1, ps.y2) && yi < std::max (ps.y1, ps.y2)) ) {
      // will need only relative angles, therefore add 2pi to all
      //std::cout << yi<< "\t" <<ps.y1<< "\t" <<ps.y2 << std::endl;
      ang.first = atan2(yi,xi)+2.*M_PI;
      // make small offset along the circle with increase of Phi
      // check that this point is on the same side as the pad center relative to
      // this pad edge
      double xi_oft = r0*cos(ang.first+phi_oft);
      double yi_oft = r0*sin(ang.first+phi_oft);
      if ( (a*xi_oft-yi_oft+b < 0 && a*m_pcX-m_pcY+b < 0) || 
           (a*xi_oft-yi_oft+b > 0 && a*m_pcX-m_pcY+b > 0)) {
        ang.second = true;
      } else ang.second = false;

      return true;
    } else return false;
}
