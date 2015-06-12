/**
* @file BeamCalFitShower.cpp
* @brief Implementation of shower finding methods
* @author Andrey sapronov <sapronov@cern.ch>
* @version 0.0.1
* @date 2015-03-31
*/

#include "BeamCalBkg.hh"
#include "BeamCalFitShower.hh"
#include "BeamCalGeoCached.hh"

// ROOT
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <cmath>

using std::vector;
using std::pair;

//=================================================================//
//                    BeamCalFitShower methods                     //
//=================================================================//

BeamCalFitShower::BeamCalFitShower(vector<EdepProfile_t*> &vep, 
		    const BCPadEnergies::BeamCalSide_t bc_side) :
                                  m_vep(vep),
				  m_BCside(bc_side),
				  m_rhom(9.3),
				  m_enTowerLimit(0.)
{
  // hardcode now, make better later
  // m_rhom = 9.3; // Moliere radius (rho_M)
  m_flagUncorr = false;
}

double BeamCalFitShower::fitShower(double &theta, double &phi, double &en_shwr, double &chi2)
{
  // select spot pads around most likely shower
  m_spotPads.clear();
  vector<int> pad_list;
  if ( this->selectSpotPads(pad_list) < 0 ){
    this->deleteSpotPads();

    theta = 0., phi = 0., en_shwr = 0., chi2 = 0.;
    return -1.;
  }

  // calculate inverse covariance matrix for spot pads
  /*
  if ( m_BCbackground->getPadsCovariance(pad_list, m_covInv, m_BCside) < 0 ){
    std::cout << "Falling back to uncorrelated errors in chi2 definition.\n";
    m_flagUncorr = true;
  }
  */
  // not using covariance for now
  m_flagUncorr = true;

  // fit the shower
  ROOT::Minuit2::Minuit2Minimizer minuit ( ROOT::Minuit2::kMigrad );
 
  minuit.SetMaxFunctionCalls(1000);
  minuit.SetMaxIterations(100);
  minuit.SetTolerance(0.1);
  //minuit.SetPrintLevel(3);
 
  const int npar = 4;
  ROOT::Math::Functor f(*this,npar); 
  //ROOT::Math::Functor f(&BeamCalFitShower::showerChi2,npar); 
 
  minuit.SetFunction(f);

  double R0 = m_spotPads.at(0)->pg->m_R;
  double dR0 = m_spotPads.at(0)->pg->m_dR;
  double phi0 = m_spotPads.at(0)->pg->m_phi;
  double dphi0 = m_spotPads.at(0)->pg->m_dphi;

  // estimate initial fit parameters
  double R_shr_center(0.), phi_shr_center(0.);
  double A0(0.), sig0(0.);
  this->estimateShowerPars(R_shr_center, phi_shr_center, A0, sig0);
  
  /*
  std::cout << R_shr_center << "\t" << dR0/20.<< "\t" << R0 - 0.499*dR0<< "\t" << R0 + 0.499*dR0 << std::endl;
  std::cout << phi_shr_center<< "\t" << dphi0/20. << "\t" <<  phi0 - 0.5*dphi0 << "\t" << phi0 + 0.5*dphi0 << std::endl;
  std::cout << A0 << "\t" << A0/20. << "\t" << 0.1*A0 << "\t" << 10.*A0 << std::endl;
  std::cout << sig0<< "\t" << sig0/20.<< "\t" << 0.1*sig0<< "\t" << 10*sig0 << std::endl;
  */
 
  // set variables and their limits
  // make radius and phi limits a bit smaller than the pad size
  minuit.SetLimitedVariable(0,"R", R_shr_center, dR0/20., R0 - 0.499*dR0, R0 + 0.499*dR0);
  minuit.SetLimitedVariable(1,"phi", phi_shr_center, dphi0/20.,  phi0 - 0.5*dphi0, phi0 + 0.5*dphi0);
  minuit.SetLimitedVariable(2,"A", A0 , A0/20., 0.1*A0, 10.*A0);
  minuit.SetLimitedVariable(3,"sig", sig0, sig0/20., 0.1*sig0, 10*sig0);
 
  minuit.Minimize(); 

  // get the minimisation result
  const double *result=minuit.X();
  chi2 = minuit.MinValue();
  /*
  if ( chi2/m_spotPads.size() < 0.1 ) {
    this->deleteSpotPads();
    m_flagUncorr = false;

    return -1;
  }
  */
  //std::cout << result[0]<< "\t" <<result[1]<< "\t" <<result[2]<< "\t" <<result[3] << std::endl;

  // call our chi2 function to calculate energies corresponding to minimum
  (*this)(result);  // wat?

  theta = result[0]/m_BCG->getLayerZDistanceToIP(m_startLayer)*1000.;
  phi = result[1]/M_PI*180.;
  en_shwr = accumulate(m_spotEint.begin(), m_spotEint.end(), 0.0);

  this->deleteSpotPads();
  m_flagUncorr = false;

  // calculate probability that this is our shower
  double prob = TMath::Prob(chi2, m_spotPads.size());

  return prob;
}

int BeamCalFitShower::selectSpotPads(vector<int> &pad_ids)
{
  const double DEGRAD = M_PI/180.;

  // find tower with largest tower chi2, provided it is > 5000.
  vector<EdepProfile_t*>::iterator it_center = m_vep.begin();
  //double max_en = (*it_center)->en;
  double max_chi2 = (*it_center)->chi2;
  double max_en(0.);

  vector<EdepProfile_t*>::iterator it_ep = m_vep.begin();
  for (;it_ep!=m_vep.end(); it_ep++){
    double chi2 = (*it_ep)->chi2;
    double en_tower = (*it_ep)->en - (*it_ep)->bg;
    if ( chi2 > 150. && en_tower > 0.7*m_enTowerLimit && max_chi2 < chi2 ){
      it_center = it_ep;
      max_chi2 = chi2;
      max_en = en_tower;
    }
  }
    //std::cout << max_en<< "\t" << std::endl;

  // no pads with high enough chi2 found:
  if ( it_center == m_vep.begin() ) return -1;

  // create entry for central pad
  double pext[6]; // tsk-tsk-tsk
  m_BCG->getPadExtentsById((*it_center)->id, pext);
  double dphi = pext[3]-pext[2];
  if (dphi < 0 ) dphi+= 360.; // in case pad extents over the -X axis
  (*it_center)->pg = new PadGeometry(pext[4], pext[5]*DEGRAD, 
    pext[1]-pext[0], dphi*DEGRAD);
  (*it_center)->pg->m_isCentral = true;
  m_spotPads.push_back(*it_center);
  pad_ids.push_back((*it_center)->id);

  // collect pads within some Moliere radii of the maximum
  for ( it_ep = m_vep.begin() ;it_ep!=m_vep.end(); it_ep++){
    double pad_dist = m_BCG->getPadsDistance((*it_center)->id, (*it_ep)->id);
    double en_tower = (*it_ep)->en - (*it_ep)->bg;
    if ( pad_dist < 2.0*m_rhom && pad_dist != 0. && en_tower > 0.1*m_enTowerLimit) {
      // if the pad is in the spot, assign its geometry too
      m_BCG->getPadExtentsById((*it_ep)->id, pext);
      double dphi = pext[3]-pext[2];
      if (dphi < 0 ) dphi+= 360.; // in case pad extents over the -X axis
      (*it_ep)->pg = new PadGeometry(pext[4], pext[5]*DEGRAD, pext[1]-pext[0], dphi*DEGRAD);
      (*it_ep)->pg->m_isCentral = false;
      m_spotPads.push_back(*it_ep);
      pad_ids.push_back((*it_ep)->id);
    }
  }

  vector<EdepProfile_t*>::iterator it_sp;
  std::cout << "\t" ;
  for ( it_sp = m_spotPads.begin() ;it_sp!=m_spotPads.end(); it_sp++){
    std::cout << (*it_sp)->en - (*it_sp)->bg<< "\t"  ;
    //std::cout << (*it_sp)->en - (*it_sp)->bg<< "\t"  ;
  }
  std::cout  << std::endl;
  
  
  return m_spotPads.size();
}

//double BeamCalFitShower::showerChi2(double *par)
// this is strangest thing I've ever coded
double BeamCalFitShower::operator()(const double *par)
{
  const int np = m_spotPads.size();
  m_spotEint.assign(np,0.);

  vector<EdepProfile_t*>::iterator it_ep = m_spotPads.begin();
  for (;it_ep != m_spotPads.end(); it_ep++){
    (*it_ep)->pg->setLocalCoords(par[0], par[1]);
  }

  double dr(0.05*m_rhom);
  double ra(0), rb(dr);
  double rmax(3*m_rhom);
  double ga = par[2]*TMath::Gaus(0.,0.,par[3]);
  double gm(0.), gb(0.);
  // simpson integration steps
  while (rb<rmax){
    double rm = (ra+rb)/2.;
    gm=par[2]*TMath::Gaus(rm,0.,par[3]);
    gb=par[2]*TMath::Gaus(rb,0.,par[3]);
    // loop over spot pads
    for (it_ep = m_spotPads.begin();it_ep != m_spotPads.end(); it_ep++){
      // approximate arc at the central point
      double arc = (*it_ep)->pg->getArcWithin(rm);
      m_spotEint.at(it_ep - m_spotPads.begin()) += 
        arc*dr/6.*(ga+4*gm+gb);
    }
    ga = gb;
    ra = rb;
    rb += dr;
  }
  
  // calculate chi2 for integral and actual deposition
  // this piece implements convolution with covariance matrix:
  // chi2 = (Edep-Eint)^T x V x (Edep-Eint)
  // where Edep is a linear algebra vector of energy depositions in spot pads,
  // Eint is the same of energy integral,
  // V is the covariance matrix for the spot pads
  double chi2 = 0.;
  if ( !m_flagUncorr ) {
    for ( int i = 0 ; i < np; i++){
      double ExV(0.);
      for ( int j = 0 ; j < np; j++){
        double d_Ej = m_spotPads[j]->en - m_spotPads[j]->bg - m_spotEint[j];
        ExV += d_Ej*m_covInv[j*np+i];
        //std::cout <<m_spotPads[i]->id<< "\t" <<m_spotPads[j]->id<< "\t" << m_covInv[j*np+i] << std::endl;
      }
      //std::cout <<  std::endl;
      double d_Ei = m_spotPads[i]->en - m_spotPads[i]->bg - m_spotEint[i];
      chi2 += ExV*d_Ei;
    }
  } else {
    for (it_ep = m_spotPads.begin();it_ep != m_spotPads.end(); it_ep++){
      chi2 += pow(((*it_ep)->en - (*it_ep)->bg - m_spotEint.at(it_ep-m_spotPads.begin()))/(*it_ep)->er,2);
    }
  }


  /*
  std::cout << chi2<< "\t" ;
  vector<double>::iterator it_ei = m_spotEint.begin();
  for (;it_ei != m_spotEint.end(); it_ei++){
    std::cout<< "\t"  << *it_ei;
  }
  std::cout << std::endl;
  */
//  std::cout << par[0]<< "\t" <<par[1]<< "\t" <<par[2]<< "\t" <<par[3] << std::endl;

  return chi2;
}

  
void BeamCalFitShower::deleteSpotPads()
{
  // delete pad geometries from spot pads
  vector<EdepProfile_t*>::iterator it_sp;
  for ( it_sp = m_spotPads.begin() ;it_sp!=m_spotPads.end(); it_sp++){
    delete (*it_sp)->pg;
  }
 
  // exclude towers of the found shower from next search
  vector<EdepProfile_t*>::iterator it_ep = m_vep.end()-1;
  for (;it_ep >= m_vep.begin(); it_ep--){
    for (it_sp = m_spotPads.begin();it_sp != m_spotPads.end(); it_sp++){
      if ( (*it_ep)->id == (*it_sp)->id ) m_vep.erase(it_ep);
    }
  }
}

void BeamCalFitShower::estimateShowerPars(double &rc, double &phic, double &A0, double &sig0)
{
  rc = 0.;
  phic = 0.;
  double esum(0.);
  vector<EdepProfile_t*>::iterator it_ep = m_spotPads.begin();
  for (it_ep = m_spotPads.begin();it_ep != m_spotPads.end(); it_ep++){
    double esh = (*it_ep)->en - (*it_ep)->bg;
    rc +=  esh * (*it_ep)->pg->m_R;
    phic += esh * (*it_ep)->pg->m_phi;
    esum += esh;
  }

  if (0. == esum ){
    std::cout << "Warning in BeamCalFitShower: unable to estimate shower center, \n" 
              << "will use hottest pad center coordinates as a fit starting point." << std::endl;
    rc = m_spotPads.at(0)->pg->m_R;
    phic = m_spotPads.at(0)->pg->m_phi;
  } else {
    rc /= esum;
    phic /= esum;
  }

  EdepProfile_t* epc = m_spotPads.at(0); // central pad profile
  sig0 = epc->pg->m_dR*0.66*esum/(epc->en - epc->bg);
  A0 = esum/sqrt(2*M_PI)/pow(sig0,2)/3.;
}

//=================================================================//
//                      PadGeometry methods                        //
//=================================================================//

PadGeometry::PadGeometry(double R, double phi, double dR,  double dphi): 
                         m_R(R),
                         m_phi(phi),
		         m_dR(dR),
		         m_dphi(dphi),
		         m_isCentral(false)
{
  m_sides.clear();
  for (int iside = 0; iside <4 ; iside++) m_sides.push_back(new PadSide_t);
}

PadGeometry::~PadGeometry()
{
  while (m_sides.size() != 0 ){
    delete m_sides.back();
    m_sides.pop_back();
  }
}

void PadGeometry::setLocalCoords(const double &R, const double &Phi)
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
      m_sides.at(iside)->x1 = x;
      m_sides.at(iside)->y1 = y;
      if (iside == 0){
        m_sides.at(m_sides.size()-1)->x2 = x;
        m_sides.at(m_sides.size()-1)->y2 = y;
      } else {
        m_sides.at(iside-1)->x2 = x;
        m_sides.at(iside-1)->y2 = y;
      }

      iside++;
    }
  }
  
  // calculate line parameters y = a*x + b for every side
  vector<PadSide_t*>::iterator it_ps = m_sides.begin();
  for (;it_ps != m_sides.end(); it_ps++){
    if ( abs((*it_ps)->x1 - (*it_ps)->x2 ) < 1.e-10)  (*it_ps)->x1+= 1.e-8;
    if ( abs((*it_ps)->y1 - (*it_ps)->y2 ) < 1.e-10)  (*it_ps)->y1+= 1.e-8;
    double x1 = (*it_ps)->x1, y1 = (*it_ps)->y1;
    double x2 = (*it_ps)->x2, y2 = (*it_ps)->y2;
    // make small offset, if x-corrdinate accidentally coinside
    (*it_ps)->a = (y2-y1)/(x2-x1);
    (*it_ps)->b = (x2*y1-x1*y2)/(x2-x1);
//    std::cout << (*it_ps)->a<< "\t" <<(*it_ps)->b<< std::endl;
  }

  // also shift pad center to the local coords:
  m_pcX = m_R*cos(m_phi) - R*cos(Phi);
  m_pcY = m_R*sin(m_phi) - R*sin(Phi);
}

double PadGeometry::getArcWithin(double &r0)
{
  // intersections with pad sides
  vector< pair<double, bool> > isects;
  //std::cout << "---------------------" << std::endl;

  // loop over pad sides
  vector<PadSide_t*>::iterator it_ps = m_sides.begin();
  for (;it_ps != m_sides.end(); it_ps++){
    // line parameters
    PadSide_t ps = *(*it_ps);
    double a = ps.a, b = ps.b;

    // check if determinant is > 0. otherwise there is no intersection
    double det = r0*r0*(1+a*a)-b*b;
    if ( det <= 0. ) continue;

    // check if both ends are within the circle, and skip if yes
    double r1 = sqrt(pow(ps.x1,2)+pow(ps.y1,2));
    double r2 = sqrt(pow(ps.x2,2)+pow(ps.y2,2));
    if ( r1 <= r0 - 0.01 && r2 <= r0 - 0.01 ) continue;
    if ( r1 >= r0 + 0.01 && r2 >= r0 + 0.01 ) continue;

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

  if ( 0 == isects.size() && m_isCentral ) return 2.*M_PI*r0;
  if ( 0 == isects.size() ) return 0.;

  // every pad must have even number of intersections
  if ( 0 != isects.size() % 2 ) {
    std::cout << "Warning in BeamCalFitShower: pad shower integration algorithm misbehaved." << std::endl;
    
    /*
    for (it_ps = m_sides.begin();it_ps != m_sides.end(); it_ps++){
      PadSide_t ps = *(*it_ps);
      std::cout << ps.x1<< "\t" <<ps.y1<< "\t" <<ps.x2<< "\t" <<ps.y2<< "\t" <<ps.a<< "\t" <<ps.b<< "\t" <<r0 << std::endl;
    }
    
    std::cout << isects.size();
    vector<pair<double,bool> >::iterator it_is = isects.begin();
    for (; it_is!=isects.end(); it_is++){
      std::cout  << "\t" <<it_is->first-2*M_PI<< "/" << it_is->second ;
    }
    std::cout << std::endl;
    std::cout << m_sides[0]->a<< "\t" <<m_sides[0]->b << std::endl;
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

bool PadGeometry::arcOpenClose(PadSide_t &ps, double &xi, 
                               double &yi, double &r0, pair<double,bool> &ang)
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
