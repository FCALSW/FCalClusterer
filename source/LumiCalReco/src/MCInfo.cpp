#include "MCInfo.h"
#include "GlobalMethodsClass.h"

#include <EVENT/MCParticle.h>

#include <iomanip>


MCInfo MCInfo::getMCParticleInfo(EVENT::MCParticle *particle, GlobalMethodsClass& gmc) {

  const double LcalZstart = gmc.GlobalParamD[GlobalMethodsClass::ZStart];
  const double Rmin = gmc.GlobalParamD[GlobalMethodsClass::RMin];
  const double Rmax = gmc.GlobalParamD[GlobalMethodsClass::RMax];
  const double thmin = gmc.GlobalParamD[GlobalMethodsClass::ThetaMin];
  const double thmax = gmc.GlobalParamD[GlobalMethodsClass::ThetaMax];

  MCInfo mcp;
  // only primary particles wanted
  if( particle->isCreatedInSimulation() ) return mcp; //<--- is primary ?

  int pdg = particle->getPDG();
  // skip neutrinos
  if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) return mcp; //<--- detectable ?

  // energy above min
  const double engy = particle->getEnergy();
  if( engy < gmc.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV] ) return mcp;  //<--- Energy >  Emin ?

  // check if within Lcal acceptance
  const double *endPoint = particle->getEndpoint();
  // did it reach at least LCal face
  if( fabs( endPoint[2] ) < LcalZstart && fabs( endPoint[2] )>0. ) return mcp; //<--- reached LCal ?

  const double* pp = particle->getMomentum();
  const int sign = int ( pp[2] > 0 ? 1 : -1 );

  // boost to lumical system
  const double pxloc = - gmc.GlobalParamD[GlobalMethodsClass::BetaGamma]*engy + 
    gmc.GlobalParamD[GlobalMethodsClass::Gamma]*pp[0];
  const double *vx = particle->getVertex();

  // particle position at LCal face ( local )
  double begX = vx[0] + double(sign)*LcalZstart*tan(pxloc/(pp[2]));
  double begY = vx[1] + double(sign)*LcalZstart*tan(pp[1]/(pp[2]));
  double rt = sqrt( begX*begX  + begY*begY );
  if( rt < Rmin || Rmax < rt  ) return mcp; //<--- within geo acceptance ?

  //  polar angle ( local )
  double theta = atan( sqrt( pxloc*pxloc + pp[1]*pp[1] ) / fabs(pp[2]));
  if( fabs(theta) < thmin  || thmax < fabs(theta)  ) return mcp; //<--- within theta range ?

  const double phi = atan2( begY, begX );

  mcp = MCInfo( particle, engy, theta, phi, begX, begY, pdg, sign );
  mcp.pp[0] = pp[0];
  mcp.pp[1] = pp[1];
  mcp.pp[2] = pp[2];

  return mcp;

}


std::ostream& operator<<(std::ostream & o, const MCInfo& rhs) {
  o  << "MCInfo: "
     << "  Energy "              << std::setw(10) << rhs.engy
     << " PDG " << std::setw(5) << rhs.pdg
     << "  pos(theta,phi) =  ( " << std::setw(10) << rhs.theta << " , " << std::setw(10) << rhs.phi << " )"
     << " Start (X, Y) = ( "
     << std::setw(10) << rhs.x
     << std::setw(10) << rhs.y
     << " ) "
     << "  Momentum (X,Y,Z) =  ( "
     << std::setw(10) << rhs.pp[0] << " , "
     << std::setw(10) << rhs.pp[1] << " , "
     << std::setw(10) << rhs.pp[2] << " )";
  return o;


}
