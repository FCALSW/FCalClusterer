#include "MCInfo.h"
#include "GlobalMethodsClass.h"

#include <EVENT/MCParticle.h>

#include <iomanip>

SMCInfo MCInfo::getMCParticleInfo(EVENT::MCParticle* particle, GlobalMethodsClass& gmc) {
  const double LcalZstart = gmc.GlobalParamD[GlobalMethodsClass::ZStart];
  const double Rmin = gmc.GlobalParamD[GlobalMethodsClass::RMin];
  const double Rmax = gmc.GlobalParamD[GlobalMethodsClass::RMax];
  const double thmin = gmc.GlobalParamD[GlobalMethodsClass::ThetaMin];
  const double thmax = gmc.GlobalParamD[GlobalMethodsClass::ThetaMax];

  auto mcp = std::make_shared<MCInfo>();
  // only primary particles wanted
  if( particle->isCreatedInSimulation() and not particle->getParents().empty() ) return mcp; //<--- is primary ?

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

  mcp        = std::make_shared<MCInfo>(particle, engy, theta, phi, begX, begY, LcalZstart, pdg, sign);
  mcp->pp[0] = pp[0];
  mcp->pp[1] = pp[1];
  mcp->pp[2] = pp[2];

  mcp->m_vtxX = double(particle->getVertex()[0]);
  mcp->m_vtxY = double(particle->getVertex()[1]);
  mcp->m_vtxZ = double(particle->getVertex()[2]);

  mcp->m_endX = double(particle->getEndpoint()[0]);
  mcp->m_endY = double(particle->getEndpoint()[1]);
  mcp->m_endZ = double(particle->getEndpoint()[2]);

  // find the original parent of the particle
  EVENT::MCParticle* mcParticleParent = particle;
  mcp->m_parentID                     = particle->id();
  while (1) {
    if (!(mcParticleParent->isCreatedInSimulation()) or mcParticleParent->getParents().empty())
      break;
    mcParticleParent = mcParticleParent->getParents()[0];
    mcp->m_parentID  = (int)mcParticleParent->id();
  }

  // correction for the polar angle for cases where the particle does not
  // come from the IP (according to a projection on the face of LumiCal)
  if (fabs(mcp->m_vtxZ) > _VERY_SMALL_NUMBER) {
    double MCr = tan(mcp->theta) * (gmc.GlobalParamD[GlobalMethodsClass::ZStart] - fabs(mcp->m_vtxZ));
    mcp->theta = atan(MCr / gmc.GlobalParamD[GlobalMethodsClass::ZStart]);
  }

  return mcp;

}


std::ostream& operator<<(std::ostream & o, const MCInfo& rhs) {
  // clang-format off
  o  << "MCInfo: "
     << "Energy " << std::setw(10) << rhs.engy
     << " PDG " << std::setw(5) << rhs.pdg
     << "  Local Pos(theta,phi) =  (" << std::setw(10) << rhs.theta << " , " << std::setw(10) << rhs.phi << " )"
     << "  Start (X, Y) = ( "
     << std::setw(13) << rhs.mcPosition[0] << ", "
     << std::setw(13) << rhs.mcPosition[1]
     << " ) "
     << " Global Momentum (X,Y,Z) =  ( "
     << std::setw(10) << rhs.pp[0] << " , "
     << std::setw(10) << rhs.pp[1] << " , "
     << std::setw(10) << rhs.pp[2] << " )";
  // clang-format on
  return o;

}
