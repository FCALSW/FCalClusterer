#ifndef MCInfo_h
#define MCInfo_h 1

#include <cstddef>

class GlobalMethodsClass;
namespace EVENT{
  class MCParticle;
}

class MCInfo{
public:
  EVENT::MCParticle* mcp;
  double engy;
  double theta;
  double phi;
  double x;
  double y;
  int pdg;
  int sign;
  double pp[3];

  MCInfo(): mcp(NULL), engy(0.0), theta(0.0), phi(0.0),
	    x(0.0), y(0.0), pdg(0), sign(0) {
    pp[0] = 0.0;
    pp[1] = 0.0;
    pp[2] = 0.0;
  }

  MCInfo( EVENT::MCParticle *m, double e, double t, double p,
	  double X, double Y, int P, int S):
    mcp(m), engy(e), theta(t), phi(p), x(X), y(Y), pdg(P), sign(S) {
    pp[0] = 0.0;
    pp[1] = 0.0;
    pp[2] = 0.0;
  }

  static MCInfo getMCParticleInfo( EVENT::MCParticle *particle, GlobalMethodsClass& gmc);

};


#endif // MCInfo_h
