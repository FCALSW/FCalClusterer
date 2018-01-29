#ifndef MCInfo_h
#define MCInfo_h 1

#include <memory>
#include <ostream>

class GlobalMethodsClass;
namespace EVENT{
  class MCParticle;
}

class MCInfo;

using SMCInfo = std::shared_ptr<MCInfo>;

class MCInfo{
public:
  EVENT::MCParticle* mcp           = nullptr;
  double             engy          = 0.0;
  double             theta         = 0.0;
  double             phi           = 0.0;
  int                pdg           = 0;
  int                sign          = 0;
  double             pp[3]         = {0.0, 0.0, 0.0};
  double             mcPosition[3] = {0.0, 0.0, 0.0};  /// position at LumiCal Front Face

  int m_id             = 0;
  int m_parentID       = 0;
  int m_numMCDaughters = 0;

  double m_vtxX = 0.0, m_vtxY = 0.0, m_vtxZ = 0.0, m_endX = 0.0, m_endY = 0.0, m_endZ = 0.0;

  MCInfo() {}

  MCInfo(EVENT::MCParticle* m, double e, double t, double p, double X, double Y, double Z, int P, int S)
      : mcp(m), engy(e), theta(t), phi(p), pdg(P), sign(S), mcPosition{X, Y, Z} {}

  static SMCInfo getMCParticleInfo(EVENT::MCParticle* particle, GlobalMethodsClass& gmc);

  friend std::ostream& operator<<(std::ostream& o, const MCInfo& rhs);

  double x() const { return mcPosition[0]; }
  double y() const { return mcPosition[1]; }
};


#endif // MCInfo_h
