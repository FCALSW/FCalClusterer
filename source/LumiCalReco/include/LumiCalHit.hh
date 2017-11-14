#ifndef LUMICALHIT_HH
#define LUMICALHIT_HH 1

#include "ProjectionInfo.hh"

#include <memory>
#include <vector>

#include <EVENT/CalorimeterHit.h>

/** Class to hold position and energy information of CalorimeterHits in the LumiCal

 * To use of projected hits, which derive from multiple CalorimeterHits this
 * wrapper is necessary to assign the original CalorimeterHits to the created
 * cluster
 */

class LumiCalHit {
public:
  LumiCalHit() {}
  LumiCalHit(int cellIDProjection, ProjectionInfo const& projection)
      : m_cellID0(cellIDProjection),
        m_cellID1(projection.getCellIdHitZ()),
        m_energy(projection.getEnergy()),
        m_position{projection.getPosition()[0], projection.getPosition()[1], projection.getPosition()[2]},
        m_caloHits(projection.getCaloHits()) {}

  ~LumiCalHit() {}

private:
  int    m_cellID0     = 0;
  int    m_cellID1     = 0;
  double m_energy      = 0.0;
  double m_position[3] = {0.0, 0.0, 0.0};
  /// original CalorimeterHits (in the global coordinate system)
  VecLCIOCalHit m_caloHits{};

public:
  void addHit(EVENT::CalorimeterHit* hit) { m_caloHits.insert(hit); }

  auto getHits() const -> decltype(m_caloHits) const & { return (m_caloHits); }
  auto beginHits() const -> decltype(m_caloHits.begin()) { return m_caloHits.begin(); }
  auto endHits() const -> decltype(m_caloHits.end()) { return m_caloHits.end(); }

  void setEnergy(double e) { m_energy = e; }
  void setCellID0(int i) { m_cellID0 = i; }
  void setCellID1(int i) { m_cellID1 = i; }
  void setPosition(double const* pos) {
    m_position[0] = pos[0];
    m_position[1] = pos[1];
    m_position[2] = pos[2];
  }  //fixme memcopy??

  double        getEnergy() const { return m_energy; }
  int           getCellID0() const { return m_cellID0; }
  int           getCellID1() const { return m_cellID1; }
  const double* getPosition() const { return m_position; }
};

#endif  // LUMICALHIT_HH
