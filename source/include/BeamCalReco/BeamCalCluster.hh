#ifndef BeamCalCluster_HH
#define BeamCalCluster_HH 1

class BCPadEnergies;

#include <map>
#include <iostream>

class BeamCalCluster{
  
public:
  BeamCalCluster():
    m_energy(0),
    m_padIndexInLayer(-1),
    m_clusterPads(),
    m_towerEnergies(),
    m_averagePhi(-9999),
    m_averageRing(-9999),
    m_averageTheta(-9999),
    m_averageZ(-9999),
    m_isReconstructedElectron(false)
  {}

  inline void setPhi(double phi)     { m_averagePhi = phi; };
  inline void setTheta(double theta) { m_averageTheta = theta; };
  inline void setRing(double ring)   { m_averageRing = ring; };
  inline void setZ(double Z)         { m_averageZ = Z; };

  inline double getPhi()   const { return m_averagePhi; };
  inline double getTheta() const { return m_averageTheta; };
  inline double getRing()  const { return m_averageRing; };
  inline double getZ()     const { return m_averageZ; };

  inline void setPadIndexInLayer(int padIn)  { m_padIndexInLayer = padIn; }
  inline int getPadIndexInLayer() const { return m_padIndexInLayer; }

  inline int getNPads() const { return m_clusterPads.size(); }
  inline double getEnergy() const { return m_energy; }

  void addPad(int padIndex, double energy);
  void addPads(const BCPadEnergies& bcp);
  void getBCPad(BCPadEnergies& bcp) const;


  friend std::ostream& operator<<(std::ostream& o, const BeamCalCluster& bcc);

private:
  /// Energy in cluster
  double m_energy;
  /// layer with energy deposit?
  int m_padIndexInLayer;
  /// global padIndex and energy pad
  std::map<int, double> m_clusterPads;
  /// global towerIndex and energy in tower
  std::map<int, double> m_towerEnergies;
  /// energy weighted azimuthal angle
  double m_averagePhi;//energyWeightedAzimuth
  /// energy weighted ring number
  double m_averageRing;
  //energy weighted Polar angle
  double m_averageTheta;

  double m_averageZ;//Should we do this?
  bool m_isReconstructedElectron;

};



#endif
