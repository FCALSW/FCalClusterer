#ifndef BCPCUTS_HH
#define BCPCUTS_HH 1

#include <vector>

#include "BeamCalCluster.hh"


/////////////////////////////////////////////////////////////////////////////////////////
// Just a container class for the cuts used in the BCPadEnergies clustering algorithms //
/////////////////////////////////////////////////////////////////////////////////////////

class BCPCuts {

public:


  BCPCuts():
    m_startingRings{0, 7},
    m_requiredRemainingEnergy{0.5, 0.2},
    m_requiredClusterEnergy{3.0, 2.0},
    m_minimumTowerSize(4),
    m_startLookingInLayer( 10 ),
    m_NShowerCountingLayers( 3 ),
    m_useConstPadCuts( false ),
    m_padSigmaCut( 0.0 ),
    m_logWeighting(-1)
    ,m_maxPadDistance(60)
  {
  }

  BCPCuts( std::vector<float> rings, 
	   std::vector<float> pads, 
	   std::vector<float> clusters, 
	   int minimumTowerSize,
	   int startLayer,
	   int countingLayers,
	   bool usePadCuts,
	   double sigmaCut, double logWeighting, double maxPadDistance):
    m_startingRings(rings),
    m_requiredRemainingEnergy(pads),
    m_requiredClusterEnergy(clusters),
    m_minimumTowerSize(minimumTowerSize),
    m_startLookingInLayer( startLayer),
    m_NShowerCountingLayers( countingLayers),
    m_useConstPadCuts(usePadCuts),
    m_padSigmaCut( sigmaCut ),
    m_logWeighting(logWeighting)
    ,m_maxPadDistance(maxPadDistance)
  {}

  bool isPadAboveThreshold(int padRing, double padEnergy) const;
  bool isClusterAboveThreshold(BeamCalCluster const& bcc) const;
  int getMinimumTowerSize() const { return m_minimumTowerSize; }
  int getStartingLayer() const { return m_startLookingInLayer; }
  int getCountingLayers() const { return m_NShowerCountingLayers;}

  bool useConstPadCuts() const { return m_useConstPadCuts; }

  double getPadSigmaCut() const { return m_padSigmaCut; }
  double getLogWeighting() const { return m_logWeighting; }

  BCPCuts& setSigmaCut(double cut) { m_padSigmaCut = cut; return *this;}
  BCPCuts& setStartLayer(int layer) { m_startLookingInLayer = layer; return *this;}
  BCPCuts& setMinimumTowerSize( int tSize) { m_minimumTowerSize = tSize; return *this; }
  BCPCuts& setLogWeighting(double logWeighting) { m_logWeighting = logWeighting; return *this; }

  inline float getMinPadEnergy() const { return m_requiredRemainingEnergy[0]; }
  inline double getMaxPadDistance() const { return m_maxPadDistance; }

private:

  std::vector<float> m_startingRings;
  std::vector<float> m_requiredRemainingEnergy;
  std::vector<float> m_requiredClusterEnergy;

  int m_minimumTowerSize;
  int m_startLookingInLayer;
  int m_NShowerCountingLayers;
  bool m_useConstPadCuts;

  double m_padSigmaCut;
  double m_logWeighting;
  double m_maxPadDistance;
};

#endif
