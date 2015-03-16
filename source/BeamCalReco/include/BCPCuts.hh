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
    m_startingRings(0),
    m_requiredRemainingEnergy(0),
    m_requiredClusterEnergy(0),
    m_minimumTowerSize(4),
    m_startLookingInLayer( 10 ),
    m_useConstPadCuts( false ),
    m_padSigmaCut( 0.0 )
  {
    m_startingRings.push_back(0);  m_requiredRemainingEnergy.push_back(0.2);  m_requiredClusterEnergy.push_back(3.0);
    m_startingRings.push_back(7);  m_requiredRemainingEnergy.push_back(0.5);  m_requiredClusterEnergy.push_back(2.0);
  }
    // m_firstOuterRing({7.0}),
    // m_innerEnergyPad(0.5),
    // m_outerEnergyPad(0.2),
    // m_innerEnergyCluster(3.0),
    // m_outerEnergyCluster(2.0),


  BCPCuts( std::vector<float> rings, 
	   std::vector<float> pads, 
	   std::vector<float> clusters, 
	   int minimumTowerSize,
	   int startLayer,
	   bool usePadCuts,
	   double sigmaCut):
    m_startingRings(rings),
    m_requiredRemainingEnergy(pads),
    m_requiredClusterEnergy(clusters),
    m_minimumTowerSize(minimumTowerSize),
    m_startLookingInLayer( startLayer),
    m_useConstPadCuts(usePadCuts),
    m_padSigmaCut( sigmaCut )
  {}

  bool isPadAboveThreshold(int padRing, double padEnergy) const;
  bool isClusterAboveThreshold(BeamCalCluster const& bcc) const;
  int getMinimumTowerSize() const { return m_minimumTowerSize; }
  int getStartingLayer() const { return m_startLookingInLayer; }

  bool useConstPadCuts() const { return m_useConstPadCuts; }

  double getPadSigmaCut() const { return m_padSigmaCut; }


  BCPCuts& setSigmaCut(double cut) { m_padSigmaCut = cut; return *this;}
  BCPCuts& setStartLayer(int layer) { m_startLookingInLayer = layer; return *this;}
  BCPCuts& setMinimumTowerSize( int tSize) { m_minimumTowerSize = tSize; return *this; }

  inline float getMinPadEnergy() const { return m_requiredRemainingEnergy[0]; }

private:

  std::vector<float> m_startingRings;
  std::vector<float> m_requiredRemainingEnergy;
  std::vector<float> m_requiredClusterEnergy;


  // double m_firstOuterRing;
  // double m_innerEnergyPad;
  // double m_outerEnergyPad;
  // double m_innerEnergyCluster;
  // double m_outerEnergyCluster;
  int m_minimumTowerSize;
  int m_startLookingInLayer;
  bool m_useConstPadCuts;

  double m_padSigmaCut;

};

#endif
