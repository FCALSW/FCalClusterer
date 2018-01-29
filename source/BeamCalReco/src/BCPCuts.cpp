#include "BCPCuts.hh"
#include "BeamCalCluster.hh"

bool BCPCuts::isPadAboveThreshold(int padRing, double padEnergy) const {
  for (int i = int(m_startingRings.size())-1; i >= 0; --i) {

    if( padRing >= m_startingRings[i] ) {
      //	std::cout << "Pad Ring " << padRing << " cutRing " << m_startingRings[i] << std::endl;
      if( padEnergy >= m_requiredRemainingEnergy[i] /* GeV */ ) 
      { 
	return true; 
      } else {
	return false;
      }
    }//

  }//for all cuts

  return false;

}//isPadAboveThreshold

bool BCPCuts::isClusterAboveThreshold(BeamCalCluster const& bcc) const {
  for (int i = int(m_startingRings.size())-1; i >= 0; --i) {

    if( bcc.getRing() >= m_startingRings[i] ) {
      //      std::cout << "\nPad Ring " << bcc.getRing() << " cutRing " << m_startingRings[i] << "\n";
      if (bcc.getEnergy() > m_requiredClusterEnergy[i] ) {
	return true;
      } else {
	return false;
      }
      
    }//found the right cut

  }//for all entries

  return false;

}//isClusterAboveThreshold
