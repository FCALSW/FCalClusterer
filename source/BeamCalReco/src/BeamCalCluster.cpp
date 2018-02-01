#include "BeamCalCluster.hh"
#include "BCPadEnergies.hh"

#include <iomanip>
#include <utility>

//void BeamCalCluster::addPads(const BCPadEnergies& /*bcp*/){}

void BeamCalCluster::addPad(int padIndex, double energy) {
  m_clusterPads[padIndex] += energy; //if not existing, create new entry in map, using default c'tor for double, i.e. 0.0
  m_energy += energy;
}

///Add energies in pads from THIS to rhs bcp object
void BeamCalCluster::getBCPad(BCPadEnergies& bcp) const {
  for (std::map<int , double>::const_iterator it = m_clusterPads.begin(); it != m_clusterPads.end(); ++it) {
    bcp.addEnergy(it->first, it->second);
  }
}//getBCPad


std::ostream& operator<<(std::ostream& o, const BeamCalCluster& bcc) {
  o << "Npads"          << std::setw(4)  << bcc.m_clusterPads.size() 
    << "  padID"        << std::setw(5)  << bcc.m_padIndexInLayer
    << "  Energy"       << std::setw(12) << bcc.m_energy << " GeV."
    << "  Phi [Deg]"    << std::setw(8)  << bcc.m_averagePhi
    << "  Ring"         << std::setw(8)  << bcc.m_averageRing
    << "  Theta [mrad]" << std::setw(8)  << bcc.getTheta();
  return o;
}
