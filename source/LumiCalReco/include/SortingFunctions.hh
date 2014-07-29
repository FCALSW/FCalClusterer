#ifndef SortingFunctions_hh
#define SortingFunctions_hh 1

#include <IMPL/CalorimeterHitImpl.h>

#include <vector>


/* =========================================================================
   Auxiliary functions
   ----------------------------------------------------------------------------
   Iftach Sadeh - ???? 2007
   ============================================================================ */


/* --------------------------------------------------------------------------
   sorting of clusterCM[id] (for a cluster with Id 'id') with respect to the cluster CM energy
   -------------------------------------------------------------------------- */
//in descending order (highest energy is first)
inline bool clusterCMEnergyCmpDesc( std::vector<double> a, std::vector<double> b ) {
  return a[0] > b[0];
}
//in ascending order (lowest energy is first)
inline bool clusterCMEnergyCmpAsc( std::vector<double> a, std::vector<double> b ) {
  return a[0] < b[0];
}


/* --------------------------------------------------------------------------
   sorting of hits with respect to their energies
   -------------------------------------------------------------------------- */
//in descending order (highest energy is first)
inline bool HitEnergyCmpDesc( IMPL::CalorimeterHitImpl* a, IMPL::CalorimeterHitImpl* b ) {
  return a->getEnergy() > b->getEnergy();
}

//in ascending order (lowest energy is first)
inline bool HitEnergyCmpAsc( IMPL::CalorimeterHitImpl* a, IMPL::CalorimeterHitImpl* b ) {
  return a->getEnergy() < b->getEnergy();
}


/* --------------------------------------------------------------------------
   sorting of hits with respect to their distance from the CM of their cluster
   -------------------------------------------------------------------------- */
//in ascending order (shortest distance is first)
template <int POS>
inline bool HitDistanceCMCmpAsc(  std::vector<double> const&a, std::vector<double> const&b ) {
  return a[POS] < b[POS];
}


// template <class Lhs, class Rhs> inline bool compareByValue(const Lhs& lhs, const Rhs& rhs) {
//   return lhs.second < rhs.second;
// }
template <class T> inline bool compareByValue(const T& lhs, const T& rhs) {
  return lhs.second < rhs.second;
}

#endif // SortingFunctions_hh
