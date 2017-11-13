#ifndef SortingFunctions_hh
#define SortingFunctions_hh 1

/* =========================================================================
   Auxiliary functions
   ----------------------------------------------------------------------------
   Iftach Sadeh - ???? 2007
   ============================================================================ */


/* --------------------------------------------------------------------------
   sorting of clusterCM[id] (for a cluster with Id 'id') with respect to the cluster CM energy
   -------------------------------------------------------------------------- */
//in descending order (highest energy is first)
template<class T>
inline bool clusterCMEnergyCmpDesc(T const& a, T const& b) {
  return a[0] > b[0];
}
//in ascending order (lowest energy is first)
template<class T>
inline bool clusterCMEnergyCmpAsc(T const& a, T const& b) {
  return a[0] < b[0];
}


/* --------------------------------------------------------------------------
   sorting of hits with respect to their energies
   -------------------------------------------------------------------------- */
//in descending order (highest energy is first)
template<class T>
inline bool HitEnergyCmpDesc(T const& a, T const& b) {
  return a->getEnergy() > b->getEnergy();
}

//in ascending order (lowest energy is first)
template<typename T>
inline bool HitEnergyCmpAsc(T const& a, T const& b) {
  return a->getEnergy() < b->getEnergy();
}


/* --------------------------------------------------------------------------
   sorting of hits with respect to their distance from the CM of their cluster
   -------------------------------------------------------------------------- */
//in ascending order (shortest distance is first)
template <class T, int POS>
inline bool HitDistanceCMCmpAsc(T const& a, T const& b) {
  return a[POS] < b[POS];
}


// template <class Lhs, class Rhs> inline bool compareByValue(const Lhs& lhs, const Rhs& rhs) {
//   return lhs.second < rhs.second;
// }
template <class T> inline bool compareByValue(const T& lhs, const T& rhs) {
  return lhs.second < rhs.second;
}

#endif // SortingFunctions_hh
