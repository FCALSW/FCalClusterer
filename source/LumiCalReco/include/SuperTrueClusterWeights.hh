#ifndef SuperTrueClusterWeights_hh
#define SuperTrueClusterWeights_hh 1

#include "LCCluster.hh"
#include <string>


/* --------------------------------------------------------------------------
   class and sort rule for computing weights required for assignement of
   reconstructed to true clusters
   -------------------------------------------------------------------------- */
class SuperTrueClusterWeights {

public:
  SuperTrueClusterWeights(int superClusterIdNow,
			  int trueClusterIdNow,
			  LCCluster superClusterCM,
			  LCCluster trueClusterCM);

  double distance2D(double *pos1, double *pos2);
  void	 setWeight(std::string weightMethod);
  void	 setWeight(std::string weightMethod, double minSeparationDistance, double minClusterEngyGeV);

  int	 superClusterId, trueClusterId;
  double distance, deltaEngy, minEngy, weight;

  static inline bool Compare(SuperTrueClusterWeights * a, SuperTrueClusterWeights * b) {
    return a->weight < b->weight;
  }

};


#endif // SuperTrueClusterWeights_hh
