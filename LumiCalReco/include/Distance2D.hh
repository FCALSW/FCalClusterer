#ifndef Distance2D_hh
#define Distance2D_hh 1


#include <cmath>


namespace LCHelper {



  /* --------------------------------------------------------------------------
     calculate the distance between two poins in 2D (in Cartesian coordinates)
     (first index of arrays is X and the second is Y coordinates {or the other way around})
     -------------------------------------------------------------------------- */
  template <class T, class U>
  static double distance2D(const T *vec1, const U *vec2) {
    const double diff0 = vec1[0]-vec2[0];
    const double diff1 = vec1[1]-vec2[1];
    const double distance = sqrt( diff0*diff0 + diff1*diff1 );
    //assert (distance >= 0);
    return distance;
  }

}//namespace

#endif // Distance2D_hh
