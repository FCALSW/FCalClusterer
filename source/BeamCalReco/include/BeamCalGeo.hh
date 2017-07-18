#ifndef BeamCalGeo_hh
#define BeamCalGeo_hh 1

#include <stdexcept>
#include <vector>


/// ABC for BeamCalGeometry information
class BeamCalGeo {

public:

  virtual ~BeamCalGeo() {};

  virtual int getPadsPerBeamCal() const;
  virtual int getPadsPerLayer() const;

  virtual void getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const;
  virtual int getPadIndex(int layer, int ring, int pad) const;
  virtual int getLayer(int padIndex) const;
  int getRing(int padIndex) const;
  int getLocalPad(int padIndex) const;

  bool arePadsNeighbours(int padIndex1, int padIndex2, bool mustBeInSameLayer = false) const;

  virtual double                getBCInnerRadius()   const = 0;
  virtual double                getBCOuterRadius()   const = 0;
  virtual int                   getBCLayers()        const = 0;
  virtual int                   getBCRings()         const = 0;
  virtual std::vector<double> const&  getPhiSegmentation() const = 0;
  virtual std::vector<double> const&  getRadSegmentation() const = 0;
  virtual std::vector<int>    const&  getNSegments()       const = 0;

  virtual double                getCutout()          const = 0;
  virtual double                getBCZDistanceToIP() const = 0;
  virtual double                getLayerZDistanceToIP(const int lr) const = 0;
  virtual double                getDeadAngle()       const = 0;

  virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const = 0;
  virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const = 0;

  virtual int getPadsInRing( int ring ) const;
  virtual int getSymmetryFold() const = 0;

  double getPadsDistance(int globalPadIndex1, int globalPadIndex2) const;
  void   getPadExtentsById(int globalPadIndex, double *extents) const;

  void   getPadExtents(int cylinder, int sector, double *extents) const;
  double getPadMiddlePhi(int cylinder, int sector) const;
  double getPadMiddleTheta(int layer, int cylinder, int sector) const;
  double getPadMiddleR(int cylinder, int sector) const;
  double getThetaFromRing(int layer, double averageRing) const;
  double getThetaFromRing(int layer, int ring) const;
  double getThetaFromRing(int ring) const;

  double getPadPhi(int ring, int pad) const;
  double getPadPhi(int globalPandIndex) const;

protected:

  //  virtual void countNumberOfPadsInRing() = 0;

  //these are also available from BeamCal Geometry
  // int m_nLayers;
  // int m_nRings;
  //std::vector<int> m_Segments;
  // int m_firstFullRing;
  // double m_fullKeyholeCutoutAngle;
  // int m_nPadsPerLayer;
  // int m_nPadsPerBeamCal;
  // std::vector<int> m_PadsBeforeRing;
  // double mradCrossingAngle;

};//Class


#endif // BeamCalGeo_hh
