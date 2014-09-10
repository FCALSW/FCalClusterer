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

  void getLayerRingPad(int padIndex, int& layer, int& ring, int& pad) const;
  int getPadIndex(int layer, int ring, int pad) const throw(std::out_of_range);
  int getLayer(int padIndex) const;
  int getRing(int padIndex) const;
  int getLocalPad(int padIndex) const;


  bool arePadsNeighbours(int padIndex1, int padIndex2, bool mustBeInSameLayer = false) const;
  double getTransversalDistancePads(int padIndex1, int padIndex2) const;
  double getTransversalDistancePadToPoint(int padIndex, double rho, double phi) const;

  virtual double                getBCInnerRadius()   const = 0;
  virtual double                getBCOuterRadius()   const = 0;
  virtual int                   getBCLayers()        const = 0;
  virtual int                   getBCRings()         const = 0;
  virtual std::vector<double> const&  getPhiSegmentation() const = 0;
  virtual std::vector<double> const&  getRadSegmentation() const = 0;
  virtual std::vector<int>    const&  getNSegments()       const = 0;
  virtual double                getCutout()          const = 0;
  virtual double                getBCZDistanceToIP() const = 0;
  virtual double                getDeadAngle()       const = 0;

  virtual int                   getFirstFullRing()   const;
  virtual double                getFullKeyHoleCutoutAngle() const = 0;
  virtual int                   getPadsBeforeRing( int ring ) const;
  virtual double                getCrossingAngle()   const = 0;

  virtual int getPadsInRing( int ring ) const;
  virtual int getSymmetryFold() const = 0;


  void   getPadExtents(int cylinder, int sector, double *extents) const;
  double getPadMiddlePhi(int cylinder, int sector) const;
  double getPadMiddleTheta(int cylinder, int sector) const;
  double getPadMiddleR(int cylinder, int sector) const;
  double getThetaFromRing(double averageRing) const;
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
