#ifndef BCPADENERGIES_HH
#define BCPADENERGIES_HH

#include <map>
#include <stdexcept>
#include <vector>

class BeamCalGeo;
class BeamCalCluster;
class BCPCuts;  
class TH1D;

class BCPadEnergies{

public:

  typedef std::map<int, int> TowerIndexList;
  typedef std::vector<int> PadIndexList;
  typedef std::vector<BeamCalCluster> BeamCalClusterList;

  enum BeamCalSide_t { kUnknown = -1, kLeft = 0 , kRight = 1};

  friend void ClusterNextToNearestNeighbourTowers(BCPadEnergies const& testPads, PadIndexList& myPadIndices, const BCPCuts &cuts, BeamCalClusterList& BeamCalClusters, bool DetailedPrintout=false);

  BCPadEnergies(const BeamCalGeo& bcg, BeamCalSide_t side = kUnknown);
  BCPadEnergies(const BeamCalGeo* bcg, BeamCalSide_t side = kUnknown);
  BCPadEnergies(const BCPadEnergies &bcp);
  BCPadEnergies(const BCPadEnergies *bcp);

  void setEnergy(int layer, int ring, int pad, double energy);
  void setEnergy(int padIndex, double energy);

  void setEnergies(const std::vector<double> &energies);
  void setEnergies(const BCPadEnergies &bcp);
  void addEnergies(const std::vector<double> &energies);
  void addEnergies(const BCPadEnergies &bcp);
  void subtractEnergies(const std::vector<double> &energies);
  void subtractEnergies(const BCPadEnergies &bcp);

  void subtractEnergiesWithCheck(const BCPadEnergies &bcp, const BCPadEnergies &sigma);
  void addEnergiesWithCheck(const BCPadEnergies &bcp, const BCPadEnergies &sigma);

  void addEnergy(int layer, int ring, int pad, double energy);
  void addEnergy(int padIndex, double energy);

  void resetEnergies();
  void scaleEnergies(double factor);

  double getEnergy(int layer, int ring, int pad) const;
  double getEnergy(int padIndex) const;

  double getTotalEnergy() const;

  std::vector<double>* getEnergies();
 
  //Here be our reconstruction functions and algorithms?helper
  BeamCalCluster lookForClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const ;
  BeamCalCluster lookForAlignedClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const ;
  BeamCalCluster lookForNeighbouringClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const ;
  BeamCalCluster lookForNeighbouringClustersOverWithVeto(const BCPadEnergies &background, const BCPCuts &cuts) const ;
  BeamCalClusterList lookForNeighbouringClustersOverWithVetoAndCheck(const BCPadEnergies &background, const BCPadEnergies &backgroundSigma, const BCPCuts &cuts) const ;
  BeamCalClusterList listOfNeighbouringClustersOverWithVeto(const BCPadEnergies &background, const BCPadEnergies &backgroundSigma, const BCPCuts &cuts) const ;

  BeamCalClusterList lookForNeighbouringClustersOverSigma( const BCPadEnergies &backgroundSigma, const BCPCuts &cuts, bool detailedPrintout = false) const;

  // Functions used for the calibration of the particle distinction algorithm
  TowerIndexList* getTopAndNeighbourTowers(double threshold) const;
  TowerIndexList* getTopAndNNNeighbourTowers(double threshold) const;
  int maxDepositTower() const;
  TowerIndexList* getMaxTowerAndWithinRadius(double radius) const;
  void getGlobalCM(double &z, double &rho, double &phi); // z expressed in layers, relative to the front of BeamCal
  TowerIndexList* getTowersWithinRadiusFromPoint(double rho, double phi, double radius) const;
  TH1D* longitudinalProfile() const;
  TH1D* longitudinalProfile(PadIndexList*) const;
  TH1D* longitudinalProfile(TowerIndexList*) const;


  inline void setSide(BeamCalSide_t side) { m_side = side; }
  inline BeamCalSide_t getSide() const { return m_side; }

private:
  //here is information pertinent to the object
  std::vector<double> m_PadEnergies;
  BeamCalSide_t m_side;

  //Reconstruction functions
  PadIndexList getPadsAboveThreshold(double threshold) const;
  PadIndexList getPadsAboveThresholds(const BCPadEnergies& testPads, const BCPCuts& cuts) const;
  PadIndexList getPadsAboveThresholds(const BCPCuts& cuts) const;
  PadIndexList getPadsAboveSigma(const BCPadEnergies& sigmas, const BCPCuts& cuts) const;
  BeamCalCluster getClusterFromPads(const PadIndexList& myPadIndices) const;
  BeamCalCluster getClusterFromAcceptedPads(const BCPadEnergies& testPads, const PadIndexList& myPadIndices, const BCPCuts& cuts) const;


  static TowerIndexList getTowersFromPads( BeamCalGeo const& geo, const PadIndexList& myPadIndices);
  //towerNumber is cellId in Layer i.e. gloadPadID % m_nPadsPerLayer
  static void removeTowerFromPads ( BeamCalGeo const& geo, PadIndexList& myPadIndices, int towerNumber );
  static void removeTowersFromPads( BeamCalGeo const& geo, PadIndexList& myPadIndices, TowerIndexList towerNumbers);

  static PadIndexList getPadsFromTowers ( BeamCalGeo const& geo,
					  PadIndexList const& allPads, 
					  TowerIndexList const& towers);


  std::string streamPad(int padId) const;

public:
  const BeamCalGeo& m_BCG;

}; // class BCPadEnergies

#endif //BCPADENERGIES_hh
