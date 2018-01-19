#ifndef BeamCalClusterReco_h
#define BeamCalClusterReco_h 1

#include <lcio.h>
#include <marlin/Processor.h>

#include <string>
#include <vector>
#include <map>

class TChain;
class TEfficiency;
class TFile;
class TH1;
class TRandom3;
class TString;

class BCPCuts;
class BCPadEnergies;
class BCRecoObject;
class BeamCal;
class BeamCalGeo;
class BeamCalBkg;

namespace EVENT {
  class CalorimeterHit;
}

class BeamCalClusterReco : public marlin::Processor {
  
 public:

  //Helper class to store information on MCParticles in the BeamCal
  class OriginalMC {
  public:
    double m_theta; 
    double m_phi; 
    double m_energy;
    bool m_wasFound;
    OriginalMC(double t, double p, double e): m_theta(t), m_phi(p), m_energy(e), m_wasFound(false) {}
    OriginalMC(): m_theta(0.0), m_phi(0.0), m_energy(0.0), m_wasFound(false) {}

  };

  
  virtual Processor*  newProcessor() { return new BeamCalClusterReco ; }
  
  
  BeamCalClusterReco() ;

  BeamCalClusterReco(const BeamCalClusterReco&) = delete;
  BeamCalClusterReco& operator=(const BeamCalClusterReco&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string m_colNameMC;
  std::string m_colNameBCal;
  std::string m_bgMethodName;
  std::string m_hitsOutColName="BeamCal_Hits";
  std::vector<std::string> m_files;

  int m_nEvt ;
  int m_specialEvent;
  int m_nBXtoOverlay;
  int m_eventSide;
  int m_minimumTowerSize;
  int m_startLookingInLayer;
  int m_NShowerCountingLayers;

  bool m_usePadCuts;
  bool m_useChi2Selection;
  bool m_createEfficienyFile;

  double m_sigmaCut;
  double m_TowerChi2ndfLimit;
  double m_calibrationFactor;

  std::vector<float> m_startingRings;
  std::vector<float> m_requiredRemainingEnergy;
  std::vector<float> m_requiredClusterEnergy;

  //  int _eventid, _nMCP, _MCNumber, bchits, pdg, hitIn, cellID0;

  BeamCalGeo *m_BCG;
  BCPCuts* m_bcpCuts;
  BeamCalBkg *m_BCbackground;

  TEfficiency *m_totalEfficiency, *m_thetaEfficieny, *m_phiEfficiency, *m_twoDEfficiency;
  TEfficiency *m_phiFake, *m_thetaFake;
  std::vector<TH1*> m_checkPlots;
  std::vector<OriginalMC> m_originalParticles;
  bool m_MCinBeamCal;
  

private:

  void findOriginalMCParticles(LCEvent *evt);
  void fillEfficiencyObjects(const std::vector<BCRecoObject*>& RecoedObjects);
  void readSignalHits(LCEvent* evt, LCCollection* colBCal, BCPadEnergies& padEnergiesLeft, BCPadEnergies& padEnergiesRight,
                      double& depositedEnergy, double& maxDeposit, int& maxLayer);
  LCCollection* createCaloHitCollection(LCCollection* simCaloHitCollection) const;

  void printBeamCalEventDisplay(BCPadEnergies& padEnergies_left, BCPadEnergies& padEnergies_right,
				int maxLayer, double maxDeposit, double depositedEnergy,
				const std::vector<BCRecoObject*> & RecoedObjects) const;

  std::vector<BCRecoObject*> FindClusters(const BCPadEnergies& signalPads, const BCPadEnergies& backgroundPads, const BCPadEnergies& backgroundSigma, const TString& title);
  std::vector<BCRecoObject*> FindClustersChi2(const BCPadEnergies& signalPads, const BCPadEnergies& backgroundPads, const BCPadEnergies& backgroundSigma, const TString& title);

  void DrawElectronMarkers ( const std::vector<BCRecoObject*> & RecoedObjects ) const;
  void DrawLineMarkers ( const std::vector<BCRecoObject*> & RecoedObjects ) const;

  std::string m_BCalClusterColName;
  std::string m_BCalRPColName;
  std::string m_EfficiencyFileName;
  std::string m_detectorName = "";

  ///Creates the BeamCalGeometry either from DD4hep if compiled with DD4hep and
  ///the geometry is available or from GearFile in all other cases
  BeamCalGeo* getBeamCalGeo();
  bool m_usingDD4HEP;
  std::map<int, std::vector<EVENT::CalorimeterHit*>> m_caloHitMap{};
} ;




#endif



