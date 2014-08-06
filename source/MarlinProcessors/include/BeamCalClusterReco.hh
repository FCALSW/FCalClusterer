#ifndef BeamCalClusterReco_h
#define BeamCalClusterReco_h 1

#include <string>
#include <vector>

#include <lcio.h>
#include <marlin/Processor.h>
#include "showerTemplates.hh"

class TChain;
class TEfficiency;
class TFile;
class TH1;
class TH3D;
class TRandom3;
class TString;

class BCPCuts;
class BCPadEnergies;
class BCRecoObject;
class BeamCal;
class BeamCalGeo;

class BeamCalClusterReco : public marlin::Processor , protected ProfileTester {
  
 public:
  
  virtual Processor*  newProcessor() { return new BeamCalClusterReco ; }
  
  
  BeamCalClusterReco() ;

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
  std::vector<std::string> m_files;

  int m_nEvt ;
  int m_specialEvent;
  int m_nBXtoOverlay;
  int m_eventSide;
  int m_minimumTowerSize;
  int m_startLookingInLayer;
  bool m_usePadCuts;
  bool m_createEfficienyFile;
  double m_sigmaCut;
  // SL: Parameters for the particle-type distinction
  double m_EMcorrelThreshold; // Correlation threshold
  double m_eFactor;				// Energy calibration coefficient
  double m_p0a;					// Energy dependence of the parameter a
  double m_p1a;					// Energy dependence of the parameter a
  double m_p0b;					// Energy dependence of the parameter b
  double m_p1b;					// Energy dependence of the parameter b

  std::vector<float> m_startingRings;
  std::vector<float> m_requiredRemainingEnergy;
  std::vector<float> m_requiredClusterEnergy;

  std::vector<double> *m_BeamCalDepositsLeft;
  std::vector<double> *m_BeamCalDepositsRight;
  BCPadEnergies* m_BeamCalAverageLeft;
  BCPadEnergies* m_BeamCalAverageRight;

  BCPadEnergies* m_BeamCalErrorsLeft;
  BCPadEnergies* m_BeamCalErrorsRight;

  //  int _eventid, _nMCP, _MCNumber, bchits, pdg, hitIn, cellID0;

  TRandom3 *m_random3;
  TChain* m_backgroundBX;

  BeamCalGeo *m_BCG;
  BCPCuts* m_bcpCuts;

  TEfficiency *m_thetaEfficieny, *m_phiEfficiency, *m_twoDEfficiency;

  

private:

  void FindOriginalMCParticle(LCEvent *evt, double& impactTheta, double& impactPhi, bool& notFoundAFake);

  void printBeamCalEventDisplay(BCPadEnergies& padEnergies_left, BCPadEnergies& padEnergies_right,
				int maxLayer, double maxDeposit, double depositedEnergy,
				const std::vector<BCRecoObject*> & RecoedObjects) const;

  BCPadEnergies* getBeamCalErrors(const BCPadEnergies *averages, const std::vector<BCPadEnergies*> singles,
				  int numberForAverage );

  std::vector<BCRecoObject*> FindClusters(const BCPadEnergies& signalPads, const BCPadEnergies& backgroundPads, const BCPadEnergies& backgroundSigma, const TString& title);

  void DrawElectronMarkers ( const std::vector<BCRecoObject*> & RecoedObjects ) const;
  void DrawLineMarkers ( const std::vector<BCRecoObject*> & RecoedObjects ) const;

  std::string m_BCalClusterColName;
  std::string m_BCalRPColName;
  std::string m_EfficiencyFileName;

  BeamCalClusterReco(const BeamCalClusterReco&);
  BeamCalClusterReco& operator=(const BeamCalClusterReco&);

  // SL: Added for the purpose of particle-type discrimination
  // Called within FindClusters to generate a histogram with the shower profile
  static TH3D* MakeBeamCalHisto(const BCPadEnergies *bcpads, TString title="temp");

} ;




#endif



