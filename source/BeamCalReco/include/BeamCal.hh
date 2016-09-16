#ifndef BEAMCAL_HH
#define BEAMCAL_HH 1

#include <TROOT.h>
#include <TH3D.h>
#include <TMath.h>
#include <TCanvas.h>

class TCrown;
class TH1D;
class TH2F;
class TH2D;
class TPad;


class BCPadEnergies;
class BeamCalGeo;

#include <EVENT/SimCalorimeterHit.h>

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearMgr.h>

#include <vector>
#include <map>


class BeamCal {

 public:

  // TH3D* createBeamCalHisto();
  // void fillBCHisto(lcio::LCCollection *colBCal, TH3D *histo);
  // void drawBeamCal(TH3D *histo, int layer=-1);
  // void rebinHisto(TH3D *histo);

  BeamCal(const BeamCal&);
  BeamCal& operator=(const BeamCal&);
  BeamCal(BeamCalGeo const& geo);

  TH3D* getBeamCalHistogram(TString title);

  ~BeamCal();

  //  void SetupBeamCal(); 
  void BeamCalDraw(TPad* pad, TH2F* frame);
  void BeamCalDraw(TPad* pad, TH2D* BeamCalEnergy, TH2F *frame);
  void BeamCalDraw(TPad* pad, TH2F* frame, Int_t layerNumber);
  TH1D* BeamCalDrawLayers(TPad *pad, Option_t* options="");
  TH1D* BeamCalDrawRadial(TPad *pad, Option_t* options="");
  TH1D* BeamCalDrawRadial(TPad *pad, Int_t LayerNumber, Option_t* options="");

  void DrawPhiDistributions(TPad *pad, Int_t layer, Option_t* options="");
  void DrawPhiDistributions(TCanvas *canv, Int_t layer, Option_t* options="");

  void BeamCalDraw(TCanvas* canvas, TH2F *frame){ BeamCalDraw((TPad*)(canvas->GetPad(0)), frame); }
  void BeamCalDraw(TCanvas* canvas, TH2D* BeamCalEnergy, TH2F *frame){ BeamCalDraw((TPad*)canvas->GetPad(0), BeamCalEnergy, frame); }
  void BeamCalDraw(TCanvas *canvas, TH2F *frame, Int_t layerNumber) {BeamCalDraw((TPad*)canvas->GetPad(0), frame, layerNumber); }
  TH1D* BeamCalDrawLayers(TCanvas *canvas) { return BeamCalDrawLayers((TPad*)canvas->GetPad(0)); }

  void FillPadFlux(TH1D *BCPadFlux);


  //Geometry functions
  // void SetMaxCutout(Double_t cutout){  m_MinCutoutRadius = cutout; }
  // inline Double_t GetMaxCutout() const { return m_MinCutoutRadius; }
  // inline Double_t GetRadMax() const {  return m_MaxRadius; }
  // inline Double_t GetRadMin() const {  return m_MinRadius; }



  void SetAxisMax(Double_t max) { m_AxisMaximum = max; }
  void SetAxisMin(Double_t min) { m_AxisMinimum = min; }
  void SetLogz    (Int_t logz)  { m_bLogZ = logz; }
  inline Double_t GetAxisMin() const { return m_AxisMinimum; }
  inline Double_t GetAxisMax() const { return m_AxisMaximum; }
  inline Int_t    GetLogz()    const { return m_bLogZ; }

  inline void SetBeamCalHisto(const TH3D *bchisto, bool) {
    if(m_h3BeamCalHisto == bchisto) return;
    if(m_h3BeamCalHisto) delete m_h3BeamCalHisto;
    //      m_h3BeamCalHisto = new TH3D(*bchisto); 
    m_h3BeamCalHisto	     = dynamic_cast<TH3D*>(bchisto->Clone(bchisto->GetName()+TString("Clone"))); 
  }

  void SetBeamCalHisto(const BCPadEnergies *bchisto, TString title="temp");
  void SetBeamCalHisto(const BCPadEnergies *bcpads, const BCPadEnergies *bcErrors, TString title="temp");

  inline TH3D* GetBeamCalHisto(){ return m_h3BeamCalHisto; }
  int GetMaxSegmentation();

  inline bool GetNormalizeByArea() { return m_NormalizeByArea; }
  inline void SetNormalizeByArea(bool norm) { m_NormalizeByArea = norm; }

  inline bool GetNormalizePerYear() { return m_NormalizePerYear; }
  inline void SetNormalizePerYear(bool norm) { m_NormalizePerYear = norm; }
  
  // inline int getNLayers() const { return m_nLayers; }
  // inline std::vector<int> getSegmentVec() const { return m_VecPhiSegments; }
  // inline double getDeadAngle() const { return m_DeadAngle; }


 private:
  BeamCalGeo const& m_BCG;
  Int_t m_bLogZ;
  Int_t m_nSymmetryfold;
  Int_t m_nFullfold;
  // std::vector<int> m_VecPhiSegments;
  // std::vector<double> m_VecRadSegments;
  std::map <int, std::map <int, TCrown> > m_MapCrowns;
  //  Double_t  m_MinRadius, m_MaxRadius, m_MinCutoutRadius;
  //Double_t m_BeamCalDist;
  Double_t m_AxisMinimum, m_AxisMaximum;
  // Double_t m_DeadAngle;
  // Double_t m_UnDeadAngle;
  bool m_NormalizeByArea, m_NormalizePerYear;
  TH3D *m_h3BeamCalHisto;

  Double_t GetPadArea(Int_t cylinder, Int_t sector) const ;
  Double_t GetCylinderArea(Int_t cylinder) const;




};


#endif
