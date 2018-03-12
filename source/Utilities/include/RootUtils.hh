#ifndef RootUtils_hh
#define RootUtils_hh 1

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>

#include <vector>
#include <string>

class TEfficiency;
class TGraph;
class TGraphErrors;
class TH1;
class TLegend;
class TCanvas;
class TPaletteAxis;

namespace RootUtils {

  void SetStyle(Int_t DefaultFont = 42);

  class Colors {
  private:
    std::vector<Color_t> m_colors;
  public:
    Colors();
    Color_t GetColor();
    void SetColors();
  };//Class Colors

  void SetAllColors( TH1* object, Color_t color );
  void SetAllColors( TEfficiency* object, Color_t color );
  void SetAllColors( TGraph* object, Color_t color );
  void SetAllColors(TGraphErrors* object, Color_t color=kBlack);
  void SetAllColors(TGraphErrors& object, Color_t color=kBlack){ SetAllColors(&object, color); }

  TLegend* BuildLegend(TCanvas& canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2);
  TLegend* BuildLegend(TCanvas& canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString& header);
  TLegend* BuildLegend(TCanvas* canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2);
  TLegend* BuildLegend(TCanvas* canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString& Header);
  TLegend* Build2DLegend(TCanvas* canvas);
  void AddLegendHeader(TLegend* leg, std::vector<std::string> const& lines);
  void SetEntryOptions(TLegend* leg, TString const& options="lpf");

  double AutoSetYRange(TCanvas& canv, double maxScale);


  class PDFFile {
  public:
    PDFFile();
    PDFFile(TString filename);
    ~PDFFile();
    void AddCanvas (TCanvas* canvas, Option_t* title="");
    void AddCanvas (TCanvas& canvas, Option_t* title="") { AddCanvas(&canvas, title); }
    void CloseFile();
  private:
    Bool_t fileOpened;
    TString fileName;
    static int fileNumber;
  };//class

  TPaletteAxis* MovePaletteHorizontally (TH1 *histo, Double_t step);

}//namespace

#endif // RootUtils_hh
