#include "RootUtils.hh"

#include <Rtypes.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TList.h>
#include <TObject.h>
#include <TPaletteAxis.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cstring>
#include <iostream>


namespace RootUtils {

  Int_t _DefaultFont = 42;

  /// implement class Colors
  Colors::Colors(): m_colors( std::vector<Color_t>(0) ) {
  }

  void Colors::SetColors() {
    m_colors.push_back(kRed-10);
    m_colors.push_back(kViolet+2);
    m_colors.push_back(kOrange);
    m_colors.push_back(kYellow+4);
    m_colors.push_back(kPink+6);
    m_colors.push_back(kRed+2);
    m_colors.push_back(kBlue);
    m_colors.push_back(kCyan+1);
    m_colors.push_back(kRed-7);
    m_colors.push_back(kGreen+2);
    m_colors.push_back(kBlack);
  }

  Color_t Colors::GetColor() {
    Color_t myColor = kBlack;
    if(m_colors.size() > 0 ) {
      myColor = m_colors.back();
      m_colors.pop_back();
      return myColor;
    } else {
      SetColors();
      return GetColor();
    }
  }
  /// end class Colors


  void SetAllColors( TH1* object, Color_t color ) {
    object->SetLineColor(color);
    object->SetMarkerColor(color);
  }

  void SetAllColors( TEfficiency* object, Color_t color ) {
    object->Paint("");
    //SetAllColors(object->GetPaintedGraph(), color);
    object->GetPaintedGraph()->SetLineColor(color);
    object->GetPaintedGraph()->SetMarkerColor(color);
    object->SetLineColor(color);
    object->SetMarkerColor(color);
  }

  void SetAllColors( TGraph* object, Color_t color ) {
    object->SetLineColor(color);
    object->SetMarkerColor(color);
    object->SetFillColor(kWhite);
  }

  void SetAllColors( TGraphErrors* object, Color_t color ) {
    object->SetLineColor(color);
    object->SetMarkerColor(color);
    object->SetFillColor(kWhite);
  }

  int PDFFile::fileNumber = 0;

  //////////////////////////////////////////////////////////////////////
  //PDFFile CLASS started
  //////////////////////////////////////////////////////////////////////
  PDFFile::~PDFFile() {
    if (fileOpened) {
      PDFFile::CloseFile();
    }
  }

  PDFFile::PDFFile() : fileOpened(false), fileName("") {
    ++fileNumber;
    fileName = Form("myPDFFile_%i.pdf", fileNumber);
  }

  PDFFile::PDFFile(TString filename) : fileOpened(false), fileName(filename) { ++fileNumber; }

  void PDFFile::AddCanvas(TCanvas* canv, Option_t* title) {
    TString myTitle(title);
    myTitle.ReplaceAll(".root", "");
    if (!fileOpened) {
      if (strlen(title)) {
        canv->Print(fileName + TString("("), "pdf,Title:" + myTitle);
      } else {
        std::cout << "No Title" << std::endl;
        canv->Print(fileName + TString("("), "pdf");
      }
      fileOpened = true;
    } else {
      if (strlen(title)) {
        canv->Print(fileName, "pdf,Title:" + myTitle);
      } else {
        canv->Print(fileName, "pdf");
      }
    }
  }  //PDFFile::AddCanvas

  void PDFFile::CloseFile() {
    TCanvas c2("emptyPDFCanv", "arg", 100, 100);
    c2.cd();
    c2.Print(fileName + TString("]"), "pdf");
    fileOpened = false;
  }
  //////////////////////////////////////////////////////////////////////
  //PDFFile CLASS finished
  //////////////////////////////////////////////////////////////////////

  void AddLegendHeader(TLegend* leg, std::vector<std::string> const& lines) {
    int counter(0);
    for (std::vector<std::string>::const_iterator lineIt = lines.begin(); lineIt != lines.end(); ++lineIt) {
      TLegendEntry* newEntry = new TLegendEntry((TObject*)nullptr, (*lineIt).c_str(), "h");
      leg->GetListOfPrimitives()->AddBefore(leg->GetListOfPrimitives()->At(counter), newEntry);
      ++counter;
    }
  }

  void SetEntryOptions(TLegend* leg, TString const& options) {
    TList* list = leg->GetListOfPrimitives();
    for (int i = 0; i < list->GetSize(); ++i) {
      static_cast<TLegendEntry*>(list->At(i))->SetOption(options);
    }
  }

  //  myRoot.BuildLegend(c1, 0.7,0.8,0.9,0.9);
  TLegend* BuildLegend(TCanvas& canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
    return BuildLegend(&canvas, x1, y1, x2, y2);
  }

  TLegend* BuildLegend(TCanvas& canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString& Header) {
    return BuildLegend(&canvas, x1, y1, x2, y2, Header);
  }

  TLegend* BuildLegend(TCanvas* canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
    TLegend* leg = (TLegend*)canvas->BuildLegend(x1, y1, x2, y2);
    if (leg) {
      leg->SetTextFont(_DefaultFont);
      leg->SetFillColor(kWhite);
      leg->SetTextSize(0.06);
      leg->SetMargin(0.15);
      leg->SetBorderSize(0.0);
    }
    return leg;
  }

  TLegend* BuildLegend(TCanvas* canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString& Header) {
    TLegend* leg = (TLegend*)canvas->BuildLegend(x1, y1, x2, y2, Header);
    if (leg) {
      leg->SetTextFont(_DefaultFont);
      leg->SetFillColor(kWhite);
      leg->SetTextSize(0.06);
      leg->SetMargin(0.15);
      leg->SetBorderSize(0.0);
    }
    return leg;
  }

  TLegend* Build2DLegend(TCanvas* canvas) { return BuildLegend(canvas, 0.4, 0.90, 0.98, 0.95); }

  void SetStyle(Int_t DefaultFont) {
    _DefaultFont = DefaultFont;
    gROOT->SetStyle("Plain"); /*Default white background for all plots*/

    /* set bkg color of all to 10: white, but not 0*/

    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetStatColor(10);
    gStyle->SetPadColor(10);

    gStyle->SetFillColor(10);

    gStyle->SetTitleFillColor(0);

    /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
    gStyle->SetPaperSize(20, 26);
    /* No yellow border around histogram*/
    gStyle->SetDrawBorder(0);
    /* remove border of canvas*/
    gStyle->SetCanvasBorderMode(0);
    /* remove border of pads*/
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetLegendBorderSize(0);  //AS:: Only change

    /* default text size*/
    gStyle->SetTextSize(0.05);
    gStyle->SetTitleSize(0.06, "xyz");  //AS changed from 0.07 for now, to be consistent with other plots!
    //    gStyle->SetTitleSize(0.07,"xyz"); //AS changed from 0.07 for now, to be consistent with other plots!
    gStyle->SetLabelSize(0.06, "xyz");

    /* title offset: distance between given text and axis, here x,y,z*/

    gStyle->SetLabelOffset(0.015, "xyz");
    gStyle->SetTitleOffset(1.1, "yz");
    gStyle->SetTitleOffset(1.0, "x");

    //AS AS AS
    gStyle->SetLabelOffset(0.01, "x");   //AS Customized //Default 0.015
    gStyle->SetLabelOffset(0.01, "yz");  //AS Customized
    gStyle->SetLabelOffset(0.005, "y");  //AS Customized

    // Use visible font for all text
    gStyle->SetTitleFont(DefaultFont);
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetStatFont(DefaultFont);
    gStyle->SetStatFontSize(0.07);
    gStyle->SetTextFont(DefaultFont);
    gStyle->SetLabelFont(DefaultFont, "xyz");
    gStyle->SetTitleFont(DefaultFont, "xyz");
    gStyle->SetTitleBorderSize(0);
    gStyle->SetStatBorderSize(1);

    /* big marker points*/
    gStyle->SetMarkerStyle(1);
    gStyle->SetLineWidth(2);
    gStyle->SetMarkerSize(1.2);
    /*set palette in 2d histogram to nice and colorful one*/
    gStyle->SetPalette(kBird, 0);

    /*No title for histograms*/
    gStyle->SetOptTitle(0);
    /* show the errors on the stat box */
    gStyle->SetOptStat(0);
    /* show errors on fitted parameters*/
    gStyle->SetOptFit(0);
    /* number of decimals used for errors*/
    gStyle->SetEndErrorSize(5);

    /* set line width to 2 by default so that histograms are visible when printed small
       idea: emphasize the data, not the frame around*/
    gStyle->SetHistLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFuncWidth(2);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetFuncColor(kRed);
    gStyle->SetLabelColor(kBlack, "xyz");

    //set the margins
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadLeftMargin(0.17);

    //set the number of divisions to show
    gStyle->SetNdivisions(506, "xy");

    //turn off xy grids
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);

    //set the tick mark style
    //    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);

    gStyle->SetCanvasDefW(800);
    gStyle->SetCanvasDefH(700);

    gROOT->ForceStyle();
  }

  double AutoSetYRange(TCanvas& canv, double maxScale) {
    TList* list      = canv.GetListOfPrimitives();
    double maximum   = 0;
    int    firstHist = -1;
    //    int isCanvasLogY = canv.GetLogy();
    for (int iPrims = 0; iPrims <= list->LastIndex(); ++iPrims) {
      TH1* hist = dynamic_cast<TH1*>(list->At(iPrims));
      if (hist) {
        //Remember histo to set maximum of, which is the first one drawn
        if (firstHist == -1) {
          firstHist = iPrims;
        }
        if (hist->GetMaximum() > maximum) {
          maximum = hist->GetMaximum();
        }
      }
    }

    if (firstHist != -1) {
      dynamic_cast<TH1*>(list->At(firstHist))->SetMaximum(maximum * maxScale);
      return maximum * maxScale;
    } else {
      std::cout << __func__ << " No Histograms found" << std::endl;
      return -1;
    }
  }

  TPaletteAxis* MovePaletteHorizontally(TH1 *histo, Double_t step) {
    TPaletteAxis *palette = (TPaletteAxis*)histo->GetListOfFunctions()->FindObject("palette");
    if (palette) {
      palette->SetX1NDC(palette->GetX1NDC()-step);
      palette->SetX2NDC(palette->GetX2NDC()-step);
      palette->GetAxis()->SetLabelOffset(0.01);
    } else {
      std::cout << "Palette not found!"  << std::endl;
    }
    return palette;
  }

}//namespace
