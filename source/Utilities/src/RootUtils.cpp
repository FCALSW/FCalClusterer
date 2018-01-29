#include "RootUtils.hh"

#include <Rtypes.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>

namespace RootUtils {

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

}//namespace
