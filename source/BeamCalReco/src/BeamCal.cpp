#include "BeamCal.hh"
#include "RootUtils.hh"
#include "BCPadEnergies.hh"
#include "BCUtilities.hh"

#include "BeamCalGeo.hh"

#include <gear/GEAR.h>
#include <gear/GearMgr.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>

#include <TCanvas.h>
#include <TCrown.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPave.h>
#include <TStyle.h>

#include <algorithm>
#include <iomanip>
#include <iostream>


BeamCal::BeamCal(BeamCalGeo const& geo):
  m_BCG(geo),
  m_bLogZ(0),
  m_nSymmetryfold(8),
  m_nFullfold(m_nSymmetryfold + 1),
  m_MapCrowns(),
  m_AxisMinimum(1e-3),
  m_AxisMaximum(-1),
  m_NormalizeByArea(false),
  m_NormalizePerYear(false),
  m_h3BeamCalHisto(NULL)
{

  for( int ring=0; ring< m_BCG.getBCRings(); ring++){
    const int padsInThisRing = m_BCG.getPadsInRing(ring);

    for(int phiPad=0; phiPad < padsInThisRing; phiPad++){
      double extents[6];
      m_BCG.getPadExtents(ring, phiPad, extents);
      m_MapCrowns[ring][phiPad] = ( extents[2] < extents[3] ) ?
		   TCrown(0,0,extents[0],extents[1],extents[2],extents[3]):
		   TCrown(0,0,extents[0],extents[1],extents[2],extents[3]+360.0);

    }//all sectors    
  }//all rings

}

BeamCal::~BeamCal() {
  delete m_h3BeamCalHisto;
}

void BeamCal::BeamCalDraw(TPad *pad, TH2D *BeamCalEnergy, TH2F *frame){

  pad->cd();
  pad->SetLogz(m_bLogZ);
  Double_t maxenergy = 0;
  Double_t minenergy = m_AxisMinimum;
  if(m_AxisMaximum > 0) {
    maxenergy=m_AxisMaximum;
    // cout << "FRA " << minenergy << "   to   " << maxenergy << endl;
    frame->SetBinContent(1, minenergy);
    frame->SetBinContent(2, maxenergy);
    frame->SetAxisRange(minenergy,maxenergy,"Z");
    frame->Draw("colz");
    frame->SetContour(99);
    gPad->Update();
  } else {
    BeamCalEnergy->SetContour(99);
    maxenergy = 1.0001*BeamCalEnergy->GetMaximum();
    BeamCalEnergy->SetAxisRange(minenergy, maxenergy, "Z");
    BeamCalEnergy->Draw("colz");
    pad->Update();
  }

  for( int i=0; i<m_BCG.getBCRings(); i++){
    for(int k=0; k<m_BCG.getPadsInRing(i); k++){
      Double_t energy = BeamCalEnergy->GetBinContent(k+1, i+1);
      if(m_NormalizeByArea) {
	energy = energy/GetPadArea(i, k);
      }
      if(m_NormalizePerYear) {
	energy = energy* BCUtil::BXperYear();
      }

      Color_t color=kWhite;
      //----------------------------------------
      if (m_bLogZ==1) {
        const Double_t wlmin   = TMath::Log10(minenergy);
        const Double_t wlmax   = TMath::Log10(maxenergy);
        const Int_t    ncolors = gStyle->GetNumberOfColors();
        const Int_t    ndivz   = 99;
        const Double_t scale   = ndivz / (wlmax - wlmin);
        const Double_t zc      = std::max(std::min(TMath::Log10(energy), wlmax), wlmin);
        const Int_t    color1  = std::min(Int_t(0.01 + (zc - wlmin) * scale), 98);
        const Int_t    color2  = Int_t((color1 + 0.99) * Double_t(ncolors) / Double_t(ndivz));
        color                  = gStyle->GetColorPalette(color2);
        //----------------------------------------
      } else {
	Double_t zc = energy;
	if(energy > 0) {
	  if (zc < minenergy) zc = minenergy;
	  if (zc > maxenergy) zc = 0.9999*maxenergy;
	}
	Int_t temp = (int)(zc / maxenergy * 100);
	if (temp > 99) temp = 99;
	color = gStyle->GetColorPalette(temp);
      }

   //----------------------------------------
      m_MapCrowns[i][k].SetFillColor(color);
      m_MapCrowns[i][k].SetLineColor(kBlack);
      m_MapCrowns[i][k].SetLineWidth(0);
      if(energy == 0) {
	m_MapCrowns[i][k].SetFillColor(kWhite);
      }
    }
  }
  
  pad->Clear();
  pad->SetLogz(m_bLogZ);
  frame->SetAxisRange(minenergy, maxenergy, "Z");
  frame->SetBinContent(1, minenergy);
  frame->SetBinContent(2, maxenergy);
  frame->GetXaxis()->SetTitle("X' [mm]");
  frame->GetYaxis()->SetTitle("Y [mm]");
  
  frame->SetContour(99);
  TH1* frameCopy = frame->DrawCopy("colz");
  pad->Update();
  TPaletteAxis* palette  = (TPaletteAxis*)frameCopy->GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->GetAxis()->SetTitle("Deposited Energy per Pad [GeV]");
    //Move Palette Horizontally 
    const double step = 0.055;
    palette->SetX1NDC(palette->GetX1NDC()-step);
    palette->SetX2NDC(palette->GetX2NDC()-step);
    palette->GetAxis()->SetTitleOffset(1.45);
    palette->GetAxis()->SetLabelOffset(0.01);
  }  else {
    std::cout << "No Palette!!!"  << std::endl;
  }
  pad->Update();

  for(int i=0; i<m_BCG.getBCRings(); i++){
    for(int k=0; k<m_BCG.getPadsInRing(i); k++){
      m_MapCrowns[i][k].DrawClone();
    }
  }

}//BeamCalDraw



void BeamCal::BeamCalDraw(TPad *pad, TH2F *frame){
  TH2D *BeamCalEnergy = (TH2D*)m_h3BeamCalHisto->Project3D("yz");
  BeamCalDraw(pad, BeamCalEnergy, frame);
}

void BeamCal::BeamCalDraw(TPad *pad, TH2F *frame, Int_t layerNumber){
  //m_h3BeamCalHisto->SetAxisRange(10.5, 37.5, "X");
  m_h3BeamCalHisto->SetAxisRange(layerNumber-0.5, layerNumber+0.5, "X");
  TH2D *BeamCalEnergy = (TH2D*)m_h3BeamCalHisto->Project3D("yz");
  BeamCalDraw(pad, BeamCalEnergy, frame);
  //We Reset the Axis Range of the Histogram
  m_h3BeamCalHisto->SetAxisRange(0, m_BCG.getBCLayers(),"X");
}


TH1D* BeamCal::BeamCalDrawLayers(TPad *pad, Option_t* options){
  TH1D *BeamCalEnergy = (TH1D*)m_h3BeamCalHisto->Project3D("x");
  BeamCalEnergy->SetXTitle("Layer");
  BeamCalEnergy->SetYTitle("Energy [GeV/BX]");
  BeamCalEnergy->SetLineWidth(3);
  if(m_NormalizePerYear){
    BeamCalEnergy->SetYTitle("Energy [GeV/yr]");
    BeamCalEnergy->Scale(BCUtil::BXperYear());
  }
  //  std::cout << "BC E Max " <<   BeamCalEnergy->GetMaximum() << std::endl;
  if(m_NormalizeByArea) {
    //Area in cm2
    Double_t area = TMath::Pi()*(m_BCG.getBCOuterRadius()*m_BCG.getBCOuterRadius()-m_BCG.getBCInnerRadius()*m_BCG.getBCInnerRadius())/100.;
    std::cout << "AREA  " << area << std::endl;
    BeamCalEnergy->Scale(1./area);
    BeamCalEnergy->SetXTitle("Layer");
    BeamCalEnergy->SetYTitle("1 MeV eq. Netruon Flux #[]{1/BX/cm^{2}}");
    if(m_NormalizePerYear){
      BeamCalEnergy->SetYTitle("1 MeV eq. Neutron Flux #[]{1/yr/cm^{2}}");
    }
  }
  //  std::cout << "BC E MAx " <<   BeamCalEnergy->GetMaximum() << std::endl;
  pad->cd();
  BeamCalEnergy->Draw(options);
  return BeamCalEnergy;
}

//This returns the Area of the pad in cm^2
Double_t BeamCal::GetPadArea(Int_t cylinder, Int_t sector) const {
  Double_t area = -1;
  //  std::map<int, std::map<int, TCrown> >::iterator crownIt = ((m_MapCrowns.find(cylinder)));
  //std::map<int, TCrown>  crownIt = ((m_MapCrowns.find(cylinder)))->second;
  //Changed this from m_MapCrowns[cylinder][sector], maybe it should be find(sector) first?
  //const TCrown& crown = ((m_MapCrowns[cylinder][sector]));
  auto const& cylinderMap = m_MapCrowns.find(cylinder);
  auto crownIt = cylinderMap->second.find(sector);
  if(crownIt != cylinderMap->second.end()) {
    auto & crown = crownIt->second;
    //Area is Rout*rout - rin*rin *DeltaPhi[rad]/2 in cm2
    Double_t rin = crown.GetR1()/10;
    Double_t rout = crown.GetR2()/10;
    Double_t deltaphi = crown.GetPhimax() - crown.GetPhimin();
    area = (rout*rout - rin*rin) * (TMath::DegToRad()*deltaphi/2.);
    // std::cout << "Rin  " << rin 
    // 	    << "  Rout " << rout 
    // 	    << "  DPho " << deltaphi 
    // 	    << std::endl;
    // std::cout << "Area " << cylinder << "  " << sector << " -  " << area << std::endl;
  } else {

  }
  return area;
}


Double_t BeamCal::GetCylinderArea(Int_t cylinder) const {

  Double_t OuterRadius,InnerRadius;
  InnerRadius = m_BCG.getBCInnerRadius()+cylinder*(m_BCG.getBCOuterRadius()-m_BCG.getBCInnerRadius())/m_BCG.getBCRings();
  OuterRadius = m_BCG.getBCInnerRadius()+(cylinder+1)*(m_BCG.getBCOuterRadius()-m_BCG.getBCInnerRadius())/m_BCG.getBCRings();
  return  (OuterRadius*OuterRadius-InnerRadius*InnerRadius)*TMath::Pi()/100.;
}





void BeamCal::FillPadFlux(TH1D *BCPadFlux) {
  BCPadFlux->SetYTitle("N");
  BCPadFlux->SetXTitle("N_{1 MeV eq.}");
  BCPadFlux->SetLineWidth(3);
  Double_t factor = 1.;
  if(m_NormalizePerYear){
    BCPadFlux->SetXTitle("N_{1 MeV eq.} #[]{1/yr}");
    factor *= BCUtil::BXperYear();
  } 
  if(m_NormalizeByArea) {
    BCPadFlux->SetXTitle("1 MeV eq. Neutron Fluence #[]{1/cm^{2}}");
  }
  if(m_NormalizeByArea && m_NormalizePerYear){ 
    BCPadFlux->SetXTitle("1 MeV eq. Neutron Flux #[]{1/yr/cm^{2}}");
  }

  for( int h=1; h<=m_h3BeamCalHisto->GetNbinsX(); ++h) {
    for( int i=1; i<=m_BCG.getBCRings(); ++i) {
      for(int k=1; k<=m_BCG.getPadsInRing(i-1); ++k) {
	if(m_NormalizeByArea) {  
	  BCPadFlux->Fill(m_h3BeamCalHisto->GetBinContent(h,i,k)*factor/GetPadArea(i-1,k));
	} else {
	  BCPadFlux->Fill(m_h3BeamCalHisto->GetBinContent(h,i,k)*factor);
	}
      }
    }
  }//Layers	
  
  
}//PadFlux



TH1D* BeamCal::BeamCalDrawRadial(TPad *pad, Option_t* options){
  TH1D *BeamCalEnergy = (TH1D*)m_h3BeamCalHisto->Project3D("y");
  BeamCalEnergy->SetXTitle("Cylinder");
  BeamCalEnergy->SetYTitle("Energy [1/BX]");
  BeamCalEnergy->SetLineWidth(3);
  if(m_NormalizePerYear){
    BeamCalEnergy->Scale(BCUtil::BXperYear());
    BeamCalEnergy->SetYTitle("Energy [1/yr]");
  }
  if(m_NormalizeByArea) {
    //Area in cm2
    BeamCalEnergy->SetYTitle("Flux #[]{kGy/yr}");
    for(Int_t i = 1; i <= BeamCalEnergy->GetNbinsX(); ++i)
      {
	BeamCalEnergy->SetBinContent(i, BeamCalEnergy->GetBinContent(i)/GetCylinderArea(i-1));
      }
  }

  pad->cd();
  BeamCalEnergy->Draw(options);
  return BeamCalEnergy;
}


TH1D* BeamCal::BeamCalDrawRadial(TPad *pad, Int_t layerNumber, Option_t* options){
  m_h3BeamCalHisto->SetAxisRange(layerNumber-0.5, layerNumber+0.5, "X");
  TH1D* hPointer = BeamCalDrawRadial(pad, options);
  m_h3BeamCalHisto->SetAxisRange(0, 40,"X");
  return hPointer;
}

void BeamCal::DrawPhiDistributions(TCanvas *canv, Int_t layer, Option_t* options){
  DrawPhiDistributions((TPad*)canv->GetPad(0), layer, options);
}

void BeamCal::DrawPhiDistributions(TPad *pad, Int_t layer, Option_t* options){
  std::vector<TH1D> phiHistos;
  
  //Create the histograms
  RootUtils::Colors mycols;
  for (Int_t i = 0; i < m_BCG.getBCRings(); ++i) {
    Int_t bins = m_BCG.getPadsInRing(i);
    Double_t extents[6];
    this->m_BCG.getPadExtents(i, 1, extents);
    if( not (extents[0] > this->m_BCG.getCutout()) ){
#pragma message "FixMe number of bins in drawphidistributions"
      //calculate bins for the limited range
      const double degreesDeadAngle = m_BCG.getDeadAngle()*TMath::RadToDeg();
      double limitedBinSize =  (360.0-degreesDeadAngle) / double(m_BCG.getPadsInRing(i));
      std::vector<double> binning;
      //leave a gap between 175 and 185..., so zero to

      for (int l = 0; l <= bins/2; ++l) {
      	binning.push_back(limitedBinSize*l);
      	binning.push_back( 180.0 + degreesDeadAngle/2.0 + limitedBinSize*l );
      }

      std::sort(binning.begin(),binning.end());

      // for (std::vector<double>::iterator it = binning.begin(); it != binning.end(); ++it) {
      // 	std::cout << *it  << std::endl;
      // }

      phiHistos.push_back(TH1D(Form("h1EnergyPhi_%i_%i", layer+1, i+1),Form("Ring %i", i+1),
			       bins+1, &binning[0]  ));


    } else {
      phiHistos.push_back(TH1D(Form("h1EnergyPhi_%i_%i", layer+1, i+1),Form("Ring %i", i+1),
			       bins, 0, 360));
    }

    if(strstr(options,"errors")) {
      phiHistos.back().Sumw2();
    }
      
    RootUtils::SetAllColors( &(phiHistos[i]), mycols.GetColor() );
  }

  phiHistos[0].SetXTitle("Angle [deg]");
  phiHistos[0].SetYTitle("Deposited Energy [GeV]");
  
  //Now fill them with the distributions
  Double_t maximum = -1;
  for(int j = 1; j <= m_BCG.getBCRings(); ++j) { //cylinder
    for(int k = 1; k <= m_BCG.getPadsInRing(j-1); ++k) { //sector
      const Double_t binContent = m_h3BeamCalHisto->GetBinContent(layer, j, k);
      const Double_t binError   = m_h3BeamCalHisto->GetBinError(layer, j, k);

      if(binContent > 0 ) {
	Double_t phiBin = this->m_BCG.getPadMiddlePhi(j-1, k-1);
	// if( j == 9 ) {
	//   int binnumber = phiHistos[j-1].FindBin(phiBin);
	//   std::cout << phiBin
	// 	    << "  " << phiHistos[j-1].FindBin(phiBin) 
	// 	    << "  "  << phiHistos[j-1].GetBinLowEdge(binnumber) 
	// 	    << "  " << binContent
	// 	    << std::endl;
	// }
	//	phiHistos[j-1].SetBinContent(k, binContent);
	phiHistos[j-1].Fill(phiBin, binContent);
	// if(j == 1) {
	//   std::cout << "DrawPhi " << std::setw(4) << k << std::setw(15) << phiBin << std::endl;
	// }
	if(strstr(options, "errors") ) {
	  phiHistos[j-1].SetBinError( phiHistos[j-1].FindBin(phiBin), binError);
	}
	if( binContent > maximum) maximum = binContent;

	// if ( j == 1 && k == 32  ) {
	//   std::cout << "Content in 0:32 " << binContent  << std::endl;
	//   std::cout << phiHistos[j-1].GetNbinsX()  << std::endl;
	//   std::cout << phiHistos[j-1].GetBinLowEdge(33)  << std::endl;
	// }

	// std::cout << this->GetPadMiddlePhi(j, k)  
	// 	  << std::setw(20) << binContent
	// 	  << std::endl;
      }//if there is a deposit
    }//sectors
  }//cylinders

  //Now Draw them all
  pad->cd();

  TH1F frame("BCPhiframe","BCPhiFrame",10, 0, 360);
  frame.SetTitle(Form("Layer %i", layer));
  frame.SetXTitle("Azimuthal Angle #phi [deg]");
  frame.SetYTitle("Deposited Energy [GeV]");
  frame.SetMaximum(maximum*1.05);
  frame.SetMinimum(0.0);
  frame.SetAxisRange(0,360,"X");
  if(!strstr(options, "same") ) {//if options contains same, dont draw it
    frame.DrawCopy();
  }

  for (std::vector<TH1D>::iterator it = phiHistos.begin(); it != phiHistos.end(); ++it) {
    (it)->SetLineWidth(1);
    if(strstr(options, "l2")) {
      (it)->SetLineWidth(2);
    }
    (it)->SetMarkerStyle(kDot);
    (it)->SetFillColor(kWhite);
    if( strstr(options, "dotted") ) {
       (it)->SetLineStyle(kDotted);
    } else if ( strstr(options,"dashed") ) {
       (it)->SetLineStyle(kDashed);
    }

    (it)->DrawCopy("same");
    
    // std::cout << __func__ << " Entries " 
    // 	      << std::setw(15) << it->GetEntries() 
    // 	      << std::setw(15) << it->Integral() 
    // 	      << std::endl;
  }

}

void BeamCal::SetBeamCalHisto(const BCPadEnergies *bcpads, TString title){
  if(m_h3BeamCalHisto) delete m_h3BeamCalHisto;
  m_h3BeamCalHisto = getBeamCalHistogram(title);

  int layer=-1, cylinder=-1, sector=-1;
  double energy;
  for (int i = 0; i < m_BCG.getPadsPerBeamCal() ;++i) {
    try {
      m_BCG.getLayerRingPad(i, layer, cylinder, sector);
      energy = bcpads->getEnergy(i);
      m_h3BeamCalHisto->Fill(layer, cylinder, sector, energy);
    }  catch(std::out_of_range &e) {
      std::cout << "ERROR, out of range! " << e.what() 
		<< std::setw(5) << layer
		<< std::setw(5) << cylinder
		<< std::setw(5) << sector
		<<std::endl;
    }catch (std::logic_error &e) {
      std::cout << "ERROR, getLayerRingPad is faulty! " << e.what() 
		<< std::setw(5) << layer
		<< std::setw(5) << cylinder
		<< std::setw(5) << sector
		<<std::endl;
    }
    

  }//for all the pads
  
}//SetBeamCalHisto


void BeamCal::SetBeamCalHisto(const BCPadEnergies *bcpads, const BCPadEnergies *bcErrors, TString title){
  if(m_h3BeamCalHisto) delete m_h3BeamCalHisto;
  m_h3BeamCalHisto = getBeamCalHistogram(title);

  int layer=-1, cylinder=-1, sector=-1;
  for (int i = 0; i < m_BCG.getPadsPerBeamCal() ;++i) {
    try {
      double energy, error;
      m_BCG.getLayerRingPad(i, layer, cylinder, sector);
      energy = bcpads->getEnergy(i);
      error =  bcErrors->getEnergy(i);
      m_h3BeamCalHisto->Fill(layer, cylinder, sector, energy);
      Int_t binNumber = m_h3BeamCalHisto->FindBin(layer,cylinder,sector);
      m_h3BeamCalHisto->SetBinError(binNumber, error);
    } catch(std::out_of_range &e) {
      std::cout << "ERROR, out of range! " << e.what() 
		<< std::setw(5) << layer
		<< std::setw(5) << cylinder
		<< std::setw(5) << sector
		<<std::endl;
    } catch (std::logic_error &e) {
      std::cout << "ERROR, getLayerRingPad is faulty! " << e.what() 
		<< std::setw(5) << layer
		<< std::setw(5) << cylinder
		<< std::setw(5) << sector
		<<std::endl;
    }
    

  }//for all the pads
  
}//SetBeamCalHisto


TH3D* BeamCal::getBeamCalHistogram(TString title){
  //#pragma "FIXME: Fix the getBeamCalHistogram Function"

  // std::cout << "Creating new Histogram with "
  // 	    << std::setw(10) << m_nLayers + 1
  // 	    << std::setw(10) << m_BCG.getBCRings()+1
  // 	    << std::setw(10) << 160
  // 	    << std::endl;
TH3D *histo = new TH3D(title, title,
		       m_BCG.getBCLayers()+1, 0.5, m_BCG.getBCLayers() + 1.5, //layers start at 1
		       m_BCG.getBCRings(),  -0.5, m_BCG.getBCRings()-0.5,      //cylinder, start at 0
		       m_BCG.getPadsInRing( m_BCG.getBCRings()-1), -0.5, m_BCG.getPadsInRing( m_BCG.getBCRings()-1)-0.5 //sectors // start at 0
		  );
 histo->SetDirectory(0);
 return histo;
}


 BeamCal::BeamCal(const BeamCal& rhs):
   m_BCG(rhs.m_BCG),
   m_bLogZ(rhs.m_bLogZ),
   m_nSymmetryfold(rhs.m_nSymmetryfold),
   m_nFullfold(rhs.m_nFullfold),
   m_MapCrowns(rhs.m_MapCrowns),
   m_AxisMinimum(rhs.m_AxisMinimum),
   m_AxisMaximum(rhs.m_AxisMaximum),
   m_NormalizeByArea(rhs.m_NormalizeByArea),
   m_NormalizePerYear(rhs.m_NormalizePerYear),
   m_h3BeamCalHisto(NULL) {
   this->SetBeamCalHisto(rhs.m_h3BeamCalHisto, true);
 }


 BeamCal& BeamCal::operator=(const BeamCal& rhs) {
   if( this == &rhs ) { return *this; } //self-assignment

   m_bLogZ	      = rhs.m_bLogZ;
   m_nSymmetryfold    = rhs.m_nSymmetryfold;
   m_nFullfold	      = rhs.m_nFullfold;
   m_MapCrowns	      = rhs.m_MapCrowns;
   m_AxisMinimum      = rhs.m_AxisMinimum;
   m_AxisMaximum      = rhs.m_AxisMaximum;
   m_NormalizeByArea  = rhs.m_NormalizeByArea;
   m_NormalizePerYear = rhs.m_NormalizePerYear;


   m_h3BeamCalHisto=NULL;
   this->SetBeamCalHisto(rhs.m_h3BeamCalHisto, true);

   return *this;

 }//assignment op


int BeamCal::GetMaxSegmentation() {  return m_BCG.getPadsInRing( m_BCG.getBCRings()-1 ); }
