/****************************************************************************
 *
 * longitudinals.cpp
 *
 * Create longitudinal projections of multiple BeamCal histograms
 * and store them in a file.
 *
 *  Strahinja Lukic, June 2014
 *
 ****************************************************************************/

#include "BeamCalGeo.hh"
#include "BeamCalGeoGear.hh"
#include "BeamCalGeoCached.hh"
#include "BeamCal.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"

//GEAR
#include <gearxml/GearXML.h>
#include <gear/GearMgr.h>

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStyle.h>

#include <string>
#include <iostream>
#include <sstream>


int main (int argn, char **argc) {

  if ( argn < 3 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "longitudinals <GearFile1> <PRootFile> "  << std::endl;
    std::exit(1);
  } 

//  const double ares = 0.03; //Constant relating the width of the deposit distribution in a cell to the square root of the deposition
  const double e0 = 1.e-3; // Effective "unit energy" to express the cell deposit distribution as Poissonian

  //Create the gearmanager from the gearfile
  std::string gearFile ( argc[1] );
  gear::GearXML gearXML( gearFile ) ;
  gear::GearMgr* gearMgr = gearXML.createGearMgr() ;
  BeamCalGeo* geoCache = new BeamCalGeoCached (gearMgr);
  //  BeamCalGeo* geoGear = new BeamCalGeoGear (gearMgr);
  std::vector<BCPadEnergies> *signalBeamCals = NULL;

  std::string rootFile ( argc[2]);
  int wrong(0);
  try {
    	signalBeamCals = BCUtil::ReadMultiRootFile(rootFile, geoCache);
    } catch (std::invalid_argument &e) {
      ++wrong;
  }
  if (!signalBeamCals) {
      std::cerr << "This file has no BCPadEnergies as a tree" << std::endl;
      return 1;
  }

  TString outName("longitudinals."); outName += rootFile;
  TFile histograms(outName.Data(), "RECREATE");

  BeamCal bc (gearMgr);

  int nEntries = signalBeamCals->size();

  bc.SetBeamCalHisto( &(signalBeamCals->at(0)) );
  TH1D *BeamCalEnergy = (TH1D*)bc.GetBeamCalHisto()->Project3D("x");
  TH1D *BCaverage = new TH1D(*BeamCalEnergy);
  BCaverage->Reset();
  int nAdded=0;
  if(BeamCalEnergy->Integral() > 0.001)
	  { BeamCalEnergy->Write("Layers_0"); BCaverage->Add(BeamCalEnergy); nAdded++;}

// Loop to create the average profile

  for( int iEvt=1; iEvt<nEntries; iEvt++)
  {
//	  std::cout << "Sum integral: " << BCaverage->Integral() << std::endl;
	  bc.SetBeamCalHisto( &(signalBeamCals->at(iEvt)) );
	  BeamCalEnergy = (TH1D*)bc.GetBeamCalHisto()->Project3D("x");
	  std::stringstream hname; hname << "Layers_" << iEvt;
	  if(BeamCalEnergy->Integral() > 1.e-6)
	  {
		  BeamCalEnergy->Write(hname.str().c_str());
		  BCaverage->Add(BeamCalEnergy);
		  nAdded++;
	  }
  }

  BCaverage->Scale(1./double(nAdded));
  TH1D *BCrms = new TH1D(*BCaverage);
  BCrms->Reset();
  TH1D BCtemp(*BCrms);
  TH1D BCsq(*BCrms);

// Loop to create the RMS of the individual bins
  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
	  bc.SetBeamCalHisto( &(signalBeamCals->at(iEvt)) );
	  BeamCalEnergy = (TH1D*)bc.GetBeamCalHisto()->Project3D("x");
	  if(BeamCalEnergy->Integral() > 1.e-6)
	  {
		  BCtemp.Add(BeamCalEnergy, BCaverage, 1., -1.);
		  BCsq.Multiply(&BCtemp,&BCtemp);
		  BCrms->Add(&BCsq);
	  }
  }
  BCrms->Scale(1./double(nAdded));

// Set the errors of the average and the individual profiles.
  //Calculate the bin factor
  TH1D *haRes = new TH1D(*BCrms);
  for( int ibin=1; ibin<=BCaverage->GetNbinsX(); ibin++)
  {
	  if(BCaverage->GetBinContent(ibin) > 0)
		  haRes->SetBinContent(ibin, sqrt(BCrms->GetBinContent(ibin) / BCaverage->GetBinContent(ibin)));
	  BCaverage->SetBinError(ibin, sqrt(BCrms->GetBinContent(ibin) / double(nAdded)));
  }

  BCaverage->Write("Layers_average");
  haRes->Write("ares");

  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
    std::stringstream hname; hname << "Layers_" << iEvt;
    histograms.GetObject(hname.str().c_str(), BeamCalEnergy);
    if(BeamCalEnergy)
    {
/*	  for( int ibin=1; ibin<=BeamCalEnergy->GetNbinsX(); ibin++)
		  BeamCalEnergy->SetBinError(ibin, /*haRes->GetBinContent(ibin)*/
//				  ares * sqrt(BeamCalEnergy->GetBinContent(ibin)));
      BeamCalEnergy->Scale(1./e0);
      BeamCalEnergy->Sumw2();
      BeamCalEnergy->Scale(e0);
	  BeamCalEnergy->Write(hname.str().c_str(), TObject::kOverwrite);
    }
  }

//  histograms.Write();


  return 0;
}
